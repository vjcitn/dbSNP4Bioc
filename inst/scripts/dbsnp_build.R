## dbsnp_build.R -- one chromosome of an NCBI dbSNP VCF -> slim, sorted Parquet
##
## Needs on PATH: bash, tabix, cut.  R packages: DBI, duckdb, jsonlite, digest.
## The VCF must have its .tbi next to it (NCBI ships one).
##
##   source("dbsnp_build.R")
##   dbsnp_chrom_to_parquet("GCF_000001405.40.gz", chrom = "17",
##                          outdir = "dbsnp_out", release = "b157", build = "GRCh38",
##                          vcf_md5 = "<md5 from NCBI .md5 file>")
##   # -> dbsnp_out/b157/GRCh38/{manifest.json, by_pos/17.parquet, by_rsid/17.parquet}
##
## Rows: SNVs only (REF and every ALT a single A/C/G/T), one row per ALT allele.
##   chrom  VARCHAR   "17" (no "chr")
##   pos    UINTEGER  1-based
##   rs     UINTEGER  numeric part of the rsID (rsIDs can exceed 2^31: never int32)
##   ref    VARCHAR   VCF REF
##   alt    VARCHAR   one ALT allele
##   n_alt  UTINYINT  number of ALT alleles on the record (bi-allelic <=> n_alt == 1)
## by_pos is sorted (pos, rs, alt); by_rsid is sorted (rs, pos, alt) with small row
## groups so zone maps can prune small remote queries.

.dbsnp_converter_version <- "0.1.0"

.q <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
.chrom_label <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

.sha256 <- function(f) {
  if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE))
    return(unname(tools::sha256sum(f)))
  digest::digest(file = f, algo = "sha256", serialize = FALSE)
}

.nc_regex <- function(chrom) {
  chrom <- toupper(.chrom_label(chrom))
  n <- if (chrom == "X") 23L else if (chrom == "Y") 24L else suppressWarnings(as.integer(chrom))
  if (is.na(n) || n < 1L || n > 24L) stop("chrom must be 1-22, X or Y")
  sprintf("^NC_%06d\\.", n)
}

## NCBI VCFs use RefSeq accessions (NC_000017.11 on GRCh38, NC_000017.10 on GRCh37);
## find the one in the tabix index instead of hard-coding versions.
dbsnp_contig <- function(vcf, chrom, tabix = "tabix") {
  contigs <- system2(tabix, c("-l", shQuote(vcf)), stdout = TRUE)
  hit <- grep(.nc_regex(chrom), contigs, value = TRUE)
  if (length(hit) != 1L)
    stop("expected exactly one NC_ contig for chrom ", chrom, "; found: ",
         paste(hit, collapse = ", "))
  hit
}

.sql_snv_table <- function(tsv, chrom) {
  sprintf("
CREATE TABLE snv AS
SELECT chrom, pos, rs, ref, alt, n_alt FROM (
  SELECT %s::VARCHAR AS chrom, pos,
         TRY_CAST(substr(id, 3) AS UINTEGER) AS rs,
         ref,
         unnest(string_split(alts, ',')) AS alt,
         CAST(len(string_split(alts, ',')) AS UTINYINT) AS n_alt
  FROM (
    SELECT pos, id, ref, alt AS alts
    FROM read_csv(%s, delim = '\\t', header = false, quote = '', escape = '',
                  auto_detect = false,
                  columns = {'chrom': 'VARCHAR', 'pos': 'UINTEGER', 'id': 'VARCHAR',
                             'ref': 'VARCHAR', 'alt': 'VARCHAR'})
    WHERE starts_with(id, 'rs')
      AND ref IN ('A','C','G','T')
      AND regexp_full_match(alt, '[ACGT](,[ACGT])*')
  )
) WHERE rs IS NOT NULL", .q(chrom), .q(tsv))
}

.sql_copy <- function(order_by, path, row_group_size)
  sprintf("COPY (SELECT * FROM snv ORDER BY %s) TO %s (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE %d)",
          order_by, .q(path), as.integer(row_group_size))

.update_manifest <- function(base, chrom, files, release, build, vcf, vcf_md5) {
  mf <- file.path(base, "manifest.json")
  m <- if (file.exists(mf)) jsonlite::read_json(mf, simplifyVector = FALSE) else list()
  src <- Filter(Negate(is.null), list(file = basename(vcf),
                                      md5 = if (is.na(vcf_md5)) NULL else vcf_md5))
  m <- utils::modifyList(m, list(
    schema_version = 1L, release = release, build = build, source = src,
    converter = list(version = .dbsnp_converter_version,
      filter = "SNV only: REF and every ALT a single A/C/G/T; one row per ALT allele; n_alt = ALT alleles on the VCF record"),
    columns = list(chrom = "VARCHAR", pos = "UINTEGER (1-based)", rs = "UINTEGER (numeric part of rsID)",
                   ref = "VARCHAR", alt = "VARCHAR", n_alt = "UTINYINT"),
    sort = list(by_pos = "pos, rs, alt", by_rsid = "rs, pos, alt"),
    created = format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")))
  if (is.null(m$files)) m$files <- list()
  m$files[[chrom]] <- files
  jsonlite::write_json(m, mf, auto_unbox = TRUE, pretty = TRUE)
  invisible(m)
}

dbsnp_chrom_to_parquet <- function(vcf, chrom, outdir, release = "b157", build = "GRCh38",
                                   vcf_md5 = NA_character_, tmpdir = tempdir(),
                                   tabix = "tabix", rg_pos = 122880L, rg_rsid = 16384L,
                                   threads = NULL, memory_limit = NULL, keep_tsv = FALSE) {
  stopifnot(file.exists(vcf))
  chrom <- .chrom_label(chrom)
  contig <- dbsnp_contig(vcf, chrom, tabix)
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
  tsv <- file.path(tmpdir, sprintf("dbsnp_%s_%s_chr%s.tsv", release, build, chrom))

  ## 1. region extract, first five VCF columns only (INFO is what makes the VCF huge)
  cmd <- sprintf("%s %s %s | cut -f1-5 > %s",
                 shQuote(tabix), shQuote(vcf), shQuote(contig), shQuote(tsv))
  if (system2("bash", c("-o", "pipefail", "-c", shQuote(cmd))) != 0L)
    stop("tabix extraction failed")
  if (!isTRUE(file.size(tsv) > 0)) stop("extraction produced no records for ", contig)
  if (!keep_tsv) on.exit(unlink(tsv), add = TRUE)

  ## 2. filter + explode + write both sort orders
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  if (!is.null(threads)) DBI::dbExecute(con, sprintf("SET threads = %d", as.integer(threads)))
  if (!is.null(memory_limit)) DBI::dbExecute(con, sprintf("SET memory_limit = %s", .q(memory_limit)))
  DBI::dbExecute(con, sprintf("SET temp_directory = %s", .q(file.path(tmpdir, "duckdb_tmp"))))
  DBI::dbExecute(con, .sql_snv_table(tsv, chrom))

  base <- file.path(outdir, release, build)
  specs <- list(by_pos  = list(order = "pos, rs, alt", rg = rg_pos),
                by_rsid = list(order = "rs, pos, alt", rg = rg_rsid))
  files <- lapply(names(specs), function(tb) {
    dir.create(file.path(base, tb), recursive = TRUE, showWarnings = FALSE)
    rel <- file.path(tb, paste0(chrom, ".parquet"))
    f <- file.path(base, rel)
    DBI::dbExecute(con, .sql_copy(specs[[tb]]$order, f, specs[[tb]]$rg))
    n <- DBI::dbGetQuery(con, sprintf("SELECT count(*) AS n FROM read_parquet(%s)", .q(f)))$n
    list(path = rel, bytes = file.size(f), sha256 = .sha256(f), rows = as.numeric(n),
         row_group_size = specs[[tb]]$rg)
  })
  names(files) <- names(specs)
  .update_manifest(base, chrom, files, release, build, vcf, vcf_md5)
  message("wrote ", files$by_pos$rows, " rows per table to ", base)
  invisible(files)
}

## ---- publishing to OSN (S3 API; bucket must allow anonymous GET) -------------------
## aws s3 sync dbsnp_out/ s3://<bucket>/dbsnp/ --endpoint-url https://<pod>.osn.xsede.org
##   (or rclone with an s3 remote of provider "Other").  Upload manifest.json LAST.
## Check:  curl -sI -H 'Range: bytes=0-99' https://<pod>.osn.xsede.org/<bucket>/dbsnp/b157/GRCh38/by_pos/17.parquet
##   -> must say 206 Partial Content; that is what lets DuckDB read small remote queries.
