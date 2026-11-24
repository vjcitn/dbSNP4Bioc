## Convert one chromosome of an NCBI dbSNP VCF to slim, sorted Parquet, and
## cache it with BiocFileCache.
##
## Rows kept: SNVs only (REF and every ALT a single A/C/G/T), one row per ALT
## allele.
##   chrom  VARCHAR   "17" (no "chr")
##   pos    UINTEGER  1-based
##   rs     UINTEGER  numeric part of the rsID (rsIDs can exceed 2^31: never int32)
##   ref    VARCHAR   VCF REF
##   alt    VARCHAR   one ALT allele
##   n_alt  UTINYINT  number of ALT alleles on the record (bi-allelic <=> n_alt == 1)
##
## by_pos is sorted (pos, rs, alt); by_rsid is sorted (rs, pos, alt), with
## small row groups so Parquet zone maps can prune queries.

.dbsnp_converter_version <- "0.1.0"

#' Default BiocFileCache for dbSNP4Bioc
#'
#' A thin wrapper around [BiocFileCache::BiocFileCache()] used as the default
#' cache for every function in this package, so callers who want a
#' non-default cache location only have to pass `cache` once.
#'
#' @param cache_dir Optional cache directory; defaults to the package's
#'   BiocFileCache default location.
#' @return A `BiocFileCache` object.
#' @export
dbsnp_cache <- function(cache_dir = NULL) {
  if (is.null(cache_dir)) BiocFileCache::BiocFileCache(ask = FALSE)
  else BiocFileCache::BiocFileCache(cache_dir, ask = FALSE)
}

#' Locate the dbSNP VCF named by `DBSNP_GCF_PATH`
#'
#' NCBI distributes dbSNP as a single bgzipped, tabix-indexed VCF per
#' assembly (e.g. `GCF_000001405.40.gz` with `GCF_000001405.40.gz.tbi` next
#' to it). This package expects the environment variable `DBSNP_GCF_PATH` to
#' point at that `.gz` file.
#'
#' @param path Path to use instead of the environment variable (mainly for
#'   testing).
#' @return The validated path, as a string.
#' @export
dbsnp_gcf_path <- function(path = Sys.getenv("DBSNP_GCF_PATH")) {
  if (!nzchar(path))
    stop("DBSNP_GCF_PATH is not set; point it at the bgzipped, tabix-indexed ",
         "dbSNP VCF (e.g. GCF_000001405.40.gz)", call. = FALSE)
  if (!file.exists(path)) stop("DBSNP_GCF_PATH does not exist: ", path, call. = FALSE)
  if (!file.exists(paste0(path, ".tbi")))
    stop("no tabix index found next to DBSNP_GCF_PATH: ", path, ".tbi", call. = FALSE)
  path
}

## NCBI VCFs use RefSeq accessions (NC_000017.11 on GRCh38, NC_000017.10 on
## GRCh37); find the one actually in the tabix index instead of hard-coding
## versions.

#' Resolve a chromosome to its RefSeq contig name in a dbSNP VCF
#'
#' @param vcf Path to the bgzipped, tabix-indexed dbSNP VCF.
#' @param chrom Chromosome label: `"1"`-`"22"`, `"X"`, or `"Y"` (a leading
#'   `"chr"` is stripped).
#' @param tabix Name or path of the `tabix` executable.
#' @return The matching `NC_######.#` contig name, as a string.
#' @export
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

.update_manifest <- function(mf, chrom, files, release, build, vcf) {
  m <- if (file.exists(mf)) jsonlite::read_json(mf, simplifyVector = FALSE) else list()
  m <- utils::modifyList(m, list(
    schema_version = 1L, release = release, build = build,
    source = list(file = basename(vcf)),
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

#' Convert one dbSNP chromosome to Parquet and cache it
#'
#' Extracts a single chromosome from the dbSNP VCF named by
#' [dbsnp_gcf_path()] (or the `vcf` argument), keeps SNV records only, writes
#' two sorted Parquet tables (`by_pos`, `by_rsid`), and stores them plus an
#' updated `manifest.json` in a [BiocFileCache::BiocFileCache()] cache. Later
#' calls to [dbsnp_tbl()] read straight from this cache.
#'
#' @param chrom Chromosome to build: `"1"`-`"22"`, `"X"`, or `"Y"`.
#' @param release dbSNP release label, used only to namespace the cache
#'   (e.g. `"b157"`).
#' @param build Genome build label, used only to namespace the cache (e.g.
#'   `"GRCh38"`).
#' @param vcf Path to the bgzipped, tabix-indexed dbSNP VCF. Defaults to
#'   [dbsnp_gcf_path()], i.e. the `DBSNP_GCF_PATH` environment variable.
#' @param cache A [BiocFileCache::BiocFileCache()] object; defaults to
#'   [dbsnp_cache()].
#' @param tabix Name or path of the `tabix` executable.
#' @param rg_pos,rg_rsid Parquet row-group sizes for the `by_pos` and
#'   `by_rsid` tables; smaller row groups let zone maps prune more.
#' @param threads,memory_limit Optional duckdb `SET threads` /
#'   `SET memory_limit` values.
#' @param tmpdir Scratch directory for the intermediate TSV and Parquet
#'   files before they are copied into the cache.
#' @return Invisibly, a data frame with columns `chrom`, `table`, `path`
#'   describing the newly cached files (as returned by [dbsnp_cached()]).
#' @export
dbsnp_build_chrom <- function(chrom, release = "b157", build = "GRCh38",
                              vcf = dbsnp_gcf_path(), cache = dbsnp_cache(),
                              tabix = "tabix", rg_pos = 122880L, rg_rsid = 16384L,
                              threads = NULL, memory_limit = NULL,
                              tmpdir = tempfile("dbsnp4bioc_build_")) {
  chrom <- .chrom_label(chrom)
  contig <- dbsnp_contig(vcf, chrom, tabix)
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tsv <- file.path(tmpdir, sprintf("dbsnp_%s_%s_chr%s.tsv", release, build, chrom))
  cmd <- sprintf("%s %s %s | cut -f1-5 > %s",
                 shQuote(tabix), shQuote(vcf), shQuote(contig), shQuote(tsv))
  if (system2("bash", c("-o", "pipefail", "-c", shQuote(cmd))) != 0L)
    stop("tabix extraction failed")
  if (!isTRUE(file.size(tsv) > 0)) stop("extraction produced no records for ", contig)

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  if (!is.null(threads)) DBI::dbExecute(con, sprintf("SET threads = %d", as.integer(threads)))
  if (!is.null(memory_limit)) DBI::dbExecute(con, sprintf("SET memory_limit = %s", .q(memory_limit)))
  DBI::dbExecute(con, sprintf("SET temp_directory = %s", .q(file.path(tmpdir, "duckdb_tmp"))))
  DBI::dbExecute(con, .sql_snv_table(tsv, chrom))

  specs <- list(by_pos  = list(order = "pos, rs, alt", rg = rg_pos),
                by_rsid = list(order = "rs, pos, alt", rg = rg_rsid))
  files <- lapply(names(specs), function(tb) {
    rel <- file.path(tb, paste0(chrom, ".parquet"))
    f <- file.path(tmpdir, rel)
    dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
    DBI::dbExecute(con, .sql_copy(specs[[tb]]$order, f, specs[[tb]]$rg))
    n <- DBI::dbGetQuery(con, sprintf("SELECT count(*) AS n FROM read_parquet(%s)", .q(f)))$n
    cached <- .bfc_put(cache, .dbsnp_rname(release, build, rel), f)
    list(path = rel, bytes = file.size(f), sha256 = .sha256(f), rows = as.numeric(n),
         row_group_size = specs[[tb]]$rg, cached_path = cached)
  })
  names(files) <- names(specs)

  mf <- file.path(tmpdir, "manifest.json")
  mf_rname <- .dbsnp_rname(release, build, "manifest.json")
  mf_rid <- .bfc_find(cache, mf_rname)
  if (!is.na(mf_rid)) file.copy(BiocFileCache::bfcpath(cache, mf_rid), mf, overwrite = TRUE)
  .update_manifest(mf, chrom, lapply(files, function(x) x[c("path", "bytes", "sha256", "rows", "row_group_size")]),
                    release, build, vcf)
  .bfc_put(cache, mf_rname, mf)

  message("cached ", files$by_pos$rows, " rows for chr", chrom, " (", release, "/", build, ")")
  have <- dbsnp_cached(release, build, cache)
  invisible(have[have$chrom == chrom, , drop = FALSE])
}

#' List what is cached for a dbSNP release/build
#'
#' @param release dbSNP release label (e.g. `"b157"`).
#' @param build Genome build label (e.g. `"GRCh38"`).
#' @param cache A [BiocFileCache::BiocFileCache()] object; defaults to
#'   [dbsnp_cache()].
#' @return A data frame with columns `chrom`, `table` (`"by_pos"` or
#'   `"by_rsid"`), and `path` (the cached file's local path).
#' @export
dbsnp_cached <- function(release = "b157", build = "GRCh38", cache = dbsnp_cache()) {
  info <- BiocFileCache::bfcinfo(cache)
  pre <- sprintf("dbsnp4bioc/%s/%s/", release, build)
  info <- info[startsWith(info$rname, pre) & grepl("\\.parquet$", info$rname), , drop = FALSE]
  if (!nrow(info))
    return(data.frame(chrom = character(), table = character(), path = character()))
  rel <- substring(info$rname, nchar(pre) + 1L)
  data.frame(chrom = sub("\\.parquet$", "", basename(rel)), table = dirname(rel),
             path = unname(vapply(info$rid, function(r) BiocFileCache::bfcpath(cache, r), "")),
             stringsAsFactors = FALSE)
}
