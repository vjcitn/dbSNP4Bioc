## dbsnp_client.R -- query slim dbSNP Parquet over HTTP, or from a BiocFileCache copy
##
## Needs: DBI, duckdb, jsonlite, BiocFileCache, digest (digest only for checksum fallback).
##
##   options(dbsnp.base_url = "https://<pod>.osn.xsede.org/<bucket>/dbsnp")
##
##   ## small queries straight over HTTP (range requests; nothing downloaded to disk)
##   db <- dbsnp_open("b157", "GRCh38", chroms = "17", source = "remote")
##   rsid_to_pos(db, c("rs1042522", "rs28934578"))
##   pos_to_rsid(db, "17", 7676154, a1 = "G", a2 = "C", match = "either_order")
##   dbsnp_close(db)
##
##   ## bulk / offline: copy chr17 into the BiocFileCache once, then no network needed
##   dbsnp_fetch_chrom("17", "b157", "GRCh38")
##   db <- dbsnp_open("b157", "GRCh38", source = "cache")

.q <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
.chrom_label <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
.dbsnp_tables <- c("by_pos", "by_rsid")

.sha256 <- function(f) {
  if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE))
    return(unname(tools::sha256sum(f)))
  digest::digest(file = f, algo = "sha256", serialize = FALSE)
}

.dbsnp_base <- function(base_url) {
  if (is.null(base_url) || !nzchar(base_url))
    stop("set options(dbsnp.base_url = ...) or pass base_url")
  sub("/+$", "", base_url)
}
.dbsnp_url <- function(base_url, release, build, rel)
  paste(.dbsnp_base(base_url), release, build, rel, sep = "/")
.dbsnp_rname <- function(release, build, rel) sprintf("dbsnp/%s/%s/%s", release, build, rel)

## ---- BiocFileCache helpers -------------------------------------------------------
.bfc_find <- function(bfc, rname) {
  q <- BiocFileCache::bfcquery(bfc, rname, field = "rname", exact = TRUE)
  if (nrow(q)) q$rid[[1L]] else NA_character_
}
## cached path (no network if already present); download otherwise
.bfc_get <- function(bfc, rname, url, refresh = FALSE) {
  rid <- .bfc_find(bfc, rname)
  if (is.na(rid)) return(unname(BiocFileCache::bfcrpath(bfc, rnames = rname, fpath = url)))
  if (refresh) BiocFileCache::bfcdownload(bfc, rid, ask = FALSE)
  unname(BiocFileCache::bfcpath(bfc, rid))
}

## manifest: fetched and re-cached when online; read from the cache when offline
dbsnp_manifest <- function(release = "b157", build = "GRCh38",
                           base_url = getOption("dbsnp.base_url"),
                           cache = BiocFileCache::BiocFileCache(ask = FALSE),
                           online = TRUE) {
  rname <- .dbsnp_rname(release, build, "manifest.json")
  path <- NULL
  if (online) {
    url <- .dbsnp_url(base_url, release, build, "manifest.json")
    path <- tryCatch(.bfc_get(cache, rname, url, refresh = TRUE), error = function(e) NULL)
  }
  if (is.null(path)) {
    rid <- .bfc_find(cache, rname)
    if (is.na(rid)) return(NULL)
    path <- unname(BiocFileCache::bfcpath(cache, rid))
  }
  jsonlite::read_json(path, simplifyVector = FALSE)
}

## Copy one or more chromosomes (both tables) into the BiocFileCache and verify sha256.
dbsnp_fetch_chrom <- function(chrom, release = "b157", build = "GRCh38",
                              tables = .dbsnp_tables,
                              base_url = getOption("dbsnp.base_url"),
                              cache = BiocFileCache::BiocFileCache(ask = FALSE),
                              verify = TRUE) {
  chrom <- .chrom_label(chrom)
  man <- dbsnp_manifest(release, build, base_url, cache)
  if (is.null(man)) stop("cannot read manifest for ", release, "/", build)
  out <- list()
  for (ch in chrom) {
    if (is.null(man$files[[ch]])) stop("chromosome ", ch, " not in manifest")
    for (tb in tables) {
      ent <- man$files[[ch]][[tb]]
      rname <- .dbsnp_rname(release, build, ent$path)
      path <- .bfc_get(cache, rname, .dbsnp_url(base_url, release, build, ent$path))
      if (verify && !identical(.sha256(path), ent$sha256)) {
        BiocFileCache::bfcremove(cache, .bfc_find(cache, rname))
        stop("sha256 mismatch for ", ent$path, "; removed from cache, retry")
      }
      out[[paste(ch, tb)]] <- path
    }
  }
  invisible(out)
}

dbsnp_cached <- function(release = "b157", build = "GRCh38",
                         cache = BiocFileCache::BiocFileCache(ask = FALSE)) {
  info <- BiocFileCache::bfcinfo(cache)
  pre <- sprintf("dbsnp/%s/%s/", release, build)
  info <- info[startsWith(info$rname, pre) & grepl("\\.parquet$", info$rname), , drop = FALSE]
  if (!nrow(info)) return(data.frame(chrom = character(), table = character(), path = character()))
  rel <- substring(info$rname, nchar(pre) + 1L)
  data.frame(chrom = sub("\\.parquet$", "", basename(rel)), table = dirname(rel),
             path = unname(vapply(info$rid, function(r) BiocFileCache::bfcpath(cache, r), "")),
             stringsAsFactors = FALSE)
}

## ---- handle ----------------------------------------------------------------------
## source = "cache": only BiocFileCache copies (never touches the network)
## source = "remote": always HTTP;  "auto": cache when present, else HTTP
dbsnp_open <- function(release = "b157", build = "GRCh38", chroms = NULL,
                       source = c("auto", "cache", "remote"),
                       base_url = getOption("dbsnp.base_url"),
                       cache = BiocFileCache::BiocFileCache(ask = FALSE),
                       memory_limit = NULL, threads = NULL) {
  source <- match.arg(source)
  cached <- dbsnp_cached(release, build, cache)
  files <- cached[0, ]; files$is_url <- logical()
  if (source != "remote" && nrow(cached)) { files <- cbind(cached, is_url = FALSE) }
  if (source != "cache") {
    man <- dbsnp_manifest(release, build, base_url, cache, online = TRUE)
    if (!is.null(man)) {
      for (ch in names(man$files)) for (tb in .dbsnp_tables) {
        have <- source != "remote" && any(files$chrom == ch & files$table == tb)
        if (!have) files <- rbind(files, data.frame(
          chrom = ch, table = tb, stringsAsFactors = FALSE,
          path = .dbsnp_url(base_url, release, build, man$files[[ch]][[tb]]$path),
          is_url = TRUE))
      }
    } else if (source == "remote") stop("cannot reach manifest")
  }
  if (!is.null(chroms)) files <- files[files$chrom %in% .chrom_label(chroms), , drop = FALSE]
  if (!nrow(files)) stop("no dbSNP files for ", release, "/", build, " (source = ", source, ")")

  con <- DBI::dbConnect(duckdb::duckdb())
  if (!is.null(threads)) DBI::dbExecute(con, sprintf("SET threads = %d", as.integer(threads)))
  if (!is.null(memory_limit)) DBI::dbExecute(con, sprintf("SET memory_limit = %s", .q(memory_limit)))
  if (any(files$is_url)) {
    try(DBI::dbExecute(con, "INSTALL httpfs"), silent = TRUE)   # needs network once
    DBI::dbExecute(con, "LOAD httpfs")
  }
  structure(list(con = con, files = files, release = release, build = build),
            class = "dbsnp_handle")
}
dbsnp_close <- function(db) invisible(DBI::dbDisconnect(db$con, shutdown = TRUE))
print.dbsnp_handle <- function(x, ...) {
  cat("<dbsnp_handle>", x$release, x$build, "\n")
  print(table(chrom = x$files$chrom, source = ifelse(x$files$is_url, "remote", "cache")) / 2)
  invisible(x)
}

.flist <- function(db, tb, chroms = NULL) {
  f <- db$files[db$files$table == tb, , drop = FALSE]
  if (!is.null(chroms)) f <- f[f$chrom %in% chroms, , drop = FALSE]
  if (!nrow(f)) stop("no ", tb, " files for requested chromosomes")
  paste0("[", paste(.q(f$path), collapse = ", "), "]")
}

.finish <- function(df) {
  df$RefSNP_id <- sprintf("rs%.0f", as.numeric(df$rs))
  df$alleles <- paste(df$ref, df$alt, sep = "/")
  df$seqnames <- df$chrom
  df[, intersect(c("query_index", "seqnames", "pos", "RefSNP_id", "ref", "alt", "alleles", "n_alt"),
                 names(df)), drop = FALSE]
}

## ---- rsID -> position/alleles ----------------------------------------------------
.parse_rsid <- function(x) {
  x <- sub("^rs", "", trimws(as.character(x)), ignore.case = TRUE)
  suppressWarnings(as.numeric(ifelse(grepl("^[0-9]+$", x), x, NA)))   # double: no int32 overflow
}

.sql_rsid <- function(files, qids, in_list_max) {
  if (length(qids) <= in_list_max)
    return(sprintf(paste0("SELECT chrom, pos, rs, ref, alt, n_alt FROM read_parquet(%s) ",
                          "WHERE rs IN (%s) ORDER BY rs, pos"),
                   files, paste(sprintf("%.0f", qids), collapse = ",")))
  sprintf(paste0("SELECT p.chrom, p.pos, p.rs, p.ref, p.alt, p.n_alt FROM read_parquet(%s) p ",
                 "JOIN (SELECT DISTINCT CAST(rs AS UINTEGER) AS rs FROM q_rs_df) q ON p.rs = q.rs"),
          files)
}

rsid_to_pos <- function(db, rsids, chroms = NULL, in_list_max = 2000L) {
  qids <- unique(.parse_rsid(rsids)); qids <- sort(qids[!is.na(qids)])
  if (!length(qids)) return(.finish(data.frame(chrom = character(), pos = numeric(), rs = numeric(),
                                               ref = character(), alt = character(), n_alt = integer())))
  if (length(qids) > in_list_max) {
    duckdb::duckdb_register(db$con, "q_rs_df", data.frame(rs = qids))
    on.exit(duckdb::duckdb_unregister(db$con, "q_rs_df"), add = TRUE)
  }
  .finish(DBI::dbGetQuery(db$con, .sql_rsid(.flist(db, "by_rsid", chroms), qids, in_list_max)))
}

## ---- position (+ alleles) -> rsID -------------------------------------------------
## match = "any"          position only
##         "either_order" {ref,alt} == {a1,a2}
##         "strand_flip"  as either_order, also against complemented a1/a2
##                        (A/T and C/G SNPs are strand-ambiguous: check before trusting)
## Returns one row per (query, matching dbSNP row); query_index is the row in your input.
.comp <- function(x) chartr("ACGT", "TGCA", toupper(x))

.sql_pos <- function(files, match, posvals, in_list_max) {
  cond <- switch(match,
    any = "TRUE",
    either_order = "((p.ref = q.a1 AND p.alt = q.a2) OR (p.ref = q.a2 AND p.alt = q.a1))",
    strand_flip = paste0("((p.ref = q.a1 AND p.alt = q.a2) OR (p.ref = q.a2 AND p.alt = q.a1) OR ",
                         "(p.ref = q.c1 AND p.alt = q.c2) OR (p.ref = q.c2 AND p.alt = q.c1))"))
  prune <- if (length(posvals) <= in_list_max)
    sprintf(" AND p.pos IN (%s)", paste(sprintf("%.0f", posvals), collapse = ",")) else ""
  sprintf(paste0("SELECT q.qi AS query_index, p.chrom, p.pos, p.rs, p.ref, p.alt, p.n_alt ",
                 "FROM read_parquet(%s) p JOIN q_pos_df q ON p.chrom = q.chrom AND p.pos = q.pos ",
                 "WHERE %s%s ORDER BY q.qi, p.rs"), files, cond, prune)
}

pos_to_rsid <- function(db, chrom, pos, a1 = NULL, a2 = NULL,
                        match = c("any", "either_order", "strand_flip"), in_list_max = 2000L) {
  match <- match.arg(match)
  n <- length(pos)
  chrom <- rep_len(.chrom_label(chrom), n)
  q <- data.frame(qi = seq_len(n), chrom = chrom, pos = as.numeric(pos), stringsAsFactors = FALSE)
  if (match != "any") {
    if (is.null(a1) || is.null(a2)) stop("a1 and a2 are required for match = '", match, "'")
    q$a1 <- toupper(as.character(a1)); q$a2 <- toupper(as.character(a2))
    q$c1 <- .comp(q$a1); q$c2 <- .comp(q$a2)
  }
  q <- q[!is.na(q$pos), , drop = FALSE]
  duckdb::duckdb_register(db$con, "q_pos_df", q)
  on.exit(duckdb::duckdb_unregister(db$con, "q_pos_df"), add = TRUE)
  chroms <- unique(q$chrom)
  sql <- .sql_pos(.flist(db, "by_pos", chroms), match, unique(q$pos), in_list_max)
  .finish(DBI::dbGetQuery(db$con, sql))
}

.sql_range <- function(files, chrom, start, end, biallelic_only = FALSE)
  sprintf(paste0("SELECT chrom, pos, rs, ref, alt, n_alt FROM read_parquet(%s) ",
                 "WHERE chrom = %s AND pos BETWEEN %.0f AND %.0f%s ORDER BY pos, rs, alt"),
          files, .q(chrom), start, end, if (biallelic_only) " AND n_alt = 1" else "")

snps_in_range <- function(db, chrom, start, end, biallelic_only = FALSE) {
  chrom <- .chrom_label(chrom)
  stopifnot(length(chrom) == 1L, length(start) == 1L, length(end) == 1L, start <= end)
  .finish(DBI::dbGetQuery(db$con,
            .sql_range(.flist(db, "by_pos", chrom), chrom, start, end, biallelic_only)))
}
