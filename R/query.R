#' A tidyverse-compatible view onto the cached dbSNP Parquet files
#'
#' Returns a lazy `dplyr`/`dbplyr` table backed by a `duckdb` connection
#' reading the Parquet files that [dbsnp_build_chrom()] cached, so that
#' ordinary `dplyr::filter()`, `dplyr::select()`, `dplyr::mutate()`,
#' `dplyr::arrange()`, and `dplyr::collect()` calls work directly against
#' the cache. Nothing is read into memory until `collect()` (or another
#' action, such as `print()` or `pull()`) is called: duckdb pushes the
#' translated SQL predicate down into the Parquet scan.
#'
#' @param table Which sort order to read: `"by_pos"` (sorted by position;
#'   best for range queries) or `"by_rsid"` (sorted by rsID; best for rsID
#'   lookups).
#' @param chrom Optional chromosome, or vector of chromosomes, to restrict
#'   to (e.g. `"17"`). Defaults to every cached chromosome for the given
#'   `table`.
#' @param release dbSNP release label used when caching (e.g. `"b157"`).
#' @param build Genome build label used when caching (e.g. `"GRCh38"`).
#' @param cache A [BiocFileCache::BiocFileCache()] object; defaults to
#'   [dbsnp_cache()].
#' @param con An existing duckdb [DBI::dbConnect()] connection to reuse;
#'   defaults to opening a fresh in-process duckdb connection. Close it with
#'   [dbsnp_disconnect()] when done (or let it be garbage-collected).
#' @return A `tbl_duckdb_connection` / `tbl_sql` object: a lazy `dplyr`
#'   table with columns `chrom`, `pos`, `rs`, `ref`, `alt`, `n_alt`.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' chr17 <- dbsnp_tbl("by_pos", chrom = "17")
#' chr17 |>
#'   filter(pos >= 38077296, pos <= 38083854, n_alt == 1) |>
#'   select(pos, rs, ref, alt) |>
#'   collect()
#' dbsnp_disconnect(chr17)
#' }
#' @export
dbsnp_tbl <- function(table = c("by_pos", "by_rsid"), chrom = NULL,
                      release = "b157", build = "GRCh38",
                      cache = dbsnp_cache(), con = NULL) {
  table <- match.arg(table)
  files <- dbsnp_cached(release, build, cache)
  files <- files[files$table == table, , drop = FALSE]
  if (!is.null(chrom)) files <- files[files$chrom %in% .chrom_label(chrom), , drop = FALSE]
  if (!nrow(files))
    stop("no cached '", table, "' files for ", release, "/", build,
         if (!is.null(chrom)) paste0(", chrom ", paste(chrom, collapse = ", ")) else "",
         "; build with dbsnp_build_chrom() first", call. = FALSE)

  if (is.null(con)) con <- DBI::dbConnect(duckdb::duckdb())
  flist <- paste(vapply(files$path, .q, ""), collapse = ", ")
  sql <- sprintf("SELECT * FROM read_parquet([%s])", flist)
  dplyr::tbl(con, dplyr::sql(sql))
}

#' Close the connection behind a [dbsnp_tbl()] result
#'
#' @param tbl A lazy table returned by [dbsnp_tbl()] (or any `tbl_sql`).
#' @return `NULL`, invisibly.
#' @export
dbsnp_disconnect <- function(tbl) {
  DBI::dbDisconnect(dbplyr::remote_con(tbl), shutdown = TRUE)
  invisible(NULL)
}

#' Look up rsIDs by genomic position
#'
#' A `dplyr`-based convenience wrapper over [dbsnp_tbl()]: filters the
#' `by_pos` table to one chromosome and position range. For exact-position
#' lookups across many positions, prefer building your own
#' [dplyr::semi_join()] against [dbsnp_tbl()] directly.
#'
#' @inheritParams dbsnp_tbl
#' @param chrom A single chromosome (e.g. `"17"`).
#' @param start,end 1-based inclusive position range.
#' @param biallelic_only If `TRUE`, keep only records with a single ALT
#'   allele (`n_alt == 1`).
#' @return A tibble with columns `chrom`, `pos`, `rs`, `ref`, `alt`, `n_alt`.
#' @importFrom rlang .data
#' @export
dbsnp_in_range <- function(chrom, start, end, biallelic_only = FALSE,
                           release = "b157", build = "GRCh38",
                           cache = dbsnp_cache()) {
  stopifnot(length(chrom) == 1L, length(start) == 1L, length(end) == 1L, start <= end)
  chrom <- .chrom_label(chrom)
  tbl <- dbsnp_tbl("by_pos", chrom = chrom, release = release, build = build, cache = cache)
  on.exit(dbsnp_disconnect(tbl), add = TRUE)
  out <- dplyr::filter(tbl, .data$chrom == !!chrom, .data$pos >= !!start, .data$pos <= !!end)
  if (biallelic_only) out <- dplyr::filter(out, .data$n_alt == 1L)
  dplyr::collect(dplyr::arrange(out, .data$pos, .data$rs, .data$alt))
}

.parse_rsid <- function(x) {
  x <- sub("^rs", "", trimws(as.character(x)), ignore.case = TRUE)
  suppressWarnings(as.numeric(ifelse(grepl("^[0-9]+$", x), x, NA)))   # double: rs ids can exceed 2^31
}

#' Look up one or more rsIDs across every cached chromosome
#'
#' A `dplyr`-based convenience wrapper over [dbsnp_tbl()]: searches the
#' `by_rsid` table across *all* chromosomes currently cached for the given
#' `release`/`build` (not just one), since an rsID does not tell you in
#' advance which chromosome it lives on.
#'
#' @inheritParams dbsnp_tbl
#' @param rsid One or more rsIDs, with or without the `"rs"` prefix (e.g.
#'   `"rs1042522"` or `1042522`).
#' @return A tibble with columns `chrom`, `pos`, `rs`, `ref`, `alt`, `n_alt`,
#'   `RefSNP_id` (one row per matching dbSNP record; a multi-allelic rsID
#'   returns one row per ALT allele). rsIDs with no match in the cache are
#'   simply absent from the result.
#' @importFrom rlang .data
#' @export
dbsnp_by_rsid <- function(rsid, release = "b157", build = "GRCh38", cache = dbsnp_cache()) {
  qids <- unique(.parse_rsid(rsid))
  qids <- qids[!is.na(qids)]
  if (!length(qids)) stop("no valid rsIDs in 'rsid'", call. = FALSE)

  tbl <- dbsnp_tbl("by_rsid", release = release, build = build, cache = cache)
  on.exit(dbsnp_disconnect(tbl), add = TRUE)
  out <- dplyr::collect(dplyr::arrange(dplyr::filter(tbl, .data$rs %in% !!qids),
                                       .data$rs, .data$chrom, .data$pos))
  out$RefSNP_id <- sprintf("rs%.0f", out$rs)
  out
}
