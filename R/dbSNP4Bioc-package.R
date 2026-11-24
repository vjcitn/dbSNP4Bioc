#' dbSNP4Bioc: build, cache, and query dbSNP VCF data as Parquet
#'
#' Converts a single chromosome of an NCBI dbSNP VCF (as distributed,
#' bgzipped and tabix-indexed) into a slim, sorted Parquet representation
#' with [dbsnp_build_chrom()], caches it with
#' [BiocFileCache::BiocFileCache()], and exposes the cached resource as a
#' lazy `dplyr`/`dbplyr` table with [dbsnp_tbl()] so that rsID and
#' genomic-position lookups can be expressed with ordinary
#' `dplyr::filter()` / `dplyr::select()` pipelines.
#'
#' @section Workflow:
#' \enumerate{
#'   \item Set the environment variable `DBSNP_GCF_PATH` to the bgzipped,
#'     tabix-indexed dbSNP VCF (e.g. `GCF_000001405.40.gz`, with
#'     `GCF_000001405.40.gz.tbi` next to it).
#'   \item Call [dbsnp_build_chrom()] once per chromosome of interest; this
#'     populates the local BiocFileCache.
#'   \item Call [dbsnp_tbl()] (or the convenience wrapper
#'     [dbsnp_in_range()]) to query the cache with `dplyr` verbs.
#' }
#'
#' @keywords internal
"_PACKAGE"
