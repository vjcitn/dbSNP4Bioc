## Internal helpers shared by the build and query layers. Not exported.

.q <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")

.chrom_label <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

.dbsnp_tables <- c("by_pos", "by_rsid")

.dbsnp_rname <- function(release, build, rel)
  sprintf("dbsnp4bioc/%s/%s/%s", release, build, rel)

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

## BiocFileCache: find an existing cache entry by rname, NA if absent
.bfc_find <- function(cache, rname) {
  hit <- BiocFileCache::bfcquery(cache, rname, field = "rname", exact = TRUE)
  if (nrow(hit)) hit$rid[[1L]] else NA_character_
}

## add (or, if rname already present, overwrite in place) a local file; returns cached path
.bfc_put <- function(cache, rname, path) {
  rid <- .bfc_find(cache, rname)
  if (is.na(rid))
    return(unname(BiocFileCache::bfcadd(cache, rname = rname, fpath = path,
                                        rtype = "local", action = "copy")))
  BiocFileCache::bfcupdate(cache, rid, fpath = path, ask = FALSE)
  unname(BiocFileCache::bfcpath(cache, rid))
}
