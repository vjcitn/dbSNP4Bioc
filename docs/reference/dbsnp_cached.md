# List what is cached for a dbSNP release/build

List what is cached for a dbSNP release/build

## Usage

``` r
dbsnp_cached(release = "b157", build = "GRCh38", cache = dbsnp_cache())
```

## Arguments

- release:

  dbSNP release label (e.g. `"b157"`).

- build:

  Genome build label (e.g. `"GRCh38"`).

- cache:

  A
  [`BiocFileCache::BiocFileCache()`](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
  object; defaults to
  [`dbsnp_cache()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_cache.md).

## Value

A data frame with columns `chrom`, `table` (`"by_pos"` or `"by_rsid"`),
and `path` (the cached file's local path).
