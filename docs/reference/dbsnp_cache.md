# Default BiocFileCache for dbSNP4Bioc

A thin wrapper around
[`BiocFileCache::BiocFileCache()`](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
used as the default cache for every function in this package, so callers
who want a non-default cache location only have to pass `cache` once.

## Usage

``` r
dbsnp_cache(cache_dir = NULL)
```

## Arguments

- cache_dir:

  Optional cache directory; defaults to the package's BiocFileCache
  default location.

## Value

A `BiocFileCache` object.
