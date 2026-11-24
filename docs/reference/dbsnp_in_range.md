# Look up rsIDs by genomic position

A `dplyr`-based convenience wrapper over
[`dbsnp_tbl()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_tbl.md):
filters the `by_pos` table to one chromosome and position range. For
exact-position lookups across many positions, prefer building your own
[`dplyr::semi_join()`](https://dplyr.tidyverse.org/reference/filter-joins.html)
against
[`dbsnp_tbl()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_tbl.md)
directly.

## Usage

``` r
dbsnp_in_range(
  chrom,
  start,
  end,
  biallelic_only = FALSE,
  release = "b157",
  build = "GRCh38",
  cache = dbsnp_cache()
)
```

## Arguments

- chrom:

  A single chromosome (e.g. `"17"`).

- start, end:

  1-based inclusive position range.

- biallelic_only:

  If `TRUE`, keep only records with a single ALT allele (`n_alt == 1`).

- release:

  dbSNP release label used when caching (e.g. `"b157"`).

- build:

  Genome build label used when caching (e.g. `"GRCh38"`).

- cache:

  A
  [`BiocFileCache::BiocFileCache()`](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
  object; defaults to
  [`dbsnp_cache()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_cache.md).

## Value

A tibble with columns `chrom`, `pos`, `rs`, `ref`, `alt`, `n_alt`.
