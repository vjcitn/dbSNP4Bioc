# A tidyverse-compatible view onto the cached dbSNP Parquet files

Returns a lazy `dplyr`/`dbplyr` table backed by a `duckdb` connection
reading the Parquet files that
[`dbsnp_build_chrom()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_build_chrom.md)
cached, so that ordinary
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html),
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html),
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html),
[`dplyr::arrange()`](https://dplyr.tidyverse.org/reference/arrange.html),
and
[`dplyr::collect()`](https://dplyr.tidyverse.org/reference/compute.html)
calls work directly against the cache. Nothing is read into memory until
[`collect()`](https://dplyr.tidyverse.org/reference/compute.html) (or
another action, such as [`print()`](https://rdrr.io/r/base/print.html)
or [`pull()`](https://dplyr.tidyverse.org/reference/pull.html)) is
called: duckdb pushes the translated SQL predicate down into the Parquet
scan.

## Usage

``` r
dbsnp_tbl(
  table = c("by_pos", "by_rsid"),
  chrom = NULL,
  release = "b157",
  build = "GRCh38",
  cache = dbsnp_cache(),
  con = NULL
)
```

## Arguments

- table:

  Which sort order to read: `"by_pos"` (sorted by position; best for
  range queries) or `"by_rsid"` (sorted by rsID; best for rsID lookups).

- chrom:

  Optional chromosome, or vector of chromosomes, to restrict to (e.g.
  `"17"`). Defaults to every cached chromosome for the given `table`.

- release:

  dbSNP release label used when caching (e.g. `"b157"`).

- build:

  Genome build label used when caching (e.g. `"GRCh38"`).

- cache:

  A
  [`BiocFileCache::BiocFileCache()`](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
  object; defaults to
  [`dbsnp_cache()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_cache.md).

- con:

  An existing duckdb
  [`DBI::dbConnect()`](https://dbi.r-dbi.org/reference/dbConnect.html)
  connection to reuse; defaults to opening a fresh in-process duckdb
  connection. Close it with
  [`dbsnp_disconnect()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_disconnect.md)
  when done (or let it be garbage-collected).

## Value

A `tbl_duckdb_connection` / `tbl_sql` object: a lazy `dplyr` table with
columns `chrom`, `pos`, `rs`, `ref`, `alt`, `n_alt`.

## Examples

``` r
if (FALSE) { # \dontrun{
library(dplyr)
chr17 <- dbsnp_tbl("by_pos", chrom = "17")
chr17 |>
  filter(pos >= 38077296, pos <= 38083854, n_alt == 1) |>
  select(pos, rs, ref, alt) |>
  collect()
dbsnp_disconnect(chr17)
} # }
```
