# dbSNP4Bioc

Build, cache, and query dbSNP VCF data as Parquet.

`dbSNP4Bioc` converts a single chromosome of an NCBI dbSNP VCF (as
distributed: bgzipped and tabix-indexed) into a slim, sorted Parquet
representation, caches it locally with `BiocFileCache`, and exposes the
cached resource as a lazy `dplyr`/`dbplyr` table backed by `duckdb`, so
that rsID and genomic-position lookups can be written with ordinary
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html) /
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
pipelines.

## Setup

Point `DBSNP_GCF_PATH` at NCBI’s bgzipped, tabix-indexed dbSNP VCF (the
`.gz.tbi` index must sit next to it):

``` sh
export DBSNP_GCF_PATH=/path/to/GCF_000001405.40.gz
```

## Build and cache a chromosome

``` r

library(dbSNP4Bioc)

dbsnp_build_chrom("17", release = "b157", build = "GRCh38")
dbsnp_cached("b157", "GRCh38")   # see what is cached
```

This extracts chromosome 17 with `tabix`, keeps SNV records only (REF
and every ALT a single A/C/G/T; one row per ALT allele), and writes two
Parquet tables into the cache: `by_pos` (sorted by position) and
`by_rsid` (sorted by rsID).

## Query the cache with dplyr

``` r

library(dplyr)

chr17 <- dbsnp_tbl("by_pos", chrom = "17", release = "b157", build = "GRCh38")

chr17 |>
  filter(pos >= 38077296, pos <= 38083854, n_alt == 1) |>
  select(pos, rs, ref, alt) |>
  collect()

dbsnp_disconnect(chr17)
```

A convenience wrapper is also provided for position-range lookups:

``` r

dbsnp_in_range("17", 38077296, 38083854, biallelic_only = TRUE,
               release = "b157", build = "GRCh38")
```

## Reference

See `inst/scripts/` for the original standalone scripts this package was
built from, including notes on fetching the VCF from NCBI and an
end-to-end example querying SNPs around a gene with
`TxDb.Hsapiens.UCSC.hg38.knownGene`.
