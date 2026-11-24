# Locate the dbSNP VCF named by `DBSNP_GCF_PATH`

NCBI distributes dbSNP as a single bgzipped, tabix-indexed VCF per
assembly (e.g. `GCF_000001405.40.gz` with `GCF_000001405.40.gz.tbi` next
to it). This package expects the environment variable `DBSNP_GCF_PATH`
to point at that `.gz` file.

## Usage

``` r
dbsnp_gcf_path(path = Sys.getenv("DBSNP_GCF_PATH"))
```

## Arguments

- path:

  Path to use instead of the environment variable (mainly for testing).

## Value

The validated path, as a string.
