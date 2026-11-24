# dbSNP4Bioc: build, cache, and query dbSNP VCF data as Parquet

Converts a single chromosome of an NCBI dbSNP VCF (as distributed,
bgzipped and tabix-indexed) into a slim, sorted Parquet representation
with
[`dbsnp_build_chrom()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_build_chrom.md),
caches it with
[`BiocFileCache::BiocFileCache()`](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html),
and exposes the cached resource as a lazy `dplyr`/`dbplyr` table with
[`dbsnp_tbl()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_tbl.md)
so that rsID and genomic-position lookups can be expressed with ordinary
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html) /
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
pipelines.

## Workflow

1.  Set the environment variable `DBSNP_GCF_PATH` to the bgzipped,
    tabix-indexed dbSNP VCF (e.g. `GCF_000001405.40.gz`, with
    `GCF_000001405.40.gz.tbi` next to it).

2.  Call
    [`dbsnp_build_chrom()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_build_chrom.md)
    once per chromosome of interest; this populates the local
    BiocFileCache.

3.  Call
    [`dbsnp_tbl()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_tbl.md)
    (or the convenience wrapper
    [`dbsnp_in_range()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_in_range.md))
    to query the cache with `dplyr` verbs.

## See also

Useful links:

- <https://github.com/vjcitn/dbSNP4Bioc>

- Report bugs at <https://github.com/vjcitn/dbSNP4Bioc/issues>

## Author

**Maintainer**: Vince Carey <stvjc@channing.harvard.edu>

Authors:

- Vince Carey <stvjc@channing.harvard.edu>
