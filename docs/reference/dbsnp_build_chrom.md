# Convert one dbSNP chromosome to Parquet and cache it

Extracts a single chromosome from the dbSNP VCF named by
[`dbsnp_gcf_path()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_gcf_path.md)
(or the `vcf` argument), keeps SNV records only, writes two sorted
Parquet tables (`by_pos`, `by_rsid`), and stores them plus an updated
`manifest.json` in a
[`BiocFileCache::BiocFileCache()`](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
cache. Later calls to
[`dbsnp_tbl()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_tbl.md)
read straight from this cache.

## Usage

``` r
dbsnp_build_chrom(
  chrom,
  release = "b157",
  build = "GRCh38",
  vcf = dbsnp_gcf_path(),
  cache = dbsnp_cache(),
  tabix = "tabix",
  rg_pos = 122880L,
  rg_rsid = 16384L,
  threads = NULL,
  memory_limit = NULL,
  tmpdir = tempfile("dbsnp4bioc_build_")
)
```

## Arguments

- chrom:

  Chromosome to build: `"1"`-`"22"`, `"X"`, or `"Y"`.

- release:

  dbSNP release label, used only to namespace the cache (e.g. `"b157"`).

- build:

  Genome build label, used only to namespace the cache (e.g.
  `"GRCh38"`).

- vcf:

  Path to the bgzipped, tabix-indexed dbSNP VCF. Defaults to
  [`dbsnp_gcf_path()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_gcf_path.md),
  i.e. the `DBSNP_GCF_PATH` environment variable.

- cache:

  A
  [`BiocFileCache::BiocFileCache()`](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
  object; defaults to
  [`dbsnp_cache()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_cache.md).

- tabix:

  Name or path of the `tabix` executable.

- rg_pos, rg_rsid:

  Parquet row-group sizes for the `by_pos` and `by_rsid` tables; smaller
  row groups let zone maps prune more.

- threads, memory_limit:

  Optional duckdb `SET threads` / `SET memory_limit` values.

- tmpdir:

  Scratch directory for the intermediate TSV and Parquet files before
  they are copied into the cache.

## Value

Invisibly, a data frame with columns `chrom`, `table`, `path` describing
the newly cached files (as returned by
[`dbsnp_cached()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_cached.md)).
