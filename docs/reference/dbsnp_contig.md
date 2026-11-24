# Resolve a chromosome to its RefSeq contig name in a dbSNP VCF

Resolve a chromosome to its RefSeq contig name in a dbSNP VCF

## Usage

``` r
dbsnp_contig(vcf, chrom, tabix = "tabix")
```

## Arguments

- vcf:

  Path to the bgzipped, tabix-indexed dbSNP VCF.

- chrom:

  Chromosome label: `"1"`-`"22"`, `"X"`, or `"Y"` (a leading `"chr"` is
  stripped).

- tabix:

  Name or path of the `tabix` executable.

## Value

The matching `NC_######.#` contig name, as a string.
