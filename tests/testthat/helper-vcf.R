## Build a tiny bgzipped, tabix-indexed VCF for chr17, for tests only.
make_test_vcf <- function(dir) {
  vcf <- file.path(dir, "test.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=NC_000017.11,length=83257441>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "NC_000017.11\t100\trs1\tA\tG\t.\t.\t.",
    "NC_000017.11\t200\trs2\tC\tT,G\t.\t.\t.",
    "NC_000017.11\t300\trs3\tA\tAT\t.\t.\t.",   # indel: dropped
    "NC_000017.11\t400\trs4\tN\tA\t.\t.\t.",    # non-ACGT REF: dropped
    "NC_000017.11\t500\trs5\tG\tC\t.\t.\t."
  ), vcf)
  gz <- paste0(vcf, ".gz")
  system2("bgzip", c("-f", shQuote(vcf)))
  system2("tabix", c("-p", "vcf", shQuote(gz)))
  gz
}
