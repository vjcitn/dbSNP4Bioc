skip_if_not(nzchar(Sys.which("tabix")) && nzchar(Sys.which("bgzip")),
            "tabix/bgzip not available")

test_that("build, cache, and tidyverse query round-trip", {
  tmp <- withr::local_tempdir()
  vcf <- make_test_vcf(tmp)
  cache <- dbsnp_cache(withr::local_tempdir())

  expect_identical(dbsnp_contig(vcf, "chr17"), "NC_000017.11")

  res <- dbsnp_build_chrom("17", release = "test", build = "GRCh38",
                           vcf = vcf, cache = cache)
  expect_setequal(res$table, c("by_pos", "by_rsid"))

  cached <- dbsnp_cached("test", "GRCh38", cache)
  expect_setequal(cached$table, c("by_pos", "by_rsid"))
  expect_true(all(file.exists(cached$path)))

  by_pos <- dplyr::collect(dbsnp_tbl("by_pos", chrom = "17", release = "test",
                                     build = "GRCh38", cache = cache))
  ## rs3 (indel) and rs4 (N ref) are dropped; rs2 explodes into two ALT rows
  expect_equal(nrow(by_pos), 4L)
  expect_setequal(by_pos$rs, c(1, 2, 2, 5))
  expect_false(3 %in% by_pos$rs)
  expect_false(4 %in% by_pos$rs)

  rng <- dbsnp_in_range("17", 1, 250, release = "test", build = "GRCh38", cache = cache)
  expect_equal(sort(rng$pos), c(100, 200, 200))

  rng_bi <- dbsnp_in_range("17", 1, 1000, biallelic_only = TRUE,
                           release = "test", build = "GRCh38", cache = cache)
  expect_true(all(rng_bi$n_alt == 1L))
  expect_setequal(rng_bi$pos, c(100, 500))
})

test_that("dbsnp_gcf_path validates the environment variable", {
  withr::local_envvar(DBSNP_GCF_PATH = "")
  expect_error(dbsnp_gcf_path(), "DBSNP_GCF_PATH is not set")

  tmp <- withr::local_tempfile()
  withr::local_envvar(DBSNP_GCF_PATH = tmp)
  expect_error(dbsnp_gcf_path(), "does not exist")
})
