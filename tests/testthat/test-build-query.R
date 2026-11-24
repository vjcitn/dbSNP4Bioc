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
  expect_true(all(by_pos$release == "test"))
  expect_true(all(by_pos$build == "GRCh38"))

  rng <- dbsnp_in_range("17", 1, 250, release = "test", build = "GRCh38", cache = cache)
  expect_equal(sort(rng$pos), c(100, 200, 200))
  expect_true(all(rng$release == "test" & rng$build == "GRCh38"))

  rng_bi <- dbsnp_in_range("17", 1, 1000, biallelic_only = TRUE,
                           release = "test", build = "GRCh38", cache = cache)
  expect_true(all(rng_bi$n_alt == 1L))
  expect_setequal(rng_bi$pos, c(100, 500))
})

test_that("building a second, different chromosome does not error on the shared manifest", {
  ## Regression test: the second dbsnp_build_chrom() call for a given
  ## release/build reuses (overwrites in place) the cached manifest.json
  ## entry. bfcupdate(fpath = ...) errors on rtype = "local" resources
  ## ("resource rtype is not 'web'"); .bfc_put() must use file.copy() instead.
  tmp <- withr::local_tempdir()
  vcf <- make_test_vcf_2chrom(tmp)
  cache <- dbsnp_cache(withr::local_tempdir())

  dbsnp_build_chrom("chr17", release = "b157", build = "GRCh38", vcf = vcf, cache = cache)
  expect_no_error(
    dbsnp_build_chrom("chr22", release = "b157", build = "GRCh38", vcf = vcf, cache = cache)
  )

  cached <- dbsnp_cached("b157", "GRCh38", cache)
  expect_setequal(cached$chrom, c("17", "22"))
  expect_true(all(file.exists(cached$path)))
})

test_that("dbsnp_by_rsid finds a SNP regardless of which cached chromosome it is on", {
  tmp <- withr::local_tempdir()
  vcf <- make_test_vcf_2chrom(tmp)
  cache <- dbsnp_cache(withr::local_tempdir())

  dbsnp_build_chrom("chr17", release = "b157", build = "GRCh38", vcf = vcf, cache = cache)
  dbsnp_build_chrom("chr22", release = "b157", build = "GRCh38", vcf = vcf, cache = cache)

  hit17 <- dbsnp_by_rsid("rs1", release = "b157", build = "GRCh38", cache = cache)
  expect_equal(hit17$chrom, "17")
  expect_equal(hit17$RefSNP_id, "rs1")

  hit22 <- dbsnp_by_rsid(20, release = "b157", build = "GRCh38", cache = cache)
  expect_equal(hit22$chrom, "22")

  both <- dbsnp_by_rsid(c("rs1", "rs20"), release = "b157", build = "GRCh38", cache = cache)
  expect_setequal(both$chrom, c("17", "22"))
  expect_true(all(both$release == "b157" & both$build == "GRCh38"))

  none <- dbsnp_by_rsid("rs999999", release = "b157", build = "GRCh38", cache = cache)
  expect_equal(nrow(none), 0L)
})

test_that("dbsnp_gcf_path validates the environment variable", {
  withr::local_envvar(DBSNP_GCF_PATH = "")
  expect_error(dbsnp_gcf_path(), "DBSNP_GCF_PATH is not set")

  tmp <- withr::local_tempfile()
  withr::local_envvar(DBSNP_GCF_PATH = tmp)
  expect_error(dbsnp_gcf_path(), "does not exist")
})
