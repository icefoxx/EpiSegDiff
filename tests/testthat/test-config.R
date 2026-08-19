test_that("example.yaml parses and validates", {
  expect_error(read_episegdiff_config(tempfile()), "does not exist")
  ex <- system.file("extdata", "example.yaml", package = "EpiSegDiff")
  if (!nzchar(ex) || !file.exists(ex)) {
    ex <- file.path(testthat::test_path(), "..", "..", "inst", "extdata", "example.yaml")
  }
  if (!file.exists(ex)) skip("example.yaml not installed")
  cfg <- read_episegdiff_config(ex)
  expect_equal(cfg$preprocessing$coordinate_system, "bed_0based")
  expect_equal(cfg$differential$score, "standardized_effect")
  expect_false(identical(cfg$assay$type, "methylation") &&
                 identical(cfg$testing$method, "edgeR_ql"))
})

test_that("default_config validates", {
  cfg <- default_config()
  cfg$samples <- list(
    list(sample_id = "c1", condition = "control", replicate = 1,
         assay = "cutnrun", target = "x", bam_or_methylation_path = "a.bam"),
    list(sample_id = "t1", condition = "treatment", replicate = 1,
         assay = "cutnrun", target = "x", bam_or_methylation_path = "b.bam")
  )
  expect_silent(validate_config(cfg))
})

test_that("methylation cannot use edgeR_ql", {
  cfg <- default_config()
  cfg$assay$type <- "methylation"
  cfg$testing$method <- "edgeR_ql"
  cfg$samples <- list(
    list(sample_id = "c1", condition = "control", replicate = 1,
         assay = "methylation", target = "cpg", bam_or_methylation_path = "a.tsv"),
    list(sample_id = "t1", condition = "treatment", replicate = 1,
         assay = "methylation", target = "cpg", bam_or_methylation_path = "b.tsv")
  )
  expect_error(validate_config(cfg), "cannot use testing.method = edgeR_ql")
})

test_that("sample IDs must be unique and two conditions required", {
  expect_error(
    validate_sample_metadata(data.frame(
      sample_id = c("a", "a"), condition = c("control", "treatment"),
      replicate = 1:2, assay = "cutnrun", target = "x",
      bam_or_methylation_path = c("a", "b"), stringsAsFactors = FALSE
    )),
    "unique"
  )
  expect_error(
    validate_sample_metadata(data.frame(
      sample_id = "a", condition = "control", replicate = 1,
      assay = "cutnrun", target = "x", bam_or_methylation_path = "a",
      stringsAsFactors = FALSE
    )),
    "two conditions"
  )
})
