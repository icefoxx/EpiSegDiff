test_that("merge_significant_bins joins adjacent same-sign calls", {
  bed <- data.frame(chrom = "1", start = 0:9 * 200L, end = 1:10 * 200L)
  effect <- c(rep(0, 3), 2, 2, 2, 0, -1, -1, 0)
  fdr <- c(rep(1, 3), 0.01, 0.01, 0.01, 1, 0.01, 0.01, 1)
  m <- merge_significant_bins(bed, effect, fdr, rep(TRUE, 10), 0.05)
  expect_equal(nrow(m), 2L)
  expect_equal(m$direction, c("GAIN", "LOSS"))
  expect_equal(m$end - m$start, c(600, 400))
})

test_that("score_calls_vs_truth is perfect on identity", {
  truth <- data.frame(start_bin = 5L, end_bin = 10L, delta = 1.5)
  calls <- data.frame(chrom = "1", start = 800L, end = 2000L,
                      effect = 1.5, direction = "GAIN", level1 = "GAIN")
  sc <- score_calls_vs_truth(calls, truth, 200L)
  expect_equal(sc$n_detected, 1)
  expect_equal(sc$precision, 1)
  expect_equal(sc$recall, 1)
  expect_equal(sc$dir_acc, 1)
  expect_equal(sc$n_merged, 0)
})
