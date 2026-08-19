test_that("transition_truth uses 20% ramps", {
  z <- transition_truth(1L, 20L, 200L)
  expect_equal(z$zone, c("left_transition", "core", "right_transition"))
  expect_equal(z$end[1] - z$start[1], 4L * 200L)
})

test_that("score_transition_zones is 0 on empty calls", {
  ev <- data.frame(start_bin = 10L, end_bin = 30L, delta = 1.5, shape = "gradual")
  empty <- data.frame(chrom = character(), start = integer(), end = integer(),
                      effect = numeric(), direction = character(),
                      stringsAsFactors = FALSE)
  s <- score_transition_zones(empty, 10L, 30L, 200L, 8000L)
  expect_equal(s$core_recall, 0)
  expect_equal(s$n_call, 0L)
})
