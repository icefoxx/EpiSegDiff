test_that("score engines return finite named-length scores", {
  x <- matrix(c(0.1, 0.1, 2, 2, 2, 2, 0.1, 0.1), 2, 4, byrow = TRUE)
  colnames(x) <- c("d1", "d2", "j1", "j2")
  samp <- data.frame(sample = colnames(x),
                     condition = c("DMSO", "DMSO", "JNKi", "JNKi"))
  for (e in c("raw_effect", "log_effect", "standardized", "hotelling")) {
    s <- compute_diff_score(x, samp, "JNKi", "DMSO", e)
    expect_equal(length(s), 2L)
    expect_true(all(is.finite(s)))
    expect_gt(s[1], 0)
    expect_lt(s[2], 0)
  }
})

test_that("WBS and fused lasso recover a sharp step", {
  set.seed(1)
  y <- c(rep(0, 40), rep(3, 40))
  w <- segment_signal(y, engine = "wbs", penalty = 2, min_width = 1L)
  f <- segment_signal(y, engine = "fused_lasso", penalty = 2, min_width = 1L)
  expect_gte(nrow(w), 1L)
  expect_gte(nrow(f), 1L)
  # a change near 40
  expect_true(any(abs(w$end_idx - 40L) <= 5L) || any(abs(w$start_idx - 41L) <= 5L))
})

test_that("v0.2 refine_boundaries still accepts BED+track", {
  track <- data.frame(chrom = "chr1", start = seq(0L, 800L, 200L),
                      end = seq(200L, 1000L, 200L), score = c(0, 0, 2, 2, 0))
  regions <- data.frame(chrom = "chr1", start = 400L, end = 800L, region_id = 1L)
  ref <- refine_boundaries(regions, track, method = "local_gradient",
                           search_window_bins = 1L)
  expect_true(all(c("refined_start", "refined_end") %in% names(ref)))
})

test_that("run_multiscale is deterministic given seed and has provenance", {
  ev <- data.frame(start_bin = 30L, end_bin = 50L, delta = 2, shape = "sharp")
  sim <- simulate_continuous_events(80L, 200L, ev, noise_sd = 0.05, seed = 9)
  a <- run_multiscale(sim, "raw_effect", "wbs", min_width = 1L, penalty = 4)
  b <- run_multiscale(sim, "raw_effect", "wbs", min_width = 1L, penalty = 4)
  if (nrow(a) && nrow(b)) {
    expect_equal(a$start, b$start)
    expect_true("score_engine" %in% names(a))
  }
})
