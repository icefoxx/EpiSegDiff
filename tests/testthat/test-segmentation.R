test_that("constant signal has no internal changepoint", {
  y <- rep(0, 40)
  seg <- segment_effect(y, penalty = 5, min_width = 2, loss = "gaussian", engine = "op")
  expect_equal(nrow(seg), 1L)
  seg2 <- segment_effect(y, penalty = 5, min_width = 2, loss = "gaussian", engine = "pelt")
  expect_equal(nrow(seg2), 1L)
})

test_that("one step produces one boundary near the true change", {
  y <- c(rep(0, 25), rep(5, 25))
  seg <- segment_effect(y, penalty = 2, min_width = 2, loss = "gaussian", engine = "op")
  expect_equal(nrow(seg), 2L)
  expect_true(abs(seg$end_idx[1] - 25L) <= 1L)
})

test_that("two steps produce the expected number of segments", {
  y <- c(rep(0, 20), rep(4, 20), rep(-3, 20))
  seg <- segment_effect(y, penalty = 2, min_width = 2, loss = "gaussian", engine = "op")
  expect_equal(nrow(seg), 3L)
})

test_that("NA gaps are handled by splitting CVRs, not interpolating", {
  track <- data.frame(
    chrom = "chr1",
    start = seq(0L, 90L, by = 10L),
    end = seq(10L, 100L, by = 10L),
    score = c(0, 0, 0, NA, NA, 4, 4, 4, 4, 4),
    valid = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  cvrs <- contiguous_valid_regions(track)
  expect_equal(nrow(cvrs), 2L)
  segs <- segment_track(track, penalty = 1, min_width = 1, loss = "gaussian", engine = "op")
  expect_true(all(segs$end <= 30 | segs$start >= 50))
})

test_that("zero-length and sub-min-width series emit no segments", {
  expect_equal(nrow(segment_effect(numeric(), penalty = 1)), 0L)
  one <- segment_effect(0.5, penalty = 1, min_width = 5)
  expect_equal(nrow(one), 0L)
})

test_that("emitted segments have at least min_width bins", {
  set.seed(1)
  y <- c(rep(0, 4), 8, 8, rep(0, 5), rep(3, 6))
  seg <- segment_effect(y, penalty = 2, min_width = 3L, engine = "op")
  expect_true(nrow(seg) >= 1L)
  expect_true(all(seg$n_bins >= 3L))
})

test_that("Huber reduces an isolated-outlier changepoint relative to Gaussian", {
  set.seed(1)
  y <- c(rep(0, 30), 20, rep(0, 30))
  g <- segment_effect(y, penalty = 8, min_width = 1, loss = "gaussian", engine = "op")
  h <- segment_effect(y, penalty = 8, min_width = 1, loss = "huber", engine = "op")
  expect_gte(nrow(g), nrow(h))
  expect_true(nrow(h) <= 3L)
})

test_that("penalty path is deterministic and stores all scales", {
  y <- c(rep(0, 20), rep(3, 20))
  p1 <- segment_penalty_path(y, penalties = c(1, 5, 20), min_width = 2,
                             loss = "gaussian", engine = "op")
  p2 <- segment_penalty_path(y, penalties = c(1, 5, 20), min_width = 2,
                             loss = "gaussian", engine = "op")
  expect_equal(p1, p2)
  expect_true(all(c(1, 5, 20) %in% p1$penalty))
})

test_that("PELT and OP agree on a small step series", {
  y <- c(rep(0, 15), rep(6, 15))
  a <- segment_effect(y, penalty = 3, min_width = 2, engine = "op")
  b <- segment_effect(y, penalty = 3, min_width = 2, engine = "pelt")
  expect_equal(a$start_idx, b$start_idx)
  expect_equal(a$end_idx, b$end_idx)
})
