test_that("Gaussian FPOP matches PELT on tractable series", {
  set.seed(42)
  cases <- list(
    rep(0, 40),
    c(rep(0, 25), rep(5, 25)),
    c(rep(0, 20), rep(4, 20), rep(-3, 20)),
    c(rep(1, 15), rep(-1, 15), rep(2, 20))
  )
  for (y in cases) {
    a <- segment_effect(y, penalty = 3, min_width = 3L, engine = "pelt")
    b <- segment_effect(y, penalty = 3, min_width = 3L, engine = "fpop")
    expect_equal(a$start_idx, b$start_idx)
    expect_equal(a$end_idx, b$end_idx)
  }
  for (i in 1:8) {
    y <- c(rnorm(30, 0, 0.2), rnorm(30, 2, 0.2), rnorm(30, 0, 0.2))
    a <- segment_effect(y, penalty = 8, min_width = 3L, engine = "pelt")
    b <- segment_effect(y, penalty = 8, min_width = 3L, engine = "fpop")
    expect_equal(a$start_idx, b$start_idx)
    expect_equal(a$end_idx, b$end_idx)
  }
})

test_that("segments never cross chromosomes", {
  track <- data.frame(
    chrom = c(rep("1", 20), rep("2", 20)),
    start = c((0:19) * 10L, (0:19) * 10L),
    end = c((1:20) * 10L, (1:20) * 10L),
    score = c(rep(0, 10), rep(4, 10), rep(0, 10), rep(-4, 10)),
    valid = TRUE
  )
  segs <- segment_genome(track, penalty = 2, min_width = 3L, engine = "pelt")$segments
  expect_true(all(segs$chrom %in% c("1", "2")))
  expect_false(any(segs$chrom == "1" & segs$end > 200))
  expect_false(any(segs$chrom == "2" & segs$start < 0 & segs$end > 200 & segs$chrom == "1"))
  expect_equal(length(unique(segs$chrom)), 2L)
})

test_that("genuine gaps become independent blocks", {
  track <- data.frame(
    chrom = "chr1",
    start = seq(0L, 190L, by = 10L),
    end = seq(10L, 200L, by = 10L),
    score = c(rep(0, 6), rep(NA, 4), rep(5, 10)),
    valid = c(rep(TRUE, 6), rep(FALSE, 4), rep(TRUE, 10))
  )
  cvrs <- contiguous_valid_regions(track)
  expect_equal(nrow(cvrs), 2L)
  segs <- segment_chromosome(track, penalty = 2, min_width = 3L, engine = "pelt")
  expect_true(all(segs$end <= 60 | segs$start >= 100))
})

test_that("incompatible checkpoints are rejected and compatible ones resume", {
  track <- data.frame(
    chrom = "1",
    start = (0:29) * 10L,
    end = (1:30) * 10L,
    score = c(rep(0, 15), rep(5, 15)),
    valid = TRUE
  )
  td <- tempfile("ck")
  dir.create(td)
  id <- list(cache_id = "abc", genome = "toy", bin_size = 10L,
             summary_fun = "mean", transform = "log2_pseudocount",
             pseudocount = 0.01, filter = "none",
             segmentation = "gaussian_pelt", penalty = 4, min_bins = 3L)
  a <- segment_genome(track, penalty = 4, min_width = 3L, engine = "pelt",
                      checkpoint_dir = td, resume = TRUE, identity = id)
  expect_false(a$chromosome_runtime$resumed[1])
  b <- segment_genome(track, penalty = 4, min_width = 3L, engine = "pelt",
                      checkpoint_dir = td, resume = TRUE, identity = id)
  expect_true(b$chromosome_runtime$resumed[1])
  expect_equal(a$segments$start, b$segments$start)
  bad <- id
  bad$penalty <- 99
  expect_error(
    segment_genome(track, penalty = 99, min_width = 3L, engine = "pelt",
                   checkpoint_dir = td, resume = TRUE, identity = bad),
    "Incompatible"
  )
})

test_that("workers=1 and workers=2 agree on scientific output", {
  track <- data.frame(
    chrom = c(rep("1", 25), rep("2", 25)),
    start = c((0:24) * 10L, (0:24) * 10L),
    end = c((1:25) * 10L, (1:25) * 10L),
    score = c(rep(0, 12), rep(6, 13), rep(0, 12), rep(-6, 13)),
    valid = TRUE
  )
  id <- list(cache_id = "w", genome = "toy", bin_size = 10L,
             summary_fun = "mean", transform = "log2_pseudocount",
             pseudocount = 0.01, filter = "none",
             segmentation = "gaussian_pelt", penalty = 3, min_bins = 3L)
  s1 <- segment_genome(track, penalty = 3, min_width = 3L, engine = "pelt",
                       workers = 1L, identity = id)$segments
  s2 <- segment_genome(track, penalty = 3, min_width = 3L, engine = "pelt",
                       workers = 2L, identity = id)$segments
  s1 <- s1[order(s1$chrom, s1$start), ]
  s2 <- s2[order(s2$chrom, s2$start), ]
  expect_equal(s1$chrom, s2$chrom)
  expect_equal(s1$start, s2$start)
  expect_equal(s1$end, s2$end)
})

test_that("filter=none synthetic integration recovers both directions", {
  set.seed(7)
  n <- 80L
  bed <- data.frame(chrom = "1", start = (seq_len(n) - 1L) * 50L,
                    end = seq_len(n) * 50L)
  sig <- cbind(
    d1 = c(rep(4, 30), rep(0.1, 50)),
    d2 = c(rep(4.1, 30), rep(0.1, 50)),
    j1 = c(rep(0.1, 30), rep(4, 50)),
    j2 = c(rep(0.1, 30), rep(4.2, 50))
  )
  samp <- data.frame(sample = colnames(sig),
                     condition = c("DMSO", "DMSO", "JNKi", "JNKi"))
  fit <- analyze_continuous_matrix(
    bed, sig, samp, "JNKi", "DMSO", bin_size = 50L,
    filter_method = "none", penalty = 8, min_bins = 3L, engine = "fpop"
  )
  der <- fit$tested[fit$tested$level1 %in% c("GAIN", "LOSS"), ]
  expect_true(any(fit$tested$effect > 0))
  expect_true(any(fit$tested$effect < 0))
  expect_true(all((fit$tested$end - fit$tested$start) >= 3L * 50L - 1L))
  expect_equal(fit$inference_mode, "replicated")
})
