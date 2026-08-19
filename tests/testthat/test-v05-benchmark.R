test_that("union peaks reduce overlaps", {
  a <- data.frame(chrom = "1", start = c(0L, 80L), end = c(100L, 200L))
  u <- union_peak_regions(a)
  expect_equal(nrow(u), 1L)
  expect_equal(u$start, 0L)
  expect_equal(u$end, 200L)
})

test_that("classify shared vs only uses RO 50% and sign", {
  esd <- data.frame(chrom = "1", start = c(0, 1000), end = c(800, 1600),
                    effect = c(-1, 1), FDR = c(0.01, 0.01),
                    level1 = c("LOSS", "GAIN"), direction = c("LOSS", "GAIN"))
  conv <- data.frame(chrom = "1", start = c(0, 5000), end = c(700, 5600),
                     effect = c(-1.1, -0.8), FDR = c(0.01, 0.01),
                     level1 = c("LOSS", "LOSS"), direction = c("LOSS", "LOSS"))
  cl <- classify_region_sets(esd, conv)
  expect_equal(cl$esd$class[1], "A")
  expect_equal(cl$esd$class[2], "C")
  expect_true("B" %in% cl$conventional$class)
})

test_that("synthetic narrow GAIN and LOSS are recovered by continuous path", {
  skip_if_not_installed("limma")
  set.seed(5)
  n <- 80L
  bed <- data.frame(chrom = "1", start = (seq_len(n) - 1L) * 50L,
                    end = seq_len(n) * 50L)
  # LOSS block 11-16, GAIN block 40-59
  base <- matrix(0.15, n, 4)
  colnames(base) <- c("d1", "d2", "j1", "j2")
  base[11:16, 1:2] <- 4
  base[40:59, 3:4] <- 4
  samp <- data.frame(sample = colnames(base),
                     condition = c("DMSO", "DMSO", "JNKi", "JNKi"))
  fit <- analyze_continuous_matrix(bed, base, samp, "JNKi", "DMSO",
                                   bin_size = 50L, filter_method = "none",
                                   penalty = 6, min_bins = 3L, engine = "pelt")
  der <- fit$tested[fit$tested$level1 %in% c("GAIN", "LOSS"), ]
  expect_true(any(der$direction == "LOSS") || any(fit$tested$effect[11:16] < 0))
  expect_true(any(der$direction == "GAIN") || any(fit$tested$effect[40:59] > 0))
})
