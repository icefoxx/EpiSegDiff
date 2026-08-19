test_that("known synthetic fold changes are recovered with the correct sign", {
  skip_if_no_edger()
  ev <- data.frame(event_class = "GAIN", start_bin = 80L, end_bin = 120L, effect = 2)
  sim <- simulate_count_assay(n_bins = 200, n_rep = 3, depth = 40, dispersion = 0.05,
                              background_sd = 0.05, events = ev, seed = 7)
  norm <- normalize_counts(sim$counts, method = "library_size")
  eff <- estimate_differential_effect(sim$counts, sim$samples, norm$offsets,
                                      score = "logFC", valid = sim$bins$valid)
  fg <- median(eff$logFC[80:120], na.rm = TRUE)
  bg <- median(eff$logFC[1:40], na.rm = TRUE)
  expect_gt(fg, 1)
  expect_gt(fg - bg, 1)
})

test_that("null simulations are centered near zero", {
  skip_if_no_edger()
  sim <- simulate_null_counts(n_bins = 150, n_rep = 3, depth = 25,
                              dispersion = 0.08, background_sd = 0.05, seed = 3)
  norm <- normalize_counts(sim$counts, method = "library_size")
  eff <- estimate_differential_effect(sim$counts, sim$samples, norm$offsets,
                                      score = "logFC")
  expect_lt(abs(median(eff$logFC, na.rm = TRUE)), 0.25)
  expect_lt(abs(mean(eff$z[is.finite(eff$z)])), 0.4)
})

test_that("replicate column reordering does not change estimates", {
  skip_if_no_edger()
  ev <- data.frame(event_class = "GAIN", start_bin = 30L, end_bin = 50L, effect = 1.8)
  sim <- simulate_count_assay(n_bins = 80, n_rep = 3, depth = 30, dispersion = 0.05,
                              background_sd = 0.02, events = ev, seed = 11)
  norm <- normalize_counts(sim$counts, method = "library_size")
  e1 <- estimate_differential_effect(sim$counts, sim$samples, norm$offsets, score = "logFC")
  ord <- c(3, 1, 2, 6, 4, 5)
  e2 <- estimate_differential_effect(sim$counts[, ord], sim$samples[ord, ],
                                     norm$offsets[ord], score = "logFC")
  expect_equal(e1$logFC, e2$logFC, tolerance = 1e-6)
})

test_that("swapping condition labels reverses the sign of beta", {
  skip_if_no_edger()
  ev <- data.frame(event_class = "GAIN", start_bin = 40L, end_bin = 70L, effect = 2)
  sim <- simulate_count_assay(n_bins = 100, n_rep = 3, depth = 30, dispersion = 0.05,
                              background_sd = 0.02, events = ev, seed = 4)
  norm <- normalize_counts(sim$counts, method = "library_size")
  e1 <- estimate_differential_effect(sim$counts, sim$samples, norm$offsets, score = "logFC")
  swapped <- sim$samples
  swapped$condition <- ifelse(swapped$condition == "control", "treatment", "control")
  e2 <- estimate_differential_effect(sim$counts, swapped, norm$offsets, score = "logFC")
  expect_equal(e1$logFC, -e2$logFC, tolerance = 1e-6)
})

test_that("effect track does not interpolate invalid bins", {
  skip_if_no_edger()
  sim <- simulate_null_counts(n_bins = 40, n_rep = 2, depth = 15, seed = 1,
                              mask_bins = 20:22)
  norm <- normalize_counts(sim$counts, method = "library_size")
  eff <- estimate_differential_effect(sim$counts, sim$samples, norm$offsets,
                                      valid = sim$bins$valid)
  tr <- write_effect_track(sim$bins, eff)
  expect_false(any(tr$valid[20:22]))
  expect_true(all(is.na(tr$score[20:22])))
})
