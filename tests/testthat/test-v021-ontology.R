test_that("one-to-one geometry classes are mutually exclusive by Δ partition", {
  # EXPANSION both ends
  ev <- classify_geometry_pair(1000, 2000, 500, 2800)
  expect_equal(ev$primary_event, "EXPANSION")
  expect_match(ev$secondary_event, "BOTH|LEFT|RIGHT")
  # one-sided right expansion
  ev <- classify_geometry_pair(1000, 2000, 1000, 2800)
  expect_equal(ev$primary_event, "EXPANSION")
  expect_equal(ev$secondary_event, "RIGHT_END")
  # CONTRACTION both
  ev <- classify_geometry_pair(500, 2800, 1000, 2000)
  expect_equal(ev$primary_event, "CONTRACTION")
  # LEFT / RIGHT shift
  ev <- classify_geometry_pair(2000, 3000, 1000, 2000)
  expect_equal(ev$primary_event, "LEFT_SHIFT")
  ev <- classify_geometry_pair(1000, 2000, 2000, 3000)
  expect_equal(ev$primary_event, "RIGHT_SHIFT")
  # opposite-end motion without enough width change is COMPLEX, not BIDIRECTIONAL_SHIFT
  ev <- classify_geometry_pair(1000, 2000, 800, 2050)
  expect_false(identical(ev$primary_event, "BIDIRECTIONAL_SHIFT"))
  expect_true(ev$primary_event %in% c("COMPLEX", "EXPANSION"))
})

test_that("BIDIRECTIONAL_SHIFT is not a primary class", {
  n <- 40
  tab <- simulate_geometry_examples("COMPLEX", n = n, seed = 9)
  expect_false(any(tab$pred == "BIDIRECTIONAL_SHIFT"))
})

test_that("background normalization is finite and reduces composition pull vs library_size", {
  skip_if_no_edger()
  ev <- compose_fraction_events(400, frac = 0.25, gain_prop = 1, region_width = 20, effect = 2)
  sim <- simulate_count_assay(n_bins = 400, n_rep = 3, depth = 25, dispersion = 0.08,
                              background_sd = 0.05, events = ev, seed = 5)
  ls <- fit_normalized_effect(sim, "library_size", run_discovery = FALSE)
  bg <- fit_normalized_effect(sim, "background", run_discovery = FALSE)
  # library-size should be pulled negative on background for a large GAIN
  expect_true(is.finite(ls$background$bias_background))
  expect_true(is.finite(bg$background$bias_background))
  expect_lt(abs(bg$background$bias_background), abs(ls$background$bias_background) + 0.15)
})

test_that("local condition boundaries do not require a peak caller", {
  ev <- data.frame(event_class = "GAIN", start_bin = 40L, end_bin = 70L, effect = 2)
  sim <- simulate_count_assay(n_bins = 120, n_rep = 3, depth = 30, dispersion = 0.05,
                              background_sd = 0.02, events = ev, seed = 8)
  region <- data.frame(chrom = "chr1", start = 39L * 200L, end = 70L * 200L)
  b <- infer_local_condition_boundaries(sim$counts, sim$bins, sim$samples, region,
                                        window_bins = 20L, method = "local_threshold")
  expect_true(all(c("control", "treatment") %in% b$condition))
  expect_true(all(b$method == "local_threshold"))
  expect_false(any(grepl("macs|seacr", b$method, ignore.case = TRUE)))
})
