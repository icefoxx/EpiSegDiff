test_that("pathb additive returns broad/local components", {
  ev <- data.frame(start_bin = c(40L, 60L), end_bin = c(80L, 61L),
                   delta = c(1.2, 2), shape = "sharp")
  sim <- simulate_continuous_events(120L, 200L, ev, noise_sd = 0.05, seed = 3)
  fit <- run_pathb(sim, "A_additive")
  expect_true(nrow(fit$regions) >= 1L)
  expect_true(all(c("component_type", "region_id") %in% names(fit$regions)))
})

test_that("nested hierarchy scorer distinguishes merge vs both", {
  truth <- data.frame(start_bin = c(10L, 20L), end_bin = c(40L, 21L), delta = c(1, 1))
  flat <- data.frame(chrom = "1", start = 1800L, end = 8200L, effect = 1,
                     direction = "GAIN", level1 = "GAIN", component_type = "flat",
                     region_id = 1L, parent_id = NA_integer_, level = 1L,
                     scale_bp = 6400L, support = 1, stringsAsFactors = FALSE)
  s1 <- score_nested_hierarchy(flat, truth, 200L)
  expect_equal(s1$p_merged, 1)
  expect_equal(s1$p_both, 0)
  hier <- rbind(
    data.frame(chrom = "1", start = 1800L, end = 8200L, effect = 1.2,
               direction = "GAIN", level1 = "GAIN", component_type = "broad",
               region_id = 1L, parent_id = NA_integer_, level = 1L,
               scale_bp = 6400L, support = 1.2, stringsAsFactors = FALSE),
    data.frame(chrom = "1", start = 4000L, end = 4200L, effect = 1.5,
               direction = "GAIN", level1 = "GAIN", component_type = "local",
               region_id = 2L, parent_id = 1L, level = 2L,
               scale_bp = 200L, support = 1.5, stringsAsFactors = FALSE)
  )
  s2 <- score_nested_hierarchy(hier, truth, 200L)
  expect_equal(s2$p_both, 1)
  expect_equal(s2$p_merged, 0)
})
