test_that("noiseless positive block has zero left and right error (symmetric walk)", {
  positions <- list(c(3L, 7L), c(5L, 9L), c(2L, 6L), c(4L, 10L))
  for (pr in positions) {
    st <- make_step_series(13L, pr[1], pr[2], height = 2, bin_size = 100L)
    got <- end_walk_series(st$y, st$start, st$end, st$true_L, st$true_R,
                           search_radius = 4L, version = "symmetric")
    expect_equal(got$L, st$true_L)
    expect_equal(got$R, st$true_R)
  }
})

test_that("noiseless negative block is not direction-dependent", {
  st <- make_step_series(13L, 5L, 9L, height = -2, bin_size = 100L)
  got <- end_walk_series(st$y, st$start, st$end, st$true_L, st$true_R,
                         search_radius = 4L, version = "symmetric")
  expect_equal(got$L, st$true_L)
  expect_equal(got$R, st$true_R)
})

test_that("reverse-signal invariance holds for the symmetric walker", {
  st <- make_step_series(20L, 6L, 12L, height = 2, bin_size = 50L)
  fwd <- end_walk_series(st$y, st$start, st$end, st$true_L, st$true_R,
                         search_radius = 5L, version = "symmetric")
  rev <- reverse_series(st$y, st$start, st$end)
  true_rev <- reverse_interval(st$true_L, st$true_R, rev$N)
  cand_rev <- reverse_interval(st$true_L, st$true_R, rev$N)
  bwd <- end_walk_series(rev$y, rev$start, rev$end, cand_rev$L, cand_rev$R,
                         search_radius = 5L, version = "symmetric")
  expect_equal(fwd$L - st$true_L, true_rev$R - bwd$R)
  expect_equal(st$true_R - fwd$R, bwd$L - true_rev$L)
  expect_equal(bwd$L, true_rev$L)
  expect_equal(bwd$R, true_rev$R)
})

test_that("Level I calls do not change when topology is disabled", {
  skip_if_no_edger()
  ev <- example_topology_event("GAIN", 40L)
  sim <- simulate_count_assay(n_bins = 120, n_rep = 3, depth = 25, dispersion = 0.05,
                              background_sd = 0.04, events = ev, seed = 22)
  cfg <- default_config()
  cfg$samples <- sim$samples
  cfg$normalization$method <- "library_size"
  cfg$segmentation$penalties <- c(8, 16)
  cfg$segmentation$selected_penalty <- 8
  cfg$output$write_bed <- FALSE
  cfg$topology$enabled <- TRUE
  cfg$output$directory <- tempfile("on")
  a <- run_episegdiff(cfg, count_matrix = sim)
  cfg$topology$enabled <- FALSE
  cfg$output$directory <- tempfile("off")
  b <- run_episegdiff(cfg, count_matrix = sim)
  cols <- intersect(c("region_id", "chrom", "start", "end", "logFC", "FDR"),
                    names(a$tested))
  expect_equal(a$tested[, cols], b$tested[, cols])
  expect_equal(nrow(b$events), 0L)
})
