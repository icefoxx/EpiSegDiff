test_that("boundary refinement moves toward a sharp gradient", {
  track <- data.frame(
    chrom = "chr1",
    start = seq(0L, 1900L, by = 100L),
    end = seq(100L, 2000L, by = 100L),
    score = c(rep(0, 8), 5, 5, 5, 5, 5, rep(0, 7)),
    valid = TRUE
  )
  regions <- data.frame(region_id = 1L, chrom = "chr1", start = 600L, end = 1600L)
  ref <- refine_boundaries(regions, track, method = "local_gradient",
                           search_window_bins = 4L)
  expect_equal(ref$coarse_start, 600L)
  expect_true(abs(ref$refined_start - 800L) <= 100L)
})

test_that("methylation is never sent through the NB count model", {
  meth <- data.frame(
    chrom = rep("chr1", 12),
    position = rep(c(10L, 20L, 30L), each = 4),
    methylated_count = c(1, 2, 18, 19, 10, 11, 10, 12, 2, 1, 2, 3),
    total_count = rep(20L, 12),
    sample = rep(c("c1", "c2", "t1", "t2"), 3),
    stringsAsFactors = FALSE
  )
  samples <- data.frame(
    sample_id = c("c1", "c2", "t1", "t2"),
    condition = c("control", "control", "treatment", "treatment"),
    replicate = c(1, 2, 1, 2),
    assay = "methylation", target = "cpg",
    bam_or_methylation_path = NA_character_,
    stringsAsFactors = FALSE
  )
  meth <- read_methylation_table(meth, min_coverage = 5)
  tr <- estimate_methylation_effect(meth, samples, score = "delta_beta")
  expect_true(all(c("effect", "uncertainty", "score", "valid") %in% names(tr)))
  expect_true(any(is.finite(tr$effect)))
  # low coverage invalid
  meth2 <- meth
  meth2$total_count <- 1L
  meth2$methylated_count <- 0L
  meth2 <- read_methylation_table(meth2, min_coverage = 5)
  expect_false(any(meth2$valid))
})

test_that("missing CpGs stay missing, not zero methylation", {
  df <- data.frame(chrom = "chr1", position = 10L, methylated_count = 3L,
                   total_count = 10L)
  out <- read_methylation_table(df)
  expect_equal(nrow(out), 1L)
  expect_false(any(out$position == 11L))
})

test_that("in-memory count pipeline is deterministic under a fixed seed", {
  skip_if_no_edger()
  ev <- data.frame(event_class = "GAIN", start_bin = 60L, end_bin = 90L, effect = 2)
  sim <- simulate_count_assay(n_bins = 150, n_rep = 3, depth = 30, dispersion = 0.05,
                              background_sd = 0.05, events = ev, seed = 21)
  cfg <- default_config()
  cfg$samples <- sim$samples
  cfg$output$directory <- tempfile("esd")
  cfg$normalization$method <- "library_size"
  cfg$segmentation$penalties <- c(8, 20)
  cfg$segmentation$selected_penalty <- 8
  cfg$segmentation$engine <- "pelt"
  cfg$output$write_bigwig <- FALSE
  r1 <- run_episegdiff(cfg, count_matrix = sim)
  r2 <- run_episegdiff(cfg, count_matrix = sim)
  expect_equal(r1$track$score, r2$track$score)
  expect_equal(r1$candidates$start, r2$candidates$start)
  expect_true(file.exists(file.path(cfg$output$directory, "run_manifest.tsv")))
  expect_true(file.exists(file.path(cfg$output$directory, "normalization_diagnostics.tsv")))
  expect_true(file.exists(file.path(cfg$output$directory, "config.resolved.yaml")))
})

test_that("ATAC insertion helper builds 1-bp intervals from a fragment table", {
  # Adapter path that does not need a BAM: insertion counting mode.
  ins <- data.frame(chrom = "chr1", start = c(10L, 25L, 25L), end = c(11L, 26L, 26L))
  bins <- build_bins(c(chr1 = 50L), bin_size = 10L)
  ct <- count_fragments_into_bins(ins, bins, "insertion")
  expect_equal(sum(ct), 3L)
  expect_equal(ct[bins$start == 10L], 1L)
  expect_equal(ct[bins$start == 20L], 2L)
})

test_that("empirical FDR helper treats null discoveries as false", {
  tested <- data.frame(FDR = c(0.001, 0.2, 0.04), start = c(0, 10, 20), end = c(5, 15, 30))
  tab <- empirical_fdr(tested, nominal = c(0.01, 0.05))
  expect_equal(tab$n_discoveries[1], 1L)
  expect_equal(tab$region_level_empirical_fdr[1], 1)
  expect_equal(tab$n_discoveries[2], 2L)
})
