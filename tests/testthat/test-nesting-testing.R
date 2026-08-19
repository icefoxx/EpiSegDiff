test_that("nesting flags partial overlap and does not force a tree", {
  segs <- data.frame(
    segment_id = 1:3,
    chrom = "chr1",
    start = c(0L, 20L, 10L),
    end = c(100L, 80L, 50L),
    penalty = c(40, 10, 10),
    effect_mean = c(0.2, 1.0, 0.9),
    stringsAsFactors = FALSE
  )
  nest <- analyze_nesting(segs)
  expect_equal(nest$nesting_status[1], "root")
  expect_true(all(nest$nesting_status[2:3] %in% c("nested", "root")))
  # overlapping children 2 and 3
  segs2 <- data.frame(
    segment_id = 1:2, chrom = "chr1",
    start = c(0L, 40L), end = c(80L, 120L),
    penalty = c(5, 5), effect_mean = c(1, 1)
  )
  nest2 <- analyze_nesting(segs2)
  expect_true(any(nest2$conflict_status == "partial_overlap"))
})

test_that("candidate strategies are distinct and configurable", {
  y <- c(rep(0, 30), rep(4, 30), rep(0, 30))
  path <- segment_penalty_path(y, penalties = c(2, 20), min_width = 2, engine = "op")
  path$chrom <- "chr1"
  # convert idx to dummy genomic
  path$start <- path$start_idx
  path$end <- path$end_idx
  a <- select_candidate_segments(path, "single_selected_penalty", selected_penalty = 20)
  b <- select_candidate_segments(analyze_nesting(path), "multi-scale_consensus")
  expect_true(all(abs(a$penalty - 20) < 1e-9))
  expect_true(nrow(b) >= 1L)
})

test_that("region recount uses raw counts not the score", {
  bins <- build_bins(c(chr1 = 400L), 100L)
  counts <- matrix(c(1L, 2L, 3L, 4L, 10L, 20L, 30L, 40L), nrow = 4)
  colnames(counts) <- c("s1", "s2")
  regions <- data.frame(chrom = "chr1", start = 100L, end = 300L, region_id = 1L)
  rec <- EpiSegDiff:::.recount_from_count_matrix(regions, bins, counts)
  expect_equal(rec$s1, 2L + 3L)
  expect_equal(rec$s2, 20L + 30L)
})

test_that("edgeR QL testing returns the required columns and a FDR caveat", {
  skip_if_no_edger()
  ev <- data.frame(event_class = "GAIN", start_bin = 50L, end_bin = 80L, effect = 2)
  sim <- simulate_count_assay(n_bins = 120, n_rep = 3, depth = 25, dispersion = 0.05,
                              events = ev, seed = 2)
  regions <- data.frame(chrom = "chr1", start = 49 * 200L, end = 80 * 200L,
                        region_id = 1L, stringsAsFactors = FALSE)
  rec <- EpiSegDiff:::.recount_from_count_matrix(regions, sim$bins, sim$counts)
  norm <- normalize_counts(sim$counts, method = "library_size")
  tt <- test_regions_edger(rec, sim$samples, offsets = norm$offsets)
  expect_true(all(c("region_id", "logFC", "logCPM", "statistic", "pvalue", "FDR") %in% names(tt)))
  expect_match(tt$fdr_interpretation[1], "not guaranteed")
  expect_gt(tt$logFC[1], 0)
})
