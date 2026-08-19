make_rle <- function(values, lengths) {
  if (requireNamespace("S4Vectors", quietly = TRUE)) {
    S4Vectors::Rle(values, lengths)
  } else {
    methods::new("Rle", values = values, lengths = as.integer(lengths))
  }
}

test_that("summary_fun mean is length-weighted and max/area are honest", {
  if (!requireNamespace("IRanges", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE)) {
    skip("IRanges/GenomicRanges")
  }
  cov <- make_rle(c(1, 3), c(4, 6))
  bins <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
  GenomeInfoDb::seqinfo(bins) <- GenomeInfoDb::Seqinfo("chr1", 10)
  cl <- IRanges::RleList(chr1 = cov)
  mu <- summarize_rle_over_bins(bins, cl, "mean")
  mx <- summarize_rle_over_bins(bins, cl, "max")
  ar <- summarize_rle_over_bins(bins, cl, "area")
  expect_equal(mu, 2.2, tolerance = 1e-8)
  expect_equal(mx, 3)
  expect_equal(ar, 22)
})

test_that("summary_fun handles constant, edge, uncovered, terminal, nonfinite", {
  if (!requireNamespace("IRanges", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE)) {
    skip("IRanges/GenomicRanges")
  }
  # constant
  cl <- IRanges::RleList(chr1 = make_rle(2, 10))
  bins <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
  expect_equal(summarize_rle_over_bins(bins, cl, "mean"), 2)
  expect_equal(summarize_rle_over_bins(bins, cl, "max"), 2)
  expect_equal(summarize_rle_over_bins(bins, cl, "area"), 20)

  # partial overlap at bin edge: 1x4 + 3x6, bins [1,5] and [6,10]
  cl <- IRanges::RleList(chr1 = make_rle(c(1, 3), c(4, 6)))
  bins <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 6), c(5, 10)))
  expect_equal(summarize_rle_over_bins(bins, cl, "mean"), c(1.4, 3), tolerance = 1e-8)
  expect_equal(summarize_rle_over_bins(bins, cl, "max"), c(3, 3))
  expect_equal(summarize_rle_over_bins(bins, cl, "area"), c(7, 15))

  # uncovered bases encoded as 0 in a full-length Rle
  cl <- IRanges::RleList(chr1 = make_rle(c(0, 5, 0), c(3, 2, 5)))
  bins <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
  expect_equal(summarize_rle_over_bins(bins, cl, "mean"), 1)
  expect_equal(summarize_rle_over_bins(bins, cl, "area"), 10)
  expect_equal(summarize_rle_over_bins(bins, cl, "max"), 5)

  # terminal short chromosome bin
  cl <- IRanges::RleList(chr1 = make_rle(1, 25))
  bins <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 11, 21), c(10, 20, 25)))
  expect_equal(summarize_rle_over_bins(bins, cl, "mean"), c(1, 1, 1))
  expect_equal(summarize_rle_over_bins(bins, cl, "area"), c(10, 10, 5))

  # non-finite values become 0
  cl <- IRanges::RleList(chr1 = make_rle(c(NaN, Inf), c(5, 5)))
  bins <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
  expect_equal(summarize_rle_over_bins(bins, cl, "mean"), 0)
})

test_that("unsupported summary_fun errors instead of silent mean", {
  if (!requireNamespace("IRanges", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE)) {
    skip("IRanges/GenomicRanges")
  }
  expect_error(summarize_rle_over_bins(
    GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 2)),
    IRanges::RleList(chr1 = make_rle(1, 2)),
    "median"
  ), "not implemented")
  expect_error(read_bigwig_assay(
    data.frame(sample = "a", condition = "c", replicate = 1, file = tempfile()),
    summary_fun = "median"
  ), regexp = "should be one of|not implemented|arg")
})

test_that("coarsen_bed_matrix averages linear-scale consecutive bins", {
  bed <- data.frame(chrom = "1", start = c(0L, 100L, 200L, 300L),
                    end = c(100L, 200L, 300L, 400L))
  sig <- cbind(a = c(1, 3, 5, 7), b = c(2, 4, 6, 8))
  out <- coarsen_bed_matrix(bed, sig, factor = 2L)
  expect_equal(out$bed$start, c(0L, 200L))
  expect_equal(out$bed$end, c(200L, 400L))
  expect_equal(out$signal[, "a"], c(2, 6))
  expect_equal(out$signal[, "b"], c(3, 7))
})

test_that("compare_continuous_runs recovers self-concordance", {
  d <- data.frame(chrom = "1", start = c(0, 1000), end = c(800, 1800),
                  effect = c(-1.2, 1.5), log2FC = c(-1.2, 1.5),
                  FDR = c(0.01, 0.01), level1 = c("LOSS", "GAIN"),
                  direction = c("LOSS", "GAIN"))
  c1 <- compare_continuous_runs(d, d, "a", "a", strong_abs_effect = 1)
  expect_equal(c1$effect_pearson, 1)
  expect_equal(c1$sign_concordance, 1)
  expect_equal(c1$DER_bp_jaccard, 1)
})

test_that("reverse contrast flips GAIN/LOSS signs", {
  x <- rbind(c(0, 0, 2, 2))
  colnames(x) <- c("d1", "d2", "j1", "j2")
  samp <- data.frame(sample = colnames(x),
                     condition = c("DMSO", "DMSO", "JNKi", "JNKi"))
  a <- estimate_continuous_effect(x, samp, "JNKi", "DMSO", "replicated")
  b <- estimate_continuous_effect(x, samp, "DMSO", "JNKi", "replicated")
  expect_equal(a$effect, -b$effect, tolerance = 1e-6)
  expect_gt(a$effect, 0)
})

test_that("continuous estimate never reports edgeR as the model", {
  x <- rbind(c(0, 0, 2, 2))
  colnames(x) <- c("d1", "d2", "j1", "j2")
  samp <- data.frame(sample = colnames(x),
                     condition = c("DMSO", "DMSO", "JNKi", "JNKi"))
  a <- estimate_continuous_effect(x, samp, "JNKi", "DMSO", "replicated")
  expect_equal(a$model, "limma")
  expect_false(identical(a$model, "edgeR"))
})
