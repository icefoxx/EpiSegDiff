test_that("assay contract rejects mixed families", {
  expect_error(validate_assay_contract("count", "bigwig"), "requires")
  expect_error(validate_assay_contract("continuous", "bam"), "requires")
  expect_silent(validate_assay_contract("continuous", "bigwig"))
})

test_that("log2 pseudocount transform is monotonic and finite", {
  x <- matrix(c(0, 0.01, 1, 10), 2, 2)
  y <- transform_continuous(x, "log2_pseudocount", 0.01)
  expect_true(all(is.finite(y)))
  expect_equal(y[1, 1], log2(0.01))
  expect_lt(y[1, 1], y[2, 1])
})

test_that("contrast direction JNKi - DMSO is GAIN when JNKi is higher", {
  x <- rbind(c(0, 0, 2, 2), c(2, 2, 0, 0))
  colnames(x) <- c("d1", "d2", "j1", "j2")
  samp <- data.frame(sample = colnames(x),
                     condition = c("DMSO", "DMSO", "JNKi", "JNKi"))
  eff <- estimate_continuous_effect(x, samp, "JNKi", "DMSO", "replicated")
  expect_gt(eff$effect[1], 0)
  expect_lt(eff$effect[2], 0)
  expect_match(eff$contrast, "JNKi - DMSO")
})

test_that("unreplicated mode cannot emit FDR", {
  x <- matrix(rnorm(20), 10, 2)
  colnames(x) <- c("c", "t")
  samp <- data.frame(sample = c("c", "t"), condition = c("DMSO", "JNKi"))
  expect_error(infer_capability_level(samp, allow_descriptive = FALSE), "Insufficient")
  expect_equal(infer_capability_level(samp, allow_descriptive = TRUE), "descriptive")
  eff <- estimate_continuous_effect(x, samp, "JNKi", "DMSO", "descriptive")
  expect_true(all(is.na(eff$pvalue)))
  expect_equal(eff$inference_mode, "descriptive")
})

test_that("common grid coordinates are 0-based half-open after conversion", {
  skip_if_no_edger()
  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) skip("GenomeInfoDb")
  si <- GenomeInfoDb::Seqinfo("chr1", 1000)
  bins <- build_genomic_grid(si, bin_size = 100L)
  bed <- granges_to_bed_frame(bins)
  expect_equal(bed$start[1], 0L)
  expect_equal(bed$end[1], 100L)
  rec <- reconstruct_interval(bed, data.frame(chrom = "chr1", start = 300L, end = 700L))
  expect_equal(rec$start, 300L)
  expect_equal(rec$end, 700L)
})

test_that("primary schema validation works", {
  expect_error(validate_primary_schema(data.frame(chrom = "1")), "missing")
  df <- standardize_primary_results(data.frame(
    chrom = "1", start = 0L, end = 200L, logFC = -1, pvalue = 0.01, FDR = 0.02,
    direction = "LOSS", level1 = "LOSS"
  ), meta = list(assay_family = "continuous", input_type = "bigwig",
                 contrast = "JNKi - DMSO", inference_mode = "replicated"))
  expect_silent(validate_primary_schema(df))
  expect_equal(df$log2FC, -1)
})

test_that("min_bins is respected in genomic width on uniform bins", {
  n <- 30L
  bs <- 200L
  y <- c(rep(0, 10), rep(4, 10), rep(0, 10))
  pos_start <- as.integer((seq_len(n) - 1L) * bs)
  seg <- segment_effect(y, positions = pos_start, penalty = 2, min_width = 3L,
                        engine = "op")
  expect_true(all(seg$n_bins >= 3L))
  # width from start of first bin to start of last bin is (n_bins-1)*bs;
  # genomic [start,end) using end = start_last + bs is n_bins * bs
  widths <- (seg$end_idx - seg$start_idx + 1L) * bs
  expect_true(all(widths >= 3L * bs))
})
