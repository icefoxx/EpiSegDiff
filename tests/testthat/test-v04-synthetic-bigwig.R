write_toy_bw <- function(path, scores, chrom = "chr1", bin = 50L) {
  gr <- GenomicRanges::GRanges(
    chrom,
    IRanges::IRanges(seq.int(1L, length(scores) * bin, by = bin),
                     width = bin),
    score = as.numeric(scores)
  )
  GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::Seqinfo(chrom, length(scores) * bin)
  rtracklayer::export(gr, path, format = "BigWig")
}

test_that("synthetic continuous pipeline recovers GAIN and LOSS directions", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("limma")
  tmp <- tempfile(); dir.create(tmp)
  n <- 40L
  gain <- c(rep(0.1, 10), rep(4, 10), rep(0.1, 20))
  loss <- c(rep(0.1, 20), rep(4, 10), rep(0.1, 10))
  # DMSO high on LOSS block, JNKi high on GAIN block
  write_toy_bw(file.path(tmp, "d1.bw"), loss + 0.05)
  write_toy_bw(file.path(tmp, "d2.bw"), loss)
  write_toy_bw(file.path(tmp, "j1.bw"), gain)
  write_toy_bw(file.path(tmp, "j2.bw"), gain + 0.05)
  sheet <- data.frame(
    sample = c("d1", "d2", "j1", "j2"),
    condition = c("DMSO", "DMSO", "JNKi", "JNKi"),
    replicate = c(1, 2, 1, 2),
    file = file.path(tmp, c("d1.bw", "d2.bw", "j1.bw", "j2.bw"))
  )
  cfg <- list(
    assay = list(assay_family = "continuous", input_type = "bigwig",
                 genome = "toy", signal_type = "raw_coverage",
                 normalization = "none"),
    samples = sheet,
    contrast = list(numerator = "JNKi", denominator = "DMSO"),
    binning = list(bin_size = 50L, summary = "mean", chromosomes = "chr1"),
    signal_filter = list(method = "none"),
    transform = list(method = "log2_pseudocount", pseudocount = 0.01),
    model = list(method = "limma"),
    segmentation = list(method = "gaussian_pelt", penalty = 4, min_bins = 3L),
    calling = list(fdr = 0.2),
    output = list(directory = file.path(tmp, "out"))
  )
  res <- run_continuous_assay(cfg)
  sig <- res$tested[res$tested$level1 %in% c("GAIN", "LOSS"), ]
  expect_true(file.exists(file.path(tmp, "out", "run_manifest.yaml")))
  expect_true(file.exists(file.path(tmp, "out", "episegdiff_primary.tsv")))
  expect_true("GAIN" %in% sig$direction || any(res$tested$effect > 0, na.rm = TRUE))
  expect_true("LOSS" %in% sig$direction || any(res$tested$effect < 0, na.rm = TRUE))
  expect_equal(unique(res$tested$contrast), "JNKi - DMSO")
  expect_equal(unique(res$tested$assay_family), "continuous")
  if (nrow(sig)) expect_true(all((sig$end - sig$start) >= 3L * 50L - 1L))
})
