#!/usr/bin/env Rscript
# GSE110378 H3K27me3 Control vs Kdm6b-KO. Filter = quantile 0.75 (predeclared).
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("results/v0.5/datasetB", recursive = TRUE, showWarnings = FALSE)
root <- "data/GSE110378"
files <- c(
  Control_R1 = file.path(root, "GSM2989282_HSC_Control_1_H3K27me3.bw"),
  Control_R2 = file.path(root, "GSM2989283_HSC_Control_2_H3K27me3.bw"),
  KO_R1 = file.path(root, "GSM2989284_6bKO_1_H3K27me3.bw"),
  KO_R2 = file.path(root, "GSM2989285_6bKO_2_H3K27me3.bw")
)
if (any(!file.exists(files))) stop("Missing Dataset B BigWigs")
sheet <- data.frame(
  sample = names(files), condition = c("Control", "Control", "KO", "KO"),
  replicate = c(1, 2, 1, 2), file = unname(files), stringsAsFactors = FALSE
)
message("import 200 bp")
assay <- read_bigwig_assay(sheet, genome = "mm10", bin_size = 200L, summary_fun = "mean")
bed <- granges_to_bed_frame(assay$bins)
fit <- analyze_continuous_matrix(
  bed, assay$signal, sheet, "KO", "Control", bin_size = 200L,
  filter_method = "quantile", filter_quantile = 0.75,
  penalty = 24, min_bins = 3L, engine = "pelt"
)
utils::write.table(fit$tested, "results/v0.5/datasetB/primary_calls_episegdiff.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
# POST HOC conventional: hiddenDomains peaks not in GEO TAR.
# Use contiguous quantile-valid bins as published-peak substitute.
keep <- filter_signal_bins(assay$signal, "quantile", quantile = 0.75)
islands <- contiguous_valid_regions(data.frame(
  chrom = bed$chrom, start = bed$start, end = bed$end, valid = keep
))
conv_bed <- data.frame(chrom = islands$chrom, start = islands$start, end = islands$end)
rsig <- summarize_matrix_over_intervals(conv_bed, bed, assay$signal)
conv <- test_fixed_regions_limma(conv_bed, rsig, sheet, "KO", "Control")
conv <- standardize_primary_results(conv, meta = list(
  assay_family = "continuous", input_type = "bigwig",
  contrast = "KO - Control", inference_mode = "replicated",
  signal_filter_method = "quantile_islands_POSTHOC",
  transform_method = "log2_pseudocount", segmentation_method = "fixed_islands"
))
utils::write.table(conv, "results/v0.5/datasetB/primary_calls_conventional.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
cl <- classify_region_sets(fit$tested, conv)
ov <- overlap_geometry_metrics(fit$tested, conv)
ov$n_opposite_1bp <- cl$n_opposite_overlaps
utils::write.table(ov, "results/v0.5/datasetB/overlap_metrics.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
yaml::write_yaml(list(
  dataset = "GSE110378",
  GAIN = fit$n_GAIN, LOSS = fit$n_LOSS, n_tested = fit$n_tested,
  n_valid = fit$n_valid, n_input = fit$n_input_bins,
  conv_GAIN = sum(conv$level1 == "GAIN"), conv_LOSS = sum(conv$level1 == "LOSS"),
  conventional_note = "POST HOC islands; hiddenDomains peaks not deposited",
  contrast = "KO - Control", expected = "GAIN of H3K27me3"
), "results/v0.5/datasetB/run_manifest.yaml")
message("Dataset B ESD GAIN/LOSS=", fit$n_GAIN, "/", fit$n_LOSS,
        " conv=", sum(conv$level1 == "GAIN"), "/", sum(conv$level1 == "LOSS"))
