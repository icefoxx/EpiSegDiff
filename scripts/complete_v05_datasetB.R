#!/usr/bin/env Rscript
# Efficient Dataset B completion: cache 200 bp matrix, island conventional, compare.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("results/v0.5/datasetB", recursive = TRUE, showWarnings = FALSE)
dir.create("cache", showWarnings = FALSE)
rds <- "cache/gse110378_200bp_signal_matrix.rds"
root <- "data/GSE110378"
files <- c(
  Control_R1 = file.path(root, "GSM2989282_HSC_Control_1_H3K27me3.bw"),
  Control_R2 = file.path(root, "GSM2989283_HSC_Control_2_H3K27me3.bw"),
  KO_R1 = file.path(root, "GSM2989284_6bKO_1_H3K27me3.bw"),
  KO_R2 = file.path(root, "GSM2989285_6bKO_2_H3K27me3.bw")
)
sheet <- data.frame(
  sample = names(files), condition = c("Control", "Control", "KO", "KO"),
  replicate = c(1L, 2L, 1L, 2L), file = unname(files), stringsAsFactors = FALSE
)
if (file.exists(rds)) {
  message("load cache")
  cache <- readRDS(rds)
  bed <- cache$bed
  signal <- cache$signal
} else {
  message("import 200 bp and cache")
  assay <- read_bigwig_assay(sheet, genome = "mm10", bin_size = 200L,
                             summary_fun = "mean")
  bed <- granges_to_bed_frame(assay$bins)
  signal <- assay$signal
  saveRDS(list(bed = bed, signal = signal, samples = sheet), rds)
  rm(assay)
  gc()
}

esd_path <- "results/v0.5/datasetB/primary_calls_episegdiff.tsv"
if (file.exists(esd_path)) {
  message("reuse ESD primary")
  tested <- utils::read.table(esd_path, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE)
} else {
  message("ESD quantile 0.75")
  fit <- analyze_continuous_matrix(
    bed, signal, sheet, "KO", "Control", bin_size = 200L,
    filter_method = "quantile", filter_quantile = 0.75,
    penalty = 24, min_bins = 3L, engine = "pelt"
  )
  tested <- fit$tested
  utils::write.table(tested, esd_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Conventional islands: contiguous quantile-valid bins, min 3 bins (POST HOC; no hiddenDomains BED)
message("conventional islands")
keep <- filter_signal_bins(signal, "quantile", quantile = 0.75)
islands <- list()
chroms <- unique(as.character(bed$chrom))
for (ch in chroms) {
  ii <- which(as.character(bed$chrom) == ch)
  v <- keep[ii]
  if (!any(v)) next
  r <- rle(v)
  end <- cumsum(r$lengths)
  start <- end - r$lengths + 1L
  # POST HOC tractability: H3K27me3 islands of <15 bins (3 kb) explode
  # to >300k regions. hiddenDomains-like broad domains are longer.
  ok <- r$values & r$lengths >= 15L
  if (!any(ok)) next
  for (k in which(ok)) {
    a <- ii[start[k]]
    b <- ii[end[k]]
    islands[[length(islands) + 1L]] <- data.frame(
      chrom = ch, start = bed$start[a], end = bed$end[b],
      a = a, b = b, stringsAsFactors = FALSE
    )
  }
}
conv_bed <- do.call(rbind, islands)
message("n islands=", nrow(conv_bed))
rsig <- matrix(NA_real_, nrow(conv_bed), ncol(signal),
               dimnames = list(NULL, colnames(signal)))
for (i in seq_len(nrow(conv_bed))) {
  rsig[i, ] <- colMeans(signal[conv_bed$a[i]:conv_bed$b[i], , drop = FALSE])
}
# faster inner: rewrite with .colMeans groups of equal length if many
conv <- test_fixed_regions_limma(conv_bed[, 1:3], rsig, sheet, "KO", "Control")
conv <- standardize_primary_results(conv, meta = list(
  assay_family = "continuous", input_type = "bigwig",
  contrast = "KO - Control", inference_mode = "replicated",
  signal_filter_method = "quantile_islands_POSTHOC",
  transform_method = "log2_pseudocount", segmentation_method = "fixed_islands"
))
utils::write.table(conv, "results/v0.5/datasetB/primary_calls_conventional.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
cl <- classify_region_sets(tested, conv)
ov <- overlap_geometry_metrics(tested, conv)
ov$n_opposite_1bp <- cl$n_opposite_overlaps
utils::write.table(ov, "results/v0.5/datasetB/overlap_metrics.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
ed <- tested[tested$level1 %in% c("GAIN", "LOSS"), ]
yaml::write_yaml(list(
  dataset = "GSE110378",
  GAIN = sum(tested$level1 == "GAIN"),
  LOSS = sum(tested$level1 == "LOSS"),
  n_tested = nrow(tested),
  conv_GAIN = sum(conv$level1 == "GAIN"),
  conv_LOSS = sum(conv$level1 == "LOSS"),
  n_islands = nrow(conv_bed),
  conventional_note = "POST HOC quantile islands; hiddenDomains not deposited",
  contrast = "KO - Control",
  expected = "GAIN of H3K27me3"
), "results/v0.5/datasetB/run_manifest.yaml")
message("B ESD GAIN/LOSS=", sum(tested$level1 == "GAIN"), "/",
        sum(tested$level1 == "LOSS"),
        " conv=", sum(conv$level1 == "GAIN"), "/", sum(conv$level1 == "LOSS"))
