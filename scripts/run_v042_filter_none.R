#!/usr/bin/env Rscript
# Full GSE280947 200 bp filter=none run using chromosome-wise FPOP/PELT.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(20272)
dir.create("results/v0.4.2/filter_none/checkpoints", recursive = TRUE, showWarnings = FALSE)

rds <- "cache/gse280947_100bp_signal_matrix.rds"
cache <- readRDS(rds)
man <- yaml::read_yaml("cache/gse280947_100bp_signal_matrix_manifest.yaml")
m200 <- coarsen_bed_matrix(cache$bed, cache$signal, factor = 2L, from_bin = 100L)
sheet <- cache$samples
workers <- as.integer(Sys.getenv("EPISEGDIFF_WORKERS", "1"))

ident <- list(
  cache_id = paste(unlist(man$checksums), collapse = ","),
  genome = "mm10",
  bin_size = 200L,
  summary_fun = "mean",
  transform = "log2_pseudocount",
  pseudocount = 0.01,
  filter = "none",
  segmentation = "gaussian_pelt",
  penalty = 24,
  min_bins = 3L,
  engine = "pelt"
)

rss0 <- process_rss_mb()
t0 <- proc.time()[["elapsed"]]
fit <- analyze_continuous_matrix(
  m200$bed, m200$signal, sheet, "JNKi", "DMSO",
  bin_size = 200L, pseudocount = 0.01, filter_method = "none",
  penalty = 24, min_bins = 3L, fdr_threshold = 0.05,
  checkpoint_dir = "results/v0.4.2/filter_none/checkpoints",
  resume = TRUE, workers = workers, engine = "pelt",
  identity = ident
)
elapsed <- proc.time()[["elapsed"]] - t0
rss1 <- process_rss_mb()

tested <- fit$tested
utils::write.table(tested, "results/v0.4.2/filter_none/primary.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(fit$chromosome_runtime,
                   "results/v0.4.2/filter_none/chromosome_runtime.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

der <- tested[tested$level1 %in% c("GAIN", "LOSS"), ]
w <- if (nrow(der)) der$end - der$start else numeric()
yaml::write_yaml(list(
  software = paste("EpiSegDiff", as.character(utils::packageVersion("EpiSegDiff"))),
  config = "v0.4.2 filter=none frozen parameters",
  assay_family = "continuous",
  input_type = "bigwig",
  genome = "mm10",
  bin_size = 200L,
  summary_fun = "mean",
  transform = "log2_pseudocount",
  pseudocount = 0.01,
  filter = "none",
  contrast = "JNKi - DMSO",
  inference = "limma",
  segmentation = "gaussian_pelt",
  penalty = 24,
  min_bins = 3L,
  engine = "pelt",
  workers = workers,
  n_input_bins = fit$n_input_bins,
  n_valid = fit$n_valid,
  n_tested = fit$n_tested,
  GAIN = fit$n_GAIN,
  LOSS = fit$n_LOSS,
  total_DER = fit$n_GAIN + fit$n_LOSS,
  min_width = if (length(w)) min(w) else NA,
  median_width = if (length(w)) stats::median(w) else NA,
  p90_width = if (length(w)) as.numeric(stats::quantile(w, 0.9)) else NA,
  runtime_sec = elapsed,
  peak_rss_mb = max(c(rss0, rss1, fit$peak_rss_mb), na.rm = TRUE),
  resumed_chromosomes = sum(fit$chromosome_runtime$resumed),
  failed_chromosomes = sum(fit$chromosome_runtime$status != "PASS"),
  cache = rds
), "results/v0.4.2/filter_none/run_manifest.yaml")
message("filter=none GAIN=", fit$n_GAIN, " LOSS=", fit$n_LOSS,
        " tested=", fit$n_tested, " sec=", round(elapsed, 1))
