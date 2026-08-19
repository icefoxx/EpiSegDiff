#!/usr/bin/env Rscript
# Complete v0.4.1 artefacts after filter=none was killed (C++ PELT).
# Does not overwrite the 12-cell count table unless missing.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(20271)
dir.create("results/v0.4.1", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

cache <- readRDS("cache/gse280947_100bp_signal_matrix.rds")
sheet <- cache$samples
num <- "JNKi"; den <- "DMSO"
strong_abs <- 1
penalty <- 24

m100 <- list(bed = cache$bed, signal = cache$signal)
m200 <- coarsen_bed_matrix(cache$bed, cache$signal, factor = 2L, from_bin = 100L)
m500 <- coarsen_bed_matrix(cache$bed, cache$signal, factor = 5L, from_bin = 100L)
mats <- list(`100` = m100, `200` = m200, `500` = m500)
message("coarsened")

ref_tab <- utils::read.table("results/v0.4.1/reference_primary.tsv",
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
conc <- list()
for (bs in c(100, 200, 500)) {
  for (pc in c(0.001, 0.01, 0.05, 0.1)) {
    key <- sprintf("b%s_p%s", bs, pc)
    if (bs == 200 && pc == 0.01) next
    message("recompute ", key, " for concordance")
    fit <- analyze_continuous_matrix(
      mats[[as.character(bs)]]$bed, mats[[as.character(bs)]]$signal,
      sheet, num, den, bin_size = bs, pseudocount = pc,
      filter_method = "absolute", filter_threshold = 0.6,
      penalty = penalty, min_bins = 3L, fdr_threshold = 0.05
    )
    conc[[length(conc) + 1L]] <- compare_continuous_runs(
      ref_tab, fit$tested, "ref200", key, strong_abs
    )
    utils::write.table(do.call(rbind, conc), "results/v0.4.1/run_concordance.tsv",
                       sep = "\t", quote = FALSE, row.names = FALSE)
    rm(fit)
    gc(verbose = FALSE)
  }
}

# Record filter=none as NOT_RUN with the measured reason
grid <- utils::read.table("results/v0.4.1/combined_sensitivity_grid.tsv",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ref_row <- grid[grid$config_id == "b200_p0.01_absolute_pen24", ]
none_row <- data.frame(
  config_id = "b200_p0.01_none_pen24",
  bin_size = 200, pseudocount = 0.01, summary_fun = "mean",
  filter_method = "none", filter_parameter = NA_real_,
  tested_segments = NA, GAIN = NA, LOSS = NA, total_DER = NA,
  median_width = NA, p90_width = NA, median_abs_effect = NA,
  n_valid_bins = ref_row$n_input_bins[1],
  n_input_bins = ref_row$n_input_bins[1],
  runtime_sec = NA,
  status = "NOT_RUN",
  failure_reason = paste(
    "Killed after >30 min at ~100% CPU in Gaussian PELT on unfiltered",
    "genome-wide z-scores (13.6M bins; one CVR per chromosome).",
    "setTimeLimit cannot interrupt pelt_gaussian_cpp.",
    "Absolute 0.6 retains 91390 / 13627618 bins (0.67%)."
  ),
  stringsAsFactors = FALSE
)
utils::write.table(rbind(none_row, ref_row), "results/v0.4.1/filter_sensitivity.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# Limma diagnostics on 200 bp linear matrix
xf <- transform_continuous(m200$signal, "log2_pseudocount", 0.01)
vd <- limma_variance_diagnostics(xf, sheet, num, den, signal_linear = m200$signal)
utils::write.table(vd, "results/v0.4.1/limma_variance_diagnostics.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
message("limma class=", vd$classification)

bin_tab <- grid[grid$pseudocount == 0.01, ]
pseudo_tab <- grid[grid$bin_size == 200, ]
filt_tab <- rbind(none_row, ref_row)
cc <- do.call(rbind, conc)

pdf("figures/v0.4.1_gain_loss_by_configuration.pdf", width = 9, height = 4)
graphics::barplot(rbind(grid$GAIN, grid$LOSS), beside = TRUE,
                  names.arg = grid$config_id, las = 2, cex.names = 0.45,
                  legend.text = c("GAIN", "LOSS"), main = "GAIN/LOSS by config")
dev.off()
pdf("figures/v0.4.1_width_by_bin_size.pdf", width = 6, height = 4)
graphics::barplot(bin_tab$median_width, names.arg = bin_tab$bin_size,
                  ylab = "median DER width", main = "Width by bin size")
dev.off()
pdf("figures/v0.4.1_der_jaccard_heatmap.pdf", width = 7, height = 4)
graphics::plot(seq_len(nrow(cc)), cc$DER_bp_jaccard, ylim = c(0, 1),
               xlab = "alternative", ylab = "bp Jaccard vs ref",
               main = "DER Jaccard", pch = 16)
graphics::text(seq_len(nrow(cc)), cc$DER_bp_jaccard, cc$alternative_id, pos = 3, cex = 0.5)
dev.off()
pdf("figures/v0.4.1_strong_region_recovery.pdf", width = 7, height = 4)
graphics::plot(seq_len(nrow(cc)), cc$strong_reference_DER_recovery, ylim = c(0, 1),
               xlab = "alternative", ylab = "strong ref recovery", pch = 16)
dev.off()
pdf("figures/v0.4.1_effect_correlation_heatmap.pdf", width = 7, height = 4)
graphics::plot(seq_len(nrow(cc)), cc$effect_pearson, ylim = c(-1, 1),
               xlab = "alternative", ylab = "effect Pearson", pch = 16)
dev.off()
pdf("figures/v0.4.1_filter_comparison.pdf", width = 5, height = 4)
vals <- rbind(c(none_row$GAIN, ref_row$GAIN), c(none_row$LOSS, ref_row$LOSS))
vals[is.na(vals)] <- 0
graphics::barplot(vals, beside = TRUE, names.arg = c("none (NOT_RUN)", "absolute 0.6"),
                  legend.text = c("GAIN", "LOSS"), main = "Filter sensitivity")
dev.off()
pdf("figures/v0.4.1_pseudocount_low_signal_sensitivity.pdf", width = 6, height = 4)
graphics::plot(pseudo_tab$pseudocount, pseudo_tab$LOSS, log = "x",
               xlab = "pseudocount", ylab = "LOSS", pch = 16, main = "LOSS vs pseudocount")
dev.off()

source("scripts/write_v041_reports.R")
message("finish_v041 done")
