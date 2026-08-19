#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(20271)
dir.create("results/v0.4.1", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

rds <- "cache/gse280947_100bp_signal_matrix.rds"
if (!file.exists(rds)) {
  stop("Build cache first: Rscript scripts/build_gse280947_100bp_cache.R")
}
cache <- readRDS(rds)
sheet <- cache$samples
num <- "JNKi"; den <- "DMSO"
strong_abs <- 1
penalty <- 24
filter_none_limit_sec <- as.numeric(Sys.getenv("EPISEGDIFF_FILTER_NONE_LIMIT", "7200"))

empty_row <- function(bin_size, pseudo, filter, fpar, status, reason) {
  data.frame(
    config_id = sprintf("b%s_p%s_%s_pen%s", bin_size, pseudo, filter, penalty),
    bin_size = bin_size, pseudocount = pseudo, summary_fun = "mean",
    filter_method = filter, filter_parameter = fpar,
    tested_segments = NA_integer_, GAIN = NA_integer_, LOSS = NA_integer_,
    total_DER = NA_integer_,
    median_width = NA_real_, p90_width = NA_real_, median_abs_effect = NA_real_,
    n_valid_bins = NA_integer_, n_input_bins = NA_integer_,
    runtime_sec = NA_real_, status = status, failure_reason = reason,
    stringsAsFactors = FALSE
  )
}

summarize_run <- function(fit, bin_size, pseudo, filter, fpar, status = "PASS",
                          reason = "") {
  sig <- fit$tested
  der <- sig[sig$level1 %in% c("GAIN", "LOSS"), ]
  w <- if (nrow(der)) der$end - der$start else numeric()
  data.frame(
    config_id = sprintf("b%s_p%s_%s_pen%s", bin_size, pseudo, filter, penalty),
    bin_size = bin_size, pseudocount = pseudo, summary_fun = "mean",
    filter_method = filter, filter_parameter = fpar,
    tested_segments = fit$n_tested,
    GAIN = fit$n_GAIN, LOSS = fit$n_LOSS,
    total_DER = fit$n_GAIN + fit$n_LOSS,
    median_width = if (length(w)) stats::median(w) else NA_real_,
    p90_width = if (length(w)) as.numeric(stats::quantile(w, 0.9)) else NA_real_,
    median_abs_effect = if (nrow(der)) stats::median(abs(der$effect), na.rm = TRUE) else NA_real_,
    n_valid_bins = fit$n_valid, n_input_bins = fit$n_input_bins,
    runtime_sec = fit$runtime_sec, status = status, failure_reason = reason,
    stringsAsFactors = FALSE
  )
}

run_one <- function(mat, bin_size, pseudo, filter, fpar, time_limit = Inf) {
  t0 <- proc.time()[["elapsed"]]
  tryCatch({
    if (is.finite(time_limit)) setTimeLimit(elapsed = time_limit, transient = TRUE)
    fit <- analyze_continuous_matrix(
      mat$bed, mat$signal, sheet, num, den, bin_size = bin_size,
      pseudocount = pseudo, filter_method = filter,
      filter_threshold = if (filter == "absolute") fpar else 0.6,
      filter_quantile = if (filter == "quantile") fpar else 0.75,
      penalty = penalty, min_bins = 3L, fdr_threshold = 0.05
    )
    setTimeLimit(elapsed = Inf)
    list(ok = TRUE, fit = fit, row = summarize_run(fit, bin_size, pseudo, filter, fpar))
  }, error = function(e) {
    setTimeLimit(elapsed = Inf)
    msg <- conditionMessage(e)
    st <- if (grepl("time.?limit|reached elapsed", msg, ignore.case = TRUE)) "NOT_RUN" else "FAILED"
    if (identical(st, "NOT_RUN")) {
      msg <- paste0("elapsed time limit ", time_limit, "s: ", msg)
    }
    list(ok = FALSE, fit = NULL,
         row = empty_row(bin_size, pseudo, filter, fpar, st, msg))
  })
}

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

message("Coarsen 100 -> 200 and 500 (linear scale)")
m100 <- list(bed = cache$bed, signal = cache$signal)
m200 <- coarsen_bed_matrix(cache$bed, cache$signal, factor = 2L, from_bin = 100L)
m500 <- coarsen_bed_matrix(cache$bed, cache$signal, factor = 5L, from_bin = 100L)
mats <- list(`100` = m100, `200` = m200, `500` = m500)
message("n bins 100/200/500 = ", nrow(m100$bed), "/", nrow(m200$bed), "/", nrow(m500$bed))

# 12-cell grid + reference (filter absolute 0.6)
grid_rows <- list()
fits <- list()
conc <- list()
gi <- 0L
for (bs in c(100, 200, 500)) {
  for (pc in c(0.001, 0.01, 0.05, 0.1)) {
    key <- sprintf("b%s_p%s", bs, pc)
    message("grid ", key)
    r <- run_one(mats[[as.character(bs)]], bs, pc, "absolute", 0.6)
    gi <- gi + 1L
    grid_rows[[gi]] <- r$row
    if (r$ok) fits[[key]] <- r$fit
    write_tsv(do.call(rbind, grid_rows), "results/v0.4.1/combined_sensitivity_grid.tsv")
  }
}
grid_tab <- do.call(rbind, grid_rows)
write_tsv(grid_tab, "results/v0.4.1/combined_sensitivity_grid.tsv")

ref_key <- "b200_p0.01"
if (is.null(fits[[ref_key]])) {
  stop("Reference cell 200bp / 0.01 failed: ",
       grid_tab$failure_reason[grid_tab$config_id == "b200_p0.01_absolute_pen24"])
}
ref_fit <- fits[[ref_key]]
ref_tab <- ref_fit$tested
write_tsv(ref_tab, "results/v0.4.1/reference_primary.tsv")
min_w <- if (nrow(ref_tab)) min(ref_tab$end - ref_tab$start) else NA_real_
yaml::write_yaml(list(
  config = "config/validation/v0.4.1_gse280947_reference.yaml",
  GAIN = ref_fit$n_GAIN, LOSS = ref_fit$n_LOSS,
  min_width = min_w,
  n_tested = ref_fit$n_tested, runtime_sec = ref_fit$runtime_sec,
  source = "cached_100bp_coarsened_to_200",
  n_input_bins = ref_fit$n_input_bins,
  n_valid_bins = ref_fit$n_valid,
  strong_abs_effect = strong_abs,
  contrast = "JNKi - DMSO",
  summary_fun = "mean",
  filter = "absolute 0.6",
  v04_expected = list(GAIN = 2, LOSS = 688, min_width = 600)
), "results/v0.4.1/reference_run_manifest.yaml")
message("Reference GAIN=", ref_fit$n_GAIN, " LOSS=", ref_fit$n_LOSS,
        " min_width=", min_w)

bin_tab <- grid_tab[grid_tab$pseudocount == 0.01, ]
write_tsv(bin_tab, "results/v0.4.1/bin_size_sensitivity.tsv")
pseudo_tab <- grid_tab[grid_tab$bin_size == 200, ]
write_tsv(pseudo_tab, "results/v0.4.1/pseudocount_sensitivity.tsv")

for (key in names(fits)) {
  if (identical(key, ref_key)) next
  conc[[length(conc) + 1L]] <- compare_continuous_runs(
    ref_tab, fits[[key]]$tested, "ref200", key, strong_abs
  )
}

# Pseudocount x signal-quartile on 200 bp bin-level effects
mx200 <- do.call(pmax, c(as.list(as.data.frame(m200$signal)), list(na.rm = TRUE)))
qcut <- tryCatch(as.integer(cut(mx200, breaks = stats::quantile(mx200, 0:4 / 4, na.rm = TRUE),
                                include.lowest = TRUE, labels = FALSE)),
                 error = function(e) rep(NA_integer_, length(mx200)))
ref_eff <- fits[[ref_key]]$bin_effect$effect
q_rows <- list()
for (pc in c(0.001, 0.01, 0.05, 0.1)) {
  key <- sprintf("b200_p%s", pc)
  if (is.null(fits[[key]])) next
  alt_eff <- fits[[key]]$bin_effect$effect
  for (q in 1:4) {
    keep <- qcut == q & is.finite(ref_eff) & is.finite(alt_eff)
    if (!any(keep)) next
    q_rows[[length(q_rows) + 1L]] <- data.frame(
      pseudocount = pc, signal_quartile = q, n = sum(keep),
      effect_pearson = suppressWarnings(stats::cor(ref_eff[keep], alt_eff[keep])),
      median_abs_effect_diff = stats::median(abs(ref_eff[keep] - alt_eff[keep])),
      sign_flip_rate = mean(sign(ref_eff[keep]) != sign(alt_eff[keep]) &
                              ref_eff[keep] != 0 & alt_eff[keep] != 0),
      stringsAsFactors = FALSE
    )
  }
}
if (length(q_rows)) {
  write_tsv(do.call(rbind, q_rows), "results/v0.4.1/pseudocount_signal_quartile.tsv")
}

# Filter=none last so a hang cannot drop the 12-cell grid
message("filter none (genome-wide 200 bp), limit ", filter_none_limit_sec, "s")
filt_none <- empty_row(200, 0.01, "none", NA, "NOT_RUN",
                       "pending genome-wide execution")
write_tsv(rbind(filt_none, grid_tab[grid_tab$config_id == "b200_p0.01_absolute_pen24", ]),
          "results/v0.4.1/filter_sensitivity.tsv")
r0 <- run_one(m200, 200, 0.01, "none", NA, time_limit = filter_none_limit_sec)
filt_tab <- rbind(r0$row, grid_tab[grid_tab$config_id == "b200_p0.01_absolute_pen24", ])
write_tsv(filt_tab, "results/v0.4.1/filter_sensitivity.tsv")
if (r0$ok) {
  conc[[length(conc) + 1L]] <- compare_continuous_runs(
    ref_tab, r0$fit$tested, "ref200", "filter_none", strong_abs
  )
}

if (length(conc)) {
  write_tsv(do.call(rbind, conc), "results/v0.4.1/run_concordance.tsv")
}

# Limma variance on reference transformed 200
xf <- transform_continuous(m200$signal, "log2_pseudocount", 0.01)
valid200 <- filter_signal_bins(m200$signal, method = "absolute", threshold = 0.6)
# bin-level FDR is not stored; use |effect| of tested DERs projected by overlap later
vd <- limma_variance_diagnostics(xf, sheet, num, den, signal_linear = m200$signal)
write_tsv(vd, "results/v0.4.1/limma_variance_diagnostics.tsv")
message("limma classification: ", vd$classification)

# Plots
pdf("figures/v0.4.1_gain_loss_by_configuration.pdf", width = 9, height = 4)
graphics::barplot(rbind(grid_tab$GAIN, grid_tab$LOSS), beside = TRUE,
                  names.arg = grid_tab$config_id, las = 2, cex.names = 0.45,
                  legend.text = c("GAIN", "LOSS"), main = "GAIN/LOSS by config")
dev.off()
pdf("figures/v0.4.1_width_by_bin_size.pdf", width = 6, height = 4)
graphics::barplot(bin_tab$median_width, names.arg = bin_tab$bin_size,
                  ylab = "median DER width", main = "Width by bin size")
dev.off()
if (file.exists("results/v0.4.1/run_concordance.tsv")) {
  cc <- utils::read.table("results/v0.4.1/run_concordance.tsv", header = TRUE, sep = "\t")
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
}
pdf("figures/v0.4.1_filter_comparison.pdf", width = 5, height = 4)
graphics::barplot(rbind(filt_tab$GAIN, filt_tab$LOSS), beside = TRUE,
                  names.arg = filt_tab$filter_method, legend.text = c("GAIN", "LOSS"),
                  main = "Filter sensitivity")
dev.off()
pdf("figures/v0.4.1_pseudocount_low_signal_sensitivity.pdf", width = 6, height = 4)
graphics::plot(pseudo_tab$pseudocount, pseudo_tab$LOSS, log = "x",
               xlab = "pseudocount", ylab = "LOSS", pch = 16, main = "LOSS vs pseudocount")
dev.off()

invisible(valid200)
message("v0.4.1 validation tables written")
