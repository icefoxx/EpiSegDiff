#!/usr/bin/env Rscript
# Task C — local condition boundaries without external peak callers.

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
  else library(EpiSegDiff)
})

set.seed(20265)
out_path <- "benchmark/simulation/boundary_validation.tsv"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

scales <- list(
  narrow = list(width = 8L, n_bins = 200L, bin_size = 50L),
  broad = list(width = 40L, n_bins = 300L, bin_size = 200L)
)
effects <- c(1.0, 1.8)
n_rep <- c(3L, 5L)
depths <- c(12, 30)
methods_norm <- c("TMM", "background", "library_size")
n_sim <- 6L
pipe <- default_config()
pipe$normalization$background_definition <- 0.5

rows <- list()
k <- 0L
for (sc in names(scales)) {
  for (fx in effects) {
    for (nr in n_rep) {
      for (dp in depths) {
        for (nm in methods_norm) {
          for (s in seq_len(n_sim)) {
            w <- scales[[sc]]$width
            ev <- data.frame(
              event_class = "GAIN", start_bin = 40L, end_bin = 40L + w - 1L,
              effect = fx,
              control_start_bin = 40L, control_end_bin = 40L + w - 1L,
              treatment_start_bin = 40L, treatment_end_bin = 40L + w - 1L
            )
            sim <- simulate_count_assay(
              n_bins = scales[[sc]]$n_bins, bin_size = scales[[sc]]$bin_size,
              n_rep = nr, depth = dp, dispersion = 0.08, background_sd = 0.05,
              events = ev, seed = 4000 + k
            )
            region <- data.frame(
              chrom = "chr1",
              start = sim$bins$start[40],
              end = sim$bins$end[40 + w - 1L]
            )
            # coarse window = candidate ± 2 widths (not a peak call)
            b <- tryCatch(
              infer_local_condition_boundaries(
                sim$counts, sim$bins, sim$samples, region,
                window_bins = max(8L, 2L * w),
                method = "end_walk",
                size_factors = normalize_counts(
                  sim$counts, method = if (nm == "background") "background" else nm,
                  background_definition = 0.5
                )$size_factors
              ),
              error = function(e) e
            )
            k <- k + 1L
            if (inherits(b, "error")) {
              rows[[k]] <- data.frame(
                sim_id = k, scale = sc, effect = fx, n_rep = nr, depth = dp,
                method = nm, error = b$message, stringsAsFactors = FALSE
              )
              next
            }
            tr <- b[b$condition == "treatment", ]
            truth_L <- region$start
            truth_R <- region$end
            dL <- if (nrow(tr)) abs(tr$L[1] - truth_L) else NA
            dR <- if (nrow(tr)) abs(tr$R[1] - truth_R) else NA
            rows[[k]] <- data.frame(
              sim_id = k, scale = sc, effect = fx, n_rep = nr, depth = dp,
              method = nm, truth_L = truth_L, truth_R = truth_R,
              pred_L = if (nrow(tr)) tr$L[1] else NA, pred_R = if (nrow(tr)) tr$R[1] else NA,
              left_abs_err = dL, right_abs_err = dR,
              peak_caller = "none",
              boundary_method = "end_walk",
              error = NA_character_,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
}
tab <- do.call(rbind, rows)
utils::write.table(tab, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
ok <- tab[!is.na(tab$left_abs_err), ]
message("Wrote ", out_path, " n=", nrow(tab),
        " median_left=", median(ok$left_abs_err),
        " p90_left=", quantile(ok$left_abs_err, 0.9),
        " median_right=", median(ok$right_abs_err))
