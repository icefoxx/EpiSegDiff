#!/usr/bin/env Rscript
# Empirical null FDR grid (Phase 6). Writes reports/Phase2_FDR_Calibration_Report.md
# material under benchmark/simulation/.

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(dirname(dirname(normalizePath(commandArgs(trailingOnly = FALSE)[1]))))
  } else {
    library(EpiSegDiff)
  }
})

# When sourced after load_all:
if (!"simulate_null_counts" %in% getNamespaceExports("EpiSegDiff") &&
    exists("simulate_null_counts", mode = "function", inherits = TRUE)) {
  # in-source
}

args_seed <- 1L
grid <- expand.grid(
  n_rep = c(2L, 3L, 5L),
  dispersion = c(0.05, 0.15, 0.40),
  depth = c(8, 20, 50),
  KEEP.OUT.ATTRS = FALSE
)
names(grid) <- c("n_rep", "dispersion", "depth")
grid$disp_label <- ifelse(grid$dispersion < 0.1, "low",
                          ifelse(grid$dispersion < 0.25, "medium", "high"))
grid$depth_label <- ifelse(grid$depth < 12, "low",
                           ifelse(grid$depth < 35, "medium", "high"))

cfg <- EpiSegDiff:::default_config()
cfg$normalization$method <- "library_size"
cfg$segmentation$penalties <- c(10, 20)
cfg$segmentation$selected_penalty <- 20
cfg$segmentation$engine <- "pelt"
cfg$segmentation$min_width <- 3L

rows <- list()
for (i in seq_len(nrow(grid))) {
  sim <- EpiSegDiff:::simulate_null_counts(
    n_bins = 200L, n_rep = grid$n_rep[i], depth = grid$depth[i],
    dispersion = grid$dispersion[i], background_sd = 0.1,
    seed = args_seed + i
  )
  fit <- tryCatch(
    EpiSegDiff:::run_simulated_discovery(sim, cfg),
    error = function(e) list(error = e$message, tested = EpiSegDiff:::.empty_test_table())
  )
  if (!is.null(fit$error)) {
    rows[[i]] <- data.frame(grid[i, ], error = fit$error, stringsAsFactors = FALSE)
    next
  }
  efd <- EpiSegDiff:::empirical_fdr(fit$tested, nominal = c(0.01, 0.05, 0.10),
                                    genome_width = 200 * 200)
  efd$n_rep <- grid$n_rep[i]
  efd$dispersion <- grid$dispersion[i]
  efd$depth <- grid$depth[i]
  efd$disp_label <- grid$disp_label[i]
  efd$depth_label <- grid$depth_label[i]
  efd$runtime_sec <- fit$runtime_sec
  efd$error <- NA_character_
  rows[[i]] <- efd
}

out <- do.call(rbind, rows)
dir.create("benchmark/simulation", recursive = TRUE, showWarnings = FALSE)
utils::write.table(out, "benchmark/simulation/null_fdr_grid.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote benchmark/simulation/null_fdr_grid.tsv")
