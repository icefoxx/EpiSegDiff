#!/usr/bin/env Rscript
# Small in-package simulation benchmark (Phase 9). Writes
# benchmark/simulation/truth.tsv and results.tsv.

devtools::load_all(".", quiet = TRUE)

set.seed(1)
classes <- c("GAIN", "LOSS", "NO_CHANGE")
truth_rows <- list()
res_rows <- list()
cfg <- default_config()
cfg$normalization$method <- "library_size"
cfg$segmentation$penalties <- c(8, 16)
cfg$segmentation$selected_penalty <- 8
cfg$segmentation$engine <- "pelt"

k <- 0L
for (cl in classes) {
  for (rep in 1:5) {
    k <- k + 1L
    ev <- example_topology_event(cl, offset = 80L)
    sim <- simulate_count_assay(
      n_bins = 200L, n_rep = 3L, depth = 25, dispersion = 0.08,
      background_sd = 0.08, events = ev, seed = 100 + k
    )
    tr <- sim$truth
    tr$sim_id <- k
    truth_rows[[k]] <- tr
    fit <- run_simulated_discovery(sim, cfg)
    pred <- fit$tested
    if (nrow(pred) && "FDR" %in% names(pred)) {
      pred <- pred[!is.na(pred$FDR) & pred$FDR <= 0.05, ]
    }
    truth_iv <- tr[!tr$event_class %in% c("NO_CHANGE", "NULL"), c("chrom", "start", "end"), drop = FALSE]
    if (!nrow(pred)) {
      pred <- data.frame(chrom = character(), start = integer(), end = integer())
    }
    met <- compute_region_metrics(pred, truth_iv)
    bmet <- compute_boundary_metrics(pred, truth_iv, met$matches)
    res_rows[[k]] <- data.frame(
      sim_id = k, event_class = cl,
      precision = met$precision, recall = met$recall, f1 = met$f1,
      mean_iou = met$mean_iou, bp_precision = met$bp_precision,
      bp_recall = met$bp_recall, start_mae = bmet$start_mae,
      end_mae = bmet$end_mae, hausdorff = bmet$hausdorff,
      runtime_sec = fit$runtime_sec,
      stringsAsFactors = FALSE
    )
  }
}

dir.create("benchmark/simulation", recursive = TRUE, showWarnings = FALSE)
utils::write.table(do.call(rbind, truth_rows), "benchmark/simulation/truth.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(do.call(rbind, res_rows), "benchmark/simulation/results.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote benchmark/simulation/{truth,results}.tsv")
