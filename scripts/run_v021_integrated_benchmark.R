#!/usr/bin/env Rscript
# Task E — run only after A–D reports exist.

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
  else library(EpiSegDiff)
})

spec <- yaml::read_yaml("config/validation/integrated_benchmark.yaml")
set.seed(as.integer(spec$seed))
out_path <- spec$output
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
classes <- unlist(spec$classes)
methods <- unlist(spec$methods)
n_sim <- as.integer(spec$n_sim)
n_bins <- as.integer(spec$n_bins)

pipe <- default_config()
pipe$normalization$background_definition <- 0.5
pipe$segmentation$penalties <- c(8, 16)
pipe$segmentation$selected_penalty <- 8
pipe$segmentation$engine <- "pelt"

rows <- list()
k <- 0L
for (cl in classes) {
  for (sc in unlist(spec$scales)) {
    for (nm in methods) {
      for (s in seq_len(n_sim)) {
        off <- if (identical(sc, "narrow")) 60L else 40L
        ev <- example_topology_event(cl, offset = off)
        if (identical(sc, "broad") && cl %in% c("GAIN", "LOSS", "EXPANSION",
                                               "CONTRACTION", "LEFT_SHIFT",
                                               "RIGHT_SHIFT")) {
          ev$end_bin <- ev$start_bin + 50L
          if ("control_end_bin" %in% names(ev) && !is.na(ev$control_end_bin[1])) {
            ev$control_end_bin <- ev$control_start_bin + 40L
            ev$treatment_end_bin <- ev$treatment_start_bin + 55L
          }
        }
        sim <- simulate_count_assay(
          n_bins = n_bins, n_rep = 3L, depth = 22, dispersion = 0.1,
          background_sd = 0.07, events = ev, seed = spec$seed + k
        )
        pipe$normalization$method <- nm
        fit <- tryCatch(run_simulated_discovery(sim, pipe), error = function(e) e)
        k <- k + 1L
        if (inherits(fit, "error")) {
          rows[[k]] <- data.frame(sim_id = k, truth = cl, scale = sc, method = nm,
                                  error = fit$message, stringsAsFactors = FALSE)
          next
        }
        disc <- summarize_discovery_vs_truth(fit$tested, sim$truth, q = 0.05)
        # local domains from counts, not MACS
        C <- build_condition_domains(sim$counts, sim$bins, sim$samples, "control",
                                     method = "condition_mean_threshold",
                                     mean_quantile = 0.85,
                                     min_width = 3L)
        Tm <- build_condition_domains(sim$counts, sim$bins, sim$samples, "treatment",
                                      method = "condition_mean_threshold",
                                      mean_quantile = 0.85,
                                      min_width = 3L)
        sig <- if (!is.null(fit$tested) && nrow(fit$tested)) {
          x <- fit$tested
          if (!"direction" %in% names(x) && "logFC" %in% names(x)) {
            x$direction <- ifelse(x$logFC > 0, "GAIN", "LOSS")
          }
          x[!is.na(x$FDR) & x$FDR <= 0.05, ]
        } else NULL
        evs <- classify_boundary_events(
          match_domains(C, Tm), C, Tm, differential = sig
        )
        pred_cls <- if (nrow(evs)) evs$primary_event[1] else "NONE"
        rows[[k]] <- data.frame(
          sim_id = k, truth = cl, scale = sc, method = nm,
          precision = disc$precision, recall = disc$recall, f1 = disc$f1,
          mean_iou = disc$mean_iou, bp_precision = disc$bp_precision,
          bp_recall = disc$bp_recall, start_mae = disc$start_mae,
          end_mae = disc$end_mae, hausdorff = disc$hausdorff,
          pred_event = pred_cls,
          n_events = nrow(evs),
          complex_flag = as.integer(any(evs$primary_event == "COMPLEX")),
          peak_caller = "none",
          error = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
tab <- do.call(rbind, rows)
utils::write.table(tab, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", out_path, " n=", nrow(tab))
