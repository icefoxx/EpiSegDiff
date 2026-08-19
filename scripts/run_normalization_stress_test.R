#!/usr/bin/env Rscript
# Task A — normalization composition-bias stress test (reduced balanced design).

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
  else library(EpiSegDiff)
})

cfg_path <- "config/validation/normalization_stress.yaml"
spec <- yaml::read_yaml(cfg_path)
set.seed(as.integer(spec$seed))
out_path <- spec$output
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

methods <- unlist(spec$methods)
fractions <- as.numeric(spec$fractions)
gain_props <- as.numeric(spec$gain_props)
n_bins <- as.integer(spec$n_bins)
P <- spec$primary
n_sim <- as.integer(spec$n_sim_primary)
disc_fracs <- as.numeric(spec$run_discovery_fractions)

pipe <- default_config()
pipe$normalization$background_definition <- 0.5
pipe$segmentation$penalties <- c(10, 20)
pipe$segmentation$selected_penalty <- 10
pipe$segmentation$engine <- "pelt"
pipe$normalization$method <- "library_size"

rows <- list()
rid <- 0L
add_row <- function(...) {
  rid <<- rid + 1L
  rows[[rid]] <<- data.frame(sim_id = rid, ..., stringsAsFactors = FALSE)
}

run_cell <- function(frac, gp, method, n_rep, depth, disp, effect, width, nrep, tag) {
  for (s in seq_len(nrep)) {
    ev <- compose_fraction_events(n_bins, frac, gp, width, effect)
    sim <- simulate_count_assay(
      n_bins = n_bins, n_rep = n_rep, depth = depth, dispersion = disp,
      background_sd = 0.08, events = ev, seed = spec$seed + rid + s
    )
    actual_frac <- length(sim$event_idx) / n_bins
    do_disc <- frac %in% disc_fracs && abs(gp - 0.5) < 1e-9 ||
      frac %in% disc_fracs && abs(gp - 1) < 1e-9
    fit <- tryCatch(
      fit_normalized_effect(sim, method, run_discovery = do_disc, config = pipe),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      add_row(regime = tag, method = method, frac = frac, gain_prop = gp,
              n_rep = n_rep, depth = depth, dispersion = disp, effect = effect,
              region_width = width, actual_frac = actual_frac,
              bias_background = NA, mae_background = NA, bias_z = NA,
              n_pred = NA, precision = NA, recall = NA, f1 = NA,
              n_false_regions = NA, start_mae = NA, end_mae = NA,
              false_boundary_count = NA, error = fit$message)
      next
    }
    disc <- if (!is.null(fit$tested)) {
      summarize_discovery_vs_truth(fit$tested, sim$truth, q = 0.05,
                                   bin_size = 200L)
    } else NULL
    add_row(regime = tag, method = method, frac = frac, gain_prop = gp,
            n_rep = n_rep, depth = depth, dispersion = disp, effect = effect,
            region_width = width, actual_frac = actual_frac,
            bias_background = fit$background$bias_background,
            mae_background = fit$background$mae_background,
            bias_z = fit$background$bias_z,
            n_pred = if (is.null(disc)) NA else disc$n_pred,
            precision = if (is.null(disc)) NA else disc$precision,
            recall = if (is.null(disc)) NA else disc$recall,
            f1 = if (is.null(disc)) NA else disc$f1,
            n_false_regions = if (is.null(disc)) NA else disc$n_false_regions,
            start_mae = if (is.null(disc)) NA else disc$start_mae,
            end_mae = if (is.null(disc)) NA else disc$end_mae,
            false_boundary_count = if (is.null(disc)) NA else disc$false_boundary_count,
            error = NA_character_)
  }
}

# Primary grid
for (frac in fractions) {
  gps <- if (frac == 0) 1 else gain_props
  for (gp in gps) {
    for (method in methods) {
      run_cell(frac, gp, method, P$n_rep, P$depth, P$dispersion, P$effect,
               P$region_width, n_sim, "primary")
    }
  }
}

# Satellites around frac=0.10, 50:50
A <- spec$satellite_anchor
ns <- as.integer(spec$n_sim_satellite)
for (method in methods) {
  for (nr in as.integer(unlist(A$n_rep))) {
    run_cell(A$fraction, A$gain_prop, method, nr, P$depth, P$dispersion,
             P$effect, P$region_width, ns, "sat_nrep")
  }
  for (dp in as.numeric(unlist(A$depth))) {
    run_cell(A$fraction, A$gain_prop, method, P$n_rep, dp, P$dispersion,
             P$effect, P$region_width, ns, "sat_depth")
  }
  for (ds in as.numeric(unlist(A$dispersion))) {
    run_cell(A$fraction, A$gain_prop, method, P$n_rep, P$depth, ds,
             P$effect, P$region_width, ns, "sat_disp")
  }
  for (fx in as.numeric(unlist(A$effect))) {
    run_cell(A$fraction, A$gain_prop, method, P$n_rep, P$depth, P$dispersion,
             fx, P$region_width, ns, "sat_effect")
  }
  for (rw in as.integer(unlist(A$region_width))) {
    run_cell(A$fraction, A$gain_prop, method, P$n_rep, P$depth, P$dispersion,
             P$effect, rw, ns, "sat_width")
  }
}

tab <- do.call(rbind, rows)
utils::write.table(tab, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", out_path, " n=", nrow(tab))
