#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
  else library(EpiSegDiff)
})

spec <- yaml::read_yaml("config/validation/null_calibration.yaml")
set.seed(as.integer(spec$seed))
out_path <- spec$output
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
qs <- as.numeric(spec$q)
n_bins <- as.integer(spec$n_bins)
n_reg <- as.integer(spec$n_per_regime)

pipe <- default_config()
pipe$normalization$method <- "TMM"
pipe$normalization$background_definition <- 0.5
pipe$segmentation$penalties <- c(10, 20)
pipe$segmentation$selected_penalty <- 10
pipe$segmentation$engine <- "pelt"

rows <- list()
rid <- 0L
.empty_cal_row <- function() {
  data.frame(
    sim_id = NA_integer_, family = NA_character_, regime = NA_character_,
    nominal_q = NA_real_, n_discoveries = NA_integer_,
    any_discovery = NA_integer_, bp_discovered = NA_real_,
    false_boundary_count = NA_integer_, genome_width = NA_real_,
    effect = NA_real_, error = NA_character_,
    stringsAsFactors = FALSE
  )
}

push <- function(...) {
  rid <<- rid + 1L
  base <- .empty_cal_row()
  extra <- list(...)
  for (nm in names(extra)) base[[nm]] <- extra[[nm]]
  base$sim_id <- rid
  rows[[rid]] <<- base
}

summarize_qs <- function(tested, family, regime, extra = list()) {
  gw <- n_bins * 200
  for (q in qs) {
    efd <- empirical_fdr(tested, nominal = q, genome_width = gw)
    args <- utils::modifyList(list(
      family = family, regime = regime, nominal_q = q,
      n_discoveries = efd$n_discoveries[1],
      any_discovery = as.integer(efd$n_discoveries[1] > 0),
      bp_discovered = efd$bp_discovered[1],
      false_boundary_count = efd$false_boundary_count[1],
      genome_width = gw, effect = NA_real_, error = NA_character_
    ), extra)
    do.call(push, args)
  }
}

# D1 strict null
for (rg in spec$strict_null_regimes) {
  for (s in seq_len(n_reg)) {
    sim <- simulate_null_counts(
      n_bins = n_bins, n_rep = rg$n_rep, depth = rg$depth,
      dispersion = rg$dispersion, background_sd = rg$background_sd,
      seed = spec$seed + rid + s
    )
    fit <- tryCatch(run_simulated_discovery(sim, pipe), error = function(e) e)
    if (inherits(fit, "error")) {
      push(family = "strict_null", regime = rg$name, nominal_q = NA,
           n_discoveries = NA, any_discovery = NA, bp_discovered = NA,
           false_boundary_count = NA, genome_width = n_bins * 200,
           error = fit$message)
      next
    }
    summarize_qs(fit$tested, "strict_null", rg$name, list(error = NA_character_))
  }
}

# D2 near-null
for (fx in as.numeric(unlist(spec$near_null$effects))) {
  for (s in seq_len(as.integer(spec$near_null$n_sim))) {
    ev <- data.frame(event_class = "GAIN", start_bin = 80L,
                     end_bin = 80L + as.integer(spec$near_null$region_width) - 1L,
                     effect = fx)
    sim <- simulate_count_assay(
      n_bins = n_bins, n_rep = spec$near_null$n_rep, depth = spec$near_null$depth,
      dispersion = spec$near_null$dispersion, background_sd = 0.08,
      events = ev, seed = spec$seed + 8000 + rid + s
    )
    fit <- tryCatch(run_simulated_discovery(sim, pipe), error = function(e) e)
    if (inherits(fit, "error")) next
    summarize_qs(fit$tested, "near_null", sprintf("log2FC_%s", fx),
                 list(effect = fx, error = NA_character_))
  }
}

# D3 heterogeneous backgrounds (true null effect)
for (s in seq_len(as.integer(spec$heterogeneous$n_sim))) {
  kind <- c("drift", "high_bg", "coverage", "mask")[((s - 1L) %% 4L) + 1L]
  extra <- list(drift_amp = 0, high_bg_blocks = NULL, coverage_scale = NULL,
                mask_bins = NULL)
  if (kind == "drift") extra$drift_amp <- 0.8
  if (kind == "high_bg") {
    extra$high_bg_blocks <- data.frame(start_bin = 30L, end_bin = 70L, log_boost = 1.4)
  }
  if (kind == "coverage") extra$coverage_scale <- pmax(0.25, 1 + 0.6 * sin(seq_len(n_bins) / 12))
  if (kind == "mask") extra$mask_bins <- 100:115
  sim <- do.call(simulate_null_counts, c(list(
    n_bins = n_bins, n_rep = spec$heterogeneous$n_rep,
    depth = spec$heterogeneous$depth, dispersion = spec$heterogeneous$dispersion,
    background_sd = 0.12, seed = spec$seed + 9000 + s
  ), extra[names(extra) %in% c("mask_bins")]))
  if (kind == "drift" || kind == "high_bg" || kind == "coverage") {
    sim <- simulate_count_assay(
      n_bins = n_bins, n_rep = spec$heterogeneous$n_rep,
      depth = spec$heterogeneous$depth, dispersion = spec$heterogeneous$dispersion,
      background_sd = 0.12, events = NULL, seed = spec$seed + 9000 + s,
      drift_amp = extra$drift_amp,
      high_bg_blocks = extra$high_bg_blocks,
      coverage_scale = extra$coverage_scale
    )
  }
  fit <- tryCatch(run_simulated_discovery(sim, pipe), error = function(e) e)
  if (inherits(fit, "error")) next
  summarize_qs(fit$tested, "heterogeneous_null", kind, list(error = NA_character_))
}

# D5 sample-split vs same-data on a subset of strict nulls
for (s in seq_len(as.integer(spec$split_n_sim))) {
  sim <- simulate_null_counts(n_bins = n_bins, n_rep = 4, depth = 20,
                              dispersion = 0.1, background_sd = 0.1,
                              seed = spec$seed + 12000 + s)
  same <- tryCatch(run_simulated_discovery(sim, pipe), error = function(e) NULL)
  spl <- tryCatch(sample_split_discovery(sim, pipe), error = function(e) NULL)
  if (!is.null(same)) {
    summarize_qs(same$tested, "selection_same_data", "null_4rep",
                 list(error = NA_character_))
  }
  if (!is.null(spl)) {
    summarize_qs(spl$tested, "selection_sample_split", "null_4rep",
                 list(error = NA_character_))
  }
}

tab <- do.call(rbind, rows)
utils::write.table(tab, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", out_path, " n=", nrow(tab))
