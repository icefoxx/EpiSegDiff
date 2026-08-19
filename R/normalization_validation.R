#' Plant GAIN/LOSS blocks covering a target genome fraction
#'
#' Used by the v0.2.1 composition-bias stress test. Size factors are
#' *not* assigned arbitrarily; composition changes come only from
#' planted counts.
#'
#' @param n_bins Total bins.
#' @param frac Fraction of bins that are differential (0–1).
#' @param gain_prop Fraction of differential bins that are GAIN (rest LOSS).
#' @param region_width Width of each planted block in bins.
#' @param effect Log-fold planted on GAIN (LOSS uses `-effect`).
#' @return An events table for [simulate_count_assay()].
#' @export
compose_fraction_events <- function(n_bins, frac, gain_prop = 1,
                                    region_width = 20L, effect = 1.5) {
  n_bins <- as.integer(n_bins)
  region_width <- as.integer(region_width)
  target <- as.integer(round(frac * n_bins))
  if (target <= 0L || frac <= 0) {
    return(data.frame(event_class = "NO_CHANGE", start_bin = 1L, end_bin = 1L,
                      effect = 0, stringsAsFactors = FALSE))
  }
  n_gain <- as.integer(round(target * gain_prop))
  n_loss <- max(0L, target - n_gain)
  cursor <- 5L
  ev <- list()
  add_blocks <- function(n_need, cls, fx) {
    left <- n_need
    while (left > 0L && cursor + 2L < n_bins - 5L) {
      w <- min(region_width, left, n_bins - cursor - 5L)
      if (w < 2L) break
      ev[[length(ev) + 1L]] <<- data.frame(
        event_class = cls, start_bin = cursor, end_bin = cursor + w - 1L,
        effect = fx, stringsAsFactors = FALSE
      )
      left <- left - w
      cursor <<- cursor + w + 8L
    }
  }
  add_blocks(n_gain, "GAIN", effect)
  add_blocks(n_loss, "LOSS", effect)
  if (!length(ev)) {
    return(data.frame(event_class = "NO_CHANGE", start_bin = 1L, end_bin = 1L,
                      effect = 0, stringsAsFactors = FALSE))
  }
  do.call(rbind, ev)
}

#' Background-bin effect metrics
#'
#' @param logFC,z Per-bin estimates.
#' @param background_idx Integer indices of true-background bins.
#' @export
background_effect_metrics <- function(logFC, z, background_idx) {
  b <- as.integer(background_idx)
  b <- b[b >= 1L & b <= length(logFC)]
  fb <- logFC[b]
  zb <- z[b]
  list(
    n_background = length(b),
    bias_background = stats::median(fb, na.rm = TRUE),
    mae_background = stats::median(abs(fb), na.rm = TRUE),
    bias_z = stats::median(zb, na.rm = TRUE),
    mae_z = stats::median(abs(zb), na.rm = TRUE)
  )
}

#' Fit one simulated cohort under a named normalization
#'
#' @param sim Output of [simulate_count_assay()].
#' @param method Normalization method.
#' @param run_discovery If TRUE, run segmentation + testing.
#' @param config Pipeline config (defaults used for missing keys).
#' @export
fit_normalized_effect <- function(sim, method, run_discovery = FALSE,
                                  config = NULL) {
  if (is.null(config)) config <- default_config()
  ref <- if (length(sim$spikein_idx)) sim$spikein_idx else NULL
  norm <- normalize_counts(
    sim$counts, method = method,
    extreme_factor_ratio = config$normalization$extreme_factor_ratio,
    spikein_usage = identical(method, "spikein"),
    background_definition = config$normalization$background_definition %||% 0.5,
    reference_bins = ref
  )
  valid <- sim$bins$valid
  if (length(sim$spikein_idx)) valid[sim$spikein_idx] <- FALSE
  eff <- estimate_differential_effect(
    sim$counts, sim$samples, offsets = norm$offsets,
    score = config$differential$score,
    control_level = config$differential$control_level,
    treatment_level = config$differential$treatment_level,
    valid = valid
  )
  bg <- background_effect_metrics(eff$logFC, eff$z, sim$background_idx)
  out <- list(norm = norm, effect = eff, background = bg, tested = NULL,
              candidates = NULL)
  if (!isTRUE(run_discovery)) return(out)
  cfg <- config
  cfg$normalization$method <- method
  fit <- run_simulated_discovery(sim, cfg)
  out$tested <- fit$tested
  out$candidates <- fit$candidates
  out$track <- fit$track
  out$runtime_sec <- fit$runtime_sec
  out
}

#' Summarize false discovery against planted truth intervals
#' @export
summarize_discovery_vs_truth <- function(tested, truth, q = 0.05,
                                         bin_size = 200L) {
  if (is.null(tested) || !nrow(tested) || !"FDR" %in% names(tested)) {
    pred <- data.frame(chrom = character(), start = integer(), end = integer())
  } else {
    pred <- tested[!is.na(tested$FDR) & tested$FDR <= q, , drop = FALSE]
    if (!all(c("chrom", "start", "end") %in% names(pred))) {
      pred <- data.frame(chrom = character(), start = integer(), end = integer())
    }
  }
  tr <- truth
  if (nrow(tr) && "event_class" %in% names(tr)) {
    tr <- tr[!tr$event_class %in% c("NO_CHANGE", "NULL"), ,
             drop = FALSE]
  }
  if (!nrow(tr)) {
    tr <- data.frame(chrom = character(), start = integer(), end = integer())
  }
  met <- compute_region_metrics(pred, tr)
  bmet <- compute_boundary_metrics(pred, tr, met$matches)
  n_false <- if (!nrow(tr)) nrow(pred) else {
    if (!nrow(pred)) 0L else sum(!seq_len(nrow(pred)) %in% met$matches$pred_i[met$matches$matched])
  }
  list(
    n_pred = nrow(pred), n_truth = nrow(tr),
    precision = met$precision, recall = met$recall, f1 = met$f1,
    mean_iou = met$mean_iou,
    n_false_regions = n_false,
    bp_precision = met$bp_precision, bp_recall = met$bp_recall,
    start_mae = bmet$start_mae, end_mae = bmet$end_mae,
    start_p90 = bmet$start_p90, end_p90 = bmet$end_p90,
    hausdorff = bmet$hausdorff,
    false_boundary_count = 2L * n_false
  )
}

.binomial_ci <- function(x, n, level = 0.95) {
  if (n <= 0) return(c(NA_real_, NA_real_))
  as.numeric(stats::binom.test(x, n, conf.level = level)$conf.int)
}
