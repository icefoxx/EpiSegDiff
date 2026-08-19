#' Closure-analysis helpers (v0.7 review / pre-v0.8).
#'
#' These functions score simulations. They do not implement a
#' hierarchical production architecture.

#' MAD-threshold calls for one score × engine (same rule as v0.7 factorial)
#' @export
closure_calls <- function(sim, score_engine, seg_engine, min_width,
                          penalty = 8) {
  sc <- compute_diff_score(sim$signal, sim$samples, "JNKi", "DMSO",
                           score_engine)
  seg <- segment_signal(sc, positions = sim$bed$start, engine = seg_engine,
                        penalty = penalty, min_width = min_width)
  thr <- 0.15 + 0.35 * stats::mad(sc, constant = 1)
  if (!is.finite(thr) || thr <= 0) thr <- 0.2
  .segments_to_calls(seg, sim$bed, thresh = thr)
}

#' Transition-zone truth for a gradual planted block
#'
#' Uses the same 20% end-ramp rule as `simulate_continuous_events`.
#' @export
transition_truth <- function(start_bin, end_bin, bin_size) {
  w <- as.integer(end_bin - start_bin + 1L)
  t1 <- max(1L, floor(0.2 * w))
  if (w < 5L) {
    return(data.frame(
      zone = "core",
      start = as.integer((start_bin - 1L) * bin_size),
      end = as.integer(end_bin * bin_size),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    zone = c("left_transition", "core", "right_transition"),
    start = as.integer(c(start_bin - 1L, start_bin + t1 - 1L,
                         end_bin - t1) * bin_size),
    end = as.integer(c(start_bin + t1 - 1L, end_bin - t1,
                       end_bin) * bin_size),
    stringsAsFactors = FALSE
  )
}

#' Score calls against core / transition / background
#' @export
score_transition_zones <- function(calls, start_bin, end_bin, bin_size,
                                   series_end) {
  z <- transition_truth(start_bin, end_bin, bin_size)
  zones <- GenomicRanges::GRanges(
    "1", IRanges::IRanges(z$start + 1L, z$end)
  )
  names(zones) <- z$zone
  full <- GenomicRanges::GRanges(
    "1", IRanges::IRanges(as.integer((start_bin - 1L) * bin_size + 1L),
                          as.integer(end_bin * bin_size))
  )
  bg <- GenomicRanges::GRanges("1", IRanges::IRanges(1L, as.integer(series_end)))
  bg <- GenomicRanges::setdiff(bg, full)
  if (!nrow(calls)) {
    return(data.frame(
      core_recall = 0, trans_recall = 0, trans_overlap = 0,
      over_extension = 0, under_extension = 1,
      midpoint_err = NA_real_, n_call = 0L,
      stringsAsFactors = FALSE
    ))
  }
  cgr <- GenomicRanges::GRanges(
    calls$chrom, IRanges::IRanges(as.integer(calls$start) + 1L,
                                  as.integer(calls$end))
  )
  ov_frac <- function(a, b) {
    if (!length(a) || !length(b) || sum(GenomicRanges::width(a)) == 0)
      return(0)
    inter <- GenomicRanges::intersect(a, b)
    sum(as.numeric(GenomicRanges::width(inter))) /
      sum(as.numeric(GenomicRanges::width(a)))
  }
  core <- zones["core"]
  trans <- zones[names(zones) != "core"]
  core_rec <- ov_frac(core, cgr)
  trans_rec <- ov_frac(trans, cgr)
  pred <- sum(as.numeric(GenomicRanges::width(cgr)))
  bg_bp <- if (length(bg)) {
    sum(as.numeric(GenomicRanges::width(GenomicRanges::intersect(cgr, bg))))
  } else 0
  over <- if (pred > 0) bg_bp / pred else 0
  under <- 1 - core_rec
  # midpoint of planted interval vs nearest call
  true_mid <- ((start_bin - 1L) * bin_size + end_bin * bin_size) / 2
  call_mid <- (calls$start + calls$end) / 2
  mid_err <- min(abs(call_mid - true_mid))
  data.frame(
    core_recall = core_rec,
    trans_recall = trans_rec,
    trans_overlap = trans_rec,
    over_extension = over,
    under_extension = under,
    midpoint_err = mid_err,
    n_call = nrow(calls),
    stringsAsFactors = FALSE
  )
}

#' Detection / conditional-boundary / unconditional F1
#' @export
decompose_detection_boundary <- function(scored) {
  det <- is.finite(scored$n_detected) & scored$n_detected > 0
  data.frame(
    n = nrow(scored),
    p_detect = mean(det),
    f1_uncond = mean(scored$f1, na.rm = TRUE),
    recall_uncond = mean(scored$recall, na.rm = TRUE),
    precision_uncond = mean(scored$precision, na.rm = TRUE),
    med_bound_err_cond = if (any(det & is.finite(scored$med_bound_err)))
      stats::median(scored$med_bound_err[det], na.rm = TRUE) else NA_real_,
    mean_bound_err_cond = if (any(det & is.finite(scored$med_bound_err)))
      mean(scored$med_bound_err[det], na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}
