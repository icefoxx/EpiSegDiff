#' v0.8.1 same-direction nested instrumentation for B_hierarchy
#' @name v081_nested

#' Simulate nested broad + local with controllable geometry
#' @export
simulate_nested_pair <- function(broad_width_bp, local_width_bp, B, L,
                                 res = 200L, noise_sd = 0.2, seed = 1L,
                                 position = c("center", "p25", "p10",
                                              "near_left", "near_right"),
                                 boundary_gap_bp = NA_integer_) {
  position <- match.arg(position)
  nb <- max(1L, as.integer(broad_width_bp / res))
  nl <- max(1L, as.integer(local_width_bp / res))
  if (nl > nb) nl <- nb
  pad <- max(40L, as.integer(nb / 2) + 20L)
  n_bins <- nb + 2L * pad
  b0 <- pad + 1L
  b1 <- b0 + nb - 1L
  # place local
  if (!is.na(boundary_gap_bp)) {
    gap <- max(0L, as.integer(boundary_gap_bp / res))
    l0 <- b0 + gap
    l1 <- l0 + nl - 1L
    if (l1 > b1) {
      l1 <- b1
      l0 <- max(b0, l1 - nl + 1L)
    }
  } else if (position == "center") {
    l0 <- b0 + floor((nb - nl) / 2)
    l1 <- l0 + nl - 1L
  } else if (position == "p25") {
    l0 <- b0 + floor(0.25 * (nb - nl))
    l1 <- l0 + nl - 1L
  } else if (position == "p10") {
    l0 <- b0 + floor(0.10 * (nb - nl))
    l1 <- l0 + nl - 1L
  } else if (position == "near_left") {
    l0 <- b0
    l1 <- l0 + nl - 1L
  } else {
    l1 <- b1
    l0 <- l1 - nl + 1L
  }
  ev <- data.frame(start_bin = c(b0, l0), end_bin = c(b1, l1),
                   delta = c(B, L), shape = "sharp")
  sim <- simulate_continuous_events(n_bins, res, ev, noise_sd = noise_sd,
                                    seed = seed)
  dist_left <- (l0 - b0) * res
  dist_right <- (b1 - l1) * res
  sim$meta <- list(
    B = B, L = L, broad_width_bp = broad_width_bp, local_width_bp = local_width_bp,
    res = res, noise_sd = noise_sd, position = position,
    b0 = b0, b1 = b1, l0 = l0, l1 = l1,
    dist_boundary_bp = min(dist_left, dist_right),
    width_ratio = local_width_bp / broad_width_bp
  )
  sim
}

#' Parent residual constructors
.parent_residual <- function(score, ii, method = c("mean", "median", "trim10")) {
  method <- match.arg(method)
  x <- score[ii]
  mu <- switch(method,
               mean = mean(x),
               median = stats::median(x),
               trim10 = mean(x, trim = 0.10))
  list(resid = score[ii] - mu, parent_effect = mu, method = method)
}

#' Instrumented B_hierarchy with stagewise child loss traces
#'
#' @param parent_method `mean`, `median`, or `trim10` residualization.
#' @param child_engine `fused_lasso` or `pelt`.
#' @param child_thr_mult multiplier on MAD threshold for children.
#' @param child_penalty fused/pelt penalty for children.
#' @param reconcile if TRUE, drop children whose |effect| < 0.5*|parent|.
#' @export
pathb_hierarchy_trace <- function(score, bed, truth = NULL, bin_size = 200L,
                                  w_B = 10L, penalty = 8,
                                  parent_method = c("mean", "median", "trim10"),
                                  child_engine = c("fused_lasso", "pelt"),
                                  child_thr_mult = 0.5,
                                  child_penalty = NULL,
                                  reconcile = FALSE,
                                  allow_orphan_local = FALSE) {
  parent_method <- match.arg(parent_method)
  child_engine <- match.arg(child_engine)
  if (is.null(child_penalty)) child_penalty <- max(2, penalty / 2)
  score <- as.numeric(score)
  thr_p <- .pathb_mad_thr(score)
  segP <- segment_signal(score, positions = bed$start, engine = "pelt",
                         penalty = penalty, min_width = w_B)
  parents <- .pathb_seg_to_regions(segP, bed, score, "broad", 1L, thr = thr_p)
  rows <- list()
  rid <- 0L
  traces <- list()
  # truth local interval in genomic coords if provided
  truth_local <- NULL
  truth_L <- NA_real_
  if (!is.null(truth) && nrow(truth) >= 2L) {
    truth_L <- truth$delta[2]
    truth_local <- list(
      start = as.integer((truth$start_bin[2] - 1L) * bin_size),
      end = as.integer(truth$end_bin[2] * bin_size)
    )
  }
  if (nrow(parents)) {
    for (i in seq_len(nrow(parents))) {
      rid <- rid + 1L
      pr <- parents[i, , drop = FALSE]
      pr$region_id <- rid
      pr$parent_id <- NA_integer_
      rows[[length(rows) + 1L]] <- pr
      ii <- which(bed$start >= pr$start & bed$end <= pr$end)
      if (length(ii) < 4L) next
      bed_i <- bed[ii, , drop = FALSE]
      prs <- .parent_residual(score, ii, parent_method)
      resid <- prs$resid
      # local contrast before/after on truth local bins if overlapping
      local_before <- NA_real_
      local_after <- NA_real_
      absorption <- NA_real_
      residual_snr <- NA_real_
      if (!is.null(truth_local)) {
        jj <- which(bed_i$start < truth_local$end & bed_i$end > truth_local$start)
        if (length(jj)) {
          local_before <- mean(score[ii][jj])
          local_after <- mean(resid[jj])
          if (is.finite(truth_L) && abs(truth_L) > 1e-8)
            absorption <- 1 - local_after / truth_L
          residual_snr <- abs(local_after) /
            max(stats::mad(resid, constant = 1), 1e-6)
        }
      }
      # raw child segments (no threshold)
      segC <- segment_signal(resid, positions = bed_i$start, engine = child_engine,
                             penalty = child_penalty,
                             min_width = if (child_engine == "pelt") 1L else 1L)
      thr_c <- max(0.05, child_thr_mult * .pathb_mad_thr(resid))
      raw_n <- nrow(segC)
      # thresholded
      ch <- .pathb_seg_to_regions(segC, bed_i, resid, "local", 2L,
                                  parent_id = rid, thr = thr_c)
      thr_n <- nrow(ch)
      # optional reconciliation: require child |effect| >= 0.5 |parent|
      removed_recon <- 0L
      if (reconcile && nrow(ch)) {
        keep <- abs(ch$effect) >= 0.5 * abs(pr$effect)
        removed_recon <- sum(!keep)
        ch <- ch[keep, , drop = FALSE]
      }
      final_n <- nrow(ch)
      # did any child overlap truth local?
      child_hit_raw <- FALSE
      child_hit_thr <- FALSE
      child_hit_final <- FALSE
      .hit_local <- function(reg) {
        if (is.null(truth_local) || !nrow(reg)) return(FALSE)
        ok <- !is.na(reg$chrom) & is.finite(reg$start) & is.finite(reg$end) &
          reg$end > reg$start
        reg <- reg[ok, , drop = FALSE]
        if (!nrow(reg)) return(FALSE)
        g <- GenomicRanges::GRanges(
          as.character(reg$chrom),
          IRanges::IRanges(as.integer(reg$start) + 1L, as.integer(reg$end))
        )
        tgr <- GenomicRanges::GRanges(
          as.character(bed$chrom[1]),
          IRanges::IRanges(truth_local$start + 1L, truth_local$end)
        )
        any(GenomicRanges::countOverlaps(g, tgr) > 0)
      }
      if (!is.null(truth_local) && raw_n) {
        raw_reg <- .pathb_seg_to_regions(segC, bed_i, resid, "local", 2L,
                                         parent_id = rid, thr = 0)
        child_hit_raw <- .hit_local(raw_reg)
      }
      ch_thr <- .pathb_seg_to_regions(segC, bed_i, resid, "local", 2L,
                                      parent_id = rid, thr = thr_c)
      child_hit_thr <- .hit_local(ch_thr)
      child_hit_final <- .hit_local(ch)
      traces[[length(traces) + 1L]] <- data.frame(
        parent_id = rid,
        parent_effect = pr$effect,
        parent_method = parent_method,
        child_engine = child_engine,
        child_thr = thr_c,
        local_signal_before = local_before,
        residual_local_effect = local_after,
        absorption = absorption,
        residual_snr = residual_snr,
        n_raw_child_segs = raw_n,
        n_thr_child = thr_n,
        n_final_child = final_n,
        n_removed_recon = removed_recon,
        child_hit_raw = child_hit_raw,
        child_hit_thr = child_hit_thr,
        child_hit_final = child_hit_final,
        stringsAsFactors = FALSE
      )
      if (nrow(ch)) {
        for (j in seq_len(nrow(ch))) {
          rid <- rid + 1L
          ch$region_id[j] <- rid
          ch$parent_id[j] <- pr$region_id
          rows[[length(rows) + 1L]] <- ch[j, , drop = FALSE]
        }
      }
    }
  }
  if (allow_orphan_local && !nrow(parents)) {
    segL <- segment_signal(score, positions = bed$start, engine = "fused_lasso",
                           penalty = penalty, min_width = 1L)
    loc <- .pathb_seg_to_regions(segL, bed, score, "local", 1L, thr = thr_p)
    if (nrow(loc)) {
      for (j in seq_len(nrow(loc))) {
        rid <- rid + 1L
        loc$region_id[j] <- rid
        loc$parent_id[j] <- NA_integer_
        rows[[length(rows) + 1L]] <- loc[j, , drop = FALSE]
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else .empty_pathb_regions()
  tr <- if (length(traces)) do.call(rbind, traces) else data.frame()
  list(regions = out, trace = tr, architecture = "B_hierarchy_trace",
       params = list(w_B = w_B, parent_method = parent_method,
                     child_engine = child_engine, child_thr_mult = child_thr_mult,
                     child_penalty = child_penalty, reconcile = reconcile))
}

#' Evaluate nested hierarchy + absorption summary for one sim
#' @export
eval_nested_trace <- function(sim, config = list()) {
  defaults <- list(w_B = 10L, penalty = 8, parent_method = "mean",
                   child_engine = "fused_lasso", child_thr_mult = 0.5,
                   child_penalty = NULL, reconcile = FALSE)
  cfg <- modifyList(defaults, config)
  sc <- compute_diff_score(sim$signal, sim$samples, "JNKi", "DMSO", "raw_effect")
  fit <- pathb_hierarchy_trace(
    sc, sim$bed, truth = sim$truth, bin_size = sim$bin_size,
    w_B = cfg$w_B, penalty = cfg$penalty,
    parent_method = cfg$parent_method, child_engine = cfg$child_engine,
    child_thr_mult = cfg$child_thr_mult, child_penalty = cfg$child_penalty,
    reconcile = cfg$reconcile, allow_orphan_local = FALSE
  )
  nh <- score_nested_hierarchy(fit$regions, sim$truth, sim$bin_size)
  sc0 <- score_calls_vs_truth(fit$regions, sim$truth, sim$bin_size)
  # parent boundary error vs broad truth
  broad <- fit$regions[fit$regions$component_type == "broad", , drop = FALSE]
  parent_bound_err <- NA_real_
  if (nrow(broad) && nrow(sim$truth)) {
    t0 <- (sim$truth$start_bin[1] - 1L) * sim$bin_size
    t1 <- sim$truth$end_bin[1] * sim$bin_size
    parent_bound_err <- min(abs(broad$start - t0) + abs(broad$end - t1))
  }
  absorp <- if (nrow(fit$trace)) mean(fit$trace$absorption, na.rm = TRUE) else NA
  rsnr <- if (nrow(fit$trace)) mean(fit$trace$residual_snr, na.rm = TRUE) else NA
  p_cand <- if (nrow(fit$trace)) mean(fit$trace$child_hit_raw) else 0
  p_pass <- if (nrow(fit$trace)) {
    mean(fit$trace$child_hit_thr[fit$trace$child_hit_raw], na.rm = TRUE)
  } else NA
  if (length(p_pass) != 1 || is.nan(p_pass)) p_pass <- NA_real_
  p_ret <- if (nrow(fit$trace)) {
    mean(fit$trace$child_hit_final[fit$trace$child_hit_thr], na.rm = TRUE)
  } else NA
  if (length(p_ret) != 1 || is.nan(p_ret)) p_ret <- NA_real_
  meta <- sim$meta
  data.frame(
    B = meta$B, L = meta$L, L_over_B = meta$L / meta$B,
    broad_width_bp = meta$broad_width_bp, local_width_bp = meta$local_width_bp,
    width_ratio = meta$width_ratio, res = meta$res, noise_sd = meta$noise_sd,
    position = meta$position, dist_boundary_bp = meta$dist_boundary_bp,
    parent_method = cfg$parent_method, child_engine = cfg$child_engine,
    child_thr_mult = cfg$child_thr_mult, reconcile = cfg$reconcile,
    p_parent = nh$p_parent, p_child = nh$p_child, p_both = nh$p_both,
    p_merged = nh$p_merged, false_child = nh$false_child,
    f1 = sc0$f1, n_regions = nrow(fit$regions),
    n_broad = sum(fit$regions$component_type == "broad"),
    n_local = sum(fit$regions$component_type == "local"),
    absorption = absorp, residual_snr = rsnr,
    parent_bound_err = parent_bound_err,
    p_candidate_exists = p_cand,
    p_pass_threshold = p_pass,
    p_retained = p_ret,
    stringsAsFactors = FALSE
  )
}
