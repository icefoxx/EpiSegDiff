#' Refine candidate region boundaries
#'
#' Methods:
#' - `local_gradient`: move each end to the local maximum of
#'   `|diff(score)|` inside a window.
#' - `local_likelihood`: maximize a two-mean Gaussian contrast in the
#'   window (step at the boundary).
#' - `cross_scale`: snap to the nearest nested finer-scale boundary if
#'   present, else `local_gradient`.
#'
#' @param regions data.frame with chrom, start, end, region_id.
#' @param track Effect track from [write_effect_track()].
#' @param method Refinement method.
#' @param search_window_bins Half-window in bins.
#' @param finer_segments Optional finer-scale segments for `cross_scale`.
#' @return data.frame with coarse/refined coordinates and support.
#' @export
refine_boundaries <- function(regions, track,
                              method = c("local_gradient", "local_likelihood",
                                         "cross_scale"),
                              search_window_bins = 5L,
                              finer_segments = NULL) {
  method <- match.arg(method)
  regions <- .check_bed_frame(regions, "regions")
  if (!"region_id" %in% names(regions)) regions$region_id <- seq_len(nrow(regions))
  if (nrow(regions) == 0L) {
    return(.empty_refined())
  }
  w <- as.integer(search_window_bins)
  rows <- lapply(seq_len(nrow(regions)), function(i) {
    .refine_one(regions[i, ], track, method, w, finer_segments)
  })
  do.call(rbind, rows)
}

.empty_refined <- function() {
  data.frame(
    region_id = integer(),
    coarse_start = integer(), coarse_end = integer(),
    refined_start = integer(), refined_end = integer(),
    start_support = numeric(), end_support = numeric(),
    method = character(),
    stringsAsFactors = FALSE
  )
}

.refine_one <- function(reg, track, method, w, finer) {
  chrom <- as.character(reg$chrom)
  idx <- which(track$chrom == chrom)
  if (!length(idx)) {
    return(data.frame(
      region_id = reg$region_id, coarse_start = reg$start, coarse_end = reg$end,
      refined_start = reg$start, refined_end = reg$end,
      start_support = NA_real_, end_support = NA_real_, method = method,
      stringsAsFactors = FALSE
    ))
  }
  tr <- track[idx, , drop = FALSE]
  s0 <- which.min(abs(tr$start - reg$start))
  e0 <- which.min(abs(tr$end - reg$end))
  if (e0 < s0) {
    tmp <- s0; s0 <- e0; e0 <- tmp
  }

  if (method == "cross_scale" && !is.null(finer) && nrow(finer)) {
    snapped <- .snap_cross_scale(reg, finer)
    if (!is.null(snapped)) {
      return(data.frame(
        region_id = reg$region_id,
        coarse_start = reg$start, coarse_end = reg$end,
        refined_start = snapped$start, refined_end = snapped$end,
        start_support = snapped$start_support, end_support = snapped$end_support,
        method = "cross_scale",
        stringsAsFactors = FALSE
      ))
    }
    method_use <- "local_gradient"
  } else {
    method_use <- method
  }

  if (method_use == "local_likelihood") {
    rs <- .best_step(tr, s0, w, side = "start")
    re <- .best_step(tr, e0, w, side = "end")
  } else {
    rs <- .best_gradient(tr, s0, w)
    re <- .best_gradient(tr, e0, w)
  }
  rs_i <- min(rs$index, re$index)
  re_i <- max(rs$index, re$index)
  data.frame(
    region_id = reg$region_id,
    coarse_start = reg$start, coarse_end = reg$end,
    refined_start = tr$start[rs_i], refined_end = tr$end[re_i],
    start_support = rs$support, end_support = re$support,
    method = method_use,
    stringsAsFactors = FALSE
  )
}

.best_gradient <- function(tr, i0, w) {
  n <- nrow(tr)
  lo <- max(2L, i0 - w)
  hi <- min(n, i0 + w)
  if (hi < lo) {
    return(list(index = i0, support = NA_real_))
  }
  d <- abs(diff(tr$score))
  # gradient at boundary k is |score[k]-score[k-1]|, associated with start of bin k
  sl <- max(1L, lo - 1L)
  sh <- min(length(d), hi)
  if (sl > sh) return(list(index = i0, support = NA_real_))
  win <- d[sl:sh]
  if (!any(is.finite(win))) return(list(index = i0, support = NA_real_))
  j <- which.max(replace(win, !is.finite(win), -Inf))
  idx <- sl + j # bin index of the right side of the jump
  idx <- min(max(idx, 1L), n)
  list(index = idx, support = as.numeric(win[j]))
}

.best_step <- function(tr, i0, w, side) {
  n <- nrow(tr)
  lo <- max(2L, i0 - w)
  hi <- min(n - 1L, i0 + w)
  if (hi < lo) return(list(index = i0, support = NA_real_))
  y <- tr$score
  best <- -Inf
  best_i <- i0
  for (k in lo:hi) {
    left <- y[max(1L, k - w): (k - 1L)]
    right <- y[k:min(n, k + w)]
    left <- left[is.finite(left)]
    right <- right[is.finite(right)]
    if (length(left) < 1L || length(right) < 1L) next
    # two-mean contrast
    stat <- length(left) * length(right) / (length(left) + length(right)) *
      (mean(left) - mean(right))^2
    if (is.finite(stat) && stat > best) {
      best <- stat
      best_i <- k
    }
  }
  list(index = best_i, support = if (is.finite(best)) best else NA_real_)
}

.snap_cross_scale <- function(reg, finer) {
  f <- finer[as.character(finer$chrom) == as.character(reg$chrom), , drop = FALSE]
  if (!nrow(f)) return(NULL)
  # nearest starts/ends among nested or overlapping finer segments
  ov <- .interval_intersection(reg$start, reg$end, f$start, f$end) > 0L
  if (!any(ov)) return(NULL)
  f <- f[ov, , drop = FALSE]
  s <- f$start[which.min(abs(f$start - reg$start))]
  e <- f$end[which.min(abs(f$end - reg$end))]
  if (e <= s) return(NULL)
  list(start = s, end = e, start_support = 1, end_support = 1)
}

#' Infer condition-specific boundaries near a differential candidate
#'
#' Uses only control/treatment signal in a local window. Does **not**
#' call MACS, SEACR, or any external peak caller. Peak callers may be
#' compared later as diagnostics only.
#'
#' @param counts Bin × sample matrix.
#' @param bins Bin table.
#' @param samples Sample metadata.
#' @param region One interval (`chrom`, `start`, `end`).
#' @param window_bins Half-window in bins around the candidate.
#' @param method `local_threshold` or `local_changepoint`.
#' @param size_factors Optional size factors.
#' @param k MAD multiplier for `local_threshold`.
#' @return data.frame with `condition`, `L`, `R`, `method`, `support`.
#' @export
infer_local_condition_boundaries <- function(counts, bins, samples, region,
                                             window_bins = 15L,
                                             method = c("local_threshold",
                                                        "local_changepoint",
                                                        "end_walk"),
                                             size_factors = NULL,
                                             k = 2,
                                             search_radius_bins = 8L,
                                             walk_version = c("symmetric", "legacy")) {
  method <- match.arg(method)
  walk_version <- match.arg(walk_version)
  region <- .check_bed_frame(as.data.frame(region), "region")
  if (nrow(region) != 1L) stop("region must be a single interval.", call. = FALSE)
  if (is.null(size_factors)) {
    sf <- colSums(counts)
    sf <- sf / exp(mean(log(pmax(sf, 1))))
  } else {
    sf <- as.numeric(size_factors)
  }
  chrom <- as.character(region$chrom[1])
  idx <- which(bins$chrom == chrom)
  if (!length(idx)) {
    return(data.frame(condition = character(), L = integer(), R = integer(),
                      method = character(), support = numeric(),
                      stringsAsFactors = FALSE))
  }
  loc <- bins[idx, , drop = FALSE]
  center_s <- which.min(abs(loc$start - region$start[1]))
  center_e <- which.min(abs(loc$end - region$end[1]))
  w <- as.integer(window_bins)
  lo <- max(1L, min(center_s, center_e) - w)
  hi <- min(nrow(loc), max(center_s, center_e) + w)
  win <- loc[lo:hi, , drop = FALSE]
  gidx <- idx[lo:hi]
  conds <- unique(as.character(samples$condition))
  rows <- lapply(conds, function(cn) {
    keep <- samples$condition == cn
    mat <- sweep(counts[gidx, keep, drop = FALSE], 2, sf[keep], "/")
    mu <- rowMeans(mat)
    if (method == "end_walk") {
      got <- end_walk_series(
        mu, win$start, win$end, region$start[1], region$end[1],
        search_radius = as.integer(search_radius_bins),
        version = walk_version
      )
      L <- got$L
      R <- got$R
      supp <- NA_real_
    } else if (method == "local_changepoint") {
      if (sum(is.finite(mu)) < 4L) {
        return(data.frame(condition = cn, L = region$start[1], R = region$end[1],
                          method = method, support = NA_real_,
                          stringsAsFactors = FALSE))
      }
      sc <- mu
      sc[!is.finite(sc)] <- median(mu, na.rm = TRUE)
      seg <- segment_effect(sc, positions = win$start, penalty = 4,
                            min_width = 2L, loss = "gaussian", engine = "pelt")
      if (!nrow(seg)) {
        L <- win$start[1]; R <- win$end[nrow(win)]; supp <- NA_real_
      } else {
        j <- which.max(seg$effect_mean)
        L <- seg$start[j]
        last <- lo + seg$end_idx[j] - 1L
        last <- min(max(last, 1L), nrow(win))
        R <- win$end[last]
        supp <- seg$effect_mean[j]
      }
    } else {
      tail_n <- max(2L, floor(nrow(win) / 5))
      bg <- c(mu[seq_len(tail_n)], mu[seq.int(nrow(win) - tail_n + 1L, nrow(win))])
      bg <- bg[is.finite(bg)]
      b0 <- if (length(bg)) stats::median(bg) else stats::median(mu, na.rm = TRUE)
      s0 <- if (length(bg) >= 3L) stats::mad(bg, constant = 1.4826) else 0
      thr <- b0 + k * max(s0, 1e-6)
      on <- is.finite(mu) & mu >= thr
      if (!any(on)) {
        L <- region$start[1]; R <- region$end[1]; supp <- NA_real_
      } else {
        first <- which(on)[1]
        last <- tail(which(on), 1)
        L <- win$start[first]
        R <- win$end[last]
        supp <- mean(mu[on])
      }
    }
    data.frame(condition = cn, L = as.integer(L), R = as.integer(R),
               method = method, support = as.numeric(supp),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
