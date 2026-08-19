#' Walk left and right ends on a 1-D series
#'
#' `version = "legacy"` is the v0.2.1 implementation (kept for audit).
#' `version = "symmetric"` is the v0.2.2 fix: each end searches only
#' its local edge, and jump index `k` (between bins `k` and `k+1`)
#' maps to `L = start[k+1]` and `R = end[k]`.
#'
#' @param y Numeric series (one value per bin).
#' @param start,end Integer bin coordinates (0-based half-open).
#' @param cand_L,cand_R Candidate interval in the same coordinates.
#' @param search_radius Max bins to search from each candidate end.
#' @param version `legacy` or `symmetric`.
#' @return list with `L`, `R`, `left_edge`, `right_edge`.
#' @export
end_walk_series <- function(y, start, end, cand_L, cand_R,
                            search_radius = 5L, version = c("symmetric", "legacy")) {
  version <- match.arg(version)
  y <- as.numeric(y)
  start <- as.integer(start)
  end <- as.integer(end)
  n <- length(y)
  if (length(start) != n || length(end) != n) {
    stop("y, start, and end must have the same length.", call. = FALSE)
  }
  if (n == 0L) return(list(L = cand_L, R = cand_R, left_edge = NA, right_edge = NA))
  iL <- which.min(abs(start - as.integer(cand_L)))
  iR <- which.min(abs(end - as.integer(cand_R)))
  if (iR < iL) {
    tmp <- iL; iL <- iR; iR <- tmp
  }
  r <- as.integer(search_radius)
  if (n == 1L) return(list(L = start[1], R = end[1], left_edge = NA, right_edge = NA))
  d <- abs(diff(y))
  # d[k] lives on the edge between bin k and bin k+1 (k = 1..n-1)

  if (version == "legacy") {
    pick <- function(i0) {
      lo2 <- max(1L, i0 - r)
      hi2 <- min(length(d), i0 + r)
      if (hi2 < lo2 || !any(is.finite(d[lo2:hi2]))) return(i0)
      lo2 + which.max(replace(d[lo2:hi2], !is.finite(d[lo2:hi2]), -Inf)) - 1L
    }
    is <- pick(iL)
    ie <- pick(iR)
    if (ie < is) {
      tmp <- is; is <- ie; ie <- tmp
    }
    return(list(
      L = start[max(1L, is)],
      R = end[min(n, ie)],
      left_edge = is, right_edge = ie, version = "legacy"
    ))
  }

  # Symmetric: restrict each search to its own half, break ties toward
  # the candidate end, map edges to [start, end) consistently.
  mid_edge <- max(1L, min(length(d), as.integer(floor((iL + iR) / 2))))
  pick_edge <- function(i0, lo_lim, hi_lim, prefer = c("left", "right")) {
    prefer <- match.arg(prefer)
    lo2 <- max(lo_lim, i0 - r)
    hi2 <- min(hi_lim, i0 + r)
    if (hi2 < lo2) return(NA_integer_)
    sl <- d[lo2:hi2]
    sl[!is.finite(sl)] <- -Inf
    if (!any(is.finite(sl))) return(NA_integer_)
    mx <- max(sl)
    hits <- which(sl == mx)
    j <- if (prefer == "left") hits[1] else hits[length(hits)]
    lo2 + j - 1L
  }
  eL <- pick_edge(iL, 1L, mid_edge, prefer = "right")
  eR <- pick_edge(iR, mid_edge, length(d), prefer = "left")
  if (is.na(eL)) eL <- max(1L, iL - 1L)
  if (is.na(eR)) eR <- min(length(d), iR)
  L <- start[min(n, eL + 1L)]
  R <- end[min(n, eR)]
  if (R <= L) {
    L <- start[iL]
    R <- end[iR]
  }
  list(L = L, R = R, left_edge = eL, right_edge = eR, version = "symmetric")
}

#' Local-gradient refinement on a series (same mapping as symmetric walk)
#' @export
gradient_series <- function(y, start, end, cand_L, cand_R, search_radius = 5L) {
  end_walk_series(y, start, end, cand_L, cand_R, search_radius, version = "symmetric")
}

#' Local two-mean likelihood refinement on a series
#' @export
likelihood_series <- function(y, start, end, cand_L, cand_R, search_radius = 5L) {
  n <- length(y)
  iL <- which.min(abs(start - as.integer(cand_L)))
  iR <- which.min(abs(end - as.integer(cand_R)))
  r <- as.integer(search_radius)
  best_step <- function(i0) {
    lo <- max(2L, i0 - r)
    hi <- min(n - 1L, i0 + r)
    if (hi < lo) return(i0)
    best <- -Inf
    best_i <- i0
    for (k in lo:hi) {
      left <- y[max(1L, k - r):(k - 1L)]
      right <- y[k:min(n, k + r)]
      left <- left[is.finite(left)]
      right <- right[is.finite(right)]
      if (!length(left) || !length(right)) next
      stat <- (length(left) * length(right) / (length(left) + length(right))) *
        (mean(left) - mean(right))^2
      if (is.finite(stat) && stat > best) {
        best <- stat
        best_i <- k
      }
    }
    best_i
  }
  a <- best_step(iL)
  b <- best_step(iR)
  if (b < a) {
    tmp <- a; a <- b; b <- tmp
  }
  list(L = start[a], R = end[b], left_edge = a, right_edge = b, version = "likelihood")
}

#' Noiseless step-function series with a single block
#' @export
make_step_series <- function(n = 13L, L_idx = 5L, R_idx = 9L, height = 2,
                             bin_size = 100L) {
  y <- numeric(n)
  y[L_idx:R_idx] <- height
  start <- as.integer((seq_len(n) - 1L) * bin_size)
  end <- start + as.integer(bin_size)
  list(
    y = y, start = start, end = end,
    true_L = start[L_idx], true_R = end[R_idx],
    L_idx = L_idx, R_idx = R_idx, n = n, bin_size = bin_size
  )
}

#' Map an interval through series reversal
#'
#' Forward `[L,R)` on a length-`N` coordinate axis (0-based, last end = N)
#' becomes `[N-R, N-L)` after reversing bins.
#' @export
reverse_interval <- function(L, R, N) {
  list(L = as.integer(N) - as.integer(R), R = as.integer(N) - as.integer(L))
}

#' Reverse a binned series
#' @export
reverse_series <- function(y, start, end) {
  n <- length(y)
  width <- end - start
  N <- end[n]
  yrev <- rev(y)
  # reversed bins still 0-based half-open covering [0, N)
  start_rev <- N - rev(end)
  end_rev <- N - rev(start)
  list(y = yrev, start = as.integer(start_rev), end = as.integer(end_rev), N = N)
}

#' Reconstruct a region from overlapping bins
#' @export
reconstruct_interval <- function(bins, region) {
  hit <- bins$chrom == region$chrom & bins$start < region$end & bins$end > region$start
  if (!any(hit)) {
    return(list(start = NA_integer_, end = NA_integer_, n_bins = 0L))
  }
  list(start = min(bins$start[hit]), end = max(bins$end[hit]), n_bins = sum(hit))
}
