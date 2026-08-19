#' Segment a differential-effect score
#'
#' Operates on a single contiguous numeric series. Callers must split
#' chromosomes / CVRs first. NA values are not allowed inside `score`.
#'
#' @param score Numeric vector (no NA).
#' @param positions Optional integer positions (bin index or genomic start).
#' @param penalty Non-negative penalty \(\lambda\).
#' @param min_width Minimum segment length in bins.
#' @param loss `gaussian`, `huber`, or `winsorized_gaussian`.
#' @param engine `pelt`, `fpop`, or `op`.
#' @param huber_k Huber threshold.
#' @param winsor_k MAD multiples for winsorization.
#' @return data.frame of segments with `start_idx`, `end_idx`, `start`,
#'   `end`, `effect_mean`, `n_bins`, `penalty`, `loss`.
#' @export
segment_effect <- function(score, positions = NULL, penalty,
                           min_width = 1L,
                           loss = c("gaussian", "huber", "winsorized_gaussian"),
                           engine = c("pelt", "fpop", "op"),
                           huber_k = 1.345,
                           winsor_k = 3) {
  loss <- match.arg(loss)
  engine <- match.arg(engine)
  score <- as.numeric(score)
  n <- length(score)
  if (!is.null(positions) && length(positions) != n) {
    stop("positions and score have different lengths.", call. = FALSE)
  }
  if (is.null(positions)) positions <- seq_len(n)
  if (anyNA(score)) {
    stop("score contains NA; split into contiguous valid regions first.", call. = FALSE)
  }
  if (!is.finite(penalty) || penalty < 0) {
    stop("penalty must be finite and non-negative.", call. = FALSE)
  }
  min_width <- as.integer(min_width)
  if (is.na(min_width) || min_width < 1L) {
    stop("min_width must be a positive integer.", call. = FALSE)
  }

  if (n == 0L) {
    return(.empty_segments())
  }
  # A series shorter than min_width is not an emitted segment.
  if (n < min_width) {
    return(.empty_segments())
  }

  y <- score
  if (loss == "winsorized_gaussian") {
    y <- .winsorize_mad(y, winsor_k)
    loss_fit <- "gaussian"
  } else {
    loss_fit <- loss
  }

  if (loss_fit == "gaussian" && engine == "fpop") {
    cps <- .fpop_gaussian(y, penalty, min_width)
  } else if (loss_fit == "gaussian" && engine == "pelt") {
    cps <- .pelt_gaussian(y, penalty, min_width)
  } else if (loss_fit == "gaussian" && engine == "op") {
    cps <- .op_gaussian(y, penalty, min_width)
  } else if (loss_fit == "huber") {
    if (engine == "fpop") {
      warning("Huber FPOP is not implemented; using PELT with Huber cost.",
              call. = FALSE)
    }
    if (engine == "op") {
      cps <- .op_huber(y, penalty, min_width, huber_k)
    } else {
      cps <- .pelt_huber(y, penalty, min_width, huber_k)
    }
  } else {
    stop("Unsupported loss/engine combination.", call. = FALSE)
  }

  seg <- .segments_from_cps(cps, score, positions, penalty, loss)
  .enforce_min_bins(seg, score, positions, min_width, penalty, loss)
}

# Merge or drop segments with n_bins < min_width so emitted widths
# match the documented minimum. Short last-segments were previously
# allowed by PELT (`len < min_width && t < n`).
.enforce_min_bins <- function(seg, score, positions, min_width, penalty, loss) {
  if (!nrow(seg) || min_width <= 1L) return(seg)
  changed <- TRUE
  while (changed && nrow(seg)) {
    changed <- FALSE
    short <- which(seg$n_bins < min_width)
    if (!length(short)) break
    i <- short[length(short)]
    if (i > 1L) {
      seg$end_idx[i - 1L] <- seg$end_idx[i]
      seg$end[i - 1L] <- seg$end[i]
      idx <- seg$start_idx[i - 1L]:seg$end_idx[i - 1L]
      seg$n_bins[i - 1L] <- length(idx)
      seg$effect_mean[i - 1L] <- mean(score[idx])
      seg <- seg[-i, , drop = FALSE]
      changed <- TRUE
    } else if (nrow(seg) > 1L) {
      seg$start_idx[2L] <- seg$start_idx[1L]
      seg$start[2L] <- seg$start[1L]
      idx <- seg$start_idx[2L]:seg$end_idx[2L]
      seg$n_bins[2L] <- length(idx)
      seg$effect_mean[2L] <- mean(score[idx])
      seg <- seg[-1L, , drop = FALSE]
      changed <- TRUE
    } else {
      return(.empty_segments())
    }
  }
  if (nrow(seg)) rownames(seg) <- NULL
  seg
}

.empty_segments <- function() {
  data.frame(
    start_idx = integer(), end_idx = integer(),
    start = integer(), end = integer(),
    effect_mean = numeric(), n_bins = integer(),
    penalty = numeric(), loss = character(),
    stringsAsFactors = FALSE
  )
}

.segments_from_cps <- function(cps, score, positions, penalty, loss) {
  n <- length(score)
  ends <- c(as.integer(cps), n)
  starts <- c(1L, as.integer(cps) + 1L)
  data.frame(
    start_idx = starts,
    end_idx = ends,
    start = as.integer(positions[starts]),
    end = as.integer(positions[ends]),
    effect_mean = vapply(seq_along(starts), function(i) {
      mean(score[starts[i]:ends[i]])
    }, numeric(1)),
    n_bins = ends - starts + 1L,
    penalty = penalty,
    loss = loss,
    stringsAsFactors = FALSE
  )
}

.winsorize_mad <- function(y, k) {
  m <- stats::median(y)
  s <- stats::mad(y, constant = 1.4826)
  if (!is.finite(s) || s == 0) return(y)
  lo <- m - k * s
  hi <- m + k * s
  pmin(pmax(y, lo), hi)
}

# --- Gaussian costs via prefix sums ---

.gaussian_cost_factory <- function(y) {
  S <- .prefix_sum(y)
  Q <- .prefix_sum(y * y)
  function(s, t) {
    n <- t - s + 1L
    ss <- S[t + 1L] - S[s]
    qq <- Q[t + 1L] - Q[s]
    qq - ss * ss / n
  }
}

.pelt_gaussian <- function(y, penalty, min_width) {
  n <- length(y)
  if (isTRUE(getOption("EpiSegDiff.use_cpp", TRUE))) {
    out <- tryCatch(
      pelt_gaussian_cpp(y, penalty, as.integer(min_width)),
      error = function(e) NULL
    )
    if (!is.null(out)) return(as.integer(out$changepoints))
  }
  .pelt_generic(y, penalty, min_width, .gaussian_cost_factory(y))
}

.op_gaussian <- function(y, penalty, min_width) {
  .op_generic(y, penalty, min_width, .gaussian_cost_factory(y))
}

.fpop_gaussian <- function(y, penalty, min_width) {
  # Functional pruning for 1-D Gaussian mean. Falls back to PELT if the
  # C++ symbol is unavailable. This is not a custom NB-FPOP engine.
  if (isTRUE(getOption("EpiSegDiff.use_cpp", TRUE))) {
    out <- tryCatch(
      fpop_gaussian_cpp(y, penalty, as.integer(min_width)),
      error = function(e) NULL
    )
    if (!is.null(out)) return(as.integer(out$changepoints))
  }
  .pelt_gaussian(y, penalty, min_width)
}

# --- Huber ---

.huber_mean <- function(y, k = 1.345, maxit = 50L, tol = 1e-8) {
  mu <- stats::median(y)
  s <- stats::mad(y, constant = 1.4826)
  if (!is.finite(s) || s == 0) s <- 1
  for (i in seq_len(maxit)) {
    r <- (y - mu) / s
    w <- ifelse(abs(r) <= k, 1, k / pmax(abs(r), 1e-12))
    mu_new <- stats::weighted.mean(y, w)
    if (abs(mu_new - mu) < tol) {
      mu <- mu_new
      break
    }
    mu <- mu_new
  }
  mu
}

.huber_cost <- function(y, k = 1.345) {
  mu <- .huber_mean(y, k)
  r <- y - mu
  sum(ifelse(abs(r) <= k, 0.5 * r * r, k * abs(r) - 0.5 * k * k))
}

.pelt_huber <- function(y, penalty, min_width, k) {
  cost <- function(s, t) .huber_cost(y[s:t], k)
  .pelt_generic(y, penalty, min_width, cost)
}

.op_huber <- function(y, penalty, min_width, k) {
  cost <- function(s, t) .huber_cost(y[s:t], k)
  .op_generic(y, penalty, min_width, cost)
}

# --- Generic PELT / OP ---

.pelt_generic <- function(y, penalty, min_width, cost) {
  n <- length(y)
  Fval <- rep(Inf, n + 1L)
  Fval[1] <- -penalty
  last <- integer(n)
  R <- 0L
  for (t in seq_len(n)) {
    best <- Inf
    best_tau <- 0L
    for (tau in R) {
      len <- t - tau
      if (len < min_width && t < n) next
      if (len < 1L) next
      val <- Fval[tau + 1L] + cost(tau + 1L, t) + penalty
      if (val < best) {
        best <- val
        best_tau <- tau
      }
    }
    if (!is.finite(best)) {
      best <- Fval[1] + cost(1L, t) + penalty
      best_tau <- 0L
    }
    Fval[t + 1L] <- best
    last[t] <- best_tau
    Rnew <- integer()
    for (tau in R) {
      len <- t - tau
      if (len < 1L) next
      if (Fval[tau + 1L] + cost(tau + 1L, t) <= Fval[t + 1L]) {
        Rnew <- c(Rnew, tau)
      }
    }
    R <- c(Rnew, t)
  }
  .backtrack_cps(last, n)
}

.op_generic <- function(y, penalty, min_width, cost) {
  n <- length(y)
  Fval <- rep(Inf, n + 1L)
  Fval[1] <- -penalty
  last <- integer(n)
  for (t in seq_len(n)) {
    best <- Inf
    best_tau <- 0L
    max_tau <- t - min_width
    if (t < min_width) {
      Fval[t + 1L] <- Fval[1] + cost(1L, t) + penalty
      last[t] <- 0L
      next
    }
    for (tau in 0:max_tau) {
      val <- Fval[tau + 1L] + cost(tau + 1L, t) + penalty
      if (val < best) {
        best <- val
        best_tau <- tau
      }
    }
    Fval[t + 1L] <- best
    last[t] <- best_tau
  }
  .backtrack_cps(last, n)
}

.backtrack_cps <- function(last, n) {
  cps <- integer()
  curr <- n
  while (curr > 0L) {
    tau <- last[curr]
    if (tau > 0L) cps <- c(tau, cps)
    curr <- tau
    if (length(cps) > n) stop("Changepoint backtrack failed.", call. = FALSE)
  }
  as.integer(cps)
}

#' Multi-penalty segmentation
#'
#' Stores all candidate segments. Does not claim they form a hierarchy.
#'
#' @param score Numeric vector.
#' @param penalties Numeric penalty grid.
#' @param positions Optional positions.
#' @param ... Passed to [segment_effect()].
#' @return data.frame with chrom-less segments plus `penalty`, `segment_id`.
#' @export
segment_penalty_path <- function(score, penalties, positions = NULL, ...) {
  penalties <- sort(unique(as.numeric(penalties)))
  parts <- lapply(penalties, function(lam) {
    seg <- segment_effect(score, positions = positions, penalty = lam, ...)
    if (nrow(seg) == 0L) return(seg)
    seg
  })
  out <- do.call(rbind, parts)
  if (is.null(out) || nrow(out) == 0L) {
    out <- .empty_segments()
    out$segment_id <- integer()
    return(out)
  }
  out$segment_id <- seq_len(nrow(out))
  out
}

#' Segment an effect track over chromosomes and CVRs
#'
#' @param track data.frame from [write_effect_track()].
#' @param penalty Scalar penalty.
#' @param ... Passed to [segment_effect()].
#' @return Segment table with `chrom`, genomic `start`/`end`.
#' @keywords internal
segment_track <- function(track, penalty, min_width = 3L, ...) {
  bins <- data.frame(
    chrom = track$chrom, start = track$start, end = track$end,
    valid = track$valid & is.finite(track$score),
    stringsAsFactors = FALSE
  )
  cvrs <- contiguous_valid_regions(bins)
  if (nrow(cvrs) == 0L) {
    out <- .empty_segments()
    out$chrom <- character()
    out$cvr_id <- integer()
    return(out)
  }
  parts <- lapply(seq_len(nrow(cvrs)), function(i) {
    a <- cvrs$bin_start_idx[i]
    b <- cvrs$bin_end_idx[i]
    sc <- track$score[a:b]
    pos <- track$start[a:b]
    seg <- segment_effect(sc, positions = pos, penalty = penalty,
                          min_width = min_width, ...)
    if (nrow(seg) == 0L) return(NULL)
    seg$chrom <- cvrs$chrom[i]
    seg$cvr_id <- cvrs$cvr_id[i]
    # genomic end is the end of the last bin in the segment
    last_bins <- a + seg$end_idx - 1L
    first_bins <- a + seg$start_idx - 1L
    seg$start <- track$start[first_bins]
    seg$end <- track$end[last_bins]
    seg
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) {
    out <- .empty_segments()
    out$chrom <- character()
    out$cvr_id <- integer()
    return(out)
  }
  do.call(rbind, parts)
}

#' Penalty path over a full effect track
#' @keywords internal
segment_track_penalty_path <- function(track, penalties, min_width = 3L, ...) {
  penalties <- sort(unique(as.numeric(penalties)))
  parts <- lapply(penalties, function(lam) {
    seg <- segment_track(track, penalty = lam, min_width = min_width, ...)
    if (nrow(seg) == 0L) return(NULL)
    seg
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) {
    out <- .empty_segments()
    out$chrom <- character()
    out$cvr_id <- integer()
    out$segment_id <- integer()
    return(out)
  }
  out <- do.call(rbind, parts)
  out$segment_id <- seq_len(nrow(out))
  out
}
