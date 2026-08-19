#' Differential score engines (v0.7)
#'
#' @param engine One of `raw_effect`, `log_effect`, `moderated_t`,
#'   `standardized`, `signed_pval`, `hotelling`.
#' @export
compute_diff_score <- function(signal, samples, numerator, denominator,
                               engine = c("raw_effect", "log_effect",
                                          "moderated_t", "standardized",
                                          "signed_pval", "hotelling"),
                               pseudocount = 0.01) {
  engine <- match.arg(engine)
  a <- signal[, samples$condition == denominator, drop = FALSE]
  b <- signal[, samples$condition == numerator, drop = FALSE]
  ma <- rowMeans(a)
  mb <- rowMeans(b)
  na <- ncol(a)
  nb <- ncol(b)
  va <- if (na > 1L) apply(a, 1L, stats::var) else rep(0, nrow(signal))
  vb <- if (nb > 1L) apply(b, 1L, stats::var) else rep(0, nrow(signal))
  sp2 <- ((na - 1) * va + (nb - 1) * vb) / pmax(na + nb - 2, 1)
  raw <- mb - ma
  out <- switch(
    engine,
    raw_effect = raw,
    log_effect = log2(pmax(mb, 0) + pseudocount) - log2(pmax(ma, 0) + pseudocount),
    moderated_t = {
      xf <- transform_continuous(signal, "log2_pseudocount", pseudocount)
      estimate_continuous_effect(xf, samples, numerator, denominator, "replicated")$t
    },
    standardized = {
      # variance-standardized effect: not the v0.4 se-floor z
      tau2 <- stats::mad(raw, constant = 1)^2
      raw / sqrt(pmax(sp2, 0) + tau2 + 1e-8)
    },
    signed_pval = {
      xf <- transform_continuous(signal, "log2_pseudocount", pseudocount)
      ef <- estimate_continuous_effect(xf, samples, numerator, denominator, "replicated")
      sign(ef$effect) * (-log10(pmax(ef$pvalue, 1e-16)))
    },
    hotelling = {
      # 1-D Hotelling T^2 with sign of the mean difference
      t2 <- (na * nb / pmax(na + nb, 1)) * (raw^2) / pmax(sp2, 1e-8)
      sign(raw) * sqrt(pmax(t2, 0))
    }
  )
  attr(out, "score_engine") <- engine
  out
}

#' 1-D total-variation / fused-lasso denoise (dual projected gradient)
.tv_denoise <- function(y, lambda, niter = 80L) {
  y <- as.numeric(y)
  n <- length(y)
  if (n < 2L) return(y)
  z <- numeric(n - 1L)
  step <- 0.25
  for (it in seq_len(niter)) {
    g <- c(z, 0) - c(0, z)
    x <- y - g
    d <- diff(x)
    z <- pmin(pmax(z + step * d, -lambda), lambda)
  }
  g <- c(z, 0) - c(0, z)
  y - g
}

.cps_from_jumps <- function(x, min_jump = NULL) {
  d <- diff(as.numeric(x))
  if (is.null(min_jump)) min_jump <- stats::mad(d, constant = 1) * 2
  if (!is.finite(min_jump) || min_jump <= 0) min_jump <- 1e-4
  which(abs(d) > min_jump)
}

.cusum_best <- function(y) {
  y <- as.numeric(y)
  n <- length(y)
  if (n < 4L) return(list(k = NA_integer_, stat = 0))
  m <- mean(y)
  cs <- cumsum(y - m)
  cs <- cs[-n]
  k <- which.max(abs(cs))
  list(k = k, stat = abs(cs[k]) / sqrt(n))
}

.binseg_cps <- function(y, penalty, min_width = 1L) {
  rec <- function(s, e, acc) {
    if (e - s + 1L < 2L * min_width + 1L) return(acc)
    sub <- y[s:e]
    cs <- .cusum_best(sub)
    if (!is.finite(cs$stat) || cs$stat^2 < penalty) return(acc)
    cp <- s + cs$k - 1L
    if (cp - s + 1L < min_width || e - cp < min_width) return(acc)
    acc <- c(acc, cp)
    acc <- rec(s, cp, acc)
    rec(cp + 1L, e, acc)
  }
  sort(unique(rec(1L, length(y), integer())))
}

.wbs_cps <- function(y, penalty, min_width = 1L, M = 200L) {
  n <- length(y)
  cands <- integer()
  stats <- numeric()
  for (m in seq_len(M)) {
    if (n < 2L * min_width + 3L) break
    maxs <- n - 2L * min_width - 1L
    if (maxs < 1L) break
    s <- sample.int(maxs, 1L)
    rem <- n - s - 2L * min_width
    if (rem < 1L) next
    e <- s + 2L * min_width + sample.int(rem, 1L) - 1L
    e <- min(max(e, s + 2L), n)
    if (e - s + 1L < 2L * min_width + 1L) next
    cs <- .cusum_best(y[s:e])
    if (is.finite(cs$stat) && cs$stat^2 >= penalty) {
      cands <- c(cands, s + cs$k - 1L)
      stats <- c(stats, cs$stat)
    }
  }
  if (!length(cands)) return(integer())
  o <- order(stats, decreasing = TRUE)
  keep <- integer()
  for (i in o) {
    cp <- cands[i]
    if (all(abs(keep - cp) >= min_width)) keep <- c(keep, cp)
  }
  sort(keep)
}

.hmm3_cps <- function(y, min_width = 1L) {
  y <- as.numeric(y)
  n <- length(y)
  s <- stats::mad(y, constant = 1.4826)
  if (!is.finite(s) || s <= 0) s <- stats::sd(y)
  if (!is.finite(s) || s <= 0) s <- 1
  mu <- c(-2, 0, 2) * s
  logp <- sapply(mu, function(m) -(y - m)^2 / (2 * s^2))
  # Viterbi, stay-biased transitions
  stay <- log(0.98)
  sw <- log(0.01)
  dp <- matrix(-Inf, n, 3L)
  ptr <- matrix(1L, n, 3L)
  dp[1, ] <- logp[1, ]
  for (i in 2:n) {
    for (k in 1:3) {
      trans <- dp[i - 1L, ] + ifelse(1:3 == k, stay, sw)
      ptr[i, k] <- which.max(trans)
      dp[i, k] <- max(trans) + logp[i, k]
    }
  }
  st <- integer(n)
  st[n] <- which.max(dp[n, ])
  for (i in n:2) st[i - 1L] <- ptr[i, st[i]]
  cps <- which(st[-1] != st[-n])
  cps[cps >= min_width & (n - cps) >= min_width]
}

.hsmm_cps <- function(y, min_width = 3L) {
  # duration-aware: HMM then drop short sojourns by merging to neighbor
  cps <- .hmm3_cps(y, min_width = 1L)
  if (!length(cps)) return(integer())
  ends <- c(cps, length(y))
  starts <- c(1L, cps + 1L)
  keep <- (ends - starts + 1L) >= min_width
  if (all(keep)) return(cps)
  # merge short pieces: drop their left boundary
  drop <- which(!keep)
  drop_cp <- unique(c(starts[drop] - 1L, ends[drop]))
  drop_cp <- drop_cp[drop_cp > 0 & drop_cp < length(y)]
  setdiff(cps, drop_cp)
}

#' Unified segmentation interface
#'
#' @param engine `pelt`, `binseg`, `wbs`, `fpop`, `fused_lasso`, `hmm`, `hsmm`.
#' @export
segment_signal <- function(score, positions = NULL, engine = "pelt",
                           penalty = 8, min_width = 1L, lambda = NULL) {
  engine <- match.arg(engine, c("pelt", "binseg", "wbs", "fpop",
                                "fused_lasso", "hmm", "hsmm"))
  score <- as.numeric(score)
  n <- length(score)
  if (is.null(positions)) positions <- seq_len(n)
  min_width <- as.integer(min_width)
  if (n < 2L) return(.empty_segments())
  if (engine %in% c("pelt", "fpop")) {
    return(segment_effect(score, positions, penalty = penalty,
                          min_width = min_width, engine = engine))
  }
  cps <- switch(
    engine,
    binseg = .binseg_cps(score, penalty, min_width),
    wbs = .wbs_cps(score, penalty, min_width),
    fused_lasso = {
      lam <- if (is.null(lambda)) stats::mad(diff(score), constant = 1) * 2 else lambda
      xh <- .tv_denoise(score, lam)
      .cps_from_jumps(xh)
    },
    hmm = .hmm3_cps(score, min_width),
    hsmm = .hsmm_cps(score, max(min_width, 3L))
  )
  cps <- sort(unique(as.integer(cps)))
  cps <- cps[cps >= 1L & cps < n]
  seg <- .segments_from_cps(cps, score, positions, penalty, engine)
  if (engine %in% c("hmm", "hsmm", "binseg", "wbs") && min_width > 1L) {
    seg <- .enforce_min_bins(seg, score, positions, min_width, penalty, engine)
  }
  seg$segmentation_engine <- engine
  seg
}

#' Calls from a piecewise-constant fit: keep pieces whose mean exceeds `thresh`
.segments_to_calls <- function(seg, bed, thresh = 0, fdr = 0.05) {
  if (!nrow(seg)) {
    return(data.frame(chrom = character(), start = integer(), end = integer(),
                      effect = numeric(), direction = character(),
                      level1 = character(), FDR = numeric(),
                      n_bins = integer(), estimated_scale_bp = integer(),
                      stringsAsFactors = FALSE))
  }
  i1 <- seg$start_idx
  i2 <- seg$end_idx
  out <- data.frame(
    chrom = bed$chrom[i1],
    start = bed$start[i1],
    end = bed$end[i2],
    effect = seg$effect_mean,
    direction = ifelse(seg$effect_mean > 0, "GAIN", "LOSS"),
    level1 = ifelse(abs(seg$effect_mean) > thresh,
                    ifelse(seg$effect_mean > 0, "GAIN", "LOSS"), "STABLE"),
    FDR = ifelse(abs(seg$effect_mean) > thresh, 0.01, 1),
    n_bins = seg$n_bins,
    estimated_scale_bp = bed$end[i2] - bed$start[i1],
    stringsAsFactors = FALSE
  )
  out[out$level1 %in% c("GAIN", "LOSS"), , drop = FALSE]
}

#' Optional local CUSUM refinement of segment edges (v0.7)
#'
#' Distinct from [refine_boundaries()], which operates on BED regions
#' and an effect track.
#' @export
refine_segment_edges <- function(score, seg, window = 5L) {
  if (!nrow(seg)) return(seg)
  n <- length(score)
  window <- as.integer(window)
  for (i in seq_len(nrow(seg))) {
    if (i == 1L) next
    cp <- seg$start_idx[i] - 1L
    lo <- max(1L, cp - window)
    hi <- min(n, cp + window)
    if (hi - lo < 3L) next
    cs <- .cusum_best(score[lo:hi])
    if (is.finite(cs$k)) {
      newcp <- lo + cs$k - 1L
      seg$start_idx[i] <- newcp + 1L
      seg$end_idx[i - 1L] <- newcp
    }
  }
  ok <- seg$end_idx >= seg$start_idx
  seg[ok, , drop = FALSE]
}
