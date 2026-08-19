#' Simulate a 4-sample continuous series with planted events
#' @export
simulate_continuous_events <- function(n_bins, bin_size = 200L, events,
                                       baseline = 0.12, noise_sd = 0.2,
                                       seed = 1L, background = "flat") {
  set.seed(as.integer(seed))
  n_bins <- as.integer(n_bins)
  mat <- matrix(baseline, n_bins, 4L)
  colnames(mat) <- c("d1", "d2", "j1", "j2")
  truth <- events
  for (k in seq_len(nrow(events))) {
    a <- events$start_bin[k]; b <- events$end_bin[k]
    delta <- events$delta[k]
    shape <- if ("shape" %in% names(events)) events$shape[k] else "sharp"
    w <- b - a + 1L
    ramp <- if (identical(shape, "gradual") && w >= 5L) {
      t1 <- max(1L, floor(0.2 * w))
      wt <- c(seq(0, 1, length.out = t1),
              rep(1, max(0L, w - 2L * t1)),
              seq(1, 0, length.out = t1))
      if (length(wt) < w) wt <- c(wt, rep(0, w - length(wt)))
      wt[seq_len(w)]
    } else {
      rep(1, w)
    }
    if (delta > 0) mat[a:b, 3:4] <- mat[a:b, 3:4] + delta * ramp
    else mat[a:b, 1:2] <- mat[a:b, 1:2] + abs(delta) * ramp
  }
  if (identical(background, "hetero")) {
    noise_sd <- noise_sd * (0.5 + seq_len(n_bins) / n_bins)
  } else if (identical(background, "elevated")) {
    mid <- seq.int(floor(n_bins * 0.4), floor(n_bins * 0.6))
    mat[mid, ] <- mat[mid, ] + 0.4
  } else if (identical(background, "zeros")) {
    z <- seq.int(1L, floor(n_bins * 0.35))
    mat[z, ] <- 0
  }
  mat <- pmax(mat + matrix(stats::rnorm(n_bins * 4L, 0, noise_sd), n_bins, 4L), 0)
  bed <- data.frame(chrom = "1",
                    start = as.integer((seq_len(n_bins) - 1L) * bin_size),
                    end = as.integer(seq_len(n_bins) * bin_size),
                    stringsAsFactors = FALSE)
  samp <- data.frame(sample = colnames(mat),
                     condition = c("DMSO", "DMSO", "JNKi", "JNKi"))
  list(bed = bed, signal = mat, samples = samp, truth = truth,
       bin_size = as.integer(bin_size))
}

#' Merge adjacent same-sign FDR bins into windows
#' @export
merge_significant_bins <- function(bed, effect, fdr, valid, alpha = 0.05) {
  call <- valid & is.finite(fdr) & fdr <= alpha & is.finite(effect) & effect != 0
  if (!any(call)) {
    return(data.frame(chrom = character(), start = integer(), end = integer(),
                      effect = numeric(), direction = character(),
                      level1 = character(), FDR = numeric(),
                      stringsAsFactors = FALSE))
  }
  sgn <- sign(effect)
  br <- c(TRUE, call[-1] != call[-length(call)] |
            (call[-1] & sgn[-1] != sgn[-length(sgn)]) |
            bed$chrom[-1] != bed$chrom[-length(bed$chrom)] |
            (call[-1] & bed$start[-1] != bed$end[-length(bed$end)]))
  id <- cumsum(br)
  keep <- call
  out <- do.call(rbind, lapply(split(which(keep), id[keep]), function(ii) {
    data.frame(chrom = bed$chrom[ii[1]], start = bed$start[ii[1]],
               end = bed$end[ii[length(ii)]],
               effect = mean(effect[ii]), 
               direction = if (mean(effect[ii]) > 0) "GAIN" else "LOSS",
               level1 = if (mean(effect[ii]) > 0) "GAIN" else "LOSS",
               FDR = min(fdr[ii], na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  out
}

#' Score calls against planted truth events
#' @export
score_calls_vs_truth <- function(calls, truth, bin_size) {
  if (!nrow(truth)) {
    return(data.frame(n_truth = 0, n_call = nrow(calls),
                      n_detected = 0, n_merged = 0, n_split = 0,
                      bp_tp = 0, bp_pred = 0, bp_tru = 0,
                      precision = NA, recall = NA, f1 = NA,
                      dir_acc = NA, opp_err = NA,
                      med_bound_err = NA, p95_bound_err = NA,
                      mean_iou = NA, one_call_frac = NA,
                      stringsAsFactors = FALSE))
  }
  tgr <- GenomicRanges::GRanges(
    "1", IRanges::IRanges(as.integer((truth$start_bin - 1L) * bin_size + 1L),
                          as.integer(truth$end_bin * bin_size))
  )
  if (!nrow(calls)) {
    tru_bp <- sum(GenomicRanges::width(tgr))
    return(data.frame(n_truth = nrow(truth), n_call = 0, n_detected = 0,
                      n_merged = 0, n_split = 0, bp_tp = 0, bp_pred = 0,
                      bp_tru = tru_bp, precision = 0, recall = 0, f1 = 0,
                      dir_acc = NA, opp_err = NA, med_bound_err = NA,
                      p95_bound_err = NA, mean_iou = NA, one_call_frac = 0,
                      stringsAsFactors = FALSE))
  }
  cgr <- GenomicRanges::GRanges(
    calls$chrom, IRanges::IRanges(as.integer(calls$start) + 1L,
                                  as.integer(calls$end))
  )
  ov <- GenomicRanges::findOverlaps(tgr, cgr)
  n_det <- sum(GenomicRanges::countOverlaps(tgr, cgr) > 0)
  n_split <- sum(GenomicRanges::countOverlaps(tgr, cgr) >= 2L)
  n_merged <- sum(GenomicRanges::countOverlaps(cgr, tgr) >= 2L)
  inter <- if (length(cgr)) GenomicRanges::intersect(cgr, tgr) else
    GenomicRanges::GRanges()
  tp <- sum(as.numeric(GenomicRanges::width(inter)))
  pred <- sum(as.numeric(GenomicRanges::width(cgr)))
  tru <- sum(as.numeric(GenomicRanges::width(tgr)))
  prec <- if (pred) tp / pred else 0
  rec <- if (tru) tp / tru else 0
  f1 <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  # direction of best-IoU call
  dir_ok <- opp <- 0
  berr <- iou <- numeric()
  for (i in seq_len(nrow(truth))) {
    hits <- S4Vectors::subjectHits(ov)[S4Vectors::queryHits(ov) == i]
    want <- if (truth$delta[i] > 0) "GAIN" else "LOSS"
    if (!length(hits)) next
    ww <- vapply(hits, function(j) {
      sum(as.numeric(GenomicRanges::width(
        GenomicRanges::intersect(tgr[i], cgr[j]))))
    }, numeric(1))
    j <- hits[which.max(ww / pmax(as.numeric(GenomicRanges::width(cgr[hits])), 1))]
    got <- if ("direction" %in% names(calls)) calls$direction[j] else
      ifelse(calls$effect[j] > 0, "GAIN", "LOSS")
    if (identical(got, want)) dir_ok <- dir_ok + 1L else opp <- opp + 1L
    berr <- c(berr, abs(calls$start[j] - (truth$start_bin[i] - 1L) * bin_size) +
                abs(calls$end[j] - truth$end_bin[i] * bin_size))
    uni <- sum(as.numeric(GenomicRanges::width(
      GenomicRanges::reduce(c(tgr[i], cgr[j])))))
    iou <- c(iou, if (uni > 0) max(ww) / uni else 0)
  }
  data.frame(
    n_truth = nrow(truth), n_call = nrow(calls), n_detected = n_det,
    n_merged = n_merged, n_split = n_split,
    bp_tp = tp, bp_pred = pred, bp_tru = tru,
    precision = prec, recall = rec, f1 = f1,
    dir_acc = if ((dir_ok + opp) > 0) dir_ok / (dir_ok + opp) else NA,
    opp_err = if ((dir_ok + opp) > 0) opp / (dir_ok + opp) else NA,
    med_bound_err = if (length(berr)) stats::median(berr) else NA,
    p95_bound_err = if (length(berr)) as.numeric(stats::quantile(berr, 0.95)) else NA,
    mean_iou = if (length(iou)) mean(iou) else NA,
    one_call_frac = mean(GenomicRanges::countOverlaps(tgr, cgr) == 1L),
    stringsAsFactors = FALSE
  )
}

#' Run ESD + window-aggregate on a simulated series
#' @export
run_two_methods_on_sim <- function(sim, alpha = 0.05, penalty = 8,
                                   min_bins = 3L) {
  fit <- analyze_continuous_matrix(
    sim$bed, sim$signal, sim$samples, "JNKi", "DMSO",
    bin_size = sim$bin_size, filter_method = "none",
    penalty = penalty, min_bins = min_bins, engine = "pelt",
    fdr_threshold = alpha
  )
  esd <- fit$tested[fit$tested$level1 %in% c("GAIN", "LOSS"), ]
  xf <- transform_continuous(sim$signal, "log2_pseudocount", 0.01)
  eff <- estimate_continuous_effect(xf, sim$samples, "JNKi", "DMSO", "replicated")
  fdr <- p.adjust(eff$pvalue, "BH")
  win <- merge_significant_bins(sim$bed, eff$effect, fdr, rep(TRUE, nrow(sim$bed)),
                                alpha)
  list(esd = esd, window = win, fit = fit)
}
