#' Infer DERs from segments using limma on region means
#' @export
infer_regions <- function(seg, signal, bed, samples, numerator, denominator,
                          alpha = 0.05, mode = "same_data") {
  if (!nrow(seg)) {
    return(data.frame(chrom = character(), start = integer(), end = integer(),
                      effect = numeric(), direction = character(),
                      level1 = character(), FDR = numeric(),
                      inference_mode = character(), stringsAsFactors = FALSE))
  }
  track <- data.frame(chrom = bed$chrom, start = bed$start, end = bed$end,
                      score = 0, valid = TRUE)
  segs <- data.frame(chrom = bed$chrom[seg$start_idx],
                     start = bed$start[seg$start_idx],
                     end = bed$end[seg$end_idx],
                     region_id = seq_len(nrow(seg)))
  xf <- transform_continuous(signal, "log2_pseudocount", 0.01)
  inf <- tryCatch(infer_capability_level(samples), error = function(e) "descriptive")
  tested <- test_regions_limma(segs, xf, track, samples, numerator, denominator,
                               inference_mode = inf, fdr_threshold = alpha)
  if (identical(inf, "descriptive")) {
    tested$level1 <- ifelse(abs(tested$effect) >= 0.5,
                            ifelse(tested$effect > 0, "GAIN", "LOSS"), "STABLE")
    tested$FDR <- ifelse(tested$level1 == "STABLE", 1, 0.05)
  }
  tested <- tested[tested$level1 %in% c("GAIN", "LOSS"), , drop = FALSE]
  if (nrow(tested)) tested$inference_mode <- mode
  tested
}

#' Multi-scale candidate union with overlap resolution
#' @export
union_resolve_calls <- function(call_list) {
  allc <- Filter(function(d) is.data.frame(d) && nrow(d) > 0, call_list)
  if (!length(allc)) {
    return(data.frame(chrom = character(), start = integer(), end = integer(),
                      effect = numeric(), direction = character(),
                      level1 = character(), FDR = numeric(),
                      stringsAsFactors = FALSE))
  }
  d <- do.call(rbind, lapply(seq_along(allc), function(i) {
    x <- allc[[i]]
    x$source_scale <- names(allc)[i]
    x
  }))
  d <- d[order(-abs(d$effect), d$FDR), ]
  keep <- rep(TRUE, nrow(d))
  if (nrow(d) > 1L) {
    g <- GenomicRanges::GRanges(d$chrom, IRanges::IRanges(d$start + 1L, d$end))
    for (i in seq_len(nrow(d))) {
      if (!keep[i]) next
      ov <- GenomicRanges::findOverlaps(g[i], g)
      jj <- S4Vectors::subjectHits(ov)
      jj <- jj[jj > i & keep[jj]]
      if (length(jj)) keep[jj] <- FALSE
    }
  }
  d[keep, , drop = FALSE]
}

#' Coarse-to-fine: segment coarse score, refine edges on fine score
#' @export
hierarchical_refine <- function(fine_score, fine_bed, coarse_factor,
                                engine, penalty, min_width) {
  n <- length(fine_score)
  g <- ceiling(seq_len(n) / coarse_factor)
  cs <- as.numeric(tapply(fine_score, g, mean))
  cpos <- tapply(fine_bed$start, g, min)
  segc <- segment_signal(cs, positions = as.integer(cpos), engine = engine,
                         penalty = penalty, min_width = max(1L, min_width))
  if (!nrow(segc)) return(segc)
  # map coarse idx to fine idx
  fine_start <- tapply(seq_len(n), g, min)
  fine_end <- tapply(seq_len(n), g, max)
  segf <- segc
  segf$start_idx <- as.integer(fine_start[segc$start_idx])
  segf$end_idx <- as.integer(fine_end[segc$end_idx])
  refine_segment_edges(fine_score, segf, window = max(3L, coarse_factor))
}

#' Route a 1-D score to a production candidate (no truth used)
#'
#' High first-difference roughness relative to the score MAD selects
#' fused-lasso / standardized (narrow-sensitive); otherwise PELT / raw.
#' Threshold frozen before held-out evaluation.
#' @export
route_production_engine <- function(score, roughness_cut = 1.2) {
  score <- as.numeric(score)
  m <- stats::mad(score, constant = 1)
  md <- stats::mad(diff(score), constant = 1)
  rough <- if (is.finite(m) && m > 0 && is.finite(md)) md / m else 0
  if (isTRUE(rough > roughness_cut)) {
    list(score_engine = "standardized", segmentation_engine = "fused_lasso",
         min_width = 1L, roughness = rough, route = "narrow")
  } else {
    list(score_engine = "raw_effect", segmentation_engine = "pelt",
         min_width = 3L, roughness = rough, route = "broad")
  }
}

#' Compare several score × segmentation engines on one simulation
#' @export
compare_engines <- function(sim, score_engines, seg_engines, min_width = 1L,
                            penalty = 8) {
  rows <- list()
  for (sc in score_engines) {
    for (se in seg_engines) {
      mw <- if (se %in% c("pelt", "fpop", "hmm", "hsmm")) max(min_width, 3L) else min_width
      row <- tryCatch({
        score <- compute_diff_score(sim$signal, sim$samples, "JNKi", "DMSO", sc)
        seg <- segment_signal(score, positions = sim$bed$start, engine = se,
                              penalty = penalty, min_width = mw)
        thr <- 0.15 + 0.35 * stats::mad(score, constant = 1)
        if (!is.finite(thr) || thr <= 0) thr <- 0.2
        calls <- .segments_to_calls(seg, sim$bed, thresh = thr)
        sc0 <- score_calls_vs_truth(calls, sim$truth, sim$bin_size)
        cbind(score_engine = sc, segmentation_engine = se, min_width = mw, sc0)
      }, error = function(e) NULL)
      if (!is.null(row)) rows[[length(rows) + 1L]] <- row
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

#' Run one engine configuration on a simulated series
#' @export
run_multiscale <- function(sim, score_engine = "raw_effect",
                           segmentation_engine = "pelt",
                           scale_mode = "single",
                           min_width = 1L, penalty = 8,
                           inference_mode = "same_data",
                           alpha = 0.05) {
  num <- "JNKi"; den <- "DMSO"
  sc <- compute_diff_score(sim$signal, sim$samples, num, den, score_engine)
  if (inference_mode == "pooled") {
    sc <- rowMeans(sim$signal)
  } else if (inference_mode == "split") {
    # even samples for score
    ev <- sim$samples$sample[c(TRUE, FALSE)]
    sc <- compute_diff_score(sim$signal[, ev, drop = FALSE],
                             sim$samples[sim$samples$sample %in% ev, ],
                             num, den, "raw_effect")
  }
  if (scale_mode == "coarse_fine") {
    seg <- hierarchical_refine(sc, sim$bed, 5L, segmentation_engine, penalty, min_width)
  } else {
    seg <- segment_signal(sc, positions = sim$bed$start, engine = segmentation_engine,
                          penalty = penalty, min_width = min_width)
  }
  if (inference_mode == "split") {
    odd <- sim$samples$sample[c(FALSE, TRUE)]
    tes <- infer_regions(seg, sim$signal[, odd, drop = FALSE],
                         sim$bed, sim$samples[sim$samples$sample %in% odd, ],
                         num, den, alpha, "split")
  } else {
    tes <- infer_regions(seg, sim$signal, sim$bed, sim$samples, num, den,
                         alpha, inference_mode)
  }
  if (nrow(tes)) {
    tes$score_engine <- score_engine
    tes$segmentation_engine <- segmentation_engine
    tes$scale_mode <- scale_mode
    tes$inference_mode <- inference_mode
  }
  tes
}

#' Cross-fit: A segments / B tests and reverse; keep regions passing either held-out test
#' @export
run_crossfit <- function(sim, score_engine = "raw_effect",
                         segmentation_engine = "wbs",
                         penalty = 8, min_width = 1L, alpha = 0.05) {
  s1 <- sim$samples$sample[c(1, 3)]
  s2 <- sim$samples$sample[c(2, 4)]
  go <- function(train, test) {
    sc <- compute_diff_score(sim$signal[, train, drop = FALSE],
                             sim$samples[sim$samples$sample %in% train, ],
                             "JNKi", "DMSO", score_engine)
    seg <- segment_signal(sc, positions = sim$bed$start,
                          engine = segmentation_engine, penalty = penalty,
                          min_width = min_width)
    infer_regions(seg, sim$signal[, test, drop = FALSE], sim$bed,
                  sim$samples[sim$samples$sample %in% test, ],
                  "JNKi", "DMSO", alpha, "crossfit")
  }
  a <- tryCatch(go(s1, s2), error = function(e) NULL)
  b <- tryCatch(go(s2, s1), error = function(e) NULL)
  union_resolve_calls(list(a = a, b = b))
}
