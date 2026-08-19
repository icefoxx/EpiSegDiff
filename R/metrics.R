#' Region discovery metrics against interval truth
#'
#' Shared matcher for EpiSegDiff and competitors.
#'
#' @param pred data.frame with chrom, start, end.
#' @param truth data.frame with chrom, start, end.
#' @param min_reciprocal_overlap Matching threshold.
#' @return list of scalar metrics and the match table.
#' @export
compute_region_metrics <- function(pred, truth, min_reciprocal_overlap = 0.2) {
  if (nrow(pred) == 0L && nrow(truth) == 0L) {
    return(list(precision = 1, recall = 1, f1 = 1, mean_iou = NA_real_,
                bp_precision = 1, bp_recall = 1, n_pred = 0L, n_truth = 0L,
                matches = .empty_match_table()))
  }
  matches <- .greedy_interval_match(pred, truth, min_reciprocal_overlap)
  tp <- sum(matches$matched)
  prec <- if (nrow(pred) == 0L) 0 else tp / nrow(pred)
  rec <- if (nrow(truth) == 0L) 0 else tp / nrow(truth)
  f1 <- if (prec + rec == 0) 0 else 2 * prec * rec / (prec + rec)
  iou <- if (any(matches$matched)) mean(matches$iou[matches$matched]) else NA_real_
  bp <- .bp_overlap(pred, truth)
  list(
    precision = prec, recall = rec, f1 = f1, mean_iou = iou,
    bp_precision = bp$precision, bp_recall = bp$recall,
    n_pred = nrow(pred), n_truth = nrow(truth),
    matches = matches
  )
}

.empty_match_table <- function() {
  data.frame(
    pred_i = integer(), truth_j = integer(), reciprocal_overlap = numeric(),
    iou = numeric(), matched = logical(), stringsAsFactors = FALSE
  )
}

.greedy_interval_match <- function(pred, truth, min_ro) {
  if (nrow(pred) == 0L || nrow(truth) == 0L) {
    return(.empty_match_table())
  }
  pairs <- list()
  for (i in seq_len(nrow(pred))) {
    for (j in seq_len(nrow(truth))) {
      if (as.character(pred$chrom[i]) != as.character(truth$chrom[j])) next
      inter <- .interval_intersection(pred$start[i], pred$end[i],
                                      truth$start[j], truth$end[j])
      if (inter <= 0L) next
      ro <- .reciprocal_overlap(pred$start[i], pred$end[i],
                                truth$start[j], truth$end[j])
      union <- (pred$end[i] - pred$start[i]) + (truth$end[j] - truth$start[j]) - inter
      iou <- if (union > 0) inter / union else 0
      if (ro >= min_ro) {
        pairs[[length(pairs) + 1L]] <- data.frame(
          pred_i = i, truth_j = j, reciprocal_overlap = ro, iou = iou,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(pairs)) return(.empty_match_table())
  cand <- do.call(rbind, pairs)
  cand <- cand[order(-cand$reciprocal_overlap, -cand$iou, cand$pred_i, cand$truth_j), ]
  used_p <- used_t <- integer()
  keep <- logical(nrow(cand))
  for (k in seq_len(nrow(cand))) {
    if (cand$pred_i[k] %in% used_p || cand$truth_j[k] %in% used_t) next
    keep[k] <- TRUE
    used_p <- c(used_p, cand$pred_i[k])
    used_t <- c(used_t, cand$truth_j[k])
  }
  cand$matched <- keep
  cand
}

.bp_overlap <- function(pred, truth) {
  # Per-chromosome discretized union via interval merge of intersections.
  pred_bp <- if (nrow(pred)) sum(pmax(0L, pred$end - pred$start)) else 0
  truth_bp <- if (nrow(truth)) sum(pmax(0L, truth$end - truth$start)) else 0
  inter_bp <- 0
  chrs <- union(unique(as.character(pred$chrom)), unique(as.character(truth$chrom)))
  for (chr in chrs) {
    p <- pred[as.character(pred$chrom) == chr, , drop = FALSE]
    t <- truth[as.character(truth$chrom) == chr, , drop = FALSE]
    if (!nrow(p) || !nrow(t)) next
    for (i in seq_len(nrow(p))) {
      inter_bp <- inter_bp + sum(.interval_intersection(
        p$start[i], p$end[i], t$start, t$end
      ))
    }
  }
  # Intersection counted per pred×truth pair can double-count; this is an
  # upper bound used as a documented approximation unless intervals are disjoint.
  list(
    precision = if (pred_bp == 0) 0 else min(1, inter_bp / pred_bp),
    recall = if (truth_bp == 0) 0 else min(1, inter_bp / truth_bp)
  )
}

#' Boundary error metrics for matched intervals
#'
#' @param matches Match table from [compute_region_metrics()] (must include
#'   coordinates if `pred`/`truth` are supplied).
#' @param pred,truth Interval tables.
#' @return list of MAE / median / p90 / Hausdorff.
#' @export
compute_boundary_metrics <- function(pred, truth, matches = NULL) {
  if (is.null(matches)) {
    matches <- compute_region_metrics(pred, truth)$matches
  }
  m <- matches[matches$matched, , drop = FALSE]
  if (nrow(m) == 0L) {
    return(list(
      start_mae = NA_real_, end_mae = NA_real_,
      start_median = NA_real_, end_median = NA_real_,
      start_p90 = NA_real_, end_p90 = NA_real_,
      hausdorff = NA_real_, n_matched = 0L
    ))
  }
  ds <- abs(pred$start[m$pred_i] - truth$start[m$truth_j])
  de <- abs(pred$end[m$pred_i] - truth$end[m$truth_j])
  pred_b <- c(pred$start, pred$end)
  truth_b <- c(truth$start, truth$end)
  haus <- max(
    max(vapply(pred_b, function(x) min(abs(x - truth_b)), numeric(1))),
    max(vapply(truth_b, function(x) min(abs(x - pred_b)), numeric(1)))
  )
  q90 <- function(x) as.numeric(stats::quantile(x, 0.9, names = FALSE, type = 7))
  list(
    start_mae = mean(ds), end_mae = mean(de),
    start_median = stats::median(ds), end_median = stats::median(de),
    start_p90 = q90(ds), end_p90 = q90(de),
    hausdorff = haus, n_matched = nrow(m)
  )
}

#' Event-class confusion matrix and macro-F1
#'
#' @param truth_class Character vector.
#' @param pred_class Character vector (same length).
#' @param classes Optional class set.
#' @export
event_class_metrics <- function(truth_class, pred_class, classes = NULL) {
  truth_class <- as.character(truth_class)
  pred_class <- as.character(pred_class)
  if (length(truth_class) != length(pred_class)) {
    stop("truth_class and pred_class must have the same length.", call. = FALSE)
  }
  if (is.null(classes)) {
    classes <- sort(unique(c(truth_class, pred_class)))
  }
  tab <- table(
    factor(truth_class, levels = classes),
    factor(pred_class, levels = classes)
  )
  f1 <- vapply(classes, function(cl) {
    tp <- tab[cl, cl]
    fp <- sum(tab[, cl]) - tp
    fn <- sum(tab[cl, ]) - tp
    prec <- if (tp + fp == 0) 0 else tp / (tp + fp)
    rec <- if (tp + fn == 0) 0 else tp / (tp + fn)
    if (prec + rec == 0) 0 else 2 * prec * rec / (prec + rec)
  }, numeric(1))
  acc <- if (!length(truth_class)) NA_real_ else mean(truth_class == pred_class)
  list(
    confusion = tab,
    per_class_f1 = f1,
    macro_f1 = mean(f1),
    accuracy = acc,
    n = length(truth_class)
  )
}
