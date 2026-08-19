#' Evaluate nesting among multi-penalty segments
#'
#' Does not force a tree. Partial overlaps are flagged as conflicts.
#'
#' @param segments data.frame with `segment_id`, `chrom`, `start`, `end`,
#'   and typically `penalty`.
#' @return The segments plus `parent_id`, `nesting_status`, `conflict_status`.
#' @export
analyze_nesting <- function(segments) {
  need <- c("segment_id", "chrom", "start", "end")
  if (!all(need %in% names(segments))) {
    stop("segments must contain segment_id, chrom, start, end.", call. = FALSE)
  }
  n <- nrow(segments)
  if (n == 0L) {
    segments$parent_id <- integer()
    segments$nesting_status <- character()
    segments$conflict_status <- character()
    return(segments)
  }
  # Prefer larger (coarser) intervals as parents. If penalty is present,
  # higher penalty is coarser and preferred.
  width <- segments$end - segments$start
  coarse <- if ("penalty" %in% names(segments)) {
    order(-segments$penalty, -width, segments$segment_id)
  } else {
    order(-width, segments$segment_id)
  }

  parent <- rep(NA_integer_, n)
  status <- rep("root", n)
  conflict <- rep("none", n)

  # Map segment_id to row
  id <- as.integer(segments$segment_id)
  for (ii in seq_along(coarse)) {
    i <- coarse[ii]
    cand_parent <- NA_integer_
    cand_width <- Inf
    for (jj in seq_len(ii - 1L)) {
      j <- coarse[jj]
      if (segments$chrom[i] != segments$chrom[j]) next
      rel <- .nesting_relation(segments$start[i], segments$end[i],
                               segments$start[j], segments$end[j])
      if (rel == "equal") {
        # Same interval at another scale: parent is the coarser equal.
        if (is.na(parent[i])) {
          parent[i] <- id[j]
          status[i] <- "equal_scale_duplicate"
        }
      } else if (rel == "i_in_j") {
        wj <- width[j]
        if (wj < cand_width) {
          cand_width <- wj
          cand_parent <- id[j]
        }
      } else if (rel == "partial") {
        conflict[i] <- "partial_overlap"
        conflict[j] <- "partial_overlap"
      }
    }
    if (!is.na(cand_parent) && status[i] != "equal_scale_duplicate") {
      parent[i] <- cand_parent
      status[i] <- "nested"
    }
  }

  # Any partial overlap remains flagged; we do not drop or rewrite intervals.
  segments$parent_id <- parent
  segments$nesting_status <- status
  segments$conflict_status <- conflict
  segments
}

.nesting_relation <- function(is, ie, js, je) {
  if (is == js && ie == je) return("equal")
  if (is >= js && ie <= je) return("i_in_j")
  if (js >= is && je <= ie) return("j_in_i")
  inter <- .interval_intersection(is, ie, js, je)
  if (inter > 0L) return("partial")
  "disjoint"
}

#' Select candidate intervals from a penalty path
#'
#' Two strategies: a single selected penalty, or a multi-scale consensus
#' that keeps coarsest intervals and nested children that change sign or
#' exceed a mean-effect gap.
#'
#' @param segments Output of a penalty-path segmentation (with `penalty`).
#' @param strategy `single_selected_penalty` or `multi-scale_consensus`.
#' @param selected_penalty Penalty used by the single-penalty strategy.
#' @param effect_gap Minimum |mean effect| difference to keep a nested child.
#' @return Subset of segments labeled with `candidate_source`.
#' @export
select_candidate_segments <- function(segments,
                                      strategy = c("single_selected_penalty",
                                                   "multi-scale_consensus"),
                                      selected_penalty = NULL,
                                      effect_gap = 0.5) {
  strategy <- match.arg(strategy)
  if (nrow(segments) == 0L) {
    segments$candidate_source <- character()
    return(segments)
  }
  if (strategy == "single_selected_penalty") {
    if (is.null(selected_penalty)) {
      stop("selected_penalty is required for single_selected_penalty.", call. = FALSE)
    }
    keep <- abs(segments$penalty - selected_penalty) < 1e-9
    if (!any(keep)) {
      # nearest penalty
      pen <- unique(segments$penalty)
      chosen <- pen[which.min(abs(pen - selected_penalty))]
      warning(sprintf("selected_penalty %s not in path; using nearest %s.",
                      selected_penalty, chosen), call. = FALSE)
      keep <- abs(segments$penalty - chosen) < 1e-9
    }
    out <- segments[keep, , drop = FALSE]
    out$candidate_source <- sprintf("single_selected_penalty:%s", selected_penalty)
    return(out)
  }

  nest <- if (!"parent_id" %in% names(segments)) analyze_nesting(segments) else segments
  # Consensus: keep roots at the coarsest penalty that still has ≥1 segment,
  # plus nested children whose effect_mean differs from the parent by effect_gap.
  max_pen <- max(nest$penalty)
  keep <- nest$penalty == max_pen
  if ("effect_mean" %in% names(nest) && "segment_id" %in% names(nest)) {
    em <- stats::setNames(nest$effect_mean, nest$segment_id)
    for (i in seq_len(nrow(nest))) {
      if (is.na(nest$parent_id[i])) next
      pid <- as.character(nest$parent_id[i])
      if (!pid %in% names(em)) next
      if (is.finite(nest$effect_mean[i]) && is.finite(em[[pid]]) &&
          abs(nest$effect_mean[i] - em[[pid]]) >= effect_gap) {
        keep[i] <- TRUE
      }
    }
  }
  # Do not drop conflict-flagged intervals silently; keep them for inspection
  # but mark source.
  out <- nest[keep, , drop = FALSE]
  src <- rep("multi-scale_consensus", nrow(out))
  if ("conflict_status" %in% names(out)) {
    src[out$conflict_status != "none"] <- "multi-scale_consensus:conflict_kept"
  }
  out$candidate_source <- src
  out
}
