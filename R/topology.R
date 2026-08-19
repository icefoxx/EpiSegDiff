#' Build condition-specific supported domains
#'
#' Not a differential test. Used only for geometry classification.
#'
#' @param counts Bin × sample count matrix.
#' @param bins Bin table.
#' @param samples Sample metadata.
#' @param condition Condition name to represent.
#' @param method `condition_mean_threshold` or `per_condition_segmentation`.
#' @param mean_quantile Quantile threshold for support.
#' @param size_factors Optional size factors for normalization.
#' @param min_width Minimum domain width in bins.
#' @param penalty Penalty if `per_condition_segmentation`.
#' @return Interval data.frame with `domain_id`, `condition`.
#' @export
build_condition_domains <- function(counts, bins, samples, condition,
                                    method = c("condition_mean_threshold",
                                               "per_condition_segmentation"),
                                    mean_quantile = 0.8,
                                    size_factors = NULL,
                                    min_width = 3L,
                                    penalty = 20) {
  method <- match.arg(method)
  keep <- samples$condition == condition
  if (!any(keep)) {
    stop("No samples for condition '", condition, "'.", call. = FALSE)
  }
  mat <- counts[, keep, drop = FALSE]
  if (is.null(size_factors)) {
    sf <- colSums(mat)
    sf <- sf / exp(mean(log(pmax(sf, 1))))
  } else {
    sf <- size_factors[keep]
  }
  normed <- sweep(mat, 2, sf, "/")
  mu <- rowMeans(normed)
  if (method == "condition_mean_threshold") {
    thr <- as.numeric(stats::quantile(mu[bins$valid %||% TRUE],
                                      probs = mean_quantile, na.rm = TRUE,
                                      names = FALSE))
    supported <- is.finite(mu) & mu >= thr & (bins$valid %||% TRUE)
  } else {
    sc <- log1p(mu)
    sc[!(bins$valid %||% TRUE)] <- NA_real_
    track <- data.frame(
      chrom = bins$chrom, start = bins$start, end = bins$end,
      score = sc, valid = is.finite(sc), stringsAsFactors = FALSE
    )
    seg <- segment_track(track, penalty = penalty, min_width = min_width,
                         loss = "gaussian", engine = "pelt")
    supported <- rep(FALSE, nrow(bins))
    if (nrow(seg)) {
      # keep segments with above-median mean
      keep_seg <- seg$effect_mean >= stats::median(sc, na.rm = TRUE)
      for (i in which(keep_seg)) {
        hit <- bins$chrom == seg$chrom[i] & bins$start >= seg$start[i] &
          bins$end <= seg$end[i]
        supported[hit] <- TRUE
      }
    }
  }
  tmp <- bins
  tmp$valid <- supported
  cvrs <- contiguous_valid_regions(tmp)
  if (nrow(cvrs) == 0L) {
    return(data.frame(
      domain_id = integer(), chrom = character(), start = integer(),
      end = integer(), condition = character(), method = character(),
      stringsAsFactors = FALSE
    ))
  }
  cvrs <- cvrs[cvrs$n_bins >= as.integer(min_width), , drop = FALSE]
  if (nrow(cvrs) == 0L) {
    return(data.frame(
      domain_id = integer(), chrom = character(), start = integer(),
      end = integer(), condition = character(), method = character(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    domain_id = seq_len(nrow(cvrs)),
    chrom = cvrs$chrom,
    start = cvrs$start,
    end = cvrs$end,
    condition = condition,
    method = method,
    stringsAsFactors = FALSE
  )
}

#' Match control and treatment domains
#'
#' @param control_domains,treatment_domains Interval tables.
#' @param min_reciprocal_overlap,max_match_distance,high_reciprocal_overlap
#'   Thresholds from `topology.*`.
#' @return list with `edges` and `components`.
#' @export
match_domains <- function(control_domains, treatment_domains,
                          min_reciprocal_overlap = 0.2,
                          max_match_distance = 5000L,
                          high_reciprocal_overlap = 0.5) {
  C <- control_domains
  T <- treatment_domains
  if (!"domain_id" %in% names(C) && nrow(C)) C$domain_id <- seq_len(nrow(C))
  if (!"domain_id" %in% names(T) && nrow(T)) T$domain_id <- seq_len(nrow(T))
  edges <- .empty_match_edges()
  if (nrow(C) == 0L && nrow(T) == 0L) {
    return(list(edges = edges, components = list()))
  }
  if (nrow(C) && nrow(T)) {
    cand <- list()
    for (i in seq_len(nrow(C))) {
      for (j in seq_len(nrow(T))) {
        if (as.character(C$chrom[i]) != as.character(T$chrom[j])) next
        ro <- .reciprocal_overlap(C$start[i], C$end[i], T$start[j], T$end[j])
        dist <- .gap_distance(C$start[i], C$end[i], T$start[j], T$end[j])
        wr <- (T$end[j] - T$start[j]) / max(C$end[i] - C$start[i], 1)
        if (ro >= min_reciprocal_overlap ||
            (ro == 0 && dist <= max_match_distance)) {
          cand[[length(cand) + 1L]] <- data.frame(
            c_id = C$domain_id[i], t_id = T$domain_id[j],
            chrom = C$chrom[i],
            c_start = C$start[i], c_end = C$end[i],
            t_start = T$start[j], t_end = T$end[j],
            reciprocal_overlap = ro, distance = dist, width_ratio = wr,
            left_shift = T$start[j] - C$start[i],
            right_shift = T$end[j] - C$end[i],
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (length(cand)) {
      edges <- do.call(rbind, cand)
      edges <- edges[order(-edges$reciprocal_overlap, edges$distance,
                           abs(log(pmax(edges$width_ratio, 1e-12))),
                           edges$c_start, edges$t_start), ]
      edges <- .prune_exclusive_edges(edges, high_reciprocal_overlap)
    }
  }
  comps <- .matching_components(C, T, edges)
  list(edges = edges, components = comps)
}

.empty_match_edges <- function() {
  data.frame(
    c_id = integer(), t_id = integer(), chrom = character(),
    c_start = integer(), c_end = integer(),
    t_start = integer(), t_end = integer(),
    reciprocal_overlap = numeric(), distance = integer(),
    width_ratio = numeric(), left_shift = integer(), right_shift = integer(),
    stringsAsFactors = FALSE
  )
}

.prune_exclusive_edges <- function(edges, high_ro) {
  if (nrow(edges) <= 1L) return(edges)
  keep <- rep(TRUE, nrow(edges))
  for (k in seq_len(nrow(edges))) {
    if (!keep[k]) next
    if (edges$reciprocal_overlap[k] < high_ro) next
    # Drop later lower-overlap exclusive competitors on the same C or T
    # only when this edge is a strong 1-1 lock. Split/merge need multi edges,
    # so we only drop a competing edge if it is also high-overlap with a
    # *different* partner (true conflict), not when multiple weak partners exist.
    same_c <- which(keep & edges$c_id == edges$c_id[k] & seq_len(nrow(edges)) != k)
    same_t <- which(keep & edges$t_id == edges$t_id[k] & seq_len(nrow(edges)) != k)
    for (u in c(same_c, same_t)) {
      if (edges$reciprocal_overlap[u] >= high_ro &&
          edges$reciprocal_overlap[u] < edges$reciprocal_overlap[k] - 1e-12) {
        keep[u] <- FALSE
      }
    }
  }
  edges[keep, , drop = FALSE]
}

.matching_components <- function(C, T, edges) {
  c_ids <- if (nrow(C)) C$domain_id else integer()
  t_ids <- if (nrow(T)) T$domain_id else integer()
  used_c <- integer(); used_t <- integer()
  comps <- list()
  # unmatched
  if (nrow(edges) == 0L) {
    for (id in c_ids) {
      comps[[length(comps) + 1L]] <- list(
        type = "control_unmatched", c_ids = id, t_ids = integer(),
        edges = edges[0, ]
      )
    }
    for (id in t_ids) {
      comps[[length(comps) + 1L]] <- list(
        type = "treatment_unmatched", c_ids = integer(), t_ids = id,
        edges = edges[0, ]
      )
    }
    return(comps)
  }
  # connected components on bipartite edges
  n_e <- nrow(edges)
  seen <- rep(FALSE, n_e)
  for (s in seq_len(n_e)) {
    if (seen[s]) next
    stack <- s
    memb <- integer()
    while (length(stack)) {
      k <- stack[1]; stack <- stack[-1]
      if (seen[k]) next
      seen[k] <- TRUE
      memb <- c(memb, k)
      nbr <- which(edges$c_id == edges$c_id[k] | edges$t_id == edges$t_id[k])
      stack <- c(stack, nbr[!seen[nbr]])
    }
    ed <- edges[sort(unique(memb)), , drop = FALSE]
    cid <- unique(ed$c_id); tid <- unique(ed$t_id)
    used_c <- c(used_c, cid); used_t <- c(used_t, tid)
    type <- if (length(cid) == 1L && length(tid) == 1L) "one_to_one"
    else if (length(cid) == 1L && length(tid) >= 2L) "one_to_many"
    else if (length(cid) >= 2L && length(tid) == 1L) "many_to_one"
    else "many_to_many"
    comps[[length(comps) + 1L]] <- list(
      type = type, c_ids = cid, t_ids = tid, edges = ed
    )
  }
  for (id in setdiff(c_ids, used_c)) {
    comps[[length(comps) + 1L]] <- list(
      type = "control_unmatched", c_ids = id, t_ids = integer(),
      edges = edges[0, ]
    )
  }
  for (id in setdiff(t_ids, used_t)) {
    comps[[length(comps) + 1L]] <- list(
      type = "treatment_unmatched", c_ids = integer(), t_ids = id,
      edges = edges[0, ]
    )
  }
  comps
}

#' Classify a single matching component
#'
#' @param component One element of [match_domains()]$components.
#' @param control_domains,treatment_domains Domain tables.
#' @param differential Optional tested differential intervals with
#'   `direction` in `GAIN`/`LOSS`.
#' @param thresholds Named list of topology thresholds.
#' @return A one-row event data.frame.
#' @export
classify_boundary_event <- function(component,
                                    control_domains,
                                    treatment_domains,
                                    differential = NULL,
                                    thresholds = NULL) {
  th <- .topology_thresholds(thresholds)
  C <- control_domains
  T <- treatment_domains
  type <- component$type
  ed <- component$edges

  if (type == "many_to_many") {
    return(.event_row("COMPLEX", NA, NA, NA, NA, NA, NA, NA, 0, 0, NA, NA,
                      "RULE0_MANY_TO_MANY", th))
  }

  if (type == "one_to_many") {
    ev <- .try_split(component, C, T, th)
    if (!is.null(ev)) return(ev)
    return(.event_row("COMPLEX", NA, NA, NA, NA, NA, NA, NA, 0, 0, NA, NA,
                      "RULE1_SPLIT_FAILED", th))
  }
  if (type == "many_to_one") {
    ev <- .try_merge(component, C, T, th)
    if (!is.null(ev)) return(ev)
    return(.event_row("COMPLEX", NA, NA, NA, NA, NA, NA, NA, 0, 0, NA, NA,
                      "RULE2_MERGE_FAILED", th))
  }

  if (type == "treatment_unmatched") {
    trow <- T[T$domain_id == component$t_ids[1], ]
    dir <- .overlap_direction(trow, differential)
    cls <- if (identical(dir, "LOSS")) "COMPLEX" else "GAIN"
    return(.event_row(cls, trow$chrom, trow$start, trow$end,
                      NA, NA, trow$start, trow$end, 0, 0, NA, NA,
                      "RULE8_UNMATCHED_TREATMENT", th))
  }
  if (type == "control_unmatched") {
    crow <- C[C$domain_id == component$c_ids[1], ]
    dir <- .overlap_direction(crow, differential)
    cls <- if (identical(dir, "GAIN")) "COMPLEX" else "LOSS"
    return(.event_row(cls, crow$chrom, crow$start, crow$end,
                      crow$start, crow$end, NA, NA, 0, 0, NA, NA,
                      "RULE8_UNMATCHED_CONTROL", th))
  }

  # one_to_one: Level II geometry is a partition of (Δ_L, Δ_R, R_W).
  # BIDIRECTIONAL_SHIFT is not a primary class (canonical opposite-end
  # motion is EXPANSION or CONTRACTION).
  e1 <- ed[1, ]
  ev <- .classify_one_to_one(e1, differential, th)
  if (!is.null(ev)) return(ev)
  .event_row("COMPLEX", e1$chrom, min(e1$c_start, e1$t_start),
             max(e1$c_end, e1$t_end), e1$c_start, e1$c_end, e1$t_start, e1$t_end,
             e1$left_shift, e1$right_shift, e1$width_ratio, e1$reciprocal_overlap,
             "RULE9_COMPLEX", th)
}

#' Classify all matched components
#'
#' @inheritParams classify_boundary_event
#' @param matching Output of [match_domains()].
#' @export
classify_boundary_events <- function(matching, control_domains, treatment_domains,
                                     differential = NULL, thresholds = NULL) {
  if (!length(matching$components)) {
    return(.empty_events())
  }
  rows <- lapply(seq_along(matching$components), function(i) {
    ev <- classify_boundary_event(
      matching$components[[i]], control_domains, treatment_domains,
      differential = differential, thresholds = thresholds
    )
    ev$event_id <- i
    ev
  })
  do.call(rbind, rows)
}

.topology_thresholds <- function(thresholds) {
  d <- list(
    min_reciprocal_overlap = 0.2,
    high_reciprocal_overlap = 0.5,
    max_match_distance = 5000,
    min_width_ratio_change = 1.25,
    min_boundary_shift = 100,
    min_internal_gap_width = 200,
    similar_width_ratio = 1.25,
    max_overlap_for_shift = 0.85
  )
  if (is.null(thresholds)) return(d)
  utils::modifyList(d, thresholds)
}

.empty_events <- function() {
  data.frame(
    event_id = integer(), chrom = character(), start = integer(), end = integer(),
    primary_event = character(), secondary_event = character(),
    control_start = integer(), control_end = integer(),
    treatment_start = integer(), treatment_end = integer(),
    boundary_shift_left = numeric(), boundary_shift_right = numeric(),
    width_ratio = numeric(), matching_score = numeric(),
    classification_rule = character(),
    stringsAsFactors = FALSE
  )
}

.event_row <- function(primary, chrom, start, end, cs, ce, ts, te,
                       ls, rs, wr, ms, rule, th, secondary = NA_character_) {
  data.frame(
    event_id = NA_integer_,
    chrom = as.character(chrom %||% NA),
    start = as.integer(start %||% NA),
    end = as.integer(end %||% NA),
    primary_event = primary,
    secondary_event = secondary,
    control_start = as.integer(cs %||% NA),
    control_end = as.integer(ce %||% NA),
    treatment_start = as.integer(ts %||% NA),
    treatment_end = as.integer(te %||% NA),
    boundary_shift_left = as.numeric(ls %||% NA),
    boundary_shift_right = as.numeric(rs %||% NA),
    width_ratio = as.numeric(wr %||% NA),
    matching_score = as.numeric(ms %||% NA),
    classification_rule = sprintf(
      "%s;min_ro=%s;high_ro=%s;min_wr=%s;min_shift=%s;min_gap=%s;sim_wr=%s",
      rule, th$min_reciprocal_overlap, th$high_reciprocal_overlap,
      th$min_width_ratio_change, th$min_boundary_shift,
      th$min_internal_gap_width, th$similar_width_ratio
    ),
    stringsAsFactors = FALSE
  )
}

.try_split <- function(component, C, T, th) {
  trows <- T[T$domain_id %in% component$t_ids, ]
  trows <- trows[order(trows$start), ]
  if (nrow(trows) < 2L) return(NULL)
  crow <- C[C$domain_id == component$c_ids[1], ]
  gaps <- trows$start[-1] - trows$end[-nrow(trows)]
  if (any(gaps < th$min_internal_gap_width)) return(NULL)
  # each T should overlap C
  if (nrow(crow) != 1L) return(NULL)
  ro <- .reciprocal_overlap(crow$start, crow$end, trows$start, trows$end)
  cont <- .interval_intersection(trows$start, trows$end, crow$start, crow$end) /
    pmax(trows$end - trows$start, 1)
  ok <- pmax(ro, cont)
  if (!length(ok) || any(!is.finite(ok)) || !all(ok >= th$min_reciprocal_overlap)) {
    return(NULL)
  }
  tw <- sum(trows$end - trows$start)
  cw <- crow$end - crow$start
  if (tw > cw * th$similar_width_ratio) return(NULL)
  .event_row("SPLIT", crow$chrom, crow$start, crow$end,
             crow$start, crow$end, min(trows$start), max(trows$end),
             0, 0, tw / max(cw, 1), mean(ro), "RULE1_SPLIT", th)
}

.try_merge <- function(component, C, T, th) {
  crows <- C[C$domain_id %in% component$c_ids, ]
  crows <- crows[order(crows$start), ]
  if (nrow(crows) < 2L) return(NULL)
  trow <- T[T$domain_id == component$t_ids[1], ]
  gaps <- crows$start[-1] - crows$end[-nrow(crows)]
  if (any(gaps < th$min_internal_gap_width)) return(NULL)
  if (nrow(trow) != 1L) return(NULL)
  ro <- .reciprocal_overlap(trow$start, trow$end, crows$start, crows$end)
  cont <- .interval_intersection(crows$start, crows$end, trow$start, trow$end) /
    pmax(crows$end - crows$start, 1)
  ok <- pmax(ro, cont)
  if (!length(ok) || any(!is.finite(ok)) || !all(ok >= th$min_reciprocal_overlap)) {
    return(NULL)
  }
  .event_row("MERGE", trow$chrom, trow$start, trow$end,
             min(crows$start), max(crows$end), trow$start, trow$end,
             0, 0, (trow$end - trow$start) / max(sum(crows$end - crows$start), 1),
             mean(ro), "RULE2_MERGE", th)
}

.classify_one_to_one <- function(e1, differential, th) {
  dL <- as.numeric(e1$left_shift)
  dR <- as.numeric(e1$right_shift)
  rw <- as.numeric(e1$width_ratio)
  if (!is.finite(rw) || rw <= 0) rw <- NA_real_
  delta <- as.numeric(th$min_boundary_shift)
  rho <- as.numeric(th$min_width_ratio_change)
  rsim <- as.numeric(th$similar_width_ratio)

  left_out <- is.finite(dL) && dL < -delta
  left_in <- is.finite(dL) && dL > delta
  left_st <- is.finite(dL) && abs(dL) <= delta
  right_out <- is.finite(dR) && dR > delta
  right_in <- is.finite(dR) && dR < -delta
  right_st <- is.finite(dR) && abs(dR) <= delta
  wider <- is.finite(rw) && rw >= rho
  narrower <- is.finite(rw) && rw <= (1 / rho)
  similar <- is.finite(rw) && rw >= (1 / rsim) && rw <= rsim

  geo_row <- function(cls, rule, sec = NA_character_) {
    .event_row(cls, e1$chrom, min(e1$c_start, e1$t_start),
               max(e1$c_end, e1$t_end), e1$c_start, e1$c_end,
               e1$t_start, e1$t_end, dL, dR, rw, e1$reciprocal_overlap,
               rule, th, sec)
  }

  # Internal / intensity only when outer geometry is stable (non-overlapping
  # with expansion/contraction/shift, which require |Δ| > δ on at least one end).
  if (left_st && right_st && similar) {
    ev <- .try_internal(e1, differential, th)
    if (!is.null(ev)) return(ev)
    ev <- .try_intensity(e1, differential, th)
    if (!is.null(ev)) return(ev)
    return(geo_row("COMPLEX", "RULE_STABLE_UNCLASSIFIED"))
  }

  # 3×3 partition of (left state × right state). Opposite-end motion is
  # never labeled BIDIRECTIONAL_SHIFT.
  if (left_out && right_out) {
    if (wider) return(geo_row("EXPANSION", "RULE3_EXPANSION_BOTH", "BOTH_ENDS"))
    return(geo_row("COMPLEX", "RULE_BOTH_OUT_NOT_WIDER"))
  }
  if (left_in && right_in) {
    if (narrower) return(geo_row("CONTRACTION", "RULE4_CONTRACTION_BOTH", "BOTH_ENDS"))
    return(geo_row("COMPLEX", "RULE_BOTH_IN_NOT_NARROWER"))
  }
  if (left_out && right_in) {
    if (similar) return(geo_row("LEFT_SHIFT", "RULE5_LEFT_SHIFT"))
    return(geo_row("COMPLEX", "RULE_LEFT_MOVE_WIDTH_NOT_SIMILAR"))
  }
  if (left_in && right_out) {
    if (similar) return(geo_row("RIGHT_SHIFT", "RULE5_RIGHT_SHIFT"))
    return(geo_row("COMPLEX", "RULE_RIGHT_MOVE_WIDTH_NOT_SIMILAR"))
  }
  if (left_out && right_st) {
    if (wider) return(geo_row("EXPANSION", "RULE3_EXPANSION_LEFT", "LEFT_END"))
    return(geo_row("COMPLEX", "RULE_LEFT_OUT_NOT_WIDER"))
  }
  if (left_st && right_out) {
    if (wider) return(geo_row("EXPANSION", "RULE3_EXPANSION_RIGHT", "RIGHT_END"))
    return(geo_row("COMPLEX", "RULE_RIGHT_OUT_NOT_WIDER"))
  }
  if (left_in && right_st) {
    if (narrower) return(geo_row("CONTRACTION", "RULE4_CONTRACTION_LEFT", "LEFT_END"))
    return(geo_row("COMPLEX", "RULE_LEFT_IN_NOT_NARROWER"))
  }
  if (left_st && right_in) {
    if (narrower) return(geo_row("CONTRACTION", "RULE4_CONTRACTION_RIGHT", "RIGHT_END"))
    return(geo_row("COMPLEX", "RULE_RIGHT_IN_NOT_NARROWER"))
  }
  geo_row("COMPLEX", "RULE9_COMPLEX")
}

.try_internal <- function(e1, differential, th) {
  r <- th$similar_width_ratio
  if (e1$reciprocal_overlap < th$high_reciprocal_overlap) return(NULL)
  if (!(e1$width_ratio >= 1 / r && e1$width_ratio <= r)) return(NULL)
  if (is.null(differential) || nrow(differential) == 0L) return(NULL)
  inter_s <- max(e1$c_start, e1$t_start)
  inter_e <- min(e1$c_end, e1$t_end)
  d <- differential[as.character(differential$chrom) == as.character(e1$chrom), ]
  if (!nrow(d) || !"direction" %in% names(d)) return(NULL)
  for (i in seq_len(nrow(d))) {
    cont <- .interval_intersection(d$start[i], d$end[i], inter_s, inter_e) /
      max(d$end[i] - d$start[i], 1)
    covers_left <- abs(d$start[i] - min(e1$c_start, e1$t_start)) <= th$min_boundary_shift
    covers_right <- abs(d$end[i] - max(e1$c_end, e1$t_end)) <= th$min_boundary_shift
    if (cont >= th$high_reciprocal_overlap && !covers_left && !covers_right) {
      if (d$direction[i] == "GAIN") {
        return(.event_row("INTERNAL_GAIN", e1$chrom, d$start[i], d$end[i],
                          e1$c_start, e1$c_end, e1$t_start, e1$t_end,
                          e1$left_shift, e1$right_shift, e1$width_ratio,
                          e1$reciprocal_overlap, "RULE7_INTERNAL_GAIN", th))
      }
      if (d$direction[i] == "LOSS") {
        return(.event_row("INTERNAL_LOSS", e1$chrom, d$start[i], d$end[i],
                          e1$c_start, e1$c_end, e1$t_start, e1$t_end,
                          e1$left_shift, e1$right_shift, e1$width_ratio,
                          e1$reciprocal_overlap, "RULE7_INTERNAL_LOSS", th))
      }
    }
  }
  NULL
}

.try_intensity <- function(e1, differential, th) {
  r <- th$similar_width_ratio
  stable <- e1$reciprocal_overlap >= th$high_reciprocal_overlap &&
    e1$width_ratio >= 1 / r && e1$width_ratio <= r &&
    abs(e1$left_shift) < th$min_boundary_shift &&
    abs(e1$right_shift) < th$min_boundary_shift
  if (!stable) return(NULL)
  dir <- "GAIN"
  if (!is.null(differential) && nrow(differential) && "direction" %in% names(differential)) {
    d <- differential[as.character(differential$chrom) == as.character(e1$chrom), ]
    ov <- .interval_intersection(e1$c_start, e1$c_end, d$start, d$end) > 0
    if (any(ov)) dir <- d$direction[which(ov)[1]]
  }
  if (!dir %in% c("GAIN", "LOSS")) dir <- "GAIN"
  .event_row(dir, e1$chrom, min(e1$c_start, e1$t_start), max(e1$c_end, e1$t_end),
             e1$c_start, e1$c_end, e1$t_start, e1$t_end,
             e1$left_shift, e1$right_shift, e1$width_ratio, e1$reciprocal_overlap,
             if (dir == "GAIN") "RULE8_GAIN" else "RULE8_LOSS", th)
}

.overlap_direction <- function(interval, differential) {
  if (is.null(differential) || !nrow(differential) || !"direction" %in% names(differential)) {
    return(NA_character_)
  }
  d <- differential[as.character(differential$chrom) == as.character(interval$chrom), ]
  if (!nrow(d)) return(NA_character_)
  ov <- .interval_intersection(interval$start, interval$end, d$start, d$end) > 0
  if (!any(ov)) return(NA_character_)
  d$direction[which(ov)[1]]
}
