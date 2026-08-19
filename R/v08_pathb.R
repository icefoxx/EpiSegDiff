#' Path B hierarchical / multiresolution prototypes (v0.8)
#'
#' Candidates A–D are research prototypes. Flat engines remain controls.
#' @name v08_pathb

.empty_pathb_regions <- function() {
  data.frame(
    chrom = character(), start = integer(), end = integer(),
    effect = numeric(), direction = character(), level1 = character(),
    component_type = character(), region_id = integer(),
    parent_id = integer(), level = integer(),
    scale_bp = integer(), support = numeric(),
    stringsAsFactors = FALSE
  )
}

.pathb_mad_thr <- function(x) {
  thr <- 0.15 + 0.35 * stats::mad(as.numeric(x), constant = 1)
  if (!is.finite(thr) || thr <= 0) 0.2 else thr
}

.pathb_seg_to_regions <- function(seg, bed, score, component_type, level,
                                  parent_id = NA_integer_, thr = NULL) {
  if (is.null(thr)) thr <- .pathb_mad_thr(score)
  if (!nrow(seg)) return(.empty_pathb_regions())
  i1 <- seg$start_idx
  i2 <- seg$end_idx
  eff <- seg$effect_mean
  keep <- is.finite(eff) & abs(eff) > thr
  if (!any(keep)) return(.empty_pathb_regions())
  i1 <- i1[keep]; i2 <- i2[keep]; eff <- eff[keep]
  data.frame(
    chrom = bed$chrom[i1],
    start = bed$start[i1],
    end = bed$end[i2],
    effect = eff,
    direction = ifelse(eff > 0, "GAIN", "LOSS"),
    level1 = ifelse(eff > 0, "GAIN", "LOSS"),
    component_type = component_type,
    region_id = NA_integer_,
    parent_id = parent_id,
    level = as.integer(level),
    scale_bp = as.integer(bed$end[i2] - bed$start[i1]),
    support = abs(eff),
    stringsAsFactors = FALSE
  )
}

#' Soft-threshold a 1-D series
.soft_thresh <- function(x, tau) {
  sign(x) * pmax(abs(x) - tau, 0)
}

#' Candidate A: sequential TV broad + residual local
#' @export
pathb_additive <- function(score, bed, lambda_B_mult = 8, w_B = 10L,
                           penalty = 8) {
  score <- as.numeric(score)
  n <- length(score)
  md <- stats::mad(diff(score), constant = 1)
  if (!is.finite(md) || md <= 0) md <- 1e-3
  B <- .tv_denoise(score, lambda = lambda_B_mult * md, niter = 40L)
  L <- score - B
  segB <- segment_signal(B, positions = bed$start, engine = "pelt",
                         penalty = penalty, min_width = w_B)
  segL <- segment_signal(L, positions = bed$start, engine = "fused_lasso",
                         penalty = penalty, min_width = 1L)
  thrB <- .pathb_mad_thr(B)
  thrL <- max(0.1, 0.5 * .pathb_mad_thr(L))
  broad <- .pathb_seg_to_regions(segB, bed, B, "broad", 1L, thr = thrB)
  local <- .pathb_seg_to_regions(segL, bed, L, "local", 2L, thr = thrL)
  out <- rbind(broad, local)
  if (nrow(out)) out$region_id <- seq_len(nrow(out))
  list(regions = out, B = B, L = L, architecture = "A_additive",
       params = list(lambda_B_mult = lambda_B_mult, w_B = w_B))
}

#' Candidate B: hierarchical PELT parents + residual children
#' @export
pathb_hierarchy <- function(score, bed, w_B = 10L, penalty = 8,
                            allow_orphan_local = TRUE) {
  score <- as.numeric(score)
  segP <- segment_signal(score, positions = bed$start, engine = "pelt",
                         penalty = penalty, min_width = w_B)
  thr <- .pathb_mad_thr(score)
  parents <- .pathb_seg_to_regions(segP, bed, score, "broad", 1L, thr = thr)
  rows <- list()
  rid <- 0L
  if (nrow(parents)) {
    for (i in seq_len(nrow(parents))) {
      rid <- rid + 1L
      pr <- parents[i, , drop = FALSE]
      pr$region_id <- rid
      pr$parent_id <- NA_integer_
      rows[[length(rows) + 1L]] <- pr
      ii <- which(bed$start >= pr$start & bed$end <= pr$end)
      if (length(ii) < 4L) next
      resid <- score[ii] - mean(score[ii])
      bed_i <- bed[ii, , drop = FALSE]
      segC <- segment_signal(resid, positions = bed_i$start, engine = "fused_lasso",
                             penalty = max(2, penalty / 2), min_width = 1L)
      ch <- .pathb_seg_to_regions(segC, bed_i, resid, "local", 2L,
                                  parent_id = rid,
                                  thr = max(0.1, 0.5 * .pathb_mad_thr(resid)))
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
  # orphan locals only when no broad parent (isolated narrow)
  if (allow_orphan_local && !nrow(parents)) {
    segL <- segment_signal(score, positions = bed$start, engine = "fused_lasso",
                           penalty = penalty, min_width = 1L)
    loc <- .pathb_seg_to_regions(segL, bed, score, "local", 1L, thr = thr)
    if (nrow(loc)) {
      for (j in seq_len(nrow(loc))) {
        rid <- rid + 1L
        loc$region_id[j] <- rid
        loc$parent_id[j] <- NA_integer_
        loc$component_type[j] <- "local"
        loc$level[j] <- 1L
        rows[[length(rows) + 1L]] <- loc[j, , drop = FALSE]
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else .empty_pathb_regions()
  list(regions = out, B = score, L = NULL, architecture = "B_hierarchy",
       params = list(w_B = w_B, allow_orphan_local = allow_orphan_local))
}

#' Candidate C: coarse / medium / fine multiresolution
#' @export
pathb_multires <- function(score, bed, k_c = 5L, k_m = 2L, penalty = 8) {
  score <- as.numeric(score)
  n <- length(score)
  # coarse
  g <- ceiling(seq_len(n) / k_c)
  cs <- as.numeric(tapply(score, g, mean))
  cpos <- as.integer(tapply(bed$start, g, min))
  cbed <- data.frame(chrom = bed$chrom[match(cpos, bed$start)],
                     start = cpos,
                     end = as.integer(tapply(bed$end, g, max)))
  segC <- segment_signal(cs, positions = cpos, engine = "pelt",
                         penalty = penalty, min_width = max(2L, as.integer(10 / k_c)))
  broad <- .pathb_seg_to_regions(segC, cbed, cs, "broad", 1L, thr = .pathb_mad_thr(cs))
  # map coarse fit back
  B <- cs[g]
  Lmed <- score - B
  gm <- ceiling(seq_len(n) / k_m)
  ms <- as.numeric(tapply(Lmed, gm, mean))
  mpos <- as.integer(tapply(bed$start, gm, min))
  mbed <- data.frame(chrom = bed$chrom[match(mpos, bed$start)],
                     start = mpos,
                     end = as.integer(tapply(bed$end, gm, max)))
  segM <- segment_signal(ms, positions = mpos, engine = "pelt",
                         penalty = penalty, min_width = 2L)
  med <- .pathb_seg_to_regions(segM, mbed, ms, "intermediate", 2L,
                               thr = .pathb_mad_thr(ms))
  Lf <- Lmed - ms[gm]
  segF <- segment_signal(Lf, positions = bed$start, engine = "fused_lasso",
                         penalty = penalty, min_width = 1L)
  loc <- .pathb_seg_to_regions(segF, bed, Lf, "local", 3L, thr = .pathb_mad_thr(Lf))
  out <- rbind(broad, med, loc)
  if (nrow(out)) out$region_id <- seq_len(nrow(out))
  list(regions = out, B = B, L = Lf, architecture = "C_multires",
       params = list(k_c = k_c, k_m = k_m))
}

#' Candidate D: joint alternating dual-TV
#' @export
pathb_joint <- function(score, bed, lambda_B_mult = 8, lambda_L_mult = 1,
                        w_B = 10L, niter = 8L, penalty = 8) {
  score <- as.numeric(score)
  md <- stats::mad(diff(score), constant = 1)
  if (!is.finite(md) || md <= 0) md <- 1e-3
  lamB <- lambda_B_mult * md
  lamL <- lambda_L_mult * md
  L <- numeric(length(score))
  B <- score
  for (it in seq_len(niter)) {
    B <- .tv_denoise(score - L, lambda = lamB, niter = 20L)
    L <- .tv_denoise(score - B, lambda = lamL, niter = 20L)
  }
  segB <- segment_signal(B, positions = bed$start, engine = "pelt",
                         penalty = penalty, min_width = w_B)
  segL <- segment_signal(L, positions = bed$start, engine = "fused_lasso",
                         penalty = penalty, min_width = 1L)
  thrB <- .pathb_mad_thr(B)
  thrL <- max(0.1, 0.5 * .pathb_mad_thr(L))
  broad <- .pathb_seg_to_regions(segB, bed, B, "broad", 1L, thr = thrB)
  local <- .pathb_seg_to_regions(segL, bed, L, "local", 2L, thr = thrL)
  out <- rbind(broad, local)
  if (nrow(out)) out$region_id <- seq_len(nrow(out))
  list(regions = out, B = B, L = L, architecture = "D_joint",
       params = list(lambda_B_mult = lambda_B_mult, lambda_L_mult = lambda_L_mult,
                     w_B = w_B, niter = niter))
}

#' Flat control wrappers returning Path B-like region tables
#' @export
pathb_control <- function(score, bed, which = c("pelt", "fused", "std_pelt"),
                          penalty = 8) {
  which <- match.arg(which)
  if (which == "pelt") {
    seg <- segment_signal(score, positions = bed$start, engine = "pelt",
                          penalty = penalty, min_width = 3L)
    reg <- .pathb_seg_to_regions(seg, bed, score, "flat", 1L)
  } else if (which == "fused") {
    seg <- segment_signal(score, positions = bed$start, engine = "fused_lasso",
                          penalty = penalty, min_width = 1L)
    reg <- .pathb_seg_to_regions(seg, bed, score, "flat", 1L)
  } else {
    # caller must pass standardized score
    seg <- segment_signal(score, positions = bed$start, engine = "pelt",
                          penalty = penalty, min_width = 1L)
    reg <- .pathb_seg_to_regions(seg, bed, score, "flat", 1L)
  }
  if (nrow(reg)) reg$region_id <- seq_len(nrow(reg))
  list(regions = reg, B = score, L = NULL, architecture = paste0("CTL_", which),
       params = list())
}

#' Run a named Path B / control architecture on a simulation
#' @export
run_pathb <- function(sim, architecture = "A_additive", score_engine = "raw_effect",
                      ...) {
  sc <- compute_diff_score(sim$signal, sim$samples, "JNKi", "DMSO", score_engine)
  bed <- sim$bed
  fit <- switch(
    architecture,
    A_additive = pathb_additive(sc, bed, ...),
    B_hierarchy = pathb_hierarchy(sc, bed, ...),
    C_multires = pathb_multires(sc, bed, ...),
    D_joint = pathb_joint(sc, bed, ...),
    CTL_pelt = pathb_control(sc, bed, "pelt"),
    CTL_fused = pathb_control(sc, bed, "fused"),
    CTL_std = {
      sc2 <- compute_diff_score(sim$signal, sim$samples, "JNKi", "DMSO",
                                "standardized")
      pathb_control(sc2, bed, "std_pelt")
    },
    CTL_union = {
      a <- pathb_control(sc, bed, "fused")$regions
      b <- pathb_control(sc, bed, "pelt")$regions
      if (nrow(a) && !"FDR" %in% names(a)) a$FDR <- 0.01
      if (nrow(b) && !"FDR" %in% names(b)) b$FDR <- 0.01
      u <- tryCatch(union_resolve_calls(list(a = a, b = b)),
                    error = function(e) rbind(a, b))
      if (nrow(u)) {
        u$component_type <- "flat"
        u$region_id <- seq_len(nrow(u))
        u$parent_id <- NA_integer_
        u$level <- 1L
        if (!"scale_bp" %in% names(u))
          u$scale_bp <- as.integer(u$end - u$start)
        if (!"support" %in% names(u)) u$support <- abs(u$effect)
      }
      list(regions = if (nrow(u)) u else .empty_pathb_regions(),
           B = sc, L = NULL, architecture = "CTL_union", params = list())
    },
    stop("Unknown architecture: ", architecture, call. = FALSE)
  )
  fit$score_engine <- score_engine
  fit
}

#' Nested / hierarchy metrics vs two-event truth (broad then local)
#' @export
score_nested_hierarchy <- function(regions, truth, bin_size) {
  if (nrow(truth) < 2L) {
    return(data.frame(
      p_parent = NA, p_child = NA, p_both = NA, p_merged = NA,
      p_child_only = NA, p_parent_only = NA, false_child = NA,
      n_broad_call = NA, n_local_call = NA, n_levels = NA,
      stringsAsFactors = FALSE
    ))
  }
  # truth row1 = broad, row2 = local by convention
  tgr <- GenomicRanges::GRanges(
    "1", IRanges::IRanges(as.integer((truth$start_bin - 1L) * bin_size + 1L),
                          as.integer(truth$end_bin * bin_size))
  )
  if (!nrow(regions)) {
    return(data.frame(
      p_parent = 0, p_child = 0, p_both = 0, p_merged = 0,
      p_child_only = 0, p_parent_only = 0, false_child = 0,
      n_broad_call = 0, n_local_call = 0, n_levels = 0,
      stringsAsFactors = FALSE
    ))
  }
  broad <- regions[regions$component_type %in% c("broad", "flat"), , drop = FALSE]
  local <- regions[regions$component_type %in% c("local", "intermediate"), ,
                   drop = FALSE]
  # flat controls: use all regions for both roles
  if (!nrow(broad)) broad <- regions
  if (!nrow(local) && all(regions$component_type == "flat")) local <- regions
  to_gr <- function(d) {
    if (!nrow(d)) return(GenomicRanges::GRanges())
    GenomicRanges::GRanges(d$chrom, IRanges::IRanges(d$start + 1L, d$end))
  }
  bg <- to_gr(broad)
  lg <- to_gr(local)
  ag <- to_gr(regions)
  has_scale <- any(regions$component_type %in% c("broad", "local", "intermediate"))
  parent_hit <- length(bg) > 0 && GenomicRanges::countOverlaps(tgr[1], bg) > 0
  child_hit <- length(lg) > 0 && GenomicRanges::countOverlaps(tgr[2], lg) > 0
  # Flat controls cannot assign scale; a call overlapping both truths is a merge.
  flat_merge <- FALSE
  if (!has_scale && length(ag)) {
    ov1 <- GenomicRanges::countOverlaps(ag, tgr[1]) > 0
    ov2 <- GenomicRanges::countOverlaps(ag, tgr[2]) > 0
    flat_merge <- any(ov1 & ov2)
    parent_hit <- any(ov1)
    child_hit <- any(ov2)
  }
  # Path B: parent may legally overlap embedded local truth; absorption =
  # parent without a separate local/child component.
  absorbed <- has_scale && parent_hit && !child_hit
  merged <- if (has_scale) absorbed else flat_merge
  # Separate recovery only counts when scale-typed components both hit.
  both_ok <- has_scale && parent_hit && child_hit
  false_child <- if (has_scale && length(lg)) {
    mean(GenomicRanges::countOverlaps(lg, tgr[2]) == 0)
  } else {
    0
  }
  data.frame(
    p_parent = as.numeric(parent_hit),
    p_child = as.numeric(child_hit),
    p_both = as.numeric(both_ok),
    p_merged = as.numeric(merged),
    p_child_only = as.numeric(child_hit && !parent_hit),
    p_parent_only = as.numeric(parent_hit && !child_hit),
    false_child = false_child,
    n_broad_call = nrow(broad),
    n_local_call = nrow(local),
    n_levels = if (nrow(regions)) max(regions$level, na.rm = TRUE) else 0,
    stringsAsFactors = FALSE
  )
}

#' Component leakage: truth broad-only or local-only
#' @export
score_component_leakage <- function(fit, mode = c("broad_only", "local_only")) {
  mode <- match.arg(mode)
  reg <- fit$regions
  n_b <- sum(reg$component_type == "broad")
  n_l <- sum(reg$component_type %in% c("local", "intermediate"))
  n_f <- sum(reg$component_type == "flat")
  if (mode == "broad_only") {
    # false local rate
    data.frame(mode = mode, n_broad = n_b, n_local = n_l, n_flat = n_f,
               false_local = n_l > 0, n_false_local = n_l,
               false_broad = FALSE, n_false_broad = 0L,
               stringsAsFactors = FALSE)
  } else {
    data.frame(mode = mode, n_broad = n_b, n_local = n_l, n_flat = n_f,
               false_local = FALSE, n_false_local = 0L,
               false_broad = n_b > 0, n_false_broad = n_b,
               stringsAsFactors = FALSE)
  }
}
