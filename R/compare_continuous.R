#' Compare two continuous-assay primary tables
#'
#' Effects are compared on overlapping called+tested intervals when
#' possible. Called-region-only metrics are secondary.
#'
#' @param reference,alternative data.frames with chrom/start/end/effect
#'   and optional FDR/level1/log2FC.
#' @param reference_id,alternative_id Labels.
#' @param strong_abs_effect Predeclared |effect| for strong DERs.
#' @param fdr_threshold FDR cutoff for called DERs.
#' @export
compare_continuous_runs <- function(reference, alternative,
                                    reference_id = "ref",
                                    alternative_id = "alt",
                                    strong_abs_effect = 1,
                                    fdr_threshold = 0.05) {
  .need_pkg("GenomicRanges", "compare_continuous_runs()")
  ref <- reference
  alt <- alternative
  if (!"effect" %in% names(ref) && "log2FC" %in% names(ref)) ref$effect <- ref$log2FC
  if (!"effect" %in% names(alt) && "log2FC" %in% names(alt)) alt$effect <- alt$log2FC
  called <- function(d) {
    if ("level1" %in% names(d)) {
      d[d$level1 %in% c("GAIN", "LOSS"), , drop = FALSE]
    } else if ("FDR" %in% names(d)) {
      d[!is.na(d$FDR) & d$FDR <= fdr_threshold, , drop = FALSE]
    } else {
      d
    }
  }
  to_gr <- function(d) {
    if (!nrow(d)) return(GenomicRanges::GRanges())
    GenomicRanges::GRanges(d$chrom, IRanges::IRanges(as.integer(d$start) + 1L, as.integer(d$end)))
  }
  overlap_stats <- function(a, b) {
    ga <- to_gr(a)
    gb <- to_gr(b)
    ov <- if (length(ga) && length(gb)) GenomicRanges::findOverlaps(ga, gb) else NULL
    if (is.null(ov) || !length(ov)) {
      return(list(pear = NA_real_, spear = NA_real_, signc = NA_real_,
                  signc_hi = NA_real_, med_b = NA_real_, p90_b = NA_real_))
    }
    e1 <- a$effect[S4Vectors::queryHits(ov)]
    e2 <- b$effect[S4Vectors::subjectHits(ov)]
    hi <- is.finite(e1) & abs(e1) >= strong_abs_effect
    dL <- abs(a$start[S4Vectors::queryHits(ov)] - b$start[S4Vectors::subjectHits(ov)])
    dR <- abs(a$end[S4Vectors::queryHits(ov)] - b$end[S4Vectors::subjectHits(ov)])
    list(
      pear = suppressWarnings(stats::cor(e1, e2, use = "complete", method = "pearson")),
      spear = suppressWarnings(stats::cor(e1, e2, use = "complete", method = "spearman")),
      signc = mean(sign(e1) == sign(e2), na.rm = TRUE),
      signc_hi = if (any(hi)) mean(sign(e1[hi]) == sign(e2[hi]), na.rm = TRUE) else NA_real_,
      med_b = stats::median(c(dL, dR), na.rm = TRUE),
      p90_b = as.numeric(stats::quantile(c(dL, dR), 0.9, na.rm = TRUE))
    )
  }
  rc <- called(ref)
  ac <- called(alt)
  all_s <- overlap_stats(ref, alt)
  der_s <- overlap_stats(rc, ac)
  rg <- to_gr(rc)
  ag <- to_gr(ac)
  recov <- if (!length(rg)) NA_real_ else mean(GenomicRanges::countOverlaps(rg, ag) > 0)
  strong <- if (nrow(rc)) rc[is.finite(rc$effect) & abs(rc$effect) >= strong_abs_effect, ] else rc
  sg <- to_gr(strong)
  srec <- if (!length(sg)) NA_real_ else mean(GenomicRanges::countOverlaps(sg, ag) > 0)
  jac <- .bp_jaccard(rg, ag)
  data.frame(
    reference_id = reference_id, alternative_id = alternative_id,
    effect_pearson = all_s$pear, effect_spearman = all_s$spear,
    sign_concordance = all_s$signc, sign_concordance_strong = all_s$signc_hi,
    called_effect_pearson = der_s$pear, called_sign_concordance = der_s$signc,
    DER_bp_jaccard = jac,
    reference_DER_recovery = recov,
    strong_reference_DER_recovery = srec,
    median_boundary_displacement = der_s$med_b,
    p90_boundary_displacement = der_s$p90_b,
    GAIN_reference = sum(rc$level1 == "GAIN" | rc$direction == "GAIN", na.rm = TRUE),
    GAIN_alternative = sum(ac$level1 == "GAIN" | ac$direction == "GAIN", na.rm = TRUE),
    LOSS_reference = sum(rc$level1 == "LOSS" | rc$direction == "LOSS", na.rm = TRUE),
    LOSS_alternative = sum(ac$level1 == "LOSS" | ac$direction == "LOSS", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

.bp_jaccard <- function(a, b) {
  if (!length(a) && !length(b)) return(1)
  if (!length(a) || !length(b)) return(0)
  inter <- sum(as.numeric(GenomicRanges::width(GenomicRanges::intersect(a, b))))
  uni <- sum(as.numeric(GenomicRanges::width(GenomicRanges::reduce(c(a, b)))))
  if (uni <= 0) return(0)
  inter / uni
}
