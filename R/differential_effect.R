#' Estimate a per-bin differential-effect track
#'
#' Fits a replicate-aware edgeR model and extracts logFC, SE, and a
#' configurable segmentation score. Dispersion is estimated on bins
#' before any segmentation.
#'
#' @param counts Integer matrix (bins × samples).
#' @param samples Sample metadata (data.frame).
#' @param offsets Numeric sample offsets from [normalize_counts()].
#' @param score `logFC`, `standardized_effect`, or `shrunken_effect`.
#' @param control_level Reference condition.
#' @param treatment_level Treatment condition.
#' @param prior_count edgeR prior count used only where edgeR requires it.
#' @param valid Optional logical vector; invalid bins get NA scores.
#' @return A list with `effect` data.frame-ready vectors and fit diagnostics.
#' @export
estimate_differential_effect <- function(counts, samples, offsets,
                                         score = c("logFC", "standardized_effect",
                                                   "shrunken_effect"),
                                         control_level = "control",
                                         treatment_level = "treatment",
                                         prior_count = 0.125,
                                         valid = NULL) {
  score <- match.arg(score)
  .need_pkg("edgeR", "estimate_differential_effect()")
  if (ncol(counts) != nrow(samples)) {
    stop("counts columns and samples rows differ.", call. = FALSE)
  }
  if (!is.null(colnames(counts)) && "sample_id" %in% names(samples)) {
    if (!identical(as.character(colnames(counts)), as.character(samples$sample_id))) {
      # Reorder counts to metadata if names match as a set.
      if (setequal(colnames(counts), samples$sample_id)) {
        counts <- counts[, samples$sample_id, drop = FALSE]
      } else {
        stop("count matrix column names do not match sample_id.", call. = FALSE)
      }
    }
  }

  cond <- .condition_factor(samples, control_level, treatment_level)
  design <- stats::model.matrix(~ cond)
  colnames(design) <- c("Intercept", "treatment")

  y <- edgeR::DGEList(counts = counts, samples = data.frame(condition = cond))
  # Attach offsets: edgeR offset is added on the natural log mean scale.
  off <- matrix(as.numeric(offsets), nrow = nrow(counts), ncol = ncol(counts),
                byrow = TRUE)
  y <- edgeR::scaleOffset(y, offset = off)

  # estimateDisp fails on all-zero matrices; handle explicitly.
  if (all(counts == 0L)) {
    n <- nrow(counts)
    return(.empty_effect(n, score, valid))
  }

  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)

  coef_idx <- 2L
  # edgeR 4.x QL objects do not always expose stdev.unscaled. Recover SE
  # from the QL F statistic: t = sign(logFC)*sqrt(F), SE = logFC / t.
  tt <- tryCatch(edgeR::glmQLFTest(fit, coef = coef_idx), error = function(e) NULL)
  if (is.null(tt) || is.null(tt$table)) {
    stop("edgeR QL test failed; cannot form per-bin effects.", call. = FALSE)
  }
  beta <- as.numeric(tt$table$logFC)
  fstat <- as.numeric(tt$table$F)
  tstat <- sign(beta) * sqrt(pmax(fstat, 0))
  se <- ifelse(is.finite(tstat) & abs(tstat) > 1e-8, beta / tstat, NA_real_)
  if (!is.null(fit$coefficients) && ncol(as.matrix(fit$coefficients)) >= coef_idx) {
    beta_fit <- as.numeric(fit$coefficients[, coef_idx])
    if (length(beta_fit) == length(beta)) {
      # Prefer the fitted coefficient when the test table was shrunk/aliased.
      beta <- beta_fit
      se <- ifelse(is.finite(tstat) & abs(tstat) > 1e-8, beta / tstat, se)
    }
  }
  z <- beta / se
  z[!is.finite(z)] <- NA_real_

  shrunken <- .shrink_effects(beta, se)

  score_vec <- switch(
    score,
    logFC = beta,
    standardized_effect = z,
    shrunken_effect = shrunken
  )

  if (!is.null(valid)) {
    score_vec[!valid] <- NA_real_
    beta[!valid] <- NA_real_
    se[!valid] <- NA_real_
    z[!valid] <- NA_real_
    shrunken[!valid] <- NA_real_
  }

  # Per-bin QL test (diagnostic; not the region-level confirmatory test)
  pval <- as.numeric(tt$table$PValue)
  if (!is.null(valid)) pval[!valid] <- NA_real_

  list(
    logFC = beta,
    se = se,
    z = z,
    shrunken = shrunken,
    score = score_vec,
    score_name = score,
    pvalue_bin = pval,
    dispersion_common = y$common.dispersion %||% NA_real_,
    dispersion_trended = y$trended.dispersion,
    dispersion_tagwise = y$tagwise.dispersion,
    AveLogCPM = edgeR::aveLogCPM(y),
    design = design,
    control_level = control_level,
    treatment_level = treatment_level,
    model = "edgeR_ql"
  )
}

.empty_effect <- function(n, score, valid) {
  na <- rep(NA_real_, n)
  if (!is.null(valid)) {
    # keep NAs
  }
  list(
    logFC = na, se = na, z = na, shrunken = na, score = na,
    score_name = score, pvalue_bin = na,
    dispersion_common = NA_real_,
    dispersion_trended = na, dispersion_tagwise = na,
    AveLogCPM = na, design = NULL,
    control_level = NA_character_, treatment_level = NA_character_,
    model = "edgeR_ql"
  )
}

.shrink_effects <- function(beta, se) {
  ok <- is.finite(beta) & is.finite(se) & se > 0
  out <- beta
  if (sum(ok) < 5L) return(out)
  # Simple EB: tau^2 = max(0, var(beta) - mean(se^2))
  tau2 <- max(0, stats::var(beta[ok]) - mean(se[ok]^2))
  w <- tau2 / (tau2 + se^2)
  w[!is.finite(w)] <- 0
  out[ok] <- beta[ok] * w[ok]
  out
}

#' Bind bins and effect estimates into a track table
#'
#' @param bins Bin table.
#' @param effect Output of [estimate_differential_effect()].
#' @return data.frame with chrom, start, end, effect, uncertainty, score, valid.
#' @export
write_effect_track <- function(bins, effect) {
  data.frame(
    chrom = bins$chrom,
    start = bins$start,
    end = bins$end,
    effect = as.numeric(effect$logFC),
    uncertainty = as.numeric(effect$se),
    score = as.numeric(effect$score),
    valid = as.logical(bins$valid) & is.finite(effect$score),
    stringsAsFactors = FALSE
  )
}
