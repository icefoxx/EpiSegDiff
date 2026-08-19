#' Per-bin continuous differential effect (limma when replicated)
#'
#' @param x Transformed bins × samples matrix.
#' @param samples Sample table with `condition`.
#' @param numerator,denominator Contrast: numerator − denominator.
#' @param inference_mode `replicated` or `descriptive`.
#' @export
estimate_continuous_effect <- function(x, samples, numerator, denominator,
                                       inference_mode = c("replicated", "descriptive")) {
  inference_mode <- match.arg(inference_mode)
  if (!numerator %in% samples$condition) {
    stop("Contrast numerator '", numerator, "' is not in samples.", call. = FALSE)
  }
  if (!denominator %in% samples$condition) {
    stop("Contrast denominator '", denominator, "' is not in samples.", call. = FALSE)
  }
  cond <- factor(samples$condition, levels = c(denominator, numerator))
  if (nlevels(droplevels(cond)) < 2L) {
    stop("Contrast is not estimable.", call. = FALSE)
  }
  a <- x[, cond == denominator, drop = FALSE]
  b <- x[, cond == numerator, drop = FALSE]
  mean_a <- rowMeans(a)
  mean_b <- rowMeans(b)
  effect <- mean_b - mean_a
  if (identical(inference_mode, "descriptive")) {
    return(list(
      effect = effect, se = rep(NA_real_, nrow(x)),
      t = rep(NA_real_, nrow(x)), pvalue = rep(NA_real_, nrow(x)),
      mean_a = mean_a, mean_b = mean_b,
      n_a = ncol(a), n_b = ncol(b),
      inference_mode = "descriptive",
      contrast = paste(numerator, "-", denominator),
      model = "none"
    ))
  }
  .need_pkg("limma", "estimate_continuous_effect()")
  if (ncol(a) < 2L || ncol(b) < 2L) {
    stop("Replicated mode needs >=2 samples per contrast condition.", call. = FALSE)
  }
  design <- stats::model.matrix(~ cond)
  fit <- limma::lmFit(x, design)
  fit <- tryCatch(
    limma::eBayes(fit, trend = nrow(x) >= 10L, robust = nrow(x) >= 10L),
    error = function(e) limma::eBayes(fit, trend = FALSE, robust = FALSE)
  )
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  list(
    effect = as.numeric(tt$logFC),
    se = as.numeric(tt$logFC / tt$t),
    t = as.numeric(tt$t),
    pvalue = as.numeric(tt$P.Value),
    mean_a = mean_a, mean_b = mean_b,
    n_a = ncol(a), n_b = ncol(b),
    inference_mode = "replicated",
    contrast = paste(numerator, "-", denominator),
    model = "limma"
  )
}

#' Test candidate regions with limma on region-mean transformed signal
#' @export
test_regions_limma <- function(intervals, transformed, track, samples,
                               numerator, denominator,
                               inference_mode = "replicated",
                               fdr_threshold = 0.05) {
  cond <- factor(samples$condition, levels = c(denominator, numerator))
  rmat <- .region_means_matrix(intervals, transformed, track)
  if (identical(inference_mode, "descriptive")) {
    ea <- rowMeans(rmat[, cond == denominator, drop = FALSE])
    eb <- rowMeans(rmat[, cond == numerator, drop = FALSE])
    logFC <- eb - ea
    p <- FDR <- rep(NA_real_, nrow(intervals))
    stat <- se <- rep(NA_real_, nrow(intervals))
    model <- "descriptive"
  } else {
    .need_pkg("limma", "test_regions_limma()")
    design <- stats::model.matrix(~ cond)
    fit <- limma::lmFit(rmat, design)
    fit <- tryCatch(
      limma::eBayes(fit, trend = nrow(rmat) >= 10L, robust = nrow(rmat) >= 10L),
      error = function(e) limma::eBayes(fit, trend = FALSE, robust = FALSE)
    )
    tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
    logFC <- as.numeric(tt$logFC)
    stat <- as.numeric(tt$t)
    se <- logFC / stat
    p <- as.numeric(tt$P.Value)
    FDR <- as.numeric(tt$adj.P.Val)
    ea <- rowMeans(rmat[, cond == denominator, drop = FALSE])
    eb <- rowMeans(rmat[, cond == numerator, drop = FALSE])
    model <- "limma"
  }
  dir <- ifelse(logFC > 0, "GAIN", ifelse(logFC < 0, "LOSS", "STABLE"))
  called <- if (identical(inference_mode, "descriptive")) {
    rep(FALSE, length(logFC))
  } else {
    !is.na(FDR) & FDR <= fdr_threshold
  }
  data.frame(
    chrom = intervals$chrom, start = intervals$start, end = intervals$end,
    region_id = if ("region_id" %in% names(intervals)) intervals$region_id else seq_len(nrow(intervals)),
    effect = logFC, log2FC = logFC, p_value = p, FDR = FDR,
    statistic = stat, se = se,
    direction = dir,
    level1 = ifelse(called, dir, "STABLE"),
    mean_group1 = ea, mean_group2 = eb,
    n_group1 = sum(cond == denominator), n_group2 = sum(cond == numerator),
    contrast = paste(numerator, "-", denominator),
    inference_mode = inference_mode,
    testing_method = model,
    stringsAsFactors = FALSE
  )
}

.region_means_matrix <- function(intervals, x, track) {
  .need_pkg("GenomicRanges", ".region_means_matrix")
  gr_i <- GenomicRanges::GRanges(
    intervals$chrom,
    IRanges::IRanges(as.integer(intervals$start) + 1L, as.integer(intervals$end))
  )
  gr_b <- GenomicRanges::GRanges(
    track$chrom,
    IRanges::IRanges(as.integer(track$start) + 1L, as.integer(track$end))
  )
  ov <- GenomicRanges::findOverlaps(gr_i, gr_b)
  out <- matrix(NA_real_, nrow(intervals), ncol(x), dimnames = list(NULL, colnames(x)))
  if (!length(ov)) return(out)
  spl <- split(S4Vectors::subjectHits(ov), S4Vectors::queryHits(ov))
  for (nm in names(spl)) {
    i <- as.integer(nm)
    out[i, ] <- colMeans(x[spl[[nm]], , drop = FALSE])
  }
  out
}
