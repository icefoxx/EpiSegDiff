#' Analyze a pre-binned continuous matrix (no BigWig re-import)
#'
#' @param bed 0-based bin table (`chrom`,`start`,`end`).
#' @param signal bins × samples linear-scale matrix.
#' @param samples Sample sheet (`sample`,`condition`).
#' @export
analyze_continuous_matrix <- function(bed, signal, samples,
                                      numerator, denominator,
                                      bin_size,
                                      transform_method = "log2_pseudocount",
                                      pseudocount = 0.01,
                                      filter_method = "absolute",
                                      filter_threshold = 0.6,
                                      filter_quantile = 0.75,
                                      penalty = 24,
                                      min_bins = 3L,
                                      fdr_threshold = 0.05,
                                      allow_descriptive = FALSE,
                                      checkpoint_dir = NULL,
                                      resume = TRUE,
                                      workers = 1L,
                                      engine = "auto",
                                      identity = NULL) {
  t0 <- proc.time()[["elapsed"]]
  rss0 <- process_rss_mb()
  colnames(signal) <- samples$sample
  mode <- infer_capability_level(samples, allow_descriptive = allow_descriptive)
  filt <- filter_signal_bins(signal, method = filter_method,
                             threshold = filter_threshold,
                             quantile = filter_quantile)
  xf <- transform_continuous(signal, transform_method, pseudocount)
  eff <- estimate_continuous_effect(xf, samples, numerator, denominator, mode)
  score <- if (identical(mode, "descriptive")) {
    eff$effect
  } else {
    z <- eff$effect / pmax(abs(eff$se), 1e-6)
    z[!is.finite(z)] <- 0
    z
  }
  track <- data.frame(
    chrom = bed$chrom, start = as.integer(bed$start), end = as.integer(bed$end),
    score = as.numeric(score),
    valid = as.logical(filt) & is.finite(score),
    stringsAsFactors = FALSE
  )
  ident <- identity %||% list(
    genome = NA_character_,
    bin_size = bin_size,
    summary_fun = "mean",
    transform = transform_method,
    pseudocount = pseudocount,
    filter = filter_method,
    segmentation = "gaussian_pelt",
    penalty = penalty,
    min_bins = min_bins,
    engine = engine
  )
  sg <- segment_genome(
    track, penalty = penalty, min_width = as.integer(min_bins),
    checkpoint_dir = checkpoint_dir, resume = resume,
    workers = as.integer(workers), identity = ident,
    engine = engine, loss = "gaussian"
  )
  segs <- sg$segments
  if (!nrow(segs)) {
    tested <- standardize_primary_results(data.frame())
  } else {
    segs$region_id <- seq_len(nrow(segs))
    tested <- test_regions_limma(segs, xf, track, samples, numerator, denominator,
                                 inference_mode = mode, fdr_threshold = fdr_threshold)
    tested <- standardize_primary_results(tested, meta = list(
      assay_family = "continuous", input_type = "bigwig",
      contrast = paste(numerator, "-", denominator),
      inference_mode = mode,
      signal_filter_method = filter_method,
      transform_method = transform_method,
      segmentation_method = "gaussian_pelt"
    ))
  }
  list(
    tested = tested,
    n_input_bins = nrow(bed),
    n_valid = sum(track$valid),
    n_tested = nrow(tested),
    n_GAIN = sum(tested$level1 == "GAIN", na.rm = TRUE),
    n_LOSS = sum(tested$level1 == "LOSS", na.rm = TRUE),
    runtime_sec = proc.time()[["elapsed"]] - t0,
    bin_effect = data.frame(chrom = bed$chrom, start = bed$start, end = bed$end,
                            effect = eff$effect, stringsAsFactors = FALSE),
    inference_mode = mode,
    chromosome_runtime = sg$chromosome_runtime,
    peak_rss_mb = {
      rss1 <- process_rss_mb()
      if (is.finite(rss0) && is.finite(rss1)) max(rss0, rss1) else rss1
    }
  )
}

#' Limma small-variance diagnostics on a transformed matrix
#'
#' Characterizes (does not suppress) limma offset/zero-variance warnings.
#' Production fitting is unchanged; this is a diagnostic only.
#' @export
limma_variance_diagnostics <- function(transformed, samples, numerator, denominator,
                                       fdr = NULL, signal_linear = NULL,
                                       der_mask = NULL) {
  cond <- factor(samples$condition, levels = c(denominator, numerator))
  mu <- rowMeans(transformed, na.rm = TRUE)
  v <- rowSums((transformed - mu)^2, na.rm = TRUE) / pmax(rowSums(is.finite(transformed)) - 1, 1)
  tiny <- is.finite(v) & v < 1e-8
  out <- data.frame(
    n_units = nrow(transformed),
    n_tiny_var = sum(tiny, na.rm = TRUE),
    frac_tiny_var = mean(tiny, na.rm = TRUE),
    median_var = stats::median(v, na.rm = TRUE),
    q01_var = as.numeric(stats::quantile(v, 0.01, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  if (!is.null(signal_linear) && nrow(signal_linear) == length(tiny)) {
    mx <- do.call(pmax, c(as.list(as.data.frame(signal_linear)), list(na.rm = TRUE)))
    out$median_signal_tiny <- if (any(tiny)) stats::median(mx[tiny], na.rm = TRUE) else NA_real_
    out$median_signal_ok <- if (any(!tiny)) stats::median(mx[!tiny], na.rm = TRUE) else NA_real_
  } else {
    out$median_signal_tiny <- NA_real_
    out$median_signal_ok <- NA_real_
  }
  if (!is.null(fdr) && length(fdr) == length(tiny)) {
    out$frac_sig_among_tiny <- if (any(tiny)) mean(fdr[tiny] <= 0.05, na.rm = TRUE) else NA_real_
    out$frac_sig_among_ok <- if (any(!tiny)) mean(fdr[!tiny] <= 0.05, na.rm = TRUE) else NA_real_
  } else {
    out$frac_sig_among_tiny <- NA_real_
    out$frac_sig_among_ok <- NA_real_
  }
  if (!is.null(der_mask) && length(der_mask) == length(tiny)) {
    out$frac_tiny_among_der <- mean(tiny[der_mask], na.rm = TRUE)
    out$frac_tiny_among_nonder <- mean(tiny[!der_mask], na.rm = TRUE)
  } else {
    out$frac_tiny_among_der <- NA_real_
    out$frac_tiny_among_nonder <- NA_real_
  }
  warns <- character()
  effect_def <- effect_rob <- NULL
  if (requireNamespace("limma", quietly = TRUE) &&
      nlevels(droplevels(cond)) == 2L &&
      all(table(cond) >= 2L)) {
    design <- stats::model.matrix(~ cond)
    fit <- limma::lmFit(transformed, design)
    withCallingHandlers({
      fit_def <- limma::eBayes(fit, trend = FALSE, robust = FALSE)
      fit_rob <- tryCatch(
        limma::eBayes(fit, trend = TRUE, robust = TRUE),
        error = function(e) limma::eBayes(fit, trend = FALSE, robust = FALSE)
      )
      effect_def <- as.numeric(limma::topTable(fit_def, coef = 2, number = Inf,
                                               sort.by = "none")$logFC)
      effect_rob <- as.numeric(limma::topTable(fit_rob, coef = 2, number = Inf,
                                               sort.by = "none")$logFC)
    }, warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  }
  out$n_warnings_captured <- length(warns)
  out$warning_text <- if (length(warns)) paste(unique(warns), collapse = " | ") else ""
  if (!is.null(effect_def) && !is.null(effect_rob)) {
    out$effect_cor_default_vs_robust <- suppressWarnings(
      stats::cor(effect_def, effect_rob, use = "complete")
    )
    out$median_abs_effect_diff_robust <- stats::median(abs(effect_def - effect_rob),
                                                       na.rm = TRUE)
  } else {
    out$effect_cor_default_vs_robust <- NA_real_
    out$median_abs_effect_diff_robust <- NA_real_
  }
  influential <- isTRUE(out$frac_tiny_var > 0.25) &&
    isTRUE(out$frac_sig_among_tiny > 2 * out$frac_sig_among_ok) &&
    isTRUE(is.finite(out$frac_sig_among_ok))
  material <- isTRUE(out$frac_tiny_among_der > 0.2)
  out$classification <- if (material) {
    "materially influential"
  } else if (influential) {
    "potentially influential"
  } else {
    "harmless / expected"
  }
  out
}
