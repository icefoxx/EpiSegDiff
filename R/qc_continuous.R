#' Per-sample and replicate QC for a continuous matrix
#' @export
qc_continuous_signal <- function(signal, samples) {
  qs <- t(apply(signal, 2, function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(c(min = NA, q25 = NA, median = NA, mean = NA,
                             q75 = NA, max = NA, zero_frac = NA, nonfinite_frac = 1))
    q <- stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
    c(min = q[1], q25 = q[2], median = q[3], mean = mean(v), q75 = q[4], max = q[5],
      zero_frac = mean(v == 0), nonfinite_frac = mean(!is.finite(signal[, 1])) * NA)
  }))
  nf <- colMeans(!is.finite(signal))
  per <- data.frame(
    sample = colnames(signal),
    condition = samples$condition[match(colnames(signal), samples$sample)],
    min = qs[, "min"], q25 = qs[, "q25"], median = qs[, "median"],
    mean = qs[, "mean"], q75 = qs[, "q75"], max = qs[, "max"],
    zero_frac = qs[, "zero_frac"],
    nonfinite_frac = nf,
    stringsAsFactors = FALSE
  )
  pearson <- stats::cor(signal, use = "pairwise.complete.obs", method = "pearson")
  spearman <- stats::cor(signal, use = "pairwise.complete.obs", method = "spearman")
  list(per_sample = per, pearson = pearson, spearman = spearman)
}

#' Write continuous-assay QC plots
#' @export
plot_continuous_qc <- function(signal, samples, qc, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  pdf <- function(name, w = 6, h = 5) {
    grDevices::pdf(file.path(directory, name), width = w, height = h)
  }
  pdf("signal_distribution.pdf")
  graphics::boxplot(as.data.frame(signal), outline = FALSE, las = 2,
                    main = "Signal distribution", ylab = "score")
  grDevices::dev.off()
  # first two samples of the same condition if possible
  conds <- unique(samples$condition)
  if (length(conds)) {
    ids <- samples$sample[samples$condition == conds[1]]
    if (length(ids) >= 2L && all(ids[1:2] %in% colnames(signal))) {
      pdf("replicate_scatter.pdf")
      x <- signal[, ids[1]]; y <- signal[, ids[2]]
      ok <- is.finite(x) & is.finite(y)
      x <- x[ok]; y <- y[ok]
      if (length(x) > 40000L) {
        keep <- sample.int(length(x), 40000L)
        x <- x[keep]; y <- y[keep]
      }
      graphics::plot(x, y, pch = 16, cex = 0.25, col = grDevices::rgb(0, 0, 0, 0.15),
                     xlab = ids[1], ylab = ids[2],
                     main = paste("Replicate scatter", conds[1], "(subsampled)"))
      grDevices::dev.off()
    }
  }
  pdf("correlation_heatmap.pdf")
  graphics::image(qc$pearson, axes = FALSE, main = "Pearson correlation")
  graphics::axis(1, at = seq(0, 1, length.out = ncol(qc$pearson)),
                 labels = colnames(qc$pearson), las = 2, cex.axis = 0.7)
  graphics::axis(2, at = seq(0, 1, length.out = ncol(qc$pearson)),
                 labels = colnames(qc$pearson), las = 2, cex.axis = 0.7)
  grDevices::dev.off()
  pdf("condition_signal.pdf")
  cond <- samples$condition[match(colnames(signal), samples$sample)]
  graphics::boxplot(split(colMeans(signal), cond), ylab = "genome-wide mean",
                    main = "Condition-level signal")
  grDevices::dev.off()
  invisible(directory)
}
