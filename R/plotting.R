#' Plot a differential-effect track
#'
#' @param track Effect track.
#' @param chrom Chromosome to plot.
#' @param segments Optional segment table.
#' @param file Optional PDF/PNG path.
#' @export
plot_effect_track <- function(track, chrom = NULL, segments = NULL, file = NULL) {
  if (!is.null(file)) {
    if (grepl("\\.png$", file, ignore.case = TRUE)) {
      grDevices::png(file, width = 1200, height = 400)
    } else {
      grDevices::pdf(file, width = 10, height = 3.5)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  if (!is.null(chrom)) track <- track[track$chrom == chrom, , drop = FALSE]
  if (!nrow(track)) {
    plot.new()
    return(invisible(NULL))
  }
  x <- (track$start + track$end) / 2
  y <- track$score
  y[!track$valid] <- NA
  graphics::plot(x, y, type = "l", xlab = "position", ylab = "score",
                 main = paste(unique(track$chrom), collapse = ","))
  graphics::abline(h = 0, lty = 2)
  if (!is.null(segments) && nrow(segments)) {
    if (!is.null(chrom)) segments <- segments[segments$chrom == chrom, ]
    graphics::rect(segments$start, min(y, na.rm = TRUE), segments$end,
                   max(y, na.rm = TRUE), border = NA, col = grDevices::rgb(1, 0, 0, 0.08))
  }
  invisible(NULL)
}

#' Plot segment means
#' @param segments Segment table.
#' @param file Optional path.
#' @export
plot_segments <- function(segments, file = NULL) {
  if (!is.null(file)) {
    grDevices::pdf(file, width = 8, height = 4)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  if (!nrow(segments)) {
    plot.new()
    return(invisible(NULL))
  }
  graphics::plot(segments$start, segments$effect_mean, pch = 16,
                 xlab = "start", ylab = "mean effect")
  invisible(NULL)
}

#' Plot bin-level dispersion versus abundance
#' @param effect Output of [estimate_differential_effect()].
#' @param file Optional path.
#' @export
plot_dispersion <- function(effect, file = NULL) {
  if (!is.null(file)) {
    grDevices::pdf(file, width = 5, height = 4)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  x <- effect$AveLogCPM
  y <- effect$dispersion_tagwise
  if (is.null(x) || is.null(y)) {
    plot.new()
    return(invisible(NULL))
  }
  graphics::plot(x, y, pch = 16, cex = 0.4, xlab = "AveLogCPM",
                 ylab = "tagwise dispersion")
  if (!is.null(effect$dispersion_trended)) {
    o <- order(x)
    graphics::lines(x[o], effect$dispersion_trended[o], col = "red")
  }
  invisible(NULL)
}
