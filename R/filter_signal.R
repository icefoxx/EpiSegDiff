#' Filter bins by continuous signal
#'
#' @param signal Bins × samples matrix (untransformed scores).
#' @param method `none`, `absolute`, `quantile`, `control_based`, `auto`.
#' @param threshold Absolute threshold (`absolute`).
#' @param quantile Quantile of per-bin max (`quantile`).
#' @param control_signal Optional control matrix (`control_based`).
#' @param control_mult Multiplier of control median for `control_based`.
#' @return Logical keep vector (TRUE = retain for segmentation).
#' @export
filter_signal_bins <- function(signal, method = c("none", "absolute", "quantile",
                                                  "control_based", "auto"),
                               threshold = 0.6, quantile = 0.75,
                               control_signal = NULL, control_mult = 2) {
  method <- match.arg(method)
  mx <- do.call(pmax, c(as.list(as.data.frame(signal)), list(na.rm = TRUE)))
  keep <- switch(
    method,
    none = rep(TRUE, nrow(signal)),
    absolute = is.finite(mx) & mx >= as.numeric(threshold),
    quantile = {
      q <- as.numeric(stats::quantile(mx[is.finite(mx)], probs = quantile,
                                      names = FALSE, na.rm = TRUE))
      is.finite(mx) & mx >= q
    },
    control_based = {
      if (is.null(control_signal)) {
        stop("control_based filter requires control_signal.", call. = FALSE)
      }
      cmx <- do.call(pmax, c(as.list(as.data.frame(control_signal)), list(na.rm = TRUE)))
      thr <- stats::median(cmx[is.finite(cmx)], na.rm = TRUE) * as.numeric(control_mult)
      is.finite(mx) & mx >= thr
    },
    auto = {
      # Deterministic: keep bins above the 75th percentile of per-bin max.
      q <- as.numeric(stats::quantile(mx[is.finite(mx)], 0.75, names = FALSE, na.rm = TRUE))
      is.finite(mx) & mx >= q
    }
  )
  attr(keep, "filter_method") <- method
  attr(keep, "filter_threshold") <- switch(
    method,
    none = NA_real_,
    absolute = as.numeric(threshold),
    quantile = as.numeric(quantile),
    control_based = as.numeric(control_mult),
    auto = 0.75
  )
  keep
}
