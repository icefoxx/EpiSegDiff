#' Transform a continuous signal matrix
#'
#' @param x Numeric matrix.
#' @param method `log2_pseudocount`, `log1p`, or `none`.
#' @param pseudocount Added before log2 when `log2_pseudocount`.
#' @export
transform_continuous <- function(x, method = c("log2_pseudocount", "log1p", "none"),
                                 pseudocount = 0.01) {
  method <- match.arg(method)
  x <- as.matrix(x)
  out <- switch(
    method,
    none = x,
    log1p = log1p(pmax(x, 0)),
    log2_pseudocount = {
      if (!is.finite(pseudocount) || pseudocount < 0) {
        stop("pseudocount must be finite and non-negative.", call. = FALSE)
      }
      log2(pmax(x, 0) + pseudocount)
    }
  )
  if (any(!is.finite(out))) {
    nbad <- sum(!is.finite(out))
    stop("Transformation produced ", nbad, " non-finite values.", call. = FALSE)
  }
  attr(out, "transform_method") <- method
  attr(out, "pseudocount") <- if (method == "log2_pseudocount") pseudocount else NA_real_
  out
}
