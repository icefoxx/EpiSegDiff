#' Estimate normalization size factors / offsets
#'
#' Offsets are stored separately from visualization tracks.
#'
#' @param counts Integer matrix (bins × samples).
#' @param method `library_size`, `TMM`, `deseq2`, `spikein`, or `background`.
#' @param extreme_factor_ratio Warn if max/min factor exceeds this.
#' @param spikein_usage Must be TRUE to use `spikein`.
#' @param background_definition Required for `background`.
#' @param lib_sizes Optional library sizes; default column sums.
#' @param reference_bins Optional integer row indices used as the
#'   composition reference for `spikein` or as an override for
#'   `background`. These must be biologically/simulationally justified
#'   (e.g. dedicated spike-in bins), never arbitrary factors.
#' @return A list with `method`, `size_factors`, `offsets`, `raw_library_size`,
#'   `effective_library_size`, `warnings`.
#' @export
normalize_counts <- function(counts,
                             method = c("library_size", "TMM", "deseq2",
                                        "spikein", "background"),
                             extreme_factor_ratio = 5,
                             spikein_usage = FALSE,
                             background_definition = NULL,
                             lib_sizes = NULL,
                             reference_bins = NULL) {
  method <- match.arg(method)
  if (is.null(dim(counts))) {
    stop("counts must be a matrix.", call. = FALSE)
  }
  if (is.null(lib_sizes)) {
    lib_sizes <- colSums(counts)
  }
  lib_sizes <- as.numeric(lib_sizes)
  names(lib_sizes) <- colnames(counts)
  if (any(!is.finite(lib_sizes)) || any(lib_sizes < 0)) {
    stop("Library sizes must be finite and non-negative.", call. = FALSE)
  }
  if (any(lib_sizes == 0)) {
    stop("At least one sample has library size 0; cannot normalize.", call. = FALSE)
  }

  warns <- character()
  if (method == "library_size") {
    gm <- exp(mean(log(lib_sizes)))
    size_factors <- lib_sizes / gm
  } else if (method == "TMM") {
    .need_pkg("edgeR", "TMM normalization")
    y <- edgeR::DGEList(counts = counts)
    y <- edgeR::calcNormFactors(y, method = "TMM")
    # edgeR effective library size = lib.size * norm.factors
    size_factors <- (y$samples$lib.size * y$samples$norm.factors) /
      exp(mean(log(y$samples$lib.size * y$samples$norm.factors)))
    names(size_factors) <- colnames(counts)
  } else if (method == "deseq2") {
    size_factors <- .deseq2_size_factors(counts)
  } else if (method == "spikein") {
    if (!isTRUE(spikein_usage)) {
      stop("normalization.method = spikein requires spikein_usage = TRUE.",
           call. = FALSE)
    }
    if (is.null(reference_bins) || !length(reference_bins)) {
      stop("spikein normalization requires reference_bins (dedicated spike-in / control bins). Do not pass arbitrary size factors.",
           call. = FALSE)
    }
    size_factors <- .deseq2_size_factors(counts[as.integer(reference_bins), , drop = FALSE])
  } else if (method == "background") {
    size_factors <- .background_size_factors(
      counts,
      background_definition = background_definition,
      reference_bins = reference_bins
    )
  }

  if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
    stop("Normalization produced NA/Inf/non-positive factors.", call. = FALSE)
  }
  ratio <- max(size_factors) / min(size_factors)
  if (is.finite(extreme_factor_ratio) && ratio > extreme_factor_ratio) {
    msg <- sprintf(
      "Extreme size-factor ratio %.3f exceeds threshold %.3f. Global biological shifts can invalidate TMM/library-size normalization.",
      ratio, extreme_factor_ratio
    )
    warning(msg, call. = FALSE)
    warns <- c(warns, msg)
  }

  list(
    method = method,
    size_factors = as.numeric(size_factors),
    offsets = log(as.numeric(size_factors)),
    raw_library_size = lib_sizes,
    effective_library_size = lib_sizes * (as.numeric(size_factors) / (lib_sizes / exp(mean(log(lib_sizes))))),
    sample = colnames(counts),
    spikein_usage = isTRUE(spikein_usage),
    background_definition = background_definition,
    warnings = warns
  )
}

.background_size_factors <- function(counts, background_definition = NULL,
                                     reference_bins = NULL) {
  if (!is.null(reference_bins) && length(reference_bins)) {
    idx <- as.integer(reference_bins)
    idx <- idx[idx >= 1L & idx <= nrow(counts)]
    if (!length(idx)) {
      stop("background reference_bins are out of range.", call. = FALSE)
    }
    return(.deseq2_size_factors(counts[idx, , drop = FALSE]))
  }
  q <- 0.5
  if (is.null(background_definition)) {
    q <- 0.5
  } else if (is.numeric(background_definition) && length(background_definition) == 1L) {
    q <- as.numeric(background_definition)
  } else if (is.character(background_definition) &&
             background_definition %in% c("lower_half", "lower_quantile")) {
    q <- 0.5
  } else if (is.list(background_definition) && !is.null(background_definition$quantile)) {
    q <- as.numeric(background_definition$quantile)
  } else if (is.character(background_definition) && grepl("^[0-9.]+$", background_definition)) {
    q <- as.numeric(background_definition)
  } else {
    stop("background_definition must be a quantile in (0,1], 'lower_half', or list(quantile=...).",
         call. = FALSE)
  }
  if (!is.finite(q) || q <= 0 || q > 1) {
    stop("background quantile must be in (0, 1].", call. = FALSE)
  }
  row_mu <- rowMeans(as.matrix(counts))
  thr <- as.numeric(stats::quantile(row_mu, probs = q, na.rm = TRUE, names = FALSE))
  keep <- is.finite(row_mu) & row_mu <= thr
  if (sum(keep) < 10L) {
    warning("Fewer than 10 background bins; falling back to all-bin median-of-ratios.",
            call. = FALSE)
    return(.deseq2_size_factors(counts))
  }
  .deseq2_size_factors(counts[keep, , drop = FALSE])
}

.deseq2_size_factors <- function(counts) {
  # Median-of-ratios without requiring DESeq2. Geometric mean per bin
  # over samples; skip bins with a zero in any sample.
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  geo <- exp(rowMeans(log(counts)))
  # zeros → geo 0 → exclude
  keep <- is.finite(geo) & geo > 0
  if (!any(keep)) {
    stop("DESeq2-style size factors need at least one bin positive in all samples.",
         call. = FALSE)
  }
  ratios <- sweep(counts[keep, , drop = FALSE], 1, geo[keep], "/")
  sf <- apply(ratios, 2, stats::median)
  sf / exp(mean(log(sf)))
}

#' Write normalization diagnostics TSV
#'
#' @param norm Output of [normalize_counts()].
#' @param path Output path.
#' @return The diagnostic data.frame.
#' @export
write_normalization_diagnostics <- function(norm, path) {
  df <- data.frame(
    sample = norm$sample,
    raw_library_size = as.numeric(norm$raw_library_size),
    effective_library_size = as.numeric(norm$effective_library_size),
    normalization_factor = as.numeric(norm$size_factors),
    method = norm$method,
    stringsAsFactors = FALSE
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  df
}
