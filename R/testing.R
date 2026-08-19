#' Region-level edgeR QL test
#'
#' BH-adjusted values are a nominal multiple-testing adjustment after
#' adaptive region discovery. They are not a guaranteed selective FDR.
#'
#' @param counts Integer matrix (regions × samples) or a data.frame from
#'   [recount_regions()] with sample columns.
#' @param samples Sample metadata.
#' @param offsets Optional sample offsets (log size factors).
#' @param control_level,treatment_level Condition levels.
#' @param fdr_method Passed to [stats::p.adjust()] (default `"BH"`).
#' @param candidate_source,segmentation_penalty,testing_method Provenance.
#' @return data.frame of test results.
#' @export
test_regions_edger <- function(counts, samples, offsets = NULL,
                               control_level = "control",
                               treatment_level = "treatment",
                               fdr_method = "BH",
                               candidate_source = NA_character_,
                               segmentation_penalty = NA_real_,
                               testing_method = "edgeR_ql") {
  .need_pkg("edgeR", "test_regions_edger()")
  meta <- NULL
  if (is.data.frame(counts) && "region_id" %in% names(counts)) {
    sample_ids <- samples$sample_id
    missing <- setdiff(sample_ids, names(counts))
    if (length(missing)) {
      stop("recount table missing sample columns: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    meta <- counts[, setdiff(names(counts), sample_ids), drop = FALSE]
    counts <- as.matrix(counts[, sample_ids, drop = FALSE])
    storage.mode(counts) <- "integer"
  }
  counts <- as.matrix(counts)
  if (ncol(counts) != nrow(samples)) {
    stop("counts columns and samples rows differ.", call. = FALSE)
  }

  cond <- .condition_factor(samples, control_level, treatment_level)
  design <- stats::model.matrix(~ cond)
  y <- edgeR::DGEList(counts = counts, samples = data.frame(condition = cond))
  if (!is.null(offsets)) {
    off <- matrix(as.numeric(offsets), nrow = nrow(counts), ncol = ncol(counts),
                  byrow = TRUE)
    y <- edgeR::scaleOffset(y, offset = off)
  } else {
    y <- edgeR::calcNormFactors(y)
  }

  nreg <- nrow(counts)
  if (nreg == 0L) {
    return(.empty_test_table())
  }
  if (nreg == 1L) {
    # estimateDisp needs care with 1 row; still fit QL if possible.
    y$common.dispersion <- 0.05
    y$trended.dispersion <- 0.05
    y$tagwise.dispersion <- 0.05
  } else {
    y <- tryCatch(
      edgeR::estimateDisp(y, design, robust = TRUE),
      error = function(e) {
        y$common.dispersion <- 0.05
        y$trended.dispersion <- rep(0.05, nreg)
        y$tagwise.dispersion <- rep(0.05, nreg)
        y
      }
    )
  }

  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  tt <- edgeR::glmQLFTest(fit, coef = 2)
  tab <- tt$table
  p <- as.numeric(tab$PValue)
  fdr <- stats::p.adjust(p, method = fdr_method)

  out <- data.frame(
    region_id = if (!is.null(meta) && "region_id" %in% names(meta)) meta$region_id else seq_len(nreg),
    logFC = as.numeric(tab$logFC),
    logCPM = as.numeric(tab$logCPM),
    statistic = as.numeric(tab$F),
    pvalue = p,
    FDR = fdr,
    testing_method = testing_method,
    fdr_method = fdr_method,
    candidate_source = candidate_source,
    segmentation_penalty = segmentation_penalty,
    fdr_interpretation = "nominal BH after adaptive discovery; not guaranteed selective FDR",
    stringsAsFactors = FALSE
  )
  if (!is.null(meta)) {
    extra <- setdiff(names(meta), names(out))
    out <- cbind(out, meta[, extra, drop = FALSE])
  }
  out
}

.empty_test_table <- function() {
  data.frame(
    region_id = integer(), logFC = numeric(), logCPM = numeric(),
    statistic = numeric(), pvalue = numeric(), FDR = numeric(),
    testing_method = character(), fdr_method = character(),
    candidate_source = character(), segmentation_penalty = numeric(),
    fdr_interpretation = character(),
    stringsAsFactors = FALSE
  )
}

#' Dispatch region testing by assay/method
#'
#' @param counts Count or methylation region table.
#' @param samples Sample metadata.
#' @param method Testing method.
#' @param ... Passed to the method-specific tester.
#' @export
test_regions <- function(counts, samples, method = c("edgeR_ql", "dss_style",
                                                     "beta_binomial"),
                         ...) {
  method <- match.arg(method)
  if (method == "edgeR_ql") {
    return(test_regions_edger(counts, samples, testing_method = method, ...))
  }
  if (method %in% c("dss_style", "beta_binomial")) {
    return(test_dmrs(counts, samples, method = method, ...))
  }
}
