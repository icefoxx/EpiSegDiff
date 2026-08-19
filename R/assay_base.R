.allowed_assay_families <- c("count", "continuous")
.allowed_input_types <- c("bam", "bigwig")
.allowed_signal_types <- c(
  "raw_coverage", "CPM", "RPKM", "RPGC", "fold_enrichment",
  "log2_ratio", "spikein_normalized", "other", "unknown"
)

#' Validate assay family / input type pairing
#'
#' @param assay_family `count` or `continuous`.
#' @param input_type `bam` or `bigwig`.
#' @export
validate_assay_contract <- function(assay_family, input_type) {
  if (!assay_family %in% .allowed_assay_families) {
    stop("assay_family must be 'count' or 'continuous'.", call. = FALSE)
  }
  if (!input_type %in% .allowed_input_types) {
    stop("input_type must be 'bam' or 'bigwig'.", call. = FALSE)
  }
  if (identical(assay_family, "count") && !identical(input_type, "bam")) {
    stop("count assay requires input_type = 'bam'.", call. = FALSE)
  }
  if (identical(assay_family, "continuous") && !identical(input_type, "bigwig")) {
    stop("continuous assay requires input_type = 'bigwig'.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Read a BigWig sample sheet
#'
#' Required columns: `sample` (or `sample_id`), `condition`, `replicate`, `file`.
#' @param path TSV path or data.frame.
#' @export
read_bigwig_sample_sheet <- function(path) {
  df <- if (is.data.frame(path)) path else {
    if (!file.exists(path)) stop("Sample sheet does not exist: ", path, call. = FALSE)
    utils::read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      check.names = FALSE, comment.char = "")
  }
  if ("sample_id" %in% names(df) && !"sample" %in% names(df)) {
    df$sample <- df$sample_id
  }
  need <- c("sample", "condition", "replicate", "file")
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    stop("BigWig sample sheet missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(df$sample)) stop("Duplicated sample IDs.", call. = FALSE)
  if (any(!nzchar(as.character(df$sample)))) stop("Empty sample IDs.", call. = FALSE)
  df$sample <- as.character(df$sample)
  df$condition <- as.character(df$condition)
  df$replicate <- as.integer(df$replicate)
  df$file <- as.character(df$file)
  df
}

#' Infer statistical capability level from replication
#' @export
infer_capability_level <- function(samples, allow_descriptive = FALSE) {
  tab <- table(samples$condition)
  if (length(tab) < 2L) {
    stop("At least two conditions are required.", call. = FALSE)
  }
  if (all(tab >= 2L)) return("replicated")
  if (isTRUE(allow_descriptive)) return("descriptive")
  stop("Insufficient replicates for inferential mode. Set allow_descriptive=TRUE for descriptive mode.",
       call. = FALSE)
}
