.sample_required_cols <- c(
  "sample_id", "condition", "replicate", "assay", "target",
  "bam_or_methylation_path"
)

.sample_optional_cols <- c(
  "control_path", "spikein_status", "library_metadata"
)

#' Read a sample metadata table
#'
#' @param path TSV/CSV path.
#' @return A data.frame of samples.
#' @export
read_sample_table <- function(path) {
  if (!file.exists(path)) {
    stop("Sample table does not exist: ", path, call. = FALSE)
  }
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  df <- utils::read.table(
    path, header = TRUE, sep = sep, stringsAsFactors = FALSE,
    check.names = FALSE, quote = "", comment.char = ""
  )
  validate_sample_metadata(df)
}

#' Validate sample metadata
#'
#' Accepts a data.frame or a list of named lists (YAML `samples:`).
#'
#' @param samples Sample metadata.
#' @return A data.frame with required columns.
#' @export
validate_sample_metadata <- function(samples) {
  if (is.null(samples) || (is.list(samples) && length(samples) == 0L)) {
    stop("At least one sample is required.", call. = FALSE)
  }
  if (is.list(samples) && !is.data.frame(samples)) {
    if (!is.null(names(samples)) && all(c("sample_id", "condition") %in% names(samples))) {
      samples <- as.data.frame(samples, stringsAsFactors = FALSE)
    } else {
      rows <- lapply(samples, function(x) {
        if (!is.list(x)) stop("Each samples[] entry must be a named list.", call. = FALSE)
        as.data.frame(lapply(x, function(v) if (is.null(v)) NA else v),
                      stringsAsFactors = FALSE)
      })
      samples <- do.call(rbind, rows)
    }
  }
  if (!is.data.frame(samples)) {
    stop("samples must be a data.frame or a list of sample records.", call. = FALSE)
  }
  rownames(samples) <- NULL

  missing <- setdiff(.sample_required_cols, names(samples))
  if (length(missing)) {
    stop("Sample metadata missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  for (opt in .sample_optional_cols) {
    if (!opt %in% names(samples)) samples[[opt]] <- NA
  }

  samples$sample_id <- as.character(samples$sample_id)
  samples$condition <- as.character(samples$condition)
  samples$replicate <- as.integer(samples$replicate)
  samples$assay <- as.character(samples$assay)
  samples$target <- as.character(samples$target)
  samples$bam_or_methylation_path <- as.character(samples$bam_or_methylation_path)

  if (any(!nzchar(samples$sample_id) | is.na(samples$sample_id))) {
    stop("sample_id must be non-empty.", call. = FALSE)
  }
  if (anyDuplicated(samples$sample_id)) {
    stop("sample_id values must be unique.", call. = FALSE)
  }
  if (any(is.na(samples$condition) | !nzchar(samples$condition))) {
    stop("condition must be non-empty for every sample.", call. = FALSE)
  }
  if (any(is.na(samples$replicate) | samples$replicate < 1L)) {
    stop("replicate must be a positive integer.", call. = FALSE)
  }
  allowed_assay <- c("cutnrun", "chipseq", "atac", "methylation")
  if (any(!samples$assay %in% allowed_assay)) {
    stop("assay must be one of: ", paste(allowed_assay, collapse = ", "), call. = FALSE)
  }
  if (length(unique(samples$assay)) > 1L) {
    stop("A single run cannot mix assay types. Split the cohort.", call. = FALSE)
  }
  if (length(unique(samples$target)) > 1L) {
    warning("Multiple targets in one run; this is allowed but unusual.", call. = FALSE)
  }
  if (length(unique(samples$condition)) < 2L) {
    stop("At least two conditions are required for a differential run.", call. = FALSE)
  }
  if (length(unique(samples$condition)) > 2L) {
    stop("The v0.2 MVP supports exactly two conditions.", call. = FALSE)
  }

  samples
}

.condition_factor <- function(samples, control_level, treatment_level) {
  cond <- as.character(samples$condition)
  if (!control_level %in% cond) {
    stop("control_level '", control_level, "' is not present in samples.", call. = FALSE)
  }
  if (!treatment_level %in% cond) {
    stop("treatment_level '", treatment_level, "' is not present in samples.", call. = FALSE)
  }
  factor(cond, levels = c(control_level, treatment_level))
}
