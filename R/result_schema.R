.primary_schema_cols <- c(
  "chrom", "start", "end", "region_id", "assay_family", "input_type",
  "contrast", "effect", "log2FC", "p_value", "FDR", "direction",
  "mean_group1", "mean_group2", "n_group1", "n_group2",
  "inference_mode", "signal_filter_method", "transform_method",
  "segmentation_method", "level1"
)

#' Standardize a Level I result table
#' @export
standardize_primary_results <- function(tested, meta = list()) {
  if (!nrow(tested)) {
    out <- as.data.frame(setNames(replicate(length(.primary_schema_cols),
                                            logical(0), simplify = FALSE),
                                  .primary_schema_cols))
    return(out)
  }
  out <- tested
  rename <- c(pvalue = "p_value", logFC = "log2FC")
  for (nm in names(rename)) {
    if (nm %in% names(out) && !rename[[nm]] %in% names(out)) {
      out[[rename[[nm]]]] <- out[[nm]]
    }
  }
  if (!"effect" %in% names(out) && "log2FC" %in% names(out)) out$effect <- out$log2FC
  if (!"direction" %in% names(out) && "level1" %in% names(out)) {
    out$direction <- out$level1
  }
  for (k in names(meta)) {
    if (!k %in% names(out)) out[[k]] <- meta[[k]]
  }
  for (col in .primary_schema_cols) {
    if (!col %in% names(out)) out[[col]] <- NA
  }
  out[, unique(c(.primary_schema_cols, setdiff(names(out), .primary_schema_cols))),
      drop = FALSE]
}

#' Validate required primary columns
#' @export
validate_primary_schema <- function(df) {
  missing <- setdiff(.primary_schema_cols, names(df))
  if (length(missing)) {
    stop("Primary result missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}
