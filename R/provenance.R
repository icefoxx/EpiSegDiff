#' Collect software and package versions
#'
#' @return A named character vector.
#' @export
session_software_versions <- function() {
  si <- utils::sessionInfo()
  pkgs <- c(
    "EpiSegDiff", "yaml", "Rcpp", "edgeR", "limma", "GenomicRanges",
    "GenomicAlignments", "Rsamtools", "IRanges", "rtracklayer",
    "DESeq2", "data.table", "testthat"
  )
  vers <- vapply(pkgs, function(p) {
    if (p == "EpiSegDiff") {
      desc <- tryCatch(utils::packageDescription("EpiSegDiff"), error = function(e) NULL)
      if (is.null(desc) || is.na(desc$Version)) return("in-source")
      return(as.character(desc$Version))
    }
    if (requireNamespace(p, quietly = TRUE)) {
      as.character(utils::packageVersion(p))
    } else {
      "not-installed"
    }
  }, character(1))
  c(
    R = paste(si$R.version$major, si$R.version$minor, sep = "."),
    platform = si$platform,
    running = si$running %||% NA_character_,
    vers
  )
}

#' Start a provenance environment for a run
#'
#' @param config Validated config.
#' @param command_line Optional command string.
#' @return An environment used by [finalize_run_provenance()].
#' @export
start_run_provenance <- function(config, command_line = NULL) {
  env <- new.env(parent = emptyenv())
  env$t0 <- proc.time()[["elapsed"]]
  env$started_at <- as.character(Sys.time())
  env$config <- config
  env$seed <- as.integer(config$reproducibility$seed)
  set.seed(env$seed)
  env$command_line <- command_line %||% paste(commandArgs(trailingOnly = FALSE), collapse = " ")
  env$warnings <- character()
  env$errors <- character()
  env$input_ids <- character()
  withCallingHandlers(
    NULL,
    warning = function(w) NULL
  )
  env
}

#' Record input identifiers on a provenance environment
#' @keywords internal
.record_input_ids <- function(prov, paths, labels = NULL) {
  paths <- as.character(paths)
  if (is.null(labels)) labels <- paths
  ids <- vapply(paths, file_stable_id, character(1))
  names(ids) <- labels
  prov$input_ids <- c(prov$input_ids, ids)
  invisible(ids)
}

#' Write a key/value run manifest
#'
#' @param path Output TSV path.
#' @param entries Named character/numeric/logical vector or list.
#' @export
write_run_manifest <- function(path, entries) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  keys <- names(entries)
  if (is.null(keys)) stop("entries must be named.", call. = FALSE)
  vals <- vapply(entries, function(v) {
    if (is.null(v)) return("null")
    if (length(v) != 1L) return(paste(as.character(v), collapse = ";"))
    as.character(v)
  }, character(1))
  df <- data.frame(key = keys, value = vals, stringsAsFactors = FALSE)
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(df)
}

#' Finalize provenance and write `run_manifest.tsv`
#'
#' @param prov Environment from [start_run_provenance()].
#' @param output_dir Output directory.
#' @param extra Extra named entries.
#' @return The manifest data.frame.
#' @export
finalize_run_provenance <- function(prov, output_dir, extra = list()) {
  elapsed <- proc.time()[["elapsed"]] - prov$t0
  versions <- session_software_versions()
  entries <- c(
    list(
      started_at = prov$started_at,
      finished_at = as.character(Sys.time()),
      runtime_sec = elapsed,
      seed = prov$seed,
      command_line = prov$command_line,
      config_path = prov$config$`_config_path` %||% NA_character_,
      assay_type = prov$config$assay$type,
      assay_target = prov$config$assay$target,
      genome = prov$config$assay$genome,
      bin_size = prov$config$preprocessing$bin_size,
      counting_mode = prov$config$preprocessing$counting_mode,
      normalization_method = prov$config$normalization$method,
      differential_score = prov$config$differential$score,
      differential_model = prov$config$differential$model,
      dispersion_strategy = prov$config$differential$dispersion,
      segmentation_loss = prov$config$segmentation$loss,
      segmentation_engine = prov$config$segmentation$engine,
      candidate_strategy = prov$config$segmentation$candidate_strategy,
      testing_method = prov$config$testing$method,
      fdr_method = prov$config$testing$fdr_method,
      fdr_note = "BH after adaptive discovery is not guaranteed region-level FDR",
      boundary_refinement = prov$config$boundary$refinement_method,
      warnings = paste(prov$warnings, collapse = " | "),
      errors = paste(prov$errors, collapse = " | ")
    ),
    as.list(versions),
    as.list(prov$input_ids),
    extra
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_config(prov$config, file.path(output_dir, "config.resolved.yaml"))
  write_run_manifest(file.path(output_dir, "run_manifest.tsv"), entries)
}
