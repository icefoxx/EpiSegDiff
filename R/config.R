#' Default EpiSegDiff configuration
#'
#' Every scientific choice is explicit. See `config/example.yaml`.
#'
#' @return A nested list.
#' @export
default_config <- function() {
  list(
    assay = list(
      type = "cutnrun",
      assay_family = "count",
      input_type = "bam",
      target = "unknown",
      genome = "unknown"
    ),
    samples = list(),
    preprocessing = list(
      bin_size = 200L,
      chrom_sizes = NULL,
      blacklist = NULL,
      extra_exclude = NULL,
      counting_mode = "fragment_midpoint",
      fragment_channels = list(),
      paired = TRUE,
      mapq_min = 0L,
      chromosomes = NULL,
      mask_low_coverage = FALSE,
      min_count_sum = 0L,
      coordinate_system = "bed_0based"
    ),
    normalization = list(
      method = "TMM",
      spikein_usage = FALSE,
      background_definition = NULL,
      extreme_factor_ratio = 5
    ),
    differential = list(
      model = "edgeR_ql",
      design = "~ condition",
      control_level = "control",
      treatment_level = "treatment",
      score = "standardized_effect",
      dispersion = "edgeR_estimateDisp",
      prior_count = 0.125
    ),
    segmentation = list(
      loss = "gaussian",
      engine = "pelt",
      penalties = c(5, 10, 20, 40, 80),
      min_width = 3L,
      huber_k = 1.345,
      winsor_k = 3,
      candidate_strategy = "single_selected_penalty",
      selected_penalty = 20
    ),
    testing = list(
      method = "edgeR_ql",
      fdr_method = "BH",
      fdr_threshold = 0.05,
      effect_threshold = 0,
      report_empirical_fdr = TRUE
    ),
    boundary = list(
      refinement_method = "local_gradient",
      search_window_bins = 5L
    ),
    topology = list(
      enabled = TRUE,
      domain_method = "condition_mean_threshold",
      domain_mean_quantile = 0.8,
      min_reciprocal_overlap = 0.2,
      high_reciprocal_overlap = 0.5,
      max_match_distance = 5000L,
      min_width_ratio_change = 1.25,
      min_boundary_shift = 100L,
      min_internal_gap_width = 200L,
      similar_width_ratio = 1.25,
      max_overlap_for_shift = 0.85
    ),
    output = list(
      directory = "results",
      write_bigwig = FALSE,
      write_bed = TRUE,
      write_tsv = TRUE
    ),
    reproducibility = list(
      seed = 1L
    ),
    atac = list(
      tn5_plus_offset = 4L,
      tn5_minus_offset = -5L,
      fragment_channel_breaks = list(
        short = c(0, 100),
        mono = c(180, 247),
        di = c(315, 473)
      )
    ),
    methylation = list(
      min_coverage = 5L,
      model = "beta_binomial_approx",
      score = "standardized_effect"
    )
  )
}

#' Read a YAML configuration and apply defaults
#'
#' @param path Path to a YAML file.
#' @return A validated config list.
#' @export
read_episegdiff_config <- function(path) {
  if (!file.exists(path)) {
    stop("Configuration file does not exist: ", path, call. = FALSE)
  }
  raw <- yaml::read_yaml(path)
  cfg <- .merge_config(default_config(), raw)
  cfg$`_config_path` <- normalizePath(path, winslash = "/", mustWork = FALSE)
  validate_config(cfg)
}

#' Write a configuration list to YAML
#'
#' @param config Config list.
#' @param path Output path.
#' @export
write_config <- function(config, path) {
  drop <- grepl("^_", names(config))
  yaml::write_yaml(config[!drop], path)
  invisible(path)
}

.merge_config <- function(default, override) {
  if (is.null(override)) return(default)
  if (!is.list(override) || is.null(names(override))) return(override)
  if (!is.list(default)) return(override)
  out <- default
  for (nm in names(override)) {
    if (nm %in% names(default) && is.list(default[[nm]]) && is.list(override[[nm]]) &&
        !nm %in% c("samples", "fragment_channels", "penalties")) {
      out[[nm]] <- .merge_config(default[[nm]], override[[nm]])
    } else {
      out[[nm]] <- override[[nm]]
    }
  }
  out
}

#' Validate a configuration list
#'
#' @param config Config list (typically from [read_episegdiff_config()]).
#' @return The config, invisibly, after checks.
#' @export
validate_config <- function(config) {
  req <- c("assay", "samples", "preprocessing", "normalization",
           "differential", "segmentation", "testing", "boundary",
           "topology", "output", "reproducibility")
  missing <- setdiff(req, names(config))
  if (length(missing)) {
    stop("Config missing top-level keys: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  at <- config$assay$type
  if (!at %in% c("cutnrun", "chipseq", "atac", "methylation")) {
    stop("assay.type must be cutnrun, chipseq, atac, or methylation.", call. = FALSE)
  }

  cm <- config$preprocessing$counting_mode
  if (!cm %in% c("fragment_midpoint", "fragment_overlap", "insertion")) {
    stop("preprocessing.counting_mode is not recognized.", call. = FALSE)
  }
  if (!identical(config$preprocessing$coordinate_system, "bed_0based")) {
    stop("Only coordinate_system = 'bed_0based' is supported.", call. = FALSE)
  }
  bs <- as.integer(config$preprocessing$bin_size)
  if (is.na(bs) || bs < 1L) {
    stop("preprocessing.bin_size must be a positive integer.", call. = FALSE)
  }
  config$preprocessing$bin_size <- bs

  nm <- config$normalization$method
  if (!nm %in% c("library_size", "TMM", "deseq2", "spikein", "background")) {
    stop("normalization.method is not recognized.", call. = FALSE)
  }
  if (nm %in% c("spikein", "background")) {
    # Allowed in schema; implementation errors unless fully specified later.
  }

  sc <- config$differential$score
  if (!sc %in% c("logFC", "standardized_effect", "shrunken_effect")) {
    stop("differential.score must be logFC, standardized_effect, or shrunken_effect.",
         call. = FALSE)
  }

  loss <- config$segmentation$loss
  if (!loss %in% c("gaussian", "huber", "winsorized_gaussian")) {
    stop("segmentation.loss is not recognized.", call. = FALSE)
  }
  if (!config$segmentation$engine %in% c("pelt", "fpop", "op")) {
    stop("segmentation.engine must be pelt, fpop, or op.", call. = FALSE)
  }
  if (!config$segmentation$candidate_strategy %in%
      c("single_selected_penalty", "multi-scale_consensus")) {
    stop("segmentation.candidate_strategy is not recognized.", call. = FALSE)
  }
  config$segmentation$penalties <- sort(unique(as.numeric(config$segmentation$penalties)))
  if (any(!is.finite(config$segmentation$penalties)) ||
      any(config$segmentation$penalties < 0)) {
    stop("segmentation.penalties must be finite and non-negative.", call. = FALSE)
  }

  if (!config$testing$method %in% c("edgeR_ql", "dss_style", "beta_binomial")) {
    stop("testing.method is not recognized.", call. = FALSE)
  }
  if (identical(at, "methylation") && identical(config$testing$method, "edgeR_ql")) {
    stop("Methylation assays cannot use testing.method = edgeR_ql.", call. = FALSE)
  }
  if (!identical(at, "methylation") && config$testing$method %in% c("dss_style", "beta_binomial")) {
    stop("Count assays cannot use a methylation testing method.", call. = FALSE)
  }

  if (!config$boundary$refinement_method %in%
      c("local_gradient", "local_likelihood", "cross_scale")) {
    stop("boundary.refinement_method is not recognized.", call. = FALSE)
  }

  seed <- as.integer(config$reproducibility$seed)
  if (is.na(seed)) stop("reproducibility.seed must be an integer.", call. = FALSE)
  config$reproducibility$seed <- seed

  if (length(config$samples)) {
    config$samples <- validate_sample_metadata(config$samples)
  }

  config
}
