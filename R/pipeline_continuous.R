#' Run the native continuous BigWig assay
#'
#' Uses limma on transformed continuous signal. Does not call edgeR.
#'
#' @param config A v0.4-style list or YAML path (see `configs/gse280947_bigwig.yaml`).
#' @export
run_continuous_assay <- function(config) {
  if (is.character(config) && length(config) == 1L) {
    config <- read_continuous_config(config)
  }
  validate_assay_contract(config$assay$assay_family, config$assay$input_type)
  outdir <- config$output$directory
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  t0 <- proc.time()[["elapsed"]]

  sheet <- if (is.character(config$samples)) {
    read_bigwig_sample_sheet(config$samples)
  } else {
    read_bigwig_sample_sheet(config$samples)
  }
  num <- config$contrast$numerator
  den <- config$contrast$denominator
  allow_desc <- isTRUE(config$model$allow_descriptive)
  mode <- infer_capability_level(sheet, allow_descriptive = allow_desc)

  assay <- read_bigwig_assay(
    sheet,
    genome = config$assay$genome,
    chromosomes = config$binning$chromosomes,
    bin_size = config$binning$bin_size,
    summary_fun = config$binning$summary %||% "mean"
  )
  sig <- assay$signal
  if (all(sig == 0) || !any(is.finite(sig))) {
    stop("All-zero or empty continuous signal.", call. = FALSE)
  }

  filt <- filter_signal_bins(
    sig,
    method = config$signal_filter$method %||% "none",
    threshold = config$signal_filter$threshold %||% 0.6,
    quantile = config$signal_filter$quantile %||% 0.75,
    control_signal = config$signal_filter$control_signal
  )
  xf <- transform_continuous(
    sig,
    method = config$transform$method %||% "log2_pseudocount",
    pseudocount = config$transform$pseudocount %||% 0.01
  )
  bed <- granges_to_bed_frame(assay$bins)
  eff <- estimate_continuous_effect(
    xf, sheet, numerator = num, denominator = den, inference_mode = mode
  )
  score <- if (identical(mode, "descriptive")) {
    eff$effect
  } else {
    z <- eff$effect / pmax(abs(eff$se), 1e-6)
    z[!is.finite(z)] <- 0
    z
  }
  track <- data.frame(
    chrom = bed$chrom, start = bed$start, end = bed$end,
    score = as.numeric(score), valid = as.logical(filt) & is.finite(score),
    stringsAsFactors = FALSE
  )
  segs <- segment_track(
    track,
    penalty = as.numeric(config$segmentation$penalty %||% 24),
    min_width = as.integer(config$segmentation$min_bins %||% 3L),
    loss = "gaussian", engine = "pelt"
  )
  if (!nrow(segs)) {
    tested <- standardize_primary_results(
      data.frame(),
      meta = list(assay_family = "continuous", input_type = "bigwig",
                  contrast = paste(num, "-", den), inference_mode = mode)
    )
  } else {
    segs$region_id <- seq_len(nrow(segs))
    tested <- test_regions_limma(
      segs, xf, track, sheet, num, den,
      inference_mode = mode,
      fdr_threshold = config$calling$fdr %||% 0.05
    )
    tested <- standardize_primary_results(tested, meta = list(
      assay_family = "continuous",
      input_type = "bigwig",
      contrast = paste(num, "-", den),
      inference_mode = mode,
      signal_filter_method = config$signal_filter$method %||% "none",
      transform_method = config$transform$method %||% "log2_pseudocount",
      segmentation_method = "gaussian_pelt"
    ))
  }
  validate_primary_schema(tested)

  qc <- qc_continuous_signal(sig, data.frame(sample = sheet$sample,
                                             condition = sheet$condition))
  utils::write.table(qc$per_sample, file.path(outdir, "sample_qc.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  tryCatch(plot_continuous_qc(sig, data.frame(sample = sheet$sample,
                                              condition = sheet$condition),
                              qc, file.path(outdir, "qc_plots")),
           error = function(e) warning("QC plots failed: ", e$message, call. = FALSE))

  utils::write.table(tested, file.path(outdir, "episegdiff_primary.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  sigrows <- tested[!is.na(tested$level1) & tested$level1 %in% c("GAIN", "LOSS"), ]
  if (nrow(sigrows)) {
    utils::write.table(
      sigrows[, c("chrom", "start", "end", "region_id", "direction")],
      file.path(outdir, "episegdiff_primary.bed"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
  }
  yaml::write_yaml(config, file.path(outdir, "analysis_parameters.yaml"))
  writeLines(utils::capture.output(utils::sessionInfo()),
             file.path(outdir, "session_info.txt"))
  sig_type <- config$assay$signal_type %||% "unknown"
  if (identical(sig_type, "unknown")) {
    warning("signal_type is 'unknown'; recorded in provenance.", call. = FALSE)
  }
  man <- list(
    software = paste0("EpiSegDiff ", as.character(utils::packageVersion("EpiSegDiff"))),
    date = as.character(Sys.time()),
    assay_family = "continuous",
    input_type = "bigwig",
    genome = config$assay$genome,
    signal_type = sig_type,
    normalization = config$assay$normalization %||% NA,
    bin_size = config$binning$bin_size,
    transform = config$transform$method,
    pseudocount = config$transform$pseudocount %||% NA,
    filter_method = config$signal_filter$method,
    filter_threshold = config$signal_filter$threshold %||% NA,
    model = if (mode == "replicated") "limma" else "none",
    contrast = paste(num, "-", den),
    segmentation_method = "gaussian_pelt",
    penalty = config$segmentation$penalty,
    min_bins = config$segmentation$min_bins,
    fdr_threshold = config$calling$fdr,
    inference_mode = mode,
    n_bins = nrow(track),
    n_valid = sum(track$valid),
    n_segments = nrow(segs),
    n_GAIN = sum(tested$level1 == "GAIN", na.rm = TRUE),
    n_LOSS = sum(tested$level1 == "LOSS", na.rm = TRUE),
    runtime_sec = proc.time()[["elapsed"]] - t0,
    inputs = paste(sheet$file, collapse = ";")
  )
  yaml::write_yaml(man, file.path(outdir, "run_manifest.yaml"))
  list(assay = assay, track = track, segments = segs, tested = tested,
       qc = qc, config = config, manifest = man)
}

#' Read a v0.4 continuous-assay YAML config
#' @export
read_continuous_config <- function(path) {
  raw <- yaml::read_yaml(path)
  if (is.null(raw$assay$assay_family)) raw$assay$assay_family <- "continuous"
  if (is.null(raw$assay$input_type)) raw$assay$input_type <- "bigwig"
  validate_assay_contract(raw$assay$assay_family, raw$assay$input_type)
  if (!is.null(raw$assay$signal_type) &&
      !raw$assay$signal_type %in% .allowed_signal_types) {
    stop("Unknown signal_type: ", raw$assay$signal_type, call. = FALSE)
  }
  raw$binning$bin_size <- as.integer(raw$binning$bin_size %||% 200L)
  raw$segmentation$penalty <- as.numeric(raw$segmentation$penalty %||% 24)
  raw$segmentation$min_bins <- as.integer(raw$segmentation$min_bins %||% 3L)
  raw$calling$fdr <- as.numeric(raw$calling$fdr %||% 0.05)
  raw$signal_filter$method <- raw$signal_filter$method %||% "none"
  raw$transform$method <- raw$transform$method %||% "log2_pseudocount"
  raw$transform$pseudocount <- raw$transform$pseudocount %||% 0.01
  raw$output$directory <- raw$output$directory %||% "results"
  raw
}
