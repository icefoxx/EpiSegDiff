#' Run the EpiSegDiff v0.2 pipeline
#'
#' Count-assay path: BAM/fragments → counts → normalization →
#' differential-effect score → multi-scale segmentation → recount →
#' edgeR QL testing → boundary refinement → topology classification.
#'
#' Methylation path uses the methylation adapter and never the NB count
#' model.
#'
#' @param config Config list or YAML path.
#' @param count_matrix Optional prebuilt list from [build_count_matrix()]
#'   or [simulate_count_assay()] (`counts`, `bins`, optional `samples`).
#' @param command_line Optional provenance string.
#' @return A list with all intermediate objects. Side effect: writes
#'   outputs under `config$output$directory`.
#' @export
run_episegdiff <- function(config, count_matrix = NULL, command_line = NULL) {
  if (is.character(config) && length(config) == 1L) {
    config <- read_episegdiff_config(config)
  } else {
    config <- validate_config(.merge_config(default_config(), config))
  }
  outdir <- config$output$directory
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  prov <- start_run_provenance(config, command_line = command_line)

  handler <- function(w) {
    prov$warnings <- c(prov$warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }

  result <- withCallingHandlers({
    if (identical(config$assay$type, "methylation")) {
      .run_methylation_pipeline(config, count_matrix, prov)
    } else {
      .run_count_pipeline(config, count_matrix, prov)
    }
  }, warning = handler)

  extra <- list(
    n_bins = if (!is.null(result$bins)) nrow(result$bins) else NA,
    n_candidates = if (!is.null(result$candidates)) nrow(result$candidates) else NA,
    n_tested = if (!is.null(result$tested)) nrow(result$tested) else NA,
    n_sig = if (!is.null(result$tested) && "FDR" %in% names(result$tested))
      sum(result$tested$FDR <= config$testing$fdr_threshold, na.rm = TRUE) else NA
  )
  finalize_run_provenance(prov, outdir, extra = extra)
  result$provenance <- prov
  result$config <- config
  result
}

.run_count_pipeline <- function(config, count_matrix, prov) {
  samples <- if (length(config$samples)) {
    validate_sample_metadata(config$samples)
  } else if (!is.null(count_matrix$samples)) {
    count_matrix$samples
  } else {
    stop("No samples in config and no count_matrix$samples.", call. = FALSE)
  }

  if (is.null(count_matrix)) {
    paths <- samples$bam_or_methylation_path
    names(paths) <- samples$sample_id
    .record_input_ids(prov, paths, samples$sample_id)
    genome <- config$preprocessing$chrom_sizes
    if (is.null(genome)) {
      stop("preprocessing.chrom_sizes is required when counting from BAMs.",
           call. = FALSE)
    }
    cm <- build_count_matrix(
      paths, genome,
      bin_size = config$preprocessing$bin_size,
      blacklist = config$preprocessing$blacklist,
      counting_mode = config$preprocessing$counting_mode,
      sample_ids = samples$sample_id,
      chromosomes = config$preprocessing$chromosomes,
      paired = config$preprocessing$paired,
      mapq_min = config$preprocessing$mapq_min,
      extra_exclude = config$preprocessing$extra_exclude,
      mask_low_coverage = config$preprocessing$mask_low_coverage,
      min_count_sum = config$preprocessing$min_count_sum
    )
  } else {
    cm <- count_matrix
    if (is.null(cm$bins) || is.null(cm$counts)) {
      stop("count_matrix must contain bins and counts.", call. = FALSE)
    }
  }

  counts <- cm$counts
  bins <- cm$bins
  if (!is.null(colnames(counts))) {
    counts <- counts[, samples$sample_id, drop = FALSE]
  }

  norm <- normalize_counts(
    counts,
    method = config$normalization$method,
    extreme_factor_ratio = config$normalization$extreme_factor_ratio,
    spikein_usage = config$normalization$spikein_usage,
    background_definition = config$normalization$background_definition
  )
  write_normalization_diagnostics(norm, file.path(config$output$directory,
                                                  "normalization_diagnostics.tsv"))

  effect <- estimate_differential_effect(
    counts, samples, offsets = norm$offsets,
    score = config$differential$score,
    control_level = config$differential$control_level,
    treatment_level = config$differential$treatment_level,
    prior_count = config$differential$prior_count,
    valid = bins$valid
  )
  track <- write_effect_track(bins, effect)
  utils::write.table(track, file.path(config$output$directory, "effect_track.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  if (isTRUE(config$output$write_bigwig)) {
    tryCatch(
      .write_bigwig_track(track, file.path(config$output$directory, "effect.bw")),
      error = function(e) warning("BigWig write failed: ", e$message, call. = FALSE)
    )
  }

  segs_all <- segment_track_penalty_path(
    track,
    penalties = config$segmentation$penalties,
    min_width = as.integer(config$segmentation$min_width),
    loss = config$segmentation$loss,
    engine = config$segmentation$engine,
    huber_k = config$segmentation$huber_k,
    winsor_k = config$segmentation$winsor_k
  )
  if (nrow(segs_all)) segs_all <- analyze_nesting(segs_all)
  cand <- select_candidate_segments(
    segs_all,
    strategy = config$segmentation$candidate_strategy,
    selected_penalty = config$segmentation$selected_penalty
  )

  if (nrow(cand) == 0L) {
    tested <- .empty_test_table()
    refined <- .empty_refined()
    events <- .empty_events()
  } else {
    recount <- .recount_from_count_matrix(cand, bins, counts)
    if (!is.null(count_matrix) && is.null(count_matrix$bam_files)) {
      # in-memory path already recounted from bins
    }
    tested <- test_regions_edger(
      recount, samples, offsets = norm$offsets,
      control_level = config$differential$control_level,
      treatment_level = config$differential$treatment_level,
      fdr_method = config$testing$fdr_method,
      candidate_source = cand$candidate_source[1],
      segmentation_penalty = config$segmentation$selected_penalty
    )
    tested$direction <- ifelse(tested$logFC > 0, "GAIN",
                               ifelse(tested$logFC < 0, "LOSS", "STABLE"))
    tested$level1 <- ifelse(
      !is.na(tested$FDR) & tested$FDR <= config$testing$fdr_threshold &
        abs(tested$logFC) >= config$testing$effect_threshold,
      tested$direction, "STABLE"
    )

    refined <- refine_boundaries(
      tested, track,
      method = config$boundary$refinement_method,
      search_window_bins = config$boundary$search_window_bins,
      finer_segments = segs_all
    )

    sig <- tested[!is.na(tested$FDR) &
                    tested$FDR <= config$testing$fdr_threshold, , drop = FALSE]
    if (isTRUE(config$topology$enabled %||% TRUE)) {
      ctrl_dom <- build_condition_domains(
        counts, bins, samples, config$differential$control_level,
        method = config$topology$domain_method,
        mean_quantile = config$topology$domain_mean_quantile,
        size_factors = norm$size_factors,
        min_width = config$segmentation$min_width
      )
      trt_dom <- build_condition_domains(
        counts, bins, samples, config$differential$treatment_level,
        method = config$topology$domain_method,
        mean_quantile = config$topology$domain_mean_quantile,
        size_factors = norm$size_factors,
        min_width = config$segmentation$min_width
      )
      matching <- match_domains(
        ctrl_dom, trt_dom,
        min_reciprocal_overlap = config$topology$min_reciprocal_overlap,
        max_match_distance = config$topology$max_match_distance,
        high_reciprocal_overlap = config$topology$high_reciprocal_overlap
      )
      events <- classify_boundary_events(
        matching, ctrl_dom, trt_dom,
        differential = sig,
        thresholds = config$topology
      )
    } else {
      events <- .empty_events()
    }
  }

  write_standard_outputs(
    config$output$directory,
    segments = if (nrow(cand)) cand else segs_all,
    tested = tested,
    events = events,
    write_bed = config$output$write_bed
  )

  list(
    samples = samples, bins = bins, counts = counts, qc = cm$qc,
    norm = norm, effect = effect, track = track,
    segments = segs_all, candidates = cand,
    tested = tested, refined = refined, events = events
  )
}

.run_methylation_pipeline <- function(config, count_matrix, prov) {
  samples <- validate_sample_metadata(config$samples)
  if (is.null(count_matrix)) {
    paths <- samples$bam_or_methylation_path
    names(paths) <- samples$sample_id
    .record_input_ids(prov, paths, samples$sample_id)
    tabs <- lapply(seq_along(paths), function(i) {
      x <- read_methylation_table(paths[[i]],
                                  min_coverage = config$methylation$min_coverage)
      x$sample <- samples$sample_id[i]
      x
    })
    meth <- do.call(rbind, tabs)
  } else {
    meth <- count_matrix$methylation
    if (is.null(meth)) stop("count_matrix$methylation is required.", call. = FALSE)
  }

  track <- estimate_methylation_effect(
    meth, samples,
    control_level = config$differential$control_level,
    treatment_level = config$differential$treatment_level,
    score = if (config$methylation$score %in% c("delta_beta", "standardized_effect"))
      config$methylation$score else "standardized_effect",
    min_coverage = config$methylation$min_coverage
  )
  utils::write.table(track, file.path(config$output$directory, "effect_track.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  segs_all <- segment_track_penalty_path(
    track,
    penalties = config$segmentation$penalties,
    min_width = as.integer(config$segmentation$min_width),
    loss = config$segmentation$loss,
    engine = config$segmentation$engine
  )
  cand <- select_candidate_segments(
    segs_all,
    strategy = config$segmentation$candidate_strategy,
    selected_penalty = config$segmentation$selected_penalty
  )
  # Region M/N reconstruction
  tested <- .empty_test_table()
  if (nrow(cand)) {
    mn <- .methylation_region_counts(cand, meth, samples)
    tested <- test_dmrs(
      mn, samples, method = if (config$testing$method == "edgeR_ql")
        "beta_binomial" else config$testing$method,
      control_level = config$differential$control_level,
      treatment_level = config$differential$treatment_level,
      fdr_method = config$testing$fdr_method
    )
    tested$chrom <- cand$chrom
    tested$start <- cand$start
    tested$end <- cand$end
  }
  write_standard_outputs(config$output$directory, segments = cand,
                         tested = tested, events = .empty_events(),
                         write_bed = config$output$write_bed)
  list(
    samples = samples, methylation = meth, track = track,
    segments = segs_all, candidates = cand, tested = tested,
    events = .empty_events()
  )
}

.methylation_region_counts <- function(regions, meth, samples) {
  out <- regions
  if (!"region_id" %in% names(out)) out$region_id <- seq_len(nrow(out))
  for (sid in samples$sample_id) {
    out[[paste0("meth_", sid)]] <- 0L
    out[[paste0("cov_", sid)]] <- 0L
  }
  if (!"sample" %in% names(meth) && "sample_id" %in% names(meth)) {
    meth$sample <- meth$sample_id
  }
  for (i in seq_len(nrow(regions))) {
    hit <- meth$chrom == regions$chrom[i] &
      meth$position >= regions$start[i] & meth$position < regions$end[i]
    if (!any(hit)) next
    sub <- meth[hit, ]
    for (sid in samples$sample_id) {
      ss <- sub[sub$sample == sid, ]
      out[i, paste0("meth_", sid)] <- sum(ss$methylated_count, na.rm = TRUE)
      out[i, paste0("cov_", sid)] <- sum(ss$total_count, na.rm = TRUE)
    }
  }
  out
}
