#' Empirical FDR from a null run
#'
#' @param tested data.frame with `FDR` or `pvalue` and interval columns.
#' @param nominal Vector of nominal q thresholds.
#' @param genome_width Total width of the analyzed genome (bp); used for
#'   bp-level FDR on a null (all discovered bp are false).
#' @return data.frame of empirical rates. On a true-null simulation every
#'   discovery is a false discovery.
#' @export
empirical_fdr <- function(tested, nominal = c(0.01, 0.05, 0.10),
                          genome_width = NA_real_) {
  if (!"FDR" %in% names(tested)) {
    stop("tested must contain an FDR column.", call. = FALSE)
  }
  has_iv <- all(c("start", "end") %in% names(tested))
  rows <- lapply(nominal, function(q) {
    hit <- !is.na(tested$FDR) & tested$FDR <= q
    n_hit <- sum(hit)
    bp <- if (has_iv && n_hit) sum(tested$end[hit] - tested$start[hit]) else 0
    data.frame(
      nominal_q = q,
      n_discoveries = n_hit,
      region_level_empirical_fdr = if (n_hit == 0L) 0 else 1,
      bp_discovered = bp,
      bp_level_empirical_fdr = if (bp == 0) 0 else 1,
      false_boundary_count = 2L * n_hit,
      genome_width = genome_width,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Permute condition labels
#'
#' @param samples Sample metadata.
#' @param seed Optional seed.
#' @return A copy of `samples` with shuffled `condition`.
#' @export
permute_condition_labels <- function(samples, seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  out <- samples
  out$condition <- sample(out$condition)
  if (length(unique(out$condition)) < 2L) {
    warning("Permutation produced a single condition; retrying with a new shuffle.",
            call. = FALSE)
    for (i in 1:20) {
      out$condition <- sample(samples$condition)
      if (length(unique(out$condition)) >= 2L) break
    }
  }
  out
}

#' Run the discovery+testing pipeline on a simulated cohort (in memory)
#'
#' Helper for FDR grids and unit tests. Does not read BAMs.
#'
#' @param sim Output of [simulate_count_assay()].
#' @param config Config list (defaults used for missing keys).
#' @return list with track, segments, tested, timing.
#' @keywords internal
run_simulated_discovery <- function(sim, config = default_config()) {
  t0 <- proc.time()[["elapsed"]]
  norm <- normalize_counts(sim$counts, method = config$normalization$method,
                           extreme_factor_ratio = config$normalization$extreme_factor_ratio)
  eff <- estimate_differential_effect(
    sim$counts, sim$samples, offsets = norm$offsets,
    score = config$differential$score,
    control_level = config$differential$control_level,
    treatment_level = config$differential$treatment_level,
    prior_count = config$differential$prior_count,
    valid = sim$bins$valid
  )
  track <- write_effect_track(sim$bins, eff)
  segs <- segment_track_penalty_path(
    track,
    penalties = config$segmentation$penalties,
    min_width = as.integer(config$segmentation$min_width),
    loss = config$segmentation$loss,
    engine = config$segmentation$engine,
    huber_k = config$segmentation$huber_k,
    winsor_k = config$segmentation$winsor_k
  )
  cand <- select_candidate_segments(
    segs,
    strategy = config$segmentation$candidate_strategy,
    selected_penalty = config$segmentation$selected_penalty
  )
  if (nrow(cand) == 0L) {
    tested <- .empty_test_table()
  } else {
    # Recount from the same binned counts by summing bins overlapping each region.
    recount <- .recount_from_count_matrix(cand, sim$bins, sim$counts)
    tested <- test_regions_edger(
      recount, sim$samples, offsets = norm$offsets,
      control_level = config$differential$control_level,
      treatment_level = config$differential$treatment_level,
      fdr_method = config$testing$fdr_method,
      candidate_source = cand$candidate_source[1],
      segmentation_penalty = config$segmentation$selected_penalty
    )
    if (!"start" %in% names(tested) && "start" %in% names(recount)) {
      tested$start <- recount$start
      tested$end <- recount$end
      tested$chrom <- recount$chrom
    }
  }
  list(
    norm = norm, effect = eff, track = track, segments = segs,
    candidates = cand, tested = tested,
    runtime_sec = proc.time()[["elapsed"]] - t0
  )
}

.recount_from_count_matrix <- function(regions, bins, counts) {
  regions <- .check_bed_frame(regions, "regions")
  if (!"region_id" %in% names(regions)) regions$region_id <- seq_len(nrow(regions))
  out <- regions
  for (j in seq_len(ncol(counts))) {
    nm <- colnames(counts)[j]
    out[[nm]] <- 0L
  }
  for (i in seq_len(nrow(regions))) {
    hit <- bins$chrom == regions$chrom[i] &
      bins$start < regions$end[i] & bins$end > regions$start[i]
    if (any(hit)) {
      out[i, colnames(counts)] <- as.list(as.integer(colSums(counts[hit, , drop = FALSE])))
    }
  }
  out
}
