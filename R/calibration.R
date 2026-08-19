#' Wilson / binomial interval for a count
#'
#' @param x Successes.
#' @param n Trials.
#' @param level Confidence level.
#' @export
binomial_interval <- function(x, n, level = 0.95) {
  x <- as.integer(x); n <- as.integer(n)
  if (is.na(n) || n <= 0L) {
    return(data.frame(est = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  est <- x / n
  ci <- tryCatch(
    as.numeric(stats::binom.test(x, n, conf.level = level)$conf.int),
    error = function(e) c(NA_real_, NA_real_)
  )
  data.frame(est = est, lo = ci[1], hi = ci[2])
}

#' Sample-split discovery / testing
#'
#' Discover on half the replicates (balanced by condition), test on the
#' held-out half. Measurement experiment; not the default pipeline.
#'
#' @param sim Simulated cohort with ≥2 replicates per condition.
#' @param config Config list.
#' @export
sample_split_discovery <- function(sim, config = default_config()) {
  samples <- sim$samples
  conds <- unique(samples$condition)
  disc_ids <- unlist(lapply(conds, function(cn) {
    ids <- samples$sample_id[samples$condition == cn]
    if (length(ids) < 2L) {
      stop("sample_split_discovery requires ≥2 replicates per condition.",
           call. = FALSE)
    }
    ids[seq_len(max(1L, floor(length(ids) / 2)))]
  }))
  test_ids <- setdiff(samples$sample_id, disc_ids)
  if (!length(test_ids) || !length(disc_ids)) {
    stop("Could not split samples.", call. = FALSE)
  }

  sim_d <- sim
  sim_d$counts <- sim$counts[, disc_ids, drop = FALSE]
  sim_d$samples <- samples[match(disc_ids, samples$sample_id), , drop = FALSE]
  sim_t <- sim
  sim_t$counts <- sim$counts[, test_ids, drop = FALSE]
  sim_t$samples <- samples[match(test_ids, samples$sample_id), , drop = FALSE]

  fit_d <- run_simulated_discovery(sim_d, config)
  cand <- fit_d$candidates
  if (is.null(cand) || !nrow(cand)) {
    return(list(discovery = fit_d, tested = .empty_test_table(),
                n_discover = length(disc_ids), n_test = length(test_ids)))
  }
  rec <- .recount_from_count_matrix(cand, sim_t$bins, sim_t$counts)
  norm_t <- normalize_counts(sim_t$counts, method = config$normalization$method,
                             extreme_factor_ratio = config$normalization$extreme_factor_ratio)
  tested <- test_regions_edger(
    rec, sim_t$samples, offsets = norm_t$offsets,
    control_level = config$differential$control_level,
    treatment_level = config$differential$treatment_level,
    fdr_method = config$testing$fdr_method,
    candidate_source = "sample_split",
    testing_method = "edgeR_ql_sample_split"
  )
  if ("start" %in% names(rec) && !"start" %in% names(tested)) {
    tested$chrom <- rec$chrom
    tested$start <- rec$start
    tested$end <- rec$end
  }
  list(discovery = fit_d, tested = tested,
       n_discover = length(disc_ids), n_test = length(test_ids))
}

#' Permutation-calibrated region p-values
#'
#' Discovers on the original labels, then compares each region QL
#' statistic to the same statistic after permuting condition labels
#' `B` times. Measurement experiment.
#'
#' @param sim Simulated cohort.
#' @param config Config.
#' @param B Number of permutations.
#' @export
permutation_calibrate_regions <- function(sim, config = default_config(), B = 20L) {
  fit <- run_simulated_discovery(sim, config)
  tested <- fit$tested
  if (is.null(tested) || !nrow(tested) || !"statistic" %in% names(tested)) {
    return(list(fit = fit, tested = tested, B = B))
  }
  stat0 <- tested$statistic
  null_ge <- integer(length(stat0))
  for (b in seq_len(as.integer(B))) {
    perm <- sim
    perm$samples <- permute_condition_labels(sim$samples, seed = 1000L + b)
    if (length(unique(perm$samples$condition)) < 2L) next
    fb <- tryCatch(run_simulated_discovery(perm, config), error = function(e) NULL)
    if (is.null(fb) || is.null(fb$tested) || !nrow(fb$tested)) next
    # compare original stats to the permutation's region stats (pooled)
    null_ge <- null_ge + as.integer(stat0 <= max(fb$tested$statistic, na.rm = TRUE))
  }
  tested$perm_pvalue <- (1 + null_ge) / (1 + B)
  list(fit = fit, tested = tested, B = B)
}
