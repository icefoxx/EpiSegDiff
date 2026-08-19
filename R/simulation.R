#' Simulate a two-condition count-assay cohort
#'
#' Generates NB counts on a single chromosome with optional planted
#' ground-truth events. Used for unit tests, FDR calibration, and
#' benchmarks.
#'
#' @param n_bins Number of bins.
#' @param bin_size Bin width (bp).
#' @param chrom Chromosome name.
#' @param n_rep Replicates per condition.
#' @param depth Mean background count per bin (before size factors).
#' @param dispersion NB dispersion \(\phi\) (`variance = mu + phi * mu^2`).
#' @param background_sd Extra log-normal background heterogeneity.
#' @param size_factors Optional sample size factors (length `2 * n_rep`).
#' @param events data.frame of planted events (see Details) or NULL.
#' @param mask_bins Optional integer indices to mark invalid.
#' @param seed Optional seed.
#' @return A list with `counts`, `bins`, `samples`, `truth`, `params`.
#'
#' @details
#' `events` columns: `event_class`, `start_bin`, `end_bin` (1-based,
#' inclusive), optional `effect` (logFC), optional geometry fields
#' `control_start_bin`, `control_end_bin`, `treatment_start_bin`,
#' `treatment_end_bin`.
#'
#' @export
simulate_count_assay <- function(n_bins = 400L, bin_size = 200L, chrom = "chr1",
                                 n_rep = 3L, depth = 20, dispersion = 0.1,
                                 background_sd = 0.2, size_factors = NULL,
                                 events = NULL, mask_bins = NULL, seed = NULL,
                                 drift_amp = 0, high_bg_blocks = NULL,
                                 coverage_scale = NULL, spikein_n_bins = 0L) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  n_bins <- as.integer(n_bins)
  n_rep <- as.integer(n_rep)
  n_s <- 2L * n_rep
  sample_ids <- c(paste0("ctrl_r", seq_len(n_rep)), paste0("trt_r", seq_len(n_rep)))
  samples <- data.frame(
    sample_id = sample_ids,
    condition = c(rep("control", n_rep), rep("treatment", n_rep)),
    replicate = c(seq_len(n_rep), seq_len(n_rep)),
    assay = "cutnrun",
    target = "sim",
    bam_or_methylation_path = NA_character_,
    control_path = NA,
    spikein_status = "none",
    library_metadata = "simulated",
    stringsAsFactors = FALSE
  )
  if (is.null(size_factors)) size_factors <- rep(1, n_s)
  if (length(size_factors) != n_s) {
    stop("size_factors must have length 2 * n_rep.", call. = FALSE)
  }

  alpha <- log(depth) + stats::rnorm(n_bins, 0, background_sd)
  if (is.finite(drift_amp) && drift_amp != 0) {
    x <- seq_len(n_bins) / n_bins
    alpha <- alpha + drift_amp * sin(2 * pi * x)
  }
  if (!is.null(high_bg_blocks) && nrow(as.data.frame(high_bg_blocks))) {
    hb <- as.data.frame(high_bg_blocks)
    for (i in seq_len(nrow(hb))) {
      a <- max(1L, as.integer(hb$start_bin[i]))
      b <- min(n_bins, as.integer(hb$end_bin[i]))
      boost <- if ("log_boost" %in% names(hb)) hb$log_boost[i] else 1.2
      alpha[a:b] <- alpha[a:b] + boost
    }
  }
  beta <- numeric(n_bins)
  spikein_idx <- integer()
  if (as.integer(spikein_n_bins) > 0L) {
    spikein_idx <- seq.int(n_bins - as.integer(spikein_n_bins) + 1L, n_bins)
    alpha[spikein_idx] <- log(depth) + 1.5
    beta[spikein_idx] <- 0
  }
  truth <- .empty_truth()

  if (!is.null(events) && nrow(events) > 0L) {
    for (i in seq_len(nrow(events))) {
      ev <- events[i, , drop = FALSE]
      planted <- .plant_event(beta, ev, n_bins)
      beta <- planted$beta
      truth <- rbind(truth, planted$truth)
    }
  }

  mu <- matrix(NA_real_, n_bins, n_s)
  X <- c(rep(0, n_rep), rep(1, n_rep))
  if (is.null(coverage_scale)) coverage_scale <- rep(1, n_bins)
  if (length(coverage_scale) != n_bins) {
    stop("coverage_scale must have length n_bins.", call. = FALSE)
  }
  for (j in seq_len(n_s)) {
    mu[, j] <- exp(log(size_factors[j]) + alpha + beta * X[j]) * coverage_scale
  }
  counts <- matrix(0L, n_bins, n_s, dimnames = list(NULL, sample_ids))
  for (j in seq_len(n_s)) {
    counts[, j] <- stats::rnbinom(n_bins, mu = mu[, j], size = 1 / max(dispersion, 1e-8))
  }

  starts <- as.integer((seq_len(n_bins) - 1L) * bin_size)
  bins <- data.frame(
    chrom = chrom,
    start = starts,
    end = starts + as.integer(bin_size),
    bin_id = seq_len(n_bins),
    valid = TRUE,
    stringsAsFactors = FALSE
  )
  if (!is.null(mask_bins)) {
    mask_bins <- as.integer(mask_bins)
    bins$valid[mask_bins] <- FALSE
    counts[mask_bins, ] <- 0L
  }

  if (nrow(truth) > 0L) {
    truth$chrom <- chrom
    truth$start <- bins$start[truth$start_bin]
    truth$end <- bins$end[truth$end_bin]
    if (!"control_start_bin" %in% names(truth) || all(is.na(truth$control_start_bin))) {
      truth$control_start <- truth$start
      truth$control_end <- truth$end
    } else {
      truth$control_start <- ifelse(is.na(truth$control_start_bin), NA,
                                    bins$start[pmax(1L, pmin(n_bins, truth$control_start_bin))])
      truth$control_end <- ifelse(is.na(truth$control_end_bin), NA,
                                  bins$end[pmax(1L, pmin(n_bins, truth$control_end_bin))])
    }
    if (!"treatment_start_bin" %in% names(truth) || all(is.na(truth$treatment_start_bin))) {
      truth$treatment_start <- truth$start
      truth$treatment_end <- truth$end
    } else {
      truth$treatment_start <- ifelse(is.na(truth$treatment_start_bin), NA,
                                      bins$start[pmax(1L, pmin(n_bins, truth$treatment_start_bin))])
      truth$treatment_end <- ifelse(is.na(truth$treatment_end_bin), NA,
                                    bins$end[pmax(1L, pmin(n_bins, truth$treatment_end_bin))])
    }
  }

  list(
    counts = counts,
    bins = bins,
    samples = samples,
    truth = truth,
    beta = beta,
    alpha = alpha,
    size_factors = as.numeric(size_factors),
    spikein_idx = spikein_idx,
    background_idx = setdiff(which(abs(beta) < 1e-8), spikein_idx),
    event_idx = which(abs(beta) >= 1e-8),
    params = list(
      n_bins = n_bins, bin_size = bin_size, chrom = chrom, n_rep = n_rep,
      depth = depth, dispersion = dispersion, background_sd = background_sd,
      seed = seed %||% NA_integer_,
      drift_amp = drift_amp
    )
  )
}

.empty_truth <- function() {
  data.frame(
    event_class = character(),
    start_bin = integer(), end_bin = integer(),
    effect_size = numeric(),
    control_start_bin = integer(), control_end_bin = integer(),
    treatment_start_bin = integer(), treatment_end_bin = integer(),
    stringsAsFactors = FALSE
  )
}

.plant_event <- function(beta, ev, n_bins) {
  cls <- as.character(ev$event_class[1])
  s <- as.integer(ev$start_bin[1])
  e <- as.integer(ev$end_bin[1])
  fx <- if ("effect" %in% names(ev) && is.finite(ev$effect[1])) ev$effect[1] else 1.5
  cs <- if ("control_start_bin" %in% names(ev)) ev$control_start_bin[1] else s
  ce <- if ("control_end_bin" %in% names(ev)) ev$control_end_bin[1] else e
  ts <- if ("treatment_start_bin" %in% names(ev)) ev$treatment_start_bin[1] else s
  te <- if ("treatment_end_bin" %in% names(ev)) ev$treatment_end_bin[1] else e

  clamp <- function(x) as.integer(pmax(1L, pmin(n_bins, x)))
  s <- clamp(s); e <- clamp(e); if (e < s) e <- s
  cs <- clamp(cs); ce <- clamp(ce); if (ce < cs) ce <- cs
  ts <- clamp(ts); te <- clamp(te); if (te < ts) te <- ts

  # Baseline support is encoded as extra signal via a large positive alpha
  # locally by adding to beta on treatment-only or control-only intervals.
  if (cls == "GAIN") {
    beta[s:e] <- beta[s:e] + fx
  } else if (cls == "LOSS") {
    beta[s:e] <- beta[s:e] - abs(fx)
  } else if (cls == "EXPANSION") {
    beta[ts:te] <- beta[ts:te] + fx
    # control-only core already present as baseline: reduce beta on core so
    # shared core is closer to zero effect and flanks are gain.
    if (cs >= ts && ce <= te) {
      beta[cs:ce] <- beta[cs:ce] - fx * 0.15
    }
  } else if (cls == "CONTRACTION") {
    beta[cs:ce] <- beta[cs:ce] - abs(fx) * 0.2
    # treatment-only loss on flanks
    if (ts > cs) beta[cs:(ts - 1L)] <- beta[cs:(ts - 1L)] - abs(fx)
    if (te < ce) beta[(te + 1L):ce] <- beta[(te + 1L):ce] - abs(fx)
  } else if (cls == "LEFT_SHIFT") {
    beta[ts:te] <- beta[ts:te] + fx
    beta[cs:ce] <- beta[cs:ce] - fx * 0.3
  } else if (cls == "RIGHT_SHIFT") {
    beta[ts:te] <- beta[ts:te] + fx
    beta[cs:ce] <- beta[cs:ce] - fx * 0.3
  } else if (cls == "SPLIT") {
    # loss in the internal gap, gain on the two products
    beta[cs:ce] <- beta[cs:ce] - abs(fx) * 0.2
    g1e <- if ("gap_start_bin" %in% names(ev)) ev$gap_start_bin[1] else
      cs + (ce - cs) %/% 3L
    g2s <- if ("gap_end_bin" %in% names(ev)) ev$gap_end_bin[1] else
      ce - (ce - cs) %/% 3L
    g1e <- clamp(g1e); g2s <- clamp(g2s)
    if (g2s > g1e) beta[g1e:g2s] <- beta[g1e:g2s] - abs(fx)
  } else if (cls == "MERGE") {
    beta[ts:te] <- beta[ts:te] + fx
  } else if (cls == "INTERNAL_GAIN") {
    beta[s:e] <- beta[s:e] + fx
  } else if (cls == "INTERNAL_LOSS") {
    beta[s:e] <- beta[s:e] - abs(fx)
  } else if (cls %in% c("NO_CHANGE", "NULL")) {
    # no beta change; "NULL" is accepted only as a deprecated alias
  } else {
    stop("Unknown event_class: ", cls, call. = FALSE)
  }

  truth <- data.frame(
    event_class = cls,
    start_bin = s, end_bin = e,
    effect_size = fx,
    control_start_bin = cs, control_end_bin = ce,
    treatment_start_bin = ts, treatment_end_bin = te,
    stringsAsFactors = FALSE
  )
  list(beta = beta, truth = truth)
}

#' Simulate a null (no biological difference) count cohort
#'
#' @inheritParams simulate_count_assay
#' @export
simulate_null_counts <- function(n_bins = 400L, bin_size = 200L, chrom = "chr1",
                                 n_rep = 3L, depth = 20, dispersion = 0.1,
                                 background_sd = 0.2, size_factors = NULL,
                                 mask_bins = NULL, seed = NULL) {
  simulate_count_assay(
    n_bins = n_bins, bin_size = bin_size, chrom = chrom, n_rep = n_rep,
    depth = depth, dispersion = dispersion, background_sd = background_sd,
    size_factors = size_factors, mask_bins = mask_bins,
    events = data.frame(event_class = "NO_CHANGE", start_bin = 1L, end_bin = 1L,
                        effect = 0, stringsAsFactors = FALSE),
    seed = seed
  )
}

#' Build a standard event-class example set for topology tests
#'
#' Returns event tables on a 400-bin canvas, one class at a time, placed
#' far from chromosome ends.
#'
#' @param event_class One of the locked topology classes.
#' @param offset Starting bin for the motif.
#' @export
example_topology_event <- function(event_class, offset = 50L) {
  o <- as.integer(offset)
  if (identical(event_class, "NULL")) event_class <- "NO_CHANGE"
  switch(
    event_class,
    GAIN = data.frame(event_class = "GAIN", start_bin = o, end_bin = o + 20L,
                      effect = 1.6, stringsAsFactors = FALSE),
    LOSS = data.frame(event_class = "LOSS", start_bin = o, end_bin = o + 20L,
                      effect = 1.6, stringsAsFactors = FALSE),
    EXPANSION = data.frame(
      event_class = "EXPANSION", start_bin = o, end_bin = o + 40L, effect = 1.4,
      control_start_bin = o + 10L, control_end_bin = o + 30L,
      treatment_start_bin = o, treatment_end_bin = o + 40L,
      stringsAsFactors = FALSE
    ),
    CONTRACTION = data.frame(
      event_class = "CONTRACTION", start_bin = o, end_bin = o + 40L, effect = 1.4,
      control_start_bin = o, control_end_bin = o + 40L,
      treatment_start_bin = o + 10L, treatment_end_bin = o + 30L,
      stringsAsFactors = FALSE
    ),
    LEFT_SHIFT = data.frame(
      event_class = "LEFT_SHIFT", start_bin = o, end_bin = o + 30L, effect = 1.4,
      control_start_bin = o + 10L, control_end_bin = o + 30L,
      treatment_start_bin = o, treatment_end_bin = o + 20L,
      stringsAsFactors = FALSE
    ),
    RIGHT_SHIFT = data.frame(
      event_class = "RIGHT_SHIFT", start_bin = o, end_bin = o + 30L, effect = 1.4,
      control_start_bin = o, control_end_bin = o + 20L,
      treatment_start_bin = o + 10L, treatment_end_bin = o + 30L,
      stringsAsFactors = FALSE
    ),
    SPLIT = data.frame(
      event_class = "SPLIT", start_bin = o, end_bin = o + 50L, effect = 1.4,
      control_start_bin = o, control_end_bin = o + 50L,
      treatment_start_bin = o, treatment_end_bin = o + 50L,
      gap_start_bin = o + 18L, gap_end_bin = o + 32L,
      stringsAsFactors = FALSE
    ),
    MERGE = data.frame(
      event_class = "MERGE", start_bin = o, end_bin = o + 50L, effect = 1.4,
      control_start_bin = o, control_end_bin = o + 50L,
      treatment_start_bin = o, treatment_end_bin = o + 50L,
      stringsAsFactors = FALSE
    ),
    INTERNAL_GAIN = data.frame(
      event_class = "INTERNAL_GAIN", start_bin = o + 10L, end_bin = o + 20L,
      effect = 1.5,
      control_start_bin = o, control_end_bin = o + 40L,
      treatment_start_bin = o, treatment_end_bin = o + 40L,
      stringsAsFactors = FALSE
    ),
    INTERNAL_LOSS = data.frame(
      event_class = "INTERNAL_LOSS", start_bin = o + 10L, end_bin = o + 20L,
      effect = 1.5,
      control_start_bin = o, control_end_bin = o + 40L,
      treatment_start_bin = o, treatment_end_bin = o + 40L,
      stringsAsFactors = FALSE
    ),
    NO_CHANGE = data.frame(event_class = "NO_CHANGE", start_bin = o, end_bin = o,
                           effect = 0, stringsAsFactors = FALSE),
    stop("Unknown event_class: ", event_class, call. = FALSE)
  )
}
