#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("results/v0.7", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
write_tsv <- function(x, p) utils::write.table(x, p, sep = "\t", quote = FALSE,
                                               row.names = FALSE)
nrep <- 50L
bs <- 200L
scores <- c("raw_effect", "log_effect", "standardized", "hotelling")
segs <- c("pelt", "binseg", "wbs", "fused_lasso", "hmm", "hsmm")
widths <- c(200, 500, 2000, 10000, 50000)

eval_one <- function(sim, score_e, seg_e, min_width, penalty = 8,
                     scale_mode = "single", inf = "same_data") {
  # Synthetic detection uses |segment mean of score| (frozen 0.35*MAD+0.15),
  # not per-eval limma — limma is reserved for inference-architecture cells.
  calls <- tryCatch({
    sc <- compute_diff_score(sim$signal, sim$samples, "JNKi", "DMSO", score_e)
    if (identical(inf, "pooled")) sc <- as.numeric(rowMeans(sim$signal))
    if (identical(scale_mode, "coarse_fine")) {
      seg <- hierarchical_refine(sc, sim$bed, 5L, seg_e, penalty, min_width)
    } else {
      seg <- segment_signal(sc, positions = sim$bed$start, engine = seg_e,
                            penalty = penalty, min_width = min_width)
    }
    thr <- 0.15 + 0.35 * stats::mad(sc, constant = 1)
    if (!is.finite(thr) || thr <= 0) thr <- 0.2
    .segments_to_calls(seg, sim$bed, thresh = thr)
  }, error = function(e) NULL)
  if (is.null(calls)) return(NULL)
  sc <- score_calls_vs_truth(calls, sim$truth, sim$bin_size)
  cbind(score_engine = score_e, segmentation_engine = seg_e,
        min_width = min_width, inference = inf, scale_mode = scale_mode,
        n_call = nrow(calls), sc)
}

# ---- score engines (PELT min1 + v0.6 PELT min3) ----
message("score engines")
sr <- list(); si <- 0L
for (w in widths) {
  nb <- max(1L, as.integer(w / bs))
  n_bins <- max(400L, nb * 4L + 80L)
  a <- floor(n_bins / 2) - floor(nb / 2)
  ev0 <- data.frame(start_bin = a, end_bin = a + nb - 1L, delta = 1.5, shape = "sharp")
  for (sc in scores) {
    for (mw in c(1L, 3L)) {
      for (r in seq_len(nrep)) {
        sim <- simulate_continuous_events(n_bins, bs, ev0, noise_sd = 0.2,
                                          seed = 70000L + w + mw * 17 + r)
        row <- eval_one(sim, sc, "pelt", mw)
        if (is.null(row)) next
        si <- si + 1L
        sr[[si]] <- cbind(width = w, rep = r, row)
      }
    }
  }
}
score_tab <- do.call(rbind, sr)
write_tsv(score_tab, "results/v0.7/score_engine_benchmark.tsv")
message("  score rows ", nrow(score_tab))

# ---- segmentation engines on raw_effect ----
message("segmentation engines")
er <- list(); ei <- 0L
for (w in widths) {
  nb <- max(1L, as.integer(w / bs))
  n_bins <- max(400L, nb * 4L + 80L)
  a <- floor(n_bins / 2) - floor(nb / 2)
  ev0 <- data.frame(start_bin = a, end_bin = a + nb - 1L, delta = 1.5, shape = "sharp")
  for (se in segs) {
    mw <- if (se %in% c("fused_lasso", "wbs", "binseg")) 1L else 3L
    for (r in seq_len(nrep)) {
      sim <- simulate_continuous_events(n_bins, bs, ev0, noise_sd = 0.2,
                                        seed = 71000L + w + match(se, segs) * 13 + r)
      row <- eval_one(sim, "raw_effect", se, mw)
      if (is.null(row)) next
      ei <- ei + 1L
      er[[ei]] <- cbind(width = w, rep = r, row)
    }
  }
}
seg_tab <- do.call(rbind, er)
write_tsv(seg_tab, "results/v0.7/segmentation_engine_benchmark.tsv")
message("  seg rows ", nrow(seg_tab))

# ---- nested + gradual + merge ----
message("nested / merge / gradual")
nr <- list()
for (r in seq_len(nrep)) {
  # broad 100 bins + narrow 2 bins inside
  ev <- data.frame(start_bin = c(80L, 120L), end_bin = c(180L, 121L),
                   delta = c(1.2, 3.0), shape = c("sharp", "sharp"))
  sim <- simulate_continuous_events(260L, bs, ev, noise_sd = 0.15, seed = 72000L + r)
  for (se in c("pelt", "wbs", "fused_lasso")) {
    calls <- tryCatch(run_multiscale(sim, "raw_effect", se, "single",
                                     if (se == "pelt") 3L else 1L, 8),
                      error = function(e) NULL)
    if (is.null(calls)) next
    sc <- score_calls_vs_truth(calls, ev, bs)
    nr[[length(nr) + 1L]] <- cbind(rep = r, engine = se, geometry = "nested", sc)
  }
  # union multi-scale
  c1 <- tryCatch(run_multiscale(sim, "raw_effect", "wbs", "single", 1L, 8),
                 error = function(e) NULL)
  c2 <- tryCatch(run_multiscale(sim, "raw_effect", "pelt", "coarse_fine", 1L, 8),
                 error = function(e) NULL)
  u <- union_resolve_calls(list(narrow = c1, coarse = c2))
  scu <- score_calls_vs_truth(u, ev, bs)
  nr[[length(nr) + 1L]] <- cbind(rep = r, engine = "union_wbs_pelt",
                                 geometry = "nested", scu)
}
write_tsv(do.call(rbind, nr), "results/v0.7/nested_event_metrics.tsv")

# merge/split same-dir gap
mr <- list()
for (gap in c(0L, 2L, 10L)) {
  for (r in seq_len(nrep)) {
    ev <- data.frame(start_bin = c(100L, 110L + gap),
                     end_bin = c(109L, 119L + gap),
                     delta = c(1.5, 1.5), shape = "sharp")
    sim <- simulate_continuous_events(250L, bs, ev, noise_sd = 0.2,
                                      seed = 73000L + gap * 50 + r)
    for (se in c("pelt", "wbs", "fused_lasso")) {
      row <- eval_one(sim, "raw_effect", se, if (se == "pelt") 3L else 1L)
      if (is.null(row)) next
      mr[[length(mr) + 1L]] <- cbind(gap_bins = gap, rep = r, row)
    }
  }
}
write_tsv(do.call(rbind, mr), "results/v0.7/merge_split_metrics.tsv")

# gradual vs sharp at 5 kb
gr <- list()
for (shp in c("sharp", "gradual")) {
  ev <- data.frame(start_bin = 80L, end_bin = 104L, delta = 1.5, shape = shp)
  for (r in seq_len(nrep)) {
    sim <- simulate_continuous_events(200L, bs, ev, noise_sd = 0.2,
                                      seed = 74000L + nchar(shp) + r)
    for (se in c("pelt", "wbs", "fused_lasso")) {
      row <- eval_one(sim, "raw_effect", se, 1L)
      if (is.null(row)) next
      gr[[length(gr) + 1L]] <- cbind(shape = shp, rep = r, row)
    }
  }
}
write_tsv(do.call(rbind, gr), "results/v0.7/boundary_refinement.tsv")

# ---- calibration + inference ----
message("calibration / inference")
cal <- list()
for (inf in c("same_data", "pooled", "split")) {
  for (se in c("pelt", "wbs", "fused_lasso")) {
    for (r in seq_len(nrep)) {
      ev <- data.frame(start_bin = integer(), end_bin = integer(),
                       delta = numeric(), shape = character())
      sim <- simulate_continuous_events(400L, bs, ev, noise_sd = 0.2,
                                        seed = 75000L + match(inf, c("same_data", "pooled", "split")) * 80 + r)
      calls <- tryCatch(run_multiscale(sim, "raw_effect", se, "single",
                                       if (se == "pelt") 3L else 1L, 8, inf),
                        error = function(e) NULL)
      n <- if (is.null(calls)) NA_integer_ else nrow(calls)
      cal[[length(cal) + 1L]] <- data.frame(
        inference = inf, engine = se, rep = r, n_der = n,
        any_der = isTRUE(n > 0), stringsAsFactors = FALSE
      )
    }
  }
}
# crossfit null
for (r in seq_len(nrep)) {
  ev <- data.frame(start_bin = integer(), end_bin = integer(),
                   delta = numeric(), shape = character())
  sim <- simulate_continuous_events(400L, bs, ev, noise_sd = 0.2, seed = 76000L + r)
  calls <- tryCatch(run_crossfit(sim, "raw_effect", "wbs", 8, 1L),
                    error = function(e) NULL)
  n <- if (is.null(calls)) 0L else nrow(calls)
  cal[[length(cal) + 1L]] <- data.frame(
    inference = "crossfit", engine = "wbs", rep = r, n_der = n,
    any_der = n > 0, stringsAsFactors = FALSE
  )
}
calt <- do.call(rbind, cal)
write_tsv(calt, "results/v0.7/calibration_metrics.tsv")

# power vs FDP: same cells with planted vs null already
pw <- list()
for (inf in c("same_data", "pooled", "split")) {
  recs <- numeric()
  for (r in seq_len(nrep)) {
    ev <- data.frame(start_bin = 150L, end_bin = 159L, delta = 1.5, shape = "sharp")
    sim <- simulate_continuous_events(320L, bs, ev, noise_sd = 0.2, seed = 77000L + r)
    calls <- tryCatch(run_multiscale(sim, "raw_effect", "wbs", "single", 1L, 8, inf),
                      error = function(e) NULL)
    if (is.null(calls)) next
    recs <- c(recs, score_calls_vs_truth(calls, ev, bs)$recall)
  }
  fdp <- mean(calt$any_der[calt$inference == inf & calt$engine == "wbs"])
  pw[[length(pw) + 1L]] <- data.frame(inference = inf, engine = "wbs",
                                      mean_recall = mean(recs, na.rm = TRUE),
                                      p_any_false = fdp, stringsAsFactors = FALSE)
}
write_tsv(do.call(rbind, pw), "results/v0.7/power_calibration_frontier.tsv")

# oracle = best F1 across engines per width (from seg_tab)
or <- aggregate(f1 ~ width + segmentation_engine, seg_tab, mean, na.rm = TRUE)
best <- do.call(rbind, lapply(split(or, or$width), function(d) d[which.max(d$f1), ]))
best$role <- "oracle_per_width"
write_tsv(best, "results/v0.7/oracle_comparison.tsv")

# narrow / broad slices
write_tsv(seg_tab[seg_tab$width %in% c(200, 500), ], "results/v0.7/narrow_event_metrics.tsv")
write_tsv(seg_tab[seg_tab$width %in% c(10000, 50000), ], "results/v0.7/broad_event_metrics.tsv")

# selection table
sel <- aggregate(cbind(f1, med_bound_err, n_merged, n_split) ~
                   segmentation_engine + width, seg_tab, mean, na.rm = TRUE)
write_tsv(sel, "results/v0.7/engine_selection.tsv")
write_tsv(data.frame(step = "synthetic_v07", nrep = nrep, status = "PASS"),
          "results/v0.7/runtime_memory.tsv")
message("v0.7 synthetic written")
