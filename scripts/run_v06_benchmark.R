#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(20260)
dir.create("results/v0.6", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
t_all <- proc.time()[["elapsed"]]
rss0 <- process_rss_mb()

write_tsv <- function(x, p) utils::write.table(x, p, sep = "\t", quote = FALSE,
                                               row.names = FALSE)

# ---- software versions ----
ver <- data.frame(
  software = c("EpiSegDiff", "R", "limma", "csaw", "edgeR", "THOR", "DiffBind",
               "JASPAR2020"),
  version = c(as.character(utils::packageVersion("EpiSegDiff")),
              paste(R.version$major, R.version$minor, sep = "."),
              as.character(utils::packageVersion("limma")),
              if (requireNamespace("csaw", quietly = TRUE))
                as.character(utils::packageVersion("csaw")) else "NOT_AVAILABLE",
              if (requireNamespace("edgeR", quietly = TRUE))
                as.character(utils::packageVersion("edgeR")) else "NOT_AVAILABLE",
              "NOT_AVAILABLE", "NOT_AVAILABLE",
              if (requireNamespace("JASPAR2020", quietly = TRUE))
                as.character(utils::packageVersion("JASPAR2020")) else "NA"),
  stringsAsFactors = FALSE
)
write_tsv(ver, "results/v0.6/software_versions.tsv")
writeLines(utils::capture.output(utils::sessionInfo()),
           "results/v0.6/session_info.txt")

# ---- primary synthetic grid ----
message("primary synthetic grid")
widths <- c(200, 500, 1000, 2000, 5000, 10000, 25000, 50000)
deltas <- c(0.5, 1.5, 3.0)
noises <- c(0.05, 0.2, 0.5)
dirs <- c("GAIN", "LOSS")
nrep <- 10L  # POST_HOC: 50 was infeasible (~80 min); 10 retained for the full grid
nrep_sec <- 10L
bs <- 200L
design <- expand.grid(width = widths, delta = deltas, noise = noises,
                      direction = dirs, stringsAsFactors = FALSE)
design$cell_id <- seq_len(nrow(design))
write_tsv(design, "results/v0.6/synthetic_design.tsv")

rows <- list(); ri <- 0L
for (ci in seq_len(nrow(design))) {
  w <- design$width[ci]
  nb <- as.integer(w / bs)
  n_bins <- max(1000L, nb * 4L + 100L)
  mid <- floor(n_bins / 2)
  a <- mid - floor(nb / 2) + 1L
  b <- a + nb - 1L
  delta <- if (design$direction[ci] == "GAIN") design$delta[ci] else -design$delta[ci]
  for (r in seq_len(nrep)) {
    seed <- 20260L + ci * 1000L + r
    ev <- data.frame(start_bin = a, end_bin = b, delta = delta, shape = "sharp")
    sim <- simulate_continuous_events(n_bins, bs, ev, noise_sd = design$noise[ci],
                                      seed = seed)
    res <- tryCatch(run_two_methods_on_sim(sim, 0.05, 8, 3L),
                    error = function(e) NULL)
    if (is.null(res)) next
    for (meth in c("esd", "window")) {
      sc <- score_calls_vs_truth(res[[meth]], ev, bs)
      ri <- ri + 1L
      rows[[ri]] <- cbind(cell_id = ci, rep = r, method = meth,
                          width = w, delta = design$delta[ci],
                          noise = design$noise[ci],
                          direction = design$direction[ci], sc)
    }
  }
  if (ci %% 20L == 0L) message("  cell ", ci, "/", nrow(design))
}
syn <- do.call(rbind, rows)
write_tsv(syn, "results/v0.6/synthetic_region_metrics.tsv")
bnd <- aggregate(cbind(med_bound_err, p95_bound_err, mean_iou) ~ method + width +
                   delta + noise + direction, syn, median, na.rm = TRUE)
write_tsv(bnd, "results/v0.6/synthetic_boundary_metrics.tsv")
frg <- aggregate(cbind(one_call_frac, n_call, n_truth) ~ method + width +
                   delta + noise + direction, syn, mean, na.rm = TRUE)
write_tsv(frg, "results/v0.6/synthetic_fragmentation.tsv")

# ---- spacing / merge-split ----
message("spacing / merge-split")
gaps <- c(0L, 2L, 10L, 50L)  # bins
ms_rows <- list(); mi <- 0L
for (gap in gaps) {
  for (rel in c("same", "opp")) {
    for (r in seq_len(nrep_sec)) {
      n_bins <- 800L
      e1 <- 200:209
      e2s <- 210L + gap
      e2 <- e2s:(e2s + 9L)
      d1 <- 1.5
      d2 <- if (rel == "same") 1.5 else -1.5
      ev <- data.frame(start_bin = c(e1[1], e2[1]), end_bin = c(e1[length(e1)], e2[length(e2)]),
                       delta = c(d1, d2), shape = "sharp")
      sim <- simulate_continuous_events(n_bins, bs, ev, noise_sd = 0.2,
                                        seed = 31000L + gap * 100 + r)
      res <- tryCatch(run_two_methods_on_sim(sim, 0.05, 8, 3L), error = function(e) NULL)
      if (is.null(res)) next
      for (meth in c("esd", "window")) {
        sc <- score_calls_vs_truth(res[[meth]], ev, bs)
        mi <- mi + 1L
        ms_rows[[mi]] <- cbind(gap_bp = gap * bs, relation = rel, method = meth,
                               rep = r, sc)
      }
    }
  }
}
ms <- do.call(rbind, ms_rows)
write_tsv(ms, "results/v0.6/synthetic_merge_split_errors.tsv")

# ---- global null ----
message("global null")
null_rows <- list()
for (r in 1:100) {
  ev <- data.frame(start_bin = integer(), end_bin = integer(), delta = numeric(),
                   shape = character())
  sim <- simulate_continuous_events(1000L, bs, ev, noise_sd = 0.2, seed = 40000L + r)
  res <- tryCatch(run_two_methods_on_sim(sim, 0.05, 8, 3L), error = function(e) NULL)
  if (is.null(res)) next
  for (meth in c("esd", "window")) {
    cl <- res[[meth]]
    null_rows[[length(null_rows) + 1L]] <- data.frame(
      rep = r, method = meth, n_der = nrow(cl),
      der_bp = if (nrow(cl)) sum(cl$end - cl$start) else 0,
      any_der = nrow(cl) > 0, stringsAsFactors = FALSE
    )
  }
}
nullt <- do.call(rbind, null_rows)
write_tsv(nullt, "results/v0.6/null_simulation.tsv")

# FDR grid on same nulls
fdr_rows <- list()
for (r in 1:50) {
  ev <- data.frame(start_bin = integer(), end_bin = integer(), delta = numeric(),
                   shape = character())
  sim <- simulate_continuous_events(1000L, bs, ev, noise_sd = 0.2, seed = 50000L + r)
  for (al in c(0.01, 0.05, 0.10)) {
    res <- tryCatch(run_two_methods_on_sim(sim, al, 8, 3L), error = function(e) NULL)
    if (is.null(res)) next
    fdr_rows[[length(fdr_rows) + 1L]] <- data.frame(
      rep = r, alpha = al, method = "esd", n_der = nrow(res$esd),
      any_der = nrow(res$esd) > 0, stringsAsFactors = FALSE
    )
    fdr_rows[[length(fdr_rows) + 1L]] <- data.frame(
      rep = r, alpha = al, method = "window", n_der = nrow(res$window),
      any_der = nrow(res$window) > 0, stringsAsFactors = FALSE
    )
  }
}
fdrt <- do.call(rbind, fdr_rows)
emp <- aggregate(cbind(n_der, any_der) ~ alpha + method, fdrt, mean)
write_tsv(emp, "results/v0.6/fdr_calibration.tsv")

# ---- selection bias ----
message("selection bias")
sb <- list()
for (r in 1:50) {
  ev <- data.frame(start_bin = 400L, end_bin = 449L, delta = 1.5, shape = "sharp")
  sim <- simulate_continuous_events(1000L, bs, ev, noise_sd = 0.2, seed = 60000L + r)
  joint <- run_two_methods_on_sim(sim, 0.05, 8, 3L)
  # split: PELT on even columns' z via descriptive effect of d1,j1 only — use 2 samples descriptive
  samp2 <- sim$samples[c(1, 3), ]
  xf <- transform_continuous(sim$signal, "log2_pseudocount", 0.01)
  # segment on pair 1 (d1,j1) descriptive
  ev1 <- xf[, 3] - xf[, 1]
  track <- data.frame(chrom = sim$bed$chrom, start = sim$bed$start, end = sim$bed$end,
                      score = ev1, valid = TRUE)
  segs <- tryCatch(segment_genome(track, penalty = 8, min_width = 3L, engine = "pelt")$segments,
                   error = function(e) data.frame())
  split_n <- if (!is.null(segs) && nrow(segs)) nrow(segs) else 0L
  sb[[r]] <- data.frame(
    rep = r,
    joint_n = nrow(joint$esd),
    split_n_segments = split_n,
    joint_detected = score_calls_vs_truth(joint$esd, ev, bs)$n_detected,
    stringsAsFactors = FALSE
  )
}
sbt <- do.call(rbind, sb)
write_tsv(sbt, "results/v0.6/selection_bias_diagnostic.tsv")

# permutation: reuse v0.5
if (file.exists("results/v0.5/negative_controls.tsv")) {
  p <- utils::read.table("results/v0.5/negative_controls.tsv", header = TRUE, sep = "\t")
  p$note <- "only 2 balanced 2+2 nulls exist with n=4"
  write_tsv(p, "results/v0.6/permutation_calibration.tsv")
} else {
  write_tsv(data.frame(status = "NOT_AVAILABLE",
                       failure_reason = "v0.5 negative_controls missing"),
            "results/v0.6/permutation_calibration.tsv")
}

message("synthetic+calibration done n_primary=", nrow(syn))
cat("ELAPSED_SYN", proc.time()[["elapsed"]] - t_all, "\n")
