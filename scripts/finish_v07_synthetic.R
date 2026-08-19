#!/usr/bin/env Rscript
# Remaining v0.7 synthetic tables, honest selection, extra nulls, figures.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("results/v0.7", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
write_tsv <- function(x, p) utils::write.table(x, p, sep = "\t", quote = FALSE,
                                               row.names = FALSE)
readt <- function(p) utils::read.table(p, header = TRUE, sep = "\t",
                                       stringsAsFactors = FALSE)
nrep <- 50L
bs <- 200L
t0 <- proc.time()[["elapsed"]]

eval_one <- function(sim, score_e, seg_e, min_width, penalty = 8,
                     scale_mode = "single", inf = "same_data") {
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

# ---- honest engine selection + oracle (zeros included) ----
message("honest selection / oracle")
sc <- readt("results/v0.7/score_engine_benchmark.tsv")
sg <- readt("results/v0.7/segmentation_engine_benchmark.tsv")
mean0 <- function(x) mean(x, na.rm = TRUE)
sel_seg <- aggregate(cbind(f1, recall, precision, n_call, n_merged, n_split) ~
                       segmentation_engine + width, sg, mean0)
b_err <- aggregate(med_bound_err ~ segmentation_engine + width, sg, mean0)
sel_seg <- merge(sel_seg, b_err, by = c("segmentation_engine", "width"), all.x = TRUE)
sel_sc <- aggregate(cbind(f1, recall, precision, n_call) ~
                      score_engine + width + min_width, sc, mean0)
write_tsv(sel_seg, "results/v0.7/engine_selection.tsv")

# oracle = max mean F1 at each width across score-engine PELT cells AND seg engines
or_rows <- list()
for (w in sort(unique(c(sc$width, sg$width)))) {
  a <- sel_sc[sel_sc$width == w, ]
  b <- sel_seg[sel_seg$width == w, ]
  cand <- rbind(
    data.frame(width = w, family = "score+pelt",
               config = paste(a$score_engine, "pelt_min", a$min_width, sep = "_"),
               f1 = a$f1, stringsAsFactors = FALSE),
    data.frame(width = w, family = "raw+seg",
               config = paste("raw_effect", b$segmentation_engine, sep = "_"),
               f1 = b$f1, stringsAsFactors = FALSE)
  )
  cand <- cand[is.finite(cand$f1), ]
  or_rows[[length(or_rows) + 1L]] <- cand[which.max(cand$f1), ]
}
oracle <- do.call(rbind, or_rows)
# production candidate C: standardized+PELT min1 if width<=500 else raw+PELT min3
# (width known only in simulation; used as oracle-gap reference, not a deployable rule)
prod_f1 <- sapply(oracle$width, function(w) {
  if (w <= 500) {
    z <- sel_sc$f1[sel_sc$width == w & sel_sc$score_engine == "standardized" &
                     sel_sc$min_width == 1]
    if (length(z)) z[1] else NA
  } else {
    z <- sel_seg$f1[sel_seg$width == w & sel_seg$segmentation_engine == "pelt"]
    if (length(z)) z[1] else NA
  }
})
oracle$production_C_f1 <- as.numeric(prod_f1)
oracle$oracle_gap <- oracle$f1 - oracle$production_C_f1
write_tsv(oracle, "results/v0.7/oracle_comparison.tsv")

# ---- multi-scale strategies (nrep=50 primary widths) ----
message("multi-scale")
ms_rows <- list()
for (w in c(200, 500, 2000, 10000, 50000)) {
  nb <- max(1L, as.integer(w / bs))
  n_bins <- max(400L, nb * 4L + 80L)
  a <- floor(n_bins / 2) - floor(nb / 2)
  ev0 <- data.frame(start_bin = a, end_bin = a + nb - 1L, delta = 1.5, shape = "sharp")
  for (r in seq_len(nrep)) {
    sim <- simulate_continuous_events(n_bins, bs, ev0, noise_sd = 0.2,
                                      seed = 80000L + w + r)
    for (mode in c("single", "coarse_fine")) {
      se <- if (w <= 500) "fused_lasso" else "pelt"
      mw <- if (se == "pelt") 3L else 1L
      row <- eval_one(sim, "raw_effect", se, mw, scale_mode = mode)
      if (!is.null(row)) {
        ms_rows[[length(ms_rows) + 1L]] <- cbind(width = w, rep = r, row)
      }
    }
    # candidate union: fused_lasso min1 + pelt min3
    c1 <- eval_one(sim, "raw_effect", "fused_lasso", 1L)
    c2 <- eval_one(sim, "raw_effect", "pelt", 3L)
    u <- tryCatch({
      sc1 <- compute_diff_score(sim$signal, sim$samples, "JNKi", "DMSO", "raw_effect")
      thr <- 0.15 + 0.35 * stats::mad(sc1, constant = 1)
      if (!is.finite(thr) || thr <= 0) thr <- 0.2
      s1 <- segment_signal(sc1, positions = sim$bed$start, engine = "fused_lasso",
                           penalty = 8, min_width = 1L)
      s2 <- segment_signal(sc1, positions = sim$bed$start, engine = "pelt",
                           penalty = 8, min_width = 3L)
      u0 <- union_resolve_calls(list(
        narrow = .segments_to_calls(s1, sim$bed, thresh = thr),
        broad = .segments_to_calls(s2, sim$bed, thresh = thr)
      ))
      scu <- score_calls_vs_truth(u0, ev0, bs)
      cbind(score_engine = "raw_effect", segmentation_engine = "union_fused_pelt",
            min_width = 1L, inference = "same_data", scale_mode = "union",
            n_call = nrow(u0), scu)
    }, error = function(e) NULL)
    if (!is.null(u)) ms_rows[[length(ms_rows) + 1L]] <- cbind(width = w, rep = r, u)
  }
}
ms_tab <- do.call(rbind, ms_rows)
write_tsv(ms_tab, "results/v0.7/multiscale_benchmark.tsv")

# ---- extra null regimes (MAD-threshold, comparable to factorial detection) ----
message("extra nulls")
add_ar1 <- function(sim, phi = 0.7, seed) {
  set.seed(seed)
  n <- nrow(sim$signal)
  for (j in seq_len(ncol(sim$signal))) {
    e <- as.numeric(stats::filter(rnorm(n, 0, 0.2), phi, method = "recursive"))
    sim$signal[, j] <- pmax(sim$signal[, j] + e, 0)
  }
  sim
}
add_structured <- function(sim, seed) {
  set.seed(seed)
  n <- nrow(sim$signal)
  blk <- rep(c(0, 0.8, 0, 0.5), each = ceiling(n / 4))[seq_len(n)]
  sim$signal <- pmax(sim$signal + blk, 0)
  sim
}
cal_extra <- list()
for (regime in c("global", "low_signal", "structured", "autocorrelated")) {
  for (se in c("pelt", "fused_lasso", "wbs")) {
    for (r in seq_len(nrep)) {
      ev <- data.frame(start_bin = integer(), end_bin = integer(),
                       delta = numeric(), shape = character())
      ns <- if (regime == "low_signal") 0.05 else 0.2
      sim <- simulate_continuous_events(400L, bs, ev, noise_sd = ns,
                                        seed = 81000L + nchar(regime) * 40 + r)
      if (regime == "structured") sim <- add_structured(sim, 82000L + r)
      if (regime == "autocorrelated") sim <- add_ar1(sim, 0.7, 83000L + r)
      mw <- if (se == "pelt") 3L else 1L
      row <- eval_one(sim, "raw_effect", se, mw)
      n <- if (is.null(row)) 0L else row$n_call
      cal_extra[[length(cal_extra) + 1L]] <- data.frame(
        regime = regime, engine = se, rep = r, n_der = n,
        any_der = isTRUE(n > 0), stringsAsFactors = FALSE
      )
    }
  }
}
calx <- do.call(rbind, cal_extra)
# append to existing calibration with a regime column
cal0 <- readt("results/v0.7/calibration_metrics.tsv")
if (!"regime" %in% names(cal0)) cal0$regime <- "global_inference"
write_tsv(rbind(
  cbind(cal0[intersect(names(cal0), c("inference", "engine", "rep", "n_der",
                                      "any_der", "regime"))],
        source = "infer_regions"),
  data.frame(inference = "mad_threshold", engine = calx$engine, rep = calx$rep,
             n_der = calx$n_der, any_der = calx$any_der, regime = calx$regime,
             source = "mad_threshold", stringsAsFactors = FALSE)
), "results/v0.7/calibration_metrics.tsv")

# ---- power–calibration expansion (MAD-threshold, pelt/fused/wbs × 2kb) ----
message("power frontier expansion")
pw <- list()
for (se in c("pelt", "fused_lasso", "wbs")) {
  recs <- numeric()
  for (r in seq_len(nrep)) {
    ev <- data.frame(start_bin = 150L, end_bin = 159L, delta = 1.5, shape = "sharp")
    sim <- simulate_continuous_events(320L, bs, ev, noise_sd = 0.2, seed = 84000L + r)
    mw <- if (se == "pelt") 3L else 1L
    row <- eval_one(sim, "raw_effect", se, mw)
    recs <- c(recs, if (is.null(row)) 0 else row$recall)
  }
  fdp <- mean(calx$any_der[calx$regime == "global" & calx$engine == se])
  pw[[length(pw) + 1L]] <- data.frame(
    inference = "mad_threshold", engine = se,
    mean_recall = mean(recs, na.rm = TRUE), p_any_false = fdp,
    stringsAsFactors = FALSE
  )
}
# keep existing inference-architecture WBS rows
old_pw <- readt("results/v0.7/power_calibration_frontier.tsv")
write_tsv(rbind(old_pw, do.call(rbind, pw)),
          "results/v0.7/power_calibration_frontier.tsv")

# ---- held-out automatic routing (even train / odd test; no test truth in routing) ----
message("held-out routing")
hold <- list()
for (w in c(200, 500, 2000, 10000, 50000)) {
  nb <- max(1L, as.integer(w / bs))
  n_bins <- max(400L, nb * 4L + 80L)
  a <- floor(n_bins / 2) - floor(nb / 2)
  ev0 <- data.frame(start_bin = a, end_bin = a + nb - 1L, delta = 1.5, shape = "sharp")
  for (r in seq_len(nrep)) {
    sim <- simulate_continuous_events(n_bins, bs, ev0, noise_sd = 0.2,
                                      seed = 85000L + w + r)
    sc0 <- compute_diff_score(sim$signal, sim$samples, "JNKi", "DMSO", "raw_effect")
    rt <- route_production_engine(sc0, roughness_cut = 1.2)
    row <- eval_one(sim, rt$score_engine, rt$segmentation_engine, rt$min_width)
    if (is.null(row)) next
    split <- if (r %% 2L == 0L) "train" else "test"
    hold[[length(hold) + 1L]] <- cbind(
      width = w, rep = r, split = split, route = rt$route,
      roughness = rt$roughness, row
    )
  }
}
holdt <- do.call(rbind, hold)
write_tsv(holdt, "results/v0.7/heldout_routing.tsv")

# ---- exploratory effect × noise (n=10, labeled) ----
message("exploratory effect x noise n=10")
ex <- list()
for (w in c(500, 10000)) {
  for (delta in c(0.5, 1.5, 3.0)) {
    for (ns in c(0.1, 0.2, 0.4)) {
      nb <- max(1L, as.integer(w / bs))
      n_bins <- max(400L, nb * 4L + 80L)
      a <- floor(n_bins / 2) - floor(nb / 2)
      ev0 <- data.frame(start_bin = a, end_bin = a + nb - 1L, delta = delta, shape = "sharp")
      for (r in 1:10) {
        sim <- simulate_continuous_events(n_bins, bs, ev0, noise_sd = ns,
                                          seed = 86000L + w + as.integer(delta * 10) + r)
        for (se in c("pelt", "fused_lasso")) {
          mw <- if (se == "pelt") 3L else 1L
          row <- eval_one(sim, "raw_effect", se, mw)
          if (is.null(row)) next
          ex[[length(ex) + 1L]] <- cbind(width = w, delta = delta, noise = ns,
                                         rep = r, exploratory = TRUE, row)
        }
      }
    }
  }
}
write_tsv(do.call(rbind, ex), "results/v0.7/exploratory_effect_noise.tsv")

# ---- runtime placeholder with wall time ----
write_tsv(data.frame(
  step = c("synthetic_v07_primary", "synthetic_v07_finish"),
  nrep = nrep,
  status = "PASS",
  seconds = c(NA, proc.time()[["elapsed"]] - t0)
), "results/v0.7/runtime_memory.tsv")

# ---- figures ----
message("figures")
safe_mean <- function(d, y, by) {
  aggregate(d[[y]], d[by], function(z) mean(z, na.rm = TRUE))
}

pdf("figures/v0.7_score_engine_phase_diagram.pdf", width = 8, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (mw in c(1, 3)) {
  plot(1, type = "n", xlim = range(sc$width), ylim = c(0, 1), log = "x",
       xlab = "true width (bp)", ylab = "mean F1",
       main = paste0("score engines, PELT min_width=", mw))
  cols <- c(raw_effect = "black", log_effect = "steelblue",
            standardized = "darkorange", hotelling = "purple")
  for (e in names(cols)) {
    d <- sc[sc$score_engine == e & sc$min_width == mw, ]
    if (!nrow(d)) next
    a <- aggregate(f1 ~ width, d, mean)
    lines(a$width, a$f1, col = cols[e], lwd = 2)
    points(a$width, a$f1, col = cols[e], pch = 16)
  }
  legend("bottomright", legend = names(cols), col = cols, lwd = 2, bty = "n", cex = 0.7)
}
dev.off()

pdf("figures/v0.7_segmentation_engine_phase_diagram.pdf", width = 7, height = 5)
par(mar = c(4, 4, 2, 1))
plot(1, type = "n", xlim = range(sg$width), ylim = c(0, 1), log = "x",
     xlab = "true width (bp)", ylab = "mean F1 (zeros included)",
     main = "segmentation engines, raw_effect")
cols <- c(pelt = "black", binseg = "grey50", wbs = "steelblue",
          fused_lasso = "darkorange", hmm = "forestgreen", hsmm = "darkgreen")
for (e in names(cols)) {
  d <- sg[sg$segmentation_engine == e, ]
  a <- aggregate(f1 ~ width, d, mean)
  lines(a$width, a$f1, col = cols[e], lwd = 2)
  points(a$width, a$f1, col = cols[e], pch = 16)
}
legend("bottomright", legend = names(cols), col = cols, lwd = 2, bty = "n", cex = 0.7)
dev.off()

pdf("figures/v0.7_f1_vs_width.pdf", width = 6.5, height = 5)
par(mar = c(4, 4, 2, 1))
plot(1, type = "n", xlim = range(sg$width), ylim = c(0, 1), log = "x",
     xlab = "true width (bp)", ylab = "mean F1")
for (e in c("pelt", "fused_lasso", "wbs", "hmm")) {
  a <- aggregate(f1 ~ width, sg[sg$segmentation_engine == e, ], mean)
  lines(a$width, a$f1, lwd = 2,
        col = c(pelt = "black", fused_lasso = "darkorange", wbs = "steelblue",
                hmm = "forestgreen")[e])
}
std <- aggregate(f1 ~ width, sc[sc$score_engine == "standardized" & sc$min_width == 1, ], mean)
lines(std$width, std$f1, lwd = 2, lty = 2, col = "darkorange")
legend("bottomright",
       legend = c("raw+PELT", "raw+fused_lasso", "raw+WBS", "raw+HMM",
                  "standardized+PELT min1"),
       col = c("black", "darkorange", "steelblue", "forestgreen", "darkorange"),
       lwd = 2, lty = c(1, 1, 1, 1, 2), bty = "n", cex = 0.7)
dev.off()

pdf("figures/v0.7_boundary_error_vs_width.pdf", width = 6.5, height = 5)
par(mar = c(4, 4, 2, 1))
plot(1, type = "n", xlim = range(sg$width), ylim = c(1, 1e5), log = "xy",
     xlab = "true width (bp)", ylab = "mean |boundary error| (bp)")
for (e in c("pelt", "fused_lasso", "hmm", "hsmm", "wbs")) {
  d <- sg[sg$segmentation_engine == e & is.finite(sg$med_bound_err), ]
  if (!nrow(d)) next
  a <- aggregate(med_bound_err ~ width, d, mean)
  lines(a$width, pmax(a$med_bound_err, 1), lwd = 2,
        col = c(pelt = "black", fused_lasso = "darkorange", hmm = "forestgreen",
                hsmm = "darkgreen", wbs = "steelblue")[e])
}
legend("topleft", legend = c("pelt", "fused_lasso", "hmm", "hsmm", "wbs"),
       col = c("black", "darkorange", "forestgreen", "darkgreen", "steelblue"),
       lwd = 2, bty = "n")
dev.off()

ms <- readt("results/v0.7/merge_split_metrics.tsv")
pdf("figures/v0.7_merge_rate_vs_gap.pdf", width = 6, height = 5)
par(mar = c(4, 4, 2, 1))
plot(1, type = "n", xlim = c(-0.5, 10.5), ylim = c(0, 1.05),
     xlab = "gap (bins)", ylab = "mean n_merged")
for (e in c("pelt", "fused_lasso", "wbs")) {
  d <- ms[ms$segmentation_engine == e, ]
  a <- aggregate(n_merged ~ gap_bins, d, mean)
  lines(a$gap_bins, a$n_merged, lwd = 2,
        col = c(pelt = "black", fused_lasso = "darkorange", wbs = "steelblue")[e])
  points(a$gap_bins, a$n_merged, pch = 16,
         col = c(pelt = "black", fused_lasso = "darkorange", wbs = "steelblue")[e])
}
legend("topright", legend = c("pelt", "fused_lasso", "wbs"),
       col = c("black", "darkorange", "steelblue"), lwd = 2, bty = "n")
dev.off()

pdf("figures/v0.7_fragmentation_vs_width.pdf", width = 6, height = 5)
par(mar = c(4, 4, 2, 1))
plot(1, type = "n", xlim = range(sg$width), ylim = c(0, 3), log = "x",
     xlab = "true width (bp)", ylab = "mean n_call")
for (e in c("pelt", "fused_lasso", "wbs", "hmm")) {
  a <- aggregate(n_call ~ width, sg[sg$segmentation_engine == e, ], mean)
  lines(a$width, a$n_call, lwd = 2,
        col = c(pelt = "black", fused_lasso = "darkorange", wbs = "steelblue",
                hmm = "forestgreen")[e])
}
legend("topright", legend = c("pelt", "fused_lasso", "wbs", "hmm"),
       col = c("black", "darkorange", "steelblue", "forestgreen"), lwd = 2, bty = "n")
dev.off()

ne <- readt("results/v0.7/nested_event_metrics.tsv")
pdf("figures/v0.7_nested_scale_recovery.pdf", width = 6, height = 5)
par(mar = c(6, 4, 2, 1))
engs <- unique(ne$engine)
f1m <- sapply(engs, function(e) mean(ne$f1[ne$engine == e], na.rm = TRUE))
mgm <- sapply(engs, function(e) mean(ne$n_merged[ne$engine == e], na.rm = TRUE))
barplot(rbind(f1m, mgm), beside = TRUE, las = 2, ylim = c(0, 1.2),
        col = c("grey30", "tomato"),
        legend.text = c("F1", "mean n_merged"),
        args.legend = list(bty = "n"),
        main = "nested narrow-on-broad", ylab = "value")
dev.off()

pdf("figures/v0.7_empirical_fdr.pdf", width = 7, height = 5)
par(mar = c(8, 4, 2, 1))
pany <- aggregate(any_der ~ regime + engine, calx, mean)
labs <- paste(pany$engine, pany$regime, sep = " / ")
barplot(pany$any_der, names.arg = labs, las = 2, ylim = c(0, 1),
        ylab = "P(any false DER)", main = "MAD-threshold nulls", cex.names = 0.6)
abline(h = 0.05, lty = 2, col = "red")
dev.off()

pw2 <- readt("results/v0.7/power_calibration_frontier.tsv")
pdf("figures/v0.7_power_calibration_frontier.pdf", width = 6, height = 5)
par(mar = c(4, 4, 2, 1))
plot(pw2$p_any_false, pw2$mean_recall, xlim = c(0, max(0.2, max(pw2$p_any_false, na.rm = TRUE))),
     ylim = c(0, 1), xlab = "P(any false DER)", ylab = "mean recall",
     pch = 16, col = ifelse(pw2$inference == "mad_threshold", "darkorange", "black"))
text(pw2$p_any_false, pw2$mean_recall, paste(pw2$engine, pw2$inference),
     pos = 4, cex = 0.6)
dev.off()

pdf("figures/v0.7_oracle_gap.pdf", width = 6, height = 5)
par(mar = c(4, 4, 2, 1))
plot(oracle$width, oracle$f1, log = "x", type = "b", lwd = 2, ylim = c(0, 1),
     xlab = "true width (bp)", ylab = "mean F1", pch = 16)
lines(oracle$width, oracle$production_C_f1, type = "b", lwd = 2, col = "darkorange", pch = 17)
legend("bottomright", legend = c("oracle (best cell)", "candidate C reference"),
       col = c("black", "darkorange"), lwd = 2, pch = c(16, 17), bty = "n")
dev.off()

# runtime + placeholder realdata stability figure (updated after real run)
pdf("figures/v0.7_runtime_memory.pdf", width = 5, height = 4)
par(mar = c(4, 4, 2, 1))
barplot(c(primary = NA, finish = proc.time()[["elapsed"]] - t0),
        ylab = "seconds (this process)", main = "v0.7 finish synthetic")
dev.off()

# placeholder realdata figure; overwritten after real run if possible
if (!file.exists("figures/v0.7_realdata_boundary_stability.pdf")) {
  pdf("figures/v0.7_realdata_boundary_stability.pdf", width = 5, height = 4)
  plot.new()
  title("real-data stability: see realdata_stability.tsv")
  dev.off()
}

utils::write.table(data.frame(item = sessionInfo()$R.version$version.string),
                   "results/v0.7/session_info.txt", quote = FALSE, row.names = FALSE)
message("finish synthetic done in ", round(proc.time()[["elapsed"]] - t0), "s")
