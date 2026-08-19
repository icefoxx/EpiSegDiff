#!/usr/bin/env Rscript
# v0.7-review closure R1–R5. Does not implement a v0.8 architecture.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("results/v0.8", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
write_tsv <- function(x, p) utils::write.table(x, p, sep = "\t", quote = FALSE,
                                               row.names = FALSE)
nrep <- 50L
t0 <- proc.time()[["elapsed"]]

configs <- list(
  raw_pelt_min3 = list(score = "raw_effect", seg = "pelt", mw = 3L),
  raw_pelt_min1 = list(score = "raw_effect", seg = "pelt", mw = 1L),
  std_pelt_min1 = list(score = "standardized", seg = "pelt", mw = 1L),
  raw_fused_min1 = list(score = "raw_effect", seg = "fused_lasso", mw = 1L)
)

eval_cfg <- function(sim, cfg, penalty = 8) {
  calls <- tryCatch(
    closure_calls(sim, cfg$score, cfg$seg, cfg$mw, penalty),
    error = function(e) NULL
  )
  if (is.null(calls)) return(NULL)
  sc <- score_calls_vs_truth(calls, sim$truth, sim$bin_size)
  list(calls = calls, sc = sc)
}

# ---- R1 resolution × narrow-width ----
message("R1 resolution x width")
r1 <- list(); i1 <- 0L
for (res in c(50L, 100L, 200L)) {
  for (w in c(100L, 200L, 300L, 500L, 750L, 1000L, 2000L)) {
    nb <- max(1L, as.integer(w / res))
    span <- max(40000L, as.integer(w * 4L + 8000L))
    n_bins <- as.integer(span / res)
    a <- floor(n_bins / 2) - floor(nb / 2)
    ev <- data.frame(start_bin = a, end_bin = a + nb - 1L, delta = 1.5,
                     shape = "sharp")
    for (cn in names(configs)) {
      cfg <- configs[[cn]]
      for (r in seq_len(nrep)) {
        sim <- simulate_continuous_events(n_bins, res, ev, noise_sd = 0.2,
                                          seed = 90000L + res + w +
                                            match(cn, names(configs)) * 19 + r)
        out <- eval_cfg(sim, cfg)
        if (is.null(out)) next
        i1 <- i1 + 1L
        r1[[i1]] <- cbind(resolution = res, width = w, config = cn,
                          score_engine = cfg$score,
                          segmentation_engine = cfg$seg,
                          min_width = cfg$mw, n_event_bins = nb,
                          rep = r, out$sc)
      }
    }
  }
}
r1t <- do.call(rbind, r1)
write_tsv(r1t, "results/v0.8/r1_resolution_width.tsv")
message("  R1 rows ", nrow(r1t))

# ---- R2 same-direction gap curve ----
message("R2 gap curve")
r2 <- list()
for (w in c(2000L, 10000L)) {
  for (gap in c(0L, 1L, 2L, 3L, 5L, 10L, 20L)) {
    for (delta in c(1.5)) {
      res <- 200L
      nb <- max(1L, as.integer(w / res))
      n_bins <- max(400L, 2L * nb + gap + 80L)
      a1 <- 40L
      e1 <- a1 + nb - 1L
      a2 <- e1 + 1L + gap
      e2 <- a2 + nb - 1L
      ev <- data.frame(start_bin = c(a1, a2), end_bin = c(e1, e2),
                       delta = c(delta, delta), shape = "sharp")
      for (cn in c("raw_pelt_min3", "raw_fused_min1")) {
        cfg <- configs[[cn]]
        for (r in seq_len(nrep)) {
          sim <- simulate_continuous_events(n_bins, res, ev, noise_sd = 0.2,
                                            seed = 91000L + w + gap * 17 + r +
                                              match(cn, names(configs)) * 3)
          out <- eval_cfg(sim, cfg)
          if (is.null(out)) next
          r2[[length(r2) + 1L]] <- cbind(
            width = w, gap_bins = gap, gap_bp = gap * res, delta = delta,
            config = cn, score_engine = cfg$score,
            segmentation_engine = cfg$seg, rep = r, out$sc
          )
        }
      }
    }
  }
}
r2t <- do.call(rbind, r2)
write_tsv(r2t, "results/v0.8/r2_gap_curve.tsv")
message("  R2 rows ", nrow(r2t))

# ---- R3 nested B × L ----
message("R3 nested B x L")
r3 <- list()
res <- 200L
broad_w <- 100L  # 20 kb
local_w <- 2L    # 400 bp
n_bins <- 220L
b0 <- 50L
b1 <- b0 + broad_w - 1L
l0 <- 98L
l1 <- l0 + local_w - 1L
for (B in c(0.6, 1.2, 1.8)) {
  for (L in c(-1.8, -1.2, -0.6, 0.6, 1.2, 1.8)) {
    ev <- data.frame(start_bin = c(b0, l0), end_bin = c(b1, l1),
                     delta = c(B, L), shape = "sharp")
    for (cn in c("raw_pelt_min3", "raw_fused_min1", "std_pelt_min1")) {
      cfg <- configs[[cn]]
      for (r in seq_len(nrep)) {
        sim <- simulate_continuous_events(n_bins, res, ev, noise_sd = 0.2,
                                          seed = 92000L + as.integer(B * 10) +
                                            as.integer(L * 10) + r +
                                            match(cn, names(configs)) * 7)
        out <- eval_cfg(sim, cfg)
        if (is.null(out)) next
        r3[[length(r3) + 1L]] <- cbind(
          B = B, L = L, relation = if (B * L > 0) "same" else "opposite",
          config = cn, score_engine = cfg$score,
          segmentation_engine = cfg$seg, rep = r, out$sc
        )
      }
    }
  }
}
r3t <- do.call(rbind, r3)
write_tsv(r3t, "results/v0.8/r3_nested_BL.tsv")
message("  R3 rows ", nrow(r3t))

# ---- R4 transition-zone gradual ----
message("R4 transition zones")
r4 <- list()
for (w in c(5000L, 10000L, 25000L)) {
  res <- 200L
  nb <- max(5L, as.integer(w / res))
  n_bins <- max(200L, nb * 4L + 40L)
  a <- floor(n_bins / 2) - floor(nb / 2)
  ev <- data.frame(start_bin = a, end_bin = a + nb - 1L, delta = 1.5,
                   shape = "gradual")
  series_end <- n_bins * res
  for (cn in c("raw_pelt_min3", "raw_fused_min1")) {
    cfg <- configs[[cn]]
    for (r in seq_len(nrep)) {
      sim <- simulate_continuous_events(n_bins, res, ev, noise_sd = 0.2,
                                        seed = 93000L + w + r +
                                          match(cn, names(configs)))
      out <- eval_cfg(sim, cfg)
      if (is.null(out)) next
      tz <- score_transition_zones(out$calls, a, a + nb - 1L, res, series_end)
      r4[[length(r4) + 1L]] <- cbind(
        width = w, config = cn, score_engine = cfg$score,
        segmentation_engine = cfg$seg, rep = r, out$sc, tz
      )
    }
  }
}
r4t <- do.call(rbind, r4)
write_tsv(r4t, "results/v0.8/r4_transition_zones.tsv")
message("  R4 rows ", nrow(r4t))

# ---- R5 detection vs conditional boundary ----
message("R5 decomposition")
decomp_by <- function(d, keys) {
  spl <- split(d, d[keys], drop = TRUE)
  do.call(rbind, lapply(spl, function(x) {
    cbind(x[1, keys, drop = FALSE], decompose_detection_boundary(x))
  }))
}
r5_r1 <- decomp_by(r1t, c("resolution", "width", "config"))
r5_r1$source <- "R1"
r5_r2 <- decomp_by(r2t, c("width", "gap_bins", "config"))
r5_r2$source <- "R2"
r5_r3 <- decomp_by(r3t, c("B", "L", "config"))
r5_r3$source <- "R3"
r5_r4 <- decomp_by(r4t, c("width", "config"))
r5_r4$source <- "R4"
# pad to common columns
pad <- function(x) {
  need <- c("source", "resolution", "width", "config", "gap_bins", "B", "L",
            "n", "p_detect", "f1_uncond", "recall_uncond", "precision_uncond",
            "med_bound_err_cond", "mean_bound_err_cond")
  for (nm in need) if (!nm %in% names(x)) x[[nm]] <- NA
  x[need]
}
r5t <- rbind(pad(r5_r1), pad(r5_r2), pad(r5_r3), pad(r5_r4))
write_tsv(r5t, "results/v0.8/r5_detection_boundary.tsv")

# also reuse v0.7 seg table
if (file.exists("results/v0.7/segmentation_engine_benchmark.tsv")) {
  v7 <- utils::read.table("results/v0.7/segmentation_engine_benchmark.tsv",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  v7$source <- "v0.7_seg"
  v7$config <- paste(v7$score_engine, v7$segmentation_engine, sep = "_")
  v7d <- decomp_by(v7, c("width", "config"))
  v7d$source <- "v0.7_seg"
  write_tsv(pad(v7d), "results/v0.8/r5_v07_reuse.tsv")
}

write_tsv(data.frame(
  step = c("R1", "R2", "R3", "R4", "R5"),
  nrep = nrep,
  rows = c(nrow(r1t), nrow(r2t), nrow(r3t), nrow(r4t), nrow(r5t)),
  seconds = proc.time()[["elapsed"]] - t0,
  status = "PASS"
), "results/v0.8/closure_runtime.tsv")

# ---- figures ----
message("figures")
mean0 <- function(x) mean(x, na.rm = TRUE)

pdf("figures/v0.8_r1_f1_resolution_width.pdf", width = 9, height = 4.2)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
cols <- c(raw_pelt_min3 = "black", raw_pelt_min1 = "grey40",
          std_pelt_min1 = "darkorange", raw_fused_min1 = "steelblue")
for (res in c(50, 100, 200)) {
  plot(1, type = "n", xlim = c(80, 2200), ylim = c(0, 1), log = "x",
       xlab = "true width (bp)", ylab = "mean F1",
       main = paste0("resolution ", res, " bp"))
  for (cn in names(cols)) {
    d <- r1t[r1t$resolution == res & r1t$config == cn, ]
    a <- aggregate(f1 ~ width, d, mean0)
    lines(a$width, a$f1, col = cols[cn], lwd = 2)
    points(a$width, a$f1, col = cols[cn], pch = 16)
  }
  if (res == 50)
    legend("bottomright", legend = names(cols), col = cols, lwd = 2,
           bty = "n", cex = 0.6)
}
dev.off()

pdf("figures/v0.8_r2_merge_vs_gap.pdf", width = 7, height = 4.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (w in c(2000, 10000)) {
  plot(1, type = "n", xlim = c(-0.5, 20.5), ylim = c(0, 1.05),
       xlab = "gap (bins @ 200 bp)", ylab = "mean n_merged",
       main = paste0("event width ", w, " bp"))
  for (cn in c("raw_pelt_min3", "raw_fused_min1")) {
    d <- r2t[r2t$width == w & r2t$config == cn, ]
    a <- aggregate(n_merged ~ gap_bins, d, mean0)
    lines(a$gap_bins, a$n_merged, lwd = 2,
          col = if (cn == "raw_pelt_min3") "black" else "steelblue")
    points(a$gap_bins, a$n_merged, pch = 16,
           col = if (cn == "raw_pelt_min3") "black" else "steelblue")
  }
  legend("topright", legend = c("PELT min3", "fused"),
         col = c("black", "steelblue"), lwd = 2, bty = "n")
}
dev.off()

pdf("figures/v0.8_r3_nested_merge.pdf", width = 7.5, height = 4.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (rel in c("same", "opposite")) {
  d0 <- r3t[r3t$relation == rel, ]
  plot(1, type = "n", xlim = range(d0$L), ylim = c(0, 1.05),
       xlab = "local additional L", ylab = "mean n_merged",
       main = paste0("nested ", rel, "-direction"))
  for (cn in unique(d0$config)) {
    a <- aggregate(n_merged ~ L, d0[d0$config == cn, ], mean0)
    lines(a$L, a$n_merged, lwd = 2)
    points(a$L, a$n_merged, pch = 16)
  }
}
dev.off()

pdf("figures/v0.8_r4_transition.pdf", width = 7, height = 4.5)
par(mar = c(4, 4, 2, 1))
plot(1, type = "n", xlim = c(4000, 26000), ylim = c(0, 1), log = "x",
     xlab = "planted width (bp)", ylab = "mean core recall / over-extension")
for (cn in c("raw_pelt_min3", "raw_fused_min1")) {
  d <- r4t[r4t$config == cn, ]
  a <- aggregate(cbind(core_recall, over_extension) ~ width, d, mean0)
  col <- if (cn == "raw_pelt_min3") "black" else "steelblue"
  lines(a$width, a$core_recall, col = col, lwd = 2)
  lines(a$width, a$over_extension, col = col, lwd = 2, lty = 2)
}
legend("right",
       legend = c("PELT core", "PELT over-ext", "fused core", "fused over-ext"),
       col = c("black", "black", "steelblue", "steelblue"),
       lty = c(1, 2, 1, 2), lwd = 2, bty = "n", cex = 0.7)
dev.off()

pdf("figures/v0.8_r5_pdetect_vs_error.pdf", width = 6, height = 5)
par(mar = c(4, 4, 2, 1))
d <- r5_r1
plot(d$p_detect, pmax(d$med_bound_err_cond, 1), log = "y",
     xlab = "P(detected)", ylab = "median |boundary error| | detected (bp)",
     pch = 16, col = as.integer(factor(d$config)),
     main = "R1: detection vs conditional error")
legend("topright", legend = unique(d$config),
       col = seq_along(unique(d$config)), pch = 16, bty = "n", cex = 0.6)
dev.off()

message("closure done in ", round(proc.time()[["elapsed"]] - t0), "s")
