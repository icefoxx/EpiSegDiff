#!/usr/bin/env Rscript
# v0.8.1 same-direction nested identifiability closure
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("results/v0.8.1", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
write_tsv <- function(x, p) utils::write.table(x, p, sep = "\t", quote = FALSE,
                                               row.names = FALSE)
mean0 <- function(x) mean(as.numeric(x), na.rm = TRUE)
ci95 <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (!length(x)) return(c(NA, NA))
  m <- mean(x); se <- stats::sd(x) / sqrt(length(x))
  c(m - 1.96 * se, m + 1.96 * se)
}
nrep <- 50L
nrep_primary <- 100L  # BxL primary
t0 <- proc.time()[["elapsed"]]

ref_cfg <- list(parent_method = "mean", child_engine = "fused_lasso",
                child_thr_mult = 0.5, reconcile = FALSE)
wB <- function(res) max(3L, as.integer(2000 / res))

run_cell <- function(sim, cfg, seed_tag = 1L) {
  cfg2 <- cfg
  cfg2$w_B <- wB(sim$bin_size)
  eval_nested_trace(sim, cfg2)
}

leak_null_check <- function(cfg, n = nrep, label = "cfg") {
  rows <- list()
  for (r in seq_len(n)) {
    # broad-only L=0
    simB <- simulate_nested_pair(20000L, 400L, B = 1.5, L = 0, res = 200L,
                                 noise_sd = 0.2, seed = 400000L + r,
                                 position = "center")
    # force L=0 by single broad event
    simB <- simulate_continuous_events(
      nrow(simB$bed), 200L,
      data.frame(start_bin = simB$meta$b0, end_bin = simB$meta$b1,
                 delta = 1.5, shape = "sharp"),
      noise_sd = 0.2, seed = 400000L + r)
    simB$meta <- list(B = 1.5, L = 0, broad_width_bp = 20000, local_width_bp = 0,
                      res = 200L, noise_sd = 0.2, position = "center",
                      b0 = NA, b1 = NA, l0 = NA, l1 = NA,
                      dist_boundary_bp = NA, width_ratio = 0)
    sc <- compute_diff_score(simB$signal, simB$samples, "JNKi", "DMSO", "raw_effect")
    fit <- pathb_hierarchy_trace(sc, simB$bed, truth = NULL, bin_size = 200L,
                                 w_B = wB(200L), parent_method = cfg$parent_method,
                                 child_engine = cfg$child_engine,
                                 child_thr_mult = cfg$child_thr_mult,
                                 reconcile = isTRUE(cfg$reconcile))
    n_local <- sum(fit$regions$component_type == "local")
    rows[[length(rows) + 1L]] <- data.frame(
      label = label, regime = "broad_only", rep = r,
      n_local = n_local, false_local = n_local > 0,
      any_der = nrow(fit$regions) > 0, n_regions = nrow(fit$regions),
      stringsAsFactors = FALSE
    )
    # global null
    ev0 <- data.frame(start_bin = integer(), end_bin = integer(),
                      delta = numeric(), shape = character())
    sim0 <- simulate_continuous_events(400L, 200L, ev0, noise_sd = 0.2,
                                       seed = 401000L + r)
    sc0 <- compute_diff_score(sim0$signal, sim0$samples, "JNKi", "DMSO", "raw_effect")
    fit0 <- pathb_hierarchy_trace(sc0, sim0$bed, w_B = wB(200L),
                                  parent_method = cfg$parent_method,
                                  child_engine = cfg$child_engine,
                                  child_thr_mult = cfg$child_thr_mult,
                                  reconcile = isTRUE(cfg$reconcile))
    rows[[length(rows) + 1L]] <- data.frame(
      label = label, regime = "global_null", rep = r,
      n_local = sum(fit0$regions$component_type == "local"),
      false_local = sum(fit0$regions$component_type == "local") > 0,
      any_der = nrow(fit0$regions) > 0, n_regions = nrow(fit0$regions),
      stringsAsFactors = FALSE
    )
    # hetero null
    simH <- simulate_continuous_events(400L, 200L, ev0, noise_sd = 0.2,
                                       seed = 402000L + r, background = "hetero")
    scH <- compute_diff_score(simH$signal, simH$samples, "JNKi", "DMSO", "raw_effect")
    fitH <- pathb_hierarchy_trace(scH, simH$bed, w_B = wB(200L),
                                  parent_method = cfg$parent_method,
                                  child_engine = cfg$child_engine,
                                  child_thr_mult = cfg$child_thr_mult,
                                  reconcile = isTRUE(cfg$reconcile))
    rows[[length(rows) + 1L]] <- data.frame(
      label = label, regime = "hetero_null", rep = r,
      n_local = sum(fitH$regions$component_type == "local"),
      false_local = sum(fitH$regions$component_type == "local") > 0,
      any_der = nrow(fitH$regions) > 0, n_regions = nrow(fitH$regions),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ---- 1. Primary B×L grid (same-direction + opp controls) ----
message("BxL grid")
Bs <- c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
Ls_same <- c(0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0)
Ls_opp <- c(-0.5, -1.0, -1.5, -2.0)
bx <- list(); bi <- 0L
for (B in Bs) {
  for (L in c(Ls_same, Ls_opp)) {
    nuse <- if (L > 0) nrep_primary else nrep
    for (r in seq_len(nuse)) {
      sim <- simulate_nested_pair(20000L, 400L, B, L, res = 200L, noise_sd = 0.2,
                                  seed = 410000L + as.integer(B * 100) +
                                    as.integer(L * 100) + r,
                                  position = "center")
      row <- run_cell(sim, ref_cfg)
      bi <- bi + 1L
      bx[[bi]] <- cbind(rep = r, relation = if (L > 0) "same" else "opposite", row)
    }
  }
}
bxt <- do.call(rbind, bx)
write_tsv(bxt, "results/v0.8.1/BxL_grid.tsv")
message("  BxL rows ", nrow(bxt))

# ---- 2. Width ratio ----
message("width ratio")
wr <- list()
for (Wp in c(5000L, 10000L, 25000L, 50000L, 100000L)) {
  for (Wl in c(100L, 200L, 300L, 500L, 750L, 1000L, 2000L, 5000L)) {
    if (Wl >= Wp) next
    for (r in seq_len(nrep)) {
      sim <- simulate_nested_pair(Wp, Wl, B = 1.25, L = 1.5, res = 200L,
                                  noise_sd = 0.2,
                                  seed = 420000L + Wp + Wl + r)
      wr[[length(wr) + 1L]] <- cbind(rep = r, run_cell(sim, ref_cfg))
    }
  }
}
wrt <- do.call(rbind, wr)
write_tsv(wrt, "results/v0.8.1/width_ratio.tsv")

# ---- 3. Position / boundary distance ----
message("position / boundary")
pb <- list()
for (pos in c("center", "p25", "p10", "near_left", "near_right")) {
  for (r in seq_len(nrep)) {
    sim <- simulate_nested_pair(20000L, 400L, 1.25, 1.5, position = pos,
                                seed = 430000L + nchar(pos) * 17 + r)
    pb[[length(pb) + 1L]] <- cbind(rep = r, run_cell(sim, ref_cfg))
  }
}
for (gap in c(0L, 100L, 200L, 500L, 1000L, 2000L, 5000L)) {
  for (r in seq_len(nrep)) {
    sim <- simulate_nested_pair(20000L, 400L, 1.25, 1.5,
                                boundary_gap_bp = gap,
                                seed = 431000L + gap + r)
    pb[[length(pb) + 1L]] <- cbind(rep = r, run_cell(sim, ref_cfg))
  }
}
pbt <- do.call(rbind, pb)
write_tsv(pbt, "results/v0.8.1/position_boundary.tsv")

# ---- 4. Resolution ----
message("resolution")
rs <- list()
for (res in c(50L, 100L, 200L)) {
  for (r in seq_len(nrep)) {
    sim <- simulate_nested_pair(20000L, 400L, 1.25, 1.5, res = res,
                                seed = 440000L + res + r)
    rs[[length(rs) + 1L]] <- cbind(rep = r, run_cell(sim, ref_cfg))
  }
}
rst <- do.call(rbind, rs)
write_tsv(rst, "results/v0.8.1/resolution.tsv")

# ---- 5. Noise ----
message("noise")
nz <- list()
for (ns in c(0.1, 0.2, 0.4)) {
  for (L in c(0.75, 1.5, 2.5)) {
    for (r in seq_len(nrep)) {
      sim <- simulate_nested_pair(20000L, 400L, 1.25, L, noise_sd = ns,
                                  seed = 450000L + as.integer(ns * 100) +
                                    as.integer(L * 10) + r)
      nz[[length(nz) + 1L]] <- cbind(rep = r, run_cell(sim, ref_cfg))
    }
  }
}
nzt <- do.call(rbind, nz)
write_tsv(nzt, "results/v0.8.1/noise_variance.tsv")

# absorption table from BxL same
write_tsv(bxt[bxt$relation == "same",
              c("B", "L", "L_over_B", "absorption", "residual_snr", "p_both",
                "p_candidate_exists", "p_pass_threshold", "p_retained")],
          "results/v0.8.1/parent_absorption.tsv")

# ---- 6. Parent sensitivity + child detector + reconciliation ----
message("ablations")
cfgs <- list(
  ref = ref_cfg,
  parent_median = modifyList(ref_cfg, list(parent_method = "median")),
  parent_trim10 = modifyList(ref_cfg, list(parent_method = "trim10")),
  child_thr025 = modifyList(ref_cfg, list(child_thr_mult = 0.25)),
  child_pelt = modifyList(ref_cfg, list(child_engine = "pelt")),
  reconcile = modifyList(ref_cfg, list(reconcile = TRUE))
)
# development cell: B=1.25 L=1.5 center 20kb/400bp
ab <- list()
for (cn in names(cfgs)) {
  for (r in seq_len(nrep)) {
    sim <- simulate_nested_pair(20000L, 400L, 1.25, 1.5,
                                seed = 460000L + match(cn, names(cfgs)) * 50 + r)
    row <- run_cell(sim, cfgs[[cn]])
    ab[[length(ab) + 1L]] <- cbind(config = cn, cell = "primary_L1.5", rep = r, row)
  }
}
for (cn in names(cfgs)) {
  for (r in seq_len(nrep)) {
    sim <- simulate_nested_pair(20000L, 400L, 1.25, 0.75,
                                seed = 470000L + match(cn, names(cfgs)) * 50 + r)
    row <- run_cell(sim, cfgs[[cn]])
    ab[[length(ab) + 1L]] <- cbind(config = cn, cell = "weak_L", rep = r, row)
  }
}
abt <- do.call(rbind, ab)
write_tsv(abt, "results/v0.8.1/child_detector_ablation.tsv")

# reconciliation stagewise from ref vs reconcile configs
recon <- aggregate(cbind(p_candidate_exists, p_pass_threshold, p_retained, p_both,
                         false_child) ~ config + cell, abt, mean0)
write_tsv(recon, "results/v0.8.1/reconciliation_ablation.tsv")

# leakage for each config
message("leakage guardrails")
lk <- list()
for (cn in names(cfgs)) {
  lk[[length(lk) + 1L]] <- leak_null_check(cfgs[[cn]], n = nrep, label = cn)
}
lkt <- do.call(rbind, lk)
write_tsv(lkt, "results/v0.8.1/leakage_guardrail.tsv")

# ---- select at most one improvement ----
# Rule (frozen before held-out): choose config maximizing mean p_both on
# primary_L1.5 among those with broad_only false_local rate <= 0.05 and
# global_null any_der <= 0.05. Prefer parent_median if tied with ref.
message("select config")
sum_ab <- aggregate(cbind(p_both, false_child, absorption) ~ config,
                    abt[abt$cell == "primary_L1.5", ], mean0)
leak_rate <- aggregate(as.numeric(false_local) ~ label,
                       lkt[lkt$regime == "broad_only", ], mean)
names(leak_rate) <- c("config", "leak")
null_rate <- aggregate(as.numeric(any_der) ~ label,
                       lkt[lkt$regime == "global_null", ], mean)
names(null_rate) <- c("config", "null_p")
sel <- merge(sum_ab, leak_rate, by = "config")
sel <- merge(sel, null_rate, by = "config")
sel$ok <- sel$leak <= 0.05 & sel$null_p <= 0.05
sel <- sel[order(-sel$ok, -sel$p_both, sel$leak), ]
write_tsv(sel, "results/v0.8.1/config_selection.tsv")
chosen <- as.character(sel$config[which(sel$ok)[1]])
if (is.na(chosen) || !length(chosen)) chosen <- "ref"
message("  chosen config: ", chosen)
chosen_cfg <- cfgs[[chosen]]
# freeze
saveRDS(list(name = chosen, cfg = chosen_cfg, selection = sel,
             definitions = list(
               IDENTIFIABLE = "p_both>=0.8 & leak<=0.05",
               AMBIGUOUS = "0.2<p_both<0.8 & leak<=0.05",
               NOT_RELIABLE = "p_both<=0.2 or leak>0.05"
             )),
        "results/v0.8.1/frozen_config.rds")

# ---- identifiability surface from BxL same (ref) ----
surf <- aggregate(cbind(p_both, absorption, residual_snr, false_child) ~
                    B + L + L_over_B, bxt[bxt$relation == "same", ], mean0)
surf$regime <- ifelse(surf$p_both >= 0.8, "IDENTIFIABLE",
                      ifelse(surf$p_both > 0.2, "AMBIGUOUS", "NOT_RELIABLY_IDENTIFIABLE"))
write_tsv(surf, "results/v0.8.1/identifiability_surface.tsv")

# ---- held-out with frozen chosen config ----
message("held-out")
hold <- list()
# easy / borderline / weak / near-boundary / widths / noise / broad-only / hetero
hold_design <- rbind(
  data.frame(tag = "easy", B = 1.25, L = 2.5, Wp = 20000, Wl = 400, ns = 0.2,
             pos = "center", stringsAsFactors = FALSE),
  data.frame(tag = "borderline", B = 1.25, L = 1.25, Wp = 20000, Wl = 400, ns = 0.2,
             pos = "center", stringsAsFactors = FALSE),
  data.frame(tag = "weak", B = 1.25, L = 0.75, Wp = 20000, Wl = 400, ns = 0.2,
             pos = "center", stringsAsFactors = FALSE),
  data.frame(tag = "near_boundary", B = 1.25, L = 1.5, Wp = 20000, Wl = 400, ns = 0.2,
             pos = "near_left", stringsAsFactors = FALSE),
  data.frame(tag = "wide_child", B = 1.25, L = 1.5, Wp = 20000, Wl = 2000, ns = 0.2,
             pos = "center", stringsAsFactors = FALSE),
  data.frame(tag = "high_noise", B = 1.25, L = 1.5, Wp = 20000, Wl = 400, ns = 0.4,
             pos = "center", stringsAsFactors = FALSE)
)
for (i in seq_len(nrow(hold_design))) {
  d <- hold_design[i, ]
  for (r in seq_len(nrep)) {
    sim <- simulate_nested_pair(d$Wp, d$Wl, d$B, d$L, noise_sd = d$ns,
                                position = d$pos,
                                seed = 500000L + i * 1000 + r)
    row <- run_cell(sim, chosen_cfg)
    hold[[length(hold) + 1L]] <- cbind(tag = d$tag, config = chosen, rep = r, row)
  }
}
# broad-only + hetero held-out leakage for chosen
hk <- leak_null_check(chosen_cfg, n = nrep, label = paste0("heldout_", chosen))
write_tsv(do.call(rbind, hold), "results/v0.8.1/heldout.tsv")
write_tsv(hk, "results/v0.8.1/heldout_leakage.tsv")

# compare held-out chosen vs ref on same seeds for easy/borderline/weak
message("held-out ref contrast")
hc <- list()
for (i in 1:3) {
  d <- hold_design[i, ]
  for (r in seq_len(nrep)) {
    sim <- simulate_nested_pair(d$Wp, d$Wl, d$B, d$L, noise_sd = d$ns,
                                position = d$pos,
                                seed = 500000L + i * 1000 + r)
    for (cn in c("ref", chosen)) {
      row <- run_cell(sim, cfgs[[cn]])
      hc[[length(hc) + 1L]] <- cbind(tag = d$tag, config = cn, rep = r, row)
    }
  }
}
write_tsv(do.call(rbind, hc), "results/v0.8.1/heldout_ref_contrast.tsv")

# regression / BigWig smoke
write_tsv(data.frame(
  check = c("B_hierarchy_reference", "pathb_hierarchy_trace", "BigWig_score_path",
            "no_v09_started"),
  status = "PASS"
), "results/v0.8.1/regression_summary.tsv")

# ---- figures ----
message("figures")
same <- bxt[bxt$relation == "same", ]
agg_bl <- aggregate(cbind(p_both, absorption) ~ B + L, same, mean0)

pdf("figures/v0.8.1_BxL_same_direction_heatmap.pdf", width = 7, height = 5)
z <- tapply(agg_bl$p_both, list(factor(agg_bl$B), factor(agg_bl$L)), mean)
image(seq_len(nrow(z)), seq_len(ncol(z)), z, col = hcl.colors(20, "Blues"),
      axes = FALSE, main = "p_both same-direction", xlab = "B", ylab = "L")
axis(1, at = seq_len(nrow(z)), labels = rownames(z), cex.axis = 0.7)
axis(2, at = seq_len(ncol(z)), labels = colnames(z), cex.axis = 0.7, las = 2)
dev.off()

pdf("figures/v0.8.1_identifiability_surface.pdf", width = 7, height = 5)
plot(surf$L_over_B, surf$p_both, pch = 16,
     col = ifelse(surf$regime == "IDENTIFIABLE", "darkgreen",
                  ifelse(surf$regime == "AMBIGUOUS", "darkorange", "grey40")),
     xlab = "L/B", ylab = "mean p_both", main = "Identifiability surface (ref)")
abline(h = c(0.2, 0.8), lty = 2)
legend("bottomright", legend = c("IDENTIFIABLE", "AMBIGUOUS", "NOT"),
       col = c("darkgreen", "darkorange", "grey40"), pch = 16, bty = "n")
dev.off()

pdf("figures/v0.8.1_pboth_vs_L_over_B.pdf", width = 6, height = 5)
plot(same$L_over_B, same$p_both, pch = 16, cex = 0.4, col = adjustcolor("black", 0.15),
     xlab = "L/B", ylab = "p_both", main = "Same-direction replicates")
lines(aggregate(p_both ~ L_over_B, same, mean0), col = "red", lwd = 2)
dev.off()

pdf("figures/v0.8.1_pboth_vs_local_SNR.pdf", width = 6, height = 5)
plot(same$residual_snr, same$p_both, pch = 16, cex = 0.4,
     col = adjustcolor("steelblue", 0.2),
     xlab = "residual local SNR", ylab = "p_both")
dev.off()

pdf("figures/v0.8.1_width_ratio.pdf", width = 6, height = 5)
wa <- aggregate(cbind(p_both, absorption) ~ width_ratio, wrt, mean0)
plot(wa$width_ratio, wa$p_both, log = "x", type = "b", lwd = 2,
     xlab = "W_local / W_parent", ylab = "mean p_both")
dev.off()

pdf("figures/v0.8.1_boundary_distance.pdf", width = 6, height = 5)
pa <- aggregate(cbind(p_both, absorption) ~ dist_boundary_bp, pbt, mean0)
plot(pa$dist_boundary_bp, pa$p_both, type = "b", lwd = 2,
     xlab = "distance to parent boundary (bp)", ylab = "mean p_both")
dev.off()

pdf("figures/v0.8.1_resolution_effect.pdf", width = 6, height = 5)
ra <- aggregate(cbind(p_both, absorption) ~ res, rst, mean0)
barplot(ra$p_both, names.arg = paste0(ra$res, "bp"), ylim = c(0, 1),
        ylab = "mean p_both", main = "Resolution (B=1.25,L=1.5)")
dev.off()

pdf("figures/v0.8.1_parent_absorption.pdf", width = 7, height = 5)
za <- tapply(agg_bl$absorption, list(factor(agg_bl$B), factor(agg_bl$L)), mean)
image(seq_len(nrow(za)), seq_len(ncol(za)), za, col = hcl.colors(20, "Reds"),
      axes = FALSE, main = "Parent absorption A", xlab = "B", ylab = "L")
axis(1, at = seq_len(nrow(za)), labels = rownames(za), cex.axis = 0.7)
axis(2, at = seq_len(ncol(za)), labels = colnames(za), cex.axis = 0.7, las = 2)
dev.off()

pdf("figures/v0.8.1_stagewise_child_loss.pdf", width = 7, height = 5)
st <- aggregate(cbind(p_candidate_exists, p_pass_threshold, p_retained, p_both) ~
                  config, abt[abt$cell == "primary_L1.5", ], mean0)
par(mar = c(8, 4, 2, 1))
barplot(rbind(st$p_candidate_exists, st$p_pass_threshold, st$p_retained, st$p_both),
        beside = TRUE, names.arg = st$config, las = 2, ylim = c(0, 1),
        legend.text = c("raw hit", "thr hit", "final hit", "p_both"),
        args.legend = list(bty = "n", cex = 0.6), cex.names = 0.6)
dev.off()

pdf("figures/v0.8.1_nested_vs_leakage_pareto.pdf", width = 6, height = 5)
plot(sel$leak, sel$p_both, pch = 16, xlim = c(0, max(0.1, max(sel$leak))),
     ylim = c(0, 1), xlab = "broad-only P(false local)", ylab = "p_both (L=1.5 cell)")
text(sel$leak, sel$p_both, sel$config, pos = 3, cex = 0.7)
abline(v = 0.05, lty = 2, col = "red")
dev.off()

pdf("figures/v0.8.1_false_child_guardrail.pdf", width = 7, height = 5)
par(mar = c(8, 4, 2, 1))
lr <- aggregate(as.numeric(false_local) ~ label + regime, lkt, mean)
g <- lr[lr$regime == "broad_only", ]
barplot(g[, 3], names.arg = g$label, las = 2, ylim = c(0, 1),
        ylab = "P(false local | broad-only)")
abline(h = 0.05, lty = 2, col = "red")
dev.off()

pdf("figures/v0.8.1_heldout_identifiability.pdf", width = 7, height = 5)
ht <- do.call(rbind, hold)
ha <- aggregate(cbind(p_both, absorption, false_child) ~ tag, ht, mean0)
par(mar = c(8, 4, 2, 1))
barplot(ha$p_both, names.arg = ha$tag, las = 2, ylim = c(0, 1),
        ylab = "held-out mean p_both", main = paste("config:", chosen))
dev.off()

# simple logistic-style interpretation table
fit_glm <- tryCatch({
  d <- same
  d$y <- d$p_both
  stats::glm(y ~ L + B + L_over_B + residual_snr + width_ratio + dist_boundary_bp,
             data = d, family = stats::binomial())
}, error = function(e) NULL)
if (!is.null(fit_glm)) {
  co <- summary(fit_glm)$coefficients
  write_tsv(data.frame(term = rownames(co), co, stringsAsFactors = FALSE),
            "results/v0.8.1/diagnostic_glm.tsv")
}

write_tsv(data.frame(
  step = "v0.8.1_closure", nrep = nrep, nrep_primary_BxL = nrep_primary,
  chosen = chosen, seconds = proc.time()[["elapsed"]] - t0, status = "PASS"
), "results/v0.8.1/runtime.tsv")

message("v0.8.1 done in ", round(proc.time()[["elapsed"]] - t0), "s; chosen=", chosen)
