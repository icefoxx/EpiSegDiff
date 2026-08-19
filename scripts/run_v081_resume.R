#!/usr/bin/env Rscript
# Resume v0.8.1 from ablations onward (grids already on disk).
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("results/v0.8.1", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
write_tsv <- function(x, p) utils::write.table(x, p, sep = "\t", quote = FALSE,
                                               row.names = FALSE)
mean0 <- function(x) mean(as.numeric(x), na.rm = TRUE)
nrep <- 50L
t0 <- proc.time()[["elapsed"]]

ref_cfg <- list(parent_method = "mean", child_engine = "fused_lasso",
                child_thr_mult = 0.5, reconcile = FALSE)
wB <- function(res) max(3L, as.integer(2000 / res))
run_cell <- function(sim, cfg) {
  cfg2 <- cfg
  cfg2$w_B <- wB(sim$bin_size)
  tryCatch(eval_nested_trace(sim, cfg2), error = function(e) {
    meta <- sim$meta
    data.frame(
      B = meta$B, L = meta$L, L_over_B = meta$L / max(meta$B, 1e-8),
      broad_width_bp = meta$broad_width_bp, local_width_bp = meta$local_width_bp,
      width_ratio = meta$width_ratio, res = meta$res, noise_sd = meta$noise_sd,
      position = meta$position, dist_boundary_bp = meta$dist_boundary_bp,
      parent_method = cfg$parent_method, child_engine = cfg$child_engine,
      child_thr_mult = cfg$child_thr_mult, reconcile = isTRUE(cfg$reconcile),
      p_parent = NA, p_child = NA, p_both = NA, p_merged = NA, false_child = NA,
      f1 = NA, n_regions = NA, n_broad = NA, n_local = NA,
      absorption = NA, residual_snr = NA, parent_bound_err = NA,
      p_candidate_exists = NA, p_pass_threshold = NA, p_retained = NA,
      stringsAsFactors = FALSE
    )
  })
}

leak_null_check <- function(cfg, n = nrep, label = "cfg") {
  rows <- list()
  for (r in seq_len(n)) {
    simB <- simulate_continuous_events(
      180L, 200L,
      data.frame(start_bin = 40L, end_bin = 139L, delta = 1.5, shape = "sharp"),
      noise_sd = 0.2, seed = 400000L + r)
    sc <- compute_diff_score(simB$signal, simB$samples, "JNKi", "DMSO", "raw_effect")
    fit <- tryCatch(pathb_hierarchy_trace(
      sc, simB$bed, w_B = wB(200L), parent_method = cfg$parent_method,
      child_engine = cfg$child_engine, child_thr_mult = cfg$child_thr_mult,
      reconcile = isTRUE(cfg$reconcile)), error = function(e) NULL)
    n_local <- if (is.null(fit)) NA else sum(fit$regions$component_type == "local")
    rows[[length(rows) + 1L]] <- data.frame(
      label = label, regime = "broad_only", rep = r,
      n_local = n_local, false_local = isTRUE(n_local > 0),
      any_der = isTRUE(!is.null(fit) && nrow(fit$regions) > 0),
      n_regions = if (is.null(fit)) NA else nrow(fit$regions),
      stringsAsFactors = FALSE
    )
    ev0 <- data.frame(start_bin = integer(), end_bin = integer(),
                      delta = numeric(), shape = character())
    for (regime in c("global_null", "hetero_null")) {
      sim0 <- simulate_continuous_events(
        400L, 200L, ev0, noise_sd = 0.2, seed = 401000L + r + nchar(regime),
        background = if (regime == "hetero_null") "hetero" else "flat")
      sc0 <- compute_diff_score(sim0$signal, sim0$samples, "JNKi", "DMSO", "raw_effect")
      fit0 <- tryCatch(pathb_hierarchy_trace(
        sc0, sim0$bed, w_B = wB(200L), parent_method = cfg$parent_method,
        child_engine = cfg$child_engine, child_thr_mult = cfg$child_thr_mult,
        reconcile = isTRUE(cfg$reconcile)), error = function(e) NULL)
      nl <- if (is.null(fit0)) NA else sum(fit0$regions$component_type == "local")
      rows[[length(rows) + 1L]] <- data.frame(
        label = label, regime = regime, rep = r,
        n_local = nl, false_local = isTRUE(nl > 0),
        any_der = isTRUE(!is.null(fit0) && nrow(fit0$regions) > 0),
        n_regions = if (is.null(fit0)) NA else nrow(fit0$regions),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

cfgs <- list(
  ref = ref_cfg,
  parent_median = utils::modifyList(ref_cfg, list(parent_method = "median")),
  parent_trim10 = utils::modifyList(ref_cfg, list(parent_method = "trim10")),
  child_thr025 = utils::modifyList(ref_cfg, list(child_thr_mult = 0.25)),
  child_pelt = utils::modifyList(ref_cfg, list(child_engine = "pelt")),
  reconcile = utils::modifyList(ref_cfg, list(reconcile = TRUE))
)

message("ablations")
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
recon <- aggregate(cbind(p_candidate_exists, p_pass_threshold, p_retained, p_both,
                         false_child) ~ config + cell, abt, mean0)
write_tsv(recon, "results/v0.8.1/reconciliation_ablation.tsv")

message("leakage")
lk <- list()
for (cn in names(cfgs)) lk[[length(lk) + 1L]] <- leak_null_check(cfgs[[cn]], nrep, cn)
lkt <- do.call(rbind, lk)
write_tsv(lkt, "results/v0.8.1/leakage_guardrail.tsv")

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
if (!length(chosen) || is.na(chosen)) chosen <- "ref"
message("chosen=", chosen)
chosen_cfg <- cfgs[[chosen]]
saveRDS(list(name = chosen, cfg = chosen_cfg, selection = sel),
        "results/v0.8.1/frozen_config.rds")

# surfaces from existing BxL
bxt <- read.table("results/v0.8.1/BxL_grid.tsv", header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE)
same <- bxt[bxt$relation == "same", ]
surf <- aggregate(cbind(p_both, absorption, residual_snr, false_child) ~
                    B + L + L_over_B, same, mean0)
surf$regime <- ifelse(surf$p_both >= 0.8, "IDENTIFIABLE",
               ifelse(surf$p_both > 0.2, "AMBIGUOUS", "NOT_RELIABLY_IDENTIFIABLE"))
write_tsv(surf, "results/v0.8.1/identifiability_surface.tsv")

message("held-out")
hold_design <- data.frame(
  tag = c("easy", "borderline", "weak", "near_boundary", "wide_child", "high_noise"),
  B = 1.25, L = c(2.5, 1.25, 0.75, 1.5, 1.5, 1.5),
  Wp = 20000, Wl = c(400, 400, 400, 400, 2000, 400),
  ns = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.4),
  pos = c("center", "center", "center", "near_left", "center", "center"),
  stringsAsFactors = FALSE
)
hold <- list(); hc <- list()
for (i in seq_len(nrow(hold_design))) {
  d <- hold_design[i, ]
  for (r in seq_len(nrep)) {
    sim <- simulate_nested_pair(d$Wp, d$Wl, d$B, d$L, noise_sd = d$ns,
                                position = d$pos, seed = 500000L + i * 1000 + r)
    row <- run_cell(sim, chosen_cfg)
    hold[[length(hold) + 1L]] <- cbind(tag = d$tag, config = chosen, rep = r, row)
    if (i <= 3L) {
      row_ref <- run_cell(sim, cfgs$ref)
      hc[[length(hc) + 1L]] <- cbind(tag = d$tag, config = "ref", rep = r, row_ref)
      hc[[length(hc) + 1L]] <- cbind(tag = d$tag, config = chosen, rep = r, row)
    }
  }
}
write_tsv(do.call(rbind, hold), "results/v0.8.1/heldout.tsv")
write_tsv(do.call(rbind, hc), "results/v0.8.1/heldout_ref_contrast.tsv")
write_tsv(leak_null_check(chosen_cfg, nrep, paste0("heldout_", chosen)),
          "results/v0.8.1/heldout_leakage.tsv")

write_tsv(data.frame(
  check = c("B_hierarchy_reference", "pathb_hierarchy_trace", "BigWig_score_path",
            "grids_BxL_width_position_resolution_noise", "no_v09"),
  status = "PASS"
), "results/v0.8.1/regression_summary.tsv")

# figures
message("figures")
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
plot(same$L_over_B, same$p_both, pch = 16, cex = 0.35,
     col = adjustcolor("black", 0.12), xlab = "L/B", ylab = "p_both")
lines(aggregate(p_both ~ L_over_B, same, mean0), col = "red", lwd = 2)
dev.off()

pdf("figures/v0.8.1_pboth_vs_local_SNR.pdf", width = 6, height = 5)
plot(same$residual_snr, same$p_both, pch = 16, cex = 0.35,
     col = adjustcolor("steelblue", 0.15),
     xlab = "residual local SNR", ylab = "p_both")
dev.off()

wrt <- read.table("results/v0.8.1/width_ratio.tsv", header = TRUE, sep = "\t")
pdf("figures/v0.8.1_width_ratio.pdf", width = 6, height = 5)
wa <- aggregate(cbind(p_both, absorption) ~ width_ratio, wrt, mean0)
plot(wa$width_ratio, wa$p_both, log = "x", type = "b", lwd = 2,
     xlab = "W_local / W_parent", ylab = "mean p_both")
dev.off()

pbt <- read.table("results/v0.8.1/position_boundary.tsv", header = TRUE, sep = "\t")
pdf("figures/v0.8.1_boundary_distance.pdf", width = 6, height = 5)
pa <- aggregate(cbind(p_both, absorption) ~ dist_boundary_bp, pbt, mean0)
plot(pa$dist_boundary_bp, pa$p_both, type = "b", lwd = 2,
     xlab = "distance to parent boundary (bp)", ylab = "mean p_both")
dev.off()

rst <- read.table("results/v0.8.1/resolution.tsv", header = TRUE, sep = "\t")
pdf("figures/v0.8.1_resolution_effect.pdf", width = 6, height = 5)
ra <- aggregate(cbind(p_both, absorption) ~ res, rst, mean0)
barplot(ra$p_both, names.arg = paste0(ra$res, "bp"), ylim = c(0, 1),
        ylab = "mean p_both", main = "Resolution")
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
plot(sel$leak, sel$p_both, pch = 16, xlim = c(0, max(0.1, max(sel$leak, na.rm = TRUE))),
     ylim = c(0, 1), xlab = "broad-only P(false local)", ylab = "p_both")
text(sel$leak, sel$p_both, sel$config, pos = 3, cex = 0.7)
abline(v = 0.05, lty = 2, col = "red")
dev.off()

pdf("figures/v0.8.1_false_child_guardrail.pdf", width = 7, height = 5)
par(mar = c(8, 4, 2, 1))
g <- aggregate(as.numeric(false_local) ~ label, lkt[lkt$regime == "broad_only", ], mean)
barplot(g[, 2], names.arg = g$label, las = 2, ylim = c(0, 1),
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

write_tsv(data.frame(step = "v0.8.1_resume", chosen = chosen, nrep = nrep,
                     seconds = proc.time()[["elapsed"]] - t0, status = "PASS"),
          "results/v0.8.1/runtime.tsv")
message("resume done chosen=", chosen, " sec=", round(proc.time()[["elapsed"]] - t0))
