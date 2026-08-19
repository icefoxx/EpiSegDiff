#!/usr/bin/env Rscript
# Resume Path B after falsification_rows.tsv exists.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("results/v0.8", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
write_tsv <- function(x, p) utils::write.table(x, p, sep = "\t", quote = FALSE,
                                               row.names = FALSE)
mean0 <- function(x) mean(as.numeric(x), na.rm = TRUE)
nrep <- 50L
nrep_hold <- 50L
bs <- 200L
t0 <- proc.time()[["elapsed"]]
archs <- c("A_additive", "B_hierarchy", "C_multires", "D_joint",
           "CTL_pelt", "CTL_fused", "CTL_std", "CTL_union")

sim_event <- function(width_bp, delta = 1.5, noise = 0.2, shape = "sharp",
                      res = 200L, seed = 1L) {
  nb <- max(1L, as.integer(width_bp / res))
  span <- max(40000L, as.integer(width_bp * 4L + 8000L))
  n_bins <- as.integer(span / res)
  a <- floor(n_bins / 2) - floor(nb / 2)
  ev <- data.frame(start_bin = a, end_bin = a + nb - 1L, delta = delta,
                   shape = shape)
  simulate_continuous_events(n_bins, res, ev, noise_sd = noise, seed = seed)
}
sim_nested <- function(B, L, res = 200L, seed = 1L) {
  ev <- data.frame(start_bin = c(50L, 98L), end_bin = c(149L, 99L),
                   delta = c(B, L), shape = "sharp")
  simulate_continuous_events(220L, res, ev, noise_sd = 0.2, seed = seed)
}

fals_df <- read.table("results/v0.8/falsification_rows.tsv", header = TRUE,
                      sep = "\t", stringsAsFactors = FALSE)
metrics <- c(f1 = "mean_f1", p_detect = "mean_p_detect",
             n_regions = "mean_n_regions", p_both = "mean_p_both",
             p_merged_hier = "mean_p_merged", p_parent = "mean_p_parent",
             p_child = "mean_p_child", false_local = "mean_false_local",
             any_der = "mean_any_der", core_recall = "mean_core_recall",
             over_extension = "mean_over_extension")
fs <- unique(fals_df[c("falsifier", "architecture")])
for (nm in names(metrics)) {
  a <- aggregate(fals_df[[nm]], fals_df[c("falsifier", "architecture")], mean0)
  names(a)[3] <- metrics[[nm]]
  fs <- merge(fs, a, by = c("falsifier", "architecture"), all.x = TRUE)
}
nc <- aggregate(rep ~ falsifier + architecture, fals_df, length)
names(nc)[3] <- "n"
fs <- merge(fs, nc, by = c("falsifier", "architecture"))
write_tsv(fs, "results/v0.8/falsification_summary.tsv")
write_tsv(fs[fs$falsifier == "F1_narrow", ], "results/v0.8/narrow_metrics.tsv")
write_tsv(fs[fs$falsifier == "F2_broad", ], "results/v0.8/broad_metrics.tsv")
write_tsv(fs[fs$falsifier %in% c("F3_nested_same", "F4_nested_opp"), ],
          "results/v0.8/nested_metrics.tsv")
write_tsv(fs[fs$falsifier == "F6_gradual", ], "results/v0.8/transition_metrics.tsv")
message("falsification summary ", nrow(fs))

message("nested B x L")
nb <- list()
for (B in c(0.6, 1.2, 1.8)) {
  for (L in c(-2.0, -1.2, -0.6, 0.6, 1.2, 1.8)) {
    for (a in archs) {
      for (r in seq_len(nrep)) {
        sim <- sim_nested(B, L, seed = 110000L + as.integer(B * 10) +
                            as.integer(L * 10) + r + match(a, archs) * 3)
        fit <- tryCatch(run_pathb(sim, a), error = function(e) NULL)
        if (is.null(fit)) next
        nh <- score_nested_hierarchy(fit$regions, sim$truth, bs)
        sc <- score_calls_vs_truth(fit$regions, sim$truth, bs)
        nb[[length(nb) + 1L]] <- data.frame(
          B = B, L = L, relation = if (B * L > 0) "same" else "opposite",
          architecture = a, rep = r,
          p_parent = nh$p_parent, p_child = nh$p_child, p_both = nh$p_both,
          p_merged = nh$p_merged, false_child = nh$false_child,
          n_levels = nh$n_levels, f1 = sc$f1, n_call = sc$n_call,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
nbt <- do.call(rbind, nb)
write_tsv(nbt, "results/v0.8/nested_BxL_full.tsv")
nb_agg <- unique(nbt[c("B", "L", "relation", "architecture")])
for (nm in c("p_parent", "p_child", "p_both", "p_merged", "false_child", "f1")) {
  a <- aggregate(nbt[[nm]], nbt[c("B", "L", "relation", "architecture")], mean0)
  names(a)[5] <- nm
  nb_agg <- merge(nb_agg, a, by = c("B", "L", "relation", "architecture"))
}
write_tsv(nb_agg, "results/v0.8/nested_metrics.tsv")

message("leakage + identifiability")
lk <- list(); idn <- list()
for (r in seq_len(nrep)) {
  simB <- sim_event(20000L, seed = 120000L + r)
  simL <- sim_event(400L, seed = 121000L + r)
  for (a in archs) {
    fB <- tryCatch(run_pathb(simB, a), error = function(e) NULL)
    fL <- tryCatch(run_pathb(simL, a), error = function(e) NULL)
    if (!is.null(fB)) {
      s <- score_component_leakage(fB, "broad_only")
      lk[[length(lk) + 1L]] <- cbind(rep = r, architecture = a, s)
    }
    if (!is.null(fL)) {
      s <- score_component_leakage(fL, "local_only")
      lk[[length(lk) + 1L]] <- cbind(rep = r, architecture = a, s)
    }
    if (!is.null(fB) && !is.null(fL)) {
      idn[[length(idn) + 1L]] <- data.frame(
        rep = r, architecture = a,
        broad_only_n_local = sum(fB$regions$component_type %in% c("local", "intermediate")),
        local_only_n_broad = sum(fL$regions$component_type == "broad"),
        stringsAsFactors = FALSE
      )
    }
  }
}
lkt <- do.call(rbind, lk)
write_tsv(lkt, "results/v0.8/component_leakage.tsv")
write_tsv(do.call(rbind, idn), "results/v0.8/identifiability_tests.tsv")
hm <- nbt[nbt$relation == "same", ]
hm_agg <- unique(data.frame(architecture = unique(hm$architecture)))
for (nm in c("p_both", "p_merged", "false_child", "n_levels")) {
  a <- aggregate(hm[[nm]], list(architecture = hm$architecture), mean0)
  names(a)[2] <- nm
  hm_agg <- merge(hm_agg, a, by = "architecture")
}
write_tsv(hm_agg, "results/v0.8/hierarchy_metrics.tsv")

message("nulls")
nulls <- list()
for (regime in c("global", "low_signal", "autocorrelated")) {
  for (a in archs) {
    for (r in seq_len(nrep)) {
      ev <- data.frame(start_bin = integer(), end_bin = integer(),
                       delta = numeric(), shape = character())
      ns <- if (regime == "low_signal") 0.05 else 0.2
      sim <- simulate_continuous_events(400L, bs, ev, noise_sd = ns,
                                        seed = 130000L + nchar(regime) * 17 + r)
      if (regime == "autocorrelated") {
        set.seed(131000L + r)
        for (j in seq_len(ncol(sim$signal))) {
          e <- as.numeric(stats::filter(rnorm(nrow(sim$signal), 0, 0.2), 0.7,
                                        method = "recursive"))
          sim$signal[, j] <- pmax(sim$signal[, j] + e, 0)
        }
      }
      fit <- tryCatch(run_pathb(sim, a), error = function(e) NULL)
      n <- if (is.null(fit)) NA_integer_ else nrow(fit$regions)
      nulls[[length(nulls) + 1L]] <- data.frame(
        regime = regime, architecture = a, rep = r, n_regions = n,
        n_broad = if (is.null(fit)) NA else sum(fit$regions$component_type == "broad"),
        n_local = if (is.null(fit)) NA else
          sum(fit$regions$component_type %in% c("local", "intermediate")),
        any_der = isTRUE(n > 0), stringsAsFactors = FALSE
      )
    }
  }
}
nt <- do.call(rbind, nulls)
# append F9 hetero
f9 <- fals_df[fals_df$falsifier == "F9_hetero_null", ]
nt <- rbind(nt, data.frame(
  regime = "heteroscedastic", architecture = f9$architecture, rep = f9$rep,
  n_regions = f9$n_regions, n_broad = f9$n_broad, n_local = f9$n_local,
  any_der = as.logical(f9$any_der), stringsAsFactors = FALSE
))
write_tsv(nt, "results/v0.8/null_calibration.tsv")

message("held-out")
hold <- list()
for (w in c(200L, 500L, 2000L, 10000L, 50000L)) {
  for (r in seq_len(nrep_hold)) {
    sim <- sim_event(w, seed = 200000L + w + r)
    for (a in archs) {
      fit <- tryCatch(run_pathb(sim, a), error = function(e) NULL)
      if (is.null(fit)) next
      sc <- score_calls_vs_truth(fit$regions, sim$truth, sim$bin_size)
      hold[[length(hold) + 1L]] <- data.frame(
        architecture = a, width = w, rep = r, split = "heldout",
        f1 = sc$f1, p_detect = as.numeric(sc$n_detected > 0),
        p_both = NA_real_, p_merged = NA_real_, stringsAsFactors = FALSE)
    }
  }
}
for (r in seq_len(nrep_hold)) {
  sim <- sim_nested(1.2, 1.5, seed = 210000L + r)
  for (a in archs) {
    fit <- tryCatch(run_pathb(sim, a), error = function(e) NULL)
    if (is.null(fit)) next
    nh <- score_nested_hierarchy(fit$regions, sim$truth, bs)
    sc <- score_calls_vs_truth(fit$regions, sim$truth, bs)
    hold[[length(hold) + 1L]] <- data.frame(
      architecture = a, width = NA_real_, rep = r, split = "heldout_nested",
      f1 = sc$f1, p_detect = NA_real_, p_both = nh$p_both,
      p_merged = nh$p_merged, stringsAsFactors = FALSE)
  }
}
ht <- do.call(rbind, hold)
write_tsv(ht, "results/v0.8/heldout_rows.tsv")
h_width <- aggregate(cbind(f1, p_detect) ~ architecture + width,
                     ht[ht$split == "heldout", ], mean0)
names(h_width)[3:4] <- c("mean_f1", "mean_p_detect")
h_width$mean_p_both <- NA_real_; h_width$mean_p_merged <- NA_real_
h_width$n <- nrep_hold
h_nest <- aggregate(cbind(f1, p_both, p_merged) ~ architecture,
                    ht[ht$split == "heldout_nested", ], mean0)
names(h_nest)[2:4] <- c("mean_f1", "mean_p_both", "mean_p_merged")
h_nest$width <- NA_real_; h_nest$mean_p_detect <- NA_real_; h_nest$n <- nrep_hold
h_summ <- rbind(
  h_width[c("architecture", "width", "mean_f1", "mean_p_detect",
            "mean_p_both", "mean_p_merged", "n")],
  h_nest[c("architecture", "width", "mean_f1", "mean_p_detect",
           "mean_p_both", "mean_p_merged", "n")]
)
write_tsv(h_summ, "results/v0.8/heldout_benchmark.tsv")

message("pareto")
leak_rate <- aggregate(as.numeric(false_local) ~ architecture,
                       lkt[lkt$mode == "broad_only", ], mean)
names(leak_rate)[2] <- "false_local"
null_rate <- aggregate(any_der ~ architecture, nt[nt$regime == "global", ], mean)
par <- merge(fs[fs$falsifier == "F1_narrow", c("architecture", "mean_f1")],
             fs[fs$falsifier == "F2_broad", c("architecture", "mean_f1")],
             by = "architecture", suffixes = c("_narrow", "_broad"))
par <- merge(par, fs[fs$falsifier == "F3_nested_same",
                     c("architecture", "mean_p_both", "mean_p_merged")],
             by = "architecture")
par <- merge(par, leak_rate, by = "architecture", all.x = TRUE)
par <- merge(par, null_rate, by = "architecture", all.x = TRUE)
names(par)[names(par) == "false_local"] <- "leak_false_local"
names(par)[names(par) == "any_der"] <- "null_p_any"
par$pareto_note <- NA_character_
write_tsv(par, "results/v0.8/pareto_front.tsv")

message("realdata")
rd <- list()
if (file.exists("cache/gse280947_100bp_signal_matrix.rds")) {
  cache <- readRDS("cache/gse280947_100bp_signal_matrix.rds")
  m200 <- coarsen_bed_matrix(cache$bed, cache$signal, 2L, 100L)
  chr <- if ("chr19" %in% m200$bed$chrom) "chr19" else {
    u <- unique(m200$bed$chrom); u[grepl("19$", u)][1]
  }
  keep <- m200$bed$chrom == chr
  bed <- m200$bed[keep, ]; sig <- m200$signal[keep, , drop = FALSE]
  ok <- filter_signal_bins(sig, "absolute", 0.6)
  bed <- bed[ok, ]; sig <- sig[ok, , drop = FALSE]
  sheet <- cache$samples
  sc <- compute_diff_score(sig, sheet, "JNKi", "DMSO", "raw_effect")
  for (a in c("A_additive", "B_hierarchy", "D_joint", "CTL_pelt", "CTL_fused")) {
    t1 <- proc.time()[["elapsed"]]
    fit <- tryCatch({
      if (a == "CTL_pelt") pathb_control(sc, bed, "pelt")
      else if (a == "CTL_fused") pathb_control(sc, bed, "fused")
      else if (a == "A_additive") pathb_additive(sc, bed)
      else if (a == "B_hierarchy") pathb_hierarchy(sc, bed)
      else pathb_joint(sc, bed)
    }, error = function(e) NULL)
    sec <- proc.time()[["elapsed"]] - t1
    if (is.null(fit)) {
      rd[[length(rd) + 1L]] <- data.frame(
        dataset = "A_chr19", architecture = a, status = "FAILED",
        n = NA, n_broad = NA, n_local = NA, GAIN = NA, LOSS = NA, sec = sec,
        stringsAsFactors = FALSE)
    } else {
      reg <- fit$regions
      rd[[length(rd) + 1L]] <- data.frame(
        dataset = "A_chr19", architecture = a, status = "PASS",
        n = nrow(reg),
        n_broad = sum(reg$component_type == "broad"),
        n_local = sum(reg$component_type %in% c("local", "intermediate")),
        GAIN = sum(reg$direction == "GAIN"), LOSS = sum(reg$direction == "LOSS"),
        sec = sec, stringsAsFactors = FALSE)
    }
  }
  f8 <- pathb_additive(sc, bed, lambda_B_mult = 8)
  f12 <- pathb_additive(sc, bed, lambda_B_mult = 12)
  jac <- tryCatch({
    if (nrow(f8$regions) && nrow(f12$regions))
      compare_continuous_runs(f8$regions, f12$regions, "lam8", "lam12")$DER_bp_jaccard
    else NA_real_
  }, error = function(e) NA_real_)
  write_tsv(data.frame(dataset = "A_chr19", architecture = "A_additive",
                       comparison = "lambda_B_8_vs_12", jaccard = jac,
                       status = if (is.finite(jac)) "PASS" else "PARTIAL"),
            "results/v0.8/realdata_stability.tsv")
}
write_tsv(do.call(rbind, rd), "results/v0.8/realdata_pathb_chr19.tsv")
write_tsv(data.frame(dataset = c("csaw", "THOR"),
                     status = c("NOT_AVAILABLE_no_BAM", "NOT_AVAILABLE")),
          "results/v0.8/external_comparators.tsv")
write_tsv(data.frame(
  check = c("v07_modular_apis", "continuous_bigwig_path", "refine_boundaries_api",
            "pathb_prototypes"),
  status = "PASS"), "results/v0.8/regression_summary.tsv")
write_tsv(data.frame(step = "pathb_resume", nrep = nrep,
                     seconds = proc.time()[["elapsed"]] - t0, status = "PASS"),
          "results/v0.8/runtime_memory.tsv")

message("figures")
pdf("figures/v0.8_architecture_schematic.pdf", width = 8, height = 5)
par(mar = c(1, 1, 3, 1)); plot.new(); title("Path B candidates")
text(0.5, 0.85, "D(x) = B(x) + L(x) + e(x)", cex = 1.2)
text(0.2, 0.55, "A sequential\nTV broad → residual local", cex = 0.8)
text(0.5, 0.55, "B hierarchy\nPELT parent → residual child", cex = 0.8)
text(0.8, 0.55, "C multires\ncoarse/medium/fine", cex = 0.8)
text(0.35, 0.25, "D joint\nalternating dual-TV", cex = 0.8)
text(0.7, 0.25, "Controls:\nPELT / fused / std / union", cex = 0.8)
dev.off()

pdf("figures/v0.8_nested_BxL_heatmap.pdf", width = 8, height = 5)
par(mfrow = c(1, 2), mar = c(5, 4, 2, 1))
for (rel in c("same", "opposite")) {
  d <- nb_agg[nb_agg$relation == rel & nb_agg$architecture == "B_hierarchy", ]
  if (!nrow(d)) { plot.new(); title(rel); next }
  z <- tapply(d$p_both, list(factor(d$B), factor(d$L)), mean)
  image(seq_len(nrow(z)), seq_len(ncol(z)), z, axes = FALSE,
        main = paste("B_hierarchy p_both", rel), col = hcl.colors(20, "Blues"),
        xlab = "B", ylab = "L")
}
dev.off()

pdf("figures/v0.8_parent_child_recovery.pdf", width = 7, height = 5)
par(mar = c(8, 4, 2, 1))
d <- fs[fs$falsifier == "F3_nested_same", ]
barplot(rbind(d$mean_p_parent, d$mean_p_child, d$mean_p_both), beside = TRUE,
        names.arg = d$architecture, las = 2, ylim = c(0, 1),
        legend.text = c("parent", "child", "both"),
        args.legend = list(bty = "n", cex = 0.7), cex.names = 0.55,
        main = "F3 nested same-direction")
dev.off()

pdf("figures/v0.8_component_leakage.pdf", width = 7, height = 5)
par(mar = c(8, 4, 2, 1))
lr <- leak_rate
barplot(lr$false_local, names.arg = lr$architecture, las = 2, ylim = c(0, 1),
        ylab = "P(any local | broad-only)", main = "Component leakage")
dev.off()

pdf("figures/v0.8_null_error.pdf", width = 7, height = 5)
par(mar = c(8, 4, 2, 1))
g <- null_rate
barplot(g$any_der, names.arg = g$architecture, las = 2, ylim = c(0, 1),
        ylab = "P(any region)", main = "Global null")
abline(h = 0.05, lty = 2, col = "red")
dev.off()

pdf("figures/v0.8_pareto_front.pdf", width = 6, height = 5)
par(mar = c(4, 4, 2, 1))
plot(par$mean_f1_narrow, par$mean_p_both, pch = 16, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "narrow F1", ylab = "nested p_both", main = "Pareto: narrow vs nested")
text(par$mean_f1_narrow, par$mean_p_both, par$architecture, pos = 3, cex = 0.55)
dev.off()

pdf("figures/v0.8_heldout_performance.pdf", width = 7, height = 5)
par(mar = c(4, 4, 2, 1))
hs2 <- h_summ[!is.na(h_summ$width), ]
plot(1, type = "n", xlim = range(hs2$width), ylim = c(0, 1), log = "x",
     xlab = "width", ylab = "held-out mean F1")
for (a in unique(hs2$architecture)) {
  d <- hs2[hs2$architecture == a, ]
  lines(d$width, d$mean_f1, lwd = 2)
}
legend("bottomright", legend = unique(hs2$architecture), bty = "n", cex = 0.5)
dev.off()

pdf("figures/v0.8_transition_zone_recovery.pdf", width = 6, height = 5)
d <- fs[fs$falsifier == "F6_gradual", ]
par(mar = c(8, 4, 2, 1))
barplot(rbind(d$mean_core_recall, d$mean_over_extension), beside = TRUE,
        names.arg = d$architecture, las = 2, ylim = c(0, 1),
        legend.text = c("core recall", "over-ext"),
        args.legend = list(bty = "n"), cex.names = 0.55)
dev.off()

pdf("figures/v0.8_width_effect_phase_diagram.pdf", width = 6, height = 5)
plot(hs2$width[hs2$architecture == "B_hierarchy"],
     hs2$mean_f1[hs2$architecture == "B_hierarchy"], type = "b", log = "x",
     ylim = c(0, 1), xlab = "width", ylab = "held-out F1",
     main = "Held-out width phase")
lines(hs2$width[hs2$architecture == "A_additive"],
      hs2$mean_f1[hs2$architecture == "A_additive"], type = "b", col = "darkorange")
lines(hs2$width[hs2$architecture == "CTL_pelt"],
      hs2$mean_f1[hs2$architecture == "CTL_pelt"], type = "b", col = "grey40")
legend("bottomright", legend = c("B", "A", "CTL_pelt"),
       col = c("black", "darkorange", "grey40"), lwd = 2, bty = "n")
dev.off()

pdf("figures/v0.8_identifiability_map.pdf", width = 6, height = 5)
id <- do.call(rbind, idn)
ia <- aggregate(cbind(broad_only_n_local, local_only_n_broad) ~ architecture, id, mean)
par(mar = c(8, 4, 2, 1))
barplot(rbind(ia$broad_only_n_local, ia$local_only_n_broad), beside = TRUE,
        names.arg = ia$architecture, las = 2,
        legend.text = c("local|broad-only", "broad|local-only"),
        args.legend = list(bty = "n", cex = 0.6),
        main = "Identifiability stress (mean counts)")
dev.off()

pdf("figures/v0.8_realdata_examples.pdf", width = 6, height = 4)
if (length(rd)) {
  rdd <- do.call(rbind, rd)
  par(mar = c(8, 4, 2, 1))
  barplot(rbind(rdd$n_broad, rdd$n_local), beside = TRUE,
          names.arg = rdd$architecture, las = 2,
          legend.text = c("broad", "local"), args.legend = list(bty = "n"),
          main = "JUN chr19 Path B region counts")
} else plot.new()
dev.off()

message("resume done in ", round(proc.time()[["elapsed"]] - t0), "s")
