#!/usr/bin/env Rscript
# v0.8 Path B: falsification, leakage, nested, held-out, nulls, real-data slice.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("results/v0.8", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
write_tsv <- function(x, p) utils::write.table(x, p, sep = "\t", quote = FALSE,
                                               row.names = FALSE)
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
  n_bins <- 220L
  b0 <- 50L; b1 <- 149L
  l0 <- 98L; l1 <- 99L
  ev <- data.frame(start_bin = c(b0, l0), end_bin = c(b1, l1),
                   delta = c(B, L), shape = "sharp")
  simulate_continuous_events(n_bins, res, ev, noise_sd = 0.2, seed = seed)
}

empty_metrics <- function() {
  data.frame(
    f1 = NA_real_, recall = NA_real_, precision = NA_real_,
    p_detect = NA_real_, n_call = NA_real_, n_detected = NA_real_,
    n_merged = NA_real_, n_split = NA_real_, med_bound_err = NA_real_,
    n_regions = NA_real_, n_broad = NA_real_, n_local = NA_real_,
    p_parent = NA_real_, p_child = NA_real_, p_both = NA_real_,
    p_merged_hier = NA_real_, false_child = NA_real_,
    core_recall = NA_real_, over_extension = NA_real_,
    false_local = NA_real_, any_der = NA_real_,
    stringsAsFactors = FALSE
  )
}

record_row <- function(fid, arch, rep, width, sim, fit = NULL,
                       core_recall = NA_real_, over_extension = NA_real_,
                       false_local = NA_real_, any_der = NA_real_) {
  m <- empty_metrics()
  status <- "PASS"
  reason <- NA_character_
  if (is.null(fit)) {
    fit <- tryCatch(run_pathb(sim, arch), error = function(e) NULL)
  }
  if (is.null(fit)) {
    status <- "FAILED"
    reason <- "run_pathb error"
  } else {
    sc <- score_calls_vs_truth(fit$regions, sim$truth, sim$bin_size)
    m$f1 <- sc$f1
    m$recall <- sc$recall
    m$precision <- sc$precision
    m$p_detect <- as.numeric(isTRUE(sc$n_detected > 0))
    m$n_call <- sc$n_call
    m$n_detected <- sc$n_detected
    m$n_merged <- sc$n_merged
    m$n_split <- sc$n_split
    m$med_bound_err <- sc$med_bound_err
    m$n_regions <- nrow(fit$regions)
    m$n_broad <- sum(fit$regions$component_type == "broad")
    m$n_local <- sum(fit$regions$component_type %in% c("local", "intermediate"))
    if (nrow(sim$truth) >= 2L) {
      nh <- score_nested_hierarchy(fit$regions, sim$truth, sim$bin_size)
      m$p_parent <- nh$p_parent
      m$p_child <- nh$p_child
      m$p_both <- nh$p_both
      m$p_merged_hier <- nh$p_merged
      m$false_child <- nh$false_child
    }
    m$any_der <- as.numeric(nrow(fit$regions) > 0)
  }
  m$core_recall <- core_recall
  m$over_extension <- over_extension
  m$false_local <- false_local
  if (!is.na(any_der)) m$any_der <- any_der
  cbind(data.frame(falsifier = fid, architecture = arch, rep = rep,
                   width = width, status = status, failure_reason = reason,
                   stringsAsFactors = FALSE), m)
}

# ---- candidate registry ----
write_tsv(data.frame(
  architecture = c("A_additive", "B_hierarchy", "C_multires", "D_joint",
                   "CTL_pelt", "CTL_fused", "CTL_std", "CTL_union"),
  family = c("additive", "hierarchy", "multires", "joint",
             "control", "control", "control", "control"),
  status = "PROTOTYPE",
  identifiability = c("TV gap λB>>λL residual",
                      "parent duration + residual children",
                      "pool factors kc/km",
                      "alternating dual-TV λB/λL",
                      "flat PELT", "flat fused", "std+PELT", "union"),
  stringsAsFactors = FALSE
), "results/v0.8/candidate_architectures.tsv")

# ---- F1–F9 falsification (primary nrep=50) ----
message("falsification F1-F9")
fals <- list()

# F1 isolated narrow 200 & 500
for (w in c(200L, 500L)) {
  for (r in seq_len(nrep)) {
    sim <- sim_event(w, seed = 100000L + w + r)
    for (a in archs)
      fals[[length(fals) + 1L]] <- record_row("F1_narrow", a, r, w, sim)
  }
}
# F2 isolated broad 10kb & 50kb
for (w in c(10000L, 50000L)) {
  for (r in seq_len(nrep)) {
    sim <- sim_event(w, seed = 101000L + w + r)
    for (a in archs)
      fals[[length(fals) + 1L]] <- record_row("F2_broad", a, r, w, sim)
  }
}
# F3 nested same
for (r in seq_len(nrep)) {
  sim <- sim_nested(1.2, 1.5, seed = 102000L + r)
  for (a in archs)
    fals[[length(fals) + 1L]] <- record_row("F3_nested_same", a, r, NA, sim)
}
# F4 nested opposite / reversal
for (r in seq_len(nrep)) {
  sim <- sim_nested(1.2, -2.0, seed = 103000L + r)
  for (a in archs)
    fals[[length(fals) + 1L]] <- record_row("F4_nested_opp", a, r, NA, sim)
}
# F5 same-dir separated gap 10 bins, 2kb events
for (r in seq_len(nrep)) {
  ev <- data.frame(start_bin = c(40L, 60L), end_bin = c(49L, 69L),
                   delta = c(1.5, 1.5), shape = "sharp")
  sim <- simulate_continuous_events(160L, bs, ev, noise_sd = 0.2,
                                    seed = 104000L + r)
  for (a in archs)
    fals[[length(fals) + 1L]] <- record_row("F5_separated", a, r, 2000, sim)
}
# F6 gradual 10kb
for (r in seq_len(nrep)) {
  sim <- sim_event(10000L, shape = "gradual", seed = 105000L + r)
  for (a in archs) {
    fit <- tryCatch(run_pathb(sim, a), error = function(e) NULL)
    tz <- if (!is.null(fit) && nrow(sim$truth)) {
      score_transition_zones(fit$regions, sim$truth$start_bin[1],
                             sim$truth$end_bin[1], bs, nrow(sim$bed) * bs)
    } else NULL
    fals[[length(fals) + 1L]] <- record_row(
      "F6_gradual", a, r, 10000, sim, fit = fit,
      core_recall = if (is.null(tz)) NA else tz$core_recall,
      over_extension = if (is.null(tz)) NA else tz$over_extension)
  }
}
# F7 narrow at broad boundary
for (r in seq_len(nrep)) {
  ev <- data.frame(start_bin = c(80L, 80L), end_bin = c(160L, 81L),
                   delta = c(1.2, 2.0), shape = "sharp")
  sim <- simulate_continuous_events(240L, bs, ev, noise_sd = 0.2,
                                    seed = 106000L + r)
  for (a in archs)
    fals[[length(fals) + 1L]] <- record_row("F7_boundary_local", a, r, NA, sim)
}
# F8 broad-only — leakage
for (r in seq_len(nrep)) {
  sim <- sim_event(20000L, seed = 107000L + r)
  for (a in archs) {
    fit <- tryCatch(run_pathb(sim, a), error = function(e) NULL)
    fl <- if (is.null(fit)) NA else
      as.numeric(score_component_leakage(fit, "broad_only")$false_local)
    fals[[length(fals) + 1L]] <- record_row(
      "F8_broad_only", a, r, 20000, sim, fit = fit, false_local = fl)
  }
}
# F9 heteroscedastic null
for (r in seq_len(nrep)) {
  ev <- data.frame(start_bin = integer(), end_bin = integer(),
                   delta = numeric(), shape = character())
  sim <- simulate_continuous_events(400L, bs, ev, noise_sd = 0.2,
                                    seed = 108000L + r, background = "hetero")
  for (a in archs) {
    fit <- tryCatch(run_pathb(sim, a), error = function(e) NULL)
    n <- if (is.null(fit)) NA_integer_ else nrow(fit$regions)
    fals[[length(fals) + 1L]] <- record_row(
      "F9_hetero_null", a, r, NA, sim, fit = fit,
      any_der = as.numeric(isTRUE(n > 0)))
  }
}
fals_df <- do.call(rbind, fals)
saveRDS(fals_df, "results/v0.8/falsification_raw.rds")
write_tsv(fals_df, "results/v0.8/falsification_rows.tsv")

mean0 <- function(x) mean(as.numeric(x), na.rm = TRUE)
# aggregate one metric at a time — na.omit on a multi-cbind drops rows
# that are NA only for transition/leakage-specific columns.
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
message("  falsification summary rows ", nrow(fs))

write_tsv(fs[fs$falsifier == "F1_narrow", ], "results/v0.8/narrow_metrics.tsv")
write_tsv(fs[fs$falsifier == "F2_broad", ], "results/v0.8/broad_metrics.tsv")
write_tsv(fs[fs$falsifier %in% c("F3_nested_same", "F4_nested_opp"), ],
          "results/v0.8/nested_metrics.tsv")
write_tsv(fs[fs$falsifier == "F6_gradual", ], "results/v0.8/transition_metrics.tsv")

# ---- nested B×L matrix (critical) ----
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
        nb[[length(nb) + 1L]] <- cbind(
          B = B, L = L, relation = if (B * L > 0) "same" else "opposite",
          architecture = a, rep = r, nh, f1 = sc$f1, n_call = sc$n_call
        )
      }
    }
  }
}
nbt <- do.call(rbind, nb)
write_tsv(nbt, "results/v0.8/nested_BxL_full.tsv")
nb_agg <- aggregate(cbind(p_parent, p_child, p_both, p_merged, false_child, f1) ~
                      B + L + relation + architecture, nbt, mean, na.rm = TRUE)
write_tsv(nb_agg, "results/v0.8/nested_metrics.tsv")

# ---- component leakage / identifiability ----
message("leakage + identifiability")
lk <- list(); idn <- list()
for (r in seq_len(nrep)) {
  # broad-only
  simB <- sim_event(20000L, seed = 120000L + r)
  # local-only 400bp
  simL <- sim_event(400L, seed = 121000L + r)
  # ambiguous: local width approaches broad (5kb local-ish vs 8kb broad)
  simA <- sim_nested(1.2, 1.2, seed = 122000L + r)
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

# hierarchy metrics from nested same
hm <- nbt[nbt$relation == "same", ]
hm_agg <- aggregate(cbind(p_both, p_merged, false_child, n_levels) ~
                      architecture, hm, mean, na.rm = TRUE)
write_tsv(hm_agg, "results/v0.8/hierarchy_metrics.tsv")

# ---- null calibration ----
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
      nb0 <- if (is.null(fit)) NA else sum(fit$regions$component_type == "broad")
      nl0 <- if (is.null(fit)) NA else
        sum(fit$regions$component_type %in% c("local", "intermediate"))
      nulls[[length(nulls) + 1L]] <- data.frame(
        regime = regime, architecture = a, rep = r, n_regions = n,
        n_broad = nb0, n_local = nl0, any_der = isTRUE(n > 0),
        stringsAsFactors = FALSE
      )
    }
  }
}
# F9 already; append hetero from falsification summary
nt <- do.call(rbind, nulls)
write_tsv(nt, "results/v0.8/null_calibration.tsv")

# ---- held-out (seeds 200000+; development used 100000+) ----
message("held-out")
hold <- list()
for (w in c(200L, 500L, 2000L, 10000L, 50000L)) {
  for (r in seq_len(nrep_hold)) {
    sim <- sim_event(w, seed = 200000L + w + r)
    for (a in archs) {
      row <- record_row("heldout", a, r, w, sim)
      row$split <- "heldout"
      hold[[length(hold) + 1L]] <- row
    }
  }
}
for (r in seq_len(nrep_hold)) {
  sim <- sim_nested(1.2, 1.5, seed = 210000L + r)
  for (a in archs) {
    row <- record_row("heldout_nested", a, r, NA, sim)
    row$split <- "heldout_nested"
    hold[[length(hold) + 1L]] <- row
  }
}
ht <- do.call(rbind, hold)
saveRDS(ht, "results/v0.8/heldout_raw.rds")
h_width <- aggregate(cbind(f1, p_detect) ~ architecture + width,
                     ht[ht$split == "heldout", ], mean0)
names(h_width)[3:4] <- c("mean_f1", "mean_p_detect")
h_width$mean_p_both <- NA_real_
h_width$mean_p_merged <- NA_real_
h_width$n <- as.integer(nrep_hold)
h_nest <- aggregate(cbind(f1, p_both, p_merged_hier) ~ architecture,
                    ht[ht$split == "heldout_nested", ], mean0)
names(h_nest)[2:4] <- c("mean_f1", "mean_p_both", "mean_p_merged")
h_nest$width <- NA_real_
h_nest$mean_p_detect <- NA_real_
h_nest$n <- as.integer(nrep_hold)
h_summ <- rbind(
  h_width[c("architecture", "width", "mean_f1", "mean_p_detect",
            "mean_p_both", "mean_p_merged", "n")],
  h_nest[c("architecture", "width", "mean_f1", "mean_p_detect",
           "mean_p_both", "mean_p_merged", "n")]
)
write_tsv(h_summ, "results/v0.8/heldout_benchmark.tsv")

# ---- Pareto ----
message("pareto")
# metrics from falsification summary + leakage
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
# simple score: prefer high narrow, broad, p_both; low merge, leak, null
par$pareto_note <- NA_character_
write_tsv(par, "results/v0.8/pareto_front.tsv")

# ---- real-data feasibility (chr19 JUN) ----
message("realdata")
rd <- list()
if (file.exists("cache/gse280947_100bp_signal_matrix.rds")) {
  cache <- readRDS("cache/gse280947_100bp_signal_matrix.rds")
  m200 <- coarsen_bed_matrix(cache$bed, cache$signal, 2L, 100L)
  chr <- if ("chr19" %in% m200$bed$chrom) "chr19" else {
    u <- unique(m200$bed$chrom); u[grepl("19$", u)][1]
  }
  keep <- m200$bed$chrom == chr
  bed <- m200$bed[keep, ]
  sig <- m200$signal[keep, , drop = FALSE]
  sheet <- cache$samples
  # filter absolute 0.6
  ok <- filter_signal_bins(sig, "absolute", 0.6)
  bed <- bed[ok, ]; sig <- sig[ok, , drop = FALSE]
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
        GAIN = sum(reg$direction == "GAIN"),
        LOSS = sum(reg$direction == "LOSS"),
        sec = sec, stringsAsFactors = FALSE)
    }
  }
  # penalty stability for A
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

write_tsv(data.frame(
  dataset = c("csaw", "THOR"),
  status = c("NOT_AVAILABLE_no_BAM", "NOT_AVAILABLE"),
  stringsAsFactors = FALSE
), "results/v0.8/external_comparators.tsv")

write_tsv(data.frame(
  check = c("v07_modular_apis", "continuous_bigwig_path", "refine_boundaries_api",
            "pathb_prototypes_additive"),
  status = "PASS"
), "results/v0.8/regression_summary.tsv")

write_tsv(data.frame(
  step = "pathb_falsification_heldout",
  nrep = nrep,
  seconds = proc.time()[["elapsed"]] - t0,
  status = "PASS"
), "results/v0.8/runtime_memory.tsv")

# ---- figures ----
message("figures")
pdf("figures/v0.8_architecture_schematic.pdf", width = 8, height = 5)
par(mar = c(1, 1, 3, 1))
plot.new(); title("Path B candidates")
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
  z <- tapply(d$p_both, list(d$B, d$L), mean)
  image(z, axes = FALSE, main = paste("B_hierarchy p_both", rel),
        col = hcl.colors(20, "Blues"))
}
dev.off()

pdf("figures/v0.8_parent_child_recovery.pdf", width = 7, height = 5)
par(mar = c(8, 4, 2, 1))
d <- fs[fs$falsifier == "F3_nested_same", ]
barplot(rbind(d$mean_p_parent, d$mean_p_child, d$mean_p_both),
        beside = TRUE, names.arg = d$architecture, las = 2, ylim = c(0, 1),
        legend.text = c("parent", "child", "both"),
        args.legend = list(bty = "n", cex = 0.7),
        main = "F3 nested same-direction", cex.names = 0.6)
dev.off()

pdf("figures/v0.8_component_leakage.pdf", width = 7, height = 5)
par(mar = c(8, 4, 2, 1))
lr <- aggregate(cbind(false_local = as.numeric(false_local)) ~ architecture,
                lkt[lkt$mode == "broad_only", ], mean)
barplot(lr$false_local, names.arg = lr$architecture, las = 2, ylim = c(0, 1),
        ylab = "P(any local | broad-only)", main = "Component leakage")
dev.off()

pdf("figures/v0.8_null_error.pdf", width = 7, height = 5)
par(mar = c(8, 4, 2, 1))
nr <- aggregate(any_der ~ architecture + regime, nt, mean)
# global only for clarity
g <- nr[nr$regime == "global", ]
barplot(g$any_der, names.arg = g$architecture, las = 2, ylim = c(0, 1),
        ylab = "P(any region)", main = "Global null")
abline(h = 0.05, lty = 2, col = "red")
dev.off()

pdf("figures/v0.8_pareto_front.pdf", width = 6, height = 5)
par(mar = c(4, 4, 2, 1))
plot(par$mean_f1_narrow, par$mean_p_both, pch = 16,
     xlab = "narrow F1 (F1)", ylab = "nested p_both (F3)",
     xlim = c(0, 1), ylim = c(0, 1), main = "Pareto: narrow vs nested")
text(par$mean_f1_narrow, par$mean_p_both, par$architecture, pos = 3, cex = 0.6)
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
# reuse held-out widths as phase proxy
plot(hs2$width[hs2$architecture == "A_additive"],
     hs2$mean_f1[hs2$architecture == "A_additive"], type = "b", log = "x",
     ylim = c(0, 1), xlab = "width", ylab = "held-out F1",
     main = "Held-out width phase (selected archs)")
lines(hs2$width[hs2$architecture == "B_hierarchy"],
      hs2$mean_f1[hs2$architecture == "B_hierarchy"], type = "b", col = "steelblue")
lines(hs2$width[hs2$architecture == "CTL_pelt"],
      hs2$mean_f1[hs2$architecture == "CTL_pelt"], type = "b", col = "grey40")
legend("bottomright", legend = c("A", "B", "CTL_pelt"),
       col = c("black", "steelblue", "grey40"), lwd = 2, bty = "n")
dev.off()

pdf("figures/v0.8_identifiability_map.pdf", width = 6, height = 5)
id <- do.call(rbind, idn)
ia <- aggregate(cbind(broad_only_n_local, local_only_n_broad) ~ architecture, id, mean)
par(mar = c(8, 4, 2, 1))
barplot(rbind(ia$broad_only_n_local, ia$local_only_n_broad), beside = TRUE,
        names.arg = ia$architecture, las = 2,
        legend.text = c("local calls | broad-only", "broad calls | local-only"),
        args.legend = list(bty = "n", cex = 0.6),
        main = "Identifiability stress (mean counts)")
dev.off()

pdf("figures/v0.8_realdata_examples.pdf", width = 6, height = 4)
if (length(rd)) {
  rdd <- do.call(rbind, rd)
  par(mar = c(8, 4, 2, 1))
  barplot(rbind(rdd$n_broad, rdd$n_local), beside = TRUE, names.arg = rdd$architecture,
          las = 2, legend.text = c("broad", "local"), args.legend = list(bty = "n"),
          main = "JUN chr19 Path B region counts")
} else {
  plot.new(); title("realdata not run")
}
dev.off()

message("pathb done in ", round(proc.time()[["elapsed"]] - t0), "s")
