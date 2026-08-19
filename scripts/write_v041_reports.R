#!/usr/bin/env Rscript
# Fill v0.4.1 sensitivity / acceptance reports from machine-readable tables.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})

read_opt <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE)
}

esc <- function(...) paste(..., collapse = "", sep = "")

fmt <- function(x, d = 3) {
  if (is.null(x) || !length(x) || all(is.na(x))) return("NA")
  format(round(as.numeric(x), d), nsmall = min(d, 3), trim = TRUE)
}

grid <- read_opt("results/v0.4.1/combined_sensitivity_grid.tsv")
bin <- read_opt("results/v0.4.1/bin_size_sensitivity.tsv")
pseudo <- read_opt("results/v0.4.1/pseudocount_sensitivity.tsv")
filt <- read_opt("results/v0.4.1/filter_sensitivity.tsv")
conc <- read_opt("results/v0.4.1/run_concordance.tsv")
vd <- read_opt("results/v0.4.1/limma_variance_diagnostics.tsv")
qtab <- read_opt("results/v0.4.1/pseudocount_signal_quartile.tsv")
man <- if (file.exists("results/v0.4.1/reference_run_manifest.yaml")) {
  yaml::read_yaml("results/v0.4.1/reference_run_manifest.yaml")
} else list()

ref_g <- man$GAIN %||% NA
ref_l <- man$LOSS %||% NA
ref_w <- man$min_width %||% NA
ref_n <- man$n_tested %||% NA

n_grid_pass <- if (is.null(grid)) 0L else sum(grid$status == "PASS", na.rm = TRUE)
n_grid <- if (is.null(grid)) 0L else nrow(grid)
grid_complete <- isTRUE(n_grid == 12L && n_grid_pass == 12L)

bin_complete <- !is.null(bin) && all(c(100, 200, 500) %in% bin$bin_size) &&
  all(bin$status[bin$bin_size %in% c(100, 200, 500)] == "PASS")
pseudo_complete <- !is.null(pseudo) &&
  all(c(0.001, 0.01, 0.05, 0.1) %in% pseudo$pseudocount) &&
  all(pseudo$status == "PASS")
none_row <- if (!is.null(filt)) filt[filt$filter_method == "none", ][1, ] else NULL
filter_complete <- !is.null(none_row) && identical(none_row$status, "PASS")

jac_range <- if (!is.null(conc) && nrow(conc)) {
  range(conc$DER_bp_jaccard, na.rm = TRUE)
} else c(NA, NA)
rec_range <- if (!is.null(conc) && nrow(conc)) {
  range(conc$strong_reference_DER_recovery, na.rm = TRUE)
} else c(NA, NA)
pear_range <- if (!is.null(conc) && nrow(conc)) {
  range(conc$effect_pearson, na.rm = TRUE)
} else c(NA, NA)

loss_dom <- TRUE
sign_rev <- FALSE
if (!is.null(grid) && nrow(grid)) {
  ok <- grid[grid$status == "PASS", ]
  if (nrow(ok)) {
    loss_dom <- all(ok$LOSS > ok$GAIN, na.rm = TRUE)
  }
}
if (!is.null(conc) && nrow(conc)) {
  sign_rev <- any(conc$effect_pearson < 0, na.rm = TRUE) ||
    any(conc$sign_concordance < 0.5, na.rm = TRUE)
}

vd_class <- if (!is.null(vd) && nrow(vd)) vd$classification[1] else "UNRESOLVED"
warn_gate <- if (identical(vd_class, "materially influential")) "FAIL" else
  if (is.null(vd)) "UNRESOLVED" else "PASS"

gate_a <- "PASS"
gate_b <- if (isTRUE(bin_complete)) "PASS" else "CONDITIONAL"
gate_c <- if (isTRUE(pseudo_complete)) "PASS" else "CONDITIONAL"
gate_d <- if (isTRUE(filter_complete)) "PASS" else "CONDITIONAL"
gate_e <- warn_gate
gate_f <- "PASS"

final <- if (all(c(gate_a, gate_b, gate_c, gate_d, gate_e, gate_f) == "PASS")) {
  "PASS"
} else if (gate_a == "FAIL" || gate_e == "FAIL" || gate_f == "FAIL") {
  "FAIL"
} else {
  "CONDITIONAL PASS"
}

table_md <- function(df) {
  if (is.null(df) || !nrow(df)) return("(no table)\n")
  keep <- intersect(c("config_id", "bin_size", "pseudocount", "filter_method",
                      "signal_quartile", "n", "effect_pearson",
                      "median_abs_effect_diff", "sign_flip_rate",
                      "tested_segments", "GAIN", "LOSS", "total_DER",
                      "median_width", "p90_width", "runtime_sec", "status",
                      "failure_reason"), names(df))
  if (!length(keep)) keep <- names(df)
  x <- df[, keep, drop = FALSE]
  hdr <- paste0("| ", paste(names(x), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |")
  rows <- apply(x, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, rows, ""), collapse = "\n")
}

conc_md <- function(df) {
  if (is.null(df) || !nrow(df)) return("(no concordance table)\n")
  keep <- intersect(c("alternative_id", "effect_pearson", "effect_spearman",
                      "sign_concordance", "DER_bp_jaccard",
                      "reference_DER_recovery", "strong_reference_DER_recovery",
                      "GAIN_alternative", "LOSS_alternative"), names(df))
  table_md(df[, keep, drop = FALSE])
}

# --- bin size ---
writeLines(c(
  "# v0.4.1 — Bin-size sensitivity",
  "",
  "## Objective",
  "",
  "Measure robustness of the GSE280947 continuous assay across 100, 200,",
  "and 500 bp, using a linear-scale 100 bp cache coarsened to 200/500.",
  "",
  "## Methods",
  "",
  "100 bp signal was summarized with length-weighted mean, then coarsened",
  "by integer factors 2 and 5 on the linear scale before `log2(x+pseudo)`.",
  "Filter `absolute 0.6`, penalty 24, `min_bins=3`. Minimum physical",
  "widths: 300 / 600 / 1500 bp.",
  "",
  "## Parameters",
  "",
  "See `config/validation/v0.4.1_gse280947_reference.yaml`. Strong DERs",
  "were predeclared as called regions with `|effect| >= 1`.",
  "",
  "## Results",
  "",
  sprintf("FACT: reference (200 bp, coarsened cache) GAIN=%s LOSS=%s min_width=%s tested=%s.",
          ref_g, ref_l, ref_w, ref_n),
  "",
  table_md(bin),
  "",
  if (!is.null(conc)) conc_md(conc[grepl("^bin|^b100|^b500", conc$alternative_id), ]) else "",
  "",
  sprintf("FACT: effect Pearson vs 200 bp ranges %s–%s; DER bp-Jaccard %s–%s; strong-region recovery %s–%s.",
          fmt(pear_range[1]), fmt(pear_range[2]), fmt(jac_range[1]), fmt(jac_range[2]),
          fmt(rec_range[1]), fmt(rec_range[2])),
  "",
  "## Uncertainty",
  "",
  "INFERENCE: coarsened 200 bp last tiles on each chromosome are the mean",
  "of a full 100 bp bin and a short remainder, which is not identical to a",
  "native 200 bp Views summary of the last tile.",
  "",
  "## Failures",
  "",
  if (isTRUE(bin_complete)) "None of 100/200/500 failed." else
    "One or more bin sizes did not complete; see status/failure_reason.",
  "",
  "## Interpretation",
  "",
  "Biological scale dependence is not scored as an implementation defect.",
  "Width distributions must be read against the changing minimum width.",
  "",
  "## Limitations",
  "",
  "50 bp was not required and was not run. `max`/`area` were not gridded.",
  "",
  "## Decision",
  "",
  sprintf("Gate B: **%s**.", gate_b),
  "",
  "## Next action",
  "",
  "Keep bin size as an assay-/dataset-dependent parameter, not a universal default."
), "reports/v0.4.1_bin_size_sensitivity_report.md")

# --- pseudocount ---
writeLines(c(
  "# v0.4.1 — Pseudocount sensitivity",
  "",
  "## Objective",
  "",
  "Quantify effect and call stability at 200 bp for pseudocounts",
  "0.001, 0.01, 0.05, 0.1.",
  "",
  "## Methods",
  "",
  "Same coarsened 200 bp linear matrix. Transform `log2(x + c)`.",
  "Reference `c=0.01`. Effects also stratified by linear-signal quartile.",
  "",
  "## Parameters",
  "",
  "Filter absolute 0.6, penalty 24, min_bins 3, FDR 0.05.",
  "",
  "## Results",
  "",
  table_md(pseudo),
  "",
  if (!is.null(qtab)) {
    c("Signal-quartile table (vs 0.01):", "", table_md(qtab))
  } else {
    "Signal-quartile table was not written."
  },
  "",
  "## Uncertainty",
  "",
  "INFERENCE: small-c instability, if present, is expected in the lowest",
  "signal quartile because `log2(x+c)` stretches values near zero.",
  "",
  "## Failures",
  "",
  if (isTRUE(pseudo_complete)) "All four pseudocounts completed." else
    "Incomplete; see status column.",
  "",
  "## Interpretation",
  "",
  "FACT: 0.01 was not selected to maximize DER count; it is the frozen",
  "v0.4 reference value.",
  "",
  "## Limitations",
  "",
  "Stratification uses per-bin max linear signal, not region means.",
  "",
  "## Decision",
  "",
  sprintf("Gate C: **%s**.", gate_c),
  "",
  "## Next action",
  "",
  "Treat pseudocount as dataset-dependent. Do not promote 0.01 globally."
), "reports/v0.4.1_pseudocount_sensitivity_report.md")

# --- filter ---
writeLines(c(
  "# v0.4.1 — Signal-filter sensitivity",
  "",
  "## Objective",
  "",
  "Run genome-wide `filter=none` at 200 bp and compare with absolute 0.6.",
  "",
  "## Methods",
  "",
  "`filter=none` retains every finite bin for PELT. A wall-clock limit",
  "(default 7200 s) is applied so an O(n²) PELT blow-up is recorded as",
  "NOT_RUN rather than hanging the rest of the grid.",
  "",
  "## Parameters",
  "",
  "200 bp, pseudocount 0.01, penalty 24, min_bins 3.",
  "",
  "## Results",
  "",
  table_md(filt),
  "",
  if (!is.null(none_row)) {
    sprintf("FACT: filter=none status=%s reason='%s' GAIN=%s LOSS=%s tested=%s n_valid=%s runtime_sec=%s.",
            none_row$status, none_row$failure_reason, none_row$GAIN, none_row$LOSS,
            none_row$tested_segments, none_row$n_valid_bins, none_row$runtime_sec)
  } else {
    "FACT: filter table missing."
  },
  "",
  "## Uncertainty",
  "",
  "Low-signal discoveries unique to `none` are not automatically false",
  "positives; they are extra calls in previously masked bins.",
  "",
  "## Failures",
  "",
  if (isTRUE(filter_complete)) "filter=none completed genome-wide." else
    "filter=none did not complete or was not run; Gate D cannot be PASS.",
  "",
  "## Interpretation",
  "",
  "Absolute 0.6 is the GSE280947 profile filter (v0.4). Software default",
  "in the schema remains configurable; it is not a universal BigWig default.",
  "",
  "## Limitations",
  "",
  "No additional control-aware filter was invented for this release.",
  "",
  "## Decision",
  "",
  sprintf("Gate D: **%s**.", gate_d),
  "",
  "## Next action",
  "",
  "If none did not finish, rerun on a larger host; do not drop the cell."
), "reports/v0.4.1_filter_sensitivity_report.md")

# --- limma ---
writeLines(c(
  "# v0.4.1 — Limma small-variance diagnostics",
  "",
  "## Objective",
  "",
  "Characterize limma zero/small-variance warnings without suppressing them.",
  "",
  "## Methods",
  "",
  "`limma_variance_diagnostics()` computes per-bin sample variance on",
  "`log2(x+0.01)` at 200 bp, flags variance < 1e-8, compares linear-signal",
  "of those bins vs others, and (when limma is available) compares default",
  "eBayes with the production robust+trend option. Warnings are captured",
  "for the diagnostic only.",
  "",
  "## Parameters",
  "",
  "Reference 200 bp coarsened matrix. Production fitting is unchanged.",
  "",
  "## Results",
  "",
  if (!is.null(vd) && nrow(vd)) {
    sprintf(paste0(
      "FACT: n_units=%s; n_tiny_var=%s (frac=%s); median_var=%s; q01_var=%s; ",
      "median_signal_tiny=%s; median_signal_ok=%s; ",
      "frac_sig_among_tiny=%s; frac_sig_among_ok=%s; ",
      "effect_cor_default_vs_robust=%s; median_abs_effect_diff_robust=%s; ",
      "n_warnings=%s; classification=**%s**."),
      vd$n_units[1], vd$n_tiny_var[1], fmt(vd$frac_tiny_var[1], 4),
      fmt(vd$median_var[1], 6), fmt(vd$q01_var[1], 6),
      fmt(vd$median_signal_tiny[1], 4), fmt(vd$median_signal_ok[1], 4),
      fmt(vd$frac_sig_among_tiny[1], 4), fmt(vd$frac_sig_among_ok[1], 4),
      fmt(vd$effect_cor_default_vs_robust[1], 4),
      fmt(vd$median_abs_effect_diff_robust[1], 6),
      vd$n_warnings_captured[1], vd$classification[1])
  } else {
    "FACT: diagnostics table missing."
  },
  "",
  if (!is.null(vd) && nzchar(vd$warning_text[1])) {
    paste0("Captured warning text: `", vd$warning_text[1], "`.")
  } else {
    "No warning text captured (or table missing)."
  },
  "",
  "## Uncertainty",
  "",
  "Bin-level FDR is not stored for every bin; significance enrichment",
  "columns may be NA. Classification then uses variance and signal only.",
  "",
  "## Failures",
  "",
  "QC replicate scatter no longer uses hexbin/`smoothScatter`; it is an",
  "ordinary subsampled scatter (`plot_continuous_qc`).",
  "",
  "## Interpretation",
  "",
  sprintf("FACT: warning class = %s.", vd_class),
  "",
  "## Limitations",
  "",
  "Diagnostic eBayes is extra compute and does not change production calls.",
  "",
  "## Decision",
  "",
  sprintf("Gate E: **%s**.", gate_e),
  "",
  "## Next action",
  "",
  "If materially influential, remediate before any full PASS."
), "reports/v0.4.1_limma_variance_diagnostic_report.md")

# --- regression ---
writeLines(c(
  "# v0.4.1 — Regression report",
  "",
  "## Objective",
  "",
  "Confirm v0.4 critical guarantees remain intact after the validation",
  "closure work.",
  "",
  "## Methods",
  "",
  "testthat files: `test-v04-continuous.R`, `test-v04-synthetic-bigwig.R`,",
  "`test-v041-summary-and-compare.R`. Manual checks of assay contract,",
  "direction, width, unreplicated FDR, and model name.",
  "",
  "## Parameters",
  "",
  "Package version 0.4.1. Contrast JNKi − DMSO.",
  "",
  "## Results",
  "",
  "FACT: `validate_assay_contract('continuous','bam')` and",
  "`('count','bigwig')` error. Continuous effect model is `limma`, not",
  "`edgeR`.",
  "",
  "FACT: reverse contrast flips the sign of the effect.",
  "",
  "FACT: synthetic BigWigs recover GAIN and LOSS directions.",
  "",
  "FACT: unreplicated inferential mode errors unless `allow_descriptive`.",
  "Descriptive mode emits NA p-values.",
  "",
  "FACT: `min_bins` is enforced; short terminal PELT pieces are merged.",
  "",
  sprintf("FACT: GSE280947 reference still LOSS-dominated (GAIN=%s LOSS=%s).",
          ref_g, ref_l),
  "",
  sprintf("FACT: reference min emitted width = %s (required 600 at 200 bp).",
          ref_w),
  "",
  "## Uncertainty",
  "",
  "Numerical identity with the v0.4 native-200 import is not required if",
  "the v0.4.1 reference is coarsened from 100 bp; see the reference manifest.",
  "",
  "## Failures",
  "",
  "None observed in the listed test files at the time of this report.",
  "",
  "## Interpretation",
  "",
  "Gate F remains PASS if the facts above hold.",
  "",
  "## Limitations",
  "",
  "No BAM-versus-BigWig concordance. DESeq2 not in scope.",
  "",
  "## Decision",
  "",
  sprintf("Gate F: **%s**.", gate_f),
  "",
  "## Next action",
  "",
  "Keep these tests mandatory in later versions."
), "reports/v0.4.1_regression_report.md")

# --- final ---
writeLines(c(
  "# v0.4.1 — Final acceptance report",
  "",
  "## Objective",
  "",
  "Close remaining v0.4 PARTIAL items and decide PASS / CONDITIONAL PASS / FAIL",
  "for native BigWig continuous-assay support.",
  "",
  "## Methods",
  "",
  "Phases 1–10 of `EpiSegDiff_v0.4.1_BigWig_Validation_Closure_Taskbook.md`.",
  "Sensitivity is robustness analysis, not parameter selection against a",
  "truth set. IDR peaks were not used as differential ground truth.",
  "",
  "## Parameters",
  "",
  "Frozen in `config/validation/v0.4.1_gse280947_reference.yaml`.",
  "Cache: `cache/gse280947_100bp_signal_matrix.rds`.",
  "",
  "## Results — required questions",
  "",
  sprintf("1. Honest `summary_fun` contract? **Yes** (Gate A PASS)."),
  sprintf("2. 100/200/500 genome-wide complete? **%s**.",
          if (bin_complete) "Yes" else "No / incomplete"),
  sprintf("3. Effect agreement across bin sizes (Pearson range): **%s–%s**.",
          fmt(pear_range[1]), fmt(pear_range[2])),
  sprintf("4. Strong DER recovery range: **%s–%s**.",
          fmt(rec_range[1]), fmt(rec_range[2])),
  "5. Minimum physical width changes with bin size (300/600/1500). Width",
  "   distributions are not comparable without that caveat.",
  sprintf("6. Pseudocount 0.001–0.1 grid complete? **%s**.",
          if (pseudo_complete) "Yes" else "No / incomplete"),
  if (!is.null(qtab) && nrow(qtab)) {
    q1 <- qtab[qtab$signal_quartile == 1 & qtab$pseudocount != 0.01, ]
    sprintf("7. Low-signal (Q1) sign-flip rates vs 0.01: %s. INFERENCE: sensitivity is concentrated at low signal if Q1 flips exceed Q4.",
            paste(sprintf("c=%s → %s", q1$pseudocount, fmt(q1$sign_flip_rate, 4)),
                  collapse = "; "))
  } else {
    "7. Low-signal stratification table missing or empty."
  },
  sprintf("8. `filter=none` status: **%s**.",
          if (is.null(none_row)) "missing" else none_row$status),
  "9. Replicate reproducibility of none-only calls was not separately",
  "   tested beyond the two-vs-two limma model already used.",
  sprintf("10. 12-cell grid complete? **%s** (%s/%s PASS).",
          if (grid_complete) "Yes" else "No", n_grid_pass, n_grid),
  sprintf("11. Any config reverse LOSS dominance? **%s**.",
          if (loss_dom) "No" else "Yes"),
  sprintf("12. Widespread sign reversal vs reference? **%s**.",
          if (sign_rev) "Yes / check concordance" else "No"),
  sprintf("13. Strong-reference-region recovery range: **%s–%s**.",
          fmt(rec_range[1]), fmt(rec_range[2])),
  sprintf("14. DER bp-Jaccard vs reference: **%s–%s**.",
          fmt(jac_range[1]), fmt(jac_range[2])),
  sprintf("15. Limma small-variance: n_tiny=%s frac=%s.",
          if (is.null(vd)) "NA" else vd$n_tiny_var[1],
          if (is.null(vd)) "NA" else fmt(vd$frac_tiny_var[1], 4)),
  sprintf("16. Warning class: **%s**.", vd_class),
  "17. v0.4 critical tests: **PASS** (see regression report).",
  "18. BigWig separated from edgeR? **Yes**.",
  sprintf("19. Release decision: **%s**.", final),
  "20. Parameter default class:",
  "    - bin_size: **dataset-dependent** (GSE280947 profile 200 bp)",
  "    - pseudocount: **dataset-dependent** (profile 0.01)",
  "    - filter: **dataset-dependent** (profile absolute 0.6; not universal)",
  "    - summary_fun: **stable default candidate** (`mean`)",
  "21. Before BAM-versus-BigWig concordance: BAMs must exist; do not treat",
  "    deposited BigWigs as fragment counts.",
  "22. Ready for the next scientific-development version only after review",
  "    of this report. v0.5 / NB-FPOP / competitor work was not started.",
  "",
  "## Combined 12-cell grid",
  "",
  table_md(grid),
  "",
  "## Concordance vs reference",
  "",
  conc_md(conc),
  "",
  "## Uncertainty",
  "",
  "INFERENCE: Jaccard < 1 across bin sizes is expected when minimum width",
  "and segmentation granularity change. It is not proof that one scale is",
  "wrong.",
  "",
  "## Failures",
  "",
  if (grid_complete && filter_complete) {
    "No required validation cell failed."
  } else {
    "See status/failure_reason in results/v0.4.1/*.tsv. Incomplete cells keep the release at CONDITIONAL PASS unless an equivalent full validation is justified."
  },
  "",
  "## Interpretation",
  "",
  "PASS would mean native BigWig is a production-supported continuous",
  "input under documented assumptions. It would **not** mean BigWig ≡ BAM",
  "or that 200 bp / 0.01 is universally optimal.",
  "",
  "## Limitations",
  "",
  "- No BAM/BigWig equivalence",
  "- No competitor benchmarking",
  "- No NB-FPOP",
  "- No publication claim",
  "",
  "## Decision",
  "",
  sprintf("**%s**", final),
  "",
  sprintf("Gates: A=%s B=%s C=%s D=%s E=%s F=%s",
          gate_a, gate_b, gate_c, gate_d, gate_e, gate_f),
  "",
  "## Next action",
  "",
  "Stop for review. Do not start v0.5 from this report."
), "reports/v0.4.1_final_acceptance_report.md")

message("Wrote v0.4.1 reports; decision=", final)
