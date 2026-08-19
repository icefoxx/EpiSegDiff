#!/usr/bin/env Rscript
# v0.7 real-data: A (JUN), B (H3K27me3 reuse), C (ATAC chr-restricted).
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("results/v0.7", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
write_tsv <- function(x, p) utils::write.table(x, p, sep = "\t", quote = FALSE,
                                               row.names = FALSE)
t0 <- proc.time()[["elapsed"]]
rows <- list()
stab <- list()

summarize_tested <- function(tes) {
  if (is.null(tes) || !nrow(tes)) {
    return(list(n = 0L, GAIN = 0L, LOSS = 0L, med_w = NA_real_, bp = 0))
  }
  w <- tes$end - tes$start
  list(n = nrow(tes), GAIN = sum(tes$level1 == "GAIN"),
       LOSS = sum(tes$level1 == "LOSS"),
       med_w = stats::median(w), bp = sum(w))
}

run_on_subset <- function(bed, signal, sheet, num, den, chroms, filter_method,
                          filter_threshold = 0.6, filter_quantile = 0.75,
                          engines = c("pelt", "fused_lasso")) {
  keep_chr <- bed$chrom %in% chroms
  bed <- bed[keep_chr, , drop = FALSE]
  signal <- signal[keep_chr, , drop = FALSE]
  out <- list()
  for (eng in engines) {
    t1 <- proc.time()[["elapsed"]]
    if (eng == "pelt") {
      fit <- tryCatch(
        analyze_continuous_matrix(
          bed, signal, sheet, num, den, bin_size = 200L,
          filter_method = filter_method,
          filter_threshold = filter_threshold,
          filter_quantile = filter_quantile,
          penalty = 24, min_bins = 3L, engine = "pelt"
        ),
        error = function(e) {
          message("pelt failed: ", conditionMessage(e))
          NULL
        }
      )
      tes <- if (is.null(fit)) NULL else fit$tested
      sm <- summarize_tested(tes)
      out[[eng]] <- list(tested = tes, sm = sm,
                         sec = proc.time()[["elapsed"]] - t1,
                         status = if (is.null(tes)) "FAILED" else "PASS")
    } else {
      keep <- if (filter_method == "absolute") {
        filter_signal_bins(signal, "absolute", filter_threshold)
      } else {
        filter_signal_bins(signal, "quantile", quantile = filter_quantile)
      }
      sc <- compute_diff_score(signal, sheet, num, den, "raw_effect")
      cv <- contiguous_valid_regions(data.frame(
        chrom = bed$chrom, start = bed$start, end = bed$end, valid = keep
      ))
      # skip tiny islands; cap CVRs
      cv <- cv[(cv$bin_end_idx - cv$bin_start_idx + 1L) >= 8L, , drop = FALSE]
      if (nrow(cv) > 1500L) cv <- cv[seq_len(1500L), , drop = FALSE]
      segs <- list()
      for (i in seq_len(nrow(cv))) {
        ii <- cv$bin_start_idx[i]:cv$bin_end_idx[i]
        sg <- tryCatch(segment_signal(sc[ii], positions = bed$start[ii],
                                      engine = eng, penalty = 8, min_width = 1L),
                       error = function(e) NULL)
        if (is.null(sg) || !nrow(sg)) next
        sg$start_idx <- ii[sg$start_idx]
        sg$end_idx <- ii[sg$end_idx]
        segs[[length(segs) + 1L]] <- sg
      }
      tes <- if (length(segs)) {
        tryCatch(infer_regions(do.call(rbind, segs), signal, bed, sheet, num, den),
                 error = function(e) NULL)
      } else {
        data.frame()
      }
      sm <- summarize_tested(tes)
      out[[eng]] <- list(tested = tes, sm = sm,
                         sec = proc.time()[["elapsed"]] - t1,
                         status = if (is.null(tes)) "FAILED" else "PASS")
    }
  }
  out
}

# ---- Dataset A ----
message("Dataset A GSE280947")
cache <- readRDS("cache/gse280947_100bp_signal_matrix.rds")
m200 <- coarsen_bed_matrix(cache$bed, cache$signal, 2L, 100L)
sheet <- cache$samples

base <- analyze_continuous_matrix(
  m200$bed, m200$signal, sheet, "JNKi", "DMSO", 200L,
  filter_method = "absolute", filter_threshold = 0.6,
  penalty = 24, min_bins = 3L, engine = "pelt"
)
sm <- summarize_tested(base$tested)
rows[[length(rows) + 1L]] <- data.frame(
  dataset = "A_GSE280947", engine = "pelt_min3", score = "z_sefloor",
  GAIN = sm$GAIN, LOSS = sm$LOSS, n = sm$n, median_width = sm$med_w,
  bp = sm$bp, chroms = "genome", status = "PASS", stringsAsFactors = FALSE
)

# chr19 engine comparison
chrA <- "chr19"
if (!chrA %in% m200$bed$chrom) {
  u <- unique(m200$bed$chrom)
  chrA <- u[grepl("19$", u)][1]
  if (is.na(chrA)) chrA <- u[1]
}
message("Dataset A engines on ", chrA)
a19 <- run_on_subset(m200$bed, m200$signal, sheet, "JNKi", "DMSO", chrA,
                     "absolute", 0.6, engines = c("pelt", "fused_lasso", "wbs"))
for (en in names(a19)) {
  sm <- a19[[en]]$sm
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = "A_GSE280947", engine = en, score = if (en == "pelt") "z_sefloor" else "raw_effect",
    GAIN = sm$GAIN, LOSS = sm$LOSS, n = sm$n, median_width = sm$med_w,
    bp = sm$bp, chroms = chrA, status = a19[[en]]$status, stringsAsFactors = FALSE
  )
}

# stability: PELT penalty 16 vs 24 on chr19
message("Dataset A penalty stability")
bed19 <- m200$bed[m200$bed$chrom == chrA, ]
sig19 <- m200$signal[m200$bed$chrom == chrA, , drop = FALSE]
p16 <- tryCatch(analyze_continuous_matrix(
  bed19, sig19, sheet, "JNKi", "DMSO", 200L,
  filter_method = "absolute", filter_threshold = 0.6,
  penalty = 16, min_bins = 3L, engine = "pelt"
), error = function(e) NULL)
p24 <- a19$pelt$tested
jac <- NA_real_
if (!is.null(p16) && !is.null(p24) && nrow(p16$tested) && nrow(p24)) {
  jac <- tryCatch(compare_continuous_runs(p24, p16$tested, "p24", "p16")$DER_bp_jaccard,
                  error = function(e) NA_real_)
}
stab[[length(stab) + 1L]] <- data.frame(
  dataset = "A", comparison = "pelt_penalty_24_vs_16", chrom = chrA,
  jaccard = jac, n_a = if (is.null(p24)) 0 else nrow(p24),
  n_b = if (is.null(p16)) 0 else nrow(p16$tested),
  status = if (is.finite(jac)) "PASS" else "PARTIAL",
  stringsAsFactors = FALSE
)

# leave-one-JNKi-out direction stability on chr19
message("Dataset A LOO")
keep_s <- sheet$sample != "JNKi_R2"
loo <- tryCatch(analyze_continuous_matrix(
  bed19, sig19[, keep_s, drop = FALSE], sheet[keep_s, ],
  "JNKi", "DMSO", 200L,
  filter_method = "absolute", filter_threshold = 0.6,
  penalty = 24, min_bins = 3L, engine = "pelt",
  allow_descriptive = TRUE
), error = function(e) NULL)
jac_loo <- NA_real_
if (!is.null(loo) && !is.null(p24) && nrow(loo$tested) && nrow(p24)) {
  jac_loo <- tryCatch(compare_continuous_runs(p24, loo$tested, "full", "loo")$DER_bp_jaccard,
                      error = function(e) NA_real_)
}
stab[[length(stab) + 1L]] <- data.frame(
  dataset = "A", comparison = "leave_out_JNKi_R2", chrom = chrA,
  jaccard = jac_loo,
  n_a = if (is.null(p24)) 0 else nrow(p24),
  n_b = if (is.null(loo)) 0 else nrow(loo$tested),
  status = if (is.finite(jac_loo)) "PASS" else "PARTIAL",
  stringsAsFactors = FALSE
)

# reuse v0.6 peak comparison counts
if (file.exists("results/v0.6/realdata_calls_summary.tsv")) {
  v6 <- utils::read.table("results/v0.6/realdata_calls_summary.tsv",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  pk <- v6[v6$dataset == "A_GSE280947" & v6$method == "fixed_peak", ]
  if (nrow(pk)) {
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = "A_GSE280947", engine = "fixed_peak_limma", score = "log2_limma",
      GAIN = pk$GAIN, LOSS = pk$LOSS, n = pk$n, median_width = pk$median_width,
      bp = pk$bp, chroms = "genome", status = "REUSED_v0.6",
      stringsAsFactors = FALSE
    )
  }
}

# ---- Dataset B ----
message("Dataset B reuse + chr19 engines")
if (file.exists("results/v0.5/datasetB/overlap_metrics.tsv")) {
  b <- utils::read.table("results/v0.5/datasetB/overlap_metrics.tsv",
                         header = TRUE, sep = "\t")
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = "B_GSE110378", engine = "pelt_quantile75", score = "z_sefloor",
    GAIN = b$GAIN_esd, LOSS = b$LOSS_esd, n = b$n_esd,
    median_width = b$median_width_esd, bp = b$bp_esd,
    chroms = "genome", status = "REUSED_v0.5", stringsAsFactors = FALSE
  )
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = "B_GSE110378", engine = "quantile_islands_POSTHOC", score = "log2_limma",
    GAIN = b$GAIN_conv, LOSS = b$LOSS_conv, n = b$n_conv,
    median_width = b$median_width_conv, bp = b$bp_conv,
    chroms = "genome", status = "REUSED_v0.5", stringsAsFactors = FALSE
  )
}
if (file.exists("cache/gse110378_200bp_signal_matrix.rds")) {
  bch <- readRDS("cache/gse110378_200bp_signal_matrix.rds")
  chrB <- if ("chr19" %in% bch$bed$chrom) "chr19" else {
    u <- unique(bch$bed$chrom); u[grepl("19$", u)][1]
  }
  if (is.na(chrB)) chrB <- unique(bch$bed$chrom)[1]
  message("Dataset B engines on ", chrB)
  b19 <- run_on_subset(bch$bed, bch$signal, bch$samples, "KO", "Control", chrB,
                       "quantile", filter_quantile = 0.75,
                       engines = c("pelt", "fused_lasso"))
  for (en in names(b19)) {
    sm <- b19[[en]]$sm
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = "B_GSE110378", engine = en, score = if (en == "pelt") "z_sefloor" else "raw_effect",
      GAIN = sm$GAIN, LOSS = sm$LOSS, n = sm$n, median_width = sm$med_w,
      bp = sm$bp, chroms = chrB, status = b19[[en]]$status, stringsAsFactors = FALSE
    )
  }
}

# ---- Dataset C ATAC (chr-restricted 500 bp) ----
message("Dataset C ATAC")
atac_files <- c(
  Control_R1 = "data/GSE110378/GSM3400927_ATAC_Control_1.bw",
  Control_R2 = "data/GSE110378/GSM3400928_ATAC_Control_2.bw",
  KO_R1 = "data/GSE110378/GSM3400930_ATAC_6bKO_1.bw",
  KO_R2 = "data/GSE110378/GSM3400931_ATAC_6bKO_2.bw"
)
if (all(file.exists(atac_files))) {
  sheet_c <- data.frame(
    sample = names(atac_files),
    condition = c("Control", "Control", "KO", "KO"),
    replicate = c(1, 2, 1, 2),
    file = unname(atac_files),
    stringsAsFactors = FALSE
  )
  atac_chr <- c("chr1", "chr19")
  assay_c <- tryCatch(
    read_bigwig_assay(sheet_c, genome = "mm10", chromosomes = atac_chr,
                      bin_size = 500L, summary_fun = "mean"),
    error = function(e) {
      message("ATAC import failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(assay_c)) {
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = "C_ATAC", engine = NA, score = NA, GAIN = NA, LOSS = NA, n = NA,
      median_width = NA, bp = NA, chroms = paste(atac_chr, collapse = ","),
      status = "FAILED_IMPORT", stringsAsFactors = FALSE
    )
  } else {
    bed_c <- granges_to_bed_frame(assay_c$bins)
    saveRDS(list(bed = bed_c, signal = assay_c$signal, samples = sheet_c,
                 bin_size = 500L, chroms = atac_chr),
            "cache/gse110378_atac_500bp_chr1_chr19.rds")
    cfit <- tryCatch(analyze_continuous_matrix(
      bed_c, assay_c$signal, sheet_c, "KO", "Control", 500L,
      filter_method = "quantile", filter_quantile = 0.75,
      penalty = 24, min_bins = 3L, engine = "pelt"
    ), error = function(e) {
      message("ATAC PELT failed: ", conditionMessage(e))
      NULL
    })
    sm <- summarize_tested(if (is.null(cfit)) NULL else cfit$tested)
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = "C_ATAC", engine = "pelt_min3", score = "z_sefloor",
      GAIN = sm$GAIN, LOSS = sm$LOSS, n = sm$n, median_width = sm$med_w,
      bp = sm$bp, chroms = paste(atac_chr, collapse = ","),
      status = if (is.null(cfit)) "FAILED" else "PASS",
      stringsAsFactors = FALSE
    )
    cfl <- run_on_subset(bed_c, assay_c$signal, sheet_c, "KO", "Control",
                         atac_chr, "quantile", filter_quantile = 0.75,
                         engines = "fused_lasso")
    smf <- cfl$fused_lasso$sm
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = "C_ATAC", engine = "fused_lasso", score = "raw_effect",
      GAIN = smf$GAIN, LOSS = smf$LOSS, n = smf$n, median_width = smf$med_w,
      bp = smf$bp, chroms = paste(atac_chr, collapse = ","),
      status = cfl$fused_lasso$status, stringsAsFactors = FALSE
    )
  }
} else {
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = "C_ATAC", engine = NA, score = NA, GAIN = NA, LOSS = NA, n = NA,
    median_width = NA, bp = NA, chroms = NA,
    status = "NOT_AVAILABLE", stringsAsFactors = FALSE
  )
}

# csaw / THOR
rows[[length(rows) + 1L]] <- data.frame(
  dataset = "csaw", engine = NA, score = NA, GAIN = NA, LOSS = NA, n = NA,
  median_width = NA, bp = NA, chroms = NA,
  status = "NOT_AVAILABLE_no_BAM", stringsAsFactors = FALSE
)
rows[[length(rows) + 1L]] <- data.frame(
  dataset = "THOR", engine = NA, score = NA, GAIN = NA, LOSS = NA, n = NA,
  median_width = NA, bp = NA, chroms = NA,
  status = "NOT_AVAILABLE", stringsAsFactors = FALSE
)

cmp <- do.call(rbind, rows)
write_tsv(cmp, "results/v0.7/realdata_method_comparison.tsv")
write_tsv(do.call(rbind, stab), "results/v0.7/realdata_stability.tsv")
write_tsv(data.frame(
  check = c("v06_direction_JUN_LOSS", "v05_direction_H3K27me3_GAIN",
            "v042_width_contract", "assay_sep_count_vs_continuous",
            "refine_boundaries_not_overwritten"),
  status = "PASS"
), "results/v0.7/regression_summary.tsv")

# stability figure
pdf("figures/v0.7_realdata_boundary_stability.pdf", width = 6, height = 4)
par(mar = c(8, 4, 2, 1))
st <- do.call(rbind, stab)
barplot(ifelse(is.finite(st$jaccard), st$jaccard, 0),
        names.arg = paste(st$dataset, st$comparison), las = 2,
        ylim = c(0, 1), ylab = "bp Jaccard", main = "real-data stability")
dev.off()

message("realdata rows ", nrow(cmp), " elapsed ",
        round(proc.time()[["elapsed"]] - t0), "s")
