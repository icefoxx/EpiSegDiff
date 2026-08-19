#!/usr/bin/env Rscript
# GSE280947 Level I proxy from deposited spike-in+CPM bigWigs.
# NOT fragment counts. NOT the NB/edgeR count model.
# Contrast: JNKi − DMSO. GAIN = higher JUN after inhibitor.

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
})

t_start <- proc.time()[["elapsed"]]
set.seed(20271)
dir.create("results", showWarnings = FALSE)
dir.create("data/GSE280947/episegdiff", recursive = TRUE, showWarnings = FALSE)

base <- "data/GSE280947/conventional"
jun_bw <- c(
  dmso_jun_r1 = file.path(base, "GSM8607876_DMSO1_cJun.bigWig"),
  dmso_jun_r2 = file.path(base, "GSM8607877_DMSO2_cJun.bigWig"),
  jnki_jun_r1 = file.path(base, "GSM8607882_JNKI1_cJun.bigWig"),
  jnki_jun_r2 = file.path(base, "GSM8607883_JNKI2_cJun.bigWig")
)
igg_bw <- c(
  dmso_igg_r1 = file.path(base, "GSM8607874_DMSO1_Neg.bigWig"),
  dmso_igg_r2 = file.path(base, "GSM8607875_DMSO2_Neg.bigWig"),
  jnki_igg_r1 = file.path(base, "GSM8607880_JNKI1_Neg.bigWig"),
  jnki_igg_r2 = file.path(base, "GSM8607881_JNKI2_Neg.bigWig")
)
stopifnot(all(file.exists(jun_bw)), all(file.exists(igg_bw)))

samples <- data.frame(
  sample_id = names(jun_bw),
  condition = c("DMSO", "DMSO", "JNKi", "JNKi"),
  replicate = c(1L, 2L, 1L, 2L),
  stringsAsFactors = FALSE
)
cond <- factor(samples$condition, levels = c("DMSO", "JNKi"))

bin_size <- 200L
std_chr <- c(as.character(1:19), "X", "Y")
si <- seqinfo(BigWigFile(jun_bw[[1]]))
si <- si[intersect(std_chr, seqnames(si))]
bins <- GenomicRanges::tileGenome(si, tilewidth = bin_size,
                                  cut.last.tile.in.chrom = TRUE)
if (inherits(bins, "GRangesList")) {
  bins <- unlist(bins, use.names = FALSE)
}
message("bins ", length(bins), " bin_size ", bin_size)

import_binned <- function(paths, bins) {
  mat <- matrix(0, length(bins), length(paths),
                dimnames = list(NULL, names(paths)))
  for (j in seq_along(paths)) {
    message("import ", names(paths)[j])
    cov <- rtracklayer::import(paths[[j]], as = "RleList")
    keep <- intersect(seqlevels(bins), names(cov))
    bins_j <- bins[as.character(seqnames(bins)) %in% keep]
    seqlevels(bins_j) <- keep
    cov <- cov[keep]
    ba <- GenomicRanges::binnedAverage(bins_j, cov, "score")
    idx <- as.character(seqnames(bins)) %in% keep
    mat[idx, j] <- as.numeric(S4Vectors::mcols(ba)$score)
    rm(cov, ba)
    gc(verbose = FALSE)
  }
  mat
}

score <- import_binned(jun_bw, bins)
igg <- import_binned(igg_bw, bins)

# QC: IDR peak enrichment and replicate correlation
read_idr <- function(path) {
  x <- utils::read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  GenomicRanges::GRanges(x[[1]], IRanges::IRanges(as.integer(x[[2]]) + 1L, as.integer(x[[3]])))
}
idr_d <- read_idr("data/GSE280947/conventional/idr_dmso_jun.bed")
idr_j <- read_idr("data/GSE280947/conventional/idr_jnki_jun.bed")
idr_u <- GenomicRanges::reduce(c(idr_d, idr_j))
seqlevelsStyle(bins) <- seqlevelsStyle(idr_d)[1]
hits_d <- unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(bins, idr_d)))
hits_u <- unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(bins, idr_u)))

qc <- data.frame(
  metric = c(
    "cor_dmso_jun_reps_idr", "cor_jnki_jun_reps_idr",
    "mean_dmso_jun_in_dmsoIDR", "mean_dmso_igg_in_dmsoIDR",
    "mean_jnki_jun_in_dmsoIDR", "mean_jnki_igg_in_dmsoIDR",
    "jun_over_igg_dmso_idr", "jun_over_igg_jnki_in_dmsoIDR"
  ),
  value = c(
    stats::cor(score[hits_d, 1], score[hits_d, 2], use = "complete"),
    stats::cor(score[hits_d, 3], score[hits_d, 4], use = "complete"),
    mean(score[hits_d, 1:2]), mean(igg[hits_d, 1:2]),
    mean(score[hits_d, 3:4]), mean(igg[hits_d, 3:4]),
    mean(score[hits_d, 1:2]) / max(mean(igg[hits_d, 1:2]), 1e-8),
    mean(score[hits_d, 3:4]) / max(mean(igg[hits_d, 3:4]), 1e-8)
  )
)
utils::write.table(qc, "data/GSE280947/qc/bigwig_proxy_qc.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
print(qc)

# Per-bin differential (JNKi − DMSO) on as-deposited scores
pseudo <- 0.01
logsc <- log2(score + pseudo)
m0 <- rowMeans(logsc[, cond == "DMSO", drop = FALSE])
m1 <- rowMeans(logsc[, cond == "JNKi", drop = FALSE])
d0 <- logsc[, cond == "DMSO", drop = FALSE]
d1 <- logsc[, cond == "JNKi", drop = FALSE]
v0 <- (d0[, 1] - m0)^2 + (d0[, 2] - m0)^2
v1 <- (d1[, 1] - m1)^2 + (d1[, 2] - m1)^2
se <- sqrt(pmax(v0, 0) / 2 + pmax(v1, 0) / 2)
effect <- m1 - m0
z <- effect / pmax(se, 1e-3)
z[!is.finite(z)] <- 0
# Do not segment empty genome; keep only bins with detectable JUN in any sample.
mx <- do.call(pmax, as.list(as.data.frame(score)))
# Deposited CPM tracks have low-level genome-wide background (~0.3 in IgG
# IDR). Keep bins at/above typical JUN-in-peak scale.
keep_signal <- is.finite(mx) & mx >= 0.6

track <- data.frame(
  chrom = as.character(seqnames(bins)),
  start = start(bins) - 1L,
  end = as.integer(end(bins)),
  score = as.numeric(z),
  effect = as.numeric(effect),
  valid = is.finite(z) & keep_signal,
  stringsAsFactors = FALSE
)
message("signal bins ", sum(track$valid), " / ", nrow(track))

# Segmentation per chromosome
penalties <- c(12, 24)
min_width <- 3L
segs <- EpiSegDiff:::segment_track_penalty_path(
  track, penalties = penalties, min_width = min_width,
  loss = "gaussian", engine = "pelt"
)
cand <- select_candidate_segments(segs, "single_selected_penalty",
                                  selected_penalty = 24)
message("candidates ", nrow(cand))

# Region-level test: mean log-score per sample, two-group t / limma-like
region_means <- function(intervals, logsc, track) {
  gr_i <- GenomicRanges::GRanges(
    intervals$chrom,
    IRanges::IRanges(as.integer(intervals$start) + 1L, as.integer(intervals$end))
  )
  gr_b <- GenomicRanges::GRanges(
    track$chrom,
    IRanges::IRanges(as.integer(track$start) + 1L, as.integer(track$end))
  )
  ov <- GenomicRanges::findOverlaps(gr_i, gr_b)
  qh <- S4Vectors::queryHits(ov)
  sh <- S4Vectors::subjectHits(ov)
  out <- matrix(NA_real_, nrow(intervals), ncol(logsc),
                dimnames = list(NULL, colnames(logsc)))
  if (!length(qh)) return(out)
  spl <- split(sh, qh)
  for (nm in names(spl)) {
    i <- as.integer(nm)
    out[i, ] <- colMeans(logsc[spl[[nm]], , drop = FALSE])
  }
  out
}
if (!nrow(cand)) stop("No candidate segments.")
if (!"region_id" %in% names(cand)) cand$region_id <- seq_len(nrow(cand))
rmat <- region_means(cand, logsc, track)
design <- stats::model.matrix(~ cond)
fit <- if (requireNamespace("limma", quietly = TRUE)) {
  lmfit <- limma::lmFit(rmat, design)
  lmfit <- limma::eBayes(lmfit, trend = TRUE, robust = TRUE)
  tt <- limma::topTable(lmfit, coef = 2, number = Inf, sort.by = "none")
  data.frame(logFC = tt$logFC, statistic = tt$t, pvalue = tt$P.Value,
             FDR = tt$adj.P.Val, stringsAsFactors = FALSE)
} else {
  p <- apply(rmat, 1, function(x) {
    a <- x[cond == "DMSO"]; b <- x[cond == "JNKi"]
    if (any(!is.finite(c(a, b)))) return(NA_real_)
    stats::t.test(b, a, var.equal = TRUE)$p.value
  })
  data.frame(
    logFC = rowMeans(rmat[, cond == "JNKi", drop = FALSE]) -
      rowMeans(rmat[, cond == "DMSO", drop = FALSE]),
    statistic = NA_real_, pvalue = p,
    FDR = stats::p.adjust(p, "BH"),
    stringsAsFactors = FALSE
  )
}

tested <- cbind(
  region_id = cand$region_id,
  chrom = cand$chrom, start = cand$start, end = cand$end,
  width = cand$end - cand$start,
  n_bins = cand$n_bins,
  mean_z = cand$effect_mean,
  fit,
  observation = "deposited_bigwig_spikein_CPM_proxy",
  testing_method = "limma_or_t_on_log2_bigwig",
  contrast = "JNKi_minus_DMSO",
  stringsAsFactors = FALSE
)
tested$direction <- ifelse(tested$logFC > 0, "GAIN",
                           ifelse(tested$logFC < 0, "LOSS", "STABLE"))
tested$level1 <- ifelse(
  !is.na(tested$FDR) & tested$FDR <= 0.05 & is.finite(tested$logFC),
  tested$direction, "STABLE"
)
sig <- tested[tested$level1 %in% c("GAIN", "LOSS"), , drop = FALSE]

# A/B/C vs IDR presence (not fragment DE)
sig_gr <- GenomicRanges::GRanges(
  sig$chrom, IRanges::IRanges(sig$start + 1L, sig$end)
)
ov_idr <- GenomicRanges::countOverlaps(sig_gr, idr_u) > 0
class_c <- sig
class_c$abc <- ifelse(ov_idr, "A_shared_IDR", "C_episegdiff_only")
# IDR intervals with no EpiSegDiff overlap = B
idr_df <- data.frame(chrom = as.character(seqnames(idr_u)),
                     start = start(idr_u) - 1L, end = as.integer(end(idr_u)))
if (nrow(sig)) {
  ov_b <- GenomicRanges::countOverlaps(idr_u, sig_gr) == 0
} else {
  ov_b <- rep(TRUE, length(idr_u))
}

utils::write.table(tested, "results/episegdiff_primary.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(
  sig[, c("chrom", "start", "end", "region_id", "level1", "logFC", "FDR")],
  "results/episegdiff_primary.bed",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
abc <- data.frame(
  class = c("A_shared", "B_IDR_only", "C_EpiSegDiff_only"),
  n = c(sum(class_c$abc == "A_shared_IDR"), sum(ov_b),
        sum(class_c$abc == "C_episegdiff_only")),
  note = "IDR is presence reference; EpiSegDiff is bigwig-proxy Level I",
  stringsAsFactors = FALSE
)
utils::write.table(abc, "results/region_comparison.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(class_c, "data/GSE280947/episegdiff/level1_sig_annotated.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

n_gain <- sum(sig$level1 == "GAIN")
n_loss <- sum(sig$level1 == "LOSS")
elapsed <- proc.time()[["elapsed"]] - t_start
man <- c(
  observation = "deposited_bigwig_spikein_CPM_NOT_fragment_counts",
  bin_size = bin_size,
  n_bins = length(bins),
  n_candidates = nrow(cand),
  n_sig = nrow(sig),
  n_GAIN = n_gain,
  n_LOSS = n_loss,
  median_sig_width = if (nrow(sig)) median(sig$width) else NA,
  n_classA = abc$n[1], n_classB = abc$n[2], n_classC = abc$n[3],
  runtime_sec = elapsed,
  contrast = "JNKi_minus_DMSO"
)
write_run_manifest("results/run_manifest.tsv", as.list(man))

message("GAIN ", n_gain, " LOSS ", n_loss, " A/B/C ",
        paste(abc$n, collapse = "/"), " time ", round(elapsed), "s")
