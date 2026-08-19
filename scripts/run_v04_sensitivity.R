#!/usr/bin/env Rscript
# Compact GSE280947 sensitivity: one import, then transform/filter/penalty grid.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(20271)
cfg <- read_continuous_config("configs/gse280947_bigwig.yaml")
sheet <- read_bigwig_sample_sheet(cfg$samples)
assay <- read_bigwig_assay(sheet, genome = "mm10", bin_size = 200L)
sig <- assay$signal
bed <- granges_to_bed_frame(assay$bins)
num <- "JNKi"; den <- "DMSO"
ref <- utils::read.table("results/gse280947_v04/episegdiff_primary.tsv",
                         header = TRUE, sep = "\t")
ref_sig <- ref[ref$level1 %in% c("GAIN", "LOSS"), ]
ref_gr <- GenomicRanges::GRanges(ref_sig$chrom,
                                 IRanges::IRanges(ref_sig$start + 1L, ref_sig$end))

jaccard <- function(a, b) {
  if (!length(a) && !length(b)) return(1)
  if (!length(a) || !length(b)) return(0)
  inter <- sum(GenomicRanges::countOverlaps(a, b) > 0)
  union <- length(a) + length(b) - inter
  inter / max(union, 1)
}

rows <- list(); k <- 0L
# Do not run filter=none genome-wide (13M bins). Compare adaptive vs
# the GSE280947 profile threshold 0.6.
grid <- expand.grid(
  pseudo = c(0.001, 0.01, 0.05, 0.1),
  filter = c("absolute"),
  penalty = c(24),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
grid <- rbind(grid,
              data.frame(pseudo = 0.01, filter = "quantile", penalty = 24),
              data.frame(pseudo = 0.01, filter = "auto", penalty = 24),
              data.frame(pseudo = 0.01, filter = "absolute", penalty = 12),
              data.frame(pseudo = 0.01, filter = "absolute", penalty = 48))

for (i in seq_len(nrow(grid))) {
  t0 <- proc.time()[["elapsed"]]
  filt <- filter_signal_bins(sig, method = grid$filter[i],
                             threshold = 0.6, quantile = 0.75)
  xf <- transform_continuous(sig, "log2_pseudocount", grid$pseudo[i])
  eff <- estimate_continuous_effect(xf, sheet, num, den, "replicated")
  z <- eff$effect / pmax(abs(eff$se), 1e-6)
  z[!is.finite(z)] <- 0
  track <- data.frame(chrom = bed$chrom, start = bed$start, end = bed$end,
                      score = z, valid = filt & is.finite(z))
  segs <- EpiSegDiff:::segment_track(track, penalty = grid$penalty[i],
                                     min_width = 3L, loss = "gaussian", engine = "pelt")
  if (!nrow(segs)) {
    ng <- nl <- 0L; medw <- p90 <- NA; jac <- NA; signc <- NA
  } else {
    segs$region_id <- seq_len(nrow(segs))
    tested <- test_regions_limma(segs, xf, track, sheet, num, den, "replicated", 0.05)
    sigc <- tested[tested$level1 %in% c("GAIN", "LOSS"), ]
    ng <- sum(sigc$level1 == "GAIN"); nl <- sum(sigc$level1 == "LOSS")
    w <- sigc$end - sigc$start
    medw <- if (nrow(sigc)) stats::median(w) else NA
    p90 <- if (nrow(sigc)) as.numeric(stats::quantile(w, 0.9)) else NA
    if (nrow(sigc)) {
      g <- GenomicRanges::GRanges(sigc$chrom, IRanges::IRanges(sigc$start + 1L, sigc$end))
      jac <- jaccard(g, ref_gr)
      # sign concordance on overlapping ref
      ov <- GenomicRanges::findOverlaps(g, ref_gr)
      if (length(ov)) {
        signc <- mean(sign(sigc$log2FC[S4Vectors::queryHits(ov)]) ==
                        sign(ref_sig$log2FC[S4Vectors::subjectHits(ov)]))
      } else signc <- NA
    } else {
      jac <- 0; signc <- NA
    }
  }
  k <- k + 1L
  rows[[k]] <- data.frame(
    bin_size = 200, pseudocount = grid$pseudo[i], filter = grid$filter[i],
    penalty = grid$penalty[i], n_tested = if (exists("tested")) nrow(tested) else 0,
    n_GAIN = ng, n_LOSS = nl, median_width = medw, p90_width = p90,
    jaccard_vs_ref = jac, sign_concordance = signc,
    runtime_sec = proc.time()[["elapsed"]] - t0
  )
  rm(tested)
}
# 500 bp coarsen, reference transform/filter
cs <- coarsen_signal_grid(assay$bins, sig, factor = 2L)
bed5 <- granges_to_bed_frame(cs$bins)
filt <- filter_signal_bins(cs$signal, "absolute", 0.6)
xf <- transform_continuous(cs$signal, "log2_pseudocount", 0.01)
eff <- estimate_continuous_effect(xf, sheet, num, den, "replicated")
z <- eff$effect / pmax(abs(eff$se), 1e-6); z[!is.finite(z)] <- 0
track <- data.frame(chrom = bed5$chrom, start = bed5$start, end = bed5$end,
                    score = z, valid = filt & is.finite(z))
t0 <- proc.time()[["elapsed"]]
segs <- EpiSegDiff:::segment_track(track, 24, 3L, loss = "gaussian", engine = "pelt")
segs$region_id <- seq_len(nrow(segs))
tested <- test_regions_limma(segs, xf, track, sheet, num, den, "replicated", 0.05)
sigc <- tested[tested$level1 %in% c("GAIN", "LOSS"), ]
rows[[length(rows) + 1L]] <- data.frame(
  bin_size = 400, # 200*2 coarsen
  pseudocount = 0.01, filter = "absolute", penalty = 24,
  n_tested = nrow(tested),
  n_GAIN = sum(sigc$level1 == "GAIN"), n_LOSS = sum(sigc$level1 == "LOSS"),
  median_width = median(sigc$end - sigc$start),
  p90_width = as.numeric(quantile(sigc$end - sigc$start, 0.9)),
  jaccard_vs_ref = NA, sign_concordance = NA,
  runtime_sec = proc.time()[["elapsed"]] - t0
)
tab <- do.call(rbind, rows)
dir.create("results/gse280947_v04", showWarnings = FALSE)
utils::write.table(tab, "results/gse280947_v04/sensitivity.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
print(tab)
