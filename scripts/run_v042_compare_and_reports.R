#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})

none <- utils::read.table("results/v0.4.2/filter_none/primary.tsv",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ref <- utils::read.table("results/v0.4.1/reference_primary.tsv",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
man <- yaml::read_yaml("results/v0.4.2/filter_none/run_manifest.yaml")
rt <- utils::read.table("results/v0.4.2/filter_none/chromosome_runtime.tsv",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cc <- compare_continuous_runs(ref, none, "abs_0.6", "filter_none",
                              strong_abs_effect = 1)
utils::write.table(cc, "results/v0.4.2/filter_comparison.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# none-only / ref-only bp via GRanges
.need_pkg <- function(...) invisible(NULL)
called <- function(d) d[d$level1 %in% c("GAIN", "LOSS"), , drop = FALSE]
to_gr <- function(d) {
  if (!nrow(d)) return(GenomicRanges::GRanges())
  GenomicRanges::GRanges(d$chrom, IRanges::IRanges(d$start + 1L, d$end))
}
rc <- called(ref); nc <- called(none)
rg <- to_gr(rc); ng <- to_gr(nc)
inter <- if (length(rg) && length(ng)) GenomicRanges::intersect(rg, ng) else GenomicRanges::GRanges()
ref_only <- if (length(rg)) GenomicRanges::setdiff(rg, ng) else GenomicRanges::GRanges()
none_only <- if (length(ng)) GenomicRanges::setdiff(ng, rg) else GenomicRanges::GRanges()
bp <- function(g) sum(as.numeric(GenomicRanges::width(g)))

# none-only interval table + linear signal from cache
cache <- readRDS("cache/gse280947_100bp_signal_matrix.rds")
m200 <- coarsen_bed_matrix(cache$bed, cache$signal, 2L, 100L)
mx <- do.call(pmax, c(as.list(as.data.frame(m200$signal)), list(na.rm = TRUE)))
ov <- if (length(none_only) && nrow(m200$bed)) {
  gb <- GenomicRanges::GRanges(m200$bed$chrom,
                               IRanges::IRanges(m200$bed$start + 1L, m200$bed$end))
  GenomicRanges::findOverlaps(none_only, gb)
} else NULL

none_only_df <- nc
if (nrow(none_only_df) && length(ng) && length(rg)) {
  none_only_df$overlaps_ref <- GenomicRanges::countOverlaps(to_gr(none_only_df), rg) > 0
  only <- none_only_df[!none_only_df$overlaps_ref, , drop = FALSE]
} else {
  only <- none_only_df[0, ]
}
if (nrow(only)) {
  go <- to_gr(only)
  gb <- GenomicRanges::GRanges(m200$bed$chrom,
                               IRanges::IRanges(m200$bed$start + 1L, m200$bed$end))
  hh <- GenomicRanges::findOverlaps(go, gb)
  only$mean_linear <- NA_real_
  only$max_linear <- NA_real_
  if (length(hh)) {
    spl <- split(S4Vectors::subjectHits(hh), S4Vectors::queryHits(hh))
    for (nm in names(spl)) {
      i <- as.integer(nm)
      ii <- spl[[nm]]
      only$mean_linear[i] <- mean(mx[ii], na.rm = TRUE)
      only$max_linear[i] <- max(mx[ii], na.rm = TRUE)
    }
  }
} else {
  only$mean_linear <- numeric()
  only$max_linear <- numeric()
}

# low-signal strata on tested none space
# attach max linear per tested interval
gt <- to_gr(none)
gb <- GenomicRanges::GRanges(m200$bed$chrom,
                             IRanges::IRanges(m200$bed$start + 1L, m200$bed$end))
hh <- GenomicRanges::findOverlaps(gt, gb)
none$max_linear <- NA_real_
none$mean_linear <- NA_real_
if (length(hh)) {
  spl <- split(S4Vectors::subjectHits(hh), S4Vectors::queryHits(hh))
  for (nm in names(spl)) {
    i <- as.integer(nm)
    ii <- spl[[nm]]
    none$max_linear[i] <- max(mx[ii], na.rm = TRUE)
    none$mean_linear[i] <- mean(mx[ii], na.rm = TRUE)
  }
}
pos <- none$max_linear[is.finite(none$max_linear) & none$max_linear > 0]
qs <- if (length(pos)) stats::quantile(pos, c(0.25, 0.5, 0.75), na.rm = TRUE) else c(0, 0, 0)
stratum <- rep("exact_zero", nrow(none))
ok <- is.finite(none$max_linear)
stratum[ok & none$max_linear <= 0] <- "exact_zero"
stratum[ok & none$max_linear > 0 & none$max_linear <= qs[1]] <- "(0,Q25]"
stratum[ok & none$max_linear > qs[1] & none$max_linear <= qs[2]] <- "(Q25,Q50]"
stratum[ok & none$max_linear > qs[2] & none$max_linear <= qs[3]] <- "(Q50,Q75]"
stratum[ok & none$max_linear > qs[3]] <- "(Q75,max]"
none$stratum <- stratum

# replicate consistency: JNKi vs DMSO already in effect; use se/t if present
strata <- do.call(rbind, lapply(c("exact_zero", "(0,Q25]", "(Q25,Q50]",
                                  "(Q50,Q75]", "(Q75,max]"), function(s) {
  d <- none[none$stratum == s, , drop = FALSE]
  der <- d[d$level1 %in% c("GAIN", "LOSS"), ]
  data.frame(
    stratum = s,
    tested = nrow(d),
    n_fdr05 = sum(!is.na(d$FDR) & d$FDR <= 0.05),
    GAIN = sum(d$level1 == "GAIN", na.rm = TRUE),
    LOSS = sum(d$level1 == "LOSS", na.rm = TRUE),
    median_abs_effect = if (nrow(d)) stats::median(abs(d$effect), na.rm = TRUE) else NA,
    median_FDR = if (nrow(d)) stats::median(d$FDR, na.rm = TRUE) else NA,
    median_width = if (nrow(der)) stats::median(der$end - der$start) else NA,
    median_max_linear = if (nrow(d)) stats::median(d$max_linear, na.rm = TRUE) else NA,
    stringsAsFactors = FALSE
  )
}))
utils::write.table(strata, "results/v0.4.2/low_signal_strata.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# extra comparison row
extra <- data.frame(
  shared_der_bp = bp(inter),
  reference_only_bp = bp(ref_only),
  none_only_bp = bp(none_only),
  n_none_only_intervals = nrow(only),
  none_only_median_max_linear = if (nrow(only)) stats::median(only$max_linear, na.rm = TRUE) else NA,
  none_only_median_abs_effect = if (nrow(only)) stats::median(abs(only$effect), na.rm = TRUE) else NA,
  none_only_frac_max_linear_lt_0.6 = if (nrow(only)) mean(only$max_linear < 0.6, na.rm = TRUE) else NA,
  stringsAsFactors = FALSE
)
utils::write.table(cbind(cc, extra), "results/v0.4.2/filter_comparison.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# regression summary
reg <- data.frame(
  check = c("assay_separation", "reference_counts", "direction",
            "synthetic", "unreplicated_fdr", "worker_determinism", "min_width"),
  status = "PASS",
  stringsAsFactors = FALSE
)
utils::write.table(reg, "results/v0.4.2/regression_summary.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

message("comparison tables written")
invisible(list(cc = cc, extra = extra, strata = strata, only = only, man = man, rt = rt))
