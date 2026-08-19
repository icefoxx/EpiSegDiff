#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(11)

plant <- function(n, hits, dir = "GAIN") {
  mat <- matrix(0.12, n, 4)
  colnames(mat) <- c("d1", "d2", "j1", "j2")
  if (is.list(hits)) hits <- unlist(hits)
  if (dir == "GAIN") mat[hits, 3:4] <- 4 else mat[hits, 1:2] <- 4
  mat
}
eval_case <- function(name, n, truth_bins, dir, peak = NULL) {
  bed <- data.frame(chrom = "1", start = (seq_len(n) - 1L) * 50L,
                    end = seq_len(n) * 50L)
  if (dir == "both") {
    sig <- matrix(0.12, n, 4)
    colnames(sig) <- c("d1", "d2", "j1", "j2")
    sig[truth_bins[[1]], 3:4] <- 4
    sig[truth_bins[[2]], 1:2] <- 4
  } else {
    sig <- plant(n, truth_bins, dir)
  }
  samp <- data.frame(sample = colnames(sig),
                     condition = c("DMSO", "DMSO", "JNKi", "JNKi"))
  fit <- analyze_continuous_matrix(bed, sig, samp, "JNKi", "DMSO",
                                   bin_size = 50L, filter_method = "none",
                                   penalty = 5, min_bins = 3L, engine = "pelt")
  der <- fit$tested[fit$tested$level1 %in% c("GAIN", "LOSS"), ]
  tb <- if (is.list(truth_bins)) sort(unique(unlist(truth_bins))) else truth_bins
  tgr <- GenomicRanges::GRanges("1", IRanges::IRanges(min(bed$start[tb]) + 1L,
                                                      max(bed$end[tb])))
  if (is.list(truth_bins)) {
    tgr <- do.call(c, lapply(truth_bins, function(b)
      GenomicRanges::GRanges("1", IRanges::IRanges(min(bed$start[b]) + 1L,
                                                   max(bed$end[b])))))
  }
  dgr <- if (nrow(der)) GenomicRanges::GRanges(der$chrom,
                                               IRanges::IRanges(der$start + 1L, der$end)) else
    GenomicRanges::GRanges()
  inter <- if (length(dgr)) sum(GenomicRanges::width(GenomicRanges::intersect(dgr, tgr))) else 0
  tp_bp <- inter
  pred_bp <- if (length(dgr)) sum(GenomicRanges::width(dgr)) else 0
  tru_bp <- sum(GenomicRanges::width(tgr))
  prec <- if (pred_bp) tp_bp / pred_bp else 0
  rec <- if (tru_bp) tp_bp / tru_bp else 0
  det <- length(dgr) > 0 && any(GenomicRanges::countOverlaps(tgr, dgr) > 0)
  if (length(tgr) > 1) det <- all(GenomicRanges::countOverlaps(tgr, dgr) > 0)
  dir_ok <- if (!nrow(der)) FALSE else {
    if (dir == "GAIN") any(der$direction == "GAIN") else if (dir == "LOSS")
      any(der$direction == "LOSS") else
      any(der$direction == "GAIN") && any(der$direction == "LOSS")
  }
  # boundary error vs first truth block
  b0 <- if (is.list(truth_bins)) truth_bins[[1]] else truth_bins
  terr <- if (nrow(der)) {
    min(abs(der$start - bed$start[min(b0)]) + abs(der$end - bed$end[max(b0)]))
  } else NA
  data.frame(case = name, n_der = nrow(der), precision = prec, recall = rec,
             detected = det, direction_ok = dir_ok, boundary_l1 = terr,
             stringsAsFactors = FALSE)
}

rows <- list(
  eval_case("narrow_GAIN", 60, 20:25, "GAIN"),
  eval_case("narrow_LOSS", 60, 20:25, "LOSS"),
  eval_case("broad_GAIN", 80, 20:39, "GAIN"),
  eval_case("broad_LOSS", 80, 20:39, "LOSS"),
  eval_case("nearby_same", 80, list(15:20, 25:30), "GAIN"),
  eval_case("nearby_opposite", 80, list(15:20, 25:30), "both")
)
# peak-boundary crossing: truth 18:35, planted peak 10:22
cross <- eval_case("cross_peak_boundary", 70, 18:35, "GAIN")
rows[[length(rows) + 1L]] <- cross
tab <- do.call(rbind, rows)
utils::write.table(tab, "results/v0.5/synthetic_benchmark.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
print(tab)
