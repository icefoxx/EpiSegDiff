#!/usr/bin/env Rscript
# Deposited MACS2/IDR peaks as a conventional *presence* reference.
# This is not DiffBind/edgeR on fragments.

dmso_f <- "data/GSE280947/conventional/GSM8607877_DMSO2_cJun_sorted.bed_idr.txt.gz"
jnki_f <- "data/GSE280947/conventional/GSM8607883_JNKI2_cJun_sorted.bed_idr.txt.gz"
read_idr <- function(path) {
  x <- utils::read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  data.frame(chrom = as.character(x[[1]]), start = as.integer(x[[2]]),
             end = as.integer(x[[3]]), stringsAsFactors = FALSE)
}
d <- read_idr(dmso_f)
j <- read_idr(jnki_f)
# crude reciprocal-overlap match
ro <- function(a, b, thr = 0.2) {
  n <- 0L
  for (i in seq_len(nrow(a))) {
    same <- b$chrom == a$chrom[i]
    if (!any(same)) next
    bb <- b[same, ]
    inter <- pmax(0L, pmin(a$end[i], bb$end) - pmax(a$start[i], bb$start))
    rec <- pmin(inter / pmax(a$end[i] - a$start[i], 1), inter / pmax(bb$end - bb$start, 1))
    if (any(rec >= thr)) n <- n + 1L
  }
  n
}
n_shared_from_d <- ro(d, j)
n_shared_from_j <- ro(j, d)
dir.create("results", showWarnings = FALSE)
out <- data.frame(
  reference = "authors_MACS2_IDR_presence",
  n_dmso_idr = nrow(d),
  n_jnki_idr = nrow(j),
  n_dmso_with_jnki_overlap = n_shared_from_d,
  n_jnki_with_dmso_overlap = n_shared_from_j,
  n_dmso_only_approx = nrow(d) - n_shared_from_d,
  n_jnki_only_approx = nrow(j) - n_shared_from_j,
  median_width_dmso = median(d$end - d$start),
  median_width_jnki = median(j$end - j$start),
  note = "presence-absence of IDR peaks, not fragment DE; not EpiSegDiff",
  stringsAsFactors = FALSE
)
utils::write.table(out, "results/conventional_differential.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(d, "data/GSE280947/conventional/idr_dmso_jun.bed",
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
utils::write.table(j, "data/GSE280947/conventional/idr_jnki_jun.bed",
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
print(out)
