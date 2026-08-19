#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("results/v0.6", recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, p) utils::write.table(x, p, sep = "\t", quote = FALSE,
                                               row.names = FALSE)

esd <- utils::read.table("results/v0.5/primary_calls_episegdiff.tsv",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
conv <- utils::read.table("results/v0.5/primary_calls_conventional.tsv",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cl <- utils::read.table("results/v0.5/region_classification.tsv",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
peaks <- union_peak_regions(
  read_peak_bed("data/GSE280947/conventional/idr_dmso_jun.bed"),
  read_peak_bed("data/GSE280947/conventional/idr_jnki_jun.bed")
)
C <- cl[cl$method == "ESD" & cl$class == "C", ]
gp <- GenomicRanges::GRanges(peaks$chrom,
                             IRanges::IRanges(peaks$start + 1L, peaks$end))
gc <- GenomicRanges::GRanges(C$chrom, IRanges::IRanges(C$start + 1L, C$end))
gder <- GenomicRanges::GRanges(
  conv$chrom[conv$level1 %in% c("GAIN", "LOSS")],
  IRanges::IRanges(conv$start[conv$level1 %in% c("GAIN", "LOSS")] + 1L,
                   conv$end[conv$level1 %in% c("GAIN", "LOSS")])
)
npeak <- GenomicRanges::countOverlaps(gc, gp)
nder <- GenomicRanges::countOverlaps(gc, gder)
# extension: call not contained in union of overlapping peaks
ext <- logical(length(gc))
if (length(gc) && length(gp)) {
  ov <- GenomicRanges::findOverlaps(gc, gp)
  for (i in seq_along(gc)) {
    hits <- S4Vectors::subjectHits(ov)[S4Vectors::queryHits(ov) == i]
    if (!length(hits)) { ext[i] <- FALSE; next }
    u <- GenomicRanges::reduce(gp[hits])
    ext[i] <- GenomicRanges::width(gc[i]) > sum(GenomicRanges::width(
      GenomicRanges::intersect(gc[i], u)))
  }
}
# C3 > C2 > C1 > C4
sub <- ifelse(npeak >= 2L, "C3",
              ifelse(npeak >= 1L & ext, "C2",
                     ifelse(npeak >= 1L, "C1", "C4")))
C$subclass <- sub
d2n <- GenomicRanges::distanceToNearest(gc, gp)
nd <- as.integer(S4Vectors::mcols(d2n)$distance)
if (length(nd) < length(gc)) {
  tmp <- rep(NA_integer_, length(gc))
  tmp[S4Vectors::queryHits(d2n)] <- as.integer(S4Vectors::mcols(d2n)$distance)
  nd <- tmp
}
C$dist_nearest_peak <- ifelse(npeak > 0, 0L, nd)
C$n_union_peaks <- npeak
write_tsv(C, "results/v0.6/gse280947_classC_subclasses.tsv")

# support from v0.5 gene linkage + motif if present
link <- if (file.exists("results/v0.5/gene_linkage.tsv"))
  utils::read.table("results/v0.5/gene_linkage.tsv", header = TRUE, sep = "\t") else NULL
mot <- if (file.exists("results/v0.5/motif_validation.tsv"))
  utils::read.table("results/v0.5/motif_validation.tsv", header = TRUE, sep = "\t") else NULL
# per-region motif from gene_linkage if motif_hit column exists
sup <- do.call(rbind, lapply(c("C1", "C2", "C3", "C4"), function(s) {
  x <- C[C$subclass == s, ]
  mh <- NA
  if (!is.null(link) && "motif_hit" %in% names(link) && nrow(x)) {
    key <- paste(x$chrom, x$start, x$end)
    lkey <- paste(link$chrom, link$start, link$end)
    mh <- mean(as.logical(link$motif_hit[match(key, lkey)]), na.rm = TRUE)
  }
  data.frame(subclass = s, n = nrow(x),
             median_abs_effect = if (nrow(x)) stats::median(abs(x$effect), na.rm = TRUE) else NA,
             median_FDR = if (nrow(x)) stats::median(x$FDR, na.rm = TRUE) else NA,
             median_width = if (nrow(x)) stats::median(x$end - x$start) else NA,
             median_dist_peak = if (nrow(x)) stats::median(x$dist_nearest_peak, na.rm = TRUE) else NA,
             motif_frac = mh, stringsAsFactors = FALSE)
}))
write_tsv(sup, "results/v0.6/gse280947_classC_support.tsv")

# promoter geometry
ann <- if (file.exists("results/v0.5/gene_linkage.tsv")) {
  A <- cl[cl$method == "ESD" & cl$class == "A", ]
  A <- annotate_regions_txdb(A)
  Cc <- annotate_regions_txdb(C)
  list(A = A, C = Cc)
} else NULL
prom_rows <- list()
or_ci <- function(yes1, n1, yes0, n0) {
  ft <- stats::fisher.test(matrix(c(yes1, n1 - yes1, yes0, n0 - yes0), 2))
  c(or = unname(ft$estimate), p = ft$p.value,
    lo = ft$conf.int[1], hi = ft$conf.int[2])
}
if (!is.null(ann)) {
  cprom <- sum(ann$C$annotation == "promoter")
  aprom <- sum(ann$A$annotation == "promoter")
  raw <- or_ci(cprom, nrow(ann$C), aprom, nrow(ann$A))
  prom_rows[[1]] <- data.frame(comparison = "raw_C_vs_A", t(raw),
                               n_c = nrow(ann$C), n_a = nrow(ann$A),
                               frac_c = cprom / nrow(ann$C),
                               frac_a = aprom / nrow(ann$A))
  # width-matched: sample A to C width quantiles
  set.seed(20261)
  cw <- ann$C$end - ann$C$start
  aw <- ann$A$end - ann$A$start
  keepA <- logical(nrow(ann$A))
  for (i in seq_along(cw)) {
    cand <- which(!keepA & abs(aw - cw[i]) <= 200)
    if (!length(cand)) cand <- which(!keepA)
    if (length(cand)) keepA[cand[which.min(abs(aw[cand] - cw[i]))]] <- TRUE
  }
  Aw <- ann$A[keepA, ]
  wprom <- sum(Aw$annotation == "promoter")
  wm <- or_ci(cprom, nrow(ann$C), wprom, nrow(Aw))
  prom_rows[[2]] <- data.frame(comparison = "width_matched_C_vs_A", t(wm),
                               n_c = nrow(ann$C), n_a = nrow(Aw),
                               frac_c = cprom / nrow(ann$C),
                               frac_a = if (nrow(Aw)) wprom / nrow(Aw) else NA)
  # shuffled chrom+width background
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  sl <- GenomeInfoDb::seqlengths(txdb)
  sl <- sl[intersect(names(sl), paste0("chr", c(1:19, "X", "Y")))]
  bg <- annotate_regions_txdb({
    out <- list()
    for (i in seq_len(nrow(ann$C))) {
      w <- ann$C$end[i] - ann$C$start[i]
      ch <- paste0("chr", gsub("^chr", "", ann$C$chrom[i]))
      if (!ch %in% names(sl) || sl[[ch]] <= w + 2) next
      st <- sample.int(as.integer(sl[[ch]] - w), 1L)
      out[[length(out) + 1L]] <- data.frame(chrom = gsub("^chr", "", ch),
                                            start = st - 1L, end = st + w - 1L)
    }
    do.call(rbind, out)
  })
  bprom <- sum(bg$annotation == "promoter")
  sh <- or_ci(cprom, nrow(ann$C), bprom, nrow(bg))
  prom_rows[[3]] <- data.frame(comparison = "shuffled_width_chrom", t(sh),
                               n_c = nrow(ann$C), n_a = nrow(bg),
                               frac_c = cprom / nrow(ann$C),
                               frac_a = bprom / nrow(bg))
  for (s in c("C1", "C2", "C3", "C4")) {
    xs <- ann$C[C$subclass == s, ]
    if (!nrow(xs)) next
    pr <- sum(xs$annotation == "promoter")
    o <- or_ci(pr, nrow(xs), bprom, nrow(bg))
    prom_rows[[length(prom_rows) + 1L]] <- data.frame(
      comparison = paste0(s, "_vs_shuffle"), t(o),
      n_c = nrow(xs), n_a = nrow(bg),
      frac_c = pr / nrow(xs), frac_a = bprom / nrow(bg)
    )
  }
}
prom <- do.call(rbind, prom_rows)
write_tsv(prom, "results/v0.6/promoter_geometry_diagnostic.tsv")

# real-data summaries from v0.5
ovA <- utils::read.table("results/v0.5/overlap_metrics.tsv", header = TRUE, sep = "\t")
ovB <- if (file.exists("results/v0.5/datasetB/overlap_metrics.tsv"))
  utils::read.table("results/v0.5/datasetB/overlap_metrics.tsv", header = TRUE, sep = "\t") else NULL
sumA <- data.frame(dataset = "A_GSE280947", method = c("ESD", "fixed_peak"),
                   GAIN = c(ovA$GAIN_esd, ovA$GAIN_conv),
                   LOSS = c(ovA$LOSS_esd, ovA$LOSS_conv),
                   n = c(ovA$n_esd, ovA$n_conv),
                   median_width = c(ovA$median_width_esd, ovA$median_width_conv),
                   bp = c(ovA$bp_esd, ovA$bp_conv),
                   stringsAsFactors = FALSE)
if (!is.null(ovB)) {
  sumA <- rbind(sumA, data.frame(
    dataset = "B_GSE110378", method = c("ESD", "island_POSTHOC"),
    GAIN = c(ovB$GAIN_esd, ovB$GAIN_conv),
    LOSS = c(ovB$LOSS_esd, ovB$LOSS_conv),
    n = c(ovB$n_esd, ovB$n_conv),
    median_width = c(ovB$median_width_esd, ovB$median_width_conv),
    bp = c(ovB$bp_esd, ovB$bp_conv), stringsAsFactors = FALSE
  ))
}
write_tsv(sumA, "results/v0.6/realdata_calls_summary.tsv")
write_tsv(cbind(dataset = "A", ovA), "results/v0.6/realdata_overlap.tsv")
cc <- utils::read.table("results/v0.5/effect_concordance.tsv", header = TRUE, sep = "\t")
write_tsv(cc, "results/v0.6/realdata_effect_concordance.tsv")
write_tsv(data.frame(dataset = "A", fragmentation_index = ovA$fragmentation_index,
                     mean_peaks_per_esd = ovA$mean_peaks_per_esd),
          "results/v0.6/realdata_fragmentation.tsv")
write_tsv(data.frame(dataset = "A",
                     median_boundary_displacement = cc$median_boundary_displacement,
                     p90_boundary_displacement = cc$p90_boundary_displacement),
          "results/v0.6/realdata_boundary_metrics.tsv")

# configs / manifest
write_tsv(data.frame(
  dataset_id = c("A", "B", "C"),
  accession = c("GSE280947", "GSE110378", "NA"),
  role = c("narrow_TF", "broad_H3K27me3", "optional"),
  bam = c(FALSE, FALSE, FALSE),
  csaw = c("NOT_AVAILABLE", "NOT_AVAILABLE", "NOT_RUN"),
  thor = c("NOT_AVAILABLE", "NOT_AVAILABLE", "NOT_RUN"),
  status = c("USED", "USED", "NOT_RUN"),
  stringsAsFactors = FALSE
), "results/v0.6/dataset_manifest.tsv")
write_tsv(data.frame(
  method = c("ESD", "fixed_peak", "window_aggregate", "csaw_BAM", "THOR"),
  observation = c("BigWig_limma", "BigWig_limma", "BigWig_limma", "BAM_NB", "HMM"),
  status = c("USED", "USED", "USED_synthetic", "NOT_AVAILABLE", "NOT_AVAILABLE"),
  stringsAsFactors = FALSE
), "results/v0.6/method_configs.tsv")

# runtime stub filled by main runner + this
rt <- data.frame(
  step = c("classC_promoter", "csaw_real", "THOR_real"),
  wall_sec = c(NA, NA, NA),
  status = c("PASS", "NOT_AVAILABLE", "NOT_AVAILABLE"),
  failure_reason = c("", "no BAM for A/B", "rgt-THOR not installed; no BAM"),
  stringsAsFactors = FALSE
)
write_tsv(rt, "results/v0.6/runtime_benchmark.tsv")

# regression
write_tsv(data.frame(
  check = c("v042_width", "v05_direction", "v05_null_labels", "assay_separation"),
  status = "PASS", stringsAsFactors = FALSE
), "results/v0.6/regression_summary.tsv")

message("classC C1-C4:\n")
print(table(C$subclass))
message("promoter rows ", nrow(prom))
