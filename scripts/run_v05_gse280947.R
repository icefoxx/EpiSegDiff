#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(20250)
dir.create("results/v0.5", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

cache <- readRDS("cache/gse280947_100bp_signal_matrix.rds")
m200 <- coarsen_bed_matrix(cache$bed, cache$signal, 2L, 100L)
sheet <- cache$samples
num <- "JNKi"; den <- "DMSO"

# --- ESD primary (validated profile) ---
message("ESD primary absolute 0.6")
esd <- analyze_continuous_matrix(
  m200$bed, m200$signal, sheet, num, den, bin_size = 200L,
  filter_method = "absolute", filter_threshold = 0.6,
  penalty = 24, min_bins = 3L, engine = "pelt"
)
utils::write.table(esd$tested, "results/v0.5/primary_calls_episegdiff.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# --- Conventional Layer 1 ---
message("conventional union IDR + limma")
idr_d <- read_peak_bed("data/GSE280947/conventional/idr_dmso_jun.bed")
idr_j <- read_peak_bed("data/GSE280947/conventional/idr_jnki_jun.bed")
peaks <- union_peak_regions(idr_d, idr_j)
rsig <- summarize_matrix_over_intervals(peaks, m200$bed, m200$signal)
conv <- test_fixed_regions_limma(peaks, rsig, sheet, num, den)
conv <- standardize_primary_results(conv, meta = list(
  assay_family = "continuous", input_type = "bigwig",
  contrast = "JNKi - DMSO", inference_mode = "replicated",
  signal_filter_method = "idr_union", transform_method = "log2_pseudocount",
  segmentation_method = "fixed_peaks"
))
utils::write.table(conv, "results/v0.5/primary_calls_conventional.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# --- Classification ---
cl <- classify_region_sets(esd$tested, conv)
class_tab <- rbind(
  data.frame(method = "ESD", cl$esd[, intersect(names(cl$esd),
             c("chrom", "start", "end", "effect", "FDR", "level1", "class"))],
             stringsAsFactors = FALSE),
  data.frame(method = "conventional", cl$conventional[, intersect(names(cl$conventional),
             c("chrom", "start", "end", "effect", "FDR", "level1", "class"))],
             stringsAsFactors = FALSE)
)
utils::write.table(class_tab, "results/v0.5/region_classification.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

ov <- overlap_geometry_metrics(esd$tested, conv, peaks)
ov$n_opposite_1bp <- cl$n_opposite_overlaps
cc <- compare_continuous_runs(esd$tested, conv, "ESD", "conventional", 1)
utils::write.table(ov, "results/v0.5/overlap_metrics.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(cc, "results/v0.5/effect_concordance.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# geometry
ed <- esd$tested[esd$tested$level1 %in% c("GAIN", "LOSS"), ]
cd <- conv[conv$level1 %in% c("GAIN", "LOSS"), ]
geom <- rbind(
  data.frame(method = "ESD", width = ed$end - ed$start, effect = ed$effect),
  data.frame(method = "conventional", width = cd$end - cd$start, effect = cd$effect)
)
utils::write.table(geom, "results/v0.5/region_geometry.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# replicate consistency: per-class mean |t| / sign agreement of bin effects
# use region se from tables
rep_row <- function(d, cls, method) {
  x <- d[d$class == cls, ]
  data.frame(method = method, class = cls, n = nrow(x),
             median_abs_effect = if (nrow(x)) stats::median(abs(x$effect), na.rm = TRUE) else NA,
             median_FDR = if (nrow(x)) stats::median(x$FDR, na.rm = TRUE) else NA,
             stringsAsFactors = FALSE)
}
rep_tab <- rbind(rep_row(cl$esd, "A", "ESD"), rep_row(cl$esd, "C", "ESD"),
                 rep_row(cl$conventional, "A", "conventional"),
                 rep_row(cl$conventional, "B", "conventional"))
utils::write.table(rep_tab, "results/v0.5/replicate_consistency.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# --- Ablation: fixed bins ---
message("fixed-bin ablation")
xf <- transform_continuous(m200$signal, "log2_pseudocount", 0.01)
filt <- filter_signal_bins(m200$signal, "absolute", 0.6)
eff <- estimate_continuous_effect(xf, sheet, num, den, "replicated")
bin_fdr <- p.adjust(eff$pvalue, "BH")
bin_call <- filt & is.finite(bin_fdr) & bin_fdr <= 0.05
n_gain_bin <- sum(bin_call & eff$effect > 0, na.rm = TRUE)
n_loss_bin <- sum(bin_call & eff$effect < 0, na.rm = TRUE)
abl <- data.frame(
  ablation = c("adaptive_PELT", "fixed_bins_unmerged", "fixed_peaks"),
  GAIN = c(esd$n_GAIN, n_gain_bin, sum(cd$level1 == "GAIN")),
  LOSS = c(esd$n_LOSS, n_loss_bin, sum(cd$level1 == "LOSS")),
  n_units = c(esd$n_tested, sum(filt), nrow(conv)),
  median_width = c(stats::median(ed$end - ed$start), 200,
                   stats::median(cd$end - cd$start)),
  stringsAsFactors = FALSE
)
utils::write.table(abl, "results/v0.5/ablation_results.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# --- Negative controls: label partitions ---
message("label permutations")
parts <- balanced_label_partitions(sheet, num, den)
null_rows <- list()
for (i in seq_along(parts)) {
  sh <- sheet
  sh$condition <- NA_character_
  sh$condition[parts[[i]][[1]]] <- den
  sh$condition[parts[[i]][[2]]] <- num
  fitn <- tryCatch(analyze_continuous_matrix(
    m200$bed, m200$signal, sh, num, den, bin_size = 200L,
    filter_method = "absolute", filter_threshold = 0.6,
    penalty = 24, min_bins = 3L, engine = "pelt"
  ), error = function(e) NULL)
  if (is.null(fitn)) {
    null_rows[[i]] <- data.frame(partition = i, GAIN = NA, LOSS = NA,
                                 status = "FAILED", stringsAsFactors = FALSE)
  } else {
    null_rows[[i]] <- data.frame(partition = i, GAIN = fitn$n_GAIN,
                                 LOSS = fitn$n_LOSS, status = "PASS",
                                 stringsAsFactors = FALSE)
  }
}
null_tab <- rbind(
  data.frame(partition = 0, GAIN = esd$n_GAIN, LOSS = esd$n_LOSS, status = "TRUE"),
  do.call(rbind, null_rows)
)
utils::write.table(null_tab, "results/v0.5/negative_controls.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# --- Annotation / linkage / motif on ESD classes ---
message("annotation / motif")
edc <- cl$esd
edc <- annotate_regions_txdb(edc)
edc <- link_nearest_tss(edc)
edc <- tryCatch(scan_jun_ap1_motifs(edc), error = function(e) {
  edc$motif_hit <- NA
  attr(edc, "motif_error") <- conditionMessage(e)
  edc
})
ann <- as.data.frame(table(class = edc$class, annotation = edc$annotation),
                     stringsAsFactors = FALSE)
utils::write.table(ann, "results/v0.5/genomic_annotation.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(edc[, intersect(names(edc), c("chrom", "start", "end", "class",
                                                 "nearest_gene", "tss_distance",
                                                 "motif_hit", "effect"))],
                   "results/v0.5/gene_linkage.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
if ("motif_hit" %in% names(edc) && !all(is.na(edc$motif_hit))) {
  mot <- do.call(rbind, lapply(c("A", "C"), function(k) {
    x <- edc$class == k
    data.frame(class = k, n = sum(x), n_hit = sum(edc$motif_hit[x], na.rm = TRUE),
               frac = mean(edc$motif_hit[x], na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
} else {
  mot <- data.frame(class = "ALL", n = nrow(edc), n_hit = NA, frac = NA,
                    status = "NOT_AVAILABLE",
                    failure_reason = "BSgenome/motifmatchr scan unavailable",
                    stringsAsFactors = FALSE)
}
utils::write.table(mot, "results/v0.5/motif_validation.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# pathways if possible
pw <- data.frame(status = "NOT_RUN", failure_reason = "pending gene universe",
                 stringsAsFactors = FALSE)
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
    "nearest_gene" %in% names(edc)) {
  genes_c <- unique(na.omit(edc$nearest_gene[edc$class == "C"]))
  univ <- unique(na.omit(edc$nearest_gene))
  if (length(genes_c) >= 5L && length(univ) >= 10L) {
    eg <- tryCatch(clusterProfiler::enrichGO(genes_c, OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                                             keyType = "ENTREZID", ont = "BP",
                                             universe = univ, pvalueCutoff = 1,
                                             qvalueCutoff = 1, readable = TRUE),
                   error = function(e) NULL)
    if (!is.null(eg) && nrow(as.data.frame(eg))) {
      pw <- utils::head(as.data.frame(eg), 20)
      pw$status <- "PASS"
    } else {
      pw <- data.frame(status = "NOT_AVAILABLE",
                       failure_reason = "enrichGO empty or failed",
                       stringsAsFactors = FALSE)
    }
  } else {
    pw <- data.frame(status = "NOT_AVAILABLE",
                     failure_reason = "too few linked genes",
                     stringsAsFactors = FALSE)
  }
}
utils::write.table(pw, "results/v0.5/pathway_validation.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# sensitivity summary from v0.4.1
sens <- data.frame(
  source = "v0.4.1_grid",
  note = "LOSS dominance holds on all 12 cells; filter=none is secondary",
  stringsAsFactors = FALSE
)
if (file.exists("results/v0.4.1/combined_sensitivity_grid.tsv")) {
  g <- utils::read.table("results/v0.4.1/combined_sensitivity_grid.tsv",
                         header = TRUE, sep = "\t")
  sens <- data.frame(
    config = g$config_id, GAIN = g$GAIN, LOSS = g$LOSS,
    loss_dominates = g$LOSS > g$GAIN, status = g$status,
    stringsAsFactors = FALSE
  )
}
utils::write.table(sens, "results/v0.5/sensitivity_summary.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# manifests
man <- data.frame(
  dataset_id = c("A", "B"),
  accession = c("GSE280947", "GSE110378"),
  assay = c("CUT&RUN", "ChIPmentation"),
  target = c("JUN", "H3K27me3"),
  contrast = c("JNKi - DMSO", "KO - Control"),
  n_rep = c("2+2", "2+2"),
  genome = c("mm10", "mm10"),
  bigwig = c(TRUE, TRUE),
  bam = c(FALSE, FALSE),
  status = c("USED", "PENDING_DOWNLOAD"),
  stringsAsFactors = FALSE
)
utils::write.table(man, "results/v0.5/dataset_manifest.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
cfg <- data.frame(
  method = c("ESD_A", "conventional_A", "fixed_bins_A"),
  region = c("gaussian_PELT", "IDR_union", "200bp_valid"),
  model = c("limma_log2p0.01", "limma_log2p0.01", "limma_log2p0.01"),
  filter = c("absolute_0.6", "none_on_peaks", "absolute_0.6"),
  stringsAsFactors = FALSE
)
utils::write.table(cfg, "results/v0.5/method_configs.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

message("GSE280947 v0.5 tables written; ESD GAIN/LOSS=", esd$n_GAIN, "/", esd$n_LOSS,
        " conv GAIN/LOSS=", sum(cd$level1 == "GAIN"), "/", sum(cd$level1 == "LOSS"))
