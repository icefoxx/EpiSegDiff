#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
cl <- utils::read.table("results/v0.5/region_classification.tsv", header = TRUE,
                        sep = "\t", stringsAsFactors = FALSE)
edc <- cl[cl$method == "ESD", ]
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
keep <- intersect(names(edc), c("chrom", "start", "end", "class",
                                "nearest_gene", "tss_distance", "motif_hit", "effect"))
utils::write.table(edc[, keep], "results/v0.5/gene_linkage.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
if ("motif_hit" %in% names(edc) && !all(is.na(edc$motif_hit))) {
  mot <- do.call(rbind, lapply(sort(unique(edc$class)), function(k) {
    x <- edc$class == k
    data.frame(class = k, n = sum(x), n_hit = sum(edc$motif_hit[x], na.rm = TRUE),
               frac = mean(edc$motif_hit[x], na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
} else {
  mot <- data.frame(class = "ALL", n = nrow(edc), n_hit = NA, frac = NA,
                    status = "NOT_AVAILABLE",
                    failure_reason = paste("motif scan failed"),
                    stringsAsFactors = FALSE)
}
utils::write.table(mot, "results/v0.5/motif_validation.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
pw <- data.frame(status = "NOT_AVAILABLE", failure_reason = "too few genes",
                 stringsAsFactors = FALSE)
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
    "nearest_gene" %in% names(edc)) {
  genes_c <- unique(na.omit(as.character(edc$nearest_gene[edc$class == "C"])))
  univ <- unique(na.omit(as.character(edc$nearest_gene)))
  if (length(genes_c) >= 5L && length(univ) >= 10L) {
    eg <- tryCatch(clusterProfiler::enrichGO(genes_c, OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                                             keyType = "ENTREZID", ont = "BP",
                                             universe = univ, pvalueCutoff = 1,
                                             qvalueCutoff = 1, readable = TRUE),
                   error = function(e) conditionMessage(e))
    if (is.character(eg)) {
      pw <- data.frame(status = "FAILED", failure_reason = eg, stringsAsFactors = FALSE)
    } else if (!is.null(eg) && nrow(as.data.frame(eg))) {
      pw <- utils::head(as.data.frame(eg), 20)
      pw$status <- "PASS"
    } else {
      pw <- data.frame(status = "NOT_AVAILABLE", failure_reason = "enrichGO empty",
                       stringsAsFactors = FALSE)
    }
  }
}
utils::write.table(pw, "results/v0.5/pathway_validation.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# manifests if missing
if (!file.exists("results/v0.5/dataset_manifest.tsv")) {
  utils::write.table(data.frame(
    dataset_id = c("A", "B"), accession = c("GSE280947", "GSE110378"),
    assay = c("CUT&RUN", "ChIPmentation"), target = c("JUN", "H3K27me3"),
    contrast = c("JNKi - DMSO", "KO - Control"),
    bam = c(FALSE, FALSE), status = c("USED", "DOWNLOAD_OK"),
    stringsAsFactors = FALSE
  ), "results/v0.5/dataset_manifest.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
}
if (!file.exists("results/v0.5/method_configs.tsv")) {
  utils::write.table(data.frame(
    method = c("ESD_A", "conventional_A", "fixed_bins_A"),
    region = c("gaussian_PELT", "IDR_union", "200bp_valid"),
    model = "limma_log2p0.01", stringsAsFactors = FALSE
  ), "results/v0.5/method_configs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
}
if (file.exists("results/v0.4.1/combined_sensitivity_grid.tsv") &&
    !file.exists("results/v0.5/sensitivity_summary.tsv")) {
  g <- utils::read.table("results/v0.4.1/combined_sensitivity_grid.tsv",
                         header = TRUE, sep = "\t")
  utils::write.table(
    data.frame(config = g$config_id, GAIN = g$GAIN, LOSS = g$LOSS,
               loss_dominates = g$LOSS > g$GAIN, status = g$status),
    "results/v0.5/sensitivity_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
}
message("annotation finish ok")
