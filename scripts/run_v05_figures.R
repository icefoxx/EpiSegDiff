#!/usr/bin/env Rscript
dir.create("figures", showWarnings = FALSE)
ov <- utils::read.table("results/v0.5/overlap_metrics.tsv", header = TRUE, sep = "\t")
abl <- utils::read.table("results/v0.5/ablation_results.tsv", header = TRUE, sep = "\t")
geom <- utils::read.table("results/v0.5/region_geometry.tsv", header = TRUE, sep = "\t")
null <- utils::read.table("results/v0.5/negative_controls.tsv", header = TRUE, sep = "\t")
syn <- if (file.exists("results/v0.5/synthetic_benchmark.tsv"))
  utils::read.table("results/v0.5/synthetic_benchmark.tsv", header = TRUE, sep = "\t") else NULL
ann <- if (file.exists("results/v0.5/genomic_annotation.tsv"))
  utils::read.table("results/v0.5/genomic_annotation.tsv", header = TRUE, sep = "\t") else NULL
mot <- if (file.exists("results/v0.5/motif_validation.tsv"))
  utils::read.table("results/v0.5/motif_validation.tsv", header = TRUE, sep = "\t") else NULL
sens <- if (file.exists("results/v0.5/sensitivity_summary.tsv"))
  utils::read.table("results/v0.5/sensitivity_summary.tsv", header = TRUE, sep = "\t") else NULL

pdf("figures/v0.5_method_overview.pdf", 6, 3)
plot.new(); title("Layer 1: same BigWig+limma, different regions")
text(0.5, 0.6, "ESD: PELT segments\nConv: IDR union peaks", cex = 1)
dev.off()

pdf("figures/v0.5_gain_loss_comparison.pdf", 6, 4)
m <- rbind(c(ov$GAIN_esd, ov$GAIN_conv), c(ov$LOSS_esd, ov$LOSS_conv))
barplot(m, beside = TRUE, names.arg = c("ESD", "conventional"),
        legend.text = c("GAIN", "LOSS"), main = "GSE280947 primary calls")
dev.off()

pdf("figures/v0.5_region_overlap.pdf", 6, 4)
barplot(c(ov$bp_shared, ov$bp_esd - ov$bp_shared, ov$bp_conv - ov$bp_shared) / 1e3,
        names.arg = c("shared", "ESD-only", "conv-only"),
        ylab = "kb", main = "DER bp overlap")
dev.off()

pdf("figures/v0.5_effect_concordance.pdf", 5, 4)
plot(1, cc <- 1, type = "n", axes = FALSE, xlab = "", ylab = "",
     main = sprintf("sign concordance=%.3f  Pearson=%.3f", 0.984, 0.864))
dev.off()
ccf <- utils::read.table("results/v0.5/effect_concordance.tsv", header = TRUE, sep = "\t")
pdf("figures/v0.5_effect_concordance.pdf", 5, 4)
barplot(c(ccf$effect_pearson, ccf$sign_concordance), ylim = c(0, 1),
        names.arg = c("Pearson", "sign"), main = "ESD vs conventional effects")
dev.off()

pdf("figures/v0.5_region_widths.pdf", 6, 4)
boxplot(width ~ method, data = geom, outline = FALSE, ylab = "width (bp)",
        main = "Called region widths")
dev.off()

pdf("figures/v0.5_fragmentation.pdf", 5, 4)
barplot(c(ov$fragmentation_index, ov$mean_peaks_per_esd),
        names.arg = c("frag_index", "peaks/ESD"),
        main = "Fragmentation (predeclared)")
dev.off()

if (!is.null(mot) && "frac" %in% names(mot) && !all(is.na(mot$frac))) {
  pdf("figures/v0.5_motif_support.pdf", 5, 4)
  barplot(mot$frac, names.arg = mot$class, ylim = c(0, 1),
          main = "JUN/AP-1 motif fraction")
  dev.off()
} else {
  pdf("figures/v0.5_motif_support.pdf", 5, 3)
  plot.new(); title("Motif support: see motif_validation.tsv")
  dev.off()
}

if (!is.null(ann) && nrow(ann)) {
  pdf("figures/v0.5_annotation_comparison.pdf", 6, 4)
  tab <- xtabs(Freq ~ annotation + class, data = ann)
  barplot(tab, beside = TRUE, legend.text = TRUE, main = "Annotation by class")
  dev.off()
} else {
  pdf("figures/v0.5_annotation_comparison.pdf", 5, 3)
  plot.new(); title("annotation table missing")
  dev.off()
}

pdf("figures/v0.5_pathway_comparison.pdf", 5, 3)
plot.new(); title("GO BP: no BH-significant terms (q~0.96)")
dev.off()

pdf("figures/v0.5_negative_control.pdf", 5, 4)
barplot(rbind(null$GAIN, null$LOSS), beside = TRUE, names.arg = null$partition,
        legend.text = c("GAIN", "LOSS"), main = "True (0) vs scrambled labels")
dev.off()

if (!is.null(syn)) {
  pdf("figures/v0.5_synthetic_boundary_accuracy.pdf", 7, 4)
  barplot(rbind(syn$precision, syn$recall), beside = TRUE, names.arg = syn$case,
          las = 2, cex.names = 0.6, legend.text = c("precision", "recall"),
          main = "Synthetic bp accuracy")
  dev.off()
} else {
  pdf("figures/v0.5_synthetic_boundary_accuracy.pdf", 5, 3)
  plot.new(); title("synthetic pending")
  dev.off()
}

pdf("figures/v0.5_ablation.pdf", 6, 4)
barplot(rbind(abl$GAIN, abl$LOSS), beside = TRUE, names.arg = abl$ablation,
        las = 2, cex.names = 0.7, legend.text = c("GAIN", "LOSS"),
        main = "Ablation GAIN/LOSS")
dev.off()

message("figures written")
