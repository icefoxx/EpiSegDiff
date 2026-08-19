#!/usr/bin/env Rscript
dir.create("figures", showWarnings = FALSE)
pdfw <- function(f, w = 6, h = 4, expr) { pdf(f, w, h); expr; dev.off() }

pdfw("figures/v0.6_method_landscape.pdf", 7, 3, {
  plot.new(); title("ESD change-point vs peak vs window-aggregate vs csaw/THOR")
  text(0.5, 0.5, "csaw/THOR real BAM: NOT_AVAILABLE\nLayer-1: same BigWig+limma")
})

if (file.exists("results/v0.6/synthetic_region_metrics.tsv")) {
  syn <- read.table("results/v0.6/synthetic_region_metrics.tsv", header = TRUE, sep = "\t")
  agg <- aggregate(f1 ~ method + width, syn, mean, na.rm = TRUE)
  pdfw("figures/v0.6_synthetic_power.pdf", 7, 4, {
    plot(agg$width[agg$method == "esd"], agg$f1[agg$method == "esd"],
         log = "x", type = "b", ylim = c(0, 1), xlab = "width", ylab = "mean F1",
         main = "Synthetic F1 by width")
    lines(agg$width[agg$method == "window"], agg$f1[agg$method == "window"],
          type = "b", col = 2)
    legend("bottomright", c("ESD", "window"), col = 1:2, lty = 1)
  })
  if (file.exists("results/v0.6/synthetic_boundary_metrics.tsv")) {
    b <- read.table("results/v0.6/synthetic_boundary_metrics.tsv", header = TRUE, sep = "\t")
    pdfw("figures/v0.6_boundary_error_by_width.pdf", 7, 4, {
      plot(b$width[b$method == "esd"], b$med_bound_err[b$method == "esd"],
           log = "x", xlab = "width", ylab = "median bound err",
           main = "Boundary error")
      points(b$width[b$method == "window"], b$med_bound_err[b$method == "window"], col = 2)
    })
  }
  if (file.exists("results/v0.6/synthetic_fragmentation.tsv")) {
    f <- read.table("results/v0.6/synthetic_fragmentation.tsv", header = TRUE, sep = "\t")
    pdfw("figures/v0.6_fragmentation_by_width.pdf", 7, 4, {
      plot(f$width[f$method == "esd"], f$one_call_frac[f$method == "esd"],
           log = "x", ylim = c(0, 1), xlab = "width", ylab = "one-call fraction")
      points(f$width[f$method == "window"], f$one_call_frac[f$method == "window"], col = 2)
    })
  }
} else {
  for (fn in c("figures/v0.6_synthetic_power.pdf",
               "figures/v0.6_boundary_error_by_width.pdf",
               "figures/v0.6_fragmentation_by_width.pdf")) {
    pdfw(fn, 5, 3, { plot.new(); title("synthetic pending") })
  }
}

if (file.exists("results/v0.6/synthetic_merge_split_errors.tsv")) {
  ms <- read.table("results/v0.6/synthetic_merge_split_errors.tsv", header = TRUE, sep = "\t")
  pdfw("figures/v0.6_merge_split_error.pdf", 6, 4, {
    boxplot(n_merged ~ method + gap_bp, data = ms, las = 2, cex.axis = 0.6,
            main = "Merge errors")
  })
} else pdfw("figures/v0.6_merge_split_error.pdf", 5, 3, { plot.new(); title("pending") })

pdfw("figures/v0.6_effect_noise_phase_diagram.pdf", 5, 3, {
  plot.new(); title("See synthetic_region_metrics.tsv")
})

if (file.exists("results/v0.6/null_simulation.tsv")) {
  nl <- read.table("results/v0.6/null_simulation.tsv", header = TRUE, sep = "\t")
  pdfw("figures/v0.6_null_pvalue_calibration.pdf", 6, 4, {
    boxplot(n_der ~ method, data = nl, main = "Null DER counts")
  })
} else pdfw("figures/v0.6_null_pvalue_calibration.pdf", 5, 3, { plot.new(); title("pending") })

if (file.exists("results/v0.6/fdr_calibration.tsv")) {
  fd <- read.table("results/v0.6/fdr_calibration.tsv", header = TRUE, sep = "\t")
  pdfw("figures/v0.6_nominal_vs_empirical_fdr.pdf", 6, 4, {
    plot(fd$alpha, fd$any_der, xlab = "nominal FDR", ylab = "P(any DER)",
         ylim = c(0, 1), main = "Null P(any call)")
    abline(0, 1, lty = 2)
  })
} else pdfw("figures/v0.6_nominal_vs_empirical_fdr.pdf", 5, 3, { plot.new(); title("pending") })

rd <- read.table("results/v0.6/realdata_calls_summary.tsv", header = TRUE, sep = "\t")
pdfw("figures/v0.6_realdata_gain_loss.pdf", 7, 4, {
  barplot(rbind(rd$GAIN, rd$LOSS), beside = TRUE, names.arg = paste(rd$dataset, rd$method),
          las = 2, cex.names = 0.5, legend.text = c("GAIN", "LOSS"))
})
pdfw("figures/v0.6_realdata_domain_fragmentation.pdf", 5, 4, {
  fr <- read.table("results/v0.6/realdata_fragmentation.tsv", header = TRUE, sep = "\t")
  barplot(fr$fragmentation_index, names.arg = fr$dataset, main = "peaks per ESD DER")
})

cs <- read.table("results/v0.6/gse280947_classC_support.tsv", header = TRUE, sep = "\t")
pdfw("figures/v0.6_classC_decomposition.pdf", 5, 4, {
  barplot(cs$n, names.arg = cs$subclass, main = "Class C subclasses")
})
pr <- read.table("results/v0.6/promoter_geometry_diagnostic.tsv", header = TRUE, sep = "\t")
pdfw("figures/v0.6_promoter_geometry_control.pdf", 7, 4, {
  barplot(rbind(pr$frac_c, pr$frac_a), beside = TRUE, names.arg = pr$comparison,
          las = 2, cex.names = 0.45, legend.text = c("C", "control"),
          main = "Promoter fraction")
})
pdfw("figures/v0.6_runtime_memory.pdf", 5, 3, {
  plot.new(); title("See runtime_benchmark.tsv")
})
message("v0.6 figures written")
