#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(20269)
dir.create("benchmark/simulation", showWarnings = FALSE, recursive = TRUE)
rows <- list(); k <- 0L
pipe <- default_config()
pipe$normalization$method <- "TMM"
pipe$segmentation$penalties <- c(8, 16)
pipe$segmentation$selected_penalty <- 8

for (w_bp in c(10000L, 20000L)) {
  for (bs in c(100L, 200L)) {
    wbin <- as.integer(w_bp / bs)
    n_bins <- wbin + 40L
    for (s in 1:8) {
      L_idx <- 15L
      R_idx <- L_idx + wbin - 1L
      ev <- data.frame(event_class = "GAIN", start_bin = L_idx, end_bin = R_idx,
                       effect = 1.6)
      sim <- simulate_count_assay(
        n_bins = n_bins, bin_size = bs, n_rep = 3, depth = 25,
        dispersion = 0.08, background_sd = 0.05, events = ev, seed = 5000 + k
      )
      true_L <- sim$bins$start[L_idx]
      true_R <- sim$bins$end[R_idx]
      # Oracle candidate = exact truth
      region_exact <- data.frame(chrom = "chr1", start = true_L, end = true_R)
      # Displaced / truncated candidate (upstream discovery error)
      region_off <- data.frame(
        chrom = "chr1",
        start = sim$bins$start[max(1L, L_idx + 3L)],
        end = sim$bins$end[min(n_bins, R_idx - 3L)]
      )
      sf_true <- sim$size_factors
      sf_est <- normalize_counts(sim$counts, method = "TMM")$size_factors
      for (cand_name in c("oracle_exact", "truncated")) {
        reg <- if (cand_name == "oracle_exact") region_exact else region_off
        for (sf_name in c("true_sf", "tmm_sf")) {
          sf <- if (sf_name == "true_sf") sf_true else sf_est
          for (ver in c("legacy", "symmetric")) {
            # treatment mean series on all bins
            keep <- sim$samples$condition == "treatment"
            mu <- rowMeans(sweep(sim$counts[, keep, drop = FALSE], 2, sf[keep], "/"))
            got <- end_walk_series(mu, sim$bins$start, sim$bins$end,
                                   reg$start, reg$end,
                                   search_radius = if (ver == "legacy") 40L else 8L,
                                   version = ver)
            k <- k + 1L
            rows[[k]] <- data.frame(
              sim_id = k, width_bp = w_bp, bin_size = bs,
              candidate = cand_name, norm = sf_name, method = ver,
              abs_L = abs(got$L - true_L), abs_R = abs(got$R - true_R),
              signed_L = got$L - true_L, signed_R = got$R - true_R,
              hausdorff = max(abs(got$L - true_L), abs(got$R - true_R)),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
}
tab <- do.call(rbind, rows)
utils::write.table(tab, "benchmark/simulation/broad_boundary_oracle_tests.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
print(aggregate(cbind(abs_L, abs_R, hausdorff) ~ candidate + norm + method, tab, median))
