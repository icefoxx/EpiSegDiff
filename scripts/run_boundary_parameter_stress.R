#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(20267)
dir.create("benchmark/simulation", showWarnings = FALSE, recursive = TRUE)
spec <- yaml::read_yaml("config/validation/boundary_audit.yaml")
rows <- list(); k <- 0L

push <- function(...) {
  k <<- k + 1L
  rows[[k]] <<- data.frame(sim_id = k, ..., stringsAsFactors = FALSE)
}

# Noiseless + light-noise series; not a full NB genome for every cell.
for (bs in as.integer(unlist(spec$stress$bin_sizes))) {
  for (width_bp in c(as.integer(unlist(spec$stress$narrow_bp)),
                     as.integer(unlist(spec$stress$broad_bp)))) {
    if (width_bp / bs < 3) next
    if (width_bp / bs > 80 && bs < 50) next
    n_bins <- as.integer(min(400L, max(40L, 6L * ceiling(width_bp / bs))))
    wbin <- as.integer(round(width_bp / bs))
    for (h in as.numeric(unlist(spec$stress$heights))) {
      for (rad in as.integer(unlist(spec$stress$radii))) {
        for (padx in as.numeric(unlist(spec$stress$paddings))) {
          for (ver in c("legacy", "symmetric", "likelihood")) {
            for (s in seq_len(as.integer(spec$stress$n_sim))) {
              L_idx <- 8L
              R_idx <- L_idx + wbin - 1L
              if (R_idx >= n_bins - 4L) next
              y <- numeric(n_bins)
              y[L_idx:R_idx] <- h
              if (h < 3.5) y <- y + stats::rnorm(n_bins, 0, 0.15)
              start <- as.integer((seq_len(n_bins) - 1L) * bs)
              end <- start + bs
              true_L <- start[L_idx]; true_R <- end[R_idx]
              pad_bins <- as.integer(round(padx * wbin))
              cand_L <- start[max(1L, L_idx - pad_bins)]
              cand_R <- end[min(n_bins, R_idx + pad_bins)]
              got <- if (ver == "likelihood") {
                likelihood_series(y, start, end, cand_L, cand_R, rad)
              } else {
                end_walk_series(y, start, end, cand_L, cand_R, rad, version = ver)
              }
              eL <- got$L - true_L
              eR <- got$R - true_R
              push(bin_size = bs, width_bp = width_bp, height = h,
                   search_radius = rad, padding_widths = padx,
                   method = ver, scale = if (width_bp < 5000) "narrow" else "broad",
                   true_L = true_L, true_R = true_R,
                   hat_L = got$L, hat_R = got$R,
                   signed_L = eL, signed_R = eR,
                   abs_L = abs(eL), abs_R = abs(eR),
                   A = abs(eL) - abs(eR),
                   A_norm = (abs(eL) - abs(eR)) / (abs(eL) + abs(eR) + 1e-6),
                   hausdorff = max(abs(eL), abs(eR)))
            }
          }
        }
      }
    }
  }
}
tab <- do.call(rbind, rows)
utils::write.table(tab, "benchmark/simulation/boundary_parameter_stress.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
message("stress n=", nrow(tab))
print(aggregate(cbind(abs_L, abs_R, A) ~ method + scale, tab, median))
