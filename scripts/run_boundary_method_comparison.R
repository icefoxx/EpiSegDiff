#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(20268)
dir.create("benchmark/simulation", showWarnings = FALSE, recursive = TRUE)
methods <- c("legacy", "symmetric", "likelihood")
rows <- list(); k <- 0L
t0_all <- proc.time()[["elapsed"]]
for (scale in c("narrow", "broad")) {
  n <- if (scale == "narrow") 40L else 80L
  w <- if (scale == "narrow") 8L else 25L
  bs <- if (scale == "narrow") 50L else 200L
  for (s in 1:30) {
    L_idx <- 10L
    R_idx <- L_idx + w - 1L
    y <- numeric(n)
    y[L_idx:R_idx] <- 2 + stats::rnorm(w, 0, 0.2)
    y <- y + stats::rnorm(n, 0, 0.12)
    start <- as.integer((seq_len(n) - 1L) * bs)
    end <- start + bs
    true_L <- start[L_idx]; true_R <- end[R_idx]
    cand_L <- start[max(1L, L_idx - 2L)]
    cand_R <- end[min(n, R_idx + 2L)]
    # reverse pair
    rev <- reverse_series(y, start, end)
    tr <- reverse_interval(true_L, true_R, rev$N)
    cr <- reverse_interval(cand_L, cand_R, rev$N)
    for (m in methods) {
      t0 <- proc.time()[["elapsed"]]
      got <- if (m == "likelihood") {
        likelihood_series(y, start, end, cand_L, cand_R, 6L)
      } else end_walk_series(y, start, end, cand_L, cand_R, 6L, m)
      bwd <- if (m == "likelihood") {
        likelihood_series(rev$y, rev$start, rev$end, cr$L, cr$R, 6L)
      } else end_walk_series(rev$y, rev$start, rev$end, cr$L, cr$R, 6L, m)
      rt <- proc.time()[["elapsed"]] - t0
      k <- k + 1L
      eL <- got$L - true_L; eR <- got$R - true_R
      rows[[k]] <- data.frame(
        sim_id = s, scale = scale, method = m, runtime_sec = rt,
        abs_L = abs(eL), abs_R = abs(eR),
        signed_L = eL, signed_R = eR,
        p90_pair = max(abs(eL), abs(eR)),
        hausdorff = max(abs(eL), abs(eR)),
        inv_diff = abs(abs(eL) - abs(bwd$R - tr$R)) + abs(abs(eR) - abs(bwd$L - tr$L)),
        stringsAsFactors = FALSE
      )
    }
  }
}
tab <- do.call(rbind, rows)
utils::write.table(tab, "benchmark/simulation/boundary_method_comparison.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
print(aggregate(cbind(abs_L, abs_R, hausdorff, inv_diff, runtime_sec) ~ method + scale,
                tab, median))
