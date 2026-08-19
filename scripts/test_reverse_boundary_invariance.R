#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(20266)
dir.create("benchmark/simulation", showWarnings = FALSE, recursive = TRUE)
versions <- c("legacy", "symmetric")
rows <- list()
k <- 0L
for (ver in versions) {
  for (s in 1:40) {
    n <- 30L
    bs <- 50L
    w <- sample(4:10, 1)
    L_idx <- sample(4:(n - w - 3L), 1)
    R_idx <- L_idx + w - 1L
    h <- sample(c(-3, -2, 2, 3), 1)
    st <- make_step_series(n, L_idx, R_idx, height = h, bin_size = bs)
    # optional 1-bin candidate padding
    pad <- sample(0:2, 1)
    cand_L <- st$start[max(1L, L_idx - pad)]
    cand_R <- st$end[min(n, R_idx + pad)]
    fwd <- end_walk_series(st$y, st$start, st$end, cand_L, cand_R,
                           search_radius = 8L, version = ver)
    rev <- reverse_series(st$y, st$start, st$end)
    tr <- reverse_interval(st$true_L, st$true_R, rev$N)
    cr <- reverse_interval(cand_L, cand_R, rev$N)
    bwd <- end_walk_series(rev$y, rev$start, rev$end, cr$L, cr$R,
                           search_radius = 8L, version = ver)
    eLf <- fwd$L - st$true_L
    eRf <- fwd$R - st$true_R
    eLb <- bwd$L - tr$L
    eRb <- bwd$R - tr$R
    # expected: eLf ≈ -eRb and eRf ≈ -eLb (signed), abs: |eLf|≈|eRb|
    k <- k + 1L
    rows[[k]] <- data.frame(
      simulation_id = s, orientation = "forward", version = ver,
      true_left = st$true_L, true_right = st$true_R,
      estimated_left = fwd$L, estimated_right = fwd$R,
      left_error = abs(eLf), right_error = abs(eRf),
      mapped_reverse_left_error = abs(eRb),
      mapped_reverse_right_error = abs(eLb),
      invariance_difference = abs(abs(eLf) - abs(eRb)) + abs(abs(eRf) - abs(eLb)),
      signed_left = eLf, signed_right = eRf,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
    rows[[k]] <- data.frame(
      simulation_id = s, orientation = "reverse", version = ver,
      true_left = tr$L, true_right = tr$R,
      estimated_left = bwd$L, estimated_right = bwd$R,
      left_error = abs(eLb), right_error = abs(eRb),
      mapped_reverse_left_error = abs(eRf),
      mapped_reverse_right_error = abs(eLf),
      invariance_difference = abs(abs(eLb) - abs(eRf)) + abs(abs(eRb) - abs(eLf)),
      signed_left = eLb, signed_right = eRb,
      stringsAsFactors = FALSE
    )
  }
}
tab <- do.call(rbind, rows)
utils::write.table(tab, "benchmark/simulation/boundary_reverse_invariance.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
for (ver in versions) {
  d <- tab[tab$version == ver & tab$orientation == "forward", ]
  cat(sprintf("%s  mean|A|=%.1f  mean inv.diff=%.1f  median Lerr=%.1f  median Rerr=%.1f\n",
              ver, mean(d$left_error - d$right_error), mean(d$invariance_difference),
              median(d$left_error), median(d$right_error)))
}
