#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(20270)
dir.create("benchmark/simulation", showWarnings = FALSE, recursive = TRUE)
rows <- list()
for (s in 1:8) {
  ev <- example_topology_event(if (s %% 2L) "GAIN" else "LOSS", 50L)
  sim <- simulate_count_assay(n_bins = 160, n_rep = 3, depth = 24,
                              dispersion = 0.06, background_sd = 0.05,
                              events = ev, seed = 7000 + s)
  cfg <- default_config()
  cfg$samples <- sim$samples
  cfg$normalization$method <- "TMM"
  cfg$segmentation$penalties <- c(8, 16)
  cfg$segmentation$selected_penalty <- 8
  cfg$output$write_bed <- FALSE
  cfg$topology$enabled <- TRUE
  cfg$output$directory <- tempfile("l1on")
  a <- run_episegdiff(cfg, count_matrix = sim)
  cfg$topology$enabled <- FALSE
  cfg$output$directory <- tempfile("l1off")
  b <- run_episegdiff(cfg, count_matrix = sim)
  cols <- c("region_id", "logFC", "pvalue", "FDR")
  same <- isTRUE(all.equal(a$tested[, cols], b$tested[, cols], tolerance = 1e-10))
  same_iv <- isTRUE(all.equal(a$tested[, c("chrom", "start", "end")],
                              b$tested[, c("chrom", "start", "end")]))
  rows[[s]] <- data.frame(
    sim_id = s, event = ev$event_class[1],
    n_tested_on = nrow(a$tested), n_tested_off = nrow(b$tested),
    n_events_on = nrow(a$events), n_events_off = nrow(b$events),
    stats_identical = same, intervals_identical = same_iv,
    stringsAsFactors = FALSE
  )
}
tab <- do.call(rbind, rows)
utils::write.table(tab, "benchmark/simulation/level1_independence_test.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
print(tab)
message("all intervals identical: ", all(tab$intervals_identical),
        " all stats identical: ", all(tab$stats_identical))
