#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
  else library(EpiSegDiff)
})

spec <- yaml::read_yaml("config/validation/topology_validation.yaml")
n <- as.integer(spec$n_per_class)
seed <- as.integer(spec$seed)
classes <- unlist(spec$classes)
geom <- c("EXPANSION", "CONTRACTION", "LEFT_SHIFT", "RIGHT_SHIFT",
          "INTERNAL_GAIN", "INTERNAL_LOSS", "COMPLEX")
matchc <- c("GAIN", "LOSS", "SPLIT", "MERGE")

parts <- list()
for (i in seq_along(classes)) {
  cl <- classes[i]
  if (cl %in% geom) {
    raw <- simulate_geometry_examples(cl, n = n, seed = seed + i)
    parts[[i]] <- data.frame(
      example_id = raw$example_id, truth = raw$truth, pred = raw$pred,
      rule = raw$rule, dL = raw$dL, dR = raw$dR, width_ratio = raw$width_ratio,
      borderline = raw$borderline, class_family = "geometry",
      stringsAsFactors = FALSE
    )
  } else if (cl %in% matchc) {
    raw <- simulate_match_class_examples(cl, n = n, seed = seed + 100 + i)
    parts[[i]] <- data.frame(
      example_id = raw$example_id, truth = raw$truth, pred = raw$pred,
      rule = raw$rule, dL = NA_real_, dR = NA_real_, width_ratio = NA_real_,
      borderline = NA, class_family = "matching",
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unknown class ", cl)
  }
}
tab <- do.call(rbind, parts)
dir.create("benchmark/simulation", recursive = TRUE, showWarnings = FALSE)
utils::write.table(tab, spec$output, sep = "\t", quote = FALSE, row.names = FALSE)
met <- event_class_metrics(tab$truth, tab$pred, classes = classes)
message("Wrote ", spec$output, " n=", nrow(tab),
        " macroF1=", round(met$macro_f1, 3),
        " COMPLEX_rate=", mean(tab$pred == "COMPLEX"))
saveRDS(met, "benchmark/simulation/topology_validation_metrics.rds")
