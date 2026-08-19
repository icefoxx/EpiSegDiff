#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
cfg_path <- NULL
i <- 1L
while (i <= length(args)) {
  if (args[i] %in% c("--config", "-c") && i < length(args)) {
    cfg_path <- args[i + 1L]; i <- i + 2L
  } else {
    i <- i + 1L
  }
}
if (is.null(cfg_path)) stop("Usage: Rscript scripts/run_episegdiff.R --config <yaml>")
if (file.exists("DESCRIPTION")) {
  suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(EpiSegDiff))
}
raw <- yaml::read_yaml(cfg_path)
if (!is.null(raw$assay$assay_family) && identical(raw$assay$assay_family, "continuous")) {
  res <- run_continuous_assay(cfg_path)
  message("continuous run GAIN=", res$manifest$n_GAIN,
          " LOSS=", res$manifest$n_LOSS)
} else {
  res <- run_episegdiff(cfg_path)
  message("count run complete")
}
invisible(res)
