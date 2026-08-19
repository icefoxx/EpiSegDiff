#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("cache", showWarnings = FALSE)
sheet <- read_bigwig_sample_sheet("configs/gse280947_samples.tsv")
ids <- vapply(sheet$file, file_stable_id, character(1))
man_path <- "cache/gse280947_100bp_signal_matrix_manifest.yaml"
rds_path <- "cache/gse280947_100bp_signal_matrix.rds"
if (file.exists(rds_path) && file.exists(man_path)) {
  old <- yaml::read_yaml(man_path)
  old_ids <- unname(unlist(old$checksums))
  same_ids <- identical(as.character(old_ids), as.character(unname(ids)))
  same_bin <- as.integer(old$bin_size) == 100L
  same_fun <- identical(old$summary_fun, "mean")
  same_genome <- identical(old$genome, "mm10")
  if (same_ids && same_bin && same_fun && same_genome) {
    message("Cache valid; skip import")
    quit(save = "no")
  }
  message("Cache stale; rebuilding")
}
t0 <- proc.time()[["elapsed"]]
message("Importing 4 BigWigs at 100 bp (linear mean)...")
assay <- read_bigwig_assay(sheet, genome = "mm10", bin_size = 100L, summary_fun = "mean")
bed <- granges_to_bed_frame(assay$bins)
cache <- list(bed = bed, signal = assay$signal, samples = sheet,
              bin_size = 100L, summary_fun = "mean", genome = "mm10",
              valid = apply(assay$signal, 1, function(z) any(is.finite(z))))
saveRDS(cache, rds_path)
yaml::write_yaml(list(
  path = rds_path,
  bin_size = 100L,
  summary_fun = "mean",
  genome = "mm10",
  n_bins = nrow(bed),
  samples = sheet$sample,
  checksums = as.list(ids),
  software = as.character(utils::packageVersion("EpiSegDiff")),
  built = as.character(Sys.time()),
  runtime_sec = proc.time()[["elapsed"]] - t0,
  coordinate_convention = "BED_0_based_half_open"
), man_path)
message("Wrote ", rds_path, " n_bins=", nrow(bed),
        " runtime_sec=", round(proc.time()[["elapsed"]] - t0, 1))
