#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("benchmark/simulation", showWarnings = FALSE, recursive = TRUE)
bins <- build_bins(c(chr1 = 1000L), 100L)
rec <- reconstruct_interval(bins, data.frame(chrom = "chr1", start = 300L, end = 700L))
track <- data.frame(
  chrom = "chr1", start = bins$start, end = bins$end,
  score = ifelse(bins$start >= 300L & bins$end <= 700L, 2, 0), valid = TRUE
)
segs <- segment_track(track, penalty = 0.5, min_width = 1L, engine = "op")
hit <- segs[segs$effect_mean > 1, ]
cvr <- contiguous_valid_regions(bins)
out <- data.frame(
  check = c("bin_starts", "overlap_starts", "overlap_ends",
            "reconstruct_start", "reconstruct_end",
            "segment_start", "segment_end", "cvr_start", "cvr_end"),
  observed = c(paste(bins$start, collapse = ","),
               paste(bins$start[bins$start < 700 & bins$end > 300], collapse = ","),
               paste(bins$end[bins$start < 700 & bins$end > 300], collapse = ","),
               rec$start, rec$end,
               min(hit$start), max(hit$end),
               min(cvr$start), max(cvr$end)),
  expected = c("0,100,200,300,400,500,600,700,800,900",
               "300,400,500,600", "400,500,600,700",
               "300", "700", "300", "700", "0", "1000"),
  stringsAsFactors = FALSE
)
out$pass <- out$observed == out$expected
utils::write.table(out, "benchmark/simulation/boundary_coordinate_tests.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
message("coordinate tests pass=", all(out$pass))
print(out)
