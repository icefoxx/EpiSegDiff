# Shared test helpers

skip_if_no_edger <- function() {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    skip("edgeR not installed")
  }
}

skip_if_no_rsamtools <- function() {
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    skip("Rsamtools not installed")
  }
}

# Minimal paired-end SAM (coordinate-sorted) for chr1 of length 1000.
# Two proper pairs:
#   pair A: fragment [100, 250)  mid = 175
#   pair B: fragment [400, 500)  mid = 450
write_tiny_paired_sam <- function(path, extra_unmapped = FALSE, extra_unpaired = FALSE) {
  lines <- c(
    "@HD\tVN:1.0\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:1000",
    "A\t99\tchr1\t101\t60\t50M\t=\t201\t150\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    "A\t147\tchr1\t201\t60\t50M\t=\t101\t-150\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    "B\t99\tchr1\t401\t60\t50M\t=\t451\t100\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    "B\t147\tchr1\t451\t60\t50M\t=\t401\t-100\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
  )
  if (extra_unmapped) {
    lines <- c(lines, "U\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII")
  }
  if (extra_unpaired) {
    lines <- c(lines,
               "S\t0\tchr1\t10\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII")
  }
  writeLines(lines, path)
  path
}

write_empty_mapped_sam <- function(path) {
  lines <- c(
    "@HD\tVN:1.0\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:1000"
  )
  writeLines(lines, path)
  path
}

sam_to_bam <- function(sam, dest_prefix) {
  skip_if_no_rsamtools()
  Rsamtools::asBam(sam, dest_prefix, overwrite = TRUE)
}

tiny_chrom_sizes <- function() {
  c(chr1 = 1000L)
}
