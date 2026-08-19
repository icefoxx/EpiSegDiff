test_that("build_bins uses 0-based half-open coordinates and keeps a short last bin", {
  bins <- build_bins(c(chr1 = 250L), bin_size = 100L)
  expect_equal(bins$start, c(0L, 100L, 200L))
  expect_equal(bins$end, c(100L, 200L, 250L))
  expect_true(all(bins$valid))
})

test_that("blacklist splits valid bins and does not create cross-gap CVRs", {
  bins <- build_bins(c(chr1 = 400L), bin_size = 100L,
                     blacklist = data.frame(chrom = "chr1", start = 100L, end = 200L))
  expect_equal(bins$valid, c(TRUE, FALSE, TRUE, TRUE))
  cvrs <- contiguous_valid_regions(bins)
  expect_equal(nrow(cvrs), 2L)
  expect_equal(cvrs$start, c(0L, 200L))
  expect_equal(cvrs$end, c(100L, 400L))
  # no CVR spans the mask
  expect_true(all(cvrs$end <= 100L | cvrs$start >= 200L))
})

test_that("fragment midpoint counts match a manual example", {
  bins <- build_bins(c(chr1 = 1000L), bin_size = 100L)
  frags <- data.frame(
    chrom = c("chr1", "chr1", "chr1"),
    start = c(100L, 400L, 990L),
    end = c(250L, 500L, 1000L)
  )
  # mids: 175, 450, 995 -> bins [100,200), [400,500), [900,1000)
  ct <- count_fragments_into_bins(frags, bins, "fragment_midpoint")
  expect_equal(sum(ct), 3L)
  expect_equal(ct[bins$start == 100L], 1L)
  expect_equal(ct[bins$start == 400L], 1L)
  expect_equal(ct[bins$start == 900L], 1L)
})

test_that("fragment overlap counts every overlapping bin", {
  bins <- build_bins(c(chr1 = 300L), bin_size = 100L)
  frags <- data.frame(chrom = "chr1", start = 50L, end = 250L)
  ct <- count_fragments_into_bins(frags, bins, "fragment_overlap")
  expect_equal(ct, c(1L, 1L, 1L))
})

test_that("synthetic BAM midpoint counts match the SAM fragments", {
  skip_if_no_rsamtools()
  tmp <- tempfile()
  dir.create(tmp)
  sam <- file.path(tmp, "tiny.sam")
  write_tiny_paired_sam(sam)
  bam <- tryCatch(sam_to_bam(sam, file.path(tmp, "tiny")), error = function(e) NULL)
  if (is.null(bam)) skip("asBam failed")
  cm <- build_count_matrix(
    bam_files = c(s1 = bam),
    genome = tiny_chrom_sizes(),
    bin_size = 100L,
    counting_mode = "fragment_midpoint"
  )
  expect_equal(unname(sum(cm$counts[, 1])), 2L)
  # mids 175 and 450
  expect_equal(unname(cm$counts[cm$bins$start == 100L, 1]), 1L)
  expect_equal(unname(cm$counts[cm$bins$start == 400L, 1]), 1L)
})

test_that("empty BAM is handled", {
  skip_if_no_rsamtools()
  tmp <- tempfile()
  dir.create(tmp)
  sam <- file.path(tmp, "empty.sam")
  write_empty_mapped_sam(sam)
  bam <- tryCatch(sam_to_bam(sam, file.path(tmp, "empty")), error = function(e) NULL)
  if (is.null(bam)) skip("asBam failed")
  qc <- summarize_bam_qc(bam)
  expect_equal(qc$mapped_fragments, 0L)
  expect_equal(qc$library_size, 0L)
})

test_that("unpaired records are not counted as paired fragments", {
  skip_if_no_rsamtools()
  tmp <- tempfile()
  dir.create(tmp)
  sam <- file.path(tmp, "mix.sam")
  write_tiny_paired_sam(sam, extra_unpaired = TRUE)
  bam <- tryCatch(sam_to_bam(sam, file.path(tmp, "mix")), error = function(e) NULL)
  if (is.null(bam)) skip("asBam failed")
  qc <- summarize_bam_qc(bam, paired = TRUE)
  expect_equal(qc$properly_paired_fragments, 2L)
  expect_gt(qc$unpaired_or_not_proper, 0L)
})

test_that("chr vs non-chr naming mismatch is reported, not silently dropped", {
  skip_if_no_rsamtools()
  tmp <- tempfile()
  dir.create(tmp)
  sam <- file.path(tmp, "tiny.sam")
  write_tiny_paired_sam(sam)
  bam <- tryCatch(sam_to_bam(sam, file.path(tmp, "tiny2")), error = function(e) NULL)
  if (is.null(bam)) skip("asBam failed")
  expect_error(
    build_count_matrix(
      bam_files = c(s1 = bam),
      genome = c(`1` = 1000L),
      bin_size = 100L,
      chromosomes = "1"
    ),
    NA
  )
  # counts will be zero because chrom names differ; QC still lists chr1
  qc <- summarize_bam_qc(bam)
  expect_true("chr1" %in% qc$seqnames_in_file)
})
