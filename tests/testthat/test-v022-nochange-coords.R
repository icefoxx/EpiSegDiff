test_that("NO_CHANGE survives YAML read/write and is not YAML-null", {
  p <- tempfile(fileext = ".yaml")
  yaml::write_yaml(list(classes = c("GAIN", "NO_CHANGE", "LOSS")), p)
  raw <- paste(readLines(p), collapse = "\n")
  expect_match(raw, "NO_CHANGE")
  expect_false(grepl("(^|\\n)\\s*- null\\s*(\\n|$)", raw, ignore.case = TRUE))
  back <- yaml::read_yaml(p)
  expect_true("NO_CHANGE" %in% unlist(back$classes))
  expect_false(any(vapply(back$classes, is.null, logical(1))))
})

test_that("NO_CHANGE survives TSV export/import", {
  p <- tempfile(fileext = ".tsv")
  df <- data.frame(event_class = c("GAIN", "NO_CHANGE", "LOSS"),
                   start = c(0L, 10L, 20L), end = c(5L, 15L, 25L))
  utils::write.table(df, p, sep = "\t", quote = FALSE, row.names = FALSE)
  got <- utils::read.table(p, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expect_equal(got$event_class[2], "NO_CHANGE")
  expect_false(any(is.na(got$event_class)))
})

test_that("example_topology_event and simulator emit NO_CHANGE not YAML-null", {
  ev <- example_topology_event("NO_CHANGE", 10L)
  expect_equal(ev$event_class, "NO_CHANGE")
  # deprecated alias still accepted
  ev2 <- example_topology_event("NULL", 10L)
  expect_equal(ev2$event_class, "NO_CHANGE")
  sim <- simulate_null_counts(n_bins = 20, n_rep = 2, depth = 10, seed = 1)
  expect_true(all(sim$truth$event_class %in% c("NO_CHANGE")))
})

test_that("toy chromosome [300,700) reconstructs exactly", {
  bins <- build_bins(c(chr1 = 1000L), bin_size = 100L)
  expect_equal(bins$start, seq.int(0L, 900L, by = 100L))
  expect_equal(bins$end, seq.int(100L, 1000L, by = 100L))
  rec <- reconstruct_interval(bins, data.frame(chrom = "chr1", start = 300L, end = 700L))
  expect_equal(rec$start, 300L)
  expect_equal(rec$end, 700L)
  expect_equal(rec$n_bins, 4L)
  ov <- bins$start < 700L & bins$end > 300L
  expect_equal(bins$start[ov], c(300L, 400L, 500L, 600L))
  expect_equal(bins$end[ov], c(400L, 500L, 600L, 700L))
})

test_that("segment and BED export keep [300,700) exact", {
  track <- data.frame(
    chrom = "chr1",
    start = seq.int(0L, 900L, by = 100L),
    end = seq.int(100L, 1000L, by = 100L),
    score = c(rep(0, 3), rep(2, 4), rep(0, 3)),
    valid = TRUE
  )
  segs <- segment_track(track, penalty = 0.5, min_width = 1L, engine = "op")
  hit <- segs[segs$effect_mean > 1, ]
  expect_equal(min(hit$start), 300L)
  expect_equal(max(hit$end), 700L)
  tmp <- tempfile()
  dir.create(tmp)
  write_standard_outputs(tmp, segments = data.frame(
    chrom = "chr1", start = 300L, end = 700L, segment_id = 1L,
    penalty = 1, effect_mean = 2, candidate_source = "toy"
  ), write_bed = TRUE)
  bed <- utils::read.table(file.path(tmp, "differential_segments.bed"),
                           header = FALSE, sep = "\t")
  expect_equal(as.integer(bed$V2[1]), 300L)
  expect_equal(as.integer(bed$V3[1]), 700L)
})
