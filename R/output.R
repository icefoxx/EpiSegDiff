#' Write standard EpiSegDiff output tables
#'
#' @param output_dir Directory.
#' @param segments Candidate / multi-scale segments.
#' @param tested Region-level tests.
#' @param events Boundary events.
#' @param write_bed Write BED files.
#' @return Invisibly, the paths written.
#' @export
write_standard_outputs <- function(output_dir, segments = NULL, tested = NULL,
                                   events = NULL, write_bed = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- character()
  if (!is.null(segments) && nrow(segments)) {
    tab <- .format_differential_segments(segments)
    p <- file.path(output_dir, "differential_segments.tsv")
    utils::write.table(tab, p, sep = "\t", quote = FALSE, row.names = FALSE)
    paths <- c(paths, p)
    if (isTRUE(write_bed)) {
      pb <- file.path(output_dir, "differential_segments.bed")
      .write_bed(tab, pb)
      paths <- c(paths, pb)
    }
  }
  if (!is.null(tested) && nrow(tested)) {
    p <- file.path(output_dir, "tested_regions.tsv")
    utils::write.table(tested, p, sep = "\t", quote = FALSE, row.names = FALSE)
    paths <- c(paths, p)
    if (isTRUE(write_bed) && all(c("chrom", "start", "end") %in% names(tested))) {
      .write_bed(tested, file.path(output_dir, "tested_regions.bed"))
    }
  }
  if (!is.null(events) && nrow(events)) {
    p <- file.path(output_dir, "boundary_events.tsv")
    utils::write.table(events, p, sep = "\t", quote = FALSE, row.names = FALSE)
    paths <- c(paths, p)
    if (isTRUE(write_bed) && all(c("chrom", "start", "end") %in% names(events))) {
      .write_bed(events, file.path(output_dir, "boundary_events.bed"),
                 name_col = "primary_event")
    }
  }
  invisible(paths)
}

.format_differential_segments <- function(segments) {
  if (!"segment_id" %in% names(segments)) segments$segment_id <- seq_len(nrow(segments))
  if (!"direction" %in% names(segments)) {
    segments$direction <- ifelse(!is.finite(segments$effect_mean), "NA",
                                 ifelse(segments$effect_mean > 0, "GAIN",
                                        ifelse(segments$effect_mean < 0, "LOSS", "STABLE")))
  }
  if (!"candidate_source" %in% names(segments)) {
    segments$candidate_source <- NA_character_
  }
  if (!"scale" %in% names(segments)) {
    segments$scale <- if ("penalty" %in% names(segments)) segments$penalty else NA
  }
  keep <- c("chrom", "start", "end", "segment_id", "scale", "effect_mean",
            "direction", "candidate_source")
  extra <- setdiff(names(segments), keep)
  segments[, c(keep, extra), drop = FALSE]
}

.write_bed <- function(df, path, name_col = NULL) {
  name <- if (!is.null(name_col) && name_col %in% names(df)) df[[name_col]]
  else if ("region_id" %in% names(df)) df$region_id
  else if ("segment_id" %in% names(df)) df$segment_id
  else seq_len(nrow(df))
  bed <- data.frame(
    chrom = df$chrom, start = df$start, end = df$end, name = name,
    stringsAsFactors = FALSE
  )
  utils::write.table(bed, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = FALSE)
}

.write_bigwig_track <- function(track, path) {
  .need_pkg("rtracklayer", "BigWig export")
  .need_pkg("GenomicRanges", "BigWig export")
  .need_pkg("IRanges", "BigWig export")
  ok <- track$valid & is.finite(track$score)
  if (!any(ok)) {
    warning("No finite valid scores to write to BigWig: ", path, call. = FALSE)
    return(invisible(NULL))
  }
  gr <- GenomicRanges::GRanges(
    seqnames = track$chrom[ok],
    ranges = IRanges::IRanges(start = track$start[ok] + 1L, end = track$end[ok]),
    score = track$score[ok]
  )
  rtracklayer::export(gr, path, format = "BigWig")
  invisible(path)
}
