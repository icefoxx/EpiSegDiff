#' Summarize BAM / SAM fragment QC
#'
#' Collects mapped fragments, properly paired fragments, duplicate
#' status, chromosome coverage, fragment-size distribution, and library
#' size. Does not call peaks.
#'
#' @param bam_path Path to a BAM or SAM file.
#' @param paired Logical; if TRUE, count proper pairs as fragments.
#' @param mapq_min Minimum mapping quality for the counted read/pair.
#' @return A list with scalar QC fields and per-chromosome tables.
#' @export
summarize_bam_qc <- function(bam_path, paired = TRUE, mapq_min = 0L) {
  if (!file.exists(bam_path)) {
    stop("BAM/SAM does not exist: ", bam_path, call. = FALSE)
  }
  .need_pkg("Rsamtools", "summarize_bam_qc()")
  .need_pkg("GenomicAlignments", "summarize_bam_qc()")
  .need_pkg("GenomicRanges", "summarize_bam_qc()")

  bf <- Rsamtools::BamFile(bam_path)
  hdr <- tryCatch(Rsamtools::scanBamHeader(bf), error = function(e) NULL)
  chrom_lengths <- integer()
  if (!is.null(hdr)) {
    tg <- hdr$targets
    if (!is.null(tg)) chrom_lengths <- as.integer(tg)
  }

  flag <- Rsamtools::scanBamFlag(
    isUnmappedQuery = FALSE,
    isSecondaryAlignment = FALSE
  )
  what <- c("qname", "flag", "rname", "pos", "mapq", "isize", "strand")
  param <- Rsamtools::ScanBamParam(flag = flag, what = what)
  recs <- Rsamtools::scanBam(bf, param = param)[[1]]

  n_mapped <- length(recs$qname)
  if (n_mapped == 0L) {
    return(list(
      path = bam_path,
      mapped_records = 0L,
      mapped_fragments = 0L,
      properly_paired_fragments = 0L,
      unpaired_or_not_proper = 0L,
      duplicate_records = 0L,
      library_size = 0L,
      fragment_size_summary = .empty_fragment_size_summary(),
      chromosome_coverage = data.frame(
        chrom = character(), n_records = integer(), n_fragments = integer(),
        stringsAsFactors = FALSE
      ),
      chrom_lengths = chrom_lengths,
      seqnames_in_file = character(),
      warnings = "BAM contains no mapped primary records"
    ))
  }

  flags <- recs$flag
  is_dup <- bitwAnd(flags, 1024L) != 0L
  is_paired <- bitwAnd(flags, 1L) != 0L
  is_proper <- bitwAnd(flags, 2L) != 0L
  is_first <- bitwAnd(flags, 64L) != 0L
  mapq <- recs$mapq
  mapq[is.na(mapq)] <- 0L
  pass_mapq <- mapq >= as.integer(mapq_min)

  chrom <- as.character(recs$rname)
  seqnames_in_file <- sort(unique(chrom))

  # Naming mismatch hint (chr1 vs 1)
  naming_warning <- NULL
  if (length(chrom_lengths)) {
    hdr_names <- names(chrom_lengths)
    if (length(intersect(seqnames_in_file, hdr_names)) == 0L &&
        length(seqnames_in_file) > 0L) {
      naming_warning <- sprintf(
        "No overlap between record seqnames (%s) and header targets (%s).",
        paste(utils::head(seqnames_in_file, 5), collapse = ","),
        paste(utils::head(hdr_names, 5), collapse = ",")
      )
    }
  }

  if (isTRUE(paired)) {
    keep <- is_paired & is_proper & is_first & pass_mapq
    frag_len <- abs(recs$isize[keep])
    frag_chrom <- chrom[keep]
    mapped_fragments <- sum(keep, na.rm = TRUE)
    properly_paired_fragments <- mapped_fragments
    unpaired_or_not_proper <- sum(is_paired & !is_proper & is_first, na.rm = TRUE) +
      sum(!is_paired, na.rm = TRUE)
  } else {
    keep <- pass_mapq
    frag_len <- rep(NA_real_, sum(keep))
    frag_chrom <- chrom[keep]
    mapped_fragments <- sum(keep)
    properly_paired_fragments <- sum(is_paired & is_proper & is_first & pass_mapq)
    unpaired_or_not_proper <- sum(!is_paired | !is_proper)
  }

  chr_tab <- as.data.frame(table(chrom), stringsAsFactors = FALSE)
  names(chr_tab) <- c("chrom", "n_records")
  frag_tab <- as.data.frame(table(frag_chrom), stringsAsFactors = FALSE)
  names(frag_tab) <- c("chrom", "n_fragments")
  chromosome_coverage <- merge(chr_tab, frag_tab, by = "chrom", all = TRUE)
  chromosome_coverage$n_records[is.na(chromosome_coverage$n_records)] <- 0L
  chromosome_coverage$n_fragments[is.na(chromosome_coverage$n_fragments)] <- 0L

  warns <- character()
  if (!is.null(naming_warning)) warns <- c(warns, naming_warning)
  if (any(grepl("^(chr)?[Mm][Tt]?$|^chrM$|^MT$", seqnames_in_file))) {
    warns <- c(warns, "Mitochondrial contig present; not dropped unless excluded by user.")
  }
  if (any(grepl("^(chr)?[XY]$", seqnames_in_file))) {
    warns <- c(warns, "Sex chromosomes present; not dropped unless excluded by user.")
  }

  list(
    path = bam_path,
    mapped_records = as.integer(n_mapped),
    mapped_fragments = as.integer(mapped_fragments),
    properly_paired_fragments = as.integer(properly_paired_fragments),
    unpaired_or_not_proper = as.integer(unpaired_or_not_proper),
    duplicate_records = as.integer(sum(is_dup)),
    library_size = as.integer(mapped_fragments),
    fragment_size_summary = .summarize_fragment_sizes(frag_len),
    chromosome_coverage = chromosome_coverage,
    chrom_lengths = chrom_lengths,
    seqnames_in_file = seqnames_in_file,
    warnings = if (length(warns)) paste(warns, collapse = " ") else NA_character_
  )
}

.empty_fragment_size_summary <- function() {
  data.frame(
    n = 0L, min = NA_real_, q25 = NA_real_, median = NA_real_,
    q75 = NA_real_, max = NA_real_, mean = NA_real_,
    stringsAsFactors = FALSE
  )
}

.summarize_fragment_sizes <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) return(.empty_fragment_size_summary())
  qs <- stats::quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
  data.frame(
    n = length(x), min = qs[1], q25 = qs[2], median = qs[3],
    q75 = qs[4], max = qs[5], mean = mean(x),
    stringsAsFactors = FALSE
  )
}

#' Extract fragments from a BAM as a 0-based data.frame
#'
#' @param bam_path BAM/SAM path.
#' @param paired Proper-pair fragments if TRUE; otherwise single-end reads.
#' @param mapq_min Minimum MAPQ.
#' @param counting_mode Unused for extraction; retained for a stable API.
#' @return data.frame with `chrom`, `start`, `end` (0-based half-open).
#' @export
extract_fragments <- function(bam_path, paired = TRUE, mapq_min = 0L,
                              counting_mode = "fragment_midpoint") {
  .need_pkg("Rsamtools", "extract_fragments()")
  .need_pkg("GenomicAlignments", "extract_fragments()")
  .need_pkg("GenomicRanges", "extract_fragments()")

  if (isTRUE(paired)) {
    flag <- Rsamtools::scanBamFlag(
      isPaired = TRUE,
      isProperPair = TRUE,
      isUnmappedQuery = FALSE,
      hasUnmappedMate = FALSE,
      isSecondaryAlignment = FALSE,
      isFirstMateRead = TRUE
    )
    param <- Rsamtools::ScanBamParam(
      flag = flag,
      what = c("rname", "pos", "qwidth", "isize", "mapq", "strand"),
      mapqFilter = as.integer(mapq_min)
    )
    recs <- Rsamtools::scanBam(bam_path, param = param)[[1]]
    if (!length(recs$pos)) {
      return(data.frame(chrom = character(), start = integer(), end = integer(),
                        stringsAsFactors = FALSE))
    }
    pos <- as.integer(recs$pos)             # 1-based BAM
    isize <- as.integer(recs$isize)
    qwidth <- as.integer(recs$qwidth)
    # Proper-pair fragment on the reference: start = min(pos, mate), width = |isize|
    # When isize is 0 or NA, fall back to the sequenced read.
    use_is <- is.finite(isize) & isize != 0L
    frag_start1 <- pos
    frag_end1 <- pos + qwidth - 1L
    frag_start1[use_is] <- ifelse(isize[use_is] > 0,
                                  pos[use_is],
                                  pos[use_is] + isize[use_is])
    frag_end1[use_is] <- frag_start1[use_is] + abs(isize[use_is]) - 1L
    start0 <- pmax(0L, frag_start1 - 1L)
    end0 <- pmax(start0, frag_end1)
    data.frame(
      chrom = as.character(recs$rname),
      start = as.integer(start0),
      end = as.integer(end0),
      stringsAsFactors = FALSE
    )
  } else {
    flag <- Rsamtools::scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE
    )
    param <- Rsamtools::ScanBamParam(
      flag = flag,
      what = c("rname", "pos", "qwidth", "mapq"),
      mapqFilter = as.integer(mapq_min)
    )
    recs <- Rsamtools::scanBam(bam_path, param = param)[[1]]
    if (!length(recs$pos)) {
      return(data.frame(chrom = character(), start = integer(), end = integer(),
                        stringsAsFactors = FALSE))
    }
    start0 <- pmax(0L, as.integer(recs$pos) - 1L)
    end0 <- start0 + as.integer(recs$qwidth)
    data.frame(
      chrom = as.character(recs$rname),
      start = as.integer(start0),
      end = as.integer(end0),
      stringsAsFactors = FALSE
    )
  }
}
