#' Extract Tn5 insertion positions from an ATAC BAM
#'
#' Offsets are configurable and are not applied to CUT&RUN/ChIP-seq.
#'
#' @param bam_path BAM path.
#' @param plus_offset Added to plus-strand 5' ends (default 4).
#' @param minus_offset Added to minus-strand 5' ends (default -5).
#' @param mapq_min Minimum MAPQ.
#' @return data.frame of 1-bp insertions (`chrom`, `start`, `end`) 0-based.
#' @export
extract_tn5_insertions <- function(bam_path, plus_offset = 4L, minus_offset = -5L,
                                   mapq_min = 0L) {
  .need_pkg("Rsamtools", "extract_tn5_insertions()")
  flag <- Rsamtools::scanBamFlag(
    isUnmappedQuery = FALSE,
    isSecondaryAlignment = FALSE
  )
  param <- Rsamtools::ScanBamParam(
    flag = flag,
    what = c("rname", "pos", "qwidth", "strand", "mapq"),
    mapqFilter = as.integer(mapq_min)
  )
  recs <- Rsamtools::scanBam(bam_path, param = param)[[1]]
  if (!length(recs$pos)) {
    return(data.frame(chrom = character(), start = integer(), end = integer(),
                      strand = character(), stringsAsFactors = FALSE))
  }
  strand <- as.character(recs$strand)
  pos1 <- as.integer(recs$pos)
  qwidth <- as.integer(recs$qwidth)
  # 5' end in 1-based BAM coordinates
  five1 <- ifelse(strand == "-", pos1 + qwidth - 1L, pos1)
  five1[strand != "-"] <- five1[strand != "-"] + as.integer(plus_offset)
  five1[strand == "-"] <- five1[strand == "-"] + as.integer(minus_offset)
  start0 <- pmax(0L, as.integer(five1) - 1L)
  data.frame(
    chrom = as.character(recs$rname),
    start = start0,
    end = start0 + 1L,
    strand = strand,
    stringsAsFactors = FALSE
  )
}

#' Build an ATAC insertion-count matrix
#'
#' Optional fragment channels are applied only when requested.
#'
#' @param bam_files Named BAM paths.
#' @param genome Chrom sizes.
#' @param bin_size Bin width.
#' @param plus_offset,minus_offset Tn5 offsets.
#' @param fragment_channels Named list of `c(min, max)` lengths, or empty.
#' @param ... Passed to [build_count_matrix()] for masks / chromosomes.
#' @export
build_atac_count_matrix <- function(bam_files, genome, bin_size,
                                    plus_offset = 4L, minus_offset = -5L,
                                    fragment_channels = list(),
                                    mapq_min = 0L, ...) {
  if (length(fragment_channels)) {
    warning("Fragment channels are optional ATAC summaries; the primary matrix is still all insertions.",
            call. = FALSE)
  }
  insertions <- lapply(bam_files, function(p) {
    if (is.data.frame(p)) return(p)
    extract_tn5_insertions(p, plus_offset = plus_offset,
                           minus_offset = minus_offset, mapq_min = mapq_min)
  })
  names(insertions) <- names(bam_files) %||% paste0("sample_", seq_along(bam_files))
  primary <- build_count_matrix(
    insertions, genome, bin_size,
    counting_mode = "insertion",
    sample_ids = names(insertions),
    paired = FALSE,
    ...
  )
  primary$assay <- "atac"
  primary$tn5_plus_offset <- as.integer(plus_offset)
  primary$tn5_minus_offset <- as.integer(minus_offset)
  primary$fragment_channels <- fragment_channels
  primary
}
