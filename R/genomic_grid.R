#' Build a common genomic tile grid
#'
#' Coordinates on the returned GRanges are 1-based inclusive (Bioconductor).
#' Downstream TSV/BED conversion uses 0-based half-open `[start, end)`.
#'
#' @param seqinfo A `Seqinfo` or named integer chrom sizes.
#' @param bin_size Bin width in bp.
#' @param chromosomes Optional restriction.
#' @export
build_genomic_grid <- function(seqinfo, bin_size = 200L, chromosomes = NULL) {
  bin_size <- as.integer(bin_size)
  if (is.na(bin_size) || bin_size < 1L) {
    stop("bin_size must be a positive integer.", call. = FALSE)
  }
  if (is.numeric(seqinfo) && !is.null(names(seqinfo))) {
    .need_pkg("GenomeInfoDb", "build_genomic_grid()")
    seqinfo <- GenomeInfoDb::Seqinfo(names(seqinfo), seqlengths = as.integer(seqinfo))
  }
  if (!is.null(chromosomes)) {
    keep <- intersect(as.character(chromosomes), as.character(GenomeInfoDb::seqnames(seqinfo)))
    if (!length(keep)) stop("No chromosomes remain after restriction.", call. = FALSE)
    seqinfo <- seqinfo[keep]
  }
  bins <- GenomicRanges::tileGenome(seqinfo, tilewidth = bin_size,
                                    cut.last.tile.in.chrom = TRUE)
  if (inherits(bins, "GRangesList")) bins <- unlist(bins, use.names = FALSE)
  bins
}

#' Convert a GRanges grid to a 0-based BED-style track frame
#' @export
granges_to_bed_frame <- function(gr) {
  data.frame(
    chrom = as.character(GenomicRanges::seqnames(gr)),
    start = as.integer(GenomicRanges::start(gr)) - 1L,
    end = as.integer(GenomicRanges::end(gr)),
    stringsAsFactors = FALSE
  )
}

#' Coarsen a binned matrix by an integer factor
#'
#' Linear-scale means only. Do not pass already log-transformed values.
#' @export
coarsen_signal_grid <- function(bins, signal, factor = 2L) {
  factor <- as.integer(factor)
  if (factor <= 1L) {
    return(list(bins = bins, signal = signal))
  }
  bed <- granges_to_bed_frame(bins)
  out <- coarsen_bed_matrix(bed, signal, factor = factor, from_bin = NA_integer_)
  new_bins <- GenomicRanges::GRanges(
    out$bed$chrom,
    IRanges::IRanges(out$bed$start + 1L, out$bed$end)
  )
  list(bins = new_bins, signal = out$signal)
}

#' Coarsen a BED-style linear signal matrix by consecutive bins
#'
#' Groups `factor` successive bins within each chromosome and takes the
#' unweighted mean of those bin summaries (already length-weighted inside
#' each fine bin). The last incomplete group on a chromosome is averaged
#' over the remaining bins.
#'
#' @param bed data.frame with `chrom`,`start`,`end` (0-based half-open).
#' @param signal bins × samples numeric matrix (linear scale).
#' @param factor Integer coarsening factor (2 → 200 bp from 100 bp).
#' @param from_bin Unused; retained for callers/manifests.
#' @export
coarsen_bed_matrix <- function(bed, signal, factor = 2L, from_bin = NA_integer_) {
  factor <- as.integer(factor)
  if (factor <= 1L) {
    return(list(bed = bed, signal = signal))
  }
  if (nrow(bed) != nrow(signal)) {
    stop("bed and signal have different numbers of rows.", call. = FALSE)
  }
  chroms <- unique(as.character(bed$chrom))
  beds <- vector("list", length(chroms))
  mats <- vector("list", length(chroms))
  nc <- ncol(signal)
  cn <- colnames(signal)
  for (k in seq_along(chroms)) {
    ii <- which(as.character(bed$chrom) == chroms[k])
    n <- length(ii)
    n_full <- n %/% factor
    rem <- n - n_full * factor
    ng <- n_full + as.integer(rem > 0L)
    new_start <- integer(ng)
    new_end <- integer(ng)
    new_m <- matrix(NA_real_, ng, nc, dimnames = list(NULL, cn))
    if (n_full) {
      src <- ii[seq_len(n_full * factor)]
      new_start[seq_len(n_full)] <- bed$start[ii[seq(1L, n_full * factor, by = factor)]]
      new_end[seq_len(n_full)] <- bed$end[ii[seq(factor, n_full * factor, by = factor)]]
      for (j in seq_len(nc)) {
        new_m[seq_len(n_full), j] <- .colMeans(signal[src, j], factor, n_full)
      }
    }
    if (rem) {
      last <- n_full + 1L
      src <- ii[(n_full * factor + 1L):n]
      new_start[last] <- bed$start[src[1]]
      new_end[last] <- bed$end[src[length(src)]]
      new_m[last, ] <- colMeans(signal[src, , drop = FALSE], na.rm = TRUE)
    }
    beds[[k]] <- data.frame(chrom = chroms[k], start = new_start, end = new_end,
                            stringsAsFactors = FALSE)
    mats[[k]] <- new_m
  }
  list(bed = do.call(rbind, beds), signal = do.call(rbind, mats))
}
