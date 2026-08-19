#' Build genomic bins
#'
#' Coordinates are BED 0-based half-open `[start, end)`. The last bin of
#' each chromosome is kept even if shorter than `bin_size`.
#'
#' @param chrom_sizes Named integer vector of chromosome lengths, or a
#'   two-column table path (`chrom`, `length`).
#' @param bin_size Bin width in bp.
#' @param chromosomes Optional character vector to restrict chromosomes.
#' @param blacklist Optional BED data.frame or path.
#' @param extra_exclude Optional additional BED data.frame or path.
#' @return data.frame with `chrom`, `start`, `end`, `bin_id`, `valid`.
#' @export
build_bins <- function(chrom_sizes, bin_size, chromosomes = NULL,
                       blacklist = NULL, extra_exclude = NULL) {
  bin_size <- as.integer(bin_size)
  if (is.na(bin_size) || bin_size < 1L) {
    stop("bin_size must be a positive integer.", call. = FALSE)
  }
  sizes <- .read_chrom_sizes(chrom_sizes)
  if (!is.null(chromosomes)) {
    chromosomes <- as.character(chromosomes)
    missing <- setdiff(chromosomes, names(sizes))
    if (length(missing)) {
      stop("Requested chromosomes absent from chrom_sizes: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    sizes <- sizes[chromosomes]
  }
  rows <- vector("list", length(sizes))
  for (i in seq_along(sizes)) {
    chr <- names(sizes)[i]
    L <- as.integer(sizes[[i]])
    if (is.na(L) || L < 1L) {
      rows[[i]] <- data.frame(
        chrom = character(), start = integer(), end = integer(),
        stringsAsFactors = FALSE
      )
      next
    }
    starts <- seq.int(0L, L - 1L, by = bin_size)
    ends <- pmin(starts + bin_size, L)
    rows[[i]] <- data.frame(
      chrom = chr, start = as.integer(starts), end = as.integer(ends),
      stringsAsFactors = FALSE
    )
  }
  bins <- do.call(rbind, rows)
  if (is.null(bins) || nrow(bins) == 0L) {
    bins <- data.frame(
      chrom = character(), start = integer(), end = integer(),
      bin_id = integer(), valid = logical(), stringsAsFactors = FALSE
    )
    return(bins)
  }
  bins$bin_id <- seq_len(nrow(bins))
  bins$valid <- TRUE

  mask <- .bind_beds(blacklist, extra_exclude)
  if (!is.null(mask) && nrow(mask) > 0L) {
    bins$valid <- !.bin_overlaps_mask(bins, mask)
  }
  bins
}

.read_chrom_sizes <- function(chrom_sizes) {
  if (is.numeric(chrom_sizes) && !is.null(names(chrom_sizes))) {
    v <- as.integer(chrom_sizes)
    names(v) <- names(chrom_sizes)
    return(v)
  }
  if (is.data.frame(chrom_sizes)) {
    if (!all(c("chrom", "length") %in% names(chrom_sizes))) {
      stop("chrom_sizes data.frame needs columns chrom, length.", call. = FALSE)
    }
    v <- as.integer(chrom_sizes$length)
    names(v) <- as.character(chrom_sizes$chrom)
    return(v)
  }
  if (is.character(chrom_sizes) && length(chrom_sizes) == 1L) {
    if (!file.exists(chrom_sizes)) {
      stop("chrom_sizes file does not exist: ", chrom_sizes, call. = FALSE)
    }
    tab <- utils::read.table(chrom_sizes, header = FALSE, sep = "\t",
                             stringsAsFactors = FALSE)
    if (ncol(tab) < 2L) stop("chrom_sizes file must have chrom and length.", call. = FALSE)
    v <- as.integer(tab[[2]])
    names(v) <- as.character(tab[[1]])
    return(v)
  }
  stop("chrom_sizes must be a named integer vector, data.frame, or file path.",
       call. = FALSE)
}

.bind_beds <- function(...) {
  items <- list(...)
  frames <- list()
  for (x in items) {
    if (is.null(x)) next
    if (is.character(x) && length(x) == 1L) {
      if (!nzchar(x) || !file.exists(x)) {
        if (nzchar(x)) {
          stop("BED path does not exist: ", x, call. = FALSE)
        }
        next
      }
      tab <- utils::read.table(x, header = FALSE, sep = "\t",
                               stringsAsFactors = FALSE, comment.char = "#")
      x <- data.frame(chrom = tab[[1]], start = as.integer(tab[[2]]),
                      end = as.integer(tab[[3]]), stringsAsFactors = FALSE)
    }
    if (is.data.frame(x)) {
      frames[[length(frames) + 1L]] <- .check_bed_frame(x, "mask BED")
    }
  }
  if (!length(frames)) return(NULL)
  do.call(rbind, frames)
}

.bin_overlaps_mask <- function(bins, mask) {
  hit <- rep(FALSE, nrow(bins))
  for (chr in unique(bins$chrom)) {
    bi <- which(bins$chrom == chr)
    mi <- which(mask$chrom == chr)
    if (!length(mi)) next
    for (j in bi) {
      if (any(.interval_intersection(bins$start[j], bins$end[j],
                                     mask$start[mi], mask$end[mi]) > 0L)) {
        hit[j] <- TRUE
      }
    }
  }
  hit
}

#' Count fragments into bins
#'
#' @param fragments data.frame with `chrom`, `start`, `end` (0-based).
#' @param bins Output of [build_bins()].
#' @param counting_mode `fragment_midpoint` or `fragment_overlap`.
#' @return Integer vector of counts, one per bin (same order as `bins`).
#' @export
count_fragments_into_bins <- function(fragments, bins,
                                      counting_mode = c("fragment_midpoint",
                                                        "fragment_overlap",
                                                        "insertion")) {
  counting_mode <- match.arg(counting_mode)
  counts <- integer(nrow(bins))
  if (nrow(bins) == 0L || nrow(fragments) == 0L) return(counts)

  fragments <- .check_bed_frame(fragments, "fragments")
  if (counting_mode %in% c("fragment_midpoint", "insertion")) {
    mid <- fragments$start + (fragments$end - fragments$start) %/% 2L
    fragments$start <- mid
    fragments$end <- mid + 1L
  }

  for (chr in intersect(unique(bins$chrom), unique(fragments$chrom))) {
    bi <- which(bins$chrom == chr)
    fi <- which(fragments$chrom == chr)
    bstart <- bins$start[bi]
    bend <- bins$end[bi]
    # Assume bins on a chromosome are contiguous and sorted.
    bs <- bstart[1]
    width <- bend[1] - bstart[1]
    # Variable last-bin width: map using findInterval on starts.
    idx <- findInterval(fragments$start[fi], bstart, left.open = FALSE)
    in_range <- idx >= 1L & idx <= length(bi)
    # Confirm the fragment actually intersects the chosen bin.
    if (counting_mode == "fragment_overlap") {
      # Each fragment may overlap multiple bins.
      for (k in fi) {
        ov <- .interval_intersection(fragments$start[k], fragments$end[k],
                                     bstart, bend) > 0L
        if (any(ov)) counts[bi[ov]] <- counts[bi[ov]] + 1L
      }
    } else {
      ok <- in_range
      if (any(ok)) {
        pos <- fragments$start[fi][ok]
        ii <- idx[ok]
        inside <- pos >= bstart[ii] & pos < bend[ii]
        if (any(inside)) {
          tab <- table(ii[inside])
          counts[bi[as.integer(names(tab))]] <-
            counts[bi[as.integer(names(tab))]] + as.integer(tab)
        }
      }
    }
  }
  counts
}

#' Build a binned fragment-count matrix
#'
#' @param bam_files Named character vector of BAM paths, or a list of
#'   fragment data.frames (for tests / simulation).
#' @param genome Named chrom sizes (see [build_bins()]).
#' @param bin_size Bin width in bp.
#' @param blacklist Optional BED path or data.frame.
#' @param counting_mode `fragment_midpoint` or `fragment_overlap`.
#' @param sample_ids Sample names; default `names(bam_files)`.
#' @param chromosomes Optional chromosome restriction.
#' @param paired Passed to [extract_fragments()].
#' @param mapq_min Minimum MAPQ.
#' @param extra_exclude Optional extra mask BED.
#' @param mask_low_coverage If TRUE, bins with count-sum < `min_count_sum`
#'   are marked invalid (can create artificial CVR edges).
#' @param min_count_sum Threshold used only if `mask_low_coverage`.
#' @return A list with `bins`, `counts` (integer matrix), and `qc`.
#' @export
build_count_matrix <- function(bam_files, genome, bin_size, blacklist = NULL,
                               counting_mode = c("fragment_midpoint",
                                                 "fragment_overlap",
                                                 "insertion"),
                               sample_ids = NULL,
                               chromosomes = NULL,
                               paired = TRUE,
                               mapq_min = 0L,
                               extra_exclude = NULL,
                               mask_low_coverage = FALSE,
                               min_count_sum = 0L) {
  counting_mode <- match.arg(counting_mode)
  if (is.null(sample_ids)) {
    sample_ids <- names(bam_files)
  }
  if (is.null(sample_ids) || any(!nzchar(sample_ids))) {
    sample_ids <- paste0("sample_", seq_along(bam_files))
  }
  if (length(sample_ids) != length(bam_files)) {
    stop("sample_ids and bam_files have different lengths.", call. = FALSE)
  }

  bins <- build_bins(genome, bin_size, chromosomes = chromosomes,
                     blacklist = blacklist, extra_exclude = extra_exclude)
  mat <- matrix(0L, nrow = nrow(bins), ncol = length(sample_ids),
                dimnames = list(NULL, sample_ids))
  qc <- vector("list", length(sample_ids))
  names(qc) <- sample_ids

  for (j in seq_along(bam_files)) {
    item <- bam_files[[j]]
    if (is.data.frame(item)) {
      frags <- .check_bed_frame(item, sample_ids[j])
      qc[[j]] <- list(
        path = NA_character_,
        mapped_fragments = nrow(frags),
        library_size = nrow(frags),
        source = "fragment_data_frame"
      )
    } else {
      qc[[j]] <- summarize_bam_qc(item, paired = paired, mapq_min = mapq_min)
      frags <- extract_fragments(item, paired = paired, mapq_min = mapq_min,
                                 counting_mode = counting_mode)
    }
    mat[, j] <- count_fragments_into_bins(frags, bins, counting_mode)
  }

  if (isTRUE(mask_low_coverage)) {
    rs <- rowSums(mat)
    bins$valid <- bins$valid & (rs >= as.integer(min_count_sum))
  }

  list(
    bins = bins,
    counts = mat,
    qc = qc,
    counting_mode = counting_mode,
    bin_size = as.integer(bin_size),
    coordinate_system = "bed_0based"
  )
}

#' Split bins into contiguous valid regions
#'
#' @param bins Bin table with `chrom`, `start`, `end`, `valid`.
#' @return data.frame of CVRs: `cvr_id`, `chrom`, `start`, `end`,
#'   `bin_start_idx`, `bin_end_idx`, `n_bins`.
#' @export
contiguous_valid_regions <- function(bins) {
  if (nrow(bins) == 0L) {
    return(data.frame(
      cvr_id = integer(), chrom = character(), start = integer(),
      end = integer(), bin_start_idx = integer(), bin_end_idx = integer(),
      n_bins = integer(), stringsAsFactors = FALSE
    ))
  }
  valid <- as.logical(bins$valid)
  valid[is.na(valid)] <- FALSE
  run_break <- c(TRUE, bins$chrom[-1] != bins$chrom[-nrow(bins)] |
                   (valid[-1] != valid[-nrow(bins)]) |
                   (valid[-1] & (bins$start[-1] != bins$end[-nrow(bins)])))
  run_id <- cumsum(run_break)
  keep <- valid
  if (!any(keep)) {
    return(data.frame(
      cvr_id = integer(), chrom = character(), start = integer(),
      end = integer(), bin_start_idx = integer(), bin_end_idx = integer(),
      n_bins = integer(), stringsAsFactors = FALSE
    ))
  }
  idx <- split(seq_len(nrow(bins))[keep], run_id[keep])
  rows <- lapply(seq_along(idx), function(i) {
    ii <- idx[[i]]
    data.frame(
      cvr_id = i,
      chrom = bins$chrom[ii[1]],
      start = bins$start[ii[1]],
      end = bins$end[ii[length(ii)]],
      bin_start_idx = ii[1],
      bin_end_idx = ii[length(ii)],
      n_bins = length(ii),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
