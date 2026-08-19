#' Recount raw fragments in candidate intervals
#'
#' Does not reuse the differential score as a count.
#'
#' @param regions data.frame with chrom, start, end (and optional region_id).
#' @param bam_files Named BAM paths or fragment data.frames.
#' @param sample_ids Sample names.
#' @param paired,mapq_min Passed to [extract_fragments()].
#' @param counting_mode `fragment_midpoint` or `fragment_overlap`.
#' @param discovery_scale Optional label stored on the output.
#' @return data.frame with region metadata and one count column per sample.
#' @export
recount_regions <- function(regions, bam_files, sample_ids = NULL,
                            paired = TRUE, mapq_min = 0L,
                            counting_mode = "fragment_overlap",
                            discovery_scale = NA) {
  regions <- .check_bed_frame(regions, "regions")
  if (nrow(regions) == 0L) {
    return(regions)
  }
  if (is.null(sample_ids)) sample_ids <- names(bam_files)
  if (is.null(sample_ids)) sample_ids <- paste0("sample_", seq_along(bam_files))
  if (!"region_id" %in% names(regions)) {
    regions$region_id <- seq_len(nrow(regions))
  }
  regions$discovery_scale <- discovery_scale

  for (j in seq_along(bam_files)) {
    item <- bam_files[[j]]
    if (is.data.frame(item)) {
      frags <- .check_bed_frame(item, sample_ids[j])
    } else {
      frags <- extract_fragments(item, paired = paired, mapq_min = mapq_min)
    }
    regions[[sample_ids[j]]] <- .count_frags_in_regions(frags, regions, counting_mode)
  }
  regions
}

.count_frags_in_regions <- function(frags, regions, counting_mode) {
  counts <- integer(nrow(regions))
  if (nrow(frags) == 0L) return(counts)
  if (counting_mode == "fragment_midpoint" || counting_mode == "insertion") {
    mid <- frags$start + (frags$end - frags$start) %/% 2L
    frags$start <- mid
    frags$end <- mid + 1L
  }
  for (chr in intersect(unique(frags$chrom), unique(regions$chrom))) {
    fi <- which(frags$chrom == chr)
    ri <- which(regions$chrom == chr)
    for (r in ri) {
      ov <- .interval_intersection(frags$start[fi], frags$end[fi],
                                   regions$start[r], regions$end[r]) > 0L
      counts[r] <- counts[r] + sum(ov)
    }
  }
  counts
}
