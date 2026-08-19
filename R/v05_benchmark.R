#' Read a 3-column BED (0-based) with optional extra columns
#' @export
read_peak_bed <- function(path) {
  x <- utils::read.table(path, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE, comment.char = "#")
  names(x)[1:3] <- c("chrom", "start", "end")
  x$chrom <- as.character(x$chrom)
  x$start <- as.integer(x[[2]])
  x$end <- as.integer(x[[3]])
  x[, 1:3]
}

#' Union of peak BEDs (reduce overlaps)
#' @export
union_peak_regions <- function(...) {
  beds <- list(...)
  d <- do.call(rbind, beds)
  .need_pkg("GenomicRanges", "union_peak_regions()")
  g <- GenomicRanges::reduce(GenomicRanges::GRanges(
    d$chrom, IRanges::IRanges(d$start + 1L, d$end)
  ))
  data.frame(chrom = as.character(GenomicRanges::seqnames(g)),
             start = as.integer(GenomicRanges::start(g)) - 1L,
             end = as.integer(GenomicRanges::end(g)),
             stringsAsFactors = FALSE)
}

#' Mean of a binned linear matrix over arbitrary intervals
#' @export
summarize_matrix_over_intervals <- function(intervals, bed, signal) {
  .need_pkg("GenomicRanges", "summarize_matrix_over_intervals()")
  gi <- GenomicRanges::GRanges(intervals$chrom,
                               IRanges::IRanges(as.integer(intervals$start) + 1L,
                                                as.integer(intervals$end)))
  gb <- GenomicRanges::GRanges(bed$chrom,
                               IRanges::IRanges(as.integer(bed$start) + 1L,
                                                as.integer(bed$end)))
  ov <- GenomicRanges::findOverlaps(gi, gb)
  out <- matrix(NA_real_, nrow(intervals), ncol(signal),
                dimnames = list(NULL, colnames(signal)))
  if (!length(ov)) return(out)
  spl <- split(S4Vectors::subjectHits(ov), S4Vectors::queryHits(ov))
  for (nm in names(spl)) {
    i <- as.integer(nm)
    out[i, ] <- colMeans(signal[spl[[nm]], , drop = FALSE], na.rm = TRUE)
  }
  out
}

#' Limma on fixed regions given a linear signal matrix over those regions
#' @export
test_fixed_regions_limma <- function(intervals, region_signal, samples,
                                     numerator, denominator,
                                     transform_method = "log2_pseudocount",
                                     pseudocount = 0.01,
                                     fdr_threshold = 0.05) {
  xf <- transform_continuous(region_signal, transform_method, pseudocount)
  dummy <- data.frame(chrom = intervals$chrom,
                      start = intervals$start, end = intervals$end,
                      score = 0, valid = TRUE, stringsAsFactors = FALSE)
  segs <- dummy
  segs$region_id <- seq_len(nrow(segs))
  test_regions_limma(segs, xf, dummy, samples, numerator, denominator,
                     inference_mode = "replicated",
                     fdr_threshold = fdr_threshold)
}

.to_gr <- function(d) {
  if (!nrow(d)) return(GenomicRanges::GRanges())
  GenomicRanges::GRanges(d$chrom,
                         IRanges::IRanges(as.integer(d$start) + 1L,
                                          as.integer(d$end)))
}

.called <- function(d, fdr = 0.05) {
  if ("level1" %in% names(d)) {
    d[d$level1 %in% c("GAIN", "LOSS"), , drop = FALSE]
  } else {
    d[!is.na(d$FDR) & d$FDR <= fdr, , drop = FALSE]
  }
}

#' Classify ESD vs conventional DERs (predeclared A/B/C)
#' @export
classify_region_sets <- function(esd, conv, fdr = 0.05) {
  .need_pkg("GenomicRanges", "classify_region_sets()")
  e <- .called(esd, fdr)
  cnv <- .called(conv, fdr)
  if (!"effect" %in% names(e) && "log2FC" %in% names(e)) e$effect <- e$log2FC
  if (!"effect" %in% names(cnv) && "log2FC" %in% names(cnv)) cnv$effect <- cnv$log2FC
  ge <- .to_gr(e)
  gc <- .to_gr(cnv)
  class_e <- rep("C", nrow(e))
  class_c <- rep("B", nrow(cnv))
  opp_n <- 0L
  if (length(ge) && length(gc)) {
    ov <- GenomicRanges::findOverlaps(ge, gc)
    if (length(ov)) {
      qi <- S4Vectors::queryHits(ov)
      si <- S4Vectors::subjectHits(ov)
      we <- GenomicRanges::width(ge)[qi]
      wc <- GenomicRanges::width(gc)[si]
      inter <- GenomicRanges::width(GenomicRanges::pintersect(ge[qi], gc[si]))
      ro <- pmin(inter / we, inter / wc)
      same <- sign(e$effect[qi]) == sign(cnv$effect[si])
      opp_n <- sum(!same, na.rm = TRUE)
      hit <- same & ro >= 0.5
      class_e[unique(qi[hit])] <- "A"
      class_c[unique(si[hit])] <- "A"
    }
  }
  e$class <- class_e
  cnv$class <- class_c
  list(esd = e, conventional = cnv, n_opposite_overlaps = opp_n)
}

#' Overlap and geometry metrics between two called sets
#' @export
overlap_geometry_metrics <- function(esd, conv, peaks = NULL) {
  .need_pkg("GenomicRanges", "overlap_geometry_metrics()")
  e <- .called(esd)
  cnv <- .called(conv)
  ge <- .to_gr(e)
  gc <- .to_gr(cnv)
  bp <- function(g) sum(as.numeric(GenomicRanges::width(g)))
  inter <- if (length(ge) && length(gc)) GenomicRanges::intersect(ge, gc) else GenomicRanges::GRanges()
  jac <- if (!length(ge) && !length(gc)) 1 else if (!length(ge) || !length(gc)) 0 else
    bp(inter) / bp(GenomicRanges::reduce(c(ge, gc)))
  n_e_ov <- if (length(ge)) mean(GenomicRanges::countOverlaps(ge, gc) > 0) else NA
  n_c_ov <- if (length(gc)) mean(GenomicRanges::countOverlaps(gc, ge) > 0) else NA
  peaks_per <- NA_real_
  frag <- NA_real_
  if (!is.null(peaks) && nrow(peaks) && length(ge)) {
    gp <- .to_gr(peaks)
    co <- GenomicRanges::countOverlaps(ge, gp)
    peaks_per <- if (length(co)) mean(co) else NA
    frag <- sum(co) / max(length(ge), 1L)
  }
  data.frame(
    n_esd = nrow(e), n_conv = nrow(cnv),
    GAIN_esd = sum(e$level1 == "GAIN" | e$direction == "GAIN", na.rm = TRUE),
    LOSS_esd = sum(e$level1 == "LOSS" | e$direction == "LOSS", na.rm = TRUE),
    GAIN_conv = sum(cnv$level1 == "GAIN" | cnv$direction == "GAIN", na.rm = TRUE),
    LOSS_conv = sum(cnv$level1 == "LOSS" | cnv$direction == "LOSS", na.rm = TRUE),
    bp_esd = bp(ge), bp_conv = bp(gc), bp_shared = bp(inter),
    DER_bp_jaccard = jac,
    frac_esd_any_overlap = n_e_ov,
    frac_conv_any_overlap = n_c_ov,
    median_width_esd = if (nrow(e)) stats::median(e$end - e$start) else NA,
    median_width_conv = if (nrow(cnv)) stats::median(cnv$end - cnv$start) else NA,
    fragmentation_index = frag,
    mean_peaks_per_esd = peaks_per,
    stringsAsFactors = FALSE
  )
}

#' Annotate intervals against a TxDb hierarchy
#' @export
annotate_regions_txdb <- function(intervals, txdb = NULL) {
  if (is.null(txdb)) {
    if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)) {
      intervals$annotation <- NA_character_
      return(intervals)
    }
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  }
  d <- intervals
  d$chrom_ucsc <- ifelse(grepl("^chr", d$chrom), d$chrom, paste0("chr", d$chrom))
  gr <- GenomicRanges::GRanges(d$chrom_ucsc,
                               IRanges::IRanges(d$start + 1L, d$end))
  prom <- suppressWarnings(GenomicFeatures::promoters(txdb, upstream = 2000L, downstream = 2000L))
  ex <- GenomicFeatures::exons(txdb)
  gn <- GenomicFeatures::genes(txdb)
  is_prom <- GenomicRanges::countOverlaps(gr, prom) > 0
  is_ex <- GenomicRanges::countOverlaps(gr, ex) > 0
  is_gene <- GenomicRanges::countOverlaps(gr, gn) > 0
  d$annotation <- ifelse(is_prom, "promoter",
                         ifelse(is_ex, "exon",
                                ifelse(is_gene, "intron", "intergenic")))
  d
}

#' Nearest TSS linkage
#' @export
link_nearest_tss <- function(intervals, txdb = NULL) {
  if (is.null(txdb)) {
    if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)) {
      intervals$nearest_gene <- NA
      intervals$tss_distance <- NA
      return(intervals)
    }
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  }
  d <- intervals
  d$chrom_ucsc <- ifelse(grepl("^chr", d$chrom), d$chrom, paste0("chr", d$chrom))
  gr <- GenomicRanges::GRanges(d$chrom_ucsc,
                               IRanges::IRanges(d$start + 1L, d$end))
  tss <- GenomicRanges::resize(GenomicFeatures::genes(txdb), width = 1L, fix = "start")
  n <- GenomicRanges::nearest(gr, tss)
  d$nearest_gene <- NA_character_
  d$tss_distance <- NA_integer_
  ok <- !is.na(n)
  if (any(ok)) {
    d$nearest_gene[ok] <- as.character(names(tss)[n[ok]])
    if (!length(names(tss))) {
      d$nearest_gene[ok] <- as.character(GenomicRanges::mcols(tss)$gene_id[n[ok]])
    }
    d$tss_distance[ok] <- GenomicRanges::distance(gr[ok], tss[n[ok]])
  }
  d
}

#' JUN/AP-1 motif scan (predeclared JASPAR IDs)
#' @export
scan_jun_ap1_motifs <- function(intervals, genome = NULL, p = 1e-4) {
  if (!requireNamespace("motifmatchr", quietly = TRUE) ||
      !requireNamespace("JASPAR2020", quietly = TRUE) ||
      !requireNamespace("TFBSTools", quietly = TRUE)) {
    intervals$motif_hit <- NA
    return(intervals)
  }
  if (is.null(genome)) {
    if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
      intervals$motif_hit <- NA
      return(intervals)
    }
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  }
  d <- intervals
  d$chrom_ucsc <- ifelse(grepl("^chr", d$chrom), d$chrom, paste0("chr", d$chrom))
  gr <- GenomicRanges::GRanges(d$chrom_ucsc,
                               IRanges::IRanges(d$start + 1L, d$end))
  mats <- list(
    tryCatch(TFBSTools::getMatrixByID(JASPAR2020::JASPAR2020, "MA0488.1"),
             error = function(e) NULL),
    tryCatch(TFBSTools::getMatrixByID(JASPAR2020::JASPAR2020, "MA0099.3"),
             error = function(e) NULL)
  )
  mats <- Filter(Negate(is.null), mats)
  if (!length(mats)) {
    d$motif_hit <- NA
    attr(d, "motif_error") <- "JASPAR IDs MA0488.1/MA0099.3 not found"
    return(d)
  }
  pwm <- do.call(TFBSTools::PWMatrixList, lapply(mats, TFBSTools::toPWM))
  names(pwm) <- vapply(mats, function(x) x@ID, "")
  hits <- motifmatchr::matchMotifs(pwm, gr, genome = genome, p.cutoff = p)
  mm <- motifmatchr::motifMatches(hits)
  d$motif_hit <- as.logical(Matrix::rowSums(mm) > 0)
  d
}

#' Balanced 2+2 condition permutations excluding true and reverse
#' @export
balanced_label_partitions <- function(samples, numerator, denominator) {
  idx <- split(seq_len(nrow(samples)), samples$condition)
  a <- idx[[denominator]]
  b <- idx[[numerator]]
  all_i <- c(a, b)
  comb <- utils::combn(all_i, 2L, simplify = FALSE)
  true <- list(sort(a), sort(b))
  revr <- list(sort(b), sort(a))
  out <- list()
  for (g1 in comb) {
    g2 <- setdiff(all_i, g1)
    if (length(g2) != 2L) next
    key <- list(sort(g1), sort(g2))
    if (identical(key, true) || identical(key, revr)) next
    # canonical order
    dup <- any(vapply(out, function(z) identical(z, key) ||
                        identical(z, list(key[[2]], key[[1]])), logical(1)))
    if (dup) next
    out[[length(out) + 1L]] <- key
  }
  out
}
