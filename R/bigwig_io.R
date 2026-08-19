#' Validate a BigWig file header and basic signal stats
#'
#' @param path BigWig path.
#' @param genome Optional genome label (recorded only).
#' @param chromosomes Optional requested seqnames.
#' @param seqname_map Optional named character vector mapping file names
#'   to requested names.
#' @return A list of validation fields.
#' @export
validate_bigwig <- function(path, genome = NULL, chromosomes = NULL,
                            seqname_map = NULL) {
  if (!file.exists(path)) stop("BigWig does not exist: ", path, call. = FALSE)
  .need_pkg("rtracklayer", "validate_bigwig()")
  bw <- tryCatch(rtracklayer::BigWigFile(path), error = function(e) {
    stop("BigWig header is not readable: ", path, " (", e$message, ")", call. = FALSE)
  })
  si <- tryCatch(GenomeInfoDb::seqinfo(bw), error = function(e) {
    stop("Could not read BigWig seqinfo: ", path, call. = FALSE)
  })
  seqs <- as.character(GenomeInfoDb::seqnames(si))
  if (!length(seqs)) stop("BigWig has no chromosomes: ", path, call. = FALSE)
  style <- if (any(grepl("^chr", seqs, ignore.case = TRUE))) "UCSC" else "NCBI/Ensembl"
  if (!is.null(chromosomes)) {
    req <- as.character(chromosomes)
    if (!is.null(seqname_map)) {
      mapped <- unname(seqname_map[req])
      req_file <- ifelse(is.na(mapped), req, mapped)
    } else {
      req_file <- req
    }
    if (!any(req_file %in% seqs) && !any(req %in% seqs)) {
      stop("Requested chromosomes are incompatible with ", path,
           " (file has e.g. ", paste(utils::head(seqs, 4), collapse = ","), ").",
           call. = FALSE)
    }
  }
  # Light scan: first seq, up to 1e5 scores
  rng <- tryCatch({
    gr <- GenomicRanges::GRanges(
      seqs[1], IRanges::IRanges(1L, min(as.integer(GenomeInfoDb::seqlengths(si)[1]), 1e6L))
    )
    rtracklayer::import(bw, which = gr, as = "NumericList")[[1]]
  }, error = function(e) numeric())
  nfin <- sum(is.finite(rng))
  list(
    path = path,
    n_seq = length(seqs),
    seqnames = seqs,
    seqname_style = style,
    genome = genome %||% NA_character_,
    empty = length(rng) == 0L,
    finite_fraction = if (!length(rng)) NA_real_ else nfin / length(rng),
    min = if (nfin) min(rng[is.finite(rng)]) else NA_real_,
    max = if (nfin) max(rng[is.finite(rng)]) else NA_real_,
    q = if (nfin) stats::quantile(rng[is.finite(rng)], c(0, 0.25, 0.5, 0.75, 1),
                                  names = FALSE) else rep(NA_real_, 5)
  )
}

#' Summarize BigWigs onto a common genomic grid
#'
#' @param sample_sheet Output of [read_bigwig_sample_sheet()].
#' @param genome Genome label.
#' @param chromosomes Chromosomes to keep (NULL = canonical autosomes+X+Y if possible).
#' @param bin_size Bin width in bp.
#' @param summary_fun `mean` (length-weighted), `max`, or `area`
#'   (`sum(x * overlap_length)`). `median` and `sum` are rejected.
#' @param seqname_map Optional name map.
#' @return A list: `bins` (GRanges), `signal` (bins × samples), `samples`, `qc`.
#' @export
read_bigwig_assay <- function(sample_sheet, genome = NA_character_,
                              chromosomes = NULL, bin_size = 200L,
                              summary_fun = c("mean", "max", "area"),
                              seqname_map = NULL) {
  summary_fun <- match.arg(summary_fun)
  .need_pkg("rtracklayer", "read_bigwig_assay()")
  .need_pkg("GenomicRanges", "read_bigwig_assay()")
  .need_pkg("IRanges", "read_bigwig_assay()")
  sample_sheet <- read_bigwig_sample_sheet(sample_sheet)
  val <- lapply(sample_sheet$file, validate_bigwig, genome = genome,
                chromosomes = chromosomes, seqname_map = seqname_map)
  names(val) <- sample_sheet$sample

  si <- GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(sample_sheet$file[1]))
  std <- if (is.null(chromosomes)) {
    cand <- c(as.character(1:22), "X", "Y", paste0("chr", c(1:22, "X", "Y")))
    intersect(cand, as.character(GenomeInfoDb::seqnames(si)))
  } else {
    as.character(chromosomes)
  }
  if (!length(std)) {
    stop("No requested chromosomes found in the first BigWig.", call. = FALSE)
  }
  bins <- build_genomic_grid(si, bin_size = bin_size, chromosomes = std)
  mat <- matrix(0, length(bins), nrow(sample_sheet),
                dimnames = list(NULL, sample_sheet$sample))
  for (j in seq_len(nrow(sample_sheet))) {
    cov <- rtracklayer::import(sample_sheet$file[j], as = "RleList")
    keep <- intersect(as.character(GenomeInfoDb::seqlevels(bins)), names(cov))
    bins_j <- bins[as.character(GenomicRanges::seqnames(bins)) %in% keep]
    GenomeInfoDb::seqlevels(bins_j) <- keep
    vals <- summarize_rle_over_bins(bins_j, cov[keep], summary_fun)
    idx <- as.character(GenomicRanges::seqnames(bins)) %in% keep
    mat[idx, j] <- vals
    rm(cov, bins_j)
  }
  list(
    bins = bins,
    signal = mat,
    samples = sample_sheet,
    genome = genome,
    bin_size = as.integer(bin_size),
    summary_fun = summary_fun,
    validation = val,
    assay_family = "continuous",
    input_type = "bigwig"
  )
}

#' Length-weighted summary of an RleList over bins
#'
#' `mean` is coverage-length weighted. `area` is `sum(x_k * L_k)`.
#' `max` is the maximum value overlapping the bin. Uncovered bases are 0
#' in a typical whole-chromosome Rle import.
#'
#' @param bins GRanges.
#' @param cov Named RleList.
#' @param fun `mean`, `max`, or `area`.
#' @export
summarize_rle_over_bins <- function(bins, cov, fun = c("mean", "max", "area")) {
  if (length(fun) == 1L && !fun %in% c("mean", "max", "area")) {
    stop("summary_fun='", fun, "' is not implemented. Use mean, max, or area.",
         call. = FALSE)
  }
  fun <- match.arg(fun)
  out <- numeric(length(bins))
  chrs <- intersect(as.character(unique(GenomicRanges::seqnames(bins))), names(cov))
  for (chr in chrs) {
    ii <- which(as.character(GenomicRanges::seqnames(bins)) == chr)
    if (!length(ii)) next
    ir <- GenomicRanges::ranges(bins[ii])
    v <- IRanges::Views(cov[[chr]], ir)
    out[ii] <- switch(
      fun,
      mean = as.numeric(IRanges::viewMeans(v, na.rm = TRUE)),
      max = as.numeric(IRanges::viewMaxs(v, na.rm = TRUE)),
      area = as.numeric(IRanges::viewSums(v, na.rm = TRUE))
    )
  }
  out[!is.finite(out)] <- 0
  out
}
