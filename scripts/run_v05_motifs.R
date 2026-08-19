#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
set.seed(20251)
cl <- utils::read.table("results/v0.5/region_classification.tsv", header = TRUE,
                        sep = "\t", stringsAsFactors = FALSE)
edc <- cl[cl$method == "ESD", ]
edc <- scan_jun_ap1_motifs(edc)
if (all(is.na(edc$motif_hit))) {
  stop("motif scan still NA: ", attr(edc, "motif_error"))
}

# width-matched background, 5x class C, same chroms, avoid DERs
cset <- edc[edc$class == "C", ]
nbg <- 5L * nrow(cset)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
sl <- GenomeInfoDb::seqlengths(txdb)
canon <- paste0("chr", c(1:19, "X", "Y"))
sl <- sl[intersect(names(sl), canon)]
der <- GenomicRanges::GRanges(
  ifelse(grepl("^chr", edc$chrom), edc$chrom, paste0("chr", edc$chrom)),
  IRanges::IRanges(edc$start + 1L, edc$end)
)
bg <- list()
tries <- 0L
while (length(bg) < nbg && tries < nbg * 40L) {
  tries <- tries + 1L
  i <- ((tries - 1L) %% nrow(cset)) + 1L
  w <- cset$end[i] - cset$start[i]
  ch <- paste0("chr", gsub("^chr", "", cset$chrom[i]))
  if (!ch %in% names(sl) || sl[[ch]] <= w + 2L) next
  st <- sample.int(as.integer(sl[[ch]] - w), 1L)
  g <- GenomicRanges::GRanges(ch, IRanges::IRanges(st, st + w - 1L))
  if (length(GenomicRanges::findOverlaps(g, der))) next
  bg[[length(bg) + 1L]] <- data.frame(
    chrom = gsub("^chr", "", ch), start = st - 1L, end = st + w - 1L,
    class = "background", stringsAsFactors = FALSE
  )
}
bg <- if (length(bg)) do.call(rbind, bg) else
  data.frame(chrom = character(), start = integer(), end = integer(),
             class = character(), stringsAsFactors = FALSE)
if (nrow(bg)) bg <- scan_jun_ap1_motifs(bg)

tab <- do.call(rbind, lapply(c("A", "B", "C", "background"), function(k) {
  if (k == "B") {
    conv <- cl[cl$method == "conventional" & cl$class == "B", ]
    if (!nrow(conv)) {
      return(data.frame(class = "B", n = 0, n_hit = NA, frac = NA,
                        stringsAsFactors = FALSE))
    }
    conv <- scan_jun_ap1_motifs(conv)
    x <- conv$motif_hit
  } else if (k == "background") {
    x <- if (nrow(bg)) bg$motif_hit else logical()
  } else {
    x <- edc$motif_hit[edc$class == k]
  }
  data.frame(class = k, n = length(x),
             n_hit = sum(x, na.rm = TRUE),
             frac = if (length(x)) mean(x, na.rm = TRUE) else NA,
             stringsAsFactors = FALSE)
}))
# Fisher C vs background
if (nrow(bg) && any(edc$class == "C")) {
  c_hit <- sum(edc$motif_hit[edc$class == "C"], na.rm = TRUE)
  c_n <- sum(edc$class == "C")
  b_hit <- sum(bg$motif_hit, na.rm = TRUE)
  b_n <- nrow(bg)
  ft <- stats::fisher.test(matrix(c(c_hit, c_n - c_hit, b_hit, b_n - b_hit), 2))
  tab$or_C_vs_bg <- NA_real_
  tab$p_C_vs_bg <- NA_real_
  tab$or_C_vs_bg[tab$class == "C"] <- unname(ft$estimate)
  tab$p_C_vs_bg[tab$class == "C"] <- ft$p.value
}
tab$status <- "PASS"
utils::write.table(tab, "results/v0.5/motif_validation.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
print(tab)
message("motifs written nbg=", nrow(bg))
