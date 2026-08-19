#!/usr/bin/env Rscript
# Phase 1: profile current vs new engines on real GSE280947 chromosomes.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})
dir.create("results/v0.4.2", recursive = TRUE, showWarnings = FALSE)

rds <- "cache/gse280947_100bp_signal_matrix.rds"
if (!file.exists(rds)) stop("Need 100 bp cache")
cache <- readRDS(rds)
m200 <- coarsen_bed_matrix(cache$bed, cache$signal, factor = 2L, from_bin = 100L)
sheet <- cache$samples
xf <- transform_continuous(m200$signal, "log2_pseudocount", 0.01)
eff <- estimate_continuous_effect(xf, sheet, "JNKi", "DMSO", "replicated")
z <- eff$effect / pmax(abs(eff$se), 1e-6)
z[!is.finite(z)] <- 0
track <- data.frame(
  chrom = m200$bed$chrom, start = m200$bed$start, end = m200$bed$end,
  score = z, valid = TRUE, stringsAsFactors = FALSE
)
chroms <- unique(track$chrom)
sizes <- tapply(seq_len(nrow(track)), track$chrom, length)
ord <- names(sort(sizes))
pick <- unique(c(ord[1], ord[max(1, length(ord) %/% 2)], ord[length(ord)]))
# also a 25k prefix of the largest chrom for PELT
large <- ord[length(ord)]

rows <- list()
time_one <- function(y, engine, label, chrom, n_bins) {
  t0 <- proc.time()[["elapsed"]]
  rss0 <- process_rss_mb()
  st <- "PASS"
  nseg <- NA_integer_
  msg <- ""
  maxc <- NA_integer_
  tryCatch({
    setTimeLimit(elapsed = 90, transient = TRUE)
    if (engine == "pelt_cpp") {
      out <- pelt_gaussian_cpp(y, 24, 3L)
      nseg <- length(out$changepoints) + 1L
      maxc <- if (!is.null(out$max_candidates)) out$max_candidates else NA_integer_
    } else if (engine == "fpop_cpp") {
      out <- fpop_gaussian_cpp(y, 24, 3L)
      nseg <- length(out$changepoints) + 1L
      maxc <- if (!is.null(out$max_candidates)) out$max_candidates else NA_integer_
    } else {
      seg <- segment_effect(y, penalty = 24, min_width = 3L, engine = engine)
      nseg <- nrow(seg)
    }
    setTimeLimit(elapsed = Inf)
  }, error = function(e) {
    setTimeLimit(elapsed = Inf)
    st <<- if (grepl("time.?limit|elapsed", conditionMessage(e), ignore.case = TRUE))
      "TIMEOUT" else "FAILED"
    msg <<- conditionMessage(e)
  })
  data.frame(
    label = label, chrom = chrom, engine = engine, n_bins = n_bins,
    n_segments = nseg, max_candidates = maxc,
    runtime_sec = proc.time()[["elapsed"]] - t0,
    rss_mb = process_rss_mb(), rss_mb_before = rss0,
    status = st, failure_reason = msg, stringsAsFactors = FALSE
  )
}

for (ch in pick) {
  sub <- track$score[track$chrom == ch]
  message("profile ", ch, " n=", length(sub))
  rows[[length(rows) + 1L]] <- time_one(sub, "fpop_cpp", "full_chrom", ch, length(sub))
  # PELT only if reasonably small
  if (length(sub) <= 80000L) {
    rows[[length(rows) + 1L]] <- time_one(sub, "pelt_cpp", "full_chrom", ch, length(sub))
  } else {
    rows[[length(rows) + 1L]] <- data.frame(
      label = "full_chrom", chrom = ch, engine = "pelt_cpp", n_bins = length(sub),
      n_segments = NA, max_candidates = NA, runtime_sec = NA,
      rss_mb = NA, rss_mb_before = NA, status = "NOT_RUN",
      failure_reason = "n>80000; skipped to avoid v0.4.1 hang",
      stringsAsFactors = FALSE
    )
  }
}

# prefixes of the largest chrom
sub_large <- track$score[track$chrom == large]
for (n in c(5000L, 20000L, 50000L)) {
  y <- sub_large[seq_len(min(n, length(sub_large)))]
  rows[[length(rows) + 1L]] <- time_one(y, "pelt_cpp", paste0("prefix", n), large, length(y))
  rows[[length(rows) + 1L]] <- time_one(y, "fpop_cpp", paste0("prefix", n), large, length(y))
}

tab <- do.call(rbind, rows)
utils::write.table(tab, "results/v0.4.2/phase1_profile.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
print(tab)
message("wrote results/v0.4.2/phase1_profile.tsv")
