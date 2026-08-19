#' Read a methylation count-pair table
#'
#' Required columns: chrom, position or start, methylated_count, total_count,
#' sample. `condition` / `replicate` may be joined from sample metadata.
#'
#' @param path TSV path or a data.frame.
#' @param format `tabular`, `bismark_coverage`, or `bedgraph`.
#' @param min_coverage Sites below this become invalid (not imputed).
#' @return data.frame with the locked methylation representation.
#' @export
read_methylation_table <- function(path, format = c("tabular", "bismark_coverage",
                                                    "bedgraph"),
                                   min_coverage = 0L) {
  format <- match.arg(format)
  if (is.data.frame(path)) {
    df <- path
  } else {
    df <- utils::read.table(path, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, comment.char = "")
  }
  if (format == "bismark_coverage") {
    # chrom start end pct meth unmeth  (start may be 1-based; we store 0-based)
    nms <- tolower(names(df))
    if (ncol(df) < 6L && !all(c("chr", "start", "meth", "unmeth") %in% nms)) {
      names(df)[1:6] <- c("chrom", "start", "end", "pct", "meth", "unmeth")
    }
    df$chrom <- df[[1]]
    start <- as.integer(df[[2]])
    # If typical Bismark 1-based, convert when start>=1 and end==start
    df$position <- pmax(0L, start - 1L)
    df$methylated_count <- as.integer(df[[5]])
    df$total_count <- as.integer(df[[5]]) + as.integer(df[[6]])
  } else if (format == "bedgraph") {
    names(df)[1:4] <- c("chrom", "start", "end", "methylation_fraction")
    df$position <- as.integer(df$start)
    if (!"methylated_count" %in% names(df)) df$methylated_count <- NA_integer_
    if (!"total_count" %in% names(df)) df$total_count <- NA_integer_
  } else {
    if (!"position" %in% names(df)) {
      if ("start" %in% names(df)) df$position <- as.integer(df$start) else
        stop("tabular methylation input needs position or start.", call. = FALSE)
    }
  }

  need <- c("chrom", "position", "methylated_count", "total_count")
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    stop("Methylation table missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  df$chrom <- as.character(df$chrom)
  df$position <- as.integer(df$position)
  df$methylated_count <- as.integer(df$methylated_count)
  df$total_count <- as.integer(df$total_count)
  if (any(df$methylated_count > df$total_count, na.rm = TRUE)) {
    stop("methylated_count exceeds total_count.", call. = FALSE)
  }
  df$coverage <- df$total_count
  df$methylation_fraction <- ifelse(df$total_count > 0,
                                    df$methylated_count / df$total_count, NA_real_)
  df$valid <- is.finite(df$total_count) & df$total_count >= as.integer(min_coverage)
  # Missing CpGs are absent rows, never zero methylation.
  df
}

#' Estimate a per-site methylation differential score
#'
#' Uses a method-of-moments delta-beta and a beta-binomial-inspired SE.
#' Does **not** call edgeR on methylated counts.
#'
#' @param meth Long table from [read_methylation_table()] with `sample`.
#' @param samples Sample metadata.
#' @param control_level,treatment_level Condition names.
#' @param score `delta_beta` or `standardized_effect`.
#' @param min_coverage Coverage filter.
#' @return Site-level track (`chrom`, `start`, `end`, `effect`, `uncertainty`,
#'   `score`, `valid`).
#' @export
estimate_methylation_effect <- function(meth, samples,
                                        control_level = "control",
                                        treatment_level = "treatment",
                                        score = c("delta_beta", "standardized_effect"),
                                        min_coverage = 5L) {
  score <- match.arg(score)
  if (!"sample" %in% names(meth) && "sample_id" %in% names(meth)) {
    meth$sample <- meth$sample_id
  }
  if (!"sample" %in% names(meth)) {
    stop("Methylation table needs a sample column.", call. = FALSE)
  }
  meth <- merge(meth, samples[, c("sample_id", "condition", "replicate")],
                by.x = "sample", by.y = "sample_id", all.x = TRUE)
  if (any(is.na(meth$condition))) {
    stop("Some methylation rows have sample IDs absent from metadata.", call. = FALSE)
  }
  meth$valid <- meth$valid & meth$total_count >= as.integer(min_coverage)

  key <- paste(meth$chrom, meth$position, sep = ":")
  sites <- unique(key)
  out <- data.frame(
    chrom = sub(":.*$", "", sites),
    start = as.integer(sub("^.*:", "", sites)),
    end = as.integer(sub("^.*:", "", sites)) + 1L,
    effect = NA_real_,
    uncertainty = NA_real_,
    score = NA_real_,
    valid = FALSE,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(sites)) {
    rows <- meth[key == sites[i] & meth$valid, , drop = FALSE]
    ctrl <- rows[rows$condition == control_level, ]
    trt <- rows[rows$condition == treatment_level, ]
    if (!nrow(ctrl) || !nrow(trt)) next
    p_c <- sum(ctrl$methylated_count) / sum(ctrl$total_count)
    p_t <- sum(trt$methylated_count) / sum(trt$total_count)
    if (!is.finite(p_c) || !is.finite(p_t)) next
    db <- p_t - p_c
    # Method-of-moments SE: binomial as a lower bound; inflate by a crude
    # overdispersion using replicate p if ≥2 replicates.
    se_c <- .meth_se(ctrl)
    se_t <- .meth_se(trt)
    se <- sqrt(se_c^2 + se_t^2)
    z <- if (is.finite(se) && se > 0) db / se else NA_real_
    out$effect[i] <- db
    out$uncertainty[i] <- se
    out$score[i] <- if (score == "delta_beta") db else z
    out$valid[i] <- is.finite(out$score[i])
  }
  out[order(out$chrom, out$start), ]
}

.meth_se <- function(rows) {
  n <- sum(rows$total_count)
  if (n <= 1) return(NA_real_)
  p <- sum(rows$methylated_count) / n
  se_bin <- sqrt(p * (1 - p) / n)
  if (nrow(rows) >= 2L) {
    pr <- rows$methylated_count / pmax(rows$total_count, 1L)
    se_rep <- stats::sd(pr) / sqrt(nrow(rows))
    if (is.finite(se_rep)) return(max(se_bin, se_rep))
  }
  se_bin
}

#' Test candidate DMRs with a methylation-specific model
#'
#' `beta_binomial` uses a per-region two-group beta-binomial GLM
#' approximation (logit + overdispersed binomial). `dss_style` uses the
#' same statistic and records the method name; DSS itself is optional.
#'
#' @param regions Region table with per-sample `M`/`N` columns named
#'   `meth_<sample>` and `cov_<sample>`, **or** a list `methylated`,
#'   `coverage` matrices.
#' @param samples Sample metadata.
#' @param method `beta_binomial` or `dss_style`.
#' @param control_level,treatment_level Condition names.
#' @param fdr_method p-adjust method.
#' @export
test_dmrs <- function(regions, samples, method = c("beta_binomial", "dss_style"),
                      control_level = "control", treatment_level = "treatment",
                      fdr_method = "BH", ...) {
  method <- match.arg(method)
  if (is.list(regions) && !is.data.frame(regions) &&
      all(c("methylated", "coverage") %in% names(regions))) {
    M <- as.matrix(regions$methylated)
    N <- as.matrix(regions$coverage)
    meta <- regions$meta %||% data.frame(region_id = seq_len(nrow(M)))
  } else if (is.data.frame(regions)) {
    sids <- samples$sample_id
    mcols <- paste0("meth_", sids)
    ncols <- paste0("cov_", sids)
    if (!all(mcols %in% names(regions)) || !all(ncols %in% names(regions))) {
      stop("DMR table needs meth_<sample> and cov_<sample> columns.", call. = FALSE)
    }
    M <- as.matrix(regions[, mcols, drop = FALSE])
    N <- as.matrix(regions[, ncols, drop = FALSE])
    colnames(M) <- colnames(N) <- sids
    meta <- regions
  } else {
    stop("Unsupported regions argument for test_dmrs().", call. = FALSE)
  }

  cond <- .condition_factor(samples, control_level, treatment_level)
  nreg <- nrow(M)
  pval <- logFC <- stat <- rep(NA_real_, nreg)
  for (i in seq_len(nreg)) {
    fit <- tryCatch(
      .bb_region_test(M[i, ], N[i, ], cond),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    pval[i] <- fit$pvalue
    logFC[i] <- fit$delta_beta
    stat[i] <- fit$statistic
  }
  data.frame(
    region_id = if ("region_id" %in% names(meta)) meta$region_id else seq_len(nreg),
    logFC_or_delta = logFC,
    statistic = stat,
    pvalue = pval,
    FDR = stats::p.adjust(pval, method = fdr_method),
    testing_method = method,
    fdr_method = fdr_method,
    fdr_interpretation = "nominal BH after adaptive discovery; not guaranteed selective FDR",
    stringsAsFactors = FALSE
  )
}

.bb_region_test <- function(m, n, cond) {
  keep <- is.finite(m) & is.finite(n) & n > 0
  if (sum(keep) < 3L || length(unique(cond[keep])) < 2L) {
    stop("insufficient")
  }
  m <- m[keep]; n <- n[keep]; cond <- droplevels(cond[keep])
  p <- m / n
  # Weighted two-group z-test on p with extra-binomial variance.
  g0 <- cond == levels(cond)[1]
  g1 <- !g0
  p0 <- stats::weighted.mean(p[g0], n[g0])
  p1 <- stats::weighted.mean(p[g1], n[g1])
  # moment overdispersion
  phi <- 0
  if (sum(g0) > 1L) phi <- max(phi, stats::var(p[g0]))
  if (sum(g1) > 1L) phi <- max(phi, stats::var(p[g1]))
  v0 <- p0 * (1 - p0) / sum(n[g0]) + phi / max(sum(g0), 1)
  v1 <- p1 * (1 - p1) / sum(n[g1]) + phi / max(sum(g1), 1)
  se <- sqrt(max(v0 + v1, 1e-12))
  z <- (p1 - p0) / se
  pvalue <- 2 * stats::pnorm(-abs(z))
  list(delta_beta = p1 - p0, statistic = z, pvalue = pvalue)
}
