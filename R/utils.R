# Internal utilities. Not a scientific API.

.need_pkg <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for %s.", pkg, why), call. = FALSE)
  }
}

.is_finite_number <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

.assert_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("'%s' must be a non-missing logical scalar.", name), call. = FALSE)
  }
}

.null_or_path <- function(x) {
  is.null(x) || (is.character(x) && length(x) == 1L)
}

`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

.empty_df <- function(cols) {
  out <- data.frame(matrix(ncol = length(cols), nrow = 0), check.names = FALSE)
  names(out) <- cols
  out
}

.as_0based_int <- function(x) {
  as.integer(x)
}

# BED-style overlap: [start, end)
.recycle4 <- function(as, ae, bs, be) {
  n <- max(length(as), length(ae), length(bs), length(be))
  list(
    as = rep(as, length.out = n),
    ae = rep(ae, length.out = n),
    bs = rep(bs, length.out = n),
    be = rep(be, length.out = n)
  )
}

.interval_intersection <- function(as, ae, bs, be) {
  x <- .recycle4(as, ae, bs, be)
  pmax(0L, pmin(x$ae, x$be) - pmax(x$as, x$bs))
}

.reciprocal_overlap <- function(as, ae, bs, be) {
  x <- .recycle4(as, ae, bs, be)
  wa <- x$ae - x$as
  wb <- x$be - x$bs
  inter <- pmax(0L, pmin(x$ae, x$be) - pmax(x$as, x$bs))
  ro <- rep(0, length(wa))
  ok <- wa > 0 & wb > 0
  ro[ok] <- pmin(inter[ok] / wa[ok], inter[ok] / wb[ok])
  ro
}

.gap_distance <- function(as, ae, bs, be) {
  inter <- .interval_intersection(as, ae, bs, be)
  d <- pmin(abs(as - be), abs(bs - ae))
  d[inter > 0] <- 0L
  as.integer(d)
}

.warn_once_env <- new.env(parent = emptyenv())

.warn_once <- function(key, msg) {
  if (isTRUE(.warn_once_env[[key]])) return(invisible())
  .warn_once_env[[key]] <- TRUE
  warning(msg, call. = FALSE)
}

.digest_chr <- function(x) {
  bytes <- charToRaw(paste(x, collapse = "\n"))
  sprintf("%x", sum(as.integer(bytes)) + length(bytes) * 31L)
}

# Stable identifier: sha256 if digest is available, else size+mtime+path.
file_stable_id <- function(path) {
  if (is.null(path) || !nzchar(as.character(path)[1])) {
    return(NA_character_)
  }
  path <- as.character(path)[1]
  if (!file.exists(path)) {
    return(paste0("missing:", path))
  }
  info <- file.info(path)
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256"))
  }
  sprintf("path=%s;size=%s;mtime=%s",
          normalizePath(path, winslash = "/", mustWork = FALSE),
          as.character(info$size),
          as.character(as.numeric(info$mtime)))
}

.sort_intervals <- function(df) {
  if (nrow(df) == 0L) return(df)
  df[order(df$chrom, df$start, df$end), , drop = FALSE]
}

.chrom_levels <- function(x) {
  as.character(x)
}

.check_bed_frame <- function(df, name) {
  need <- c("chrom", "start", "end")
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    stop(sprintf("%s is missing columns: %s", name, paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  df$chrom <- as.character(df$chrom)
  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)
  if (any(df$end < df$start, na.rm = TRUE)) {
    stop(sprintf("%s has end < start.", name), call. = FALSE)
  }
  df
}

# Prefix sums with a leading 0 so sums[s:t] = P[t] - P[s-1] using 1-based s,t.
.prefix_sum <- function(x) {
  c(0, cumsum(as.numeric(x)))
}
