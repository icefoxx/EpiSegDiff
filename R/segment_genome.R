#' Process RSS of the current R session in megabytes
#' @export
process_rss_mb <- function() {
  pid <- Sys.getpid()
  out <- tryCatch(
    suppressWarnings(system(sprintf("ps -o rss= -p %s", pid), intern = TRUE)),
    error = function(e) NA_character_
  )
  mb <- suppressWarnings(as.numeric(out) / 1024)
  if (!length(mb) || !is.finite(mb[1])) return(NA_real_)
  mb[1]
}

.checkpoint_identity <- function(identity) {
  list(
    cache_id = identity$cache_id %||% NA_character_,
    genome = identity$genome %||% NA_character_,
    bin_size = as.integer(identity$bin_size %||% NA_integer_),
    summary_fun = identity$summary_fun %||% NA_character_,
    transform = identity$transform %||% NA_character_,
    pseudocount = as.numeric(identity$pseudocount %||% NA_real_),
    filter = identity$filter %||% NA_character_,
    segmentation = identity$segmentation %||% "gaussian_pelt",
    penalty = as.numeric(identity$penalty %||% NA_real_),
    min_bins = as.integer(identity$min_bins %||% NA_integer_),
    engine = identity$engine %||% NA_character_,
    software = as.character(identity$software %||%
                              utils::packageVersion("EpiSegDiff"))
  )
}

.identity_compatible <- function(saved, want) {
  keys <- c("cache_id", "genome", "bin_size", "summary_fun", "transform",
            "pseudocount", "filter", "segmentation", "penalty", "min_bins")
  for (k in keys) {
    a <- saved[[k]]; b <- want[[k]]
    if (is.null(a) || is.null(b)) return(FALSE)
    if (is.numeric(a) && is.numeric(b)) {
      if (!isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-12))) {
        return(FALSE)
      }
    } else if (!identical(as.character(a), as.character(b))) {
      return(FALSE)
    }
  }
  TRUE
}

.safe_chrom_token <- function(chrom) {
  gsub("[^A-Za-z0-9._-]+", "_", as.character(chrom))
}

#' Split a bin track by chromosome without concatenating genomes
#' @export
split_track_by_chromosome <- function(track) {
  chroms <- unique(as.character(track$chrom))
  out <- lapply(chroms, function(ch) {
    track[as.character(track$chrom) == ch, , drop = FALSE]
  })
  names(out) <- chroms
  out
}

#' Segment one chromosome through genuine contiguous valid blocks
#'
#' Blocks are defined only by invalid/non-finite scores or coordinate
#' discontinuities. No arbitrary computational chunks.
#' @export
segment_chromosome <- function(track, penalty, min_width = 3L, ...) {
  if (!nrow(track)) {
    out <- .empty_segments()
    out$chrom <- character()
    out$cvr_id <- integer()
    return(out)
  }
  chrom <- unique(as.character(track$chrom))
  if (length(chrom) != 1L) {
    stop("segment_chromosome() requires a single chromosome.", call. = FALSE)
  }
  segment_track(track, penalty = penalty, min_width = min_width, ...)
}

.read_chrom_checkpoint <- function(path, want) {
  if (!file.exists(path)) return(NULL)
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj) || is.null(obj$identity)) {
    stop("Corrupt checkpoint: ", path, call. = FALSE)
  }
  if (!.identity_compatible(obj$identity, want)) {
    stop("Incompatible checkpoint ", path, call. = FALSE)
  }
  obj
}

.write_chrom_checkpoint <- function(path, obj) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

#' Segment a genome chromosome-by-chromosome with optional checkpoints
#'
#' @param track data.frame with `chrom`,`start`,`end`,`score`,`valid`.
#' @param penalty,min_width Passed to [segment_effect()].
#' @param checkpoint_dir Directory for `chr*.rds` checkpoints, or NULL.
#' @param resume Reuse compatible completed chromosome checkpoints.
#' @param workers Number of chromosome-level workers (1 = sequential).
#' @param identity Named list recorded in every checkpoint.
#' @param engine `pelt`, `fpop`, or `auto` (fpop when a block has
#'   >= `auto_n` bins).
#' @param auto_n Size threshold for `engine="auto"`.
#' @export
segment_genome <- function(track, penalty, min_width = 3L,
                           checkpoint_dir = NULL, resume = TRUE,
                           workers = 1L, identity = list(),
                           engine = c("auto", "pelt", "fpop"),
                           auto_n = 5000000L, ...) {
  engine <- match.arg(engine)
  workers <- max(1L, as.integer(workers))
  want <- .checkpoint_identity(identity)
  want$engine <- if (engine == "auto") "auto" else engine
  chroms <- unique(as.character(track$chrom))
  extra <- list(...)

  run_one <- function(ch) {
    t0 <- proc.time()[["elapsed"]]
    rss0 <- process_rss_mb()
    tok <- .safe_chrom_token(ch)
    cpath <- if (!is.null(checkpoint_dir)) {
      file.path(checkpoint_dir, paste0("chr", tok, ".rds"))
    } else {
      NULL
    }
    if (isTRUE(resume) && !is.null(cpath)) {
      hit <- tryCatch(.read_chrom_checkpoint(cpath, want), error = function(e) {
        if (grepl("Incompatible|Corrupt", conditionMessage(e))) stop(e)
        NULL
      })
      if (!is.null(hit)) {
        hit$runtime$resumed <- TRUE
        return(hit)
      }
    }
    sub <- track[as.character(track$chrom) == ch, , drop = FALSE]
    n <- nrow(sub)
    n_valid <- sum(sub$valid & is.finite(sub$score))
    use_engine <- engine
    if (identical(engine, "auto")) {
      use_engine <- if (n_valid >= as.integer(auto_n)) "fpop" else "pelt"
    }
    args <- c(list(track = sub, penalty = penalty, min_width = min_width,
                   loss = extra$loss %||% "gaussian", engine = use_engine),
              extra[setdiff(names(extra), c("loss", "engine"))])
    segs <- do.call(segment_chromosome, args)
    rss1 <- process_rss_mb()
    obj <- list(
      chrom = ch,
      segments = segs,
      identity = want,
      runtime = data.frame(
        chrom = ch,
        n_bins = n,
        n_valid = n_valid,
        n_segments = nrow(segs),
        engine = use_engine,
        runtime_sec = proc.time()[["elapsed"]] - t0,
        rss_mb_before = rss0,
        rss_mb_after = rss1,
        resumed = FALSE,
        status = "PASS",
        stringsAsFactors = FALSE
      )
    )
    if (!is.null(cpath)) .write_chrom_checkpoint(cpath, obj)
    rm(sub)
    obj
  }

  if (workers == 1L) {
    parts <- lapply(chroms, run_one)
  } else {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(workers)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, varlist = c("track", "penalty", "min_width",
                                              "checkpoint_dir", "resume",
                                              "want", "engine", "auto_n",
                                              "extra", "run_one"),
                              envir = environment())
      parallel::clusterEvalQ(cl, {
        if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
        NULL
      })
      parts <- parallel::parLapply(cl, chroms, run_one)
    } else {
      old <- getOption("mc.cores", 1L)
      options(mc.cores = 1L)
      on.exit(options(mc.cores = old), add = TRUE)
      parts <- parallel::mclapply(chroms, run_one, mc.cores = workers,
                                  mc.preschedule = FALSE)
      err <- vapply(parts, inherits, logical(1), what = "try-error")
      if (any(err)) {
        stop("Worker failed: ", parts[[which(err)[1]]], call. = FALSE)
      }
    }
  }
  names(parts) <- chroms
  segs <- Filter(function(x) nrow(x$segments) > 0, parts)
  combined <- if (!length(segs)) {
    out <- .empty_segments()
    out$chrom <- character()
    out$cvr_id <- integer()
    out
  } else {
    do.call(rbind, lapply(segs, `[[`, "segments"))
  }
  rt <- do.call(rbind, lapply(parts, `[[`, "runtime"))
  rownames(rt) <- NULL
  if (nrow(combined)) rownames(combined) <- NULL
  list(segments = combined, chromosome_runtime = rt, identity = want)
}

#' Choose a Gaussian engine that preserves PELT semantics
#' @export
choose_gaussian_engine <- function(n, engine = "auto", auto_n = 5000000L) {
  engine <- match.arg(engine, c("auto", "pelt", "fpop"))
  if (identical(engine, "auto")) {
    if (as.integer(n) >= as.integer(auto_n)) "fpop" else "pelt"
  } else {
    engine
  }
}
