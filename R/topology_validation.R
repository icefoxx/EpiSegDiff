#' Geometry-only classification from two intervals
#'
#' Applies the v0.2.1 non-overlapping Δ_L / Δ_R / R_W partition.
#'
#' @param c_start,c_end,t_start,t_end Integer coordinates (0-based).
#' @param chrom Chromosome label.
#' @param differential Optional internal D table.
#' @param thresholds Topology thresholds.
#' @export
classify_geometry_pair <- function(c_start, c_end, t_start, t_end,
                                   chrom = "chr1", differential = NULL,
                                   thresholds = NULL) {
  th <- .topology_thresholds(thresholds)
  e1 <- data.frame(
    chrom = chrom,
    c_start = as.integer(c_start), c_end = as.integer(c_end),
    t_start = as.integer(t_start), t_end = as.integer(t_end),
    left_shift = as.integer(t_start) - as.integer(c_start),
    right_shift = as.integer(t_end) - as.integer(c_end),
    width_ratio = (as.integer(t_end) - as.integer(t_start)) /
      max(as.integer(c_end) - as.integer(c_start), 1),
    reciprocal_overlap = .reciprocal_overlap(c_start, c_end, t_start, t_end),
    distance = .gap_distance(c_start, c_end, t_start, t_end),
    stringsAsFactors = FALSE
  )
  .classify_one_to_one(e1, differential, th)
}

#' Generate labeled 1-1 geometry examples (including borderline)
#'
#' @param event_class Primary class or `COMPLEX`.
#' @param n Number of examples.
#' @param seed RNG seed.
#' @export
simulate_geometry_examples <- function(event_class, n = 200L, seed = 1L) {
  set.seed(as.integer(seed))
  n <- as.integer(n)
  delta <- 100
  rows <- vector("list", n)
  for (i in seq_len(n)) {
    wc <- sample(c(800L, 1200L, 2000L, 4000L), 1)
    lc <- sample(seq.int(2000L, 20000L, by = 50L), 1)
    rc <- lc + wc
    borderline <- (i %% 5L == 0L)
    dmag <- if (borderline) delta + sample(1:40, 1) else delta + sample(150:800, 1)
    if (event_class == "EXPANSION") {
      sided <- sample(c("both", "left", "right"), 1)
      need <- as.integer(ceiling(wc * 0.30))
      ext <- max(dmag, need)
      if (sided == "both") {
        lt <- lc - max(as.integer(ceiling(ext / 2)), delta + 20L)
        rt <- rc + max(as.integer(ceiling(ext / 2)), delta + 20L)
      } else if (sided == "left") {
        lt <- lc - ext; rt <- rc + sample(0:20, 1)
      } else {
        lt <- lc + sample(0:20, 1); rt <- rc + ext
      }
    } else if (event_class == "CONTRACTION") {
      sided <- sample(c("both", "left", "right"), 1)
      need <- as.integer(ceiling(wc * 0.25))
      shrink <- min(max(dmag, need), floor(wc * 0.40))
      if (sided == "both") {
        half <- max(as.integer(ceiling(shrink / 2)), delta + 20L)
        lt <- lc + half; rt <- rc - half
      } else if (sided == "left") {
        lt <- lc + shrink; rt <- rc - sample(0:20, 1)
      } else {
        lt <- lc + sample(0:20, 1); rt <- rc - shrink
      }
    } else if (event_class == "LEFT_SHIFT") {
      lt <- lc - dmag; rt <- rc - dmag
    } else if (event_class == "RIGHT_SHIFT") {
      lt <- lc + dmag; rt <- rc + dmag
    } else if (event_class == "INTERNAL_GAIN" || event_class == "INTERNAL_LOSS") {
      lt <- lc; rt <- rc
    } else if (event_class == "COMPLEX") {
      # Both ends move outward but width ratio stays below ρ (not expansion).
      tiny <- as.integer(max(delta + 5L, floor(wc * 0.05)))
      lt <- lc - tiny
      rt <- rc + tiny
      # R_W = 1 + 2*tiny/wc ≈ 1.10 < 1.25
    } else if (event_class == "GAIN") {
      # unmatched treatment is generated at the domain-table level
      lt <- lc; rt <- rc
    } else if (event_class == "LOSS") {
      lt <- lc; rt <- rc
    } else {
      stop("Unsupported event_class for geometry examples: ", event_class,
           call. = FALSE)
    }
    D <- NULL
    if (event_class %in% c("INTERNAL_GAIN", "INTERNAL_LOSS")) {
      ds <- lc + as.integer(round(0.3 * wc))
      de <- lc + as.integer(round(0.55 * wc))
      D <- data.frame(chrom = "chr1", start = ds, end = de,
                      direction = if (event_class == "INTERNAL_GAIN") "GAIN" else "LOSS")
    }
    pred <- classify_geometry_pair(lc, rc, lt, rt, differential = D)
    rows[[i]] <- data.frame(
      example_id = i,
      truth = event_class,
      pred = pred$primary_event[1],
      rule = pred$classification_rule[1],
      c_start = lc, c_end = rc, t_start = lt, t_end = rt,
      dL = lt - lc, dR = rt - rc,
      width_ratio = (rt - lt) / max(rc - lc, 1),
      borderline = borderline,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

#' Domain-table examples for SPLIT / MERGE / unmatched GAIN / LOSS
#' @export
simulate_match_class_examples <- function(event_class, n = 200L, seed = 2L) {
  set.seed(as.integer(seed))
  th <- .topology_thresholds(NULL)
  rows <- vector("list", n)
  for (i in seq_len(n)) {
    base <- 1000L + i * 50L
    if (event_class == "SPLIT") {
      C <- data.frame(domain_id = 1L, chrom = "chr1",
                      start = base, end = base + 5000L)
      gap <- if (i %% 5L == 0L) 220L else 400L + sample(0:200, 1)
      T <- data.frame(domain_id = 1:2, chrom = "chr1",
                      start = c(base, base + 2000L + gap),
                      end = c(base + 1800L, base + 5000L))
    } else if (event_class == "MERGE") {
      T <- data.frame(domain_id = 1L, chrom = "chr1",
                      start = base, end = base + 5000L)
      gap <- if (i %% 5L == 0L) 220L else 400L + sample(0:200, 1)
      C <- data.frame(domain_id = 1:2, chrom = "chr1",
                      start = c(base, base + 2000L + gap),
                      end = c(base + 1800L, base + 5000L))
    } else if (event_class == "GAIN") {
      C <- data.frame(domain_id = integer(), chrom = character(),
                      start = integer(), end = integer())
      T <- data.frame(domain_id = 1L, chrom = "chr1",
                      start = base, end = base + 1500L)
    } else if (event_class == "LOSS") {
      T <- data.frame(domain_id = integer(), chrom = character(),
                      start = integer(), end = integer())
      C <- data.frame(domain_id = 1L, chrom = "chr1",
                      start = base, end = base + 1500L)
    } else {
      stop("Use simulate_geometry_examples() for ", event_class, call. = FALSE)
    }
    m <- match_domains(C, T, th$min_reciprocal_overlap, th$max_match_distance,
                       th$high_reciprocal_overlap)
    ev <- classify_boundary_events(m, C, T, thresholds = th)
    rows[[i]] <- data.frame(
      example_id = i, truth = event_class,
      pred = ev$primary_event[1],
      rule = ev$classification_rule[1],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}
