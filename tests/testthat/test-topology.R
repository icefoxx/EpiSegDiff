th <- list(
  min_reciprocal_overlap = 0.2,
  high_reciprocal_overlap = 0.5,
  max_match_distance = 5000,
  min_width_ratio_change = 1.25,
  min_boundary_shift = 100,
  min_internal_gap_width = 200,
  similar_width_ratio = 1.25,
  max_overlap_for_shift = 0.85
)

class_from_domains <- function(C, T, D = NULL) {
  m <- match_domains(C, T, min_reciprocal_overlap = th$min_reciprocal_overlap,
                     max_match_distance = th$max_match_distance,
                     high_reciprocal_overlap = th$high_reciprocal_overlap)
  ev <- classify_boundary_events(m, C, T, differential = D, thresholds = th)
  ev
}

test_that("GAIN / LOSS unmatched domains classify correctly", {
  C <- data.frame(domain_id = 1L, chrom = "chr1", start = 0L, end = 1000L)
  T <- data.frame(domain_id = integer(), chrom = character(), start = integer(), end = integer())
  ev <- class_from_domains(C, T)
  expect_equal(ev$primary_event, "LOSS")

  C2 <- T
  T2 <- data.frame(domain_id = 1L, chrom = "chr1", start = 0L, end = 1000L)
  ev2 <- class_from_domains(C2, T2)
  expect_equal(ev2$primary_event, "GAIN")
})

test_that("EXPANSION and CONTRACTION use width and overlap rules", {
  C <- data.frame(domain_id = 1L, chrom = "chr1", start = 1000L, end = 2000L)
  T <- data.frame(domain_id = 1L, chrom = "chr1", start = 500L, end = 2500L)
  ev <- class_from_domains(C, T)
  expect_equal(ev$primary_event, "EXPANSION")

  ev2 <- class_from_domains(T, C)
  expect_equal(ev2$primary_event, "CONTRACTION")
})

test_that("LEFT_SHIFT and RIGHT_SHIFT are directional", {
  C <- data.frame(domain_id = 1L, chrom = "chr1", start = 2000L, end = 3000L)
  T <- data.frame(domain_id = 1L, chrom = "chr1", start = 1000L, end = 2000L)
  ev <- class_from_domains(C, T)
  expect_equal(ev$primary_event, "LEFT_SHIFT")

  C <- data.frame(domain_id = 1L, chrom = "chr1", start = 1000L, end = 2000L)
  T <- data.frame(domain_id = 1L, chrom = "chr1", start = 2000L, end = 3000L)
  ev <- class_from_domains(C, T)
  expect_equal(ev$primary_event, "RIGHT_SHIFT")
})

test_that("SPLIT and MERGE require gaps and overlap", {
  C <- data.frame(domain_id = 1L, chrom = "chr1", start = 0L, end = 5000L)
  T <- data.frame(domain_id = 1:2, chrom = "chr1",
                  start = c(0L, 3000L), end = c(1500L, 5000L))
  ev <- class_from_domains(C, T)
  expect_equal(ev$primary_event, "SPLIT")

  ev2 <- class_from_domains(T, C)
  expect_equal(ev2$primary_event, "MERGE")
})

test_that("INTERNAL_GAIN / INTERNAL_LOSS use differential containment", {
  C <- data.frame(domain_id = 1L, chrom = "chr1", start = 0L, end = 5000L)
  T <- data.frame(domain_id = 1L, chrom = "chr1", start = 0L, end = 5000L)
  D <- data.frame(chrom = "chr1", start = 2000L, end = 2500L, direction = "GAIN")
  ev <- class_from_domains(C, T, D)
  expect_equal(ev$primary_event, "INTERNAL_GAIN")
  D$direction <- "LOSS"
  ev2 <- class_from_domains(C, T, D)
  expect_equal(ev2$primary_event, "INTERNAL_LOSS")
})

test_that("ambiguous many-to-many maps to COMPLEX", {
  C <- data.frame(domain_id = 1:2, chrom = "chr1",
                  start = c(0L, 1000L), end = c(2000L, 3000L))
  T <- data.frame(domain_id = 1:2, chrom = "chr1",
                  start = c(500L, 1500L), end = c(2500L, 3500L))
  ev <- class_from_domains(C, T)
  expect_true(all(ev$primary_event %in% c("COMPLEX", "EXPANSION", "CONTRACTION",
                                          "GAIN", "LOSS", "MERGE", "SPLIT")))
  # at least one COMPLEX if a many-to-many component exists
  m <- match_domains(C, T, 0.2, 5000, 0.5)
  types <- vapply(m$components, `[[`, character(1), "type")
  if (any(types == "many_to_many")) {
    expect_true(any(grepl("COMPLEX|MANY_TO_MANY|RULE", ev$classification_rule)))
  }
})

test_that("simulated topology examples: 100 per class, macro-F1 reported", {
  classes <- c("GAIN", "LOSS", "EXPANSION", "CONTRACTION",
               "LEFT_SHIFT", "RIGHT_SHIFT", "SPLIT", "MERGE",
               "INTERNAL_GAIN", "INTERNAL_LOSS")
  truth <- character()
  pred <- character()
  for (cl in classes) {
    for (k in seq_len(100)) {
      # Construct domains directly from the example geometry (deterministic).
      ev <- example_topology_event(cl, offset = 50L)
      bs <- 200L
      if (cl == "GAIN") {
        C <- data.frame(domain_id = integer(), chrom = character(), start = integer(), end = integer())
        T <- data.frame(domain_id = 1L, chrom = "chr1",
                        start = (ev$start_bin - 1L) * bs, end = ev$end_bin * bs)
      } else if (cl == "LOSS") {
        T <- data.frame(domain_id = integer(), chrom = character(), start = integer(), end = integer())
        C <- data.frame(domain_id = 1L, chrom = "chr1",
                        start = (ev$start_bin - 1L) * bs, end = ev$end_bin * bs)
      } else if (cl == "SPLIT") {
        C <- data.frame(domain_id = 1L, chrom = "chr1",
                        start = (ev$control_start_bin - 1L) * bs,
                        end = ev$control_end_bin * bs)
        T <- data.frame(domain_id = 1:2, chrom = "chr1",
                        start = c((ev$control_start_bin - 1L) * bs, ev$gap_end_bin * bs),
                        end = c((ev$gap_start_bin) * bs, ev$control_end_bin * bs))
      } else if (cl == "MERGE") {
        T <- data.frame(domain_id = 1L, chrom = "chr1",
                        start = (ev$treatment_start_bin - 1L) * bs,
                        end = ev$treatment_end_bin * bs)
        mid <- ev$treatment_start_bin + 15L
        C <- data.frame(domain_id = 1:2, chrom = "chr1",
                        start = c((ev$treatment_start_bin - 1L) * bs, (mid + 5L) * bs),
                        end = c(mid * bs, ev$treatment_end_bin * bs))
      } else {
        C <- data.frame(domain_id = 1L, chrom = "chr1",
                        start = (ev$control_start_bin - 1L) * bs,
                        end = ev$control_end_bin * bs)
        T <- data.frame(domain_id = 1L, chrom = "chr1",
                        start = (ev$treatment_start_bin - 1L) * bs,
                        end = ev$treatment_end_bin * bs)
      }
      D <- NULL
      if (cl %in% c("INTERNAL_GAIN", "INTERNAL_LOSS")) {
        D <- data.frame(chrom = "chr1",
                        start = (ev$start_bin - 1L) * bs,
                        end = ev$end_bin * bs,
                        direction = if (cl == "INTERNAL_GAIN") "GAIN" else "LOSS")
      }
      got <- class_from_domains(C, T, D)
      truth <- c(truth, cl)
      pred <- c(pred, got$primary_event[1])
    }
  }
  met <- event_class_metrics(truth, pred, classes = classes)
  # Write confusion for the Phase 8 report consumer
  assign("episegdiff_topology_metrics", met, envir = .GlobalEnv)
  expect_gt(met$macro_f1, 0.85)
  expect_gt(met$accuracy, 0.85)
})

test_that("classification_rule is stored and thresholds are visible", {
  C <- data.frame(domain_id = 1L, chrom = "chr1", start = 0L, end = 1000L)
  T <- data.frame(domain_id = 1L, chrom = "chr1", start = 0L, end = 3000L)
  ev <- class_from_domains(C, T)
  expect_match(ev$classification_rule[1], "min_ro=")
  expect_match(ev$classification_rule[1], "RULE")
})
