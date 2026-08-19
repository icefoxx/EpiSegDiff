test_that("library_size factors are reproducible and finite", {
  set.seed(1)
  counts <- matrix(c(10, 20, 40, 80, 5, 10, 20, 40), nrow = 4)
  colnames(counts) <- c("a", "b")
  n1 <- normalize_counts(counts, method = "library_size")
  n2 <- normalize_counts(counts, method = "library_size")
  expect_equal(n1$size_factors, n2$size_factors)
  expect_true(all(is.finite(n1$size_factors)))
  expect_true(all(n1$size_factors > 0))
})

test_that("TMM factors are finite when edgeR is present", {
  skip_if_no_edger()
  set.seed(1)
  counts <- matrix(rnbinom(200, mu = 30, size = 10), nrow = 100, ncol = 2)
  colnames(counts) <- c("a", "b")
  n <- normalize_counts(counts, method = "TMM")
  expect_true(all(is.finite(n$size_factors)))
  expect_true(all(n$size_factors > 0))
})

test_that("zero library size is an error, not a silent NA", {
  counts <- matrix(c(0, 0, 5, 5), nrow = 2)
  colnames(counts) <- c("a", "b")
  expect_error(normalize_counts(counts, method = "library_size"), "library size 0")
})

test_that("extreme factor ratio warns", {
  counts <- matrix(c(1, 1, 1, 1000, 1000, 1000), nrow = 3)
  colnames(counts) <- c("a", "b")
  expect_warning(normalize_counts(counts, method = "library_size",
                                  extreme_factor_ratio = 5),
                 "Extreme size-factor")
})

test_that("spikein without reference bins errors; background uses low-count bins", {
  counts <- matrix(1:4, nrow = 2)
  colnames(counts) <- c("a", "b")
  expect_error(normalize_counts(counts, method = "spikein", spikein_usage = TRUE),
               "reference_bins")
  set.seed(1)
  big <- matrix(rnbinom(400, mu = 20, size = 10), nrow = 200, ncol = 2)
  colnames(big) <- c("a", "b")
  n <- normalize_counts(big, method = "background", background_definition = 0.5)
  expect_true(all(is.finite(n$size_factors)))
  expect_true(all(n$size_factors > 0))
})

test_that("diagnostics TSV has the required columns", {
  counts <- matrix(c(10, 20, 30, 40), nrow = 2)
  colnames(counts) <- c("s1", "s2")
  n <- normalize_counts(counts, method = "library_size")
  p <- tempfile(fileext = ".tsv")
  df <- write_normalization_diagnostics(n, p)
  expect_true(all(c("sample", "raw_library_size", "effective_library_size",
                    "normalization_factor", "method") %in% names(df)))
  expect_true(file.exists(p))
})
