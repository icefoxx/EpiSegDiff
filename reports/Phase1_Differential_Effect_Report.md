# Phase 1 / Gate 1 — Differential-effect report

**Objective.** Show that a replicate-aware per-bin NB (edgeR QL) fit
recovers known synthetic fold-changes, is centered near zero on null
data, is invariant to sample order, and reverses sign when condition
labels are swapped.

**Inputs.** In-memory NB simulations (`simulate_count_assay`,
`simulate_null_counts`). No real BAM cohort in this checkout.

**Methods.** `normalize_counts(method = "library_size"| "TMM")` then
`estimate_differential_effect()`. Score default:
`standardized_effect = logFC / SE`, with SE recovered from the QL F
statistic (`t = sign(logFC)*sqrt(F)`). Dispersion is estimated on bins
before any segmentation.

**Parameters.** See `config/example.yaml`. Tests use `seed` recorded in
each `simulate_*` call.

**Tests.** `tests/testthat/test-differential-effect.R`,
`test-normalization.R`, `test-count-matrix.R`.

**Results.**

- Known GAIN (`effect = 2` on bins 80–120, 200 bins, 3 vs 3, depth 40):
  median logFC in the planted interval > 1 and exceeds the background
  median by > 1.
- Null median \(Z\) is near 0 (test threshold `|mean(Z)| < 0.4`).
- Reordering replicates does not change `logFC` beyond `1e-6`.
- Swapping control/treatment negates `logFC` within `1e-6`.
- Masked bins remain `valid = FALSE` and `score = NA` (no interpolation).
- Library-size and TMM factors are finite, reproducible, and warn on
  extreme ratios. Spike-in / background methods error instead of
  silently pretending to work.
- Synthetic paired BAM: two proper-pair fragments at midpoints 175 and
  450 are counted into the expected 100-bp bins.

**Failures.** None in the current `testthat` suite.

**Limitations.**

- Library-size normalization on a genome with a large GAIN induces a
  composition effect: background logFC is pulled slightly opposite the
  event. This is documented in `docs/known_limitations.md` and is why
  recovery tests compare event vs background rather than requiring
  background exactly 0.
- `stdev.unscaled` is not present in edgeR 4.0 QL objects; SE is derived
  from the QL F statistic. This is explicit, not hidden.
- No real-chromosome BigWig was written in this session
  (`write_bigwig: false` unless rtracklayer export is requested).
- Dispersion plots are produced by `plot_dispersion()` when a user
  supplies an output path; they are not committed as binary artifacts.

**Files created.**

- `R/bam_qc.R`, `R/count_matrix.R`, `R/normalization.R`,
  `R/differential_effect.R`
- `reports/normalization_diagnostics.tsv`

**Decision: PASS**

**Next action.** Segmentation, recount, and empirical FDR (Phases 3–6).
The null-score distribution in unit tests does not show strong
systematic bias after allowing for the documented composition effect.
