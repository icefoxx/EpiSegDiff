# Phase 3 — Simulation benchmark (count assay)

**Objective.** Recover planted GAIN/LOSS/NULL events with the shared
interval matcher and report region, boundary, and runtime metrics.

**Inputs.** `simulate_count_assay()` with `example_topology_event()`
motifs. 200 bins × 200 bp, 3 vs 3 replicates, depth 25, \(\phi=0.08\),
5 replicates per class.

**Methods.** Same discovery path as Gate 2. Predictions = tested
regions with `FDR <= 0.05`. Matcher:
`compute_region_metrics(min_reciprocal_overlap = 0.2)`.

**Parameters.** `selected_penalty = 8`, `library_size` offsets,
Gaussian PELT.

**Tests.** Unit segmentation/topology tests plus this grid.

**Results.** Written to `benchmark/simulation/results.tsv` and
`benchmark/simulation/truth.tsv`.

| Class | Mean interval F1 | Mean recall | Mean precision | Notes |
|-------|------------------|-------------|----------------|-------|
| GAIN | 0.41 | 1.00 | 0.20–0.33 | Truth interval matched (IoU = 1); extra significant flanks |
| LOSS | 0.50 | 1.00 | 0.33 | Same pattern |
| NULL | 1.00* | 1.00* | 1.00* | No discoveries; empty-vs-empty scores as 1 |

\*Empty prediction vs empty truth is defined as perfect in
`compute_region_metrics()`. Combined with Gate 2 (0/27 null
discoveries at \(q=0.05\)), this is conservative rather than a
vacuous claim.

Matched start/end MAE on the planted interval was 0 bp (bin-aligned
truth). Hausdorff distance is large (~20 kb) **because extra predicted
intervals exist**, not because the matched event is mis-placed.

**Failures.** None hidden. Extra flank discoveries are reported. They
are consistent with library-size composition bias (treatment libraries
absorb the GAIN, so background looks like a weak LOSS and is sometimes
called).

**Limitations.**

- Width classes 200–500 bp / 1–5 kb / 20–100 kb were **not** fully
  crossed in this first grid (only the example motif).
- EXPANSION / SHIFT / SPLIT / MERGE **geometry recovery** is unit-tested
  on domains (Phase 8: 100 examples/class, macro-F1 > 0.85) but not yet
  end-to-end through BAM-scale NB simulation.
- csaw / MACS3+DiffBind / THOR were **not executed**. Logged as
  `not-attempted-in-this-checkout` in
  `benchmark/competitors/runs.tsv`.
- Runtime is ~0.03–0.05 s per 200-bin toy chromosome; not a
  genome-wide scaling measurement. No complexity claim.

**Files created.**

- `benchmark/simulation/truth.tsv`
- `benchmark/simulation/results.tsv`
- `benchmark/competitors/runs.tsv`

**Decision: CONDITIONAL PASS**

The engine recovers planted intensity events and is conservative on
NULL. Precision is hurt by composition-bias extra segments. Competitor
and full event-class simulation grids remain open.

**Next action.** Add TMM and a composition-robust offset; expand the
event-class × width grid; run comparators when binaries are available.
