# Phase 2 / Gate 2 — Empirical FDR calibration

**Objective.** Measure null behavior of the **full** adaptive
discovery + testing path. Answer: does it exhibit acceptable empirical
null behavior?

**Inputs.** `simulate_null_counts()` over

```text
replicates = 2, 3, 5
dispersion = 0.05, 0.15, 0.40   # low, medium, high
depth      = 8, 20, 50          # low, medium, high
```

180 bins, `library_size` offsets, Gaussian PELT, selected penalty 20,
edgeR QL + BH.

**Methods.** `run_simulated_discovery()` then `empirical_fdr()` at
nominal \(q \in \{0.01, 0.05, 0.10\}\). On a true-null simulation every
discovery is a false discovery.

**Parameters.** Seed offset `10 + grid_row`. Config:
`default_config()` with the overrides above.

**Tests.** Unit: `empirical_fdr()` treats null discoveries as FDR = 1.
Grid: `benchmark/simulation/null_fdr_grid.tsv`.

**Results.**

| Question | Answer |
|----------|--------|
| Null runs with any discovery at \(q=0.05\) | **0 / 27** |
| Mean \(n\) discoveries at \(q=0.01/0.05/0.10\) | **0 / 0 / 0** |
| Failed runs | **0** |

A separate 200-bin null with selected penalty 8 produced **1 candidate
region and 0 significant tests** at \(q=0.05\). A matched GAIN
simulation produced **3 significant regions** (the planted interval plus
two composition-effect flanks).

**Does the complete adaptive procedure exhibit acceptable empirical
null behavior?**

**Yes, on this short-chromosome NB grid it is conservative** (no
null discoveries at the nominal thresholds). It is **not
anti-conservative** in these 27 cells.

This is **not** a claim that BH after segmentation is a guaranteed
region-level 5% FDR on real genomes. Longer chromosomes, weaker
penalties, TMM under global shifts, and real-data label permutations
remain open.

**Failures.** None in the grid. Label permutation on real BAM data was
**not** run (no real cohort in the checkout).

**Limitations.**

- 180–200 bins is not genome-wide. False-boundary rate on chr1-scale
  series is unmeasured.
- Conservativeness may indicate under-segmentation at \(\lambda=20\)
  rather than perfect type-I control.
- Selection-bias caveat remains in every `FDR` column
  (`fdr_interpretation`).

**Files created.** `benchmark/simulation/null_fdr_grid.tsv`

**Decision: CONDITIONAL PASS**

Acceptable null behavior on the implemented simulator grid. Do not
advance to “FDR guaranteed” language. Real-data permutation and
genome-scale nulls are still required before any manuscript calibration
claim.

**Next action.** Simulation power / boundary metrics (Phase 9) and
honest reporting of extra composition-bias segments.
