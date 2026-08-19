# Phase 4 — Methylation adapter

**Objective.** Provide a methylation path that does **not** use the
count-assay NB / edgeR likelihood, then reuse the common segmentation
API on a delta-beta (or standardized) score.

**Inputs.** Tabular `methylated_count` / `total_count` test tables.
No WGBS cohort and no DSS/dmrseq/metilene binaries in this checkout.

**Methods.**

- `read_methylation_table()` preserves M, N, beta, coverage.
- Sites with `total_count < min_coverage` are invalid, not imputed.
- Missing CpGs are absent rows, never beta = 0.
- `estimate_methylation_effect()` is a method-of-moments delta-beta
  with a replicate-inflated SE.
- `test_dmrs()` is an overdispersed two-group test on region M/N
  (`beta_binomial` / `dss_style` label). Config validation **rejects**
  `testing.method = edgeR_ql` when `assay.type = methylation`.

**Parameters.** `methylation.min_coverage = 5`,
`methylation.score = standardized_effect`.

**Tests.** `tests/testthat/test-boundary-methylation-pipeline.R`,
`test-config.R` (methylation × edgeR forbidden).

**Results.**

- Low-coverage sites are invalid.
- Missing cytosines stay missing.
- A 3-site × 4-sample toy table yields finite delta-beta.
- Pipeline dispatch uses `.run_methylation_pipeline()` and never
  `estimate_differential_effect()` / edgeR on M counts.

**Failures.** DSS, dmrseq, and metilene were **not run**. Recorded in
`benchmark/competitors/runs.tsv` as `not-attempted-in-this-checkout`.
No real WGBS/RRBS dataset.

**Limitations.** The beta-binomial tester is a moment approximation, not
DSS. ONT adapter is out of scope. This phase is an adapter + unit tests,
not a methylation methods paper.

**Files created.** `R/methylation.R`

**Decision: CONDITIONAL PASS**

The observation-model separation required by the spec is implemented
and tested. Comparator and real-data methylation validation are open.

**Next action.** Import a public RRBS toy set and attempt DSS/dmrseq
when those packages are installed.
