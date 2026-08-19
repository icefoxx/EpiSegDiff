# Final acceptance report — EpiSegDiff v0.2 prototype

**Objective.** State honestly what the taskbook’s final gate can and
cannot claim after this implementation pass.

**Inputs.** Specification v0.2, taskbook v0.2, this repository, local
R 4.3.2 + edgeR 4.0.16, `testthat` suite (all passing).

**Methods.** Design lock first; then count-assay MVP; Gaussian/Huber
PELT; recount + edgeR QL; empirical null grid; topology predicates;
ATAC and methylation adapters. No custom NB-FPOP.

**Parameters.** All scientific choices live in `config/example.yaml`
and `default_config()`.

## Core

| Criterion | Status |
|-----------|--------|
| Replicate-aware differential-effect estimation | **PASS** (synthetic) |
| Differential segmentation on synthetic truth | **PASS** |
| Narrow and broad regimes recovered | **PARTIAL** (toy widths only) |
| Region-level testing | **PASS** (edgeR QL) |
| Full run reproducible under fixed seed | **PASS** (in-memory pipeline test) |

## Statistics

| Criterion | Status |
|-----------|--------|
| Null behavior measured | **PASS** (27-cell grid + unit tests) |
| Empirical FDR reported | **PASS** (conservative on the grid) |
| Selection-bias limitation documented | **PASS** |
| No unsupported FDR guarantee | **PASS** |

## Boundary dynamics

| Criterion | Status |
|-----------|--------|
| Boundary MAE reported | **PARTIAL** (matched MAE 0 on bin-aligned GAIN/LOSS; Hausdorff inflated by extra calls) |
| Expansion/contraction/shift classified | **PASS** (100/class unit geometry) |
| Split/merge tested | **PASS** (unit geometry) |
| Complex not forced | **PASS** (many-to-many → COMPLEX) |

## Benchmark

| Criterion | Status |
|-----------|--------|
| csaw | **NOT RUN** (logged) |
| MACS3+DiffBind | **NOT RUN** (logged) |
| THOR | **NOT RUN** (logged) |
| Real narrow dataset | **NOT RUN** (no data bundled) |
| Real broad dataset | **NOT RUN** (no data bundled) |

## Methylation

| Criterion | Status |
|-----------|--------|
| Separate observation model | **PASS** |
| ≥2 DMR comparators run | **FAIL / not run** (logged) |
| Not analyzed with count NB | **PASS** |

## Engineering

| Criterion | Status |
|-----------|--------|
| Tests pass | **PASS** (`testthat::test_local(".")`) |
| Configuration saved | **PASS** (`config.resolved.yaml`, `run_manifest.tsv`) |
| Version manifest saved | **PASS** |
| No hidden manual steps | **PASS** |
| Custom NB-FPOP absent | **PASS** (blocked; design note only) |

## Failures

1. No real epigenomic BAM/WGBS in the workspace — Phases 11 and 13
   real-data sections cannot be closed.
2. Competitor software was not executed; rows remain in the log rather
   than being dropped.
3. Library-size composition bias produces extra significant flanks
   around strong GAIN/LOSS. Documented, not swept aside.
4. Genome-wide runtime/memory unmeasured.

## Limitations

See `docs/known_limitations.md`. v0.2 is a research prototype.

## Files created

Package tree under `/Users/qzhang10/workstation/development/EpiSegDiff`
including `R/`, `src/` (Gaussian PELT / Huber / intervals only),
`tests/`, `docs/`, `reports/`, `config/`, `benchmark/`.

Local test command:

```bash
Rscript -e 'testthat::test_local(".")'
```

**Decision: CONDITIONAL PASS** of the **software prototype and design
lock**. The taskbook’s full publication-level acceptance gate is **not**
met (real data + competitors).

**Next action.**

1. Point `config/example.yaml` at a real paired-end CUT&RUN/ChIP cohort.
2. Run `scripts/run_null_fdr_grid.R` at chromosome scale.
3. Attempt csaw / DiffBind / THOR and keep failed logs.
4. Only then revisit `docs/NB_FPOP_Experimental_Design.md` §7.
