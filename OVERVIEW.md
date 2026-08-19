# EpiSegDiff v0.2 — Overall overview

**Package:** EpiSegDiff 0.4.0  
**Role:** research prototype for replicate-aware, multi-scale *differential* epigenomic segmentation and boundary dynamics  
**Assays:** CUT&RUN / ChIP-seq first; ATAC and DNA methylation as separate adapters  
**Status:** v0.2.2 boundary audit — **PASS** for one narrow Level I real-data pilot  
**Date of this write-up:** 2026-08-17

This file is the single entry point. Detailed science lives in `docs/`. Gate decisions live in `reports/`. Runnable defaults live in `config/example.yaml`.

---

## 1. What this project is for

Standard differential epigenomics often does:

```text
peak calling → union/consensus intervals → count matrix → DESeq2 / edgeR / DiffBind
```

That tests *predefined* intervals. Results inherit the peak caller’s boundaries.

EpiSegDiff instead does:

```text
assay-specific observations
    → replicate-aware differential-effect estimate
    → multi-scale segmentation of that effect
    → recount raw observations in candidate intervals
    → statistical testing
    → empirical null calibration
    → boundary refinement
    → geometry / topology classification
```

The scientific target is **the spatial structure of the condition effect along the genome**, not “places where raw coverage happens to change.”

A coverage boundary that is identical in both conditions is *not* a differential boundary. Shared structure belongs in the baseline \(\alpha_i\). Only structure in \(\beta_i\) is a differential candidate.

This prototype does **not** claim to be the first changepoint, DMR, or FPOP method.

---

## 2. Two result layers

### Level I — intensity

| Label | Meaning |
|-------|---------|
| `GAIN` | treatment higher than control |
| `LOSS` | treatment lower than control |
| `STABLE` | not called after testing |

### Level II — geometry / topology

| Label | Meaning |
|-------|---------|
| `EXPANSION` | treatment domain wider, still overlapping control |
| `CONTRACTION` | treatment domain narrower |
| `LEFT_SHIFT` / `RIGHT_SHIFT` | similar width, translated |
| `BIDIRECTIONAL_SHIFT` | both ends move opposite ways without a clean width change |
| `INTERNAL_GAIN` / `INTERNAL_LOSS` | intensity change inside a shared domain |
| `SPLIT` / `MERGE` | one-to-many / many-to-one domain mapping |
| `COMPLEX` | no stable mapping; **never forced** into a prettier class |

Every Level II label is a predicate on intervals and numeric thresholds. Thresholds are configuration, not visual judgement. The rule string is stored on every event.

---

## 3. Locked statistical model (count assays)

For bin \(i\), condition \(k\), replicate \(r\):

\[
Y_{ikr} \sim \mathrm{NB}(\mu_{ikr}, \phi)
\]

\[
\log \mu_{ikr} = o_{kr} + \alpha_i + \beta_i X_k
\]

| Symbol | Meaning |
|--------|---------|
| \(Y_{ikr}\) | fragment or insertion count |
| \(o_{kr}\) | sample normalization offset |
| \(\alpha_i\) | baseline genomic signal (shared) |
| \(X_k\) | 0 = control, 1 = treatment |
| \(\beta_i\) | **condition effect** (log fold-change) |
| \(\phi\) | dispersion, estimated **once on bins**, not per discovered segment |

Default segmentation score:

\[
Z_i = \frac{\hat\beta_i}{\mathrm{SE}(\hat\beta_i)}
\]

Configurable alternatives: `logFC` (\(\hat\beta_i\)) and `shrunken_effect`.

**Testing unit** is the candidate *interval* after recount of raw fragments, not the binned score. Region tests use edgeR QL. BH q-values after adaptive discovery are **nominal**. They are not claimed to be guaranteed region-level FDR.

DNA methylation does **not** use this NB model. It uses methylated/total count pairs and a beta-binomial-style score. The shared layer is segmentation + topology, not the likelihood.

A custom Negative-Binomial FPOP engine is **blocked** until this baseline defines what accuracy problem it must solve. See `docs/NB_FPOP_Experimental_Design.md`.

Full lock: `docs/statistical_model.md`, `docs/assay_models.md`, `docs/topology_definitions.md`.

---

## 4. Pipeline (what actually runs)

```text
BAM / fragment table / methylation table
        │
        ▼
binning + masks  →  contiguous valid regions (CVRs)
        │
        ▼
normalization offsets   (library_size | TMM | optional DESeq2-style)
        │
        ▼
edgeR QL per-bin fit    →  effect, SE, score track
        │
        ▼
Gaussian / Huber / winsorized PELT on each CVR
        │
        ▼
penalty path  →  optional nesting analysis  →  candidate selection
        │
        ▼
recount raw fragments in candidates
        │
        ▼
edgeR QL region test   (or methylation-specific test)
        │
        ▼
boundary refinement
        │
        ▼
condition-specific domains  →  match  →  topology class
        │
        ▼
BED / TSV / manifest / resolved config
```

High-level call:

```r
devtools::load_all(".")
cfg <- read_episegdiff_config("config/example.yaml")
# For real data, set sample BAM paths, chrom_sizes, output$directory
# result <- run_episegdiff(cfg)

# In-memory (simulations / tests):
# result <- run_episegdiff(cfg, count_matrix = sim)
```

Every run writes:

- `config.resolved.yaml`
- `run_manifest.tsv` (versions, seed, command, runtime, warnings, input IDs)
- `normalization_diagnostics.tsv`
- `effect_track.tsv`
- `differential_segments.tsv` (+ BED if enabled)
- `tested_regions.tsv`
- `boundary_events.tsv`

---

## 5. What is implemented

| Area | Implementation |
|------|----------------|
| Config / metadata / provenance | YAML schema, sample table validation, run manifest |
| BAM QC | mapped fragments, proper pairs, duplicates, chrom coverage, fragment sizes |
| Count matrix | `fragment_midpoint`, `fragment_overlap`; 0-based half-open bins; blacklist splits CVRs |
| Normalization | `library_size`, `TMM`, DESeq2-style median-of-ratios; spike-in/background **error** (not fake) |
| Differential effect | edgeR QL; SE from QL F statistic (edgeR 4.x has no `stdev.unscaled`) |
| Segmentation | PELT (C++ Gaussian), OP for tests, Huber and winsorized Gaussian |
| Multi-scale | penalty path; nesting **evaluated**, not forced into a tree |
| Candidate selection | `single_selected_penalty` or `multi-scale_consensus` |
| Testing | `test_regions_edger()` with explicit FDR caveat on every row |
| Null / FDR | simulator + 27-cell null grid |
| Boundaries | `local_gradient`, `local_likelihood`, `cross_scale` |
| Topology | deterministic matching + ordered rules |
| ATAC adapter | configurable Tn5 offsets; insertion bins; optional fragment channels |
| Methylation adapter | M/N tables, delta-beta score, DMR test; edgeR QL rejected by config |
| Simulator | planted GAIN/LOSS/geometry events, null genomes |

C++ in `src/` is **only** Gaussian PELT, Huber helpers, and interval ops. `src/experimental/` is reserved and empty of engine code.

---

## 6. Current status (honest)

Local tests: passing (`testthat::test_local(".")`, including v0.2.1 ontology tests).

```bash
Rscript -e 'testthat::test_local(".")'
```

| Gate | Decision | One-line reason |
|------|----------|-----------------|
| Phase 0 setup | **PASS** | package skeleton, config, no NB-FPOP |
| Phase 0 design lock | **PASS** | \(\beta_i\) / \(Z_i\) defined; methylation separated; topology algorithmic |
| Phase 1 differential effect | **PASS** | synthetic recovery, sign reversal, order invariance |
| Phase 2 empirical FDR | **CONDITIONAL PASS** | 0/27 null runs called at \(q=0.05\) (conservative, short chroms) |
| Phase 3 simulation | **CONDITIONAL PASS** | GAIN/LOSS recalled; extra flanks from library-size composition bias |
| Phase 4 methylation | **CONDITIONAL PASS** | adapter tested; DSS/dmrseq/metilene not run |
| Final taskbook acceptance (v0.2) | **CONDITIONAL PASS** | prototype yes; real data + competitors no |
| v0.2.1 Gate A normalization | **PASS** | TMM default; `background` for large one-sided shifts |
| v0.2.1 Gate B topology | **PASS** | non-overlapping \(\Delta_L,\Delta_R,R_W\); no BIDIRECTIONAL_SHIFT |
| v0.2.1 Gate C boundaries | **CONDITIONAL PASS** | peak-free; narrow OK; broad MAE large |
| v0.2.1 Gate D calibration | **CONDITIONAL PASS** | \(P(\mathrm{any}\mid q=0.05)\approx 0.05\) on ~200×3 nulls |
| v0.2.1 final | **CONDITIONAL PASS** | scoped narrow-mark pilot only |

**Not done (logged, not silently dropped):**

- Real TF / H3K27ac / H3K27me3 / WGBS cohorts
- csaw, MACS3+DiffBind, THOR, DSS, dmrseq, metilene
- Genome-wide runtime and memory
- Custom NB-FPOP

Composition bias is real: library-size offsets on a large GAIN make treatment libraries look “deeper,” so background can look like a weak LOSS and extra flanks can be called. That is reported in `reports/Phase3_Simulation_Benchmark_Report.md`, not hidden.

---

## 7. How to use it

### Install / load

```r
install.packages(c("yaml", "Rcpp", "testthat"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "edgeR", "limma", "GenomicRanges", "GenomicAlignments",
  "Rsamtools", "IRanges", "rtracklayer"
))
devtools::load_all(".")
```

### Minimal in-memory example

```r
library(EpiSegDiff)

ev <- example_topology_event("GAIN", offset = 80L)
sim <- simulate_count_assay(
  n_bins = 200, n_rep = 3, depth = 25, dispersion = 0.08,
  events = ev, seed = 1
)

cfg <- default_config()
cfg$samples <- sim$samples
cfg$output$directory <- "results/demo"
cfg$normalization$method <- "library_size"
cfg$segmentation$selected_penalty <- 8
cfg$segmentation$penalties <- c(8, 16)

res <- run_episegdiff(cfg, count_matrix = sim)
head(res$tested)
head(res$events)
```

### Real BAM run

Edit `config/example.yaml`:

1. `samples[].bam_or_methylation_path`
2. `preprocessing.chrom_sizes` (chrom + length TSV)
3. `preprocessing.blacklist` if you have one
4. `output$directory`
5. explicit scientific keys: `bin_size`, `normalization.method`, `differential.score`, `segmentation.loss`, `segmentation.penalties`

Then:

```r
cfg <- read_episegdiff_config("config/example.yaml")
res <- run_episegdiff(cfg)
```

### Scripts

| Script | Purpose |
|--------|---------|
| `scripts/run_null_fdr_grid.R` | null replicate × depth × dispersion grid |
| `scripts/run_simulation_benchmark.R` | small GAIN/LOSS/NULL recovery table |

---

## 8. Decisions that must stay explicit

These are never “hidden defaults.” They are always in the config and the manifest:

```text
bin size
normalization method
mask / blacklist
dispersion strategy
differential score          (logFC | standardized_effect | shrunken_effect)
segmentation loss           (gaussian | huber | winsorized_gaussian)
penalty / penalty grid
minimum region width
testing method
FDR method                  (nominal BH; not selective FDR)
boundary refinement rule
topology thresholds
candidate / merge strategy
```

Coordinates are **BED 0-based half-open** `[start, end)`.

---

## 9. Repository map

```text
EpiSegDiff/
├── OVERVIEW.md                 ← this file
├── README.md                   install / test commands
├── config/example.yaml         full explicit config
├── R/                          package code
├── src/                        Gaussian PELT, Huber, intervals
│   └── experimental/           reserved; no NB-FPOP
├── docs/
│   ├── statistical_model.md
│   ├── assay_models.md
│   ├── topology_definitions.md
│   ├── benchmark_plan.md
│   ├── known_limitations.md
│   └── NB_FPOP_Experimental_Design.md
├── reports/                    phase gates (PASS / CONDITIONAL / FAIL)
├── tests/testthat/             109 tests
├── benchmark/simulation/       null grid + toy truth/results
├── benchmark/competitors/      failed/not-run log (kept on purpose)
├── scripts/
└── vignettes/EpiSegDiff-mvp.Rmd
```

### Reports

| File | Contents |
|------|----------|
| `reports/Phase0_Project_Setup_Report.md` | repo + test harness |
| `reports/Phase0_Design_Lock_Report.md` | locked model + unresolved decisions |
| `reports/Phase1_Differential_Effect_Report.md` | effect-track validation |
| `reports/Phase2_FDR_Calibration_Report.md` | empirical null |
| `reports/Phase3_Simulation_Benchmark_Report.md` | recovery metrics |
| `reports/Phase4_Methylation_Validation_Report.md` | methylation adapter |
| `reports/Final_Acceptance_Report.md` | taskbook checklist |

---

## 10. What to do next (recommended order)

Do **not** start with C++ NB-FPOP.

1. Point the example config at one real paired-end CUT&RUN or ChIP-seq cohort (narrow mark).
2. Repeat on one broad mark (H3K27me3 / H3K9me3).
3. Re-run the null grid on chromosome-scale series; keep anti-conservative cells if they appear.
4. Attempt csaw / MACS3+DiffBind / THOR and **leave failed runs in the log**.
5. Only then fill `docs/NB_FPOP_Experimental_Design.md` §7 (what accuracy problem remains).
6. Methylation: public RRBS/WGBS toy + DSS/dmrseq if installable.

---

## 11. Positioning sentence (for talks / notes)

> EpiSegDiff is a replicate-aware, multi-scale framework for detecting differential epigenomic regions and characterizing their boundary and structural dynamics without requiring predefined peak boundaries. Assay-specific likelihoods are preserved; the shared layer is differential-effect segmentation and quantitative topology. v0.2 is a prototype whose novelty and FDR claims must be measured, not assumed.
