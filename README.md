# EpiSegDiff v0.2

Replicate-aware, multi-scale differential epigenomic segmentation and
boundary-dynamics analysis.

EpiSegDiff detects **differential epigenomic structure** without requiring
predefined peaks, then characterizes both **intensity** (gain/loss) and
**geometry/topology** (expansion, contraction, shift, split, merge).

This is a research prototype. It does not claim to be the first method of
any broad class. Novelty hypotheses must be tested empirically.

**Start here:** [OVERVIEW.md](OVERVIEW.md) — science, pipeline, status, and how to run.

## Scientific target

For count assays the observation model is Negative Binomial. Segmentation
targets the genomic **condition effect** \(\beta_i\) (or a standardized
form \(Z_i = \hat\beta_i / \mathrm{SE}(\hat\beta_i)\)), not raw coverage
changepoints shared by all conditions.

DNA methylation uses a separate beta-binomial observation model. The
shared layer is segmentation and boundary classification, not the
likelihood.

A custom Negative-Binomial FPOP engine is **not** part of the v0.2 core.
See `docs/known_limitations.md`.

## Install

```r
# From the repository root, with Bioconductor dependencies available:
install.packages(c("yaml", "Rcpp", "testthat"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "edgeR", "limma", "GenomicRanges", "GenomicAlignments",
  "Rsamtools", "IRanges", "rtracklayer", "GenomeInfoDb"
))
devtools::install(".")
```

Load without installing:

```r
devtools::load_all(".")
```

## Run

```r
library(EpiSegDiff)
cfg <- read_episegdiff_config("config/example.yaml")
# result <- run_episegdiff(cfg)
```

Every analysis parameter is configurable. Hidden scientific defaults are
not permitted; see `config/example.yaml` and `docs/statistical_model.md`.

## Tests

Local test command (documented acceptance requirement):

```bash
Rscript -e 'testthat::test_local(".")'
```

Alternatively:

```r
devtools::load_all(".")
testthat::test_dir("tests/testthat")
```

Package check (optional, requires Suggested packages):

```bash
R CMD build .
R CMD check EpiSegDiff_0.2.0.tar.gz --no-manual
```

## Pipeline (count-assay MVP)

```text
BAM / fragment counts
  → binning + masking
  → normalization offsets
  → replicate-aware differential-effect track
  → Gaussian / Huber multi-scale segmentation
  → candidate recount
  → edgeR QL region testing
  → empirical null calibration
  → boundary refinement
  → topology classification
```

## Documentation

| File | Contents |
|------|----------|
| `docs/statistical_model.md` | Count-assay model, scores, testing unit |
| `docs/assay_models.md` | CUT&RUN, ChIP-seq, ATAC, methylation |
| `docs/topology_definitions.md` | Quantitative event rules |
| `docs/benchmark_plan.md` | Simulation, comparators, metrics |
| `docs/known_limitations.md` | Explicit non-goals and caveats |
| `docs/NB_FPOP_Experimental_Design.md` | Blocked experimental branch |

## Reproducibility

Each `run_episegdiff()` writes:

- the resolved configuration
- software / package versions
- random seed
- input file hashes or stable identifiers
- command line
- runtime
- warnings and errors

BH-adjusted q-values after adaptive region discovery are **not** claimed
to be guaranteed region-level FDR. Empirical null calibration is required.

## Development stages

See `EpiSegDiff_Project_Specification_v0.2` (design) and the agent
taskbook. Phase 14 (experimental NB-FPOP) is blocked until the MVP
benchmarks define the accuracy problem a custom likelihood must solve.
