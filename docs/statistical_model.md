# Statistical model (count assays)

**Status:** Phase 0 design lock for EpiSegDiff v0.2  
**Scope:** CUT&RUN, ChIP-seq, and ATAC *count* observations  
**Not in this document:** DNA methylation (see `assay_models.md`)

This document locks the symbols, the segmentation target, dispersion
strategy, and the testing unit. Implementation must not silently change
this model.

---

## 1. Observation

For genomic bin \(i\), condition \(k\), and biological replicate \(r\):

\[
Y_{ikr} \sim \mathrm{NB}(\mu_{ikr}, \phi)
\]

| Symbol | Meaning |
|--------|---------|
| \(Y_{ikr}\) | Observed fragment (or insertion) count |
| \(\mu_{ikr}\) | Mean count |
| \(\phi\) | Dispersion (edgeR parameterization; variance \(\mu + \phi\mu^2\)) |
| \(i\) | Genomic bin index within a chromosome / contiguous valid region |
| \(k\) | Biological condition (`control`, `treatment` in the two-group MVP) |
| \(r\) | Biological replicate index within condition |

The Negative Binomial is an **observation model for count assays**. It is
not used for methylation proportions.

### What is counted

Configured by `preprocessing.counting_mode`:

| Mode | Observation |
|------|-------------|
| `fragment_midpoint` | Midpoint of a properly paired fragment falls in the bin |
| `fragment_overlap` | Fragment interval overlaps the bin (one count per overlapping fragment) |
| `insertion` | ATAC Tn5 insertion after optional offset (ATAC adapter only) |

Fragment-length channels are optional and assay-aware. They are never
forced from ATAC onto CUT&RUN/ChIP-seq.

---

## 2. Mean model

\[
\log \mu_{ikr} = o_{kr} + \alpha_i + \beta_i X_k
\]

| Symbol | Meaning |
|--------|---------|
| \(o_{kr}\) | Sample-level normalization offset (log scale) |
| \(\alpha_i\) | Baseline genomic signal (shared across conditions) |
| \(X_k\) | Design variable. Two-group MVP: \(X_{\mathrm{control}}=0\), \(X_{\mathrm{treatment}}=1\) |
| \(\beta_i\) | Genomic **condition effect** (log fold-change of treatment vs control) |

### Replicate handling

Replicates share \(\alpha_i\) and \(\beta_i\). Residual replicate variation
is absorbed by the NB (or QL) variance. The design matrix is

```text
~ condition
```

with `control` as the reference level unless the user sets another
`differential.control_level`. Sample order must not change estimates
beyond floating-point noise; this is a required test.

The two-group MVP does not include batch covariates. Additional columns
in `differential.design` are allowed only when the user supplies them
explicitly.

---

## 3. Normalization offset

Normalization is **separated from segmentation**.

\[
o_{kr} = \log s_{kr}
\]

where \(s_{kr}\) is a size factor:

| `normalization.method` | \(s_{kr}\) |
|------------------------|------------|
| `library_size` | Library size / geometric mean of library sizes |
| `TMM` | Library size \(\times\) TMM norm.factor (edgeR) |
| `deseq2` | DESeq2-style median-of-ratios (optional; not required for MVP) |
| `spikein` | Spike-in derived factor (only if `spikein_usage: true`) |
| `background` | Experimental; requires `background_definition` |

TMM / RPGC are **not** assumed valid under global biological shifts.
Every run records `normalization_method`, `size_factors`,
`spikein_usage`, `background_definition`, and `excluded_regions`.

Offsets are stored as a sample vector (or sample \(\times\) bin matrix if
a future method requires it). Visualization tracks (CPM/RPGC) are
derived separately and are **not** used as testing counts.

---

## 4. Dispersion

Dispersion is estimated **once, before segmentation**, on the binned
count matrix after optional user-specified filters.

Locked strategy for v0.2 (`differential.dispersion: edgeR_estimateDisp`):

1. Build an edgeR `DGEList` from raw integer counts.
2. Attach normalization offsets / norm.factors.
3. Estimate common, trended, and tagwise dispersion with
   `edgeR::estimateDisp` under the condition design.
4. Use quasi-likelihood (`glmQLFit`) for region-level testing.

This avoids the circularity

```text
need segments to estimate dispersion
and
need dispersion to identify segments
```

**There is no segment-specific dispersion in the MVP.** A custom NB
changepoint cost that estimates \(\phi\) per candidate segment is an
experimental branch (Phase 14) and is blocked until this baseline is
validated.

---

## 5. Differential-effect series (segmentation input)

After the local / pooled NB (edgeR) fit, each valid bin yields:

```text
beta_hat[i]     # \hat\beta_i  (logFC)
SE_beta[i]      # estimated standard error
Z[i]            # beta_hat[i] / SE_beta[i]
```

The **segmentation target is the condition effect**, not raw coverage.

Configurable score (`differential.score`):

| Name | Definition |
|------|------------|
| `logFC` | \(D_i = \hat\beta_i\) |
| `standardized_effect` | \(Z_i = \hat\beta_i / \mathrm{SE}(\hat\beta_i)\) **(default)** |
| `shrunken_effect` | variance-weighted shrink of \(\hat\beta_i\) toward 0: \(\tilde\beta_i = \hat\beta_i \cdot \tau^2 / (\tau^2 + \mathrm{SE}_i^2)\), with \(\tau^2\) estimated from the empirical distribution of \(\hat\beta\) on valid bins |

Default is `standardized_effect` because it down-weights high-variance
bins (low coverage) without using a segment-specific dispersion.

The Gaussian / Huber segmentation engines treat the score as a 1-D
numeric series with approximately constant residual scale after the
user's choice of score. The implicit Gaussian variance is 1 unless the
user pre-scales the series. Penalty values are in **the same units as
the chosen cost** (RSS for Gaussian).

---

## 6. Zero-count and missing-bin behavior

| Situation | Matrix | `valid` | Segmentation |
|-----------|--------|---------|--------------|
| All samples zero in a bin that is not masked | count 0 retained | `TRUE` unless `mask_low_coverage` | Effect may be NA; NA scores are excluded from the CVR (treated as a gap) |
| Blacklist / user exclude / assembly gap overlap | count may be computed or set 0 | `FALSE` | Splits contiguous valid regions (CVRs) |
| Chromosome absent in a BAM | 0 for that sample | still `TRUE` if not masked | May inflate apparent loss; warned in QC |
| `mask_low_coverage: true` and count-sum \(<\) `min_count_sum` | retained | `FALSE` | Splits CVRs. This **can** create artificial boundaries; it is opt-in |

**No interpolation across invalid bins.** BigWig export must not fill
masked bins.

Contiguous valid regions (CVRs) are the independent segmentation units
within each chromosome.

---

## 7. Segmentation (MVP)

Within each CVR, minimize a penalized change-in-mean cost on the
differential score \(D_i\):

\[
\sum_{j=1}^{m} C(D_{s_j:e_j}) + \lambda \, m
\]

| Loss | Segment cost \(C\) |
|------|--------------------|
| `gaussian` | residual sum of squares about the segment mean |
| `huber` | Huber loss about the Huber M-estimator of the mean |
| `winsorized_gaussian` | RSS after MAD-winsorization of the CVR |

Engines: PELT (default), optional FPOP (Gaussian), and exact optimal
partitioning (`op`) for small \(n\) and unit tests.

A penalty path \(\lambda_1 < \cdots < \lambda_K\) produces candidate
segments at multiple scales. These candidates are **not** assumed to
form a tree; nesting is evaluated separately (`analyze_nesting`).

A custom NB-FPOP / vector-valued / beta-binomial changepoint cost is
**out of scope** for the MVP.

---

## 8. Testing unit

Segmentation proposes candidate intervals. It is **not** the confirmatory
test.

1. Recount **raw fragments / insertions** in each candidate interval
   (not the binned score, not CPM).
2. Fit a replicate-aware edgeR QL model on the region \(\times\) sample
   count matrix.
3. Report `logFC`, `logCPM`, QL F statistic, nominal p-value, and
   BH-adjusted `FDR`.

The testing unit is the **candidate interval**, after selection.

**Critical limitation:** BH adjustment after adaptive region discovery is
**not** guaranteed region-level FDR. Empirical null experiments (mock vs
mock, label permutation, simulated null genomes) are required. See
`docs/benchmark_plan.md`.

---

## 9. What must not be segmented

Do not treat a coverage boundary that is **shared identically by all
conditions** as a differential boundary. Shared baseline structure lives
in \(\alpha_i\). Only structure in \(\beta_i\) (or \(Z_i\)) is a
differential candidate.

---

## 10. Locked decisions (v0.2 MVP)

1. Two biological conditions.
2. Segmentation input = configured differential score, default \(Z_i\).
3. Dispersion = edgeR `estimateDisp` on bins, not on discovered segments.
4. Coordinates = BED 0-based half-open \([start, end)\).
5. Last incomplete bin of a chromosome is **kept** if length \(\ge 1\).
6. Methylation never enters this NB model.
