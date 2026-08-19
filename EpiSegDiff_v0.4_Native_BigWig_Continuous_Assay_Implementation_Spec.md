# EpiSegDiff v0.4 — Native BigWig / Continuous-Assay Implementation Specification

**Status:** Agent-executable development specification  
**Version:** v0.4  
**Primary objective:** Promote BigWig from an exploratory proxy to a first-class continuous-signal input while preserving a statistically distinct count-assay path.  
**Reference pilot:** GSE280947 v0.3 BigWig proxy addendum.

---

## 1. Motivation and evidence from v0.3

The v0.3 GSE280947 pilot established the following facts:

- No BAM files were available in the checkout.
- Eight deposited spike-in+CPM BigWigs could be imported with `rtracklayer` and analyzed on a common mm10 grid.
- Replicate agreement was high: DMSO JUN r=0.887 and JNKi JUN r=0.871.
- The continuous-signal pilot produced 630 LOSS and 1 GAIN regions at FDR <= 0.05.
- Mean JUN signal decreased from 0.947 (DMSO) to 0.639 (JNKi), consistent with the deposited IDR peak contraction (4960 to 1632).
- The pilot used limma on `log2(score + 0.01)` region means. It did **not** use fragment counting or the NB/edgeR model.

**Decision for v0.4:** BigWig becomes a native continuous-signal assay input. This does **not** establish statistical equivalence between BigWig and BAM/count assays.

---

## 2. Core design principle

EpiSegDiff must separate **input representation**, **observation model**, and **segmentation**.

```text
                         EpiSegDiff
                             |
                    Assay abstraction
                    /               \
             count assay       continuous assay
                 |                   |
                BAM                BigWig
                 |                   |
        fragment counts      normalized signal
                 |                   |
          count model       continuous model
           (NB/edgeR)       (v0.4: limma path)
                    \               /
                     standardized effects
                             |
                        segmentation
                             |
                    differential regions
                       GAIN / LOSS
```

Do **not** convert BigWig values to pseudo-counts merely to run edgeR.

DESeq2 is **not** a required v0.4 dependency or implementation target.

---

## 3. Scope

### 3.1 Required in v0.4

1. Native BigWig input.
2. Explicit assay abstraction: `count` versus `continuous`.
3. BigWig metadata/provenance validation.
4. Common genomic grid construction.
5. Replicate QC.
6. Continuous-signal differential modeling.
7. Configurable signal filtering.
8. Configurable transformation.
9. Segmentation and GAIN/LOSS calling.
10. Standardized output schema shared as far as possible with the count path.
11. Parameter-sensitivity analysis.
12. Explicit statistical capability levels depending on replication.
13. Verification of the v0.3 width/min-width inconsistency.
14. Machine-readable run manifest plus human-readable report.
15. Regression tests preserving the v0.3 GSE280947 result direction.

### 3.2 Not required for v0.4 acceptance

- Proof that BigWig and BAM produce equivalent DER sets.
- Replacing the existing count/NB workflow.
- DESeq2 integration.
- Competitor benchmarking based only on the current BigWig pilot.
- Universal normalization of arbitrary BigWigs.
- DNA methylation-specific statistical models.

These can be addressed in later versions.

---

## 4. Input contract

### 4.1 Supported assay families

```yaml
assay_family: count | continuous
input_type: bam | bigwig
```

For v0.4:

- `count + bam`: preserve existing count workflow.
- `continuous + bigwig`: new native workflow.
- Reject incompatible combinations unless explicitly implemented later.

### 4.2 BigWig sample sheet

Minimum required fields:

```tsv
sample	condition	replicate	file
DMSO_R1	DMSO	1	path/DMSO_R1.bw
DMSO_R2	DMSO	2	path/DMSO_R2.bw
JNKi_R1	JNKi	1	path/JNKi_R1.bw
JNKi_R2	JNKi	2	path/JNKi_R2.bw
```

Optional fields:

```text
batch
control_sample
assay
signal_type
normalization
source
GEO_accession
```

### 4.3 Required continuous-signal metadata

At minimum record:

```yaml
assay_family: continuous
input_type: bigwig
genome: mm10
signal_type: CPM
normalization: spike_in_plus_CPM
```

Allowed `signal_type` vocabulary should initially include:

```text
raw_coverage
CPM
RPKM
RPGC
fold_enrichment
log2_ratio
spikein_normalized
other
unknown
```

`unknown` must generate a warning and be recorded in output provenance. Do not silently infer the signal semantics from filename alone.

### 4.4 File validation

For every BigWig:

- file exists and is readable;
- BigWig header is valid;
- chromosome naming convention is identified;
- requested genome/chromosome set is compatible;
- file is non-empty;
- finite-value fraction is reported;
- signal range and selected quantiles are reported;
- chromosome coverage is reported.

Fail early for incompatible genome/chromosome mappings unless the user explicitly supplies a mapping.

---

## 5. Statistical capability levels

### Level I — replicated continuous inference

Requirements:

- >=2 biological replicates per tested condition (minimum supported);
- valid design matrix;
- estimable contrast.

Capabilities:

- effect size;
- statistical testing;
- multiple-testing correction;
- GAIN/LOSS classification;
- segmentation;
- FDR-bearing output.

The report must still warn when replicate number is very small.

### Level II — unreplicated/descriptive continuous comparison

Example:

```text
Control.bw
Treatment.bw
```

Capabilities:

- effect map;
- segmentation;
- candidate GAIN/LOSS regions;
- descriptive ranking.

Prohibited:

- fabricated replicate variance;
- ordinary replicate-based FDR claims;
- labeling descriptive candidates as statistically significant DERs.

Required output field:

```text
inference_mode = descriptive
```

### Level III — control-aware continuous analysis

Examples:

```text
TF + IgG
ChIP + Input
```

Support control information for masking/background assessment. Do not assume subtraction or ratio is universally valid; make the selected operation explicit in configuration and provenance.

---

## 6. BigWig ingestion and common grid

### 6.1 Reader

Implement a reusable BigWig ingestion module rather than embedding `rtracklayer` calls directly in the analysis script.

Suggested interface:

```r
read_bigwig_assay(
    sample_sheet,
    genome,
    chromosomes = NULL,
    bin_size = 200,
    summary_fun = "mean"
)
```

### 6.2 Common genomic grid

All samples in a contrast must be summarized over exactly the same genomic intervals.

Default candidate:

```text
bin_size = 200 bp
```

But 200 bp is **not yet validated as a universal default**. It must remain configurable.

At minimum test:

```text
100 bp
200 bp
500 bp
```

50 bp may be supported as an opt-in high-resolution mode but is not required for the default genome-wide acceptance test.

### 6.3 Summary statistic

Default candidate:

```text
mean signal per bin
```

Architecture must allow later addition of max, integral/area, median, or other summaries without rewriting the assay layer.

---

## 7. QC requirements

Produce both machine-readable QC tables and plots.

### Required metrics

Per sample:

- min/max;
- mean/median;
- zero fraction;
- non-finite fraction;
- quantiles;
- chromosome coverage.

Across replicates:

- Pearson correlation;
- Spearman correlation;
- sample-sample correlation matrix;
- optional PCA/MDS when enough samples exist.

### Required plots

1. Signal distribution plot.
2. Replicate scatterplot or hexbin equivalent.
3. Correlation heatmap.
4. Condition-level signal comparison.

Do not hard-code a universal replicate-correlation pass threshold in v0.4. Report the values and flag obviously problematic samples using configurable QC policy.

---

## 8. Signal filtering

The v0.3 threshold `max JUN score >= 0.6` was dataset-specific and must **not** become the universal default.

Implement:

```yaml
signal_filter:
  method: none | absolute | quantile | control_based | auto
```

### Absolute

```yaml
method: absolute
threshold: 0.6
```

Only when explicitly requested.

### Quantile

Example:

```yaml
method: quantile
quantile: 0.75
```

### Control-based

Use control-track distribution to establish a signal/background criterion. The exact rule must be reported.

### Auto

If implemented, the rule must be deterministic, documented, and written to the run manifest. No opaque heuristic.

The default should be conservative. If `auto` is not sufficiently validated, use `none` as the software default and require explicit filtering in analysis profiles.

---

## 9. Transformation

v0.3 used:

```r
log2(score + 0.01)
```

This is acceptable as a pilot but the pseudocount is not yet validated.

Implement configurable transformation:

```yaml
transform:
  method: log2_pseudocount | log1p | none
  pseudocount: 0.01
```

### Required sensitivity test

For GSE280947 compare at least:

```text
0.001
0.01
0.05
0.1
```

Report:

- number of tested regions;
- number GAIN/LOSS at selected FDR;
- effect-size correlation;
- DER overlap/Jaccard;
- sign concordance;
- rank correlation.

Do not select a pseudocount solely because it maximizes the number of significant regions.

---

## 10. Continuous differential model

For v0.4, retain limma as the reference implementation for replicated continuous signal.

Conceptually:

```text
BigWig
  -> common bins/regions
  -> transformed continuous matrix
  -> design matrix
  -> contrast
  -> limma
  -> effect + P + FDR
```

Required outputs per tested unit:

```text
region/bin identifier
chrom
start
end
mean_condition_A
mean_condition_B
log2FC or clearly defined effect
SE/statistic where available
p_value
FDR
contrast
inference_mode
```

Do not call edgeR on normalized BigWig values.

Do not add DESeq2 merely for symmetry with the count path.

---

## 11. Segmentation

### 11.1 Current reference

v0.3 used:

```text
Gaussian PELT
penalty = 24
min width = 3 bins
bin size = 200 bp
```

### 11.2 Mandatory v0.4 blocker investigation

The v0.3 report states:

```text
bin size = 200 bp
min width = 3 bins = nominally 600 bp
median reported LOSS width = 400 bp
```

This apparent inconsistency must be explained before v0.4 acceptance.

Agent must trace:

1. definition of `min.width` in the segmentation implementation;
2. whether width is measured before or after masking;
3. whether post-segmentation trimming occurs;
4. whether intervals are merged/split after segmentation;
5. BED coordinate conversion;
6. chromosome-edge behavior;
7. whether the reported width corresponds to segment width or final tested DER width.

Create a minimal synthetic regression test where the expected number of bins and genomic width are known exactly.

**Acceptance requirement:** observed output widths must be mathematically consistent with the documented algorithm, or the documentation/implementation must be corrected.

### 11.3 Parameterization

Expose at least:

```yaml
segmentation:
  method: gaussian_pelt
  penalty: 24
  min_bins: 3
```

Do not claim penalty 24 is universally optimal.

---

## 12. GAIN / LOSS definition

Direction must be tied explicitly to the contrast.

For:

```text
contrast = JNKi - DMSO
```

then:

```text
GAIN = effect > 0
LOSS = effect < 0
```

Every result table/report must print the contrast definition.

Never emit GAIN/LOSS without preserving the condition ordering.

---

## 13. Standardized result schema

Primary result file:

```text
results/episegdiff_primary.tsv
```

Minimum columns:

```text
chrom
start
end
region_id
assay_family
input_type
contrast
effect
log2FC
p_value
FDR
direction
mean_group1
mean_group2
n_group1
n_group2
inference_mode
signal_filter_method
transform_method
segmentation_method
```

If a field is not statistically meaningful in descriptive mode, write `NA`; do not manufacture values.

---

## 14. Provenance and reproducibility

Every run must write:

```text
run_manifest.yaml
session_info.txt
sample_qc.tsv
analysis_parameters.yaml
```

Manifest should include:

- software version/git commit;
- date;
- input files;
- checksums where practical;
- genome;
- signal type;
- normalization metadata;
- bin size;
- transformation;
- pseudocount;
- filter method and threshold;
- model;
- design/contrast;
- segmentation parameters;
- FDR threshold;
- inference mode.

A result should be reproducible from the manifest plus input files.

---

## 15. GSE280947 regression dataset

Use GSE280947 as the primary v0.4 continuous-assay regression dataset.

The objective is **not** to hard-code exact region counts. Instead verify stable qualitative and quantitative behavior.

### Required checks

1. All eight deposited BigWigs import successfully.
2. DMSO and JNKi JUN replicate correlations remain high and close to the v0.3 values (0.887 and 0.871), allowing small implementation-dependent differences.
3. Mean JUN signal remains lower in JNKi than DMSO at the reference signal set.
4. LOSS strongly dominates GAIN under the reference configuration.
5. Direction convention is correct for `JNKi - DMSO`.
6. The pipeline produces deterministic output for a fixed environment/configuration.
7. Width/min-width behavior passes the new regression test.

Do not require exactly 630 LOSS / 1 GAIN if bug fixes or explicitly documented parameter changes legitimately alter boundaries. Any substantial change must be explained quantitatively.

---

## 16. Sensitivity analysis

Run a compact parameter grid on GSE280947.

### Bin size

```text
100, 200, 500 bp
```

### Pseudocount

```text
0.001, 0.01, 0.05, 0.1
```

### Signal filtering

At minimum:

```text
none
reference absolute threshold (0.6; GSE280947 only)
one data-adaptive/control-aware option if implemented
```

### Report

For each configuration record:

- runtime;
- memory if readily measurable;
- number of tested units;
- number of GAIN/LOSS;
- median/p90 width;
- effect correlation versus reference;
- sign concordance;
- overlap/Jaccard of called DERs.

The goal is robustness characterization, not parameter optimization against an assumed truth set.

---

## 17. IDR comparison: interpretation rules

The v0.3 classes were:

```text
Class A: shared between proxy DER and deposited IDR union = 257
Class B: IDR-only = 5144
Class C: proxy-only = 374
```

Retain this comparison only as contextual validation.

Do **not** treat deposited IDR peaks as differential ground truth.

Do **not** report `A/(A+B)` as sensitivity of the differential method.

IDR indicates reproducible peak presence; differential testing addresses a different question.

---

## 18. BAM versus BigWig relationship

v0.4 must use the following terminology:

### Count assay

```text
BAM -> fragment counts -> NB/count model
```

### Continuous assay

```text
BigWig -> continuous signal -> continuous model
```

BigWig is a **first-class input**, not globally a "proxy".

Use the term `coverage proxy` only when specifically asking whether BigWig analysis reproduces a fragment-count/BAM result.

No claim of BAM/BigWig equivalence is permitted until tested on a dataset where both representations are available.

---

## 19. Software architecture

Recommended modules:

```text
R/
  assay_base.R
  assay_count.R
  assay_bigwig.R
  bigwig_io.R
  genomic_grid.R
  qc_continuous.R
  transform_continuous.R
  filter_signal.R
  model_continuous_limma.R
  segmentation.R
  result_schema.R
  provenance.R
  reporting.R
```

Existing repository conventions take precedence if they differ; avoid unnecessary restructuring.

### Assay object concept

A continuous assay should carry at least:

```text
sample metadata
genomic intervals
signal matrix
genome
signal type
normalization provenance
transformation state
filter state
```

Downstream segmentation should consume a standardized statistical/effect object rather than directly depending on BigWig file handling.

---

## 20. CLI/config behavior

Target usage should permit a configuration-driven run such as:

```bash
Rscript scripts/run_episegdiff.R \
  --config configs/gse280947_bigwig.yaml
```

Example config:

```yaml
version: "0.4"

assay:
  assay_family: continuous
  input_type: bigwig
  genome: mm10
  signal_type: CPM
  normalization: spike_in_plus_CPM

samples: configs/gse280947_samples.tsv

contrast:
  numerator: JNKi
  denominator: DMSO

binning:
  bin_size: 200
  summary: mean

signal_filter:
  method: absolute
  threshold: 0.6

transform:
  method: log2_pseudocount
  pseudocount: 0.01

model:
  method: limma

segmentation:
  method: gaussian_pelt
  penalty: 24
  min_bins: 3

calling:
  fdr: 0.05

output:
  directory: results/gse280947_v04
```

The GSE280947 threshold of 0.6 is a dataset profile value, **not** a package-wide default.

---

## 21. Error handling

Pipeline must fail with informative messages for:

- missing BigWig;
- malformed sample sheet;
- duplicated sample IDs;
- condition absent from sample metadata;
- contrast not estimable;
- incompatible chromosome/genome naming;
- all-zero or effectively empty signal;
- non-finite transformation results;
- insufficient replicates for requested inferential mode.

For an unreplicated comparison, the pipeline may automatically downgrade to descriptive mode only if this behavior is explicit in config/report; otherwise fail and ask for descriptive mode.

---

## 22. Tests

### Unit tests

At minimum:

1. BigWig metadata parsing.
2. Common-grid coordinate correctness.
3. Mean signal summarization.
4. Transformation/pseudocount behavior.
5. Directionality of contrast.
6. GAIN/LOSS classification.
7. FDR column handling.
8. Unreplicated-mode safeguards.
9. Width/min-bin consistency.
10. BED coordinate output.
11. Result-schema validation.
12. Deterministic configuration parsing.

### Synthetic integration test

Create synthetic continuous tracks containing known:

- unchanged background;
- one GAIN block;
- one LOSS block;
- known segment boundaries;
- replicated noise.

Test whether the pipeline recovers direction and approximately correct boundaries without requiring unrealistic exact p-values.

### Real-data regression test

GSE280947 as described above.

---

## 23. Required report

Generate:

```text
reports/v0.4_Native_BigWig_Continuous_Assay_Acceptance_Report.md
```

Required sections:

1. Executive summary.
2. Implementation completed.
3. Deviations from specification.
4. Input validation.
5. GSE280947 QC.
6. Statistical model and contrast.
7. Segmentation results.
8. Width/min-width investigation.
9. Sensitivity analysis.
10. Regression-test results.
11. Failures/warnings.
12. Known limitations.
13. Acceptance matrix.
14. Final recommendation: PASS / CONDITIONAL PASS / FAIL.
15. Exact next actions.

Use explicit labels where appropriate:

```text
FACT
INFERENCE
UNRESOLVED
```

Do not convert an inference into a fact.

---

## 24. Acceptance criteria

### Critical — all required for PASS

- [ ] BigWig works as a native `continuous` assay input.
- [ ] BigWig is no longer globally labeled as an exploratory proxy.
- [ ] Count and continuous assay paths are explicitly separated.
- [ ] Normalized BigWig signal is not passed to edgeR as pseudo-count data.
- [ ] Sample and signal provenance is retained.
- [ ] Replicated continuous analysis produces valid effect/P/FDR output.
- [ ] Unreplicated data cannot produce unsupported replicate-based FDR claims.
- [ ] Contrast direction is explicit and tested.
- [ ] Width/min-width inconsistency from v0.3 is resolved and regression-tested.
- [ ] GSE280947 regression analysis completes.
- [ ] Reproducibility manifest is generated.
- [ ] Core unit/integration tests pass.

### Important — required unless technically justified

- [ ] Configurable bin size.
- [ ] Configurable transformation/pseudocount.
- [ ] Configurable signal filtering.
- [ ] Sensitivity analysis completed.
- [ ] QC tables and plots generated.
- [ ] Standardized result schema implemented.

### Explicitly not required

- [ ] Exact reproduction of 630 LOSS / 1 GAIN after justified bug fixes.
- [ ] DESeq2 support.
- [ ] BigWig/BAM equivalence claim.
- [ ] Competitor benchmark.

---

## 25. Agent execution order

Execute in this order. Do not skip the blocker investigation.

### Phase 0 — repository audit

- Inspect current v0.3 code and tests.
- Identify count-specific assumptions.
- Locate `run_gse280947_bigwig_proxy.R` and all generated outputs.
- Record baseline behavior before modification.

**Deliverable:** `reports/v0.4_phase0_repository_audit.md`

### Phase 1 — investigate width inconsistency

- Reproduce the 200 bp / min 3 bins / median 400 bp observation.
- Trace coordinates through segmentation and post-processing.
- Build synthetic width regression test.
- Fix implementation or documentation.

**Gate:** Do not proceed to final acceptance with this unresolved.

**Deliverable:** `reports/v0.4_phase1_width_investigation.md`

### Phase 2 — assay abstraction + BigWig ingestion

- Implement assay family metadata.
- Implement native BigWig reader.
- Implement common grid.
- Implement input validation and provenance.

**Deliverable:** `reports/v0.4_phase2_bigwig_ingestion.md`

### Phase 3 — continuous statistical path

- Implement configurable transformation.
- Implement filtering.
- Implement replicated limma path.
- Implement descriptive unreplicated safeguards.
- Standardize result output.

**Deliverable:** `reports/v0.4_phase3_continuous_model.md`

### Phase 4 — segmentation integration

- Connect standardized continuous statistics/effects to segmentation.
- Preserve explicit contrast direction.
- Validate GAIN/LOSS.

**Deliverable:** `reports/v0.4_phase4_segmentation.md`

### Phase 5 — QC + reporting

- Generate QC metrics and plots.
- Generate manifests/session information.
- Ensure reports identify assay family, signal type, normalization, and inference mode.

**Deliverable:** `reports/v0.4_phase5_qc_reporting.md`

### Phase 6 — sensitivity and regression

- Run GSE280947 reference configuration.
- Run parameter sensitivity grid.
- Run synthetic integration test.
- Compare with v0.3 baseline.

**Deliverable:** `reports/v0.4_phase6_regression_sensitivity.md`

### Phase 7 — final acceptance

- Run full tests.
- Summarize unresolved issues.
- Produce final acceptance report.

**Deliverable:** `reports/v0.4_Native_BigWig_Continuous_Assay_Acceptance_Report.md`

---

## 26. Stop conditions

Stop and report rather than silently working around the issue if:

1. v0.3 width inconsistency cannot be reproduced or explained;
2. chromosome/genome mapping is ambiguous;
3. BigWig signal semantics are unknown in a way that changes interpretation;
4. the design matrix is non-estimable;
5. inferential FDR is requested without sufficient replication;
6. implementation requires treating normalized continuous signal as raw counts;
7. GSE280947 direction reverses relative to the v0.3 baseline without a clear explanation;
8. a bug fix materially changes conclusions and has not been documented.

---

## 27. Definition of done

v0.4 is complete when a user can supply a valid BigWig sample sheet/config and EpiSegDiff can:

```text
validate input
   -> construct common signal assay
   -> QC replicates
   -> transform/filter explicitly
   -> perform appropriate continuous inference
   -> segment differential signal
   -> classify GAIN/LOSS
   -> export standardized results
   -> preserve complete provenance
   -> generate an acceptance-quality report
```

without pretending that continuous normalized BigWig values are fragment counts.

The release must state clearly:

> **EpiSegDiff v0.4 supports BigWig as a first-class continuous-signal input. Count assays and continuous-signal assays use distinct observation models. Statistical equivalence between BigWig and BAM/count analyses is not assumed and must be evaluated separately when matched data are available.**

---

## 28. Recommended post-v0.4 direction (not part of acceptance)

After v0.4 passes, prioritize a matched BAM + BigWig benchmark to quantify:

- effect-size correlation;
- direction concordance;
- DER overlap/Jaccard;
- rank concordance;
- boundary concordance;
- robustness to normalization;
- conditions under which continuous-signal and count-assay conclusions diverge.

Only after that benchmark should the project make evidence-based statements about when deposited BigWigs can substitute for BAM-derived counts for a given scientific objective.
