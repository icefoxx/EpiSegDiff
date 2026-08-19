# Benchmark plan

**Status:** Phase 0 design lock

---

## 1. Principles

- The same simulated truth is used for every method.
- Comparator versions, commands, parameters, and preprocessing are
  recorded. Failed competitor runs stay in the log.
- No manual removal of unfavorable cases.
- BH q-values after adaptive discovery are not treated as guaranteed
  FDR. Empirical nulls are required.
- No complexity claims before measurement.

---

## 2. Count-assay simulator

Implemented by `simulate_count_assay()`.

### Factors

```text
replicates: 2, 3, 5
depth: low, medium, high
dispersion: low, medium, high
SNR / effect size
region width: 200–500 bp, 1–5 kb, 20–100 kb, optional >100 kb
background heterogeneity
mappability masks
global normalization shifts
event class
```

### Ground-truth event classes

```text
GAIN LOSS EXPANSION CONTRACTION
LEFT_SHIFT RIGHT_SHIFT
SPLIT MERGE
INTERNAL_GAIN INTERNAL_LOSS
NULL
```

Truth table (`benchmark/simulation/truth.tsv`):

```text
chrom start end event_class effect_size control_start control_end
treatment_start treatment_end sim_id
```

---

## 3. Metrics

### Regions

- interval precision, recall, F1
- IoU / Jaccard (matched pairs)
- base-pair precision / recall

Matching: a predicted interval matches a truth interval if
`reciprocal_overlap >= 0.2` (configurable). One-to-one greedy by
overlap. The **same matcher** is used for EpiSegDiff and competitors.

### Boundaries

For matched start (end) coordinates \(b, \hat b\):

- MAE, median absolute error, 90th percentile
- Hausdorff distance between boundary sets

### Event class

- per-class accuracy
- macro-F1
- confusion matrix

### Calibration

At nominal \(q \in \{0.01, 0.05, 0.10\}\):

- region-level empirical FDR
- base-pair-level empirical FDR
- false boundary count / rate

Null designs: simulated no-difference genomes, mock vs mock, condition
label permutation, replicate splits.

### Robustness

Replicate dropout, downsampling, normalization, dispersion, bin size,
penalty.

### Performance

Runtime, peak memory, chromosome scaling. Measure, do not claim.

---

## 4. Comparators (count)

| Method | Status in v0.2 |
|--------|----------------|
| csaw | attempt; record if package missing |
| MACS3 + DiffBind | attempt; external tools |
| THOR | attempt; external tool |
| broad-domain baseline | only after a written suitability note |

DiffSegR may be added as a conceptual comparator when technically
meaningful.

### Competitor log fields

```text
software version, command, parameters, input transformation,
runtime, exit status, output parsing rule
```

---

## 5. Comparators (methylation)

DSS, dmrseq, metilene; MethCP if installable. Methylation never uses
the count-assay NB likelihood.

---

## 6. Real data (not bundled)

At least:

- Dataset A: narrow / local (TF CUT&RUN/ChIP, ATAC, or H3K27ac)
- Dataset B: broad (H3K27me3 or H3K9me3)
- Dataset C: WGBS/RRBS two conditions

Validation questions are those in the project specification §19.
No biological claim from EpiSegDiff-unique regions alone.

Required analyses: replicate concordance, downsampling, bin-size,
penalty, normalization sensitivity; overlap with a standard differential
method; method-unique regions.

---

## 7. Deliverables

```text
benchmark/simulation/truth.tsv
benchmark/simulation/results.tsv
reports/Phase2_FDR_Calibration_Report.md
reports/Phase3_Simulation_Benchmark_Report.md
reports/Phase4_Methylation_Validation_Report.md
```
