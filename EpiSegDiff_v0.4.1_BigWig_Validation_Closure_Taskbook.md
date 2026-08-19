# EpiSegDiff v0.4.1 — Native BigWig Continuous-Assay Validation Closure

**Status:** Agent-executable validation/closure specification  
**Target version:** v0.4.1  
**Baseline:** EpiSegDiff v0.4.0  
**Release type:** Validation/closure release; no major new scientific feature  
**Primary dataset:** GSE280947 JUN continuous BigWig assay  
**Primary goal:** Close all remaining PARTIAL / UNRESOLVED items from v0.4 and determine whether native BigWig support can move from CONDITIONAL PASS to PASS.

---

## 1. Baseline facts from v0.4

The v0.4 acceptance report established:

- BigWig is a first-class **continuous** assay input.
- Count/BAM and continuous/BigWig paths are separated.
- BigWig values are not passed to edgeR.
- GSE280947 reference run uses limma on `log2(score + 0.01)`.
- Contrast is `JNKi - DMSO`.
- Reference result: **2 GAIN / 688 LOSS** at FDR <= 0.05.
- Width enforcement was fixed; emitted intervals are >=600 bp.
- The former 400 bp median-width inconsistency was traced to short CVRs / short terminal PELT pieces and fixed.
- Synthetic BigWigs recover both GAIN and LOSS.
- Unreplicated data cannot produce replicate-based FDR unless explicitly run in descriptive mode.

The remaining v0.4 closure issues are:

1. bin-size sensitivity grid incomplete;
2. pseudocount sensitivity grid incomplete;
3. 500 bp coarsening incomplete;
4. `filter=none` not profiled genome-wide;
5. `summary_fun` accepts non-mean choices but still computes mean;
6. small-variance limma warnings require characterization;
7. QC hexbin warning requires cleanup or documented harmlessness;
8. 200 bp / 0.01 must not be treated as universal defaults until robustness is measured.

This release exists to close these items.

---

# 2. Scope and non-scope

## 2.1 Required in v0.4.1

1. Complete the GSE280947 bin-size sensitivity analysis.
2. Complete the pseudocount sensitivity analysis.
3. Complete 500 bp coarsening from a cached finer-resolution matrix.
4. Run and characterize `filter=none`.
5. Fix or restrict `summary_fun` so the API never claims behavior it does not implement.
6. Quantify stability relative to the v0.4 reference configuration.
7. Characterize the limma small-variance warnings.
8. Fix or suppress only demonstrably harmless plotting warnings.
9. Produce machine-readable sensitivity tables.
10. Produce a final v0.4.1 acceptance report with explicit PASS / CONDITIONAL PASS / FAIL gates.
11. Preserve all v0.4 regression tests and assay-separation guarantees.

## 2.2 Explicitly not required

Do **not** expand this release into:

- BAM-versus-BigWig equivalence testing;
- DESeq2 integration;
- NB-FPOP;
- competitor benchmarking;
- DNA methylation-specific modeling;
- new topology classes;
- publication claims;
- optimization of parameters to maximize the number of DERs.

If BAMs later become available, BAM-versus-BigWig concordance belongs to a subsequent release.

---

# 3. Core scientific rule

Parameter sensitivity is a **robustness analysis**, not parameter selection against an assumed truth set.

Do not choose bin size, pseudocount, signal filter, or segmentation penalty because a setting yields more significant regions.

The primary question is whether the scientific conclusions and strongest differential regions are stable across reasonable parameter choices.

---

# 4. Freeze the v0.4 reference configuration

Create:

```text
config/validation/v0.4.1_gse280947_reference.yaml
```

It must reproduce the v0.4 reference settings as closely as possible:

```yaml
assay_family: continuous
input_type: bigwig
genome: mm10
bin_size: 200
summary_fun: mean
transform:
  method: log2
  pseudocount: 0.01
signal_filter:
  method: absolute
  threshold: 0.6
contrast: JNKi - DMSO
inference:
  method: limma
segmentation:
  method: gaussian_pelt
  min_bins: 3
```

Do not silently change any reference parameter.

Create:

```text
results/v0.4.1/reference_run_manifest.yaml
results/v0.4.1/reference_primary.tsv
```

The reference result should be close to:

```text
2 GAIN
688 LOSS
minimum emitted width 600 bp
```

Exact identity is preferred but not mandatory if implementation changes are documented. If the reference run changes substantially, stop and explain why before interpreting sensitivity.

---

# 5. Task A — Fix `summary_fun` contract

## A1. Current problem

v0.4 accepts values other than `summary_fun = mean` but still computes the mean. This is an API correctness defect.

## A2. Required resolution

### Preferred approach

Implement at least:

```text
mean
max
sum/area
```

Optionally implement `median` only if run-length weighting is correct.

For a genomic bin `B` and piecewise BigWig coverage values, `mean` must be coverage-length weighted:

\[
\bar{x}_B = \frac{\sum_k x_k L_k}{\sum_k L_k}
\]

where `L_k` is overlap length of BigWig interval `k` with bin `B`.

Do not use an unweighted average of RLE runs when widths differ.

For `sum`, define whether the quantity is actually signal area `sum(x_k * L_k)`; rename to `area` if needed.

### Acceptable minimal alternative

If only `mean` is scientifically supported in v0.4.1:

- restrict allowed values to `mean`;
- error explicitly for every other value;
- update schema/config/help/docs.

Silent fallback to mean is prohibited.

## A3. Tests

Use hand-computable synthetic BigWigs covering:

1. constant signal;
2. unequal run lengths;
3. partial overlap at bin edge;
4. missing / uncovered bases;
5. terminal chromosome bin;
6. non-finite handling.

---

# 6. Task B — Build reusable cached 100 bp signal matrix

Create:

```text
cache/gse280947_100bp_signal_matrix.*
cache/gse280947_100bp_signal_matrix_manifest.yaml
```

Preserve:

```text
chrom
start
end
sample columns
valid-bin status
genome
source BigWig checksums
summary_fun
software version
```

Invalidate cache if source checksums, genome, chromosome set, summary statistic, or coordinate convention changes.

Do not silently reuse stale cache.

---

# 7. Task C — Bin-size sensitivity

## C1. Required grid

Run genome-wide:

```text
100 bp
200 bp
500 bp
```

200 bp is the reference. 50 bp remains optional.

## C2. Coarsening

500 bp must be genuinely implemented and analyzed. Prefer:

```text
BigWig
  ↓
linear-scale 100 bp cache
  ↓
coarsen to target resolution
  ↓
transform
  ↓
model
```

Do not average already log-transformed values unless specifically justified.

## C3. Minimum physical segment width

With `min_bins = 3`, report:

```text
100 bp -> 300 bp
200 bp -> 600 bp
500 bp -> 1500 bp
```

This must be considered when comparing width distributions and boundary displacement.

## C4. Required metrics

For every bin size record:

```text
bin_size
number_input_bins
number_tested_segments
GAIN
LOSS
total_DER
median_width
p90_width
median_abs_effect
runtime_sec
peak_memory_mb if available
```

Relative to 200 bp compute:

```text
effect correlation
sign concordance
DER overlap
bp-level Jaccard
reciprocal-overlap distribution
strong-region recovery
boundary displacement
```

Define strong regions before viewing sensitivity results using a deterministic rule.

---

# 8. Task D — Pseudocount sensitivity

At 200 bp run:

```text
0.001
0.01
0.05
0.1
```

Reference is 0.01.

For each value report:

```text
pseudocount
tested_segments
GAIN
LOSS
total_DER
median_width
p90_width
median_effect
effect quantiles
runtime_sec
```

Relative to 0.01 compute:

```text
effect Pearson correlation
effect Spearman correlation
sign concordance
DER bp Jaccard
strong-region recovery
number of sign flips
number of significance-only changes
```

Stratify pseudocount sensitivity by reference linear-scale signal quartile and report effect correlation, median absolute effect difference, and sign-flip rate per stratum.

Do not claim 0.01 is optimal merely because it yields more calls.

---

# 9. Task E — Signal-filter sensitivity

Required profiles:

```text
A. none
B. absolute threshold = 0.6
```

If an existing deterministic control-aware/adaptive method is already production-ready, add one such profile. Do not invent an opaque `auto` rule solely for this task.

`filter=none` must run genome-wide.

Record:

```text
number of bins entering segmentation
number of segments tested
GAIN
LOSS
effect distribution
FDR distribution
runtime
memory
```

Relative to absolute 0.6 compute shared/reference-only/none-only DER bp, direction agreement, and effect correlation on shared evaluable regions.

Low-signal discoveries unique to `none` must be characterized, not automatically called false positives.

---

# 10. Task F — Combined sensitivity matrix

After one-dimensional checks, run:

```text
bin_size ∈ {100, 200, 500}
pseudocount ∈ {0.001, 0.01, 0.05, 0.1}
```

Total: **12 configurations**.

Hold constant:

```text
summary_fun = mean
signal_filter = absolute 0.6
```

If runtime is prohibitive, use the cached 100 bp matrix. If full-genome 12-cell execution still cannot complete, record every unrun cell and reason; overall status remains CONDITIONAL PASS unless an equivalent full validation is scientifically justified.

---

# 11. Task G — Reusable stability metrics

Implement:

```r
compare_continuous_runs(reference, alternative)
```

Return at least:

```text
reference_id
alternative_id
effect_pearson
effect_spearman
sign_concordance
DER_bp_jaccard
reference_DER_recovery
strong_reference_DER_recovery
median_boundary_displacement
p90_boundary_displacement
GAIN_reference
GAIN_alternative
LOSS_reference
LOSS_alternative
```

Do not correlate only significant regions. Define a common evaluable domain and also provide called-region-only metrics as secondary summaries.

For sign concordance, report both all matched units and units above a predeclared absolute-effect threshold.

---

# 12. Task H — Characterize limma small-variance warnings

Do not simply suppress the warning.

Quantify:

```text
number and fraction of affected units
signal distribution of affected units
effect distribution
significant-DER enrichment
low-signal enrichment
```

If compatible, compare default limma empirical Bayes with a documented robust option as a sensitivity analysis. Do not change production behavior solely to eliminate a warning.

Classify the warning as:

```text
harmless / expected
potentially influential
materially influential
```

with quantitative evidence.

If materially influential, v0.4.1 cannot receive full PASS without remediation.

---

# 13. Task I — QC plotting warning

Resolve the current hexbin `grid too coarse` warning by choosing a data-appropriate plot resolution or using ordinary scatter for small datasets.

Do not globally suppress warnings.

Acceptance:

- QC plot generation completes without misleading warnings;
- numerical QC remains unchanged;
- unrelated warnings remain visible.

---

# 14. Task J — Regression suite

All v0.4 critical guarantees remain mandatory.

## J1. Assay separation

`continuous + bigwig` must never dispatch to edgeR/count inference. `count + bam` must preserve the existing count pathway. DESeq2 remains out of scope.

## J2. Direction

For `JNKi - DMSO`:

```text
effect > 0 -> GAIN
effect < 0 -> LOSS
```

Reverse contrast and verify sign reversal.

## J3. Width

For every output configuration:

```text
width >= bin_size * min_bins
```

unless a later explicitly documented boundary-refinement stage changes the rule. No short terminal PELT piece may bypass enforcement.

## J4. Replication safety

With one sample per condition:

- inferential/FDR mode must fail;
- descriptive mode may run;
- inferential fields must be `NA` where undefined.

## J5. Determinism

Same input, config, version, and seed must reproduce identical primary tables within documented numerical tolerance.

---

# 15. Required outputs

Create:

```text
results/v0.4.1/bin_size_sensitivity.tsv
results/v0.4.1/pseudocount_sensitivity.tsv
results/v0.4.1/filter_sensitivity.tsv
results/v0.4.1/combined_sensitivity_grid.tsv
results/v0.4.1/run_concordance.tsv
results/v0.4.1/limma_variance_diagnostics.tsv
results/v0.4.1/reference_primary.tsv
results/v0.4.1/reference_run_manifest.yaml
```

Each sensitivity row must retain:

```text
config_id
bin_size
pseudocount
summary_fun
filter_method
filter_parameter
tested_segments
GAIN
LOSS
median_width
p90_width
runtime_sec
status
failure_reason
```

Do not delete failed configurations. Use `PASS`, `FAILED`, or `NOT_RUN` with explicit reasons.

---

# 16. Required reports

Generate:

```text
reports/v0.4.1_summary_fun_contract_report.md
reports/v0.4.1_bin_size_sensitivity_report.md
reports/v0.4.1_pseudocount_sensitivity_report.md
reports/v0.4.1_filter_sensitivity_report.md
reports/v0.4.1_limma_variance_diagnostic_report.md
reports/v0.4.1_regression_report.md
reports/v0.4.1_final_acceptance_report.md
```

Each report must contain:

```text
Objective
Methods
Parameters
Results
Uncertainty
Failures
Interpretation
Limitations
Decision
Next action
```

Label statements as appropriate:

```text
FACT
INFERENCE
UNRESOLVED
```

---

# 17. Required plots

At minimum:

```text
figures/v0.4.1_gain_loss_by_configuration.*
figures/v0.4.1_effect_correlation_heatmap.*
figures/v0.4.1_der_jaccard_heatmap.*
figures/v0.4.1_strong_region_recovery.*
figures/v0.4.1_width_by_bin_size.*
figures/v0.4.1_pseudocount_low_signal_sensitivity.*
figures/v0.4.1_filter_comparison.*
```

Plots are supplemental; machine-readable metrics are required.

---

# 18. Acceptance gates

## Gate A — API correctness

PASS only if `summary_fun` does what it says or unsupported values fail explicitly. Silent fallback to mean = **FAIL**.

## Gate B — Bin-size closure

PASS only if 100, 200, and 500 bp complete genome-wide, 500 bp is genuinely implemented, minimum-width behavior is correct, and concordance metrics are reported.

Biological instability across scales is not automatically an implementation failure; it means the parameter is context-dependent.

## Gate C — Pseudocount closure

PASS only if all four pseudocount values complete and low-signal sensitivity is quantified.

## Gate D — Filter closure

PASS only if `filter=none` runs genome-wide and is quantitatively compared with the 0.6 reference profile.

## Gate E — Statistical warning closure

PASS only if limma small-variance warnings are quantitatively characterized. Material influence blocks full PASS until remediated.

## Gate F — Regression

PASS only if:

- GSE280947 still shows strong LOSS dominance under reference settings;
- direction is correct;
- minimum width is respected;
- synthetic GAIN/LOSS tests pass;
- unreplicated FDR remains forbidden;
- BigWig never enters edgeR;
- all prior v0.4 critical tests remain passing.

---

# 19. Final release decision

Allowed decisions:

```text
PASS
CONDITIONAL PASS
FAIL
```

## PASS

Use only if all critical gates pass, the required sensitivity grid is complete, no silent API misbehavior remains, warnings are characterized, and v0.4 regressions remain intact.

PASS means:

> Native BigWig / continuous-assay support is accepted as a production-supported EpiSegDiff input modality under documented statistical assumptions.

PASS does **not** mean BigWig is equivalent to BAM fragment-count analysis, nor that 200 bp / 0.01 is universally optimal.

## CONDITIONAL PASS

Use if core implementation remains valid but required validation cells remain incomplete for documented resource reasons or a parameter-sensitive regime remains unresolved.

## FAIL

Use if any of the following occurs:

- `summary_fun` silently misrepresents computation;
- BigWig values reach edgeR/count inference;
- contrast direction is wrong;
- width enforcement regresses;
- unreplicated FDR is fabricated;
- sensitivity reveals a material statistical defect invalidating the current continuous model;
- reference behavior cannot be reproduced without explanation.

---

# 20. Default-selection rule after validation

Do not automatically promote the GSE280947 reference profile to a universal software default.

Classify each parameter as:

```text
stable default candidate
assay-dependent
dataset-dependent
user-required
```

At minimum classify:

```text
bin_size
pseudocount
filter
summary_fun
```

A setting may remain the default for the GSE280947 profile without being a universal BigWig default.

---

# 21. Required final acceptance questions

`reports/v0.4.1_final_acceptance_report.md` must answer:

1. Does `summary_fun` now have an honest implementation contract?
2. Do 100, 200, and 500 bp runs all complete genome-wide?
3. How strongly do effect estimates agree across bin sizes?
4. How stable are strong DERs across bin sizes?
5. How much instability is explained by changing minimum physical segment width?
6. How sensitive are effects and calls to pseudocount 0.001–0.1?
7. Is pseudocount sensitivity concentrated in low-signal regions?
8. What changes under `filter=none`?
9. Are low-signal `filter=none` discoveries reproducible across replicates?
10. Does the full 12-cell bin-size × pseudocount grid complete?
11. Do any configurations reverse the dominant GSE280947 LOSS conclusion?
12. Does any reasonable configuration produce widespread sign reversal?
13. What is the strong-reference-region recovery range?
14. What is the DER bp-Jaccard range relative to reference?
15. What do the limma small-variance warnings affect?
16. Are those warnings harmless, potentially influential, or materially influential?
17. Are all v0.4 critical regression tests still passing?
18. Is BigWig fully separated from edgeR/count inference?
19. Is v0.4.1 PASS, CONDITIONAL PASS, or FAIL?
20. Which parameters are justified as defaults and under what assumptions?
21. What remains necessary before BAM-versus-BigWig concordance testing?
22. Is the project ready for the next scientific-development version?

---

# 22. Execution order

Execute strictly:

```text
Phase 1  Freeze reference configuration
   ↓
Phase 2  Fix summary_fun contract
   ↓
Phase 3  Build/cache 100 bp matrix
   ↓
Phase 4  Bin-size sensitivity
   ↓
Phase 5  Pseudocount sensitivity
   ↓
Phase 6  filter=none / filter sensitivity
   ↓
Phase 7  12-cell combined grid
   ↓
Phase 8  limma warning characterization
   ↓
Phase 9  QC warning cleanup
   ↓
Phase 10 full regression suite
   ↓
Phase 11 final acceptance report
```

If an earlier phase exposes a correctness bug, fix it before interpreting later sensitivity results.

---

# 23. Resource strategy

Because v0.4 left validation incomplete primarily for computational-cost reasons:

1. cache linear-scale 100 bp summarized signal once;
2. reuse it for coarser grids when mathematically valid;
3. avoid repeated BigWig import;
4. record runtime and memory where available;
5. permit HPC execution;
6. do not silently reduce the biological validation grid.

If the full 12-cell genome-wide grid cannot run, record every unrun cell and resource reason. A reduced subset may be supplemental evidence, but final status remains CONDITIONAL PASS unless an equivalent full validation is justified.

---

# 24. Prohibited behavior

The agent must not:

- select parameters because they maximize DER count;
- use deposited IDR peaks as differential ground truth;
- claim BigWig/BAM equivalence;
- globally call BigWig a coverage proxy;
- transform normalized BigWig values into fake counts for edgeR;
- introduce DESeq2 merely because it is common;
- suppress statistical warnings globally;
- omit failed sensitivity runs;
- silently substitute mean for another `summary_fun`;
- loosen width constraints to recover more regions;
- interpret low Jaccard alone as proof one setting is wrong;
- claim publication readiness from this closure release.

---

# 25. Agent start instruction

Execute **EpiSegDiff v0.4.1 validation closure only**.

Use v0.4.0 as baseline.

Start by reproducing the v0.4 reference run and fixing the `summary_fun` contract.

Do not begin v0.5 feature development.  
Do not implement NB-FPOP.  
Do not perform competitor benchmarking.  
Do not claim BAM/BigWig equivalence.

After the sensitivity grid and regression suite are complete, generate:

```text
reports/v0.4.1_final_acceptance_report.md
```

and stop for review.
