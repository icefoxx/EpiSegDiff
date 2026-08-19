# EpiSegDiff v0.2.1 --- Statistical Validation and Normalization Task

**Purpose:** Resolve the main statistical risks identified after the
v0.2 MVP before real-data benchmarking or experimental NB-FPOP
development.

**Execution rule:** Complete this as a validation release. Do **not**
implement NB-FPOP in this phase.

## 1. Objectives

1.  Quantify and reduce normalization-induced differential-effect bias.
2.  Replace ambiguous topology definitions with mathematically
    non-overlapping rules.
3.  Expand null and near-null statistical calibration.
4.  Verify topology/boundary inference does not depend on externally
    predefined peak boundaries.

Passing unit tests alone is insufficient.

## 2. Required execution order

``` text
A. Normalization stress test
        ↓
B. Topology ontology correction
        ↓
C. Boundary-inference independence validation
        ↓
D. Expanded null / near-null calibration
        ↓
E. Integrated simulation benchmark
        ↓
F. v0.2.1 acceptance report
```

Do not begin real-data competitor benchmarking before the final gate.

# 3. Task A --- Normalization Stress Test

## A1. Problem

The current simulation indicates that library-size normalization can
transform:

``` text
true:
STABLE | GAIN | STABLE
```

into:

``` text
estimated:
weak LOSS | GAIN | weak LOSS
```

when a large treatment-specific gain changes library composition.

This can create false flanks, inaccurate boundaries, false
expansion/contraction, and incorrect topology labels.

## A2. Methods

Benchmark at least:

-   `library_size`
-   `TMM`
-   DESeq2-style median-of-ratios
-   a defensible reference/background normalization if implementable

Spike-in normalization may be included as a separate coherent simulation
branch. Do not simulate it by assigning arbitrary size factors.

## A3. Differential-genome fraction

Simulate:

``` text
0%
1%
5%
10%
25%
50%
```

For non-null cases vary:

``` text
GAIN:LOSS
50:50
70:30
90:10
100:0
```

Also vary replicates (2,3,5), depth, dispersion, effect size, and region
scale.

If the full Cartesian product is too expensive, document a balanced
reduced design.

## A4. Metrics

For true background bins (B):

\[
Bias\_{background}=`\operatorname{median}`{=tex}\_{i`\in `{=tex}B}(`\hat`{=tex}`\beta`{=tex}\_i)
\]

\[
MAE\_{background}=`\operatorname{median}`{=tex}\_{i`\in `{=tex}B}\|`\hat`{=tex}`\beta`{=tex}\_i\|
\]

and:

\[ Bias_Z=`\operatorname{median}`{=tex}\_{i`\in `{=tex}B}(Z_i) \]

Report:

-   background effect bias and MAE
-   background Z bias
-   false differential regions/bp
-   false boundary count/rate
-   left/right boundary MAE
-   true-region recall and precision
-   topology classification accuracy

Create diagnostics of differential fraction versus background bias,
false-boundary rate, and boundary MAE, stratified by normalization and
gain/loss imbalance.

## A5. Gate

Do not choose normalization solely by power. Prefer the strategy that
jointly minimizes background bias, false boundaries, and boundary error
while retaining acceptable power.

Produce a validity-envelope table for each normalization strategy. Do
not fill PASS/FAIL cells until simulations provide evidence.

# 4. Task B --- Correct the Topology Ontology

For matched domains define:

\[ `\Delta`{=tex}\_L=L_T-L_C,`\qquad `{=tex}`\Delta`{=tex}\_R=R_T-R_C \]

and:

\[ W_C=R_C-L_C,`\qquad `{=tex}W_T=R_T-L_T,`\qquad `{=tex}R_W=W_T/W_C \]

Use configurable displacement (`\delta`{=tex}) and width thresholds.

## Required geometry

**RIGHT_SHIFT:** (`\Delta`{=tex}\_L\>`\delta`{=tex}) and
(`\Delta`{=tex}\_R\>`\delta`{=tex}), with approximately preserved width.

**LEFT_SHIFT:** (`\Delta`{=tex}\_L\<-`\delta`{=tex}) and
(`\Delta`{=tex}\_R\<-`\delta`{=tex}), with approximately preserved
width.

**EXPANSION:** canonical case (`\Delta`{=tex}\_L\<-`\delta`{=tex}) and
(`\Delta`{=tex}\_R\>`\delta`{=tex}), with increased width. Explicitly
define one-sided expansion.

**CONTRACTION:** canonical case (`\Delta`{=tex}\_L\>`\delta`{=tex}) and
(`\Delta`{=tex}\_R\<-`\delta`{=tex}), with reduced width. Explicitly
define one-sided contraction.

## BIDIRECTIONAL_SHIFT

Preferred action: remove it as a primary class.

Retain it only if a mathematically non-overlapping, biologically
interpretable definition can be established that is not expansion,
contraction, or shift.

## Event hierarchy

Separate Level I intensity (`GAIN/LOSS/STABLE`) from Level II
geometry/topology.

Recommended Level II decision order:

``` text
one-to-many / many-to-one
    → SPLIT / MERGE
internal differential structure
    → INTERNAL_GAIN / INTERNAL_LOSS
one-to-one geometry
    → EXPANSION / CONTRACTION / LEFT_SHIFT / RIGHT_SHIFT
otherwise
    → COMPLEX
```

Generate at least 200 synthetic examples per primary topology class,
including threshold-borderline examples.

Report confusion matrix, macro-F1, per-class precision/recall, and
COMPLEX rate.

# 5. Task C --- Boundary-Inference Independence

Differential discovery must remain:

\[ `\hat`{=tex}`\beta`{=tex}(x)`\rightarrow `{=tex}segmentation \]

and must not require MACS/SEACR or another peak caller.

Within a local window around each differential candidate, infer
condition-specific boundaries:

\[ B_C=(L_C,R_C),`\qquad `{=tex}B_T=(L_T,R_T) \]

from control/treatment signals using a documented local method such as
local changepoints, local likelihood, data-driven signal/background
thresholds, replicate consensus, or cross-scale consensus.

External peak callers may be used **only as diagnostic comparators**,
not mandatory inputs.

For simulated truth report:

-   left-boundary MAE
-   right-boundary MAE
-   median and 90th-percentile error
-   Hausdorff distance
-   support-interval coverage if implemented

Stratify by narrow/broad signal, effect strength, replicate number,
depth, and normalization method.

# 6. Task D --- Expanded Null Calibration

## D1. Strict null

Target at least 500 independent simulations per primary null regime. If
computational cost prevents this, use at least 200 and report binomial
confidence intervals.

Zero discoveries must not be interpreted as proof of FDR control.

## D2. Near-null

Add weak structured effects, for example:

``` text
log2FC = ±0.1
log2FC = ±0.2
log2FC = ±0.3
```

Test whether segmentation aggregates weak structured deviations into
apparently strong regions.

## D3. Spatially heterogeneous backgrounds

Include:

-   slow background drift
-   local high-background blocks
-   coverage heterogeneity
-   masked gaps

Clearly distinguish strict null from structured-background simulations.

## D4. Outputs

At nominal q = 0.01, 0.05, 0.10 report:

-   probability of any false discovery
-   mean/median false-region count
-   false-bp fraction
-   false-boundary count
-   empirical region FDR when defined
-   empirical bp-level FDR
-   uncertainty intervals

## D5. Selection-bias experiment

Compare adaptive discovery + same-data testing against at least one
exploratory calibration strategy, e.g. sample splitting or
permutation-calibrated region statistics.

The goal is measurement, not necessarily replacement of the default
pipeline.

# 7. Task E --- Integrated Boundary-Dynamics Simulation

After A-D, rerun the full pipeline on:

``` text
GAIN
LOSS
EXPANSION
CONTRACTION
LEFT_SHIFT
RIGHT_SHIFT
SPLIT
MERGE
INTERNAL_GAIN
INTERNAL_LOSS
NULL
```

Do not include `BIDIRECTIONAL_SHIFT` unless Task B establishes a
non-overlapping definition.

Report region precision/recall/F1, bp precision/recall, IoU, left/right
boundary MAE, Hausdorff distance, event precision/recall/F1, COMPLEX
rate, and false-topology rate.

Stratify at minimum by normalization method, region scale, effect size,
and replicate number.

# 8. Required Code/Script Changes

Update or create:

``` text
R/normalization_validation.R
R/topology.R
R/boundary_refinement.R
R/simulation.R
R/calibration.R

scripts/run_normalization_stress_test.R
scripts/run_topology_validation.R
scripts/run_extended_null_calibration.R
scripts/run_boundary_independence_test.R
scripts/run_v021_integrated_benchmark.R
```

Create validation configs under:

``` text
config/validation/
```

for normalization stress, null calibration, topology validation, and
integrated benchmarking.

# 9. Required Output Tables

Create:

-   `normalization_stress_results.tsv`
-   `topology_validation.tsv`
-   `null_calibration.tsv`
-   `boundary_validation.tsv`

Each row must retain simulation ID, parameter regime, method, truth,
prediction, and relevant metrics.

# 10. Required Reports

Generate:

``` text
reports/v0.2.1_Normalization_Stress_Test_Report.md
reports/v0.2.1_Topology_Ontology_Report.md
reports/v0.2.1_Boundary_Independence_Report.md
reports/v0.2.1_Extended_FDR_Calibration_Report.md
reports/v0.2.1_Integrated_Simulation_Report.md
reports/v0.2.1_Final_Acceptance_Report.md
```

Every report must contain:

``` text
Objective
Hypothesis
Methods
Simulation design
Parameters
Results
Uncertainty
Failures
Interpretation
Limitations
Decision
Next action
```

Explicitly label:

-   **FACT:** direct experimental result
-   **INFERENCE:** interpretation supported by results
-   **UNRESOLVED:** not established by current evidence

# 11. Acceptance Gates

## Gate A --- Normalization

PASS only if background bias is quantified, composition-bias failure
regimes are identified, and at least one strategy performs acceptably in
the intended primary regime.

If none is robust to large global shifts, document the limitation.

## Gate B --- Topology

PASS only if primary classes are mathematically non-overlapping,
expansion/contraction cannot be mislabeled as generic bidirectional
shift by definition, borderline cases are handled, and COMPLEX remains
available.

## Gate C --- Boundary inference

PASS only if simulated boundaries are quantitatively recoverable,
topology does not require external peaks, and boundary errors are
reported by event scale.

## Gate D --- Calibration

PASS or CONDITIONAL PASS only after hundreds of null simulations,
uncertainty intervals, spatially heterogeneous backgrounds, near-null
scenarios, and quantification of adaptive-selection behavior.

# 12. Final Decision

Allowed decisions:

``` text
PASS
CONDITIONAL PASS
FAIL
```

**PASS:** proceed to real narrow/local assay → real broad mark →
competitor benchmarking.

**CONDITIONAL PASS:** proceed only to an explicitly scoped real-data
pilot.

**FAIL:** fix the statistical failure before large real-data
benchmarking.

# 13. Prohibited in v0.2.1

Do not:

-   implement production NB-FPOP
-   optimize C++ before statistical validation
-   claim publication readiness
-   claim selective FDR control from BH alone
-   claim "first" without literature verification
-   tune simulations only to favorable regimes
-   remove failed simulation cells
-   require peak calls for differential discovery
-   force ambiguous topology into named classes
-   choose normalization solely because it yields more discoveries

# 14. Final Acceptance Questions

`reports/v0.2.1_Final_Acceptance_Report.md` must answer:

1.  Which normalization method should be default, and under what
    assumptions?
2.  Under what differential-genome fraction or gain/loss imbalance does
    it fail?
3.  Is background differential-effect bias acceptably centered?
4.  Does segmentation create systematic false flanks?
5.  What are median and 90th-percentile boundary errors?
6.  Are EXPANSION, CONTRACTION, LEFT_SHIFT, and RIGHT_SHIFT
    distinguishable?
7.  Was `BIDIRECTIONAL_SHIFT` removed or rigorously redefined?
8.  How often is `COMPLEX` returned?
9.  Does topology inference require external peak boundaries?
10. What happens under strict-null, near-null, and spatially
    heterogeneous simulations?
11. Is same-data discovery/testing anti-conservative, conservative, or
    approximately calibrated?
12. Is EpiSegDiff ready for real-data benchmarking?
13. What unresolved accuracy problem, if any, would justify future
    NB-FPOP development?

# 15. Agent Start Instruction

Execute **EpiSegDiff v0.2.1 statistical validation only**.

Start with Task A. Do not implement NB-FPOP and do not begin competitor
benchmarking.

Use the current v0.2 code as baseline. Preserve all existing tests and
add validation tests without weakening previous acceptance criteria.

After Tasks A-D, stop and generate their reports before running Task E.

After Task E, generate `reports/v0.2.1_Final_Acceptance_Report.md` and
stop.

Do not proceed to real-data analysis until that report gives PASS or an
explicitly scoped CONDITIONAL PASS.
