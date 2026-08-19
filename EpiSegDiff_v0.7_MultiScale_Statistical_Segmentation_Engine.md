# EpiSegDiff v0.7 — Multi-Scale Statistical Segmentation Engine

**Status:** Agent-executable method-development specification  
**Baseline:** v0.6 Methodological Benchmark — CONDITIONAL PASS  
**Release type:** Core statistical-method redesign and benchmarking  
**Primary goal:** Build and evaluate a scale-adaptive differential segmentation framework that works across narrow and broad epigenomic events, rather than optimizing only a wide-domain regime.

---

## 1. Motivation

v0.6 established that the current `200 bp + limma-derived score + Gaussian PELT + min_bins=3` implementation is a useful baseline, not a universal final method.

Observed limitations include:

- poor performance for narrow 200–500 bp synthetic events;
- improving F1 as event width increases;
- lower synthetic boundary error from window aggregation than the current EpiSegDiff baseline;
- non-trivial same-direction adjacent-event merging;
- incomplete gradual-boundary characterization;
- anti-conservative behavior under unfiltered synthetic nulls;
- unresolved selection bias from using the same data for segmentation and testing.

Therefore v0.7 must not merely tune a new PELT penalty.

The central methodological objective is:

> **Develop a multi-scale differential segmentation framework that can recover narrow local events and broad domains using one coherent architecture, while maintaining accurate boundaries, appropriate event separation, and acceptable statistical calibration.**

---

## 2. Core design principles

Separate the method into four independently testable layers:

```text
Observation representation
        ↓
Differential statistic
        ↓
Segmentation / change-point engine
        ↓
Inference strategy
        ↓
Final differential domains
```

Do not attribute failure of one combination to the entire EpiSegDiff concept.

All outputs must retain provenance:

```text
score_engine
segmentation_engine
scale_mode
scale
boundary_refinement
inference_mode
parameter_set
```

---

## 3. Primary hypotheses

### H1 — Scale adaptivity

The framework should recover both:

```text
narrow: ~100–1000 bp
broad: multi-kb to >50 kb
```

without requiring users to manually switch to unrelated methods.

### H2 — Boundary accuracy

Boundary error should improve relative to the v0.6 baseline.

### H3 — Event separation

Nearby events should remain separate when supported by intervening signal.

### H4 — Broad coherence

Broad events should not be fragmented into many artificial subregions.

### H5 — Statistical calibration

Improved detection must not come from uncontrolled null false positives.

### H6 — Assay-geometry robustness

The framework should work across TF-like, accessibility-like, and broad histone-like signal geometries.

---

## 4. Modular architecture

Implement conceptually:

```text
INPUT
  │
  ├── count-derived signal
  └── continuous BigWig signal
          ↓
MULTI-SCALE REPRESENTATION
          ↓
DIFFERENTIAL SCORE ENGINE
          ↓
SEGMENTATION ENGINE
          ↓
BOUNDARY REFINEMENT
          ↓
INFERENCE ENGINE
          ↓
GAIN / LOSS DOMAINS
```

Suggested interfaces:

```r
compute_diff_score(...)
segment_signal(...)
refine_boundaries(...)
infer_regions(...)
run_multiscale(...)
compare_engines(...)
```

Do not hard-code one preferred engine before benchmarking.

---

## 5. Differential statistic benchmark

Implement and compare at least:

### 5.1 Raw effect

```text
mean(Treatment) - mean(Control)
```

Purpose: preserve biological magnitude and reduce tiny-SE amplification.

### 5.2 Log-effect

For positive continuous signals:

```text
log2(Treatment + c) - log2(Control + c)
```

Pseudocount must remain explicit.

### 5.3 Moderated t statistic

Use the replicate-aware moderated statistic directly.

### 5.4 Standardized effect

Implement a variance-standardized effect distinct from the current ad hoc SE-floor z-score.

Document the exact formula.

### 5.5 Signed significance

Exploratory only:

```text
sign(effect) * -log10(p)
```

Must not become production default without calibration.

### 5.6 Replicate-aware multivariate score

Implement at least one score that preserves replicate structure more directly, e.g.:

- condition-centroid difference with pooled covariance;
- Hotelling-like statistic where valid;
- multivariate Gaussian segment cost;
- replicate-aware likelihood ratio.

Do not add complex scores without synthetic validation.

---

## 6. Segmentation/change-point engines

Use one common interface.

Required engines:

### 6.1 Gaussian PELT

Retain v0.6 baseline.

### 6.2 Binary Segmentation

Simple classical baseline.

### 6.3 Wild Binary Segmentation / WBS2

High priority for narrow/sharp events.

### 6.4 FPOP

Include if computationally practical. If real-data complexity is prohibitive, retain synthetic-only with explicit status.

### 6.5 Fused lasso / total variation

High priority.

Conceptual objective:

\[
\min_{\beta}
\sum_i (y_i-\beta_i)^2
+
\lambda \sum_i |\beta_i-\beta_{i-1}|
\]

or a justified weighted variant.

Important because short and long regions can coexist without a rigid fixed `min_bins` floor.

### 6.6 HMM

Three-state model where appropriate:

```text
LOSS
NEUTRAL
GAIN
```

### 6.7 HSMM / duration-aware model

High priority if feasible.

Purpose: allow different state-duration scales explicitly.

### 6.8 Optional exact DP

Segment-neighborhood or another exact dynamic-programming baseline may be added.

---

## 7. Multi-scale signal representation

Do not force one universal bin size.

Support at minimum:

```text
50 bp
100 bp
200 bp
500 bp
1 kb
```

Optionally for broad domains:

```text
2 kb
5 kb
```

Prefer deriving coarse scales from the highest-resolution **linear-scale** signal rather than averaging transformed values.

---

## 8. Multi-scale strategies

Implement at least three.

### Strategy A — Single-scale baseline

Run each scale independently.

Purpose:

> establish best achievable single-scale behavior per event class.

### Strategy B — Candidate-union multi-scale

```text
50/100 bp → narrow candidates
200/500 bp → medium candidates
1–5 kb → broad candidates
       ↓
candidate union
       ↓
conflict resolution
```

### Strategy C — Hierarchical coarse-to-fine

```text
coarse scale
→ broad candidate
→ local high-resolution refinement
→ final boundaries
```

### Strategy D — Consensus evidence

Optional but desirable.

Combine scale-specific support with a predeclared rule rather than selecting the best scale after seeing truth.

---

## 9. Boundary refinement

Boundary accuracy is a first-class target.

Evaluate optional refinement methods:

- local high-resolution change-point search;
- local likelihood scan;
- local fused-lasso refinement;
- derivative/gradient-based transition search;
- transition-profile fitting.

Refinement must not use synthetic truth.

The same logic must work on real data.

---

## 10. Event-scale output

Every final DER should include:

```text
estimated_scale_bp
source_scale
boundary_confidence
score_engine
segmentation_engine
inference_engine
```

Allow eventual categories:

```text
narrow
intermediate
broad
```

but thresholds must be explicit or learned from validated rules.

---

## 11. Synthetic benchmark redesign

Create a factorial benchmark:

\[
width \times effect \times noise \times boundary\_type
\]

### Width

```text
100 bp
200 bp
500 bp
1 kb
2 kb
5 kb
10 kb
25 kb
50 kb
100 kb
```

### Effect

```text
weak
moderate
strong
```

Define numerically.

### Noise

```text
low
moderate
high
```

### Boundary type

At minimum:

```text
sharp
linear ramp
shoulder
peak-on-domain
```

### Event relationships

Include:

```text
isolated
same-direction adjacent
same-direction separated
opposite-direction adjacent
opposite-direction separated
nested narrow-on-broad
```

Nested events are a critical test of genuine multi-scale capability.

---

## 12. Biological signal geometries

Do not simulate only piecewise-flat blocks.

Required:

### Narrow TF-like
Localized sharp event.

### ATAC-like
Narrow/intermediate event with shoulders.

### Broad histone domain
Multi-kb plateau or gradual domain.

### Peak-on-domain
Narrow event superimposed on a broad differential background.

### Multi-peak domain
Several local peaks sharing one broad condition shift.

---

## 13. Simulation replication

Use:

```text
>=50 replicates per primary cell
```

Prefer 100 where tractable.

If the full factorial grid is too large:

1. retain high replication for primary cells;
2. use lower replication only for explicitly exploratory cells.

Do not silently fall back to n=10.

Store all seeds.

---

## 14. Core metrics

For each combination measure:

### Detection

```text
event sensitivity
event precision
bp recall
bp precision
F1
```

### Direction

```text
GAIN/LOSS accuracy
opposite-direction error
```

### Boundary

```text
left error
right error
median absolute boundary error
p90/p95 error
IoU
```

### Fragmentation

```text
calls per true event
P(single-region reconstruction)
```

### Merging

```text
P(two true events merged)
```

### Splitting

```text
P(one true event split)
```

### Scale

```text
estimated scale vs true width
```

### Calibration

```text
empirical FDP
P(any false DER)
nominal vs empirical FDR
```

---

## 15. Performance phase diagrams

Generate phase diagrams over:

```text
width × effect
width × noise
width × boundary type
effect × noise
```

For each engine/configuration show at least:

- F1;
- boundary error;
- merge rate;
- empirical FDP.

The goal is not one global winner.

The goal is to identify which method is valid in which regime.

---

## 16. Oracle benchmark

For every simulation cell calculate the best result achievable across:

```text
score engines
segmentation engines
single scales
```

This oracle is not deployable.

It answers:

> Is poor performance caused by bad engine/scale choice, or by a deeper formulation problem?

Compare the final multi-scale strategy against the oracle.

---

## 17. Automatic scale/engine selection

Only after the benchmark is complete, evaluate automatic selection.

Candidate inputs:

```text
local signal width
autocorrelation
effect smoothness
replicate variance
candidate-boundary density
```

Candidate outputs:

```text
preferred scale
preferred segmentation engine
boundary refinement depth
```

Train/tune on one simulation subset and evaluate on held-out simulations.

Do not select using the held-out truth.

---

## 18. Inference strategies

Compare at least:

### I1 — Same-data baseline

```text
all data
→ differential score
→ segmentation
→ inference on same data
```

### I2 — Pooled-signal boundaries

Ignore condition labels for geometry:

```text
pooled signal
→ candidate boundaries
→ condition test within fixed regions
```

### I3 — Sample splitting

```text
subset A → segmentation
subset B → testing
```

and reverse.

### I4 — Cross-fitting

```text
A defines / B tests
B defines / A tests
```

Combine evidence with a predeclared rule.

High priority.

### I5 — Permutation-derived inference

Optional calibration reference:

run the full segmentation pipeline under label permutation.

---

## 19. Calibration benchmark

Repeat null testing for every serious production candidate.

Required nulls:

```text
global null
structured null
low-signal null
spatially autocorrelated null
```

Nominal FDR:

```text
0.01
0.05
0.10
```

Report:

```text
mean false regions
false bp
P(any false region)
empirical FDP
p-value histogram
FDR distribution
```

Materially anti-conservative candidates cannot become production defaults.

---

## 20. Power–calibration frontier

For every inference architecture plot:

```text
power
vs
empirical FDP
```

This should reveal whether:

- same-data is powerful but anti-conservative;
- sample split is calibrated but weak;
- pooled-boundary improves validity;
- cross-fitting improves the compromise.

---

## 21. Real-data validation

Use at least three regimes where possible.

### A — Narrow TF-like
GSE280947 JUN or another replicated TF dataset.

### B — Broad histone-like
Replicated H3K27me3/H3K36me3 or comparable broad-mark data.

### C — Intermediate/open chromatin
ATAC-seq or similar.

Compare:

```text
single-scale EpiSegDiff
multi-scale EpiSegDiff
fixed peaks
window aggregation
csaw where BAM exists
THOR where compatible
```

Do not create proxy comparisons for unavailable methods.

---

## 22. Narrow-event acceptance requirement

The final recommended architecture must materially improve the v0.6 deficiency at:

```text
200–500 bp
```

without catastrophic loss of broad-domain performance.

Do not pass this gate merely by globally reducing `min_bins` if false positives or splitting increase substantially.

---

## 23. Broad-event acceptance requirement

The final method must retain useful performance for:

```text
5–100 kb
```

with:

- high bp recall;
- acceptable boundary error;
- low fragmentation;
- acceptable separation of neighboring domains.

---

## 24. Nested multi-scale requirement

Critical test:

```text
broad differential domain
+
embedded narrow event
```

The final architecture must not be forced to choose only one scale.

Where truth contains both components, acceptable behavior should preserve both without artificial fragmentation.

---

## 25. Boundary-confidence estimation

Evaluate at least one boundary-confidence concept:

- bootstrap boundary frequency;
- replicate resampling;
- local likelihood profile;
- change-point support score.

On synthetic data test whether higher confidence predicts lower boundary error.

Do not expose a user-facing confidence score unless validated.

---

## 26. Engine-selection framework

Do not choose the final production engine from one metric.

Use Pareto analysis across:

```text
narrow F1
broad F1
boundary error
merge rate
split rate
empirical FDP
runtime
memory
noise robustness
```

A weighted ranking is allowed only if weights are frozen before final comparison.

Prefer Pareto-optimal solutions.

---

## 27. Production architecture candidates

Compare at least:

### Candidate A — Single universal engine
One segmentation family + multi-scale representation.

### Candidate B — Engine ensemble
Multiple engines propose candidates; consensus/selection yields final calls.

### Candidate C — Scale-specialized routing

Example:

```text
narrow → WBS/fused lasso
medium → PELT/fused lasso
broad → HSMM/PELT
```

with automatic routing.

Do not accept Candidate C merely because it is more complex.

Complexity must produce measurable gain.

---

## 28. Real-data stability

For final candidates assess:

```text
replicate leave-one-out
bootstrap
parameter perturbation
scale perturbation
```

Report:

```text
region recovery
boundary variance
direction stability
```

High synthetic accuracy with unstable real-data boundaries is not sufficient.

---

## 29. Existing-method comparison

v0.7 remains method-development focused, but where compatible BAMs exist include:

```text
csaw
THOR where possible
fixed-peak count workflow
```

Purpose:

> ensure the redesigned multi-scale engine remains competitive across narrow and broad regimes.

Do not claim superiority unless directly supported.

---

## 30. Required outputs

```text
results/v0.7/score_engine_benchmark.tsv
results/v0.7/segmentation_engine_benchmark.tsv
results/v0.7/multiscale_benchmark.tsv
results/v0.7/boundary_refinement.tsv
results/v0.7/narrow_event_metrics.tsv
results/v0.7/broad_event_metrics.tsv
results/v0.7/nested_event_metrics.tsv
results/v0.7/merge_split_metrics.tsv
results/v0.7/calibration_metrics.tsv
results/v0.7/power_calibration_frontier.tsv
results/v0.7/oracle_comparison.tsv
results/v0.7/engine_selection.tsv
results/v0.7/realdata_stability.tsv
results/v0.7/realdata_method_comparison.tsv
results/v0.7/runtime_memory.tsv
results/v0.7/regression_summary.tsv
```

Keep failed/unavailable configurations with reasons.

---

## 31. Required figures

```text
figures/v0.7_score_engine_phase_diagram.*
figures/v0.7_segmentation_engine_phase_diagram.*
figures/v0.7_boundary_error_vs_width.*
figures/v0.7_f1_vs_width.*
figures/v0.7_merge_rate_vs_gap.*
figures/v0.7_fragmentation_vs_width.*
figures/v0.7_nested_scale_recovery.*
figures/v0.7_empirical_fdr.*
figures/v0.7_power_calibration_frontier.*
figures/v0.7_oracle_gap.*
figures/v0.7_realdata_boundary_stability.*
figures/v0.7_runtime_memory.*
```

---

## 32. Required reports

```text
reports/v0.7_architecture_design.md
reports/v0.7_score_engine_benchmark.md
reports/v0.7_segmentation_engine_benchmark.md
reports/v0.7_multiscale_benchmark.md
reports/v0.7_boundary_refinement.md
reports/v0.7_calibration_report.md
reports/v0.7_realdata_narrow_benchmark.md
reports/v0.7_realdata_broad_benchmark.md
reports/v0.7_realdata_intermediate_benchmark.md
reports/v0.7_engine_selection_report.md
reports/v0.7_regression_report.md
reports/v0.7_final_methodological_acceptance_report.md
```

Every report must distinguish:

```text
FACT
INFERENCE
SPECULATION
UNRESOLVED
```

---

## 33. Acceptance gates

### Gate A — Modular correctness
All score/segmentation/inference engines produce reproducible and traceable outputs.

### Gate B — Narrow-event recovery
Final architecture materially improves 200–500 bp performance relative to v0.6.

### Gate C — Broad-event retention
Narrow improvement does not destroy 5–100 kb performance.

### Gate D — Boundary accuracy
Boundary error is competitive across scales and improved where v0.6 was weak.

### Gate E — Merge/split control
Broad coherence is not achieved by excessive merging.

### Gate F — True multi-scale behavior
Nested narrow-on-broad events demonstrate more than one useful scale can be represented.

### Gate G — Calibration — CRITICAL
Production candidates show acceptable empirical error under predeclared null regimes.

### Gate H — Selection-bias mitigation — CRITICAL
At least one inference strategy improves calibration relative to the v0.6 same-data unfiltered baseline without unacceptable power loss.

### Gate I — Real-data stability
Final calls/boundaries remain stable under replicate and parameter perturbation.

### Gate J — Existing-method competitiveness
Where compatible data exist, benchmark against csaw and/or THOR.

### Gate K — Regression
All validated v0.4–v0.6 infrastructure guarantees remain intact.

---

## 34. Decision rules

Allowed:

```text
PASS
CONDITIONAL PASS
FAIL
```

### PASS

Requires:

- useful narrow-event recovery;
- retained broad-event performance;
- genuine multi-scale behavior;
- acceptable boundary accuracy;
- controlled merging/splitting;
- acceptable calibration;
- defensible production architecture.

PASS does not require one engine to dominate every regime.

### CONDITIONAL PASS

Use if multi-scale methodology is validated but one critical issue, such as calibration or real-data narrow validation, remains incomplete.

### FAIL

Use if:

- narrow recovery remains fundamentally poor;
- broad performance collapses after narrow optimization;
- multi-scale output is unstable;
- calibration remains materially anti-conservative;
- automatic selector overfits simulations;
- held-out performance collapses;
- complexity adds no measurable benefit;
- critical regressions fail.

---

## 35. Required final questions

The final acceptance report must answer:

1. Which differential statistic is best for narrow events?
2. Which is best for broad events?
3. Which score is best calibrated?
4. Which segmentation engine gives best narrow boundary recovery?
5. Which gives best broad-domain coherence?
6. Which has the lowest merge rate?
7. Which has the lowest split rate?
8. Does fused lasso materially improve the baseline?
9. Does WBS/WBS2 improve narrow sharp events?
10. Does HMM/HSMM help broad/duration-dependent events?
11. Does FPOP add practical value?
12. What is the best single scale at each true width?
13. How close does multi-scale come to the oracle?
14. Does coarse-to-fine improve boundaries?
15. Does candidate-union improve detection?
16. Can nested narrow-on-broad events be recovered?
17. How does boundary error vary with effect and noise?
18. How do gradual boundaries affect each engine?
19. What is empirical FDP for each inference architecture?
20. Does pooled-boundary inference improve calibration?
21. How much power does sample splitting lose?
22. Does cross-fitting improve the power/calibration tradeoff?
23. Is same-data inference still anti-conservative?
24. Which production architecture is Pareto-optimal?
25. Is automatic scale selection reliable on held-out simulations?
26. Does final EpiSegDiff improve 200–500 bp performance?
27. Does it retain 5–100 kb performance?
28. Are real-data boundaries stable?
29. Does the final method remain directionally coherent on JUN?
30. Does it remain coherent on broad histone marks?
31. Does it work on intermediate/open chromatin?
32. How does it compare with csaw where BAM is available?
33. How does it compare with THOR where compatible?
34. Which method components are retained for production?
35. Which components are rejected?
36. What claims are now supported?
37. What failure regimes remain?
38. Is v0.7 PASS, CONDITIONAL PASS, or FAIL?
39. Is the core statistical engine ready for manuscript-focused benchmarking?

---

## 36. Execution order

```text
Phase 1
modularize score / segmentation / inference interfaces
    ↓
Phase 2
build full factorial synthetic benchmark
    ↓
Phase 3
benchmark differential statistics
    ↓
Phase 4
benchmark segmentation engines
    ↓
Phase 5
benchmark single-scale behavior
    ↓
Phase 6
implement multi-scale strategies
    ↓
Phase 7
boundary refinement
    ↓
Phase 8
nested-scale benchmark
    ↓
Phase 9
null calibration
    ↓
Phase 10
pooled-boundary / sample-split / cross-fit inference
    ↓
Phase 11
power-calibration frontier
    ↓
Phase 12
held-out simulation evaluation
    ↓
Phase 13
real narrow benchmark
    ↓
Phase 14
real broad benchmark
    ↓
Phase 15
real intermediate benchmark
    ↓
Phase 16
csaw/THOR compatibility benchmark
    ↓
Phase 17
engine selection / Pareto analysis
    ↓
Phase 18
regression suite
    ↓
final acceptance
```

---

## 37. Stop conditions

Stop and report if:

- synthetic benchmark reveals a fundamental correctness defect;
- calibration is materially invalid across all inference strategies;
- multi-scale merging introduces uncontrolled event fusion;
- held-out simulation performance collapses;
- complexity increases without measurable benefit;
- v0.4–v0.6 validated infrastructure regresses.

Do not hide failed engines. Rejected methods are useful results.

---

## 38. Agent instruction

Implement **EpiSegDiff v0.7 Multi-Scale Statistical Segmentation Engine**.

Treat the v0.6 `200 bp + PELT + min_bins=3` workflow as a baseline only.

Systematically benchmark:

1. multiple differential statistics;
2. multiple segmentation/change-point engines;
3. multiple genomic scales;
4. boundary refinement;
5. narrow, broad, gradual, adjacent, and nested geometries;
6. merge/split behavior;
7. multiple inference architectures;
8. empirical calibration;
9. held-out simulation generalization;
10. real narrow, intermediate, and broad assays.

High-priority methods:

```text
effect-based score
moderated t
PELT
WBS/WBS2
fused lasso / total variation
HMM/HSMM
FPOP where practical
```

High-priority inference strategies:

```text
same-data baseline
pooled-signal boundaries
sample splitting
cross-fitting
permutation reference
```

Do not tune against one dataset.

Do not redefine EpiSegDiff as a wide-only method.

Do not sacrifice broad performance merely to repair narrow performance.

Do not claim FDR calibration without empirical evidence.

Do not select the final engine using training simulations only.

Finally generate:

```text
reports/v0.7_final_methodological_acceptance_report.md
```

Then stop for formal review.
