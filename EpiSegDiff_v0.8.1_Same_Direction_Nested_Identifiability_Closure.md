# EpiSegDiff v0.8.1 — Same-Direction Nested Identifiability Closure

**Status:** Agent-executable focused methodological closure  
**Basis:** `v0.8_final_architecture_acceptance.md`  
**v0.8 verdict:** `PATH_B_CONDITIONAL_PASS`  
**Reference implementation:** `B_hierarchy`  
**Purpose:** Determine why same-direction local changes are absorbed by broad parent domains, characterize the empirical identifiability boundary, and test whether nested recovery can improve without sacrificing leakage/null control.  
**Production freeze:** **NOT AUTHORIZED**  
**v0.9:** **NOT AUTHORIZED**

## 1. Evidence basis

The v0.8 acceptance report established approximately:

```text
B_hierarchy:
same-direction nested p_both = 0.16
opposite-direction nested p_both = 0.40
broad→local leakage = 0.02
false_child ≈ 0.03
broad performance ≈ 1.00
held-out same-direction nested ≈ 0.10
```

Candidate C had stronger nested recovery but unacceptable leakage (~0.34). Candidate A had low leakage but nested `p_both=0`. Therefore the immediate question is:

> Can same-direction nested structure become identifiable without converting broad structure or noise into false local children?

## 2. Working model and competing explanations

For broad parent effect B and embedded local increment L:

`D(x) = B(x) + L(x) + error(x)`.

Same-direction children change magnitude without changing sign, so the current parent model may absorb them. This is an inference, not an established fact.

Test at least:

- H1: fundamental signal-identifiability limit;
- H2: parent estimator over-absorbs local signal;
- H3: residual contains child but child penalty/threshold is too conservative;
- H4: observation resolution dilutes local contrast;
- H5: parent-boundary proximity causes confounding;
- H6: noise/replicate variance drives the apparent |L|≈1.2–1.8 threshold;
- H7: reconciliation removes valid children after initial detection.

Do not assume H1 until H2–H7 are tested.

## 3. Instrument the hierarchy pipeline

Separate:

```text
parent estimation
→ residual/local contrast construction
→ child candidate detection
→ child threshold/support
→ parent-child reconciliation
→ final hierarchy
```

For every simulated event retain:

```text
truth_B, truth_L
broad_width, local_width
local_position, distance_to_parent_boundary
resolution, noise, replicate_variance
parent_detected_initial, parent_boundary_initial, parent_effect_estimate
local_signal_before_parent_fit, local_signal_after_parent_fit
residual_local_effect, residual_SNR
child_candidate_generated, child_candidate_boundary
child_candidate_effect, child_support
child_removed_by_threshold
child_removed_by_reconciliation
child_merged_into_parent
parent_final, child_final, hierarchy_correct
```

The objective is to locate where same-direction children disappear.

## 4. Primary B × L identifiability grid

Suggested broad effects:

```text
B = 0.5, 0.75, 1.0, 1.25, 1.5, 2.0
```

Suggested local increments:

```text
L = 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0
```

Primary focus is B>0, L>0. Retain negative-L attenuation/reversal controls.

Use >=100 replicates per primary cell where tractable; minimum 50. Do not use n=10 for headline conclusions.

## 5. Width-ratio grid

Parent widths:

```text
5 kb, 10 kb, 25 kb, 50 kb, 100 kb
```

Child widths:

```text
100 bp, 200 bp, 300 bp, 500 bp, 750 bp, 1 kb, 2 kb, 5 kb
```

Analyze performance against child width, parent width, and `W_local/W_parent`.

## 6. Child-position / boundary-distance benchmark

Place children at:

```text
center
25% position
10% position
near left boundary
near right boundary
```

Where possible parameterize physical distance to nearest parent boundary:

```text
0, 100 bp, 200 bp, 500 bp, 1 kb, 2 kb, 5 kb
```

Distinguish nested-identifiability failure from boundary-confounding failure.

## 7. Resolution benchmark

For nested truth test:

```text
50 bp
100 bp
200 bp
```

Optionally use 25 bp for selected difficult cells. Do not infer nested resolution-independence from the earlier isolated-event R1 experiment.

## 8. Noise and replicate variance

Test at least low/moderate/high noise. If possible vary between-replicate SD separately from within-region noise.

The apparent |L| threshold must be reported conditional on SNR/noise, not as a universal constant.

Calculate interpretable contrasts such as:

`C1 = |L| / |B|`

and a precisely defined local SNR such as:

`C2 = |L| / sigma_local`.

Determine whether recovery is better predicted by absolute L, L/B, local SNR, width ratio, or boundary distance.

## 9. Parent-over-absorption experiment — CRITICAL

Quantify how much local increment remains after parent estimation.

Conceptually define an absorption quantity:

`A = 1 - estimated_residual_local_effect / truth_L`.

Interpret cautiously:

```text
A ≈ 0  local signal preserved
A ≈ 1  local signal absorbed
A outside [0,1]  estimator distortion/overshoot
```

Use an operational definition appropriate to the implementation.

Plot absorption against B, L, widths, position, resolution and noise.

## 10. Controlled parent-model sensitivity

Test the current parent estimator plus only a small predeclared set (2–4) of controlled alternatives/regularization strengths.

Ask whether a smoother/coarser parent:

- preserves more child residual;
- damages parent boundaries;
- biases broad effects;
- creates edge residuals.

No unrestricted hyperparameter search.

## 11. Child-detector sensitivity

On the same residual compare, where available:

```text
current B_hierarchy child detector
fused-based child detector
standardized/local-contrast detector
```

Determine whether the child is absent from the residual or present but missed by detection.

## 12. Reconciliation ablation

Measure hierarchy at:

```text
A. raw child candidates
B. thresholded candidates
C. final reconciled hierarchy
```

Report:

```text
P(candidate exists)
P(candidate passes threshold | exists)
P(child retained | passes)
P(final parent+child correct)
```

If reconciliation removes substantial true children, fix reconciliation before changing the decomposition objective.

## 13. Leakage and null guardrails — CRITICAL

Every proposed improvement must simultaneously run:

### Broad-only truth (L=0)

Measure:

```text
false child probability
broad→local leakage
false child count per parent
false local bp
```

### Global null

Measure:

```text
P(any false)
false regions
false hierarchy
```

### Heteroscedastic null

Ensure variance changes are not interpreted as local children.

The v0.8 B reference is approximately:

```text
same-direction p_both = 0.16
broad→local leakage = 0.02
false_child ≈ 0.03
```

A provisional engineering guardrail is broad→local leakage <=0.05. This is not a formal statistical guarantee. Report the full leakage curve.

Do not reproduce Candidate C's failure mode: high nested recovery with leakage ≈0.34.

## 14. Primary success view

Always evaluate:

`P(correct parent AND correct child)` versus `false-child/leakage`.

Do not optimize p_both alone.

Also report:

```text
p_parent
p_child
p_both
p_merged
p_parent_only
p_child_only
p_false_child
broad_to_local_leakage
local_to_broad_leakage
parent_boundary_error
child_boundary_error | detected
child_detection_probability
residual_effect_retention
runtime
```

Include confidence intervals where practical and document the interval method.

## 15. Empirical identifiability surface

Estimate:

`P(both correct | B,L,W_parent,W_local,noise,position,resolution)`.

Do not seek one universal effect threshold. Characterize regimes such as:

```text
IDENTIFIABLE
AMBIGUOUS
NOT RELIABLY IDENTIFIABLE
```

Provisional research definitions may be:

```text
IDENTIFIABLE: p_both >= 0.8 and leakage guardrail satisfied
AMBIGUOUS: 0.2 < p_both < 0.8
NOT RELIABLY IDENTIFIABLE: p_both <= 0.2
```

Freeze definitions before held-out evaluation. Do not describe them as biological laws.

If useful, evaluate an `AMBIGUOUS_CHILD` output rather than forcing false certainty.

## 16. Diagnostic statistical modeling

Use simple interpretable models to explain `P(both correct)` using:

```text
L
B
L/B
local SNR
child width
parent width
width ratio
boundary distance
resolution
noise
```

Purpose: mechanism interpretation, not predictive leaderboard optimization.

## 17. Held-out validation — CRITICAL

After diagnostics:

1. select at most one modified B_hierarchy configuration;
2. freeze parent estimator;
3. freeze child detector;
4. freeze thresholds/reconciliation;
5. evaluate untouched seeds/cells.

Held-out tests must include:

```text
easy same-direction nested
borderline same-direction nested
weak child
near-boundary child
multiple width ratios
moderate/high noise
broad-only truth
heteroscedastic null
```

No post-hoc retuning.

## 18. Required figures

```text
figures/v0.8.1_BxL_same_direction_heatmap.*
figures/v0.8.1_identifiability_surface.*
figures/v0.8.1_pboth_vs_L_over_B.*
figures/v0.8.1_pboth_vs_local_SNR.*
figures/v0.8.1_width_ratio.*
figures/v0.8.1_boundary_distance.*
figures/v0.8.1_resolution_effect.*
figures/v0.8.1_parent_absorption.*
figures/v0.8.1_stagewise_child_loss.*
figures/v0.8.1_nested_vs_leakage_pareto.*
figures/v0.8.1_false_child_guardrail.*
figures/v0.8.1_heldout_identifiability.*
```

## 19. Required tables

```text
results/v0.8.1/BxL_grid.tsv
results/v0.8.1/width_ratio.tsv
results/v0.8.1/position_boundary.tsv
results/v0.8.1/resolution.tsv
results/v0.8.1/noise_variance.tsv
results/v0.8.1/parent_absorption.tsv
results/v0.8.1/child_detector_ablation.tsv
results/v0.8.1/reconciliation_ablation.tsv
results/v0.8.1/leakage_guardrail.tsv
results/v0.8.1/identifiability_surface.tsv
results/v0.8.1/heldout.tsv
results/v0.8.1/regression_summary.tsv
```

## 20. Required reports

```text
reports/v0.8.1_same_direction_failure_localization.md
reports/v0.8.1_parent_absorption_analysis.md
reports/v0.8.1_child_detector_analysis.md
reports/v0.8.1_reconciliation_analysis.md
reports/v0.8.1_identifiability_report.md
reports/v0.8.1_leakage_guardrail_report.md
reports/v0.8.1_heldout_validation.md
reports/v0.8.1_final_closure_decision.md
```

Every report must label conclusions as FACT, INFERENCE, SPECULATION, or UNRESOLVED.

## 21. Critical gates

- **Gate A — Failure localization:** identify whether loss occurs during parent absorption, child detection, thresholding, or reconciliation.
- **Gate B — Same-direction improvement — CRITICAL:** meaningful improvement over the v0.8 B reference in defined identifiable regimes.
- **Gate C — Leakage — CRITICAL:** gain must not arise from broad→local over-decomposition.
- **Gate D — Null — CRITICAL:** global and heteroscedastic null behavior remains acceptable.
- **Gate E — Broad integrity:** parent broad performance/boundaries remain strong.
- **Gate F — Width robustness:** gain is not limited to one hand-picked width.
- **Gate G — Noise robustness:** gain persists outside easiest noise.
- **Gate H — Boundary-position robustness:** characterize/improve near-boundary behavior.
- **Gate I — Held-out validation — CRITICAL:** gain persists after freezing.
- **Gate J — Regression/BigWig:** validated infrastructure and BigWig path remain functional.

## 22. Decision states

Allowed:

```text
IDENTIFIABILITY_CLOSURE_PASS
IDENTIFIABILITY_CLOSURE_PARTIAL
IDENTIFIABILITY_CLOSURE_FAIL
```

### PASS

Requires meaningful held-out same-direction improvement, controlled leakage/null behavior, preserved broad integrity, and a characterized empirical identifiable regime.

### PARTIAL

Mechanism is understood and some leakage-safe improvement exists, but weak-child/near-boundary/other important regimes remain unresolved.

### FAIL

Use if p_both cannot improve without substantial leakage, gains disappear held-out, parent integrity deteriorates, or no stable identifiable regime exists.

FAIL is scientifically informative and may justify abandoning the current B_hierarchy formulation.

## 23. Interpretation of possible outcomes

1. **Parent absorption dominates** → redesign parent estimation/orthogonalization.
2. **Residual retains child but detector fails** → improve local detection.
3. **Reconciliation removes children** → redesign reconciliation.
4. **Recovery depends mainly on SNR/width ratio** → define an empirical identifiability regime and allow ambiguous outputs.
5. **Recovery/leakage trade-off is irreducible** → reconsider the B_hierarchy formulation rather than tuning indefinitely.

## 24. Final questions

The final report must answer:

1. Where do same-direction children disappear?
2. How much local signal is absorbed by the parent?
3. Does absorption increase with B?
4. Does it decrease with L?
5. Does width ratio predict absorption?
6. Does boundary proximity predict failure?
7. Does finer resolution help?
8. Is absolute L, L/B, or local SNR more predictive?
9. Does replicate variance shift the identifiability boundary?
10. Are child candidates generated then removed?
11. Does reconciliation materially reduce p_both?
12. Can parent estimation change without damaging broad boundaries?
13. Which child detector performs best on identical residuals?
14. Can p_both improve while leakage remains controlled?
15. What is the nested-recovery/leakage Pareto frontier?
16. What happens under broad-only truth?
17. What happens under heteroscedastic null?
18. What happens near parent boundaries?
19. Does improvement persist across child widths?
20. Does it persist at moderate/high noise?
21. Does it persist held-out?
22. Is there a stable identifiable regime?
23. Is there an ambiguous regime where the method should abstain?
24. Is |L|≈1.2–1.8 reproducible across conditions?
25. If not, what better describes the threshold?
26. Does B_hierarchy remain the best leakage-safe Path B reference?
27. Is the hierarchy worth further development?
28. Should the next step be production freeze, model redesign, or abandonment?

## 25. Execution order

```text
freeze v0.8 B_hierarchy reference
→ instrument stagewise traces
→ B×L grid
→ width-ratio benchmark
→ child-position/boundary benchmark
→ resolution benchmark
→ noise/replicate benchmark
→ parent-absorption analysis
→ controlled parent sensitivity
→ controlled child-detector sensitivity
→ reconciliation ablation
→ leakage/null guardrails
→ select at most one improvement
→ freeze configuration
→ held-out validation
→ regression/BigWig check
→ final closure decision
```

## 26. Agent instruction

Execute **EpiSegDiff v0.8.1 Same-Direction Nested Identifiability Closure**.

Do not redesign the entire project and do not begin v0.9.

Use `B_hierarchy` from v0.8 as the frozen reference. Determine why same-direction nested children are poorly recovered and whether recovery can improve without sacrificing:

```text
low broad→local leakage
low false-child rate
strong broad-domain recovery
acceptable null behavior
BigWig compatibility
```

Instrument loss at parent absorption, child detection, thresholding, and reconciliation.

Systematically study B, L, L/B, local SNR, parent/child widths, width ratio, child position, boundary distance, resolution, noise, and replicate variance.

Do not optimize p_both without simultaneous leakage/null measurement. Freeze any proposed improvement before held-out testing.

At completion generate:

```text
reports/v0.8.1_final_closure_decision.md
```

Then stop for formal review.

**No production architecture freeze and no v0.9 without explicit authorization.**
