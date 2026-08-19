# EpiSegDiff v0.8 --- Path B Hierarchical / Multiresolution Architecture Specification

**Status:** Agent-executable architecture and methodological research
specification\
**Authorization basis:** v0.7 strict methodological review + v0.7.1
R1--R5 closure + v0.8 architecture decision gate\
**Scientific direction:** **Path B --- hierarchical / multiresolution
redesign**\
**Path A status:** Rejected as the primary scientific plan; existing
engines remain controls/sub-solvers\
**Release objective:** Define, prototype, falsify, and compare
statistically identifiable multiscale formulations before any production
architecture is frozen.

------------------------------------------------------------------------

# 1. Decision basis

The architecture decision gate established all five redesign criteria:

1.  no single tested engine was acceptable across narrow and broad
    regimes;
2.  nested events remained systematically merged after resolution
    correction;
3.  improved same-direction separation traded against
    fragmentation/missed detection;
4.  one-scale routing was unstable on held-out simulations;
5.  gradual boundaries remained structurally difficult under
    point-change formulations.

Therefore v0.8 must not be another PELT/fused/WBS tuning release.

The central problem is now:

> **What statistical representation of genomic differential signal can
> simultaneously preserve broad domains, local narrow changes, nested
> structure, accurate boundaries, and valid inference?**

------------------------------------------------------------------------

# 2. Core hypothesis

The flat model

``` text
one differential score
        ↓
one segmentation
        ↓
one non-overlapping region set
```

may be insufficient when differential structure exists simultaneously at
multiple genomic scales.

A primary candidate representation is:

\[ D(x)=B(x)+L(x)+`\epsilon`{=tex}(x) \]

where:

-   (D(x)): observed/estimated differential signal;
-   (B(x)): broad/low-frequency differential component;
-   (L(x)): local/high-frequency differential component;
-   (`\epsilon`{=tex}(x)): residual noise.

This equation is a **candidate formulation**, not an accepted model.

v0.8 must test whether this or an equivalent
hierarchical/multiresolution representation is identifiable, stable,
calibrated, and empirically superior to flat controls.

------------------------------------------------------------------------

# 3. Non-negotiable design principle: identifiability first

Without constraints,

\[ D=B_1+L_1=B_2+L_2 \]

can have infinitely many decompositions.

Therefore no additive decomposition may be accepted without explicit
constraints defining what distinguishes broad from local structure.

Candidate constraints may involve:

``` text
smoothness
total variation
minimum/typical duration
frequency/scale
sparsity
hierarchical support
penalty structure
```

For every proposed model document:

1.  mathematical objective;
2.  assumptions;
3.  identifiability mechanism;
4.  tuning parameters;
5.  failure cases;
6.  computational complexity;
7.  inferential implications.

------------------------------------------------------------------------

# 4. Architectures that must be compared

Do not select one architecture before benchmarking.

## A. Additive broad + local decomposition

Concept:

\[ D(x)=B(x)+L(x)+`\epsilon`{=tex}(x) \]

Candidate penalized form:

\[ `\min`{=tex}\_{B,L} `\mathcal{L}`{=tex}(D,B+L) +`\lambda`{=tex}\_B
P_B(B) +`\lambda`{=tex}\_L P_L(L) \]

Possible choices:

-   smooth/trend-filtered (B) + sparse/fused (L);
-   low-frequency (B) + sparse local (L);
-   broad piecewise-constant (B) + local residual segmentation.

### Required question

Can the model recover an embedded local event without inventing local
residuals inside a broad-only truth?

------------------------------------------------------------------------

## B. Hierarchical segmentation tree

Concept:

``` text
chromosome
└── broad differential domain
    ├── subdomain
    │   └── local event
    └── local event
```

Regions are allowed to have parent-child relationships.

Required outputs should include:

``` text
region_id
parent_id
level
start
end
direction
effect
scale
support
```

### Required question

Does hierarchy recover true nested structure without generating
arbitrary deep trees?

------------------------------------------------------------------------

## C. Explicit multiresolution decomposition

Concept:

\[ D(x)=D\_{coarse}(x)+D\_{medium}(x)+D\_{fine}(x)+`\epsilon`{=tex}(x)
\]

Candidate families may include:

-   wavelet/multiresolution decomposition;
-   scale-space decomposition;
-   multiresolution fused lasso;
-   trend filtering;
-   multiscale likelihood methods.

Wavelets are **not automatically preferred**.

### Required question

Are scale-specific components stable to modest changes in resolution and
tuning?

------------------------------------------------------------------------

## D. Multiscale penalized likelihood

Develop at least one joint objective where multiple scales are estimated
simultaneously rather than routed after detection.

Conceptually:

\[ `\mathcal{J}`{=tex} = -`\log `{=tex}p(Y`\mid `{=tex}B,L)
+`\lambda`{=tex}\_B P_B(B) +`\lambda`{=tex}\_L P_L(L) \]

or a justified extension.

This candidate is especially important if it provides a more explicit
statistical model than sequential residualization.

------------------------------------------------------------------------

# 5. Existing methods remain mandatory controls

Do not delete prior engines.

Retain:

``` text
raw-effect + PELT
standardized + PELT
fused lasso
window aggregation
candidate union
```

Where practical retain WBS as an experimental control, but do not devote
major development effort to it unless implementation validation
justifies doing so.

These controls answer:

> Does Path B actually improve the failures that motivated redesign?

------------------------------------------------------------------------

# 6. Observation resolution

Separate biological event width from observation resolution.

Required base resolutions:

``` text
50 bp
100 bp
200 bp
```

Optional:

``` text
25 bp
500 bp
1 kb
```

where computationally and biologically sensible.

Do not infer that a 200 bp event is intrinsically undetectable merely
because a 200 bp bin representation fails.

------------------------------------------------------------------------

# 7. Broad/local scale definition

Do not hard-code an arbitrary biological threshold such as:

``` text
<1 kb = local
>1 kb = broad
```

unless it is explicitly a benchmark definition.

Architecture candidates must instead define scale through their
objective or multiresolution representation.

For reporting, descriptive width categories may be used:

``` text
narrow
intermediate
broad
```

but these labels must not silently determine the algorithm.

------------------------------------------------------------------------

# 8. Mandatory falsification suite

Every Path B prototype must be tested against all ten falsifiers before
promotion.

## F1 --- Isolated narrow event

Truth:

``` text
100–500 bp local differential event
```

Must recover useful detection and boundaries.

## F2 --- Isolated broad event

Truth:

``` text
5–100 kb broad differential domain
```

Must preserve broad coherence.

## F3 --- Nested same-direction event

Example:

``` text
broad +1
embedded local additional +1
```

Must recover broad baseline and local enhancement separately when the
local effect is identifiable.

## F4 --- Nested opposite-direction event

Example:

``` text
broad +1
embedded local -2
```

Must recover local reversal.

## F5 --- Same-direction separated events

Two events with a true neutral gap must not be systematically merged.

## F6 --- Broad gradual transition

Truth includes explicit transition intervals rather than only point
boundaries.

## F7 --- Narrow event at broad boundary

Must not generate duplicate/artifactual events from residualization or
edge effects.

## F8 --- Broad-only null for local component

Truth:

\[ L(x)=0 \]

Method must not invent local structure.

## F9 --- Heteroscedastic null

Local variance changes without differential mean change.

Method must remain controlled.

## F10 --- Replicate perturbation / real-data stability

Calls and hierarchical relationships must remain reasonably stable under
replicate leave-one-out/bootstrap/parameter perturbation.

------------------------------------------------------------------------

# 9. Expand falsification to component leakage

Path B introduces new error modes.

Measure:

## Broad-to-local leakage

Truth:

\[ D=B,`\quad `{=tex}L=0 \]

but:

\[ `\hat `{=tex}L`\neq0`{=tex} \]

## Local-to-broad leakage

Truth:

\[ D=L,`\quad `{=tex}B=0 \]

but:

\[ `\hat `{=tex}B`\neq0`{=tex} \]

## Scale misassignment

Correct region detected but assigned to the wrong scale/component.

## False hierarchy

Correct parent event but false child, or correct child but false parent.

## Double counting

Same differential effect reported simultaneously as broad and local
without a justified decomposition.

These metrics are mandatory.

------------------------------------------------------------------------

# 10. Synthetic truth families

Use at minimum:

``` text
isolated narrow
isolated intermediate
isolated broad
same-direction adjacent
opposite-direction adjacent
nested same-direction
nested attenuation
nested reversal
peak-on-domain
multi-peak broad domain
gradual broad transition
boundary-local event
broad-only
local-only
global null
heteroscedastic null
autocorrelated null
```

------------------------------------------------------------------------

# 11. Core width grid

Use approximately:

``` text
100 bp
200 bp
300 bp
500 bp
750 bp
1 kb
2 kb
5 kb
10 kb
25 kb
50 kb
100 kb
```

Nested simulations must independently vary:

``` text
broad width
local width
broad effect B
local effect L
local position inside broad domain
```

------------------------------------------------------------------------

# 12. Effect/noise grid

Do not select architecture using only an easy effect/noise regime.

At minimum:

``` text
effect:
weak
moderate
strong

noise:
low
moderate
high
```

Numerical values must be frozen and documented.

Use at least 50 replicates for primary benchmark cells; prefer 100 where
feasible.

Store seeds.

------------------------------------------------------------------------

# 13. Nested B × L matrix

This is a critical v0.8 benchmark.

For each broad effect (B), vary local effect (L):

``` text
same-direction enhancement
weak attenuation
strong attenuation
local reversal
```

Measure:

``` text
P(parent detected)
P(child detected)
P(parent + child both correct)
P(merged)
P(child-only)
P(parent-only)
false-child rate
boundary error
component leakage
```

Do not summarize nested performance with union F1 alone.

------------------------------------------------------------------------

# 14. Transition-zone truth

For gradual events define:

``` text
stable background
left transition interval
core differential domain
right transition interval
stable background
```

Metrics:

``` text
core-domain recall
transition-zone overlap
transition midpoint error
over-extension
under-extension
false sharp-boundary confidence
```

A method must not be rewarded for pretending a gradual transition has an
exact biological breakpoint.

------------------------------------------------------------------------

# 15. Hierarchical truth metrics

If the model emits parent/child structure, evaluate:

``` text
parent precision/recall
child precision/recall
parent-child relationship accuracy
tree edit distance or simpler declared hierarchy loss
depth overfitting
orphan child rate
false child per true parent
```

Prefer simple interpretable metrics before complex tree metrics.

------------------------------------------------------------------------

# 16. Detection and boundary must remain separate

Always report:

\[ P(`\text{detected}`{=tex}) \]

and:

\[ E\[`\text{boundary error}`{=tex}`\mid`{=tex}`\text{detected}`{=tex}\]
\]

plus unconditional reconstruction metrics.

Do not report boundary error only among successful detections without
detection probability.

------------------------------------------------------------------------

# 17. Inference architecture must be explicit

For every candidate separate:

``` text
1. candidate discovery
2. component/scale assignment
3. boundary estimation
4. effect estimation
5. hypothesis testing
6. multiplicity correction
```

No report may use "FDR" without specifying which hypothesis family is
controlled.

------------------------------------------------------------------------

# 18. Inferential targets to evaluate

At least conceptually distinguish:

## Region-level FDR

\[ E`\left[\frac{V}{R\vee1}\right]`{=tex}\]

over selected reported regions.

## Component-level FDR

False broad/local components among reported components.

## Base-pair error

False differential bp / reported differential bp.

## Family-wise probability

\[ P(`\text{any false region}`{=tex}) \]

## Hierarchical error

Examples:

``` text
false child conditional on parent
false hierarchy edge
false scale assignment
```

The production target must be chosen explicitly before final acceptance.

------------------------------------------------------------------------

# 19. Selection-bias strategies

Re-evaluate inference only after the Path B discovery model is defined.

Candidate approaches:

``` text
same-data descriptive baseline
cross-fitting
selective-inference approximation
full-pipeline permutation
bootstrap stability/support
hierarchical testing
```

Naive sample splitting should remain a control because v0.7 showed
severe power loss.

Do not assume cross-fitting is calibrated until tested.

------------------------------------------------------------------------

# 20. Null benchmark

Required null regimes:

``` text
global null
low-signal null
heteroscedastic null
autocorrelated null
broad-only truth when testing local component
local-only truth when testing broad component
```

Nominal levels where inferential p/FDR procedures exist:

``` text
0.01
0.05
0.10
```

Report:

``` text
false regions
false components
false bp
P(any false)
empirical FDP
component leakage
false hierarchy
```

------------------------------------------------------------------------

# 21. Identifiability stress tests

Explicitly test cases where broad and local components are difficult to
distinguish.

Examples:

``` text
local width approaches broad width
broad effect becomes weak
local effect becomes weak
local event near broad edge
two local events mimic a broad plateau
dense local events mimic broad change
gradual broad transition resembles multiple local changes
```

For each candidate architecture identify the empirical region where
decomposition becomes unreliable.

The method must be allowed to return:

``` text
AMBIGUOUS_SCALE
```

or an equivalent uncertainty state if scale is not identifiable.

Do not force false certainty.

------------------------------------------------------------------------

# 22. Model complexity control

Hierarchical/multiresolution models can overfit.

Track:

``` text
number of components
number of hierarchy levels
number of child events
effective degrees of freedom where available
runtime
memory
```

A more complex model must show measurable improvement on held-out truth.

------------------------------------------------------------------------

# 23. Training / tuning / test separation

Partition synthetic designs/seeds into:

``` text
development
validation
held-out test
```

Architecture/penalty choices may use development and validation sets.

Final headline metrics must come from held-out simulations.

Do not tune on the final test set.

------------------------------------------------------------------------

# 24. Model-selection strategy

Do not choose a winner from one F1 value.

Use Pareto analysis over:

``` text
narrow detection
broad detection
nested recovery
boundary accuracy
merge rate
split rate
component leakage
null false positives
hierarchy accuracy
runtime
memory
stability
```

Prefer the simplest Pareto-optimal architecture.

------------------------------------------------------------------------

# 25. Required prototype stages

## Stage 1 --- Mathematical specification

No major production coding.

For each candidate:

``` text
objective
constraints
identifiability
optimization
parameters
expected failure modes
inference compatibility
```

## Stage 2 --- Minimal synthetic prototype

Implement only enough to test F1--F9.

## Stage 3 --- Comparative falsification

Compare A/B/C/D architectures against flat controls.

## Stage 4 --- Select at most two candidates

Only candidates passing minimum falsification gates continue.

## Stage 5 --- Expanded held-out benchmark

Run the full width/effect/noise/nested grid.

## Stage 6 --- Real-data feasibility

Only after synthetic evidence supports the architecture.

## Stage 7 --- Inference/calibration development

Formalize testing and multiplicity.

## Stage 8 --- Production decision

Freeze architecture only after all gates.

------------------------------------------------------------------------

# 26. Real-data validation

After synthetic selection, evaluate at least:

## Narrow

Replicated TF-like dataset, including JUN where appropriate.

## Intermediate

Replicated ATAC/open-chromatin dataset.

## Broad

Replicated broad histone-mark dataset such as H3K27me3.

Where raw BAM exists compare against:

``` text
csaw
THOR where compatible
fixed-peak workflow
window aggregation
```

Do not claim superiority if external comparator data remain unavailable.

------------------------------------------------------------------------

# 27. BigWig input remains a required production path

Path B must preserve continuous-signal input.

Required considerations:

``` text
replicate BigWigs
common chromosome naming
common genome build
missing coverage
binning/resampling
linear vs transformed scale
normalization assumptions
```

Do not design the new model around BAM-only inputs.

Where counts are available, count-based workflows may provide additional
inference options, but BigWig must remain supported.

------------------------------------------------------------------------

# 28. Output model

Path B may require richer output than a flat BED.

Minimum flat-compatible fields:

``` text
chrom
start
end
direction
effect
qvalue_or_support
scale
component_type
region_id
parent_id
level
boundary_support
```

If hierarchical output is adopted, also export a flattened BED/TSV
representation for compatibility.

------------------------------------------------------------------------

# 29. Required tables

``` text
results/v0.8/candidate_architectures.tsv
results/v0.8/identifiability_tests.tsv
results/v0.8/falsification_summary.tsv
results/v0.8/narrow_metrics.tsv
results/v0.8/broad_metrics.tsv
results/v0.8/nested_metrics.tsv
results/v0.8/transition_metrics.tsv
results/v0.8/component_leakage.tsv
results/v0.8/hierarchy_metrics.tsv
results/v0.8/null_calibration.tsv
results/v0.8/heldout_benchmark.tsv
results/v0.8/pareto_front.tsv
results/v0.8/runtime_memory.tsv
results/v0.8/realdata_stability.tsv
results/v0.8/external_comparators.tsv
results/v0.8/regression_summary.tsv
```

Failed configurations remain in tables with explicit reasons.

------------------------------------------------------------------------

# 30. Required figures

``` text
figures/v0.8_architecture_schematic.*
figures/v0.8_identifiability_map.*
figures/v0.8_nested_BxL_heatmap.*
figures/v0.8_parent_child_recovery.*
figures/v0.8_component_leakage.*
figures/v0.8_width_effect_phase_diagram.*
figures/v0.8_transition_zone_recovery.*
figures/v0.8_null_error.*
figures/v0.8_pareto_front.*
figures/v0.8_heldout_performance.*
figures/v0.8_realdata_examples.*
```

------------------------------------------------------------------------

# 31. Required reports

``` text
reports/v0.8_candidate_model_specifications.md
reports/v0.8_identifiability_report.md
reports/v0.8_minimal_prototype_falsification.md
reports/v0.8_nested_structure_report.md
reports/v0.8_transition_boundary_report.md
reports/v0.8_component_leakage_report.md
reports/v0.8_inference_target_report.md
reports/v0.8_null_calibration_report.md
reports/v0.8_heldout_benchmark.md
reports/v0.8_realdata_validation.md
reports/v0.8_external_comparison.md
reports/v0.8_final_architecture_acceptance.md
```

Every report must distinguish:

``` text
FACT
INFERENCE
SPECULATION
UNRESOLVED
```

------------------------------------------------------------------------

# 32. Architecture gates

## Gate A --- Identifiability

Candidate has explicit broad/local/scale constraints and does not
collapse into arbitrary decomposition.

## Gate B --- Narrow recovery

Useful 100--500 bp recovery.

## Gate C --- Broad recovery

Retains coherent 5--100 kb domains.

## Gate D --- Nested recovery --- CRITICAL

Substantial improvement over flat PELT/fused/union on same-direction and
opposite-direction nested truth.

## Gate E --- Separation/coherence

Improved local separation must not cause unacceptable broad
fragmentation.

## Gate F --- Transition representation

Gradual events are represented without pretending unsupported exact
breakpoints.

## Gate G --- Component leakage --- CRITICAL

Broad-only truth does not routinely generate local calls and vice versa.

## Gate H --- Null behavior --- CRITICAL

No material uncontrolled false structure under required nulls.

## Gate I --- Held-out generalization

Benefits persist on held-out simulations.

## Gate J --- Complexity justification

Added complexity produces measurable gain.

## Gate K --- Real-data stability

Calls/components/hierarchy are reasonably stable.

## Gate L --- Inference declaration

Inferential target and multiplicity family are explicit.

## Gate M --- BigWig compatibility

New architecture works through the continuous-signal input path.

## Gate N --- Regression

Validated v0.4--v0.7 infrastructure remains intact.

------------------------------------------------------------------------

# 33. Minimum success criteria for Path B

A candidate cannot be accepted merely because nested F1 improves.

It must jointly demonstrate:

``` text
narrow recovery
broad coherence
nested recovery
acceptable separation
low component leakage
reasonable transition behavior
held-out robustness
acceptable null behavior
```

If no candidate achieves this, v0.8 must FAIL rather than forcing an
architecture.

------------------------------------------------------------------------

# 34. Stop conditions

Stop a candidate early if:

-   decomposition is non-identifiable in ordinary benchmark regimes;
-   local component is routinely hallucinated in broad-only truth;
-   broad component absorbs isolated local truth;
-   nested gain comes only from massive overcalling;
-   broad fragmentation becomes unacceptable;
-   held-out gains disappear;
-   runtime/memory becomes impractical without compensating benefit.

Negative results must remain documented.

------------------------------------------------------------------------

# 35. What v0.8 must not claim

Until demonstrated:

``` text
universal scale adaptivity
formal region-level FDR control
validated hierarchical FDR
superiority to csaw
superiority to THOR
biological novelty of nested calls
exact biological boundaries for gradual transitions
```

------------------------------------------------------------------------

# 36. What v0.8 is allowed to establish

If supported by data:

-   flat segmentation is insufficient for specified nested regimes;
-   a particular multiresolution formulation improves nested recovery;
-   broad/local decomposition is identifiable within defined regimes;
-   component leakage is controlled within tested conditions;
-   hierarchical output better represents specified synthetic truth;
-   a Path B architecture improves the narrow/broad Pareto frontier.

Claims must remain conditional on tested regimes.

------------------------------------------------------------------------

# 37. Final decision states

Allowed:

``` text
PATH_B_PASS
PATH_B_CONDITIONAL_PASS
PATH_B_FAIL
```

## PATH_B_PASS

Requires all critical gates, including nested recovery, leakage, null
behavior, held-out generalization, and BigWig compatibility.

## PATH_B_CONDITIONAL_PASS

Architecture is promising but one non-fatal validation block remains.

Do not use this state if identifiability, leakage, or null behavior
fails.

## PATH_B_FAIL

Use if no candidate provides a defensible improvement over flat
controls.

Failure is scientifically acceptable.

------------------------------------------------------------------------

# 38. Required final questions

The final architecture report must answer:

1.  Is additive broad+local decomposition identifiable?
2.  Is hierarchical segmentation more stable than additive
    decomposition?
3.  Does explicit multiresolution modeling outperform simple
    residualization?
4.  Which architecture best recovers isolated narrow events?
5.  Which best preserves isolated broad domains?
6.  Which best recovers same-direction nested events?
7.  Which best recovers nested reversals?
8.  At what B×L combinations is nested structure identifiable?
9.  At what width ratios does scale assignment become ambiguous?
10. How often does broad signal leak into the local component?
11. How often does local signal leak into the broad component?
12. Does the model create false children?
13. Does hierarchy become unnecessarily deep?
14. Does improved separation fragment broad domains?
15. How are gradual transitions represented?
16. Are transition zones more accurate than point breakpoints?
17. Does performance persist at 50/100/200 bp observation resolution?
18. Does performance persist across effect/noise regimes?
19. Does performance persist on held-out simulations?
20. What happens under global null?
21. What happens under heteroscedastic null?
22. What happens under autocorrelated null?
23. What is the false local-component rate under broad-only truth?
24. What is the false broad-component rate under local-only truth?
25. What inferential target is selected?
26. Is region-level FDR actually calibrated?
27. If not, what support statistic is reported instead?
28. Are real JUN results stable?
29. Are real ATAC results stable?
30. Are broad histone results stable?
31. Does BigWig input behave consistently across resolutions?
32. How does the selected architecture compare with raw+PELT?
33. How does it compare with fused?
34. How does it compare with window aggregation?
35. How does it compare with csaw where possible?
36. How does it compare with THOR where possible?
37. What complexity cost is introduced?
38. Which candidate architectures were rejected and why?
39. What failure regimes remain?
40. Is Path B PASS, CONDITIONAL PASS, or FAIL?

------------------------------------------------------------------------

# 39. Execution order

``` text
Phase 0
freeze v0.7 controls and benchmark definitions
    ↓
Phase 1
write mathematical candidate specifications
    ↓
Phase 2
implement minimal prototypes only
    ↓
Phase 3
run F1–F10 falsification suite
    ↓
Phase 4
reject clearly invalid candidates
    ↓
Phase 5
run identifiability + component leakage stress tests
    ↓
Phase 6
select at most two candidates
    ↓
Phase 7
full width × effect × noise × resolution benchmark
    ↓
Phase 8
nested B × L benchmark
    ↓
Phase 9
transition-zone benchmark
    ↓
Phase 10
held-out simulation evaluation
    ↓
Phase 11
define inference target
    ↓
Phase 12
null/calibration work
    ↓
Phase 13
BigWig compatibility and regression
    ↓
Phase 14
real narrow / ATAC / broad validation
    ↓
Phase 15
external comparator benchmark where available
    ↓
Phase 16
Pareto and complexity analysis
    ↓
Phase 17
final architecture acceptance
```

------------------------------------------------------------------------

# 40. Agent instruction

Implement **EpiSegDiff v0.8 Path B** as a staged methodological research
program.

Do **not** begin by committing the production codebase to wavelets,
residualization, hierarchical trees, HSMMs, or any other favored
technique.

First formalize and compare at least:

``` text
A. additive broad + local decomposition
B. hierarchical segmentation
C. explicit multiresolution decomposition
D. joint multiscale penalized likelihood
```

Existing flat methods remain mandatory controls.

The critical scientific problem is not simply better F1.

The candidate must demonstrate that it can represent:

``` text
broad structure
+
local structure
+
nested structure
```

without:

``` text
false local residuals
false broad components
excessive fragmentation
uncontrolled null calls
unstable scale assignment
```

Treat identifiability, component leakage, nested recovery, held-out
generalization, and null behavior as critical gates.

Maintain BigWig as a first-class input path.

Do not claim region-level FDR unless it is explicitly defined and
empirically calibrated.

At the end generate:

``` text
reports/v0.8_final_architecture_acceptance.md
```

Then stop for formal review.

Do not begin manuscript-superiority claims or v0.9 implementation
without review authorization.
