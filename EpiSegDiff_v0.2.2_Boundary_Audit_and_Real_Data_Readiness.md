# EpiSegDiff v0.2.2 --- Boundary Audit and Real-Data Readiness

**Purpose:** Resolve the directional boundary-error issue identified in
v0.2.1, fix the `NULL` configuration bug, and determine whether
EpiSegDiff Level I is ready for a scoped real-data pilot.

**Input evidence:** `v0.2.1_Final_Acceptance_Report.md`

**Scope:** Small validation release only.

**Do not implement NB-FPOP. Do not begin full competitor benchmarking.
Do not expand Level II topology features in this phase.**

------------------------------------------------------------------------

# 1. Motivation

v0.2.1 established that:

-   TMM is a reasonable default under sparse/moderate differential
    composition.
-   library-size normalization can generate false flanks.
-   strict-null calibration is approximately acceptable under the tested
    simulation regime.
-   Level I GAIN/LOSS detection is sufficiently promising for a
    real-data pilot.
-   topology definitions are mathematically cleaner after removal of
    `BIDIRECTIONAL_SHIFT`.
-   Level II topology inference remains unreliable because
    condition-specific boundary inference is weak.
-   boundary errors are strongly asymmetric.

Observed Task C boundary errors:

``` text
Narrow:
left  median/p90 = 350 / 350 bp
right median/p90 =   0 / 400 bp

Broad:
left  median/p90 = 7800 / 7800 bp
right median/p90 =    0 / 8000 bp
```

This asymmetry requires an implementation/algorithm audit before
interpreting the errors as an inherent statistical limitation.

------------------------------------------------------------------------

# 2. Primary Questions

This release must answer:

1.  Is the left/right boundary asymmetry caused by indexing,
    coordinates, search direction, initialization, or candidate
    anchoring?
2.  Is `end_walk` intrinsically directional?
3.  Does reversing the genomic signal reverse the error pattern as
    mathematically expected?
4.  How do bin size, search radius, padding, effect strength, and domain
    width affect boundary error?
5.  Can narrow-domain boundary accuracy be improved without changing the
    core statistical model?
6.  Is broad-domain failure an implementation problem, a
    local-boundary-method limitation, or unresolved?
7.  After the audit, is Level I ready for a real CUT&RUN/ChIP-seq pilot?

------------------------------------------------------------------------

# 3. Task A --- Fix the YAML NULL Bug

The v0.2.1 report states that the `NULL` class was dropped from Task E.

Do not use an unquoted YAML `NULL` label.

Replace the simulation truth label with:

``` text
NO_CHANGE
```

Preferred over `"NULL"` to avoid ambiguity across:

``` text
YAML null
JSON null
R NULL
Python None
NA/missing-value semantics
```

Update:

-   simulator;
-   configs;
-   truth tables;
-   benchmark scripts;
-   tests;
-   documentation.

## Acceptance

-   [ ] `NO_CHANGE` survives YAML read/write.
-   [ ] `NO_CHANGE` survives TSV export/import.
-   [ ] R and Python interoperability does not convert it to missing.
-   [ ] Task E includes the expected number of no-change simulations.
-   [ ] Regression test added.

------------------------------------------------------------------------

# 4. Task B --- Audit Coordinate and Index Conventions

EpiSegDiff uses BED-style:

\[ \[start,end) \]

Verify every boundary-related function uses the same convention.

Audit:

``` text
bin generation
candidate start/end
CVR boundaries
segment start/end
local search window
end_walk
BED export
truth simulator
boundary evaluation
```

For every module document:

``` text
coordinate origin
inclusive/exclusive convention
bin-index → genomic-coordinate conversion
```

## Required toy example

Use:

``` text
chromosome length = 1000 bp
bin size = 100 bp
true region = [300,700)
```

Expected bins:

``` text
[300,400)
[400,500)
[500,600)
[600,700)
```

Verify that all internal and exported representations return exactly:

``` text
start = 300
end   = 700
```

No ±1 bp or ±1 bin tolerance is allowed for this deterministic test.

------------------------------------------------------------------------

# 5. Task C --- Symmetric Step-Function Tests

Create deterministic/noiseless differential-effect signals.

## C1. Single positive block

Example:

``` text
0 0 0 0  +2 +2 +2 +2 +2  0 0 0 0
          ^               ^
          L               R
```

Test at multiple positions.

Expected:

\[ Error_L = Error_R = 0 \]

within bin resolution.

## C2. Negative block

Repeat with:

``` text
-2
```

to verify GAIN/LOSS sign does not affect boundary direction.

## C3. Multiple blocks

Test:

``` text
GAIN
LOSS
GAIN → LOSS
LOSS → GAIN
```

## Acceptance

If noiseless symmetric signals produce systematic left/right asymmetry,
classify the release as **FAIL** until fixed.

------------------------------------------------------------------------

# 6. Task D --- Reverse-Signal Invariance Test

This is a mandatory diagnostic.

For a signal of length (N):

\[ S=(s_1,s_2,`\ldots`{=tex},s_N) \]

construct:

\[ S\^{rev}=(s_N,`\ldots`{=tex},s_2,s_1) \]

Run the complete boundary-refinement procedure on both.

For a forward interval:

\[ \[L,R) \]

the expected reversed interval is:

\[ \[N-R,N-L) \]

after coordinate transformation.

Boundary errors should transform approximately as:

\[ Error_L\^{forward} `\approx`{=tex} Error_R\^{reverse} \]

and:

\[ Error_R\^{forward} `\approx`{=tex} Error_L\^{reverse} \]

within numerical/bin tolerance.

## Required output

``` text
boundary_reverse_invariance.tsv
```

Columns:

``` text
simulation_id
orientation
true_left
true_right
estimated_left
estimated_right
left_error
right_error
mapped_reverse_left_error
mapped_reverse_right_error
invariance_difference
```

## Acceptance

Systematic orientation-dependent discrepancies indicate algorithmic bias
and must be fixed before real-data pilot.

------------------------------------------------------------------------

# 7. Task E --- Audit `end_walk`

Inspect the current implementation specifically for:

-   left/right loop asymmetry;
-   different stopping conditions;
-   off-by-one indexing;
-   initialization from candidate start versus end;
-   threshold comparison asymmetry;
-   use of `>` versus `>=`;
-   bin-center versus bin-edge conversion;
-   padding clipping;
-   chromosome/CVR clipping;
-   candidate anchoring;
-   search-radius truncation.

Create:

``` text
docs/end_walk_algorithm.md
```

Document pseudocode separately for left and right searches.

If the algorithms are intended to be symmetric, demonstrate their
symmetry explicitly.

------------------------------------------------------------------------

# 8. Task F --- Boundary Parameter Stress Test

After fixing any implementation issue, test:

## Bin size

``` text
10 bp
25 bp
50 bp
100 bp
250 bp
```

Use only scales computationally appropriate for each simulated domain.

## Candidate padding

Example grid:

``` text
0
1× region width
2× region width
5× region width
```

or an equivalent documented grid.

## Search radius

Test multiple radii relative to expected domain width.

## Effect strength

``` text
weak
moderate
strong
```

## Replicates

``` text
2
3
5
```

## Region widths

### Narrow

``` text
200 bp
500 bp
1 kb
2 kb
```

### Broad

``` text
10 kb
20 kb
50 kb
100 kb
```

------------------------------------------------------------------------

# 9. Required Boundary Metrics

Report separately for left and right boundaries:

\[ MAE_L=`\operatorname{median}`{=tex}\|`\hat `{=tex}L-L\| \]

\[ MAE_R=`\operatorname{median}`{=tex}\|`\hat `{=tex}R-R\| \]

Also report:

``` text
median error
mean error
90th percentile
95th percentile
maximum error
signed error
Hausdorff distance
```

The signed error is important for detecting systematic inward/outward
bias.

Define:

\[ SignedError_L=`\hat `{=tex}L-L \]

\[ SignedError_R=`\hat `{=tex}R-R \]

------------------------------------------------------------------------

# 10. Boundary Asymmetry Metric

Define:

\[ A= MAE_L-MAE_R \]

and optionally normalized:

\[ A\_{norm} = `\frac{MAE_L-MAE_R}`{=tex} {MAE_L+MAE_R+`\epsilon`{=tex}}
\]

Report by:

``` text
bin size
region width
effect strength
replicate number
boundary method
orientation
```

The goal is not to force (A=0), but to determine whether observed
asymmetry is explained by noise or by implementation directionality.

------------------------------------------------------------------------

# 11. Task G --- Compare Boundary Refinement Methods

Compare current:

``` text
end_walk
```

against at least one alternative already compatible with the framework.

Preferred candidates:

``` text
local differential-gradient maximum
local likelihood refinement
cross-scale consensus
```

Do not introduce a large new modeling framework.

Compare:

``` text
left MAE
right MAE
p90 error
Hausdorff distance
orientation invariance
runtime
```

------------------------------------------------------------------------

# 12. Task H --- Broad-Domain Failure Analysis

Do not merely report that broad boundaries are inaccurate.

Determine whether the error originates from:

``` text
initial differential segmentation
candidate interval truncation
normalization
local boundary refinement
bin resolution
search radius
broad low-gradient signal
```

Perform an oracle experiment.

## Oracle candidate test

Provide the boundary-refinement method with a candidate centered around
the known true domain.

If error remains several kb:

> boundary refinement itself is the likely bottleneck.

If error collapses:

> upstream candidate discovery/anchoring is the likely bottleneck.

## Oracle normalization test

Use the true simulated normalization factors where available.

If boundary accuracy remains poor:

> normalization is not the primary explanation.

------------------------------------------------------------------------

# 13. Task I --- Level I Isolation Test

Because Level I is proposed as the primary real-data output, verify that
it does not depend on Level II success.

Run the pipeline with topology disabled:

``` text
differential effect
    ↓
segmentation
    ↓
candidate recount
    ↓
region testing
    ↓
GAIN / LOSS
```

Compare with the full pipeline.

Required check:

``` text
Level I region calls identical before/after disabling topology
```

except for explicitly documented post-processing dependencies.

If topology modifies Level I calls, document and remove unnecessary
coupling.

------------------------------------------------------------------------

# 14. Real-Data Readiness Criteria

A real-data pilot may begin if:

## Mandatory

-   [ ] YAML `NO_CHANGE` bug fixed.
-   [ ] deterministic coordinate tests pass.
-   [ ] noiseless step tests show no directional implementation bias.
-   [ ] reverse-signal invariance passes within bin tolerance.
-   [ ] `end_walk` algorithm is documented.
-   [ ] Level I is independent of topology.
-   [ ] no unresolved indexing/coordinate bug remains.

## Narrow-domain requirement

Boundary refinement should show substantially improved or explained
narrow-domain errors.

Do not impose an arbitrary publication threshold yet.

The report must provide actual:

``` text
median
p90
signed bias
```

and decide whether the remaining error is acceptable relative to bin
size and biological feature width.

## Broad-domain requirement

Broad-domain Level II does **not** need to pass for the narrow real-data
pilot.

If broad boundaries remain inaccurate, explicitly keep:

``` text
broad Level II = EXPERIMENTAL / BLOCKED
```

------------------------------------------------------------------------

# 15. Real-Data Pilot Scope if PASS

If v0.2.2 passes, the next release may analyze exactly one initial
cohort.

Preferred:

``` text
paired-end CUT&RUN or ChIP-seq
biological replicates >= 2 per condition
narrow/local signal
```

Examples:

``` text
TF CUT&RUN
TF ChIP-seq
H3K27ac
H3K4me3
```

Do not start with H3K27me3/H3K9me3 as the first pilot.

Primary output:

``` text
Level I GAIN / LOSS
```

Secondary experimental output:

``` text
Level II boundary dynamics
```

Normalization:

``` text
TMM = primary
background = sensitivity analysis
```

------------------------------------------------------------------------

# 16. Required Files

Create/update:

``` text
docs/end_walk_algorithm.md

scripts/test_coordinate_conventions.R
scripts/test_reverse_boundary_invariance.R
scripts/run_boundary_parameter_stress.R
scripts/run_boundary_method_comparison.R
scripts/run_broad_boundary_failure_analysis.R

config/validation/boundary_audit.yaml
```

------------------------------------------------------------------------

# 17. Required Result Tables

``` text
boundary_coordinate_tests.tsv
boundary_reverse_invariance.tsv
boundary_parameter_stress.tsv
boundary_method_comparison.tsv
broad_boundary_oracle_tests.tsv
level1_independence_test.tsv
```

------------------------------------------------------------------------

# 18. Required Reports

Generate:

``` text
reports/v0.2.2_Coordinate_and_Direction_Audit_Report.md
reports/v0.2.2_Boundary_Method_Comparison_Report.md
reports/v0.2.2_Broad_Domain_Failure_Report.md
reports/v0.2.2_Level1_Independence_Report.md
reports/v0.2.2_Real_Data_Readiness_Report.md
```

Every report must distinguish:

``` text
FACT
INFERENCE
UNRESOLVED
```

and contain:

``` text
Objective
Methods
Parameters
Results
Failures
Interpretation
Limitations
Decision
Next action
```

------------------------------------------------------------------------

# 19. Final Decision Rules

Allowed decisions:

``` text
PASS
CONDITIONAL PASS
FAIL
```

## PASS

Proceed to one narrow real-data pilot.

## CONDITIONAL PASS

Proceed only with:

``` text
Level I as primary
Level II explicitly experimental
```

and document the unresolved boundary limitation.

## FAIL

Do not begin real-data pilot if:

-   directional implementation bias remains;
-   coordinate conventions remain inconsistent;
-   reverse invariance fails unexplained;
-   Level I depends incorrectly on topology;
-   simulation truth cannot be reconstructed in noiseless tests.

------------------------------------------------------------------------

# 20. NB-FPOP Decision

NB-FPOP remains blocked.

Do not interpret broad-domain boundary failure as evidence that NB-FPOP
is required.

Only reconsider NB-FPOP if later evidence identifies a specific failure
caused by the observation/segmentation loss rather than:

``` text
normalization
candidate anchoring
boundary refinement
bin resolution
topology inference
```

------------------------------------------------------------------------

# 21. Agent Start Instruction

Execute **EpiSegDiff v0.2.2 Boundary Audit and Real-Data Readiness**
only.

Start by replacing the simulation label `NULL` with `NO_CHANGE`.

Then audit coordinate conventions and run deterministic symmetric
step-function tests.

Next run the mandatory reverse-signal invariance experiment before
changing the boundary algorithm.

If directional bias is found, identify and fix its cause, add a
regression test, and rerun the complete audit.

Afterward compare `end_walk` with at least one simple alternative
boundary-refinement method and perform the broad-domain oracle
experiments.

Finally verify Level I independence from Level II and generate:

``` text
reports/v0.2.2_Real_Data_Readiness_Report.md
```

Stop after this report.

Do not run real biological data, competitor benchmarking, or NB-FPOP
during v0.2.2.
