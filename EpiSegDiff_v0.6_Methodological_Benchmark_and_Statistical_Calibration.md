# EpiSegDiff v0.6 --- Methodological Benchmark and Statistical Calibration

**Status:** Agent-executable specification\
**Baseline:** v0.5 Scientific Validation --- CONDITIONAL PASS\
**Purpose:** Determine whether change-point-based differential
segmentation provides a reproducible methodological advantage over
established differential epigenomic approaches, and whether EpiSegDiff
region-level inference is statistically calibrated.

## 1. Why v0.6

v0.5 established that EpiSegDiff is technically valid and scientifically
non-trivial:

-   adaptive PELT segmentation differs from fixed-bin testing;
-   synthetic narrow/broad/boundary-crossing events were recovered
    accurately;
-   GSE280947 EpiSegDiff and fixed-peak limma effects were strongly
    directionally concordant;
-   label scrambling collapsed the signal;
-   EpiSegDiff-only regions showed modest motif support above matched
    background;
-   however, superiority over established de novo differential-region
    methods was not established.

The key question is now:

> Does
> `differential evidence → change-point segmentation → adaptive domain`
> provide measurable advantages over fixed peaks, significant-window
> aggregation, and HMM/state-based differential-region methods?

A valid result may show superiority only for particular event
geometries, equivalence, a sensitivity/domain-coherence tradeoff, or a
statistical weakness. Do not redesign the benchmark merely because
EpiSegDiff does not win.

------------------------------------------------------------------------

## 2. Working method landscape

``` text
Fixed peak / DiffBind-like
predefined regions → counts/signal → differential test

csaw
sliding-window counts → differential test → aggregate significant windows

THOR
replicated genomic read signal → HMM/state model → differential regions

EpiSegDiff
genome-wide differential evidence → change-point segmentation → adaptive domains
```

Do **not** claim that EpiSegDiff is novel simply because it is
peak-independent. The specific hypothesis is that **change-point-based
adaptive differential-domain reconstruction** provides useful behavior
distinct from fixed-region, window-aggregation, and HMM formulations.

------------------------------------------------------------------------

## 3. Primary hypotheses

**H1 Boundary recovery.** EpiSegDiff accurately reconstructs known
differential-event boundaries.

**H2 Broad domains.** For broad events, EpiSegDiff reduces artificial
fragmentation while retaining bp-level recovery.

**H3 Boundary-crossing events.** Events crossing conventional peak
boundaries are represented more completely by adaptive segmentation.

**H4 Nearby events.** Separate events remain separable when an
intervening signal supports a true boundary.

**H5 Calibration.** Under null data, EpiSegDiff must not exhibit
material uncontrolled false-positive/FDR inflation.

**H6 Tradeoff.** Better domain reconstruction must not be purchased
through unacceptable merging or false positives.

**H7 Real-data coherence.** EpiSegDiff domains should remain
replicate-consistent and biologically coherent.

**H8 Assay geometry.** Broad marks may benefit more from adaptive
domains than narrow TF binding. This is a hypothesis, not an established
fact.

------------------------------------------------------------------------

## 4. Required comparators

### 4.1 Fixed-peak comparator

Retain the validated v0.5 comparator.

For continuous BigWig:

``` text
fixed union/consensus peaks
→ summarize signal per peak
→ same continuous replicate-aware model where appropriate
→ differential fixed regions
```

Where BAMs are available:

``` text
fixed peaks → fragment counts → appropriate count model → differential peaks
```

DiffBind may be included where technically appropriate.

### 4.2 csaw --- highest-priority external comparator

Use a documented standard csaw workflow:

``` text
BAM
→ sliding-window counts
→ filtering/normalization
→ replicate-aware testing
→ aggregation / region-level inference
```

Use defensible documented settings. Do not deliberately handicap csaw. A
small predeclared sensitivity grid is allowed.

### 4.3 THOR

Use THOR where assay/design/input requirements are met.

Purpose:

> HMM/state-based differential segmentation versus change-point
> differential segmentation.

If incompatible, record `NOT_AVAILABLE` with the exact reason.

### 4.4 Optional

PePr or another established replicated differential peak method may be
secondary. Do not add methods merely to enlarge the benchmark.

------------------------------------------------------------------------

## 5. Fairness rules

For every run record:

``` text
dataset
sample IDs
condition/replicate
assay and target/mark
genome build
BAM/BigWig/peak provenance
normalization
blacklist/chromosome policy
contrast
software version
configuration
```

Use the same biological samples and contrast wherever possible.

Separate two questions:

1.  **Region-definition comparison:** keep the observation/statistical
    model as similar as possible.
2.  **End-to-end comparison:** allow BAM/count versus BigWig/continuous
    differences, but explicitly acknowledge the confounding.

Never feed normalized BigWig values to edgeR/NB.

Never tune the primary analysis to equalize call counts.

------------------------------------------------------------------------

## 6. Datasets

### Dataset A --- GSE280947 JUN

Role: narrow/localized TF-binding benchmark.

Primary EpiSegDiff profile remains the validated high-signal
configuration:

``` text
200 bp
mean
log2(score + 0.01)
absolute filter 0.6
limma
Gaussian PELT
min_bins=3
JNKi-DMSO
```

Reference: 2 GAIN / 688 LOSS.

`filter=none` remains sensitivity only.

### Dataset B --- broad-domain benchmark

Select a replicated broad-mark perturbation dataset, preferably H3K27me3
or H3K36me3, with BAM strongly preferred and BigWig preferred.

Requirements:

-   =2 biological replicates/condition;

-   clear perturbation/control;

-   adequate QC;

-   known genome build;

-   input compatible with at least EpiSegDiff + csaw and preferably
    THOR.

GSE110378 may be retained only if it satisfies the intended comparator
requirements.

### Dataset C --- optional but strongly preferred

A different signal geometry, e.g. ATAC-seq, CUT&RUN TF, or another
histone mark.

Create:

`reports/v0.6_dataset_and_comparator_selection.md`

Do not select datasets after seeing favorable EpiSegDiff results.

------------------------------------------------------------------------

## 7. Predeclared benchmark plan

Before final results, create:

`reports/v0.6_predeclared_benchmark_plan.md`

Freeze:

-   datasets;
-   comparators;
-   synthetic parameter grid;
-   event definitions;
-   overlap/matching rules;
-   primary metrics;
-   null simulations/permutations;
-   FDR thresholds;
-   software versions;
-   parameter-selection rules;
-   failure criteria.

Any later change must be marked `POST_HOC`.

------------------------------------------------------------------------

## 8. Synthetic truth benchmark

Generate replicate-aware synthetic genomic data with known truth.

Vary at minimum:

### Direction

`GAIN`, `LOSS`

### Width

`200 bp, 500 bp, 1 kb, 2 kb, 5 kb, 10 kb, 25 kb, 50 kb`

### Effect

`weak, moderate, strong`, numerically defined in the manifest.

### Replicate noise

`low, moderate, high`

### Boundary shape

`sharp, gradual`

For gradual events define the truth transition interval explicitly.

### Event spacing

Pairs with `adjacent, small-gap, medium-gap, large-gap`, including
same-direction and opposite-direction pairs.

### Peak geometry

Truth that: - matches a conventional peak; - is narrower than a peak; -
is broader than a peak; - crosses one peak boundary; - crosses multiple
peaks; - lies between peaks; - is truly peak-distal.

### Background

Include flat, heteroscedastic, locally elevated, and long zero/near-zero
backgrounds.

Do not simulate only conditions favorable to EpiSegDiff.

Use \>=50 random replicates per primary simulation cell; prefer 100
where feasible. Store all seeds.

------------------------------------------------------------------------

## 9. Synthetic metrics

For every method/cell compute:

**Detection** - region sensitivity/precision; - bp recall/precision; -
F1.

**Direction** - GAIN/LOSS accuracy; - opposite-direction error rate.

**Boundary** - left/right absolute boundary error; - mean/median/p95
boundary error; - IoU.

**Fragmentation** - calls per true event; - probability of one-call
representation.

**Merge error** - probability that two true events are incorrectly
merged.

**Split error** - probability that one true event is incorrectly split.

Call count alone is not a performance metric.

------------------------------------------------------------------------

## 10. Secondary matched-call benchmark

In addition to primary method-valid thresholds, perform a secondary
analysis at approximately matched:

-   number of calls, or
-   called genomic bp.

Label it explicitly `SECONDARY_MATCHED_CALL`.

Purpose: determine whether boundary/fragmentation differences persist
after discovery volume is controlled.

Do not use it as the primary benchmark.

------------------------------------------------------------------------

## 11. Statistical calibration --- critical gate

EpiSegDiff uses adaptive segmentation plus statistical inference. v0.6
must test whether segmentation/testing dependence creates
selection-related inflation.

### 11.1 Global null simulation

Generate no-condition-effect data preserving realistic: - signal
distribution; - replicate variance; - spatial autocorrelation; -
low/zero signal structure.

Run the **complete** pipeline including segmentation.

Record: - probability of any DER; - false DER count; - false
differential bp; - empirical FDP; - p-value/FDR distributions.

### 11.2 Real-data label permutation

Use all valid unique permutations when feasible; otherwise a predeclared
random set. Target \>=100 valid permutations when design permits.

Do not create invalid pseudo-replication.

### 11.3 Nominal calibration

At nominal FDR:

`0.01, 0.05, 0.10`

compare nominal and empirical error.

Do not claim formal FDR control merely because few null calls occur.

### 11.4 Selection-bias diagnostic

Where possible compare:

``` text
A. same data define segments and test them
B. independent/split information defines segments and tests them
```

Possible diagnostic:

``` text
replicate subset 1 → boundaries
replicate subset 2 → inference
```

and reverse.

If real replication is insufficient, use synthetic data.

Final classification must be one of:

``` text
well calibrated within tested regimes
conservative
mildly anti-conservative
materially anti-conservative
unresolved
```

with quantitative evidence.

------------------------------------------------------------------------

## 12. If calibration fails

Do not hide the result.

Explore at least one remediation concept:

-   sample splitting;
-   cross-fitting;
-   independent segmentation statistic;
-   permutation-derived region significance.

This is exploratory. Do not silently replace the production method
without separate validation.

------------------------------------------------------------------------

## 13. Real-data benchmark

For compatible datasets run:

``` text
EpiSegDiff
fixed-peak comparator
csaw
THOR where available
```

Report: - GAIN/LOSS; - total differential bp; - number/width of
regions; - effect/sign concordance; - replicate consistency; -
fragmentation; - boundary displacement.

For broad marks additionally report: - calls per coherent broad
domain; - fraction of domain covered; - internal gaps.

For narrow TF data additionally report: - summit proximity; - motif
support; - peak-boundary extension; - truly peak-distal fraction.

Do not directly interpret method-specific FDRs as equivalent truth
probabilities when hypothesis spaces differ.

------------------------------------------------------------------------

## 14. v0.5 Class-C closure

Decompose GSE280947 EpiSegDiff-only regions into predeclared operational
subclasses:

``` text
C1 — overlaps conventional peak but fails reciprocal-overlap criterion
C2 — peak-boundary extension
C3 — between/multi-peak domain
C4 — truly peak-distal
```

For each compare: - distance to nearest peak; - distance to summit; -
JUN/AP-1 motif support; - mean/max signal; - effect/FDR; - replicate
consistency; - promoter/enhancer annotation; - nearest TSS; -
IgG/background where available.

The goal is to distinguish:

``` text
novel peak-independent biology
```

from:

``` text
adaptive boundary/domain reconstruction
```

Either conclusion is acceptable.

------------------------------------------------------------------------

## 15. Promoter-enrichment diagnostic

v0.5 observed a high promoter fraction in EpiSegDiff-only regions. Test
whether this is biological or geometric.

Required:

1.  width-matched comparison;
2.  distance-to-peak-matched comparison;
3.  chromosome/width-preserving shuffled background;
4.  C1/C2/C3/C4 separately;
5.  odds ratios and confidence intervals.

Do not interpret raw promoter percentages alone.

------------------------------------------------------------------------

## 16. Biological support

Use only predeclared evidence appropriate to each assay:

-   motif;
-   RNA-seq;
-   ATAC-seq;
-   another histone mark/TF;
-   independent cohort;
-   appropriate functional annotations.

Prioritize orthogonal evidence over agreement with another caller.

Missing evidence must be `NOT_AVAILABLE`, not inferred.

------------------------------------------------------------------------

## 17. Runtime/reproducibility

Record per method:

-   wall time;
-   CPU time where available;
-   peak RSS;
-   disk usage;
-   worker count;
-   restartability.

Preserve EpiSegDiff chromosome-wise checkpoint/restart.

Record exact versions for EpiSegDiff, R/Bioconductor, csaw, edgeR,
DiffBind if used, THOR, PePr if used, genome annotation and motif
database.

Store session/environment info, commands/configs and random seeds.

------------------------------------------------------------------------

## 18. Required outputs

``` text
results/v0.6/dataset_manifest.tsv
results/v0.6/software_versions.tsv
results/v0.6/method_configs.tsv
results/v0.6/synthetic_design.tsv
results/v0.6/synthetic_region_metrics.tsv
results/v0.6/synthetic_boundary_metrics.tsv
results/v0.6/synthetic_fragmentation.tsv
results/v0.6/synthetic_merge_split_errors.tsv
results/v0.6/null_simulation.tsv
results/v0.6/permutation_calibration.tsv
results/v0.6/fdr_calibration.tsv
results/v0.6/selection_bias_diagnostic.tsv
results/v0.6/realdata_calls_summary.tsv
results/v0.6/realdata_overlap.tsv
results/v0.6/realdata_effect_concordance.tsv
results/v0.6/realdata_fragmentation.tsv
results/v0.6/realdata_boundary_metrics.tsv
results/v0.6/gse280947_classC_subclasses.tsv
results/v0.6/gse280947_classC_support.tsv
results/v0.6/promoter_geometry_diagnostic.tsv
results/v0.6/runtime_benchmark.tsv
results/v0.6/regression_summary.tsv
```

Keep `FAILED`, `NOT_RUN`, and `NOT_AVAILABLE` rows with reasons.

------------------------------------------------------------------------

## 19. Required figures

``` text
figures/v0.6_method_landscape.*
figures/v0.6_synthetic_power.*
figures/v0.6_boundary_error_by_width.*
figures/v0.6_fragmentation_by_width.*
figures/v0.6_merge_split_error.*
figures/v0.6_effect_noise_phase_diagram.*
figures/v0.6_null_pvalue_calibration.*
figures/v0.6_nominal_vs_empirical_fdr.*
figures/v0.6_realdata_gain_loss.*
figures/v0.6_realdata_domain_fragmentation.*
figures/v0.6_classC_decomposition.*
figures/v0.6_promoter_geometry_control.*
figures/v0.6_runtime_memory.*
```

Generate genome-track panels for predeclared examples of shared events,
boundary extension, multi-peak domains, truly peak-distal candidates,
and broad-domain events. Selection rules must be documented.

------------------------------------------------------------------------

## 20. Required reports

``` text
reports/v0.6_dataset_and_comparator_selection.md
reports/v0.6_predeclared_benchmark_plan.md
reports/v0.6_comparator_installation_and_validation.md
reports/v0.6_synthetic_benchmark.md
reports/v0.6_statistical_calibration.md
reports/v0.6_gse280947_method_comparison.md
reports/v0.6_broad_domain_benchmark.md
reports/v0.6_classC_boundary_biology_closure.md
reports/v0.6_runtime_reproducibility.md
reports/v0.6_regression_report.md
reports/v0.6_final_methodological_acceptance_report.md
```

Use `FACT`, `INFERENCE`, `SPECULATION`, `UNRESOLVED`.

------------------------------------------------------------------------

## 21. Acceptance gates

### Gate A --- comparator validity

Comparator workflows must be validated and fairly configured. csaw is
mandatory where compatible BAM data exist. THOR is required where
compatible or explicitly `NOT_AVAILABLE`.

### Gate B --- synthetic direction

GAIN/LOSS accuracy must remain high across the predeclared grid.

### Gate C --- boundary performance

A claimed boundary/domain advantage must occur reproducibly in a
**predeclared** event class.

### Gate D --- fragmentation versus merging

Broad-domain gains must not come from systematic merging of biologically
separate events.

### Gate E --- null calibration --- CRITICAL

Null simulations/permutations must not reveal material uncontrolled
false-positive inflation.

Material anti-conservatism blocks full PASS.

### Gate F --- selection-bias assessment --- CRITICAL

The dependence between adaptive segmentation and subsequent testing must
be explicitly evaluated.

### Gate G --- real-data coherence

Real-data EpiSegDiff regions must remain directionally and replicate
coherent.

### Gate H --- Class-C interpretation

The v0.5 EpiSegDiff-only class must be decomposed sufficiently to
distinguish peak-distal discovery from boundary/domain effects.

### Gate I --- broad-domain benchmark

At least one compatible broad-mark dataset must be evaluated.

### Gate J --- reproducibility

Versions/configs/seeds/manifests must permit reproduction.

### Gate K --- regression

All critical v0.4.2/v0.5 guarantees remain passing.

------------------------------------------------------------------------

## 22. Decision rules

Allowed:

`PASS`, `CONDITIONAL PASS`, `FAIL`.

### PASS

Requires all critical gates, especially statistical calibration.

PASS does not require EpiSegDiff to win every metric.

A valid supported conclusion might be:

> EpiSegDiff improves reconstruction of broad and boundary-crossing
> differential domains while maintaining acceptable calibration.

or:

> EpiSegDiff primarily reduces fragmentation and improves domain
> representation rather than increasing sensitivity.

Use only the statement supported by results.

### CONDITIONAL PASS

Use if methodological distinction is supported but a major comparator,
independent dataset, or calibration issue remains incomplete.

### FAIL

Use if: - comparator fairness fails; - null false positives are
materially inflated; - segmentation/testing selection invalidates
inference; - broad-domain gains arise from uncontrolled merging; -
conclusions require post hoc tuning; - critical regressions fail.

------------------------------------------------------------------------

## 23. Claims matrix

The final report must mark every claim as `SUPPORTED`, `PARTIAL`, or
`UNSUPPORTED`:

  Claim                                     Status
  ----------------------------------------- --------
  Detects differential epigenomic regions   
  Is peak-independent                       
  Reconstructs adaptive boundaries          
  Improves broad-domain representation      
  Reduces fragmentation                     
  Improves sensitivity                      
  Improves specificity                      
  Discovers truly peak-distal biology       
  Has calibrated region-level FDR           
  Generalizes across assay types            
  Is superior to csaw                       
  Is superior to THOR                       

No ambiguous cells.

------------------------------------------------------------------------

## 24. Final acceptance questions

The final report must answer:

1.  What methodological niche remains after comparison with fixed peaks,
    csaw and THOR?
2.  Which comparisons isolate region definition from statistical-model
    differences?
3.  Were all comparators fairly configured?
4.  Which method has lowest synthetic boundary error?
5.  How does performance vary with width, effect size and noise?
6.  What happens with gradual boundaries?
7.  Which method fragments broad events most?
8.  Which method incorrectly merges nearby events most?
9.  Can opposite-direction neighbors be separated?
10. Does EpiSegDiff retain any advantage in matched-call analysis?
11. Under global null, what false DER burden occurs?
12. What is empirical error at nominal FDR 0.01/0.05/0.10?
13. Do label permutations support calibration?
14. Is segmentation/testing selection bias detectable?
15. Does split-sample/cross-fit analysis alter conclusions?
16. How do methods compare on GSE280947 where compatible?
17. How do methods compare on the broad-domain dataset?
18. Does EpiSegDiff reduce fragmentation without excessive merging?
19. What fraction of v0.5 Class C is C1/C2/C3/C4?
20. Does motif support persist in truly peak-distal C4?
21. Does promoter enrichment survive width/distance-matched controls?
22. Is promoter enrichment biological, geometric, or unresolved?
23. Do conclusions reproduce across assay geometries?
24. What are computational tradeoffs?
25. Which claims are supported, partial, or unsupported?
26. What observations could falsify the preferred interpretation?
27. What alternative explanations remain?
28. Is v0.6 PASS, CONDITIONAL PASS or FAIL?
29. Is the project ready to move from core-method development to
    manuscript-focused validation?

------------------------------------------------------------------------

## 25. Execution order

``` text
dataset/comparator selection
→ predeclared benchmark plan
→ install/validate comparators
→ synthetic benchmark
→ boundary/fragmentation/merge-split analysis
→ global-null calibration
→ permutation + selection-bias diagnostics
→ GSE280947 benchmark
→ Class-C + promoter closure
→ broad-domain benchmark
→ optional third dataset
→ runtime/reproducibility
→ v0.4.2/v0.5 regressions
→ final methodological acceptance
```

If statistical calibration shows material anti-conservatism, stop
superiority claims and prioritize diagnosis/remediation.

------------------------------------------------------------------------

## 26. Stop conditions

Stop and report if:

-   comparator installation cannot be reproduced;
-   sample/genome provenance is incompatible;
-   no valid dataset can support required csaw/THOR comparisons;
-   synthetic tests reveal a correctness defect;
-   null calibration is materially invalid;
-   direction/width/replication regressions appear;
-   benchmark fairness cannot be maintained.

Do not bury these failures.

------------------------------------------------------------------------

## 27. Agent instruction

Implement **EpiSegDiff v0.6 Methodological Benchmark and Statistical
Calibration** only.

Priority order:

1.  csaw as principal de novo/window comparator;
2.  THOR as HMM/state comparator where compatible;
3.  fixed peaks as conventional reference;
4.  synthetic truth for boundary and fragmentation claims;
5.  global-null/permutation calibration;
6.  explicit segmentation/testing selection-bias assessment;
7.  GSE280947 Class-C decomposition and promoter-geometry analysis;
8.  at least one broad-domain real-data benchmark.

Do not optimize parameters to make EpiSegDiff win.

Do not equate more/fewer calls with superiority.

Do not introduce DESeq2 merely as another comparator.

Do not start NB-FPOP unless separately specified.

Do not modify validated v0.4.2 BigWig infrastructure without regression
tests.

Finally generate:

`reports/v0.6_final_methodological_acceptance_report.md`

Then stop for formal review.
