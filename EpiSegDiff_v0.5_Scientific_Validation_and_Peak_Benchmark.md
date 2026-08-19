# EpiSegDiff v0.5 — Scientific Validation and Conventional Peak-Based Benchmark

**Status:** Agent-executable specification  
**Baseline:** v0.4.2 PASS  
**Goal:** Determine whether adaptive differential segmentation adds reproducible scientific information beyond defensible conventional peak-based analysis.

## 1. Scientific transition

v0.4.x established native BigWig infrastructure, including scalable `filter=none`. v0.5 must move from engineering validation to scientific validation.

Primary question:

> Does EpiSegDiff recover reproducible, biologically coherent differential genomic domains that fixed peak-based analysis misses, fragments, or represents less appropriately?

The benchmark must not be tuned to make EpiSegDiff win.

## 2. Hypotheses

- **H1 Direction:** known perturbation effects are recovered in biologically relevant signal domains.
- **H2 Reproducibility:** DER effects are consistent across replicates.
- **H3 Peak-independent recovery:** some valid differential regions may fall outside or cross fixed peak boundaries.
- **H4 Domain representation:** segmentation may better represent broad, shifted, or multi-peak differential domains.
- **H5 Biological support:** EpiSegDiff-specific regions should show independent support comparable to or exceeding matched background.
- **H6 Robustness:** conclusions should survive reasonable predeclared parameter changes.

Do not claim BigWig=BAM, more calls=better, EpiSegDiff-only=true, conventional-only=false, or IDR=differential truth.

## 3. Required datasets

### Dataset A — GSE280947

Use as continuity/development benchmark. Preserve the validated primary profile:

- continuous BigWig
- 200 bp
- mean
- `log2(score+0.01)`
- absolute filter 0.6
- limma
- Gaussian PELT
- `min_bins=3`
- `JNKi-DMSO`

Reference: 2 GAIN / 688 LOSS.

Use `filter=none` only as secondary sensitivity; v0.4.2 showed 4,866 GAIN / 4,608 LOSS and predominantly low-signal none-only calls. Do not treat 0.6 as universal.

### Dataset B — independent validation dataset

Select at least one independent public perturbation dataset with:
- >=2 biological replicates/condition;
- same assay target/mark across conditions;
- clear control/perturbation;
- adequate signal;
- public BigWig and preferably BAM;
- known genome build.

Preferred: CUT&RUN, CUT&Tag, ChIP-seq, or ATAC-seq.

Create `reports/v0.5_dataset_selection.md` with accession, assay, target, conditions, replicates, BAM/BigWig/peak availability, genome, expected biology, and selection/rejection reason.

Do not select based on favorable EpiSegDiff results.

## 4. Predeclared analysis plan

Before final benchmarking create `reports/v0.5_predeclared_analysis_plan.md` freezing:
- datasets;
- samples/contrast;
- EpiSegDiff config;
- comparator;
- genomic evaluation universe;
- overlap definitions;
- motif source/set;
- annotation source;
- gene-linkage rule;
- gene-set source/universe;
- negative controls;
- primary metrics;
- multiple-testing correction.

Any later change must be labeled **post hoc**.

## 5. Comparator design

### Layer 1 — primary BigWig region-definition benchmark

Hold the observation/statistical model as similar as possible:

```text
EpiSegDiff:
adaptive segmentation → regional BigWig signal → limma

Conventional:
fixed union/consensus peaks → peak BigWig signal → limma
```

This isolates adaptive region definition from fixed peak boundaries.

### Layer 2 — full conventional count workflow

Only where suitable BAMs exist:

```text
consensus/union peaks
→ fragment counts
→ normalization
→ count/NB differential model
→ GAIN/LOSS
```

Use the existing count pathway/edgeR if appropriate. Never send normalized BigWig values into edgeR.

Interpret Layer 2 cautiously because both region definition and observation model differ.

## 6. Fairness requirements

Both methods must use the same biological samples and contrast whenever possible. Record input provenance, normalization, genome build and replicates.

Do not tune thresholds to equalize call counts in the primary benchmark.

Define a common evaluable genomic domain before sensitivity-like comparisons.

Deposited IDR peaks may provide context but are not differential truth.

## 7. Region classes

Classify:
- **A shared**
- **B conventional-only**
- **C EpiSegDiff-only**

Evaluate overlap using:
- >=1 bp;
- reciprocal overlap >=50%;
- bp-level overlap.

Separate same-direction and opposite-direction overlap. Investigate opposite-direction overlaps explicitly.

## 8. Primary metrics

For matched regions report:
- effect Pearson/Spearman;
- sign concordance;
- replicate consistency;
- baseline/perturbed/max/mean signal;
- absolute effect;
- p/FDR with hypothesis-space caveat;
- width;
- total differential bp;
- nearest peak-summit distance;
- conventional peaks per EpiSegDiff domain;
- boundary displacement.

Define and report a fragmentation metric plus its raw numerator/denominator.

## 9. Synthetic truth benchmark

Extend synthetic tests to:
- narrow GAIN;
- narrow LOSS;
- broad GAIN;
- broad LOSS;
- nearby same-direction events;
- nearby opposite-direction events;
- signal event crossing a fixed peak boundary.

Because truth is known, compute:
- bp precision/recall;
- region detection rate;
- direction accuracy;
- boundary error.

This directly tests whether adaptive segmentation provides a representational advantage.

## 10. Ablation

Run, where possible with the same continuous model:

1. full adaptive EpiSegDiff;
2. fixed-width bins without adaptive segmentation;
3. fixed conventional peak regions.

Purpose: distinguish the contribution of segmentation from merely using BigWig signal.

## 11. Biological validation

### Motif support
For TF-binding datasets, compare shared, conventional-only, EpiSegDiff-only and matched background using a predeclared motif set.

Report motif-positive fraction, odds ratio, CI, adjusted p-value and motif-to-region-center distance.

For JUN, use an appropriate predeclared JUN/AP-1 motif source. Do not choose motifs after seeing results.

### Annotation
Use one annotation hierarchy for all classes: promoter, proximal/enhancer if justified, intron, exon, intergenic, other.

### Gene linkage
Start with a predeclared nearest-TSS linkage. Do not equate nearest gene with regulatory target.

### Pathways
Use the same gene universe for all region classes. Record database, test, correction and enrichment effect.

### Orthogonal evidence
Where available, prioritize independent RNA-seq, ATAC, another mark/TF, or an independent cohort. Orthogonal evidence is stronger than agreement with another peak caller.

## 12. Negative controls

Require at least one valid null analysis, preferably:
- condition-label permutation where design permits;
- matched background/null genomic regions;
- descriptive within-condition pseudo-contrast where appropriate.

Expected primary biological signal should substantially collapse under a valid null. Do not overinterpret under-replicated pseudo-contrasts.

## 13. Robustness

Use v0.4.1/v0.4.2 sensitivity knowledge rather than retuning.

For GSE280947 compare primary absolute-0.6 conclusions with `filter=none` sensitivity.

Assess whether the scientific claim changes under reasonable bin/pseudocount settings already validated.

Do not select the setting yielding the most DERs.

## 14. Independent replication

Repeat the key benchmark on Dataset B.

A strong general claim requires the central advantage to reproduce qualitatively on the independent dataset.

If it does not reproduce, report dataset-specific evidence rather than broadening the claim.

## 15. Required outputs

Results:
```text
results/v0.5/dataset_manifest.tsv
results/v0.5/method_configs.tsv
results/v0.5/primary_calls_episegdiff.tsv
results/v0.5/primary_calls_conventional.tsv
results/v0.5/region_classification.tsv
results/v0.5/overlap_metrics.tsv
results/v0.5/effect_concordance.tsv
results/v0.5/replicate_consistency.tsv
results/v0.5/region_geometry.tsv
results/v0.5/motif_validation.tsv
results/v0.5/genomic_annotation.tsv
results/v0.5/gene_linkage.tsv
results/v0.5/pathway_validation.tsv
results/v0.5/negative_controls.tsv
results/v0.5/synthetic_benchmark.tsv
results/v0.5/ablation_results.tsv
results/v0.5/sensitivity_summary.tsv
```

Keep FAILED / NOT_RUN / NOT_AVAILABLE rows with reasons.

Figures:
```text
figures/v0.5_method_overview.*
figures/v0.5_gain_loss_comparison.*
figures/v0.5_region_overlap.*
figures/v0.5_effect_concordance.*
figures/v0.5_region_widths.*
figures/v0.5_fragmentation.*
figures/v0.5_motif_support.*
figures/v0.5_annotation_comparison.*
figures/v0.5_pathway_comparison.*
figures/v0.5_negative_control.*
figures/v0.5_synthetic_boundary_accuracy.*
figures/v0.5_ablation.*
```

Representative genome-track loci must be selected by a documented rule, not visual attractiveness.

Reports:
```text
reports/v0.5_dataset_selection.md
reports/v0.5_predeclared_analysis_plan.md
reports/v0.5_conventional_comparator_validation.md
reports/v0.5_gse280947_benchmark.md
reports/v0.5_independent_dataset_benchmark.md
reports/v0.5_synthetic_validation.md
reports/v0.5_negative_control_report.md
reports/v0.5_ablation_report.md
reports/v0.5_biological_validation.md
reports/v0.5_final_scientific_acceptance_report.md
```

Use FACT / INFERENCE / UNRESOLVED labels.

## 16. Acceptance gates

### Gate A — fairness
Same samples/contrast where required; provenance documented; defensible comparator; no BigWig→NB misuse; no primary call-count tuning.

### Gate B — comparator validity
Conventional workflow passes internal QC and direction checks; region universe is correctly constructed.

### Gate C — synthetic truth
Both directions recovered; narrow/broad/adjacent/boundary-shift cases behave acceptably.

### Gate D — real-data reproducibility
Primary EpiSegDiff effects are replicate-supported and not predominantly unstable low-signal artifacts.

### Gate E — independent biological support
EpiSegDiff-specific findings are evaluated against appropriate background. PASS may include a transparent finding of no added support; do not fabricate an advantage.

### Gate F — negative controls
Primary biological signal substantially collapses under valid null analyses.

### Gate G — ablation
The contribution of adaptive segmentation is separable from simply using continuous signal.

### Gate H — robustness
Main conclusion is not reversed by reasonable predeclared sensitivity.

### Gate I — independent replication
At least one independent real dataset is analyzed before a strong general scientific claim.

### Gate J — regression
All critical v0.4.2 infrastructure guarantees remain passing.

## 17. Scientific decision

Allowed: **PASS / CONDITIONAL PASS / FAIL**.

PASS requires a fair benchmark plus evidence of a reproducible scientific contribution. The report must state exactly what contribution is supported.

CONDITIONAL PASS: valid implementation/benchmark with promising but incomplete independent validation or dataset-limited evidence.

FAIL: unfair comparator, irreproducible primary signal, persistent null signal, major synthetic failure, instability, or an apparent advantage dependent on post hoc tuning.

A scientifically negative result is acceptable and must not trigger automatic benchmark redesign.

## 18. Final questions

`reports/v0.5_final_scientific_acceptance_report.md` must answer:

1. What scientific advantage was hypothesized?
2. Was the analysis predeclared?
3. Were samples/contrasts matched?
4. What comparator was used and why?
5. What are primary calls for each method?
6. How much differential bp is shared?
7. Same- vs opposite-direction overlaps?
8. Effect concordance?
9. Replicate support in shared/B/C classes?
10. Signal-strength differences among classes?
11. Does EpiSegDiff reduce fragmentation or capture coherent broader domains?
12. Synthetic boundary accuracy?
13. Boundary-shift recovery?
14. Motif support for EpiSegDiff-specific TF regions?
15. Annotation differences?
16. Linked-gene/pathway coherence?
17. Orthogonal evidence?
18. Null-control behavior?
19. Fixed-bin ablation result?
20. Fixed-peak ablation result?
21. Robustness to reasonable parameters?
22. Effect of `filter=none` sensitivity?
23. Does Dataset B reproduce the key conclusion?
24. What would falsify the claimed advantage?
25. What alternative explanations remain?
26. PASS / CONDITIONAL PASS / FAIL?
27. What exact claim is justified?
28. What claims remain unjustified?
29. What next version is warranted?

## 19. Execution order

```text
dataset selection
→ predeclared analysis plan
→ comparator validation
→ synthetic truth
→ GSE280947 primary benchmark
→ region-class characterization
→ biological validation
→ negative controls
→ ablation
→ robustness
→ independent dataset
→ v0.4.2 regressions
→ final scientific acceptance
```

Stop before final interpretation if sample identity/genome/provenance is incompatible, comparator validation fails, synthetic tests expose a correctness defect, or v0.4.2 regressions fail.

## 20. Agent instruction

Implement **v0.5 Scientific Validation** only.

The goal is to test whether adaptive differential segmentation adds reproducible information beyond fixed peak-based analysis—not to maximize EpiSegDiff performance.

Start with dataset selection and a predeclared analysis plan.

Use the continuous fixed-peak comparator as the primary BigWig benchmark to isolate region-definition effects. Add a BAM/count comparator only when suitable BAMs genuinely exist.

Do not start NB-FPOP. Do not introduce DESeq2 without a separately justified count-assay requirement. Do not use IDR peaks as differential ground truth. Do not tune parameters to make EpiSegDiff win.

Require synthetic truth, negative controls, ablation, biological validation, robustness analysis, and at least one independent real dataset before a strong general claim.

Finally generate:
`reports/v0.5_final_scientific_acceptance_report.md`

Then stop for formal review.
