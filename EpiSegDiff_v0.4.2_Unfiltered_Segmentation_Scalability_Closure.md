# EpiSegDiff v0.4.2 — Unfiltered Continuous Segmentation Scalability Closure

**Status:** Agent-executable implementation and acceptance specification
**Baseline:** v0.4.1
**Release type:** Narrow engineering/scalability closure
**Primary dataset:** GSE280947 JUN BigWig continuous assay

## 1. Objective

Close the only remaining v0.4.1 blocker: genome-wide `filter=none` segmentation did not complete on ~13.6 million 200-bp bins.

The release must make unfiltered whole-genome continuous segmentation tractable without changing established statistical semantics. This is not a new scientific-feature release.

## 2. Frozen baseline

- BigWig is a native continuous assay.
- BigWig must never enter count/NB/edgeR inference.
- Reference: 200 bp, `summary_fun=mean`, `log2(score+0.01)`, absolute filter 0.6, limma, Gaussian PELT, `min_bins=3`, `JNKi-DMSO`.
- Reference result: 2 GAIN / 688 LOSS; minimum width 600 bp.
- v0.4.1 completed the 100/200/500-bp and 0.001/0.01/0.05/0.1 grid.
- Only `filter=none` remained NOT_RUN because the genome-wide PELT job was killed after >30 min at 100% CPU.

## 3. Scope

Required:
1. Profile the scalability failure.
2. Segment chromosomes independently.
3. Split true coordinate/validity discontinuities into independent blocks.
4. Make execution memory-bounded.
5. Add checkpoint/resume.
6. Add deterministic chromosome-level parallelism.
7. Prove equivalence to the legacy implementation on tractable domains.
8. Complete the full GSE280947 `filter=none` run.
9. Compare `none` against absolute 0.6.
10. Re-run critical regressions.
11. Issue PASS / CONDITIONAL PASS / FAIL.

Out of scope: NB-FPOP, DESeq2, BAM/BigWig equivalence, competitors, new assay types, new segmentation science, publication benchmarking, or rerunning the full v0.4.1 grid.

## 4. Scientific invariants

### 4.1 Chromosome isolation
A change point may never cross chromosomes. PELT state must be reset for every chromosome.

### 4.2 True gaps are boundaries
Missing/non-finite bins or genuine coordinate discontinuities may define independent blocks. Arbitrary computational chunks must NOT become biological boundaries.

### 4.3 Statistical semantics are frozen
Do not change bin size, pseudocount, score, model, penalty, `min_bins`, or contrast merely to make the run finish.

### 4.4 Optimization must be equivalent
On data where both old and new implementations finish, boundaries and directions must agree; numerical fields must agree within documented floating-point tolerance.

## 5. Phase 1 — Profile the current implementation

Create `reports/v0.4.2_phase1_scalability_profile.md`.

Profile several actual mm10 chromosomes spanning small/medium/large sizes. Record:
- chromosome
- number of bins / finite bins / valid blocks
- total and segmentation runtime
- peak memory
- number of segments
- status

Determine whether the current code:
- concatenates chromosomes;
- makes large intermediate copies;
- retains unnecessary genome-wide objects;
- behaves poorly on long zero/near-zero regions.

Label findings FACT / INFERENCE / UNRESOLVED.

## 6. Phase 2 — Chromosome-wise dispatcher

Implement in the existing pipeline, conceptually:

```r
segment_genome <- function(x, config) {
  out <- lapply(split_by_chromosome(x), segment_chromosome, config=config)
  combine_segments(out)
}
```

Exact implementation may differ.

Output must retain deterministic chromosome ordering and:
`chrom,start,end,n_bins,effect,direction,statistical fields,segmentation metadata`.

## 7. Phase 3 — Contiguous valid blocks

Within each chromosome identify independent blocks created by genuine:
- missing bins;
- non-finite transformed values;
- excluded genomic gaps;
- coordinate discontinuities.

Segment each valid block independently.

**Do not** split biology at arbitrary memory/compute chunk boundaries. If arbitrary chunks are necessary, implement boundary reconciliation and demonstrate equivalence.

## 8. Phase 4 — Memory-bounded execution

Preferred pattern:

```text
select chromosome
→ derive required score
→ segment
→ write/checkpoint
→ release chromosome workspace
→ next chromosome
```

Avoid redundant genome-wide copies of tables, transformed matrices, z-scores and PELT workspaces. Measure peak memory before/after.

## 9. Phase 5 — Checkpoint/resume

Create e.g.:
`results/v0.4.2/filter_none/checkpoints/chr*.rds`

Checkpoint validity must include:
- source/cache identity or checksum
- genome
- bin size
- summary function
- transform and pseudocount
- filter
- segmentation method/parameters
- software version

On restart: valid completed chromosomes are reused; missing/incompatible ones are recomputed. Never silently reuse incompatible checkpoints.

## 10. Phase 6 — Deterministic parallelism

Support configurable chromosome-level workers, e.g.:

```yaml
execution:
  workers: 4
  checkpoint: true
  resume: true
```

Requirements:
- identical scientific output for workers=1 and workers>1;
- no uncontrolled nested parallelism;
- collision-safe checkpoint writes;
- same semantics locally and on HPC.

## 11. Phase 7 — Legacy equivalence gate

Before full genome interpretation, compare old vs new segmentation on tractable chromosomes/subsets with identical parameters.

Compare:
- segment count
- chrom/start/end/n_bins
- effect
- direction
- p/FDR where applicable

Acceptance: exact deterministic boundaries where expected; numerical equality within stated tolerance; zero direction differences.

If unexplained differences occur, stop before interpreting the full run.

## 12. Phase 8 — Required full `filter=none` run

Use exactly:

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
  method: none
contrast: JNKi - DMSO
inference:
  method: limma
segmentation:
  method: gaussian_pelt
  min_bins: 3
```

Required outputs:
- `results/v0.4.2/filter_none/primary.tsv`
- `results/v0.4.2/filter_none/run_manifest.yaml`
- `results/v0.4.2/filter_none/chromosome_runtime.tsv`

Record total bins, tested segments, GAIN, LOSS, total DER, width summaries, runtime, peak memory, workers, resumed and failed chromosomes.

## 13. Phase 9 — `none` versus absolute 0.6

Compare against the frozen 2 GAIN / 688 LOSS reference.

Compute:
- shared/reference-only/none-only DER bp
- bp-level Jaccard
- direction concordance
- Pearson/Spearman effect concordance
- strong-reference recovery
- GAIN/LOSS and tested segment counts

Characterize `none-only` calls by:
- mean/max linear signal
- effect
- FDR
- width
- replicate consistency where defined
- IgG/background context where available

Do not label `none-only` calls false positives merely because they fail the 0.6 filter.

## 14. Phase 10 — Low-signal characterization

Stratify the `filter=none` tested space deterministically:

```text
exact zero
(0,Q25]
(Q25,Q50]
(Q50,Q75]
(Q75,max]
```

For each report:
- tested count
- FDR<=0.05 count
- GAIN / LOSS
- median |effect|
- median FDR
- replicate consistency where meaningful

This is descriptive. Do not invent a new threshold in v0.4.2.

## 15. Phase 11 — Scalability report

Create `reports/v0.4.2_scalability_report.md`.

Include before/after:
- whether full genome completes
- runtime
- peak memory
- restartability
- chromosome independence
- configurable workers

Also report empirical runtime and memory versus chromosome bin count. Do not claim formal asymptotic complexity unless established.

## 16. Phase 12 — Critical regressions

Re-run all critical v0.4/v0.4.1 gates:

1. `continuous+bigwig` never dispatches to edgeR/count inference.
2. Absolute-0.6 reference reproduces ~2 GAIN / 688 LOSS and min width >=600 bp.
3. `JNKi-DMSO`: effect>0=GAIN, effect<0=LOSS; reversed contrast reverses signs.
4. Synthetic GAIN and LOSS are recovered.
5. One replicate/condition cannot generate inferential FDR.
6. workers=1 and workers>1 produce equivalent primary output.
7. Width enforcement remains correct.

## 17. Required deliverables

Results:
```text
results/v0.4.2/filter_none/primary.tsv
results/v0.4.2/filter_none/run_manifest.yaml
results/v0.4.2/filter_none/chromosome_runtime.tsv
results/v0.4.2/filter_comparison.tsv
results/v0.4.2/low_signal_strata.tsv
results/v0.4.2/regression_summary.tsv
```

Reports:
```text
reports/v0.4.2_phase1_scalability_profile.md
reports/v0.4.2_legacy_equivalence_report.md
reports/v0.4.2_filter_none_report.md
reports/v0.4.2_filter_comparison_report.md
reports/v0.4.2_scalability_report.md
reports/v0.4.2_regression_report.md
reports/v0.4.2_final_acceptance_report.md
```

Tests must cover chromosome isolation, valid-block isolation, checkpoint validation/resume, worker determinism, legacy/new equivalence and `filter=none` integration.

## 18. Acceptance gates

### Gate A — Correctness
PASS only if chromosome/gap semantics are correct and optimized output matches legacy behavior on tractable cases.

### Gate B — Scalability
PASS only if the complete GSE280947 200-bp `filter=none` genome-wide run finishes.

A larger machine alone is not sufficient closure if execution remains genome-monolithic and non-restartable.

### Gate C — Restartability
PASS only if interrupted execution resumes completed valid chromosomes and rejects incompatible checkpoints.

### Gate D — Parallel determinism
PASS only if worker count does not alter scientific output.

### Gate E — Filter comparison
PASS only if `none` vs 0.6 is quantitatively compared and none-only calls are characterized.

### Gate F — Regression
PASS only if all prior critical guarantees remain intact.

## 19. Final decision

Allowed: **PASS / CONDITIONAL PASS / FAIL**.

PASS requires Gates A–F.

PASS means native BigWig continuous-assay infrastructure is production-supported under documented assumptions, including tractable unfiltered whole-genome execution.

PASS does **not** mean:
- BigWig == BAM;
- 200 bp is universally optimal;
- pseudocount 0.01 is universally optimal;
- absolute 0.6 is universally optimal.

FAIL conditions include unexplained segmentation changes, cross-chromosome segments, artificial chunk boundaries, BigWig→edgeR leakage, width/direction regression, unsafe checkpoint reuse, worker-dependent calls, or fabricated unreplicated FDR.

## 20. Final report questions

`reports/v0.4.2_final_acceptance_report.md` must answer:

1. What caused the v0.4.1 scalability failure?
2. Was the genome processed monolithically?
3. Is segmentation now chromosome-independent?
4. Are genuine contiguous blocks handled correctly?
5. Does computational chunking alter boundaries?
6. Does new reproduce legacy on tractable domains?
7. Before/after runtime and memory?
8. Is execution restartable?
9. Are incompatible checkpoints rejected?
10. Is parallelism deterministic?
11. Did the full `filter=none` run complete?
12. Runtime and peak memory?
13. GAIN/LOSS under `none`?
14. Comparison with 2 GAIN / 688 LOSS?
15. DER bp-Jaccard?
16. Direction/effect concordance?
17. How many none-only DERs are low signal?
18. Are none-only low-signal regions replicate-consistent where assessable?
19. Does removing the filter alter LOSS dominance?
20. Do all prior critical regressions pass?
21. Does any BigWig path call edgeR?
22. Final PASS / CONDITIONAL PASS / FAIL?
23. If PASS, can v0.4.x infrastructure work close and the project proceed?

## 21. Execution order

```text
profile failure
→ chromosome dispatcher
→ contiguous blocks
→ memory bounding
→ checkpoint/resume
→ deterministic parallelism
→ legacy equivalence
→ full filter=none
→ none vs 0.6
→ low-signal characterization
→ scaling report
→ regressions
→ final acceptance report
```

If legacy equivalence fails, stop before interpreting the full run.

## 22. Agent instruction

Implement **v0.4.2 scalability closure only**.

Do not start v0.5. Do not implement NB-FPOP or DESeq2. Do not benchmark competitors. Do not claim BAM/BigWig equivalence. Do not alter frozen scientific parameters merely to make `filter=none` finish.

Prefer chromosome-wise, memory-bounded, restartable segmentation with deterministic chromosome-level parallelism.

After the full GSE280947 `filter=none` analysis and all acceptance gates are complete, generate:

`reports/v0.4.2_final_acceptance_report.md`

Then stop for review.
