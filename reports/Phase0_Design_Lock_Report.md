# Phase 0 — Scientific design lock

**Objective.** Lock the count-assay statistical model, the separate
methylation observation model, and quantitative topology rules before
further engine work. The differential-effect target must be defined.

**Inputs.** Project specification v0.2; taskbook Phase 0 tasks 0.1–0.3.

**Methods.** Written specifications in `docs/`, not code-first guesses.

**Parameters.** Defaults recorded in `config/example.yaml` and
`default_config()`. All are user-visible.

## Locked model (count assays)

\[
Y_{ikr} \sim \mathrm{NB}(\mu_{ikr},\phi)
\]

\[
\log \mu_{ikr} = o_{kr} + \alpha_i + \beta_i X_k
\]

- \(o_{kr}\): normalization offset.
- \(\alpha_i\): baseline genomic signal.
- \(X_k\): 0 = control, 1 = treatment (MVP).
- \(\beta_i\): genomic condition effect (segmentation target).

**Default segmentation score:**
\(Z_i = \hat\beta_i / \mathrm{SE}(\hat\beta_i)\)
(`differential.score: standardized_effect`). Alternatives: `logFC`,
`shrunken_effect`.

**Dispersion:** `edgeR::estimateDisp` on bins **before** segmentation.
Not segment-specific. This breaks the circular “need segments for
dispersion / need dispersion for segments” loop.

**Zero counts:** retained. NA effect estimates become CVR gaps.

**Missing / masked bins:** `valid = FALSE`; no interpolation; independent
CVRs.

**Testing unit:** candidate interval after recount of **raw** fragments;
edgeR QL + BH. BH is nominal, not guaranteed selective FDR.

Methylation is a beta-binomial count-pair model and **never** enters the
NB / edgeR count likelihood.

Topology labels are predicates on intervals + configurable thresholds.
Ambiguity → `COMPLEX`.

## Acceptance

- [x] All symbols defined (`docs/statistical_model.md`).
- [x] Replicate handling defined (shared \(\alpha_i,\beta_i\); `~ condition`).
- [x] Dispersion strategy specified (pre-segmentation edgeR).
- [x] Zero-count and missing-bin behavior specified.
- [x] No circular segment-specific dispersion.
- [x] Methylation representation locked (`docs/assay_models.md`).
- [x] Topology events algorithmic (`docs/topology_definitions.md`).

## Unresolved scientific decisions

These are recorded, not hidden. They do **not** leave the
differential-effect target undefined.

1. **Whether \(Z_i\) or \(\hat\beta_i\) has better FDR/power** on real
   marks. Default is \(Z_i\); both remain configurable. Must be
   benchmarked (Phase 9).
2. **Penalty selection.** No universal \(\lambda\). Two candidate
   strategies are implemented; neither is declared superior.
3. **Condition-specific domain construction** for topology
   (`condition_mean_quantile`) is a support heuristic. A second method
   exists. Which matches visual “domains” on H3K27me3 is empirical.
4. **Global biological shifts** can invalidate TMM/library-size
   offsets. Spike-in / background methods are specified but not
   implemented in the MVP.
5. **Post-selection inference** (sample splitting, selective inference)
   is optional research and does not block the prototype.
6. **Multi-condition (>2) designs** are out of the MVP.
7. **ONT methylation** adapter is out of v0.2 core.
8. **Whether a custom NB-FPOP engine is justified** is explicitly
   unanswered until Phases 1–11 benchmarks exist
   (`docs/NB_FPOP_Experimental_Design.md`).
9. **Comparator operational status** (csaw / MACS3 / THOR / DSS / …)
   depends on the local environment; failed runs must be logged, not
   dropped.
10. **Coordinate of BAM fragment reconstruction** when `isize` is
    missing/zero: fall back to the sequenced read. Documented in
    `extract_fragments()`.

**Tests.** Document review against taskbook §2. Differential-effect
target is defined (section above).

**Results.** Design lock is sufficient to implement the MVP.

**Failures.** None that block proceeding. Items 1–10 remain empirical
or out of scope, not definitional holes.

**Limitations.** Design lock is not empirical validation.

**Files created.**

- `docs/statistical_model.md`
- `docs/assay_models.md`
- `docs/topology_definitions.md`
- `docs/benchmark_plan.md`
- `docs/known_limitations.md`
- `docs/NB_FPOP_Experimental_Design.md`

**Decision: PASS**

**Next action.** Phase 1 count-assay MVP (BAM QC, count matrix,
normalization), then Phases 2–8 in the published order.
