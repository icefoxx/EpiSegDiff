# EpiSegDiff v0.7 — Strict Methodological Review and v0.8 Decision Gate

**Review type:** Adversarial / strict methodological review  
**Reviewed artifact:** `v0.7_final_methodological_acceptance_report.md`  
**Formal verdict:** **CONDITIONAL PASS — method-development milestone only**  
**Purpose:** Determine what v0.7 actually established, distinguish parameter-level from formulation-level failures, and define the evidence required before v0.8 architecture freeze.

---

## 1. Executive verdict

v0.7 does **not** establish a universal narrow-to-broad segmentation method, validated automatic multi-scale engine, calibrated region-level FDR framework, or superiority over mature external methods.

The strongest defensible conclusion is:

> v0.7 successfully modularized EpiSegDiff and demonstrated that narrow-event recovery is possible with alternative score/segmentation configurations while strong broad-event performance is retained by the existing raw-effect + PELT configuration. However, it failed to produce a validated scale-adaptive decision rule, failed nested multi-scale recovery, retained same-direction merging and gradual-boundary limitations, and did not establish region-level FDR control.

This is substantial method-development progress, but not a final statistical architecture.

---

## 2. What v0.7 actually proved

### FACT 2.1 — Narrow-event failure is not structurally unavoidable

Reported results:

```text
standardized + PELT:
200 bp F1 = 0.47
500 bp F1 = 0.64

raw_effect + fused:
200 bp F1 = 0.62
```

Therefore the v0.6 narrow failure was at least partly configuration-dependent.

**Supported:** narrow-event recovery is achievable inside the modular framework.

**Not supported:** narrow events are now robustly solved.

### FACT 2.2 — Raw-effect PELT is the strongest demonstrated broad configuration

Reported:

```text
raw_effect:
F1 = 1 at 2–50 kb
one call
0 bp error
```

This supports PELT as the current broad baseline. It does not establish universal broad biological superiority because gradual boundaries degrade, real truth boundaries are unknown, and csaw/THOR are unavailable.

### FACT 2.3 — Fused lasso is useful but not universal

Reported:

```text
200 bp F1 = 0.62
500 bp raw-effect F1 = 0.06
50 kb F1 = 0.64 vs PELT 1
H3K27me3 chr19 = 4 GAIN / 12 LOSS
```

Therefore fused contains useful narrow-scale information but cannot replace PELT globally.

### FACT 2.4 — Current WBS implementation failed its intended role

Reported WBS F1=0 at 200–2000 bp. Its low merge rate cannot be interpreted as good separation when it does not call the events.

### FACT 2.5 — Automatic routing failed held-out validation

Reported:

```text
all 2 kb events routed "narrow"
held-out 2 kb F1 = 0.52
PELT = 1.00
500 bp test F1 = 0
```

Therefore v0.7 did not produce a validated automatic scale/engine selector.

### FACT 2.6 — True nested multi-scale recovery failed

For broad + embedded narrow truth:

```text
PELT  → merged
WBS   → merged
fused → merged
union → merged
```

This is a major failure, not a cosmetic benchmark defect.

### FACT 2.7 — Candidate union is not itself a multi-scale statistical solution

Union preserves the better caller in some isolated regimes but does not solve the 500 bp gap or nested events. It is an aggregation mechanism, not demonstrated multi-scale inference.

---

## 3. Central methodological concern

The combined results suggest that the current representation may be too restrictive because a locus is ultimately explained by one segmentation of one derived differential signal.

For nested truth, consider:

\[
D(x)=B(x)+L(x)+\epsilon(x)
\]

where \(B(x)\) is a broad differential component and \(L(x)\) is a local component.

A flat one-layer partition tries to explain \(D(x)\) with one segmentation. But the biological structure may require:

```text
broad component
+
local residual component
+
noise
```

This is an **inference from the failure pattern**, not a proven theorem. It must be tested against the alternative explanation that better parameterization/resolution is sufficient.

---

## 4. Parameter-level versus formulation-level failures

### 4.1 Likely parameter/implementation-sensitive

**WBS narrow failure:** penalty, interval sampling, thresholding, implementation, or score compatibility could matter. WBS as a general method is not falsified.

**Fused 500 bp instability:** likely sensitive to lambda, thresholding, discretization, noise, and scaling. However, broad real-data sign incoherence argues against global fused deployment.

**PELT narrow deficiency:** `200 bp × min_bins=3` imposes an obvious scale floor. Alternative score/minimum-length choices already improve narrow recovery.

### 4.2 Likely formulation-level

**Nested narrow-on-broad:** multiple segmentation families and union all fail the same conceptual test. This is the strongest evidence for a representation/objective limitation.

**Automatic routing:** selecting either a narrow or broad mode cannot naturally represent a locus containing both scales. The target may need to change from *choose one scale* to *decompose multiple scales*.

**Same-direction merging:** PELT gap-10 merge rate=0.90 versus fused=0.12 indicates a broad-coherence/local-separation tradeoff. A single global penalty may not resolve both without fragmentation.

**Gradual boundaries:** PELT F1=0.89 with ~1 kb error versus sharp F1=1/0 bp error suggests mismatch between a point-change model and transition-zone truth.

---

## 5. Strict statistical review

### 5.1 Frozen MAD detection is not region-level FDR inference

The factorial benchmark used a frozen MAD threshold, whereas `infer_regions` used limma/descriptive inference. These are not interchangeable.

`P(any false)=0` for some MAD-threshold null runs supports only the conservativeness of that tested detection rule.

It does **not** prove calibrated region-level FDR.

### 5.2 Same-data inference remains unresolved

Reported:

```text
same-data P(any false) = 0.06
pooled = 0.08
split = 0
cross-fit = 0
```

Interpretation is limited because P(any false) is not FDR, tested region sets can differ, split/cross-fit lose power, and cross-fit power was not adequately characterized.

**FACT:** severe v0.6 null inflation is reduced in tested v0.7 configurations.

**INFERENCE:** filtering/score choices improve practical null behavior.

**UNSUPPORTED:** formal region-level FDR is calibrated.

### 5.3 Pooled boundaries did not solve calibration

The tested pooled implementation produced 0.08 versus 0.06 for same-data. It cannot be claimed to improve calibration.

### 5.4 Naive sample splitting is too costly in current low-replicate settings

WBS recall reportedly falls 0.12 → 0.02. This shows that the tested naive split is poorly matched to low-replicate genomic studies; it does not invalidate sample splitting in principle.

Cross-fitting remains interesting but insufficiently validated.

---

## 6. Boundary review

Boundary accuracy must be decomposed into:

\[
P(\text{detected})
\]

and

\[
E[\text{boundary error}\mid\text{detected}]
\]

because fused can have exact 200 bp boundaries **when it fires**, while detection is only 0.62.

For gradual transitions, single-point boundary error is potentially inappropriate. Future truth should distinguish:

```text
core domain
left transition interval
right transition interval
background
```

and evaluate transition-zone overlap, core recovery, transition midpoint, and over/under-extension.

---

## 7. Why the two-mode proposal is only an engineering fallback

Current Pareto proposal:

```text
broad default:
raw + PELT min3

optional narrow overlay:
standardized + PELT min1
or fused
```

This is reasonable as a temporary benchmark/control, but not yet a satisfactory final method.

A user should not need to know the unknown event scale beforehand, and broad + embedded narrow truth cannot be represented simply by choosing one mode.

Retain this pair as a v0.8 control, not as the final architecture.

---

## 8. Recommended conceptual shift

The next question should change from:

> Which segmentation engine should be used at each scale?

to:

> Can differential signal be decomposed into scale-specific components and reconstructed hierarchically?

Conceptually:

\[
D(x)=D_{broad}(x)+D_{local}(x)+\epsilon(x)
\]

Potential workflow:

```text
high-resolution differential signal
              ↓
estimate broad component
              ↓
broad differential domains
              ↓
compute local residual
              ↓
detect narrow/intermediate events
              ↓
hierarchical reconciliation
              ↓
final region structure
```

This is fundamentally different from `PELT + fused + union`.

---

## 9. Candidate future directions

**Do not implement these yet; first perform the closure analyses in Section 12.**

### A — Broad-background residual decomposition

Estimate \(\hat B(x)\), then:

\[
R(x)=D(x)-\hat B(x)
\]

and detect local events on \(R(x)\).

Test edge artifacts, double counting, and identifiability.

### B — Multi-resolution decomposition

Candidates include wavelets, multiresolution fused lasso, trend filtering, scale-space approaches, or multiscale likelihood segmentation.

The essential property is additive/multiple-scale representation rather than one-scale routing.

### C — Hierarchical segmentation tree

Represent:

```text
chromosome
   ↓
broad domain
   ↓
subdomain
   ↓
local event
```

This raises an important design question: should output remain a flat BED-like DER set or become a hierarchical region graph?

### D — Multiscale penalized likelihood

Conceptually:

\[
\mathcal{L}(D;Y)+\lambda_B P(B)+\lambda_L P(L)
\]

with \(D=B+L\).

This is methodologically preferable to an empirical router if identifiability and calibration can be demonstrated.

---

## 10. Falsification tests for any redesigned objective

Any future architecture must pass:

1. **Isolated narrow:** recover 200–500 bp events.
2. **Isolated broad:** retain 5–100 kb coherence.
3. **Nested same-direction:** broad + local additional same-direction effect.
4. **Nested opposite-direction:** broad background plus local reversal.
5. **Same-direction separated:** do not merge across true neutral gap.
6. **Broad gradual:** represent transition interval appropriately.
7. **Narrow at broad boundary:** avoid residualization artifacts/duplicate calls.
8. **Null broad background:** smoothing must not generate local false events.
9. **Heteroscedastic null:** remain controlled under local variance changes.
10. **Replicate perturbation:** real-data calls/boundaries remain reasonably stable.

---

## 11. Benchmark weaknesses still present in v0.7

### 11.1 Sparse width grid

Primary widths `{200,500,2000,10000,50000}` leave important gaps. Future core grid should include approximately:

```text
100, 200, 300, 500, 750, 1k, 2k, 5k, 10k, 25k, 50k, 100k
```

### 11.2 Primary effect/noise regime is too dominant

Primary Δ=1.5/noise=0.2, while effect/noise exploration used n=10. Final architecture selection requires replicated effect × noise × width testing.

### 11.3 Main factorial still uses 200 bp bins

This confounds biological event width with observation resolution. A 200 bp event represented by one 200 bp bin is not equivalent to four 50 bp observations.

### 11.4 ATAC generalization is incomplete

chr1+chr19 execution establishes technical feasibility, not genome-wide biological validation.

### 11.5 External comparator gap remains

csaw and THOR are still unavailable. This blocks mature comparative/superiority claims.

---

## 12. Required closure analyses before v0.8 architecture freeze

### R1 — Resolution × narrow-width matrix

Run:

```text
resolution:
50, 100, 200 bp

event width:
100, 200, 300, 500, 750, 1000, 2000 bp
```

with the best current PELT and fused configurations.

**Purpose:** separate observation-resolution failure from segmentation-objective failure.

### R2 — Same-direction neutral-gap curve

Vary:

```text
0
1 bin
2 bins
3 bins
5 bins
10 bins
20 bins
```

across multiple widths/effects.

Estimate the separation phase boundary for PELT and fused.

### R3 — Nested broad × local effect matrix

Systematically vary broad effect \(B\) and local additional effect \(L\), including same- and opposite-direction nested events.

### R4 — Transition-zone gradual-boundary benchmark

Replace point truth for ramp events with explicit transition intervals.

### R5 — Conditional/unconditional boundary decomposition

For every engine report:

```text
detection probability
boundary error conditional on detection
overall reconstruction performance
```

These five analyses are mandatory before deciding whether v0.8 is tuning or redesign.

---

## 13. Optional reconstruction loss

Preserve all individual metrics. If a composite is needed, define only after freezing weights:

\[
L =
w_mL_{miss}
+w_bL_{boundary}
+w_gL_{merge}
+w_sL_{split}
+w_fL_{false}
\]

Prefer Pareto analysis if defensible weights cannot be specified.

Never hide component failures behind one score.

---

## 14. Statistical architecture requirement

Any future method must specify separately:

1. candidate discovery;
2. boundary estimation;
3. effect estimation;
4. hypothesis testing;
5. multiple-testing correction.

It must also state the inferential target explicitly:

```text
per-region type-I error?
FDR over selected regions?
false differential bp?
family-wise P(any false region)?
```

These are not interchangeable.

---

## 15. Manuscript readiness

### Ready for method-development reporting

- modular score/segmentation benchmarking;
- event-geometry-dependent regime map;
- failure of simplistic automatic routing;
- adaptive-inference calibration challenges;
- distinction between detection and FDR inference.

### Not ready for central claims that

- EpiSegDiff is universally scale-adaptive;
- EpiSegDiff outperforms established differential-region methods;
- region-level FDR is calibrated;
- one validated model automatically solves narrow and broad events;
- nested/hierarchical differential structure is resolved.

---

## 16. Strict gate reassessment

| Gate | v0.7 report | Strict review |
|---|---|---|
| A Modular correctness | PASS | **PASS** |
| B Narrow recovery | PASS if overlay | **PARTIAL — recoverable, not unified** |
| C Broad retention | PASS | **PASS within tested PELT regimes** |
| D Boundary accuracy | PASS ≥2 kb | **PARTIAL — sharp broad strong; gradual/narrow incomplete** |
| E Merge/split | PARTIAL | **PARTIAL / BLOCKER** |
| F True multi-scale / nested | FAIL | **FAIL — major blocker** |
| G Calibration | PARTIAL | **PARTIAL — no region-FDR proof** |
| H Selection-bias mitigation | PARTIAL | **PARTIAL — no satisfactory power/calibration solution** |
| I Real-data stability | PARTIAL | **PARTIAL** |
| J csaw/THOR | NOT_AVAILABLE | **OPEN BLOCKER for comparative claims** |
| K Regression | PASS | **PASS** |

---

## 17. Strict claims matrix

| Claim | Status |
|---|---|
| Modular framework can host multiple engines | **SUPPORTED** |
| Narrow events are recoverable in some configurations | **SUPPORTED** |
| Broad sharp synthetic events are strongly recovered by raw+PELT | **SUPPORTED** |
| One universal segmentation engine works across scales | **UNSUPPORTED** |
| Current union is genuinely multi-scale | **UNSUPPORTED** |
| Automatic routing works | **REFUTED in tested held-out benchmark** |
| Nested narrow-on-broad is recovered | **REFUTED in tested benchmark** |
| PELT is universally best for broad biology | **UNSUPPORTED** |
| Fused is universally best for narrow events | **PARTIAL / regime-dependent** |
| Current WBS is useful | **UNSUPPORTED** |
| Gradual boundaries are accurately localized | **UNSUPPORTED** |
| Region-level FDR is calibrated | **UNSUPPORTED** |
| Same-data selection bias is solved | **UNSUPPORTED** |
| Cross-fitting is production-ready | **UNSUPPORTED** |
| ATAC generalization is established | **UNSUPPORTED** |
| Superior to csaw/THOR | **UNSUPPORTED** |
| v0.7 provides a useful regime map | **SUPPORTED** |

---

## 18. Final methodological interpretation

The combined evidence supports:

> **v0.7 demonstrates that the narrow-versus-broad problem is not reliably solved by selecting one existing segmentation engine or naïvely unioning narrow and broad callers.**

The failure pattern is consistent with:

\[
\boxed{\text{single-layer segmentation} \neq \text{hierarchical multi-scale differential structure}}
\]

This remains an inference, not a theorem.

Supporting observations:

- PELT strong broad / weak narrow;
- fused useful narrow / weaker broad;
- union fails nested structure;
- automatic routing fails held-out simulations;
- same-direction neighboring events remain difficult;
- gradual boundaries expose point-change mismatch.

**Alternative explanation:** better tuning could close some gaps, especially WBS, fused 500 bp behavior, and PELT narrow behavior. It is less convincing as a complete explanation for the nested failure shared across engines.

---

## 19. Decision

### v0.7

**CONDITIONAL PASS — method-development milestone only.**

Do not freeze the narrow-overlay + broad-PELT pair as the final scientific architecture.

### v0.8

**NOT YET AUTHORIZED FOR FULL IMPLEMENTATION.**

First complete:

```text
R1 resolution × narrow-width matrix
R2 same-direction gap curve
R3 nested broad × local effect matrix
R4 transition-zone gradual-boundary benchmark
R5 conditional/unconditional boundary decomposition
```

Then choose:

**Path A — targeted engine refinement**  
Only if failures largely disappear after fair resolution/parameter correction.

**Path B — hierarchical/multiresolution objective redesign**  
If nested, same-direction separation, and scale-conflict failures persist.

---

## 20. Redesign decision criterion

Proceed to hierarchical objective redesign if, after fair parameter/resolution optimization:

1. no single engine achieves acceptable narrow and broad Pareto performance;
2. nested events remain systematically merged;
3. improving same-direction separation materially fragments broad events;
4. automatic one-scale routing remains unstable;
5. gradual-boundary performance remains structurally poor under point-change models.

If these conditions hold, stop tuning the flat segmentation architecture.

The next method should explicitly model something like:

\[
D(x)=B(x)+L(x)+\epsilon(x)
\]

or an equivalent hierarchical/multiresolution formulation.

---

## 21. Final reviewer statement

v0.7 is scientifically useful precisely because several negative results were retained rather than hidden.

The project has moved beyond:

> Can EpiSegDiff detect differential regions?

The central question is now:

> **What statistical representation of genomic differential signal can simultaneously preserve broad domains, local narrow changes, nested structure, accurate boundaries, and valid inference?**

v0.7 does not answer that question yet. It provides enough evidence to justify answering it next.

**Final strict-review verdict: CONDITIONAL PASS.**

**Immediate next action: focused closure analysis before v0.8 architecture freeze.**
