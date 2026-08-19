# Experimental NB-FPOP design note

**Status:** BLOCKED until Phases 1–11 of the taskbook pass.  
**Location for any future code:** `src/experimental/`  
**This file is a design note, not an implementation.**

A custom Negative-Binomial functional-pruning engine must not replace
the Gaussian/Huber MVP until the questions below have answers and a
benchmark shows a defined accuracy problem the custom engine solves.

---

## 1. What exact NB cost would be optimized?

Candidate (not implemented): for a segment \([s,t]\) with counts
\(Y_{ikr}\) and fixed offsets \(o_{kr}\), maximize the NB log-likelihood
under a two-group mean

\[
\log \mu_{kr} = o_{kr} + \alpha + \beta X_k
\]

with \(\alpha,\beta\) segment-constant. The segment cost is
\(-2\times\) that maximized log-likelihood.

An alternative “effect-only” cost would model the already-estimated
\(\hat\beta_i\) as NB-derived with known weights. That is closer to the
MVP and may not justify a new engine.

## 2. Is dispersion fixed, pooled, shrunken, or estimated per candidate?

**Required decision before coding:** dispersion must not be freely
estimated per candidate in a way that recreates

```text
need segments to estimate dispersion
and
need dispersion to identify segments
```

The only pruning-compatible default considered here is **dispersion
fixed** at the pre-estimated bin-level (or common) \(\hat\phi\) from
edgeR. Segment-specific \(\phi\) is out of scope even for the
experimental branch unless a separate paper-level derivation is added.

## 3. What parameters are segment-specific?

Allowed: \(\alpha,\beta\) (or a single mean if an offset-only model).  
Not allowed in the first experimental prototype: per-segment \(\phi\),
per-segment size factors.

## 4. What pruning conditions remain valid?

FPOP / PELT pruning requires a cost that satisfies the relevant
inequality (typically a non-negative per-observation cost with an
additive constant \(K\)). The NB two-parameter profile likelihood is
**not** automatically a 1-D quadratic. Functional pruning needs an
explicit derivation (parameter to prune over, envelope of piecewise
functions). Until that derivation exists, only PELT with \(K=0\) as an
empirical heuristic would be attempted — and that is **not** a
correctness proof.

## 5. Theoretical worst-case complexity

Without a proven prune, PELT/OP on NB costs is \(O(N^2)\) per
chromosome / CVR, times the cost of a 2-parameter NB MLE per candidate
segment. No \(O(N\log N)\) claim is permitted.

## 6. Expected empirical speed improvement

Unknown. The MVP already segments a 1-D score in essentially linear
time after pruning. NB-FPOP is unlikely to be faster. It is motivated
by **accuracy**, not speed.

## 7. What accuracy problem of the MVP is this intended to solve?

To be filled after Phase 3 simulation benchmarks. Hypotheses to test
first:

- Gaussian-on-\(Z\) over-segments low-count bins even after
  standardization.
- Huber still misses heavy NB tails at very low depth.
- Boundary bias when effects are large and variance is mean-dependent.

If those hypotheses fail (MVP already well calibrated), NB-FPOP should
not be written.

---

## Comparison required if implemented

```text
Gaussian effect-FPOP / PELT
robust effect segmentation
experimental NB-FPOP
```

Metrics: boundary accuracy, region power, empirical FDR, runtime,
memory.

Do not replace the stable engine until that comparison is stored under
`benchmark/`.
