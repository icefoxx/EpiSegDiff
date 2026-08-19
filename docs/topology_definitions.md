# Topology and boundary-event definitions

**Status:** Phase 0 design lock  
**Rule:** every label is computed from intervals and numeric thresholds.
No visual judgement. Ambiguous mappings are `COMPLEX`.

All thresholds live in `topology.*` of the configuration.

---

## 1. Inputs

Classification uses three interval sets on one chromosome:

| Set | Meaning | How generated |
|-----|---------|---------------|
| \(C = \{C_u\}\) | Control supported domains | `build_condition_domains()` |
| \(T = \{T_v\}\) | Treatment supported domains | same, other condition |
| \(D = \{D_w\}\) | Significant differential intensity intervals | segmentation + testing |

Condition-specific domains are **not** the primary differential test.
They exist only to interpret geometry.

### Domain construction (locked default)

`topology.domain_method: condition_mean_threshold`

1. For each condition, compute the per-bin mean of
   library-size-normalized counts across replicates.
2. A bin is “supported” if that mean is \(\ge\) the configured
   chromosome-wide quantile (`domain_mean_quantile`, default 0.8).
3. Merge consecutive supported bins into domains.
4. Drop domains shorter than `segmentation.min_width` bins.

Alternative (optional): `per_condition_segmentation` runs Gaussian PELT
on the condition-mean track. The method used is stored on every event.

---

## 2. Interval measurements

For intervals \(A=[a_s,a_e)\), \(B=[b_s,b_e)\) on the same chromosome:

```text
width(A)              = a_e - a_s
intersection          = max(0, min(a_e,b_e) - max(a_s,b_s))
reciprocal_overlap    = min( intersection/width(A), intersection/width(B) )
                       if both widths > 0 else 0
containment(A in B)   = intersection / width(A)
distance              = 0 if intersection > 0
                       else min(|a_s-b_e|, |b_s-a_e|)   # gap length
left_shift(A→B)       = b_s - a_s     # >0: B starts to the right of A
right_shift(A→B)      = b_e - a_e
width_ratio(A→B)      = width(B) / width(A)
```

Coordinates are BED 0-based half-open.

---

## 3. Configurable thresholds

| Key | Default | Role |
|-----|---------|------|
| `min_reciprocal_overlap` | 0.2 | minimum overlap to consider a candidate match |
| `high_reciprocal_overlap` | 0.5 | “substantial” overlap |
| `max_match_distance` | 5000 bp | unmatched gap still eligible as a neighbor |
| `min_width_ratio_change` | 1.25 | expansion / contraction |
| `similar_width_ratio` | 1.25 | “similar width” if ratio \(\in [1/r, r]\) |
| `max_overlap_for_shift` | 0.85 | overlap above this is not called a shift |
| `min_boundary_shift` | 100 bp | shift / expansion displacement |
| `min_internal_gap_width` | 200 bp | gap that can separate split products |

---

## 4. Matching (deterministic)

Build a bipartite candidate graph:

- Edge \(C_u\)—\(T_v\) if they share a chromosome **and**
  (`reciprocal_overlap >= min_reciprocal_overlap`
  **or** (`intersection == 0` and `distance <= max_match_distance`
  and the two intervals do not have another higher-overlap partner
  blocking; distance-only edges are used only for shift detection)).

Priority for assigning matches (stable, sorted):

1. Higher `reciprocal_overlap`
2. Smaller `distance`
3. Smaller `|log(width_ratio)|`
4. Leftmost `C` start, then leftmost `T` start (tie break)

One-to-many and many-to-one are **allowed**. Greedy assignment keeps
every edge that does not contradict a higher-priority exclusive
one-to-one if both sides already have an overlapping partner with
strictly higher overlap **and** that overlap is `>= high_reciprocal_overlap`.
Otherwise multiple partners are retained (needed for SPLIT / MERGE).

The matching type for a connected component:

| Component | Type |
|-----------|------|
| 1 C, 1 T | `one_to_one` |
| 1 C, \(n\ge 2\) T | `one_to_many` (split candidate) |
| \(n\ge 2\) C, 1 T | `many_to_one` (merge candidate) |
| \(n\ge 2\) C, \(m\ge 2\) T | `many_to_many` → `COMPLEX` unless decomposable |
| 1 C, 0 T | `control_unmatched` |
| 0 C, 1 T | `treatment_unmatched` |

---

## 5. Event classes

Every event returns **exactly one primary class** plus optional
secondary annotations. The `classification_rule` string records the
branch taken.

Conflict resolution is the **ordered rule list below**. The first
matching rule wins.

### Rule 0 — `COMPLEX` pre-checks

- many-to-many component that cannot be split into disjoint
  one-to-many / many-to-one parts
- partners on different chromosomes (should not occur)
- zero-width interval

### Rule 1 — `SPLIT`

**Required intervals:** 1 control domain \(C\), \(n\ge 2\) treatment
domains \(T_{(1)} < T_{(2)} < \cdots\) in the same component.

| Check | Requirement |
|-------|-------------|
| overlap | each \(T_{(j)}\) has `containment(T in C) >= min_reciprocal_overlap` **or** `reciprocal_overlap >= min_reciprocal_overlap` |
| distance | consecutive \(T\) are separated by a gap \(\ge\) `min_internal_gap_width` |
| width | `sum(width(T_j)) <= width(C) * similar_width_ratio` (not a huge expansion) |
| direction | not required |

If gaps are smaller than `min_internal_gap_width`, do **not** call SPLIT
(fall through; often `COMPLEX` or `INTERNAL_LOSS`).

### Rule 2 — `MERGE`

**Required intervals:** \(n\ge 2\) control domains, 1 treatment domain.

Symmetric to SPLIT with C/T swapped.

### Rule 3 — `EXPANSION`

**Required intervals:** 1–1 match \((C,T)\).

| Check | Requirement |
|-------|-------------|
| overlap | `reciprocal_overlap >= min_reciprocal_overlap` and `containment(C in T) >= high_reciprocal_overlap` |
| width | `width(T) / width(C) >= min_width_ratio_change` |
| distance | overlapping |
| direction | at least one side extends by \(\ge\) `min_boundary_shift` (`left_shift < 0` and/or `right_shift > 0`) |

Secondary: `LEFT_EXPANSION`, `RIGHT_EXPANSION`, `BIDIRECTIONAL_EXPANSION`.

### Rule 4 — `CONTRACTION`

1–1 match.

| Check | Requirement |
|-------|-------------|
| overlap | `reciprocal_overlap >= min_reciprocal_overlap` and `containment(T in C) >= high_reciprocal_overlap` |
| width | `width(C) / width(T) >= min_width_ratio_change` |
| direction | at least one side recedes by \(\ge\) `min_boundary_shift` |

### Rule 3–5 — one-to-one geometry (v0.2.1 partition)

For a 1–1 match define \(\Delta_L = L_T - L_C\), \(\Delta_R = R_T - R_C\),
\(R_W = W_T / W_C\), displacement threshold \(\delta\) =
`min_boundary_shift`, width threshold \(\rho\) =
`min_width_ratio_change`.

States: left/right *out* / *in* / *stable* by comparing \(\Delta\) to
\(\pm\delta\). These nine cells are exclusive.

| \(\Delta_L,\Delta_R\) | Width | Primary class |
|----------------------|-------|----------------|
| out, out | \(R_W \ge \rho\) | `EXPANSION` (both ends) |
| in, in | \(R_W \le 1/\rho\) | `CONTRACTION` (both ends) |
| out, in | similar | `LEFT_SHIFT` |
| in, out | similar | `RIGHT_SHIFT` |
| out, stable | \(R_W \ge \rho\) | `EXPANSION` (left) |
| stable, out | \(R_W \ge \rho\) | `EXPANSION` (right) |
| in, stable | \(R_W \le 1/\rho\) | `CONTRACTION` (left) |
| stable, in | \(R_W \le 1/\rho\) | `CONTRACTION` (right) |
| stable, stable | similar | `INTERNAL_*` or Level I `GAIN`/`LOSS` |
| any other | — | `COMPLEX` |

`BIDIRECTIONAL_SHIFT` is **not** a primary class. Opposite-end motion
is expansion or contraction when the width rule holds; otherwise
`COMPLEX`.

### Rule 5 (removed) — `LEFT_SHIFT` / `RIGHT_SHIFT` (legacy overlap form)

1–1 match.

| Check | Requirement |
|-------|-------------|
| width | `width_ratio` \(\in [1/r, r]\) with \(r=\) `similar_width_ratio` |
| overlap | `reciprocal_overlap < high_reciprocal_overlap` (reduced exact overlap) |
| distance | overlapping **or** `distance <= max_match_distance` |
| direction | both boundaries move the **same** direction by \(\ge\) `min_boundary_shift` |

- `left_shift < -min_boundary_shift` and `right_shift < -min_boundary_shift` → `LEFT_SHIFT`
- `left_shift > min_boundary_shift` and `right_shift > min_boundary_shift` → `RIGHT_SHIFT`

### Rule 6 — `BIDIRECTIONAL_SHIFT`

1–1 match, similar width, both boundaries move by
\(\ge\) `min_boundary_shift` in **opposite** directions, but the width
ratio is still within `similar_width_ratio` (so not expansion/contraction).
This is unusual; most opposite-direction moves are expansion/contraction
and already caught. Remaining cases → `BIDIRECTIONAL_SHIFT`.

### Rule 7 — `INTERNAL_GAIN` / `INTERNAL_LOSS`

Uses differential intervals \(D_w\) relative to a **shared** domain.

Shared domain: 1–1 match with `reciprocal_overlap >= high_reciprocal_overlap`
and width ratio within `similar_width_ratio` (geometry essentially
unchanged).

| Class | Requirement |
|-------|-------------|
| `INTERNAL_GAIN` | a tested GAIN interval \(D\) with `containment(D in C∩T) >= high_reciprocal_overlap`, and \(D\) does not cover either outer boundary of \(C\cup T\) within `min_boundary_shift` |
| `INTERNAL_LOSS` | same with LOSS |

If a GAIN/LOSS interval **is** the whole matched domain, use Rule 8
instead (domain-level intensity).

### Rule 8 — `GAIN` / `LOSS`

Primary intensity change **without** a required geometry interpretation.

| Class | Required intervals | Rule |
|-------|--------------------|------|
| `GAIN` | unmatched \(T\), or 1–1 match with stable geometry, overlapping a GAIN \(D\) | treatment support appears or intensity increases without expansion/shift rules |
| `LOSS` | unmatched \(C\), or 1–1 match with stable geometry, overlapping a LOSS \(D\) | control support lost or intensity decreases |

Stable geometry: `reciprocal_overlap >= high_reciprocal_overlap` and
width ratio within `similar_width_ratio` and both
`|left_shift|`, `|right_shift| < min_boundary_shift`.

Unmatched treatment domain with no overlapping LOSS in control → `GAIN`.
Unmatched control domain with no overlapping GAIN in treatment → `LOSS`.

### Rule 9 — `COMPLEX`

Everything else, including:

- conflicting geometry + intensity that fails Rules 1–8
- overlapping expansion and shift criteria that were not resolved
  (should be rare because expansion is checked first)
- many-to-many leftovers

---

## 6. Intensity-only labels (Level I)

Independent of geometry, each tested region also carries:

```text
GAIN     if logFC > 0 and q <= fdr_threshold and |logFC| >= effect_threshold
LOSS     if logFC < 0 and q <= fdr_threshold and |logFC| >= effect_threshold
STABLE   otherwise
```

Level I is always written. Level II (this document) is additional.

---

## 7. Output fields

See taskbook table `boundary_events.tsv`. `classification_rule` is the
rule identifier (`RULE_SPLIT`, `RULE_EXPANSION`, …) plus the threshold
values used, so the label is reproducible from the stored rule.

---

## 8. Acceptance properties

- Every class is a predicate on intervals + numbers.
- Thresholds are configuration keys.
- Ambiguity → `COMPLEX`, never a forced prettier label.
- Deterministic given sort order defined in §4.
