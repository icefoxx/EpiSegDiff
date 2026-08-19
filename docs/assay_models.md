# Assay observation models

**Status:** Phase 0 design lock  
**Principle:** One likelihood is not valid for every assay. Shared
software is the segmentation and boundary framework, not the
observation model.

---

## 1. Common metadata

Every sample must supply:

```text
sample_id
condition
replicate
assay
target
bam_or_methylation_path
control_path          # optional
spikein_status        # optional
library_metadata      # optional free text / JSON
```

`assay` is one of `cutnrun`, `chipseq`, `atac`, `methylation`.

---

## 2. CUT&RUN / ChIP-seq (count assays)

### Observation

Paired-end **fragment** counts per genomic bin.

| Channel | Default | Notes |
|---------|---------|-------|
| total fragments | yes | midpoint or overlap mode |
| short fragments | optional | CUT&RUN TF only; **not** a universal footprint |
| coverage / midpoint | optional | broad histone marks |

Do **not** interpret `<120 bp` as a universal TF footprint across
targets. Channel boundaries are user-configured when used.

Recommended starting assay for the MVP: paired-end CUT&RUN or ChIP-seq.

### Counting modes

- `fragment_midpoint`: bin containing \((start + end) // 2\) of the
  proper pair fragment.
- `fragment_overlap`: every bin overlapping \([frag\_start, frag\_end)\).

Unpaired / not-properly-paired records are excluded from fragment
counts and recorded in BAM QC. Duplicates are counted as present in the
BAM (flag 1024) and summarized; they are **not** silently removed unless
the user pre-filters the BAM.

### Coordinate convention

Internal and output intervals use **BED 0-based half-open**
\([start, end)\). BAM 1-based positions are converted on import.

---

## 3. ATAC-seq (count assay, separate adapter)

Preferred signal: **Tn5 insertion / cut-site counts**, not fragment
coverage.

Configurable Tn5 offsets (Buenrostro-style defaults, not mandatory):

```text
plus strand  : +4
minus strand : -5
```

Optional fragment-class channels (ATAC only):

```text
short : [0, 100]
mono  : [180, 247]
di    : [315, 473]
all insertions
```

These channels must not be forced onto CUT&RUN or ChIP-seq.

ATAC then uses the **same** differential-effect → segmentation →
testing interface as other count assays, with `counting_mode: insertion`.

---

## 4. DNA methylation (not Negative Binomial)

Methylation is a **proportion / count-pair** assay.

### Required representation

```text
chrom
position              # 0-based cytosine start; width 1
methylated_count      # M
total_count           # N  (coverage)
methylation_fraction  # M / N when N > 0, else NA
coverage              # alias of total_count
condition
replicate
sample
```

### Conceptual model (WGBS / RRBS)

\[
M_{ir} \sim \mathrm{BetaBinomial}(N_{ir}, p_{ir}, \phi)
\]

A useful two-group mean model on the methylation probability is

\[
\mathrm{logit}\,p_{ir} = \alpha_i + \beta_i X_k + \text{(optional offset)}
\]

The v0.2 adapter estimates a per-site (or per-bin) **delta-beta** and an
uncertainty, then forms a differential score for the **shared**
segmentation engine.

Locked rules:

1. Methylation **never** enters the count-assay NB / edgeR count model.
2. Sites with `total_count < methylation.min_coverage` are invalid
   (masked), not imputed.
3. Missing CpGs (no row) are absent, not zero methylation.
4. Low coverage must not be interpolated across.
5. Coordinate convention: cytosine start is 0-based; strand is retained
   if supplied, otherwise unstranded.
6. ONT modification calls are a **separate future adapter** (probability
   / count). Not implemented in v0.2 core.

### Differential methylation score

| Name | Definition |
|------|------------|
| `delta_beta` | \(\hat p_{\mathrm{trt}} - \hat p_{\mathrm{ctrl}}\) |
| `standardized_effect` | delta-beta / SE (method-of-moments or GLM) |

Region testing uses a methylation-specific model (`dss_style` or a
beta-binomial GLM approximation), never edgeR on methylated counts
alone.

### Import formats

- Tabular columns listed above (required)
- Optional: Bismark coverage (`chrom, start, end, pct, meth, unmeth`)
- Optional: bedGraph-like `chrom, start, end, beta` plus a coverage
  table if available

---

## 5. Binning and masks (all assays)

Supported bin sizes: 10, 25, 50, 100, 250, 1000, 5000, 10000 bp, plus
any user integer. No single default is biologically optimal.

Exclude or flag:

- blacklist
- low-mappability (if user supplies)
- assembly gaps (if user supplies)
- extra exclusions

Valid genomic intervals are stored as **contiguous valid regions**.
Mask edges are not differential biology.

---

## 6. What is shared vs what is not

| Layer | Shared? |
|-------|---------|
| Sample metadata schema | yes |
| Binning / CVR / masks | yes (units differ: bp bins vs CpG sites) |
| NB mean model | count assays only |
| Beta-binomial / DSS-style | methylation only |
| Gaussian/Huber segmentation of a differential score | yes |
| Topology classifier | yes |
| edgeR QL region test | count assays only |
