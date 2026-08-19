# EpiSegDiff v0.3 --- Real-Data Pilot: JUN CUT&RUN under JNK Inhibition

**Status:** Agent-executable validation plan\
**Primary dataset:** GEO **GSE280947**\
**Primary contrast:** DMSO vs JNK inhibitor (SP600125), JUN/c-Jun
CUT&RUN\
**Genome:** mouse GRCm38/mm10\
**Primary goal:** Validate EpiSegDiff **Level I GAIN/LOSS** on a real
narrow-factor CUT&RUN cohort.\
**Secondary goal:** Evaluate boundary behavior and Level II events
experimentally without making Level II a release gate.

------------------------------------------------------------------------

# 1. Why this dataset was selected

GSE280947 is selected as the first real-data pilot because it closely
matches the v0.2.2 readiness envelope.

## Evidence from GEO / publication

The dataset contains:

``` text
DMSO + IgG       n=2
DMSO + JUN       n=2
DMSO + H3K4me3   n=2

JNK inhibitor + IgG       n=2
JNK inhibitor + JUN       n=2
JNK inhibitor + H3K4me3   n=2
```

The primary JUN samples are:

``` text
GSM8607876  DMSO, JUN, replicate 1
GSM8607877  DMSO, JUN, replicate 2

GSM8607882  JNK inhibitor, JUN, replicate 1
GSM8607883  JNK inhibitor, JUN, replicate 2
```

Relevant controls:

``` text
GSM8607874  DMSO, IgG, replicate 1
GSM8607875  DMSO, IgG, replicate 2
GSM8607880  JNK inhibitor, IgG, replicate 1
GSM8607881  JNK inhibitor, IgG, replicate 2
```

Positive-control H3K4me3:

``` text
GSM8607878
GSM8607879
GSM8607884
GSM8607885
```

The experiment used Min6 mouse beta cells treated for 48 h with either
DMSO or 20 µM SP600125. CUT&RUN used a JUN/c-Jun antibody, paired-end
2×150-bp sequencing, approximately 30 million reads per sample, and E.
coli spike-in DNA. The original analysis used MACS2 narrow peaks against
IgG and IDR across the two biological replicates.

## Why it is particularly useful for EpiSegDiff

1.  JUN is a narrow transcription-factor binding assay, matching the
    first-pilot recommendation.
2.  The perturbation directly targets JNK signaling, creating a
    biologically plausible expectation of altered JUN occupancy.
3.  Two biological replicates exist per condition.
4.  Paired-end raw reads are available.
5.  IgG controls exist for both conditions.
6.  E. coli spike-in information is present, allowing a valuable
    normalization sensitivity analysis.
7.  Original MACS2/IDR peak outputs provide an external reference
    without being required by EpiSegDiff discovery.
8.  The published processing provides a reproducible conventional
    baseline.

**Important:** Published peaks/IDR are evaluation references only. They
must never be used to define EpiSegDiff candidate regions.

------------------------------------------------------------------------

# 2. Scientific questions

The pilot must answer:

1.  Can EpiSegDiff recover strong, reproducible JUN occupancy changes
    without predefined peaks?
2.  Are EpiSegDiff Level I calls stable to replicate removal?
3.  How sensitive are calls to TMM versus spike-in/reference
    normalization?
4.  How sensitive are calls to bin size and segmentation penalty?
5.  How much do EpiSegDiff calls overlap conventional MACS2/IDR-based
    signal changes?
6.  What characterizes:
    -   shared calls,
    -   conventional-only calls,
    -   EpiSegDiff-only calls?
7.  Do EpiSegDiff-only regions show reproducible raw-signal evidence?
8.  Does the v0.2.2 symmetric boundary algorithm behave plausibly on
    real narrow-factor data?
9.  Is there evidence that the current method is ready for a larger
    benchmark?

------------------------------------------------------------------------

# 3. Explicit hypotheses

## H1 --- Level I reproducibility

Strong EpiSegDiff GAIN/LOSS regions should be supported by both
biological replicates.

## H2 --- Perturbation consistency

Because the perturbation targets JNK signaling, a substantial subset of
JUN occupancy should change between DMSO and inhibitor-treated cells.

This is a biological expectation, **not** a requirement that all or most
JUN sites decrease.

SP600125 may have indirect/off-target consequences; do not use a
simplistic "all JUN must be LOSS" criterion.

## H3 --- Normalization

TMM should be acceptable if differential occupancy is not an extreme
genome-wide one-directional shift.

Spike-in-derived scaling provides an independent sensitivity analysis.

If TMM and spike-in normalization disagree substantially, investigate
rather than selecting the result with more discoveries.

## H4 --- Standard-method concordance

Strong conventional JUN peaks with clear replicate-consistent
differential signal should often be recovered by EpiSegDiff.

EpiSegDiff is not required to reproduce all conventional peaks because
its target is differential regions rather than peak presence.

------------------------------------------------------------------------

# 4. Dataset acquisition

Create:

``` text
data/GSE280947/
├── metadata/
├── fastq/
├── bam/
├── qc/
├── spikein/
├── conventional/
└── episegdiff/
```

Retrieve raw data from SRA for the eight JUN + IgG samples first.

H3K4me3 samples are optional QC/reference samples and must not delay the
primary pilot.

Create:

``` text
data/GSE280947/metadata/samples.tsv
```

with:

``` text
sample_id
gsm
condition
antibody
replicate
sra_run
fastq_1
fastq_2
bam
spikein_available
```

## Gate 0

Verify:

-   [ ] all four JUN samples exist;
-   [ ] all four IgG samples exist;
-   [ ] paired-end layout confirmed;
-   [ ] replicate/condition labels match GEO;
-   [ ] mm10/GRCm38 reference identified;
-   [ ] raw-data checksums recorded.

Generate:

``` text
reports/v0.3_GSE280947_Data_Acquisition_Report.md
```

------------------------------------------------------------------------

# 5. Alignment and preprocessing

Use a reproducible paired-end pipeline.

Recommended primary processing:

``` text
FASTQ
  ↓
adapter/quality QC
  ↓
Bowtie2 → mm10
  ↓
mapping-quality filtering
  ↓
paired-fragment BAM
  ↓
EpiSegDiff count matrix
```

Also align/quantify E. coli spike-in reads consistently with the
dataset's experimental design.

Do not copy the original authors' processed bigWigs as the source of
EpiSegDiff counts.

Processed bigWigs may be used only for cross-checking visualization.

Record:

``` text
raw read pairs
mapped read pairs
proper pairs
MAPQ-filtered pairs
duplicate fraction/status
fragment-length distribution
mouse-mapped fragments
E. coli-mapped fragments
spike-in scaling factor
```

Do not remove duplicates merely because another pipeline did so; choose
a policy, justify it, and run a sensitivity check if duplicate burden is
substantial.

------------------------------------------------------------------------

# 6. CUT&RUN QC

Generate per-sample:

``` text
fragment-size distribution
library complexity
mouse mapping rate
E. coli mapping fraction
IgG background
JUN enrichment
replicate correlation
genome-browser tracks
```

Evaluate replicate concordance before differential analysis.

## Pilot stop condition

Stop and report if:

-   one replicate is grossly discordant;
-   a sample has severe mapping failure;
-   JUN signal is not detectably enriched over IgG;
-   sample labels appear inconsistent.

Do not silently remove a replicate.

------------------------------------------------------------------------

# 7. Primary EpiSegDiff analysis

## Primary assay

``` text
DMSO JUN n=2
vs
JNK inhibitor JUN n=2
```

IgG is a QC/background reference, not a mandatory peak mask.

## Primary normalization

``` text
TMM
```

## Primary Level I output

``` text
GAIN
LOSS
```

where direction is defined explicitly as:

``` text
JNK inhibitor relative to DMSO
```

Therefore:

``` text
GAIN = higher JUN occupancy after inhibitor
LOSS = lower JUN occupancy after inhibitor
```

Do not reverse this convention in downstream plots.

------------------------------------------------------------------------

# 8. Normalization sensitivity

Run at least:

``` text
A. TMM
B. spike-in-derived normalization
```

If technically defensible, add:

``` text
C. IgG/background reference normalization
```

Do not make IgG subtraction the default observation.

For each method report:

``` text
normalization factors
GAIN count
LOSS count
total DER count
total DER bp
effect-size distribution
TMM-vs-spike-in overlap
direction agreement
```

Define pairwise region concordance using both:

``` text
interval overlap
effect-direction agreement
```

Flag regions where normalization changes the sign.

These regions require special review.

------------------------------------------------------------------------

# 9. Bin-size sensitivity

Primary:

``` text
50 bp
```

Sensitivity:

``` text
25 bp
100 bp
```

Optional:

``` text
200 bp
```

if runtime is acceptable.

For each bin size report:

``` text
DER count
GAIN/LOSS ratio
median DER width
strong-call overlap with primary run
boundary displacement
runtime
memory
```

Do not judge robustness solely by identical region coordinates.

Strong biological regions may shift by approximately one bin.

------------------------------------------------------------------------

# 10. Penalty sensitivity

Define a primary penalty from the current validated configuration.

Run a small grid around it, e.g.:

``` text
0.5×
1×
2×
```

or the appropriate parameterization for the current implementation.

Report:

``` text
segment count
tested candidate count
significant DER count
median width
strong-region recovery
fragmentation rate
merge rate
```

The goal is to identify stable conclusions, not maximize region count.

------------------------------------------------------------------------

# 11. Leave-One-Replicate-Out stability

This is mandatory.

With two replicates per condition, perform four LOO perturbations:

``` text
remove DMSO rep1
remove DMSO rep2
remove inhibitor rep1
remove inhibitor rep2
```

**Important statistical limitation:** each LOO run leaves one condition
with n=1. These runs are diagnostic stability analyses, not independent
inferential analyses and must not be interpreted as valid
replicate-aware differential tests equivalent to the full model.

If the software cannot support the asymmetric n=1/n=2 design safely,
implement a simpler diagnostic:

``` text
full-model DER
    ↓
measure per-replicate effect consistency inside each DER
```

For each full-model DER calculate:

``` text
replicate-pair direction consistency
effect consistency
LOO recovery if available
```

Define a pilot stability score:

\[ Stability(R)=
`\frac{\#\text{diagnostic runs recovering an overlapping same-direction region}}`{=tex}
{#`\text{valid diagnostic runs}`{=tex}} \]

Store:

``` text
region_id
full_effect
stability_score
direction_consistency
replicate_support
```

Do not use stability score as a replacement for FDR.

------------------------------------------------------------------------

# 12. Conventional reference analysis

Run a conventional pipeline independently:

``` text
JUN BAM
  ↓
MACS2 narrow peak calling using condition-matched IgG
  ↓
replicate consistency / IDR where appropriate
  ↓
union/consensus reference regions
  ↓
count fragments
  ↓
edgeR/DiffBind-style differential testing
```

This is the first pilot comparator, not the full competitor benchmark.

Record every parameter.

Do not optimize MACS2 specifically to favor or disfavor EpiSegDiff.

Use the authors' deposited IDR peaks as a secondary external reference
where useful.

------------------------------------------------------------------------

# 13. A/B/C region classification

Create three principal classes.

## Class A --- Shared

``` text
EpiSegDiff significant
AND
conventional differential reference significant
```

## Class B --- Conventional only

``` text
conventional significant
AND
no matching EpiSegDiff DER
```

## Class C --- EpiSegDiff only

``` text
EpiSegDiff significant
AND
no matching conventional differential region
```

Also retain:

``` text
direction disagreement
normalization-sensitive
parameter-sensitive
low-stability
```

as orthogonal flags.

Output:

``` text
results/region_comparison.tsv
```

------------------------------------------------------------------------

# 14. Analysis of Class A

Class A is the primary positive-validation set.

Report:

``` text
number
fraction of strong conventional changes recovered
effect correlation
direction agreement
boundary differences
width differences
stability
```

Generate browser-ready examples representing:

``` text
strong LOSS
strong GAIN
narrow change
wider local change
```

Do not cherry-pick only perfect examples.

Select examples by predeclared quantitative ranking.

------------------------------------------------------------------------

# 15. Analysis of Class B

For conventional-only regions determine whether failure is associated
with:

``` text
very narrow peak
weak effect
low counts
replicate inconsistency
segmentation penalty
bin size
candidate merging
normalization
FDR threshold
```

Produce a categorized failure table.

This class measures EpiSegDiff sensitivity limitations.

------------------------------------------------------------------------

# 16. Analysis of Class C

Class C is scientifically interesting but must be treated
conservatively.

For each EpiSegDiff-only region evaluate:

``` text
raw replicate signal
effect consistency
TMM/spike-in agreement
bin-size stability
penalty stability
IgG/background behavior
distance to conventional JUN peak
overlap with deposited IDR peaks
local boundary behavior
```

Classify each into:

``` text
C1 — strong reproducible non-peak-boundary differential signal
C2 — extension/boundary change around conventional peak
C3 — weak but reproducible candidate
C4 — normalization-sensitive
C5 — parameter-sensitive
C6 — likely artifact / unsupported
```

Do not describe C1--C3 as "novel biology" in this pilot.

Use:

> EpiSegDiff-specific candidate

until independent validation exists.

------------------------------------------------------------------------

# 17. Replicate stability metrics

For each DER compute:

``` text
within-condition replicate correlation
per-replicate normalized signal
condition effect consistency
LOO stability if valid
```

Rank confidence using separate dimensions rather than one opaque score.

Recommended fields:

``` text
statistical_FDR
absolute_effect
replicate_consistency
normalization_consistency
bin_consistency
penalty_consistency
```

Do not collapse them into a single final confidence score in v0.3.

------------------------------------------------------------------------

# 18. Browser-ready outputs

Generate BigWigs for:

``` text
each JUN replicate
each IgG replicate
condition mean JUN
TMM-normalized differential effect
spike-in-normalized differential effect
```

Generate BED tracks:

``` text
EpiSegDiff GAIN
EpiSegDiff LOSS
conventional differential
Class A
Class B
Class C
normalization-sensitive
```

Create an IGV/UCSC-ready track manifest.

Generate an HTML or Markdown review report with automatically selected
example loci.

------------------------------------------------------------------------

# 19. Automatic example selection

Avoid manual cherry-picking.

Select examples by rules such as:

``` text
top 5 Class A LOSS by |effect| with high stability
top 5 Class A GAIN by |effect| with high stability
top 10 Class C by FDR/effect with high cross-normalization stability
top 5 Class B by conventional effect
top 5 normalization-discordant regions
top 5 parameter-sensitive regions
```

If fewer regions exist, report all available regions.

------------------------------------------------------------------------

# 20. Level II exploratory analysis

Level II remains experimental.

For high-confidence Level I candidates only, run:

``` text
local condition-specific boundary inference
    ↓
EXPANSION
CONTRACTION
LEFT_SHIFT
RIGHT_SHIFT
SPLIT
MERGE
INTERNAL_GAIN
INTERNAL_LOSS
COMPLEX
```

Requirements:

-   use the symmetric v0.2.2 boundary method;
-   do not use conventional peaks to define EpiSegDiff boundaries;
-   report Level II separately;
-   do not use Level II success as the primary pilot PASS condition.

Quantify how many Level II calls remain stable across bin/penalty
choices.

------------------------------------------------------------------------

# 21. Primary pilot metrics

Report:

## Dataset/QC

``` text
mapping
fragment size
library complexity
replicate concordance
JUN/IgG enrichment
```

## Level I

``` text
GAIN count
LOSS count
region widths
effect sizes
FDR
```

## Robustness

``` text
TMM vs spike-in concordance
25/50/100-bp concordance
penalty concordance
replicate stability
```

## Conventional comparison

``` text
Class A count
Class B count
Class C count
direction agreement
effect correlation
```

## Boundary

``` text
boundary displacement across parameter settings
Level II stability
```

## Engineering

``` text
runtime
peak memory
disk usage
failed chromosomes
warnings
```

------------------------------------------------------------------------

# 22. No arbitrary success-by-region-count rule

Do not define success as:

``` text
EpiSegDiff finds more regions
```

or:

``` text
EpiSegDiff finds fewer regions
```

The pilot succeeds if the method produces interpretable, reproducible
differential regions with acceptable robustness and no major
statistical/engineering failure.

------------------------------------------------------------------------

# 23. Pilot PASS criteria

## Mandatory

-   [ ] all four JUN samples process successfully;
-   [ ] sample identity/QC acceptable;
-   [ ] full genome run completes reproducibly;
-   [ ] no systematic chromosome/index failure;
-   [ ] Level I GAIN/LOSS results generated;
-   [ ] TMM vs spike-in sensitivity completed;
-   [ ] bin sensitivity completed;
-   [ ] penalty sensitivity completed;
-   [ ] replicate stability evaluated;
-   [ ] conventional reference pipeline completed;
-   [ ] A/B/C classification completed;
-   [ ] browser-ready review outputs generated.

## Scientific readiness

PASS does **not** require a predetermined overlap percentage.

The report must instead determine whether:

1.  strong shared differential signals are recovered;
2.  major discrepancies can be explained;
3.  EpiSegDiff-only calls contain a reproducible supported subset;
4.  results are not dominated by normalization or parameter instability;
5.  Level I is useful enough to justify multi-dataset benchmarking.

------------------------------------------------------------------------

# 24. Stop conditions

Stop and report rather than silently modifying the analysis if:

-   sample metadata conflict with GEO;
-   one replicate is grossly discordant;
-   spike-in scaling is unusable;
-   TMM and spike-in produce globally inverted conclusions;
-   results change catastrophically across 25/50/100-bp bins;
-   penalty changes reverse most strong calls;
-   the conventional pipeline cannot be reproduced;
-   EpiSegDiff produces pathological genome-wide fragmentation;
-   memory/runtime makes genome-scale use impractical.

------------------------------------------------------------------------

# 25. Required outputs

Create:

``` text
results/
├── episegdiff_primary.tsv
├── episegdiff_primary.bed
├── normalization_sensitivity.tsv
├── bin_sensitivity.tsv
├── penalty_sensitivity.tsv
├── replicate_stability.tsv
├── conventional_differential.tsv
├── region_comparison.tsv
├── level2_experimental.tsv
└── run_manifest.tsv
```

------------------------------------------------------------------------

# 26. Required reports

Generate:

``` text
reports/v0.3_GSE280947_Data_Acquisition_Report.md
reports/v0.3_GSE280947_QC_Report.md
reports/v0.3_GSE280947_Level1_Report.md
reports/v0.3_GSE280947_Normalization_Sensitivity_Report.md
reports/v0.3_GSE280947_Parameter_Robustness_Report.md
reports/v0.3_GSE280947_Replicate_Stability_Report.md
reports/v0.3_GSE280947_Conventional_Comparison_Report.md
reports/v0.3_GSE280947_Level2_Exploratory_Report.md
reports/v0.3_GSE280947_Final_Pilot_Report.md
```

Each report must explicitly distinguish:

``` text
FACT
INFERENCE
UNRESOLVED
```

------------------------------------------------------------------------

# 27. Final Pilot Report Questions

`v0.3_GSE280947_Final_Pilot_Report.md` must answer:

1.  Did EpiSegDiff complete a reproducible genome-wide real-data run?
2.  How many GAIN and LOSS regions were found?
3.  Are strong calls replicate-supported?
4.  How concordant are TMM and spike-in normalization?
5.  Which calls change direction under normalization?
6.  How robust are strong calls to 25/50/100-bp bins?
7.  How robust are strong calls to segmentation penalty?
8.  What fraction of conventional strong differential regions are
    recovered?
9.  What are the principal causes of Class B conventional-only regions?
10. How many Class C EpiSegDiff-only regions have strong reproducibility
    support?
11. Are Class C regions enriched for peak-boundary extensions or other
    structures missed by consensus-peak analysis?
12. Does real-data boundary behavior remain symmetric/plausible?
13. Does Level II produce interpretable events, even if still
    experimental?
14. Are runtime and memory compatible with routine genome-scale
    analysis?
15. Is Level I ready for a multi-dataset competitor benchmark?
16. What specific weakness should be addressed next?
17. Is there any evidence that NB-FPOP is needed?

------------------------------------------------------------------------

# 28. Decision options

Allowed:

``` text
PASS
CONDITIONAL PASS
FAIL
```

## PASS

Proceed to:

``` text
second independent narrow dataset
    ↓
broad-mark pilot
    ↓
formal competitor benchmark
```

## CONDITIONAL PASS

Proceed to a second narrow dataset only, with specified limitations.

## FAIL

Return to the identified statistical or implementation failure before
expansion.

------------------------------------------------------------------------

# 29. NB-FPOP remains blocked

Do not implement NB-FPOP during this pilot.

Only reconsider it if real data identify a reproducible failure
specifically attributable to the current segmentation loss/model rather
than:

``` text
normalization
replicate variance
candidate generation
boundary refinement
binning
penalty selection
```

------------------------------------------------------------------------

# 30. Agent Start Instruction

Execute the **GSE280947 JUN CUT&RUN v0.3 real-data pilot**.

First retrieve and verify metadata and raw paired-end data for JUN and
IgG samples. Generate the acquisition report and stop if sample
identities or sequencing layouts do not match the specification.

Then perform QC and alignment to mm10/GRCm38, including E. coli spike-in
quantification.

Run EpiSegDiff Level I with:

``` text
contrast = JNK inhibitor vs DMSO
primary normalization = TMM
primary bin size = 50 bp
Level I = primary
Level II = experimental
```

Then perform:

``` text
spike-in normalization sensitivity
25/50/100-bp bin sensitivity
penalty sensitivity
replicate stability
independent MACS2/IDR/conventional differential reference
A/B/C region classification
browser-ready automatic example selection
```

Generate all required reports.

Stop after:

``` text
reports/v0.3_GSE280947_Final_Pilot_Report.md
```

Do not begin a second dataset, full competitor benchmark, broad-mark
validation, methylation validation, or NB-FPOP until this report is
reviewed.
