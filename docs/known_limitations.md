# Known limitations (v0.2)

EpiSegDiff v0.2 is a research prototype. This file is the explicit
non-goal and caveat list.

## Scientific

- Not claimed: novelty of FPOP, epigenomic changepoint detection, or
  DMR segmentation.
- Not claimed: guaranteed \(O(N \log N)\) for generalized engines.
- BH q-values after adaptive segmentation are **not** guaranteed
  region-level FDR.
- TMM / library-size normalization can be invalid under global
  biological shifts. Sensitivity analysis is required.
- Two-condition designs only in the MVP.
- Segmentation of \(Z_i\) uses a Gaussian / Huber working model. That
  is an approximation, not the NB sampling model.
- Condition-specific domains used for topology are a support heuristic,
  not a second differential test.
- Penalty paths are not automatically valid trees.

## Implementation (current)

- Custom NB-FPOP is **not implemented** and is blocked (Phase 14).
- Spike-in normalization requires dedicated `reference_bins` (not
  arbitrary factors). `background` (lower-quantile median-of-ratios)
  is implemented; it fails when GAIN and LOSS are both large (Task A).
- TMM is the default. Library-size is not acceptable under one-sided
  composition shifts. See `reports/v0.2.1_Normalization_Stress_Test_Report.md`.
- `BIDIRECTIONAL_SHIFT` is not a primary topology class.
- DESeq2 size factors are optional and require DESeq2.
- Real BAM datasets and competitor binaries are not bundled.
- ONT methylation adapter is not implemented.
- Multi-sample batch correction is not implemented.
- BigWig writing requires `rtracklayer`.
- Selective inference / sample splitting is documented as optional
  research and does not block the prototype.

## Operational

- Chromosomes, samples, or competitor failures are never silently
  dropped; they are logged.
- Results that depend on undocumented manual edits are invalid.
