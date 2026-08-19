# Phase 0 — Project setup report

**Objective.** Create an installable R package skeleton that matches the
v0.2 specification layout, with a working test framework, a parseable
example configuration, and no custom NB-FPOP engine.

**Inputs.**

- `/Users/qzhang10/Downloads/EpiSegDiff_Project_Specification_v0.2.md`
- `/Users/qzhang10/Downloads/EpiSegDiff_Agent_Implementation_Taskbook_v0.2.md`

**Methods.** New R package `EpiSegDiff` 0.2.0 under the workspace root.
C++ sources implement only Gaussian PELT / interval / Huber helpers.

**Parameters.** Package version `0.2.0`. License MIT.

**Tests.**

- [x] Repository layout created (`R/`, `src/`, `inst/`, `tests/`,
      `benchmark/`, `docs/`, `vignettes/`, `scripts/`, `config/`,
      `reports/`).
- [x] Example configuration `config/example.yaml` and
      `inst/extdata/example.yaml` parse via `read_episegdiff_config()`.
- [x] Unit-test framework: `testthat` edition 3.
- [x] Local test command documented in `README.md`:
      `Rscript -e 'testthat::test_local(".")'`
- [x] No NB-FPOP implementation. `src/experimental/` is reserved and
      empty of engine code. `docs/NB_FPOP_Experimental_Design.md` is a
      blocked design note only.

**Results.** Package source is present and intended to load after
`Rcpp::compileAttributes()` + compile. Bioconductor packages are
Suggests, not hard load-time imports.

**Failures.** None for the setup gate itself. Install/load is verified
in the implementation session after compilation (see later reports if
the check fails).

**Limitations.** Real BAM datasets and competitor binaries are not
bundled. CI workflow is provided but not executed in this checkout.

**Files created.** See repository tree. Key setup files:

- `DESCRIPTION`, `NAMESPACE`, `LICENSE`, `README.md`
- `config/example.yaml`, `inst/schemas/config_schema.yaml`
- `.github/workflows/R-CMD-check.yaml`

**Decision: PASS**

**Next action.** Complete Phase 0 design lock (`docs/*.md` and
`reports/Phase0_Design_Lock_Report.md`), then implement the count-assay
MVP in the taskbook order.
