#!/usr/bin/env Rscript
# GSE280947 JUN CUT&RUN Level I pilot.
# Refuses to invent counts from deposited bigWigs.

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE)
})

meta <- utils::read.table("data/GSE280947/metadata/samples.tsv",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
jun <- meta[meta$role == "primary", ]
missing <- jun$bam[!file.exists(jun$bam)]
if (length(missing)) {
  stop(
    "Primary JUN BAMs are not present. EpiSegDiff will not use deposited bigWigs as counts.\n",
    "Missing:\n  ", paste(missing, collapse = "\n  "),
    "\nSee reports/v0.3_GSE280947_Data_Acquisition_Report.md.",
    call. = FALSE
  )
}
cfg <- read_episegdiff_config("config/pilot_gse280947.yaml")
res <- run_episegdiff(cfg)
utils::write.table(res$tested, "results/episegdiff_primary.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote results/episegdiff_primary.tsv n=", nrow(res$tested))
