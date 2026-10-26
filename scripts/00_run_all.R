# ==============================================================================
# Runs the complete analysis from the repository root:
#   Rscript scripts/00_run_all.R
# Each script can also be run on its own, in this order. Logs go to results/logs/.
# Tables 1 and 2 (clinical characteristics) are produced separately by
# notebooks/Tables_1_2_clinical_characteristics.ipynb (confidential input data).
# ==============================================================================

scripts <- c("01_sample_qc.R", "02_curation_and_composition.R", "03_nasopharyngeal_analysis.R",
             "04_gut_analysis.R", "05_cross_niche_model.R", "06_niche_comparison_and_networks.R",
             "07_supplementary_tables.R")
dir.create("results/logs", recursive = TRUE, showWarnings = FALSE)
for (s in scripts) {
  message("== ", s)
  log <- file.path("results/logs", sub("\\.R$", ".log", s))
  status <- system2(file.path(R.home("bin"), "Rscript"), file.path("scripts", s), stdout = log, stderr = log)
  if (status != 0) stop(s, " failed; see ", log)
}
# Versions of the packages used (written after loading them all)
suppressPackageStartupMessages(for (p in c("phyloseq", "mia", "ANCOMBC", "vegan", "SpiecEasi", "igraph", "nlme",
                                           "emmeans", "pROC", "broom", "dplyr", "tidyr", "ggplot2", "patchwork",
                                           "cowplot", "ggpubr", "writexl")) library(p, character.only = TRUE))
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
message("Done. Figures in figures/, supplementary tables in tables/supplementary/.")
