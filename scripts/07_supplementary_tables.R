# ==============================================================================
# 07. SUPPLEMENTARY TABLES 1-5 (Excel) FROM THE RESULTS OF SCRIPTS 01-06
# ==============================================================================
#   Supplementary Table 1  ANCOM-BC2, nasopharynx (global and pairwise)
#   Supplementary Table 2  ANCOM-BC2, gut (global and pairwise)
#   Supplementary Table 3  Network modules vs. clinical and environmental variables
#   Supplementary Table 4  Age-restricted case-control comparison  [Revision 1, new]
#   Supplementary Table 5  Sample quality control and sample flow   [Revision 1, new]
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})
source("R/functions.R")
make_dirs(DIR_TAB_SUPP)

rd <- function(f) read.csv(file.path(DIR_RESULTS, f), check.names = FALSE)
sci <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, sub("e-0?", "×10−", sprintf("%.2e", p)), sprintf("%.3f", p)))
pretty_group <- function(x) gsub("-/", "−/", gsub("Wheeze-", "Wheeze−", x))

# Writes one workbook; optional notes (table legend) go to a first "Notes" sheet
write_book <- function(sheets, file, notes = NULL) {
  if (!is.null(notes)) sheets <- c(list(Notes = data.frame(Note = notes)), sheets)
  writexl::write_xlsx(sheets, file)
  message("  table written: ", file)
}

# ------------------------------------------------------------------------------
# Supplementary Tables 1-2: ANCOM-BC2. One "Global test" row per robust genus followed
# by all its pairwise contrasts (A vs. CTRL, or A vs. B), ordered by p-value.
# ------------------------------------------------------------------------------
ancom_table <- function(niche_dir) {
  gl <- rd(file.path(niche_dir, "ancombc2", "Table_ANCOMBC2_Global_Robust.csv"))
  pw <- rd(file.path(niche_dir, "ancombc2", "Table_ANCOMBC2_Pairwise_Robust_Taxa.csv"))
  bind_rows(lapply(seq_len(nrow(gl)), function(i) {
    g <- gl[i, ]
    p <- pw %>% filter(Genus == g$taxon) %>% arrange(`p-value`)
    bind_rows(
      data.frame(`Taxon / Pairwise comparison` = paste0(g$taxon, " (Global test)"), W = round(g$W, 2),
                 LFC = NA, SE = NA, `p-value` = sci(g$p_val), `q-value (FDR)` = sci(g$q_val),
                 `Passed SS` = g$passed_ss, Robust = g$diff_robust_abn, check.names = FALSE),
      data.frame(`Taxon / Pairwise comparison` = paste0("   ", pretty_group(p$Group1), " vs. ", pretty_group(p$Group2)),
                 W = round(p$`W Statistic`, 2), LFC = round(p$`Log Fold Change (LFC)`, 2), SE = round(p$`Standard Error (SE)`, 2),
                 `p-value` = sci(p$`p-value`), `q-value (FDR)` = sci(p$`q-value (FDR)`),
                 `Passed SS` = p$`Passed Sensitivity (SS)`, Robust = NA, check.names = FALSE))
  }))
}
note_ancom <- c("ANCOM-BC2 (genus level, unrarefied counts; fixed effects: age + clinical group). Global rows: multi-group test for genera with robust global differences (q < 0.05 and passing the sensitivity analysis).",
                "Pairwise rows: all pairwise contrasts for that genus (LFC of the first vs. the second group); exploratory. W: test statistic; LFC: log fold change; SE: standard error; q: FDR-adjusted p; SS: sensitivity to pseudo-count.")
write_book(list(`Supplementary Table 1` = ancom_table("03_nasopharynx")),
           file.path(DIR_TAB_SUPP, "Supplementary_Table_1.xlsx"), note_ancom)
write_book(list(`Supplementary Table 2` = ancom_table("04_gut")),
           file.path(DIR_TAB_SUPP, "Supplementary_Table_2.xlsx"), note_ancom)

# ------------------------------------------------------------------------------
# Supplementary Table 3: network modules
# ------------------------------------------------------------------------------
write_book(list(`Supplementary Table 3` = rd("06_niches_and_networks/Module_Clinical_Associations.csv")),
           file.path(DIR_TAB_SUPP, "Supplementary_Table_3.xlsx"),
           c("Standardised multivariable linear-regression coefficients (beta) and FDR-adjusted p-values for module eigen-scores vs. clinical and environmental variables.",
                                            "Module numbering corresponds to the networks built on the final sample set (Revision 1)."))

# ------------------------------------------------------------------------------
# Supplementary Table 4: NPA case-control comparison restricted to cases within the age
# range of the controls  [Revision 1 - R1.1, R1.d]. A single sheet, without notes.
# ------------------------------------------------------------------------------
age_restricted <- bind_rows(
  rd("03_nasopharynx/age_restricted/AgeRestricted_Alpha.csv") %>% mutate(Analysis = "Alpha diversity", .before = 1),
  rd("03_nasopharynx/age_restricted/AgeRestricted_PERMANOVA.csv") %>% rename(Metric = Test) %>% mutate(Analysis = "PERMANOVA", .before = 1),
  rd("03_nasopharynx/age_restricted/AgeRestricted_ANCOMBC2.csv") %>% rename(Metric = Genus) %>% mutate(Analysis = "ANCOM-BC2 (age-adjusted)", .before = 1))
num <- vapply(age_restricted, is.numeric, TRUE)
age_restricted[num] <- lapply(age_restricted[num], signif, 4)
write_book(list(`Supplementary Table 4` = age_restricted), file.path(DIR_TAB_SUPP, "Supplementary_Table_4.xlsx"))

# ------------------------------------------------------------------------------
# Supplementary Table 5: sample QC  [Revision 1 - R1.5/a, R2.1]
# ------------------------------------------------------------------------------
write_book(list(
  Criteria = data.frame(Criterion = c("(1) Sequencing depth", "(2) Anatomical-site concordance", "(3) Exclusion"),
    Definition = c(">= 10,000 high-quality reads",
                   "Flagged if the mean Bray-Curtis dissimilarity (genus level) to samples of the other anatomical site is lower than to samples of its own site AND all 5 nearest neighbours belong to the other site",
                   "Flagged sample that is also near-identical (Bray-Curtis < 0.10) to a sample from a different participant (sample mix-up), or failing criterion (1)")),
  Sample_flow       = rd("01_sample_qc/Table_Sample_Flow.csv"),
  Site_concordance  = rd("01_sample_qc/Table_QC_Site_Concordance.csv"),
  Sequencing_depth  = rd("01_sample_qc/Table_QC_Sequencing_Depth.csv"),
  Rarefaction_slope = rd("01_sample_qc/Table_QC_Rarefaction_Slope.csv") %>% mutate(slope_per_1000 = signif(slope_per_1000, 3))),
  file.path(DIR_TAB_SUPP, "Supplementary_Table_5.xlsx"))
