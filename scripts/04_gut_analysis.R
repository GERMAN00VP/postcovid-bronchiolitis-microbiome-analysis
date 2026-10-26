# ==============================================================================
# 04. GUT (STOOL) MICROBIOTA ACROSS CLINICAL-VIROLOGICAL GROUPS
# ==============================================================================
#   Alpha and beta diversity (age-adjusted)             -> Supplementary Figure 2
#   Differential abundance (ANCOM-BC2, age-adjusted)     -> Supplementary Figure 3,
#                                                           Supplementary Table 2
#   Wheezing risk models within bronchiolitis cases      -> Supplementary Figure 4
# Same methods as script 03 (see the comments there).
# ==============================================================================

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(mia)
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(broom)
  library(pROC)
  library(vegan)
})
source("R/functions.R")
set.seed(123)

out_dir <- file.path(DIR_RESULTS, "04_gut")
dirs <- file.path(out_dir, c("alpha", "beta", "ancombc2", "models"))
names(dirs) <- c("alpha", "beta", "ancombc2", "models")
make_dirs(dirs)

# ------------------------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------------------------
pseq_raw_master  <- readRDS(file.path(DIR_CURATED, "phyloseq_RAW_curated_global.rds"))
pseq_rare_master <- readRDS(file.path(DIR_CURATED, "phyloseq_RAREFIED_global.rds"))

pseq_gut      <- subset_samples(pseq_raw_master,  Sample_Type == "GUT")
pseq_gut_rare <- subset_samples(pseq_rare_master, Sample_Type == "GUT")
metadata_base <- as(sample_data(pseq_gut), "data.frame")
tse_gut       <- convertFromPhyloseq(pseq_gut)

dist_gut  <- phyloseq::distance(pseq_gut_rare, method = "bray")
meta_beta <- metadata_base[labels(dist_gut), ]
# [Revision 1 - R1.5/a] 86 stool samples: BQLHSO91GUT reinstated, HULP71GUT excluded
cat("Gut samples:", nsamples(pseq_gut), "\n"); print(table(metadata_base$Condicion_Clinica))

# ------------------------------------------------------------------------------
# 2. Alpha diversity  [Revision 1 - R1.2, R2.3: joint model with spline age]
# ------------------------------------------------------------------------------
alpha_idx <- estimate_richness(pseq_gut_rare, measures = c("Observed", "Shannon", "Simpson"))
alpha_idx$Sample_ID <- rownames(alpha_idx)
df_alpha <- metadata_base %>% select(-any_of("Sample_ID")) %>% tibble::rownames_to_column("Sample_ID") %>%
  left_join(alpha_idx, by = "Sample_ID") %>% filter(!is.na(Condicion_Clinica))

alpha_tests <- alpha_age_adjusted_tests(df_alpha, c("Observed", "Shannon", "Simpson"),
                                        group = "Condicion_Clinica", age = "Age")
write.csv(alpha_tests$global, file.path(dirs["alpha"], "Table_Alpha_AgeAdjusted_Global.csv"), row.names = FALSE)
write.csv(alpha_tests$pairwise %>% arrange(Metric, p_FDR),
          file.path(dirs["alpha"], "Table_Alpha_AgeAdjusted_Pairwise_FDR.csv"), row.names = FALSE)

df_alpha_long <- df_alpha %>% select(Sample_ID, Condicion_Clinica, Observed, Shannon, Simpson) %>%
  pivot_longer(c(Observed, Shannon, Simpson), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Observed", "Shannon", "Simpson")))
write.csv(df_alpha_long %>% group_by(Metric, Condicion_Clinica) %>%
            summarise(n = n(), Median = median(Value), Q1 = quantile(Value, .25), Q3 = quantile(Value, .75),
                      Mean = mean(Value), SD = sd(Value), .groups = "drop"),
          file.path(dirs["alpha"], "Table_Alpha_Summary.csv"), row.names = FALSE)
stat_results <- alpha_tests$pairwise %>% filter(p_FDR < 0.05) %>%
  transmute(Metric = factor(Metric, levels = c("Observed", "Shannon", "Simpson")),
            group1 = Group1, group2 = Group2, p.adj = p_FDR, p.adj.signif = as.character(Signif)) %>%
  bracket_positions(df_alpha_long)

# ------------------------------------------------------------------------------
# 3. Beta diversity  [Revision 1 - R2.3: age-adjusted, FDR-corrected pairwise PERMANOVA]
# ------------------------------------------------------------------------------
permanova <- adonis2(dist_gut ~ Age + Condicion_Clinica, data = meta_beta, permutations = 999, by = "margin")
write.csv(as.data.frame(permanova), file.path(dirs["beta"], "Table_PERMANOVA_AgeAdjusted.csv"))
bdisp <- betadisper(dist_gut, meta_beta$Condicion_Clinica)
write.csv(as.data.frame(permutest(bdisp, permutations = 999)$tab), file.path(dirs["beta"], "Table_Betadisper.csv"))
write.csv(pairwise_permanova_age(dist_gut, meta_beta), file.path(dirs["beta"], "Table_Pairwise_PERMANOVA_AgeAdjusted_FDR.csv"),
          row.names = FALSE)

var_pct <- round(100 * bdisp$eig / sum(bdisp$eig), 1)
df_pcoa <- data.frame(PCoA1 = bdisp$vectors[, 1], PCoA2 = bdisp$vectors[, 2], Clinical_Group = meta_beta$Condicion_Clinica)
save_figure(plot_diversity_panel(df_alpha_long, stat_results, df_pcoa, var_pct), 2, width = 10, height = 9,
            supplementary = TRUE)

# ------------------------------------------------------------------------------
# 4. Differential abundance: ANCOM-BC2
# ------------------------------------------------------------------------------
ancom <- ancombc2(data = tse_gut, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "Age + Condicion_Clinica", p_adj_method = "fdr", group = "Condicion_Clinica",
                  dunnet = FALSE, pairwise = TRUE, global = TRUE, alpha = 0.05)
robust_taxa <- ancom$res_global$taxon[ancom$res_global$diff_robust_abn]
write.csv(ancom$res_global %>% filter(diff_robust_abn), file.path(dirs["ancombc2"], "Table_ANCOMBC2_Global_Robust.csv"), row.names = FALSE)
ancom_pairs <- ancom_pairwise_long(ancom$res_pair, robust_taxa)
write.csv(ancom_pairs, file.path(dirs["ancombc2"], "Table_ANCOMBC2_Pairwise_Robust_Taxa.csv"), row.names = FALSE)
cat("ANCOM-BC2 robust genera:", paste(robust_taxa, collapse = ", "), "\n")

df_bact <- ancom_log_abundances(ancom)
df_taxa <- df_bact[, robust_taxa, drop = FALSE]
colnames(df_taxa) <- gsub("-", "_", colnames(df_taxa))
df_taxa$Sample_ID <- rownames(df_taxa)

# ------------------------------------------------------------------------------
# 5. Wheezing risk models within bronchiolitis cases -> Supplementary Figure 4
# ------------------------------------------------------------------------------
df_model <- metadata_base %>% select(-any_of("Sample_ID")) %>% tibble::rownames_to_column("Sample_ID") %>%
  left_join(df_taxa, by = "Sample_ID") %>% filter(Bronchiolitis == "Yes") %>%
  mutate(Wheezing_Bin = ifelse(Wheezing.treatment == "Yes", 1, 0), Age = as.numeric(Age))
rm_gut <- fit_niche_risk_models(df_model, setdiff(colnames(df_taxa), "Sample_ID"))
write.csv(rm_gut$univariate, file.path(dirs["models"], "Table_Univariate_AgeAdjusted.csv"), row.names = FALSE)
write.csv(rm_gut$multivariable, file.path(dirs["models"], "Table_Multivariable.csv"), row.names = FALSE)
cat(sprintf("Gut multivariable model AUC = %.3f (%.3f-%.3f)\n", rm_gut$auc_ci[2], rm_gut$auc_ci[1], rm_gut$auc_ci[3]))
save_figure(plot_risk_panel(rm_gut), 4, width = 11, height = 10, supplementary = TRUE)

# ------------------------------------------------------------------------------
# 6. Supplementary Figure 3: age-adjusted, bias-corrected abundances of the robust genera
# ------------------------------------------------------------------------------
supp3 <- plot_ancom_boxplots(df_bact, metadata_base, robust_taxa, ancom_pairs, ncol = 2, step_increase = 0.13)
save_figure(supp3, 3, width = 11, height = 18, supplementary = TRUE)
