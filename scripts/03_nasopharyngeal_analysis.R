# ==============================================================================
# 03. NASOPHARYNGEAL (NPA) MICROBIOTA ACROSS CLINICAL-VIROLOGICAL GROUPS
# ==============================================================================
#   Alpha and beta diversity (age-adjusted)             -> Figure 2
#   Cases restricted to the age range of the controls    -> Supplementary Table 4
#   Differential abundance (ANCOM-BC2, age-adjusted)     -> Figure 3, Supplementary Table 1
#   Wheezing risk models within bronchiolitis cases      -> Figure 4
#
# The order of the steps that draw random numbers (permutations, plot jitter, imputation
# noise) is the one used for the results reported in the manuscript; the age-restricted
# analysis runs inside with_preserved_seed() so that it does not alter it.
# ==============================================================================

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(mia)          # loaded before the tidyverse so that dplyr verbs are not masked
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

out_dir <- file.path(DIR_RESULTS, "03_nasopharynx")
dirs <- file.path(out_dir, c("alpha", "beta", "ancombc2", "models", "age_restricted"))
names(dirs) <- c("alpha", "beta", "ancombc2", "models", "age_restricted")
make_dirs(dirs)

# ------------------------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------------------------
pseq_raw_master  <- readRDS(file.path(DIR_CURATED, "phyloseq_RAW_curated_global.rds"))
pseq_rare_master <- readRDS(file.path(DIR_CURATED, "phyloseq_RAREFIED_global.rds"))

pseq_npa      <- subset_samples(pseq_raw_master,  Sample_Type == "NPA")
pseq_npa_rare <- subset_samples(pseq_rare_master, Sample_Type == "NPA")
metadata_base <- as(sample_data(pseq_npa), "data.frame")
tse_npa       <- convertFromPhyloseq(pseq_npa)

dist_npa  <- phyloseq::distance(pseq_npa_rare, method = "bray")
meta_beta <- metadata_base[labels(dist_npa), ]
cat("NPA samples:", nsamples(pseq_npa), "\n"); print(table(metadata_base$Condicion_Clinica))

# ------------------------------------------------------------------------------
# 2. Alpha diversity (rarefied ASV table)
# ------------------------------------------------------------------------------
alpha_idx <- estimate_richness(pseq_npa_rare, measures = c("Observed", "Shannon", "Simpson"))
alpha_idx$Sample_ID <- rownames(alpha_idx)
df_alpha <- metadata_base %>% select(-any_of("Sample_ID")) %>% tibble::rownames_to_column("Sample_ID") %>%
  left_join(alpha_idx, by = "Sample_ID") %>% filter(!is.na(Condicion_Clinica))

# [Revision 1 - R1.2, R2.3] Joint model  index ~ group + ns(Age, 3 df)  with FDR-corrected
# pairwise contrasts, replacing residuals of index ~ Age followed by Wilcoxon tests.
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

# Figure 2a shows observed (not residualised) values; brackets only for FDR-significant contrasts
stat_results <- alpha_tests$pairwise %>% filter(p_FDR < 0.05) %>%
  transmute(Metric = factor(Metric, levels = c("Observed", "Shannon", "Simpson")),
            group1 = Group1, group2 = Group2, p.adj = p_FDR, p.adj.signif = as.character(Signif)) %>%
  bracket_positions(df_alpha_long)

# ------------------------------------------------------------------------------
# 3. Beta diversity (Bray-Curtis on the rarefied ASV table)
# ------------------------------------------------------------------------------
permanova <- adonis2(dist_npa ~ Age + Condicion_Clinica, data = meta_beta, permutations = 999, by = "margin")
write.csv(as.data.frame(permanova), file.path(dirs["beta"], "Table_PERMANOVA_AgeAdjusted.csv"))
bdisp <- betadisper(dist_npa, meta_beta$Condicion_Clinica)
write.csv(as.data.frame(permutest(bdisp, permutations = 999)$tab), file.path(dirs["beta"], "Table_Betadisper.csv"))

# [Revision 1 - R2.3] Pairwise PERMANOVA adjusted for age and truly BH-FDR corrected
# (pairwiseAdonis::pairwise.adonis2 ignored the requested correction in the submission)
write.csv(pairwise_permanova_age(dist_npa, meta_beta), file.path(dirs["beta"], "Table_Pairwise_PERMANOVA_AgeAdjusted_FDR.csv"),
          row.names = FALSE)

var_pct <- round(100 * bdisp$eig / sum(bdisp$eig), 1)
df_pcoa <- data.frame(PCoA1 = bdisp$vectors[, 1], PCoA2 = bdisp$vectors[, 2], Clinical_Group = meta_beta$Condicion_Clinica)
save_figure(plot_diversity_panel(df_alpha_long, stat_results, df_pcoa, var_pct), 2, width = 10, height = 9)

# ------------------------------------------------------------------------------
# 4. Cases restricted to the age range of the controls (3.7-9.0 months)
#    [Revision 1 - R1.1, R1.d] -> Supplementary Table 4
# ------------------------------------------------------------------------------
with_preserved_seed({
  df_a <- df_alpha %>% mutate(Age = as.numeric(Age),
                              CTRL = factor(ifelse(Condicion_Clinica == "CTRL", "CTRL", "Bronchiolitis"),
                                            levels = c("CTRL", "Bronchiolitis")))
  meta_s <- meta_beta %>% mutate(Age = as.numeric(Age),
                                 CTRL = factor(ifelse(Condicion_Clinica == "CTRL", "CTRL", "Bronchiolitis"),
                                               levels = c("CTRL", "Bronchiolitis")))
  metrics <- c("Observed", "Shannon", "Simpson")
  rng <- range(df_a$Age[df_a$CTRL == "CTRL"])
  df_r <- df_a[df_a$Age >= rng[1] & df_a$Age <= rng[2], ]
  cat("Age-restricted subset (", round(rng[1], 2), "-", round(rng[2], 2), "months):",
      sum(df_r$CTRL == "CTRL"), "controls vs", sum(df_r$CTRL != "CTRL"), "cases\n")
  ar_alpha <- bind_rows(lapply(metrics, function(m) data.frame(
    Metric = m, n_CTRL = sum(df_r$CTRL == "CTRL"), n_cases = sum(df_r$CTRL != "CTRL"),
    Median_CTRL = median(df_r[[m]][df_r$CTRL == "CTRL"]), Median_cases = median(df_r[[m]][df_r$CTRL != "CTRL"]),
    p_Wilcoxon = wilcox.test(df_r[[m]] ~ df_r$CTRL, exact = FALSE)$p.value,
    p_lm_spline_age = summary(lm(df_r[[m]] ~ df_r$CTRL + ns(df_r$Age, df = 2)))$coefficients[2, 4])))
  ids_r <- df_r$Sample_ID
  d_r <- as.dist(as.matrix(dist_npa)[ids_r, ids_r]); m_r <- meta_s[ids_r, ]
  set.seed(123); a_r  <- adonis2(d_r ~ Age + CTRL, data = m_r, permutations = 999, by = "margin")
  set.seed(123); bd_r <- permutest(betadisper(d_r, m_r$CTRL), permutations = 999)
  ar_perm <- data.frame(Test = "PERMANOVA CTRL vs bronchiolitis (age-restricted, age-adjusted)",
                        R2 = a_r["CTRL", "R2"], p = a_r["CTRL", "Pr(>F)"], p_betadisper = bd_r$tab[1, "Pr(>F)"])
  tse_r <- convertFromPhyloseq(prune_samples(ids_r, pseq_npa))
  SummarizedExperiment::colData(tse_r)$CTRL <- m_r[colnames(tse_r), "CTRL"]
  SummarizedExperiment::colData(tse_r)$Age  <- m_r[colnames(tse_r), "Age"]
  anc_r <- ancombc2(data = tse_r, assay_name = "counts", tax_level = "Genus", fix_formula = "Age + CTRL",
                    p_adj_method = "fdr", group = "CTRL", global = FALSE, pairwise = FALSE, alpha = 0.05)
  ar_ancom <- anc_r$res %>%
    filter(taxon %in% c("Dolosigranulum", "Corynebacterium", "Haemophilus", "Burkholderia-Caballeronia-Paraburkholderia")) %>%
    transmute(Genus = taxon, LFC_Bronchiolitis_vs_CTRL = lfc_CTRLBronchiolitis, SE = se_CTRLBronchiolitis,
              p = p_CTRLBronchiolitis, q_FDR = q_CTRLBronchiolitis)
  write.csv(ar_alpha, file.path(dirs["age_restricted"], "AgeRestricted_Alpha.csv"), row.names = FALSE)
  write.csv(ar_perm,  file.path(dirs["age_restricted"], "AgeRestricted_PERMANOVA.csv"), row.names = FALSE)
  write.csv(ar_ancom, file.path(dirs["age_restricted"], "AgeRestricted_ANCOMBC2.csv"), row.names = FALSE)
  print(ar_alpha); print(ar_perm); print(ar_ancom)
})

# ------------------------------------------------------------------------------
# 5. Differential abundance: ANCOM-BC2 (genus level, unrarefied counts, age-adjusted)
# ------------------------------------------------------------------------------
ancom <- ancombc2(data = tse_npa, assay_name = "counts", tax_level = "Genus",
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
# 6. Wheezing risk models within bronchiolitis cases -> Figure 4
# ------------------------------------------------------------------------------
df_model <- metadata_base %>% select(-any_of("Sample_ID")) %>% tibble::rownames_to_column("Sample_ID") %>%
  left_join(df_taxa, by = "Sample_ID") %>% filter(Bronchiolitis == "Yes") %>%
  mutate(Wheezing_Bin = ifelse(Wheezing.treatment == "Yes", 1, 0), Age = as.numeric(Age))
rm_npa <- fit_niche_risk_models(df_model, setdiff(colnames(df_taxa), "Sample_ID"))
write.csv(rm_npa$univariate, file.path(dirs["models"], "Table_Univariate_AgeAdjusted.csv"), row.names = FALSE)
write.csv(rm_npa$multivariable, file.path(dirs["models"], "Table_Multivariable.csv"), row.names = FALSE)
cat(sprintf("NPA multivariable model AUC = %.3f (%.3f-%.3f)\n", rm_npa$auc_ci[2], rm_npa$auc_ci[1], rm_npa$auc_ci[3]))
save_figure(plot_risk_panel(rm_npa), 4, width = 11, height = 10)

# ------------------------------------------------------------------------------
# 7. Figure 3: age-adjusted, bias-corrected abundances of the robust genera
# ------------------------------------------------------------------------------
short_name <- function(x) gsub("Burkholderia-Caballeronia-Paraburkholderia", "Burkholderia", x)
fig3 <- plot_ancom_boxplots(df_bact, metadata_base, robust_taxa, ancom_pairs, display = short_name,
                            ncol = 3, step_increase = 0.11)
save_figure(fig3, 3, width = 14, height = 10)
