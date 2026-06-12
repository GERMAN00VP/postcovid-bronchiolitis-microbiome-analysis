# ==============================================================================
# SECTION 1: LIBRARIES, CONFIGURATION & ENVIRONMENT
# ==============================================================================
library(ANCOMBC)
library(dplyr)
library(tidyr)
library(phyloseq)
library(mia)
library(broom)
library(pROC)
library(ggplot2)
library(ggpubr)
library(knitr)
library(vegan)
library(pairwiseAdonis)
library(WGCNA)
library(cluster)
library(igraph)
library(pheatmap)
library(stringr)
library(rstatix)

set.seed(123)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4) 

setwd("/mnt/usb/BQL/BQL_ANALYSIS")

# Configuración estricta de la nueva jerarquía de carpetas requerida (REANALISIS_3)
base_dir        <- "REANALISIS/NPA/"
dir_alpha       <- paste0(base_dir, "alpha/")
dir_beta        <- paste0(base_dir, "beta/")
dir_abundances  <- paste0(base_dir, "abundances/")
dir_models      <- paste0(base_dir, "models/") # Guardado directo sin subcarpetas
output_dir_wgcna <- paste0(base_dir, "WGCNA/") # Directorio de WGCNA

# Creación recursiva de todos los directorios del proyecto
for(path in c(dir_alpha, dir_beta, dir_abundances, dir_models, output_dir_wgcna)){
  if(!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# ==============================================================================
# SECTION 2: GLOBAL DATA PRE-PROCESSING & SUBSETTING
# ==============================================================================
data_input_dir <- "REANALISIS_2/curated_data_global/"

pseq_raw_master  <- readRDS(paste0(data_input_dir, "phyloseq_RAW_curated_global.rds"))
pseq_rare_master <- readRDS(paste0(data_input_dir, "phyloseq_RAREFIED_global.rds"))
metadata_master  <- readRDS(paste0(data_input_dir, "metadata_curated_global.rds"))

metadata_master <- metadata_master %>%
  mutate(
    Condicion_Clinica = factor(Condicion_Clinica, 
                               levels = c("CTRL", "RSV-/Wheeze-", "RSV-/Wheeze+", "RSV+/Wheeze-", "RSV+/Wheeze+")),
    Patient_ID = as.factor(Patient_ID),
    Sample_Type = as.factor(Sample_Type)
  )

sample_data(pseq_raw_master)  <- sample_data(metadata_master)
sample_data(pseq_rare_master) <- sample_data(metadata_master)

pseq_anf      <- subset_samples(pseq_raw_master, Sample_Type == "ANF" | Sample_Type == "NPA")
pseq_anf_rare <- subset_samples(pseq_rare_master, Sample_Type == "ANF" | Sample_Type == "NPA")

metadata_base <- as(sample_data(pseq_anf), "data.frame")
tse_anf       <- convertFromPhyloseq(pseq_anf)

dist_matrix_final <- phyloseq::distance(pseq_anf_rare, method = "bray")
muestras_comunes  <- labels(dist_matrix_final)
metadata_beta     <- metadata_base[muestras_comunes, ]

# ==============================================================================
# SECTION 3: ALPHA DIVERSITY (UNIVARIATE ANALYSIS - AGE ADJUSTED)
# ==============================================================================
alpha_indices <- estimate_richness(pseq_anf_rare, measures = c("Observed", "Shannon", "Simpson"))
alpha_indices$Sample_ID <- rownames(alpha_indices)

df_alpha_global <- metadata_base %>%
  dplyr::select(-any_of("Sample_ID")) %>% 
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::left_join(alpha_indices, by = "Sample_ID") %>%
  filter(!is.na(Condicion_Clinica) & Condicion_Clinica != "NA")

for (metric in c("Observed", "Shannon", "Simpson")) {
  df_alpha_global[[paste0(metric, "_adj")]] <- residuals(
    lm(as.numeric(df_alpha_global[[metric]]) ~ as.numeric(Age), data = df_alpha_global, na.action = na.exclude)
  ) + mean(df_alpha_global[[metric]], na.rm = TRUE)
}

df_alpha_long <- df_alpha_global %>%
  select(Sample_ID, Condicion_Clinica, Observed = Observed_adj, Shannon = Shannon_adj, Simpson = Simpson_adj) %>%
  pivot_longer(cols = c(Observed, Shannon, Simpson), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Observed", "Shannon", "Simpson")))

tabla_resumen_alpha <- df_alpha_long %>%
  group_by(Metric, Condicion_Clinica) %>%
  summarise(Mean = mean(Value), SD = sd(Value), .groups = 'drop')
write.csv(tabla_resumen_alpha, paste0(dir_alpha, "Table_Alpha_Diversity_Summary.csv"), row.names = FALSE)

tabla_comparaciones_alpha <- df_alpha_long %>%
  group_by(Metric) %>%
  wilcox_test(Value ~ Condicion_Clinica) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  select(Metric, Group1 = group1, Group2 = group2, `p-raw` = p, `p-adj (FDR)` = p.adj, Significancia = p.adj.signif) %>%
  arrange(Metric, `p-adj (FDR)`)

write.csv(tabla_comparaciones_alpha, paste0(dir_alpha, "Table_Alpha_Diversity_Wilcoxon_Pairwise_FDR.csv"), row.names = FALSE)

stat_results <- df_alpha_long %>%
  group_by(Metric) %>%
  wilcox_test(Value ~ Condicion_Clinica) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  filter(p.adj < 0.05)

plot_alpha <- ggplot(df_alpha_long, aes(x = Condicion_Clinica, y = Value, fill = Condicion_Clinica)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.3, color = "black", size = 0.9) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  theme_bw() + scale_fill_brewer(palette = "Set2") +
  labs(title = "Nasopharyngeal Microbiota Alpha Diversity Across Clinical Conditions",
       subtitle = "Age-adjusted residuals | Pairwise comparisons (Wilcoxon FDR)",
       x = "Clinical Group (RSV & Recurrent Wheezing Status)", y = "Alpha Diversity Index Value (Age-Corrected)") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        text = element_text(size = 12), plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, face = "italic"),
        strip.background = element_rect(fill = "#2c3e50"), strip.text = element_text(color = "white", face = "bold", size = 11),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.28)))

if(nrow(stat_results) > 0) {
  stat_results <- stat_results %>% add_y_position(step.increase = 0.12)
  plot_alpha <- plot_alpha + stat_pvalue_manual(stat_results, label = "p.adj.signif", hide.ns = TRUE, step.increase = 0.09, tip.length = 0.01)
}
ggsave(paste0(dir_alpha, "Plot_Alpha_Diversity_Clean.png"), plot = plot_alpha, width = 13, height = 6, dpi = 300)

# ==============================================================================
# SECTION 4: BETA DIVERSITY (MULTIVARIATE SPATIAL DYNAMICS)
# ==============================================================================
permanova_results <- adonis2(dist_matrix_final ~ Age + Condicion_Clinica, data = metadata_beta, permutations = 999, by = "margin")
write.csv(as.data.frame(permanova_results), paste0(dir_beta, "Table_Beta_Diversity_PERMANOVA.csv"))

beta_disp_final <- betadisper(dist_matrix_final, metadata_beta$Condicion_Clinica)
disp_test <- permutest(beta_disp_final, permutations = 999)
write.csv(as.data.frame(disp_test$tab), paste0(dir_beta, "Table_Beta_Diversity_Dispersion.csv"))

post_hoc_beta <- pairwise.adonis2(dist_matrix_final ~ Condicion_Clinica, data = metadata_beta, p.adjust.m = "fdr")
tabla_post_hoc_beta <- data.frame(Comparison = character(), F_Model = numeric(), R2 = numeric(), p_value_FDR = numeric())

for(pair in names(post_hoc_beta)[names(post_hoc_beta) != "parent_call"]) {
  res_pair <- post_hoc_beta[[pair]]
  tabla_post_hoc_beta <- rbind(tabla_post_hoc_beta, data.frame(Comparison = pair, F_Model = res_pair$F[1], R2 = res_pair$R2[1], p_value_FDR = res_pair$`Pr(>F)`[1]))
}
write.csv(tabla_post_hoc_beta, paste0(dir_beta, "Table_Beta_Diversity_Pairwise_PostHoc.csv"), row.names = FALSE)

porcentaje_variacion <- round(100 * (beta_disp_final$eig / sum(beta_disp_final$eig)), 1)
df_pcoa <- data.frame(PCoA1 = beta_disp_final$vectors[, 1], PCoA2 = beta_disp_final$vectors[, 2], Clinical_Group = metadata_beta$Condicion_Clinica)

plot_pcoa <- ggplot(df_pcoa, aes(x = PCoA1, y = PCoA2, color = Clinical_Group, fill = Clinical_Group)) +
  geom_point(size = 2.5, alpha = 0.6) + stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95, size = 0.5) +
  theme_bw() + scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2") +
  labs(title = "Principal Coordinates Analysis (PCoA) of Nasopharyngeal Microbiota",
       subtitle = paste0("PERMANOVA p = ", format.pval(permanova_results$`Pr(>F)`[2], digits = 3), " | Betadisper p = ", format.pval(disp_test$tab$`Pr(>F)`[1], digits = 3)),
       x = paste0("PCoA 1 (", porcentaje_variacion[1], "%)"), y = paste0("PCoA 2 (", porcentaje_variacion[2], "%)"), color = "Clinical Group", fill = "Clinical Group") +
  theme(text = element_text(size = 12), plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"), panel.grid.minor = element_blank())
ggsave(paste0(dir_beta, "Plot_Beta_Diversity_PCoA.png"), plot = plot_pcoa, width = 8.5, height = 6, dpi = 300)

# ==============================================================================
# SECTION 5: COMPOSITIONAL ANALYSIS (ANCOM-BC2) WITH BIAS-CORRECTED JITTER IMPUTATION
# ==============================================================================
output = ancombc2(
  data = tse_anf, assay_name = "counts", tax_level = "Genus",         
  fix_formula = "Age + Condicion_Clinica", p_adj_method = "fdr", group = "Condicion_Clinica",        
  dunnet = FALSE, pairwise = TRUE, global = TRUE, alpha = 0.05
)

df_global <- output$res_global
tabla_ancom_significativa <- df_global %>%
  filter(diff_abn == TRUE) %>%
  mutate(W = round(W, 3), p_val = format.pval(p_val, digits = 4), q_val = format.pval(q_val, digits = 4)) %>%
  select(Genus = taxon, `W Statistic` = W, `p-value` = p_val, `q-value (FDR)` = q_val) %>%
  arrange(`q-value (FDR)`)

write.csv(tabla_ancom_significativa, paste0(dir_abundances, "Table_ANCOM_Global_Significant.csv"), row.names = FALSE)

matriz_ancom_log <- as.data.frame(output$bias_correct_log_table)

# Imputación avanzada agregando ruido basal normal/uniforme para evitar enlaces artificiales en redes
matriz_ancom_log_imputada <- as.data.frame(lapply(matriz_ancom_log, function(x) {
  if (any(is.na(x))) {
    min_val <- min(x, na.rm = TRUE)
    if (is.infinite(min_val)) min_val <- -10 
    n_nas <- sum(is.na(x))
    if (n_nas > 0) {
      ruido_basal <- rnorm(n_nas, mean = min_val, sd = 0.01)
      x[is.na(x)] <- ruido_basal
    }
  }
  return(x)
}))
rownames(matriz_ancom_log_imputada) <- rownames(matriz_ancom_log)

df_bacterias <- as.data.frame(t(matriz_ancom_log_imputada))
bacterias_ancom <- df_global$taxon[df_global$diff_abn == TRUE]
df_bacterias_filtrado <- df_bacterias[, bacterias_ancom, drop = FALSE]
colnames(df_bacterias_filtrado) <- gsub("-", "_", colnames(df_bacterias_filtrado))

df_bacterias_filtrado$Sample_ID <- rownames(df_bacterias_filtrado)

df_final_modelo <- metadata_base %>%
  dplyr::select(-any_of("Sample_ID")) %>%
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::left_join(df_bacterias_filtrado, by = "Sample_ID")

# ==============================================================================
# SECTION 6: RISK MODELS - GLOBAL BRONCHIOLITIS SUBSET (GUARDADO DIRECTO EN MODELS)
# ==============================================================================
print("--- RUNNING RISK MODELS FOR TOTAL BRONCHIOLITIS ---")

df_model_bql <- df_final_modelo[df_final_modelo$Bronchiolitis == "Yes", ]
df_model_bql$Wheezing_Bin <- ifelse(df_model_bql$Wheezing.treatment == "Yes", 1, 0)
df_model_bql$Age <- as.numeric(df_model_bql$Age)

bacterias_presentes <- setdiff(intersect(colnames(df_bacterias_filtrado), colnames(df_model_bql)), "Sample_ID")

# 6.1 Univariate Model
lista_resultados_global <- list()
for (bac in bacterias_presentes) {
  if(length(unique(df_model_bql[[bac]])) < 2) next
  if(length(unique(df_model_bql$Wheezing_Bin)) < 2) next
  
  formula_u <- as.formula(paste0("Wheezing_Bin ~ Age + `", bac, "`"))
  modelo_u  <- tryCatch(glm(formula_u, data = df_model_bql, family = binomial), error = function(e) NULL)
  
  if(!is.null(modelo_u)) {
    res <- tryCatch(tidy(modelo_u, exponentiate = TRUE, conf.int = TRUE), error = function(e) tidy(modelo_u, exponentiate = TRUE, conf.int = FALSE))
    lista_resultados_global[[bac]] <- res %>% filter(term == bac) %>% mutate(Bacteria = bac)
  }
}
tabla_univariada_global <- bind_rows(lista_resultados_global)

if(nrow(tabla_univariada_global) > 0) {
  tabla_univariada_global <- tabla_univariada_global %>%
    select(Bacteria, `Odds Ratio (OR)` = estimate, `Lower CI (2.5%)` = conf.low, `Upper CI (97.5%)` = conf.high, `p-value` = p.value)
  write.csv(tabla_univariada_global, paste0(dir_models, "Table_Univariate_Risk.csv"), row.names = FALSE)
  
  bacterias_sig_uni_global <- tabla_univariada_global %>% filter(`p-value` < 0.05) %>% pull(Bacteria)
  
  # Plot Univariate
  tabla_univariada_global$Bacteria <- factor(tabla_univariada_global$Bacteria, levels = rev(tabla_univariada_global$Bacteria))
  plot_forest_uni_g <- ggplot(tabla_univariada_global, aes(x = `Odds Ratio (OR)`, y = Bacteria)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 0.8) +
    geom_errorbarh(aes(xmin = `Lower CI (2.5%)`, xmax = `Upper CI (97.5%)`), height = 0.2, color = "#34495e", size = 0.8) +
    geom_point(size = 3.5, color = "#2980b9") + theme_minimal() +
    labs(title = "Forest Plot: Univariate Predictors (Total Bronchiolitis)", subtitle = "Adjusted by Age", x = "Odds Ratio", y = "") +
    scale_x_log10() + theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle = element_text(hjust=0.5))
  ggsave(paste0(dir_models, "Plot_ForestPlot_Univariate.png"), plot = plot_forest_uni_g, width = 8, height = 6, dpi = 300)
}

# 6.2 Multivariate Model & ROC Validation
if(exists("bacterias_sig_uni_global") && length(bacterias_sig_uni_global) > 0) {
  formula_m_g <- as.formula(paste("Wheezing_Bin ~ Age +", paste(paste0("`", bacterias_sig_uni_global, "`"), collapse = " + ")))
  modelo_m_g   <- glm(formula_m_g, data = df_model_bql, family = binomial, na.action = na.exclude)
  
  tabla_m_g <- tidy(modelo_m_g, exponentiate = TRUE, conf.int = TRUE) %>% filter(term != "(Intercept)") %>%
    select(term, `Odds Ratio (OR)` = estimate, `Lower CI (2.5%)` = conf.low, `Upper CI (97.5%)` = conf.high, `p-value` = p.value)
  write.csv(tabla_m_g, paste0(dir_models, "Table_Multivariate_Model_Definitive.csv"), row.names = FALSE)
  
  df_model_bql$probabilidades <- predict(modelo_m_g, type = "response")
  df_roc_g <- df_model_bql %>% filter(!is.na(probabilidades) & !is.na(Wheezing_Bin))
  objeto_roc_g <- roc(df_roc_g$Wheezing_Bin, df_roc_g$probabilidades, levels = c(0, 1), direction = "<", quiet = TRUE)
  auc_ci_g <- ci.auc(objeto_roc_g)
  
  png(filename = paste0(dir_models, "Plot_ROC_Curve_Multivariate.png"), width = 2000, height = 1800, res = 300)
  plot(objeto_roc_g, col = "#1f77b4", lwd = 3, main = paste0("ROC Curve (Total Bronchiolitis)\nAUC = ", round(auc_ci_g[2], 3)), identity.col = "grey", identity.lty = 2)
  grid()
  dev.off()
}

# ==============================================================================
# SECTION 7: ADVANCED VISUALIZATIONS (AGE-ADJUSTED MEGABOXPLOTS)
# ==============================================================================
print("--- RENDERING CLEAN AGE-ADJUSTED MEGABOXPLOTS ---")

df_bacterias_cambio <- df_bacterias
for(b in bacterias_ancom) {
  col_original <- colnames(df_bacterias_cambio)[grep(paste0("^", gsub("-", ".", b)), colnames(df_bacterias_cambio))]
  if(length(col_original) > 0) colnames(df_bacterias_cambio)[colnames(df_bacterias_cambio) == col_original] <- b
}
df_bacterias_cambio$Sample_ID <- rownames(df_bacterias_cambio)

metadata_base$Clinical_Group <- metadata_base$Condicion_Clinica

df_base_ajuste <- metadata_base %>%
  dplyr::select(-any_of("Sample_ID")) %>% tibble::rownames_to_column("Sample_ID") %>%
  dplyr::left_join(df_bacterias_cambio, by = "Sample_ID") %>% mutate(Age = as.numeric(Age))

bacterias_procesar <- gsub("Burkholderia-Caballeronia-Paraburkholderia", "Burkholderia", bacterias_ancom)
colnames(df_base_ajuste) <- gsub("Burkholderia-Caballeronia-Paraburkholderia", "Burkholderia", colnames(df_base_ajuste))

bacterias_procesar <- make.names(bacterias_procesar)
colnames(df_base_ajuste) <- make.names(colnames(df_base_ajuste))

for (bac in bacterias_procesar) {
  formula_res <- as.formula(paste0("`", bac, "` ~ Age"))
  fit_res     <- lm(formula_res, data = df_base_ajuste, na.action = na.exclude)
  df_base_ajuste[[paste0(bac, "_adj")]] <- residuals(fit_res) + mean(df_base_ajuste[[bac]], na.rm = TRUE)
}

columnas_adj <- paste0(bacterias_procesar, "_adj")
df_long_adj <- df_base_ajuste %>%
  select(Sample_ID, Clinical_Group, Age, all_of(columnas_adj)) %>%
  pivot_longer(cols = all_of(columnas_adj), names_to = "Bacteria", values_to = "Abundancia_Adjusted") %>%
  mutate(Bacteria = gsub("_adj$", "", Bacteria), Bacteria = gsub("\\.", "-", Bacteria))

plot_boxplots_clean <- ggplot(df_long_adj, aes(x = Clinical_Group, y = Abundancia_Adjusted, fill = Clinical_Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) + geom_jitter(width = 0.15, alpha = 0.3, color = "black", size = 0.9) +
  facet_wrap(~ Bacteria, scales = "free_y", ncol = 3) + theme_bw() + scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "Age-Adjusted Abundance (Residuals + Mean)", title = "Bacterial Abundance Dynamics by Clinical Group", subtitle = "Linear Adjustment by Age") + 
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 14, hjust = 0.5), plot.subtitle = element_text(size = 11, hjust = 0.5, face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#f8f9fa", color = "gray80"), strip.text = element_text(color = "black", face = "bold", size = 10))

ggsave(filename = paste0(dir_abundances, "Plot_MegaBoxplot_Minimal.png"), plot = plot_boxplots_clean, width = 12, height = 9, dpi = 300)
