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
library(cluster)
library(pheatmap)
library(stringr)
library(rstatix)

set.seed(123)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4) 

setwd("/mnt/usb/BQL/BQL_ANALYSIS/COMPLETE_REANALISIS_BQL/")

# Configuración estricta de la jerarquía de carpetas enfocada en el nicho intestinal (GUT)
base_dir        <- "REANALISIS/GUT/"
dir_alpha       <- paste0(base_dir, "alpha/")
dir_beta        <- paste0(base_dir, "beta/")
dir_abundances  <- paste0(base_dir, "abundances/")
dir_models      <- paste0(base_dir, "models/") # Guardado directo sin subcarpetas

# Creación recursiva de todos los directorios del proyecto para GUT
for(path in c(dir_alpha, dir_beta, dir_abundances, dir_models)){
  if(!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# ==============================================================================
# SECTION 2: GLOBAL DATA PRE-PROCESSING & SUBSETTING (GUT FOCUS)
# ==============================================================================
data_input_dir <- "REANALISIS/Global/curated_data_global/"

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

# Filtrado exclusivo del nicho intestinal (GUT / Gut)
pseq_gut      <- subset_samples(pseq_raw_master, Sample_Type == "GUT" | Sample_Type == "Gut")
pseq_gut_rare <- subset_samples(pseq_rare_master, Sample_Type == "GUT" | Sample_Type == "Gut")

metadata_base <- as(sample_data(pseq_gut), "data.frame")
tse_gut       <- convertFromPhyloseq(pseq_gut)

dist_matrix_final <- phyloseq::distance(pseq_gut_rare, method = "bray")
muestras_comunes  <- labels(dist_matrix_final)
metadata_beta     <- metadata_base[muestras_comunes, ]

# ==============================================================================
# SECTION 3: ALPHA DIVERSITY (UNIVARIATE ANALYSIS - AGE ADJUSTED)
# ==============================================================================
alpha_indices <- estimate_richness(pseq_gut_rare, measures = c("Observed", "Shannon", "Simpson"))
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

# ==============================================================================
# SECTION 3 & 4: DIVERSITY PLOTTING (INTEGRATED PREMIUM PANEL)
# ==============================================================================
cat("\n[3 & 4] Generando Panel de Figuras Unificado para GUT (Modo Paper - Proporciones Corregidas)...\n")
library(patchwork)

# --- Diseño: ALFA DIVERSIDAD (Panel Superior 2a) ---
plot_alpha <- ggplot(df_alpha_long, aes(x = Condicion_Clinica, y = Value, fill = Condicion_Clinica)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.55, color = "#2c3e50", lwd = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.25, color = "black", size = 0.8) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 11) + 
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "Alpha Diversity Value\n(Age-Corrected)") +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        strip.background = element_rect(fill = "#2c3e50", color = NA), 
        strip.text = element_text(color = "white", face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#f0f0f0")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))

if(nrow(stat_results) > 0) {
  stat_results <- stat_results %>% add_y_position(step.increase = 0.04)
  plot_alpha <- plot_alpha + stat_pvalue_manual(stat_results, label = "p.adj.signif", hide.ns = TRUE, step.increase = 0.04, tip.length = 0.01)
}

# --- Diseño: BETA DIVERSIDAD (Panel Inferior 2b) ---
plot_pcoa <- ggplot(df_pcoa, aes(x = PCoA1, y = PCoA2, color = Clinical_Group, fill = Clinical_Group)) +
  geom_point(size = 2.2, alpha = 0.65) + 
  stat_ellipse(geom = "polygon", alpha = 0.08, level = 0.95, size = 0.4, linetype = 2) +
  theme_bw(base_size = 11) + 
  scale_color_brewer(palette = "Set2", name = "Clinical Phenotype Group") + 
  scale_fill_brewer(palette = "Set2", name = "Clinical Phenotype Group") +
  labs(x = paste0("PCoA 1 (", porcentaje_variacion[1], "%)"), 
       y = paste0("PCoA 2 (", porcentaje_variacion[2], "%)")) +
  theme(axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#f0f0f0"),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        legend.position = "right")

# --- Ensamblado y guardado final (Ajuste de ratio 0.7 a 1 para balancear alturas) ---
integrated_panel <- plot_alpha / plot_pcoa + 
  plot_layout(ncol = 1, heights = c(0.7, 1)) + # CORRECCIÓN: Menos alto arriba, más espacio abajo
  plot_annotation(
   tag_levels = 'a'
  ) & 
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "#2c3e50"),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic", color = "#555555"),
    plot.tag = element_text(face = "bold", size = 14, color = "#2c3e50")
  )

ggsave(
  filename = paste0(base_dir, "Figure2_GUT_Diversity_Integrated_Panel.png"), 
  plot = integrated_panel, 
  width = 10, 
  height = 9.0, # Ajustado levemente el alto global para compactar el lienzo de dibujo
  dpi = 300
)

cat(" >> ¡Análisis de diversidad de GUT completado con éxito con proporciones balanceadas!\n")

# ==============================================================================
# SECTION 5: COMPOSITIONAL ANALYSIS (ANCOM-BC2) WITH BIAS-CORRECTED JITTER IMPUTATION
# ==============================================================================
output = ancombc2(
  data = tse_gut, assay_name = "counts", tax_level = "Genus",         
  fix_formula = "Age + Condicion_Clinica", p_adj_method = "fdr", group = "Condicion_Clinica",        
  dunnet = FALSE, pairwise = TRUE, global = TRUE, alpha = 0.05
)

# 1. Extraer las tablas nativas del objeto de salida
df_global <- output$res_global
df_pair   <- output$res_pair  # Contiene los LFC y p_val para cada par de grupos

# 2. Filtrado Global Robusto (Uso estricto de diff_robust_abn para evitar falsos positivos)
taxones_robust_globales <- df_global %>% 
  filter(diff_robust_abn == TRUE) %>% 
  pull(taxon)

# 3. Mapeo de diferencias entre grupos (Pairwise Analysis)
library(dplyr)
library(tidyr)
library(stringr)

tabla_ancom_grupos <- df_pair %>%
  # 1. Filtramos por tus taxones robustos globales
  filter(taxon %in% taxones_robust_globales) %>%
  
  # 2. Pivotamos usando un patrón adaptado exactamente a tus prefijos de ANCOM-BC2.
  # Captura los prefijos conocidos (lfc, se, W, p, q, diff, passed_ss, diff_robust) 
  # y deja todo lo que viene detrás en la columna "Comparacion"
  pivot_longer(
    cols = -taxon,
    names_to = c(".value", "Comparacion"),
    names_pattern = "^(lfc|se|W|p|q|diff|passed_ss)_(.*)"
  ) %>%
  
  # 3. Limpiamos por completo el molesto prefijo "Condicion_Clinica" de los nombres de los grupos
  mutate(Comparacion = str_remove_all(Comparacion, "Condicion_Clinica")) %>%
  
  # 4. Seleccionamos y renombramos las columnas para tu reporte final
  # Incluyendo p_value, q_value (FDR) y el control de calidad passed_ss
  select(
    Genus = taxon, 
    Comparacion_Especifica = Comparacion, 
    `Log Fold Change (LFC)` = lfc, 
    `W Statistic` = W,
    `Standard Error (SE)` = se,
    `p-value` = p,
    `q-value (FDR)` = q,
    `Passed Sensitivity (SS)` = passed_ss
  ) %>%
  
  # 5. Opcional: Si quieres quitar las filas de control basal de intercepto (ej. las que solo dicen "RSV-/Wheeze-")
  # y quedarte solo con las comparaciones reales cara a cara (las que contienen el guion bajo "_"), desenta esta linea:
  # filter(str_detect(Comparacion_Especifica, "_")) %>%
  
  # 6. Ordenamos por Genus y por significancia estadística del p-valor bruto
  # ... (todo tu select anterior igual) ...
  # JUSTO AQUÍ AÑADES ESTA LÍNEA:
  filter(!is.na(`Log Fold Change (LFC)`)) %>%
  arrange(Genus, `p-value`)

# 4. Almacenamiento de los resultados para Modelos y Plots
# Guardamos la tabla detallada que te dice exactamente qué par de grupos difiere y en qué dirección (LFC)
write.csv(tabla_ancom_grupos, paste0(dir_abundances, "Table_ANCOM_Pairwise_From_Global_Robust.csv"), row.names = FALSE)

# Guardamos también un resumen limpio de los ganadores globales para los plots
df_global_significativa <- df_global %>%
  filter(diff_robust_abn == TRUE)
write.csv(df_global_significativa, paste0(dir_abundances, "Table_ANCOM_Global_Robust_Signif.csv"), row.names = FALSE)

cat(paste0(" >> ANCOM-BC2 completado. Se detectaron ", length(taxones_robust_globales), " géneros con variación global robusta.\n"))

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
bacterias_ancom <- df_global$taxon[df_global$diff_robust_abn == TRUE]
df_bacterias_filtrado <- df_bacterias[, bacterias_ancom, drop = FALSE]
colnames(df_bacterias_filtrado) <- gsub("-", "_", colnames(df_bacterias_filtrado))

df_bacterias_filtrado$Sample_ID <- rownames(df_bacterias_filtrado)

df_final_modelo <- metadata_base %>%
  dplyr::select(-any_of("Sample_ID")) %>%
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::left_join(df_bacterias_filtrado, by = "Sample_ID")

# ==============================================================================
# SECTION 6: RISK MODELS - GLOBAL BRONCHIOLITIS SUBSET & MULTIVARIATE FOREST PLOT
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
    labs(title = "Forest Plot: Univariate Predictors (Total Bronchiolitis)", subtitle = "Adjusted by Age", x = "Odds Ratio (Log Scale)", y = "") +
    scale_x_log10() + theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle = element_text(hjust=0.5))
  ggsave(paste0(dir_models, "Plot_ForestPlot_Univariate.png"), plot = plot_forest_uni_g, width = 8, height = 6, dpi = 300)
}

# 6.2 Multivariate Model, ROC Validation & Forest Plot
if(exists("bacterias_sig_uni_global") && length(bacterias_sig_uni_global) > 0) {
  formula_m_g <- as.formula(paste("Wheezing_Bin ~ Age +", paste(paste0("`", bacterias_sig_uni_global, "`"), collapse = " + ")))
  modelo_m_g   <- glm(formula_m_g, data = df_model_bql, family = binomial, na.action = na.exclude)
  
  tabla_m_g <- tidy(modelo_m_g, exponentiate = TRUE, conf.int = TRUE) %>% 
    filter(term != "(Intercept)") %>%
    mutate(Variable = sapply(term, function(x) {
      x_clean <- gsub("`", "", x)
      if(x_clean == "Age") return("Age")
      idx <- grep(paste0("^", gsub("-", ".", x_clean)), gsub("-", ".", bacterias_ancom))
      if(length(idx) > 0) return(bacterias_ancom[idx[1]]) else return(x_clean)
    })) %>%
    select(Variable, `Odds Ratio (OR)` = estimate, `Lower CI (2.5%)` = conf.low, `Upper CI (97.5%)` = conf.high, `p-value` = p.value)
  
  write.csv(tabla_m_g, paste0(dir_models, "Table_Multivariate_Model_Definitive.csv"), row.names = FALSE)
  
  # Validación de Curva ROC
  df_model_bql$probabilidades <- predict(modelo_m_g, type = "response")
  df_roc_g <- df_model_bql %>% filter(!is.na(probabilidades) & !is.na(Wheezing_Bin))
  objeto_roc_g <- roc(df_roc_g$Wheezing_Bin, df_roc_g$probabilidades, levels = c(0, 1), direction = "<", quiet = TRUE)
  auc_ci_g <- ci.auc(objeto_roc_g)
  
  png(filename = paste0(dir_models, "Plot_ROC_Curve_Multivariate.png"), width = 2000, height = 1800, res = 300)
  plot(objeto_roc_g, col = "#1f77b4", lwd = 3, main = paste0("ROC Curve (Total Bronchiolitis)\nAUC = ", round(auc_ci_g[2], 3)), identity.col = "grey", identity.lty = 2)
  grid()
  dev.off()
  
  # Generación del Forest Plot Multivariante 
  variables_ordenadas  <- unique(tabla_m_g$Variable)
  tabla_m_g$Variable   <- factor(tabla_m_g$Variable, levels = rev(variables_ordenadas))
  
  plot_forest_multi <- ggplot(tabla_m_g, aes(x = `Odds Ratio (OR)`, y = Variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 0.8) +
    geom_errorbarh(aes(xmin = `Lower CI (2.5%)`, xmax = `Upper CI (97.5%)`), height = 0.2, color = "#2c3e50", size = 0.8) +
    geom_point(size = 4, color = "#16a085") + 
    theme_minimal() +
    labs(title    = "Forest Plot: Multivariate Predictors of Recurrent Wheezing",
         subtitle = paste0("Model Adjusted by Age (AUC = ", round(auc_ci_g[2], 3), ")"), 
         x        = "Odds Ratio (Log Scale)", y        = "") +
    scale_x_log10() + 
    theme(text = element_text(size = 13), 
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5), 
          panel.grid.minor = element_blank())
  
  ggsave(paste0(dir_models, "Plot_ForestPlot_Multivariate.png"), plot = plot_forest_multi, width = 8, height = 5, dpi = 300)
  print("--- TABLA Y FOREST PLOT MULTIVARIANTE GENERADOS CON ÉXITO ---")
}

# ==============================================================================
# SECTION 7: ADVANCED VISUALIZATIONS (2x5 GRID & FIXED ANNOTATIONS - GUT)
# ==============================================================================
print("--- RENDERING CLEAN AGE-ADJUSTED MEGABOXPLOTS (2x5 GRID) ---")

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(stringr)

# 1. Copiamos y homogeneizamos nombres en la matriz de bacterias original
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

# --- TRATAMIENTO SEGURO DE CARACTERES ---
bacterias_mapeo <- tibble(
  Original = bacterias_ancom,
  Seguro   = make.names(bacterias_ancom)
)

for(i in 1:nrow(bacterias_mapeo)){
  colnames(df_base_ajuste)[colnames(df_base_ajuste) == bacterias_mapeo$Original[i]] <- bacterias_mapeo$Seguro[i]
}

for (bac_segura in bacterias_mapeo$Seguro) {
  formula_res <- as.formula(paste0("`", bac_segura, "` ~ Age"))
  fit_res     <- lm(formula_res, data = df_base_ajuste, na.action = na.exclude)
  df_base_ajuste[[paste0(bac_segura, "_adj")]] <- residuals(fit_res) + mean(df_base_ajuste[[bac_segura]], na.rm = TRUE)
}

columnas_adj <- paste0(bacterias_mapeo$Seguro, "_adj")
df_long_adj <- df_base_ajuste %>%
  dplyr::select(Sample_ID, Clinical_Group, Age, all_of(columnas_adj)) %>%
  tidyr::pivot_longer(cols = all_of(columnas_adj), names_to = "Bacteria_Segura", values_to = "Abundancia_Adjusted") %>%
  dplyr::mutate(Bacteria_Segura = gsub("_adj$", "", Bacteria_Segura)) %>%
  dplyr::left_join(bacterias_mapeo, by = c("Bacteria_Segura" = "Seguro")) %>%
  dplyr::select(Sample_ID, Clinical_Group, Age, Bacteria = Original, Abundancia_Adjusted)

# ==============================================================================
# 1. MATRIZ DE ASTERISCOS - TRADUCCIÓN EXACTA BASADA EN EL INTERCEPTO (CTRL)
# ==============================================================================
# Extraemos el valor máximo real de cada faceta para que sirva de base de altura
df_max_valores <- df_long_adj %>%
  dplyr::group_by(Bacteria) %>%
  dplyr::summarise(max_val = max(Abundancia_Adjusted, na.rm = TRUE), .groups = "drop")

df_p_annotations <- tabla_ancom_grupos %>%
  dplyr::filter(`p-value` < 0.05) %>%
  dplyr::mutate(Bacteria = Genus) %>%  
  dplyr::mutate(
    # LÓGICA ESTRICTA: 
    # Si NO hay guion bajo, es una comparación directa del Intercepto (CTRL) vs ese grupo.
    # Si SÍ hay guion bajo, separamos de forma limpia el grupo 1 del grupo 2.
    group1 = case_when(
      str_detect(Comparacion_Especifica, "_") ~ str_split_i(Comparacion_Especifica, "_", 1),
      TRUE                                    ~ "CTRL"
    ),
    group2 = case_when(
      str_detect(Comparacion_Especifica, "_") ~ str_split_i(Comparacion_Especifica, "_", 2),
      TRUE                                    ~ Comparacion_Especifica
    )
  ) %>%
  # Filtro de seguridad por si acaso quedara alguna autocomparación huérfana
  dplyr::filter(group1 != group2) %>% 
  # Mapeo clásico de significancia por asteriscos
  dplyr::mutate(label = case_when(
    `p-value` < 0.001 ~ "***",
    `p-value` < 0.01  ~ "**",
    TRUE              ~ "*"
  )) %>%
  dplyr::select(Bacteria, group1, group2, label) %>%
  dplyr::inner_join(df_max_valores, by = "Bacteria") %>% 
  
  # ESCALONAMIENTO SEGURO: 
  # Si una bacteria (ej. Blautia) tiene varias comparaciones significativas, 
  # subimos cada bracket un "piso" (0.15 de margen incremental) para que no se pisen entre ellos.
  dplyr::group_by(Bacteria) %>%
  dplyr::mutate(idx = row_number()) %>%
  dplyr::mutate(y.position = max_val + (abs(max_val) * (0.05 + (idx - 1) * 0.15))) %>%
  dplyr::ungroup()

# ==============================================================================
# 2. GENERACIÓN DINÁMICA DE LETRAS DE PANEL (A-J)
# ==============================================================================
bacterias_ordenadas <- sort(unique(df_long_adj$Bacteria))
df_letras_paneles <- tibble(
  Bacteria = bacterias_ordenadas,
  Letra = letters[1:length(bacterias_ordenadas)] 
)

# ==============================================================================
# 3. CONSTRUCCIÓN DEL PLOT (CUADRÍCULA OPTIMIZADA DE 2 COLUMNAS x 5 FILAS)
# ==============================================================================
plot_boxplots_clean <- ggplot(df_long_adj, aes(x = Clinical_Group, y = Abundancia_Adjusted, fill = Clinical_Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) + 
  geom_jitter(width = 0.15, alpha = 0.3, color = "black", size = 1.0) +
  # CAMBIO CLAVE: Cambiado de ncol = 3 a ncol = 2
  facet_wrap(~ Bacteria, scales = "free_y", ncol = 2) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set2") +
  
  labs(
    x = NULL, 
    y = "Bias-Corrected Age-Adjusted Abundance (ANCOM-BC2 Residuals)"
  ) + 
  
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 10)),
    axis.text.y  = element_text(size = 11, color = "black"),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95"),
    strip.background = element_rect(fill = "#f8f9fa", color = "gray80"), 
    strip.text = element_text(color = "black", face = "bold.italic", size = 13) 
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) 

# Inyección de letras de panel
plot_boxplots_clean <- plot_boxplots_clean +
  geom_text(
    data = df_letras_paneles,
    aes(label = Letra),
    x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5,
    inherit.aes = FALSE, fontface = "bold", size = 6, color = "black"
  )

# Inyección de brackets y asteriscos corregidos
if (nrow(df_p_annotations) > 0) {
  plot_boxplots_clean <- plot_boxplots_clean +
    stat_pvalue_manual(
      data = df_p_annotations,
      label = "label",
      y.position = "y.position",
      step.increase = 0.13, 
      hide.ns = TRUE,
      tip.length = 0.02,
      size = 5
    )
}

# Guardamos el TIFF modificando las dimensiones (más estrecho y más alto para el layout 2x5)
ggsave(filename = paste0(dir_abundances, "Figure_Microbiome_Dynamics.tiff"), 
       plot = plot_boxplots_clean, width = 11, height = 18, dpi = 300, device = "tiff")