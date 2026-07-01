# ==============================================================================
# SECTION 1: LIBRARIES, CONFIGURATION & ENVIRONMENT (CROSS-NICHE MULTIVARIATE)
# ==============================================================================
suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pROC)
  library(broom)
  library(patchwork)
})

set.seed(123)

setwd("/mnt/usb/BQL/BQL_ANALYSIS")
base_dir       <- "COMPLETE_REANALISIS_BQL/REANALISIS/CROSS_NICHE/"
dir_abundances <- paste0(base_dir, "abundances/")
dir_models     <- paste0(base_dir, "models/")

for (d in c(dir_abundances, dir_models)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ==============================================================================
# FUNCIONES AUXILIARES NATIVAS
# ==============================================================================
extract_meta <- function(pseq_obj) {
  df <- data.frame(as(sample_data(pseq_obj), "data.frame"),
                   stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Sample_ID" %in% colnames(df)) df$Sample_ID <- rownames(df)
  rownames(df) <- df$Sample_ID
  return(df)
}

clean_taxa_names <- function(x) {
  x <- gsub("[^[:alnum:]_]", "_", x)
  x <- gsub("_{2,}", "_", x)
  x <- gsub("^_|_$", "", x)
  return(make.unique(x, sep = "_dup"))
}

rename_to_genus <- function(pseq_obj, genus_col) {
  tt <- as.data.frame(tax_table(pseq_obj), stringsAsFactors = FALSE)
  gen_labels <- tt[[genus_col]]
  higher_ranks <- rev(colnames(tt)[colnames(tt) != genus_col])
  for (i in seq_along(gen_labels)) {
    if (is.na(gen_labels[i]) || gen_labels[i] == "" ||
        grepl("^NA$|^uncultured$|^unknown$", gen_labels[i], ignore.case = TRUE)) {
      for (rk in higher_ranks) {
        val <- tt[i, rk]
        if (!is.na(val) && val != "" && !grepl("^NA$|^uncultured$|^unknown$", val, ignore.case = TRUE)) {
          gen_labels[i] <- paste0("Unclassified_", val)
          break
        }
      }
      if (is.na(gen_labels[i])) gen_labels[i] <- "Unclassified_Bacteria"
    }
  }
  taxa_names(pseq_obj) <- clean_taxa_names(gen_labels)
  rownames(tax_table(pseq_obj)) <- taxa_names(pseq_obj)
  return(pseq_obj)
}

# ==============================================================================
# 2. CARGA, AGREGACIÓN A GÉNERO Y FILTRADO PAREADO ESTRICTO
# ==============================================================================
cat("\n[2] Procesando y emparejando matriz trans-nicho...\n")
data_input_dir <- "COMPLETE_REANALISIS_BQL/REANALISIS/Global/curated_data_global/"
pseq_raw_master <- readRDS(paste0(data_input_dir, "phyloseq_RAW_curated_global.rds"))

genus_col <- colnames(tax_table(pseq_raw_master))[grepl("^[Gg]enus$", colnames(tax_table(pseq_raw_master)))][1]

pseq_genus <- tax_glom(pseq_raw_master, taxrank = genus_col, NArm = FALSE) %>% 
  rename_to_genus(genus_col)

pseq_prop <- transform_sample_counts(pseq_genus, function(x) log10((x / sum(x)) * 100 + 1))

metadata_df <- extract_meta(pseq_prop) %>%
  mutate(
    Age = as.numeric(scale(as.numeric(Age))),
    Patient_ID  = as.factor(Patient_ID),
    Sample_Type = ifelse(Sample_Type %in% c("ANF", "NPA"), "NPA", "GUT")
  )

metadata_df$Bronchiolitis <- ifelse(metadata_df$Condicion_Clinica == "CTRL", "No", "Yes")
metadata_df$Wheezing_Bin <- ifelse(
  grepl("\\+$", metadata_df$Condicion_Clinica) | metadata_df$Condicion_Clinica %in% c("RSV-/Wheeze+", "RSV+/Wheeze+"), 
  1, 0
)

df_model_bql_pre <- metadata_df %>% filter(Bronchiolitis == "Yes")

patients_with_both <- df_model_bql_pre %>% 
  group_by(Patient_ID) %>%
  filter(all(c("NPA", "GUT") %in% Sample_Type)) %>% 
  pull(Patient_ID) %>% 
  unique()

otu_mat <- as.data.frame(t(as(otu_table(pseq_prop), "matrix")))

target_npa_clean <- rownames(otu_mat)[grep("Burkholderia", rownames(otu_mat), ignore.case = TRUE)[1]]
target_gut_clean <- rownames(otu_mat)[grep("Ruminococcus_gnavus", rownames(otu_mat), ignore.case = TRUE)[1]]

vec_burkholderia <- as.numeric(otu_mat[target_npa_clean, ])
names(vec_burkholderia) <- colnames(otu_mat)

vec_ruminococcus <- as.numeric(otu_mat[target_gut_clean, ])
names(vec_ruminococcus) <- colnames(otu_mat)

df_model_bql_pre$Burkholderia_NPA <- vec_burkholderia[df_model_bql_pre$Sample_ID]
df_model_bql_pre$Ruminococcus_GUT <- vec_ruminococcus[df_model_bql_pre$Sample_ID]

df_npa_side <- df_model_bql_pre %>% 
  filter(Sample_Type == "NPA" & Patient_ID %in% patients_with_both) %>%
  select(Patient_ID, Age, Wheezing_Bin, Burkholderia_NPA, 
         Respiratory.syncytial.virus, Family.history.atopy, Passive.smoking, Cesarean.section, Previous.antibiotics)

df_gut_side <- df_model_bql_pre %>% 
  filter(Sample_Type == "GUT" & Patient_ID %in% patients_with_both) %>%
  select(Patient_ID, Ruminococcus_GUT)

df_final_modelo <- df_npa_side %>%
  dplyr::inner_join(df_gut_side, by = "Patient_ID") %>%
  drop_na(Burkholderia_NPA, Ruminococcus_GUT, Age, Wheezing_Bin,
          Respiratory.syncytial.virus, Family.history.atopy, Passive.smoking, Cesarean.section, Previous.antibiotics)

# ------------------------------------------------------------------------------
# ARMONIZACIÓN Y CONVERSIÓN FACTOR-DUMMY DE COVARIABLES CLÍNICAS (RSV INCLUIDO)
# ------------------------------------------------------------------------------
df_final_modelo <- df_final_modelo %>%
  mutate(
    RSV                  = ifelse(Respiratory.syncytial.virus == "Yes", 1, 0),
    Family.history.atopy = ifelse(Family.history.atopy == "Yes", 1, 0),
    Passive.smoking      = ifelse(Passive.smoking == "Yes", 1, 0),
    Cesarean.section     = ifelse(Cesarean.section == "Yes", 1, 0),
    Previous.antibiotics = ifelse(Previous.antibiotics == "Yes", 1, 0)
  )

cat(paste(" >> [SUCCESS] Cohorte BQL+ finalizada con", nrow(df_final_modelo), "pacientes perfectamente armonizados numéricamente.\n"))

# ==============================================================================
# 3. MODELO MULTIVARIANTE AJUSTADO CLÍNICO-MICROBIOMA (GLM)
# ==============================================================================
cat("\n[3] Ajustando modelo logístico multivariante expandido...\n")

mv_model <- glm(Wheezing_Bin ~ Burkholderia_NPA + Ruminococcus_GUT + Age + 
                  RSV + Family.history.atopy + Passive.smoking + 
                  Cesarean.section + Previous.antibiotics, 
                data = df_final_modelo, family = binomial)

# Procesamiento estricto e inequívoco de coeficientes
model_hits <- tidy(mv_model, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = case_when(
      term == "Burkholderia_NPA"        ~ "Burkholderia-\nCaballeronia-\nParaburkholderia\n(NPA)",
      term == "Ruminococcus_GUT"        ~ "[Ruminococcus]\ngnavus group\n(GUT)",
      term == "Age"                     ~ "Age (Months)",
      term == "RSV"                     ~ "RSV\n(Yes)",
      term == "Family.history.atopy"    ~ "Family History\nof Atopy (Yes)",
      term == "Passive.smoking"         ~ "Passive Smoking\n(Yes)",
      term == "Cesarean.section"        ~ "Cesarean Section\n(Yes)",
      term == "Previous.antibiotics"    ~ "Previous Antibiotics\n(Yes)",
      TRUE ~ term
    )
  )

# Guardar coeficientes numéricos en CSV de forma segura
write.csv(model_hits, file = paste0(dir_models, "Table_Wheezing_Predictors_Coefficients.csv"), row.names = FALSE)
cat(" >> [SUCCESS] Tabla de coeficientes exportada a CSV con éxito.\n")

# Forzar orden de factores para la visualización vertical
model_hits$term <- factor(model_hits$term, levels = rev(model_hits$term))

# ==============================================================================
# 4. PREPARACIÓN DINÁMICA DE LA CURVA ROC CON INTERVALOS DE CONFIANZA
# ==============================================================================
df_final_modelo$pred_prob <- predict(mv_model, type = "response")
roc_obj <- roc(df_final_modelo$Wheezing_Bin, df_final_modelo$pred_prob, quiet = TRUE)

roc_ci <- ci.se(roc_obj, specificities = seq(0, 1, l=50))
df_roc_ci <- data.frame(
  Specificity = as.numeric(rownames(roc_ci)),
  Lower = roc_ci[, 1],
  Upper = roc_ci[, 3]
)

df_roc_curve <- data.frame(
  Sensitivity = roc_obj$sensitivities,
  Specificity = roc_obj$specificities
)

real_ci <- ci.auc(roc_obj)
ci_low  <- round(real_ci[1], 3)
ci_high <- round(real_ci[3], 3)
auc_val <- round(auc(roc_obj), 3)

# ==============================================================================
# 5. GENERACIÓN AUTOMÁTICA DEL PANEL COMBINADO DE GRÁFICOS (PATCHWORK RE-ESCALADO)
# ==============================================================================
# PANEL a: Forest Plot Limpio y Ajustado en escala simétrica Log10
plot_7a <- ggplot(model_hits, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#e74c3c", size = 0.7) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15, color = "#2c3e50", size = 0.8) +
  geom_point(color = "#3498db", size = 3.5) +
  scale_x_log10(
    breaks = c(0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 50, 500),
    labels = c("0.001", "0.01", "0.1", "0.5", "1", "2", "5", "10", "50", "500")
  ) +
  theme_bw(base_size = 11) +
  labs(x = "Odds Ratio (95% Confidence Interval)", y = "") +
  theme(
    axis.text.y = element_text(face = "bold", color = "black", lineheight = 0.85),
    axis.text.x = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray95")
  )

# PANEL b: Curva ROC Dinámica sin subtítulos internos
plot_7b <- ggplot() +
  geom_ribbon(data = df_roc_ci, aes(x = 1 - Specificity, ymin = Lower, ymax = Upper), 
              fill = "#2ecc71", alpha = 0.15) +
  geom_path(data = df_roc_curve, aes(x = 1 - Specificity, y = Sensitivity), 
            color = "#2ecc71", size = 1.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60") +
  theme_bw(base_size = 11) +
  annotate("text", x = 0.95, y = 0.15, 
           label = paste0("AUC = ", auc_val, "\n95% CI: [", ci_low, " - ", ci_high, "]"), 
           size = 4.2, fontface = "bold", color = "#2c3e50", hjust = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme(panel.grid.minor = element_blank())

# ENSAMBLADO FINAL CON PATCHWORK EQUILIBRADO
supp_figure_7 <- (plot_7a / plot_7b) + 
  plot_layout(heights = c(1.8, 1)) +
  plot_annotation(tag_levels = list(c("a", "b"))) & 
  theme(plot.tag = element_text(face = "bold", size = 13))

# Guardado directo sobreescribiendo el panel defectuoso
ggsave(filename = paste0(dir_models, "Supplementary_Figure_7_Integrated_Panel.png"), 
       plot = supp_figure_7, width = 7, height = 10, dpi = 300)

cat("\n >> [PIPELINE FINISHED] Modelo corregido, CSV guardado y gráfica Premium generada.\n")
