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
})

set.seed(123)

# --- Configuración de Rutas de Trabajo (Adaptadas a tu jerarquía CROSS_NICHE) ---
setwd("/mnt/usb/BQL/BQL_ANALYSIS")
base_dir       <- "COMPLETE_REANALISIS_BQL/REANALISIS/CROSS_NICHE/"
dir_abundances <- paste0(base_dir, "abundances/")
dir_models     <- paste0(base_dir, "models/")

for (d in c(dir_abundances, dir_models)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Paleta de colores unificada oficial de tu v9
clin_colors <- c(
  "CTRL"         = "#2c3e50", 
  "RSV-/Wheeze-" = "#3498db", 
  "RSV-/Wheeze+" = "#2ecc71", 
  "RSV+/Wheeze-" = "#e67e22", 
  "RSV+/Wheeze+" = "#e74c3c"  
)

# ==============================================================================
# FUNCIONES AUXILIares NATIVAS (DE TU SCRIPT V9)
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
# 2. CARGA, AGREGACIÓN A GÉNERO Y FILTRADO PAREADO ESTRICTO (CORREGIDO)
# ==============================================================================
cat("\n[2] Procesando y emparejando matriz trans-nicho...\n")
data_input_dir <- "COMPLETE_REANALISIS_BQL/REANALISIS/Global/curated_data_global/"
pseq_raw_master <- readRDS(paste0(data_input_dir, "phyloseq_RAW_curated_global.rds"))

# Encontrar columna Genus automáticamente
genus_col <- colnames(tax_table(pseq_raw_master))[grepl("^[Gg]enus$", colnames(tax_table(pseq_raw_master)))][1]

# Consolidar a género y renombrar según tus funciones v9
pseq_genus <- tax_glom(pseq_raw_master, taxrank = genus_col, NArm = FALSE) %>% 
  rename_to_genus(genus_col)

# Transformación de Abundancia Relativa Log-estabilizada
pseq_prop <- transform_sample_counts(pseq_genus, function(x) log10((x / sum(x)) * 100 + 1))

# Extraer metadatos unificados
metadata_df <- extract_meta(pseq_prop) %>%
  mutate(
    Age = as.numeric(Age),
    Patient_ID  = as.factor(Patient_ID),
    Sample_Type = ifelse(Sample_Type %in% c("ANF", "NPA"), "NPA", "GUT")
  )

# --- NUEVO FILTRADO ESTILO SECCIÓN 6 ---
# 1. Definir que todos los que no son CTRL tienen Bronquiolitis
metadata_df$Bronchiolitis <- ifelse(metadata_df$Condicion_Clinica == "CTRL", "No", "Yes")

# 2. Definir Wheezing_Bin (1 si termina en Wheeze+, 0 si es Wheeze-)
metadata_df$Wheezing_Bin <- ifelse(
  grepl("\\+$", metadata_df$Condicion_Clinica) | metadata_df$Condicion_Clinica %in% c("RSV-/Wheeze+", "RSV+/Wheeze+"), 
  1, 0
)

# 3. Filtrar para quedarse ÚNICAMENTE con los Bronchiolitis == "Yes"
df_model_bql_pre <- metadata_df %>% filter(Bronchiolitis == "Yes")

# Identificar pacientes BQL+ estrictamente pareados (tienen NPA y GUT)
patients_with_both <- df_model_bql_pre %>% 
  group_by(Patient_ID) %>%
  filter(all(c("NPA", "GUT") %in% Sample_Type)) %>% 
  pull(Patient_ID) %>% 
  unique()

# Extraer matriz de abundancias transpuesta (Filas = Taxones, Columnas = Muestras)
otu_mat <- as.data.frame(t(as(otu_table(pseq_prop), "matrix")))

# Buscar tus dos taxones significativos
target_npa_clean <- rownames(otu_mat)[grep("Burkholderia", rownames(otu_mat), ignore.case = TRUE)[1]]
target_gut_clean <- rownames(otu_mat)[grep("Ruminococcus_gnavus", rownames(otu_mat), ignore.case = TRUE)[1]]

print(paste(">> Taxón detectado para NPA:", target_npa_clean))
print(paste(">> Taxón detectado para GUT:", target_gut_clean))

# Extraer los vectores de abundancia como numéricos indexados por el ID de la muestra
vec_burkholderia <- as.numeric(otu_mat[target_npa_clean, ])
names(vec_burkholderia) <- colnames(otu_mat)

vec_ruminococcus <- as.numeric(otu_mat[target_gut_clean, ])
names(vec_ruminococcus) <- colnames(otu_mat)

# Inyectar abundancias directas usando el subconjunto de Bronquiolitis
df_model_bql_pre$Burkholderia_NPA <- vec_burkholderia[df_model_bql_pre$Sample_ID]
df_model_bql_pre$Ruminococcus_GUT <- vec_ruminococcus[df_model_bql_pre$Sample_ID]

# Dividir por nicho para hacer el pivote/merge horizontal (formato WIDE por paciente)
df_npa_side <- df_model_bql_pre %>% 
  filter(Sample_Type == "NPA" & Patient_ID %in% patients_with_both) %>%
  select(Patient_ID, Age, Sex..Male., Wheezing_Bin, Burkholderia_NPA)

df_gut_side <- df_model_bql_pre %>% 
  filter(Sample_Type == "GUT" & Patient_ID %in% patients_with_both) %>%
  select(Patient_ID, Ruminococcus_GUT)

# Combinar en el data frame final cruzado (Cambiamos el nombre de la variable para el modelo)
df_final_modelo <- df_npa_side %>%
  dplyr::inner_join(df_gut_side, by = "Patient_ID") %>%
  drop_na(Burkholderia_NPA, Ruminococcus_GUT, Age, Wheezing_Bin)

cat(paste(" >> [SUCCESS] Cohorte BQL+ cruzada finalizada con", nrow(df_final_modelo), "pacientes perfectamente pareados.\n"))

# ==============================================================================
# 3. MODELO CONJUNTO MULTIVARIANTE Y FOREST PLOT (SECCIÓN 6 STYLE)
# ==============================================================================
cat("\n[3] Ajustando modelo logístico multivariante para Wheezing (BQL+)...\n")

# CAMBIO DE VARIABLE: Ahora predecimos Wheezing_Bin en lugar de Is_Severe
mv_model <- glm(Wheezing_Bin ~ Burkholderia_NPA + Ruminococcus_GUT + Age, 
                data = df_final_modelo, family = binomial)

# Extraer coeficientes exponenciados para obtener Odds Ratios (OR) directos
model_hits <- tidy(mv_model, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = case_when(
      term == "Burkholderia_NPA" ~ "Burkholderia-Caballeronia (NPA)",
      term == "Ruminococcus_GUT" ~ "[Ruminococcus] gnavus group (GUT)",
      term == "Age"              ~ "Edad (Años)",
      TRUE ~ term
    )
  )

# Renderizar el Forest Plot Premium (Actualizado títulos)
plot_forest <- ggplot(model_hits, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#e74c3c", size = 0.7) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15, color = "#2c3e50", size = 0.8) +
  geom_point(color = "#3498db", size = 4) + # Color azul para ir a juego con Wheezing de tu v9
  theme_bw(base_size = 13) +
  labs(
    title = "Predictors of Wheezing in Bronchiolitis Patients",
    subtitle = "Integrated Cross-Niche Paired GLM (Gut-Lung Axis)",
    x = "Odds Ratio (95% Confidence Interval)",
    y = ""
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "italic", hjust = 0.5),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray96")
  )

ggsave(filename = paste0(dir_models, "ForestPlot_Wheezing_BQL_Model.png"), 
       plot = plot_forest, width = 7.5, height = 4.5, dpi = 300)

# ==============================================================================
# 4. EVALUACIÓN PREDICTIVA: MULTIVARIATE ROC CURVE
# ==============================================================================
cat("\n[4] Calculando rendimiento diagnóstico del modelo (ROC/AUC)...\n")

# CAMBIO: Usamos df_final_modelo y Wheezing_Bin
df_final_modelo$pred_prob <- predict(mv_model, type = "response")
roc_obj <- roc(df_final_modelo$Wheezing_Bin, df_final_modelo$pred_prob, quiet = TRUE)
auc_val <- round(auc(roc_obj), 3)

df_roc_curve <- data.frame(
  Sensitivity = roc_obj$sensitivities,
  Specificity = roc_obj$specificities
)

plot_roc <- ggplot(df_roc_curve, aes(x = 1 - Specificity, y = Sensitivity)) +
  geom_path(color = "#2ecc71", size = 1.3) + # Verde premium
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60") +
  theme_bw(base_size = 13) +
  annotate("text", x = 0.65, y = 0.2, 
           label = paste0("Joint AUC = ", auc_val), 
           size = 5, fontface = "bold", color = "#2c3e50") +
  labs(
    title = "ROC Curve: Wheezing Risk Model",
    subtitle = "Combined Biomarkers: Burkholderia (NPA) + Ruminococcus (GUT) + Age",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 9, face = "italic", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

ggsave(filename = paste0(dir_models, "ROC_Wheezing_BQL_Model.png"), 
       plot = plot_roc, width = 6, height = 5.5, dpi = 300)

# Guardar coeficientes numéricos en CSV
write.csv(model_hits, file = paste0(dir_models, "Table_Wheezing_Predictors_Coefficients.csv"), row.names = FALSE)

cat("\n >> [PIPELINE FINISHED] Cross-niche risk panel for Wheezing executed successfully.\n")
library(ggplot2)
library(dplyr)
library(broom)
library(pROC)
library(patchwork)

# ==============================================================================
# 1. PREPARACIÓN DE DATOS PARA GRÁFICOS (SUPP FIG 7 - CORREGIDA)
# ==============================================================================
# Extraemos las probabilidades y generamos el objeto ROC real
df_final_modelo$pred_prob <- predict(mv_model, type = "response")
roc_obj <- roc(df_final_modelo$Wheezing_Bin, df_final_modelo$pred_prob, quiet = TRUE)

# Calcular IC de la curva ROC usando bootstrap/sensibilidad
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

# Calcular el intervalo de confianza matemático real del AUC (Método DeLong)
real_ci <- ci.auc(roc_obj)
ci_low  <- round(real_ci[1], 3)
ci_high <- round(real_ci[3], 3)
auc_val <- round(auc(roc_obj), 3)

# Formateamos las etiquetas con saltos de línea (\n) y nombres completos
tabla_plot_forest <- data.frame(
  term = c(
    "Burkholderia-\nCaballeronia-\nParaburkholderia\n(NPA)", 
    "[Ruminococcus]\ngnavus group\n(GUT)", 
    "Age (Months)"
  ),
  estimate = c(0.290828143784031, 3.80651971369552, 0.890712616318975),
  conf.low = c(0.106397583317266, 1.24562018832318, 0.760608351046915),
  conf.high = c(0.715908292116628, 14.875050817664, 1.03441135126556)
)

# Forzar el orden de los factores para que no se altere en el gráfico
tabla_plot_forest$term <- factor(tabla_plot_forest$term, levels = rev(tabla_plot_forest$term))

# ==============================================================================
# PANEL A: FOREST PLOT (SIN SUBTÍTULO, ETIQUETAS CORREGIDAS)
# ==============================================================================
plot_7a <- ggplot(tabla_plot_forest, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#e74c3c", size = 0.7) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15, color = "#2c3e50", size = 0.8) +
  geom_point(color = "#3498db", size = 4) +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 15)) +
  theme_bw(base_size = 12) +
  labs(
    title = "",
    subtitle = NULL, # Eliminado según requerimiento
    x = "Odds Ratio (95% Confidence Interval)",
    y = ""
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13, color = "black"),
    axis.text.y = element_text(face = "bold", color = "black", lineheight = 0.9), # Ajuste del interlineado de la etiqueta
    panel.grid.minor = element_blank()
  )

# ==============================================================================
# PANEL B: CURVA ROC (SENSITIVITY VS 1-SPECIFICITY, TEXTO MOVIDO Y REAJUSTADO)
# ==============================================================================
plot_7b <- ggplot() +
  geom_ribbon(data = df_roc_ci, aes(x = 1 - Specificity, ymin = Lower, ymax = Upper), 
              fill = "#2ecc71", alpha = 0.15) +
  geom_path(data = df_roc_curve, aes(x = 1 - Specificity, y = Sensitivity), 
            color = "#2ecc71", size = 1.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60") +
  theme_bw(base_size = 12) +
  # Corregida posición (x=0.95, hjust=1) para alinear a la derecha dentro de la caja y evitar cortes
  annotate("text", x = 0.95, y = 0.18, 
           label = paste0("AUC = ", auc_val, "\n95% CI: [", ci_low, " - ", ci_high, "]"), 
           size = 4.5, fontface = "bold", color = "#2c3e50", hjust = 1) +
  labs(
    title = NULL,
    subtitle = NULL, # Eliminados subtítulos internos
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

# ==============================================================================
# COMBINACIÓN INTERACTIVA CON PATCHWORK (FORMATO VERTICAL CON PLOTS ETIQUETADOS)
# ==============================================================================
# Usamos plot_annotation para añadir el marcado de secciones A y B de manera limpia
supp_figure_7 <- (plot_7a / plot_7b) + 
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = list(c("a", "b"))) & 
  theme(plot.tag = element_text(face = "bold", size = 12))

# Guardar panel integrado
dir_models <- "COMPLETE_REANALISIS_BQL/REANALISIS/CROSS_NICHE/models/"
if(!dir.exists(dir_models)) dir.create(dir_models, recursive = TRUE)

ggsave(filename = paste0(dir_models, "Supplementary_Figure_7_Integrated_Panel.png"), 
       plot = supp_figure_7, width = 6.5, height = 8.5, dpi = 300)