# ==============================================================================
# SCRIPT REFINADO v8: FLUJO DE TRABAJO REESTRUCTURADO Y OPTIMIZADO
# ==============================================================================
suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(nlme)
  library(emmeans)
  library(vegan)
  library(SpiecEasi)
  library(igraph)
  library(scales)
  library(RColorBrewer)
  library(patchwork) 
  library(ggsignif)
})

set.seed(123)

# --- Configuración de Rutas de Trabajo ---
setwd("/mnt/usb/BQL/BQL_ANALYSIS")
output_dir_global <- "COMPLETE_REANALISIS_BQL/REANALISIS/Global/"
data_input_dir    <- "COMPLETE_REANALISIS_BQL/REANALISIS/Global/curated_data_global/"
output_dir_stats  <- "COMPLETE_REANALISIS_BQL/REANALISIS/Global/ESTADISTICA_CONJUNTA/"
output_dir_nets   <- "COMPLETE_REANALISIS_BQL/REANALISIS/NETWORKS/"

for (d in c(output_dir_stats, output_dir_nets, output_dir_global)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ==============================================================================
# PALETAS DE COLORES UNIFICADAS
# ==============================================================================
clin_colors <- c(
  "CTRL"         = "#2c3e50", 
  "RSV-/Wheeze-" = "#3498db", 
  "RSV-/Wheeze+" = "#2ecc71", 
  "RSV+/Wheeze-" = "#e67e22", 
  "RSV+/Wheeze+" = "#e74c3c"  
)

niche_colors <- c(
  "NPA" = "#2980b9", 
  "GUT" = "#27ae60"  
)

# ==============================================================================
# FUNCIONES AUXILIARES
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

# ==============================================================================
# 1. CARGA Y FILTRADO PAREADO DE DATOS MAESTROS
# ==============================================================================
cat("\n[1] Cargando y preparando datos maestros...\n")
pseq_raw  <- readRDS(paste0(data_input_dir, "phyloseq_RAW_curated_global.rds"))
pseq_rare <- readRDS(paste0(data_input_dir, "phyloseq_RAREFIED_global.rds"))

genus_col <- colnames(tax_table(pseq_rare))[grepl("^[Gg]enus$", colnames(tax_table(pseq_rare)))][1]

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

pseq_rare_genus <- tax_glom(pseq_rare, taxrank = genus_col, NArm = FALSE) %>% rename_to_genus(genus_col)

metadata_df <- extract_meta(pseq_rare_genus) %>%
  mutate(
    Condicion_Clinica = factor(Condicion_Clinica, levels = names(clin_colors)),
    Patient_ID  = as.factor(Patient_ID),
    Sample_Type = factor(Sample_Type, levels = c("NPA", "GUT"))
  )

patients_with_both <- metadata_df %>% group_by(Patient_ID) %>%
  filter(all(c("NPA", "GUT") %in% Sample_Type)) %>% pull(Patient_ID) %>% unique()

pseq_rare_genus <- subset_samples(pseq_rare_genus, Patient_ID %in% patients_with_both)
metadata_df <- extract_meta(pseq_rare_genus)
sample_data(pseq_rare_genus) <- sample_data(metadata_df)


# ==============================================================================
# 2. CÁLCULO DE DIVERSIDAD ALFA Y BETA (ENFOQUE DE MÁXIMO RIGOR ESTADÍSTICO)
# ==============================================================================
cat("\n[2] Procesando Diversidades con Modelos de Efectos Mixtos y PERMANOVA Pareado...\n")

# --- 2.1 ALFA DIVERSIDAD: Modelo Lineal Mixto (LME) ---
alpha_metrics <- estimate_richness(pseq_rare_genus, measures = "Shannon")
df_alpha <- merge(data.frame(Sample_ID = rownames(alpha_metrics), Shannon = alpha_metrics$Shannon),
                  metadata_df, by = "Sample_ID")

# El modelo ideal: Ajuste simultáneo por Edad y Tipo de muestra con intercepto aleatorio por Paciente
lme_alpha <- lme(Shannon ~ Sample_Type + Age, random = ~1|Patient_ID, data = df_alpha, method = "REML")

# Extracción de estadísticos para la tabla y la corrección gráfica
summary_lme <- summary(lme_alpha)$tTable
beta_age  <- summary_lme["Age", "Value"]

# Guardar tabla de Alfa Diversidad (LME)
write.csv(summary_lme, 
          file = paste0(output_dir_stats, "LME_Alpha_Diversity_Model.csv"), 
          row.names = TRUE)

# Guardar la variable corregida por edad basada en el coeficiente del modelo mixto
df_alpha$Shannon_Corrected <- df_alpha$Shannon - (beta_age * df_alpha$Age)


# --- 2.2 BETA DIVERSIDAD: PERMANOVA Secuencial Pareado y Dispersión ---
bray_dist <- phyloseq::distance(pseq_rare_genus, method = "bray")
ord_pcoa  <- ordinate(pseq_rare_genus, method="PCoA", distance="bray")
df_pcoa   <- plot_ordination(pseq_rare_genus, ord_pcoa, justDF = TRUE)

# Estructura de permutación pareada: Restringe los cambios ÚNICAMENTE dentro del mismo paciente
perm_structure <- how(blocks = metadata_df$Patient_ID, nperm = 999)

# PERMANOVA: Ajusta primero por edad (secuencial), controlando el diseño pareado
permanova_res <- adonis2(bray_dist ~ Age + Sample_Type, 
                         data = metadata_df, 
                         permutations = perm_structure,
                         by = "terms") # Evalúa en orden: Age primero, luego Sample_Type

# Guardar tabla de PERMANOVA
write.csv(as.data.frame(permanova_res), 
          file = paste0(output_dir_stats, "PERMANOVA_Beta_Diversity_Pareado.csv"), 
          row.names = TRUE)

# Homogeneidad de la dispersión de grupos (Betadisper)
dispersion_res <- betadisper(bray_dist, group = metadata_df$Sample_Type)
permtest_disp  <- permutest(dispersion_res, permutations = perm_structure)

# Guardar tabla de Dispersión
write.csv(as.data.frame(permtest_disp$tab), 
          file = paste0(output_dir_stats, "Betadisper_Homogeneity_Test.csv"), 
          row.names = TRUE)


# --- 2.3 GRÁFICOS PARA PUBLICACIÓN (Limpios, sin títulos estadísticos) ---

p_alpha_box <- ggplot(df_alpha, aes(x = Sample_Type, y = Shannon_Corrected, fill = Sample_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1.5, color = "#2c3e50") +
  scale_fill_manual(values = niche_colors) +
  theme_bw(base_size = 12) +
  labs(x = "Niche / Sample Type", 
       y = "Shannon Index (Age-Corrected)") +
  theme(axis.title = element_text(face = "bold", size = 11),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_text(face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p_beta_pcoa <- ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, color = Sample_Type)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(fill = Sample_Type), geom = "polygon", alpha = 0.1, type = "t", linetype = 2) +
  scale_color_manual(values = niche_colors, name = "Niche Location") +
  scale_fill_manual(values = niche_colors, name = "Niche Location") +
  theme_bw(base_size = 12) +
  labs(x = paste0("PCoA 1 (", round(ord_pcoa$values$Relative_eig[1] * 100, 1), "%)"),
       y = paste0("PCoA 2 (", round(ord_pcoa$values$Relative_eig[2] * 100, 1), "%)")) +
  theme(axis.title = element_text(face = "bold", size = 11),
        axis.text = element_text(color = "black", size = 10),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# --- 2.4 Ensamblado Final con Patchwork y Guardado (300 DPI) ---
panel_joint <- p_alpha_box + p_beta_pcoa + 
  plot_layout(ncol = 2, widths = c(1, 1.2)) + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(paste0(output_dir_global, "Plot_Diversity_Simplified_Panel.png"), 
       plot = panel_joint, width = 10, height = 4.5, dpi = 300)

# ==============================================================================
# 3. REDES DE CADA NICHO (INTRA-NICHO CO-ABUNDANCIA)
# ==============================================================================
cat("\n[3] Módulo 3: Construcción de Redes Intra-Nicho por separado (Generación SpiecEasi)...\n")

prevalence_threshold <- 0.10
pseq_filt <- filter_taxa(pseq_rare_genus,
                         function(x) sum(x>0) >= (prevalence_threshold * nsamples(pseq_rare_genus)),
                         prune = TRUE)
taxa_names(pseq_filt) <- clean_taxa_names(taxa_names(pseq_filt))

pseq_npa <- subset_samples(pseq_filt, Sample_Type == "NPA") %>% prune_samples(sample_names(.), .)
pseq_gut <- subset_samples(pseq_filt, Sample_Type == "GUT") %>% prune_samples(sample_names(.), .)

remove_zero_var <- function(ps) {
  otu_m <- as(otu_table(ps), "matrix")
  keep <- if (taxa_are_rows(ps)) apply(otu_m, 1, var) > 0 else apply(otu_m, 2, var) > 0
  return(prune_taxa(keep, ps))
}
pseq_npa <- remove_zero_var(pseq_npa)
pseq_gut <- remove_zero_var(pseq_gut)

extract_otu_safe <- function(ps) {
  mat <- otu_table(ps)
  if (taxa_are_rows(ps)) mat <- t(mat)
  otu_standard <- as.matrix(mat)
  # Forzamos que las filas queden identificadas directamente por Patient_ID para emparejamiento seguro
  rownames(otu_standard) <- as.character(sample_data(ps)$Patient_ID)
  return(otu_standard)
}

otu_npa <- extract_otu_safe(pseq_npa)
otu_gut <- extract_otu_safe(pseq_gut)

arrange_rows_by_name <- function(mat) mat[order(rownames(mat)), , drop=FALSE]

# Aseguramos alineación perfecta de pacientes entre matrices intra-nicho antes de modelar
common_patients_internal <- intersect(rownames(otu_npa), rownames(otu_gut))
otu_npa <- otu_npa[common_patients_internal, , drop=FALSE] %>% arrange_rows_by_name() # orden alfabético estable
otu_gut <- otu_gut[common_patients_internal, , drop=FALSE] %>% arrange_rows_by_name()

otu_npa <- arrange_rows_by_name(otu_npa)
otu_gut <- arrange_rows_by_name(otu_gut)

run_spiec <- function(otu_mat, label) {
  cat(sprintf("  >> SpiecEasi [%s]...\n", label))
  spiec.easi(otu_mat, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50, ncores=1))
}
se_npa <- run_spiec(otu_npa, "NPA")
se_gut <- run_spiec(otu_gut, "GUT")

extract_igraph <- function(se_obj, otu_mat, label) {
  beta_raw <- getOptBeta(se_obj)
  if (is.null(beta_raw) || length(beta_raw)==0) {
    p <- ncol(otu_mat)
    beta_sym <- matrix(0, p, p, dimnames=list(colnames(otu_mat), colnames(otu_mat)))
  } else {
    beta_sym <- as.matrix(symBeta(beta_raw, mode="maxabs"))
    p <- ncol(otu_mat)
    if (nrow(beta_sym) != p || ncol(beta_sym) != p) {
      beta_new <- matrix(0, p, p, dimnames=list(colnames(otu_mat), colnames(otu_mat)))
      if (!is.null(colnames(beta_sym)) && all(colnames(beta_sym) %in% colnames(otu_mat))) {
        idx <- match(colnames(beta_sym), colnames(otu_mat))
        beta_new[idx, idx] <- beta_sym
      } else {
        m <- min(p, nrow(beta_sym))
        beta_new[1:m,1:m] <- beta_sym[1:m,1:m]
      }
      beta_sym <- beta_new
    }
    diag(beta_sym) <- 0
    rownames(beta_sym) <- colnames(beta_sym) <- colnames(otu_mat)
  }
  g <- graph_from_adjacency_matrix(beta_sym, mode="undirected", weighted=TRUE, diag=FALSE)
  V(g)$name <- colnames(otu_mat)
  V(g)$label <- colnames(otu_mat)
  return(g)
}
g_npa <- extract_igraph(se_npa, otu_npa, "NPA")
g_gut <- extract_igraph(se_gut, otu_gut, "GUT")

detect_communities <- function(g, label) {
  g_sub <- delete_vertices(g, degree(g)==0)
  if (vcount(g_sub)==0) {
    V(g)$module <- 0
    return(list(graph=g, communities=NA))
  }
  comm <- cluster_louvain(g_sub, weights=abs(E(g_sub)$weight))
  mod_map <- setNames(membership(comm), V(g_sub)$name)
  V(g)$module <- mod_map[V(g)$name]
  V(g)$module[is.na(V(g)$module)] <- 0
  cat(sprintf("  >> [%s] %d módulos detectados por Louvain.\n", label, max(comm$membership, na.rm=TRUE)))
  return(list(graph=g, communities=comm))
}
res_npa <- detect_communities(g_npa, "NPA"); g_npa <- res_npa$graph
res_gut <- detect_communities(g_gut, "GUT"); g_gut <- res_gut$graph

identify_hubs <- function(g, top_pct=0.05) {
  dg <- degree(g)
  thresh <- quantile(dg, probs=1-top_pct)
  V(g)$is_hub <- dg >= thresh
  V(g)$hub_label <- ifelse(V(g)$is_hub, V(g)$name, "")
  return(g)
}
g_npa <- identify_hubs(g_npa)
g_gut <- identify_hubs(g_gut)

hub_table <- bind_rows(
  data.frame(Network="NPA", Genus=V(g_npa)$name[V(g_npa)$is_hub], Degree=degree(g_npa)[V(g_npa)$is_hub], Module=V(g_npa)$module[V(g_npa)$is_hub]),
  data.frame(Network="GUT", Genus=V(g_gut)$name[V(g_gut)$is_hub], Degree=degree(g_gut)[V(g_gut)$is_hub], Module=V(g_gut)$module[V(g_gut)$is_hub])
) %>% arrange(Network, desc(Degree))
write.csv(hub_table, paste0(output_dir_nets, "Hub_Genera_NPA_GUT.csv"), row.names=FALSE)


# ------------------------------------------------------------------------------
# 3.1 EXTRACCIÓN DE MÓDULOS Y EXTRACCIÓN DE EIGENVECTORS POR MÓDULO (EIGEN-SCORES)
# ------------------------------------------------------------------------------
cat("\n[3.1] Calculando Eigen-Scores ponderados por la topología interna del módulo...\n")
df_analysis <- data.frame(Patient_ID = rownames(otu_npa), stringsAsFactors = FALSE)

otu_npa_log <- log1p(otu_npa)
otu_gut_log <- log1p(otu_gut)

adj_npa_full <- as.matrix(as_adjacency_matrix(g_npa, attr="weight", sparse=FALSE))
adj_gut_full <- as.matrix(as_adjacency_matrix(g_gut, attr="weight", sparse=FALSE))

# --- Cómputo Eigen-Scores para Nasofaringe (NPA) ---
npa_modules <- sort(unique(V(g_npa)$module)); npa_modules <- npa_modules[npa_modules != 0]
for(m in npa_modules) {
  genera_in_mod <- V(g_npa)$name[V(g_npa)$module == m]
  genera_in_mod <- intersect(genera_in_mod, colnames(otu_npa))
  
  if(length(genera_in_mod) > 1) {
    beta_sub <- adj_npa_full[genera_in_mod, genera_in_mod, drop=FALSE]
    ev <- eigen(beta_sub)$vectors[,1]
    names(ev) <- genera_in_mod
    
    hub_mod <- genera_in_mod[which.max(degree(g_npa, v=genera_in_mod))]
    if(ev[hub_mod] < 0) ev <- ev * -1
    
    df_analysis[[paste0("NPA_Mod", m)]] <- as.matrix(otu_npa_log[, genera_in_mod]) %*% ev
  } else if (length(genera_in_mod) == 1) {
    df_analysis[[paste0("NPA_Mod", m)]] <- otu_npa_log[, genera_in_mod]
  }
}

# --- Cómputo Eigen-Scores para Intestino (GUT) ---
gut_modules <- sort(unique(V(g_gut)$module)); gut_modules <- gut_modules[gut_modules != 0]
for(m in gut_modules) {
  genera_in_mod = V(g_gut)$name[V(g_gut)$module == m]
  genera_in_mod = intersect(genera_in_mod, colnames(otu_gut))
  
  if(length(genera_in_mod) > 1) {
    beta_sub <- adj_gut_full[genera_in_mod, genera_in_mod, drop=FALSE]
    ev <- eigen(beta_sub)$vectors[,1]
    names(ev) <- genera_in_mod
    
    hub_mod <- genera_in_mod[which.max(degree(g_gut, v=genera_in_mod))]
    if(ev[hub_mod] < 0) ev <- ev * -1
    
    df_analysis[[paste0("GUT_Mod", m)]] <- as.matrix(otu_gut_log[, genera_in_mod]) %*% ev
  } else if (length(genera_in_mod) == 1) {
    df_analysis[[paste0("GUT_Mod", m)]] <- otu_gut_log[, genera_in_mod]
  }
}


# ------------------------------------------------------------------------------
# 3.2 TRADUCCIÓN CLÍNICA, MODELADO DE REGRESIÓN MÚLTIPLE Y HEATMAP DE CORRIDO
# ------------------------------------------------------------------------------
cat("\n[3.2] Desacoplando variables clínicas y ajustando Modelos de Regresión Múltiple...\n")

raw_metadata <- data.frame(as(sample_data(pseq_rare_genus), "data.frame"), stringsAsFactors = FALSE, check.names = FALSE)

# Extracción limpia y robusta de los componentes individuales de fenotipo y co-variables
metadata_cleaned <- raw_metadata %>%
  mutate(
    # Las 3 variables principales solicitadas (Binarias fijas a niveles No/Yes para estabilizar coeficientes)
    `Bronchiolitis`        = factor(ifelse(grepl("^Yes$|1", Bronchiolitis, ignore.case=TRUE), "Yes", "No"), levels=c("No", "Yes")),
    `RSV Infection`        = factor(ifelse(grepl("^Yes$|1|RSV\\+", Respiratory.syncytial.virus) | grepl("RSV\\+", Condicion_Clinica), "Yes", "No"), levels=c("No", "Yes")),
    `Wheezing`             = factor(ifelse(grepl("Wheeze\\+", Condicion_Clinica) | grepl("^Yes$|1", Wheezing.treatment, ignore.case=TRUE) | (suppressWarnings(as.numeric(Wheezing.count)) > 0), "Yes", "No"), levels=c("No", "Yes")),
    
    # Resto de co-variables fenotípicas de contexto
    `Cesarean Delivery`    = factor(ifelse(grepl("^Yes$|1", Cesarean.section, ignore.case=TRUE), "Yes", "No"), levels=c("No", "Yes")),
    `Previous Antibiotics` = factor(ifelse(grepl("^Yes$|1", Previous.antibiotics, ignore.case=TRUE), "Yes", "No"), levels=c("No", "Yes")),
    `Family History Atopy` = factor(ifelse(grepl("^Yes$|1", Family.history.atopy, ignore.case=TRUE), "Yes", "No"), levels=c("No", "Yes")),
    `Breastfeeding`        = factor(ifelse(grepl("^Yes$|1", Breastfeeding, ignore.case=TRUE), "Yes", "No"), levels=c("No", "Yes")),
    `Age (Months)`         = as.numeric(Age)
  ) %>%
  select(Patient_ID, `Bronchiolitis`, `RSV Infection`, `Wheezing`, `Age (Months)`, `Cesarean Delivery`, `Previous Antibiotics`, `Family History Atopy`, `Breastfeeding`) %>%
  distinct(Patient_ID, .keep_all = TRUE)

df_analysis_clean_base <- df_analysis %>% select(Patient_ID, starts_with("NPA_Mod"), starts_with("GUT_Mod"))

# CORRECCIÓN: Forzamos de forma explícita el uso de dplyr para evitar conflictos con 'mia'
df_analysis_extended   <- dplyr::inner_join(df_analysis_clean_base, metadata_cleaned, by="Patient_ID")

target_modules <- c(grep("^NPA_Mod", colnames(df_analysis_extended), value=TRUE), grep("^GUT_Mod", colnames(df_analysis_extended), value=TRUE))
predictors_list <- c("Bronchiolitis", "RSV Infection", "Wheezing", "Age (Months)", "Breastfeeding", "Cesarean Delivery", "Previous Antibiotics", "Family History Atopy")

# Escalamos la variable de edad para que su coeficiente beta sea directamente interpretable y comparable frente a los binarios
df_modelling <- df_analysis_extended
df_modelling$`Age (Months)` <- as.vector(scale(df_modelling$`Age (Months)`))

regression_results <- data.frame()

if(length(target_modules) > 0) {
  for(mod in target_modules) {
    # Ajustamos un único modelo lineal múltiple por módulo con todos los predictores controlados mutuamente
    formula_mod <- as.formula(paste0("`", mod, "` ~ Bronchiolitis + `RSV Infection` + Wheezing + `Age (Months)` + Breastfeeding + `Cesarean Delivery` + `Previous Antibiotics` + `Family History Atopy`"))
    fit <- tryCatch(lm(formula_mod, data=df_modelling), error = function(e) NULL)
    
    if(!is.null(fit)) {
      coef_matrix <- summary(fit)$coefficients
      for(pred in predictors_list) {
        # Añadimos fixed = TRUE para que los paréntesis de la Edad se lean como texto literal
        matched_row <- rownames(coef_matrix)[grepl(pred, rownames(coef_matrix), fixed = TRUE)]
        if(length(matched_row) > 0) {
          eff_size <- coef_matrix[matched_row[1], "Estimate"]
          p_val    <- coef_matrix[matched_row[1], "Pr(>|t|)"]
          
          regression_results <- rbind(regression_results, data.frame(
            Module = mod, Predictor = pred, Effect_Size = eff_size, P_Value = p_val, stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  # Corrección global de FDR (Benjamini-Hochberg) cruzando todos los coeficientes calculados
  regression_results <- regression_results %>% mutate(FDR_Value = p.adjust(P_Value, method = "fdr"))
}

# Generación automática de etiquetas premium (Valor de Beta + Estrellas de significación FDR)
regression_results <- regression_results %>%
  mutate(
    Significance_Label = case_when(
      FDR_Value < 0.001 ~ "***",
      FDR_Value < 0.01  ~ "**",
      FDR_Value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Cell_Text = sprintf("%.2f%s", Effect_Size, Significance_Label),
    Predictor = factor(Predictor, levels = predictors_list) # Fijamos orden estricto de columnas
  )

# --- REPRESENTACIÓN GRÁFICA: HEATMAP CONTINUO "DE CORRIDO" ---
plot_heatmap_expanded <- ggplot(regression_results, aes(x = Predictor, y = Module, fill = Effect_Size)) +
  geom_tile(color = "white", lwd = 0.5) +
  geom_text(aes(label = Cell_Text), color = "black", size = 2.8, fontface = "bold") +
  scale_fill_gradient2(low = "#2c3e50", mid = "#f1c40f", high = "#e74c3c", midpoint = 0, 
                       name = "Standardized\nEffect Size (β)") +
  theme_minimal(base_size = 11) +
  labs(title = "Standardized Effect Sizes (β) Across Clinical Traits",
             x = "Patient Clinical Traits & Contextual Covariates", 
       y = "Structural Network Discovery Modules") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", color = "#2c3e50"),
        axis.text.y = element_text(face = "bold", color = "#2c3e50"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 9))

ggsave(paste0(output_dir_nets, "Plot_Heatmap_Expanded_Clinical_Covariates.png"), plot = plot_heatmap_expanded, width = 11, height = 7, dpi = 300)


# --- EXPORTAR LA TABLA MAESTRA DE ASOCIACIÓN POR BACTERIA EN FORMATO WIDE (UNA FILA POR MÓDULO) ---
cat("  >> Compilando tabla maestra de coeficientes resumida por Módulo...\n")

# 1. Extraemos los grados (degree) de todas las bacterias para identificar el Hub de cada módulo
deg_npa <- data.frame(Genus = V(g_npa)$name, Niche = "Nasopharynx (NPA)", Degree = degree(g_npa), stringsAsFactors = FALSE)
deg_gut <- data.frame(Genus = V(g_gut)$name, Niche = "Gut (GUT)", Degree = degree(g_gut), stringsAsFactors = FALSE)
deg_master <- rbind(deg_npa, deg_gut)

taxa_npa_status <- data.frame(Genus = V(g_npa)$name, Niche = "Nasopharynx (NPA)", Structural_Module = V(g_npa)$module, stringsAsFactors = FALSE)
taxa_gut_status <- data.frame(Genus = V(g_gut)$name, Niche = "Gut (GUT)", Structural_Module = V(g_gut)$module, stringsAsFactors = FALSE)

# 2. Combinamos la taxonomía base, eliminamos el módulo de ruido (0) y cruzamos con sus grados
master_taxa_raw <- rbind(taxa_npa_status, taxa_gut_status) %>% 
  filter(Structural_Module != 0) %>%
  dplyr::left_join(deg_master, by = c("Genus", "Niche"))

# 3. Colapsamos la taxonomía: una sola línea por módulo con sus géneros unidos por espacios y su Hub
master_modules_collapsed <- master_taxa_raw %>%
  group_by(Niche, Structural_Module) %>%
  summarise(
    Genera_In_Module = paste(Genus, collapse = " "),
    Hub_Genus        = Genus[which.max(Degree)],
    .groups = "drop"
  )

# 4. Si existen resultados de las regresiones, pivotamos y fusionamos todo
if(exists("regression_results") && nrow(regression_results) > 0) {
  clinical_p_wide <- regression_results %>%
    mutate(String_Report = sprintf("%.3f (FDR: %.4f)", Effect_Size, FDR_Value)) %>%
    select(Module, Predictor, String_Report) %>%
    distinct(Module, Predictor, .keep_all = TRUE) %>%
    mutate(Module_Num = as.numeric(gsub("NPA_Mod|GUT_Mod", "", Module)),
           Niche_Prefix = ifelse(grepl("^NPA", Module), "Nasopharynx (NPA)", "Gut (GUT)")) %>%
    pivot_wider(id_cols = c(Module_Num, Niche_Prefix), names_from = Predictor, values_from = String_Report, names_prefix = "Beta_")
  
  # Combinamos las betas clínicas con nuestra estructura de módulos colapsados
  master_clinical_report <- master_modules_collapsed %>%
    dplyr::left_join(clinical_p_wide, by = c("Structural_Module" = "Module_Num", "Niche" = "Niche_Prefix")) %>%
    mutate(Structural_Module = paste0(ifelse(Niche == "Nasopharynx (NPA)", "NPA_Mod", "GUT_Mod"), Structural_Module)) %>%
    arrange(desc(Niche), Structural_Module)
  
  # 5. Reorganizamos los elementos para que los géneros y el Hub queden al final de todo de forma estricta
  master_clinical_report <- master_clinical_report %>%
    relocate(Genera_In_Module, .after = last_col()) %>%
    relocate(Hub_Genus, .after = last_col())
  
  # Guardamos el archivo final refinado
  write.csv(master_clinical_report, paste0(output_dir_nets, "Master_Taxa_Module_Clinical_Associations_FINAL.csv"), row.names = FALSE)
  cat("  >> ¡Tabla maestra optimizada guardada con éxito! (1 fila por módulo, géneros agrupados y columna Hub al final).\n")
}


# ------------------------------------------------------------------------------
# 3.3 PLOTS DE GRAFOS INDEPENDIENTES DE CADA NICHO (NPA Y GUT, SIN MÓDULO 0)
# ------------------------------------------------------------------------------
cat("\n[3.3] Generando visualizaciones premium de Redes Intra-Nicho Disecadas (Sin Módulo 0)...\n")

log_mean_npa <- log1p(colMeans(otu_npa, na.rm = TRUE))
log_mean_gut <- log1p(colMeans(otu_gut, na.rm = TRUE))

# --- PLOT INDEPENDIENTE: NASOFARINGE (NPA) ---
g_npa_clean <- delete_vertices(g_npa, V(g_npa)[V(g_npa)$module == 0])
V(g_npa_clean)$label <- V(g_npa_clean)$name
E(g_npa_clean)$color <- ifelse(E(g_npa_clean)$weight > 0, "#2ecc71DD", "#e74c3cDD")
E(g_npa_clean)$width <- abs(E(g_npa_clean)$weight) * 5.5
V(g_npa_clean)$size <- 3 + (log_mean_npa[V(g_npa_clean)$name] * 1.2)

unique_npa_mods <- sort(unique(V(g_npa_clean)$module))
color_map_npa <- setNames(paste0(colorRampPalette(brewer.pal(8, "Set1"))(length(unique_npa_mods)), "CC"), unique_npa_mods)
V(g_npa_clean)$color_node <- color_map_npa[as.character(V(g_npa_clean)$module)]

out_npa_dissected <- paste0(output_dir_nets, "Network_Dissected_Internal_NPA_NoMod0.png")
png(out_npa_dissected, width=3200, height=2800, res=300)
layout(matrix(c(1,2), nrow=1), widths=c(3,1))
par(mar=c(2, 2, 5, 2), bg="white")
set.seed(42)
plot(g_npa_clean, layout=layout_with_fr(g_npa_clean, weights=abs(E(g_npa_clean)$weight)), 
     vertex.color=V(g_npa_clean)$color_node, vertex.size=V(g_npa_clean)$size, vertex.frame.color="#2c3e50", vertex.frame.width=0.7,
     vertex.label=V(g_npa_clean)$label, vertex.label.cex=0.50, vertex.label.font=4, vertex.label.color="#2c3e50",
     main="Nasopharynx (NPA) Connected Network (Excluding Mod0)")
par(mar=c(2, 0, 5, 1)); plot.new()
legend("topleft", legend=paste("NPA Mod", names(color_map_npa)), col=color_map_npa, pch=19, bty="n", title="Modules", cex=1.0, title.font=2)
legend("center", legend=c("Positive Alignment (+)", "Mutual Exclusion (-)"), col=c("#2ecc71", "#e74c3c"), lty=1, lwd=4, bty="n", title="Interaction Type", cex=1.0, title.font=2)
dev.off()

# --- PLOT INDEPENDIENTE: INTESTINO (GUT) ---
g_gut_clean <- delete_vertices(g_gut, V(g_gut)[V(g_gut)$module == 0])
V(g_gut_clean)$label <- V(g_gut_clean)$name
E(g_gut_clean)$color <- ifelse(E(g_gut_clean)$weight > 0, "#2ecc71DD", "#e74c3cDD")
E(g_gut_clean)$width <- abs(E(g_gut_clean)$weight) * 5.5
V(g_gut_clean)$size <- 3 + (log_mean_gut[V(g_gut_clean)$name] * 1.2)

unique_gut_mods <- sort(unique(V(g_gut_clean)$module))
color_map_gut <- setNames(paste0(colorRampPalette(brewer.pal(8, "Dark2"))(length(unique_gut_mods)), "CC"), unique_gut_mods)
V(g_gut_clean)$color_node <- color_map_gut[as.character(V(g_gut_clean)$module)]

out_gut_dissected <- paste0(output_dir_nets, "Network_Dissected_Internal_GUT_NoMod0.png")
png(out_gut_dissected, width=3200, height=2800, res=300)
layout(matrix(c(1,2), nrow=1), widths=c(3,1))
par(mar=c(2, 2, 5, 2), bg="white")
set.seed(88)
plot(g_gut_clean, layout=layout_with_fr(g_gut_clean, weights=abs(E(g_gut_clean)$weight)),
     vertex.color=V(g_gut_clean)$color_node, vertex.size=V(g_gut_clean)$size, vertex.frame.color="#2c3e50", vertex.frame.width=0.7,
     vertex.label=V(g_gut_clean)$label, vertex.label.cex=0.50, vertex.label.font=4, vertex.label.color="#2c3e50",
     main="Gut (GUT) Connected Network (Excluding Mod0)")
par(mar=c(2, 0, 5, 1)); plot.new()
legend("topleft", legend=paste("GUT Mod", names(color_map_gut)), col=color_map_gut, pch=19, bty="n", title="Modules", cex=1.0, title.font=2)
legend("center", legend=c("Positive Alignment (+)", "Mutual Exclusion (-)"), col=c("#2ecc71", "#e74c3c"), lty=1, lwd=4, bty="n", title="Interaction Type", cex=1.0, title.font=2)
dev.off()


# ==============================================================================
# 4. ANÁLISIS DE RED INTER-NICHO (CROSS-DOMAIN TALK NPA <-> GUT)
# ==============================================================================
cat("\n[4] Módulo 4: Ejecutando Análisis de Red Cross-Domain (Eje Aerodigestivo NPA <-> GUT)...\n")

# Re-confirmamos matrices pareadas estrictas
otu_npa_paired <- as.matrix(otu_npa)
otu_gut_paired <- as.matrix(otu_gut)

se_cross <- spiec.easi(list(otu_npa_paired, otu_gut_paired), method="mb", 
                       lambda.min.ratio=1e-2, nlambda=20, 
                       pulsar.params=list(rep.num=50, ncores=1))

beta_cross <- as.matrix(symBeta(getOptBeta(se_cross), mode="maxabs"))
diag(beta_cross) <- 0

n_npa <- ncol(otu_npa_paired)
n_gut <- ncol(otu_gut_paired)
total <- n_npa + n_gut

# --- EXPORTACIÓN DE TABLA DE INTERACCIONES INTER-NICHO SIGNIFICATIVAS ---
block_cross_npa_gut <- beta_cross[1:n_npa, (n_npa+1):total]
edges_cross <- which(block_cross_npa_gut != 0, arr.ind = TRUE)

if (nrow(edges_cross) > 0) {
  cross_df <- data.frame(
    Genus_NPA = colnames(otu_npa_paired)[edges_cross[,1]],
    Genus_GUT = colnames(otu_gut_paired)[edges_cross[,2]],
    Weight    = block_cross_npa_gut[edges_cross]
  ) %>% arrange(desc(abs(Weight)))
  write.csv(cross_df, paste0(output_dir_nets, "CrossDomain_Associations_NPA_GUT.csv"), row.names = FALSE)
  cat("  >> Tabla de interacciones guardada: CrossDomain_Associations_NPA_GUT.csv\n")
}

# --- CONSTRUCCIÓN Y REPRESENTACIÓN GRÁFICA DEL GRAFO CROSS-DOMAIN ---
g_cross <- graph_from_adjacency_matrix(beta_cross, mode="undirected", weighted=TRUE, diag=FALSE)
V(g_cross)$name  <- c(colnames(otu_npa_paired), colnames(otu_gut_paired))
V(g_cross)$tissue <- c(rep("NPA", n_npa), rep("GUT", n_gut))

g_cross_sub <- delete_vertices(g_cross, which(degree(g_cross) == 0))

if (vcount(g_cross_sub) > 0 && ecount(g_cross_sub) > 0) {
  v_tissue <- V(g_cross_sub)$tissue
  el <- as_edgelist(g_cross_sub, names = FALSE)
  is_cross_edge <- (v_tissue[el[,1]] != v_tissue[el[,2]])
  
  E(g_cross_sub)$weight_abs <- abs(E(g_cross_sub)$weight)
  E(g_cross_sub)$width <- ifelse(is_cross_edge, 4.0, 0.7) 
  E(g_cross_sub)$color <- ifelse(is_cross_edge,
                                 ifelse(E(g_cross_sub)$weight > 0, "#2ecc71BB", "#e74c3cBB"), 
                                 ifelse(E(g_cross_sub)$weight > 0, "#2ecc7115", "#e74c3c15")) 
  
  degree_cross_only <- adjacent_vertices(g_cross_sub, V(g_cross_sub))
  nodes_with_cross_interactions <- sapply(1:vcount(g_cross_sub), function(i) {
    any(v_tissue[degree_cross_only[[i]]] != v_tissue[i])
  })
  V(g_cross_sub)$label <- ifelse(nodes_with_cross_interactions, V(g_cross_sub)$name, "")
  
  out_path_premium <- paste0(output_dir_nets, "Network_CrossDomain_SIMPLIFIED.png")
  png(out_path_premium, width=3200, height=2800, res=300)
  layout(matrix(c(1,2), nrow=1), widths=c(3,1))
  par(mar=c(2, 2, 5, 2), bg="white")
  
  set.seed(123)
  plot(g_cross_sub, 
       layout=layout_with_fr(g_cross_sub, weights=E(g_cross_sub)$weight_abs),
       vertex.color=ifelse(V(g_cross_sub)$tissue == "NPA", "#2980b9", "#27ae60"), 
       vertex.size=6, vertex.frame.color="#2c3e50", vertex.frame.width=1,
       vertex.label=V(g_cross_sub)$label, vertex.label.cex=0.6, vertex.label.font=4, vertex.label.color="#2c3e50",
       main="Aerodigestive Axis: Inter-Niche Translocation Slopes (NPA <-> GUT)")
  
  par(mar=c(2, 0, 5, 1)); plot.new()
  legend("topleft", legend=c("Nasopharynx (NPA)", "Gut (GUT)"), col=c("#2980b9", "#27ae60"), pch=21, pt.bg=c("#2980b9", "#27ae60"), bty="n", title="Niche Source Location", cex=1.1, title.font=2)
  legend("center", legend=c("Positive Alignment (+)", "Negative Alignment (-)"), col=c("#2ecc71", "#e74c3c"), lty=1, lwd=4, bty="n", title="Inter-Niche Edge Profile", cex=1.1, title.font=2)
  dev.off()
}

cat("\n[Pipeline Finalizado de forma Exitosa. Estructura ordenada y robusta].\n")
