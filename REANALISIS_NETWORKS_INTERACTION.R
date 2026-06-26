# ==============================================================================
# SCRIPT REFINADO v9: AMPLIADO CON MULTIPLES METRICAS Y PROPORCIONES PREMIUM
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
# 2. CÁLCULO DE DIVERSIDAD ALFA Y BETA
# ==============================================================================
cat("\n[2] Procesando Diversidades con Modelos de Efectos Mixtos y PERMANOVA Pareado...\n")

# --- 2.1 ALFA DIVERSIDAD: Modelo Lineal Mixto (LME) para múltiples métricas ---
alpha_metrics <- estimate_richness(pseq_rare_genus, measures = c("Observed", "Shannon", "Simpson"))
df_alpha <- merge(data.frame(Sample_ID = rownames(alpha_metrics), 
                             Observed = alpha_metrics$Observed,
                             Shannon = alpha_metrics$Shannon, 
                             Simpson = alpha_metrics$Simpson),
                  metadata_df, by = "Sample_ID")

# Inicializar dataframe para almacenar significancias de los corchetes
signif_list <- list()

for (metric in c("Observed", "Shannon", "Simpson")) {
  # Ajuste del modelo mixto para cada métrica
  lme_form <- as.formula(paste(metric, "~ Sample_Type + Age"))
  lme_alpha <- lme(lme_form, random = ~1|Patient_ID, data = df_alpha, method = "REML")
  
  summary_lme <- summary(lme_alpha)$tTable
  beta_age  <- summary_lme["Age", "Value"]
  
  # Guardar tabla de resultados
  write.csv(summary_lme, 
            file = paste0(output_dir_stats, "LME_Alpha_Diversity_", metric, "_Model.csv"), 
            row.names = TRUE)
  
  # Corrección por edad
  df_alpha[[paste0(metric, "_Corrected")]] <- df_alpha[[metric]] - (beta_age * df_alpha$Age)
  
  # Extraer p-valor corregido para el contraste post-hoc entre NPA y GUT
  emm <- emmeans(lme_alpha, specs = pairwise ~ Sample_Type, adjust = "fdr")
  p_val <- as.data.frame(emm$contrasts)$p.value
  
  # Convertir a formato asterisco
  annot_star <- "ns"
  if(p_val < 0.05)  annot_star <- "*"
  if(p_val < 0.01)  annot_star <- "**"
  if(p_val < 0.001) annot_star <- "***"
  
  signif_list[[metric]] <- data.frame(p = p_val, annot = annot_star)
}


# --- 2.2 BETA DIVERSIDAD: PERMANOVA Secuencial Pareado y Dispersión ---
bray_dist <- phyloseq::distance(pseq_rare_genus, method = "bray")
ord_pcoa  <- ordinate(pseq_rare_genus, method="PCoA", distance="bray")
df_pcoa   <- plot_ordination(pseq_rare_genus, ord_pcoa, justDF = TRUE)

# Estructura de permutación pareada
perm_structure <- how(blocks = metadata_df$Patient_ID, nperm = 999)

# PERMANOVA
permanova_res <- adonis2(bray_dist ~ Age + Sample_Type, 
                         data = metadata_df, 
                         permutations = perm_structure,
                         by = "terms")

write.csv(as.data.frame(permanova_res), 
          file = paste0(output_dir_stats, "PERMANOVA_Beta_Diversity_Pareado.csv"), 
          row.names = TRUE)

# Homogeneidad de la dispersión de grupos
dispersion_res <- betadisper(bray_dist, group = metadata_df$Sample_Type)
permtest_disp  <- permutest(dispersion_res, permutations = perm_structure)

write.csv(as.data.frame(permtest_disp$tab), 
          file = paste0(output_dir_stats, "Betadisper_Homogeneity_Test.csv"), 
          row.names = TRUE)


# ==============================================================================
# 3. GENERACIÓN DE GRÁFICOS INDEPENDIENTES PARA PANEL PREMIUM
# ==============================================================================

# Función para automatizar los boxplots individuales de Alpha
make_alpha_plot <- function(df, y_col, label_y, signif_df) {
  p <- ggplot(df, aes_string(x = "Sample_Type", y = y_col, fill = "Sample_Type")) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.5, color = "#2c3e50") +
    geom_jitter(width = 0.1, alpha = 0.4, size = 1.5, color = "#2c3e50") +
    scale_fill_manual(values = niche_colors) +
    theme_bw(base_size = 11) +
    labs(x = NULL, y = label_y) +
    theme(axis.title.y = element_text(face = "bold", size = 10),
          axis.text = element_text(color = "black", size = 9),
          axis.text.x = element_text(face = "bold"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#f0f0f0"))
  
  # Si el contraste es significativo, añadir corchete manual seguro
  if (signif_df$annot != "ns") {
    max_y <- max(df[[y_col]], na.rm = TRUE)
    range_y <- max_y - min(df[[y_col]], na.rm = TRUE)
    p <- p + geom_signif(
      comparisons = list(c("NPA", "GUT")),
      annotations = signif_df$annot,
      y_position = max_y + (range_y * 0.05),
      tip_length = 0.02,
      vjust = 0.5,
      textsize = 4,
      color = "black"
    )
  }
  return(p)
}

# Crear las 3 gráficas individuales de Alpha
p_observed <- make_alpha_plot(df_alpha, "Observed_Corrected", "Observed Richness\n(Age-Corrected)", signif_list[["Observed"]])
p_shannon  <- make_alpha_plot(df_alpha, "Shannon_Corrected", "Shannon Index\n(Age-Corrected)", signif_list[["Shannon"]])
p_simpson  <- make_alpha_plot(df_alpha, "Simpson_Corrected", "Simpson Index\n(Age-Corrected)", signif_list[["Simpson"]])

# --- Gráfico de Beta Diversidad (Ocupará todo el ancho abajo) ---
p_beta_pcoa <- ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, color = Sample_Type)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(fill = Sample_Type), geom = "polygon", alpha = 0.1, type = "t", linetype = 2) +
  scale_color_manual(values = niche_colors, name = "Niche Location") +
  scale_fill_manual(values = niche_colors, name = "Niche Location") +
  theme_bw(base_size = 11) +
  labs(x = paste0("PCoA 1 (", round(ord_pcoa$values$Relative_eig[1] * 100, 1), "%)"),
       y = paste0("PCoA 2 (", round(ord_pcoa$values$Relative_eig[2] * 100, 1), "%)")) +
  theme(axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(color = "black", size = 9),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#f0f0f0"),
        legend.position = "right")

# ==============================================================================
# 4. ENSAMBLADO FINAL CON PATCHWORK (ESTRUCTURA 3 ARRIBA, 1 ABAJO)
# ==============================================================================
# Empaquetamos los tres de arriba en una sola fila
top_row <- p_observed + p_shannon + p_simpson + plot_layout(ncol = 3)

# Ensamblado final vertical: fila superior (alpha) / fila inferior (beta)
# Cambiamos las proporciones para que el PCoA tenga un poco más de aire (0.7 a 1)
panel_joint <- top_row / p_beta_pcoa + 
  plot_layout(heights = c(0.7, 1)) + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = "bold", size = 14, color = "#2c3e50"))

# Guardar con dimensiones optimizadas para evitar solapamientos verticales
ggsave(paste0(output_dir_global, "Figure_Supplementary_1_Diversity_Expanded_Panel.png"), 
       plot = panel_joint, width = 11, height = 8.5, dpi = 300)

cat(" >> ¡Análisis de diversidad expandido y panel guardado exitosamente!\n")

# ==============================================================================
# 3. REDES DE CADA NICHO (INTRA-NICHO CO-ABUNDANCIA OPTIMIZADA)
# ==============================================================================
cat("\n[3] Módulo 3: Construcción de Redes Intra-Nicho por separado (Generación SpiecEasi NATIVA)...\n")

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
  rownames(otu_standard) <- as.character(sample_data(ps)$Patient_ID)
  return(otu_standard)
}

otu_npa <- extract_otu_safe(pseq_npa)
otu_gut <- extract_otu_safe(pseq_gut)

arrange_rows_by_name <- function(mat) mat[order(rownames(mat)), , drop=FALSE]

common_patients_internal <- intersect(rownames(otu_npa), rownames(otu_gut))
otu_npa <- otu_npa[common_patients_internal, , drop=FALSE] %>% arrange_rows_by_name()
otu_gut <- otu_gut[common_patients_internal, , drop=FALSE] %>% arrange_rows_by_name()

run_spiec <- function(otu_mat, label) {
  cat(sprintf("  >> SpiecEasi [%s] con estimación Sparser MB...\n", label))
  spiec.easi(otu_mat, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50, ncores=1))
}
se_npa <- run_spiec(otu_npa, "NPA")
se_gut <- run_spiec(otu_gut, "GUT")

extract_igraph_robust <- function(se_obj, otu_mat) {
  refit_mat <- as.matrix(getRefit(se_obj))
  beta_raw  <- getOptBeta(se_obj)
  beta_sym  <- as.matrix(symBeta(beta_raw, mode="maxabs"))
  
  colnames(refit_mat) <- rownames(refit_mat) <- colnames(otu_mat)
  colnames(beta_sym)  <- rownames(beta_sym)  <- colnames(otu_mat)
  
  g <- graph_from_adjacency_matrix(beta_sym, mode="undirected", weighted=TRUE, diag=FALSE)
  V(g)$name  <- colnames(otu_mat)
  V(g)$label <- colnames(otu_mat)
  return(g)
}

g_npa <- extract_igraph_robust(se_npa, otu_npa)
g_gut <- extract_igraph_robust(se_gut, otu_gut)

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
  cat(sprintf("  >> [%s] %d módulos funcionales detectados por el algoritmo de Louvain.\n", label, max(comm$membership, na.rm=TRUE)))
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
# 3.1 EXTRACCIÓN DE MÓDULOS Y EXTRACCIÓN DE EIGENVECTORS (CON Z-SCORE IMPERATIVO)
# ------------------------------------------------------------------------------
cat("\n[3.1] Calculando Eigen-Scores ponderados y escalando a Z-scores estándar...\n")
df_analysis <- data.frame(Patient_ID = rownames(otu_npa), stringsAsFactors = FALSE)

otu_npa_log <- log1p(otu_npa)
otu_gut_log <- log1p(otu_gut)

adj_npa_full <- as.matrix(as_adjacency_matrix(g_npa, attr="weight", sparse=FALSE))
adj_gut_full <- as.matrix(as_adjacency_matrix(g_gut, attr="weight", sparse=FALSE))

# --- Nasofaringe (NPA) ---
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
    
    df_analysis[[paste0("NPA_Mod", m)]] <- as.vector(scale(as.matrix(otu_npa_log[, genera_in_mod]) %*% ev))
  } else if (length(genera_in_mod) == 1) {
    df_analysis[[paste0("NPA_Mod", m)]] <- as.vector(scale(otu_npa_log[, genera_in_mod]))
  }
}

# --- Intestino (GUT) ---
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
    
    df_analysis[[paste0("GUT_Mod", m)]] <- as.vector(scale(as.matrix(otu_gut_log[, genera_in_mod]) %*% ev))
  } else if (length(genera_in_mod) == 1) {
    df_analysis[[paste0("GUT_Mod", m)]] <- as.vector(scale(otu_gut_log[, genera_in_mod]))
  }
}


# ------------------------------------------------------------------------------
# 3.2 TRADUCCIÓN CLÍNICA Y MODELADO DE REGRESIÓN MÚLTIPLE (BETAS COMPARABLES)
# ------------------------------------------------------------------------------
cat("\n[3.2] Desacoplando variables clínicas y ajustando Modelos de Regresión Múltiple...\n")

raw_metadata <- data.frame(as(sample_data(pseq_rare_genus), "data.frame"), stringsAsFactors = FALSE, check.names = FALSE)

metadata_cleaned <- raw_metadata %>%
  mutate(
    `Bronchiolitis`        = factor(ifelse(grepl("^Yes$|1", Bronchiolitis, ignore.case=TRUE), "Yes", "No"), levels=c("No", "Yes")),
    `RSV Infection`        = factor(ifelse(grepl("^Yes$|1|RSV\\+", Respiratory.syncytial.virus) | grepl("RSV\\+", Condicion_Clinica), "Yes", "No"), levels=c("No", "Yes")),
    `Wheezing`             = factor(ifelse(grepl("Wheeze\\+", Condicion_Clinica) | grepl("^Yes$|1", Wheezing.treatment, ignore.case=TRUE) | (suppressWarnings(as.numeric(Wheezing.count)) > 0), "Yes", "No"), levels=c("No", "Yes")),
    `Cesarean Delivery`    = factor(ifelse(grepl("^Yes$|1", Cesarean.section, ignore.case=TRUE), "Yes", "No"), levels=c("No", "Yes")),
    `Previous Antibiotics` = factor(ifelse(grepl("^Yes$|1", Previous.antibiotics, ignore.case=TRUE), "Yes", "No"), levels=c("No", "Yes")),
    `Family History Atopy` = factor(ifelse(grepl("^Yes$|1", Family.history.atopy, ignore.case=TRUE), "Yes", "No"), levels=c("No", "Yes")),
    `Breastfeeding`        = factor(ifelse(grepl("^Yes$|1", Breastfeeding, ignore.case=TRUE), "Yes", "No"), levels=c("No", "Yes")),
    `Age (Months)`         = as.numeric(Age)
  ) %>%
  select(Patient_ID, `Bronchiolitis`, `RSV Infection`, `Wheezing`, `Age (Months)`, `Cesarean Delivery`, `Previous Antibiotics`, `Family History Atopy`, `Breastfeeding`) %>%
  distinct(Patient_ID, .keep_all = TRUE)

df_analysis_clean_base <- df_analysis %>% select(Patient_ID, starts_with("NPA_Mod"), starts_with("GUT_Mod"))
df_analysis_extended   <- dplyr::inner_join(df_analysis_clean_base, metadata_cleaned, by="Patient_ID")

target_modules  <- c(grep("^NPA_Mod", colnames(df_analysis_extended), value=TRUE), grep("^GUT_Mod", colnames(df_analysis_extended), value=TRUE))
predictors_list <- c("Bronchiolitis", "RSV Infection", "Wheezing", "Age (Months)", "Breastfeeding", "Cesarean Delivery", "Previous Antibiotics", "Family History Atopy")

df_modelling <- df_analysis_extended
df_modelling$`Age (Months)` <- as.vector(scale(df_modelling$`Age (Months)`))

regression_results <- data.frame()

if(length(target_modules) > 0) {
  for(mod in target_modules) {
    formula_mod <- as.formula(paste0("`", mod, "` ~ Bronchiolitis + `RSV Infection` + Wheezing + `Age (Months)` + Breastfeeding + `Cesarean Delivery` + `Previous Antibiotics` + `Family History Atopy`"))
    fit <- tryCatch(lm(formula_mod, data=df_modelling), error = function(e) NULL)
    
    if(!is.null(fit)) {
      coef_matrix <- summary(fit)$coefficients
      for(pred in predictors_list) {
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
  regression_results <- regression_results %>% mutate(FDR_Value = p.adjust(P_Value, method = "fdr"))
}

regression_results <- regression_results %>%
  mutate(
    Significance_Label = case_when(
      FDR_Value < 0.001 ~ "***",
      FDR_Value < 0.01  ~ "**",
      FDR_Value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Cell_Text = sprintf("%.2f%s", Effect_Size, Significance_Label),
    Predictor = factor(Predictor, levels = predictors_list)
  )

plot_heatmap_expanded <- ggplot(regression_results, aes(x = Predictor, y = Module, fill = Effect_Size)) +
  geom_tile(color = "white", lwd = 0.5) +
  geom_text(aes(label = Cell_Text), color = "black", size = 2.8, fontface = "bold") +
  scale_fill_gradient2(low = "#2c3e50", mid = "#f1c40f", high = "#e74c3c", midpoint = 0, name = "Standardized\nEffect Size (β)") +
  theme_minimal(base_size = 11) +
  labs(title = "Standardized Effect Sizes (β) Across Clinical Traits", x = "Patient Clinical Traits & Contextual Covariates", y = "Structural Network Discovery Modules") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", color = "#2c3e50"), axis.text.y = element_text(face = "bold", color = "#2c3e50"), panel.grid = element_blank())

ggsave(paste0(output_dir_nets, "Plot_Heatmap_Expanded_Clinical_Covariates.png"), plot = plot_heatmap_expanded, width = 11, height = 7, dpi = 300)


# ==============================================================================
# RETORNADO: COMPILACIÓN DE LA TABLA MAESTRA EN FORMATO WIDE
# ==============================================================================
cat("  >> Compilando tabla maestra de coeficientes resumida por Módulo...\n")

deg_npa <- data.frame(Genus = V(g_npa)$name, Niche = "Nasopharynx (NPA)", Degree = degree(g_npa), stringsAsFactors = FALSE)
deg_gut <- data.frame(Genus = V(g_gut)$name, Niche = "Gut (GUT)", Degree = degree(g_gut), stringsAsFactors = FALSE)
deg_master <- rbind(deg_npa, deg_gut)

taxa_npa_status <- data.frame(Genus = V(g_npa)$name, Niche = "Nasopharynx (NPA)", Structural_Module = V(g_npa)$module, stringsAsFactors = FALSE)
taxa_gut_status <- data.frame(Genus = V(g_gut)$name, Niche = "Gut (GUT)", Structural_Module = V(g_gut)$module, stringsAsFactors = FALSE)

master_taxa_raw <- rbind(taxa_npa_status, taxa_gut_status) %>% 
  filter(Structural_Module != 0) %>%
  dplyr::left_join(deg_master, by = c("Genus", "Niche"))

master_modules_collapsed <- master_taxa_raw %>%
  group_by(Niche, Structural_Module) %>%
  summarise(
    Genera_In_Module = paste(Genus, collapse = " "),
    Hub_Genus        = Genus[which.max(Degree)],
    .groups = "drop"
  )

if(exists("regression_results") && nrow(regression_results) > 0) {
  clinical_p_wide <- regression_results %>%
    mutate(String_Report = sprintf("%.3f (FDR: %.4f)", Effect_Size, FDR_Value)) %>%
    select(Module, Predictor, String_Report) %>%
    distinct(Module, Predictor, .keep_all = TRUE) %>%
    mutate(Module_Num = as.numeric(gsub("NPA_Mod|GUT_Mod", "", Module)),
           Niche_Prefix = ifelse(grepl("^NPA", Module), "Nasopharynx (NPA)", "Gut (GUT)")) %>%
    pivot_wider(id_cols = c(Module_Num, Niche_Prefix), names_from = Predictor, values_from = String_Report, names_prefix = "Beta_")
  
  master_clinical_report <- master_modules_collapsed %>%
    dplyr::left_join(clinical_p_wide, by = c("Structural_Module" = "Module_Num", "Niche" = "Niche_Prefix")) %>%
    mutate(Structural_Module = paste0(ifelse(Niche == "Nasopharynx (NPA)", "NPA_Mod", "GUT_Mod"), Structural_Module)) %>%
    arrange(desc(Niche), Structural_Module)
  
  master_clinical_report <- master_clinical_report %>%
    relocate(Genera_In_Module, .after = last_col()) %>%
    relocate(Hub_Genus, .after = last_col())
  
  write.csv(master_clinical_report, paste0(output_dir_nets, "Master_Taxa_Module_Clinical_Associations_FINAL.csv"), row.names = FALSE)
  cat("  >> ¡Tabla maestra de asociación por bacteria exportada con éxito!\n")
}


# ------------------------------------------------------------------------------
# 3.4 COMPARACIÓN TOPOLÓGICA DE REDES Y RIGOR MATEMÁTICO (REVISOR 2)
# ------------------------------------------------------------------------------
cat("\n[3.4] Computando métricas estructurales globales para Tabla Comparativa de Redes...\n")

compute_global_metrics <- function(g_full, label_name) {
  total_nodes         <- vcount(g_full)
  isolated_nodes      <- sum(degree(g_full) == 0)
  connected_nodes     <- total_nodes - isolated_nodes
  pct_nodes_discarded <- (isolated_nodes / total_nodes) * 100
  
  g_sub <- delete_vertices(g_full, V(g_full)[degree(g_full) == 0])
  
  edges_count  <- ecount(g_sub)
  net_density  <- edge_density(g_sub)
  cluster_coef <- transitivity(g_sub, type="global") 
  avg_path_len <- mean_distance(g_sub, directed=FALSE, weights=NA)
  deg_cent     <- centralization.degree(g_sub)$centralization
  
  pos_edges    <- sum(E(g_sub)$weight > 0)
  neg_edges    <- sum(E(g_sub)$weight < 0)
  
  return(data.frame(
    Niche = label_name,
    Total_Input_Genera = total_nodes,
    Connected_Genera = connected_nodes,
    Isolated_Genera_Discarded = isolated_nodes,
    Pct_Genera_Discarded = round(pct_nodes_discarded, 2),
    Total_Edges = edges_count,
    Positive_Alignments = pos_edges,
    Mutual_Exclusions = neg_edges,
    Network_Density = round(net_density, 4),
    Global_Clustering_Coefficient = round(cluster_coef, 4),
    Average_Path_Length = round(avg_path_len, 4),
    Degree_Centralization = round(deg_cent, 4),
    stringsAsFactors = FALSE
  ))
}

network_comparison_matrix <- rbind(
  compute_global_metrics(g_npa, "Nasopharynx (NPA)"),
  compute_global_metrics(g_gut, "Gut (GUT)")
)

write.csv(network_comparison_matrix, paste0(output_dir_nets, "Network_Mathematical_Comparison_Matrix.csv"), row.names = FALSE)
cat("  >> ¡Tabla de comparación topológica global guardada con éxito! (Network_Mathematical_Comparison_Matrix.csv).\n")


# ==============================================================================
# SECCIÓN ADICIONAL: GENERACIÓN DE FIGURA INTEGRADA MÁSTER (Figure_5) - FIJADA
# ==============================================================================
cat("\n[Figuras] Generando la Figura 5 Integrada con maximización real de espacio...\n")

library(cowplot)
library(grid)
library(RColorBrewer)

log_mean_npa <- log1p(colMeans(otu_npa, na.rm = TRUE))
log_mean_gut <- log1p(colMeans(otu_gut, na.rm = TRUE))

# --- PANEL A: Red de Nasofaringe (NPA) (Layout Original + Enlaces Gruesos) ---
g_npa_fig5 <- delete_vertices(g_npa, V(g_npa)[V(g_npa)$module == 0])
V(g_npa_fig5)$label <- gsub("([-_])", "\\1\n", V(g_npa_fig5)$name)
E(g_npa_fig5)$color <- ifelse(E(g_npa_fig5)$weight > 0, "#2ecc71DD", "#e74c3cDD")
E(g_npa_fig5)$width <- abs(E(g_npa_fig5)$weight) * 8.5
V(g_npa_fig5)$size  <- 4 + (log_mean_npa[V(g_npa_fig5)$name] * 1.5)

unique_npa_mods5 <- sort(unique(V(g_npa_fig5)$module))
n_mods <- length(unique_npa_mods5)
colores_dinamicos <- colorRampPalette(brewer.pal(min(8, n_mods), "Set1"))(n_mods)
color_map_npa5    <- setNames(paste0(colores_dinamicos, "E6"), unique_npa_mods5)
V(g_npa_fig5)$color_node <- color_map_npa5[as.character(V(g_npa_fig5)$module)]

p_network_a <- ~{
  layout(matrix(c(1,2), nrow=1), widths=c(3.8, 0.9))
  par(mar=c(0, 0, 0, 0), bg="white")
  set.seed(42)
  
  plot(g_npa_fig5, 
       layout=layout_with_fr(g_npa_fig5, weights=abs(E(g_npa_fig5)$weight)), 
       vertex.color=V(g_npa_fig5)$color_node, 
       vertex.size=V(g_npa_fig5)$size, 
       vertex.frame.color="#1a252f", 
       vertex.frame.width=0.9,
       vertex.label=V(g_npa_fig5)$label, 
       vertex.label.cex=0.85,          
       vertex.label.font=2, 
       vertex.label.color="#2c3e50",
       rescale=TRUE,                 
       xlim=c(-1, 1),            
       ylim=c(-1.0, 1.0))
  
  par(mar=c(1, 0, 1, 1))
  plot.new()
  legend(x = "center", y = 0.85, legend=paste("Mod", names(color_map_npa5)), 
         col=color_map_npa5, pch=19, bty="n", title="Modules", cex=1, title.font=2, y.intersp=1.2)
  legend(x = "top", y = 0.95, legend=c("Alignment (+)", "Exclusion (-)"), 
         col=c("#2ecc71", "#e74c3c"), lty=1, lwd=6, bty="n", title="Interactions", cex=1.1, title.font=2, y.intersp=1.2)
}

# --- PANEL B: Heatmap Clínico ---
p_heatmap_b <- ggplot(regression_results, aes(x = Predictor, y = Module, fill = Effect_Size)) +
  geom_tile(color = "white", lwd = 0.6) +
  geom_text(aes(label = Cell_Text), color = "black", size = 4.0, fontface = "bold") + 
  scale_fill_gradient2(low = "#2c3e50", mid = "#f1c40f", high = "#d35400", midpoint = 0, 
                       name = "Standardized\nEffect Size (β)") +
  theme_minimal(base_size = 14) + 
  labs(x = "Patient Clinical Traits & Contextual Covariates", 
       y = "Structural Network Discovery Modules") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", color = "#1a252f", size = 12),
        axis.text.y = element_text(face = "bold", color = "#1a252f", size = 12), 
        axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
        axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 15)),
        panel.grid = element_blank(),
        plot.title = element_blank(), 
        legend.title = element_text(face = "bold", size = 11),
        legend.position = "right",
        plot.margin = margin(t = 5, r = 5, b = 5, l = 15, unit = "pt"))

# --- COMBINACIÓN Y GUARDADO DE FIGURE_5 MÁSTER ---
figure_5_final <- plot_grid(
  p_network_a, 
  p_heatmap_b, 
  labels = c("a", "b"), 
  label_size = 20, 
  label_fontface = "bold",
  ncol = 1, 
  rel_heights = c(1.3, 1.0)
)

ggsave(
  filename = paste0(output_dir_nets, "Figure_5.png"), 
  plot = figure_5_final, 
  width = 12, 
  height = 16, 
  dpi = 300, 
  bg = "white"
)
cat("  >> [Éxito] Figure_5 guardada maximizando realmente el espacio horizontal.\n")


# ==============================================================================
# GENERACIÓN DE LA SUPPLEMENTARY FIGURE 5: RED DE INTESTINO (GUT)
# ==============================================================================
cat("\n[Figuras] Generando la Supplementary Figure 5 (Red GUT Independiente CORREGIDA)...\n")

g_gut_supp5 <- delete_vertices(g_gut, V(g_gut)[V(g_gut)$module == 0])
V(g_gut_supp5)$label <- gsub("([-_])", "\\1\n", V(g_gut_supp5)$name)
E(g_gut_supp5)$color <- ifelse(E(g_gut_supp5)$weight > 0, "#2ecc71DD", "#e74c3cDD")
E(g_gut_supp5)$width <- abs(E(g_gut_supp5)$weight) * 8.5
V(g_gut_supp5)$size  <- 5 + (log_mean_gut[V(g_gut_supp5)$name] * 1.5)

unique_gut_mods5 <- sort(unique(V(g_gut_supp5)$module))
n_mods_gut5 <- length(unique_gut_mods5)
colores_gut5 <- colorRampPalette(brewer.pal(min(8, n_mods_gut5), "Dark2"))(n_mods_gut5)
color_map_gut5 <- setNames(paste0(colores_gut5, "E6"), unique_gut_mods5)
V(g_gut_supp5)$color_node <- color_map_gut5[as.character(V(g_gut_supp5)$module)]

out_supp5 <- paste0(output_dir_nets, "Supplementary_Figure_5.png")
png(out_supp5, width = 3400, height = 2400, res = 300)

layout(matrix(c(1,2), nrow=1), widths=c(3.9, 0.9))
par(mar=c(1, 1, 1, 1), bg="white")
set.seed(88)

plot(g_gut_supp5, 
     layout=layout_with_fr(g_gut_supp5, weights=abs(E(g_gut_supp5)$weight)), 
     vertex.color=V(g_gut_supp5)$color_node, 
     vertex.size=V(g_gut_supp5)$size, 
     vertex.frame.color="#1a252f", 
     vertex.frame.width=0.9,
     vertex.label=V(g_gut_supp5)$label, 
     vertex.label.cex=0.82,          
     vertex.label.font=2, 
     vertex.label.color="#2c3e50")

par(mar=c(1, 0, 1, 1))
plot.new()
legend(x = "center", y = 0.85, legend=paste("Mod", names(color_map_gut5)), 
       col=color_map_gut5, pch=19, bty="n", title="GUT Modules", cex=1.1, title.font=2, y.intersp=1.3)
legend(x = "bottom", y = 0.15, legend=c("Alignment (+)", "Exclusion (-)"), 
       col=c("#2ecc71", "#e74c3c"), lty=1, lwd=6, bty="n", title="Interactions", cex=1.1, title.font=2, y.intersp=1.3)

dev.off()
cat("  >> [Éxito] Supplementary_Figure_5 guardada correctamente con red simetrizada.\n")
cat("\n[Pipeline Finalizado de forma Exitosa. Estructura ordenada y robusta].\n")