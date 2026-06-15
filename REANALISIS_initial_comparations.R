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
library(forcats)
library(patchwork)

set.seed(123)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4)

setwd("/mnt/usb/BQL/BQL_ANALYSIS/COMPLETE_REANALISIS_BQL/")

output_dir_global <- "REANALISIS/Global/"

if(!dir.exists(output_dir_global)) dir.create(output_dir_global, recursive = TRUE)

# ==============================================================================
# SECTION 2: GLOBAL DATA PRE-PROCESSING & PCoA (BEFORE OUTLIER REMOVAL, AT GENUS LEVEL)
# ==============================================================================
pseq_raw_init <- readRDS("data/phyloseq_object.rds")

# A. Sincronización y curación de metadatos iniciales con TODO (incluyendo outliers)
sample_data(pseq_raw_init)$Sample_ID <- sample_names(pseq_raw_init)
metadata_init <- as(sample_data(pseq_raw_init), "data.frame") %>%
  mutate(
    Sample_Type = case_when(
      str_detect(Sample_ID, "ANF$") ~ "NPA", 
      str_detect(Sample_ID, "GUT$") ~ "GUT",
      TRUE ~ Sample_Type
    ),
    Patient_ID = str_remove(Sample_ID, "(ANF|GUT)$")
  )
sample_data(pseq_raw_init) <- sample_data(metadata_init)

# B. Transformación a nivel de GÉNERO y Rarefacción manteniendo los outliers
pseq_init_genus <- tax_glom(pseq_raw_init, taxrank = "Genus", NArm = FALSE)
min_depth_init  <- min(sample_sums(pseq_init_genus))
pseq_init_rare  <- rarefy_even_depth(pseq_init_genus, sample.size = min_depth_init, rngseed = 123, replace = FALSE, trimOTUs = TRUE)

# C. Cálculo de distancias y ordenación (Nivel Género con Outliers)
dist_bray_init <- phyloseq::distance(pseq_init_rare, method = "bray")
pcoa_init      <- ordinate(pseq_init_rare, method = "PCoA", distance = dist_bray_init)
eigenvalues    <- pcoa_init$values$Relative_eig * 100

df_pcoa_diag <- data.frame(
  PCoA1 = pcoa_init$vectors[, 1],
  PCoA2 = pcoa_init$vectors[, 2],
  Sample_Type = sample_data(pseq_init_rare)$Sample_Type,
  Patient_ID = sample_data(pseq_init_rare)$Patient_ID
)

# D. Renderizado del PCoA global original a nivel de Género
plot_pcoa_diag <- ggplot(df_pcoa_diag, aes(x = PCoA1, y = PCoA2, color = Sample_Type)) +
  geom_point(size = 4.5, alpha = 0.80, stroke = 0.3) + 
  stat_ellipse(aes(fill = Sample_Type), geom = "polygon", alpha = 0.10, linewidth = 0.8, linetype = "dashed") +
  theme_bw() +
  scale_color_manual(values = c("NPA" = "#3498db", "GUT" = "#2ecc71"), name = "Sample Type") +
  scale_fill_manual(values = c("NPA" = "#3498db", "GUT" = "#2ecc71"), name = "Sample Type") +
  labs(
    title = "Global PCoA (Genus Level - Includes Outliers)",
    subtitle = "",
    x = paste0("PCoA 1 (", round(eigenvalues[1], 1), "%)"),
    y = paste0("PCoA 2 (", round(eigenvalues[2], 1), "%)")
  ) +
  theme(
    text = element_text(size = 13, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, face = "italic", hjust = 0.5, margin = margin(b = 15)),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right"
  )

ggsave(paste0(output_dir_global, "Diagnostic_PCoA_Bray_Curtis_FINAL.png"), plot = plot_pcoa_diag, width = 10, height = 8, dpi = 300)


# ==============================================================================
# ELIMINACIÓN POSTERIOR DE LOS 3 OUTLIERS (Para el resto del flujo estándar)
# ==============================================================================
pseq <- pseq_raw_init # Hereda los metadatos corregidos estructurados arriba
muestras_a_mantener <- setdiff(sample_names(pseq), c("BQLHSO91GUT", "HULP71ANF","HULP71GUT"))
pseq <- prune_samples(muestras_a_mantener, pseq)

# Recalculamos la rarefacción limpia para el entorno río abajo (sin los 3 outliers)
min_seq_depth <- min(sample_sums(pseq))
pseq_rarefied <- rarefy_even_depth(pseq, sample.size = min_seq_depth, rngseed = 123, replace = FALSE, trimOTUs = TRUE)
dist_bray_global <- phyloseq::distance(pseq_rarefied, method = "bray")

# ==============================================================================
# SECTION 3: ALIGNED DUAL BARPLOT WITH HIGH-CONTRAST VIBRANT PALETTE
# ==============================================================================
pseq_clean <- subset_samples(pseq, Sample_Type %in% c("NPA", "GUT"))

metadata_clean <- as(sample_data(pseq_clean), "data.frame") %>%
  mutate(
    RSV_Wheezing = paste0(Respiratory.syncytial.virus, "_", Wheezing.treatment),
    RSV_Wheezing = if_else(Bronchiolitis == "No", "CTRL", RSV_Wheezing),
    Condicion_Clinica = case_when(
      RSV_Wheezing == "CTRL"        ~ "CTRL",
      RSV_Wheezing == "No_No"       ~ "RSV-/Wheeze-",
      RSV_Wheezing == "No_Yes"      ~ "RSV-/Wheeze+",
      RSV_Wheezing == "Yes_No"      ~ "RSV+/Wheeze-",
      RSV_Wheezing == "Yes_Yes"     ~ "RSV+/Wheeze+",
      TRUE                          ~ RSV_Wheezing
    )
  ) %>%
  mutate(Condicion_Clinica = factor(Condicion_Clinica, 
                                    levels = c("CTRL", "RSV-/Wheeze-", "RSV-/Wheeze+", "RSV+/Wheeze-", "RSV+/Wheeze+")))

sample_data(pseq_clean) <- sample_data(metadata_clean)

all_patients_ordered <- metadata_clean %>%
  filter(Sample_Type == "NPA") %>% 
  arrange(Condicion_Clinica, Patient_ID) %>%
  pull(Patient_ID) %>%
  unique()

pseq_genus <- tax_glom(pseq_clean, taxrank = "Genus", NArm = FALSE)
pseq_rel <- transform_sample_counts(pseq_genus, function(x) (x / sum(x)) * 100)
df_long_all <- psmelt(pseq_rel)

# --- TOP 8 POR ÓRGANO ---
top_npa <- df_long_all %>% filter(Sample_Type == "NPA") %>%
  group_by(Genus) %>% summarise(m = mean(Abundance)) %>% arrange(desc(m)) %>% slice_head(n = 8) %>% pull(Genus)

top_gut <- df_long_all %>% filter(Sample_Type == "GUT") %>%
  group_by(Genus) %>% summarise(m = mean(Abundance)) %>% arrange(desc(m)) %>% slice_head(n = 8) %>% pull(Genus)

combined_top_taxa <- unique(c(top_npa, top_gut))

df_long_all <- df_long_all %>%
  mutate(Taxon_Final = if_else(Genus %in% combined_top_taxa, as.character(Genus), "Other"))

df_long_all$Patient_ID <- factor(df_long_all$Patient_ID, levels = all_patients_ordered)

taxa_levels <- c(setdiff(combined_top_taxa, "Other"), "Other")
df_long_all$Taxon_Final <- factor(df_long_all$Taxon_Final, levels = taxa_levels)

# --- PALETA CHILLONA / ALTO CONTRASTE ---
taxa_sin_other <- setdiff(combined_top_taxa, "Other")
color_count <- length(taxa_sin_other)

colores_vivos <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", 
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",             
  "#1F78B4", "#33A02C", "#FB9A99", "#FDBF6F", "#CAB2D6", "#6A3D9A",             
  "#01665E", "#8C510A", "#053061", "#67001F", "#40004B"
)

unique_colors <- head(colores_vivos, color_count)
names(unique_colors) <- taxa_sin_other
unique_colors["Other"] <- "#b3b3b3" 

# --- PANEL 1: NPA PLOT ---
df_npa_plot <- df_long_all %>% filter(Sample_Type == "NPA") %>%
  group_by(Patient_ID, Condicion_Clinica, Taxon_Final) %>% summarise(Abundance = sum(Abundance), .groups = 'drop')

p_npa <- ggplot(df_npa_plot, aes(x = Patient_ID, y = Abundance, fill = Taxon_Final)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_manual(values = unique_colors, name = "Microbiota Genus") +
  scale_x_discrete(drop = FALSE) + 
  theme_minimal() +
  labs(title = "Nasopharyngeal Aspirate Profile (NPA)", y = "Relative Abundance (%)", x = NULL) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),                  # Fuente de eje y más grande
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 10))
  ) +
  scale_y_continuous(expand = c(0, 0))

# --- PANEL 2: GUT PLOT ---
df_gut_plot <- df_long_all %>% filter(Sample_Type == "GUT") %>%
  group_by(Patient_ID, Condicion_Clinica, Taxon_Final) %>% summarise(Abundance = sum(Abundance), .groups = 'drop')

p_gut <- ggplot(df_gut_plot, aes(x = Patient_ID, y = Abundance, fill = Taxon_Final)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_manual(values = unique_colors, name = "Microbiota Genus") +
  scale_x_discrete(drop = FALSE) + 
  theme_minimal() +
  labs(title = "Gut Microbiota Profile (GUT)", y = "Relative Abundance (%)", x = NULL) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),                  # Fuente de eje y más grande
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 10))
  ) +
  scale_y_continuous(expand = c(0, 0))

# --- PANEL 3: CLINICAL METADATA TRACK ---
df_meta_bar <- df_long_all %>% 
  filter(Sample_Type == "NPA") %>% 
  select(Patient_ID, Condicion_Clinica) %>% 
  distinct()

colors_clinicos <- c("CTRL" = "#4DAF4A", "RSV-/Wheeze-" = "#FFFF33", "RSV-/Wheeze+" = "#984EA3", "RSV+/Wheeze-" = "#FF7F00", "RSV+/Wheeze+" = "#377EB8")

p_meta_bar <- ggplot(df_meta_bar, aes(x = Patient_ID, y = 1, fill = Condicion_Clinica)) +
  geom_tile() +
  scale_fill_manual(values = colors_clinicos, name = "Clinical Group") +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  labs(x = "Patients", y = NULL) +
  theme(
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
    panel.grid = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),                   # Títulos de leyendas más legibles
    legend.text = element_text(size = 13)                                    # Textos de leyendas más legibles
  ) +
  scale_y_continuous(expand = c(0, 0))

# --- ASSEMBLING WITH PATCHWORK ---
# Se mantiene tu estructura de layouts y colecta de guías tal como estaba originalmente
plot_mirror_final <- (p_npa / p_gut / p_meta_bar) + 
  plot_layout(heights = c(10, 10, 1.2), guides = "collect")

# Asignamos tags 'a' y 'b' solo a los dos paneles principales, dejando el track clínico sin letra
plot_mirror_final <- plot_mirror_final + 
  plot_annotation(tag_levels = list(c('a', 'b', ''))) & 
  theme(plot.tag = element_text(face = "bold", size = 18),
        legend.title = element_text(size = 14, face = "bold"),                # Asegura legibilidad en la leyenda colectada
        legend.text = element_text(size = 13))

ggsave(
  filename = paste0(output_dir_global, "Aligned_NPA_vs_GUT_Curated_Plot.png"),
  plot = plot_mirror_final, 
  width = 20, 
  height = 13, 
  dpi = 300
)

# ==============================================================================
# SECTION 4: EXPORT CURATED GLOBAL ENVIRONMENT FOR DOWNSTREAM ANALYSES
# ==============================================================================
data_output_dir <- paste0(output_dir_global, "curated_data_global/")
if(!dir.exists(data_output_dir)) dir.create(data_output_dir, recursive = TRUE)

cat("Construyendo metadatos clínicos finales e integrándolos en el objeto global...\n")

metadata_global_final <- as(sample_data(pseq), "data.frame") %>%
  mutate(
    RSV_Wheezing = paste0(Respiratory.syncytial.virus, "_", Wheezing.treatment),
    RSV_Wheezing = if_else(Bronchiolitis == "No", "CTRL", RSV_Wheezing),
    Condicion_Clinica = case_when(
      RSV_Wheezing == "CTRL"        ~ "CTRL",
      RSV_Wheezing == "No_No"       ~ "RSV-/Wheeze-",
      RSV_Wheezing == "No_Yes"      ~ "RSV-/Wheeze+",
      RSV_Wheezing == "Yes_No"      ~ "RSV+/Wheeze-",
      RSV_Wheezing == "Yes_Yes"     ~ "RSV+/Wheeze+",
      TRUE                          ~ RSV_Wheezing
    )
  ) %>%
  mutate(
    Condicion_Clinica = factor(Condicion_Clinica, 
                               levels = c("CTRL", "RSV-/Wheeze-", "RSV-/Wheeze+", "RSV+/Wheeze-", "RSV+/Wheeze+")),
    Patient_ID = as.factor(Patient_ID),
    Sample_Type = as.factor(Sample_Type)
  )

# Sincronizamos los metadatos completos (SIN OUTLIERS) dentro de ambos objetos antes de guardar
sample_data(pseq)          <- sample_data(metadata_global_final)
sample_data(pseq_rarefied) <- sample_data(metadata_global_final)

# Guardado de objetos limpios en disco (exactamente igual que antes)
saveRDS(pseq, file = paste0(data_output_dir, "phyloseq_RAW_curated_global.rds"))
saveRDS(pseq_rarefied, file = paste0(data_output_dir, "phyloseq_RAREFIED_global.rds"))
saveRDS(dist_bray_global, file = paste0(data_output_dir, "distance_matrix_bray_global.rds"))
saveRDS(metadata_global_final, file = paste0(data_output_dir, "metadata_curated_global.rds"))
write.csv(metadata_global_final, file = paste0(data_output_dir, "metadata_curated_global.csv"), row.names = FALSE)

cat("\n[✓] ¡PROCESO DE CURACIÓN Y ENRIQUECIMIENTO COMPLETADO!")