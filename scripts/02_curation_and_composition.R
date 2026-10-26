# ==============================================================================
# 02. DATA CURATION AND TAXONOMIC COMPOSITION (Figure 1)
# ==============================================================================
# - Diagnostic PCoA of all sequenced samples (genus level), kept as a QC output
# - Removal of the samples excluded by the pre-specified QC rule (script 01)
# - Rarefaction to the minimum library size (used only for alpha/beta diversity)
# - Clinical-virological phenotype groups and curated objects for scripts 03-06
# - Figure 1: genus-level composition of NPA and stool, aligned by infant
#
# Outputs
#   results/curated_data/*.rds|csv            curated phyloseq objects and metadata
#   results/02_curation/Diagnostic_PCoA_all_samples.png
#   figures/main/Figure_1.jpg
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})
source("R/functions.R")
set.seed(123)

out_dir <- file.path(DIR_RESULTS, "02_curation")
make_dirs(out_dir, DIR_CURATED)

# ------------------------------------------------------------------------------
# 1. Load data and harmonise sample metadata
# ------------------------------------------------------------------------------
pseq_all <- readRDS("data/phyloseq_object.rds")
sample_data(pseq_all)$Sample_ID <- sample_names(pseq_all)
md <- as(sample_data(pseq_all), "data.frame") %>%
  mutate(Sample_Type = case_when(str_detect(Sample_ID, "ANF$") ~ "NPA",
                                 str_detect(Sample_ID, "GUT$") ~ "GUT",
                                 TRUE ~ Sample_Type),
         Patient_ID = str_remove(Sample_ID, "(ANF|GUT)$"))
sample_data(pseq_all) <- sample_data(md)

# ------------------------------------------------------------------------------
# 2. Diagnostic PCoA of all sequenced samples (genus level, before exclusions)
# ------------------------------------------------------------------------------
pseq_all_genus <- tax_glom(pseq_all, taxrank = "Genus", NArm = FALSE)
pseq_all_rare  <- rarefy_even_depth(pseq_all_genus, sample.size = min(sample_sums(pseq_all_genus)),
                                    rngseed = 123, replace = FALSE, trimOTUs = TRUE)
pcoa_all <- ordinate(pseq_all_rare, method = "PCoA", distance = phyloseq::distance(pseq_all_rare, method = "bray"))
eig <- pcoa_all$values$Relative_eig * 100
df_diag <- data.frame(PCoA1 = pcoa_all$vectors[, 1], PCoA2 = pcoa_all$vectors[, 2],
                      Sample_Type = sample_data(pseq_all_rare)$Sample_Type)
p_diag <- ggplot(df_diag, aes(PCoA1, PCoA2, color = Sample_Type)) +
  geom_point(size = 4.5, alpha = 0.8, stroke = 0.3) +
  stat_ellipse(aes(fill = Sample_Type), geom = "polygon", alpha = 0.10, linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = c(NPA = "#3498db", GUT = "#2ecc71"), name = "Sample Type") +
  scale_fill_manual(values = c(NPA = "#3498db", GUT = "#2ecc71"), name = "Sample Type") +
  labs(title = "Global PCoA of all sequenced samples (genus level)",
       x = paste0("PCoA 1 (", round(eig[1], 1), "%)"), y = paste0("PCoA 2 (", round(eig[2], 1), "%)")) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "Diagnostic_PCoA_all_samples.png"), p_diag, width = 10, height = 8, dpi = 300)

# ------------------------------------------------------------------------------
# 3. Exclusions and rarefaction
# ------------------------------------------------------------------------------
# [Revision 1 - R1.5/a] Exclusions are no longer hard-coded after visual PCoA inspection
# (originally BQLHSO91GUT, HULP71ANF and HULP71GUT); they are read from the
# pre-specified QC rule of script 01 (HULP71ANF and HULP71GUT).
excluded_samples <- readLines(file.path(DIR_RESULTS, "01_sample_qc", "excluded_samples.txt"))
cat("Excluded by QC:", paste(excluded_samples, collapse = ", "), "\n")
pseq <- prune_samples(setdiff(sample_names(pseq_all), excluded_samples), pseq_all)

# Rarefied ASV table (minimum library size), used only for alpha and beta diversity;
# differential abundance and networks use unrarefied counts  [R2.1]
pseq_rarefied <- rarefy_even_depth(pseq, sample.size = min(sample_sums(pseq)), rngseed = 123,
                                   replace = FALSE, trimOTUs = TRUE)
cat("Rarefaction depth:", min(sample_sums(pseq)), "reads\n")
dist_bray_global <- phyloseq::distance(pseq_rarefied, method = "bray")

# ------------------------------------------------------------------------------
# 4. Figure 1: composition (top 8 genera per niche), columns aligned by infant
# ------------------------------------------------------------------------------
pseq_clean <- subset_samples(pseq, Sample_Type %in% c("NPA", "GUT"))
md_clean <- add_clinical_group(as(sample_data(pseq_clean), "data.frame"))
sample_data(pseq_clean) <- sample_data(md_clean)

patients_ordered <- md_clean %>% filter(Sample_Type == "NPA") %>%
  arrange(Condicion_Clinica, Patient_ID) %>% pull(Patient_ID) %>% unique()

df_long <- psmelt(transform_sample_counts(tax_glom(pseq_clean, taxrank = "Genus", NArm = FALSE),
                                          function(x) x / sum(x) * 100))
top_genera <- function(type) df_long %>% filter(Sample_Type == type) %>% group_by(Genus) %>%
  summarise(m = mean(Abundance)) %>% arrange(desc(m)) %>% slice_head(n = 8) %>% pull(Genus)
combined_top <- unique(c(top_genera("NPA"), top_genera("GUT")))
df_long <- df_long %>%
  mutate(Taxon_Final = factor(if_else(Genus %in% combined_top, as.character(Genus), "Other"),
                              levels = c(setdiff(combined_top, "Other"), "Other")),
         Patient_ID = factor(Patient_ID, levels = patients_ordered))

taxa_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF",
                 "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
                 "#1F78B4", "#33A02C", "#FB9A99", "#FDBF6F", "#CAB2D6", "#6A3D9A",
                 "#01665E", "#8C510A", "#053061", "#67001F", "#40004B")
taxa_sin_other <- setdiff(combined_top, "Other")
fill_colors <- setNames(head(taxa_colors, length(taxa_sin_other)), taxa_sin_other)
fill_colors["Other"] <- "#b3b3b3"
wrap_label <- function(x) { x <- as.character(x)
  x[x == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia-\n Caballeronia-\n Paraburkholderia"; x }

composition_panel <- function(type, title) {
  d <- df_long %>% filter(Sample_Type == type) %>%
    group_by(Patient_ID, Condicion_Clinica, Taxon_Final) %>% summarise(Abundance = sum(Abundance), .groups = "drop")
  ggplot(d, aes(x = Patient_ID, y = Abundance, fill = Taxon_Final)) +
    geom_bar(stat = "identity", width = 0.9) +
    scale_fill_manual(values = fill_colors, name = "Microbiota Genus", labels = wrap_label) +
    scale_x_discrete(drop = FALSE) + scale_y_continuous(expand = c(0, 0)) +
    theme_minimal() + labs(title = title, y = "Relative Abundance (%)", x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.title.y = element_text(size = 14, face = "bold"), panel.grid = element_blank(),
          plot.title = element_text(size = 16, face = "bold", margin = margin(b = 10)))
}
clin_colors <- c("CTRL" = "#4DAF4A", "RSV-/Wheeze-" = "#FFFF33", "RSV-/Wheeze+" = "#984EA3",
                 "RSV+/Wheeze-" = "#FF7F00", "RSV+/Wheeze+" = "#377EB8")
p_track <- ggplot(df_long %>% filter(Sample_Type == "NPA") %>% distinct(Patient_ID, Condicion_Clinica),
                  aes(x = Patient_ID, y = 1, fill = Condicion_Clinica)) +
  geom_tile() + scale_fill_manual(values = clin_colors, name = "Clinical Group") +
  scale_x_discrete(drop = FALSE) + scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() + labs(x = "Patients", y = NULL) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
        panel.grid = element_blank())

fig1 <- (composition_panel("NPA", "Nasopharyngeal Aspirate Profile (NPA)") /
           composition_panel("GUT", "Gut Microbiota Profile (GUT)") / p_track) +
  plot_layout(heights = c(10, 10, 1.2), guides = "collect") +
  plot_annotation(tag_levels = list(c("a", "b", ""))) &
  theme(plot.tag = element_text(face = "bold", size = 18),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 13))
save_figure(fig1, 1, width = 20, height = 13)

# ------------------------------------------------------------------------------
# 5. Curated objects for downstream scripts
# ------------------------------------------------------------------------------
md_final <- add_clinical_group(as(sample_data(pseq), "data.frame")) %>%
  mutate(Patient_ID = as.factor(Patient_ID), Sample_Type = as.factor(Sample_Type))
sample_data(pseq)          <- sample_data(md_final)
sample_data(pseq_rarefied) <- sample_data(md_final)

saveRDS(pseq,             file.path(DIR_CURATED, "phyloseq_RAW_curated_global.rds"))
saveRDS(pseq_rarefied,    file.path(DIR_CURATED, "phyloseq_RAREFIED_global.rds"))
saveRDS(dist_bray_global, file.path(DIR_CURATED, "distance_matrix_bray_global.rds"))
saveRDS(md_final,         file.path(DIR_CURATED, "metadata_curated_global.rds"))
write.csv(md_final,       file.path(DIR_CURATED, "metadata_curated_global.csv"), row.names = FALSE)
cat("Curated dataset:", nsamples(pseq), "samples (",
    paste(names(table(md_final$Sample_Type)), table(md_final$Sample_Type), collapse = ", "), ")\n")
