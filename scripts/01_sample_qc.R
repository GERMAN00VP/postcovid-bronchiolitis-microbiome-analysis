# ==============================================================================
# 01. SAMPLE QUALITY CONTROL, SAMPLE FLOW AND RAREFACTION CURVES
# ==============================================================================
# [Revision 1 - R1.5/a, R2.1] New step. In the submitted version samples were
# excluded after visual inspection of the PCoA, without a formal criterion. Every
# sample is now evaluated with three pre-specified quantitative criteria:
#   (1) Sequencing depth >= 10,000 high-quality reads.
#   (2) Sample identity: a sample is flagged as site-discordant when its genus-level
#       profile is closer to samples of the other anatomical site than to samples of
#       its own site (mean Bray-Curtis dissimilarity) AND its 5 nearest neighbours all
#       belong to the other site.
#   (3) Exclusion: a site-discordant sample is excluded only when it is also
#       near-identical (Bray-Curtis < 0.10) to a sample from a DIFFERENT participant,
#       i.e. there is evidence of a sample mix-up. Site discordance alone is not enough,
#       because genuine profiles (e.g. Streptococcus-dominated stool in young infants)
#       can resemble the other site.
# Result: both samples of control HULP71 are excluded (inverted, near-identical to
# samples of other participants); BQLHSO91GUT, excluded in the original submission,
# does not meet the rule and is retained.
#
# Outputs
#   results/01_sample_qc/excluded_samples.txt     list read by script 02
#   results/01_sample_qc/*.csv                    depth, site concordance, sample flow,
#                                                 rarefaction slope (Supplementary Table 5)
#   figures/supplementary/Supplementary_Figure_6.jpg   rarefaction curves  [R2.1]
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
source("R/functions.R")
set.seed(123)

out_dir <- file.path(DIR_RESULTS, "01_sample_qc")
make_dirs(out_dir)

pseq <- readRDS("data/phyloseq_object.rds")
meta <- as(sample_data(pseq), "data.frame")
meta$Sample_ID   <- rownames(meta)
meta$Sample_Type <- ifelse(grepl("ANF$", meta$Sample_ID), "NPA", "GUT")
meta$Patient_ID  <- sub("(ANF|GUT)$", "", meta$Sample_ID)
meta$Group       <- ifelse(meta$Bronchiolitis == "No", "Control", "Bronchiolitis")

# ------------------------------------------------------------------------------
# 1. Sequencing depth (criterion 1)
# ------------------------------------------------------------------------------
depth <- data.frame(Sample_ID = sample_names(pseq), Reads = sample_sums(pseq)) %>%
  left_join(meta[, c("Sample_ID", "Sample_Type", "Group")], by = "Sample_ID") %>%
  mutate(Pass_Depth = Reads >= QC_MIN_DEPTH)
write.csv(depth, file.path(out_dir, "Table_QC_Sequencing_Depth.csv"), row.names = FALSE)
cat("Depth: min =", min(depth$Reads), "| median =", median(depth$Reads),
    "| IQR =", quantile(depth$Reads, .25), "-", quantile(depth$Reads, .75),
    "| samples below", QC_MIN_DEPTH, "=", sum(!depth$Pass_Depth), "\n")

# ------------------------------------------------------------------------------
# 2. Anatomical-site concordance and sample mix-up (criteria 2 and 3)
#    genus level, relative abundances, Bray-Curtis
# ------------------------------------------------------------------------------
pg   <- tax_glom(pseq, taxrank = "Genus", NArm = FALSE)
rel  <- transform_sample_counts(pg, function(x) x / sum(x))
otu  <- as(otu_table(rel), "matrix"); if (taxa_are_rows(rel)) otu <- t(otu)
D    <- as.matrix(vegdist(otu, method = "bray"))
site <- setNames(meta$Sample_Type, meta$Sample_ID)[rownames(D)]
gen  <- as.character(tax_table(rel)[, "Genus"])

conc <- bind_rows(lapply(rownames(D), function(s) {
  d_other <- D[s, setdiff(names(site)[site != site[s]], s)]
  d_own   <- D[s, setdiff(names(site)[site == site[s]], s)]
  nn      <- names(sort(D[s, colnames(D) != s]))[1:5]
  pair    <- paste0(meta$Patient_ID[meta$Sample_ID == s], ifelse(site[s] == "NPA", "GUT", "ANF"))
  closest <- names(sort(D[s, colnames(D) != s]))[1]
  top     <- which.max(otu[s, ])
  data.frame(Sample_ID = s, Site = site[s],
             Mean_BC_own_site = round(mean(d_own), 3),
             Mean_BC_other_site = round(mean(d_other), 3),
             NN5_other_site = sum(site[nn] != site[s]),
             Closest_sample = closest, BC_closest = round(D[s, closest], 3),
             BC_to_own_paired_sample = if (pair %in% colnames(D)) round(D[s, pair], 3) else NA,
             Dominant_genus = gen[top], Dominant_rel_abund = round(otu[s, top], 3),
             N_genera = sum(otu[s, ] > 0))
}))
conc <- conc %>%
  mutate(Site_Discordant = Mean_BC_other_site < Mean_BC_own_site & NN5_other_site == 5,
         Closest_other_participant = sub("(ANF|GUT)$", "", Closest_sample) != sub("(ANF|GUT)$", "", Sample_ID),
         Mixup = Site_Discordant & BC_closest < QC_MIXUP_BC_MAX & Closest_other_participant) %>%
  left_join(depth[, c("Sample_ID", "Pass_Depth")], by = "Sample_ID") %>%
  mutate(Excluded = Mixup | !Pass_Depth) %>%
  select(-Mixup) %>%
  left_join(meta[, c("Sample_ID", "Group", "Age")], by = "Sample_ID")

EXCLUDED_SAMPLES <- conc$Sample_ID[conc$Excluded]
writeLines(EXCLUDED_SAMPLES, file.path(out_dir, "excluded_samples.txt"))
write.csv(conc, file.path(out_dir, "Table_QC_Site_Concordance.csv"), row.names = FALSE)
cat("\nSite-discordant samples:\n")
print(conc[conc$Site_Discordant, c("Sample_ID", "Group", "Age", "Mean_BC_own_site", "Mean_BC_other_site",
                                   "Closest_sample", "BC_closest", "Dominant_genus", "Dominant_rel_abund", "Excluded")],
      row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Sample flow (collected / sequenced / excluded / analysed)
# ------------------------------------------------------------------------------
flow <- meta %>%
  mutate(Analysed = !(Sample_ID %in% EXCLUDED_SAMPLES)) %>%
  group_by(Group, Sample_Type) %>%
  summarise(Sequenced = n(), Excluded = sum(!Analysed), Analysed = sum(Analysed), .groups = "drop") %>%
  left_join(meta %>% distinct(Patient_ID, Group) %>% count(Group, name = "Participants"), by = "Group") %>%
  mutate(Not_collected = Participants - Sequenced)
write.csv(flow, file.path(out_dir, "Table_Sample_Flow.csv"), row.names = FALSE)
cat("\nSample flow:\n"); print(as.data.frame(flow))

# ------------------------------------------------------------------------------
# 4. Rarefaction curves at ASV level  [R2.1] -> Supplementary Figure 6
# ------------------------------------------------------------------------------
asv   <- as(otu_table(pseq), "matrix"); if (taxa_are_rows(pseq)) asv <- t(asv)
steps <- c(seq(100, 1000, by = 300), seq(2000, 80000, by = 2000))
rc <- bind_rows(lapply(rownames(asv), function(s) {
  n <- sum(asv[s, ]); st <- c(steps[steps < n], n)
  data.frame(Sample_ID = s, Reads = st, ASVs = as.numeric(rarefy(asv[s, , drop = FALSE], st)))
})) %>%
  left_join(meta[, c("Sample_ID", "Sample_Type")], by = "Sample_ID") %>%
  mutate(Sample_Type = factor(Sample_Type, levels = c("NPA", "GUT"),
                              labels = c("Nasopharyngeal aspirate (NPA)", "Gut (stool)")))

p_rc <- ggplot(rc, aes(Reads, ASVs, group = Sample_ID)) +
  geom_line(alpha = 0.45, linewidth = 0.4, colour = "#2c3e50") +
  geom_vline(xintercept = min(depth$Reads), linetype = 2, colour = "#c0392b") +
  facet_wrap(~ Sample_Type, scales = "free_y") +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Sequencing depth (reads)", y = "Observed ASVs",
       caption = paste0("Dashed line: rarefaction depth used for diversity analyses (",
                        format(min(depth$Reads), big.mark = ","), " reads)")) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#2c3e50"),
        strip.text = element_text(colour = "white", face = "bold"))
save_figure(p_rc, 6, width = 10, height = 4.8, supplementary = TRUE)

# Plateau metric: new ASVs per 1,000 additional reads at the rarefaction depth
slope <- rc %>% group_by(Sample_ID, Sample_Type) %>% arrange(Reads) %>%
  summarise(slope_per_1000 = {
    i <- which.min(abs(Reads - min(depth$Reads))); j <- min(i + 1, n())
    if (j == i) NA else 1000 * (ASVs[j] - ASVs[i]) / (Reads[j] - Reads[i])
  }, .groups = "drop")
write.csv(slope, file.path(out_dir, "Table_QC_Rarefaction_Slope.csv"), row.names = FALSE)
cat("\nNew ASVs per 1,000 additional reads at the rarefaction depth:\n")
print(slope %>% group_by(Sample_Type) %>%
        summarise(median = median(slope_per_1000, na.rm = TRUE),
                  q1 = quantile(slope_per_1000, .25, na.rm = TRUE),
                  q3 = quantile(slope_per_1000, .75, na.rm = TRUE)))
