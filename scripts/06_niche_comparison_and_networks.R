# ==============================================================================
# 06. NPA vs. GUT COMPARISON AND CO-ABUNDANCE NETWORKS
# ==============================================================================
#   Paired comparison of the two niches (LME, blocked PERMANOVA) -> Supplementary Figure 1
#   SpiecEasi (MB) networks per niche, Louvain modules, module eigen-scores and
#   their association with clinical traits                      -> Figure 6,
#                                                                  Supplementary Figure 5,
#                                                                  Supplementary Table 3
# Infants with both samples only. Genus level, rarefied counts, prevalence >= 10%.
#
# [Revision 1 - R1.5/a] With the final sample set (BQLHSO91GUT reinstated, HULP71
# excluded) rarefaction changes for all samples, so the networks were rebuilt and
# their modules renumbered relative to the submission; the main findings are unchanged.
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
  library(RColorBrewer)
  library(patchwork)
  library(ggsignif)
  library(cowplot)
})
source("R/functions.R")
set.seed(123)

out_dir <- file.path(DIR_RESULTS, "06_niches_and_networks")
make_dirs(out_dir)

niche_colors <- c(NPA = "#2980b9", GUT = "#27ae60")

# ------------------------------------------------------------------------------
# 1. Data: genus level, infants with both samples
# ------------------------------------------------------------------------------
pseq_rare <- readRDS(file.path(DIR_CURATED, "phyloseq_RAREFIED_global.rds"))
pseq_rare_genus <- tax_glom(pseq_rare, taxrank = "Genus", NArm = FALSE) %>% rename_to_genus("Genus")
md <- extract_meta(pseq_rare_genus) %>%
  mutate(Condicion_Clinica = factor(Condicion_Clinica, levels = GROUP_LEVELS),
         Patient_ID = as.factor(Patient_ID), Sample_Type = factor(Sample_Type, levels = c("NPA", "GUT")))
both <- md %>% group_by(Patient_ID) %>% filter(all(c("NPA", "GUT") %in% Sample_Type)) %>% pull(Patient_ID) %>% unique()
pseq_rare_genus <- subset_samples(pseq_rare_genus, Patient_ID %in% both)
metadata_df <- extract_meta(pseq_rare_genus)
sample_data(pseq_rare_genus) <- sample_data(metadata_df)

# ------------------------------------------------------------------------------
# 2. Niche comparison: alpha (LME with random intercept per infant) and beta
#    (PERMANOVA with permutations blocked by infant)
# ------------------------------------------------------------------------------
alpha <- estimate_richness(pseq_rare_genus, measures = c("Observed", "Shannon", "Simpson"))
df_alpha <- merge(data.frame(Sample_ID = rownames(alpha), alpha[, c("Observed", "Shannon", "Simpson")]),
                  metadata_df, by = "Sample_ID")
signif <- list()
for (metric in c("Observed", "Shannon", "Simpson")) {
  fit <- lme(as.formula(paste(metric, "~ Sample_Type + Age")), random = ~ 1 | Patient_ID, data = df_alpha, method = "REML")
  tt  <- summary(fit)$tTable
  write.csv(tt, file.path(out_dir, paste0("LME_Alpha_", metric, ".csv")))
  df_alpha[[paste0(metric, "_Corrected")]] <- df_alpha[[metric]] - tt["Age", "Value"] * df_alpha$Age
  p <- as.data.frame(emmeans(fit, specs = pairwise ~ Sample_Type, adjust = "fdr")$contrasts)$p.value
  signif[[metric]] <- data.frame(p = p, annot = ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

bray_dist <- phyloseq::distance(pseq_rare_genus, method = "bray")
ord_pcoa  <- ordinate(pseq_rare_genus, method = "PCoA", distance = "bray")
df_pcoa   <- plot_ordination(pseq_rare_genus, ord_pcoa, justDF = TRUE)
perm_blocks <- how(blocks = metadata_df$Patient_ID, nperm = 999)
write.csv(as.data.frame(adonis2(bray_dist ~ Age + Sample_Type, data = metadata_df, permutations = perm_blocks, by = "terms")),
          file.path(out_dir, "PERMANOVA_NPA_vs_GUT_paired.csv"))
write.csv(as.data.frame(permutest(betadisper(bray_dist, group = metadata_df$Sample_Type), permutations = perm_blocks)$tab),
          file.path(out_dir, "Betadisper_NPA_vs_GUT.csv"))

alpha_box <- function(y_col, label_y, sg) {
  p <- ggplot(df_alpha, aes(x = Sample_Type, y = .data[[y_col]], fill = Sample_Type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.5, color = "#2c3e50") +
    geom_jitter(width = 0.1, alpha = 0.4, size = 1.5, color = "#2c3e50") +
    scale_fill_manual(values = niche_colors) + theme_bw(base_size = 11) + labs(x = NULL, y = label_y) +
    theme(axis.title.y = element_text(face = "bold", size = 10), axis.text = element_text(color = "black", size = 9),
          axis.text.x = element_text(face = "bold"), legend.position = "none",
          panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#f0f0f0"))
  if (sg$annot != "ns") {
    mx <- max(df_alpha[[y_col]], na.rm = TRUE); rg <- mx - min(df_alpha[[y_col]], na.rm = TRUE)
    p <- p + geom_signif(comparisons = list(c("NPA", "GUT")), annotations = sg$annot, y_position = mx + rg * 0.05,
                         tip_length = 0.02, vjust = 0.5, textsize = 4, color = "black")
  }
  p
}
p_beta <- ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, color = Sample_Type)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(fill = Sample_Type), geom = "polygon", alpha = 0.1, type = "t", linetype = 2) +
  scale_color_manual(values = niche_colors, name = "Niche Location") +
  scale_fill_manual(values = niche_colors, name = "Niche Location") + theme_bw(base_size = 11) +
  labs(x = paste0("PCoA 1 (", round(ord_pcoa$values$Relative_eig[1] * 100, 1), "%)"),
       y = paste0("PCoA 2 (", round(ord_pcoa$values$Relative_eig[2] * 100, 1), "%)")) +
  theme(axis.title = element_text(face = "bold", size = 10), axis.text = element_text(color = "black", size = 9),
        legend.title = element_text(face = "bold", size = 10), legend.text = element_text(size = 9),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#f0f0f0"))
supp1 <- (alpha_box("Observed_Corrected", "Observed Richness\n(Age-Corrected)", signif$Observed) +
            alpha_box("Shannon_Corrected", "Shannon Index\n(Age-Corrected)", signif$Shannon) +
            alpha_box("Simpson_Corrected", "Simpson Index\n(Age-Corrected)", signif$Simpson) + plot_layout(ncol = 3)) /
  p_beta + plot_layout(heights = c(0.7, 1)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14, color = "#2c3e50"))
save_figure(supp1, 1, width = 11, height = 8.5, supplementary = TRUE)

# ------------------------------------------------------------------------------
# 3. Co-abundance networks per niche (SpiecEasi, Meinshausen-Buhlmann)
# ------------------------------------------------------------------------------
pseq_filt <- filter_taxa(pseq_rare_genus, function(x) sum(x > 0) >= 0.10 * nsamples(pseq_rare_genus), prune = TRUE)
taxa_names(pseq_filt) <- clean_taxa_names(taxa_names(pseq_filt))
remove_zero_var <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  prune_taxa(if (taxa_are_rows(ps)) apply(m, 1, var) > 0 else apply(m, 2, var) > 0, ps)
}
otu_by_patient <- function(ps) {
  m <- otu_table(ps); if (taxa_are_rows(ps)) m <- t(m)
  m <- as.matrix(m); rownames(m) <- as.character(sample_data(ps)$Patient_ID); m
}
otu_npa <- otu_by_patient(remove_zero_var(subset_samples(pseq_filt, Sample_Type == "NPA")))
otu_gut <- otu_by_patient(remove_zero_var(subset_samples(pseq_filt, Sample_Type == "GUT")))
common  <- intersect(rownames(otu_npa), rownames(otu_gut))
otu_npa <- otu_npa[sort(common), , drop = FALSE]
otu_gut <- otu_gut[sort(common), , drop = FALSE]

run_spiec <- function(m) spiec.easi(m, method = "mb", lambda.min.ratio = 1e-2, nlambda = 20,
                                    pulsar.params = list(rep.num = 50, ncores = 1))
se_npa <- run_spiec(otu_npa)
se_gut <- run_spiec(otu_gut)

to_igraph <- function(se_obj, otu_mat) {
  b <- as.matrix(symBeta(getOptBeta(se_obj), mode = "maxabs"))
  colnames(b) <- rownames(b) <- colnames(otu_mat)
  g <- graph_from_adjacency_matrix(b, mode = "undirected", weighted = TRUE, diag = FALSE)
  V(g)$name <- V(g)$label <- colnames(otu_mat)
  g
}
add_modules <- function(g, label) {   # Louvain on connected genera; isolated genera = module 0
  g_sub <- delete_vertices(g, degree(g) == 0)
  comm  <- cluster_louvain(g_sub, weights = abs(E(g_sub)$weight))
  V(g)$module <- setNames(membership(comm), V(g_sub)$name)[V(g)$name]
  V(g)$module[is.na(V(g)$module)] <- 0
  cat(sprintf("  [%s] %d Louvain modules\n", label, max(comm$membership)))
  g
}
add_hubs <- function(g, top_pct = 0.05) {
  dg <- degree(g); V(g)$is_hub <- dg >= quantile(dg, probs = 1 - top_pct); g
}
g_npa <- add_hubs(add_modules(to_igraph(se_npa, otu_npa), "NPA"))
g_gut <- add_hubs(add_modules(to_igraph(se_gut, otu_gut), "GUT"))

write.csv(bind_rows(
  data.frame(Network = "NPA", Genus = V(g_npa)$name[V(g_npa)$is_hub], Degree = degree(g_npa)[V(g_npa)$is_hub], Module = V(g_npa)$module[V(g_npa)$is_hub]),
  data.frame(Network = "GUT", Genus = V(g_gut)$name[V(g_gut)$is_hub], Degree = degree(g_gut)[V(g_gut)$is_hub], Module = V(g_gut)$module[V(g_gut)$is_hub])) %>%
    arrange(Network, desc(Degree)), file.path(out_dir, "Hub_Genera.csv"), row.names = FALSE)

# [Revision 1] Final networks exported (igraph objects and signed edge lists) for reproducibility
for (nm in c("NPA", "GUT")) {
  g <- if (nm == "NPA") g_npa else g_gut
  saveRDS(g, file.path(out_dir, paste0("igraph_Genus_", nm, ".rds")))
  el <- as.data.frame(as_edgelist(g)); colnames(el) <- c("Genus_1", "Genus_2")
  el$Weight   <- E(g)$weight
  el$Module_1 <- V(g)$module[match(el$Genus_1, V(g)$name)]
  el$Module_2 <- V(g)$module[match(el$Genus_2, V(g)$name)]
  write.csv(el, file.path(out_dir, paste0("Edge_List_", nm, ".csv")), row.names = FALSE)
}

# ------------------------------------------------------------------------------
# 4. Module eigen-scores (first eigenvector of the module adjacency, oriented to the hub,
#    applied to log1p counts and z-scored) and multivariable association with clinical traits
# ------------------------------------------------------------------------------
module_scores <- function(g, otu_mat, prefix) {
  adj <- as.matrix(as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  lg  <- log1p(otu_mat); out <- list()
  for (m in setdiff(sort(unique(V(g)$module)), 0)) {
    gen <- intersect(V(g)$name[V(g)$module == m], colnames(otu_mat))
    if (length(gen) > 1) {
      ev <- eigen(adj[gen, gen, drop = FALSE])$vectors[, 1]; names(ev) <- gen
      if (ev[gen[which.max(degree(g, v = gen))]] < 0) ev <- -ev
      out[[paste0(prefix, m)]] <- as.vector(scale(as.matrix(lg[, gen]) %*% ev))
    } else if (length(gen) == 1) out[[paste0(prefix, m)]] <- as.vector(scale(lg[, gen]))
  }
  as.data.frame(out, check.names = FALSE)
}
df_scores <- cbind(Patient_ID = rownames(otu_npa), module_scores(g_npa, otu_npa, "NPA_Mod"),
                   module_scores(g_gut, otu_gut, "GUT_Mod"))

yes <- function(x) factor(ifelse(grepl("^Yes$|1", x, ignore.case = TRUE), "Yes", "No"), levels = c("No", "Yes"))
clin <- extract_meta(pseq_rare_genus) %>%
  mutate(`Bronchiolitis` = yes(Bronchiolitis),
         `RSV Infection` = factor(ifelse(grepl("^Yes$|1|RSV\\+", Respiratory.syncytial.virus) | grepl("RSV\\+", Condicion_Clinica), "Yes", "No"), levels = c("No", "Yes")),
         `Wheezing` = factor(ifelse(grepl("Wheeze\\+", Condicion_Clinica) | grepl("^Yes$|1", Wheezing.treatment, ignore.case = TRUE) |
                                      (suppressWarnings(as.numeric(Wheezing.count)) > 0), "Yes", "No"), levels = c("No", "Yes")),
         `Cesarean Delivery` = yes(Cesarean.section), `Previous Antibiotics` = yes(Previous.antibiotics),
         `Family History Atopy` = yes(Family.history.atopy), `Breastfeeding` = yes(Breastfeeding),
         `Age (Months)` = as.numeric(Age)) %>%
  select(Patient_ID, `Bronchiolitis`, `RSV Infection`, `Wheezing`, `Age (Months)`, `Cesarean Delivery`,
         `Previous Antibiotics`, `Family History Atopy`, `Breastfeeding`) %>%
  distinct(Patient_ID, .keep_all = TRUE)
df_mod <- inner_join(df_scores, clin, by = "Patient_ID")
df_mod$`Age (Months)` <- as.vector(scale(df_mod$`Age (Months)`))

predictors <- c("Bronchiolitis", "RSV Infection", "Wheezing", "Age (Months)", "Breastfeeding",
                "Cesarean Delivery", "Previous Antibiotics", "Family History Atopy")
modules <- c(grep("^NPA_Mod", colnames(df_mod), value = TRUE), grep("^GUT_Mod", colnames(df_mod), value = TRUE))
reg <- bind_rows(lapply(modules, function(mod) {
  fit <- lm(as.formula(paste0("`", mod, "` ~ Bronchiolitis + `RSV Infection` + Wheezing + `Age (Months)` + Breastfeeding + ",
                              "`Cesarean Delivery` + `Previous Antibiotics` + `Family History Atopy`")), data = df_mod)
  cm <- summary(fit)$coefficients
  bind_rows(lapply(predictors, function(p) {
    r <- rownames(cm)[grepl(p, rownames(cm), fixed = TRUE)][1]
    data.frame(Module = mod, Predictor = p, Effect_Size = cm[r, "Estimate"], P_Value = cm[r, "Pr(>|t|)"])
  }))
})) %>%
  mutate(FDR_Value = p.adjust(P_Value, method = "fdr"),
         Significance_Label = case_when(FDR_Value < 0.001 ~ "***", FDR_Value < 0.01 ~ "**", FDR_Value < 0.05 ~ "*", TRUE ~ ""),
         Cell_Text = sprintf("%.2f%s", Effect_Size, Significance_Label),
         Predictor = factor(Predictor, levels = predictors))
write.csv(reg, file.path(out_dir, "Module_Clinical_Regressions_long.csv"), row.names = FALSE)

# Wide table (one row per module: betas with FDR, genera and hub) -> Supplementary Table 3
deg <- rbind(data.frame(Genus = V(g_npa)$name, Niche = "Nasopharynx (NPA)", Degree = degree(g_npa), Structural_Module = V(g_npa)$module),
             data.frame(Genus = V(g_gut)$name, Niche = "Gut (GUT)", Degree = degree(g_gut), Structural_Module = V(g_gut)$module))
mods <- deg %>% filter(Structural_Module != 0) %>% group_by(Niche, Structural_Module) %>%
  summarise(Genera_In_Module = paste(Genus, collapse = " "), Hub_Genus = Genus[which.max(Degree)], .groups = "drop")
wide <- reg %>% mutate(String_Report = sprintf("%.3f (FDR: %.4f)", Effect_Size, FDR_Value),
                       Module_Num = as.numeric(gsub("NPA_Mod|GUT_Mod", "", Module)),
                       Niche = ifelse(grepl("^NPA", Module), "Nasopharynx (NPA)", "Gut (GUT)")) %>%
  select(Module_Num, Niche, Predictor, String_Report) %>%
  pivot_wider(id_cols = c(Module_Num, Niche), names_from = Predictor, values_from = String_Report, names_prefix = "Beta_")
master <- mods %>% left_join(wide, by = c("Structural_Module" = "Module_Num", "Niche")) %>%
  mutate(Structural_Module = paste0(ifelse(Niche == "Nasopharynx (NPA)", "NPA_Mod", "GUT_Mod"), Structural_Module)) %>%
  arrange(desc(Niche), Structural_Module) %>% relocate(Genera_In_Module, Hub_Genus, .after = last_col())
write.csv(master, file.path(out_dir, "Module_Clinical_Associations.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Global network topology
# ------------------------------------------------------------------------------
topology <- function(g, label) {
  g_sub <- delete_vertices(g, V(g)[degree(g) == 0])
  data.frame(Niche = label, Total_Input_Genera = vcount(g), Connected_Genera = vcount(g_sub),
             Isolated_Genera = vcount(g) - vcount(g_sub), Pct_Isolated = round(100 * (vcount(g) - vcount(g_sub)) / vcount(g), 2),
             Total_Edges = ecount(g_sub), Positive_Edges = sum(E(g_sub)$weight > 0), Negative_Edges = sum(E(g_sub)$weight < 0),
             Density = round(edge_density(g_sub), 4), Clustering = round(transitivity(g_sub, type = "global"), 4),
             Mean_Path_Length = round(mean_distance(g_sub, directed = FALSE, weights = NA), 4),
             Degree_Centralization = round(centr_degree(g_sub)$centralization, 4),
             Modules = length(setdiff(unique(V(g)$module), 0)))
}
write.csv(rbind(topology(g_npa, "Nasopharynx (NPA)"), topology(g_gut, "Gut (GUT)")),
          file.path(out_dir, "Network_Topology.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. Figure 6 (NPA network + module-trait heatmap) and Supplementary Figure 5 (gut network)
# ------------------------------------------------------------------------------
# [Revision 1] Network layout. The submitted figures used a weighted Fruchterman-Reingold
# layout in base igraph; with the rebuilt networks (one connected component, many weak
# edges) that layout collapsed the modules into overlapping clusters with unreadable
# labels. Modules are now placed apart (Fruchterman-Reingold on the module graph) and
# genera are laid out inside each module; labels are repelled so that none overlap.
# Colours, edge signs/widths and node sizes (mean abundance) are as in the submission.
module_layout <- function(g, seed, spread = 1.6, gap = 1.5) {
  set.seed(seed)
  mods <- sort(unique(V(g)$module))
  # genera inside each module (Kamada-Kawai), scaled on the 75th centile so that dense cores are expanded
  inner <- lapply(mods, function(m) {
    idx <- which(V(g)$module == m)
    if (length(idx) == 1) return(matrix(0, 1, 2))
    l <- layout_with_kk(induced_subgraph(g, idx), weights = NA)
    l <- sweep(l, 2, colMeans(l)); nr <- sqrt(rowSums(l^2))
    l / quantile(nr, 0.75) * length(idx)^0.65 * spread * 0.7
  })
  radius <- vapply(inner, function(l) max(sqrt(rowSums(l^2)), 0.5), 0)
  # module centres: Fruchterman-Reingold on the module graph, then pushed apart until no two modules overlap
  mg <- contract(g, match(V(g)$module, mods), vertex.attr.comb = "ignore")
  E(mg)$weight <- 1
  mg <- simplify(mg, edge.attr.comb = list(weight = "sum"))
  centres <- layout_with_fr(mg, weights = E(mg)$weight, niter = 5000)
  centres <- centres / max(dist(centres)) * sum(radius)
  repeat {
    moved <- FALSE
    for (i in seq_along(mods)) for (j in seq_along(mods)) if (i < j) {
      d <- centres[j, ] - centres[i, ]; dd <- sqrt(sum(d^2)); need <- radius[i] + radius[j] + gap
      if (dd < need - 1e-6) {
        sh <- (need - dd) / 2 * d / max(dd, 1e-6); centres[i, ] <- centres[i, ] - sh; centres[j, ] <- centres[j, ] + sh; moved <- TRUE
      }
    }
    if (!moved) break
  }
  xy <- matrix(0, vcount(g), 2)
  for (k in seq_along(mods)) xy[V(g)$module == mods[k], ] <- sweep(inner[[k]], 2, centres[k, ], "+")
  xy
}
wrap_label <- function(x) {   # long names on two lines, split at the space nearest the middle
  vapply(x, function(s) {
    sp <- gregexpr(" ", s)[[1]]
    if (nchar(s) <= 18 || sp[1] < 0) return(s)
    k <- sp[which.min(abs(sp - nchar(s) / 2))]
    paste0(substr(s, 1, k - 1), "\n", substr(s, k + 1, nchar(s)))
  }, "")
}
plot_network <- function(g, otu_mat, palette, seed, legend_title = "Modules") {
  g <- delete_vertices(g, V(g)[V(g)$module == 0])
  xy <- module_layout(g, seed)
  mods <- sort(unique(V(g)$module))
  cols <- setNames(colorRampPalette(brewer.pal(min(8, length(mods)), palette))(length(mods)), mods)
  nodes <- data.frame(x = xy[, 1], y = xy[, 2], label = wrap_label(gsub("_$", "", gsub("_", " ", V(g)$name))),
                      module = factor(V(g)$module, levels = mods),
                      abundance = log1p(colMeans(otu_mat, na.rm = TRUE))[V(g)$name])
  nodes$size <- 2.5 + 7.5 * (nodes$abundance - min(nodes$abundance)) / diff(range(nodes$abundance))   # point size (mm)
  el <- as_edgelist(g, names = FALSE)
  edges <- data.frame(x = xy[el[, 1], 1], y = xy[el[, 1], 2], xend = xy[el[, 2], 1], yend = xy[el[, 2], 2],
                      w = abs(E(g)$weight), sign = ifelse(E(g)$weight > 0, "Alignment (+)", "Exclusion (-)"))
  ggplot() +
    geom_segment(data = edges, aes(x, y, xend = xend, yend = yend, linewidth = w, colour = sign), alpha = 0.75) +
    scale_colour_manual(values = c("Alignment (+)" = "#2ecc71", "Exclusion (-)" = "#e74c3c"), name = "Interactions", drop = FALSE) +
    scale_linewidth(range = c(0.2, 2.4), guide = "none") +
    geom_point(data = nodes, aes(x, y, fill = module, size = size), shape = 21, colour = "#1a252f", stroke = 0.5) +
    scale_fill_manual(values = cols, labels = paste("Mod", mods), name = legend_title) +
    scale_size_identity() +
    ggrepel::geom_text_repel(data = nodes, aes(x, y, label = label, point.size = size), size = 3.2, fontface = "bold", family = "serif", lineheight = 0.85,
                             colour = "#2c3e50", max.overlaps = Inf, box.padding = 0.5, point.padding = 0.25, force = 3, max.iter = 50000, max.time = 5,
                             min.segment.length = 0.3, segment.colour = "grey60", seed = seed) +
    coord_cartesian(clip = "off") + theme_void(base_size = 13) +
    theme(legend.title = element_text(face = "bold"), plot.margin = margin(10, 10, 10, 10)) +
    guides(colour = guide_legend(order = 1, override.aes = list(linewidth = 2.5)),
           fill = guide_legend(order = 2, override.aes = list(size = 4.5)))
}
p_network <- plot_network(g_npa, otu_npa, "Set1", seed = 42)
p_heat <- ggplot(reg, aes(x = Predictor, y = Module, fill = Effect_Size)) +
  geom_tile(color = "white", lwd = 0.6) +
  geom_text(aes(label = Cell_Text), color = "black", size = 4.0, fontface = "bold") +
  scale_fill_gradient2(low = "#2c3e50", mid = "#f1c40f", high = "#d35400", midpoint = 0, name = "Standardized\nEffect Size (β)") +
  theme_minimal(base_size = 14) +
  labs(x = "Patient Clinical Traits & Contextual Covariates", y = "Structural Network Discovery Modules") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", color = "#1a252f", size = 12),
        axis.text.y = element_text(face = "bold", color = "#1a252f", size = 12),
        axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
        axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 15)),
        panel.grid = element_blank(), legend.title = element_text(face = "bold", size = 11),
        plot.margin = margin(5, 5, 5, 15))
fig6 <- plot_grid(p_network, p_heat, labels = c("a", "b"), label_size = 20, label_fontface = "bold",
                  ncol = 1, rel_heights = c(1.3, 1.0))
save_figure(fig6, 6, width = 12, height = 16)

save_figure(plot_network(g_gut, otu_gut, "Dark2", seed = 88, legend_title = "GUT Modules"),
            5, width = 11.3, height = 8, supplementary = TRUE)
