# ==============================================================================
# Shared functions for the bronchiolitis microbiome analysis
# ------------------------------------------------------------------------------
# Sourced by every script in scripts/. All paths are relative to the repository
# root. Comments tagged [Revision 1 - Rx.y] mark changes made in response to the
# peer review (Frontiers, revision 1): R1 = Reviewer 1, R2 = Reviewer 2, followed
# by the comment number (letters = recommendations of Reviewer 1).
# ==============================================================================

suppressPackageStartupMessages({
  library(splines)
  library(vegan)
})

# ------------------------------------------------------------------------------
# Constants and paths
# ------------------------------------------------------------------------------
GROUP_LEVELS <- c("CTRL", "RSV-/Wheeze-", "RSV-/Wheeze+", "RSV+/Wheeze-", "RSV+/Wheeze+")

# [Revision 1 - R1.2, R2.3] Age is modelled as a natural cubic spline (3 df) in the
# group comparisons of alpha diversity, instead of linear residuals.
AGE_SPLINE_DF <- 3

# [Revision 1 - R1.5/a] Pre-specified sample QC thresholds (see scripts/01_sample_qc.R)
QC_MIN_DEPTH    <- 10000   # minimum number of high-quality reads
QC_MIXUP_BC_MAX <- 0.10    # Bray-Curtis below which two samples are considered near-identical

DIR_RESULTS  <- "results"
DIR_CURATED  <- file.path(DIR_RESULTS, "curated_data")
DIR_FIG_MAIN <- file.path("figures", "main")
DIR_FIG_SUPP <- file.path("figures", "supplementary")
DIR_TAB_SUPP <- file.path("tables", "supplementary")

make_dirs <- function(...) {
  for (d in c(...)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Figures: every manuscript figure is written once, as JPEG, named by its number
# ------------------------------------------------------------------------------
# Note: plots with geom_jitter() draw random numbers when they are rendered, so each
# figure is rendered exactly once to keep the random-number stream (and therefore all
# downstream results) identical to the analysis reported in the manuscript.
save_figure <- function(plot, number, width, height, supplementary = FALSE, dpi = 300) {
  dir  <- if (supplementary) DIR_FIG_SUPP else DIR_FIG_MAIN
  name <- if (supplementary) paste0("Supplementary_Figure_", number) else paste0("Figure_", number)
  make_dirs(dir)
  file <- file.path(dir, paste0(name, ".jpg"))
  ggplot2::ggsave(file, plot = plot, width = width, height = height, dpi = dpi,
                  device = "jpeg", quality = 95, bg = "white")
  message("  figure written: ", file)
  invisible(file)
}

# Runs `expr` and then restores the random-number generator state. Used for the
# age-restricted analysis added during the revision, so that it can sit next to the
# main analysis it refers to without altering any downstream random draw.
with_preserved_seed <- function(expr) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  force(expr)
}

fmt_p <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# ------------------------------------------------------------------------------
# Metadata helpers
# ------------------------------------------------------------------------------
# Five clinical-virological phenotype groups defined by the index episode:
# healthy controls (CTRL) and bronchiolitis cases by RSV status and recurrent wheezing.
add_clinical_group <- function(md) {
  rsv_wh <- paste0(md$Respiratory.syncytial.virus, "_", md$Wheezing.treatment)
  grp <- dplyr::case_when(
    md$Bronchiolitis == "No" ~ "CTRL",
    rsv_wh == "No_No"        ~ "RSV-/Wheeze-",
    rsv_wh == "No_Yes"       ~ "RSV-/Wheeze+",
    rsv_wh == "Yes_No"       ~ "RSV+/Wheeze-",
    rsv_wh == "Yes_Yes"      ~ "RSV+/Wheeze+",
    TRUE                     ~ NA_character_
  )
  md$Condicion_Clinica <- factor(grp, levels = GROUP_LEVELS)
  md
}

extract_meta <- function(pseq_obj) {
  df <- data.frame(as(phyloseq::sample_data(pseq_obj), "data.frame"),
                   stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Sample_ID" %in% colnames(df)) df$Sample_ID <- rownames(df)
  rownames(df) <- df$Sample_ID
  df
}

clean_taxa_names <- function(x) {
  x <- gsub("[^[:alnum:]_]", "_", x)
  x <- gsub("_{2,}", "_", x)
  x <- gsub("^_|_$", "", x)
  make.unique(x, sep = "_dup")
}

# Uses the genus as taxon name; unassigned genera take the closest assigned rank
rename_to_genus <- function(pseq_obj, genus_col = "Genus") {
  tt <- as.data.frame(phyloseq::tax_table(pseq_obj), stringsAsFactors = FALSE)
  gen_labels <- tt[[genus_col]]
  higher_ranks <- rev(colnames(tt)[colnames(tt) != genus_col])
  bad <- "^NA$|^uncultured$|^unknown$"
  for (i in seq_along(gen_labels)) {
    if (is.na(gen_labels[i]) || gen_labels[i] == "" || grepl(bad, gen_labels[i], ignore.case = TRUE)) {
      for (rk in higher_ranks) {
        val <- tt[i, rk]
        if (!is.na(val) && val != "" && !grepl(bad, val, ignore.case = TRUE)) {
          gen_labels[i] <- paste0("Unclassified_", val)
          break
        }
      }
      if (is.na(gen_labels[i])) gen_labels[i] <- "Unclassified_Bacteria"
    }
  }
  phyloseq::taxa_names(pseq_obj) <- clean_taxa_names(gen_labels)
  rownames(phyloseq::tax_table(pseq_obj)) <- phyloseq::taxa_names(pseq_obj)
  pseq_obj
}

# ------------------------------------------------------------------------------
# [Revision 1 - R1.2, R2.3] Age-adjusted alpha-diversity comparisons
# ------------------------------------------------------------------------------
# For each metric fits  metric ~ group + ns(age)  (natural cubic spline) and returns
#   global   : F-test of the group term (model with vs. without group)
#   pairwise : every pairwise contrast (estimate, SE, p) with BH-FDR within metric.
# In an additive model the contrast between two groups is the difference of their
# coefficients, obtained here by re-levelling the group factor.
# Replaces the original two-step approach (residuals of metric ~ age, then Wilcoxon),
# which assumed a linear age effect and could remove group effects correlated with age.
alpha_age_adjusted_tests <- function(df, metrics, group = "Condicion_Clinica", age = "Age",
                                     spline_df = AGE_SPLINE_DF) {
  df[[age]] <- as.numeric(df[[age]])
  age_expr <- paste0("ns(", age, ", df = ", spline_df, ")")
  lv <- levels(droplevels(factor(df[[group]])))
  global <- list(); pairwise <- list()
  for (m in metrics) {
    d <- df[!is.na(df[[m]]) & !is.na(df[[age]]) & !is.na(df[[group]]), ]
    d[[group]] <- factor(d[[group]], levels = lv)
    f_full <- as.formula(paste0("`", m, "` ~ ", group, " + ", age_expr))
    f_red  <- as.formula(paste0("`", m, "` ~ ", age_expr))
    a <- anova(lm(f_red, data = d), lm(f_full, data = d))
    global[[m]] <- data.frame(Metric = m, n = nrow(d),
                              F_group = a$F[2], df_group = a$Df[2], p_group = a$`Pr(>F)`[2])
    combs <- combn(lv, 2)
    pw <- do.call(rbind, lapply(seq_len(ncol(combs)), function(i) {
      g1 <- combs[1, i]; g2 <- combs[2, i]
      d2 <- d; d2[[group]] <- relevel(d2[[group]], ref = g1)
      co <- summary(lm(f_full, data = d2))$coefficients
      term <- paste0(group, g2)
      data.frame(Metric = m, Group1 = g1, Group2 = g2,
                 Estimate_G2_minus_G1 = co[term, 1], SE = co[term, 2], p_raw = co[term, 4])
    }))
    pw$p_FDR <- p.adjust(pw$p_raw, method = "BH")
    pairwise[[m]] <- pw
  }
  pairwise <- do.call(rbind, pairwise)
  pairwise$Signif <- cut(pairwise$p_FDR, c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "ns"))
  list(global = do.call(rbind, global), pairwise = pairwise)
}

# Bracket heights for significant contrasts, computed per facet (metric)
bracket_positions <- function(sig, df_long, value = "Value", metric = "Metric", step = 0.08) {
  if (nrow(sig) == 0) return(sig)
  sig <- sig[order(sig[[metric]]), ]
  sig$y.position <- NA_real_
  for (m in unique(sig[[metric]])) {
    v <- df_long[[value]][df_long[[metric]] == m]
    rng <- diff(range(v, na.rm = TRUE)); top <- max(v, na.rm = TRUE)
    idx <- which(sig[[metric]] == m)
    sig$y.position[idx] <- top + rng * (0.06 + step * (seq_along(idx) - 1))
  }
  sig
}

# ------------------------------------------------------------------------------
# [Revision 1 - R2.3] Age-adjusted pairwise PERMANOVA with BH-FDR
# ------------------------------------------------------------------------------
# For each pair of groups: adonis2(dist ~ Age + Group, by = "margin") on the subset,
# i.e. the group effect conditional on age. In the original submission
# pairwiseAdonis::pairwise.adonis2 was used; that function ignores its `p.adjust.m`
# argument, so the p-values reported as "FDR" were in fact unadjusted.
pairwise_permanova_age <- function(dist_obj, meta, group = "Condicion_Clinica", age = "Age",
                                   permutations = 999, seed = 123) {
  m    <- as.matrix(dist_obj)
  meta <- meta[rownames(m), ]
  lv   <- levels(droplevels(factor(meta[[group]])))
  combs <- combn(lv, 2)
  out <- do.call(rbind, lapply(seq_len(ncol(combs)), function(i) {
    keep <- meta[[group]] %in% combs[, i]
    d  <- as.dist(m[keep, keep])
    md <- data.frame(Age = as.numeric(meta[[age]][keep]), Group = droplevels(factor(meta[[group]][keep])))
    set.seed(seed)
    a <- adonis2(d ~ Age + Group, data = md, permutations = permutations, by = "margin")
    data.frame(Comparison = paste(combs[1, i], "vs", combs[2, i]),
               n1 = sum(md$Group == combs[1, i]), n2 = sum(md$Group == combs[2, i]),
               F_Model = a["Group", "F"], R2 = a["Group", "R2"], p_raw = a["Group", "Pr(>F)"])
  }))
  out$p_FDR <- p.adjust(out$p_raw, method = "BH")
  out
}

# ------------------------------------------------------------------------------
# Diversity figure (alpha boxplots on observed values + PCoA), used for NPA and gut
# ------------------------------------------------------------------------------
plot_diversity_panel <- function(df_alpha_long, stat_results, df_pcoa, var_pct) {
  p_alpha <- ggplot(df_alpha_long, aes(x = Condicion_Clinica, y = Value, fill = Condicion_Clinica)) +
    geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.55, color = "#2c3e50", lwd = 0.6) +
    geom_jitter(width = 0.12, alpha = 0.25, color = "black", size = 0.8) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
    theme_bw(base_size = 11) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = "Alpha Diversity Value") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1, size = 9, color = "black"),
          axis.text.y = element_text(size = 9, color = "black"),
          strip.background = element_rect(fill = "#2c3e50", color = NA),
          strip.text = element_text(color = "white", face = "bold", size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#f0f0f0")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))
  # Brackets only for age-adjusted contrasts that survive FDR (none in the final analysis)
  if (nrow(stat_results) > 0) {
    p_alpha <- p_alpha + ggpubr::stat_pvalue_manual(stat_results, label = "p.adj.signif", y.position = "y.position",
                                                    hide.ns = TRUE, tip.length = 0.01, size = 3)
  }
  p_pcoa <- ggplot(df_pcoa, aes(x = PCoA1, y = PCoA2, color = Clinical_Group, fill = Clinical_Group)) +
    geom_point(size = 2.2, alpha = 0.65) +
    stat_ellipse(geom = "polygon", alpha = 0.08, level = 0.95, linewidth = 0.4, linetype = 2) +
    theme_bw(base_size = 11) +
    scale_color_brewer(palette = "Set2", name = "Clinical Phenotype Group") +
    scale_fill_brewer(palette = "Set2", name = "Clinical Phenotype Group") +
    labs(x = paste0("PCoA 1 (", var_pct[1], "%)"), y = paste0("PCoA 2 (", var_pct[2], "%)")) +
    theme(axis.text = element_text(size = 9, color = "black"),
          axis.title = element_text(size = 10, face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#f0f0f0"),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 9),
          legend.position = "right")
  (p_alpha / p_pcoa) +
    patchwork::plot_layout(ncol = 1, heights = c(0.7, 1)) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 14, color = "#2c3e50"))
}

# ------------------------------------------------------------------------------
# ANCOM-BC2 pairwise results for the robust globally significant genera (long format)
# ------------------------------------------------------------------------------
# Columns named "<stat>_Condicion_Clinica<A>" are contrasts of group A vs. the reference
# (CTRL); "<stat>_Condicion_Clinica<A>_Condicion_Clinica<B>" are contrasts A vs. B.
ancom_pairwise_long <- function(res_pair, taxa) {
  res_pair %>%
    dplyr::filter(taxon %in% taxa) %>%
    tidyr::pivot_longer(cols = -taxon, names_to = c(".value", "Comparacion"),
                        names_pattern = "^(lfc|se|W|p|q|diff|passed_ss)_(.*)") %>%
    dplyr::mutate(Comparacion = stringr::str_remove_all(Comparacion, "Condicion_Clinica")) %>%
    dplyr::select(Genus = taxon, Comparacion_Especifica = Comparacion,
                  `Log Fold Change (LFC)` = lfc, `W Statistic` = W, `Standard Error (SE)` = se,
                  `p-value` = p, `q-value (FDR)` = q, `Passed Sensitivity (SS)` = passed_ss) %>%
    dplyr::filter(!is.na(`Log Fold Change (LFC)`)) %>%
    dplyr::mutate(Group1 = ifelse(grepl("_", Comparacion_Especifica), sub("_.*$", "", Comparacion_Especifica), Comparacion_Especifica),
                  Group2 = ifelse(grepl("_", Comparacion_Especifica), sub("^[^_]*_", "", Comparacion_Especifica), "CTRL")) %>%
    dplyr::arrange(Genus, `p-value`)
}

# Bias-corrected log abundances from ANCOM-BC2; missing values (taxa absent from a
# sample) are imputed at the taxon minimum plus small Gaussian noise
ancom_log_abundances <- function(ancom_out) {
  m <- as.data.frame(ancom_out$bias_correct_log_table)
  imp <- as.data.frame(lapply(m, function(x) {
    if (any(is.na(x))) {
      min_val <- min(x, na.rm = TRUE)
      if (is.infinite(min_val)) min_val <- -10
      n_nas <- sum(is.na(x))
      if (n_nas > 0) x[is.na(x)] <- rnorm(n_nas, mean = min_val, sd = 0.01)
    }
    x
  }))
  rownames(imp) <- rownames(m)
  as.data.frame(t(imp))
}

# ------------------------------------------------------------------------------
# Wheezing risk models within bronchiolitis cases (one niche): univariate age-adjusted
# models for each robust ANCOM-BC2 genus, multivariable model with the genera that are
# significant univariately, ROC curve. Returns tables and the three-panel figure.
# ------------------------------------------------------------------------------
fit_niche_risk_models <- function(df_model, taxa_cols) {
  uni <- list()
  for (bac in taxa_cols) {
    if (length(unique(df_model[[bac]])) < 2 || length(unique(df_model$Wheezing_Bin)) < 2) next
    fit <- tryCatch(glm(as.formula(paste0("Wheezing_Bin ~ Age + `", bac, "`")), data = df_model, family = binomial),
                    error = function(e) NULL)
    if (is.null(fit)) next
    res <- tryCatch(broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE),
                    error = function(e) broom::tidy(fit, exponentiate = TRUE, conf.int = FALSE))
    uni[[bac]] <- res %>% dplyr::filter(gsub("`", "", term) == bac) %>% dplyr::mutate(Bacteria = bac)
  }
  uni <- dplyr::bind_rows(uni) %>%
    dplyr::select(Bacteria, `Odds Ratio (OR)` = estimate, `Lower CI (2.5%)` = conf.low,
                  `Upper CI (97.5%)` = conf.high, `p-value` = p.value)
  sig <- uni$Bacteria[uni$`p-value` < 0.05]
  out <- list(univariate = uni, selected = sig)
  if (length(sig) > 0) {
    fit_m <- glm(as.formula(paste("Wheezing_Bin ~ Age +", paste0("`", sig, "`", collapse = " + "))),
                 data = df_model, family = binomial, na.action = na.exclude)
    out$multivariable <- broom::tidy(fit_m, exponentiate = TRUE, conf.int = TRUE) %>%
      dplyr::filter(term != "(Intercept)") %>%
      dplyr::select(term, `Odds Ratio (OR)` = estimate, `Lower CI (2.5%)` = conf.low,
                    `Upper CI (97.5%)` = conf.high, `p-value` = p.value)
    pred <- predict(fit_m, type = "response")
    ok <- !is.na(pred) & !is.na(df_model$Wheezing_Bin)
    out$roc <- pROC::roc(df_model$Wheezing_Bin[ok], pred[ok], levels = c(0, 1), direction = "<", quiet = TRUE)
    out$auc_ci <- pROC::ci.auc(out$roc)
  }
  out
}

taxon_label <- function(x) gsub("_", "-\n", gsub("_group", "", gsub("`", "", x)))

plot_risk_panel <- function(rm) {
  forest <- function(d, colour) {
    ggplot(d, aes(x = `Odds Ratio (OR)`, y = Label)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "#d63031", linewidth = 0.8) +
      geom_errorbarh(aes(xmin = `Lower CI (2.5%)`, xmax = `Upper CI (97.5%)`), height = 0.2, color = "#2c3e50", linewidth = 0.8) +
      geom_point(size = 3.5, color = colour) +
      theme_bw(base_size = 12) + labs(x = "Odds Ratio (95% CI)", y = NULL) + scale_x_log10() +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#f0f0f0"),
            axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 10)),
            axis.text.y = element_text(face = "italic", color = "black", size = 12, lineheight = 0.85),
            axis.text.x = element_text(color = "black", size = 12))
  }
  u <- rm$univariate %>% dplyr::mutate(Label = taxon_label(Bacteria))
  u$Label <- factor(u$Label, levels = rev(sort(unique(u$Label))))
  m <- rm$multivariable %>% dplyr::mutate(Label = ifelse(term == "Age", "Age\n(Adjustment)", taxon_label(term)))
  m$Label <- factor(m$Label, levels = rev(sort(unique(m$Label))))
  coords <- pROC::coords(rm$roc, "all", transpose = FALSE) %>% as.data.frame() %>%
    dplyr::transmute(One_Minus_Specificity = 1 - specificity, sensitivity) %>%
    rbind(data.frame(One_Minus_Specificity = c(0, 1), sensitivity = c(0, 1))) %>%
    dplyr::arrange(One_Minus_Specificity, sensitivity) %>% dplyr::distinct()
  ci <- rm$auc_ci
  p_roc <- ggplot(coords, aes(x = One_Minus_Specificity, y = sensitivity)) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey60", linewidth = 0.6) +
    geom_path(color = "#1f77b4", linewidth = 1.3, lineend = "round") +
    scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_bw(base_size = 12) + labs(x = "1 - Specificity", y = "Sensitivity") +
    annotate("text", x = 0.75, y = 0.15, size = 5, fontface = "bold", color = "#2c3e50",
             label = sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)", ci[2], ci[1], ci[3])) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#f0f0f0"),
          axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 10)),
          axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 10)),
          axis.text = element_text(color = "black", size = 12))
  ((forest(u, "#2980b9") | forest(m, "#27ae60")) / p_roc) +
    patchwork::plot_layout(heights = c(1, 1.2)) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 18, color = "#2c3e50"))
}

# ------------------------------------------------------------------------------
# Boxplots of age-adjusted, bias-corrected abundances of the robust ANCOM-BC2 genera
# ------------------------------------------------------------------------------
# Brackets show exploratory pairwise ANCOM-BC2 contrasts with nominal p < 0.05.
# Contrast labels are mapped directly from the ANCOM-BC2 column names
# (A vs. the reference group CTRL, or A vs. B).
plot_ancom_boxplots <- function(df_bact, metadata_base, taxa, pairs, display = identity,
                                ncol = 3, step_increase = 0.11) {
  safe <- make.names(display(taxa))
  d <- metadata_base %>% dplyr::select(-dplyr::any_of("Sample_ID")) %>% tibble::rownames_to_column("Sample_ID")
  ab <- df_bact[, taxa, drop = FALSE]; colnames(ab) <- safe; ab$Sample_ID <- rownames(df_bact)
  d <- dplyr::left_join(d, ab, by = "Sample_ID") %>% dplyr::mutate(Age = as.numeric(Age))
  # Age adjustment for display: residuals of abundance ~ age, shifted by the mean
  for (s in safe) d[[paste0(s, "_adj")]] <- residuals(lm(as.formula(paste0("`", s, "` ~ Age")), data = d, na.action = na.exclude)) +
    mean(d[[s]], na.rm = TRUE)
  map <- data.frame(Safe = safe, Bacteria = display(taxa))
  long <- d %>% dplyr::select(Sample_ID, Clinical_Group = Condicion_Clinica, dplyr::all_of(paste0(safe, "_adj"))) %>%
    tidyr::pivot_longer(-c(Sample_ID, Clinical_Group), names_to = "Safe", values_to = "Abundance") %>%
    dplyr::mutate(Safe = sub("_adj$", "", Safe)) %>% dplyr::left_join(map, by = "Safe")
  maxv <- long %>% dplyr::group_by(Bacteria) %>% dplyr::summarise(max_val = max(Abundance, na.rm = TRUE), .groups = "drop")
  ann <- pairs %>% dplyr::filter(`p-value` < 0.05) %>%
    dplyr::transmute(Bacteria = display(Genus), group1 = Group1, group2 = Group2,
                     label = dplyr::case_when(`p-value` < 0.001 ~ "***", `p-value` < 0.01 ~ "**", TRUE ~ "*")) %>%
    dplyr::inner_join(maxv, by = "Bacteria") %>% dplyr::group_by(Bacteria) %>%
    dplyr::mutate(y.position = max_val + abs(max_val) * (0.05 + (dplyr::row_number() - 1) * 0.15)) %>% dplyr::ungroup()
  letters_df <- data.frame(Bacteria = sort(unique(long$Bacteria)))
  letters_df$Letra <- letters[seq_len(nrow(letters_df))]
  p <- ggplot(long, aes(x = Clinical_Group, y = Abundance, fill = Clinical_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.3, color = "black", size = 1.0) +
    facet_wrap(~ Bacteria, scales = "free_y", ncol = ncol) +
    theme_bw() + scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = "Bias-Corrected Age-Adjusted Abundance (ANCOM-BC2 Residuals)") +
    theme(legend.position = "none",
          axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 10)),
          axis.text.y = element_text(size = 11, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
          panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "gray95"),
          strip.background = element_rect(fill = "#f8f9fa", color = "gray80"),
          strip.text = element_text(color = "black", face = "bold.italic", size = 13)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    geom_text(data = letters_df, aes(label = Letra), x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5,
              inherit.aes = FALSE, fontface = "bold", size = 6, color = "black")
  if (nrow(ann) > 0) p <- p + ggpubr::stat_pvalue_manual(ann, label = "label", y.position = "y.position",
                                                        step.increase = step_increase, hide.ns = TRUE,
                                                        tip.length = 0.02, size = 5)
  p
}

# ------------------------------------------------------------------------------
# Cross-niche model data: one row per bronchiolitis case with both samples, with the
# log10 relative abundance (%) of Burkholderia-Caballeronia-Paraburkholderia in the NPA
# and of the [Ruminococcus] gnavus group in stool, plus clinical covariates.
# ------------------------------------------------------------------------------
build_cross_niche_data <- function(pseq_raw) {
  pg <- phyloseq::tax_glom(pseq_raw, taxrank = "Genus", NArm = FALSE) %>% rename_to_genus("Genus")
  pp <- phyloseq::transform_sample_counts(pg, function(x) log10((x / sum(x)) * 100 + 1))
  otu <- as(phyloseq::otu_table(pp), "matrix"); if (!phyloseq::taxa_are_rows(pp)) otu <- t(otu)
  burk <- otu[grep("Burkholderia", rownames(otu), ignore.case = TRUE)[1], ]
  rumi <- otu[grep("Ruminococcus_gnavus", rownames(otu), ignore.case = TRUE)[1], ]
  md <- extract_meta(pp) %>% dplyr::filter(Bronchiolitis == "Yes")
  both <- md %>% dplyr::group_by(Patient_ID) %>% dplyr::filter(all(c("NPA", "GUT") %in% Sample_Type)) %>%
    dplyr::pull(Patient_ID) %>% unique()
  npa <- md %>% dplyr::filter(Sample_Type == "NPA", Patient_ID %in% both) %>%
    dplyr::transmute(Patient_ID, Age_months = as.numeric(Age), Burkholderia_NPA = burk[Sample_ID],
                     Wheezing_Bin = as.integer(Wheezing.treatment == "Yes"),
                     Wheezing_count = suppressWarnings(as.numeric(Wheezing.count)),
                     Respiratory.syncytial.virus, Family.history.atopy, Passive.smoking,
                     Cesarean.section, Previous.antibiotics)
  gut <- md %>% dplyr::filter(Sample_Type == "GUT", Patient_ID %in% both) %>%
    dplyr::transmute(Patient_ID, Ruminococcus_GUT = rumi[Sample_ID])
  dplyr::inner_join(npa, gut, by = "Patient_ID") %>%
    tidyr::drop_na(Burkholderia_NPA, Ruminococcus_GUT, Age_months, Wheezing_Bin, Respiratory.syncytial.virus,
                   Family.history.atopy, Passive.smoking, Cesarean.section, Previous.antibiotics) %>%
    dplyr::mutate(RSV = as.integer(Respiratory.syncytial.virus == "Yes"),
                  Family.history.atopy = as.integer(Family.history.atopy == "Yes"),
                  Passive.smoking = as.integer(Passive.smoking == "Yes"),
                  Cesarean.section = as.integer(Cesarean.section == "Yes"),
                  Previous.antibiotics = as.integer(Previous.antibiotics == "Yes"))
}

CROSS_NICHE_COVARIATES <- "RSV + Family.history.atopy + Passive.smoking + Cesarean.section + Previous.antibiotics"
