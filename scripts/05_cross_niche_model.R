# ==============================================================================
# 05. CROSS-NICHE MODEL FOR RECURRENT WHEEZING (Figure 5)
# ==============================================================================
# Bronchiolitis cases with both samples: logistic regression of recurrent wheezing on
# nasopharyngeal Burkholderia-Caballeronia-Paraburkholderia and gut [Ruminococcus]
# gnavus group (log10 relative abundance), adjusted for age, RSV, family history of
# atopy, passive smoking, caesarean delivery and previous antibiotics.
#
# [Revision 1] The ORs given in the submitted text (0.28 and 5.59) came from an earlier
# specification with breastfeeding instead of RSV; Figure 5, the Methods and the AUC
# corresponded to this RSV-adjusted model, which is the one reported now. After the
# QC revision (BQLHSO91GUT reinstated) the model includes 79 cases (49 with wheezing).
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
source("R/functions.R")
set.seed(123)

out_dir <- file.path(DIR_RESULTS, "05_cross_niche")
make_dirs(out_dir)

# ------------------------------------------------------------------------------
# 1. Model data
# ------------------------------------------------------------------------------
pseq <- readRDS(file.path(DIR_CURATED, "phyloseq_RAW_curated_global.rds"))
df_cn <- build_cross_niche_data(pseq)
# Age standardised with the mean and SD of all analysed samples, as in the reported model
all_age <- as.numeric(sample_data(pseq)$Age)
df_cn$Age <- (df_cn$Age_months - mean(all_age)) / sd(all_age)
cat("Cross-niche cohort:", nrow(df_cn), "cases,", sum(df_cn$Wheezing_Bin), "with recurrent wheezing\n")

# ------------------------------------------------------------------------------
# 2. Multivariable model
# ------------------------------------------------------------------------------
mv_model <- glm(as.formula(paste("Wheezing_Bin ~ Burkholderia_NPA + Ruminococcus_GUT + Age +", CROSS_NICHE_COVARIATES)),
                data = df_cn, family = binomial)
term_labels <- c(Burkholderia_NPA = "Burkholderia-\nCaballeronia-\nParaburkholderia\n(NPA)",
                 Ruminococcus_GUT = "[Ruminococcus]\ngnavus group\n(GUT)",
                 Age = "Age (per SD)",    # [Revision 1] axis label: age is standardised
                 RSV = "RSV\n(Yes)", Family.history.atopy = "Family History\nof Atopy (Yes)",
                 Passive.smoking = "Passive Smoking\n(Yes)", Cesarean.section = "Cesarean Section\n(Yes)",
                 Previous.antibiotics = "Previous Antibiotics\n(Yes)")
model_hits <- tidy(mv_model, conf.int = TRUE, exponentiate = TRUE) %>% filter(term != "(Intercept)")
write.csv(model_hits, file.path(out_dir, "Table_Cross_Niche_Model.csv"), row.names = FALSE)
print(model_hits)
model_hits$label <- factor(term_labels[model_hits$term], levels = rev(term_labels[model_hits$term]))

df_cn$pred_prob <- predict(mv_model, type = "response")
roc_obj <- roc(df_cn$Wheezing_Bin, df_cn$pred_prob, quiet = TRUE)
roc_ci  <- ci.se(roc_obj, specificities = seq(0, 1, l = 50))   # bootstrap band for the plot
auc_ci  <- ci.auc(roc_obj)
cat(sprintf("Cross-niche AUC = %.3f (95%% CI %.3f-%.3f)\n", auc_ci[2], auc_ci[1], auc_ci[3]))

# ------------------------------------------------------------------------------
# 3. Figure 5
# ------------------------------------------------------------------------------
p_forest <- ggplot(model_hits, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#e74c3c", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15, color = "#2c3e50", linewidth = 0.8) +
  geom_point(color = "#3498db", size = 3.5) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 50, 500),
                labels = c("0.001", "0.01", "0.1", "0.5", "1", "2", "5", "10", "50", "500")) +
  theme_bw(base_size = 11) + labs(x = "Odds Ratio (95% Confidence Interval)", y = "") +
  theme(axis.text.y = element_text(face = "bold", color = "black", lineheight = 0.85),
        axis.text.x = element_text(color = "black"), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "gray95"))
p_roc <- ggplot() +
  geom_ribbon(data = data.frame(Specificity = as.numeric(rownames(roc_ci)), Lower = roc_ci[, 1], Upper = roc_ci[, 3]),
              aes(x = 1 - Specificity, ymin = Lower, ymax = Upper), fill = "#2ecc71", alpha = 0.15) +
  geom_path(data = data.frame(Sensitivity = roc_obj$sensitivities, Specificity = roc_obj$specificities),
            aes(x = 1 - Specificity, y = Sensitivity), color = "#2ecc71", linewidth = 1.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60") +
  annotate("text", x = 0.95, y = 0.15, hjust = 1, size = 4.2, fontface = "bold", color = "#2c3e50",
           label = sprintf("AUC = %.3f\n95%% CI: [%.3f - %.3f]", auc_ci[2], auc_ci[1], auc_ci[3])) +
  theme_bw(base_size = 11) + labs(x = "1 - Specificity", y = "Sensitivity") +
  theme(panel.grid.minor = element_blank())
fig5 <- (p_forest / p_roc) + plot_layout(heights = c(1.8, 1)) +
  plot_annotation(tag_levels = list(c("a", "b"))) & theme(plot.tag = element_text(face = "bold", size = 13))
save_figure(fig5, 5, width = 7, height = 10)

