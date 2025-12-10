# ========================================
# Script completo R: Firth + Bootstrap + ROC + Forest
# ========================================

# 1. SET WORKING DIRECTORY
setwd("/home/ucm-user/Documentos/BQL_ANALYSIS/Predictive_microbiota_analysis/RSV")

# 2. CARGAR LIBRERÍAS
library(logistf)       # Firth logistic
library(boot)          # Bootstrap
library(pROC)          # ROC curve
library(ggplot2)       # Plots
library(dplyr)

# 3. LEER DATOS
data <- read.csv("Log2_micro_RSV_data.csv", row.names = 1)

# 4. SELECCIONAR VARIABLES POR LASSO
# Suponiendo que 'data' ya está leído
X <- data[, -ncol(data)]  # Todas las columnas excepto la última
y <- data[, ncol(data)]   # La última columna como target


# 5. FIRTH LOGISTIC REGRESSION
firth_model <- logistf(formula = y ~ ., data = cbind(y, X),
                       control=logistf.control(maxit=200, maxstep=0.5))
summary(firth_model)
firth_model$converged 
# Guardar resultados Firth
firth_res <- data.frame(
  Variable = names(firth_model$coef),
  OR = exp(firth_model$coef),
  CI_low = exp(firth_model$ci.lower),
  CI_high = exp(firth_model$ci.upper),
  p_value = firth_model$prob
)
write.csv(firth_res, "Firth_ORs.csv", row.names = FALSE)

# 6. BOOTSTRAP DE ORs (1000 iteraciones)
set.seed(123)
bootstrap_fun <- function(data, indices){
  d <- data[indices, ]
  tryCatch({
    model <- logistf(y ~ ., data = d,
                     control=logistf.control(maxit=200, maxstep=0.5))
    return(exp(model$coef))  # intercepto + predictores
  }, error=function(e) rep(NA, ncol(d)))  # ahora incluye intercepto
}


boot_results <- boot(data = cbind(y, X), statistic = bootstrap_fun, R = 1000)

# Calcular median OR y 95% CI
# Mediana OR
boot_ORs <- apply(boot_results$t, 2, function(x) median(x, na.rm = TRUE))

# 95% CI bootstrap
boot_CI_low <- apply(boot_results$t, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
boot_CI_high <- apply(boot_results$t, 2, function(x) quantile(x, 0.975, na.rm = TRUE))


# MODIFICAR: Bootstrap table con p-values
bootstrap_table <- data.frame(
  Variable = names(firth_model$coef),
  OR_bootstrap = boot_ORs,
  CI_low = boot_CI_low,
  CI_high = boot_CI_high,
  p_value = firth_model$prob  # Añadir p-values de Firth
)
# AÑADIR AMBOS TESTS GLOBALES
wald_row <- data.frame(
  Variable = "Wald test (global)",
  OR_bootstrap = NA,
  CI_low = NA,
  CI_high = NA,
  p_value = 0.001536626
)

lr_row <- data.frame(
  Variable = "Likelihood ratio test (global)",
  OR_bootstrap = NA, 
  CI_low = NA,
  CI_high = NA,
  p_value = 0.2187346 
)

bootstrap_table_complete <- rbind(bootstrap_table, wald_row,lr_row)

write.csv(bootstrap_table_complete, "Firth_ORs_bootstrap.csv", row.names = FALSE)

# 7. FOREST PLOT (bootstrap ORs) - MODIFICADO PARA PNG
forest_df <- bootstrap_table %>% filter(Variable != "(Intercept)") %>%
  arrange(OR_bootstrap)

# Cambiar pdf() por png() con resolución adecuada
png("Forest_plot_bootstrap.png", width=800, height=600, res=150)
ggplot(forest_df, aes(x = OR_bootstrap, y = reorder(Variable, OR_bootstrap))) +
  geom_point(size=3) +
  geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), height=0.2) +
  geom_vline(xintercept = 1, linetype="dashed", color="red") +
  scale_x_log10() +
  xlab("Odds Ratio (log scale)") +
  ylab("Genus") +
  ggtitle("Bootstrap Firth Logistic Regression — Forest Plot") +
  theme_minimal()
dev.off()

# 8. ROC CURVE + AUC + sensibilidad/especificidad

# fitted probabilities
probs <- predict(firth_model, type="response")
roc_obj <- roc(y, probs)

# Punto óptimo (Youden)
coords_opt <- coords(roc_obj, "best", ret=c("threshold","sensitivity","specificity"), best.method="youden")
x_plot <- as.numeric(coords_opt["specificity"])  # eje X = 1 - especificidad
y_plot <- as.numeric(coords_opt["sensitivity"])

# Guardar PNG con tamaño y resolución
png("ROC_curve_Firth_large.png", width=1000, height=1000, res=150)

# Plot ROC con letras grandes
plot(roc_obj,
     main=paste0("ROC Curve (AUC = ", round(auc(roc_obj),3), ")"),
     col="blue", lwd=2,
     cex.main=2,    # tamaño título
     cex.lab=1.5,   # tamaño ejes
     cex.axis=1.3,  # tamaño ticks
     xlab="1 - Specificity", ylab="Sensitivity")

# Punto óptimo
points(x_plot, y_plot, pch=19, col="red", cex=2)

# Leyenda con tamaño mayor
legend("bottomright", legend=c(
  paste0("Sensitivity=", round(y_plot,3)),
  paste0("Specificity=", round(coords_opt["specificity"],3))
), bty="n", cex=1.5)

dev.off()


# 9. Guardar ROC metrics
roc_metrics <- data.frame(
  AUC = auc(roc_obj),
  Threshold_opt = coords_opt["threshold"],
  Sensitivity_opt = coords_opt["sensitivity"],
  Specificity_opt = coords_opt["specificity"]
)
write.csv(roc_metrics, "ROC_metrics.csv", row.names=FALSE)

# ========================================
# FIN DEL SCRIPT
# ========================================
