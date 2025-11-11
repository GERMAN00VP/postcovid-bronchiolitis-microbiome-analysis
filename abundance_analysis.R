
library(ANCOMBC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(phyloseq)

set.seed(123)

setwd("/mnt/usb/BQL/BQL_ANALYSIS")

if (file.exists("Abundances")){
} else {
    dir.create("Abundances")   
}

### Charge the data
pseq = readRDS("data/phyloseq_object.rds")

# ANF BQL COMPARISON
##########################
# Filtrar el objeto sample_data dentro del objeto phyloseq
samp_anf<- sample_data(pseq)[sample_data(pseq)$Sample_Type == "ANF", ]

# Crear un nuevo objeto phyloseq con el sample_data filtrado
pseq_anf <- phyloseq(otu_table(pseq), sample_data(samp_anf), tax_table(pseq))

pseq_anf


res <- ancombc2(data = pseq_anf, tax_level = "Genus",
                fix_formula = 'Bronchiolitis+Wheezing.treatment+Respiratory.syncytial.virus+Family.history.atopy+Breastfeeding+Cesarean.section+Age_older' , 
                rand_formula = NULL,
                p_adj_method = "hochberg", pseudo_sens = TRUE, prv_cut = 0.1, s0_perc = 0.05,
                   group = "Bronchiolitis", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 5, verbose = TRUE,
                   iter_control = list(tol = 1e-3, max_iter = 100, 
                                       verbose = TRUE),
                   mdfdr_control= list(fwer_ctrl_method="hochberg",
                                       B=1000))$res

write.csv(res,file="./Abundances/ANF_BQL_ANCOMBC2.csv")

# ANF WHEEZING COMPARISON
##########################
# Filtrar el objeto sample_data dentro del objeto phyloseq
samp_anf_wheez<- sample_data(pseq_anf)[sample_data(pseq_anf)$Bronchiolitis == "Yes", ]

# Crear un nuevo objeto phyloseq con el sample_data filtrado
pseq_anf_wheez <- phyloseq(otu_table(pseq_anf), sample_data(samp_anf_wheez), 
                           tax_table(pseq_anf))

pseq_anf_wheez


res <- ancombc2(data = pseq_anf_wheez, tax_level = "Genus",
                fix_formula = 'Wheezing.treatment + 
                    Respiratory.syncytial.virus + Family.history.atopy + 
                    Breastfeeding + Cesarean.section + 
                    Age_older' , 
                rand_formula = NULL,
                p_adj_method = "hochberg", pseudo_sens = TRUE, prv_cut = 0.1, s0_perc = 0.05,
                group = "Wheezing.treatment", struc_zero = TRUE, neg_lb = TRUE,
                alpha = 0.05, n_cl = 5, verbose = TRUE,
                iter_control = list(tol = 1e-3, max_iter = 100, 
                                    verbose = TRUE),
                mdfdr_control= list(fwer_ctrl_method="hochberg",
                                    B=1000))$res

write.csv(res,file="./Abundances/ANF_WHEEZING_ANCOMBC2.csv")


###################################################
######### GUT
###################################################

# GUT BQL COMPARISON
##########################
# Filtrar el objeto sample_data dentro del objeto phyloseq
samp_gut <- sample_data(pseq)[sample_data(pseq)$Sample_Type == "GUT", ]

# Crear un nuevo objeto phyloseq con el sample_data filtrado
pseq_gut <- phyloseq(otu_table(pseq), sample_data(samp_gut), tax_table(pseq))

pseq_gut



res <- ancombc2(data = pseq_gut, tax_level = "Genus",
                fix_formula = 'Bronchiolitis+Wheezing.treatment+Family.history.atopy+Breastfeeding+Cesarean.section+Previous.antibiotics+Age_older' , 
                rand_formula = NULL,
                p_adj_method = "hochberg", pseudo_sens = TRUE, prv_cut = 0.1, s0_perc = 0.05,
                group = "Bronchiolitis", struc_zero = TRUE, neg_lb = TRUE,
                alpha = 0.05, n_cl = 5, verbose = TRUE,
                iter_control = list(tol = 1e-3, max_iter = 100, 
                                    verbose = TRUE),
                mdfdr_control= list(fwer_ctrl_method="hochberg",
                                    B=1000))$res

write.csv(res,file="./Abundances/GUT_BQL_ANCOMBC2.csv")

# GUT WHEEZING COMPARISON
##########################
# Filtrar el objeto sample_data dentro del objeto phyloseq
samp_gut_wheez<- sample_data(pseq_gut)[sample_data(pseq_gut)$Bronchiolitis == "Yes", ]

# Crear un nuevo objeto phyloseq con el sample_data filtrado
pseq_gut_wheez <- phyloseq(otu_table(pseq_gut), sample_data(samp_gut_wheez), 
                           tax_table(pseq_gut))

pseq_gut_wheez


res <- ancombc2(data = pseq_gut_wheez, tax_level = "Genus",
                fix_formula = 'Wheezing.treatment + 
                    Respiratory.syncytial.virus + Family.history.atopy + 
                    Breastfeeding + Cesarean.section + 
                    Previous.antibiotics+Age_older' , 
                rand_formula = NULL,
                p_adj_method = "hochberg", pseudo_sens = TRUE, prv_cut = 0.1, s0_perc = 0.05,
                group = "Wheezing.treatment", struc_zero = TRUE, neg_lb = TRUE,
                alpha = 0.05, n_cl = 5, verbose = TRUE,
                iter_control = list(tol = 1e-3, max_iter = 100, 
                                    verbose = TRUE),
                mdfdr_control= list(fwer_ctrl_method="hochberg",
                                    B=1000))$res

write.csv(res,file="./Abundances/GUT_WHEEZING_ANCOMBC2.csv")
