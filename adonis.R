library(vegan)
library(dplyr)

set.seed(123)

setwd("/home/ucm-user/Documentos/BQL_ANALYSIS")




#  SAMPLE TYPE COMPARISON
###########################################

metadata <- read.csv("data/phyloseq_data/sample_data.csv",sep=",",header = TRUE,
                     row.names = 1)

data <- read.csv("Beta_div/Distances/Distance_matrix.csv",
                     sep=",",header = TRUE,skip=0, row.names = 1)

## Run ADONIS tests
adonis <- adonis2(data ~ Sample_Type+Age,
                  data = metadata, permutations = 999,
                  by = "margin" # Does the correction for all variables
)

adonis

beta_disp <-  betadisper(as.dist(data), metadata %>%  pull(Sample_Type))

permutest(beta_disp)



#  ANF BQL COMPARISON
###########################################

metadata <- read.csv("data/phyloseq_data/sample_data.csv",sep=",",header = TRUE,
                     row.names = 1)

data_anf <- read.csv("Beta_div/Distances/Distance_matrix_ANF.csv",
                 sep=",",header = TRUE,skip=0, row.names = 1)

# Subset the metadata to keep only target patients data
metadata_anf <- metadata[rownames(data_anf),]
metadata_anf
## Run ADONIS tests
adonis <- adonis2(data_anf ~ Bronchiolitis+Wheezing.treatment+Respiratory.syncytial.virus+Family.history.atopy+Breastfeeding
                  +Cesarean.section+Age,
                  data = metadata_anf, permutations = 999,set.seed(1),
                  by = "margin" # Does the correction for all variables
                  )

adonis

capture.output(adonis, file = "Beta_div/ANF_BQL_ADONIS.txt")


## Beta dispersion 
beta_disp <- betadisper(as.dist(data_anf), metadata_anf %>%  pull(Bronchiolitis))

permutest(beta_disp)
capture.output(permutest(beta_disp), 
               file = "Beta_div/ANF_BQL_betadisper.txt")


## Beta dispersion 
beta_disp <- betadisper(as.dist(data_anf), metadata_anf %>%  pull(Wheezing.treatment))

permutest(beta_disp)
capture.output(permutest(beta_disp), 
               file = "Beta_div/ANF_Wheezing_betadisper.txt")

## Beta dispersion 
beta_disp <-  betadisper(as.dist(data_anf), metadata_anf %>%  pull(Respiratory.syncytial.virus))
permutest(beta_disp)
capture.output(permutest(beta_disp), 
               file = "Beta_div/ANF_rsv_betadisper.txt")


# Not necessary if betadisp isnt significative but is informative
anosim <- anosim(data_anf, metadata_anf$Bronchiolitis)
summary(anosim)
capture.output(summary(anosim), 
               file = "Beta_div/ANF_BQL_ANOSIM.txt")

png("Beta_div/ANF_BQL_ANOSIM.png", width = 600, height = 670)

par(cex = 1.5)  # Aumenta el tamaño de todos los textos en el gráfico

# Crear el gráfico
plot(anosim, xlab = "", ylab = "Dissimilarity Rank")

dev.off()


#  ANF WHEEZING COMPARISON 
#################################################################################

data_anf_wheez <- read.csv("Beta_div/Distances/Distance_matrix_ANF_BRQ.csv"
                           ,sep=",",header = TRUE,skip=0, row.names = 1)
metadata_anf_wheez <- metadata[rownames(data_anf_wheez),]


## Run ADONIS tests
adonis <- adonis2(data_anf_wheez ~ Wheezing.treatment + Respiratory.syncytial.virus + Family.history.atopy + 
                  Breastfeeding + Cesarean.section +
                  Age_older, # Formula
                  data = metadata_anf_wheez, permutations = 999,
                  by = "margin" # Does the correction for all variables
)
adonis

capture.output(adonis, 
               file = "Beta_div/ANF_WHEEZ_ADONIS.txt")


## Beta dispersion 
beta_disp <- betadisper(as.dist(data_anf_wheez), 
                        metadata_anf_wheez %>%  pull(Wheezing.treatment))

permutest(beta_disp)
capture.output(permutest(beta_disp), 
               file = "Beta_div/ANF_WHEEZ_betadisper.txt")

# Betadisp plot
png("Beta_div/ANF_WHEEZ_Betadisper.png", 
    width = 600, height = 670)

# Personalizar los parámetros gráficos
par(cex = 1.5)  # Aumentar el tamaño del texto globalmente

# Graficar la beta dispersión
plot(beta_disp, main = "Beta Dispersion by Wheezing Status",sub="",
     xlab = "PCoA Axis 1", ylab = "PCoA Axis 2", ellipse = TRUE, hull = FALSE)


dev.off()



## Beta dispersion 
beta_disp <- betadisper(as.dist(data_anf_wheez), 
                        metadata_anf_wheez %>%  pull(Respiratory.syncytial.virus))

permutest(beta_disp)

# Not necessary if betadisp isnt significative but is informative
anosim <- anosim(data_anf_wheez, metadata_anf_wheez$Respiratory.syncytial.virus)
summary(anosim)

capture.output(summary(anosim), 
               file = "Beta_div/ANF_WHEEZ_ANOSIM.txt")

png("Beta_div/ANF_WHEEZ_ANOSIM.png", width = 600, height = 670)

par(cex = 1.5)  # Aumenta el tamaño de todos los textos en el gráfico

# Crear el gráfico
plot(anosim, xlab = "", ylab = "Dissimilarity Rank")

dev.off()

# GUT MICROBIOTA

###########

#  GUT BRONCHIOLITIS COMPARISON
##################################################################################


data_gut <- read.csv("Beta_div/Distances/Distance_matrix_GUT.csv",
                     sep=",",header = TRUE,skip=0, row.names = 1)

# Subset the metadata to keep only target patients data
metadata_gut <- metadata[rownames(data_gut),]

## Run ADONIS tests
adonis <- adonis2(data_gut ~ Bronchiolitis+Wheezing.treatment+Respiratory.syncytial.virus+Family.history.atopy+Breastfeeding
                  +Cesarean.section+ Previous.antibiotics +Age_older,
                  data = metadata_gut, permutations = 999,
                  by = "margin" # Does the correction for all variables
)

adonis

capture.output(adonis, file = "Beta_div/GUT_BQL_ADONIS.txt")


## Beta dispersion 
beta_disp <- betadisper(as.dist(data_gut), 
                        metadata_gut %>%  pull(Bronchiolitis))

permutest(beta_disp)
capture.output(permutest(beta_disp), 
               file = "Beta_div/GUT_BQL_betadisper.txt")


# Betadisp plot
png("Beta_div/GUT_BQL_Betadisp.png", 
    width = 600, height = 670)

# Personalizar los parámetros gráficos
par(cex = 1.5)  # Aumentar el tamaño del texto globalmente

# Graficar la beta dispersión
plot(beta_disp, main = "Beta Dispersion by Bronchiolitis Status",sub="",
     xlab = "PCoA Axis 1", ylab = "PCoA Axis 2", ellipse = TRUE, hull = FALSE)


dev.off()

# Not necessary if betadisp isnt significative but is informative
anosim <- anosim(data_gut, metadata_gut$Bronchiolitis)
summary(anosim)
capture.output(summary(anosim), 
               file = "Beta_div/GUT_BQL_ANOSIM.txt")

png("Beta_div/GUT_BQL_ANOSIM.png", 
    width = 800, height = 900)

par(cex = 1.5)  # Aumenta el tamaño de todos los textos en el gráfico

# Crear el gráfico
plot(anosim, xlab = "", ylab = "Dissimilarity Rank")

dev.off()


#  GUT WHEEZING COMPARISON 
#################################################################################

data_gut_wheez <- read.csv("Beta_div/Distances/Distance_matrix_GUT_BRQ.csv"
                           ,sep=",",header = TRUE,skip=0, row.names = 1)
metadata_gut_wheez <- metadata[rownames(data_gut_wheez),]


## Run ADONIS tests
adonis <- adonis2(data_gut_wheez ~ Wheezing.treatment  + 
                    Respiratory.syncytial.virus + Family.history.atopy + 
                    Breastfeeding + Cesarean.section + Previous.antibiotics+Age_older, # Formula
                  data = metadata_gut_wheez, permutations = 999,
                  by = "margin" # Does the correction for all variables
)
adonis

capture.output(adonis, 
               file = "Beta_div/GUT_WHEEZ_ADONIS.txt")


## Beta dispersion 
beta_disp <- betadisper(as.dist(data_gut_wheez), 
                        metadata_gut_wheez %>%  pull(Cesarean.section))

permutest(beta_disp)
capture.output(permutest(beta_disp), 
               file = "Beta_div/GUT_WHEEZ_betadisper.txt")

# Betadisp plot
png("Beta_div/GUT_WHEEZ_Betadisper.png", 
    width = 600, height = 670)

# Personalizar los parámetros gráficos
par(cex = 1.5)  # Aumentar el tamaño del texto globalmente

# Graficar la beta dispersión
plot(beta_disp, main = "Beta Dispersion by Wheezing Status",sub="",
     xlab = "PCoA Axis 1", ylab = "PCoA Axis 2", ellipse = TRUE, hull = FALSE)


dev.off()

# Not necessary if betadisp isnt significative but is informative
anosim <- anosim(data_gut_wheez, metadata_gut_wheez$Wheezing.treatment)
summary(anosim)
capture.output(summary(anosim), 
               file = "Beta_div/GUT_WHEEZ_ANOSIM.txt")

png("Beta_div/GUT_WHEEZ_ANOSIM.png", width = 600, height = 670)

par(cex = 1.5)  # Aumenta el tamaño de todos los textos en el gráfico

# Crear el gráfico
plot(anosim, xlab = "", ylab = "Dissimilarity Rank")

dev.off()

