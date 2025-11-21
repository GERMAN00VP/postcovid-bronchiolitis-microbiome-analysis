setwd("/home/ucm-user/Documentos/BQL_ANALYSIS")

library(phyloseq)

# Leer los datos desde los archivos
otu_table <- read.csv("data/phyloseq_data/otu_table.csv", row.names = 1, check.names = FALSE)
sample_data <- read.csv("data/phyloseq_data/sample_data.csv", row.names = 1)
tax_table <- read.csv("data/phyloseq_data/tax_table.csv", row.names = 1)

# Crear los objetos para phyloseq
otu_table_phy <- otu_table(as.matrix(otu_table), taxa_are_rows = FALSE)
sample_data_phy <- sample_data(sample_data)
tax_table_phy <- tax_table(as.matrix(tax_table))

# Crear el objeto phyloseq
physeq <- phyloseq(otu_table_phy, sample_data_phy, tax_table_phy)

# Ver el objeto phyloseq

physeq

# Guardar el objeto phyloseq en un archivo .rds
saveRDS(physeq, file = "data/phyloseq_object.rds")
