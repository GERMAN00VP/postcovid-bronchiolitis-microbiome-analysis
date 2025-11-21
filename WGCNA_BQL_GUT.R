
setwd("/home/ucm-user/Documentos/BQL_ANALYSIS/WGCNA")

#Loading WGCNA
library(WGCNA)
library(dbplyr)
library(cluster)
library(ggplot2)
library(igraph)
library(dplyr)

#Setting string not as factor
options(stringsAsFactors = FALSE)

#Enable multithread
enableWGCNAThreads(nThreads = 30)

if (file.exists("results/GUT")){
} else {
    dir.create("results/GUT",recursive=TRUE)   
}


#Reading the raw data (rows are the sample and columns the genes)
expressiondata = as.matrix(read.csv("data/GUT_BQL_counts.csv",row.names = 1))
tax <- colnames(expressiondata)


#Group data in a dendogram to check outliers
sampleTree = hclust(dist(expressiondata), method =  "average")

#If you want to save this plot in a pdf file, do not comment the line below:
png(file = "results/GUT/sampleClustering.png",  width = 2100, height = 700,res = 120) 
par(cex = 0.6)
par(mar = c(0,4,2,0))
# Extract height values from the sample tree
heights <- sampleTree$height

# Set the cutoff at the 95th percentile
cutoff <- quantile(heights, 0.95)

# Plot the dendrogram with the cutoff line
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = cutoff, col = "red")

# Cerrar el dispositivo gráfico
dev.off()


#Reading the TRAIT DATA
traitData <- as.data.frame(read.csv("data/GUT_BQL_metadata.csv",row.names = 1))

# Keep the patient Info for the linear models
metadata <- traitData

# Filter out the patient data, its non numerical and gives error
traitData <- traitData[, colnames(traitData) != "Patient", drop = FALSE]

# Convierte a matriz numérica si no es ya numérica
traitData <- as.matrix(traitData)

# Identify samples to remove (samples above the cutoff)
samples_to_remove <- c("BQLHSO97GUT","HULP71GUT")

# Filter out the identified samples
expressiondata <- expressiondata[!(rownames(expressiondata) %in% samples_to_remove), ]
traitData <- traitData[!(rownames(traitData) %in% samples_to_remove), ]

# Assign the trait colors now we've done the filtering
traitColors = numbers2colors(traitData, signed = FALSE)


#Group data in a dendogram to check outliers
sampleTree = hclust(dist(expressiondata), method =  "ward.D2")

#If you want to save this plot in a pdf file, do not comment the line below:
png(file = "results/GUT/sampleClustering_traits.png",  width = 2300, height = 1200,res = 100) 


#Plot a sample dendogram with the colors below
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = colnames(traitData), 
                    main = "Sample dendrogram and trait heatmap")

dev.off()


# Choose a set of soft threshold parameters
powers = c(c(1:20), seq(from = 22, to=30, by=2))


sft = pickSoftThreshold(expressiondata, powerVector = powers,
                        networkType = "signed",
                        verbose = 5,
                        corFnc=cor, corOptions =   list(use="p"
                          ,method = "spearman")) 


png(file = "results/GUT/softthreshold.png",    width = 1200, height = 900) 
#Plotting the results
par(mfrow = c(1,2))
cex1 = 0.9

#Index the scale free topology adjust as a function of the power soft thresholding.
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")


#This line corresponds to use a cut-off R² of h
abline(h=0.80,col="red")

#Connectivity mean as a function of soft power thresholding
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



# Cerrar el dispositivo gráfico
dev.off()

png(file = "results/GUT/softthreshold_peack.png",    width = 1200, height = 900) 
#Plotting the results
cex1 = 0.9

SOFT_PEACK = sft$fitIndices[,5]*((-sign(sft$fitIndices[,3])*sft$fitIndices[,2]))

#Connectivity mean as a function of soft power thresholding
plot(sft$fitIndices[,1], SOFT_PEACK,
     xlab="Soft Threshold (power)",ylab="Mean Connectivity*Scale Free R²", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], SOFT_PEACK, labels=powers, cex=cex1,col="red")

SOFT_PEACK

# Cerrar el dispositivo gráfico
dev.off()


expressiondata <- as.matrix(expressiondata)
expressiondata <- apply(expressiondata, 2, as.numeric)


power= 6

print(paste("R2 at selected power:",
            (sft$fitIndices[sft$fitIndices$Power==power,"SFT.R.sq"])))
print(paste("Conectivity at selected power:",
            (sft$fitIndices[sft$fitIndices$Power==power,"median.k."])))





### STEP BY STEP
adjacency = adjacency(expressiondata, power = power, type = "signed",
                      corOptions = list(use = "p",method="spearman"))

###### PLOT THE FREE SCALE 
mi_red <-simplify(graph_from_adjacency_matrix(adjacency, mode = "undirected", weighted = TRUE, diag = FALSE))

png("results/GUT/freescaleplot.png")

# Calcular la fuerza (grado ponderado)
strength_values <- strength(mi_red, mode = "all", weights = E(mi_red)$weight)

# Histograma para distribución empírica
hist_data <- hist(strength_values, breaks = 50, plot = FALSE)
freq <- hist_data$counts / sum(hist_data$counts)  # Probabilidad P(s)
centers <- hist_data$mids  # Valores de fuerza representativos


# Graficar la distribución observada en log-log
plot(centers, freq, log = "xy", pch = 19, col = "blue",
     xlab = "Fuerza (grado ponderado)", 
     ylab = "P(Fuerza)",
     main = "Distribución de Fuerza - Ajuste Ley de Potencia")

dev.off()


TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
dissTOM = 1-TOM


# PERFORM WARD.D2 CLUSTERING
geneTree = hclust(as.dist(dissTOM), method =  "ward.D2");


######################
# CLUSTERING PARAMETERS SELECTION
##################

# Definir las combinaciones de parámetros
minClusterSizes <- seq(5, 25, 5)
deepSplits <- c(2, 3, 4)

# Crear una matriz para almacenar los resultados
results <- matrix(NA, nrow = length(minClusterSizes) * length(deepSplits), ncol = 2)
row_names <- vector("character", length = length(minClusterSizes) * length(deepSplits))

# Hacer el loop para probar todas las combinaciones de deepSplit y minClusterSize
counter <- 1
for (minClusterSize in minClusterSizes) {
  for (deepSplit in deepSplits) {
    
    # Hacer el clustering con los parámetros actuales
    dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = deepSplit, 
                                 minClusterSize = minClusterSize)
    
    # Calcular el índice de silueta
    silhouette_score <- silhouette(dynamicMods, dist(dissTOM))
    avg_silhouette <- mean(silhouette_score[, 3])  # Promedio del índice de silueta
    
    # Almacenar los resultados
    results[counter, 1] <- avg_silhouette
    results[counter, 2] <- paste("dS:", deepSplit, "CS:", minClusterSize)
    row_names[counter] <- paste("dS:", deepSplit, "CS:", minClusterSize)
    
    counter <- counter + 1
  }
}

# Asignar nombres a las filas
rownames(results) <- row_names

# Convertir la matriz de resultados en un dataframe para facilitar el gráfico
results_df <- as.data.frame(results)
colnames(results_df) <- c("Silhouette_Score", "Parameters")

# Crear un gráfico de los resultados
png( "results/GUT/MODULE_CLUSTER_silhouette.png",width = 900, height = 600)

ggplot(results_df, aes(x = Parameters, y = Silhouette_Score, group = 1)) +
  geom_line() + 
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Silhouette Score para diferentes combinaciones de deepSplit y minClusterSize", 
       x = "Parámetros (deepSplit, minClusterSize)", 
       y = "Índice de Silueta") +
  theme_minimal()

dev.off()


##############################
###################


# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 10);
table(dynamicMods)
length(table(dynamicMods)) 

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
png(file = "results/GUT/MODULE_CLUSTER_COLORS.png", width = 800, height = 600);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Genera dendrogram and module colors")
dev.off()


###########
# Guardar gráfico con ajustes para publicación
png("results/GUT/MODULE_CLUSTER_COLORS_PUBLICATION.png", width = 1200, height = 1000, res = 300)

# Crear gráfico sin título, con etiquetas claras y margen para incluir la letra "A"
par(mar = c(1, 2, 4, 1))  # Márgenes (inferior, izquierda, superior, derecha)

# Gráfico
plotDendroAndColors(
  dendro = geneTree, 
  colors = dynamicColors, 
  groupLabels = "",  # Sin título
  dendroLabels = FALSE,
  hang = 0.03, 
  addGuide = TRUE, 
  guideHang = 0.05,
  main=""
)


dev.off()
##########

# Calculate eigengenes
MEList = moduleEigengenes(expressiondata, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

#Clustering the eigengenes modules
METree = hclust(as.dist(MEDiss), method = "average")

#Grouping the clusters from a cut-off
MEDissThres = 0.5
#Plotting the result

if (length(METree$labels)>2){
  
  # Dendrogram of module clusterig
  png(file="results/GUT/Clustering_module_eigengenes.png",heigh=1200,width=900,res = 200)
  
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  
  
  #Plotting a cut-off line
  abline(h=MEDissThres, col = "red")
  dev.off()
  
}



merge = mergeCloseModules(expressiondata, dynamicColors, cutHeight = MEDissThres, 
                          verbose = 3) 

## New Dendrogram plot

#Grouping module colors
mergedColors = merge$colors

#Eigengenes of new grouped modules
mergedMEs = merge$newMEs


png(file = "results/GUT/MODULE_CLUSTER_MERGED_COLORS.png", width = 800, height = 600)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors


#Building numeric labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Define numbers of genes and samples
nGenes = ncol(expressiondata);
nSamples = nrow(expressiondata);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(expressiondata, mergedColors)$eigengenes
MEs = orderMEs(MEs0)

# Calculate spearman correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, traitData, use = "p",method = "pearson");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#  Plot heatmap of module-traits relationship
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

moduleTraitCor

png("results/GUT/module_traits_correlation.png", width = 1000, height = 900,res = 100)

par(mar = c(15, 12, 5, 5),cex=1.1);

mis_labels <- c("Bronchiolitis" , "Wheezing episodes","Respiratory syncytial virus","Cesarean section",
                "Previous antibiotics","Age", "Family history atopy","Breastfeeding" )

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = round(moduleTraitCor,3),
               xLabels = mis_labels,
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = "")
dev.off()



write.csv(MEs,file="results/GUT/MEs.csv")


modNames <- substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(expressiondata, MEs, use = "p",method = "spearman"));

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
modNames



names(geneModuleMembership) = paste("MM", modNames, sep="");

names(MMPvalue) = paste("p.MM", modNames, sep="");


write.csv(geneModuleMembership,file="results/GUT/ME_MEMBERSHIP.csv")

write.csv(MMPvalue,file="results/GUT/ME_MEMBERSHIP_PVAL.csv")


adj <- TOM

adj[adj > quantile(adj, 0.7)] = 1
adj[adj != 1] = 0

colnames(adj) <- colnames(expressiondata)
rownames(adj) <- colnames(expressiondata)

network <- graph_from_adjacency_matrix(adj,mode = "undirected")
network <- simplify(network)  # removes self-loops

V(network)$color <- mergedColors
# remove unconnected nodes
#network <- delete_vertices(network, degree(network)==0)

V(network)$name 

# Guardar la red en formato GraphML
write_graph(network, file = "results/GUT/network.graphml", format = "graphml")


rownames(MEs)<- rownames(traitData)


write.csv(MEs,"results/GUT/module_eigengenes.csv")

rownames(TOM)<- tax
colnames(TOM)<- tax


write.csv(TOM,"results/GUT/ADJ_MATRIX.csv")

tax_module <- cbind(tax,mergedColors)

write.csv(tax_module,"results/GUT/Tax_modules.csv")



# Cargar la matriz de similitud TOM y la clasificación en módulos
tax_module <- read.csv("results/GUT/Tax_modules.csv", stringsAsFactors = FALSE)
TOM <- as.matrix(read.csv("results/GUT/ADJ_MATRIX.csv", row.names = 1))



analizar_redes <- function(tax_module, TOM) {
  resultados <- data.frame()
  
  # ANÁLISIS DE LA RED COMPLETA
  g_all <- simplify(graph_from_adjacency_matrix(TOM, mode = "undirected", weighted = TRUE, diag = FALSE))
  

    # Dividir los nodos en tres grupos según su módulo
  modules <- unique(tax_module$mergedColors)
  lista_matrices <- list()
  
  # ANÁLISIS DE LOS SUB-GRAFOS (MÓDULOS)
  for (mod in modules) {
    taxa_in_mod <- tax_module$tax[tax_module$mergedColors == mod]
    lista_matrices[[mod]] <- TOM[taxa_in_mod, taxa_in_mod, drop = FALSE]
  }
  
  for (modulo in names(lista_matrices)) {
    TOM_submatrix <- lista_matrices[[modulo]]
    g_mod <- simplify(graph_from_adjacency_matrix(TOM_submatrix, mode = "undirected", weighted = TRUE, diag = FALSE))
    
    # Número de nodos en el subgrafo
    nodos_mod <- vcount(g_mod)
    
    # Identificar hubs como nodos en el percentil 95 de fuerza
    fuerza_nodos_mod <- strength(g_mod, mode = "all")
    fuerza_media_mod <- mean(fuerza_nodos_mod)
    
    # Treshold para hubs en el subgrafo
    umbral_hub_mod <- quantile(fuerza_nodos_mod, 0.95)
    hubs_mod <- names(fuerza_nodos_mod[fuerza_nodos_mod >= umbral_hub_mod])

    # Agregar los resultados del módulo al dataframe
    resultados <- rbind(resultados, data.frame(
      Modulo = modulo,
      Nodos = nodos_mod,
      Fuerza_Media = fuerza_media_mod,
      Hubs = paste(hubs_mod, collapse = ", ")
    ))
  }

  
  return(resultados)
}

df <- analizar_redes(tax_module, TOM)

# Mostrar tabla final
write.csv(df,"results/GUT/Module_analysis.csv")

View(df)
