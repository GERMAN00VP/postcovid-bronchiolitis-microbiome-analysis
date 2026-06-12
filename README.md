# 🧬 Bronchiolitis Microbiome Analysis 

This repository contains the updated scripts, input data, and comprehensive outputs used for the analysis of respiratory and gut microbiota in infant bronchiolitis.

---

## 📂 Repository Structure

<details>
<summary><strong>Click to expand full tree structure</strong></summary>
```bash
.
├── data/
│   ├── phyloseq_data/
│   │   ├── otu_table.csv
│   │   ├── sample_data.csv
│   │   └── tax_table.csv
│   └── phyloseq_object.rds
│
├── REANALISIS/
│   ├── Global/             # Joint NPA vs GUT comparative analyses & multi-omics plots
│   │   ├── curated_data_global/
│   │   └── ESTADISTICA_CONJUNTA/
│   ├── GUT/                # Gut microbiota specific sub-pipeline (Alpha, Beta, ANCOM, Models)
│   │   ├── abundances/
│   │   ├── alpha/
│   │   ├── beta/
│   │   └── models/
│   ├── NPA/                # Nasopharyngeal aspirate specific sub-pipeline
│   │   ├── abundances/
│   │   ├── alpha/
│   │   ├── beta/
│   │   └── models/
│   └── NETWORKS/           # Cross-domain network interactions and clinical associations
│   
├── REANALISIS_initial_comparations.R  # Initial data parsing and global curation
├── REANALISIS_GUT.R                   # Core statistical pipeline for Gut samples
├── REANALISIS_NPA.R                   # Core statistical pipeline for NPA samples
└── REANALISIS_NETWORKS_INTERACTION.R  # Multi-domain network construction (igraph)

</details>
```
---

## 🚀 Reanalysis Workflow

The analysis pipeline is divided into four main execution blocks. All core statistical outputs, tables, and publication-ready figures are automatically exported to the REANALISIS/ directory.

### 1. Global Exploration & Curation
Initial data synchronization and comparative baseline statistics between Nasopharyngeal Aspirates (NPA) and Gut (GUT) samples.
* **Script:** REANALISIS_initial_comparations.R
* **Key Outputs (REANALISIS/Global/):**
    * Curated matrices: phyloseq_RAREFIED_global.rds, phyloseq_RAW_curated_global.rds
    * Main Panels: Plot_Alpha_Genus_Panel.png, Plot_Beta_Genus_Panel.png
    * Statistics: Joint LME ANOVA and Post-Hoc Tukey tables (ESTADISTICA_CONJUNTA/)

### 2. Gut & NPA Independent Pipelines
Parallel analysis tracks to evaluate alpha diversity, beta diversity dispersion, differential abundance (ANCOM-BC), and predictive clinical modeling (Univariate/Multivariate Logistic Regressions & ROC Curves).
* **Scripts:** REANALISIS_GUT.R and REANALISIS_NPA.R
* **Key Outputs (REANALISIS/GUT/ & REANALISIS/NPA/):**
    * alpha/: Pairwise Wilcoxon tests with FDR corrections + Clean Boxplots.
    * beta/: PERMANOVA tables and PCoA projections (Plot_Beta_Diversity_PCoA.png).
    * abundances/: Mega-Boxplots of global significant features via ANCOM.
    * models/: Univariate risk tables, definitive Multivariate Models, and Forest Plots.

### 3. Cross-Domain Network Interactions
Co-occurrence and association networks evaluating cross-domain relationships between the respiratory tract and the gut, correlated with clinical covariates.
* **Script:** REANALISIS_NETWORKS_INTERACTION.R
* **Key Outputs (REANALISIS/NETWORKS/):**
    * Network Objects: igraph_Genus_GUT.rds, igraph_Genus_NPA.rds
    * Visualizations: Dissected internal networks, simplified cross-domain maps, and expanded clinical heatmaps.
    * Data: Target Hub Genera lists and clinical association master tables.

---

## 🧱 Key Downstream Outputs

| Target Analysis | Main Figures | Main Statistical Tables |
|:---|:---|:---|
| **Global Sync** | Diagnostic_PCoA_Bray_Curtis_FINAL.png | PERMANOVA_BrayCurtis_Genus_Restringida.csv |
| **NPA / GUT Core** | Plot_MegaBoxplot_Minimal.png | Table_ANCOM_Global_Significant.csv |
| **Predictive Models** | Plot_ROC_Curve_Multivariate.png | Table_Multivariate_Model_Definitive.csv |
| **Networks** | Network_CrossDomain_SIMPLIFIED.png | Master_Taxa_Module_Clinical_Associations_FINAL.csv |

---

## 🔬 Dependencies & Data

* **R Environment:** Requires R >= 4.1
* **Core Packages:** phyloseq, vegan, ANCOMBC, igraph, ggplot2, dplyr, nlme, pROC
* **Input Data:** Processed OTU tables, taxonomic assignments, and sample metadata are fully contained within the data/phyloseq_data/ directory for immediate reproducibility.