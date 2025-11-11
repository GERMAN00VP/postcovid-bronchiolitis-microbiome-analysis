# 🧬 Post-COVID Bronchiolitis Microbiome Analysis

This repository contains the scripts and input data used for the analyses in:

**_Post-COVID Changes in Infant Bronchiolitis Microbiota and Respiratory Outcomes_**  
Vallejo-Palma G., Alcolea S., García-García ML., *et al.* (2025, manuscript in review)

---

## 📂 Repository Structure

<details>
<summary><strong>Click to expand</strong></summary>

```
BQL_ANALYSIS/
│
├── abundance_analysis.R      # Differential abundance (ANCOM-BC2) and main figure generation
├── adonis.R                  # Beta diversity (PERMANOVA, BETADISP, ANOSIM)
├── phyloseqcreation.R        # Builds the phyloseq object
├── Main_analysis.ipynb         # Exploratory data analysis notebook
│
├── WGCNA_BQL_ANF.R           # WGCNA for respiratory microbiota
├── WGCNA_BQL_GUT.R           # WGCNA for gut microbiota
│
├── utils/
│   ├── Ancova_analysis.py
│   └── class_Stastics.py
│
└── data/
    ├── phyloseq_data/
    │   ├── otu_table.csv
    │   ├── tax_table.csv
    │   └── sample_data.csv
    └── phyloseq_object.rds
```

</details>

---

## 🚀 Analysis Workflow

<details>
<summary><strong>1. Create Phyloseq Object</strong></summary>

```bash
Rscript phyloseqcreation.R
```

This creates the file:

```
data/phyloseq_object.rds
```

</details>

<details>
<summary><strong>2. Run Core Analysis (Main Figures & ANCOM-BC2)</strong></summary>

```bash
Rscript abundance_analysis.R
```

Outputs include:
- Alpha and beta diversity results
- Differential abundance tables (ANCOM-BC2)
- Manuscript figure panels
</details>

<details>
<summary><strong>3. Run Optional Downstream Analyses</strong></summary>

| Script | Purpose |
|--------|---------|
| `adonis.R` | PERMANOVA / BETADISP / ANOSIM |
| `WGCNA_BQL_ANF.R` | Network analysis (respiratory microbiota) |
| `WGCNA_BQL_GUT.R` | Network analysis (gut microbiota) |

Outputs include:
- Network modules
- Trait-module correlation heatmaps
- GraphML network files for Cytoscape/Gephi visualization

</details>

---

## 🧱 Dependencies

<details>
<summary><strong>Software Requirements</strong></summary>

| Tool | Version |
|------|---------|
| R | ≥ 4.1 |
| Python | ≥ 3.8 |

**Key R Packages:**  
`phyloseq`, `vegan`, `ANCOMBC`, `WGCNA`, `igraph`, `ggplot2`, `dplyr`
</details>

---

## 🔬 Data Availability

Raw sequencing FASTQ files will be deposited to **NCBI SRA** after publication.  
Processed feature tables and metadata are included in this repository.

---

## 📖 Citation

If you use this repository, please cite:

```
Vallejo-Palma G., Alcolea S., García-García ML., et al.
Post-COVID Changes in Infant Bronchiolitis Microbiota and Respiratory Outcomes.
2025 (manuscript in review)
```
