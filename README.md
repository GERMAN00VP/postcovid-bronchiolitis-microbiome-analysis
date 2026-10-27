# Respiratory and gut microbiota in post-COVID bronchiolitis and recurrent wheezing

Analysis code for *"Respiratory Microbiota Composition in Bronchiolitis and Its Association with Recurrent Wheezing in Infants"*. The pipeline regenerates every figure and supplementary table of the manuscript from the processed 16S rRNA data in `data/`.

## Run

From the repository root (R ≥ 4.6; package versions in `R/sessionInfo.txt`):

```bash
Rscript scripts/00_run_all.R
```

or run the scripts one by one, in order:

| Script | Content | Outputs |
|:---|:---|:---|
| `01_sample_qc.R` | Pre-specified sample QC (depth, anatomical-site concordance, sample mix-up), sample flow, rarefaction curves | Supplementary Figure 6, data for Supplementary Table 5 |
| `02_curation_and_composition.R` | Exclusions, rarefaction, clinical-virological groups, curated objects | Figure 1 |
| `03_nasopharyngeal_analysis.R` | NPA alpha/beta diversity (age-adjusted), case–control comparison restricted to the age range of the controls, ANCOM-BC2, wheezing risk models | Figures 2–4, data for Supplementary Tables 1 and 4 |
| `04_gut_analysis.R` | Same analyses for stool | Supplementary Figures 2–4, Supplementary Table 2 |
| `05_cross_niche_model.R` | Cross-niche wheezing model | Figure 5 |
| `06_niche_comparison_and_networks.R` | NPA vs. gut comparison, SpiecEasi networks, module–trait regressions | Supplementary Figure 1, Figure 6, Supplementary Figure 5, Supplementary Table 3 |
| `07_supplementary_tables.R` | Assembles the supplementary tables | `tables/supplementary/*.xlsx` |

Shared functions are in `R/functions.R`. Figures are written to `figures/main/` and `figures/supplementary/` (JPEG, named by figure number); intermediate results to `results/`.

Tables 1 and 2 (clinical characteristics) are produced by `notebooks/Tables_1_2_clinical_characteristics.ipynb` (Python 3.12: pandas, NumPy, SciPy, statsmodels, pyreadstat). Its input clinical database is confidential and is not included (`data/clinical/` is git-ignored).

---

Raw reads: NCBI SRA, BioProject [PRJNA1551058](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1551058).

---

Steps that draw random numbers (permutations, plot jitter, imputation noise, SpiecEasi, Louvain) run in a fixed order with fixed seeds; the age-restricted analysis runs inside `with_preserved_seed()` so that it does not alter the main results.