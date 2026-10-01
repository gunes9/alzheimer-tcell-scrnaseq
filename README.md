# alzheimer-pbmc-tcell-scrnaseq

# Single-Cell Pseudobulk DEG Analysis of T-Cell Subsets in Alzheimer's Disease (GSE168522)

This repository contains a reproducible bioinformatics pipeline for single-cell derived pseudobulk differential expression analysis (DEG) on CD4+ and CD8+ T-cell subsets across Alzheimer's Disease (AD) progression stages (Normal, Early AD, Late AD) using the GSE168522 dataset.

---

##  Overview & Key Findings

Due to the limited sample size (N=2 per group), standard FDR/p-value adjusted differential testing (e.g., standard `FindMarkers`) fails to capture subtle transcriptomic shifts, yielding zero significant genes. To bypass this statistical power bottleneck, this pipeline employs an effect-size filtering strategy ($\vert{}\log_2\text{FoldChange}\vert{} > 0.5$) to uncover true, literature-validated biological drivers.

Our analysis identified a distinct two-step temporal immune trajectory localized within peripheral T-cell subsets:

1. Early-Stage Antioxidant Defense Breakdown (Early AD):
   * Genes: `HBB`, `HBA2` (Down-regulated in Early AD)
  
2. Late-Stage Ribosomal Stress Response (Late AD):
   * Genes: `RPS4Y1` (Up-regulated in Late AD)


---

##  Pipeline Architecture

1. Preprocessing & Pseudobulk Aggregation: Aggregation of single-cell transcriptomic profiles into cell-type-specific pseudobulk matrices (CD4+ and CD8+ subsets).
2. Differential Expression Analysis: Application of thresholding ($\vert{}\log_2\text{FC}\vert{} > 0.5$) across chronological disease stages (`Normal` vs. `Early_AD`, `Early_AD` vs. `Late_AD`).
3. Visualization: Chronologically ordered heatmaps (`Normal` ,  `Early_AD`,  `Late_AD`) highlighting stage-specific gene dynamics.
4. Network & Pathway Integration: Functional enrichment and PPI network clustering via STRING DB and MCL algorithm. (In Progress)

---

## 📂 Repository Structure

```text
├── data/               # Input pseudobulk count matrices and metadata
├── scripts/            # R scripts for DEG analysis and visualization
├── results/            # Output heatmaps, DEG tables
├── README.md           # Project documentation
└── LICENSE             # Open-source license
