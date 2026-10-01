# ==============================================================================
# Script Name: 01_preprocessing_and_pca.R

# Description: Load pre-cleaned and annotated T-cell subsets (CD4+ and CD8+),
#              perform Seurat normalization, HVG selection, disease stage 
#              annotation, and subset-specific PCA dimensional reduction.
# ==============================================================================

# 1. Load Required Libraries
library(Seurat)
library(tidyverse)

# 2. Define File Paths (Adjust relative path if needed)
cd4_path <- "data/GSE168522_CD4_T_cell_data.txt.gz"
cd8_path <- "data/GSE168522_CD8_T_cell_data.txt.gz"

# 3. Read Raw Count Data
cd4_data <- read.table(gzfile(cd4_path), header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cd8_data <- read.table(gzfile(cd8_path), header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 4. Create Seurat Objects
cd4_seurat <- CreateSeuratObject(counts = cd4_data, project = "CD4_T")
cd8_seurat <- CreateSeuratObject(counts = cd8_data, project = "CD8_T")

# 5. Merge T-Cell Datasets
t_cells <- merge(cd4_seurat, y = cd8_seurat, add.cell.ids = c("CD4", "CD8"), project = "AD_T_cells")

# 6. Data Normalization and Feature Selection
# Note: Pre-filtering was completed by dataset authors; cell-level QC was skipped.
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2000)

# 7. Scale Data and Run Global PCA
all.genes <- rownames(t_cells)
t_cells <- ScaleData(t_cells, features = all.genes)
t_cells <- RunPCA(t_cells, features = VariableFeatures(object = t_cells), npcs = 10)

# 8. Add Disease Stage Metadata
# Sample Mapping: Normal (16, 17) | Early AD (23, 28) | Late AD (24, 27)
t_cells$disease_stage <- "Normal"
t_cells$disease_stage[grepl("23|28", colnames(t_cells))] <- "Early_AD"
t_cells$disease_stage[grepl("24|27", colnames(t_cells))] <- "Late_AD"

# 9. Subset CD4+ and CD8+ T Cells & Perform Cell-Type Specific PCA
cd4_cells <- subset(t_cells, subset = orig.ident == "CD4_T")
cd8_cells <- subset(t_cells, subset = orig.ident == "CD8_T")

# CD4+ T Cells PCA
cd4_cells <- RunPCA(cd4_cells, features = VariableFeatures(cd4_cells), npcs = 5)
pca_cd4_plot <- DimPlot(cd4_cells, reduction = "pca", group.by = "disease_stage", pt.size = 5) + 
  ggtitle("CD4+ T Cells: Disease Stages")

# CD8+ T Cells PCA
cd8_cells <- RunPCA(cd8_cells, features = VariableFeatures(cd8_cells), npcs = 5)
pca_cd8_plot <- DimPlot(cd8_cells, reduction = "pca", group.by = "disease_stage", pt.size = 5) + 
  ggtitle("CD8+ T Cells: Disease Stages")

# 10. Save Objects and PCA Visualizations
dir.create("results", showWarnings = FALSE)
ggsave("results/CD4_PCA_disease_stages.png", plot = pca_cd4_plot, width = 6, height = 5)
ggsave("results/CD8_PCA_disease_stages.png", plot = pca_cd8_plot, width = 6, height = 5)

saveRDS(t_cells, "results/t_cells_preprocessed.rds")
