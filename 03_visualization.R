# ==============================================================================
# Script Name: 03_visualization.R
# Description: Generates publication-grade heatmaps for identified AD biomarkers 
#              (HBB, HBA2, RPS4Y1) in CD4+ T cells across disease stages.
# ==============================================================================

# 1. Load Required Packages
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)
library(Seurat)

# 2. Load Preprocessed Data (if running independently)
if (!exists("cd4_cells")) {
  t_cells <- readRDS("results/t_cells_preprocessed.rds")
  cd4_cells <- subset(t_cells, subset = orig.ident == "CD4_T")
}

# 3. Define Target Biomarkers Derived from DEG Analysis
target_genes <- c("HBB", "HBA2", "RPS4Y1")

# 4. Extract Expression Matrix and Metadata for CD4+ T Cells
counts_cd4 <- GetAssayData(cd4_cells, layer = "counts")
meta_cd4   <- cd4_cells@meta.data

# Filter target gene expression matrix
mat_cd4 <- counts_cd4[target_genes, ]
colnames(mat_cd4) <- paste(meta_cd4$disease_stage, seq_len(ncol(mat_cd4)), sep = "_")

# 5. Sort Sample Columns Chronologically (Normal -> Early_AD -> Late_AD)
sample_order_cd4 <- order(factor(meta_cd4$disease_stage, levels = c("Normal", "Early_AD", "Late_AD")))
mat_cd4_sorted <- mat_cd4[, sample_order_cd4]

# Create column annotation dataframe
annotation_col_cd4 <- data.frame(Disease_Stage = meta_cd4$disease_stage[sample_order_cd4])
rownames(annotation_col_cd4) <- colnames(mat_cd4_sorted)

# Define publication-quality colors for disease stages
ann_colors <- list(
  Disease_Stage = c(Normal = "#4DAF4A", Early_AD = "#377EB8", Late_AD = "#E41A1C")
)

# 6. Render and Save Publication-Grade Heatmap (Row-wise Z-score scaling)
heatmap_file <- "results/CD4_key_biomarkers_heatmap.png"

png(filename = heatmap_file, width = 2000, height = 1500, res = 300)
pheatmap(mat_cd4_sorted,
         scale = "row",
         annotation_col = annotation_col_cd4,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,                  # Keep chronological stage ordering
         cluster_rows = TRUE,                   # Cluster similar expression dynamics
         show_colnames = TRUE,
         show_rownames = TRUE,
         main = "CD4+ T Cells: Key AD Biomarker Expression Profile",
         color = colorRampPalette(c("#313695", "#FFFFBF", "#A50026"))(100))
dev.off()

cat("\nHeatmap successfully saved to:", heatmap_file, "\n")
