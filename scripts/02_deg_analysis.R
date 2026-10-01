# ==============================================================================
# Script Name: 02_deg_analysis.R
# Description: Custom pseudobulk differential gene expression (DEG) analysis. 
#              Due to small sample size (N=2 per group), an effect-size 
#              filtering strategy (|log2FC| > 0.5) is employed to uncover 
#              biologically relevant transcriptional shifts in CD4+ and CD8+ cells.
# ==============================================================================

# 1. Load Required Libraries
library(Seurat)
library(tidyverse)

# 2. Load Preprocessed Object (if running sequentially, ensure t_cells exists)
if (!exists("t_cells")) {
  t_cells <- readRDS("results/t_cells_preprocessed.rds")
  cd4_cells <- subset(t_cells, subset = orig.ident == "CD4_T")
  cd8_cells <- subset(t_cells, subset = orig.ident == "CD8_T")
}

# ------------------------------------------------------------------------------
# 1. CD4+ T CELLS DEG CALCULATION
# ------------------------------------------------------------------------------
counts_cd4 <- GetAssayData(cd4_cells, layer = "counts")
meta_cd4 <- cd4_cells@meta.data

# Group sample indices
norm_samples_cd4  <- rownames(meta_cd4[meta_cd4$disease_stage == "Normal", ])
early_samples_cd4 <- rownames(meta_cd4[meta_cd4$disease_stage == "Early_AD", ])
late_samples_cd4  <- rownames(meta_cd4[meta_cd4$disease_stage == "Late_AD", ])

# Mean pseudobulk expression per stage
mean_norm_cd4  <- rowMeans(counts_cd4[, norm_samples_cd4])
mean_early_cd4 <- rowMeans(counts_cd4[, early_samples_cd4])
mean_late_cd4  <- rowMeans(counts_cd4[, late_samples_cd4])

# Calculate log2 Fold Change with pseudocount (+1)
log2FC_cd4_early <- log2((mean_early_cd4 + 1) / (mean_norm_cd4 + 1))
log2FC_cd4_late  <- log2((mean_late_cd4 + 1)  / (mean_norm_cd4 + 1))

# Create dataframes
deg_cd4_early <- data.frame(gene = names(log2FC_cd4_early), log2FoldChange = log2FC_cd4_early)
deg_cd4_late  <- data.frame(gene = names(log2FC_cd4_late),  log2FoldChange = log2FC_cd4_late)

# Filter for effect-size threshold (|log2FC| > 0.5)
deg_cd4_early_sig <- deg_cd4_early[abs(deg_cd4_early$log2FoldChange) > 0.5, ]
deg_cd4_late_sig  <- deg_cd4_late[abs(deg_cd4_late$log2FoldChange) > 0.5, ]

# ------------------------------------------------------------------------------
# 2. CD8+ T CELLS DEG CALCULATION
# ------------------------------------------------------------------------------
counts_cd8 <- GetAssayData(cd8_cells, layer = "counts")
meta_cd8 <- cd8_cells@meta.data

# Group sample indices
norm_samples_cd8  <- rownames(meta_cd8[meta_cd8$disease_stage == "Normal", ])
early_samples_cd8 <- rownames(meta_cd8[meta_cd8$disease_stage == "Early_AD", ])
late_samples_cd8  <- rownames(meta_cd8[meta_cd8$disease_stage == "Late_AD", ])

# Mean pseudobulk expression per stage
mean_norm_cd8  <- rowMeans(counts_cd8[, norm_samples_cd8])
mean_early_cd8 <- rowMeans(counts_cd8[, early_samples_cd8])
mean_late_cd8  <- rowMeans(counts_cd8[, late_samples_cd8])

# Calculate log2 Fold Change
log2FC_cd8_early <- log2((mean_early_cd8 + 1) / (mean_norm_cd8 + 1))
log2FC_cd8_late  <- log2((mean_late_cd8 + 1)  / (mean_norm_cd8 + 1))

# Create dataframes
deg_cd8_early <- data.frame(gene = names(log2FC_cd8_early), log2FoldChange = log2FC_cd8_early)
deg_cd8_late  <- data.frame(gene = names(log2FC_cd8_late),  log2FoldChange = log2FC_cd8_late)

# Filter for effect-size threshold
deg_cd8_early_sig <- deg_cd8_early[abs(deg_cd8_early$log2FoldChange) > 0.5, ]
deg_cd8_late_sig  <- deg_cd8_late[abs(deg_cd8_late$log2FoldChange) > 0.5, ]

# ------------------------------------------------------------------------------
# 3. PRINT SUMMARY & EXPORT CSV RESULTS
# ------------------------------------------------------------------------------
cat("\n=== DEG ANALYSIS SUMMARY (|log2FC| > 0.5) ===\n")
cat("CD4+ Early_AD vs Normal:", nrow(deg_cd4_early_sig), "genes\n")
cat("CD4+ Late_AD vs Normal :", nrow(deg_cd4_late_sig), "genes\n")
cat("CD8+ Early_AD vs Normal:", nrow(deg_cd8_early_sig), "genes\n")
cat("CD8+ Late_AD vs Normal :", nrow(deg_cd8_late_sig), "genes\n")

# Save tables into results/
write.csv(deg_cd4_early_sig, "results/DEG_CD4_Early_vs_Normal.csv", row.names = FALSE)
write.csv(deg_cd4_late_sig,  "results/DEG_CD4_Late_vs_Normal.csv",  row.names = FALSE)
write.csv(deg_cd8_early_sig, "results/DEG_CD8_Early_vs_Normal.csv", row.names = FALSE)
write.csv(deg_cd8_late_sig,  "results/DEG_CD8_Late_vs_Normal.csv",  row.names = FALSE)
