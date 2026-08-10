library(irlba)
library(data.table)
library(RSpectra)
library(magrittr)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(viridis)
library(readr)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(limma)
library(tidyverse)
library(DoubletFinder)
library(plyr)
library(rlist)
library(hdf5r)

set.seed(2024)

# ============================================================
# USER-DEFINED PATHS
# Update these variables to match your data structure
# ============================================================
analysis.path <- "/path/to/analysis/output/"
setwd(analysis.path)

# ============================================================
# STEP 1: LOAD UNINTEGRATED MERGED OBJECT WITH UNINTEGRATED PCA
# ============================================================
merged.obj <- readRDS("seurat_objects/for_integration/all_samples_merged_unintPCA_regressCC_noDoublets.RDS")

# ============================================================
# STEP 2: HARMONY INTEGRATION
# ============================================================
merged.obj <- IntegrateLayers(
  object          = merged.obj,
  method          = HarmonyIntegration,
  orig.reduction  = "unint.pca",
  dims.use        = 1:30,
  new.reduction   = "harmony",
  theta           = 0.2,
  lambda          = NULL,
  npcs            = 30L,
  verbose         = TRUE
)

# Re-join layers after integration
merged.obj[["RNA"]] <- JoinLayers(merged.obj[["RNA"]])

# ============================================================
# STEP 3: CLUSTERING AND UMAP
# ============================================================
merged.obj <- FindNeighbors(merged.obj, reduction = "harmony", dims = 1:30, graph.name = "harmony.snn")
merged.obj <- FindClusters(merged.obj,  resolution = c(0.6, 0.8, 1, 1.5, 2),              graph.name = "harmony.snn")
merged.obj <- RunUMAP(merged.obj,       dims = 1:30, reduction = "harmony", reduction.name = "harmony.umap", return.model = TRUE)

# ============================================================
# STEP 4: SAVE
# ============================================================
saveRDS(merged.obj, "seurat_objects/all_samples_Harmony_theta0.2_integrated_regressCC_noDoublets.RDS")