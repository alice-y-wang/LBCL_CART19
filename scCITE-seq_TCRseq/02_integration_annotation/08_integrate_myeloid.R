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

# Load Monocyte Object # -------
merged.mono      <- readRDS("seurat_objects/all_samples_MONO_T-Removed.RDS")
DefaultAssay(merged.mono) <- "RNA"
message("Cells in monocyte object: ", length(Cells(merged.mono)))

# filter for less expressed genes in cells 
merged.mono$tmp.sample_name <- merged.mono$sample.name

mtx   <- merged.mono[["RNA"]]$counts
sIDs  <- merged.mono$tmp.sample_name

# Fraction of cells expressing each gene, per sample
mtx_pbulk <- sapply(unique(sIDs), function(x) rowMeans(mtx[, sIDs == x] > 0))
rmaxs     <- apply(mtx_pbulk, 1, max)

# variable feature selection # ------
merged.mono[["RNA"]] <- split(merged.mono[["RNA"]], f = merged.mono$tmp.sample_name)
merged.mono          <- NormalizeData(merged.mono, normalization.method = "LogNormalize", scale.factor = 10000)
merged.mono          <- FindVariableFeatures(merged.mono, selection.method = "vst", nfeatures = 3500)

var_features <- VariableFeatures(merged.mono)

# Remove TCR alpha/beta variable genes
var_features <- var_features[!grepl("^TRBV", var_features)]
var_features <- var_features[!grepl("^TRAV", var_features)]

# Remove features expressed in <1% of cells in every sample
rare_features <- names(which(rmaxs < 0.01))
var_features  <- intersect(setdiff(rownames(mtx), rare_features), var_features)

# Remove CAR transgene
var_features <- setdiff(var_features, "scfv")

message("Final variable feature count: ", length(var_features))
VariableFeatures(merged.mono) <- var_features

# Unintegrated Scale, PCA, UMAP # ------
merged.mono <- ScaleData(
  merged.mono,
  vars.to.regress = c("G2M.Score", "S.Score", "nCount_RNA", "percent.mt", "HeatShock.Score1")
)

merged.mono <- RunPCA(merged.mono,      npcs = 50, reduction.name = "unint.mono.pca")
merged.mono <- FindNeighbors(merged.mono, dims = 1:50, reduction = "unint.mono.pca",
                             graph.name = c("unint.mono.snn", "unint.mono.nn"))
merged.mono <- FindClusters(merged.mono,  resolution = 1, graph.name = "unint.mono.snn")
merged.mono <- RunUMAP(merged.mono,      dims = 1:50, reduction = "unint.mono.pca",
                       reduction.name = "unint.mono.umap")

# save
saveRDS(merged.mono, "seurat_objects/all_samples_MONO_T-Removed_unintPCA_regressCC_3500VarFeat_50dim.RDS")

# Harmony Integration # ------

merged.obj.check <- readRDS("seurat_objects/all_samples_MONO_T-Removed_unintPCA_regressCC_3500VarFeat_50dim.RDS")

merged.obj.check <- IntegrateLayers(
  object = merged.obj.check,
  method = HarmonyIntegration,
  orig.reduction = "unint.mono.pca",
  dims.use = 1:50,
  new.reduction = "harmony.mono",
  theta = 0.2,
  lambda = NULL,
  npcs=50L,
  verbose = TRUE
)

# re-join layers after integration
merged.obj.check[["RNA"]] <- JoinLayers(merged.obj.check[["RNA"]])

merged.obj.check <- FindNeighbors(merged.obj.check, reduction = "harmony.mono", dims = 1:50, graph.name = c("harmony.nn", "harmony.snn"), prune.SNN = 1/25)
merged.obj.check <- FindClusters(merged.obj.check, resolution = c(0.8), graph.name = "harmony.snn")
merged.obj.check <- RunUMAP(merged.obj.check, dims = 1:50, reduction = "harmony.mono", reduction.name = "harmony.mono.umap",
                            min.dist = 0.1, n.neighbors = 100L, return.model=TRUE)

saveRDS(merged.obj.check, "seurat_objects/all_samples_Mono_T-Removed_Harmony.RDS")

