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


# LOAD T/NK SUBSET # ------

merged.tnk       <- readRDS("seurat_objects/all_samples_TNK_mono_dying_removed.RDS")
DefaultAssay(merged.tnk) <- "RNA"
message("Cells in T/NK object: ", length(Cells(merged.tnk)))


# PSEUDOBULK FILTER — REMOVE RARE FEATURES (<1% cells) #------

mtx   <- merged.tnk[["RNA"]]$counts
sIDs  <- merged.tnk$sample.name

# Fraction of cells expressing each gene, per sample
mtx_pbulk <- sapply(unique(sIDs), function(x) rowMeans(mtx[, sIDs == x] > 0))
rmaxs     <- apply(mtx_pbulk, 1, max)


# VARIABLE FEATURE SELECTION # -----


merged.tnk[["RNA"]] <- split(merged.tnk[["RNA"]], f = merged.tnk$sample.name)
merged.tnk          <- NormalizeData(merged.tnk, normalization.method = "LogNormalize", scale.factor = 10000)
merged.tnk          <- FindVariableFeatures(merged.tnk, selection.method = "vst", nfeatures = 3500)

var_features <- VariableFeatures(merged.tnk)

# Remove TCR alpha/beta variable genes
var_features <- var_features[!grepl("^TRBV", var_features)]
var_features <- var_features[!grepl("^TRAV", var_features)]

# Remove features expressed in <1% of cells in every sample
rare_features <- names(which(rmaxs < 0.01))
var_features  <- intersect(setdiff(rownames(mtx), rare_features), var_features)

# Remove CAR transgene
var_features <- setdiff(var_features, "scfv")

message("Final variable feature count: ", length(var_features))
VariableFeatures(merged.tnk) <- var_features


# SCALE, PCA, NEIGHBORS (unintegrated) # ------

merged.tnk <- ScaleData(
  merged.tnk,
  vars.to.regress = c("G2M.Score", "S.Score", "nCount_RNA", "percent.mt", "HeatShock.Score1")
)

merged.tnk <- RunPCA(merged.tnk,      npcs = 50, reduction.name = "unint.tnk.pca")
merged.tnk <- FindNeighbors(merged.tnk, dims = 1:50, reduction = "unint.tnk.pca",
                            graph.name = c("unint.tnk.nn", "unint.tnk.snn"))

# save 
saveRDS(merged.tnk, "seurat_objects/all_samples_TNK_mono_dying_removed_3500varFeat_unintPCA.RDS")


# HARMONY INTEGRATION run Harmony with T cells # -----

merged.obj.check <- readRDS("seurat_objects/all_samples_TNK_mono_dying_removed_3500varFeat_unintPCA.RDS")

merged.obj.check <- IntegrateLayers(
  object = merged.obj.check,
  method = HarmonyIntegration,
  orig.reduction = "unint.tnk.pca",
  dims.use = 1:50, 
  new.reduction = "harmony.tnk",
  theta = 0.2,
  lambda = NULL,
  npcs=50L,
  verbose = TRUE
)

# re-join layers after integration
merged.obj.check[["RNA"]] <- JoinLayers(merged.obj.check[["RNA"]])

merged.obj.check <- FindNeighbors(merged.obj.check, reduction = "harmony.tnk", dims = 1:50, graph.name = c("harmony.nn", "harmony.snn"), prune.SNN = 1/25) # use 50 dims for neighbors
merged.obj.check <- FindClusters(merged.obj.check, resolution = c(0.1, 0.2, 0.3, 0.4, 0.6), graph.name = "harmony.snn")
merged.obj.check <- RunUMAP(merged.obj.check, dims = 1:50, min.dist=0.1, n.neighbors = 100L, reduction = "harmony.tnk", reduction.name = "harmony.tnk.umap", return.model=TRUE)

saveRDS(merged.obj.check, "seurat_objects/all_samples_TNK_mono_dying_removed_Harmony.RDS")