library(irlba)
library(data.table)
library(RSpectra)
library(magrittr)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(viridis)
library(compiler)
library(readr)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(R.utils)
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
# PATIENT / SAMPLE METADATA TABLES
# Map each sample name to its patient ID, response status,
# and timepoint. Edit to match your cohort.
# ============================================================

# Each row: sample.name | patient.id | response | timepoint
sample.metadata <- tribble(
  ~sample.name,         ~patient.id,  ~response,       ~timepoint,
  "PT01-APH",           "PT01",       "responder",      "APH",
  "PT01-TDN",           "PT01",       "responder",      "TDN",
  "PT01-D0",            "PT01",       "responder",      "D0",
  "PT01-D7",            "PT01",       "responder",      "Peak",
  "PT01-4W",            "PT01",       "responder",      "4W",
  "PT02-APH",           "PT02",       "nonresponder",   "APH",
  "PT02-TDN",           "PT02",       "nonresponder",   "TDN",
  "PT02-D0",            "PT02",       "nonresponder",   "D0",
  "PT02-D7",            "PT02",       "nonresponder",   "Peak",
  "PT02-4W",            "PT02",       "nonresponder",   "4W"
  # Add remaining patients here following the same pattern
)

# Ordered timepoint labels for plotting
timepoint.order.map <- c(
  "APH"  = "T1-APH",
  "TDN"  = "T2-TDN",
  "D0"   = "T3-D0",
  "Peak" = "T4-Peak",
  "4W"   = "T5-4-Week"
)

# ============================================================
# STEP 1: LOAD MERGED OBJECT AND REMOVE DOUBLETS
# ============================================================

merged.obj <- readRDS("seurat_objects/all_samples_merged.RDS")
merged.obj <- subset(merged.obj, subset = Doublet_Singlet == "Singlet")

# ============================================================
# STEP 2: JOIN RNA LAYERS
# ============================================================

DefaultAssay(merged.obj) <- "RNA"
merged.obj[["RNA"]] <- JoinLayers(merged.obj[["RNA"]])
message("RNA layers joined.")

# ============================================================
# STEP 3: ADD METADATA (patient ID, response, timepoint)
# ============================================================

meta <- merged.obj@meta.data %>%
  rownames_to_column("cell.barcode") %>%
  left_join(sample.metadata, by = "sample.name") %>%
  mutate(timepoint_order = timepoint.order.map[timepoint]) %>%
  column_to_rownames("cell.barcode")

merged.obj <- AddMetaData(merged.obj, meta[, c("patient.id", "response", "timepoint", "timepoint_order")])

# ============================================================
# STEP 4: PSEUDOBULK FILTER — REMOVE RARE FEATURES (<2% cells)
# ============================================================

mtx  <- merged.obj[["RNA"]]$counts
sIDs <- merged.obj$sample.name

# Fraction of cells expressing each gene, per sample
mtx_pbulk <- sapply(unique(sIDs), function(x) rowMeans(mtx[, sIDs == x] > 0))
rmaxs     <- apply(mtx_pbulk, 1, max)

# ============================================================
# STEP 5: VARIABLE FEATURE SELECTION
# Split by sample, normalize, find variable features,
# then filter TCR genes, rare features, and CAR transgene
# ============================================================

merged.obj[["RNA"]] <- split(merged.obj[["RNA"]], f = merged.obj$sample.name)
merged.obj          <- NormalizeData(merged.obj, normalization.method = "LogNormalize", scale.factor = 10000)
merged.obj          <- FindVariableFeatures(merged.obj, selection.method = "vst", nfeatures = 2000)

var_features <- VariableFeatures(merged.obj)

# Remove TCR alpha/beta variable genes
var_features <- var_features[!grepl("^TRBV", var_features)]
var_features <- var_features[!grepl("^TRAV", var_features)]

# Remove features expressed in <2% of cells in every sample
rare_features <- names(which(rmaxs < 0.02))
var_features  <- intersect(setdiff(rownames(mtx), rare_features), var_features)

# Remove CAR transgene
var_features <- setdiff(var_features, "scfv")

message(paste0("Final variable feature count: ", length(var_features)))
VariableFeatures(merged.obj) <- var_features

# ============================================================
# STEP 6: SCALE, PCA, CLUSTERING, UMAP (unintegrated)
# ============================================================

merged.obj <- ScaleData(
  merged.obj,
  vars.to.regress = c("G2M.Score", "S.Score", "nCount_RNA", "percent.mt", "HeatShock.Score1")
)

merged.obj <- RunPCA(merged.obj,      npcs = 50, reduction.name = "unint.pca")
merged.obj <- FindNeighbors(merged.obj, dims = 1:30, reduction = "unint.pca", graph.name = "unint.snn")
merged.obj <- FindClusters(merged.obj,  resolution = 1,  graph.name = "unint.snn")
merged.obj <- RunUMAP(merged.obj,      dims = 1:50, reduction = "unint.pca", reduction.name = "unint.umap")

saveRDS(merged.obj, "seurat_objects/for_integration/all_samples_merged_unintPCA_regressCC_noDoublets.RDS")