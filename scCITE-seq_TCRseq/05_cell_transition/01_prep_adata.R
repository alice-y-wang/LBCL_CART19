################################################################################
## Convert Seurat objects to h5ad for downstream Python/scanpy analysis
##
## Assumes the .RDS objects already carry all harmonized metadata (including
## the canonical cell-type column `cell.anno`) and their Harmony/UMAP
## reductions.
################################################################################

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
library(rlist)
library(hdf5r)
library(SeuratDisk)

# ---- CONFIG -----------------------------------------------------------------
project_dir  <- "."
obj_dir      <- file.path(project_dir, "seurat_objects")

# Object filenames (shared across the pipeline)
f_tnk_all    <- "all_samples_TNK_filtered.RDS"
f_tnk_noncar <- "nonCAR_TNK_filtered.RDS"
f_tnk_car    <- "CAROnly_TNK_filtered.RDS"
f_mono       <- "all_samples_Mono_filtered.RDS"

setwd(project_dir)
set.seed(2024)

##### Conversion to adata #####
# load monocyte only # -----

mono.seurat <- readRDS(file.path(obj_dir, f_mono))

counts <- GetAssayData(mono.seurat, assay = "RNA", layer = "counts")
norm <- GetAssayData(mono.seurat, assay = "RNA", layer = "data")

meta.data <- mono.seurat@meta.data

# Select umap to add
umap.red <- mono.seurat[["harmony.mono.umap"]]
harmony.red <- mono.seurat[["harmony.mono"]]

# Create v3 seurat object
options(Seurat.object.assay.version = "v3")

dataset.to.convert <- CreateSeuratObject(counts = counts, assay="RNA", meta.data = meta.data)

# Add dimensional reduction ([[<- assigns a reduction; $<- writes a metadata column)
dataset.to.convert[["umap"]] <- umap.red
dataset.to.convert[["harmony"]] <- harmony.red

# Redo log normalization
#dataset.to.convert <- NormalizeData(dataset.to.convert, normalization.method = "LogNormalize")
dataset.to.convert[["RNA"]]$data <- norm

# Factor to character, or else your factor will be number in adata
i <- sapply(dataset.to.convert@meta.data, is.factor)
dataset.to.convert@meta.data[i] <- lapply(dataset.to.convert@meta.data[i], as.character)

# Save h5Seurat
SaveH5Seurat(dataset.to.convert, filename=file.path(obj_dir, "all_samples_Mono_toadata.h5Seurat"), overwrite = TRUE)

# Convert to h5ad
Convert(file.path(obj_dir, "all_samples_Mono_toadata.h5Seurat"), dest = "h5ad", assay="RNA", overwrite = TRUE)

# load TNK seurat obj # ------
tnk.seurat <- readRDS(file.path(obj_dir, f_tnk_all))

# SPLIT  TNK OBJ into CAR+ and CAR neg # -----

tnk.noncar <- subset(tnk.seurat, subset = CAR.exp == FALSE)

tnk.car <- subset(tnk.seurat, subset = CAR.exp == TRUE)

saveRDS(tnk.noncar, file.path(obj_dir, f_tnk_noncar))
saveRDS(tnk.car, file.path(obj_dir, f_tnk_car))

tnk.noncar <- readRDS(file.path(obj_dir, f_tnk_noncar))
tnk.car <- readRDS(file.path(obj_dir, f_tnk_car))

# Non-CAR T and NK # -------
counts <- GetAssayData(tnk.noncar, assay = "RNA", layer = "counts")
norm <- GetAssayData(tnk.noncar, assay = "RNA", layer = "data")

meta.data <- tnk.noncar@meta.data

# Select umap to add
umap.red <- tnk.noncar[["harmony.tnk.umap"]]
harmony.red <- tnk.noncar[["harmony.tnk"]]

# Create v3 seurat object
options(Seurat.object.assay.version = "v3")

dataset.to.convert <- CreateSeuratObject(counts = counts, assay="RNA", meta.data = meta.data)

# Add dimensional reduction
dataset.to.convert[["umap"]] <- umap.red
dataset.to.convert[["harmony"]] <- harmony.red

# Redo log normalization
#dataset.to.convert <- NormalizeData(dataset.to.convert, normalization.method = "LogNormalize")
dataset.to.convert[["RNA"]]$data <- norm

# Factor to character, or else your factor will be number in adata
i <- sapply(dataset.to.convert@meta.data, is.factor)
dataset.to.convert@meta.data[i] <- lapply(dataset.to.convert@meta.data[i], as.character)

# Save h5Seurat
SaveH5Seurat(dataset.to.convert, filename=file.path(obj_dir, "nonCAR_TNK_toadata.h5Seurat"), overwrite = TRUE)

# Convert to h5ad
Convert(file.path(obj_dir, "nonCAR_TNK_toadata.h5Seurat"), dest = "h5ad", assay="RNA", overwrite = TRUE)

# CARonly T and NK # ------
counts <- GetAssayData(tnk.car, assay = "RNA", layer = "counts")
norm <- GetAssayData(tnk.car, assay = "RNA", layer = "data")

meta.data <- tnk.car@meta.data

# Select umap to add
umap.red <- tnk.car[["harmony.tnk.umap"]]
harmony.red <- tnk.car[["harmony.tnk"]]

# Create v3 seurat object
options(Seurat.object.assay.version = "v3")

dataset.to.convert <- CreateSeuratObject(counts = counts, assay="RNA", meta.data = meta.data)

# Add dimensional reduction
dataset.to.convert[["umap"]] <- umap.red
dataset.to.convert[["harmony"]] <- harmony.red

# Redo log normalization
#dataset.to.convert <- NormalizeData(dataset.to.convert, normalization.method = "LogNormalize")
dataset.to.convert[["RNA"]]$data <- norm

# Factor to character, or else your factor will be number in adata
i <- sapply(dataset.to.convert@meta.data, is.factor)
dataset.to.convert@meta.data[i] <- lapply(dataset.to.convert@meta.data[i], as.character)

# Save h5Seurat
SaveH5Seurat(dataset.to.convert, filename=file.path(obj_dir, "CARonly_TNK_toadata.h5Seurat"), overwrite = TRUE)

# Convert to h5ad
Convert(file.path(obj_dir, "CARonly_TNK_toadata.h5Seurat"), dest = "h5ad", assay="RNA", overwrite = TRUE)