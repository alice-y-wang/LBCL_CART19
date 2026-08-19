
# ---- Libraries --------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(SeuratDisk)

# ---- CONFIG -----------------------------------------------------------------
analysis_dir <- "."
setwd(analysis_dir)

objects_dir  <- "seurat_objects"

dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)

# Input objects (shared across the pipeline)
f_full_dataset <- "LBCL_full_dataset_clean_new.RDS"
f_tnk_car      <- "CAROnly_TNK_filtered.RDS"
f_tnk_noncar   <- "nonCAR_TNK_filtered.RDS"

# AnnData exports
f_h5_car_myeloid    <- "CAR_myeloid_TDN_D0_LIANA.h5Seurat"
f_h5_noncar_myeloid <- "NonCAR_myeloid_TDN_D0_LIANA.h5Seurat"

# Cell types added to the T/NK compartment for the ligand-receptor analysis
myeloid_NK_celltypes <- c("CD14 Mono", "CD16 Mono", "Dendritic Cell", "NK")

set.seed(2024)

#### load objects ## ----
full.dataset <- readRDS(file.path(objects_dir, f_full_dataset))
tnk.car      <- readRDS(file.path(objects_dir, f_tnk_car))
tnk.noncar   <- readRDS(file.path(objects_dir, f_tnk_noncar))

# save the TDN-D0 CAR version # ----
CAR.cells    <- Cells(tnk.car)
nonCAR.cells <- Cells(tnk.noncar)

myeloid.seurat <- subset(full.dataset, subset = cell.anno %in% myeloid_NK_celltypes)
myeloid.cells  <- Cells(myeloid.seurat)

# combine cells
CAR_myeloid.cells    <- c(CAR.cells, myeloid.cells)
nonCAR_myeloid.cells <- c(nonCAR.cells, myeloid.cells)

# subset
CAR_myeloid.obj    <- subset(full.dataset, cells = CAR_myeloid.cells)
nonCAR_myeloid.obj <- subset(full.dataset, cells = nonCAR_myeloid.cells)

# add new timepoint
CAR_myeloid.obj$timepoint2 <- CAR_myeloid.obj$timepoint
CAR_myeloid.obj$timepoint2[CAR_myeloid.obj$timepoint %in% c("TDN", "D0")] <- "TDN-D0"

nonCAR_myeloid.obj$timepoint2 <- nonCAR_myeloid.obj$timepoint
nonCAR_myeloid.obj$timepoint2[nonCAR_myeloid.obj$timepoint %in% c("TDN", "D0")] <- "TDN-D0"

# CAR # -----
counts <- GetAssayData(CAR_myeloid.obj, assay = "RNA", layer = "counts")
norm <- GetAssayData(CAR_myeloid.obj, assay = "RNA", layer = "data")

meta.data <- CAR_myeloid.obj@meta.data

# Select umap to add
umap.red <- CAR_myeloid.obj[["harmony.umap"]]
harmony.red <- CAR_myeloid.obj[["harmony"]]

# Create v3 seurat object
options(Seurat.object.assay.version = "v3")

dataset.to.convert <- CreateSeuratObject(counts = counts, assay="RNA", meta.data = meta.data)

# Add dimensional reduction ([[<- assigns a reduction; $<- writes a metadata column)
dataset.to.convert[["umap"]] <- umap.red
dataset.to.convert[["harmony"]] <- harmony.red

# Redo log normalization
dataset.to.convert[["RNA"]]$data <- norm

# Factor to character, or else your factor will be number in adata
i <- sapply(dataset.to.convert@meta.data, is.factor)
dataset.to.convert@meta.data[i] <- lapply(dataset.to.convert@meta.data[i], as.character)

# Save h5Seurat
SaveH5Seurat(dataset.to.convert, filename=file.path(objects_dir, f_h5_car_myeloid), overwrite = TRUE)

# Convert to h5ad
Convert(file.path(objects_dir, f_h5_car_myeloid), dest = "h5ad", assay="RNA", overwrite = TRUE)

# non-CAR # -----
counts <- GetAssayData(nonCAR_myeloid.obj, assay = "RNA", layer = "counts")
norm <- GetAssayData(nonCAR_myeloid.obj, assay = "RNA", layer = "data")

meta.data <- nonCAR_myeloid.obj@meta.data

# Select umap to add
umap.red <- nonCAR_myeloid.obj[["harmony.umap"]]
harmony.red <- nonCAR_myeloid.obj[["harmony"]]

# Create v3 seurat object
options(Seurat.object.assay.version = "v3")

dataset.to.convert <- CreateSeuratObject(counts = counts, assay="RNA", meta.data = meta.data)

# Add dimensional reduction
dataset.to.convert[["umap"]] <- umap.red
dataset.to.convert[["harmony"]] <- harmony.red

# Redo log normalization
dataset.to.convert[["RNA"]]$data <- norm

# Factor to character, or else your factor will be number in adata
i <- sapply(dataset.to.convert@meta.data, is.factor)
dataset.to.convert@meta.data[i] <- lapply(dataset.to.convert@meta.data[i], as.character)

# Save h5Seurat
SaveH5Seurat(dataset.to.convert, filename=file.path(objects_dir, f_h5_noncar_myeloid), overwrite = TRUE)

# Convert to h5ad
Convert(file.path(objects_dir, f_h5_noncar_myeloid), dest = "h5ad", assay="RNA", overwrite = TRUE)
