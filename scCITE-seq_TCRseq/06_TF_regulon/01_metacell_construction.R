library(Seurat)
library(compiler)
library(Signac)
library(pheatmap)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(future)
library(tidyverse)
library(cowplot)
library(WGCNA)
library(hdWGCNA)
library(SeuratDisk)

# ---- paths (edit to your environment) --------------------------------------
analysis_dir <- "."                  # project root
objects_dir  <- "seurat_objects"     # Seurat / metacell objects
files_dir    <- "files"              # metadata tables
images_dir   <- "images"             # figure outputs

set.seed(2023)
options(future.globals.maxSize = 1000 * 1024^3)
plan("multisession", workers = 1)

setwd(analysis_dir)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# ============================================================================ #
# Load objects
# ============================================================================ #
tnk.noncar <- readRDS(file.path(objects_dir, "nonCAR_TNK_filtered.RDS"))
tnk.car    <- readRDS(file.path(objects_dir, "CAROnly_TNK_filtered.RDS"))
mono.object <- readRDS(file.path(objects_dir, "all_samples_Mono_filtered.RDS"))

# ============================================================================ #
# Call MetaCells
# ============================================================================ #

# ---- CAR T/NK --------------------------------------------------------------
# get variable genes
DefaultAssay(tnk.car) <- "RNA"
vegs <- VariableFeatures(tnk.car)

# check to make sure that at least 1% of cells express these genes.
gfreq <- rowMeans(tnk.car[["RNA"]]$counts > 0)
rgenes <- names(which(gfreq < 0.01))
vegs <- setdiff(vegs, rgenes)

# use VariableFeatures as VEGs
tnk.car <- SetupForWGCNA(
  tnk.car,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.01, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tnk_car_meta", # the name of the hdWGCNA experiment
  features = vegs
)

tnk.car$celltype_sample <- paste0(tnk.car$cell.anno, tnk.car$sample.name)

# make metacells
DefaultAssay(tnk.car) <- 'RNA'
Idents(tnk.car) <- "celltype_sample"
tnk_car_metacells_RNA <- MetacellsByGroups(
  seurat_obj = tnk.car,
  min_cells = 15,
  ident.group = "celltype_sample",
  group.by = 'celltype_sample',
  reduction = 'harmony.tnk',
  k = 15,
  max_shared = 3
)

metacell.tnk.car <- GetMetacellObject(tnk_car_metacells_RNA)
saveRDS(metacell.tnk.car, file.path(objects_dir, "Metacell_TNK_CAR.RDS"))

# ---- non-CAR T/NK ----------------------------------------------------------
# get variable genes
DefaultAssay(tnk.noncar) <- "RNA"
vegs <- VariableFeatures(tnk.noncar)

# check to make sure that at least 1% of cells express these genes.
gfreq <- rowMeans(tnk.noncar[["RNA"]]$counts > 0)
rgenes <- names(which(gfreq < 0.01))
vegs <- setdiff(vegs, rgenes)

# use VariableFeatures as VEGs
tnk.noncar <- SetupForWGCNA(
  tnk.noncar,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.01, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tnk_noncar_meta", # the name of the hdWGCNA experiment
  features = vegs
)

tnk.noncar$celltype_sample <- paste0(tnk.noncar$cell.anno, tnk.noncar$sample.name)

# make metacells
DefaultAssay(tnk.noncar) <- 'RNA'
Idents(tnk.noncar) <- "celltype_sample"
tnk_noncar_metacells_RNA <- MetacellsByGroups(
  seurat_obj = tnk.noncar,
  min_cells = 15,
  ident.group = "celltype_sample",
  group.by = 'celltype_sample',
  reduction = 'harmony.tnk',
  k = 15,
  max_shared = 3
)

metacell.tnk.noncar <- GetMetacellObject(tnk_noncar_metacells_RNA)
saveRDS(metacell.tnk.noncar, file.path(objects_dir, "Metacell_TNK_NonCAR.RDS"))

# ---- Monocytes -------------------------------------------------------------
# get variable genes
DefaultAssay(mono.object) <- "RNA"
vegs <- VariableFeatures(mono.object)

# check to make sure that at least 1% of cells express these genes.
gfreq <- rowMeans(mono.object[["RNA"]]$counts > 0)
rgenes <- names(which(gfreq < 0.01))
vegs <- setdiff(vegs, rgenes)

# use VariableFeatures as VEGs
mono.object <- SetupForWGCNA(
  mono.object,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.01, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "mono_meta", # the name of the hdWGCNA experiment
  features = vegs
)

mono.object$celltype_sample <- paste0(mono.object$cell.anno, mono.object$sample.name)

# make metacells
DefaultAssay(mono.object) <- 'RNA'
Idents(mono.object) <- "celltype_sample"
mono_metacells_RNA <- MetacellsByGroups(
  seurat_obj = mono.object,
  min_cells = 15,
  ident.group = "celltype_sample",
  group.by = 'celltype_sample',
  reduction = 'harmony.mono',
  k = 15,
  max_shared = 3
)

metacell.mono.object <- GetMetacellObject(mono_metacells_RNA)
saveRDS(metacell.mono.object, file.path(objects_dir, "Metacell_mono.RDS"))

# ============================================================================ #
# Process metacells (normalize / scale / PCA / Harmony / UMAP)
# ============================================================================ #

# ---- CAR T/NK --------------------------------------------------------------
tnk.car_obj <- readRDS(file.path(objects_dir, "Metacell_TNK_CAR.RDS"))

tnk.car_obj <- NormalizeMetacells(tnk.car_obj)
tnk.car_obj <- ScaleMetacells(tnk.car_obj, features = VariableFeatures(tnk.car_obj))
tnk.car_obj <- RunPCAMetacells(tnk.car_obj, features = VariableFeatures(tnk.car_obj))
tnk.car_obj <- RunHarmonyMetacells(tnk.car_obj, group.by.vars = "orig.ident")
tnk.car_obj <- RunUMAPMetacells(tnk.car_obj, reduction = 'harmony', dims = 1:15)

# add celltype back
df <- as.data.frame(tnk.car_obj@misc$tnk_car_meta$wgcna_metacell_obj$celltype_sample)
colnames(df) <- c("sample")

# Extract celltype name (everything before the first digit)
df <- df %>%
  separate(sample, into = c("celltype1", "celltype2"), sep = " ", extra = "merge", fill = "right")

df <- df %>%
  mutate(split1 = str_extract(celltype1, "^[^0-9]+(?=[0-9])"),
         split2 = str_extract(celltype2, "^[^0-9]+(?=[0-9])"))

df <- df %>%
  mutate(
    new_celltype = if_else(
      split1 == "CD",
      paste0(celltype1, " ", split2),
      split1
    )
  )


saveRDS(tnk.car_obj, file.path(objects_dir, "Metacell_TNK_CAR.RDS"))

# ---- non-CAR T/NK ----------------------------------------------------------
tnk.noncar_obj <- readRDS(file.path(objects_dir, "Metacell_TNK_NonCAR.RDS"))

tnk.noncar_obj <- NormalizeMetacells(tnk.noncar_obj)
tnk.noncar_obj <- ScaleMetacells(tnk.noncar_obj, features = VariableFeatures(tnk.noncar_obj))
tnk.noncar_obj <- RunPCAMetacells(tnk.noncar_obj, features = VariableFeatures(tnk.noncar_obj))
tnk.noncar_obj <- RunHarmonyMetacells(tnk.noncar_obj, group.by.vars = "orig.ident")
tnk.noncar_obj <- RunUMAPMetacells(tnk.noncar_obj, reduction = 'harmony', dims = 1:15)

# add celltype back
df <- as.data.frame(tnk.noncar_obj@misc$tnk_noncar_meta$wgcna_metacell_obj$celltype_sample)
colnames(df) <- c("sample")

# Extract celltype name (everything before the first digit)
df <- df %>%
  separate(sample, into = c("celltype1", "celltype2"), sep = " ", extra = "merge", fill = "right")

df <- df %>%
  mutate(split1 = str_extract(celltype1, "^[^0-9]+(?=[0-9])"),
         split2 = str_extract(celltype2, "^[^0-9]+(?=[0-9])"))

df <- df %>%
  mutate(
    new_celltype = if_else(
      split1 == "CD",
      paste0(celltype1, " ", split2),
      split1
    )
  )

tnk.noncar_obj@misc$tnk_noncar_meta$wgcna_metacell_obj@meta.data$new_celltype <- df$new_celltype


saveRDS(tnk.noncar_obj, file.path(objects_dir, "Metacell_TNK_NonCAR.RDS"))

# ---- Monocytes -------------------------------------------------------------
mono_obj <- readRDS(file.path(objects_dir, "Metacell_mono.RDS"))

mono_obj <- NormalizeMetacells(mono_obj)
mono_obj <- ScaleMetacells(mono_obj, features = VariableFeatures(mono_obj))
mono_obj <- RunPCAMetacells(mono_obj, features = VariableFeatures(mono_obj))
mono_obj <- RunHarmonyMetacells(mono_obj, group.by.vars = "orig.ident")
mono_obj <- RunUMAPMetacells(mono_obj, reduction = 'harmony', dims = 1:15)

# add celltype back
df <- as.data.frame(mono_obj@misc$mono_meta$wgcna_metacell_obj$celltype_sample)
colnames(df) <- c("sample")

# Extract celltype name (everything before the first digit)
celltypes <- unique(mono_obj$cell.anno)

df <- df %>%
  separate(sample, into = c("new_celltype", "sample.name"), sep = "\\|")

mono_obj@misc$mono_meta$wgcna_metacell_obj@meta.data$new_celltype <- df$new_celltype


saveRDS(mono_obj, file.path(objects_dir, "Metacell_mono.RDS"))

# ============================================================================ #
# Convert metacells into AnnData (h5ad via h5Seurat)
# ============================================================================ #

# ---- CAR T/NK --------------------------------------------------------------
metacell_tnk_car <- GetMetacellObject(tnk.car_obj)
VariableFeatures(metacell_tnk_car) <- VariableFeatures(tnk.car_obj)

counts <- GetAssayData(metacell_tnk_car, assay = "RNA", layer = "counts")

meta.data <- metacell_tnk_car@meta.data
meta.data <-
  meta.data %>%
  mutate(
    sample_id = str_remove(celltype_sample, str_c(new_celltype, collapse = "|"))
  )

# Select umap to add
umap.red <- metacell_tnk_car[["umap"]]
harmony.red <- metacell_tnk_car[["harmony"]]

# Create v3 seurat object
options(Seurat.object.assay.version = "v3")

dataset.to.convert <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = meta.data)

# Add dimensional reduction
dataset.to.convert$umap <- umap.red
dataset.to.convert$harmony <- harmony.red

# Factor to character, or else factors become numbers in adata
i <- sapply(dataset.to.convert@meta.data, is.factor)
dataset.to.convert@meta.data[i] <- lapply(dataset.to.convert@meta.data[i], as.character)

# Save h5Seurat
SaveH5Seurat(dataset.to.convert, filename = file.path(objects_dir, "Metacell_TNK_CAR.h5Seurat"), overwrite = TRUE)

# Convert to h5ad
Convert(file.path(objects_dir, "Metacell_TNK_CAR.h5Seurat"), dest = "h5ad", assay = "RNA", overwrite = TRUE)

# ---- non-CAR T/NK ----------------------------------------------------------
tnk.noncar_obj <- readRDS(file.path(objects_dir, "metacell_seurat_tnk_noncar.RDS"))
metacell_tnk.noncar <- GetMetacellObject(tnk.noncar_obj)
VariableFeatures(metacell_tnk.noncar) <- VariableFeatures(tnk.noncar_obj)

counts <- GetAssayData(metacell_tnk.noncar, assay = "RNA", layer = "counts")

meta.data <- metacell_tnk.noncar@meta.data
meta.data <-
  meta.data %>%
  mutate(
    sample_id = str_remove(celltype_sample, str_c(new_celltype, collapse = "|"))
  )

# Select umap to add
umap.red <- metacell_tnk.noncar[["umap"]]
harmony.red <- metacell_tnk.noncar[["harmony"]]

# Create v3 seurat object
options(Seurat.object.assay.version = "v3")

dataset.to.convert <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = meta.data)

# Add dimensional reduction
dataset.to.convert$umap <- umap.red
dataset.to.convert$harmony <- harmony.red

# Factor to character, or else factors become numbers in adata
i <- sapply(dataset.to.convert@meta.data, is.factor)
dataset.to.convert@meta.data[i] <- lapply(dataset.to.convert@meta.data[i], as.character)

# Save h5Seurat
SaveH5Seurat(dataset.to.convert, filename = file.path(objects_dir, "Metacell_TNK_NonCAR.h5Seurat"), overwrite = TRUE)

# Convert to h5ad
Convert(file.path(objects_dir, "Metacell_TNK_NonCAR.h5Seurat"), dest = "h5ad", assay = "RNA", overwrite = TRUE)

# ---- Monocytes -------------------------------------------------------------
metacell_mono <- GetMetacellObject(mono_obj)
VariableFeatures(metacell_mono) <- VariableFeatures(mono_obj)

counts <- GetAssayData(metacell_mono, assay = "RNA", layer = "counts")

meta.data <- metacell_mono@meta.data

# Select umap to add
umap.red <- metacell_mono[["umap"]]
harmony.red <- metacell_mono[["harmony"]]

# Create v3 seurat object
options(Seurat.object.assay.version = "v3")

dataset.to.convert <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = meta.data)

# Add dimensional reduction
dataset.to.convert$umap <- umap.red
dataset.to.convert$harmony <- harmony.red

# Factor to character, or else factors become numbers in adata
i <- sapply(dataset.to.convert@meta.data, is.factor)
dataset.to.convert@meta.data[i] <- lapply(dataset.to.convert@meta.data[i], as.character)

# Save h5Seurat
SaveH5Seurat(dataset.to.convert, filename = file.path(objects_dir, "Metacell_mono.h5Seurat"), overwrite = TRUE)

# Convert to h5ad
Convert(file.path(objects_dir, "Metacell_mono.h5Seurat"), dest = "h5ad", assay = "RNA", overwrite = TRUE)