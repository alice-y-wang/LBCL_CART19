################################################################################
## hdWGCNA metacell construction (CAR+ T/NK, CAR- T/NK, monocytes)
################################################################################

# ---- Libraries --------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(future)
library(WGCNA)
library(hdWGCNA)

# ---- CONFIG -----------------------------------------------------------------
analysis_dir <- "."                  # project root
setwd(analysis_dir)

objects_dir  <- "seurat_objects"

dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)

# Input objects
f_tnk_noncar <- "nonCAR_TNK_filtered.RDS"
f_tnk_car    <- "CAROnly_TNK_filtered.RDS"
f_mono       <- "all_samples_Mono_filtered.RDS"

# Output metacell objects
f_meta_tnk_car    <- "metacell_seurat_tnk_car.RDS"
f_meta_tnk_noncar <- "metacell_seurat_tnk_noncar.RDS"
f_meta_mono       <- "metacell_seurat_mono.RDS"

# AnnData exports
f_h5_tnk_car    <- "Metacell_TNK_CAR.h5Seurat"
f_h5_tnk_noncar <- "Metacell_TNK_NonCAR.h5Seurat"
f_h5_mono       <- "Metacell_Mono.h5Seurat"

# Minimum fraction of cells a gene must be expressed in to be kept
min_gene_fraction <- 0.01

source(utils_script)
set.seed(2023)
options(future.globals.maxSize = 1000 * 1024^3)
plan("multisession", workers = 8)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# ---- Load objects -----------------------------------------------------------
tnk.noncar  <- readRDS(file.path(objects_dir, f_tnk_noncar))
tnk.car     <- readRDS(file.path(objects_dir, f_tnk_car))
mono.object <- readRDS(file.path(objects_dir, f_mono))


# ============================================================================ #
# CAR+ T/NK
# ============================================================================ #
DefaultAssay(tnk.car) <- "RNA"
vegs <- VariableFeatures(tnk.car)

# keep only genes expressed in at least min_gene_fraction of cells
gfreq <- rowMeans(tnk.car[["RNA"]]$counts > 0)
rgenes <- names(which(gfreq < min_gene_fraction))
vegs <- setdiff(vegs, rgenes)

# use the filtered variable features as VEGs
tnk.car <- SetupForWGCNA(
  tnk.car,
  gene_select = "fraction",         # the gene selection approach
  fraction = min_gene_fraction,     # fraction of cells a gene must be expressed in
  wgcna_name = "tnk_car_meta",      # the name of the hdWGCNA experiment
  features = vegs
)

tnk.car$celltype_sample <- paste0(tnk.car$cell.anno, "|", tnk.car$sample.name)

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

saveRDS(tnk_car_metacells_RNA, file.path(objects_dir, f_meta_tnk_car))

# ============================================================================ #
# CAR- T/NK
# ============================================================================ #
DefaultAssay(tnk.noncar) <- "RNA"
vegs <- VariableFeatures(tnk.noncar)

# keep only genes expressed in at least min_gene_fraction of cells
gfreq <- rowMeans(tnk.noncar[["RNA"]]$counts > 0)
rgenes <- names(which(gfreq < min_gene_fraction))
vegs <- setdiff(vegs, rgenes)

# use the filtered variable features as VEGs
tnk.noncar <- SetupForWGCNA(
  tnk.noncar,
  gene_select = "fraction",
  fraction = min_gene_fraction,
  wgcna_name = "tnk_noncar_meta",
  features = vegs
)

tnk.noncar$celltype_sample <- paste0(tnk.noncar$cell.anno, "|", tnk.noncar$sample.name)

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

saveRDS(tnk_noncar_metacells_RNA, file.path(objects_dir, f_meta_tnk_noncar))

# ============================================================================ #
# Monocytes
# ============================================================================ #
DefaultAssay(mono.object) <- "RNA"
vegs <- VariableFeatures(mono.object)

# keep only genes expressed in at least min_gene_fraction of cells
gfreq <- rowMeans(mono.object[["RNA"]]$counts > 0)
rgenes <- names(which(gfreq < min_gene_fraction))
vegs <- setdiff(vegs, rgenes)

# use the filtered variable features as VEGs
mono.object <- SetupForWGCNA(
  mono.object,
  gene_select = "fraction",
  fraction = min_gene_fraction,
  wgcna_name = "mono_meta",
  features = vegs
)

# clusters included for finer mapping of similar cells
mono.object$celltype_sample <- paste0(mono.object$cell.anno, "|",
                                      mono.object$harmony.snn_res.0.4, "|",
                                      mono.object$sample.name)

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

saveRDS(mono_metacells_RNA, file.path(objects_dir, f_meta_mono))

# ============================================================================ #
# Process metacells: CAR+ T/NK
# ============================================================================ #
tnk.car_obj <- readRDS(file.path(objects_dir, f_meta_tnk_car))

tnk.car_obj <- NormalizeMetacells(tnk.car_obj)
tnk.car_obj <- ScaleMetacells(tnk.car_obj, features=VariableFeatures(tnk.car_obj))

# add celltype back # -----
# celltype_sample is "<cell.anno>|<sample.name>"
df <- as.data.frame(tnk.car_obj@misc$tnk_car_meta$wgcna_metacell_obj$celltype_sample)
colnames(df) <- c("sample")

df <- df %>%
  separate(sample, into = c("cell.anno", "sample.name"), sep = "\\|", extra = "merge", fill = "right")

tnk.car_obj@misc$tnk_car_meta$wgcna_metacell_obj@meta.data$cell.anno <- df$cell.anno

saveRDS(tnk.car_obj, file.path(objects_dir, f_meta_tnk_car))

# ============================================================================ #
# Process metacells: CAR- T/NK
# ============================================================================ #
tnk.noncar_obj <- readRDS(file.path(objects_dir, f_meta_tnk_noncar))

tnk.noncar_obj <- NormalizeMetacells(tnk.noncar_obj)
tnk.noncar_obj <- ScaleMetacells(tnk.noncar_obj, features=VariableFeatures(tnk.noncar_obj))

# add celltype back # -----
df <- as.data.frame(tnk.noncar_obj@misc$tnk_noncar_meta$wgcna_metacell_obj$celltype_sample)
colnames(df) <- c("sample")

df <- df %>%
  separate(sample, into = c("cell.anno", "sample.name"), sep = "\\|", extra = "merge", fill = "right")

tnk.noncar_obj@misc$tnk_noncar_meta$wgcna_metacell_obj@meta.data$cell.anno <- df$cell.anno

saveRDS(tnk.noncar_obj, file.path(objects_dir, f_meta_tnk_noncar))

# ============================================================================ #
# Process metacells: Monocytes
# ============================================================================ #
mono_obj <- readRDS(file.path(objects_dir, f_meta_mono))

mono_obj <- NormalizeMetacells(mono_obj)
mono_obj <- ScaleMetacells(mono_obj, features=VariableFeatures(mono_obj))


# add celltype back # -----
# celltype_sample is "<cell.anno>|<cluster>|<sample.name>"
df <- as.data.frame(mono_obj@misc$mono_meta$wgcna_metacell_obj$celltype_sample)
colnames(df) <- c("sample")

df <- df %>%
  separate(sample, into = c("cell.anno", "cluster", "sample.name"), sep = "\\|", extra = "merge", fill = "right")

mono_obj@misc$mono_meta$wgcna_metacell_obj@meta.data$cell.anno <- df$cell.anno

saveRDS(mono_obj, file.path(objects_dir, f_meta_mono))

# ============================================================================ #
# Convert into anndata object: CAR+ T/NK
# ============================================================================ #
metacell_tnk_car <-  GetMetacellObject(tnk.car_obj)
VariableFeatures(metacell_tnk_car) <- VariableFeatures(tnk.car_obj)

counts <- GetAssayData(metacell_tnk_car, assay = "RNA", layer = "counts")

meta.data <- metacell_tnk_car@meta.data
meta.data <-
  meta.data %>%
  mutate(
    sample_id = str_split_fixed(celltype_sample, "\\|", 2)[, 2]
  )

# Select umap to add
umap.red <- metacell_tnk_car[["umap"]]
harmony.red <- metacell_tnk_car[["harmony"]]

# Create v3 seurat object
options(Seurat.object.assay.version = "v3")

dataset.to.convert <- CreateSeuratObject(counts = counts, assay="RNA", meta.data = meta.data)

# Add dimensional reduction ([[<- assigns a reduction; $<- writes a metadata column)
dataset.to.convert[["umap"]] <- umap.red
dataset.to.convert[["harmony"]] <- harmony.red

# Factor to character, or else your factor will be number in adata
i <- sapply(dataset.to.convert@meta.data, is.factor)
dataset.to.convert@meta.data[i] <- lapply(dataset.to.convert@meta.data[i], as.character)

# Save h5Seurat
SaveH5Seurat(dataset.to.convert, filename=file.path(objects_dir, f_h5_tnk_car), overwrite = TRUE)

# Convert to h5ad
Convert(file.path(objects_dir, f_h5_tnk_car), dest = "h5ad", assay="RNA", overwrite = TRUE)

# ============================================================================ #
# Convert into anndata object: CAR- T/NK
# ============================================================================ #
metacell_tnk.noncar <-  GetMetacellObject(tnk.noncar_obj)
VariableFeatures(metacell_tnk.noncar) <- VariableFeatures(tnk.noncar_obj)

counts <- GetAssayData(metacell_tnk.noncar, assay = "RNA", layer = "counts")

meta.data <- metacell_tnk.noncar@meta.data
meta.data <-
  meta.data %>%
  mutate(
    sample_id = str_split_fixed(celltype_sample, "\\|", 2)[, 2]
  )

# Select umap to add
umap.red <- metacell_tnk.noncar[["umap"]]
harmony.red <- metacell_tnk.noncar[["harmony"]]

# Create v3 seurat object
options(Seurat.object.assay.version = "v3")

dataset.to.convert <- CreateSeuratObject(counts = counts, assay="RNA", meta.data = meta.data)

# Add dimensional reduction
dataset.to.convert[["umap"]] <- umap.red
dataset.to.convert[["harmony"]] <- harmony.red

# Factor to character, or else your factor will be number in adata
i <- sapply(dataset.to.convert@meta.data, is.factor)
dataset.to.convert@meta.data[i] <- lapply(dataset.to.convert@meta.data[i], as.character)

# Save h5Seurat
SaveH5Seurat(dataset.to.convert, filename=file.path(objects_dir, f_h5_tnk_noncar), overwrite = TRUE)

# Convert to h5ad
Convert(file.path(objects_dir, f_h5_tnk_noncar), dest = "h5ad", assay="RNA", overwrite = TRUE)

# ============================================================================ #
# Convert into anndata object: Monocytes
# ============================================================================ #
metacell_mono <-  GetMetacellObject(mono_obj)
VariableFeatures(metacell_mono) <- VariableFeatures(mono_obj)

counts <- GetAssayData(metacell_mono, assay = "RNA", layer = "counts")

meta.data <- metacell_mono@meta.data
meta.data <-
  meta.data %>%
  mutate(
    sample_id = str_split_fixed(celltype_sample, "\\|", 3)[, 3]
  )

# Select umap to add
umap.red <- metacell_mono[["umap"]]
harmony.red <- metacell_mono[["harmony"]]

# Create v3 seurat object
options(Seurat.object.assay.version = "v3")

dataset.to.convert <- CreateSeuratObject(counts = counts, assay="RNA", meta.data = meta.data)

# Add dimensional reduction
dataset.to.convert[["umap"]] <- umap.red
dataset.to.convert[["harmony"]] <- harmony.red

# Factor to character, or else your factor will be number in adata
i <- sapply(dataset.to.convert@meta.data, is.factor)
dataset.to.convert@meta.data[i] <- lapply(dataset.to.convert@meta.data[i], as.character)
