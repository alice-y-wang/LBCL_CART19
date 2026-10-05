library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

BASE_DIR   <- "/path/to/your/"

REFERENCE_SAMPLE <- "DLBCL_34774" 
SKETCH_CELLS     <- 500000         
INTEGRATION_DIMS <- 1:30
CLUSTER_RES      <- c(0.25, 0.5, 0.75)


SEURAT_DIR <- file.path(BASE_DIR, "SeuratObj")
QC_DIR     <- file.path(BASE_DIR, "QC_image", "Integration")
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# Merged per-sample CODEX objects
merged_obj_path     <- file.path(SEURAT_DIR, "DLBCL_CODEX_Merged.rds")
integrated_obj_path <- file.path(SEURAT_DIR, "DLBCL_CODEX_Integrated.rds")


# -- Clustering marker panel ---------------------------------------------------
b_cell_markers      <- c("CD19", "CD20", "CD22", "CD79A")
t_cell_markers      <- c("CD3e", "CD4", "CD8", "CD5", "CD2", "TBET", "FOXP3",
                         "CD25", "CD27", "CD45RO", "CCR7")
myeloid_markers     <- c("CD16", "CD68", "CD14", "HLA-DR", "CD15", "CCR3")
dc_markers          <- c("CD11C", "CD11B", "CD123")
nk_markers          <- c("CD56")
stromal_markers     <- c("PDPN", "Vimentin", "Collagen IV", "CD34", "NAKATPase")

CLUSTERING_MARKERS <- unique(c(b_cell_markers, t_cell_markers, myeloid_markers,
                               dc_markers, nk_markers, stromal_markers))

# ==============================================================================
# 1. Load Merged Object and Split by Sample
# ==============================================================================

obj <- readRDS(merged_obj_path)
DefaultAssay(obj) <- "CODEX"
obj <- JoinLayers(obj)

counts_mat <- LayerData(obj, assay = "CODEX", layer = "counts")
saved_meta <- obj@meta.data
obj[["CODEX"]] <- CreateAssay5Object(counts = counts_mat)
obj@meta.data  <- saved_meta
rm(counts_mat); gc()

obj[["CODEX"]] <- split(obj[["CODEX"]], f = obj$orig.ident)

# ==============================================================================
# 2. CLR Normalization, Scaling, and PCA (all cells)
# ==============================================================================

obj <- NormalizeData(obj, normalization.method = "CLR", margin = 1)
VariableFeatures(obj) <- rownames(obj)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# ==============================================================================
# 3. Sketch and RPCA Integration
# ==============================================================================

obj <- SketchData(object = obj, assay = "CODEX", ncells = SKETCH_CELLS,
                  method = "LeverageScore", sketched.assay = "sketch",
                  over.write = TRUE)

DefaultAssay(obj) <- "sketch"
VariableFeatures(obj) <- rownames(obj)
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 1)
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, npcs = 30, features = rownames(obj),
              reduction.name = "unint.pca")

reference_idx <- which(Layers(obj, search = "data") == paste0("data.", REFERENCE_SAMPLE))

obj <- IntegrateLayers(obj, method = RPCAIntegration,
                       features = rownames(obj),
                       orig = "unint.pca", new.reduction = "integrated.rpca",
                       dims = INTEGRATION_DIMS, k.anchor = 20,
                       reference = reference_idx)

obj <- ProjectIntegration(object = obj, features = rownames(obj),
                          sketched.assay = "sketch", assay = "CODEX",
                          reduction = "integrated.rpca")

# ==============================================================================
# 4. Clustering on Sketch
# ==============================================================================

DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
obj <- ScaleData(obj)

clustering_features <- intersect(CLUSTERING_MARKERS, rownames(obj))
obj <- RunPCA(obj, features = clustering_features, approx = FALSE)
obj <- FindNeighbors(obj, features = clustering_features, k.param = 30)
for (res in CLUSTER_RES) {
  obj <- FindClusters(obj, algorithm = 2, resolution = res,
                      cluster.name = paste0("CODEX_snn_res.", res))
}

# ==============================================================================
# 5. UMAP on Sketch and Projection to All Cells
# ==============================================================================

obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = INTEGRATION_DIMS,
               return.model = TRUE)

obj[["sketch"]] <- split(obj[["sketch"]], f = obj$orig.ident)

obj <- ProjectData(object = obj, assay = "CODEX",
                   full.reduction = "integrated.rpca.full",
                   sketched.assay = "sketch",
                   sketched.reduction = "integrated.rpca",
                   umap.model = "umap", dims = INTEGRATION_DIMS,
                   refdata = list(CODEX_snn_res.0.5_full  = "CODEX_snn_res.0.5"))

DefaultAssay(obj) <- "CODEX"
saveRDS(obj, integrated_obj_path)

# ==============================================================================
# 7. QC Plots
# ==============================================================================

pdf(file.path(QC_DIR, "Integration_UMAP_Clusters_Sample.pdf"),
    width = 14, height = 12)
DimPlot(obj, reduction = "full.umap", group.by = "CODEX_snn_res.0.5_full",
        label = TRUE, repel = TRUE, raster = TRUE) + coord_fixed() + NoAxes()
DimPlot(obj, reduction = "full.umap", group.by = "orig.ident",
        raster = TRUE) + coord_fixed() + NoAxes()
dev.off()
