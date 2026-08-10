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
# STEP 1: LOAD HARMONY-INTEGRATED OBJECT AND CLUSTER AT RES 0.4
# ============================================================
merged.harmony <- readRDS("seurat_objects/all_samples_Harmony_theta0.2_integrated_regressCC_noDoublets.RDS")
merged.harmony <- FindClusters(merged.harmony, resolution = c(0.4), graph.name = "harmony.snn")
# ============================================================
# STEP 2: ADD UMAP COORDINATES TO METADATA
# (used for coordinate-based annotation below)
# ============================================================
umap.coords    <- Embeddings(merged.harmony, reduction = "harmony.umap")
merged.harmony <- AddMetaData(merged.harmony, metadata = umap.coords)


# ============================================================
# STEP 3: BASIC CELL TYPE ANNOTATION
# Cluster assignments based on harmony.snn_res.0.4;
# UMAP coordinates used to resolve ambiguous clusters
# on the T cell / Monocyte boundary
# ============================================================
merged.harmony$basic.anno <- "Unassigned"

merged.harmony$basic.anno[merged.harmony$harmony.snn_res.0.4 %in% c(2, 3, 6, 9, 1)] <- "T Cell"
merged.harmony$basic.anno[merged.harmony$harmony.snn_res.0.4 %in% c(0, 4)]           <- "Monocyte"

# Resolve boundary clusters by UMAP coordinate
merged.harmony$basic.anno[merged.harmony$harmonyumap_1 < 3 &
                            merged.harmony$harmony.snn_res.0.4 %in% c(0, 4)]         <- "T Cell"
merged.harmony$basic.anno[merged.harmony$harmonyumap_1 > 3 &
                            merged.harmony$harmony.snn_res.0.4 %in% c(1, 2, 3, 6)]   <- "Monocyte"

# TDN timepoint cells are exclusively T cells
merged.harmony$basic.anno[merged.harmony$timepoint == "TDN"]                          <- "T Cell"

merged.harmony$basic.anno[merged.harmony$harmony.snn_res.0.4 %in% c(5)]              <- "NK"
merged.harmony$basic.anno[merged.harmony$harmony.snn_res.0.4 %in% c(7)]              <- "Dendritic"
merged.harmony$basic.anno[merged.harmony$harmony.snn_res.0.4 %in% c(8)]              <- "B Cell"

# ============================================================
# STEP 4: VISUALIZE
# ============================================================
pdf("images/all_samples_Harmony_Basic_Anno.pdf")
print(DimPlot(merged.harmony, reduction = "harmony.umap", group.by = "basic.anno",
              raster = FALSE, label = TRUE) + ggtitle("Basic Annotation"))
dev.off()

# ============================================================
# STEP 5: SAVE ANNOTATED OBJECT
# ============================================================
saveRDS(merged.harmony, "seurat_objects/all_samples_Harmony_theta0.2_integrated_regressCC_noDoublets_BasicAnno.RDS")

# ============================================================
# STEP 6: SUBSET BY LINEAGE FOR DOWNSTREAM ANALYSIS
# ============================================================
tcell.NK.subset <- subset(merged.harmony, subset = basic.anno %in% c("T Cell", "NK"))
bcell.subset    <- subset(merged.harmony, subset = basic.anno %in% c("B Cell"))
mono.subset     <- subset(merged.harmony, subset = basic.anno %in% c("Monocyte", "Dendritic"))

saveRDS(tcell.NK.subset, "seurat_objects/all_samples_TNK_only.RDS")
saveRDS(bcell.subset,    "seurat_objects/all_samples_Bcell_only.RDS")
saveRDS(mono.subset,     "seurat_objects/all_samples_Mono_DC_only.RDS")



# ============================================================
# STEP 7: CONVERT INTO ANNDATA OBJECTS
# ============================================================


library(SeuratDisk)

tnk.convert <- tcell.NK.subset
tnk.convert[["ADT"]] <- NULL
tnk.convert[["prediction.score.celltype.l1"]] <- NULL
tnk.convert[["prediction.score.celltype.l2"]] <- NULL
tnk.convert[["prediction.score.celltype.l3"]] <- NULL

tnk.convert@reductions[["harmony"]] <- NULL
tnk.convert@reductions[["harmony.umap"]] <- NULL
tnk.convert@reductions[["unint.pca"]] <- NULL
tnk.convert@reductions[["unint.umap"]] <- NULL
tnk.convert@reductions[["unint.tnk.pca"]] <- NULL
tnk.convert@reductions[["unint.tnk.umap"]] <- NULL

tnk.convert[["RNA"]]$scale.data <- NULL

tnk.convert[["RNA3"]] <- as(object = tnk.convert[["RNA"]], Class = "Assay")
DefaultAssay(tnk.convert) <- "RNA3"
tnk.convert[["RNA"]] <- NULL

SaveH5Seurat(tnk.convert, filename = "seurat_objects/all_samples_TNKonly.h5Seurat")
Convert("seurat_objects/all_samples_TNKonly.h5Seurat", dest = "h5ad")

# Monocytes
mono.convert <- mono.subset
mono.convert[["ADT"]] <- NULL
mono.convert[["prediction.score.celltype.l1"]] <- NULL
mono.convert[["prediction.score.celltype.l2"]] <- NULL
mono.convert[["prediction.score.celltype.l3"]] <- NULL

mono.convert@reductions[["harmony"]] <- NULL
mono.convert@reductions[["harmony.umap"]] <- NULL
mono.convert@reductions[["unint.pca"]] <- NULL
mono.convert@reductions[["unint.umap"]] <- NULL
mono.convert@reductions[["unint.mono.pca"]] <- NULL
mono.convert@reductions[["unint.mono.umap"]] <- NULL

mono.convert[["RNA"]]$scale.data <- NULL

mono.convert[["RNA3"]] <- as(object = mono.convert[["RNA"]], Class = "Assay")
DefaultAssay(mono.convert) <- "RNA3"
mono.convert[["RNA"]] <- NULL


SaveH5Seurat(mono.convert, filename = "seurat_objects/all_samples_MONOonly.h5Seurat")
Convert("seurat_objects/all_samples_MONOonly.h5Seurat", dest = "h5ad")


