library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

BASE_DIR   <- "/path/to/your/"
SCRIPT_DIR <- "/path/to/your/Scripts"
source(file.path(SCRIPT_DIR, "Utilities", "CODEX_Functions.R"))

SEURAT_DIR <- file.path(BASE_DIR, "SeuratObj")
QC_DIR     <- file.path(BASE_DIR, "QC_image", "Annotation_L1")
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

integrated_obj_path <- file.path(SEURAT_DIR, "DLBCL_CODEX_Integrated.rds")
annotated_obj_path  <- file.path(SEURAT_DIR, "DLBCL_CODEX_Annotated_L1.rds")

# -- Marker panels -------------------------------------------------------------
CLUSTERING_MARKERS <- c("CD19", "CD20", "CD22", "CD79A",
                        "CD3e", "CD4", "CD8", "CD5", "CD2", "TBET", "FOXP3",
                        "CD25", "CD27", "CD45RO", "CCR7",
                        "CD16", "CD68", "CD14", "HLA-DR", "CD15", "CCR3",
                        "CD11C", "CD11B", "CD123", "CD56",
                        "PDPN", "Vimentin", "Collagen IV", "CD34", "NAKATPase")
MYELOID_MARKERS <- c("CD80", "CD58", "PDL1", "CD47", "HVEM", "HLA-A",
                     "CD11C", "CD11B", "CD123", "CD16", "CD68", "CD14", "HLA-DR",
                     "CD163 Akoya", "CD15", "HIF1A", "Serpin B9", "KI67")
CD8T_MARKERS <- c("CD3e", "CD4", "CD8", "CD5", "CD2", "TBET", "FOXP3", "CD25",
                  "CD27", "CD45RO", "CCR7", "GATA3", "CXCR5", "CCR6", "CCR4", "CCR3",
                  "PD-1", "LAG3", "TIM3", "Vista", "CD57", "ICOS", "CD69", "GRZB",
                  "CXCR3", "TOX", "HLA-A", "KI67")
GRAN_CD4T_MARKERS <- c("CD19", "CD20", "CD22", "CD79A",
                       "CD3e", "CD4", "CD8", "CD5", "CD2", "TBET", "FOXP3", "CD25",
                       "CD27", "CD45RO", "CCR7", "GATA3",
                       "CD16", "CD68", "CD14", "HLA-DR", "CD163 Akoya", "CD15",
                       "CD11C", "CD11B", "CD123",
                       "CD80", "CD58", "PDL1", "CD47", "HVEM", "HLA-A", "KI67")

LINEAGE_MARKERS <- c("CD19", "CD20", "KI67", "CD3e", "CD4", "CD8", "FOXP3", "CD56",
                     "CD11C", "CD11B", "CD123", "CD68", "CD14", "CD16", "CD15", "HLA-DR",
                     "CD34", "Vimentin", "Collagen IV", "PDPN")

# Sub-cluster a set of cells on the integrated RPCA embedding
subcluster_ids <- function(obj, cells, res, round) {
  sub <- subset(obj, cells = cells)
  sub <- FindNeighbors(sub, reduction = "integrated.rpca.full", dims = 1:30)
  sub <- FindClusters(sub, resolution = res)
  pdf(file.path(QC_DIR, paste0(round, "_Subclusters_DotPlot.pdf")), width = 14, height = 8)
  print(DotPlot(sub, features = LINEAGE_MARKERS, scale = TRUE,
                col.min = -1.5, col.max = 1.5) + RotatedAxis() +
          scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
  dev.off()
  setNames(as.character(Idents(sub)), colnames(sub))
}

# Re-integrate one lineage and cluster it
lineage_cluster_ids <- function(obj, cells, markers, res, round) {
  sub <- subset(obj, cells = cells)
  sub <- reintegrate_sketch_rpca(sub, markers = markers, merge_small = FALSE)
  cluster_sketch_project(sub, markers = markers, res = res)
}

lineage_cluster_ids_nosketch <- function(obj, cells, markers, res, merge_samples, round) {
  reintegrate_cluster_nosketch(subset(obj, cells = cells), markers = markers,
                               res = res, merge_samples = merge_samples)
}

# Assign a label to cells
set_label <- function(ids, clusters, label) {
  cells <- names(ids)[ids %in% as.character(clusters)]
  L1[cells] <<- label
}

# ==============================================================================
# 1. Initial Annotation
# ==============================================================================

obj <- readRDS(integrated_obj_path)
DefaultAssay(obj) <- "CODEX"
obj <- JoinLayers(obj)
Idents(obj) <- "orig.ident"
obj@meta.data <- obj@meta.data[, !colnames(obj@meta.data) %in% rownames(obj[["CODEX"]])]

pdf(file.path(QC_DIR, "Integration_Clusters_DotPlot.pdf"), width = 14, height = 8)
DotPlot(obj, features = LINEAGE_MARKERS, group.by = "CODEX_snn_res.0.5_full", scale = TRUE,
        col.min = -1.5, col.max = 1.5) + RotatedAxis() +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

# Clusters annotated from their marker expression profiles
ids <- setNames(as.character(obj$CODEX_snn_res.0.5_full), colnames(obj))
L1  <- setNames(rep(NA_character_, ncol(obj)), colnames(obj))
set_label(ids, c(0, 9),               "CD8+ T-Cells")
set_label(ids, 1,                     "Dendritic Cells")
set_label(ids, c(3, 6, 7, 10, 11),    "Lymphoma")
set_label(ids, c(4, 13, 16),          "Monocytes")
set_label(ids, 5,                     "Endothelial")
set_label(ids, 12,                    "Myeloid")
set_label(ids, 14,                    "NK Cells")
set_label(ids, 15,                    "T-regulatory")
set_label(ids, c(2, 8, 17, 18, 19),   "Artifact")

# ==============================================================================
# 2. Sub-cluster Each Population to Resolve Ambiguous or Spill-Over Cells
# ==============================================================================

# Artifact
ids <- subcluster_ids(obj, names(L1)[L1 %in% "Artifact"], 0.4, "Step1_Artifact")
set_label(ids, 14,       "Dendritic Cells")
set_label(ids, c(17, 20), "Myeloid")
set_label(ids, 18,       "NK Cells")
ids <- subcluster_ids(obj, names(ids)[ids %in% c("5", "11", "13")], 0.3, "Step1_Artifact_sub")
set_label(ids, c(3, 4),  "CD8+ T-Cells")
set_label(ids, c(2, 5),  "Endothelial")

# Myeloid cells
ids <- subcluster_ids(obj, names(L1)[L1 %in% c("Monocytes", "Dendritic Cells", "Myeloid")],
                      0.4, "Step2_Myeloid")
set_label(ids, c(0, 1, 2, 4, 10, 11, 16), "Dendritic Cells")
set_label(ids, c(3, 9),                   "Monocytes")
set_label(ids, c(5, 6),                   "Monocyte/Macrophage")
set_label(ids, 12,                        "Macrophage")
set_label(ids, 14,                        "CD11B+ Myeloid")
ids <- subcluster_ids(obj, names(ids)[ids %in% c("7", "8", "13", "15", "17")], 0.4,
                      "Step2_Myeloid_sub")
set_label(ids, 0,                         "Monocytes")
set_label(ids, c(1, 2),                   "Granulocytes")
set_label(ids, 3,                         "CD8+ T-Cells")
set_label(ids, 10,                        "Lymphoma")
set_label(ids, c(4, 5, 6, 7, 8, 9, 11, 12, 13), "Artifact")

# Endothelial cells
ids <- subcluster_ids(obj, names(L1)[L1 %in% "Endothelial"], 0.4, "Step3_Endothelial")
set_label(ids, c(1, 6),  "Stromal")
set_label(ids, c(4, 11), "CD8+ T-Cells")
set_label(ids, 13,       "T-regulatory")
ids <- subcluster_ids(obj, names(ids)[ids %in% c("0", "3", "7", "12", "14")], 0.6,
                      "Step3_Endothelial_sub")
set_label(ids, c(2, 6, 10, 11),      "NK Cells")
set_label(ids, 5,                    "Stromal")
set_label(ids, 8,                    "CD4+ T-Cells")
set_label(ids, c(7, 12, 13, 14, 15), "Artifact")

# T and NK cells
ids <- subcluster_ids(obj, names(L1)[L1 %in% c("CD8+ T-Cells", "CD4+ T-Cells",
                                               "T-regulatory", "NK Cells")], 0.5, "Step4_T_NK")
set_label(ids, c(0, 3, 4, 5, 8, 9, 14, 15, 16, 18), "CD8+ T-Cells")
set_label(ids, c(6, 12), "NK Cells")
set_label(ids, 7,        "T-regulatory")
set_label(ids, 19,       "REMOVE")
ids <- subcluster_ids(obj, names(ids)[ids %in% c("1", "2", "10", "11", "13", "17", "20", "21")],
                      0.6, "Step4_T_NK_sub")
set_label(ids, c(0, 1, 2, 3, 4, 5, 6, 7, 12, 13), "CD4+ T-Cells")
set_label(ids, c(9, 10),     "CD8+ T-Cells")
set_label(ids, c(8, 11, 14), "NK Cells")
set_label(ids, c(15, 16, 17), "REMOVE")

# Lymphoma
lymphoma_cells <- names(L1)[L1 %in% "Lymphoma"]
ids  <- subcluster_ids(obj, lymphoma_cells, 0.4, "Step5_Lymphoma")
ids5 <- subcluster_ids(obj, lymphoma_cells, 0.5, "Step5_Lymphoma_res0.5")
set_label(ids, 0:8,          "Lymphoma")
set_label(ids, c(15, 16, 17), "REMOVE")
ids5 <- ids5[names(ids)[ids %in% as.character(9:14)]]
set_label(ids5, c(0, 5, 7, 8, 9),                   "Lymphoma")
set_label(ids5, c(1, 2, 3, 4, 6, 10, 11, 12, 13),   "Ambiguous")
set_label(ids5, 14,                                 "REMOVE")

# Ambiguous cells
ids <- subcluster_ids(obj, names(L1)[L1 %in% "Ambiguous"], 0.7, "Step6_Ambiguous")
set_label(ids, c(0, 9, 12, 14),    "REMOVE")
set_label(ids, c(1, 3, 6, 7, 10),  "Lymphoma")
set_label(ids, c(4, 11),           "Macrophage")
set_label(ids, 13,                 "CD8+ T-Cells")
L1[L1 %in% "Ambiguous"] <- "Artifact"

# Artifact cells
ids <- subcluster_ids(obj, names(L1)[L1 %in% "Artifact"], 0.5, "Step8_Artifact")
set_label(ids, c(3, 6, 9), "Dendritic Cells")
set_label(ids, c(11, 17),  "Granulocytes")
set_label(ids, 19,         "CD8+ T-Cells")

# NK cells
ids <- subcluster_ids(obj, names(L1)[L1 %in% "NK Cells"], 0.6, "Step9_NK")
set_label(ids, c(4, 7, 9, 10), "Artifact")
set_label(ids, c(2, 15),       "CD8+ T-Cells")
set_label(ids, 13,             "Stromal")

# ==============================================================================
# 3. Re-cluster Artifact Cells
# ==============================================================================

# Artifact cells
ids <- subcluster_ids(obj, names(L1)[L1 %in% "Artifact"], 0.4, "Artifact_rescue")
set_label(ids, c(0, 1, 2, 4, 5, 6, 8, 16), "Lymphoma")
set_label(ids, c(3, 9),                    "CD4+ T-Cells")
set_label(ids, c(10, 11, 12, 14),          "Endothelial")

ids <- subcluster_ids(obj, names(L1)[L1 %in% "Dendritic Cells"], 0.4, "DC_subcluster")
set_label(ids, 9, "REMOVE")

L1[L1 %in% c("Artifact", "NK Cells")] <- "REMOVE"

# ==============================================================================
# 4.Re-integrate and Re-cluster Lineages
# ==============================================================================

# Myeloid cells
ids <- lineage_cluster_ids(obj, names(L1)[L1 %in% c("Dendritic Cells", "Monocytes",
                                                    "Monocyte/Macrophage", "Granulocytes",
                                                    "Macrophage", "CD11B+ Myeloid")],
                           MYELOID_MARKERS, 0.5, "Myeloid_reint")
set_label(ids, c(0, 1, 2, 5),  "Dendritic Cells")
set_label(ids, c(3, 4),        "Macrophage")
set_label(ids, c(9, 10, 11),   "Monocytes")
set_label(ids, 8,              "Granulocytes")
set_label(ids, 12,             "pDC")
set_label(ids, 6,              "Lymphoma")
set_label(ids, c(7, 13, 14),   "REMOVE")

# CD8+ T cells
ids <- lineage_cluster_ids(obj, names(L1)[L1 %in% "CD8+ T-Cells"], CD8T_MARKERS, 0.4, "CD8T_reint")
set_label(ids, c(4, 5, 12), "CD4+ T-Cells")
set_label(ids, c(10, 11),   "REMOVE")

# Myeloid cells
L1[L1 %in% c("Dendritic Cells", "Macrophage", "Monocytes", "pDC")] <- "Myeloid_cell"

# Granulocytes and CD4+ T cells
gran_cells <- names(L1)[L1 %in% "Granulocytes"]
cd4t_cells <- names(L1)[L1 %in% "CD4+ T-Cells"]
ids_gran <- lineage_cluster_ids_nosketch(obj, gran_cells, GRAN_CD4T_MARKERS, 0.4,
                                         c("DLBCL_CS207490", "DLBCL_402641A", "DLBCL_SP2223235"),
                                         "Gran_reint")
ids_cd4t <- lineage_cluster_ids_nosketch(obj, cd4t_cells, GRAN_CD4T_MARKERS, 0.4,
                                         c("DLBCL_SP206535", "DLBCL_SP2223235"), "CD4T_reint")
set_label(ids_gran, c(0, 5, 7), "Lymphoma")
set_label(ids_gran, 2,          "Myeloid_cell")
set_label(ids_gran, 3,          "CD4+ T-Cells")
set_label(ids_gran, 8,          "CD8+ T-Cells")
set_label(ids_gran, 9,          "REMOVE")
set_label(ids_cd4t, c(4, 9),    "CD8+ T-Cells")
set_label(ids_cd4t, 6,          "Lymphoma")
set_label(ids_cd4t, 11,         "Myeloid_cell")

# ==============================================================================
# 5. Save Annotated Object
# ==============================================================================

obj$cell_type_final <- L1[colnames(obj)]
obj <- subset(obj, cells = colnames(obj)[!is.na(obj$cell_type_final) &
                                           obj$cell_type_final != "REMOVE"])
print(table(obj$cell_type_final))

pdf(file.path(QC_DIR, "Level1_Annotation_UMAP_DotPlot.pdf"), width = 14, height = 10)
DimPlot(obj, reduction = "full.umap", group.by = "cell_type_final", cols = L1_COLS,
        label = TRUE, repel = TRUE, raster = TRUE) + coord_fixed() + NoAxes()
DotPlot(obj, features = LINEAGE_MARKERS, group.by = "cell_type_final", scale = TRUE,
        col.min = -1.5, col.max = 1.5) + RotatedAxis() +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

saveRDS(obj, annotated_obj_path)
