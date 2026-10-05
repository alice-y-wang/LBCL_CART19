library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

BASE_DIR   <- "/path/to/your/"
SCRIPT_DIR <- "/path/to/your/Scripts"
source(file.path(SCRIPT_DIR, "Utilities", "CODEX_Functions.R"))

T_MARKERS <- c("CD3e", "CD4", "CD8", "CD5", "CD2", "TBET", "FOXP3", "CD25",
               "CD27", "CD45RO", "CCR7", "CXCR5", "CCR6", "CCR4", "CCR3",
               "PD-1", "TIM3", "Vista", "CD69", "GRZB", "CXCR3", "TOX",
               "HLA-A", "KI67", "HIF1A")
MYELOID_MARKERS <- c("CD80", "CD58", "PDL1", "HVEM", "HLA-A",
                     "CD11C", "CD11B", "CD123",
                     "CD16", "CD68", "CD14", "HLA-DR", "CD163 Akoya", "CD15",
                     "CCR6", "CCR3", "CCR4", "CXCR5",
                     "HIF1A", "Serpin B9", "KI67")

# Broad lineages
LINEAGES <- list(
  CD4T    = list(cell_type = "CD4+ T-Cells", markers = T_MARKERS,       res = 0.5),
  CD8T    = list(cell_type = "CD8+ T-Cells", markers = T_MARKERS,       res = 0.5),
  Myeloid = list(cell_type = "Myeloid_cell", markers = MYELOID_MARKERS, res = 0.4)
)

# Mixed / spill-over lineage
CONTAMINATION <- list(
  CD4Hyp   = list(lineage = "CD4T",    cluster = 8, res = 0.5),
  CD8Naive = list(lineage = "CD8T",    cluster = 1, res = 0.3),
  CD8Prol  = list(lineage = "CD8T",    cluster = 5, res = 0.2),
  GranHyp  = list(lineage = "Myeloid", cluster = 9, res = 0.5)
)


SEURAT_DIR  <- file.path(BASE_DIR, "SeuratObj")
QC_DIR      <- file.path(BASE_DIR, "QC_image", "Annotation")
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

annotated_l1_path   <- file.path(SEURAT_DIR, "DLBCL_CODEX_Annotated_L1.rds")   # from 2a
final_obj_path      <- file.path(SEURAT_DIR, "DLBCL_CODEX_Final.rds")

# ==============================================================================
# 1. Load Level-1 Annotated Object
# ==============================================================================

obj <- readRDS(annotated_l1_path)
DefaultAssay(obj) <- "CODEX"
obj <- JoinLayers(obj)
Idents(obj) <- "orig.ident"

obj@meta.data <- obj@meta.data[, !colnames(obj@meta.data) %in% rownames(obj[["CODEX"]])]

# ==============================================================================
# 2. Re-integrate and Cluster Each Lineage (CD4+ T, CD8+ T, Myeloid)
# ==============================================================================

lineage_ids <- list()
for (lin in names(LINEAGES)) {
  cfg <- LINEAGES[[lin]]
  sub <- subset(obj, cells = colnames(obj)[obj$cell_type_final == cfg$cell_type])
  sub <- reintegrate_sketch_rpca(sub, markers = cfg$markers)
  ids <- cluster_sketch_project(sub, markers = cfg$markers, res = cfg$res)
  rm(sub); gc()

  qc <- subset(obj, cells = names(ids))
  qc$cluster <- ids[colnames(qc)]
  pdf(file.path(QC_DIR, paste0(lin, "_Subclusters_DotPlot.pdf")), width = 14, height = 8)
  print(DotPlot(qc, features = cfg$markers, group.by = "cluster",
                scale = TRUE, col.min = -1.5, col.max = 1.5) + RotatedAxis() +
          scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
  dev.off()

  lineage_ids[[lin]] <- data.frame(barcode = names(ids), lineage = lin, cluster = unname(ids))
  rm(qc); gc()
}
lineage_ids <- bind_rows(lineage_ids)

# ==============================================================================
# 3. Sub-cluster Mixed Lineage Clusters
# ==============================================================================

lineage_ids$subcluster <- NA_character_

for (pop in names(CONTAMINATION)) {
  cfg <- CONTAMINATION[[pop]]
  pool <- lineage_ids$barcode[lineage_ids$lineage == cfg$lineage &
                                lineage_ids$cluster == as.character(cfg$cluster)]

  sub <- subset(obj, cells = pool)
  sub <- FindNeighbors(sub, reduction = "integrated.rpca.full", dims = 1:30)
  sub <- FindClusters(sub, resolution = cfg$res)
  res_col <- paste0("CODEX_snn_res.", cfg$res)

  pdf(file.path(QC_DIR, paste0(pop, "_Subclusters_DotPlot.pdf")), width = 16, height = 6)
  print(DotPlot(sub, features = union(T_MARKERS, MYELOID_MARKERS), group.by = res_col,
                scale = TRUE, col.min = -1.5, col.max = 1.5) + RotatedAxis() +
          scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
  dev.off()

  idx <- match(colnames(sub), lineage_ids$barcode)
  lineage_ids$subcluster[idx] <- as.character(sub@meta.data[[res_col]])
  rm(sub); gc()
}


# ==============================================================================
# 4. Final Cell Type Annotation
# ==============================================================================

meta <- data.frame(barcode = colnames(obj),
                   L1_in = as.character(obj$cell_type_final)) %>%
  left_join(lineage_ids, by = "barcode")
L2 <- meta$L1_in

set_label <- function(lin, clusters, label, subclusters = NULL) {
  idx <- which(meta$lineage == lin & meta$cluster %in% as.character(clusters))
  if (!is.null(subclusters)) idx <- idx[meta$subcluster[idx] %in% as.character(subclusters)]
  L2[idx] <<- label
}

# CD4+ T cells
set_label("CD4T", c(0, 2, 11),       "CD4T effector-exhausted")
set_label("CD4T", c(1, 3, 5, 6, 9),  "CD4 TEM")
set_label("CD4T", 4,                 "Macrophage")
set_label("CD4T", 7,                 "Lymphoma")
set_label("CD4T", 10,                "cDC")
set_label("CD4T", 8, "REMOVE",                  subclusters = 0)
set_label("CD4T", 8, "CD4T effector-exhausted", subclusters = c(1, 4))
set_label("CD4T", 8, "CD4 TEM",                 subclusters = 2)
set_label("CD4T", 8, "CD4 Naive-like",          subclusters = 3)
set_label("CD4T", 8, "Macrophage",              subclusters = 5)
set_label("CD4T", 8, "T-regulatory",            subclusters = 6)

# CD8+ T cells
set_label("CD8T", c(0, 6),           "CD8 TEM")
set_label("CD8T", c(2, 3, 10),       "CD8T effector-exhausted")
set_label("CD8T", 4,                 "Macrophage")
set_label("CD8T", 7,                 "CD4 Naive-like")
set_label("CD8T", 9,                 "Lymphoma")
set_label("CD8T", 11,                "CD14+TIMs")
set_label("CD8T", c(8, 12),          "REMOVE")
set_label("CD8T", 1, "Lymphoma",                subclusters = c(0, 1, 2, 3, 5))
set_label("CD8T", 1, "CD8T effector-exhausted", subclusters = 4)
set_label("CD8T", 1, "REMOVE",                  subclusters = 6)
set_label("CD8T", 5, "Proliferating T",         subclusters = c(0, 1, 3, 6))
set_label("CD8T", 5, "Macrophage",              subclusters = 2)
set_label("CD8T", 5, "Lymphoma",                subclusters = 4)
set_label("CD8T", 5, "Endothelial",             subclusters = 5)

# Myeloid cells
set_label("Myeloid", c(0, 1, 2, 4),  "Lymphoma")
set_label("Myeloid", c(3, 5, 6, 10), "Macrophage")
set_label("Myeloid", 7,              "CD14+TIMs")
set_label("Myeloid", 8,              "CD14+CD16+TIMs")
set_label("Myeloid", 11,             "cDC")
set_label("Myeloid", 12,             "Granulocytes")
set_label("Myeloid", 9, "Endothelial",  subclusters = 0)
set_label("Myeloid", 9, "Granulocytes", subclusters = c(1, 2, 4, 5, 6, 7, 8))
set_label("Myeloid", 9, "REMOVE",       subclusters = c(3, 9, 10))

L1 <- dplyr::recode(L2,
  "CD4 Naive-like" = "CD4+ T-Cells", "CD4 TEM" = "CD4+ T-Cells",
  "CD4T effector-exhausted" = "CD4+ T-Cells",
  "CD8 TEM" = "CD8+ T-Cells", "CD8T effector-exhausted" = "CD8+ T-Cells",
  "Proliferating T" = "CD8+ T-Cells",
  "Macrophage" = "Myeloid_cell", "CD14+TIMs" = "Myeloid_cell",
  "CD14+CD16+TIMs" = "Myeloid_cell", "cDC" = "Myeloid_cell")
L1[which(meta$lineage == "CD8T" & meta$cluster == "5" & meta$subcluster == "3")] <- "CD4+ T-Cells"

obj@meta.data[[L2_COL]] <- L2
obj@meta.data[[L1_COL]] <- L1
obj <- subset(obj, cells = colnames(obj)[L2 != "REMOVE"])

print(table(obj@meta.data[[L2_COL]]))
print(table(obj@meta.data[[L1_COL]], obj@meta.data[[L2_COL]]))

# ==============================================================================
# 5. Save Final Annotated Object
# ==============================================================================

saveRDS(obj, final_obj_path)

pdf(file.path(QC_DIR, "Final_Annotation_UMAP.pdf"), width = 16, height = 12)
DimPlot(obj, reduction = "full.umap", group.by = L1_COL, cols = L1_COLS,
        label = TRUE, repel = TRUE, raster = TRUE) + coord_fixed() + NoAxes()
DimPlot(obj, reduction = "full.umap", group.by = L2_COL, cols = L2_COLS,
        label = TRUE, repel = TRUE, raster = TRUE) + coord_fixed() + NoAxes()
dev.off()
