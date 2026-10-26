library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)

BASE_DIR   <- "/path/to/your/"
SCRIPT_DIR <- "/path/to/your/Scripts"
source(file.path(SCRIPT_DIR, "Utilities", "CODEX_Functions.R"))

T_MARKERS <- c("CD3e", "CD4", "CD8", "CD5", "CD2", "TBET", "FOXP3", "CD25", "GATA3",
               "CD27", "CD45RO", "CCR7", "CXCR5", "CCR6", "CCR4", "CCR3",
               "PD-1", "LAG3", "TIM3", "Vista", "ICOS", "CD69", "GRZB", "CXCR3",
               "TOX", "HLA-A", "KI67", "HIF1A")
MYELOID_MARKERS <- c("CD80", "CD58", "PDL1", "HVEM", "HLA-A",
                     "CD11C", "CD11B", "CD123",
                     "CD16", "CD68", "CD14", "HLA-DR", "CD163 Akoya", "CD15",
                     "CCR6", "CCR3", "CCR4", "CXCR5",
                     "HIF1A", "Serpin B9", "KI67")

T_UMAP       <- list(n.neighbors = 200, min.dist = 0.05)
MYELOID_UMAP <- list(n.neighbors = 100, min.dist = 0.01)

SEURAT_DIR <- file.path(BASE_DIR, "SeuratObj")
FIG_DIR    <- file.path(BASE_DIR, "Results", "Figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

final_obj_path   <- file.path(SEURAT_DIR, "DLBCL_CODEX_Final.rds")
tcell_obj_path   <- file.path(SEURAT_DIR, "DLBCL_CODEX_Tcells_ReIntegrated.rds")
myeloid_obj_path <- file.path(SEURAT_DIR, "DLBCL_CODEX_Myeloid_ReIntegrated.rds")

# ==============================================================================
# 1. Load Final Annotated Object
# ==============================================================================

obj <- readRDS(final_obj_path)
DefaultAssay(obj) <- "CODEX"
l2 <- as.character(obj@meta.data[[L2_COL]])

# ==============================================================================
# 2. T Cells
# ==============================================================================

t_obj <- subset(obj, cells = colnames(obj)[l2 %in% T_ORDER])
t_obj <- reintegrate_sketch_rpca(t_obj, markers = T_MARKERS)

DefaultAssay(t_obj) <- "CODEX"
t_obj <- RunUMAP(t_obj, reduction = "integrated.rpca.full", dims = 1:15,
                 reduction.name = "umap_T")
saveRDS(t_obj, tcell_obj_path)

t_draw <- c("CD4 TEM", "CD8T effector-exhausted", "CD8 TEM", "CD4T effector-exhausted",
            "T-regulatory", "Proliferating T", "CD4 Naive-like")
t_obj@meta.data[[L2_COL]] <- factor(t_obj@meta.data[[L2_COL]], levels = t_draw)

p <- DimPlot(t_obj, reduction = "umap_T", group.by = L2_COL, cols = T_COLS[t_draw],
             order = rev(t_draw), label = TRUE, repel = TRUE, label.size = 4,
             raster = TRUE) +
  ggtitle(paste0("T cells (N = ", format(ncol(t_obj), big.mark = ","), ")")) +
  coord_fixed() + NoAxes() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11))
ggsave(p, filename = file.path(FIG_DIR, "Tcell_UMAP.pdf"), width = 16, height = 12)

# ==============================================================================
# 3. Myeloid Cells
# ==============================================================================

m_obj <- subset(obj, cells = colnames(obj)[l2 %in% MYELOID_ORDER])
rm(obj); gc()

for (a in setdiff(Assays(m_obj), "CODEX")) m_obj[[a]] <- NULL
for (r in Reductions(m_obj)) m_obj[[r]] <- NULL
m_obj <- JoinLayers(m_obj)
m_obj[["CODEX"]] <- CreateAssay5Object(counts = GetAssayData(m_obj, assay = "CODEX",
                                                             layer = "counts"))
DefaultAssay(m_obj) <- "CODEX"

ob.list <- SplitObject(m_obj, split.by = "orig.ident")
rm(m_obj); gc()
ob.list <- lapply(ob.list, function(x) {
  x <- NormalizeData(x, normalization.method = "CLR", margin = 1, verbose = FALSE)
  x <- ScaleData(x, features = MYELOID_MARKERS, verbose = FALSE)
  RunPCA(x, features = MYELOID_MARKERS, approx = FALSE, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = ob.list, dims = 1:15,
                                  anchor.features = MYELOID_MARKERS,
                                  reduction = "rpca", k.anchor = 20)
m_obj <- IntegrateData(anchorset = anchors, dims = 1:15, k.weight = 30)
rm(anchors, ob.list); gc()

DefaultAssay(m_obj) <- "integrated"
m_obj <- ScaleData(m_obj, verbose = FALSE)
m_obj <- RunPCA(m_obj, features = MYELOID_MARKERS, npcs = 30, verbose = FALSE)
m_obj <- RunUMAP(m_obj, reduction = "pca", dims = 1:15,
                 reduction.name = "umap_Myeloid")
saveRDS(m_obj, myeloid_obj_path)

m_obj@meta.data[[L2_COL]] <- factor(m_obj@meta.data[[L2_COL]], levels = MYELOID_ORDER)

p <- DimPlot(m_obj, reduction = "umap_Myeloid", group.by = L2_COL, cols = MYELOID_COLS,
             order = c("Macrophage", "CD14+TIMs", "CD14+CD16+TIMs", "cDC"),
             label = TRUE, repel = TRUE, label.size = 4,
             raster = TRUE, raster.dpi = c(1024, 1024)) +
  ggtitle(paste0("Myeloid cells (N = ", format(ncol(m_obj), big.mark = ","), ")")) +
  coord_fixed() + NoAxes() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11))
ggsave(p, filename = file.path(FIG_DIR, "Fig5b_Myeloid_UMAP.pdf"), width = 14, height = 12)
