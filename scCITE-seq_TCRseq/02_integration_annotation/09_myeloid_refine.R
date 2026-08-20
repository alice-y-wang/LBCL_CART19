# =============================================================================
# Monocyte / dendritic-cell subclustering and annotation pipeline
#
# Workflow:
#   1. Load celltypist annotations and reference colour palettes.
#   2. Merge the monocyte-like cluster discovered in the TNK object
#      (cluster 5 from the TNK pipeline) into the standalone myeloid object.
#   3. Use a re-integrated Harmony object, visualize, then drop the merged-in
#      T-cell-origin cells and TDN samples -> v2.
#   4. Re-integrate (separate pipeline), iteratively subcluster CD14, CD16,
#      DC, and inflammatory populations.
#   5. Define annotations
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)        
library(Seurat)

set.seed(2024)


# CONFIG ---------------------------------------------------------------

analysis.path <- "/path/to/analysis/output/"
setwd(analysis.path)

dir.create("images",         showWarnings = FALSE, recursive = TRUE)
dir.create("seurat_objects", showWarnings = FALSE, recursive = TRUE)

# Standard plotting arguments reused throughout
umap_reduction <- "harmony.mono.umap"
fp_args <- list(
  reduction  = umap_reduction,
  min.cutoff = "q2",
  max.cutoff = "q98",
  raster     = FALSE
)


# COLOUR PALETTES ---------------------------------------------------------------

# Immune-Atlas Celltypist
majority.vote.ctypes <- c(
  "Tcm/Naive cytotoxic T cells", "Tem/Trm cytotoxic T cells",
  "Tem/Temra cytotoxic T cells", "CD16+ NK cells",
  "Proliferative germinal center B cells",
  "Tem/Effector helper T cells", "Tcm/Naive helper T cells",
  "Regulatory T cells", "MAIT cells", "pDC", "Classical monocytes",
  "Plasmablasts", "Mast cells", "DC1", "DC2", "Non-classical monocytes"
)
majority.vote.cols <- setNames(
  c("#019477", "#7DB954", "#a8fc83", "#ffe06e", "#c312c9",
    "#2874A6", "#3db1f5", "#CF9FFF", "#DD3F4E", "#f79e9e",
    "#fc8d1e", "#692727", "#ab8585", "darkred", "firebrick3", "sienna2"),
  majority.vote.ctypes
)

# Healthy-COVID-19 PBMC reference labels
covid.vote.ctypes <- c(
  "ASDC", "B_exhausted", "B_malignant", "B_naive",
  "B_non-switched_memory", "B_switched_memory",
  "C1_CD16_mono", "CD14_mono", "CD16_mono", "CD83_CD14_mono", "Mono_prolif",
  "DC1", "DC2", "DC3",
  "CD8.Naive", "CD8.TE", "CD8.EM", "CD8.Prolif",
  "CD4.Naive", "CD4.Prolif", "CD4.CM", "CD4.EM", "CD4.IL22", "CD4.Tfh", "Treg",
  "gdT", "MAIT", "HSC_CD38pos", "HSC_erythroid", "ILC1_3", "NKT",
  "NK_16hi", "NK_56hi", "NK_prolif", "pDC", "Plasmablast", "Platelets"
)
covid.cols <- setNames(
  c("#c64aff", "#732696", "#cf55cb", "#750273", "#c57fc7", "#c312c9",
    "sienna3", "#fc8d1e", "peru", "violetred", "sandybrown",
    "firebrick3", "darkred", "indianred3",
    "#019477", "#a8fc83", "#7DB954", "#076e2c",
    "#3db1f5", "#0839bf", "#3f8cbf", "#2874A6", "#114a70", "#a0dffa", "#CF9FFF",
    "#f590d7", "#ff0000", "pink", "hotpink", "gray", "wheat3",
    "#ffe06e", "#e3e332", "#c7a320", "#ab032d", "#750273", "navy"),
  covid.vote.ctypes
)

# AIFI Immune Health Atlas palettes (loaded from CSV)
load_aifi_palette <- function(csv_path, level) {
  df <- read_csv(csv_path, show_col_types = FALSE)
  setNames(df[[paste0(level, "_color")]], df[[level]])
}
aifi.l1 <- load_aifi_palette("celltypist_mapping/AIFI_colors/AIFI_L1_imm_health_atlas_type_order_colors.csv", "AIFI_L1")
aifi.l2 <- load_aifi_palette("celltypist_mapping/AIFI_colors/AIFI_L2_imm_health_atlas_type_order_colors.csv", "AIFI_L2")
aifi.l3 <- load_aifi_palette("celltypist_mapping/AIFI_colors/AIFI_L3_imm_health_atlas_type_order_colors.csv", "AIFI_L3")


# HELPER ---------------------------------------------------------------

prep_metadata <- function(labels, label_name) {
  colnames(labels)[colnames(labels) == "...1"] <- "cell.bc"
  labels <- as.data.frame(labels)
  rownames(labels) <- labels$cell.bc
  colnames(labels)[colnames(labels) == "majority_voting"] <-
    paste0(label_name, "_majority_voting")
  labels$cell.bc         <- NULL
  labels$over_clustering <- NULL
  labels
}


# LOAD CELLTYPIST ANNOTATIONS --------------------------------------------------------------

celltypist_files <- list(
  AIFI_L1       = "celltypist_mapping/AIFI_L1_MONO_celltypist_majvote_predictions.csv",
  AIFI_L2       = "celltypist_mapping/AIFI_L2_MONO_celltypist_majvote_predictions.csv",
  AIFI_L3       = "celltypist_mapping/AIFI_L3_MONO_celltypist_majvote_predictions.csv",
  Imm_Low       = "celltypist_mapping/Celltypist_Immune_Atlas_MONO_celltypist_majvote_predictions.csv",
  Healthy_Covid = "celltypist_mapping/Healthy_COVID19_PBMC_MONO_celltypist_majvote_predictions.csv"
)
celltypist_meta <- lapply(names(celltypist_files), function(nm) {
  prep_metadata(read_csv(celltypist_files[[nm]], show_col_types = FALSE), nm)
})
names(celltypist_meta) <- names(celltypist_files)


#  MERGE MONOCYTE OBJECT WITH TNK-ORIGIN CLUSTER 5 ---------------------------------------------------------------

cluster.5 <- readRDS("seurat_objects/all_samples_TNK_cluster5Mono_Only.RDS")
mono.old  <- readRDS("seurat_objects/all_samples_MONO_RPCA_wLabels.RDS")
mono.old$old_int_rpca_res.0.2 <- mono.old$int.snn_res.0.2

mono.old$is.from.tnk  <- FALSE
cluster.5$is.from.tnk <- TRUE

mono.new.merge <- merge(mono.old, y = cluster.5)
mono.new.merge[["RNA"]] <- JoinLayers(mono.new.merge[["RNA"]])

saveRDS(mono.new.merge, "seurat_objects/all_samples_MONO_merge_Cluster5.RDS")


# FIRST QC PASS ON THE HARMONY-INTEGRATED MERGED OBJECT ---------------------------------------------------------------

mono.harmony <- readRDS("seurat_objects/all_samples_Mono_Cluster5_Harmony_theta0.2_3500VF_50dimNN.RDS")

# Clustering + UMAP at multiple resolutions
mono.harmony <- FindClusters(mono.harmony,
                             resolution = c(0.2, 0.3, 0.4, 0.5),
                             graph.name = "harmony.snn")
mono.harmony <- FindNeighbors(mono.harmony,
                              reduction  = "harmony.mono",
                              dims       = 1:50,
                              graph.name = c("harmony.nn", "harmony.snn"),
                              prune.SNN  = 1/25)
mono.harmony <- FindClusters(mono.harmony,
                             resolution = c(0.2, 0.3, 0.4, 0.5),
                             graph.name = "harmony.snn")
mono.harmony <- RunUMAP(mono.harmony,
                        dims           = 1:50,
                        reduction      = "harmony.mono",
                        reduction.name = umap_reduction,
                        min.dist       = 0.1,
                        n.neighbors    = 100L,
                        return.model   = TRUE)

# --- QC PDF: metadata groups, reference labels, marker features --------------
# Helper 
plot_dim <- function(obj, group.by, cols = NULL, title = NULL, ...) {
  p <- DimPlot(obj, reduction = umap_reduction, group.by = group.by,
               cols = cols, ...)
  if (!is.null(title)) p <- p + ggtitle(title)
  print(p)
}

pdf("images/all_samples_MONO_Cluster5_Harmony_Theta0.2_3500VF_50dim.pdf",
    height = 8, width = 10)

# Metadata-coloured UMAPs
for (g in c("sample.name", "patient.id", "timepoint", "response",
            "Phase", "Doublet_Singlet")) {
  plot_dim(mono.harmony, group.by = g)
}

# Reference labels (Seurat predicted + AIFI L2/L3 + Healthy COVID)
plot_dim(mono.harmony, "predicted.celltype.l2",        cols = l2.cols)
plot_dim(mono.harmony, "AIFI_L2_majority_voting",      cols = aifi.l2,
         title = "AIFI L2")
plot_dim(mono.harmony, "AIFI_L3_majority_voting",      cols = aifi.l3,
         title = "AIFI L3")
plot_dim(mono.harmony, "Healthy_Covid_majority_voting", cols = covid.cols,
         title = "Healthy Covid Maj Vote")
plot_dim(mono.harmony, "is.from.tnk", title = "is from TNK cluster 5")

# Clustering at each tested resolution
for (res in c(0.2, 0.3, 0.4)) {
  plot_dim(mono.harmony, paste0("harmony.snn_res.", res),
           title = paste("new clustering res", res))
}

# Annotation v3 split by timepoint and response
print(DimPlot(mono.harmony, reduction = umap_reduction,
              group.by = "anno.v3.new", split.by = "timepoint_order",
              ncol = 2, cols = mono.v2.cols, raster = FALSE) +
        ggtitle("new annotations v3"))
print(DimPlot(mono.harmony, reduction = umap_reduction,
              group.by = "anno.v3.new", split.by = "response",
              ncol = 1, cols = mono.v2.cols, raster = FALSE) +
        ggtitle("new annotations v3"))

# Marker feature plots
print(do.call(FeaturePlot,
              c(list(mono.harmony, features = "scfv"),
                fp_args[c("reduction", "min.cutoff", "max.cutoff")])))
print(do.call(FeaturePlot,
              c(list(mono.harmony,
                     features = c("nCount_RNA", "nFeature_RNA",
                                  "percent.mt", "percent.ribo")),
                fp_args)) & NoAxes())
print(do.call(FeaturePlot,
              c(list(mono.harmony,
                     features = c("CD14", "FCGR3A", "HLA-DPB1", "ISG15",
                                  "CD1C", "C1QA", "ITGAX", "ITGAM", "NKG7")),
                fp_args[c("reduction", "min.cutoff", "max.cutoff")])) & NoAxes())
dev.off()

# Quick cluster-resolution check
pdf("images/all_samples_MONO_cluster5_ClusterCheck_Harmony_Theta0.2_3500VF.pdf")
plot_dim(mono.harmony, "is.from.tnk", title = "is from TNK cluster 5")
for (res in c(0.2, 0.3, 0.4, 0.5)) {
  plot_dim(mono.harmony, paste0("harmony.snn_res.", res),
           title = paste("new clustering res", res))
}
dev.off()


# v2: drop TNK-origin cells AND the TDN timepoint ---------------------------------------------------------------

merged.mono  <- readRDS("seurat_objects/all_samples_MONO_merge_Cluster5.RDS")
cells.remove <- Cells(cluster.5)

merged.mono <- subset(merged.mono, cells   = cells.remove, invert = TRUE)
merged.mono <- subset(merged.mono, subset  = timepoint == "TDN", invert = TRUE)

saveRDS(merged.mono, "seurat_objects/all_samples_MONO_merge_T-Removed.RDS")


# RE-INTEGRATED OBJECT (T-removed); cluster, project, annotate ---------------------------------------------------------------

mono.harmony <- readRDS("seurat_objects/all_samples_Mono_T-Removed_Harmony_theta0.2_3500VF_50dimNN.RDS")

mono.harmony <- FindNeighbors(mono.harmony,
                              reduction  = "harmony.mono",
                              dims       = 1:50,
                              graph.name = c("harmony.nn", "harmony.snn"),
                              prune.SNN  = 1/25)
mono.harmony <- FindClusters(mono.harmony,
                             resolution = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1),
                             graph.name = "harmony.snn")
mono.harmony <- RunUMAP(mono.harmony,
                        dims           = 1:50,
                        reduction      = "harmony.mono",
                        reduction.name = umap_reduction,
                        min.dist       = 0.1,
                        n.neighbors    = 100L,
                        return.model   = TRUE)

for (m in celltypist_meta) mono.harmony <- AddMetaData(mono.harmony, metadata = m)

# --- QC PDF (mirrors the pre-removal pass) -----------------------------------
pdf("images/all_samples_MONO_T-removed_Harmony_Theta0.2_3500VF_50dim.pdf",
    height = 8, width = 10)

for (g in c("sample.name", "patient.id", "timepoint", "response",
            "Phase", "Doublet_Singlet")) {
  plot_dim(mono.harmony, g)
}

plot_dim(mono.harmony, "predicted.celltype.l2",        cols = l2.cols)
plot_dim(mono.harmony, "AIFI_L2_majority_voting",      cols = aifi.l2,
         title = "AIFI L2")
plot_dim(mono.harmony, "AIFI_L3_majority_voting",      cols = aifi.l3,
         title = "AIFI L3")
plot_dim(mono.harmony, "Healthy_Covid_majority_voting", cols = covid.cols,
         title = "Healthy Covid Maj Vote")
plot_dim(mono.harmony, "is.from.tnk", title = "is from TNK cluster 5")
for (res in c(0.4, 0.2)) {
  plot_dim(mono.harmony, paste0("harmony.snn_res.", res),
           title = paste("new clustering res", res))
}

print(DimPlot(mono.harmony, reduction = umap_reduction,
              group.by = "anno.v3.new", split.by = "timepoint_order",
              ncol = 2, cols = mono.v2.cols, raster = FALSE) +
        ggtitle("new annotations v3"))
print(DimPlot(mono.harmony, reduction = umap_reduction,
              group.by = "anno.v3.new", split.by = "response",
              ncol = 1, cols = mono.v2.cols, raster = FALSE) +
        ggtitle("new annotations v3"))

print(do.call(FeaturePlot,
              c(list(mono.harmony, features = "scfv"),
                fp_args[c("reduction", "min.cutoff", "max.cutoff")])))
print(do.call(FeaturePlot,
              c(list(mono.harmony,
                     features = c("nCount_RNA", "nFeature_RNA",
                                  "percent.mt", "percent.ribo")),
                fp_args)) & NoAxes())
print(do.call(FeaturePlot,
              c(list(mono.harmony,
                     features = c("CD14", "FCGR3A", "HLA-DPB1", "ISG15",
                                  "CD1C", "C1QA", "ITGAX", "ITGAM", "NKG7")),
                fp_args[c("reduction", "min.cutoff", "max.cutoff")])) & NoAxes())
dev.off()

# Full clustering-resolution sweep
pdf("images/all_samples_MONO_T-removed_CheckClustering_Harmony_Theta0.2_3500VF_50dim.pdf",
    height = 8, width = 10)
for (res in c(0.2, 0.3, 0.4, 0.6, 0.8, 1)) {
  plot_dim(mono.harmony, paste0("harmony.snn_res.", res))
}
dev.off()


# Annotation ---------------------------------------------------------------
mono.harmony$final_anno_level1 <- "CD14 Mono"
mono.harmony$final_anno_level1[mono.harmony$harmony.snn_res.0.4 %in% c(2, 7)] <- "CD16 Mono"
mono.harmony$final_anno_level1[mono.harmony$harmony.snn_res.0.4 == 5]         <- "Dendritic Cell"

