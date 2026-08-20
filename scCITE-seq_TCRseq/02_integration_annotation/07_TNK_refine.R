# =============================================================================
# T/NK cell subclustering and annotation pipeline
#
# Workflow:
#   1. Load celltypist + SCimilarity annotations and add to the integrated
#      TNK Seurat object.
#   2. Visualize integration, clustering resolutions, and reference labels.
#   3. Identify and remove a monocyte-contaminant cluster -> TNK v2.
#   4. Re-integrate (separate pipeline), remove a dying-T-cell cluster -> v3.
#   5. Iteratively subcluster CD4, CD8, NK, gdT, proliferating, and MAIT
#      populations and define three levels of annotation:
#         final_anno_level1 - broad lineage
#         final_anno_level2 - cell.anno
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)        
library(Seurat)
library(data.table)      
library(RColorBrewer)     
library(pheatmap)

set.seed(2024)


# CONFIG # -----

analysis.path <- "/path/to/analysis/output/"
setwd(analysis.path)

dir.create("images",          showWarnings = FALSE, recursive = TRUE)
dir.create("seurat_objects",  showWarnings = FALSE, recursive = TRUE)

# Standard FeaturePlot parameters reused throughout
fp_args <- list(
  reduction  = "harmony.tnk.umap",
  pt.size    = 0.1,
  min.cutoff = "q01",
  max.cutoff = "q99",
  raster     = FALSE
)


# COLORS # -----

majority.vote.ctypes <- c("Tcm/Naive cytotoxic T cells", "Tem/Trm cytotoxic T cells", "Tem/Temra cytotoxic T cells", 
                          "CD16+ NK cells", 
                          "Proliferative germinal center B cells",
                          "Tem/Effector helper T cells", "Tcm/Naive helper T cells", 
                          "Regulatory T cells", 
                          "MAIT cells", "pDC", "Classical monocytes", "Plasmablasts", "Mast cells",
                          "DC1", "DC2", "Non-classical monocytes")
majority.vote.cols <- c("#019477", "#7DB954",  "#a8fc83",
                        "#ffe06e",
                        "#c312c9",
                        "#2874A6", "#3db1f5", 
                        "#CF9FFF",
                        "#DD3F4E","#f79e9e", "#fc8d1e", "#692727",  "#ab8585",
                        "darkred", "firebrick3", "sienna2"
)
names(majority.vote.cols) <- majority.vote.ctypes

covid.vote.ctypes <- c("ASDC", "B_exhausted", "B_malignant", "B_naive", "B_non-switched_memory", "B_switched_memory",
                       "C1_CD16_mono", "CD14_mono", "CD16_mono", "CD83_CD14_mono","Mono_prolif",
                       "DC1", "DC2", "DC3",
                       "CD8.Naive",  "CD8.TE", "CD8.EM" , "CD8.Prolif",  
                       "CD4.Naive", "CD4.Prolif",  'CD4.CM', "CD4.EM", "CD4.IL22", "CD4.Tfh", "Treg",
                       "gdT", "MAIT", "HSC_CD38pos", "HSC_erythroid", "ILC1_3", "NKT", 
                       "NK_16hi", "NK_56hi", "NK_prolif", "pDC", "Plasmablast", "Platelets")

covid.cols <- c(  "#c64aff", "#732696", "#cf55cb", "#750273", "#c57fc7", "#c312c9",
                  "sienna3","#fc8d1e", "peru", "violetred", "sandybrown",
                  "firebrick3", "darkred", "indianred3",
                  "#019477","#a8fc83", "#7DB954", "#076e2c",
                  "#3db1f5", "#0839bf",  "#3f8cbf", "#2874A6", "#114a70", "#a0dffa",  "#CF9FFF",
                  "#f590d7", "#ff0000", "pink", "hotpink", "gray", "wheat3",
                  "#ffe06e", "#e3e332", "#c7a320", "#ab032d",  "#750273", "navy")
names(covid.cols) <- covid.vote.ctypes


# AIFI colors 
AIFI_L1_col <- read_csv("./AIFI_L1_imm_health_atlas_type_order_colors.csv")
AIFI_L2_col <- read_csv("./AIFI_L2_imm_health_atlas_type_order_colors.csv")
AIFI_L3_col <- read_csv("./AIFI_L3_imm_health_atlas_type_order_colors.csv")

aifi.l1 <- AIFI_L1_col$AIFI_L1_color
names(aifi.l1) <- AIFI_L1_col$AIFI_L1

aifi.l2 <- AIFI_L2_col$AIFI_L2_color
names(aifi.l2) <- AIFI_L2_col$AIFI_L2

aifi.l3 <- AIFI_L3_col$AIFI_L3_color
names(aifi.l3) <- AIFI_L3_col$AIFI_L3


# HELPER: standardise celltypist prediction CSVs as Seurat metadata # -----

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


# LOAD CELLTYPIST ANNOTATIONS AND APPLY TO INTEGRATED OBJECT # -----

celltypist_files <- list(
  AIFI_L1       = "celltypist_mapping/AIFI_L1_TNK_celltypist_majvote_predictions.csv",
  AIFI_L2       = "celltypist_mapping/AIFI_L2_TNK_celltypist_majvote_predictions.csv",
  AIFI_L3       = "celltypist_mapping/AIFI_L3_TNK_celltypist_majvote_predictions.csv",
  Imm_Low       = "celltypist_mapping/Celltypist_Immune_Atlas_TNK_celltypist_majvote_predictions.csv",
  Healthy_Covid = "celltypist_mapping/Healthy_COVID19_PBMC_TNK_celltypist_majvote_predictions.csv"
)

celltypist_meta <- lapply(names(celltypist_files), function(nm) {
  prep_metadata(read_csv(celltypist_files[[nm]], show_col_types = FALSE), nm)
})
names(celltypist_meta) <- names(celltypist_files)

# Initial integration object (output of upstream harmony pipeline)
tnk.harmony <- readRDS("seurat_objects/all_samples_TNK_Only_Harmony_theta0.2_regressCC_noDoublets.RDS")

for (m in celltypist_meta) tnk.harmony <- AddMetaData(tnk.harmony, metadata = m)

# VISUALIZE # -------
# ANNOTATION FREQUENCY BAR PLOTS # -----

make_freq_barplot <- function(obj, anno_col, palette, xlab = "timepoint") {
  df <- as.data.frame(table(obj$timepoint_order, obj[[anno_col]][, 1]))
  ggplot(df, aes(fill = Var2, y = Freq, x = Var1)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = palette) +
    xlab(xlab)
}

pdf("images/all_samples_TNKonly_Timepoint_Anno_Bargraphs.pdf", width = 10, height = 6)
print(make_freq_barplot(tnk.harmony, "AIFI_L2_majority_voting",       aifi.l2))
print(make_freq_barplot(tnk.harmony, "AIFI_L3_majority_voting",       aifi.l3))
print(make_freq_barplot(tnk.harmony, "Healthy_Covid_majority_voting", covid.cols))
dev.off()


# CLUSTERING-RESOLUTION SCAN AND QC PLOTS # -----

tnk.harmony <- FindClusters(tnk.harmony,
                            resolution  = c(0.2, 0.4),
                            graph.name  = "harmony.snn",
                            algorithm   = 1)

# visualize in PDF 
pdf("images/all_samples_TNK_ONLY_Harmony_Theta0.2_integrated.pdf", height = 8, width = 10)

meta_groups <- c("sample.name", "patient.id", "timepoint", "response",
                 "Phase", "Doublet_Singlet")
for (g in meta_groups) {
  print(DimPlot(tnk.harmony, reduction = "harmony.tnk.umap", group.by = g))
}

print(DimPlot(tnk.harmony, reduction = "harmony.tnk.umap",
              group.by = "predicted.celltype.l2", cols = l2.cols, raster = FALSE))
for (s in c("timepoint_order", "patient.id", "response")) {
  ncol_i <- switch(s, "timepoint_order" = 3, "patient.id" = 4, "response" = 2)
  print(DimPlot(tnk.harmony, reduction = "harmony.tnk.umap",
                group.by = "predicted.celltype.l2", split.by = s,
                ncol = ncol_i, cols = l2.cols))
}

# multiple resolutions
for (res in c(0.2, 0.4, 1, 1.5, 2)) {
  col <- paste0("harmony.snn_res.", res)
  if (col %in% colnames(tnk.harmony@meta.data)) {
    print(DimPlot(tnk.harmony, reduction = "harmony.tnk.umap",
                  group.by = col, raster = FALSE, label = (res == 0.2)) +
            ggtitle(paste0("integrated clustering res ", res)))
  }
}

# Marker feature plots + dotplot
print(do.call(FeaturePlot, c(list(tnk.harmony, features = "scfv"),
                             fp_args[c("reduction", "min.cutoff", "max.cutoff", "raster")])))
print(do.call(FeaturePlot,
              c(list(tnk.harmony,
                     features = c("CD14", "FCGR3A", "CD4", "CD8A",
                                  "CD3E", "NCAM1", "HLA-DRA", "CCR7")),
                fp_args[c("reduction", "min.cutoff", "max.cutoff", "raster")])) & NoAxes())
print(DotPlot(tnk.harmony,
              features = c("CD14", "FCGR3A", "CD4", "CD8A",
                           "CD3E", "NCAM1", "HLA-DRA", "CCR7", "scfv")))
dev.off()


# Marker DEGs # ------

Idents(tnk.harmony) <- "harmony.snn_res.0.2"
tnk.harmony.down    <- subset(tnk.harmony, downsample = 3000)
Idents(tnk.harmony.down) <- "harmony.snn_res.0.2"

tnk.harmony.down.degs <- FindAllMarkers(
  tnk.harmony.down,
  assay             = "RNA",
  logfc.threshold   = 0.25,
  test.use          = "wilcox",
  slot              = "data",
  min.pct           = 0.01,
  verbose           = TRUE,
  random.seed       = 1,
  min.cells.feature = 15,
  min.cells.group   = 15
)

tnk.degs.filter <- tnk.harmony.down.degs %>%
  group_by(cluster) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  dplyr::filter(avg_log2FC > 0.5,
                p_val_adj  < 0.001,
                pct.1      > 0.1,
                pct.2      > 0.1)

tnk.top.10 <- tnk.degs.filter %>% slice_head(n = 10)
tnk.top.10$cluster <- factor(tnk.top.10$cluster, levels = as.character(0:5))
tnk.top.10 <- tnk.top.10 %>% arrange(cluster)


# Visualize top DEGS # -----

tnk.harmony.down <- ScaleData(tnk.harmony.down, features = tnk.top.10$gene)

degs.select <- tnk.top.10$gene
scaled.mat  <- tnk.harmony.down[["RNA"]]$scale.data[degs.select, ]

# Clip extreme values to the [5%, 95%] quantile range
up_cut  <- quantile(scaled.mat, 0.95, na.rm = TRUE)
low_cut <- min(0, quantile(scaled.mat, 0.05, na.rm = TRUE))
scaled.mat[is.na(scaled.mat)]      <- 0
scaled.mat[scaled.mat > up_cut]    <- up_cut
scaled.mat[scaled.mat < low_cut]   <- low_cut

# Column annotation = cluster identity, ordered, then downsampled for plotting
ann_column <- data.frame(cluster = tnk.harmony.down$harmony.snn_res.0.2,
                         row.names = colnames(tnk.harmony.down))
ann_column$cluster <- factor(ann_column$cluster, levels = as.character(0:5))
ann_column        <- ann_column[order(ann_column$cluster), , drop = FALSE]

set.seed(2023)
ann_column <- ann_column[sort(sample(seq_len(nrow(ann_column)), 5000)), , drop = FALSE]

gaps_row    <- which(diff(as.numeric(factor(tnk.top.10$cluster)))   != 0)
gaps_column <- which(diff(as.numeric(ann_column$cluster))           != 0)

pdf("images/all_samples_TNK_Harmony_downsample_res.0.2_wilcoxDEG.pdf",
    height = 15, width = 15)
print(pheatmap::pheatmap(
  scaled.mat[, rownames(ann_column)],
  cluster_cols   = FALSE,
  cluster_rows   = FALSE,
  show_colnames  = FALSE,
  fontsize       = 12,
  gaps_row       = gaps_row,
  gaps_col       = gaps_column,
  annotation_col = ann_column,
  breaks         = seq(-2, 2, length.out = 100),
  color          = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)
))
dev.off()

# =============================================================================
# TNK v2: remove monocyte-like cluster # --------
# =============================================================================

tnk.harmony$is.mono <- FALSE
tnk.harmony$is.mono[tnk.harmony$harmony.snn_res.0.2 == 5] <- TRUE

tnk.new   <- subset(tnk.harmony, subset = is.mono == FALSE)
cluster.5 <- subset(tnk.harmony, subset = is.mono == TRUE)

saveRDS(cluster.5, "seurat_objects/all_samples_TNK_cluster5Mono_Only.RDS")
saveRDS(tnk.new,   "seurat_objects/all_samples_TNKonly_harmony_monoremoved.RDS")


# RE-INTEGRATED, MONOCYTE-REMOVED OBJECT using same pipeline as 06_integrate_TNK_pipeline.R

tnk.theta0 <- readRDS("seurat_objects/all_samples_TNK_monoremoved_Harmony_theta0.RDS")
for (m in celltypist_meta) tnk.theta0 <- AddMetaData(tnk.theta0, metadata = m)

# Add SCimilarity annotations
load_scim <- function(path) {
  d <- read_csv(path, show_col_types = FALSE)
  d <- as.data.table(d)
  rownames(d) <- d$...1
  d$...1 <- NULL
  d
}
tnk.theta0 <- AddMetaData(tnk.theta0,
                          metadata = load_scim("scimilarity_data/all_samples_TNK_SCimilarity_predictions.csv"))
tnk.theta0 <- AddMetaData(tnk.theta0,
                          metadata = load_scim("scimilarity_data/all_samples_TNK_SCimilarity_targeted_celltype.csv"))

# =============================================================================
# TNK v3: remove dying-T-cell cluster # --------
# =============================================================================
tnk.theta0$anno.l1 <- "T cell"
tnk.theta0$anno.l1[tnk.theta0$new_clustering %in% c("5")] <- "Dying T cell"  # ~11015 cells

tnk.theta0.subset <- subset(tnk.theta0, subset = anno.l1 == "Dying T cell",
                            invert = TRUE)
saveRDS(tnk.theta0.subset, "seurat_objects/all_samples_TNK_mono_dyingremoved.RDS")

# REINTEGRATED MONOCYTE and DYING T CELL REMOVED OBJECT using same pipeline as 06_integrate_TNK_pipeline.R ------

tnk.theta <- readRDS("seurat_objects/all_samples_TNK_mono_dying_removed_Harmony.RDS")

# Iterative subclustering -------
# Helper
subcluster <- function(obj, ident_col, cluster, name, resolution,
                       algorithm = 1, graph = "harmony.snn") {
  Idents(obj) <- ident_col
  FindSubCluster(obj,
                 cluster         = cluster,
                 graph.name      = graph,
                 subcluster.name = name,
                 resolution      = resolution,
                 algorithm       = algorithm)
}

# First-pass subclustering at res 0.35
tnk.theta <- subcluster(tnk.theta, "harmony.snn_res.0.35", 7, "gdt.mait.subclustering", 0.4)
tnk.theta <- subcluster(tnk.theta, "harmony.snn_res.0.35", 0, "prolif.subclustering",   0.3)
tnk.theta <- subcluster(tnk.theta, "harmony.snn_res.0.35", 1, "CD4s.subclustering.2",   0.4)
tnk.theta <- subcluster(tnk.theta, "harmony.snn_res.0.35", 4, "NK.subclustering",       0.4)
tnk.theta <- subcluster(tnk.theta, "harmony.snn_res.0.35", 2, "CD8s.subclustering",     0.2)
tnk.theta <- subcluster(tnk.theta, "CD4s.subclustering",  "1_1", "CD4Naive.subclustering",   0.2)
tnk.theta <- subcluster(tnk.theta, "CD4Naive.subclustering", "1_0", "CD4s.1_0.subclustering", 0.2)


# Broad Lineage from sub-clusters # ------
# Confirmation DEGs at res 0.35
Idents(tnk.theta) <- "harmony.snn_res.0.35"
tnk.theta.down    <- subset(tnk.theta, downsample = 2000)

tnk.theta.DEGs <- FindAllMarkers(
  tnk.theta.down,
  assay             = "RNA",
  logfc.threshold   = 0.25,
  test.use          = "wilcox",
  slot              = "data",
  min.pct           = 0.01,
  verbose           = TRUE,
  random.seed       = 1,
  min.cells.feature = 30,
  min.cells.group   = 30
)

tnk.theta$final_anno_level1 <- "T Cell"
tnk.theta$final_anno_level1[tnk.theta$gdt.mait.subclustering %in% c("7_0", "7_2")] <- "gdT"
tnk.theta$final_anno_level1[tnk.theta$gdt.mait.subclustering == "7_1"]             <- "MAIT"
tnk.theta$final_anno_level1[tnk.theta$harmony.snn_res.0.35 == 0]                   <- "Proliferating T Cell"
tnk.theta$final_anno_level1[tnk.theta$harmony.snn_res.0.35 == 5]                   <- "Proliferating T Cell"
tnk.theta$final_anno_level1[tnk.theta$harmony.snn_res.0.35 == 6]                   <- "Treg"
tnk.theta$final_anno_level1[tnk.theta$harmony.snn_res.0.35 == 4]                   <- "NK"
tnk.theta$final_anno_level1[tnk.theta$CD4s.subclustering == "1_2"]                 <- "CD8 T"
tnk.theta$final_anno_level1[tnk.theta$CD4s.subclustering %in% c("1_0", "1_1")]     <- "CD4 T"
tnk.theta$final_anno_level1[tnk.theta$harmony.snn_res.0.35 %in% c(8)]              <- "CD4 T"
tnk.theta$final_anno_level1[tnk.theta$harmony.snn_res.0.35 %in% c(9)]              <- "Proliferating T Cell"
tnk.theta$final_anno_level1[tnk.theta$harmony.snn_res.0.35 %in% c(2, 3)]           <- "CD8 T"
tnk.theta$final_anno_level1[tnk.theta$CD8s.subclustering %in% c("2_2")]            <- "Proliferating T Cell"
tnk.theta$final_anno_level1[tnk.theta$CD4Naive.subclustering %in% c("1_1_1")]      <- "CD4 T"



# SECOND-PASS SUBCLUSTERING ON BROAD IDENTITIES # ------


tnk.theta$final_clustering_level1 <- "T Cell"
tnk.theta$final_clustering_level1[tnk.theta$gdt.mait.subclustering %in% c("7_0", "7_2")] <- 6
tnk.theta$final_clustering_level1[tnk.theta$gdt.mait.subclustering == "7_1"]             <- 7
tnk.theta$final_clustering_level1[tnk.theta$harmony.snn_res.0.35 == 0]                   <- 1
tnk.theta$final_clustering_level1[tnk.theta$harmony.snn_res.0.35 == 5]                   <- 1
tnk.theta$final_clustering_level1[tnk.theta$harmony.snn_res.0.35 == 6]                   <- 5
tnk.theta$final_clustering_level1[tnk.theta$harmony.snn_res.0.35 == 4]                   <- 4
tnk.theta$final_clustering_level1[tnk.theta$CD4s.subclustering == "1_2"]                 <- 3
tnk.theta$final_clustering_level1[tnk.theta$CD4s.subclustering %in% c("1_0", "1_1")]     <- 2
tnk.theta$final_clustering_level1[tnk.theta$harmony.snn_res.0.35 %in% c(8)]              <- 2
tnk.theta$final_clustering_level1[tnk.theta$harmony.snn_res.0.35 %in% c(9)]              <- 1
tnk.theta$final_clustering_level1[tnk.theta$harmony.snn_res.0.35 %in% c(2, 3)]           <- 3
tnk.theta$final_clustering_level1[tnk.theta$CD8s.subclustering %in% c("2_2")]            <- 1
tnk.theta$final_clustering_level1[tnk.theta$CD4Naive.subclustering %in% c("1_1_1")]      <- 2



tnk.theta <- subcluster(tnk.theta, "final_clustering_level1", 2, "CD4.subclustering",    0.4, algorithm = 2)
tnk.theta <- subcluster(tnk.theta, "final_clustering_level1", 3, "CD8.subclustering",    0.3, algorithm = 2)
tnk.theta <- subcluster(tnk.theta, "final_clustering_level1", 1, "porlif.subclustering", 0.5, algorithm = 2)
tnk.theta <- subcluster(tnk.theta, "final_clustering_level1", 4, "NK.subclustering",     0.2, algorithm = 2)
tnk.theta <- subcluster(tnk.theta, "final_clustering_level1", 6, "gdt.subclustering",    0.2, algorithm = 2)


# Cell Annotation - Refined # -------

tnk.theta$final_anno_level2 <- "T cell"
tnk.theta$final_anno_level2[tnk.theta$porlif.subclustering %in%
                              c("1_0","1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","1_9")] <- "Proliferating"
tnk.theta$final_anno_level2[tnk.theta$final_anno_level1 == "NK"]              <- "NK"
tnk.theta$final_anno_level2[tnk.theta$final_anno_level1 == "Treg"]            <- "Treg"
tnk.theta$final_anno_level2[tnk.theta$final_anno_level1 == "gdT"]             <- "gdT"
tnk.theta$final_anno_level2[tnk.theta$final_anno_level1== "MAIT"]            <- "MAIT"
tnk.theta$final_anno_level2[tnk.theta$CD4.subclustering %in%
                              c("2_5","2_2")]                                                       <- "CD4 Naive-like"
tnk.theta$final_anno_level2[tnk.theta$CD4.subclustering %in%
                              c("2_0","2_1","2_3","2_4","2_6")]                                     <- "CD4 EM-like"
tnk.theta$final_anno_level2[tnk.theta$CD8.subclustering %in%
                              c("3_0","3_1","3_2")]                                                  <- "CD8 EM"
tnk.theta$final_anno_level2[tnk.theta$CD8.subclustering == "3_3"]       <- "CD8 Naive-like"
tnk.theta$final_anno_level2[tnk.theta$c9.subclustering == "9_2"]        <- "CD8 EM"

