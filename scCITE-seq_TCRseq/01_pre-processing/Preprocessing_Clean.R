library(irlba)
library(data.table)
library(RSpectra)
library(magrittr)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(viridis)
library(compiler)
library(readr)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(R.utils)
library(limma)
library(tidyverse)
library(DoubletFinder)
library(plyr)
library(rlist)
library(hdf5r)
library(Azimuth)
library(SeuratData)
library(patchwork)

set.seed(2024)

# ============================================================
# USER-DEFINED PATHS AND SAMPLE NAMES
# Update these variables to match your data structure
# ============================================================
cellranger.path <- "/path/to/cellranger/output/"
analysis.path   <- "/path/to/analysis/output/"
setwd(analysis.path)

sample.names <- c(
  "SAMPLE1-APH",
  "SAMPLE1-TDN",
  "SAMPLE1-D7",
  "SAMPLE1-4W",
  "SAMPLE1-D0",
  "SAMPLE2-TDN"
)

cell_cycle_genes_path <- "/path/to/regev_lab_cell_cycle_genes.txt"
heat_shock_genes_path <- "/path/to/heat_shock_geneList.txt"

# ============================================================
# STEP 1: LOAD RAW DATA AND CREATE SEURAT OBJECTS
# ============================================================

for (i in seq_along(sample.names)) {
  
  TCR.df <- read_csv(paste0(cellranger.path, sample.names[i], "-TCR/outs/filtered_contig_annotations.csv"))
  TCR.df <- as.data.frame(TCR.df[c("barcode", "raw_clonotype_id", "chain", "v_gene", "j_gene", "d_gene", "c_gene", "productive")])
  TCR.df <- TCR.df[!duplicated(TCR.df$barcode), ]
  rownames(TCR.df) <- TCR.df$barcode
  
  h5.path    <- paste0(cellranger.path, sample.names[i], "/outs/filtered_feature_bc_matrix.h5")
  sample.h5  <- Read10X_h5(h5.path)
  
  # Rename ADT features (remove trailing ".1" suffix added by CellRanger)
  rownames(sample.h5[["Antibody Capture"]]) <- gsub("\\.1$", "", rownames(sample.h5[["Antibody Capture"]]))
  
  seurat.obj <- CreateSeuratObject(
    counts      = sample.h5[["Gene Expression"]],
    min.cells   = 3,
    min.features = 200,
    project     = sample.names[i]
  )
  
  seurat.obj[["ADT"]] <- CreateAssayObject(sample.h5[["Antibody Capture"]][, colnames(seurat.obj)])
  seurat.obj          <- NormalizeData(seurat.obj, assay = "ADT", normalization.method = "CLR")
  seurat.obj[["percent.mt"]]   <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
  seurat.obj          <- PercentageFeatureSet(seurat.obj, "^RP[SL]", col.name = "percent.ribo")
  
  # Add TCR metadata (fill NAs for cells without TCR)
  non.tcr.bc  <- colnames(seurat.obj)[!colnames(seurat.obj) %in% TCR.df$barcode]
  tcr.cols    <- c("barcode", "raw_clonotype_id", "chain", "v_gene", "j_gene", "d_gene", "c_gene", "productive")
  empty.df    <- data.frame(matrix(ncol = length(tcr.cols), nrow = length(non.tcr.bc)))
  colnames(empty.df)  <- tcr.cols
  rownames(empty.df)  <- non.tcr.bc
  empty.df$barcode    <- non.tcr.bc
  
  seurat.obj <- AddMetaData(seurat.obj, metadata = rbind(TCR.df, empty.df))
  seurat.obj <- RenameCells(seurat.obj, new.names = sub("-1", paste0("-", sample.names[i]), Cells(seurat.obj)))
  
  saveRDS(seurat.obj, paste0("seurat_objects/", sample.names[i], "_unfiltered_wTCR.RDS"))
}

# ============================================================
# STEP 2: DIMENSIONALITY REDUCTION ON UNFILTERED OBJECTS
# ============================================================

for (i in seq_along(sample.names)) {
  seurat.obj <- readRDS(paste0("seurat_objects/", sample.names[i], "_unfiltered_wTCR.RDS"))
  
  seurat.obj <- NormalizeData(seurat.obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(assay = "RNA", selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(assay = "RNA", verbose = TRUE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = c(0.4, 0.8), verbose = FALSE)
  
  saveRDS(seurat.obj, paste0("seurat_objects/", sample.names[i], "_unfiltered_wTCR_wReductions.RDS"))
}


# ============================================================
# STEP 3: FILTERING ROUND 1 + DOUBLET DETECTION
# ============================================================

### Dummy Parameters, refer to Supplementary Table ####
for (i in seq_along(sample.names)) {
  seurat.obj <- readRDS(paste0("seurat_objects/", sample.names[i], "_unfiltered_wTCR_wReductions.RDS"))
  
  if (grepl("TDN", sample.names[i])) {
    seurat.obj.filtered <- subset(seurat.obj, subset =
                                    nFeature_RNA > 1000 & nFeature_RNA < 6500 &
                                    percent.mt < 10 &
                                    nCount_ADT < 15000 &
                                    nCount_RNA < 30000) 
  } else {
    seurat.obj.filtered <- subset(seurat.obj, subset =
                                    nFeature_RNA > 500 & nFeature_RNA < 6000 &
                                    percent.mt < 7.5 &
                                    nCount_ADT < 20000 &
                                    nCount_RNA < 30000)
  }
  
  # Cell cycle scoring
  cycle.genes  <- fread(cell_cycle_genes_path, header = FALSE)$V1
  s.genes      <- cycle.genes[1:43]
  g2m.genes    <- cycle.genes[44:97]
  seurat.obj.filtered <- CellCycleScoring(seurat.obj.filtered, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  
  # Heat shock module score
  heat.shock.genes    <- fread(heat_shock_genes_path)$`Approved symbol`
  seurat.obj.filtered <- AddModuleScore(seurat.obj.filtered, features = list(heat.shock.genes), name = "HeatShock.Score")
  
  seurat.obj.filtered <- NormalizeData(seurat.obj.filtered, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(assay = "RNA", selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(assay = "RNA", vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score", "HeatShock.Score1"), verbose = TRUE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = c(0.2, 0.4, 0.6, 0.7, 0.8, 1), verbose = FALSE)
  
  seurat.obj.filtered <- FindDoublets(seurat.obj.filtered, print_plots_to = "images/", save_prefix = sample.names[i])
  
  saveRDS(seurat.obj.filtered, paste0("seurat_objects/", sample.names[i], "_filteredv1_wDoublets_wReductions.RDS"))
}

# ============================================================
# STEP 4: QC VISUALIZATION (FILTERED + DOUBLETS LABELED)
# ============================================================

for (i in seq_along(sample.names)) {
  seurat.obj            <- readRDS(paste0("seurat_objects/", sample.names[i], "_filteredv1_wDoublets_wReductions.RDS"))
  seurat.obj.unfiltered <- readRDS(paste0("seurat_objects/", sample.names[i], "_unfiltered_wTCR_wReductions.RDS"))
  
  n.count    <- length(Cells(seurat.obj))
  n.doublets <- table(seurat.obj$Doublet_Singlet)[["Doublet"]]
  n.filtered <- length(Cells(seurat.obj.unfiltered)) - n.count
  gg.title   <- paste0(sample.names[i], "  Filtered+Doublets  n=", n.count, "  Doublets=", n.doublets, "  Removed=", n.filtered)
  
  seurat.obj.unfiltered$filtered_cells <- !Cells(seurat.obj.unfiltered) %in% Cells(seurat.obj)
  
  pdf(paste0("images/", sample.names[i], "_QC_filteredv1_wDoublets.pdf"), width = 10, height = 6)
  print(DimPlot(seurat.obj.unfiltered, reduction = "umap", group.by = "filtered_cells", raster = FALSE) + ggtitle("Unfiltered UMAP – filtered cells highlighted"))
  print(VlnPlot(seurat.obj, features = c("nCount_ADT", "nFeature_ADT", "percent.ribo"), group.by = "orig.ident", ncol = 3, pt.size = 0))
  print(VlnPlot(seurat.obj, features = c("nCount_RNA",  "nFeature_RNA",  "percent.mt"),  group.by = "orig.ident", ncol = 3, pt.size = 0))
  print(DimPlot(seurat.obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle(gg.title))
  print(DimPlot(seurat.obj, reduction = "umap", group.by = "Doublet_Singlet") + ggtitle(gg.title))
  print(FeaturePlot(seurat.obj, features = c("nFeature_RNA", "nFeature_ADT", "nCount_RNA", "nCount_ADT")) & NoAxes())
  print(FeaturePlot(seurat.obj, features = "percent.mt")   + ggtitle(paste0(gg.title, " – percent.mt")))
  print(FeaturePlot(seurat.obj, features = "percent.ribo") + ggtitle(paste0(gg.title, " – percent.ribo")))
  print(DimPlot(seurat.obj,     reduction = "umap", group.by = "raw_clonotype_id") + ggtitle(gg.title) + NoLegend())
  print(FeaturePlot(seurat.obj, features = c("scfv", "PTPRC", "CD3E", "CD8A", "CD8B", "CD4", "CD14", "NCAM1", "FCGR3A"),
                    min.cutoff = "q1", max.cutoff = "q99", order = TRUE, ncol = 3) & NoAxes())
  print(DimPlot(seurat.obj, reduction = "umap", group.by = "Phase") + ggtitle(paste0(gg.title, " – cell cycle phase")))
  dev.off()
}

# ============================================================
# STEP 5: FILTERING ROUND 2 — REMOVE DOUBLETS
# ============================================================

for (i in seq_along(sample.names)) {
  seurat.obj <- readRDS(paste0("seurat_objects/", sample.names[i], "_filteredv1_wDoublets_wReductions.RDS"))
  
  # Keep all cells for TDN samples; remove doublets for all others
  if (grepl("TDN", sample.names[i])) {
    seurat.obj.filtered <- seurat.obj
  } else {
    seurat.obj.filtered <- subset(seurat.obj, subset = Doublet_Singlet == "Singlet")
  }
  
  vars.regress <- c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score", "HeatShock.Score1")
  
  seurat.obj.filtered <- NormalizeData(seurat.obj.filtered, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(assay = "RNA", selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(assay = "RNA", vars.to.regress = vars.regress, verbose = TRUE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = c(0.4, 0.6, 0.8, 1), verbose = FALSE)
  
  saveRDS(seurat.obj.filtered, paste0("seurat_objects/", sample.names[i], "_filtered_NOdoublets_wReductions.RDS"))
}


# ============================================================
# STEP 6: AZIMUTH CELL TYPE ANNOTATION
# ============================================================

for (i in seq_along(sample.names)) {
  seurat.obj <- readRDS(paste0("seurat_objects/", sample.names[i], "_filtered_NOdoublets_wReductions.RDS"))
  
  # Drop prior Azimuth predictions if re-running
  azimuth.cols <- grep("^predicted\\.celltype", colnames(seurat.obj@meta.data), value = TRUE)
  if (length(azimuth.cols) > 0) seurat.obj@meta.data[, azimuth.cols] <- NULL
  
  DefaultAssay(seurat.obj) <- "RNA"
  seurat.obj <- RunAzimuth(seurat.obj, reference = "pbmcref", query.modality = "RNA")
  
  saveRDS(seurat.obj, paste0("seurat_objects/", sample.names[i], "_Azimuth_mapped.RDS"))
}

