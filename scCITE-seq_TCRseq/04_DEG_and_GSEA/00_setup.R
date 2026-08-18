################################################################################
## Setup for DEG and GSEA analysis
##
## Description:
##   Loads integrated Seurat objects for the full dataset and the T/NK and
##   monocyte compartments, harmonizes metadata/annotations across objects,
##   removes a flagged sample, and runs differential expression (DEG) analyses
##   comparing response groups across cell types and timepoints.
##
## Requirements:
##   - The Seurat objects (.RDS) and metadata (.csv) referenced below.
##   - A helper script of custom utilities (see `utils_script` in CONFIG).
##
## Notes:
##   - All input/output paths are set in the CONFIG block. Edit those to match
##     your local directory layout before running.
################################################################################

# ---- Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)   # includes dplyr, ggplot2, readr, magrittr, etc.
  library(data.table)
  library(Matrix)
  library(irlba)
  library(RSpectra)
  library(Seurat)
  library(Signac)
  library(DoubletFinder)
  library(liana)
  library(limma)
  library(compiler)
  library(R.utils)
  library(hdf5r)
  library(rlist)
  library(ggrepel)
  library(ggpubr)
  library(dendsort)    # optimal leaf ordering for dendrograms
})

# ---- CONFIG -----------------------------------------------------------------
# Edit these paths to match your environment.
project_dir   <- "."                          # project root
obj_dir       <- file.path(project_dir, "seurat_objects")
files_dir     <- file.path(project_dir, "files")


# Input object filenames
f_full_dataset <- "LBCL_full_dataset_clean_new.RDS"

# need to save new objects
f_tnk_noncar   <- "nonCAR_TNK_filtered.RDS"
f_tnk_car      <- "CAROnly_TNK_filtered.RDS"
f_tnk_all      <- "all_samples_TNK_filtered.RDS"
f_mono         <- "all_samples_Mono_filtered.RDS"
f_metadata     <- "full_dataset_metadata_annotations.csv"

# Tag appended to DEG output filenames
analysis_tag <- "CRvsPD"

set.seed(2024)
setwd(project_dir)
source(utils_script)

colors <- c("#E7A75E", "white", "#4FB7C5")

# Optimal-leaf-ordering callback for pheatmap-style clustering
callback <- function(hc, ...) dendsort(hc)

# ---- Load objects -----------------------------------------------------------
full.dataset <- readRDS(file.path(obj_dir, f_full_dataset))
tnk.noncar   <- readRDS(file.path(obj_dir, f_tnk_noncar))
tnk.car      <- readRDS(file.path(obj_dir, f_tnk_car))
tnk.seurat   <- readRDS(file.path(obj_dir, f_tnk_all))
mono.object  <- readRDS(file.path(obj_dir, f_mono))

# ---- Harmonize metadata across objects --------------------------------------
general.meta <- read.csv(file.path(files_dir, f_metadata))
general.meta <- as.data.frame(general.meta)
rownames(general.meta) <- general.meta$X
general.meta$X <- NULL

meta_cols <- c(
  "basic.anno", "annos_level1", "annos_level2", "annos_level3",
  "CRS.grade", "ICANS.grade", "CAR.exp", "cell.bc",
  "final_anno_level1", "final_anno_level2", "final_anno_level3",
  "mono_anno_level1", "mono_anno_level2",
  "tcell_label_new", "tcell_label_l2", "changed_label",
  "cell.anno", "annos_level3_new", "clone_original_labels"
)
metadata.toadd <- general.meta[, meta_cols]

tnk.noncar   <- AddMetaData(tnk.noncar, metadata.toadd)
full.dataset <- AddMetaData(full.dataset, metadata.toadd)
tnk.car      <- AddMetaData(tnk.car, metadata.toadd)
tnk.seurat   <- AddMetaData(tnk.seurat, metadata.toadd)

# Plotting annotation: keep monocyte/DC level-1 labels, otherwise use level-2
add_plot_anno <- function(obj) {
  obj@meta.data <- obj@meta.data %>%
    dplyr::mutate(
      annos_level2_plot = dplyr::if_else(
        annos_level1 %in% c("CD14 Mono", "CD16 Mono", "Dendritic Cell"),
        annos_level1,
        cell.anno
      )
    )
  obj
}

full.dataset <- add_plot_anno(full.dataset)
tnk.seurat   <- add_plot_anno(tnk.seurat)
tnk.car      <- add_plot_anno(tnk.car)
tnk.noncar   <- add_plot_anno(tnk.noncar)

# Monocyte object uses level-1/level-2 annotations only (no plot anno needed)
mono.object <- AddMetaData(
  mono.object,
  metadata = general.meta[, c(
    "basic.anno", "annos_level1", "annos_level2", "annos_level3",
    "CRS.grade", "ICANS.grade", "CAR.exp", "cell.bc",
    "final_anno_level1", "final_anno_level2", "final_anno_level3",
    "mono_anno_level1", "mono_anno_level2"
  )]
)

# ---- Drop cells with no plotting annotation ---------------------------------
cells_keep <- rownames(full.dataset@meta.data)[!is.na(full.dataset@meta.data$annos_level2_plot)]
full.dataset <- subset(full.dataset, cells = cells_keep)

# ---- Derive the non-CAR (CAR-negative) subset -------------------------------
all.cells   <- Cells(full.dataset)
car.only    <- Cells(tnk.car)
noncar.only <- setdiff(all.cells, car.only)
noncar.obj  <- subset(full.dataset, cells = noncar.only)

################################################################################
## Helper functions: DEG analysis
################################################################################

#' Decide whether two groups need downsampling before DEG testing.
#'
#' @param seurat.obj Seurat object with an `idents.compare` metadata column.
#' @param idents1,idents2 The two identity labels to compare.
#' @return "skip" if either group has < 3 cells; TRUE if downsampling is
#'   warranted (group > 5000 cells or group sizes differ by > 2000);
#'   FALSE otherwise.
check_downsample <- function(seurat.obj, idents1, idents2) {
  
  n1 <- sum(seurat.obj$idents.compare == idents1)
  n2 <- sum(seurat.obj$idents.compare == idents2)
  
  print(c(idents1, idents2))
  print(c(n1, n2))
  
  if (n1 < 3 | n2 < 3) {
    print("less than 3 cells, skip")
    return("skip")
  }
  if (n1 > 5000)              return(TRUE)
  if (n2 > 5000)              return(TRUE)
  if (abs(n2 - n1) > 2000)    return(TRUE)
  return(FALSE)
}

