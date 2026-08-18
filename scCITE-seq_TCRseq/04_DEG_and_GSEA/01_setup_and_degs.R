################################################################################
## Setup + differential expression (DEG) for Lymphoma CAR-T scRNA-seq
##
## Description:
##   Loads integrated Seurat objects for the full dataset and the T/NK and
##   monocyte compartments, harmonizes metadata/annotations across objects,
##   and runs differential expression comparing response groups across cell
##   types and timepoints.
##
## Requirements:
##   - The Seurat objects (.RDS) and metadata (.csv) referenced below.
##   - A helper script of custom utilities (`utils_script` in CONFIG).
##
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
project_dir <- "."                                  # project root
setwd(project_dir)

obj_dir      <- file.path(project_dir, "seurat_objects")
files_dir    <- file.path(project_dir, "files")
results_dir  <- file.path(project_dir, "results", "deg")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Input object filenames
f_full_dataset <- "LBCL_full_dataset_clean_new.RDS"
f_tnk_noncar   <- "nonCAR_TNK_filtered.RDS" # need to save and load 
f_tnk_car      <- "CAROnly_TNK_filtered.RDS"
f_tnk_all      <- "all_samples_TNK_filtered.RDS"
f_mono         <- "all_samples_Mono_filtered.RDS"


# Comparison + date tags appended to DEG output filenames
analysis_tag <- "CRvsPD"
date_tag     <- format(Sys.Date(), "%Y%m%d")

# Shared FindMarkers / FindAllMarkers parameters (0.01 thresholds = effectively
# unfiltered, so all fold-changes are retained for downstream filtering).
LOGFC_THRESHOLD <- 0.01
MIN_PCT         <- 0.01
TEST_USE        <- "wilcox"

set.seed(2024)

colors <- c("#E7A75E", "white", "#4FB7C5")

# Optimal-leaf-ordering callback for pheatmap-style clustering
callback <- function(hc, ...) dendsort(hc)

# ---- Load objects -----------------------------------------------------------
full.dataset <- readRDS(file.path(obj_dir, f_full_dataset))
tnk.noncar   <- readRDS(file.path(obj_dir, f_tnk_noncar))
tnk.car      <- readRDS(file.path(obj_dir, f_tnk_car))
tnk.seurat   <- readRDS(file.path(obj_dir, f_tnk_all))
mono.object  <- readRDS(file.path(obj_dir, f_mono))


# ---- Derive the non-CAR (CAR-negative) subset -------------------------------
all.cells   <- Cells(full.dataset)
car.only    <- Cells(tnk.car)
noncar.only <- setdiff(all.cells, car.only)
noncar.obj  <- subset(full.dataset, cells = noncar.only)

################################################################################
## Helper functions
################################################################################

#' Decide whether two groups need downsampling before DEG testing.

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

#' Pairwise DEGs per cell type, per timepoint, between two grouping levels.

run_pairwise_timepoint_degs <- function(obj, split_col, group1, group2,
                                        out_dir, file_tag,
                                        anno_col = "cell.anno",
                                        celltypes = NULL) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  obj$idents.compare <- paste0(obj[[anno_col]][, 1], "_",
                               obj[[split_col]][, 1], "_",
                               obj$timepoint)
  Idents(obj) <- "idents.compare"
  
  if (is.null(celltypes)) celltypes <- unique(obj[[anno_col]][, 1])
  celltypes  <- celltypes[!is.na(celltypes)]
  timepoints <- unique(obj$timepoint)
  timepoints <- timepoints[!is.na(timepoints)]
  
  for (ct in celltypes) {
    for (tp in timepoints) {
      ID.1 <- paste0(ct, "_", group1, "_", tp)
      ID.2 <- paste0(ct, "_", group2, "_", tp)
      
      downsample.check <- check_downsample(obj, ID.1, ID.2)
      if (identical(downsample.check, "skip")) {
        print("skip this comparison")
        next
      }
      
      subset.obj <- obj
      if (isTRUE(downsample.check)) {
        smallest <- min(sum(obj$idents.compare == ID.1),
                        sum(obj$idents.compare == ID.2))
        subset.obj <- subset(obj, downsample = round(smallest))
        print("subsetting to smallest")
      }
      
      marker.degs <- FindMarkers(
        object          = subset.obj,
        ident.1         = ID.1,
        ident.2         = ID.2,
        test.use        = TEST_USE,
        min.pct         = MIN_PCT,
        logfc.threshold = LOGFC_THRESHOLD,
        verbose         = TRUE,
        assay           = "RNA"
      )
      
      marker.degs$geneID    <- rownames(marker.degs)
      marker.degs$celltype  <- ct
      marker.degs$timepoint <- tp
      
      ct.sub <- gsub("/", "", x = ct)
      write.csv(
        marker.degs,
        file.path(out_dir, paste0(
          ct.sub, "_", tp, "_ALLFCs_", file_tag,
          "_downsample_", downsample.check,
          "_", date_tag, ".csv")),
        row.names = FALSE
      )
    }
  }
}

#' FindAllMarkers with the project's standard parameters.

find_all_markers_std <- function(obj, out_file, ident_col = "cell.anno") {
  Idents(obj) <- ident_col
  markers <- FindAllMarkers(
    obj,
    assay             = "RNA",
    logfc.threshold   = LOGFC_THRESHOLD,
    test.use          = TEST_USE,
    slot              = "data",
    min.pct           = MIN_PCT,
    min.diff.pct      = 0.01,
    verbose           = TRUE,
    only.pos          = FALSE,
    max.cells.per.ident = Inf,
    random.seed       = 1,
    latent.vars       = NULL,
    min.cells.feature = 3,
    min.cells.group   = 3,
    mean.fxn          = NULL
  )
  write.csv(markers, out_file, row.names = FALSE)
  markers
}

################################################################################
## 1. Responder vs. non-responder DEGs (per timepoint)
################################################################################

# Non-CAR cells.
run_pairwise_timepoint_degs(
  obj       = noncar.obj,
  split_col = "response",
  group1    = "CR",
  group2    = "PD",
  out_dir   = file.path(results_dir, "noncar_all"),
  file_tag  = analysis_tag
)

# CAR cells.
run_pairwise_timepoint_degs(
  obj       = tnk.car,
  split_col = "response",
  group1    = "CR",
  group2    = "PD",
  out_dir   = file.path(results_dir, "car_all"),
  file_tag  = analysis_tag
)

################################################################################
## 2. FindAllMarkers across cell types
################################################################################

# T/NK compartment
tnk.markers <- find_all_markers_std(
  tnk.seurat,
  file.path(results_dir, "TNK_AllTime_FindAllMarkers.csv")
)

# Monocyte compartment
mono.markers <- find_all_markers_std(
  mono.object,
  file.path(results_dir, "Mono_AllTime_FindAllMarkers.csv")
)