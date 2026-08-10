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
analysis.path    <- "/path/to/analysis/output/"
integration.path <- file.path(analysis.path, "seurat_objects/for_integration")
setwd(analysis.path)

# ============================================================
# STEP 1: LIST OBJECTS TO MERGE
# Expects files named: <sample_name>_NEW_Azimuth_v5_Mapping.RDS
# ============================================================
list.objs        <- list.files(integration.path, pattern = "\\.RDS$")
list.sample.names <- sub("_NEW_Azimuth_v5_Mapping.RDS", "", list.objs)

message(paste0("Objects found: ", length(list.objs)))

# ============================================================
# STEP 2: ITERATIVELY MERGE ALL SAMPLES
# ============================================================
merged.obj              <- readRDS(file.path(integration.path, list.objs[1]))
merged.obj$sample.name  <- list.sample.names[1]
DefaultAssay(merged.obj) <- "RNA"

for (i in 2:length(list.objs)) {
  message("Merging: ", list.sample.names[i])
  
  tmp               <- readRDS(file.path(integration.path, list.objs[i]))
  tmp$sample.name   <- list.sample.names[i]
  tmp@project.name  <- list.sample.names[i]
  DefaultAssay(tmp) <- "RNA"
  
  # Drop SCT assay and associated metadata if present
  if ("SCT" %in% names(tmp@assays)) {
    tmp[["SCT"]]       <- NULL
    tmp$nCount_SCT     <- NULL
    tmp$nFeature_SCT   <- NULL
  }
  
  # Drop predicted ADT assay if present
  if ("predicted_ADT" %in% names(tmp@assays)) {
    tmp[["predicted_ADT"]] <- NULL
  }
  
  merged.obj <- merge(merged.obj, tmp)
}

# ============================================================
# STEP 3: SAVE
# ============================================================
saveRDS(merged.obj, file.path(integration.path, "all_samples_merged.RDS"))