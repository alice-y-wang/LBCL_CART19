################################################################################
## Public CAR-T reference dataset (GSE197268): filter to the published cells
## and re-integrate with Harmony
################################################################################

# ---- Libraries --------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(readr)
library(stringr)
library(Matrix)
library(Azimuth)
library(harmony)

# ---- CONFIG -----------------------------------------------------------------
analysis_dir <- "."
setwd(analysis_dir)

ref_dir <- "GSE197268_data"   # public dataset directory

f_original    <- "GSE197268_ori_Oct24.rds"          # unprocessed object
f_ref_meta    <- "Fig1_CART_global_obs.csv"         # published cell metadata
f_filtered    <- "GSE197268_filtered.RDS"           # output: cells kept
f_integrated  <- "GSE197268_Harmony_Integrated.RDS" # output: after Harmony

set.seed(2024)

# 1. Load object and published metadata # -------

original.data <- readRDS(file.path(ref_dir, f_original))

ref.meta <- read_csv(file.path(ref_dir, f_ref_meta))
ref.meta$cell_bc <- ref.meta$...1


# 2. Parse identifiers, subset to the published cells, attach metadata # -----

split_cols <- str_split_fixed(original.data$patient_id, "-", 2)
original.data$patient   <- split_cols[, 1]
original.data$timepoint <- split_cols[, 2]

split_cols <- str_split_fixed(rownames(original.data@meta.data), "_", 2)
original.data$cell_bc <- split_cols[, 2]

# keep only the cells retained in the published analysis
original.filtered <- subset(original.data, subset = cell_bc %in% ref.meta$cell_bc)

filtered.meta <- original.filtered@meta.data
filtered.meta$full_bc <- rownames(filtered.meta)

ref.meta.new <- left_join(filtered.meta, ref.meta, by = "cell_bc")
rownames(ref.meta.new) <- ref.meta.new$full_bc

original.filtered <- AddMetaData(original.filtered, metadata = ref.meta.new)

# sanity check: cells per patient x channel
count_df <- table(ref.meta.new$patient_id, ref.meta.new$channel, useNA = "ifany") %>%
  as.data.frame() %>%
  rename(patient_id = Var1, channel = Var2) %>%
  filter(Freq > 0)
print(count_df)

saveRDS(original.filtered, file.path(ref_dir, f_filtered))


# 3. Harmony integration # ------

original.filtered <- readRDS(file.path(ref_dir, f_filtered))

original.filtered[["RNA"]] <- JoinLayers(original.filtered[["RNA"]])
original.filtered[["RNA"]] <- split(original.filtered[["RNA"]], f = original.filtered$patient_id)

original.filtered <- NormalizeData(original.filtered, normalization.method = "LogNormalize", scale.factor = 10000)
original.filtered <- FindVariableFeatures(original.filtered, selection.method = "vst", nfeatures = 2000)
original.filtered <- ScaleData(original.filtered, vars.to.regress = c("total_counts", "pct_counts_mt"))
original.filtered <- RunPCA(original.filtered, npcs = 50, reduction.name = "unint.pca")

original.filtered <- IntegrateLayers(
  object = original.filtered, method = HarmonyIntegration,
  orig.reduction = "unint.pca", new.reduction = "integrated_harmony",
  verbose = FALSE
)

original.filtered <- FindNeighbors(original.filtered, reduction = "integrated_harmony", dims = 1:30)
original.filtered <- FindClusters(original.filtered, resolution = 2, cluster.name = "integrated_clusters")
original.filtered <- RunUMAP(original.filtered, reduction = "integrated_harmony", dims = 1:30, reduction.name = "umap_harmony")

saveRDS(original.filtered, file.path(ref_dir, f_integrated))