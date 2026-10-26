################################################################################
## pySCENIC regulon analysis on hdWGCNA metacells
##
## Imports consensus SCENIC loom files for the CAR+ T/NK, CAR- T/NK and monocyte
## metacell objects, attaches regulon AUC as a "TF" assay, harmonizes patient /
## response / timepoint metadata, then produces waterfall plots and AUCell
## z-score heatmaps.
##
## Cell type is always read from `cell.anno`.
## Timepoints are always one of: APH, TDN, D0, Peak, Week4.
## Response is always one of: CR, PD.
################################################################################

# ---- Libraries --------------------------------------------------------------
library(Seurat)
library(hdWGCNA)
library(loomR)
library(tidyverse)
library(pheatmap)
library(dendsort)
library(SeuratExtend)

# ---- CONFIG -----------------------------------------------------------------
analysis_dir <- "."
setwd(analysis_dir)

objects_dir  <- "seurat_objects"
files_dir    <- "files"
pyscenic_dir <- "pyscenic"
images_dir   <- file.path(pyscenic_dir, "images")
waterfall_dir <- file.path(images_dir, "waterfall_grn")

utils_script <- file.path("06_TF_regulon", "utils.R")

for (d in c(images_dir,
            file.path(waterfall_dir, "CAR"),
            file.path(waterfall_dir, "Myeloid"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Metacell objects (written by the hdWGCNA script)
f_meta_tnk_car <- "metacell_seurat_tnk_car.RDS"
f_meta_mono    <- "metacell_seurat_mono.RDS"

f_meta_tnk_noncar <- "metacell_seurat_tnk_noncar.RDS"

# Full-dataset cell metadata
f_general_meta <- "full_dataset_metadata_annos_new_Update.csv" # load from seurat object

# Consensus SCENIC loom files
f_loom_car    <- "CAR_consensus_scenic_final_loom.loom"
f_loom_noncar <- "NonCAR_consensus_scenic_final_loom.loom"
f_loom_mono   <- "Mono_consensus_scenic_final_loom.loom"

heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
waterfall_colors <- c("#E7A75E", "white", "#4FB7C5")

# optimal leaf ordering callback for pheatmap (same as 00_setup.R)
callback <- function(hc, ...) dendsort(hc)

source(utils_script)
set.seed(2024)

# TFs taken from https://humantfs.ccbr.utoronto.ca/ 
human_tfs <- readLines(file.path(files_dir, "human_tfs.txt"))

# ---- Sample / patient / response map ----------------------------------------
# Sample names per de-identified patient label (PT#); everything downstream uses PT#.
patient_samples <- list(
  PT1  = c("PT1-APH", "PT1-D0", "PT1-4W", "PT1-D7", "PT1-TDN"),
  PT2  = c("PT2-APH", "PT2-D0", "PT2-D28", "PT2-D7", "PT2-TDN"),
  PT3  = c("PT3-APH", "PT3-D0", "PT3-D28", "PT3-PD", "PT3-TDN"),
  PT4  = c("PT4-APH", "PT4-D0", "PT4-4W", "PT4-PD", "PT4-TDN"),
  PT5  = c("PT5-APH", "PT5-D0", "PT5-D28", "PT5-D7", "PT5-TDN-neg", "PT5-TDN-pos"),
  PT6  = c("PT6-APH", "PT6-D0", "PT6-D28", "PT6-D7", "PT6-TDN"),
  PT7  = c("PT7-APH", "PT7-D0", "PT7-D21", "PT7-D7", "PT7-TDN2"),
  PT8  = c("PT8-APH", "PT8-D0", "PT8-D28", "PT8-D7", "PT8-TDN"),
  PT9  = c("PT9_APH", "PT9-D0", "PT9_D7", "PT9_TDN_CARN", "PT9_TDN_CARP"),
  PT10 = c("PT10-APH", "PT10-D7", "PT10-TDN"),
  PT11 = c("PT11-APH", "PT11-D0", "PT11-D28", "PT11-PD", "PT11-TDN"),
  PT12 = c("PT12-APH", "PT12-D0", "PT12-4W", "PT12-D7", "PT12-TDN"),
  PT13 = c("PT13-APH", "PT13-D0", "PT13-D28", "PT13-PD", "PT13-TDN"),
  PT14 = c("PT14-APH", "PT14-D0", "PT14-D28", "PT14-D7", "PT14-TDN"),
  PT15 = c("PT15-APH", "PT15-D0", "PT15-D28", "PT15-D7", "PT15-TDN-neg", "PT15-TDN-pos")
)

cr_patients    <- c("PT1", "PT2", "PT3", "PT4", "PT5", "PT6", "PT7", "PT8")
pd_patients <- c("PT9", "PT10", "PT11", "PT12", "PT13", "PT14", "PT15")


# ============================================================================ #
# Load metacell objects
# ============================================================================ #
tnk.car.metacell <- readRDS(file.path(objects_dir, f_meta_tnk_car))
tnk.car.metacell <- GetMetacellObject(tnk.car.metacell)

mono.metacell <- readRDS(file.path(objects_dir, f_meta_mono))
mono.metacell <- GetMetacellObject(mono.metacell)

# ============================================================================ #
# Import CAR consensus loom file into Seurat
# ============================================================================ #
lfile <- loomR::connect(file.path(pyscenic_dir, f_loom_car), skip.validate = T)
if(is.null(tnk.car.metacell)) {} else if(length(lfile$col.attrs$CellID[]) != length(colnames(tnk.car.metacell))) {
  stop("Loom file and Seurat object have different cell numbers")
}else if(!identical(lfile$col.attrs$CellID[], colnames(tnk.car.metacell))){
  warning("Loom file and Seurat object have the same cell numbers but different cell names")
}

RegulonsAUC <- lfile$col.attrs$RegulonsAUC[]
colnames(RegulonsAUC) <- sub("_(+)", "", colnames(RegulonsAUC), fixed = T)
rownames(RegulonsAUC) <- lfile$col.attrs$CellID[]

Regulons <- lfile$row.attrs$Regulons[]
colnames(Regulons) <- sub("_(+)", "", colnames(Regulons), fixed = T)
genes <- lfile$row.attrs$Gene[]
Regulons <- apply(Regulons, 2, function(x) genes[which(x==1)])

lfile$close_all()
results <- list(RegulonsAUC = RegulonsAUC,
                Regulons = Regulons)

tnk.car.metacell@misc[["SCENIC"]] <- results
tnk.car.metacell[["TF"]] <- CreateAssayObject(data = t(RegulonsAUC))

car_tf_auc <- RegulonsAUC
car_tf_gene_list <- Regulons

# ---- CAR metadata -----------------------------------------------------------
# celltype_sample is "<cell.anno>|<sample.name>"
tnk.car.metacell$sample.name <- sapply(strsplit(tnk.car.metacell$celltype_sample, "\\|"), `[`, 2)

# patient ID
tnk.car.metacell$patient.id <- rep("none", length(Cells(tnk.car.metacell)))
for (pt in names(patient_samples)) {
  tnk.car.metacell$patient.id[tnk.car.metacell$sample.name %in% patient_samples[[pt]]] <- pt
}

# response status
tnk.car.metacell$response <- rep("none", length(Cells(tnk.car.metacell)))
tnk.car.metacell$response[tnk.car.metacell$patient.id %in% cr_patients]    <- "CR"
tnk.car.metacell$response[tnk.car.metacell$patient.id %in% pd_patients] <- "PD"

# timepoint: second field once underscores are normalized to dashes
tnk.car.metacell$timepoint_raw <- str_split_i(gsub("_", "-", tnk.car.metacell$sample.name), "-", 2)

tnk.car.metacell$timepoint <- rep("none", length(Cells(tnk.car.metacell)))
tnk.car.metacell$timepoint[tnk.car.metacell$timepoint_raw %in% c("APH")]              <- "APH"
tnk.car.metacell$timepoint[tnk.car.metacell$timepoint_raw %in% c("TDN", "TDN2")]      <- "TDN"
tnk.car.metacell$timepoint[tnk.car.metacell$timepoint_raw %in% c("D0")]               <- "D0"
tnk.car.metacell$timepoint[tnk.car.metacell$timepoint_raw %in% c("D7", "PD")]         <- "Peak"
tnk.car.metacell$timepoint[tnk.car.metacell$timepoint_raw %in% c("D21", "D28", "4W")] <- "Week4"

tnk.car.metacell$celltype_timepoint_response <- paste0(tnk.car.metacell$cell.anno, "|",
                                                       tnk.car.metacell$timepoint, "|",
                                                       tnk.car.metacell$response)

# ============================================================================ #
# Import Mono consensus loom file into Seurat
# ============================================================================ #
lfile <- loomR::connect(file.path(pyscenic_dir, f_loom_mono), skip.validate = T)
if(is.null(mono.metacell)) {} else if(length(lfile$col.attrs$CellID[]) != length(colnames(mono.metacell))) {
  stop("Loom file and Seurat object have different cell numbers")
}else if(!identical(lfile$col.attrs$CellID[], colnames(mono.metacell))){
  warning("Loom file and Seurat object have the same cell numbers but different cell names")
}

RegulonsAUC <- lfile$col.attrs$RegulonsAUC[]
colnames(RegulonsAUC) <- sub("_(+)", "", colnames(RegulonsAUC), fixed = T)
rownames(RegulonsAUC) <- lfile$col.attrs$CellID[]

Regulons <- lfile$row.attrs$Regulons[]
colnames(Regulons) <- sub("_(+)", "", colnames(Regulons), fixed = T)
genes <- lfile$row.attrs$Gene[]
Regulons <- apply(Regulons, 2, function(x) genes[which(x==1)])

lfile$close_all()
results <- list(RegulonsAUC = RegulonsAUC,
                Regulons = Regulons)

mono.metacell@misc[["SCENIC"]] <- results
mono.metacell[["TF"]] <- CreateAssayObject(data = t(RegulonsAUC))

mono_tf_auc <- RegulonsAUC
mono_tf_gene_list <- Regulons

# ---- Mono metadata ----------------------------------------------------------
# celltype_sample is "<cell.anno>|<cluster>|<sample.name>", so the sample name
# is the THIRD field here (the T/NK objects only have two fields).
mono.metacell$sample.name <- sapply(strsplit(mono.metacell$celltype_sample, "\\|"), `[`, 3)

# patient ID -- same PT# labels as the CAR object
mono.metacell$patient.id <- rep("none", length(Cells(mono.metacell)))
for (pt in names(patient_samples)) {
  mono.metacell$patient.id[mono.metacell$sample.name %in% patient_samples[[pt]]] <- pt
}

# response status
mono.metacell$response <- rep("none", length(Cells(mono.metacell)))
mono.metacell$response[mono.metacell$patient.id %in% cr_patients]    <- "CR"
mono.metacell$response[mono.metacell$patient.id %in% pd_patients] <- "PD"

# timepoint
mono.metacell$timepoint_raw <- str_split_i(gsub("_", "-", mono.metacell$sample.name), "-", 2)

mono.metacell$timepoint <- rep("none", length(Cells(mono.metacell)))
mono.metacell$timepoint[mono.metacell$timepoint_raw %in% c("APH")]              <- "APH"
mono.metacell$timepoint[mono.metacell$timepoint_raw %in% c("TDN", "TDN2")]      <- "TDN"
mono.metacell$timepoint[mono.metacell$timepoint_raw %in% c("D0")]               <- "D0"
mono.metacell$timepoint[mono.metacell$timepoint_raw %in% c("D7", "PD")]         <- "Peak"
mono.metacell$timepoint[mono.metacell$timepoint_raw %in% c("D21", "D28", "4W")] <- "Week4"

mono.metacell$celltype_timepoint_response <- paste0(mono.metacell$cell.anno, "|",
                                                    mono.metacell$timepoint, "|",
                                                    mono.metacell$response)

# sanity check: nothing should be left unassigned
print(table(tnk.car.metacell$patient.id == "none"))
print(table(tnk.car.metacell$timepoint == "none"))
print(table(mono.metacell$patient.id == "none"))
print(table(mono.metacell$timepoint == "none"))

# ============================================================================ #
# CAR T cell waterfall plots (CR vs PD)
# ============================================================================ #
tnk.idents.pairs <- list(
  c("CD4 TEM-like|TDN|CR",   "CD4 TEM-like|TDN|PD"),
  c("CD4 TEM-like|Peak|CR",  "CD4 TEM-like|Peak|PD"),
  c("CD8 TEM|TDN|CR",        "CD8 TEM|TDN|PD"),
  c("CD8 TEM|Peak|CR",       "CD8 TEM|Peak|PD"),
  c("CD8 TEM|Week4|CR",      "CD8 TEM|Week4|PD")
)

Idents(tnk.car.metacell) <- "celltype_timepoint_response"

clean_names <- gsub("\\([+-]\\)", "", rownames(tnk.car.metacell[["TF"]]))
tfs_to_subset <- intersect(clean_names, human_tfs)
tfs_to_subset <- paste0(tfs_to_subset, "(+)")

Std.matr <- Seu2Matr(
  seu = tnk.car.metacell,
  features = tfs_to_subset,
  group.by = NULL
)

for (pairs in tnk.idents.pairs){
  ident.1 <- pairs[1]
  ident.2 <- pairs[2]
  print(c(ident.1, ident.2))
  
  celltype <- str_split(ident.1, "\\|")[[1]][1]
  timepoint <- str_split(ident.1, "\\|")[[1]][2]
  
  print(c(celltype, timepoint))
  
  scores <- WaterfallPlot_Calc(
    matr = t(Std.matr$matr),
    f = Std.matr$f,
    ident.1 = ident.1,
    ident.2 = ident.2,
    exp.transform = FALSE,
    length = "logFC",
    color = "p",
    order = TRUE,
    len.threshold = 0,
    col.threshold = 0,
    log.base = "2",
    top.n = NULL
  )
  
  csub <- gsub("/", "", celltype)
  
  scores$p_value <- 10^(-abs(scores$p))
  scores$padj <- p.adjust(scores$p_value, method = "BH")
  scores$log10_padj <- sign(scores$logFC) * -log10(pmax(scores$padj, 1e-300))
  rownames(scores) <- gsub("tf_", "", scores$rank)
  rownames(scores) <- gsub("\\(\\+\\)", "", rownames(scores))
  scores$rank <- gsub("tf_", "", scores$rank)
  scores$rank <- gsub("\\(\\+\\)", "", scores$rank)
  scores$color <- NULL
  scores <- scores %>%
    dplyr::rename(color = log10_padj)
  
  write.csv(scores, file.path(waterfall_dir, "CAR",
                              paste0("CAR_", csub, "_", timepoint, "_CR_vs_PD_wilcox_scores.csv")))
  
  top_pos <- scores %>%
    dplyr::filter(padj < 0.01) %>%
    arrange(desc(logFC)) %>%
    slice_head(n = 20)
  
  top_neg <- scores %>%
    dplyr::filter(padj < 0.01) %>%
    arrange(logFC) %>%
    slice_head(n = 20) %>% arrange(desc(logFC))
  
  top_combined <- bind_rows(top_pos, top_neg)
  
  p <- WaterfallPlot_Plot_Custom(scores = top_combined, color_theme = waterfall_colors)
  pdf(file.path(waterfall_dir, "CAR",
                paste0("CAR_waterfall_", csub, "_", timepoint, "_CR_vs_PD_wilcox_scores.pdf")),
      width=7, height=5)
  print(p)
  dev.off()
}

# ============================================================================ #
# Myeloid waterfall plots (CR vs PD)
# ============================================================================ #
pairs_myeloid <- list(
  c("CD14 Mono|APH|CR",   "CD14 Mono|APH|PD"),
  c("CD14 Mono|D0|CR",    "CD14 Mono|D0|PD"),
  c("CD14 Mono|Peak|CR",  "CD14 Mono|Peak|PD"),
  c("CD14 Mono|Week4|CR", "CD14 Mono|Week4|PD"),
  
  c("CD16 Mono|APH|CR",   "CD16 Mono|APH|PD"),
  c("CD16 Mono|D0|CR",    "CD16 Mono|D0|PD"),
  c("CD16 Mono|Peak|CR",  "CD16 Mono|Peak|PD"),
  c("CD16 Mono|Week4|CR", "CD16 Mono|Week4|PD"),
  
  c("Dendritic Cell|APH|CR",  "Dendritic Cell|APH|PD"),
  c("Dendritic Cell|D0|CR",   "Dendritic Cell|D0|PD"),
  c("Dendritic Cell|Peak|CR", "Dendritic Cell|Peak|PD")
)

Idents(mono.metacell) <- "celltype_timepoint_response"

clean_names <- gsub("\\([+-]\\)", "", rownames(mono.metacell[["TF"]]))
tfs_to_subset <- intersect(clean_names, human_tfs)
tfs_to_subset <- paste0(tfs_to_subset, "(+)")

Std.matr <- Seu2Matr(
  seu = mono.metacell,
  features = tfs_to_subset,
  group.by = NULL
)

for (pairs in pairs_myeloid){
  ident.1 <- pairs[1]
  ident.2 <- pairs[2]
  print(c(ident.1, ident.2))
  
  celltype <- str_split(ident.1, "\\|")[[1]][1]
  timepoint <- str_split(ident.1, "\\|")[[1]][2]
  
  print(c(celltype, timepoint))
  
  scores <- WaterfallPlot_Calc(
    matr = t(Std.matr$matr),
    f = Std.matr$f,
    ident.1 = ident.1,
    ident.2 = ident.2,
    exp.transform = FALSE,
    length = "logFC",
    color = "p",
    order = TRUE,
    len.threshold = 0,
    col.threshold = 0,
    log.base = "2",
    top.n = NULL
  )
  
  csub <- gsub("/", "", celltype)
  
  scores$p_value <- 10^(-abs(scores$p))
  scores$padj <- p.adjust(scores$p_value, method = "BH")
  scores$log10_padj <- sign(scores$logFC) * -log10(pmax(scores$padj, 1e-300))
  rownames(scores) <- gsub("tf_", "", scores$rank)
  rownames(scores) <- gsub("\\(\\+\\)", "", rownames(scores))
  scores$rank <- gsub("tf_", "", scores$rank)
  scores$rank <- gsub("\\(\\+\\)", "", scores$rank)
  scores$color <- NULL
  scores <- scores %>%
    dplyr::rename(color = log10_padj)
  
  write.csv(scores, file.path(waterfall_dir, "Myeloid",
                              paste0("Myeloid_", csub, "_", timepoint, "_CR_vs_PD_wilcox_scores.csv")))
  
  top_pos <- scores %>%
    dplyr::filter(padj < 0.01) %>%
    arrange(desc(logFC)) %>%
    slice_head(n = 20)
  
  top_neg <- scores %>%
    dplyr::filter(padj < 0.01) %>%
    arrange(logFC) %>%
    slice_head(n = 20) %>% arrange(desc(logFC))
  
  top_combined <- bind_rows(top_pos, top_neg)
  
  p <- WaterfallPlot_Plot_Custom(scores = top_combined, color_theme = waterfall_colors)
  pdf(file.path(waterfall_dir, "Myeloid",
                paste0("Myeloid_waterfall_", csub, "_", timepoint, "_CR_vs_PD_wilcox_scores.pdf")),
      width=7, height=5)
  print(p)
  dev.off()
}


# ============================================================================ #
# AUCell heatmaps: myeloid
# ============================================================================ #
clean_names <- gsub("\\([+-]\\)", "", colnames(mono_tf_auc))
mono_tf_auc_curated <- mono_tf_auc[, clean_names %in% human_tfs]

tf_zscore_curated <- CalcStats(mono_tf_auc_curated, f = mono.metacell$celltype_timepoint_response,
                               method = "zscore", order = "p", n = 10, t = TRUE)

col_order <- c("Dendritic Cell|APH|CR",  "Dendritic Cell|APH|PD",
               "CD14 Mono|APH|CR",       "CD14 Mono|APH|PD",
               "CD16 Mono|APH|CR",       "CD16 Mono|APH|PD",
               "Dendritic Cell|D0|CR",   "Dendritic Cell|D0|PD",
               "CD14 Mono|D0|CR",        "CD14 Mono|D0|PD",
               "CD16 Mono|D0|CR",        "CD16 Mono|D0|PD",
               "Dendritic Cell|Peak|CR", "Dendritic Cell|Peak|PD",
               "CD14 Mono|Peak|CR",      "CD14 Mono|Peak|PD",
               "CD16 Mono|Peak|CR",      "CD16 Mono|Peak|PD")
col_order <- intersect(col_order, colnames(tf_zscore_curated))

# gaps between timepoint blocks
timepoint_order <- sapply(strsplit(col_order, "\\|"), `[`, 2)
gaps_col <- which(diff(as.numeric(factor(timepoint_order, levels = unique(timepoint_order)))) != 0)

ph3 <- pheatmap(tf_zscore_curated[, col_order],
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                cutree_rows = 15,
                gaps_col = gaps_col,
                breaks = seq(-2.5, 2.5, length.out = 100),
                color = heatmap_colors,
                main = "top TF regulons by Z score", clustering_callback = callback)

pdf(file.path(images_dir, "Mono_Consensus_tf_zscore_top10_APH_D0_Peak.pdf"), width=10, height=15)
print(ph3)
dev.off()

# order columns by response (CR then PD), then by timepoint
resp_order <- sapply(strsplit(colnames(tf_zscore_curated), "\\|"), `[`, 3)
tp_order   <- sapply(strsplit(colnames(tf_zscore_curated), "\\|"), `[`, 2)
col_order  <- order(resp_order, tp_order)

tf_zscore_curated.order <- tf_zscore_curated[, col_order]
celltype_order <- sapply(strsplit(colnames(tf_zscore_curated.order), "\\|"), `[`, 1)
gaps_col <- which(diff(as.numeric(factor(celltype_order))) != 0)

ph4 <- pheatmap(tf_zscore_curated.order,
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                cutree_rows = 15,
                gaps_col = gaps_col,
                breaks = seq(-3, 3, length.out = 100),
                color = heatmap_colors,
                main = "top TF regulons by Z score", clustering_callback = callback)

pdf(file.path(images_dir, "Mono_Consensus_tf_zscore_top10_orderedResponse.pdf"), width=15, height=15)
print(ph4)
dev.off()


# ============================================================================ #
# Load CAR- T/NK metacell object
# ============================================================================ #
# (the CAR+ object, its regulon AUC and metadata are already set up above)
tnk.noncar.metacell <- readRDS(file.path(objects_dir, f_meta_tnk_noncar))
tnk.noncar.metacell <- GetMetacellObject(tnk.noncar.metacell)

# ============================================================================ #
# Import NonCAR consensus loom file into Seurat
# ============================================================================ #
lfile <- loomR::connect(file.path(pyscenic_dir, f_loom_noncar), skip.validate = T)
if(is.null(tnk.noncar.metacell)) {} else if(length(lfile$col.attrs$CellID[]) != length(colnames(tnk.noncar.metacell))) {
  stop("Loom file and Seurat object have different cell numbers")
}else if(!identical(lfile$col.attrs$CellID[], colnames(tnk.noncar.metacell))){
  warning("Loom file and Seurat object have the same cell numbers but different cell names")
}

RegulonsAUC <- lfile$col.attrs$RegulonsAUC[]
colnames(RegulonsAUC) <- sub("_(+)", "", colnames(RegulonsAUC), fixed = T)
rownames(RegulonsAUC) <- lfile$col.attrs$CellID[]

Regulons <- lfile$row.attrs$Regulons[]
colnames(Regulons) <- sub("_(+)", "", colnames(Regulons), fixed = T)
genes <- lfile$row.attrs$Gene
Regulons <- apply(Regulons, 2, function(x) genes[which(x==1)])

lfile$close_all()
results <- list(RegulonsAUC = RegulonsAUC,
                Regulons = Regulons)

tnk.noncar.metacell@misc[["SCENIC"]] <- results
tnk.noncar.metacell[["TF"]] <- CreateAssayObject(data = t(RegulonsAUC))

noncar_tf_auc <- RegulonsAUC
noncar_tf_gene_list <- Regulons

# ---- CAR- metadata ----------------------------------------------------------
tnk.noncar.metacell@meta.data <-
  tnk.noncar.metacell@meta.data %>%
  mutate(
    sample.name = str_remove(celltype_sample, str_c(cell.anno, collapse = "|"))
  )


# patient ID
tnk.noncar.metacell$patient.id <- rep("none", length(Cells(tnk.noncar.metacell)))
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT4-4W', "PT4-APH", "PT4-D0", "PT4-PD", "PT4-TDN")] <- "PT4"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT9_APH', "PT9_D7", "PT9_TDN_CARN", "PT9_TDN_CARP", "PT9-D0")] <- "PT9"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT3-APH', "PT3-D0", "PT3-D28", "PT3-PD", "PT3-TDN")] <- "PT3"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT15-APH', "PT15-D0", "PT15-D28", "PT15-D7", "PT15-TDN-neg", "PT15-TDN-pos")] <- "PT15"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT10-APH', "PT10-TDN", "PT10-D7")] <- "PT10"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT5-APH', "PT5-D0", "PT5-D28", "PT5-D7", "PT5-TDN-neg", "PT5-TDN-pos")] <- "PT5"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT6-APH', "PT6-D0", "PT6-D28", "PT6-D7", "PT6-TDN")] <- "PT6"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT8-APH', "PT8-D0", "PT8-D28", "PT8-D7", "PT8-TDN")] <- "PT8"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT14-APH', "PT14-D0", "PT14-D28", "PT14-D7", "PT14-TDN")] <- "PT14"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT2-APH', "PT2-D0", "PT2-D28", "PT2-D7", "PT2-TDN")] <- "PT2"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT13-APH', "PT13-D0", "PT13-D28", "PT13-PD", "PT13-TDN")] <- "PT13"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT12-APH', "PT12-D0", "PT12-4W", "PT12-D7", "PT12-TDN")] <- "PT12"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT1-APH', "PT1-D0", "PT1-4W", "PT1-D7", "PT1-TDN")] <- "PT1"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT7-APH', "PT7-D0", "PT7-D21", "PT7-D7", "PT7-TDN2")] <- "PT7"
tnk.noncar.metacell$patient.id[tnk.noncar.metacell$sample.name %in% c('PT11-APH', "PT11-D0", "PT11-D28", "PT11-PD", "PT11-TDN")] <- "PT11"

# response status
tnk.noncar.metacell$response <- rep("none", length(Cells(tnk.noncar.metacell)))
tnk.noncar.metacell$response[tnk.noncar.metacell$patient.id %in% c("PT4", "PT3", "PT5", "PT6", "PT8", "PT2", "PT1", "PT7")] <- "CR"
tnk.noncar.metacell$response[tnk.noncar.metacell$patient.id %in% c("PT15", "PT10", "PT9", "PT14", "PT13", "PT11", "PT12")] <- "PD"

# timepoint: second dash-separated field
tnk.noncar.metacell$timepoint_raw <- str_split_i(tnk.noncar.metacell$sample.name, "-", i =2)

tnk.noncar.metacell$timepoint <- rep("none", length(Cells(tnk.noncar.metacell)))
tnk.noncar.metacell$timepoint[tnk.noncar.metacell$timepoint_raw %in% c("APH")] <- "APH"
tnk.noncar.metacell$timepoint[tnk.noncar.metacell$timepoint_raw %in% c("D0")] <- "D0"
tnk.noncar.metacell$timepoint[tnk.noncar.metacell$timepoint_raw %in% c("D21", "D28", "4W")] <- "Week4"
tnk.noncar.metacell$timepoint[tnk.noncar.metacell$timepoint_raw %in% c("D7", "PD")] <- "Peak"
tnk.noncar.metacell$timepoint[tnk.noncar.metacell$timepoint_raw %in% c("TDN", "TDN2")] <- "TDN"

# samples whose names do not split cleanly on "-"
tnk.noncar.metacell$timepoint[tnk.noncar.metacell$sample.name %in% c("PT9_APH", "PT3-APH")] <- "APH"
tnk.noncar.metacell$timepoint[tnk.noncar.metacell$sample.name %in% c("PT9_D7", "PT3-PD")] <- "Peak"
tnk.noncar.metacell$timepoint[tnk.noncar.metacell$sample.name %in% c("PT9_TDN_CARN", "PT9_TDN_CARP", "PT3-TDN")] <- "TDN"
tnk.noncar.metacell$timepoint[tnk.noncar.metacell$sample.name %in% c("PT3-D0")] <- "D0"
tnk.noncar.metacell$timepoint[tnk.noncar.metacell$sample.name %in% c("PT3-D28")] <- "Week4"

tnk.noncar.metacell$celltype_timepoint_response <- paste0(tnk.noncar.metacell$cell.anno, "|", tnk.noncar.metacell$timepoint, "|",
                                                          tnk.noncar.metacell$response)


# ============================================================================ #
# AUCell heatmaps: CD8 TEM and CD4 TEM-like, CAR- and CAR+ plotted separately
# ============================================================================ #
# Columns grouped by timepoint (CR next to PD within each timepoint)
plot_tf_heatmap <- function(obj, tf_auc, celltype, timepoints, title, out_file,
                            group_col = "celltype_timepoint_response",
                            n_top = 10, cutree = 10, width = 7, height = 12,
                            zlim = 2.5) {
  lab   <- obj@meta.data[[group_col]]
  parts <- strsplit(lab, "|", fixed = TRUE)
  ct    <- sapply(parts, `[`, 1)
  tp    <- sapply(parts, `[`, 2)
  resp  <- sapply(parts, `[`, 3)
  
  keep  <- ct == celltype & tp %in% timepoints & resp %in% c("CR", "PD")
  cells <- rownames(obj@meta.data)[keep]
  grp   <- paste0(tp[keep], "|", resp[keep])
  cat("\n", title, "\n"); print(table(grp))
  
  # curated human TFs, this object's own regulons only
  auc <- tf_auc[cells, , drop = FALSE]
  auc <- auc[, gsub("\\([+-]\\)", "", colnames(auc)) %in% human_tfs]
  
  zs <- CalcStats(auc, f = grp, method = "zscore", order = "p", n = n_top, t = TRUE)
  
  # CR next to PD within each timepoint
  col_order <- as.vector(t(outer(timepoints, c("CR", "PD"), paste, sep = "|")))
  missing_cols <- setdiff(col_order, colnames(zs))
  if (length(missing_cols) > 0) message("Groups not present: ", paste(missing_cols, collapse = ", "))
  col_order <- intersect(col_order, colnames(zs))
  
  tp_order <- sapply(strsplit(col_order, "|", fixed = TRUE), `[`, 1)
  gaps_col <- which(diff(as.numeric(factor(tp_order, levels = unique(tp_order)))) != 0)
  
  mat <- zs[, col_order]
  mat[mat >  zlim] <-  zlim   # clamp so out-of-range values get the end colors
  mat[mat < -zlim] <- -zlim
  
  ph <- pheatmap::pheatmap(mat,
                           cluster_rows = TRUE,
                           cluster_cols = FALSE,
                           cutree_rows = cutree,
                           gaps_col = gaps_col,
                           breaks = seq(-zlim, zlim, length.out = 100),
                           color = heatmap_colors,
                           main = title,
                           clustering_callback = callback)
  pdf(file.path(images_dir, out_file), width = width, height = height)
  print(ph)
  dev.off()
  
  invisible(zs)
}

# Columns grouped by response (all CR timepoints, then all PD timepoints)
plot_tf_heatmap2 <- function(obj, tf_auc, celltype, timepoints, title, out_file,
                             group_col = "celltype_timepoint_response",
                             n_top = 10, cutree = 10, width = 7, height = 12,
                             zlim=2.5) {
  lab   <- obj@meta.data[[group_col]]
  parts <- strsplit(lab, "|", fixed = TRUE)
  ct    <- sapply(parts, `[`, 1)
  tp    <- sapply(parts, `[`, 2)
  resp  <- sapply(parts, `[`, 3)
  
  keep  <- ct == celltype & tp %in% timepoints & resp %in% c("CR", "PD")
  cells <- rownames(obj@meta.data)[keep]
  grp   <- paste0(tp[keep], "|", resp[keep])
  cat("\n", title, "\n"); print(table(grp))
  
  # curated human TFs, this object's own regulons only
  auc <- tf_auc[cells, , drop = FALSE]
  auc <- auc[, gsub("\\([+-]\\)", "", colnames(auc)) %in% human_tfs]
  
  zs <- CalcStats(auc, f = grp, method = "zscore", order = "p", n = n_top, t = TRUE)
  
  # all timepoints for CR, then all timepoints for PD
  col_order <- c(paste0(timepoints, "|CR"), paste0(timepoints, "|PD"))
  missing_cols <- setdiff(col_order, colnames(zs))
  if (length(missing_cols) > 0) message("Groups not present: ", paste(missing_cols, collapse = ", "))
  col_order <- intersect(col_order, colnames(zs))
  
  # single gap between the CR and PD blocks
  resp_order <- sapply(strsplit(col_order, "|", fixed = TRUE), `[`, 2)
  gaps_col <- which(diff(as.numeric(factor(resp_order, levels = unique(resp_order)))) != 0)
  
  mat <- zs[, col_order]
  mat[mat >  zlim] <-  zlim   # clamp so out-of-range values get the end colors
  mat[mat < -zlim] <- -zlim
  
  ph <- pheatmap::pheatmap(mat,
                           cluster_rows = TRUE,
                           cluster_cols = FALSE,
                           cutree_rows = cutree,
                           gaps_col = gaps_col,
                           breaks = seq(-zlim, zlim, length.out = 100),
                           color = heatmap_colors,
                           main = title,
                           clustering_callback = callback)
  
  pdf(file.path(images_dir, out_file), width = width, height = height)
  print(ph)
  dev.off()
  
  invisible(zs)
}


# ---- all timepoints ---------------------------------------------------------
# Output names follow the myeloid heatmaps, in the same folder (images_dir):
#   <NonCAR|CAR>_<celltype>_Consensus_tf_zscore_top10_<timepoints>.pdf  (by timepoint)
#   <NonCAR|CAR>_<celltype>_Consensus_tf_zscore_top10_orderedResponse.pdf (by response)
car_tf_auc    <- tnk.car.metacell@misc[["SCENIC"]][["RegulonsAUC"]]
noncar_tf_auc <- tnk.noncar.metacell@misc[["SCENIC"]][["RegulonsAUC"]]

# sanity check: metacells x regulons, same cells as the object
dim(car_tf_auc)
all(rownames(car_tf_auc) %in% colnames(tnk.car.metacell))

tp_noncar <- c("APH", "D0", "TDN", "Peak", "Week4")
tp_car    <- c("TDN", "Peak", "Week4")

heatmap_jobs <- list(
  list(obj = tnk.noncar.metacell, auc = noncar_tf_auc, status = "NonCAR", ct = "CD8 TEM",
       tps = tp_noncar, cutree = 12, w = 8, h = 15),
  list(obj = tnk.car.metacell,    auc = car_tf_auc,    status = "CAR",    ct = "CD8 TEM",
       tps = tp_car,    cutree = 8,  w = 6, h = 12),
  list(obj = tnk.noncar.metacell, auc = noncar_tf_auc, status = "NonCAR", ct = "CD4 TEM-like",
       tps = tp_noncar, cutree = 12, w = 8, h = 15),
  list(obj = tnk.car.metacell,    auc = car_tf_auc,    status = "CAR",    ct = "CD4 TEM-like",
       tps = tp_car,    cutree = 8,  w = 6, h = 12)
)

tf_zscores <- list()
for (job in heatmap_jobs) {
  csub  <- gsub("[/ -]", "", job$ct)                      # "CD8TEM", "CD4TEMlike"
  label <- ifelse(job$status == "NonCAR", "CAR-", "CAR+")
  key   <- paste0(job$status, "_", csub)
  
  tf_zscores[[key]] <- plot_tf_heatmap(
    obj        = job$obj,
    tf_auc     = job$auc,
    celltype   = job$ct,
    timepoints = job$tps,
    title      = paste0(label, " ", job$ct, ": top TF regulons by Z score"),
    out_file   = paste0(key, "_Consensus_tf_zscore_top10_", paste(job$tps, collapse = "_"), ".pdf"),
    cutree     = job$cutree, width = job$w, height = job$h
  )
}


# ---- CAR- CD8 TEM and CD4 TEM-like: APH and D0 only -------------------------
tp_early <- c("APH", "D0")

early_jobs <- list(
  list(obj = tnk.noncar.metacell, auc = noncar_tf_auc, status = "NonCAR", ct = "CD8 TEM",
       tps = tp_early, cutree = 6, w = 5, h = 9),
  list(obj = tnk.noncar.metacell, auc = noncar_tf_auc, status = "NonCAR", ct = "CD4 TEM-like",
       tps = tp_early, cutree = 6, w = 5, h = 9)
)

# timepoint-grouped: APH|CR, APH|PD | D0|CR, D0|PD
tf_zscores_early <- list()
for (job in early_jobs) {
  csub  <- gsub("[/ -]", "", job$ct)
  label <- ifelse(job$status == "NonCAR", "CAR-", "CAR+")
  key   <- paste0(job$status, "_", csub)
  
  tf_zscores_early[[key]] <- plot_tf_heatmap(
    obj        = job$obj,
    tf_auc     = job$auc,
    celltype   = job$ct,
    timepoints = job$tps,
    title      = paste0(label, " ", job$ct, " (APH, D0): top TF regulons by Z score"),
    out_file   = paste0(key, "_Consensus_tf_zscore_top10_APH_D0.pdf"),
    cutree     = job$cutree, width = job$w, height = job$h
  )
}

