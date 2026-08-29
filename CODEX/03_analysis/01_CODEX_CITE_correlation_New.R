################################################################################
## CITE-seq (ADT) vs CODEX marker correlation
################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

options(future.globals.maxSize = 16000 * 1024^2)

# ---- CONFIG -----------------------------------------------------------------
analysis_dir <- "."
setwd(analysis_dir)

objects_dir <- "seurat_objects"
files_dir   <- file.path("results", "cite_codex_correlation", "tables")
plot_dir    <- file.path("results", "cite_codex_correlation", "plots")

for (d in c(files_dir, plot_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

f_full_dataset <- "LBCL_full_dataset_clean_new.RDS"
f_codex        <- "LBCL_CODEX_processed_clean.RDS"

# metadata column names
patient_col        <- "patient.id"        # shared by both objects
cite_celltype_col  <- "cell.anno"         # single-cell / CITE object
codex_celltype_col <- "cell_type_level2"  # CODEX object

# ---- Load objects -----------------------------------------------------------
full.dataset <- readRDS(file.path(objects_dir, f_full_dataset))

# apheresis samples only, ADT assay
aph.samples <- subset(full.dataset, subset = timepoint == "APH")
DefaultAssay(aph.samples) <- "ADT"
aph.samples[["RNA"]] <- NULL

codex.dataset <- readRDS(file.path(objects_dir, f_codex))
codex.dataset[["CODEX"]] <- JoinLayers(codex.dataset[["CODEX"]])

### aggregate through average value ####------

# ---- Features shared between the two panels ---------------------------------
cite_features      <- rownames(full.dataset[["ADT"]])
important_features <- rownames(codex.dataset[["CODEX"]])

correlative_features_cite  <- c(intersect(cite_features, important_features), c("CD11b", "CD11c", "PD1"))
correlative_features_codex <- c(intersect(cite_features, important_features), c("CD11B", "CD11C", "PD-1"))

# ---- helper to pseudobulk marker signal by patient and cell type ------------
aggregate_adt_by_group <- function(seu, patient_col = "patient.id", celltype_col = "cell.anno",
                                   assay = "ADT", slot = "data") {
  adt_mat <- GetAssayData(seu, assay = assay, layer = slot)
  meta <- seu@meta.data
  
  if(!is.null(patient_col)){
    grp <- paste(meta[[patient_col]], meta[[celltype_col]], sep = "__")
    grp_levels <- unique(grp)
    print(grp_levels)
  } else {
    grp <- paste(meta[[celltype_col]])
    grp_levels <- unique(grp)
  }
  
  pb <- sapply(grp_levels, function(g) {
    cols <- which(grp == g)
    Matrix::rowMeans(as.matrix(adt_mat[, cols, drop = FALSE]))
  })
  
  pb <- as.data.frame(pb)
  rownames(pb) <- rownames(adt_mat)
  return(pb)
}

# aggregate once
adt_pb   <- aggregate_adt_by_group(aph.samples,   patient_col = patient_col, celltype_col = cite_celltype_col)
codex_pb <- aggregate_adt_by_group(codex.dataset, patient_col = patient_col, celltype_col = codex_celltype_col,
                                   assay = "CODEX")

adt_pb_cell   <- aggregate_adt_by_group(aph.samples,   patient_col = NULL, celltype_col = cite_celltype_col)
codex_pb_cell <- aggregate_adt_by_group(codex.dataset, patient_col = NULL, celltype_col = codex_celltype_col,
                                        assay = "CODEX")

adt_pb_cell   <- adt_pb_cell[correlative_features_cite, ]
codex_pb_cell <- codex_pb_cell[correlative_features_codex, ]

adt_pb   <- adt_pb[correlative_features_cite, ]
codex_pb <- codex_pb[correlative_features_codex, ]

# CITE marker names -> CODEX spelling. toupper() already resolves CD11b/CD11c.
clean_markers <- function(x) {
  x <- toupper(x)
  x <- gsub("^PD1$", "PD-1", x)
  return(x)
}

rownames(adt_pb_cell) <- clean_markers(rownames(adt_pb_cell))
rownames(adt_pb)      <- clean_markers(rownames(adt_pb))

# ---- cell pair lists --------------------------------------------------------
cell_pairs <- list(
  c("CD8 TEM", "CD8 TEM")
)


cd8_tem_markers <- c(
  "CD8", "CD45RO", "CXCR3", "CD57", "CD27", "CD2",
  "CCR7",                          
  "CCR4", "CCR6",                  
  "TIM3", "PD1",                   
  "CD25",                          
  "CXCR5"                          
)


mono_mac_markers <- c(
  "CD14",    
  "CD16",    
  "CD11B",   
  "HLA-DR",  
  "CD11C",   
  "CD80",    
  "PDL1",    
  "PD-1",    
  "TIM3",    
  "CCR6",
  "CCR4",
  "CXCR3",
  "CCR7",
  "CXCR5",
  "CD45RO"  
)

########## Correlation for CD8EM ####################
all_results <- data.frame()
plots_list <- list()

lineage = FALSE
cell_marker = TRUE

# Helper: convert p-value to significance stars
sig_stars <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}

for (cts in cell_pairs){
  
  cite_cell  <- cts[1]
  codex_cell <- cts[2]
  
  CT.aph.data   <- as.data.frame(adt_pb_cell)   %>% dplyr::select(all_of(cite_cell))
  CT.codex.data <- as.data.frame(codex_pb_cell) %>% dplyr::select(all_of(codex_cell))
  
  colnames(CT.codex.data) <- "CODEX"
  colnames(CT.aph.data)   <- "CITE"
  CT.codex.data$celltype  <- paste0(cite_cell, "_", codex_cell)
  rownames(CT.aph.data)   <- clean_markers(rownames(CT.aph.data))
  
  print(rownames(CT.aph.data))
  CT.combined          <- cbind(CT.codex.data, CT.aph.data)
  CT.combined$celltype <- as.factor(CT.combined$celltype)
  
  # ----- Cell marker subset filter -----
  if (cell_marker) {
    marker_keep <- switch(codex_cell,
                          "CD8 TEM" = cd8_tem_markers
    )
    if (!is.null(marker_keep))
      CT.combined <- CT.combined %>% filter(rownames(CT.combined) %in% marker_keep)
  }
  
  # ----- Correlations -----
  cor_pearson  <- cor.test(CT.combined$CITE, CT.combined$CODEX, method = "pearson")
  
  # ----- Significance label strings -----
  pearson_label  <- sprintf("r = %.2f, p = %.3g %s",
                            cor_pearson$estimate,
                            cor_pearson$p.value,
                            sig_stars(cor_pearson$p.value))
  
  # ----- Flatten statistics -----
  row <- data.frame(
    celltype_pair = paste0(cite_cell, "__", codex_cell),
    pearson_cor   = cor_pearson$estimate,
    pearson_pval  = cor_pearson$p.value
  )
  all_results <- bind_rows(all_results, row)
  
  # ----- Build plots -----
  CT.combined$marker <- rownames(CT.combined)
  
  p <- ggplot(CT.combined, aes(x = CODEX, y = CITE)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_text_repel(aes(label = marker), size = 3) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
             label = pearson_label,  size = 3.5, color = "black") +
    ggtitle(paste0(cite_cell, "_CITE_vs_", codex_cell, "_CODEX")) +
    theme_minimal()
  
  plots_list[[paste0(cite_cell, "__", codex_cell)]] <- p
}

write.csv(all_results,
          file.path(files_dir, "celltype_marker_correlation.csv"))

pdf(file.path(plot_dir, "celltype_marker_correlation_plot.pdf"))
for (p in plots_list) print(p)
dev.off()

############# Correlation for Monocyte and Mono-Mac ##################

# Aggregate by cell.anno_2
aph.samples$cell.anno_2 <- aph.samples[[cite_celltype_col]][, 1]
aph.samples$cell.anno_2[aph.samples$cell.anno_2 %in% c("CD14 Mono", "CD16 Mono")] <- "Monocyte"

adt_pb_mono       <- aggregate_adt_by_group(aph.samples,
                                            patient_col  = patient_col,
                                            celltype_col = "cell.anno_2")
adt_pb_cell_mono  <- aggregate_adt_by_group(aph.samples,
                                            patient_col  = NULL,
                                            celltype_col = "cell.anno_2")

# Subset to correlative features as before
adt_pb_mono      <- adt_pb_mono[correlative_features_cite, ]
adt_pb_cell_mono <- adt_pb_cell_mono[correlative_features_cite, ]

rownames(adt_pb_cell_mono) <- clean_markers(rownames(adt_pb_cell_mono))
rownames(adt_pb_mono)      <- clean_markers(rownames(adt_pb_mono))
# so you can confirm the exact string to grep for monocytes:
print(grep("Mono", colnames(adt_pb_mono), value = TRUE))
print(grep("Mono", colnames(adt_pb_cell_mono), value = TRUE))

# Aggregate codex 
codex.dataset$cell.anno_2 <- codex.dataset[[codex_celltype_col]][, 1]
codex.dataset$cell.anno_2[codex.dataset$cell.anno_2 %in% c("CD14+CD16+TIMs",
                                                           "CD14+TIMs",
                                                           "Macrophage")] <- "Mono-Mac"

codex_pb_mono <- aggregate_adt_by_group(codex.dataset,
                                        patient_col  = patient_col,
                                        celltype_col = "cell.anno_2",
                                        assay = "CODEX")

codex_pb_cell_mono <- aggregate_adt_by_group(codex.dataset,
                                             patient_col  = NULL,
                                             celltype_col = "cell.anno_2",
                                             assay = "CODEX")

# Subset to correlative features as before
codex_pb_mono      <- codex_pb_mono[correlative_features_codex, ]
codex_pb_cell_mono <- codex_pb_cell_mono[correlative_features_codex, ]


cell_pairs_mono <- list(
  c("Monocyte", "Mono-Mac")
)

# helper 
get_markers <- function(cite_cell, codex_cell) {
  switch(codex_cell,
         "Mono-Mac" = mono_mac_markers
  )
}

# Run Correlation # ------
all_results <- data.frame()
plots_list  <- list()
lineage      <- FALSE
cell_marker  <- TRUE

for (cts in cell_pairs_mono) {
  
  cite_cell  <- cts[1]
  codex_cell <- cts[2]
  pair_key   <- paste0(cite_cell, "__", codex_cell)
  
  CT.aph.data   <- as.data.frame(adt_pb_cell_mono)   %>% dplyr::select(all_of(cite_cell))
  CT.codex.data <- as.data.frame(codex_pb_cell_mono) %>% dplyr::select(dplyr::contains(codex_cell))
  
  colnames(CT.codex.data) <- "CODEX"
  colnames(CT.aph.data)   <- "CITE"
  CT.codex.data$celltype  <- pair_key
  rownames(CT.aph.data)   <- clean_markers(rownames(CT.aph.data))
  
  CT.combined          <- cbind(CT.codex.data, CT.aph.data)
  CT.combined$celltype <- as.factor(CT.combined$celltype)
  
  if (lineage || cell_marker) {
    marker_keep <- get_markers(cite_cell, codex_cell)
    if (!is.null(marker_keep))
      CT.combined <- CT.combined %>% filter(rownames(CT.combined) %in% marker_keep)
  }
  
  if (nrow(CT.combined) == 0) {
    warning(paste("No markers remain for:", pair_key)); next
  }
  
  cor_pearson  <- cor.test(CT.combined$CITE, CT.combined$CODEX, method = "pearson")
  
  pearson_label  <- sprintf("r = %.2f, p = %.3g %s",
                            cor_pearson$estimate, cor_pearson$p.value,
                            sig_stars(cor_pearson$p.value))
  
  all_results <- bind_rows(all_results, data.frame(
    celltype_pair = pair_key,
    pearson_cor   = cor_pearson$estimate,   pearson_pval  = cor_pearson$p.value
  ))
  
  CT.combined$marker <- rownames(CT.combined)
  
  p <- ggplot(CT.combined, aes(x = CODEX, y = CITE)) +
    geom_point() + geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_text_repel(aes(label = marker), size = 3) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
             label = pearson_label,  size = 3.5, color = "black") +
    ggtitle(paste0(cite_cell, "_CITE_vs_", codex_cell, "_CODEX")) + theme_minimal()
  
  plots_list[[pair_key]] <- p
}

write.csv(all_results,
          file.path(files_dir, "celltype_marker_all_aggregates_correlation.csv"))

pdf(file.path(plot_dir, "celltype_marker_all_aggregates_correlation_plot.pdf"))
for (p in plots_list) print(p)
dev.off()