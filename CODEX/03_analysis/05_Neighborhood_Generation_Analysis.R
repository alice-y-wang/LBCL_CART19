library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(imcRtools)
library(SingleCellExperiment)
library(reshape2)
library(pheatmap)
library(ggplotify)

BASE_DIR   <- "/path/to/your/"
SCRIPT_DIR <- "/path/to/your/Scripts"
source(file.path(SCRIPT_DIR, "Utilities", "CODEX_Functions.R"))

final_obj_path   <- file.path(BASE_DIR, "SeuratObj", "DLBCL_CODEX_Final.rds")
neighborhood_dir <- file.path(BASE_DIR, "Results", "Neighborhood_Files")
FIG_DIR          <- file.path(BASE_DIR, "Results", "Figures")
dir.create(neighborhood_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

Col_order <- c("Lymphoma", "Macrophage", "CD14+TIMs", "CD14+CD16+TIMs", "cDC",
               "Granulocytes", "T-regulatory", "CD4 Naive-like", "CD4 TEM",
               "CD4T effector-exhausted", "CD8 TEM", "CD8T effector-exhausted",
               "Proliferating T", "Endothelial", "Stromal")
FDR_THRESHOLD <- 0.001

CN_REMAP      <- c("1" = 6, "2" = 10, "3" = 2, "4" = 3, "5" = 5,
                   "6" = 4, "7" = 7, "8" = 1, "9" = 8, "10" = 9)

# ==============================================================================
# 1. Parse command-line arguments
# ==============================================================================
args <- commandArgs(TRUE)
Kval <- as.numeric(args[1])
Nval <- as.numeric(args[2])
tag  <- sprintf("k%d_n%d", Kval, Nval)

# ==============================================================================
# 2. Load and prepare data
# ==============================================================================
obj <- readRDS(final_obj_path)
Idents(obj) <- L2_COL
obj <- subset(obj, idents = intersect(Col_order, unique(obj@meta.data[[L2_COL]])))
DefaultAssay(obj) <- "CODEX"
obj <- JoinLayers(obj)

keep_ct <- obj@meta.data %>%
  group_by(.data[[L2_COL]]) %>%
  summarise(sample_count = n_distinct(.data[[SAMPLE_COL]]), .groups = "drop") %>%
  filter(sample_count >= 4) %>% pull(1) %>% as.character()
obj <- subset(obj, idents = keep_ct)
Col_order <- intersect(Col_order, keep_ct)

# ==============================================================================
# 3. Build spatial graph and identify neighborhoods
# ==============================================================================
sce <- SingleCellExperiment(assays = list(logcounts = LayerData(obj, layer = "data")),
                            colData = obj@meta.data)
colnames(colData(sce))[colnames(colData(sce)) == "x.coord"] <- "Pos_X"
colnames(colData(sce))[colnames(colData(sce)) == "y.coord"] <- "Pos_Y"

sce <- buildSpatialGraph(sce, img_id = SAMPLE_COL, type = "knn", k = Kval)
sce <- aggregateNeighbors(sce, colPairName = "knn_interaction_graph",
                          aggregate_by = "metadata", count_by = L2_COL)

cn_km <- kmeans(sce$aggregatedNeighbors, centers = Nval, nstart = 50, iter.max = 500)

# ==============================================================================
# 4. Assign CN labels
# ==============================================================================
cn <- unname(setNames(cn_km$cluster, colnames(sce))[colnames(obj)])
if (Kval == 15 && Nval == 10) {
  cn <- unname(CN_REMAP[as.character(cn)])
}
obj$cn_celltypes1 <- factor(cn, levels = 1:Nval)
obj$cn_label      <- paste0("CN", cn)

cn_meta <- data.frame(barcode = colnames(obj), sample = obj@meta.data[[SAMPLE_COL]],
                      cell_type = factor(obj@meta.data[[L2_COL]], levels = Col_order),
                      cn = as.integer(cn))
write.csv(cn_meta, file.path(neighborhood_dir, paste0("CN_labels_", tag, ".csv")),
          row.names = FALSE)

# ==============================================================================
# 5. Hypergeometric enrichment test
# ==============================================================================
cell_count <- table(cn_meta$cn, cn_meta$cell_type)[as.character(1:Nval), Col_order]
cell_count <- matrix(as.numeric(cell_count), nrow = Nval,
                     dimnames = list(as.character(1:Nval), Col_order))
n_total    <- sum(cell_count)

p_value_matrix <- matrix(NA, Nval, length(Col_order), dimnames = dimnames(cell_count))
for (i in rownames(cell_count)) {
  for (j in Col_order) {
    p_value_matrix[i, j] <- phyper(cell_count[i, j] - 1, sum(cell_count[, j]),
                                   n_total - sum(cell_count[, j]), sum(cell_count[i, ]),
                                   lower.tail = FALSE)
  }
}
long_df <- reshape2::melt(p_value_matrix)
colnames(long_df) <- c("CN", "Cell", "pval")
long_df$fdr <- p.adjust(long_df$pval, method = "BH")
write.table(long_df, file.path(neighborhood_dir,
                               paste0("CN_Enrichment_Pvalue_", tag, "_BH_Adjusted.txt")),
            sep = "\t", quote = FALSE)

fdr_matrix <- matrix(long_df$fdr, Nval, dimnames = dimnames(p_value_matrix))
oe_ratio   <- cell_count / (outer(rowSums(cell_count), colSums(cell_count)) / n_total)

# ==============================================================================
# 6. Heatmap of cell type enrichment per CN with significance
# ==============================================================================
prop_mat <- sweep(cell_count, 1, rowSums(cell_count), "/")
sig_star <- ifelse(oe_ratio >= 1 & fdr_matrix < FDR_THRESHOLD, "*", "")
row_lab  <- if (Kval == 15 && Nval == 10) {
  sprintf("CN%s: %s", rownames(prop_mat), CN_LABELS[rownames(prop_mat)])
} else paste0("CN", rownames(prop_mat))

P1 <- pheatmap(prop_mat, color = colorRampPalette(c("dark blue", "white", "dark red"))(100),
               scale = "column", cluster_rows = TRUE, cluster_cols = TRUE,
               display_numbers = sig_star, number_color = "black", fontsize_number = 18,
               labels_row = row_lab,
               main = paste0(tag, " -- * FDR < ", FDR_THRESHOLD, " (enrichment)"))
ggsave(as.ggplot(P1), filename = file.path(FIG_DIR, paste0("CN_Enrichment_Heatmap_", tag, ".pdf")),
       width = 12.5, height = 8)

# ==============================================================================
# 7. Frequency plots
# ==============================================================================
cn_meta$CN <- factor(paste0("CN", cn_meta$cn), levels = paste0("CN", 1:Nval))
cn_cols    <- if (Nval == 10) CN_COLS else setNames(scales::hue_pal()(Nval), paste0("CN", 1:Nval))

cc_ct <- as.data.frame(table(CellType = cn_meta$cell_type, CN = cn_meta$CN))
pp1 <- ggplot(cc_ct, aes(x = CN, y = Freq, fill = CellType)) +
  geom_bar(position = "fill", stat = "identity") + theme_bw() +
  labs(x = "CNs", y = "Cell Type") + scale_fill_manual(values = L2_COLS) +
  RotatedAxis()
ggsave(pp1, filename = file.path(FIG_DIR, paste0("FreqPlot_CNsByCellType_", tag, ".pdf")),
       width = 9, height = 12)

cn_meta$sample_label <- paste0(PATIENT_LABELS[cn_meta$sample], ":", sub("^LBCL_", "", cn_meta$sample))
cn_meta$sample_label <- factor(cn_meta$sample_label,
                               levels = unique(cn_meta$sample_label[order(cn_meta$sample)]))
cc_sample <- as.data.frame(table(Sample = cn_meta$sample_label, CN = cn_meta$CN))
pp2 <- ggplot(cc_sample, aes(x = Sample, y = Freq, fill = CN)) +
  geom_bar(position = "fill", stat = "identity") + theme_bw() +
  labs(x = "Sample", y = "CNs") + scale_fill_manual(values = cn_cols) +
  RotatedAxis()
ggsave(pp2, filename = file.path(FIG_DIR, paste0("FreqPlot_CNsBySample_", tag, ".pdf")),
       width = 9, height = 12)

# ==============================================================================
# 8. CN proportions per sample: CR vs PD
# ==============================================================================
samp_pct <- as.data.frame(table(sample = cn_meta$sample, CN = cn_meta$CN),
                          stringsAsFactors = FALSE) %>%
  group_by(sample) %>% mutate(pct = 100 * Freq / sum(Freq)) %>% ungroup()
samp_pct$CN <- factor(samp_pct$CN, levels = paste0("CN", 1:Nval))
samp_pct$Response <- factor(RESPONSE_MAP[samp_pct$sample], levels = c("PD", "CR"))

cn_stats <- samp_pct %>% group_by(CN) %>%
  summarise(median_CR = median(pct[Response == "CR"]),
            median_PD = median(pct[Response == "PD"]),
            p_wilcox  = wilcox.test(pct[Response == "CR"], pct[Response == "PD"],
                                    exact = FALSE)$p.value,
            y = max(pct), .groups = "drop") %>%
  mutate(fdr_BH = p.adjust(p_wilcox, method = "BH"),
         label  = sprintf("p = %.3f\nFDR = %.3f", p_wilcox, fdr_BH))
write.csv(cn_stats, file.path(neighborhood_dir, paste0("CN_CR_vs_PD_Wilcoxon_", tag, ".csv")),
          row.names = FALSE)

pp3 <- ggplot(samp_pct, aes(Response, pct, fill = Response)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.9) +
  geom_text(data = cn_stats, aes(x = 1.5, y = y * 1.1, label = label),
            inherit.aes = FALSE, size = 3, lineheight = 0.95) +
  facet_wrap(~ CN, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = RESPONSE_COLS) +
  labs(x = NULL, y = "Within sample CN %") + theme_bw(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(size = 11, face = "bold"))
ggsave(pp3, filename = file.path(FIG_DIR, paste0("CN_CR_vs_PD_", tag, ".pdf")),
       width = 12, height = 7)
