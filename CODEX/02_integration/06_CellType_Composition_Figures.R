library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(grid)

BASE_DIR   <- "/path/to/your/"
SCRIPT_DIR <- "/path/to/your/Scripts"
source(file.path(SCRIPT_DIR, "Utilities", "CODEX_Functions.R"))

L1_MARKERS <- c("CD19", "CD20", "KI67", "CD11B", "CD11C", "CD68", "CD14", "CD16",
                "CD15", "CD34", "Vimentin", "Collagen IV", "CD3e", "CD4", "CD45RO",
                "CD8", "CD27", "CD25", "FOXP3", "CD2")
T_HEATMAP_MARKERS <- c("CD4", "CD8", "Vista", "CCR7", "CD2", "CD5", "HLA-A", "PD-1",
                       "CD69", "CD27", "CD25", "CXCR3", "CD45RO", "TOX", "TIM3",
                       "TBET", "GRZB", "FOXP3", "KI67")
MYELOID_HEATMAP_MARKERS <- c("CD14", "CD16", "CD68", "HLA-DR", "CD11C", "HLA-A",
                             "Serpin B9")

SEURAT_DIR  <- file.path(BASE_DIR, "SeuratObj")
FIG_DIR     <- file.path(BASE_DIR, "Results", "Figures")
RESULTS_DIR <- file.path(BASE_DIR, "Results", "CellType_Composition")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

final_obj_path <- file.path(SEURAT_DIR, "DLBCL_CODEX_Final.rds")

# ==============================================================================
# 1. Load Final Annotated Object
# ==============================================================================

obj <- readRDS(final_obj_path)
DefaultAssay(obj) <- "CODEX"
obj@meta.data <- obj@meta.data[, !colnames(obj@meta.data) %in% rownames(obj[["CODEX"]])]

l1      <- as.character(obj@meta.data[[L1_COL]])
l2      <- as.character(obj@meta.data[[L2_COL]])
samples <- as.character(obj@meta.data[[SAMPLE_COL]])
print(table(l1)); print(table(l2))

# ==============================================================================
# 2. Level-1 UMAP
# ==============================================================================

p <- DimPlot(obj, reduction = "full.umap", group.by = L1_COL, cols = L1_COLS,
             label = TRUE, repel = TRUE, label.size = 4,
             raster = TRUE) +
  ggtitle(paste0("CODEX cell types, level 1 (N = ", format(ncol(obj), big.mark = ","), ")")) +
  coord_fixed() + NoAxes() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 10))
ggsave(p, filename = file.path(FIG_DIR, "L1_UMAP.pdf"), width = 16, height = 12)

# ==============================================================================
# 3. Level-1 Marker Heatmap
# ==============================================================================

obj <- JoinLayers(obj)

avg_l1 <- AverageExpression(obj, features = L1_MARKERS, group.by = L1_COL,
                            assays = "CODEX", layer = "data")$CODEX
colnames(avg_l1) <- gsub("-", "_", colnames(avg_l1))
avg_l1 <- as.matrix(avg_l1[L1_MARKERS, gsub("-", "_", L1_ORDER)])
colnames(avg_l1) <- L1_ORDER
z_l1 <- row_zscore(avg_l1)
write.csv(z_l1, file.path(RESULTS_DIR, "L1_heatmap_zscore.csv"))

pdf(file.path(FIG_DIR, "L1_Heatmap.pdf"), width = 9, height = 12)
draw(plot_zscore_heatmap(z_l1, legend_title = "Normalized Expression"),
     heatmap_legend_side = "right", padding = unit(c(2, 15, 2, 2), "mm"))
dev.off()

# ==============================================================================
# 4. Composition Barplots
# ==============================================================================

l1_cats      <- sort(names(L1_COLS))
t_cats       <- c("CD4 Naive-like", "CD4 TEM", "CD4T effector-exhausted", "CD8 TEM",
                  "CD8T effector-exhausted", "Proliferating T", "T-regulatory")
myeloid_cats <- c("Macrophage", "CD14+TIMs", "CD14+CD16+TIMs", "cDC")

comp_l1      <- make_comp(l1, samples, l1_cats)
comp_t       <- make_comp(l2, samples, t_cats)
comp_myeloid <- make_comp(l2, samples, myeloid_cats)
write.csv(comp_l1,      file.path(RESULTS_DIR, "L1_composition.csv"),      row.names = FALSE)
write.csv(comp_t,       file.path(RESULTS_DIR, "Tcell_composition.csv"),   row.names = FALSE)
write.csv(comp_myeloid, file.path(RESULTS_DIR, "Myeloid_composition.csv"), row.names = FALSE)

pdf(file.path(FIG_DIR, "Composition_Barplots.pdf"), width = 14, height = 10)
print(plot_PDvsCR(comp_l1, L1_COLS, "Level-1 cell types", "Percentage of all cells (%)"))
print(plot_PDvsCR(comp_t, T_COLS, "T cell subsets", "Percentage of T cells (%)"))
print(plot_PDvsCR(comp_myeloid, MYELOID_COLS, "Myeloid subsets",
                  "Percentage of myeloid cells (%)"))
dev.off()

# ==============================================================================
# 5. CR vs PD Wilcoxon Rank-Sum Tests
# ==============================================================================

run_wilcoxon_block <- function(cell_types, tag, denom_lbl) {
  tab <- as.data.frame(table(cell_type = factor(l2[l2 %in% cell_types], levels = cell_types),
                             sample = samples[l2 %in% cell_types]), stringsAsFactors = FALSE)
  tab <- tab %>% group_by(sample) %>% mutate(proportion = Freq / sum(Freq)) %>% ungroup()
  tab$Response <- factor(RESPONSE_MAP[tab$sample], levels = c("PD", "CR"))

  res <- do.call(rbind, lapply(cell_types, function(ct) {
    s  <- tab[tab$cell_type == ct, ]
    pd <- s$proportion[s$Response == "PD"]
    cr <- s$proportion[s$Response == "CR"]
    wt <- tryCatch(wilcox.test(pd, cr, exact = TRUE),
                   error = function(e) wilcox.test(pd, cr, exact = FALSE))
    data.frame(cell_type = ct, mean_PD = mean(pd), mean_CR = mean(cr),
               p_value = wt$p.value)
  }))
  res$p_adj <- p.adjust(res$p_value, method = "BH")
  write.csv(res, file.path(RESULTS_DIR, paste0("Wilcoxon_", tag, ".csv")), row.names = FALSE)

  lab <- tab %>% group_by(cell_type) %>% summarise(y = 1.15 * max(proportion)) %>%
    left_join(res, by = "cell_type") %>%
    mutate(label = paste0("p = ", signif(p_value, 3), "\np.adj = ", signif(p_adj, 3)))
  tab$cell_type <- factor(tab$cell_type, levels = cell_types)
  lab$cell_type <- factor(lab$cell_type, levels = cell_types)

  p <- ggplot(tab, aes(x = Response, y = proportion, fill = Response)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(aes(color = Response), width = 0.15, size = 2, alpha = 0.8) +
    facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = RESPONSE_COLS) +
    scale_color_manual(values = RESPONSE_COLS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    geom_text(data = lab, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE,
              size = 3, fontface = "italic", lineheight = 0.9) +
    labs(x = "Response", y = paste0("Proportion (of ", denom_lbl, ")")) +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(size = 9, face = "bold"), legend.position = "bottom")
  ggsave(p, filename = file.path(FIG_DIR, paste0("Wilcoxon_", tag, ".pdf")),
         width = 12, height = 3 * ceiling(length(cell_types) / 4) + 2)
  res
}

run_wilcoxon_block(myeloid_cats, "Myeloid_within_Myeloid", "myeloid cells")
run_wilcoxon_block(t_cats, "Tcell_within_Tcell", "T cells")

# ==============================================================================
# 6. T Cell and Myeloid Marker Heatmaps
# ==============================================================================

avg_all <- AverageExpression(obj, features = T_HEATMAP_MARKERS, group.by = L2_COL,
                             assays = "CODEX", layer = "data")$CODEX
z_t <- row_zscore(as.matrix(avg_all))[T_HEATMAP_MARKERS,
                                      c("CD4 Naive-like", "CD4 TEM", "CD4T effector-exhausted",
                                        "CD8 TEM", "CD8T effector-exhausted",
                                        "T-regulatory", "Proliferating T")]
write.csv(z_t, file.path(RESULTS_DIR, "Tcell_heatmap_Z.csv"))

data_mat <- GetAssayData(obj, assay = "CODEX", layer = "data")[MYELOID_HEATMAP_MARKERS, ]
mean_all <- sapply(sort(unique(l2)), function(g)
  Matrix::rowMeans(data_mat[, which(l2 == g), drop = FALSE]))
z_m <- row_zscore(mean_all)[MYELOID_HEATMAP_MARKERS, MYELOID_ORDER]
write.csv(z_m, file.path(RESULTS_DIR, "Myeloid_heatmap_Z.csv"))

pdf(file.path(FIG_DIR, "Tcell_Heatmap.pdf"), width = 6.5, height = 10)
draw(plot_zscore_heatmap(z_t, legend_title = "Z-score"),
     heatmap_legend_side = "right", padding = unit(c(2, 15, 2, 2), "mm"))
dev.off()

pdf(file.path(FIG_DIR, "Myeloid_Heatmap.pdf"), width = 5.2, height = 4.4)
draw(plot_zscore_heatmap(z_m, legend_title = "Z-score",
                         row_fontsize = 8, col_fontsize = 8),
     heatmap_legend_side = "right", padding = unit(c(2, 15, 2, 2), "mm"))
dev.off()
