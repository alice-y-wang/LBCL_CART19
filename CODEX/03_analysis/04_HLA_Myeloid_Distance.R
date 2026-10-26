library(Seurat)
library(dplyr)
library(ggplot2)

BASE_DIR   <- "/path/to/your/"
SCRIPT_DIR <- "/path/to/your/Scripts"
source(file.path(SCRIPT_DIR, "Utilities", "CODEX_Functions.R"))

SEURAT_DIR   <- file.path(BASE_DIR, "SeuratObj")
DISTANCE_DIR <- file.path(BASE_DIR, "Distance_Analysis", "Distance_Results")
HLA_DIR      <- file.path(DISTANCE_DIR, "results", "HLA_myeloid_distance")
FIG_DIR      <- file.path(BASE_DIR, "Results", "Figures")
dir.create(HLA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

final_obj_path <- file.path(SEURAT_DIR, "DLBCL_CODEX_Final.rds")

MYELOID_TYPES <- c("Macrophage", "CD14+TIMs", "CD14+CD16+TIMs")
MARKERS       <- c("HLA-A", "HLA-DR")
TARGET        <- "CD8 TEM"
MIN_CELLS     <- 10

# ==============================================================================
# 1. HLA-A / HLA-DR Expression of Myeloid Cells
# ==============================================================================

obj  <- readRDS(final_obj_path)
meta <- obj@meta.data
mac_cells <- rownames(meta)[meta[[L2_COL]] %in% MYELOID_TYPES]

data_layers <- grep("^data", Layers(obj[["CODEX"]]), value = TRUE)
expr <- do.call(cbind, lapply(data_layers, function(l)
  LayerData(obj[["CODEX"]], layer = l)[MARKERS, , drop = FALSE]))
mac_cells <- intersect(mac_cells, colnames(expr))

hla <- data.frame(CellID   = mac_cells,
                  sample   = as.character(meta[mac_cells, SAMPLE_COL]),
                  HLA_A    = as.numeric(expr["HLA-A", mac_cells]),
                  HLA_DR   = as.numeric(expr["HLA-DR", mac_cells]))
hla$Response <- RESPONSE_MAP[hla$sample]
rm(obj, meta, expr); gc()

# ==============================================================================
# 2. High vs Low
# ==============================================================================

dist <- readRDS(file.path(DISTANCE_DIR, paste0("Distance_ByCell_", safe_name(TARGET), ".rds")))
dist <- dist[dist$source %in% MYELOID_TYPES, ]

stats <- list()
plots <- list()
for (mk in MARKERS) {
  col      <- gsub("-", "_", mk)
  med      <- median(hla[[col]])                   
  lab_low  <- paste(mk, "Low")
  lab_high <- paste(mk, "High")
  level    <- setNames(ifelse(hla[[col]] > med, lab_high, lab_low), hla$CellID)

  dd <- dist
  dd$Level    <- factor(level[dd$CellID], levels = c(lab_low, lab_high))
  dd$Response <- factor(RESPONSE_MAP[dd$sample], levels = c("CR", "PD"))
  dd <- dd[!is.na(dd$Level), ]

  dd <- dd %>% group_by(Response, Level) %>%
    filter(distance_um <= quantile(distance_um, 0.99)) %>% ungroup()

  for (rs in c("CR", "PD")) {
    lo <- dd$distance_um[dd$Response == rs & dd$Level == lab_low]
    hi <- dd$distance_um[dd$Response == rs & dd$Level == lab_high]
    alt <- if (rs == "CR") "greater" else "less"
    stats[[length(stats) + 1]] <- data.frame(
      marker = mk, Response = rs, test = ifelse(rs == "CR", "High closer", "Low closer"),
      n_low = length(lo), n_high = length(hi),
      median_low = median(lo), median_high = median(hi),
      p_value = wilcox.test(lo, hi, alternative = alt)$p.value)
  }

  plots[[mk]] <- ggplot(dd, aes(Level, distance_um, fill = Level)) +
    geom_violin(scale = "width", color = NA) +
    geom_boxplot(width = 0.2, fill = NA, outlier.shape = NA) +
    facet_wrap(~ Response) +
    scale_fill_manual(values = setNames(c("steelblue", "#DC143C"), c(lab_low, lab_high))) +
    theme_bw() + guides(fill = "none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = paste(mk, "level"), y = paste0("Min distance to ", TARGET, " (um)"))
}

stats <- do.call(rbind, stats)
stats$p_adj <- ave(stats$p_value, stats$Response, FUN = function(p) p.adjust(p, method = "BH"))
write.csv(stats, file.path(HLA_DIR, "Myeloid_HLA_Distance_CD8TEM_byResponse.csv"), row.names = FALSE)
print(stats)

pdf(file.path(FIG_DIR, "HLA_Myeloid_Distance_CD8TEM.pdf"), width = 8, height = 6)
for (mk in MARKERS) {
  s <- stats[stats$marker == mk, ]
  print(plots[[mk]] + ggtitle(paste0(mk, ": CR ", s$test[1], " p=", signif(s$p_value[1], 2),
                                     " | PD ", s$test[2], " p=", signif(s$p_value[2], 2))))
}
dev.off()
