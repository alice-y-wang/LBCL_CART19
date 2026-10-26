library(Seurat)

BASE_DIR   <- "/path/to/your/"
SCRIPT_DIR <- "/path/to/your/Scripts"
source(file.path(SCRIPT_DIR, "Utilities", "CODEX_Functions.R"))

final_obj_path <- file.path(BASE_DIR, "SeuratObj", "DLBCL_CODEX_Final.rds")
cn_labels_csv  <- file.path(BASE_DIR, "Results", "Neighborhood_Files", "CN_labels.csv")
OUTPUT_DIR     <- file.path(BASE_DIR, "Results", "Napari_Inputs")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

CELLTYPE_CODES <- setNames(1:15, c("Lymphoma", "Macrophage", "CD14+TIMs", "CD14+CD16+TIMs",
                                   "cDC", "Granulocytes", "T-regulatory", "CD4 Naive-like",
                                   "CD4 TEM", "CD4T effector-exhausted", "CD8 TEM",
                                   "CD8T effector-exhausted", "Proliferating T",
                                   "Endothelial", "Stromal"))

# ==============================================================================
# 1. Per-Cell Labels
# ==============================================================================

meta <- readRDS(final_obj_path)@meta.data
cn   <- read.csv(cn_labels_csv)

df <- data.frame(barcode        = rownames(meta),
                 sample         = meta[[SAMPLE_COL]],
                 CellID         = meta$CellID,
                 x.coord        = meta$x.coord,
                 y.coord        = meta$y.coord,
                 cell_type      = as.character(meta[[L2_COL]]),
                 cell_type_code = unname(CELLTYPE_CODES[as.character(meta[[L2_COL]])]))
df$cn_code <- cn$cn[match(df$barcode, cn$barcode)]

write.csv(data.frame(code = CELLTYPE_CODES, cell_type = names(CELLTYPE_CODES)),
          file.path(OUTPUT_DIR, "CellType_codes.csv"), row.names = FALSE)
write.csv(data.frame(code = as.integer(names(CN_LABELS)), CN = unname(CN_LABELS)),
          file.path(OUTPUT_DIR, "CN_codes.csv"), row.names = FALSE)

# ==============================================================================
# 2. Write One Table per Sample
# ==============================================================================

for (s in sort(unique(df$sample))) {
  write.csv(df[df$sample == s, c("CellID", "x.coord", "y.coord", "cell_type",
                                 "cell_type_code", "cn_code")],
            file.path(OUTPUT_DIR, paste0(s, "_labels.csv")), row.names = FALSE)
}
