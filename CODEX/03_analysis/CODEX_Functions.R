# ==============================================================================
# Helper functions for the CODEX analysis
# ==============================================================================

# -- Metadata columns of the final annotated object ----------------------------
SAMPLE_COL <- "orig.ident"
L1_COL     <- "cell_type_final_L1_reconciled"          # level 1 (8 cell types)
L2_COL     <- "cell_type_final_subtypes_myeloid5v2"    # level 2 (15 cell types)

# -- Clinical response and patient labels --------------------------------------
RESPONSE_MAP <- c(
  "DLBCL_10858"    = "PD", "DLBCL_25913A1"   = "PD", "DLBCL_34774"   = "PD",
  "DLBCL_402641A"  = "PD", "DLBCL_CS207490"  = "PD",
  "DLBCL_184419A"  = "CR", "DLBCL_31361C"    = "CR", "DLBCL_38979"   = "CR",
  "DLBCL_SP206535" = "CR", "DLBCL_SP2223235" = "CR"
)
PATIENT_LABELS <- c(
  "DLBCL_10858"    = "PT14", "DLBCL_25913A1"   = "PT10", "DLBCL_34774" = "PT9",
  "DLBCL_402641A"  = "PT12", "DLBCL_CS207490"  = "PT11",
  "DLBCL_184419A"  = "PT8",  "DLBCL_31361C"    = "PT4",  "DLBCL_38979" = "PT6",
  "DLBCL_SP206535" = "PT7",  "DLBCL_SP2223235" = "PT5"
)
SAMPLE_ORDER  <- c(sort(names(RESPONSE_MAP)[RESPONSE_MAP == "PD"]),
                   sort(names(RESPONSE_MAP)[RESPONSE_MAP == "CR"]))
RESPONSE_COLS <- c("PD" = "#E7A75E", "CR" = "#4FB7C5")

# -- Cell type palettes --------------------------------------------------------
L1_COLS <- c(
  "Lymphoma"     = "#D9D9D9", "Myeloid_cell" = "#FF9900",
  "Granulocytes" = "#E7298A", "T-regulatory" = "#984EA3",
  "Endothelial"  = "#00D6B6", "Stromal"      = "#CFC691",
  "CD8+ T-Cells" = "#19802C", "CD4+ T-Cells" = "#5391F5"
)
L1_ORDER <- c("Lymphoma", "Myeloid_cell", "Granulocytes", "Endothelial",
              "Stromal", "CD4+ T-Cells", "CD8+ T-Cells", "T-regulatory")

T_COLS <- c(
  "CD4 Naive-like"          = "#1F9CDA", "CD4 TEM"                 = "#112EB0",
  "CD4T effector-exhausted" = "#917EE6", "CD8 TEM"                 = "#5BAE03",
  "CD8T effector-exhausted" = "#6FF588", "T-regulatory"            = "#984EA3",
  "Proliferating T"         = "#808080"
)
T_ORDER <- names(T_COLS)

MYELOID_COLS <- c(
  "CD14+TIMs" = "#654321", "CD14+CD16+TIMs" = "#FF9007",
  "Macrophage" = "#F9D448", "cDC" = "#E94310"
)
MYELOID_ORDER <- names(MYELOID_COLS)

L2_COLS <- c(
  "Lymphoma" = "#D9D9D9", "Macrophage" = "#F9D448", "CD14+TIMs" = "#F2A32B",
  "CD14+CD16+TIMs" = "#BE8C06", "cDC" = "#E94310", "Granulocytes" = "#E7298A",
  T_COLS, "Endothelial" = "#00D6B6", "Stromal" = "#CFC691"
)

# -- Cellular neighborhoods ----------------------------------------------------
CN_LABELS <- c(
  "1" = "Myeloid/Tex",      "2" = "Endothelial rich",
  "3" = "Granulocyte/Stromal rich", "4" = "Lymphoma_1",
  "5" = "CD4T/cDC_1",       "6" = "CD8T/Myeloid",
  "7" = "Macrophage rich",  "8" = "TIM rich",
  "9" = "CD4T/cDC_2",       "10" = "Lymphoma_2"
)
CN_COLS <- c(
  "CN1" = "#FB8072", "CN2" = "#8DD3C7", "CN3" = "#BC80BD", "CN4" = "#D9D9D9",
  "CN5" = "#BEBADA", "CN6" = "#80B1D3", "CN7" = "#FF7F00", "CN8" = "#FDB462",
  "CN9" = "#E7298A", "CN10" = "#969696"
)

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# ==============================================================================
# Sketch-based RPCA re-integration of a cell subset
# ==============================================================================

reintegrate_sketch_rpca <- function(obj, markers, ncells = 200000, npcs = 15,
                                    reference = "DLBCL_34774", merge_small = TRUE) {
  DefaultAssay(obj) <- "CODEX"
  if ("sketch" %in% Assays(obj)) obj[["sketch"]] <- NULL
  obj <- DietSeurat(obj, assays = "CODEX", layers = c("counts", "data"),
                    dimreducs = NULL, graphs = NULL)
  if ("scale.data" %in% Layers(obj)) obj[["CODEX"]]$scale.data <- NULL

  if (merge_small && "DLBCL_SP2223235" %in% names(which(table(obj$orig.ident) < 200))) {
    obj$orig.ident[obj$orig.ident %in% c("DLBCL_SP206535", "DLBCL_SP2223235")] <-
      "DLBCL_SP206535_SP2223235"
  }

  obj <- JoinLayers(obj)
  obj[["CODEX"]] <- split(obj[["CODEX"]], f = obj$orig.ident)
  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 1)
  obj <- ScaleData(obj, features = markers, verbose = FALSE)
  obj <- RunPCA(obj, npcs = npcs, features = markers, verbose = FALSE)

  VariableFeatures(obj) <- rownames(obj)
  obj <- SketchData(object = obj, assay = "CODEX", ncells = ncells,
                    method = "LeverageScore", sketched.assay = "sketch",
                    over.write = TRUE)

  DefaultAssay(obj) <- "sketch"
  VariableFeatures(obj) <- rownames(obj)
  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 1)
  obj <- ScaleData(obj, features = markers)
  obj <- RunPCA(obj, npcs = npcs, features = markers, reduction.name = "unint.pca")

  reference_idx <- which(Layers(obj, search = "data") == paste0("data.", reference))
  obj <- IntegrateLayers(obj, method = RPCAIntegration, features = markers,
                         orig = "unint.pca", new.reduction = "integrated.rpca",
                         dims = 1:npcs, k.anchor = 20, reference = reference_idx)
  obj <- ProjectIntegration(object = obj, features = markers,
                            sketched.assay = "sketch", assay = "CODEX",
                            reduction = "integrated.rpca")
  obj
}

cluster_sketch_project <- function(obj, markers, res, npcs = 15) {
  DefaultAssay(obj) <- "sketch"
  obj <- JoinLayers(obj)
  obj <- ScaleData(obj, features = markers)
  obj <- RunPCA(obj, features = markers, approx = FALSE)
  obj <- FindNeighbors(obj, features = markers, k.param = 30)
  obj <- FindClusters(obj, algorithm = 2, resolution = res, cluster.name = "sketch_clusters")
  obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:npcs, return.model = TRUE)
  obj[["sketch"]] <- split(obj[["sketch"]], f = obj$orig.ident)
  obj <- ProjectData(object = obj, assay = "CODEX",
                     full.reduction = "integrated.rpca.full",
                     sketched.assay = "sketch", sketched.reduction = "integrated.rpca",
                     umap.model = "umap", dims = 1:npcs,
                     refdata = list(clusters_full = "sketch_clusters"))
  setNames(as.character(obj$clusters_full), colnames(obj))
}

reintegrate_cluster_nosketch <- function(obj, markers, res, merge_samples,
                                         npcs = 15, reference = "DLBCL_34774") {
  DefaultAssay(obj) <- "CODEX"
  obj <- DietSeurat(obj, assays = "CODEX", layers = c("counts", "data"),
                    dimreducs = NULL, graphs = NULL)
  obj$orig.ident[obj$orig.ident %in% merge_samples] <- paste(merge_samples, collapse = "_")
  obj <- JoinLayers(obj)
  obj[["CODEX"]] <- split(obj[["CODEX"]], f = obj$orig.ident)
  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 1)
  obj <- ScaleData(obj, features = markers, verbose = FALSE)
  obj <- RunPCA(obj, npcs = npcs, features = markers, verbose = FALSE)
  reference_idx <- which(Layers(obj, search = "data") == paste0("data.", reference))
  obj <- IntegrateLayers(obj, method = RPCAIntegration, features = markers,
                         orig.reduction = "pca", new.reduction = "integrated.rpca",
                         dims = 1:npcs, k.anchor = 20, reference = reference_idx)
  obj <- JoinLayers(obj)
  obj <- ScaleData(obj, features = markers)
  obj <- FindNeighbors(obj, features = markers, k.param = 30)
  obj <- FindClusters(obj, algorithm = 2, resolution = res)
  setNames(as.character(Idents(obj)), colnames(obj))
}

# ==============================================================================
# Plot helpers
# ==============================================================================

make_comp <- function(values, samples, keep_cats) {
  df  <- data.frame(Sample = samples, CellType = values)
  df  <- df[df$CellType %in% keep_cats, ]
  tab <- as.data.frame(table(Sample = df$Sample, CellType = df$CellType),
                       stringsAsFactors = FALSE)
  colnames(tab) <- c("Sample", "CellType", "Count")
  tab <- tab %>% group_by(Sample) %>%
    mutate(Proportion = Count / sum(Count) * 100) %>% ungroup() %>% as.data.frame()
  tab$Proportion[is.na(tab$Proportion)] <- 0
  tab$CellType <- factor(tab$CellType, levels = keep_cats)
  tab
}

plot_PDvsCR <- function(comp, cols, title, ylab) {
  comp <- comp[comp$Sample %in% SAMPLE_ORDER, ]
  comp$Sample <- factor(PATIENT_LABELS[comp$Sample],
                        levels = PATIENT_LABELS[SAMPLE_ORDER])
  npd  <- sum(RESPONSE_MAP[SAMPLE_ORDER] == "PD")
  ncr  <- sum(RESPONSE_MAP[SAMPLE_ORDER] == "CR")
  yseg <- -0.085 * 100
  ytxt <- -0.155 * 100
  ggplot(comp, aes(x = Sample, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = cols, drop = FALSE) +
    labs(title = title, x = NULL, y = ylab, fill = "Cell Type") +
    geom_vline(xintercept = npd + 0.5, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    annotate("segment", x = 0.6, xend = npd + 0.4, y = yseg, yend = yseg,
             colour = RESPONSE_COLS["PD"], linewidth = 1.4) +
    annotate("segment", x = npd + 0.6, xend = npd + ncr + 0.4, y = yseg, yend = yseg,
             colour = RESPONSE_COLS["CR"], linewidth = 1.4) +
    annotate("text", x = (npd + 1) / 2, y = ytxt, label = "PD",
             colour = RESPONSE_COLS["PD"], fontface = "bold", size = 5) +
    annotate("text", x = npd + (ncr + 1) / 2, y = ytxt, label = "CR",
             colour = RESPONSE_COLS["CR"], fontface = "bold", size = 5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title  = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 9),
          plot.margin = margin(t = 8, r = 12, b = 60, l = 12)) +
    coord_cartesian(clip = "off")
}

row_zscore <- function(m) {
  z <- t(scale(t(m)))
  z[is.na(z)] <- 0
  z[is.infinite(z)] <- 0
  z
}

plot_zscore_heatmap <- function(mz, legend_title = "Z-score", title = NULL,
                                row_fontsize = 10, col_fontsize = 11) {
  Heatmap(mz, name = "Z-score",
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          cluster_rows = FALSE, cluster_columns = FALSE,
          row_names_gp = gpar(fontsize = row_fontsize),
          column_names_gp = gpar(fontsize = col_fontsize),
          column_names_rot = 45, rect_gp = gpar(col = "grey80", lwd = 0.5),
          column_title = title,
          column_title_gp = gpar(fontsize = 12, fontface = "bold"),
          heatmap_legend_param = list(title = legend_title, at = c(-2, -1, 0, 1, 2)))
}
