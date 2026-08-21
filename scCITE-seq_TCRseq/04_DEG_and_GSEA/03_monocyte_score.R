################################################################################
## Module scoring of the interferon / TNFa / TGFb signatures in the monocyte
## compartment
################################################################################

# ---- Libraries --------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ggpubr)
library(ggrepel)

# ---- CONFIG -----------------------------------------------------------------
analysis_dir <- "."
setwd(analysis_dir)

objects_dir <- "seurat_objects"
out_root    <- file.path("results", "module_scoring_mono")
table_dir   <- file.path(out_root, "tables")
plot_dir    <- file.path(out_root, "plots")
traj_dir    <- file.path(out_root, "trajectories")

for (d in c(table_dir, plot_dir, traj_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

f_mono <- "all_samples_Mono_filtered.RDS"

celltype_col <- "cell.anno"
patient_col  <- "patient.id"

PAL        <- c(CR = "#4FB7C5", PD = "#d95f02")
timepoints <- c("APH", "D0", "Peak", "4W")   # plotting order

# signatures plotted as boxplots / violins
PLOT_SIGS <- c("CD14_TNFa", "CD14_tgfb")

# trajectory mapping: cell type -> module prefix, timepoint -> module suffix
CT_MAP    <- c("CD14 Mono" = "CD14", "CD16 Mono" = "CD16", "Dendritic Cell" = "DC")
TP_MAP    <- c(APH = "day0", D0 = "day0", Peak = "peak", `4W` = "week4")
TP_BREAKS <- c(APH = 0, D0 = 30, Peak = 37, `4W` = 58)   # modeled as days post-apheresis to get positive auc values.

set.seed(1)

# ---- Modules ----------------------------------------------------------------
modules <- list(
  CD14_peak = c(
    "MX1", "RSAD2", "IFITM1", "IFIT3", "IFITM3", "IFIT2", "ISG20", "EPSTI1",
    "IFI44L", "ISG15", "PSMB9", "SAMD9", "OASL", "CMPK2", "USP18", "OAS1",
    "LY6E", "IFI44", "SAMD9L", "LAP3", "IRF7", "DDX60", "GMPR", "RTP4",
    "PSME2", "IFIH1", "UBE2L6", "IFITM2", "IFI35", "NMI", "PLSCR1", "TRIM5",
    "TRIM21", "SELL", "PARP9", "EIF2AK2", "PNPT1", "DHX58", "IRF2", "PSMB8",
    "BST2", "GBP4", "IRF1", "B2M", "SP110", "TAP1", "TRAFD1", "PSMA3", "GBP2",
    "PSME1", "PARP14", "STAT2", "IFI30", "CASP1", "TRIM14", "HELZ2", "LPAR6",
    "IRF9", "PARP12", "NUB1", "IFIT1", "IRF4", "IDO1", "TNFSF10", "ZBP1",
    "MX2", "XAF1", "OAS3", "GZMA", "SERPING1", "OAS2", "TOR1B", "LYSMD2",
    "FCGR1A", "STAT1", "VAMP5", "FAS", "PSMB10", "HLA-DRB1", "FPR1", "GCH1",
    "APOL6", "PSMA2", "JAK2", "CD40", "FGL2", "ST8SIA4", "RBCK1", "SRI",
    "CASP4", "HERC5", "GBP1", "HLA-DQB1", "GBP5", "TRIM22", "FCGR1B",
    "HSPA1A", "RNASEL", "IFI6", "IFIT5", "TRIM34", "SP100"),
  
  CD14_day0 = c(
    "IFITM1", "RSAD2", "IFI44L", "USP18", "ISG15", "MX1", "IFIT3", "CMPK2",
    "IFI44", "ISG20", "GMPR", "IFIT2", "OASL", "LY6E", "DDX60", "IFITM2",
    "PNPT1", "EIF2AK2", "EPSTI1", "IFITM3", "DHX58", "OAS1", "HERC6", "IFI35",
    "IRF7", "PLSCR1", "TRIM14", "PARP9", "HELZ2", "LAP3", "IFIH1", "SAMD9",
    "SAMD9L", "PARP14", "PARP12", "TRIM5", "TENT5A", "MOV10", "TRIM25",
    "BST2", "SP110", "CMTR1", "UBE2L6", "SELL", "NMI", "IFIT1", "OAS3",
    "ZBP1", "CD69", "MX2", "OAS2", "SERPING1", "IRF4", "SOCS3", "TOR1B",
    "XAF1", "PIM1", "ARID5B", "TNFSF10", "GCH1", "RNF213", "SOD2", "IL15RA",
    "HLA-DRB1", "CASP3", "TDRD7", "ICAM1", "B2M", "PELI1", "PSMA3", "PFKP",
    "ADAR", "STAT3", "APOL6", "NFKBIA", "STAT2", "PSMB9", "ST8SIA4",
    "ST3GAL5", "HERC5", "IFI6", "HLA-DQB1", "TRIM22"),
  
  CD14_week4 = c(
    "USP18", "IFI44L", "IFITM1", "LY6E", "EPSTI1", "MX1", "IFIT2", "LGALS3BP",
    "IRF7", "IFITM3", "RSAD2", "IFITM2", "ISG15", "OASL", "ISG20", "CMPK2",
    "IFI44", "IFIT3", "IFI35", "PLSCR1", "EIF2AK2", "PSME2", "LAP3", "OAS1",
    "ZBP1", "IFIT1", "SERPING1", "CD40", "FCGR1A", "OAS2", "HLA-DRB1", "MX2",
    "XAF1", "UPP1", "VAMP5", "PML", "IFI6", "HLA-DQB1", "TUBB6", "GBP1",
    "HERC5", "FKBP5"),
  
  CD16_peak = c(
    "CXCL10", "RSAD2", "IFIT3", "IFI44L", "IFIT2", "MX1", "ISG20", "OASL",
    "IFITM1", "IFI44", "ISG15", "USP18", "SAMD9", "CMPK2", "EPSTI1", "OAS1",
    "SAMD9L", "IFIH1", "IFITM3", "PLSCR1", "GMPR", "EIF2AK2", "DDX60", "IRF7",
    "PARP9", "PSMB9", "LY6E", "IFI35", "UBE2L6", "LPAR6", "BST2", "PSME2",
    "LAP3", "NMI", "IFITM2", "TRIM21", "PSMA3", "IFI27", "IRF1", "SELL",
    "PARP14", "TAP1", "IFI30", "B2M", "IRF9", "SP110", "PSMB8", "IFIT1",
    "ZBP1", "XAF1", "TNFSF10", "MX2", "TOR1B", "OAS2", "OAS3", "SERPING1",
    "LYSMD2", "STAT1", "FPR1", "VAMP5", "PSMB10", "PSMA2", "FGL2", "HERC5",
    "GBP1", "IFI6", "TRIM22", "EGR1", "HSPA1A"),
  
  CD16_day0 = c(
    "RSAD2", "IFI44L", "USP18", "IFITM1", "ISG20", "ISG15", "OASL", "GMPR",
    "IFIT3", "CXCL10", "IFI44", "MX1", "CMPK2", "EPSTI1", "DHX58", "LAP3",
    "EIF2AK2", "DDX60", "PLSCR1", "IFIH1", "IFI35", "IFIT2", "PNPT1", "LY6E",
    "PARP9", "IL15", "IRF1", "HELZ2", "IRF7", "IFITM2", "SAMD9", "B2M",
    "TRIM14", "SAMD9L", "TENT5A", "PSMB9", "UBE2L6", "PARP14", "PARP12",
    "STAT2", "BST2", "IFIT1", "ZBP1", "OAS3", "SOCS3", "SERPING1", "PIM1",
    "MX2", "TOR1B", "OAS2", "XAF1", "HLA-DQA1", "HLA-DRB1", "ICAM1",
    "TNFSF10", "STAT1", "GCH1", "CASP3", "APOL6", "ST8SIA4", "HERC5", "IFI6",
    "HLA-DQB1"),
  
  DC_peak = c(
    "IFI44L", "IFIT3", "IFITM1", "MX1", "IFITM3", "ISG15", "SAMD9L", "OAS1",
    "PSMB9", "LY6E", "IFI35", "PSME2", "IFIT1", "XAF1", "TNFSF10", "STAT1",
    "HLA-DQA1", "IFI6"),
  
  DC_day0 = c(
    "RSAD2", "IFI44L", "USP18", "IFITM1", "MX1", "IFIT3", "ISG15", "ISG20",
    "IFIT2", "CMPK2", "DHX58", "OASL", "DDX60", "IFI44", "LY6E", "IFITM3",
    "EPSTI1", "IFI35", "PARP14", "IFITM2", "OAS1", "EIF2AK2", "LAP3", "IRF7",
    "IFIT1", "OAS3", "OAS2", "ZBP1", "SERPING1", "MX2", "PIM1", "XAF1",
    "HLA-DQA1", "TNFSF10", "HERC5", "IFI6", "TUBB6"),
  
  CD14_TNFa = c(
    "AREG", "SDC4", "G0S2", "IER3", "KLF9", "PER1", "JAG1", "ICOSLG", "ABCA1",
    "CLCF1", "BHLHE40", "HBEGF", "VEGFA", "KLF10", "ID2", "PFKFB3", "CD69",
    "CDKN1A", "IRS2", "IER5", "FOSL1", "BTG3", "GADD45A", "MAP2K3", "MAFF",
    "NR4A1", "RNF19B", "IL18", "SPHK1"),
  
  CD14_tgfb = c(
    "SMAD7", "FURIN", "ID1", "KLF10", "ID2", "THBS1", "CDKN1C", "ENG",
    "NCOR2", "CTNNB1", "SKIL", "ID3", "CDK9", "HIPK2", "TGIF1", "SKI",
    "SMURF1")
)


# 1. Score with AddModuleScore # ------------------

obj <- readRDS(file.path(objects_dir, f_mono))
DefaultAssay(obj) <- "RNA"

for (sig in names(modules)) {
  message("scoring ", sig)
  obj <- AddModuleScore(obj, features = modules[sig], name = paste0(sig, "_AMS"), seed = 1)
  # Seurat appends the index to the name, so rename back
  colnames(obj@meta.data)[colnames(obj@meta.data) == paste0(sig, "_AMS1")] <- paste0(sig, "_AMS")
}

meta <- obj@meta.data
meta$cell_bc <- rownames(meta)
rm(obj); gc(verbose = FALSE)

write_csv(meta, file.path(table_dir, "mono_module_scores.csv"))

sig_cols <- paste0(names(modules), "_AMS")

# ---- tidy for plotting ------------------------------------------------------
df <- meta %>%
  filter(response %in% names(PAL), timepoint %in% timepoints) %>%
  mutate(response  = factor(response, levels = names(PAL)),
         timepoint = factor(timepoint, levels = timepoints),
         celltype  = .data[[celltype_col]]) %>%
  filter(!is.na(celltype))

cat(nrow(df), "cells,", length(sig_cols), "signatures\n")

base_theme <- theme_classic(base_size = 12) +
  theme(plot.title   = element_text(face = "bold"),
        strip.text   = element_text(face = "bold"),
        axis.title   = element_text(face = "bold"),
        legend.position = "top",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))


# 2. Patient-level Wilcoxon # -------------

stats <- df %>%
  group_by(.data[[patient_col]], response, timepoint, celltype) %>%
  summarise(across(all_of(sig_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(all_of(sig_cols), names_to = "signature", values_to = "score") %>%
  group_by(signature, celltype, timepoint) %>%
  group_modify(~ {
    a <- .x$score[.x$response == "CR"]; a <- a[!is.na(a)]
    b <- .x$score[.x$response == "PD"]; b <- b[!is.na(b)]
    if (length(a) < 2 || length(b) < 2) return(tibble())
    wt <- suppressWarnings(wilcox.test(a, b))
    tibble(n_CR = length(a), n_PD = length(b),
           median_CR = median(a), median_PD = median(b),
           statistic = unname(wt$statistic), p.value = wt$p.value)
  }) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value)

write_csv(stats, file.path(table_dir, "patient_level_wilcox_celltype_x_timepoint.csv"))
print(head(stats, 20))


# TNFa and TGFb Boxplot + violin by cell.anno x timepoint # -------------

for (sig in PLOT_SIGS) {
  col <- paste0(sig, "_AMS")
  
  d <- df %>%
    select(response, timepoint, celltype, score = all_of(col)) %>%
    filter(!is.na(score))
  
  p_box <- ggplot(d, aes(response, score, fill = response)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.8, width = 0.6) +
    facet_grid(timepoint ~ celltype, scales = "free_y") +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
    scale_fill_manual(values = PAL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(title = paste0(col, " - boxplot"), x = NULL, y = "Module score") +
    base_theme
  
  p_vln <- ggplot(d, aes(response, score, fill = response)) +
    geom_violin(trim = FALSE, scale = "width", width = 0.8, alpha = 0.5,
                color = "grey30", linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
    facet_grid(timepoint ~ celltype, scales = "free_y") +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
    scale_fill_manual(values = PAL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(title = paste0(col, " - violin"), x = NULL, y = "Module score") +
    base_theme
  
  pdf(file.path(plot_dir, paste0("boxplot_violin_", col, "_celltype_x_timepoint.pdf")),
      width = 12, height = 9)
  print(p_box); print(p_vln)
  dev.off()
}


# Interferon family trajectory for celltype x timepoint # -----------


module_df <- df %>%
  select(all_of(c(patient_col, "timepoint", "response", "celltype")), all_of(sig_cols)) %>%
  filter(celltype %in% names(CT_MAP)) %>%
  mutate(time_numeric = unname(TP_BREAKS[as.character(timepoint)]))

ct <- unname(CT_MAP[as.character(module_df$celltype)])
tp <- unname(TP_MAP[as.character(module_df$timepoint)])
col <- paste0(ct, "_", tp, "_AMS")


miss <- !col %in% sig_cols
col[miss] <- sub("_week4_", "_peak_", col[miss])

M   <- as.matrix(module_df[, sig_cols, drop = FALSE])
idx <- match(col, sig_cols)
module_df$module_score <- ifelse(is.na(idx), NA_real_, M[cbind(seq_len(nrow(module_df)), idx)])
module_df <- filter(module_df, !is.na(module_score))

## median per patient x celltype x timepoint
ps <- module_df %>%
  group_by(.data[[patient_col]], response, celltype, time_numeric) %>%
  summarise(module_score = median(module_score, na.rm = TRUE), .groups = "drop")

plot_auc_traj <- function(ps, ttl, individual = TRUE) {
  p <- ggplot(ps, aes(time_numeric, module_score, color = response))
  if (individual) {
    p <- p +
      geom_line(aes(group = .data[[patient_col]]), alpha = 0.25) +
      geom_point(aes(group = .data[[patient_col]]), alpha = 0.4)
  }
  p +
    stat_summary(aes(group = response), fun = mean, geom = "line", linewidth = 1.5) +
    stat_summary(aes(group = response, fill = response), fun.data = mean_se,
                 geom = "ribbon", alpha = 0.2, color = NA) +
    scale_color_manual(values = PAL) +
    scale_fill_manual(values = PAL) +
    scale_x_continuous(
      breaks = unname(TP_BREAKS),
      labels = paste0(unname(TP_BREAKS), "\n", names(TP_BREAKS)),
      minor_breaks = NULL
    ) +
    facet_wrap(~celltype, scales = "free_y") +
    theme_classic(base_size = 14) +
    theme(panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4)) +
    labs(x = "Days post apheresis", y = "Module score", title = ttl)
}

pdf(file.path(traj_dir, "Trajectory_celltype_timepoint_AMS.pdf"), width = 15, height = 6)
print(plot_auc_traj(ps, "celltype x timepoint module", individual = TRUE))
print(plot_auc_traj(ps, "celltype x timepoint module", individual = FALSE))
dev.off()



# Interferon total normalized AUC per patient, same celltype x timepoint module# ------

# Trapezoidal area under each patient's trajectory, divided by the time span
# it covers, so patients with different numbers of timepoints stay comparable.
auc <- ps %>%
  arrange(.data[[patient_col]], celltype, time_numeric) %>%
  group_by(.data[[patient_col]], response, celltype) %>%
  filter(n() >= 2) %>%
  summarise(
    raw_auc_norm = sum(diff(time_numeric) *
                         (head(module_score, -1) + tail(module_score, -1)) / 2) /
      (max(time_numeric) - min(time_numeric)),
    n_timepoints = n(),
    .groups = "drop"
  )

write_csv(auc, file.path(table_dir, "patient_normalized_auc_celltype_timepoint.csv"))

# CR vs PD per cell type
auc_stats <- auc %>%
  group_by(celltype) %>%
  group_modify(~ {
    a <- .x$raw_auc_norm[.x$response == "CR"]; a <- a[!is.na(a)]
    b <- .x$raw_auc_norm[.x$response == "PD"]; b <- b[!is.na(b)]
    if (length(a) < 2 || length(b) < 2) return(tibble(p = NA_real_, n_CR = length(a), n_PD = length(b)))
    wt <- suppressWarnings(wilcox.test(a, b))
    tibble(p = wt$p.value, n_CR = length(a), n_PD = length(b))
  }) %>%
  ungroup()

write_csv(auc_stats, file.path(table_dir, "patient_normalized_auc_wilcox.csv"))

fmt_p <- function(p) ifelse(is.na(p), "ns", paste0("p = ", signif(p, 2)))

plot_auc_box <- function(auc, ttl, stats) {
  ann <- auc %>%
    group_by(celltype) %>%
    summarise(lo = min(raw_auc_norm, na.rm = TRUE),
              hi = max(raw_auc_norm, na.rm = TRUE), .groups = "drop") %>%
    left_join(stats, by = "celltype") %>%
    mutate(y = hi + 0.08 * (hi - lo), lab = fmt_p(p))
  
  ggplot(auc, aes(response, raw_auc_norm, fill = response)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
    geom_point(shape = 21, size = 2, stroke = 0.2, color = "black") +
    geom_text_repel(aes(label = .data[[patient_col]]), size = 2,
                    max.overlaps = 50, segment.color = "grey50") +
    geom_text(data = ann, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE, size = 3.2) +
    scale_fill_manual(values = PAL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    facet_wrap(~celltype, scales = "free_y") +
    theme_classic() +
    labs(x = "response", y = "Raw normalized AUC", title = ttl)
}

pdf(file.path(traj_dir, "AUC_boxplot_celltype_timepoint_AMS.pdf"), width = 10, height = 6)
print(plot_auc_box(auc, "celltype x timepoint module - normalized AUC", auc_stats))
dev.off()

