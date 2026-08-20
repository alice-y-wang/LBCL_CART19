################################################################################
## Module scoring of the CAR IFN and TNFA signatures in the public CAR-T
## reference dataset (GSE197268), T cell compartment
################################################################################

# ---- Libraries --------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ggpubr)
library(ggrepel)

# ---- Config -----------------------------------------------------------------
analysis_dir <- "."
setwd(analysis_dir)

ref_dir   <- "GSE197268_data"
table_dir <- file.path("results", "module_scores", "tables")
plot_dir  <- file.path("results", "module_scores", "plots")

for (d in c(table_dir, plot_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

f_ref <- "GSE197268_Harmony_Integrated.RDS"

# response labels as encoded in the public dataset
response_colors <- c(R = "#4FB7C5", NR = "#E7A75E")

set.seed(1)

# ---- Gene modules -----------------------------------------------------------
# derived from leading-edge genes in fgsea 
modules <- list(
  car_ifn = c("IFIT3", "MX1", "EPSTI1", "RSAD2", "CMPK2", "IFITM3", "USP18",
              "PARP9", "IRF7", "OAS1", "TRAFD1", "TRIM21", "IFI35", "OASL",
              "BST2", "DDX60", "PSMB9", "SAMD9L", "SAMD9", "RTP4", "IFI30",
              "ISG15", "IFITM1", "TAP1", "UBE2L6", "PSMB8", "IFIT1", "MX2",
              "STAT1", "CD38", "OAS3", "OAS2", "TNFSF10", "XAF1", "HLA-DRB1",
              "GBP1", "FANCF", "IFI6", "FKBP5", "CDK1", "IRF2", "GBP5",
              "PTPN6", "HLA-DRA"),
  
  car_tnfa = c("PER1", "NR4A2", "PFKFB3", "ZBTB10", "CDKN1A", "DUSP1", "FOSL2",
               "PDE4B", "DUSP4", "IL23A", "CD83", "KDM6B", "CD69", "JUNB",
               "NR4A3", "JUN", "BTG2", "DUSP2", "MAP3K8", "ZFP36", "REL",
               "TNFAIP3", "ZC3H12A", "NINJ1", "BCL3", "FOSB", "IER5", "KLF10",
               "BTG1", "NFKBIE", "BTG3", "RELB", "PTGER4", "NFKBIA", "NFKB1",
               "PNRC1", "GADD45B", "EIF1", "PPP1R15A", "CEBPB", "BHLHE40",
               "CCNL1", "B4GALT1", "IER2", "DUSP5", "KLF6")
)

plot_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    strip.text       = element_text(face = "bold"),
    axis.title       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

# 1. Load and subset # ------

obj <- readRDS(file.path(ref_dir, f_ref))
obj <- subset(obj, subset = cell_type %in% c("CD4 T", "CD8 T", "Infusion T"))
obj <- subset(obj, subset = response %in% names(response_colors))

split_ids     <- str_split_fixed(obj$patient_id, "-", 2)
obj$patient   <- split_ids[, 1]
obj$timepoint <- split_ids[, 2]


# 2. Score modules with AddModuleScore # ------

DefaultAssay(obj) <- "RNA"

for (sig in names(modules)) {
  message("scoring ", sig)
  obj <- AddModuleScore(obj, features = modules[sig], name = sig, seed = 1)
  # Seurat appends the cluster index to the name, so rename back
  colnames(obj@meta.data)[colnames(obj@meta.data) == paste0(sig, "1")] <- sig
}

meta <- obj@meta.data
meta$cell_bc <- rownames(meta)
rm(obj); gc(verbose = FALSE)

write_csv(meta, file.path(table_dir, "tcell_module_scores.csv"))


# 3. Patient-level Wilcoxon (R vs NR), per module and group # ------

group_vars <- c("CAR", "cell_type", "timepoint", "generic")

patient_means <- meta %>%
  dplyr::filter(generic == "Tisa-cel") %>% 
  group_by(across(all_of(c("patient", group_vars, "response")))) %>%
  summarise(across(all_of(names(modules)), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(all_of(names(modules)), names_to = "signature", values_to = "agg_expr")

wilcox_res <- patient_means %>%
  group_by(across(all_of(c("signature", group_vars)))) %>%
  group_modify(~ {
    dat  <- .x
    n_r  <- sum(dat$response == "R"  & !is.na(dat$agg_expr))
    n_nr <- sum(dat$response == "NR" & !is.na(dat$agg_expr))
    if (n_r < 2 || n_nr < 2) return(tibble())
    wt <- suppressWarnings(
      wilcox.test(dat$agg_expr[dat$response == "R"],
                  dat$agg_expr[dat$response == "NR"])
    )
    tibble(
      statistic = unname(wt$statistic),
      p.value   = wt$p.value,
      R_median  = median(dat$agg_expr[dat$response == "R"],  na.rm = TRUE),
      NR_median = median(dat$agg_expr[dat$response == "NR"], na.rm = TRUE),
      R_n       = n_r,
      NR_n      = n_nr
    )
  }) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value)

write_csv(wilcox_res, file.path(table_dir, "tcell_patient_level_wilcox.csv"))


# 4. Plots: one PDF per module (patient-level, boxplot, violin) # ------

for (sig in names(modules)) {
  
  df <- meta %>% filter(!is.na(.data[[sig]]))
  df$response <- factor(df$response, levels = names(response_colors))
  df <- df %>% filter(!is.na(response))
  df <- df %>% filter(generic == "Tisa-cel")
  if (!nrow(df)) next
  
  pd <- position_dodge(width = 0.8)
  
  p_box <- ggplot(df, aes(x = CAR, y = .data[[sig]], fill = response)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.8, width = 0.6, position = pd) +
    facet_wrap(vars(timepoint, cell_type, generic), scales = "free_x") +
    stat_compare_means(aes(group = response), method = "wilcox.test",
                       label = "p.format", hide.ns = FALSE) +
    scale_fill_manual(values = response_colors) +
    labs(title = sig, x = "CAR", y = "Module score", fill = "response") +
    plot_theme
  
  p_vln <- ggplot(df, aes(x = CAR, y = .data[[sig]], fill = response)) +
    geom_violin(trim = FALSE, alpha = 0.4, width = 0.8, scale = "width", position = pd) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8, position = pd) +
    facet_wrap(vars(timepoint, cell_type, generic), scales = "free_x") +
    stat_compare_means(aes(group = response), method = "wilcox.test",
                       label = "p.format", hide.ns = FALSE) +
    scale_fill_manual(values = response_colors) +
    labs(title = sig, x = "CAR", y = "Module score", fill = "response") +
    plot_theme
  
  agg <- df %>%
    group_by(across(all_of(c("patient", "response", "CAR", group_vars)))) %>%
    summarise(mean_module_score = mean(.data[[sig]], na.rm = TRUE), .groups = "drop")
  
  pd_pt <- position_jitterdodge(jitter.width = 0, dodge.width = 0.6, seed = 42)
  
  p_pt <- ggplot(agg, aes(x = CAR, y = mean_module_score, fill = response)) +
    geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA,
                 position = position_dodge(width = 0.6)) +
    geom_point(aes(group = response), position = pd_pt,
               size = 2, alpha = 0.8, shape = 21, stroke = 0.3) +
    geom_text_repel(aes(label = patient, group = response),
                    position = pd_pt, size = 2.5, max.overlaps = 50,
                    box.padding = 0.3, point.padding = 0.2, segment.color = "grey50") +
    facet_wrap(vars(timepoint, cell_type, generic), scales = "free_x") +
    stat_compare_means(aes(group = response), method = "wilcox.test",
                       label = "p.format", hide.ns = FALSE) +
    scale_fill_manual(values = response_colors) +
    labs(title = paste0("Patient-level ", sig), x = "CAR",
         y = "Mean module score", fill = "response") +
    plot_theme
  
  pdf(file.path(plot_dir, paste0("boxplot_tcell_", sig, ".pdf")), width = 25, height = 18)
  print(p_pt); print(p_box); print(p_vln)
  dev.off()
}