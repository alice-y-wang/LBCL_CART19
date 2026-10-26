################################################################################
## LIANA+ rank-aggregate results: Fisher's-method combination and visualization
################################################################################

# ---- Libraries --------------------------------------------------------------
library(tidyverse)
library(readr)
library(ggplot2)
library(ggrepel)
library(poolr)

# ---- CONFIG -----------------------------------------------------------------
analysis_dir <- "."
setwd(analysis_dir)

dea_dir    <- file.path("python_data", "DEA_save")     # LIANA DEA + Fisher tables
fig_dir    <- file.path("images", "liana_fishers")     # figures
files_dir  <- file.path(fig_dir, "files")              # exported top-interaction tables

for (d in c(dea_dir, fig_dir, files_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Condition names used for the Fisher output files (and as plot labels below)
conditions <- c("CAR_TDN-D0", "nonCAR_APH", "CAR_Peak", "CAR_Week4")

# Cell-type groupings (values of cell.anno)
t.cells       <- c("Proliferating", "CD8 TEM", "CD4 TEM-like")
myeloid.cells <- c("CD14 Mono", "CD16 Mono", "Dendritic Cell")

# Diverging palette shared by every interaction-statistic scale
col_low  <- "#E7A75E"
col_mid  <- "white"
col_high <- "#4FB7C5"

set.seed(2024)


# 1. Load LIANA DEA rank-aggregate tables (one per condition) # ----------------

dea.car.tdnD0    <- read_csv(file.path(dea_dir, "LIANA_DEA_CAR_TDN_D0_rank_agg_ALL.csv"))
dea.noncar.aph   <- read_csv(file.path(dea_dir, "LIANA_DEA_NonCAR_APH_rank_agg_ALL.csv"))
dea.car.peak     <- read_csv(file.path(dea_dir, "LIANA_DEA_CAR_Peak_rank_agg_ALL.csv"))
dea.car.week4    <- read_csv(file.path(dea_dir, "LIANA_DEA_CAR_Week4_rank_agg_ALL.csv"))

# same order as `conditions`
dataframes <- list(dea.car.tdnD0, dea.noncar.aph, dea.car.peak, dea.car.week4)


# 2. Fisher's method helper # ----------------

apply_fisher <- function(pvalues) {
  # Remove NA values
  pvalues <- pvalues[!is.na(pvalues)]
  
  # Check if we have valid p-values
  if (length(pvalues) == 0) {
    return(list(p = NA, statistic = NA, n = 0))
  }
  
  # poolr::fisher handles edge cases well
  result <- fisher(pvalues)
  
  return(list(
    p = result$p,
    statistic = result$statistic,
    n = length(pvalues)
  ))
}


# TRUE if gene g is one of the subunits of a complex string
has <- function(x, g) str_detect(x, paste0("(^|_)", g, "(_|$)"))


# --- target cell types allowed to use each co-receptor ---------------------
cd4_ok <- c("CD14 Mono", "CD16 Mono", "Dendritic Cell", "CD4 Naive", 
            "CD4 TEM-like", "Treg", "Proliferating")

cd8_ok <- c("CD8 TEM", "CD8 Naive", "gdT", "NK", "Proliferating")


# filter for cell surface receptors only (ligand based analysis) # -------------
surface_db <- read_csv("python_data/consensus_cell_surface_receptors_2.csv")
surface_receptors <- unique(surface_db$receptor)

surface_lig <- read_csv("python_data/consensus_extracellular_ligands.csv")
surface_ligands <- unique(surface_lig$ligand)


# 3. Combine p-values across TARGETS, per (source, ligand, receptor) # ----------------
#    -- one pass on ligand p-values, one on interaction p-values

for (i in 1:length(conditions)){
  df <- dataframes[[i]]
  fishers_results <- df %>%
    # surface/secreted ligands & receptors (all complex subunits must pass)
    filter(ligand_complex %in% surface_ligands) %>% 
    filter(receptor_complex %in% surface_receptors) %>% 
    # lineage-consistent co-receptor usage
    filter(!(has(receptor_complex, "CD8[AB]") & !target %in% cd8_ok),
           !(has(receptor_complex, "CD4")     & !target %in% cd4_ok)) %>%
    group_by(source, ligand, receptor, ligand_complex, receptor_complex) %>%
    dplyr::summarise(
      # Apply Fisher's method to ligand p-values
      fisher_result = list(apply_fisher(ligand_pvalue)),
      
      # Keep some additional useful information
      n_targets = n(),
      targets = paste(unique(target), collapse = "; "),
      mean_interaction_stat = mean(interaction_stat, na.rm = TRUE),
      min_ligand_pvalue = min(ligand_pvalue, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    # Unpack the Fisher's method results
    mutate(
      combined_pvalue = sapply(fisher_result, function(x) x$p),
      chi_square = sapply(fisher_result, function(x) x$statistic),
      n_pvalues_combined = sapply(fisher_result, function(x) x$n)
    ) %>%
    select(-fisher_result) %>%
    # Add FDR correction on combined p-values
    mutate(combined_padj = p.adjust(combined_pvalue, method = "BH")) %>%
    # Sort by combined p-value
    arrange(combined_pvalue)
  
  # Display results
  print(fishers_results)
  
  # save fishers results
  write.csv(fishers_results,
            file = file.path(dea_dir, paste0(conditions[i], "_fisher_ligand_results_strict.csv")))
  
  # Summary statistics
  cat("\n=== Summary ===\n")
  cat("Total ligand-receptor pairs analyzed:", nrow(fishers_results), "\n")
  cat("Significant at p < 0.05:", sum(fishers_results$combined_pvalue < 0.05, na.rm = TRUE), "\n")
  cat("Significant at FDR < 0.05:", sum(fishers_results$combined_padj < 0.05, na.rm = TRUE), "\n")
  
}



# 4. Combine p-values across SOURCES, per (target, ligand, receptor) # ----------------
#    -- receptor p-values

for (i in 1:length(conditions)){
  df <- dataframes[[i]]
  fishers_results <- df %>%
    # surface/secreted ligands & receptors (all complex subunits must pass)
    filter(ligand_complex %in% surface_ligands) %>% 
    filter(receptor_complex %in% surface_receptors) %>% 
    # lineage-consistent co-receptor usage
    filter(!(has(receptor_complex, "CD8[AB]") & !target %in% cd8_ok),
           !(has(receptor_complex, "CD4")     & !target %in% cd4_ok)) %>%
    group_by(target, ligand, receptor, ligand_complex, receptor_complex) %>%
    summarise(
      # Apply Fisher's method to receptor p-values
      fisher_result = list(apply_fisher(receptor_pvalue)),
      
      # Keep some additional useful information
      n_sources = n(),
      sources = paste(unique(source), collapse = "; "),
      mean_interaction_stat = mean(interaction_stat, na.rm = TRUE),
      min_receptor_pvalue = min(receptor_pvalue, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    # Unpack the Fisher's method results
    mutate(
      combined_pvalue = sapply(fisher_result, function(x) x$p),
      chi_square = sapply(fisher_result, function(x) x$statistic),
      n_pvalues_combined = sapply(fisher_result, function(x) x$n)
    ) %>%
    select(-fisher_result) %>%
    # Add FDR correction on combined p-values
    mutate(combined_padj = p.adjust(combined_pvalue, method = "BH")) %>%
    # Sort by combined p-value
    arrange(combined_pvalue)
  
  # Display results
  print(fishers_results)
  
  # save fishers results
  write.csv(fishers_results,
            file = file.path(dea_dir, paste0(conditions[i], "_fisher_receptor_results_strict.csv")))
  
  # Summary statistics
  cat("\n=== Summary ===\n")
  cat("Total ligand-receptor pairs analyzed:", nrow(fishers_results), "\n")
  cat("Significant at p < 0.05:", sum(fishers_results$combined_pvalue < 0.05, na.rm = TRUE), "\n")
  cat("Significant at FDR < 0.05:", sum(fishers_results$combined_padj < 0.05, na.rm = TRUE), "\n")
  
}

# 5. Exploratory: top ligand-based interactions per condition x cell group # ----------------

# visualize fisher results strict # --------------
list_celltypes <- list(myeloid.cells) # or list(t.cells)
for (i in 1:length(conditions)){
  
  fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[i], "_fisher_ligand_results_strict.csv")))
  
  dea.raa <- dataframes[[i]]
  
  for (c in list_celltypes){
    print(c)
    celltype <- c[1]
    if ("CD14 Mono" %in% c){
      celltype <- "Myeloid"
    } else if ("Dendritic Cell" %in% c){
      celltype <- "DC"
    } else if ("CD4 TEM-like" %in% c){
      celltype <- "All_Tcell"
    }
    top_interactions <- fishers_results %>%
      filter(source %in% c) %>%
      filter(receptor_complex %in% surface_receptors) %>%   # keep cell-surface receptors only
      filter(!is.na(combined_padj)) %>%
      arrange(combined_padj) %>%
      head(20) %>%
      mutate(
        # Create interaction label
        interaction_label = paste0(ligand_complex, "-", receptor_complex),
        # Transform p-value for visualization
        neg_log10_pvalue = -log10(combined_padj)
      )
    write_csv(top_interactions, file.path(files_dir, paste0(conditions[i], "_", celltype, "_top_interactions_surface_receptor.csv")))
    if (nrow(top_interactions) == 0) {
      message("No cell-surface receptor interactions for ", celltype, " in ", conditions[i])
      next
    }
  }
}

# 6. Curated Dotplots and Heatmaps (ligand based analysis) # --------------------------
library(tidyverse)
library(scales)
library(patchwork)

# --- Paths --------------------------------------------------------------------
base_dir      <- fig_dir
top_files_dir <- files_dir                        # top-interaction tables from section 5
new_fig_dir   <- file.path(base_dir, "figures")
dir.create(new_fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- LIANA DEA tables, keyed by the condition name used in the file names ----
dea_list <- list(
  "CAR_TDN-D0"  = dea.car.tdnD0,
  "nonCAR_APH" = dea.noncar.aph,
  "CAR_Peak"    = dea.car.peak,
  "CAR_Week4"   = dea.car.week4
)

# --- Original Fisher results (one per condition, in dea_dir) ----------------
fisher_suffix <- "_fisher_ligand_results_strict.csv"
padj_cutoff   <- 0.01

# --- Source cell groups, named exactly as in the file names ------------------
# Cell-type groupings (values of cell.anno)
t.cells       <- c("Proliferating", "CD8 TEM", "CD4 TEM-like")
myeloid.cells <- c("CD14 Mono", "CD16 Mono", "Dendritic Cell")

source_groups <- list(
  "All_Tcell" = t.cells,
  "Myeloid"   = myeloid.cells
)

# --- Cell-type border colours (same palette as before) -----------------------
mono.l1.ctypes <- c("CD16 Mono", "CD14 Mono", "Dendritic Cell")
mono.l1.cols   <- c("sienna", "darkorange1", "firebrick")
names(mono.l1.cols) <- mono.l1.ctypes
names2 <- c("CD4 TEM-like", "CD4 Naive", "CD4 Proliferating",
            "CD8 TEM", "CD8 Naive", "CD8 Proliferating",
            "gdT", "MAIT", "NK", "Treg")
final2.colors <- c("steelblue3", "lightskyblue", "navy",
                   "darkgreen", "darkolivegreen2", "seagreen4",
                   "pink1", "red", "gold", "plum2")
ctypenames.plot   <- c(mono.l1.ctypes, names2, c("removed", "B Cell", "Proliferating"))
ctype.colors.plot <- c(mono.l1.cols, final2.colors, c("brown4", "magenta", "gray"))
names(ctype.colors.plot) <- ctypenames.plot


# --- Loop over every condition x cell-type file -------------------------------
suffix <- "_top_interactions_surface_receptor.csv"
top_files <- list.files(top_files_dir, pattern = paste0(suffix, "$"))

# Plots kept as "<condition>|<celltype>" for the combined panels at the end
dot_plots  <- list()
heat_plots <- list()
heat_rows  <- list()

for (f in top_files) {
  
  # "<condition>_<celltype>_top_interactions_surface_receptor.csv"
  stem      <- str_remove(f, fixed(suffix))
  condition <- names(dea_list)[str_starts(stem, fixed(paste0(names(dea_list), "_")))]
  if (length(condition) != 1) {
    warning("Could not match a condition for ", f, " - skipping")
    next
  }
  celltype <- str_remove(stem, fixed(paste0(condition, "_")))
  message("Plotting ", condition, " | ", celltype)
  
  # 1. Ligand and receptor lists from the top-interactions file
  top_file      <- read_csv(file.path(top_files_dir, f), show_col_types = FALSE)
  ligand_list   <- unique(top_file$ligand_complex)
  receptor_list <- unique(top_file$receptor_complex)
  
  # 2. Source cells for this cell-type group
  if (!celltype %in% names(source_groups)) {
    warning("No source group named '", celltype, "' - skipping ", f)
    next
  }
  sources <- source_groups[[celltype]]
  
  # 3. Re-filter the ORIGINAL Fisher results with those lists
  fishers_results <- read_csv(file.path(dea_dir, paste0(condition, fisher_suffix)),
                              show_col_types = FALSE)
  
  top_interactions <- fishers_results %>%
    filter(source %in% sources) %>%
    filter(!is.na(combined_padj)) %>%
    filter(combined_padj < padj_cutoff) %>%
    filter(ligand_complex %in% ligand_list) %>%
    filter(receptor_complex %in% receptor_list) %>%
    arrange(combined_padj) %>%
    mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
           neg_log10_pvalue  = -log10(combined_padj))
  
  message("  ", length(ligand_list), " ligands x ", length(receptor_list),
          " receptors -> ", nrow(top_interactions), " significant interactions")
  
  if (nrow(top_interactions) == 0) {
    message("  no interactions - skipping")
    next
  }
  
  write_csv(top_interactions,
            file.path(new_fig_dir, paste0(condition, "_", str_replace_all(celltype, " ", "_"),
                                          "_refiltered_interactions.csv")))
  
  # ---- Dot plot: fill = interaction stat, border = source cell type ----------
  p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                     y = reorder(interaction_label, neg_log10_pvalue))) +
    geom_point(aes(size = scales::squish(neg_log10_pvalue, c(2, 30)),
                   fill = mean_interaction_stat,
                   colour = source),
               shape = 21, stroke = 1.5) +
    scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0,
                         limits = c(-3, 3), oob = squish,
                         name = "interaction_statistic") +
    scale_colour_manual(values = ctype.colors.plot, na.value = "black",
                        name = "Source cell type") +
    scale_size_continuous(range = c(1, 15), limits = c(2, 30), breaks = c(2, 5, 10, 20, 30),
                          name = "-log10(p-value)") +
    scale_x_continuous(limits = c(-3.3, 3.3), breaks = -3:3,
                       expand = c(0, 0), oob = squish) +
    guides(colour = guide_legend(override.aes = list(size = 5, fill = "white", stroke = 2))) +
    labs(title = paste0(condition, " | ", celltype, ": top cell-surface receptor interactions"),
         subtitle = "Ordered by combined p-value (most significant at top)",
         x = "Mean Interaction Statistic",
         y = "ligand_receptor Pair") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          legend.position = "right")
  
  # ---- Heatmap across source -> target ----------------------------------------
  top_pairs <- top_interactions %>%
    select(source, ligand_complex, receptor_complex) %>%
    distinct()
  
  interaction_order <- top_interactions %>%
    arrange(desc(combined_padj)) %>%
    pull(interaction_label) %>%
    unique()
  
  heatmap_data <- dea_list[[condition]] %>%
    inner_join(top_pairs, by = c("source", "ligand_complex", "receptor_complex")) %>%
    mutate(interaction_label = factor(paste0(ligand_complex, "-", receptor_complex),
                                      levels = interaction_order),
           source_target = paste0(source, " → ", target))
  
  p2 <- ggplot(heatmap_data, aes(x = source_target,
                                 y = interaction_label,
                                 fill = interaction_stat)) +
    geom_tile(color = "white") +
    # one panel per SENDER cell, so there is a break between sender groups
    facet_grid(cols = vars(source), scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                         midpoint = 0, name = "Interaction\nStatistic") +
    labs(title = paste0(condition, " | ", celltype, ": interaction statistics across source-target pairs"),
         subtitle = "Ordered by combined p-value",
         x = "Source → Target",
         y = "ligand receptor Pair") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 9),
          panel.grid = element_blank(),
          panel.spacing.x = unit(0.8, "lines"),
          strip.text.x = element_text(face = "bold", size = 9))
  
  key <- paste(condition, celltype, sep = "|")
  dot_plots[[key]]  <- p1
  heat_plots[[key]] <- p2
  heat_rows[[key]]  <- n_distinct(heatmap_data$interaction_label)
  
  # ---- Save -------------------------------------------------------------------
  n_cols   <- n_distinct(heatmap_data$source_target)
  n_rows   <- n_distinct(top_interactions$interaction_label)
  pdf_w    <- max(10, 3 + 0.25 * n_cols)
  pdf_h    <- max(6, 2 + 0.3 * n_rows)
  out_name <- paste0(condition, "_", str_replace_all(celltype, " ", "_"),
                     "_top_surface_receptor_interactions.pdf")
  
  pdf(file.path(new_fig_dir, out_name), width = pdf_w, height = pdf_h)
  print(p1)
  print(p2)
  dev.off()
}

# ---- Combined panels: APH, TDN-D0, Peak ---------------------------------------
panel_conditions <- c("nonCAR_APH", "CAR_TDN-D0", "CAR_Peak")
panel_celltypes  <- c("All_Tcell", "Myeloid")

# Dot plots: the three timepoints side by side, one row per cell type
row_of_dotplots <- function(celltype) {
  plots <- map(panel_conditions, function(cond) {
    p <- dot_plots[[paste(cond, celltype, sep = "|")]]
    if (is.null(p)) return(plot_spacer())
    p + labs(title = cond, subtitle = NULL)
  })
  wrap_plots(plots, nrow = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(title = paste0(celltype, " (sender): top cell-surface receptor interactions"))
}

pdf(file.path(new_fig_dir, "Combined_APH_TDN-D0_Peak_top_surface_receptor_dotplots.pdf"),
    width = 30, height = 8)
for (ct in panel_celltypes) print(row_of_dotplots(ct))
dev.off()

# Heatmaps: the three timepoints stacked, with breaks between sender cells
stack_of_heatmaps <- function(celltype) {
  keys <- paste(panel_conditions, celltype, sep = "|")
  keep <- keys[keys %in% names(heat_plots)]
  if (length(keep) == 0) return(NULL)
  plots <- map(keep, function(k) {
    heat_plots[[k]] + labs(title = str_remove(k, "\\|.*$"), subtitle = NULL)
  })
  rows <- unlist(heat_rows[keep])
  wrap_plots(plots, ncol = 1) +
    plot_layout(heights = rows) +
    plot_annotation(title = paste0(celltype, " (sender): interaction statistics across source-target pairs"))
}

for (ct in panel_celltypes) {
  p <- stack_of_heatmaps(ct)
  if (is.null(p)) next
  keys  <- paste(panel_conditions, ct, sep = "|")
  keys  <- keys[keys %in% names(heat_plots)]
  n_col <- max(map_int(keys, ~ n_distinct(heat_plots[[.x]]$data$source_target)))
  pdf(file.path(new_fig_dir, paste0("Combined_APH_TDN-D0_Peak_", ct, "_top_surface_receptor_heatmaps.pdf")),
      width = max(12, 4 + 0.3 * n_col),
      height = max(10, 3 + 0.3 * sum(unlist(heat_rows[keys])) + 2.5 * length(keys)))
  print(p)
  dev.off()
}



# 7. Concatenate Fisher ligand results across timepoints # ----------------

fishers_results_tdnd0 <- read_csv(file.path(dea_dir, paste0(conditions[1], "_fisher_ligand_results_strict.csv")))
fishers_results_aph   <- read_csv(file.path(dea_dir, paste0(conditions[2], "_fisher_ligand_results_strict.csv")))
fishers_results_peak  <- read_csv(file.path(dea_dir, paste0(conditions[3], "_fisher_ligand_results_strict.csv")))
fishers_results_week4 <- read_csv(file.path(dea_dir, paste0(conditions[4], "_fisher_ligand_results_strict.csv")))

fishers_results_tdnd0$timepoint <- "TDN-D0"
fishers_results_aph$timepoint   <- "APH"
fishers_results_peak$timepoint  <- "Peak"
fishers_results_week4$timepoint <- "Week4"

fishers_all <- do.call(rbind, list(
  fishers_results_aph,
  fishers_results_peak,
  fishers_results_tdnd0,
  fishers_results_week4
))
fishers_all$timepoint <- factor(fishers_all$timepoint, levels = c("APH", "TDN-D0", "Peak", "Week4"))


# 8. Line graph: TNF / LTA / LTB signaling over time # ---------------

fishers_all_plot <- fishers_all %>% filter(source %in% t.cells) %>%
  filter(combined_padj < 0.01) %>%
  filter(ligand_complex %in% c("TNF", "LTA", "LTB")) %>%
  filter(receptor_complex %in% c("TNFRSF1B",
                                 "TNFRSF1A", "TNFRSF1A_TNFRSF1B"))
pdf(file.path(fig_dir, "TNF_LTA_LTB_Interactions.pdf"), width = 10, height = 6)

print(
  ggplot(fishers_all_plot, aes(
    x = timepoint,
    y = mean_interaction_stat,
    group = interaction(ligand_complex, receptor_complex),
    color = interaction(ligand_complex, receptor_complex)
  )) +
    geom_line() +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_text_repel(
      data = fishers_all_plot,
      aes(label = paste0(ligand_complex, "^", receptor_complex)),
      size = 3,
      direction = "y",
      hjust = 0,
      segment.color = "grey70"
    ) +
    theme_classic() +
    labs(
      x = "Timepoint",
      y = "Interaction statistic"
    )
)

dev.off()


# 9. Exploratory: top receptor-based interactions per condition x cell group # ----------------

list_celltypes <- list(t.cells, myeloid.cells)

# RECEPTOR-based Fisher results: cell types below are the TARGET (receiving) cells
receptor_fisher_suffix <- "_fisher_receptor_results_strict.csv"

list_celltypes <- list(myeloid.cells, t.cells)

for (i in 1:length(conditions)){
  
  fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[i], receptor_fisher_suffix)))
  dea.raa <- dataframes[[i]]
  
  for (c in list_celltypes){
    print(c)
    celltype <- c[1]
    if ("CD14 Mono" %in% c){
      celltype <- "Myeloid"
    } else if ("Dendritic Cell" %in% c){
      celltype <- "DC"
    } else if ("CD4 TEM-like" %in% c){
      celltype <- "All_Tcell"
    }
    
    
    top_interactions <- fishers_results %>%
      filter(target %in% c) %>%
      filter(ligand_complex %in% surface_ligands) %>% 
      filter(receptor_complex %in% surface_receptors) %>% 
      # drop CD8A / CD8B (alone or in a complex) as receptors on CD4 T cells
      filter(!(has(receptor_complex, "CD8[AB]") & !target %in% cd8_ok),
             !(has(receptor_complex, "CD4")     & !target %in% cd4_ok)) %>%
      filter(!is.na(combined_padj)) %>%
      arrange(combined_padj) %>%
      head(20) %>%
      mutate(
        # Create interaction label
        interaction_label = paste0(ligand_complex, "-", receptor_complex),
        # Transform p-value for visualization
        neg_log10_pvalue = -log10(combined_padj)
      )
    
    if (nrow(top_interactions) == 0) {
      message("No cell-surface receptor interactions for ", celltype, " (target) in ", conditions[i])
      next
    }
    write_csv(top_interactions, file.path(files_dir, paste0(conditions[i], "_", celltype, "_top_interactions_surface_receptor_RECEPTOR.csv")))
    
    ctype.colors.plot <- c(mono.l1.cols, final2.colors, c("brown4", "magenta", "gray"))
    names(ctype.colors.plot) <- ctypenames.plot
    
  }
}

# 10. Curated Dotplots and Heatmaps, RECEPTOR-based version: the cell-type groups are the TARGET (receiving) cells, ###############-------
# and the Fisher results are the per-target RECEPTOR results.

# --- Paths --------------------------------------------------------------------
base_dir      <- fig_dir
top_files_dir <- files_dir                        # top-interaction tables from section 9
new_fig_dir   <- file.path(base_dir, "figures")
dir.create(new_fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- LIANA DEA tables, keyed by the condition name used in the file names ----
dea_list <- list(
  "CAR_TDN-D0"  = dea.car.tdnD0,
  "nonCAR_APH" = dea.noncar.aph,
  "CAR_Peak"    = dea.car.peak,
  "CAR_Week4"   = dea.car.week4
)

# --- Original Fisher RECEPTOR results (one per condition) -------------------
fisher_dir    <- dea_dir
fisher_suffix <- "_fisher_receptor_results_strict.csv"
padj_cutoff   <- 0.01

# --- Target cell groups, named exactly as in the file names ------------------
# Cell-type groupings (values of cell.anno)
t.cells       <- c("Proliferating", "CD8 TEM", "CD4 TEM-like")
myeloid.cells <- c("CD14 Mono", "CD16 Mono", "Dendritic Cell")

target_groups <- list(
  "All_Tcell" = t.cells,
  "CD8 TEM"   = "CD8 TEM",
  "Myeloid"   = myeloid.cells
)

# --- Cell-type border colours (same palette as before) -----------------------
mono.l1.ctypes <- c("CD16 Mono", "CD14 Mono", "Dendritic Cell")
mono.l1.cols   <- c("sienna", "darkorange1", "firebrick")
names(mono.l1.cols) <- mono.l1.ctypes
names2 <- c("CD4 TEM-like", "CD4 Naive", "CD4 Proliferating",
            "CD8 TEM", "CD8 Naive", "CD8 Proliferating",
            "gdT", "MAIT", "NK", "Treg")
final2.colors <- c("steelblue3", "lightskyblue", "navy",
                   "darkgreen", "darkolivegreen2", "seagreen4",
                   "pink1", "red", "gold", "plum2")
ctypenames.plot   <- c(mono.l1.ctypes, names2, c("removed", "B Cell", "Proliferating"))
ctype.colors.plot <- c(mono.l1.cols, final2.colors, c("brown4", "magenta", "gray"))
names(ctype.colors.plot) <- ctypenames.plot


# --- Loop over every condition x cell-type file -------------------------------
# Receptor-based top files written in section 9
suffix <- "_top_interactions_surface_receptor_RECEPTOR.csv"
top_files <- list.files(top_files_dir, pattern = paste0(suffix, "$"))

# Plots kept as "<condition>|<celltype>" for the combined panels at the end
dot_plots  <- list()
heat_plots <- list()
heat_rows  <- list()

for (f in top_files) {
  
  # "<condition>_<celltype>_top_interactions_surface_receptor_RECEPTOR.csv"
  stem      <- str_remove(f, fixed(suffix))
  condition <- names(dea_list)[str_starts(stem, fixed(paste0(names(dea_list), "_")))]
  if (length(condition) != 1) {
    warning("Could not match a condition for ", f, " - skipping")
    next
  }
  celltype <- str_remove(stem, fixed(paste0(condition, "_")))
  message("Plotting ", condition, " | ", celltype, " (as TARGET)")
  
  # 1. Ligand and receptor lists from the top-interactions file
  top_file      <- read_csv(file.path(top_files_dir, f), show_col_types = FALSE)
  ligand_list   <- unique(top_file$ligand_complex)
  receptor_list <- unique(top_file$receptor_complex)
  
  # 2. Target cells for this cell-type group
  if (!celltype %in% names(target_groups)) {
    warning("No target group named '", celltype, "' - skipping ", f)
    next
  }
  targets <- target_groups[[celltype]]
  
  # 3. Re-filter the ORIGINAL Fisher RECEPTOR results with those lists
  fishers_results <- read_csv(file.path(fisher_dir, paste0(condition, fisher_suffix)),
                              show_col_types = FALSE)
  
  top_interactions <- fishers_results %>%
    filter(target %in% targets) %>%
    filter(!is.na(combined_padj)) %>%
    filter(combined_padj < padj_cutoff) %>%
    filter(ligand_complex %in% ligand_list) %>%
    filter(receptor_complex %in% receptor_list) %>%
    arrange(combined_padj) %>%
    mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
           neg_log10_pvalue  = -log10(combined_padj))
  
  message("  ", length(ligand_list), " ligands x ", length(receptor_list),
          " receptors -> ", nrow(top_interactions), " significant interactions")
  
  if (nrow(top_interactions) == 0) {
    message("  no interactions - skipping")
    next
  }
  
  write_csv(top_interactions,
            file.path(new_fig_dir, paste0(condition, "_", str_replace_all(celltype, " ", "_"),
                                          "_RECEPTOR_refiltered_interactions.csv")))
  
  # ---- Dot plot: fill = interaction stat, border = target cell type ----------
  p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                     y = reorder(interaction_label, neg_log10_pvalue))) +
    geom_point(aes(size = scales::squish(neg_log10_pvalue, c(2, 30)),
                   fill = mean_interaction_stat,
                   colour = target),
               shape = 21, stroke = 1.5) +
    scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0,
                         limits = c(-3, 3), oob = squish,
                         name = "interaction_statistic") +
    scale_colour_manual(values = ctype.colors.plot, na.value = "black",
                        name = "Target cell type") +
    scale_size_continuous(range = c(1, 15), limits = c(2, 30), breaks = c(2, 5, 10, 20, 30),
                          name = "-log10(p-value)") +
    scale_x_continuous(limits = c(-3.3, 3.3), breaks = -3:3,
                       expand = c(0, 0), oob = squish) +
    guides(colour = guide_legend(override.aes = list(size = 5, fill = "white", stroke = 2))) +
    labs(title = paste0(condition, " | ", celltype, " (target): top cell-surface receptor interactions"),
         subtitle = "Receptor-based Fisher results; ordered by combined p-value (most significant at top)",
         x = "Mean Interaction Statistic",
         y = "ligand_receptor Pair") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          legend.position = "right")
  
  # ---- Heatmap across source -> target ----------------------------------------
  top_pairs <- top_interactions %>%
    select(target, ligand_complex, receptor_complex) %>%
    distinct()
  
  interaction_order <- top_interactions %>%
    arrange(desc(combined_padj)) %>%
    pull(interaction_label) %>%
    unique()
  
  heatmap_data <- dea_list[[condition]] %>%
    inner_join(top_pairs, by = c("target", "ligand_complex", "receptor_complex")) %>%
    mutate(interaction_label = factor(paste0(ligand_complex, "-", receptor_complex),
                                      levels = interaction_order),
           source_target = paste0(source, " → ", target))
  
  p2 <- ggplot(heatmap_data, aes(x = source_target,
                                 y = interaction_label,
                                 fill = interaction_stat)) +
    geom_tile(color = "white") +
    # one panel per RECEIVER cell, so there is a break between receiver groups
    facet_grid(cols = vars(target), scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                         midpoint = 0, name = "Interaction\nStatistic") +
    labs(title = paste0(condition, " | ", celltype, " (target): interaction statistics across source-target pairs"),
         subtitle = "Ordered by combined p-value",
         x = "Source → Target",
         y = "ligand receptor Pair") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 9),
          panel.grid = element_blank(),
          panel.spacing.x = unit(0.8, "lines"),
          strip.text.x = element_text(face = "bold", size = 9))
  
  key <- paste(condition, celltype, sep = "|")
  dot_plots[[key]]  <- p1
  heat_plots[[key]] <- p2
  heat_rows[[key]]  <- n_distinct(heatmap_data$interaction_label)
  
  # ---- Save -------------------------------------------------------------------
  n_cols   <- n_distinct(heatmap_data$source_target)
  n_rows   <- n_distinct(top_interactions$interaction_label)
  pdf_w    <- max(10, 3 + 0.25 * n_cols)
  pdf_h    <- max(6, 2 + 0.3 * n_rows)
  out_name <- paste0(condition, "_", str_replace_all(celltype, " ", "_"),
                     "_RECEPTOR_top_surface_receptor_interactions.pdf")
  
  pdf(file.path(new_fig_dir, out_name), width = pdf_w, height = pdf_h)
  print(p1)
  print(p2)
  dev.off()
}

# ---- Combined panels: APH, TDN-D0, Peak ---------------------------------------
panel_conditions <- c("nonCAR_APH", "CAR_TDN-D0", "CAR_Peak")
panel_celltypes  <- c("All_Tcell", "Myeloid")

# Dot plots: the three timepoints side by side, one row per cell type
row_of_dotplots <- function(celltype) {
  plots <- map(panel_conditions, function(cond) {
    p <- dot_plots[[paste(cond, celltype, sep = "|")]]
    if (is.null(p)) return(plot_spacer())
    p + labs(title = cond, subtitle = NULL)
  })
  wrap_plots(plots, nrow = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(title = paste0(celltype, " (target): top cell-surface receptor interactions"))
}

pdf(file.path(new_fig_dir, "Combined_APH_TDN-D0_Peak_RECEPTOR_top_surface_receptor_dotplots.pdf"),
    width = 30, height = 8)
for (ct in panel_celltypes) print(row_of_dotplots(ct))
dev.off()

# Heatmaps: the three timepoints stacked, with breaks between receiver cells
stack_of_heatmaps <- function(celltype) {
  keys <- paste(panel_conditions, celltype, sep = "|")
  keep <- keys[keys %in% names(heat_plots)]
  if (length(keep) == 0) return(NULL)
  plots <- map(keep, function(k) {
    heat_plots[[k]] + labs(title = str_remove(k, "\\|.*$"), subtitle = NULL)
  })
  rows <- unlist(heat_rows[keep])
  wrap_plots(plots, ncol = 1) +
    plot_layout(heights = rows) +
    plot_annotation(title = paste0(celltype, " (target): interaction statistics across source-target pairs"))
}

for (ct in panel_celltypes) {
  p <- stack_of_heatmaps(ct)
  if (is.null(p)) next
  keys  <- paste(panel_conditions, ct, sep = "|")
  keys  <- keys[keys %in% names(heat_plots)]
  n_col <- max(map_int(keys, ~ n_distinct(heat_plots[[.x]]$data$source_target)))
  pdf(file.path(new_fig_dir, paste0("Combined_APH_TDN-D0_Peak_", ct, "_RECEPTOR_top_surface_receptor_heatmaps.pdf")),
      width = max(12, 4 + 0.3 * n_col),
      height = max(10, 3 + 0.3 * sum(unlist(heat_rows[keys])) + 2.5 * length(keys)))
  print(p)
  dev.off()
}