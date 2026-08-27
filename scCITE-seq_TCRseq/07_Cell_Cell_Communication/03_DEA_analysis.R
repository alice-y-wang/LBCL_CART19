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

# Conditions, matching the names written by the Python LIANA DEA script
conditions <- c("CAR_TDN_D0", "NonCAR_APH", "CAR_Peak", "CAR_Week4")

# Cell-type groupings (values of cell.anno)
t.cells       <- c("Proliferating", "CD8 EM", "CD4 EM-like")
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


# 3. Combine p-values across TARGETS, per (source, ligand, receptor) # ----------------
#    -- one pass on ligand p-values, one on interaction p-values

for (i in 1:length(conditions)){
  df <- dataframes[[i]]
  fishers_results <- df %>%
    group_by(source, ligand, receptor, ligand_complex, receptor_complex) %>%
    summarise(
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
            file = file.path(dea_dir, paste0(conditions[i], "_fisher_ligand_results.csv")))
  
  # Summary statistics
  cat("\n=== Summary ===\n")
  cat("Total ligand-receptor pairs analyzed:", nrow(fishers_results), "\n")
  cat("Significant at p < 0.05:", sum(fishers_results$combined_pvalue < 0.05, na.rm = TRUE), "\n")
  cat("Significant at FDR < 0.05:", sum(fishers_results$combined_padj < 0.05, na.rm = TRUE), "\n")
  
  fishers_results_interaction <- df %>%
    group_by(source, ligand, receptor, ligand_complex, receptor_complex) %>%
    summarise(
      # Apply Fisher's method to interaction p-values
      fisher_result = list(apply_fisher(interaction_pvalue)),
      
      # Keep some additional useful information
      n_targets = n(),
      targets = paste(unique(target), collapse = "; "),
      mean_interaction_stat = mean(interaction_stat, na.rm = TRUE),
      min_interaction_pvalue = min(interaction_pvalue, na.rm = TRUE),
      
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
  print(fishers_results_interaction)
  
  # save fishers results
  write.csv(fishers_results_interaction,
            file = file.path(dea_dir, paste0(conditions[i], "_fisher_interaction_results.csv")))
  
  # Summary statistics
  cat("\n=== Summary ===\n")
  cat("Total ligand-receptor pairs analyzed:", nrow(fishers_results_interaction), "\n")
  cat("Significant at p < 0.05:", sum(fishers_results_interaction$combined_pvalue < 0.05, na.rm = TRUE), "\n")
  cat("Significant at FDR < 0.05:", sum(fishers_results_interaction$combined_padj < 0.05, na.rm = TRUE), "\n")
  
}


# 4. Combine p-values across SOURCES, per (target, ligand, receptor) # ----------------
#    -- receptor p-values

for (i in 1:length(conditions)){
  df <- dataframes[[i]]
  fishers_results <- df %>%
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
            file = file.path(dea_dir, paste0(conditions[i], "_fisher_receptor_results.csv")))
  
  # Summary statistics
  cat("\n=== Summary ===\n")
  cat("Total ligand-receptor pairs analyzed:", nrow(fishers_results), "\n")
  cat("Significant at p < 0.05:", sum(fishers_results$combined_pvalue < 0.05, na.rm = TRUE), "\n")
  cat("Significant at FDR < 0.05:", sum(fishers_results$combined_padj < 0.05, na.rm = TRUE), "\n")
  
}


# 5. Exploratory: top ligand-based interactions per condition x cell group # ----------------

list_celltypes <- list(myeloid.cells) # or list(t.cells)

for (i in 1:length(conditions)){
  
  fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[i], "_fisher_ligand_results.csv")))
  dea.raa <- dataframes[[i]]
  
  for (c in list_celltypes){
    print(c)
    celltype <- c[1]
    if ("CD14 Mono" %in% c){
      celltype <- "Myeloid"
    } else if ("Dendritic Cell" %in% c){
      celltype <- "DC"
    } else if ("CD4 EM-like" %in% c){
      celltype <- "All_Tcell"
    }
    top_interactions <- fishers_results %>%
      filter(source %in% c) %>%
      filter(!is.na(combined_padj)) %>%
      arrange(combined_padj) %>%
      head(20) %>%
      mutate(
        # Create interaction label
        interaction_label = paste0(ligand_complex, "-", receptor_complex),
        # Transform p-value for visualization
        neg_log10_pvalue = -log10(combined_padj)
      )
    
    p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                       y = reorder(interaction_label, neg_log10_pvalue))) +
      geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
      scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0,
                             name = "interaction_statistic") +
      scale_size_continuous(range = c(3, 10),
                            name = "-log10(p-value)") +
      labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
           subtitle = "Ordered by combined p-value (most significant at top)",
           x = "Mean Interaction Statistic",
           y = "ligand_receptor Pair") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 10),
            legend.position = "right")
    
    top_pairs <- top_interactions %>%
      select(source, ligand_complex, receptor_complex) %>%
      distinct()
    
    pdf(file.path(fig_dir, paste0(conditions[i], "_", celltype, "_fisher_ligand_images.pdf")),
        width = 10, height = 8)
    print(p1)
    dev.off()
    
  }
}


# 6. Curated ligand-based figures: T cells at APH # ----------------

fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[2], "_fisher_ligand_results.csv")))

top_interactions <- fishers_results %>%
  filter(source %in% t.cells) %>%
  filter(!is.na(combined_pvalue)) %>%
  filter(combined_padj < 0.01) %>%
  filter(ligand_complex %in% c("NUCB2", "TNF", "HLA-DRB1", "TNFSF13B")) %>%
  filter(receptor_complex %in% c("ERAP1", "VSIR", "TRAF2", "TRADD", "TNFRSF1B",
                                 "TNFRSF1A", "NOTCH1", "FAS", "ICOS", "SEMA4C", "CD4",
                                 "HLA-DPB1")) %>%
  arrange(combined_padj) %>%
  mutate(
    # Create interaction label
    interaction_label = paste0(ligand_complex, "-", receptor_complex),
    # Transform p-value for visualization
    neg_log10_pvalue = -log10(combined_padj)
  )

write.csv(top_interactions, file.path(files_dir, "APH_NonCAR_Tcell_Top_Ordered.csv"))

p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                   y = reorder(interaction_label, neg_log10_pvalue))) +
  geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
  scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3, 3),
                         name = "interaction_statistic") +
  scale_size_continuous(range = c(0, 15),
                        limits = c(2, 30),
                        name = "-log10(p-value)") +
  labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
       subtitle = "Ordered by combined p-value (most significant at top)",
       x = "Mean Interaction Statistic",
       y = "ligand_receptor Pair") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")

top_pairs <- top_interactions %>%
  select(source, ligand_complex, receptor_complex) %>%
  distinct()

# heatmap restricted to significant interactions
heatmap_data <- dea.noncar.aph %>%
  inner_join(top_pairs, by = c("source", "ligand_complex", "receptor_complex")) %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
         source_target = paste0(source, " → ", target))

interaction_order <- top_interactions %>%
  arrange(desc(combined_pvalue))  %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex)) %>%
  pull(interaction_label)

interaction_order <- unique(interaction_order)

heatmap_data <- heatmap_data %>%
  mutate(interaction_label = factor(interaction_label, levels = interaction_order))

p2 <- ggplot(heatmap_data, aes(x = source_target,
                               y = interaction_label,
                               fill = interaction_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0,
                       name = "Interaction\nStatistic") +
  labs(title = "Interaction Statistics Across Source-Target Pairs",
       subtitle = "Ordered by combined p-value",
       x = "Source → Target",
       y = "ligand receptor Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(fig_dir, "NonCAR_APH_Top_TCell_Interactions.pdf"), width = 10, height = 6)
print(p1)
print(p2)
dev.off()




# 7. Curated ligand-based figures: T cells, CAR+ TDN-D0 # ----------------

fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[1], "_fisher_ligand_results.csv")))

top_interactions <- fishers_results %>%
  filter(source %in% t.cells) %>%
  filter(!is.na(combined_pvalue)) %>%
  filter(combined_padj < 0.01) %>%
  filter(ligand_complex %in% c(
    "AGRN",
    "B2M",
    "HLA",
    "ICAM2",
    "LTB",
    "NMU",
    "SOCS2",
    "TGFB1",
    "THBS4",
    "TNFSF8"
  )
  ) %>%
  filter(receptor_complex %in% c(
    "CD3D",
    "CD47",
    "CXCR4",
    "DAG1",
    "ENG",
    "EPOR",
    "ITGAL_ITGB2",
    "ITGB1",
    "KLRD1",
    "LPP",
    "NMUR1",
    "TNFRSF1A",
    "TNFRSF8",
    "TGFBR1_TGFBR2",
    "TGFBR3"
  )) %>%
  arrange(combined_padj) %>%
  mutate(
    interaction_label = paste0(ligand_complex, "-", receptor_complex),
    neg_log10_pvalue = -log10(combined_padj)
  )

write.csv(top_interactions, file.path(files_dir, "TDND0_CAR_Tcell_Top_Ordered.csv"))

p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                   y = reorder(interaction_label, neg_log10_pvalue))) +
  geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
  scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3, 3),
                         name = "interaction_statistic") +
  scale_size_continuous(range = c(0, 10),
                        limits = c(2, 30),
                        name = "-log10(p-adj)") +
  labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
       subtitle = "Ordered by combined p-value (most significant at top)",
       x = "Mean Interaction Statistic",
       y = "ligand_receptor Pair") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")

top_pairs <- top_interactions %>%
  select(source, ligand_complex, receptor_complex) %>%
  distinct()

heatmap_data <- dea.car.tdnD0 %>%
  inner_join(top_pairs, by = c("source", "ligand_complex", "receptor_complex")) %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
         source_target = paste0(source, " → ", target)) 

interaction_order <- top_interactions %>%
  arrange(desc(combined_pvalue))  %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex)) %>%
  pull(interaction_label)

interaction_order <- unique(interaction_order)

heatmap_data <- heatmap_data %>%
  mutate(interaction_label = factor(interaction_label, levels = interaction_order))

p2 <- ggplot(heatmap_data, aes(x = source_target,
                               y = interaction_label,
                               fill = interaction_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0,
                       name = "Interaction\nStatistic") +
  labs(title = "Interaction Statistics Across Source-Target Pairs",
       subtitle = "Ordered by combined p-value",
       x = "Source → Target",
       y = "ligand receptor Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(fig_dir, "CAR_TDN_D0_Top_TCell_Interactions.pdf"), width = 10, height = 6)
print(p1)
print(p2)
dev.off()




# 8. Curated ligand-based figures: T cells, CAR+ Peak # ----------------

fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[3], "_fisher_ligand_results.csv")))

top_interactions <- fishers_results %>%
  filter(source %in% t.cells) %>%
  filter(!is.na(combined_pvalue)) %>%
  filter(combined_padj < 0.01) %>%
  filter(ligand_complex %in% c(
    "CADM1",
    "COL6A2",
    "CSF2",
    "HLA-DQA1",
    "HLA-DRB1",
    "LTA",
    "TNFSF10"
  )
  ) %>%
  filter(receptor_complex %in% c(
    "CD4",
    "CD44",
    "CRTAM",
    "IL3RA",
    "ITGA1_ITGB1",
    "ITGA3_ITGB1",
    "ITGA9_ITGB1",
    "ITGA11_ITGB1",
    "ITGAV_ITGB8",
    "LAG3",
    "RIPK1",
    "TNFRSF10A",
    "TNFRSF10B",
    "TNFRSF14",
    "TNFRSF1A",
    "TNFRSF1B"
  )) %>%
  arrange(combined_padj) %>%
  mutate(
    interaction_label = paste0(ligand_complex, "-", receptor_complex),
    neg_log10_pvalue = -log10(combined_padj)
  )

write.csv(top_interactions, file.path(files_dir, "Peak_CAR_Tcell_Top_Ordered.csv"))

p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                   y = reorder(interaction_label, neg_log10_pvalue))) +
  geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
  scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3, 3),
                         name = "interaction_statistic") +
  scale_size_continuous(range = c(0, 10),
                        limits = c(2, 30),
                        name = "-log10(p-adj)") +
  labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
       subtitle = "Ordered by combined p-value (most significant at top)",
       x = "Mean Interaction Statistic",
       y = "ligand_receptor Pair") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")

top_pairs <- top_interactions %>%
  select(source, ligand_complex, receptor_complex) %>%
  distinct()

heatmap_data <- dea.car.peak %>%
  inner_join(top_pairs, by = c("source", "ligand_complex", "receptor_complex")) %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
         source_target = paste0(source, " → ", target))

interaction_order <- top_interactions %>%
  arrange(desc(combined_pvalue))  %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex)) %>%
  pull(interaction_label)

interaction_order <- unique(interaction_order)

heatmap_data <- heatmap_data %>%
  mutate(interaction_label = factor(interaction_label, levels = interaction_order))

p2 <- ggplot(heatmap_data, aes(x = source_target,
                               y = interaction_label,
                               fill = interaction_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0,
                       name = "Interaction\nStatistic") +
  labs(title = "Interaction Statistics Across Source-Target Pairs",
       subtitle = "Ordered by combined p-value",
       x = "Source → Target",
       y = "ligand receptor Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(fig_dir, "CAR_Peak_Top_TCell_Interactions.pdf"), width = 10, height = 6)
print(p1)
print(p2)
dev.off()



# 9. Curated ligand-based figures: myeloid, APH # ----------------

fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[2], "_fisher_ligand_results.csv")))

top_interactions <- fishers_results %>%
  filter(source %in% myeloid.cells) %>%
  filter(!is.na(combined_pvalue)) %>%
  filter(combined_padj < 0.01) %>%
  filter(ligand_complex %in% c(
    "ADA",
    "FCER2",
    "HLA-DPB1",
    "HLA-DRB5",
    "HSPA8",
    "S100A8",
    "S100A9",
    "VCAN"
  )
  ) %>%
  filter(receptor_complex %in% c(
    "AGER",
    "CD36",
    "CD4",
    "CD44",
    "CD68",
    "CD69",
    "DPP4",
    "ITGA4",
    "ITGAM",
    "ITGB1",
    "ITGB2",
    "LAG3",
    "LDLR",
    "SELL",
    "TLR1",
    "TLR4"
  )) %>%
  arrange(combined_padj) %>%
  mutate(
    interaction_label = paste0(ligand_complex, "-", receptor_complex),
    neg_log10_pvalue = -log10(combined_padj)
  )

write.csv(top_interactions, file.path(files_dir, "APH_NonCAR_Myeloid_Top_Ordered.csv"))

p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                   y = reorder(interaction_label, neg_log10_pvalue))) +
  geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
  scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3, 3),
                         name = "interaction_statistic") +
  scale_size_continuous(range = c(0, 10),
                        limits = c(2, 30),
                        name = "-log10(p-value)") +
  labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
       subtitle = "Ordered by combined p-value (most significant at top)",
       x = "Mean Interaction Statistic",
       y = "ligand_receptor Pair") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")

top_pairs <- top_interactions %>%
  select(source, ligand_complex, receptor_complex) %>%
  distinct()

heatmap_data <- dea.noncar.aph %>%
  inner_join(top_pairs, by = c("source", "ligand_complex", "receptor_complex")) %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
         source_target = paste0(source, " → ", target))

interaction_order <- top_interactions %>%
  arrange(desc(combined_pvalue))  %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex)) %>%
  pull(interaction_label)

interaction_order <- unique(interaction_order)

heatmap_data <- heatmap_data %>%
  mutate(interaction_label = factor(interaction_label, levels = interaction_order))

p2 <- ggplot(heatmap_data, aes(x = source_target,
                               y = interaction_label,
                               fill = interaction_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0,
                       name = "Interaction\nStatistic") +
  labs(title = "Interaction Statistics Across Source-Target Pairs",
       subtitle = "Ordered by combined p-value",
       x = "Source → Target",
       y = "ligand receptor Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(fig_dir, "NonCAR_APH_Top_Myeloid_Interactions.pdf"), width = 10, height = 6)
print(p1)
print(p2)
dev.off()



# 10. Curated ligand-based figures: myeloid, CAR+ TDN-D0 # ----------------

fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[1], "_fisher_ligand_results.csv")))

top_interactions <- fishers_results %>%
  filter(source %in% myeloid.cells) %>%
  filter(!is.na(combined_pvalue)) %>%
  filter(combined_padj < 0.01) %>%
  filter(ligand_complex %in% c(
    "ARPC5",
    "B2M",
    "CD47",
    "HLA-DQA1",
    "HLA-DRB1",
    "IRAK4",
    "RTN4"
  )
  ) %>%
  filter(receptor_complex %in% c(
    "ADRB2",
    "CD4",
    "GJB2",
    "KLRD1",
    "LAG3",
    "LDLR",
    "RTN4R",
    "S1PR2",
    "SIRPG",
    "TLR6"
  )) %>%
  arrange(combined_padj) %>%
  mutate(
    interaction_label = paste0(ligand_complex, "-", receptor_complex),
    neg_log10_pvalue = -log10(combined_padj)
  )

write.csv(top_interactions, file.path(files_dir, "TDND0_CAR_Myeloid_Top_Ordered.csv"))

p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                   y = reorder(interaction_label, neg_log10_pvalue))) +
  geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
  scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3, 3),
                         name = "interaction_statistic") +
  scale_size_continuous(range = c(0, 10),
                        limits = c(2, 30),
                        name = "-log10(p-adj)") +
  labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
       subtitle = "Ordered by combined p-value (most significant at top)",
       x = "Mean Interaction Statistic",
       y = "ligand_receptor Pair") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")

top_pairs <- top_interactions %>%
  select(source, ligand_complex, receptor_complex) %>%
  distinct()

heatmap_data <- dea.car.tdnD0 %>%
  inner_join(top_pairs, by = c("source", "ligand_complex", "receptor_complex")) %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
         source_target = paste0(source, " → ", target)) 

interaction_order <- top_interactions %>%
  arrange(desc(combined_pvalue))  %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex)) %>%
  pull(interaction_label)

interaction_order <- unique(interaction_order)

heatmap_data <- heatmap_data %>%
  mutate(interaction_label = factor(interaction_label, levels = interaction_order))

p2 <- ggplot(heatmap_data, aes(x = source_target,
                               y = interaction_label,
                               fill = interaction_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0,
                       name = "Interaction\nStatistic") +
  labs(title = "Interaction Statistics Across Source-Target Pairs",
       subtitle = "Ordered by combined p-value",
       x = "Source → Target",
       y = "ligand receptor Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(fig_dir, "CAR_TDN_D0_Top_Myeloid_Interactions.pdf"), width = 10, height = 6)
print(p1)
print(p2)
dev.off()



# 11. Curated ligand-based figures: myeloid, CAR+ Peak # ----------------

fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[3], "_fisher_ligand_results.csv")))

top_interactions <- fishers_results %>%
  filter(source %in% myeloid.cells) %>%
  filter(!is.na(combined_pvalue)) %>%
  filter(combined_padj < 0.01) %>%
  filter(ligand_complex %in% c(
    "CALM1",
    "CCL4",
    "FCER2",
    "IL1B",
    "JAG1",
    "TGS1",
    "VEGFA"
  )
  ) %>%
  filter(receptor_complex %in% c(
    "ADRB2",
    "CCR1",
    "CCR5",
    "CCR8",
    "CD44",
    "CD46",
    "FAS",
    "ITGAM_ITGB2",
    "ITGAV",
    "ITGAV_ITGB3",
    "ITGAX_ITGB2",
    "ITGB1",
    "KCNN4",
    "NOTCH1",
    "NOTCH2",
    "NRP2",
    "PTPRA",
    "RXRA",
    "SIGIRR"
  )
  
  ) %>%
  arrange(combined_padj) %>%
  mutate(
    interaction_label = paste0(ligand_complex, "-", receptor_complex),
    neg_log10_pvalue = -log10(combined_padj)
  )

write.csv(top_interactions, file.path(files_dir, "Peak_CAR_Myeloid_Top_Ordered.csv"))

p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                   y = reorder(interaction_label, neg_log10_pvalue))) +
  geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
  scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3.2, 3),
                         name = "interaction_statistic") +
  scale_size_continuous(range = c(0, 10),
                        limits = c(2, 30),
                        name = "-log10(p-adj)") +
  labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
       subtitle = "Ordered by combined p-value (most significant at top)",
       x = "Mean Interaction Statistic",
       y = "ligand_receptor Pair") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")

top_pairs <- top_interactions %>%
  select(source, ligand_complex, receptor_complex) %>%
  distinct()

heatmap_data <- dea.car.peak %>%
  inner_join(top_pairs, by = c("source", "ligand_complex", "receptor_complex")) %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
         source_target = paste0(source, " → ", target)) 

interaction_order <- top_interactions %>%
  arrange(desc(combined_pvalue))  %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex)) %>%
  pull(interaction_label)

interaction_order <- unique(interaction_order)

heatmap_data <- heatmap_data %>%
  mutate(interaction_label = factor(interaction_label, levels = interaction_order))

p2 <- ggplot(heatmap_data, aes(x = source_target,
                               y = interaction_label,
                               fill = interaction_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0,
                       name = "Interaction\nStatistic") +
  labs(title = "Interaction Statistics Across Source-Target Pairs",
       subtitle = "Ordered by combined p-value",
       x = "Source → Target",
       y = "ligand receptor Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(fig_dir, "CAR_Peak_Top_Myeloid_Interactions.pdf"), width = 10, height = 6)
print(p1)
print(p2)
dev.off()



# 12. Concatenate Fisher ligand results across timepoints # ----------------

fishers_results_tdnd0 <- read_csv(file.path(dea_dir, paste0(conditions[1], "_fisher_ligand_results.csv")))
fishers_results_aph   <- read_csv(file.path(dea_dir, paste0(conditions[2], "_fisher_ligand_results.csv")))
fishers_results_peak  <- read_csv(file.path(dea_dir, paste0(conditions[3], "_fisher_ligand_results.csv")))
fishers_results_week4 <- read_csv(file.path(dea_dir, paste0(conditions[4], "_fisher_ligand_results.csv")))

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


# 13. Line graph: TNF / LTA / LTB signaling over time # ---------------

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


# 14. Exploratory: top receptor-based interactions per condition x cell group # ----------------

list_celltypes <- list(t.cells, myeloid.cells)

for (i in 1:length(conditions)){
  
  fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[i], "_fisher_receptor_results.csv")))
  dea.raa <- dataframes[[i]]
  
  for (c in list_celltypes){
    print(c)
    celltype <- c[1]
    if ("CD14 Mono" %in% c){
      celltype <- "All_Myeloid"
    } else if ("Dendritic Cell" %in% c){
      celltype <- "DC"
    } else if ("CD4 EM-like" %in% c){
      celltype <- "All_TCell"
    }
    top_interactions <- fishers_results %>%
      filter(target %in% c) %>%
      filter(!is.na(combined_padj)) %>%
      filter(combined_padj < 0.05) %>%
      arrange(combined_padj) %>%
      head(20) %>%
      mutate(
        interaction_label = paste0(ligand_complex, "-", receptor_complex),
        neg_log10_pvalue = -log10(combined_padj)
      )
    
    p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                       y = reorder(interaction_label, neg_log10_pvalue))) +
      geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
      scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3, 3),
                             name = "interaction_statistic") +
      scale_size_continuous(range = c(0, 10),
                            limits = c(2, 30),
                            name = "-log10(p-adj)") +
      labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
           subtitle = "Ordered by combined p-value (most significant at top)",
           x = "Mean Interaction Statistic",
           y = "ligand_receptor Pair") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 10),
            legend.position = "right")
    
    top_pairs <- top_interactions %>%
      select(ligand_complex, receptor_complex) %>%
      distinct()
    
    pdf(file.path(fig_dir, paste0(conditions[i], "_", celltype, "_fisher_receptor_images.pdf")),
        width = 15, height = 12)
    print(p1)
    dev.off()
    
  }
}


# 15. Curated receptor-based figures: T cells,  APH # ----------------

fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[2], "_fisher_receptor_results.csv")))

top_interactions <- fishers_results %>%
  filter(target %in% t.cells) %>%
  filter(!is.na(combined_pvalue)) %>%
  filter(combined_padj < 0.05) %>%
  filter(ligand_complex %in% c(
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-E",
    "HLA-F",
    "CD86",
    "AFDN",
    "GNAI2",
    "CXCL16",
    "IL23A",
    "HLA-G",
    "B2M",
    "CLEC2B"
  )) %>%
  filter(receptor_complex %in% c(
    "CD8B",
    "CTLA4",
    "EPHB6",
    "S1PR5",
    "CXCR6",
    "IL12RB2",
    "CD8A",
    "KIR2DL1",
    "KLRF1"
  )) %>%
  arrange(combined_padj) %>%
  mutate(
    interaction_label = paste0(ligand_complex, "-", receptor_complex),
    neg_log10_pvalue = -log10(combined_padj)
  )

write.csv(top_interactions, file.path(files_dir, "APH_NonCAR_Tcell_Top_Ordered_receptor_based.csv"))

p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                   y = reorder(interaction_label, neg_log10_pvalue))) +
  geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
  scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3, 3),
                         name = "interaction_statistic") +
  scale_size_continuous(range = c(0, 15),
                        limits = c(2, 30),
                        name = "-log10(p-value)") +
  labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
       subtitle = "Ordered by combined p-value (most significant at top)",
       x = "Mean Interaction Statistic",
       y = "ligand_receptor Pair") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")

top_pairs <- top_interactions %>%
  select(target, ligand_complex, receptor_complex) %>%
  distinct()

heatmap_data <- dea.noncar.aph %>%
  inner_join(top_pairs, by = c("target", "ligand_complex", "receptor_complex")) %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
         source_target = paste0(source, " → ", target))

interaction_order <- top_interactions %>%
  arrange(desc(combined_pvalue))  %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex)) %>%
  pull(interaction_label)

interaction_order <- unique(interaction_order)

heatmap_data <- heatmap_data %>%
  mutate(interaction_label = factor(interaction_label, levels = interaction_order))

p2 <- ggplot(heatmap_data, aes(x = source_target,
                               y = interaction_label,
                               fill = interaction_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0,
                       name = "Interaction\nStatistic") +
  labs(title = "Interaction Statistics Across Source-Target Pairs",
       subtitle = "Ordered by combined p-value",
       x = "Source → Target",
       y = "ligand receptor Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(fig_dir, "NonCAR_APH_Top_TCell_ReceptorBased_Interactions.pdf"), width = 10, height = 6)
print(p1)
print(p2)
dev.off()


# 16. Curated receptor-based figures: T cells, CAR+ TDN-D0 # ----------------

fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[1], "_fisher_receptor_results.csv")))

top_interactions <- fishers_results %>%
  filter(target %in% t.cells) %>%
  filter(!is.na(combined_pvalue)) %>%
  filter(combined_padj < 0.01) %>%
  filter(ligand_complex %in% c(
    "TNFSF12",
    "HMGB1",
    "ZG16B",
    "SELPLG",
    "SPON2",
    "FADD",
    "TNF",
    "TGFB1",
    "ACTR2",
    "ARPC5",
    "GNAS",
    "HSPA8",
    "POMC",
    "VEGFB",
    "CLEC2D",
    "LGALS9"
  )
  ) %>%
  filter(receptor_complex %in% c(
    "TNFRSF25",
    "CXCR4",
    "ITGAM",
    "TRADD",
    "VSIR",
    "ADRB2",
    "KLRB1",
    "SLC1A5"
  )) %>%
  arrange(combined_padj) %>%
  mutate(
    interaction_label = paste0(ligand_complex, "-", receptor_complex),
    neg_log10_pvalue = -log10(combined_padj)
  )

write.csv(top_interactions, file.path(files_dir, "TDND0_CAR_Tcell_Top_Ordered_ReceptorBased.csv"))

p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                   y = reorder(interaction_label, neg_log10_pvalue))) +
  geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
  scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3, 3),
                         name = "interaction_statistic") +
  scale_size_continuous(range = c(0, 10),
                        limits = c(2, 30),
                        name = "-log10(p-adj)") +
  labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
       subtitle = "Ordered by combined p-value (most significant at top)",
       x = "Mean Interaction Statistic",
       y = "ligand_receptor Pair") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")

top_pairs <- top_interactions %>%
  select(target, ligand_complex, receptor_complex) %>%
  distinct()

heatmap_data <- dea.car.tdnD0 %>%
  inner_join(top_pairs, by = c("target", "ligand_complex", "receptor_complex")) %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
         source_target = paste0(source, " → ", target))

interaction_order <- top_interactions %>%
  arrange(desc(combined_pvalue))  %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex)) %>%
  pull(interaction_label)

interaction_order <- unique(interaction_order)

heatmap_data <- heatmap_data %>%
  mutate(interaction_label = factor(interaction_label, levels = interaction_order))

p2 <- ggplot(heatmap_data, aes(x = source_target,
                               y = interaction_label,
                               fill = interaction_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0,
                       name = "Interaction\nStatistic") +
  labs(title = "Interaction Statistics Across Source-Target Pairs",
       subtitle = "Ordered by combined p-value",
       x = "Source → Target",
       y = "ligand receptor Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(fig_dir, "CAR_TDN_D0_Top_TCell_ReceptorBased_Interactions.pdf"), width = 10, height = 6)
print(p1)
print(p2)
dev.off()


# 17. Curated receptor-based figures: T cells, CAR+ Peak # ----------------

fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[3], "_fisher_receptor_results.csv")))

top_interactions <- fishers_results %>%
  filter(target %in% t.cells) %>%
  filter(!is.na(combined_pvalue)) %>%
  filter(combined_padj < 0.01) %>%
  filter(ligand_complex %in% c(
    "CRTAM",
    "GNAI2",
    "HMGB1",
    "TGFB1",
    "ADAM10",
    "CCL28",
    "CCL3",
    "B2M",
    "HLA-G",
    "CCL4",
    "CCL5",
    "CALM1",
    "CALM3",
    "LGALS1",
    "CLEC2D",
    "S100A8"
  )
  
  ) %>%
  filter(receptor_complex %in% c(
    "CADM1",
    "TBXA2R",
    "CXCR4",
    "CCR3",
    "KLRC1",
    "CCR8",
    "KCNQ5",
    "CD69",
    "KLRB1"
  )) %>%
  arrange(combined_padj) %>%
  mutate(
    interaction_label = paste0(ligand_complex, "-", receptor_complex),
    neg_log10_pvalue = -log10(combined_padj)
  )

write.csv(top_interactions, file.path(files_dir, "Peak_CAR_Tcell_Top_Ordered_ReceptorBased.csv"))

p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                   y = reorder(interaction_label, neg_log10_pvalue))) +
  geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
  scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3, 3),
                         name = "interaction_statistic") +
  scale_size_continuous(range = c(0, 10),
                        limits = c(2, 30),
                        name = "-log10(p-adj)") +
  labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
       subtitle = "Ordered by combined p-value (most significant at top)",
       x = "Mean Interaction Statistic",
       y = "ligand_receptor Pair") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")

top_pairs <- top_interactions %>%
  select(target, ligand_complex, receptor_complex) %>%
  distinct()

heatmap_data <- dea.car.peak %>%
  inner_join(top_pairs, by = c("target", "ligand_complex", "receptor_complex")) %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
         source_target = paste0(source, " → ", target))

interaction_order <- top_interactions %>%
  arrange(desc(combined_pvalue))  %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex)) %>%
  pull(interaction_label)

interaction_order <- unique(interaction_order)

heatmap_data <- heatmap_data %>%
  mutate(interaction_label = factor(interaction_label, levels = interaction_order))

p2 <- ggplot(heatmap_data, aes(x = source_target,
                               y = interaction_label,
                               fill = interaction_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0,
                       name = "Interaction\nStatistic") +
  labs(title = "Interaction Statistics Across Source-Target Pairs",
       subtitle = "Ordered by combined p-value",
       x = "Source → Target",
       y = "ligand receptor Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(fig_dir, "CAR_Peak_Top_TCell_ReceptorBased_Interactions.pdf"), width = 10, height = 6)
print(p1)
print(p2)
dev.off()


# 18. Curated receptor-based figures: myeloid,  APH # ----------------

fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[2], "_fisher_receptor_results.csv")))

top_interactions <- fishers_results %>%
  filter(target %in% myeloid.cells) %>%
  filter(!is.na(combined_pvalue)) %>%
  filter(combined_padj < 0.01) %>%
  filter(ligand_complex %in% c(
    "ADAM15",
    "ADAM17",
    "TNFSF13B",
    "ICAM1",
    "S100A8",
    "S100A9",
    "TNF",
    "F8",
    "ENTPD1",
    "LTA_LTB",
    "ARF1",
    "CALM1",
    "CALM3",
    "HRAS",
    "NAMPT",
    "CD59",
    "TNFSF14",
    "ADO",
    "FABP5",
    "TGS1"
  )) %>%
  filter(receptor_complex %in%  c(
    "ITGA5",
    "HLA-DPB1",
    "IL2RG",
    "CD68",
    "TRPM2",
    "ASGR2",
    "ADORA2B",
    "LTBR",
    "INSR",
    "STAB1",
    "RXRA"
  )) %>%
  arrange(combined_padj) %>%
  mutate(
    interaction_label = paste0(ligand_complex, "-", receptor_complex),
    neg_log10_pvalue = -log10(combined_padj)
  )

write.csv(top_interactions, file.path(files_dir, "APH_NonCAR_Myeloid_Top_Ordered_receptor_based.csv"))

p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                   y = reorder(interaction_label, neg_log10_pvalue))) +
  geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
  scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3, 3),
                         name = "interaction_statistic") +
  scale_size_continuous(range = c(0, 15),
                        limits = c(2, 30),
                        name = "-log10(p-value)") +
  labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
       subtitle = "Ordered by combined p-value (most significant at top)",
       x = "Mean Interaction Statistic",
       y = "ligand_receptor Pair") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")

top_pairs <- top_interactions %>%
  select(target, ligand_complex, receptor_complex) %>%
  distinct()

heatmap_data <- dea.noncar.aph %>%
  inner_join(top_pairs, by = c("target", "ligand_complex", "receptor_complex")) %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
         source_target = paste0(source, " → ", target))

interaction_order <- top_interactions %>%
  arrange(desc(combined_pvalue))  %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex)) %>%
  pull(interaction_label)

interaction_order <- unique(interaction_order)

heatmap_data <- heatmap_data %>%
  mutate(interaction_label = factor(interaction_label, levels = interaction_order))

p2 <- ggplot(heatmap_data, aes(x = source_target,
                               y = interaction_label,
                               fill = interaction_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0,
                       name = "Interaction\nStatistic") +
  labs(title = "Interaction Statistics Across Source-Target Pairs",
       subtitle = "Ordered by combined p-value",
       x = "Source → Target",
       y = "ligand receptor Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(fig_dir, "NonCAR_APH_Top_Myeloid_ReceptorBased_Interactions.pdf"), width = 10, height = 6)
print(p1)
print(p2)
dev.off()


# 19. Curated receptor-based figures: myeloid, CAR+ TDN-D0 # ----------------

fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[1], "_fisher_receptor_results.csv")))

top_interactions <- fishers_results %>%
  filter(target %in% myeloid.cells) %>%
  filter(!is.na(combined_pvalue)) %>%
  filter(combined_padj < 0.01) %>%
  filter(ligand_complex %in%  c(
    "APP",
    "TNFSF10",
    "GNAI2",
    "PTPN6",
    "NCAM1",
    "SIRPG",
    "ZG16B",
    "LGALS9",
    "HMGB1",
    "TGFB1",
    "F8",
    "FAM3C",
    "TNF",
    "NELL2"
  )
  ) %>%
  filter(receptor_complex %in% c(
    "LRP1",
    "TNFRSF10C",
    "S1PR3",
    "CD300LF",
    "ROBO3",
    "CD47",
    "CXCR4",
    "ASGR2",
    "FFAR2",
    "PTGER2"
  )) %>%
  arrange(combined_padj) %>%
  mutate(
    interaction_label = paste0(ligand_complex, "-", receptor_complex),
    neg_log10_pvalue = -log10(combined_padj)
  )

write.csv(top_interactions, file.path(files_dir, "TDND0_CAR_Myeloid_Top_Ordered_ReceptorBased.csv"))

p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                   y = reorder(interaction_label, neg_log10_pvalue))) +
  geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
  scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3, 3),
                         name = "interaction_statistic") +
  scale_size_continuous(range = c(0, 10),
                        limits = c(2, 30),
                        name = "-log10(p-adj)") +
  labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
       subtitle = "Ordered by combined p-value (most significant at top)",
       x = "Mean Interaction Statistic",
       y = "ligand_receptor Pair") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")

top_pairs <- top_interactions %>%
  select(target, ligand_complex, receptor_complex) %>%
  distinct()

heatmap_data <- dea.car.tdnD0 %>%
  inner_join(top_pairs, by = c("target", "ligand_complex", "receptor_complex")) %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
         source_target = paste0(source, " → ", target))

interaction_order <- top_interactions %>%
  arrange(desc(combined_pvalue))  %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex)) %>%
  pull(interaction_label)

interaction_order <- unique(interaction_order)

heatmap_data <- heatmap_data %>%
  mutate(interaction_label = factor(interaction_label, levels = interaction_order))

p2 <- ggplot(heatmap_data, aes(x = source_target,
                               y = interaction_label,
                               fill = interaction_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0,
                       name = "Interaction\nStatistic") +
  labs(title = "Interaction Statistics Across Source-Target Pairs",
       subtitle = "Ordered by combined p-value",
       x = "Source → Target",
       y = "ligand receptor Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(fig_dir, "CAR_TDN_D0_Top_Myeloid_ReceptorBased_Interactions.pdf"), width = 10, height = 6)
print(p1)
print(p2)
dev.off()


# 20. Curated receptor-based figures: myeloid, CAR+ Peak # ----------------

fishers_results <- read_csv(file.path(dea_dir, paste0(conditions[3], "_fisher_receptor_results.csv")))

top_interactions <- fishers_results %>%
  filter(target %in% myeloid.cells) %>%
  filter(!is.na(combined_pvalue)) %>%
  filter(combined_padj < 0.01) %>%
  filter(ligand_complex %in% c(
    "ARPC5",
    "GNAS",
    "HSPA8",
    "POMC",
    "VEGFB",
    "CALM1",
    "CALM3",
    "LGALS3",
    "TGFB1",
    "GSTO1",
    "APOA1",
    "FADD",
    "LIN7C",
    "ACTR2"
  )
  
  ) %>%
  filter(receptor_complex %in% c(
    "ADRB2",
    "PDE1B",
    "ENG",
    "RYR1",
    "ABCA1"
  )) %>%
  arrange(combined_padj) %>%
  mutate(
    interaction_label = paste0(ligand_complex, "-", receptor_complex),
    neg_log10_pvalue = -log10(combined_padj)
  )

write.csv(top_interactions, file.path(files_dir, "Peak_CAR_Myeloid_Top_Ordered_ReceptorBased.csv"))

p1 <- ggplot(top_interactions, aes(x = mean_interaction_stat,
                                   y = reorder(interaction_label, neg_log10_pvalue))) +
  geom_point(aes(size = neg_log10_pvalue, color = mean_interaction_stat)) +
  scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, limits = c(-3, 3),
                         name = "interaction_statistic") +
  scale_size_continuous(range = c(0, 10),
                        limits = c(2, 30),
                        name = "-log10(p-adj)") +
  labs(title = "Top Significant ligand_complex-receptor_complex Interactions",
       subtitle = "Ordered by combined p-value (most significant at top)",
       x = "Mean Interaction Statistic",
       y = "ligand_receptor Pair") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")

top_pairs <- top_interactions %>%
  select(target, ligand_complex, receptor_complex) %>%
  distinct()

heatmap_data <- dea.car.peak %>%
  inner_join(top_pairs, by = c("target", "ligand_complex", "receptor_complex")) %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex),
         source_target = paste0(source, " → ", target))

interaction_order <- top_interactions %>%
  arrange(desc(combined_pvalue))  %>%
  mutate(interaction_label = paste0(ligand_complex, "-", receptor_complex)) %>%
  pull(interaction_label)

interaction_order <- unique(interaction_order)

heatmap_data <- heatmap_data %>%
  mutate(interaction_label = factor(interaction_label, levels = interaction_order))

p2 <- ggplot(heatmap_data, aes(x = source_target,
                               y = interaction_label,
                               fill = interaction_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0,
                       name = "Interaction\nStatistic") +
  labs(title = "Interaction Statistics Across Source-Target Pairs",
       subtitle = "Ordered by combined p-value",
       x = "Source → Target",
       y = "ligand receptor Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(fig_dir, "CAR_Peak_Top_Myeloid_ReceptorBased_Interactions.pdf"), width = 10, height = 6)
print(p1)
print(p2)
dev.off()


# 21. Interferon-stimulated ligand sets and pathway-level temporal dynamics # ----------------

IFN_I_high <- c(
  "CXCL10","TNFSF10",   # TRAIL
  "TNFSF13B",  # BAFF
  "TNFSF13",   # APRIL
  "BST2","B2M","CD274",     # PD-L1
  "NCR3LG1",   # B7-H6
  "MICA", "MICB","ULBP1","ULBP2","ULBP3")

IFN_I_moderate <- c(
  "CCL5","CCL3","CCL4",
  "CXCL1","CXCL2","CXCL3","CXCL16",
  "IL15","IL15RA","IL10","IL32","CSF1",
  "FLT3LG","LIF","OSM","ICAM1",
  "VCAM1","CD86","CD80","CD70",
  "CD48","CD58","HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G")

IFN_II_high <- c(
  "CXCL9", "CXCL10","CXCL11","CD274", "TNFSF10",
  "ICAM1","VCAM1","B2M","HLA-A","HLA-B","HLA-C",
  "MICA","MICB","ULBP1","ULBP2","ULBP3", "IFNG")

IFN_II_moderate <- c(
  "CCL5","CXCL16","IL15","IL15RA","IL32",
  "CD86", "CD80", "CD40LG","CD70",
  "TNF","IL18","HLA-DRA","HLA-DRB1","HLA-DQA1",
  "HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA","HLA-DMB",
  "HLA-DOB")

# ---- Type I IFN -------------------------------------------------------------
fishers_type1 <- fishers_all %>% filter(ligand_complex %in% c(IFN_I_high, IFN_I_moderate)) %>%
  filter(combined_padj < 0.01) %>%
  filter(source %in% c(myeloid.cells, "CD8 EM", "CD4 EM-like", "Proliferating")) %>%
  mutate(interaction = paste(ligand_complex, receptor_complex, sep = "_"),
         sig = combined_padj < 0.01)

fishers_type1$timepoint <- factor(fishers_type1$timepoint, levels = c("APH", "TDN-D0", "Peak", "Week4"))
fishers_type1$pathway <- "IFN1"
fishers_type1$response <- rep("CR", length(fishers_type1$interaction))
fishers_type1$response[fishers_type1$mean_interaction_stat < 0] <- "PD"

IFN1_interactions <- c(
  # High confidence
  "CXCL10_CXCR3", "CXCL10_SDC4", "CXCL10_DPP4",
  "TNFSF10_TNFRSF10A", "TNFSF10_TNFRSF10B",
  "TNFSF10_TNFRSF10C", "TNFSF10_TNFRSF10D",
  "TNFSF10_RIPK1",
  "TNFSF13B_TNFRSF13C", "TNFSF13B_CD40",
  "BST2_LILRA4",
  "ULBP2_HCST_KLRK1", "MICB_HCST_KLRK1", "MICA_HCST_KLRK1",
  
  # Moderate confidence
  "CCL5_CCR1", "CCL5_CCR4", "CCL5_CCR5", "CCL5_CXCR3", "CCL5_SDC4",
  "CCL3_CCR1", "CCL3_CCR4", "CCL3_CCR5", "CCL3_CCR3",
  "CCL4_CCR1", "CCL4_CCR3", "CCL4_CCR5", "CCL4_CCR8", "CCL4_ACKR2",
  "IL15_IL2RA", "IL15_IL15RA_IL2RB_IL2RG",
  "FLT3LG_FLT3",
  "ICAM1_ITGAL_ITGB2", "ICAM1_ITGAM_ITGB2", "ICAM1_ITGAX_ITGB2",
  "IFNG_IFNGR1_IFNGR2", "IFNG_IFNGR1", "IFNG_IFNGR2"
)

fishers_type1_clean <- fishers_type1 %>% filter(interaction %in% IFN1_interactions)

# ---- Type II IFN ------------------------------------------------------------
fishers_type2 <- fishers_all %>% filter(ligand_complex %in% c(IFN_II_high, IFN_II_moderate)) %>%
  filter(combined_padj < 0.01) %>%
  filter(source %in% c(myeloid.cells, "CD8 EM", "CD4 EM-like", "Proliferating")) %>%
  mutate(interaction = paste(ligand_complex, receptor_complex, sep = "_"),
         sig = combined_padj < 0.01)

fishers_type2$timepoint <- factor(fishers_type2$timepoint, levels = c("APH", "TDN-D0", "Peak", "Week4"))
fishers_type2$pathway <- "IFN2"
fishers_type2$response <- rep("CR", length(fishers_type2$interaction))
fishers_type2$response[fishers_type2$mean_interaction_stat < 0] <- "PD"

ifng_interactions <- c("HLA-DRB1_CD4","HLA-DQA1_CD4","HLA-DRA_CD4","HLA-DPB1_CD4","HLA-DPA1_CD4","HLA-DQB1_CD4",
                       "HLA-DRB1_LAG3","HLA-DQA1_LAG3","HLA-DRA_LAG3","HLA-DPB1_LAG3","HLA-DPA1_LAG3","HLA-DQB1_LAG3",
                       "HLA-A_CD8A","HLA-A_CD8B","HLA-A_CD8B2","HLA-B_CD8A","HLA-B_CD8B","HLA-B_CD8B2","HLA-C_CD8A","HLA-C_CD8B","HLA-C_CD8B2",
                       "HLA-A_LILRB1","HLA-A_LILRB2","HLA-A_LILRA1","HLA-B_LILRB1","HLA-B_LILRB2","HLA-B_LILRA1","HLA-C_LILRB1","HLA-C_LILRB2","HLA-C_LILRA1",
                       "MICA_HCST_KLRK1","MICB_HCST_KLRK1","ULBP2_HCST_KLRK1",
                       "CXCL10_CXCR3","CXCL10_SDC4","CXCL10_DPP4","CCL5_CXCR3","CCL5_CCR5","CCL5_SDC4",
                       "TNF_TNFRSF1A","TNF_TNFRSF1B","TNF_TRADD","TNF_TRAF2","TNF_FAS",
                       "TNFSF10_TNFRSF10A","TNFSF10_TNFRSF10B","TNFSF10_TNFRSF10C","TNFSF10_TNFRSF10D","TNFSF10_RIPK1",
                       "IL18_IL18R1_IL18RAP","IL18_IL18BP",
                       "IL15_IL2RA","IL15_IL15RA_IL2RB_IL2RG",
                       "ICAM1_ITGAL_ITGB2","ICAM1_ITGAM_ITGB2","ICAM1_ITGAX_ITGB2","ICAM1_SPN")

fishers_type2_clean <- fishers_type2 %>% filter(interaction %in% ifng_interactions)

# ---- TNF --------------------------------------------------------------------
fishers_tnf <- fishers_all %>% filter(ligand_complex %in% c("TNF", "LTB", "LTA")) %>%
  filter(receptor_complex %in% c("TNFRSF1A", "TNFRSF1B", "TRADD", "TRAF2", "TRPM2", "FFAR2",
                                 "RIPK1")) %>%
  filter(source %in% t.cells) %>%
  filter(combined_padj < 0.01) %>%
  mutate(interaction = paste(ligand, receptor, sep = "_"),
         sig = combined_padj < 0.01)

fishers_tnf$timepoint <- factor(fishers_tnf$timepoint, levels = c("APH", "TDN-D0", "Peak", "Week4"))
fishers_tnf$pathway <- "TNF"
fishers_tnf$response <- rep("CR", length(fishers_tnf$interaction))
fishers_tnf$response[fishers_tnf$mean_interaction_stat < 0] <- "PD"

# ---- Interaction-level combined plot ----------------------------------------
fishers_combined_interaction <- do.call(rbind, list(fishers_type1_clean, fishers_tnf, fishers_type2_clean))

fishers_combined_interaction <- fishers_combined_interaction  %>%
  mutate(signed_logp = sign(mean_interaction_stat) * -log10(combined_padj)) %>%
  group_by(timepoint, interaction,  pathway, response, ligand_complex) %>%
  summarise(
    mean_interaction_stat = mean(mean_interaction_stat, na.rm = TRUE),
    pathway_score = mean(signed_logp, na.rm = TRUE),
    .groups = "drop"
  )

fishers_combined_interaction <- fishers_combined_interaction %>% filter(abs(pathway_score) > -log10(0.05))


# ---- Pathway-level summaries ------------------------------------------------
my_colors <- c(
  "IFN1" = "cornflowerblue",
  "IFN1_IFN2" = "dodgerblue",
  "IFN2" = "blue1",
  "TNF" = "#e31a1c",       # red
  "IFN2_TNF" = "magenta"
)


# same, with 95% confidence intervals
pathway_summary_score_ci <- fishers_combined_interaction %>%
  group_by(timepoint, pathway) %>%
  summarise(
    mean_score = mean(pathway_score, na.rm = TRUE),
    sd_score   = sd(pathway_score, na.rm = TRUE),
    n_interactions = sum(!is.na(pathway_score)),
    se_score   = sd_score / sqrt(n_interactions),
    ci_lower   = mean_score - qt(0.95, df = n_interactions - 1) * se_score,
    ci_upper   = mean_score + qt(0.95, df = n_interactions - 1) * se_score,
    .groups = "drop"
  )



# signed -log10(padj), CI ribbon
p6.ci <- ggplot(pathway_summary_score_ci,
                aes(x = timepoint,
                    y = mean_score,
                    color = pathway,
                    group = pathway)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_ribbon(aes(ymin = ci_lower,
                  ymax = ci_upper,
                  fill = pathway),
              alpha = 0.2,
              color = NA) +
  theme_classic() +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    y = "Mean signed -log10(p-adjusted) ",
    x = "Timepoint",
    title = "Temporal dynamics of signaling by combined pathway"
  )



pdf(file.path(fig_dir, "Opposing_Interactions_PathwayLevel_Linegraph.pdf"), width = 10, height = 5)
print(p6.ci)
dev.off()