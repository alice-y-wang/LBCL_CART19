################################################################################
## TCR repertoire analysis: diversity metrics, clonal expansion, and testing
################################################################################

# ---- Libraries --------------------------------------------------------------
library(tidyverse)
library(readr)
library(ggplot2)
library(ggpubr)
library(ggrastr)
library(rstatix)

# ---- CONFIG -----------------------------------------------------------------
analysis_dir <- "."
setwd(analysis_dir)

files_dir  <- "files"
tcr_dir    <- "TCR_Data"
stats_dir  <- file.path(tcr_dir, "wilcox_test")
boxplot_dir <- file.path(tcr_dir, "unpaired_boxplot")

for (d in c(tcr_dir, stats_dir, boxplot_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

source(file.path("TCR_Code", "00_functions.R"))
set.seed(11)

# Response labels. recode() leaves values that are already CR/PD unchanged.
response_levels <- c("CR", "PD")

# Timepoint order used by the TCR metadata
timepoint_levels <- c("Apheresis", "Product", "Day 0", "Peak", "Week 4")

# Plot colours
response_colors <- c("CR" = "#4FB7C5", "PD" = "#E7A75E")
car_colors      <- c("CAR+" = "steelblue2", "CAR-" = "gray")
celltype_colors <- c("CD4 T" = "skyblue2", "CD8 T" = "darkolivegreen3")


# 1. Load TCR metadata # -----------

tcrs_metadata = read_csv(file.path(files_dir, "LBCL_TCR_metadata_clean_forupload.csv")) %>%
  mutate(response = recode(as.character(response),
                           "Responder" = "CR", "Non-Responder" = "PD"),
         timepoint = factor(timepoint, levels = timepoint_levels),
         response = factor(response, levels = c("CR", "PD")))


# 2. Clone size and clonal homeostasis tables # ------


# clonal expansion
df_clonal_homeostasis = tcrs_metadata %>%
  filter(has_tcr_gex) %>%
  group_by(patient_id, response, CTaa) %>%
  summarise(clone_size = n(), .groups = "drop_last") %>%
  mutate(freq = clone_size / sum(clone_size)) %>%   # within patient
  ungroup() %>%
  mutate(bin = case_when(
    freq > 0.1 ~ "Hyperexpanded (0.1 < X <= 1)",
    freq > 0.01 ~ "Large (0.01 < X <= 0.1)",
    freq > 0.001 ~ "Medium (0.001 < X <= 0.01)",
    freq > 0.0001 ~ "Small (0.0001 < X <= 0.001)",
    freq > 0.00001 ~ "Rare (0.00001 < X <= 0.0001)",
    TRUE ~ NA_character_
  ))



# 3. Diversity metrics split by CD4 / CD8 # -------
# ---- CAR+ -------------------------------------------------------------------
meta_to_use = tcrs_metadata %>%
  filter(has_tcr_gex,
         CAR_value == 'CAR+',
         tcell_label %in% c('CD4 T', 'CD8 T'))

# Metrics: richness, shannon entropy, gini-simpson, gini index, evenness, dominance
metrics_df_CAR_CD4_CD8 = data.frame()
meta_to_use$sample2 = paste(meta_to_use$patient_id, meta_to_use$timepoint, sep="_")

for(celltype in c('CD4 T', 'CD8 T')){
  df = calculate_metrics(meta_to_use %>% filter(tcell_label == celltype), min_size = 200)
  df$celltype = celltype
  metrics_df_CAR_CD4_CD8 = rbind(metrics_df_CAR_CD4_CD8, df)
}

metrics_df_CAR_CD4_CD8 = metrics_df_CAR_CD4_CD8 %>%
  mutate(
         celltype = factor(celltype, levels = c("CD4 T", "CD8 T"))) %>%
  relocate(response, .after = patient_id)

write_csv(metrics_df_CAR_CD4_CD8, file.path(tcr_dir, "diversity_metrics_CAR_CD4_CD8_Tcells.csv"))

# ---- CAR- -------------------------------------------------------------------
meta_to_use = tcrs_metadata %>%
  filter(has_tcr_gex,
         CAR_value == 'CAR-',
         tcell_label %in% c('CD4 T', 'CD8 T'))

# Metrics: richness, shannon entropy, gini-simpson, gini index, evenness, dominance
metrics_df_noncar_CD4_CD8 = data.frame()
meta_to_use$sample2 = paste(meta_to_use$patient_id, meta_to_use$timepoint, sep="_")

for(celltype in c('CD4 T', 'CD8 T')){
  df = calculate_metrics(meta_to_use %>% filter(tcell_label == celltype), min_size = 200)
  df$celltype = celltype
  metrics_df_noncar_CD4_CD8 = rbind(metrics_df_noncar_CD4_CD8, df)
}

metrics_df_noncar_CD4_CD8 = metrics_df_noncar_CD4_CD8 %>%
  mutate(celltype = factor(celltype, levels = c("CD4 T", "CD8 T"))) %>%
  relocate(response, .after = patient_id)

write_csv(metrics_df_noncar_CD4_CD8, file.path(tcr_dir, "diversity_metrics_noncar_CD4_CD8_Tcells.csv"))


# 4. UMAP of hyperexpanded clones # -------------

df = tcrs_metadata %>%
  filter(has_tcr_gex) %>%
  group_by(patient_id, response, CTaa) %>%
  summarise(clone_size = n(), .groups = "drop_last") %>%
  mutate(freq = clone_size / sum(clone_size)) %>%   # within patient
  ungroup() %>%
  mutate(bin = case_when(
    freq > 0.1 ~ "Hyperexpanded (0.1 < X <= 1)",
    freq > 0.01 ~ "Large (0.01 < X <= 0.1)",
    freq > 0.001 ~ "Medium (0.001 < X <= 0.01)",
    freq > 0.0001 ~ "Small (0.0001 < X <= 0.001)",
    freq > 0.00001 ~ "Rare (0.00001 < X <= 0.0001)",
    TRUE ~ NA_character_
  ),
  bin = factor(bin, levels = c("Hyperexpanded (0.1 < X <= 1)", "Large (0.01 < X <= 0.1)", "Medium (0.001 < X <= 0.01)", "Small (0.0001 < X <= 0.001)", "Rare (0.00001 < X <= 0.0001)"))) %>%
  ungroup()

tcrs_plot = tcrs_metadata %>%
  filter(has_tcr_gex) %>%
  left_join(df %>% dplyr::select(patient_id, CTaa, bin, freq), by = c("patient_id", "CTaa"))

tcrs_plot$bin <- factor(
  tcrs_plot$bin,
  levels = c(
    "Rare (0.00001 < X <= 0.0001)",
    "Small (0.0001 < X <= 0.001)",
    "Medium (0.001 < X <= 0.01)",
    "Large (0.01 < X <= 0.1)",
    "Hyperexpanded (0.1 < X <= 1)"
  )
)

tcrs_plot <- tcrs_plot[order(tcrs_plot$bin), ]

# Umap coloured by singletons vs expanded clones split by timepoint
pdf(file.path(tcr_dir, "umaps_by_hyperexpanded_clones_timepoint.pdf"),
    width = 24, height = 5, onefile = TRUE)

# UMAP coordinates derived from tnk object
p <- ggplot(tcrs_plot, aes(x = harmonytnkumap_1,
                           y = harmonytnkumap_2,
                           color = bin,
                           size = bin)) +
  geom_point_rast(raster.dpi = 300) +   # rasterized points
  scale_color_brewer(palette = "GnBu", name = "Clone abundance") +
  scale_size_manual(
    name = "Clone Abundance",
    values = c(
      "Rare (0.00001 < X <= 0.0001)" = 0.4,
      "Small (0.0001 < X <= 0.001)" = 0.5,
      "Medium (0.001 < X <= 0.01)" = 0.6,
      "Large (0.01 < X <= 0.1)" = 0.8,
      "Hyperexpanded (0.1 < X <= 1)" = 1
    )
  ) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "") +
  theme_minimal() +
  facet_wrap(~timepoint, ncol = 5) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14)
  )

print(p)
dev.off()


# 5. Reload the CD4/CD8 metric tables and tag CAR status # -------

metrics_df_CD4_CD8_noncar <- read_csv(file.path(tcr_dir, "diversity_metrics_noncar_CD4_CD8_Tcells.csv"))
metrics_df_CD4_CD8_car    <- read_csv(file.path(tcr_dir, "diversity_metrics_CAR_CD4_CD8_Tcells.csv"))

metrics_df_CAR_CD4_CD8$CAR <- "CAR+"
metrics_df_noncar_CD4_CD8$CAR <- "CAR-"

df_metrics_car <- rbind(metrics_df_CAR_CD4_CD8, metrics_df_noncar_CD4_CD8)

df_metrics <- df_metrics_car


# 6. PAIRED differences in CAR at each timepoint, split by celltype and response # -------

plots_list <- list()
for (t in unique(df_metrics$timepoint)){
  
  if (t %in% c("Apheresis", "Day 0")){
    next
  }
  print(t)
  df_peak <- df_metrics %>%
    filter(timepoint == t)
  
  df_ok <- df_peak %>%
    group_by(celltype, CAR, response) %>%
    mutate(n_per_group = n()) %>%
    ungroup() %>%
    group_by(celltype, response) %>%
    filter(
      n_distinct(CAR) == 2,          # both CAR levels present
      all(n_per_group >= 2)          # each has at least 2 obs
    ) %>%
    ungroup() %>%
    select(-n_per_group)
  

  df_wide <- df_ok %>%
    filter(CAR %in% c("CAR+","CAR-")) %>%
    select(patient_id, celltype, CAR, norm_entropy, response) %>%
    pivot_wider(names_from = CAR, values_from = norm_entropy) %>%
    filter(!is.na(`CAR+`) & !is.na(`CAR-`))   # keep only paired patients
  
  df_long <- df_wide %>%
    pivot_longer(cols = c(`CAR+`, `CAR-`),
                 names_to = "CAR",
                 values_to = "norm_entropy")
  

  stat.test <- df_long %>%
    group_by(celltype, response) %>%
    wilcox_test(norm_entropy ~ CAR, paired = TRUE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
  
    left_join(
      df_long %>%
        group_by(celltype, response) %>%
        summarise(y.position = 1 + 0.02, .groups = "drop"),
      by = c("response", "celltype")
    )
  
  p2 <- ggplot(df_peak, aes(x = CAR, y = norm_entropy, fill = CAR)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = CAR),
                width = 0.15,
                size = 2,
                alpha = 0.8) +
    scale_fill_manual(values = car_colors) +
    scale_color_manual(values = car_colors) +
    stat_pvalue_manual(stat.test, label = "p.adj") +
    facet_grid(~ celltype + response) +
    theme_minimal() +
    theme(legend.position = "right") +
    ggtitle(paste0(t, " CAR vs nonCAR in CD4, CD8, and response compartments p-val adjusted"))
  
  plots_list[[length(plots_list)+1]] <- p2
  
  write.csv(as.data.frame(stat.test),
            file.path(stats_dir, paste0("paired_wilcox_test_", t, "_CAR_vs_nonCAR_celltype_and_response.csv")))
}

pdf(file.path(tcr_dir, "paired_boxplot_wilcox_CAR_vs_nonCAR_celltype_response_split.pdf"), width=8, height=6)
for (plot in plots_list){
  print(plot)
}
dev.off()


# 7. PAIRED differences in CD4 vs CD8 at each timepoint, split by response and CAR # ------

df_metrics_car <- rbind(metrics_df_CAR_CD4_CD8, metrics_df_noncar_CD4_CD8)
df_metrics <- df_metrics_car

plots_list <- list()
for (t in unique(df_metrics$timepoint)){
  if( t == "Week 4"){
    next
  }
  print(t)
  df_peak <- df_metrics %>%
    filter(timepoint == t)
  
  df_ok <- df_peak %>%
    group_by(celltype, CAR, response) %>%
    mutate(n_per_group = n()) %>%
    ungroup() %>%
    group_by(CAR, response) %>%
    filter(
      n_distinct(celltype) == 2,     # both celltypes present
      all(n_per_group >= 2)          # each has at least 2 obs
    ) %>%
    ungroup() %>%
    select(-n_per_group)
  

  df_wide <- df_ok %>%
    filter(celltype %in% c("CD4 T","CD8 T")) %>%
    select(patient_id, celltype, CAR, norm_entropy, response) %>%
    pivot_wider(names_from = celltype, values_from = norm_entropy) %>%
    filter(!is.na(`CD4 T`) & !is.na(`CD8 T`))   # keep only paired patients
  
  df_long <- df_wide %>%
    pivot_longer(cols = c(`CD4 T`, `CD8 T`),
                 names_to = "celltype",
                 values_to = "norm_entropy")
  

  stat.test <- df_long %>%
    group_by(CAR, response) %>%
    wilcox_test(norm_entropy ~ celltype, paired = TRUE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
  
    left_join(
      df_long %>%
        group_by(CAR, response) %>%
        summarise(y.position = 1 + 0.02, .groups = "drop"),
      by = c("response", "CAR")
    )
  
  p2 <- ggplot(df_peak, aes(x = celltype, y = norm_entropy, fill = celltype)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = celltype),
                width = 0.15,
                size = 2,
                alpha = 0.8) +
    scale_fill_manual(values = celltype_colors) +
    scale_color_manual(values = celltype_colors) +
    stat_pvalue_manual(stat.test, label = "p") +
    facet_grid(~ CAR + response) +
    theme_minimal() +
    theme(legend.position = "right") +
    ggtitle(paste0(t, " CD4 vs CD8 in CAR and response compartments p-val non-adjusted"))
  
  plots_list[[length(plots_list)+1]] <- p2
  
  write.csv(as.data.frame(stat.test),
            file.path(stats_dir, paste0("pairwise_wilcox_test_", t, "_CD4_vs_CD8_CAR_and_response.csv")))
}

pdf(file.path(tcr_dir, "pairwise_boxplot_wilcox_CD4_vs_CD8_CAR_response_split.pdf"), width=8, height=6)
for (plot in plots_list){
  print(plot)
}
dev.off()


# 8. PAIRED change between Product and Peak, by CAR / celltype / response # ------

df_metrics <- rbind(metrics_CD4_CD8_car, metrics_CD4_CD8_noncar)
df_metrics$timepoint <- factor(df_metrics$timepoint, levels = timepoint_levels)

stat.test <- df_metrics %>%
  filter(timepoint %in% c("Product", "Peak")) %>%
  group_by(CAR, celltype, patient_id, response) %>%
  filter(n_distinct(timepoint) == 2) %>%
  ungroup() %>%
  group_by(CAR, celltype, response) %>%
  wilcox_test(norm_entropy ~ timepoint, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "timepoint")

stat.test <- stat.test %>% dplyr::select(-groups)
write.csv(stat.test,
          file.path(stats_dir, "paired_wilcoxon_test_norm_entropy_product_vs_peak_CAR_celltype_response_split.csv"))

pdf(file.path(tcr_dir, "paired_wilcoxon_test_norm_entropy_product_vs_peak_CAR_celltype_response_split.pdf"),
    width=18, height=6)
plot_df <- df_metrics %>%
  filter(timepoint %in% c("Product","Peak"))

print(
  ggplot(plot_df,
         aes(x = timepoint, y = norm_entropy)) +
    geom_boxplot(aes(fill = timepoint),
                 outlier.shape = NA,
                 alpha = 0.7) +
    geom_line(aes(group = patient_id),
              color = "grey70",
              alpha = 0.5) +
    geom_point(aes(color = timepoint),
               size = 2) +
    stat_pvalue_manual(
      stat.test,
      label = "p.adj",
      tip.length = 0.005
    ) + facet_grid(~CAR+celltype+response) + scale_color_manual(values = c(
      "Product" = "gold",
      "Peak" = "darkolivegreen3"
    )) +
    scale_fill_manual(values = c(
      "Product" = "gold",
      "Peak" = "darkolivegreen3"
    )) + theme_classic() +
    labs(
      x = NULL,
      y = "normalized entropy"
    ) +
    theme(
      legend.position = "right"
    ) + ggtitle("paired change in entropy p-value adjusted")
)

print(
  ggplot(plot_df,
         aes(x = timepoint, y = norm_entropy)) +
    geom_boxplot(aes(fill = timepoint),
                 outlier.shape = NA,
                 alpha = 0.7) +
    geom_line(aes(group = patient_id),
              color = "grey70",
              alpha = 0.5) +
    geom_point(aes(color = timepoint),
               size = 2) +
    stat_pvalue_manual(
      stat.test,
      label = "p"
    ) + facet_grid(~ CAR+celltype+response) + scale_color_manual(values = c(
      "Product" = "gold",
      "Peak" = "darkolivegreen3"
    )) +
    scale_fill_manual(values = c(
      "Product" = "gold",
      "Peak" = "darkolivegreen3"
    )) + theme_classic() +
    labs(
      x = NULL,
      y = "normalized entropy"
    ) +
    theme(
      legend.position = "right"
    ) + ggtitle("paired change in entropy p-value unadjusted")
)
dev.off()


# 9. Proportion of large and hyperexpanded clones # ------

df_clone <-  tcrs_metadata %>%
  left_join(df_clonal_homeostasis, by = c("patient_id", "response", "CTaa")) %>%
  filter(has_tcr_gex)


df_clone <- df_clone %>% group_by(response, patient_id, timepoint, tcell_label, CAR_value, bin) %>%
  summarise(n_bin = n(), .groups = "drop") %>%
  group_by(response, patient_id, timepoint, CAR_value, tcell_label) %>%
  mutate(prop = n_bin / sum(n_bin)) %>%
  ungroup()

# ---- Proportion of large and hyperexpanded clones ---------------------------
df_ids <- df_clone %>%
  distinct(patient_id, response, timepoint, tcell_label, CAR_value)

df_expanded <- df_clone %>%
  filter(bin %in% c(
    "Large (0.01 < X <= 0.1)",
    "Hyperexpanded (0.1 < X <= 1)"
  )) %>%
  group_by(patient_id, response, timepoint, tcell_label, CAR_value) %>%
  summarise(
    prop_expanded = sum(prop, na.rm = TRUE),
    .groups = "drop"
  )

df_expanded_full <- df_ids %>%
  left_join(df_expanded,
            by = c("patient_id","response","timepoint","tcell_label", "CAR_value")) %>%
  mutate(prop_expanded = ifelse(is.na(prop_expanded), 0, prop_expanded))

# ---- paired test within response, CAR and celltype --------------------------
stat.test <- df_expanded_full %>%
  filter(timepoint %in% c("Product", "Peak"),
         tcell_label %in% c("CD4 T","CD8 T")) %>%
  group_by(tcell_label, CAR_value, response, patient_id) %>%
  filter(n_distinct(timepoint) == 2) %>%
  ungroup() %>%
  group_by(tcell_label, CAR_value, response) %>%
  wilcox_test(prop_expanded ~ timepoint, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "timepoint")

stat.test <- stat.test %>% dplyr::select(-groups)
write.csv(stat.test,
          file.path(stats_dir, "paired_wilcoxon_test_expanded_clones_proportion_product_vs_peak_CAR_celltype_response_split.csv"))

pdf(file.path(tcr_dir, "paired_wilcoxon_unadj_proportion_expanded_clones_product_vs_peak_CAR_Celltype_Response_split.pdf"),
    width=18, height=6)
plot_df <- df_expanded_full %>%
  filter(timepoint %in% c("Product","Peak"),
         tcell_label %in% c("CD4 T","CD8 T"))

print(
  ggplot(plot_df,
         aes(x = timepoint, y = prop_expanded)) +
    
    geom_boxplot(aes(fill = timepoint),
                 outlier.shape = NA,
                 alpha = 0.7) +
    
    geom_line(aes(group = patient_id),
              color = "grey70",
              alpha = 0.5) +
    
    geom_point(aes(color = response),
               size = 2) +
    
    stat_pvalue_manual(
      stat.test,
      label = "p",
      tip.length = 0.005
    ) +
    
    facet_grid(~response + tcell_label + CAR_value) +
    
    scale_color_manual(values = response_colors) +
    
    scale_fill_manual(values = c(
      "Product" = "gold",
      "Peak" = "darkolivegreen3"
    )) +
    
    theme_classic() +
    
    labs(
      x = NULL,
      y = "Proportion Expanded Clones"
    ) +
    
    theme(
      legend.position = "right"
    )
)
dev.off()


# 10. CR vs PD in CAR+ and CAR- CD4 T and CD8 T # -----

df_metrics <- rbind(metrics_CD4_CD8_noncar, metrics_CD4_CD8_car)

wilcox_results <- df_metrics %>%
  group_by(timepoint, celltype, CAR) %>%
  filter(n_distinct(response) == 2) %>%
  summarise(
    n_total = n(),
    n_CR = sum(response == "CR"),
    n_PD = sum(response == "PD"),
    p_value = wilcox.test(norm_entropy ~ response)$p.value,
    .groups = "drop"
  ) %>%
  group_by(timepoint) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

print(wilcox_results)

wilcox_plot_df <- wilcox_results %>%
  mutate(
    group1 = "CR",
    group2 = "PD",
    label = paste0("p = ", signif(p_adj, 2))
  )

y_positions <- df_metrics %>%
  group_by(timepoint, celltype, CAR) %>%
  summarise(
    y.position = max(norm_entropy, na.rm = TRUE) * 1.1,
    .groups = "drop"
  )

wilcox_plot_df <- wilcox_plot_df %>%
  left_join(y_positions, by = c("timepoint", "celltype", "CAR"))

p <- ggplot(df_metrics, aes(x = response, y = norm_entropy, fill = response)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = response),
              width = 0,
              size = 2,
              alpha = 0.8) +
  facet_grid(~timepoint+celltype+CAR, scales = "free_y") +
  labs(
    x = "Response Group",
    y = "Normalized Shannon Entropy",
    title = "Entropy differences across response groups"
  ) +
  scale_fill_manual(values = response_colors) +
  scale_color_manual(values = response_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )

p <- p +
  stat_pvalue_manual(
    wilcox_plot_df,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    inherit.aes = FALSE
  )

write.csv(wilcox_results, file.path(stats_dir, "unpaired_wilcox_CR_vs_PD_split_CAR_celltype.csv"))

pdf(file.path(boxplot_dir, "boxplot_unpaired_wilcox_CR_vs_PD_split_CAR_celltype.pdf"), width=14, height=6)
print(p)
dev.off()