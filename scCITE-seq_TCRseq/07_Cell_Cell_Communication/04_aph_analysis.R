## =============================================================================
## APH ligand-receptor interactions (Non-CAR)
## Fisher-aggregated LIANA results + rank-aggregate DEA
##
##   Section A: CD8 TEM as SOURCE  -> myeloid as TARGET   (ligand side)
##   Section B: monocytes as SOURCE -> T cells as TARGET  (receptor side)
##
## The two sections run in opposite directions, so their cell-type vectors are
## named separately rather than reassigning source_cells / target_cells.
## =============================================================================

## ---- 0. setup ----------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(patchwork)
library(igraph)
library(ggraph)

set.seed(42)   # ggraph "fr" layout is stochastic

# All I/O is relative to the project root -- set these two and nothing below
# needs a path.
data_dir <- "data"
fig_dir  <- "figs"


# Shared diverging palette (low / mid / high)
response_colors <- c("#E7A75E", "white", "#4FB7C5")

# Cell-type labels
t_cell_types  <- c("CD4 Naive", "CD4 TEM-like", "CD8 Naive", "CD8 TEM")
myeloid_types <- c("CD14 Mono", "CD16 Mono", "Dendritic Cell")
mono_types    <- c("CD14 Mono", "CD16 Mono")

tp <- "Non-CAR_APH"   # Non-CAR_APH


## ---- 1. load -----------------------------------------------------------------

# Input filenames have been shortened Rename the files on
# disk to match, or point these three lines at whatever they are called.
fishers_results_aph <- read_csv(file.path(
  data_dir, paste0(tp, "_fisher_ligand.csv")
))

fishers_results_aph_r <- read_csv(file.path(
  data_dir, paste0(tp, "_fisher_receptor.csv")
))

car.dea.raa.aph <- read_csv(file.path(
  data_dir, paste0(tp, "_liana_dea.csv")
))


## =============================================================================
## SECTION A -- IFNG / CD8 TEM as source, myeloid as target
## =============================================================================

t_source     <- c("CD8 TEM")
myeloid_targ <- myeloid_types

## ----  fisher hits ---------------------------------------------------------

# Ligand side (source = CD8 TEM, targets string contains a myeloid type)
fishers_ligand_t_to_mono <- fishers_results_aph %>%
  filter(
    source %in% t_source,
    combined_padj < 0.05
  ) %>%
  rowwise() %>%
  mutate(
    matching_sources = paste(
      myeloid_targ[
        myeloid_targ %in% str_split(targets, ";\\s*")[[1]]
      ],
      collapse = ";"
    )
  ) %>%
  ungroup() %>%
  filter(matching_sources != "")

names(fishers_ligand_t_to_mono)   # sanity check on column names

unique_ligands   <- unique(fishers_ligand_t_to_mono$ligand_complex)
unique_receptors <- unique(fishers_ligand_t_to_mono$receptor_complex)

t_source_mono_target <- car.dea.raa.aph %>%
  filter(source %in% t_source, target %in% myeloid_targ) %>%
  filter(ligand_complex %in% unique_ligands, receptor_complex %in% unique_receptors)


## ----  dot plot ------------------------------------------------------------

plot_df_lig <- fishers_ligand_t_to_mono %>%
  mutate(
    interaction    = paste0(ligand_complex, " -> ", receptor_complex),
    neg_log10_padj = -log10(pmax(combined_padj, 1e-300))
  )

lvls <- plot_df_lig %>%
  arrange(desc(combined_padj), mean_interaction_stat) %>%
  pull(interaction)

plot_df_lig <- plot_df_lig %>% mutate(interaction = factor(interaction, levels = lvls))

lim <- max(abs(plot_df_lig$mean_interaction_stat), na.rm = TRUE)

p_dots <- ggplot(plot_df_lig,
                 aes(x = source, y = interaction,
                     colour = mean_interaction_stat, size = neg_log10_padj)) +
  geom_point() +
  scale_colour_gradient2(
    low = response_colors[1], mid = "grey90", high = response_colors[3],
    midpoint = 0, limits = c(-lim, lim)
  ) +
  scale_size_continuous(range = c(2, 8)) +
  scale_y_discrete(limits = lvls, drop = FALSE) +
  labs(
    x = NULL, y = "Interaction (ligand -> receptor)",
    colour = "Mean\ninteraction stat",
    size   = "-log10(combined padj)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(face = "bold"),
    legend.position  = "left"
  )


## ---- assemble ------------------------------------------------------------



# Tile grid variant instead of the text list
tile_df <- plot_df_lig %>%
  select(interaction, targets) %>%
  separate_rows(targets, sep = "\\s*;\\s*") %>%
  filter(!is.na(targets), targets != "") %>%
  distinct() %>%
  mutate(interaction = factor(interaction, levels = lvls))

p_tiles <- ggplot(tile_df, aes(x = targets, y = interaction)) +
  geom_tile(fill = "grey30", colour = "white", linewidth = 0.4) +
  scale_y_discrete(limits = lvls, drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_blank())

pdf(file.path(fig_dir, "cd8tem_myeloid_dotplot.pdf"),
    width = 10, height = 8)
print(p_dots + p_tiles + plot_layout(widths = c(1, 1.6)))
dev.off()


## =============================================================================
## SECTION B -- T cells as RECEPTOR at APH (monocyte source -> T cell target)
## =============================================================================

mono_source <- mono_types
t_target    <- c("CD8 TEM", "CD8 Naive", "CD4 Naive", "CD4 TEM-like")

mono_targ_t_source <- car.dea.raa.aph %>%
  filter(source %in% mono_source, target %in% t_target)


fishers_r_mono_to_tcell <- fishers_results_aph_r %>%
  filter(
    target %in% t_target,
    combined_padj < 0.05
  ) %>%
  rowwise() %>%
  mutate(
    matching_sources = paste(
      mono_source[
        mono_source %in% str_split(sources, ";\\s*")[[1]]
      ],
      collapse = ";"
    )
  ) %>%
  ungroup() %>%
  filter(matching_sources != "") %>%
  select(
    target,
    ligand,
    receptor,
    matching_sources,
    n_sources,
    mean_interaction_stat,
    combined_padj
  ) %>%
  arrange(combined_padj)


plot_df_rec <- fishers_r_mono_to_tcell %>%
  separate_rows(matching_sources, sep = ";") %>%
  mutate(
    matching_sources = trimws(matching_sources),
    interaction      = paste(ligand, receptor, sep = " -> "),
    neglog10_padj    = -log10(combined_padj)
  ) %>%
  left_join(
    mono_targ_t_source %>%
      select(
        target,
        source,
        ligand,
        receptor,
        interaction_stat
      ),
    by = c(
      "target",
      "matching_sources" = "source",
      "ligand",
      "receptor"
    )
  )

interaction_order <- plot_df_rec %>%
  group_by(interaction) %>%
  summarise(best_padj = min(combined_padj), .groups = "drop") %>%
  arrange(best_padj) %>%
  pull(interaction)

plot_df_rec$interaction <- factor(
  plot_df_rec$interaction,
  levels = rev(interaction_order)
)

top_int <- plot_df_rec %>%
  group_by(interaction) %>%
  summarise(best_padj = min(combined_padj), .groups = "drop") %>%
  arrange(best_padj) %>%
  slice_head(n = 25)


## ---- dot plots -----------------------------------------------------------


pdf(file.path(fig_dir, "dotplot_mono_to_tcell_facet.pdf"),
    width = 12, height = 10)
print(
  plot_df_rec %>%
    ggplot(aes(x = matching_sources,
               y = reorder(interaction, neglog10_padj))) +
    geom_point(aes(size = neglog10_padj, color = interaction_stat)) +
    facet_wrap(~target, scales = "free_y") +
    scale_color_gradient2(
      low      = response_colors[1],
      mid      = response_colors[2],
      high     = response_colors[3],
      midpoint = 0
    ) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8))
)
dev.off()


## ---- edge networks -------------------------------------------------------

edges <- plot_df_rec %>%
  group_by(target, matching_sources) %>%
  summarise(
    score            = max(neglog10_padj),
    interaction_stat = mean(interaction_stat, na.rm = TRUE),
    .groups          = "drop"
  )

edges_pos <- edges %>% filter(interaction_stat > 0)
edges_neg <- edges %>% filter(interaction_stat < 0)

g_pos <- graph_from_data_frame(edges_pos)
g_neg <- graph_from_data_frame(edges_neg)

pdf(file.path(fig_dir, "network_mono_to_tcell.pdf"),
    width = 6, height = 6)
print(
  ggraph(g_pos, layout = "fr") +
    geom_edge_link(aes(width = score, edge_color = interaction_stat), alpha = 0.7) +
    scale_edge_color_gradient2(
      low      = response_colors[1],
      mid      = response_colors[2],
      high     = response_colors[3],
      midpoint = 0
    ) +
    geom_node_point(size = 8) +
    geom_node_text(aes(label = name), repel = TRUE) +
    theme_void()
)
print(
  ggraph(g_neg, layout = "fr") +
    geom_edge_link(aes(width = score, edge_color = interaction_stat), alpha = 0.7) +
    scale_edge_color_gradient2(
      low      = response_colors[1],
      mid      = response_colors[2],
      high     = response_colors[3],
      midpoint = 0
    ) +
    geom_node_point(size = 8) +
    geom_node_text(aes(label = name), repel = TRUE) +
    theme_void()
)
dev.off()