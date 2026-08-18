################################################################################
## Cell-type proportion analyses: responder (CR) vs non-responder (PD)
##
## Boxplots + Wilcoxon tests across compartments (CAR+ T/NK, CAR- T/NK, all
## CAR- cells, myeloid), timepoints, and CAR expression status.
##
## Assumes the full dataset .RDS already carries all harmonized metadata,
## including the canonical cell-type column `cell.anno`. 
##
## PATHS / IDENTIFIERS: everything environment- or study-specific lives in the
## CONFIG block. 
################################################################################

# --- Libraries ---------------------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(ggpubr)

# =============================================================================
# CONFIG
# =============================================================================
project_dir <- "."

obj_dir    <- file.path(project_dir, "seurat_objects")
files_dir  <- file.path(project_dir, "files")
output_dir <- file.path(project_dir, "images")

wilcox_dir  <- file.path(output_dir, "files", "wilcoxon_tables")
prop_dir    <- file.path(output_dir, "files", "proportion_tables")
boxplot_dir <- file.path(output_dir, "celltype_boxplots")

for (d in c(obj_dir, wilcox_dir, prop_dir, boxplot_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Object filenames (shared with the setup/DEG script)
f_full_dataset <- "LBCL_full_dataset_clean_new.RDS"
f_tnk_noncar   <- "nonCAR_TNK_filtered.RDS"
f_tnk_car      <- "CAROnly_TNK_filtered.RDS"
f_tnk_all      <- "all_samples_TNK_filtered.RDS"
f_mono         <- "all_samples_Mono_filtered.RDS"

# Cell-type definitions -- all keyed off `cell.anno`
tnk_celltypes <- c("CD4 EM-like", "CD4 Naive-like", "CD8 EM", "CD8 Naive-like",
                   "gdT", "MAIT", "NK", "Proliferating", "Treg")
myeloid_celltypes <- c("CD14 Mono", "CD16 Mono", "Dendritic Cell")

# Timepoints (values of `timepoint_order`)
timepoint_levels <- c("T1-APH", "T2-TDN", "T3-D0", "T4-Peak", "T5-4-Week")
tdn_label        <- "T2-TDN"   # pre-infusion timepoint, dropped from some figures

# Metadata columns pulled from every object
meta_cols <- c("patient.id", "sample.name", "response", "timepoint_order", "cell.anno")

# Plot colours
response_colours <- c("CR" = "#4FB7C5", "PD" = "#E7A75E")
carexp_colours   <- c("CAR" = "steelblue2", "nonCAR" = "gray")

# --- Sample-merge map --------------------------------------------------------
sample_merge_map <- list(
       "PT5-TDN" = c("PT5-TDN-neg",  "PT5-TDN-pos"),
       "PT15-TDN" = c("PT15-TDN-neg",  "PT15-TDN-pos"),
       "PT9-TDN" = c("PT9_TDN_CARN", "PT9_C_TDN_CARP")
   )

# =============================================================================
# LOAD AND SUBSET OBJECTS
# =============================================================================
full.dataset <- readRDS(file.path(obj_dir, f_full_dataset))

tnk.seurat  <- subset(full.dataset, subset = cell.anno %in% tnk_celltypes)
tnk.car     <- subset(full.dataset, subset = CAR.exp == TRUE  & cell.anno %in% tnk_celltypes)
tnk.noncar  <- subset(full.dataset, subset = CAR.exp == FALSE & cell.anno %in% tnk_celltypes)
mono.object <- subset(full.dataset, subset = cell.anno %in% myeloid_celltypes)

# All CAR-negative cells (not just T/NK). Derived by set difference rather than
# `CAR.exp == FALSE` so that cells where CAR.exp is NA (e.g. myeloid) are kept.
all.cells   <- Cells(full.dataset)
car.only    <- Cells(tnk.car)
noncar.only <- setdiff(all.cells, car.only)
noncar.obj  <- subset(full.dataset, cells = noncar.only)

saveRDS(tnk.seurat,  file.path(obj_dir, f_tnk_all))
saveRDS(tnk.car,     file.path(obj_dir, f_tnk_car))
saveRDS(tnk.noncar,  file.path(obj_dir, f_tnk_noncar))
saveRDS(mono.object, file.path(obj_dir, f_mono))

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Rejoin split libraries into a single sample.name2 column using the map above.
apply_sample_merge <- function(meta, merge_map) {
  meta$sample.name2 <- meta$sample.name
  for (target in names(merge_map)) {
    meta$sample.name2[meta$sample.name %in% merge_map[[target]]] <- target
  }
  meta
}

# Pull the standard metadata columns and apply the sample merge in one step.
get_prop_metadata <- function(obj, merge_map = sample_merge_map) {
  apply_sample_merge(obj@meta.data[, meta_cols], merge_map)
}

# Column-wise proportions from a celltype x sample table, joined to metadata.
get_prop_table <- function(ptable, big.df) {
  ptable  <- prop.table(ptable, margin = 2)
  prop.df <- as.data.frame(ptable)
  colnames(prop.df) <- c("celltype", "sample.name2", "prop")
  prop.df$celltype  <- as.character(prop.df$celltype)
  
  prop.df <- left_join(
    prop.df,
    big.df[, c("patient.id", "sample.name2", "response", "timepoint_order")],
    by = "sample.name2"
  )
  prop.df %>% distinct()
}

# Raw counts from a celltype x sample table, joined to metadata.
get_counts_table <- function(ptable, big.df) {
  prop.df <- as.data.frame(ptable)
  colnames(prop.df) <- c("celltype", "sample.name2", "count")
  prop.df$celltype  <- as.character(prop.df$celltype)
  
  prop.df <- left_join(
    prop.df,
    big.df[, c("patient.id", "sample.name2", "response", "timepoint_order")],
    by = "sample.name2"
  )
  prop.df %>% distinct()
}

# Build the proportion table for one compartment straight from its metadata.
build_prop_table <- function(meta) {
  get_prop_table(table(meta$cell.anno, meta$sample.name2), meta)
}

# Unpaired two-group Wilcoxon, returns a tidy data frame.
wilcox_2group <- function(d, group_var, value_var = "prop") {
  wilcox.test(
    as.formula(paste(value_var, "~", group_var)),
    data  = d,
    exact = FALSE
  ) %>% broom::tidy()
}

# Paired Wilcoxon across two levels of `group_var`, matched on `id_var`.
wilcox_2group_paired <- function(d,
                                 id_var    = "patient.id",
                                 group_var = "CAR.exp",
                                 value_var = "prop") {
  
  wide <- d %>%
    dplyr::select(dplyr::all_of(c(id_var, group_var, value_var))) %>%
    tidyr::pivot_wider(names_from  = dplyr::all_of(group_var),
                       values_from = dplyr::all_of(value_var),
                       values_fn   = mean) %>%   # guard against duplicate rows
    tidyr::drop_na()
  
  stopifnot(ncol(wide) == 3)   # id column + the two groups
  
  wilcox.test(wide[[2]], wide[[3]],
              paired = TRUE, exact = FALSE) %>%
    broom::tidy()
}

# CR vs PD Wilcoxon per (timepoint x celltype), BH-corrected within timepoint.
run_resp_vs_non <- function(prop_df) {
  prop_df %>%
    group_by(timepoint_order, celltype) %>%
    filter(n_distinct(response) == 2) %>%      # need both CR and PD
    group_modify(~ wilcox_2group(.x, "response")) %>%
    ungroup() %>%
    group_by(timepoint_order) %>%
    mutate(p_adj_BH_time = p.adjust(p.value, method = "BH")) %>%
    ungroup()
}

# Standard CR vs PD proportion boxplot, faceted by cell type.
plot_prop_boxplot <- function(props_df, title, outfile, height = 15, width = 18) {
  p <- ggplot(props_df, aes(x = timepoint_order, y = prop, fill = response)) +
    geom_boxplot(alpha = 0.7, width = 0.6) +
    geom_jitter(position = position_dodge(width = 0.6),
                shape = 21, size = 2, stroke = 0.2) +
    stat_compare_means(method = "wilcox.test", label = "p.format", hide.ns = TRUE) +
    facet_wrap(~ celltype, scales = "free_y") +
    labs(y = "Proportion of Cells", x = "Timepoint", title = title) +
    scale_fill_manual(values = response_colours) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf(outfile, height = height, width = width)
  print(p)
  dev.off()
  
  invisible(p)
}

# Single-timepoint boxplot with one Wilcoxon bracket per cell type.
plot_tdn_boxplot <- function(props_df,
                             title,
                             outfile,
                             timepoint   = tdn_label,
                             dodge_width = 0.8) {
  
  df_tdn <- props_df %>%
    filter(timepoint_order == timepoint) %>%
    mutate(
      response = factor(response, levels = names(response_colours)),
      celltype = factor(celltype)
    )
  
  # Per-celltype Wilcoxon test (need both groups present)
  wilcox_df <- df_tdn %>%
    group_by(celltype) %>%
    filter(n_distinct(response) == 2) %>%
    summarise(
      p = wilcox.test(prop ~ response, exact = FALSE)$p.value,
      .groups = "drop"
    ) %>%
    mutate(p.adj = p.adjust(p, method = "BH"))
  
  # Bracket coordinates so the p-value spans the two dodged boxes per celltype
  cell_levels <- levels(df_tdn$celltype)
  y_top       <- max(df_tdn$prop, na.rm = TRUE) * 1.05
  wilcox_plot <- wilcox_df %>%
    mutate(
      x_centre   = as.numeric(factor(celltype, levels = cell_levels)),
      xmin       = x_centre - dodge_width / 2 + 0.05,
      xmax       = x_centre + dodge_width / 2 - 0.05,
      y.position = y_top,
      label      = paste0("p=", signif(p, 3)),
      group1     = "CR",
      group2     = "PD"
    )
  
  p <- ggplot(df_tdn, aes(x = celltype, y = prop, fill = response)) +
    geom_boxplot(alpha = 0.7, width = 0.6,
                 position = position_dodge(width = dodge_width)) +
    geom_jitter(position = position_dodge(width = dodge_width),
                shape = 21, size = 2, stroke = 0.2) +
    scale_fill_manual(values = response_colours) +
    stat_pvalue_manual(
      wilcox_plot,
      label      = "label",
      xmin       = "xmin",
      xmax       = "xmax",
      y.position = "y.position",
      tip.length = 0.02
    ) +
    theme_minimal() +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    labs(title = title, x = "Cell Type", y = "Proportion")
  
  ggsave(outfile, p, width = 10, height = 6)
  
  invisible(list(plot = p, stats = wilcox_df))
}

# All pairwise timepoint comparisons within each (response x celltype), paired
# on patient.id. BH-adjusted within each (response, celltype) family.
run_paired_timepoint_wilcox <- function(props_df, min_pairs = 2) {
  
  results <- list()
  
  for (resp in unique(props_df$response)) {
    for (ct in unique(props_df$celltype)) {
      
      sub <- props_df %>%
        filter(response == resp, celltype == ct) %>%
        select(patient.id, timepoint_order, prop)
      
      tps <- unique(as.character(sub$timepoint_order))
      if (length(tps) < 2) next
      combs <- combn(tps, 2, simplify = FALSE)
      
      for (cmb in combs) {
        tp1 <- cmb[1]
        tp2 <- cmb[2]
        
        # Wide-format so each patient contributes one paired observation
        pair_data <- sub %>%
          filter(timepoint_order %in% c(tp1, tp2)) %>%
          pivot_wider(names_from  = timepoint_order,
                      values_from = prop,
                      values_fn   = mean) %>%
          drop_na(all_of(c(tp1, tp2)))
        
        if (nrow(pair_data) < min_pairs) next
        
        test <- wilcox.test(pair_data[[tp1]], pair_data[[tp2]],
                            paired = TRUE, exact = FALSE)
        
        results[[length(results) + 1]] <- data.frame(
          response   = resp,
          celltype   = ct,
          timepoint1 = tp1,
          timepoint2 = tp2,
          n_pairs    = nrow(pair_data),
          p_value    = test$p.value
        )
      }
    }
  }
  
  if (length(results) == 0) return(data.frame())
  
  bind_rows(results) %>%
    group_by(response, celltype) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}

# =============================================================================
# METADATA EXTRACTION AND PROPORTION TABLES
# =============================================================================
metadata.car        <- get_prop_metadata(tnk.car)      # CAR+ T/NK
metadata.noncar.tnk <- get_prop_metadata(tnk.noncar)   # CAR- T/NK
metadata.noncar.all <- get_prop_metadata(noncar.obj)   # all CAR- cells
metadata.mono       <- get_prop_metadata(mono.object)  # myeloid
metadata.tnk        <- get_prop_metadata(tnk.seurat)   # all T/NK (CAR +/-)

# Proportions are computed WITHIN each compartment, i.e. each column of the
# table sums to 1 across that compartment's cell types.
celltype_props_car        <- build_prop_table(metadata.car)
celltype_props_noncar_tnk <- build_prop_table(metadata.noncar.tnk)
celltype_props_noncar_all <- build_prop_table(metadata.noncar.all)
celltype_props_mono       <- build_prop_table(metadata.mono)

# =============================================================================
# 1. CR vs PD BOXPLOTS + WILCOXON, PER COMPARTMENT
# =============================================================================
compartments <- list(
  CAR = list(
    props  = celltype_props_car,
    title  = "CAR+ T/NK cell type proportions",
    pdf    = "CAR_CelltypeProp.pdf",
    csv    = "CAR_celltype_props_wilcoxon.csv",
    height = 15, width = 18
  ),
  NonCAR_All = list(
    props  = celltype_props_noncar_all,
    title  = "All CAR- cell type proportions",
    pdf    = "NonCAR_All_CelltypeProp.pdf",
    csv    = "NonCAR_All_celltype_props_wilcoxon.csv",
    height = 15, width = 22
  ),
  NonCAR_TNK = list(
    props  = celltype_props_noncar_tnk,
    title  = "CAR- T/NK cell type proportions",
    pdf    = "NonCAR_TNK_CelltypeProp.pdf",
    csv    = "NonCAR_TNK_celltype_props_wilcoxon.csv",
    height = 14, width = 16
  ),
  Myeloid = list(
    props  = celltype_props_mono,
    title  = "Myeloid cell type proportions",
    pdf    = "Myeloid_CelltypeProp.pdf",
    csv    = "Myeloid_celltype_props_wilcoxon.csv",
    height = 10, width = 14
  )
)

for (nm in names(compartments)) {
  cmp <- compartments[[nm]]
  
  plot_prop_boxplot(cmp$props, cmp$title,
                    file.path(boxplot_dir, cmp$pdf),
                    height = cmp$height, width = cmp$width)
  
  write.csv(run_resp_vs_non(cmp$props),
            file.path(wilcox_dir, cmp$csv),
            row.names = FALSE)
  
  write.csv(cmp$props,
            file.path(prop_dir, paste0(nm, "_celltype_proportions.csv")),
            row.names = FALSE)
}

# =============================================================================
# 2. TDN-TIMEPOINT BOXPLOTS (CAR+ and CAR- T/NK)
# =============================================================================
car_tdn <- plot_tdn_boxplot(
  props_df = celltype_props_car,
  title    = "CAR+ cell type proportions at TDN",
  outfile  = file.path(boxplot_dir, "TDN_all_celltypes_boxplot_CRvsPD.pdf")
)

noncar_tdn <- plot_tdn_boxplot(
  props_df = celltype_props_noncar_tnk,
  title    = "CAR- cell type proportions at TDN",
  outfile  = file.path(boxplot_dir, "nonCAR_TDN_all_celltypes_boxplot_CRvsPD.pdf")
)

# =============================================================================
# 3. PAIRED WILCOXON ACROSS TIMEPOINTS, WITHIN RESPONSE GROUP
# =============================================================================
paired_inputs <- list(
  CARonly       = celltype_props_car,
  NONCAR_T_only = celltype_props_noncar_tnk,
  NONCAR_ALL    = celltype_props_noncar_all
)

for (nm in names(paired_inputs)) {
  paired_results <- run_paired_timepoint_wilcox(paired_inputs[[nm]])
  write.csv(paired_results,
            file.path(wilcox_dir, paste0(nm, "_timepoint_paired_wilcoxon_results.csv")),
            row.names = FALSE)
}

# =============================================================================
# 4. COMBINED FIGURE: CR vs PD plus within-response timepoint comparisons
# =============================================================================

# Drop the pre-infusion timepoint for this summary figure
plot_df <- celltype_props_noncar_all %>%
  filter(timepoint_order != tdn_label)

post_tdn <- setdiff(timepoint_levels, tdn_label)
timepoint_comparisons <- combn(post_tdn, 2, simplify = FALSE)   # all 6 pairs

max_prop <- max(plot_df$prop, na.rm = TRUE)
cr_label_y <- seq(0.65, 0.85, length.out = length(timepoint_comparisons)) * max_prop
pd_label_y <- seq(0.90, 1.10, length.out = length(timepoint_comparisons)) * max_prop

p_combined <- ggplot(plot_df, aes(x = timepoint_order, y = prop, fill = response)) +
  
  geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA,
               position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75),
              shape = 21, size = 2, stroke = 0.2) +
  
  # CR vs PD at each timepoint
  stat_compare_means(
    aes(group = response),
    method = "wilcox.test", label = "p.format", hide.ns = TRUE
  ) +
  
  # Within-CR timepoint comparisons
  stat_compare_means(
    data          = subset(plot_df, response == "CR"),
    comparisons   = timepoint_comparisons,
    method        = "wilcox.test",
    label         = "p.signif",
    hide.ns       = TRUE,
    color         = response_colours[["CR"]],
    step.increase = 0.05,
    label.y       = cr_label_y
  ) +
  
  # Within-PD timepoint comparisons (offset upward)
  stat_compare_means(
    data          = subset(plot_df, response == "PD"),
    comparisons   = timepoint_comparisons,
    method        = "wilcox.test",
    label         = "p.signif",
    hide.ns       = TRUE,
    color         = response_colours[["PD"]],
    step.increase = 0.05,
    label.y       = pd_label_y
  ) +
  
  facet_wrap(~ celltype, scales = "free_y") +
  labs(y = "Proportion of Cells", x = "Timepoint",
       title = "All CAR- cell type proportions by patient") +
  scale_fill_manual(values  = response_colours) +
  scale_color_manual(values = response_colours) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 90, hjust = 1),
    strip.text      = element_text(face = "bold"),
    legend.position = "top"
  )

pdf(file.path(boxplot_dir, "NonCAR_All_CelltypeProp_withTimepointStats.pdf"),
    height = 25, width = 20)
print(p_combined)
dev.off()

# =============================================================================
# 5. CAR vs NON-CAR T CELL COMPARISON
# =============================================================================
# Both tables must have: celltype, sample.name2, prop, patient.id, response,
# timepoint_order. Proportions are within-compartment, so a CAR vs nonCAR
# contrast compares each cell type's share of its own compartment.
celltype_props_tnk <- bind_rows(
  celltype_props_car        %>% mutate(CAR.exp = "CAR"),
  celltype_props_noncar_tnk %>% mutate(CAR.exp = "nonCAR")
) %>%
  mutate(
    car.time      = paste0(timepoint_order, "_", CAR.exp),
    response.time = paste0(timepoint_order, "_", response)
  )

# --- Boxplots ----------------------------------------------------------------
bxplot_response <- ggplot(celltype_props_tnk,
                          aes(x = car.time, y = prop, fill = response)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  geom_jitter(position = position_dodge(width = 0.6),
              shape = 21, size = 2, stroke = 0.2) +
  stat_compare_means(aes(group = response),
                     method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = response_colours) +
  facet_wrap(~ celltype, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf(file.path(boxplot_dir, "boxplot_car_noncar_props.pdf"), height = 12, width = 25)
print(bxplot_response)
dev.off()

bxplot_carexp <- ggplot(celltype_props_tnk,
                        aes(x = response.time, y = prop, fill = CAR.exp)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  geom_jitter(position = position_dodge(width = 0.6),
              shape = 21, size = 2, stroke = 0.2) +
  stat_compare_means(aes(group = CAR.exp),
                     method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = carexp_colours) +
  facet_wrap(~ celltype, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf(file.path(boxplot_dir, "boxplot_car_noncar_props_fillCAR.pdf"), height = 12, width = 25)
print(bxplot_carexp)
dev.off()


# --- Wilcoxon contrasts (per timepoint x celltype) ---------------------------

# (1) CR vs PD, within each CAR / nonCAR stratum (unpaired)
resp_vs_non_by_car_time <- celltype_props_tnk %>%
  group_by(timepoint_order, celltype, CAR.exp) %>%
  filter(n_distinct(response) == 2) %>%
  group_modify(~ wilcox_2group(.x, "response")) %>%
  mutate(contrast = "Resp_vs_NonResp_within_CARexp") %>%
  ungroup()

# (2) CAR vs nonCAR within responders, paired on patient.id
car_vs_noncar_resp_time_pair <- celltype_props_tnk %>%
  filter(response == "CR") %>%
  group_by(timepoint_order, celltype) %>%
  filter(n_distinct(CAR.exp)    == 2,
         n_distinct(patient.id) > 1) %>%
  group_modify(~ wilcox_2group_paired(.x)) %>%
  mutate(contrast = "CAR_vs_nonCAR_within_responder_paired") %>%
  ungroup()

# (3) CAR vs nonCAR within non-responders, paired on patient.id
car_vs_noncar_nonresp_time_pair <- celltype_props_tnk %>%
  filter(response == "PD") %>%
  group_by(timepoint_order, celltype) %>%
  filter(n_distinct(CAR.exp)    == 2,
         n_distinct(patient.id) > 1) %>%
  group_modify(~ wilcox_2group_paired(.x)) %>%
  mutate(contrast = "CAR_vs_nonCAR_within_nonresponder_paired") %>%
  ungroup()

# --- BH adjustment -----------------------------------------------------------
# First pass: within each contrast family, across celltypes per timepoint
# (and per CAR.exp stratum for contrast 1).
adjust_within_timepoint <- function(df, extra_group = NULL) {
  grp <- c("timepoint_order", extra_group)
  df %>%
    group_by(across(all_of(grp))) %>%
    mutate(p_adj_BH_acrosstest = p.adjust(p.value, method = "BH")) %>%
    ungroup()
}

resp_vs_non_by_car_time_adj         <- adjust_within_timepoint(resp_vs_non_by_car_time, "CAR.exp")
car_vs_noncar_resp_time_pair_adj    <- adjust_within_timepoint(car_vs_noncar_resp_time_pair)
car_vs_noncar_nonresp_time_pair_adj <- adjust_within_timepoint(car_vs_noncar_nonresp_time_pair)

# Paired tests, combined
all_tests_time_pair <- bind_rows(
  car_vs_noncar_resp_time_pair_adj,
  car_vs_noncar_nonresp_time_pair_adj
)

# --- Write outputs -----------------------------------------------------------
write.csv(resp_vs_non_by_car_time_adj,
          file.path(wilcox_dir, "Tcell_wilcoxon_CR_vs_PD_comparisons_CAR_and_nonCAR.csv"),
          row.names = FALSE)

write.csv(all_tests_time_pair,
          file.path(wilcox_dir, "Tcell_wilcoxon_CARvsnonCAR_paired_comparisons.csv"),
          row.names = FALSE)

write.csv(celltype_props_tnk,
          file.path(prop_dir, "Tcell_proportions_raw_per_sample_CAR_and_nonCAR.csv"),
          row.names = FALSE)