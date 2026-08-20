################################################################################
## Tensor-cell2cell factor loadings: response x timepoint boxplots + Wilcoxon
################################################################################

# ---- Libraries --------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(readxl)
library(rstatix)
library(ggpubr)

# ---- CONFIG -----------------------------------------------------------------
analysis_dir <- "."
setwd(analysis_dir)

objects_dir <- "seurat_objects"
tensor_dir  <- file.path("python_data", "images", "Tensor", "FullDataset_bysample")
fig_dir     <- file.path("images", "tensor_factors")
stats_dir   <- file.path(fig_dir, "wilcox_test")

for (d in c(fig_dir, stats_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

f_full_dataset <- "LBCL_full_dataset_clean_new.RDS"

set.seed(2024)


# 1. Load factor loadings # -------------

Loadings <- read_excel(file.path(tensor_dir, "Loadings.xlsx")) # loadings results from tensor factorization

colnames(Loadings) <- c("sample.name", "Factor_1", "Factor_2", "Factor_3",
                        "Factor_4", "Factor_5", "Factor_6")

# sample.name is "<patient.id>_<timepoint>" (sample_timepoint2 from the Python side)
Loadings <- Loadings %>%
  tidyr::separate(sample.name, into = c("patient_id", "timepoint"),
                  sep = "_", extra = "merge", fill = "right")


# 2. attach metadata # -----------------

full.dataset <- readRDS(file.path(objects_dir, f_full_dataset))

general.meta <- full.dataset@meta.data
df2 <- general.meta %>%
  arrange(patient.id, response) %>%
  distinct(patient.id, response)

Loadings <- Loadings %>%
  left_join(df2 %>% select(patient.id, response), by = c("patient_id" = "patient.id"))

Loadings$response_timepoint <- paste0(Loadings$response, "|", Loadings$timepoint)


# 3. Pivot long and run Wilcoxon tests # ----------------

df_long <- Loadings %>%
  tidyr::pivot_longer(
    cols = starts_with("Factor_"),
    names_to = "factor",
    values_to = "value"
  )

# ---- unpaired: CR vs PD within (and across) timepoints ----------------------
comparisons_response <- list(
  c("CR|APH", "PD|APH"),
  c("CR|TDN-D0", "PD|TDN-D0"),
  c("CR|Peak", "PD|Peak"),
  c("CR|4W", "PD|4W"),
  c("PD|TDN-D0", "CR|Peak"),
  c("PD|APH", "CR|Peak")
)


comparisons_time <- list(
  c("CR|APH", "CR|TDN-D0"),
  c("CR|APH", "CR|Peak"),
  c("CR|APH", "CR|4W"),
  c("CR|TDN-D0", "CR|Peak"),
  c("CR|TDN-D0", "CR|4W"),
  c("CR|Peak", "CR|4W"),
  c("PD|APH", "PD|TDN-D0"),
  c("PD|APH", "PD|Peak"),
  c("PD|APH", "PD|4W"),
  c("PD|TDN-D0", "PD|Peak"),
  c("PD|TDN-D0", "PD|4W"),
  c("PD|Peak", "PD|4W")
)



# ---- all comparisons, unpaired ---------
comparisons <- c(comparisons_response, comparisons_time)

stats <- df_long %>%
  group_by(factor) %>%
  pairwise_wilcox_test(
    value ~ response_timepoint,
    p.adjust.method = "BH",
    comparisons = comparisons
  ) %>%
  add_y_position()


# 4. Boxplots by response x timepoint, faceted by factor # -----

response_timepoint_levels <- c("CR|APH", "PD|APH",
                               "CR|TDN-D0", "PD|TDN-D0",
                               "CR|Peak", "PD|Peak",
                               "CR|4W", "PD|4W")

df_long$response_timepoint <- factor(df_long$response_timepoint,
                                     levels = response_timepoint_levels)

response_colors <- c("#4FB7C5", "#E7A75E", "#4FB7C5", "#E7A75E",
                     "#4FB7C5", "#E7A75E", "#4FB7C5", "#E7A75E")
names(response_colors) <- response_timepoint_levels

p <- ggplot(df_long,
            aes(x = response_timepoint,
                y = value,
                fill = response_timepoint)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = response_colors) +
  geom_jitter(aes(fill = response_timepoint), shape = 21,
              width = 0, alpha = 0.5) +
  facet_wrap(~factor, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(
    stats %>% dplyr::filter(p.adj <= 0.05),
    label = "p.adj",
    tip.length = 0.01
  )

pdf(file.path(fig_dir, "factor_boxplots_padj.pdf"), width = 12, height = 10)
print(p)
dev.off()


# 5. Save statistics # ------

stats$groups <- NULL

write.csv(stats, file.path(stats_dir, "unpaired_factor_wilcox_test.csv"))

