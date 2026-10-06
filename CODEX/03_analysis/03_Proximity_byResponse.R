library(dplyr)
library(ggplot2)

BASE_DIR   <- "/path/to/your/"
SCRIPT_DIR <- "/path/to/your/Scripts"
source(file.path(SCRIPT_DIR, "Utilities", "CODEX_Functions.R"))

DISTANCE_DIR <- file.path(BASE_DIR, "Distance_Analysis", "Distance_Results")
PROX_DIR     <- file.path(DISTANCE_DIR, "results", "proximity_by_response")
FIG_DIR      <- file.path(BASE_DIR, "Results", "Figures")
dir.create(PROX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_SAMPLES <- 4   # keep source-target pairs observed in >= 4 samples

# Figure 5f: proximity of selected source cell types to CD8 TEM
FIG_TARGET  <- "CD8 TEM"
FIG_SOURCES <- c("CD14+CD16+TIMs", "Endothelial")

# ==============================================================================
# 1. Normalized Proximity Ranks
# ==============================================================================

prox_files <- list.files(DISTANCE_DIR, pattern = "^Proximity_.*\\.rds$", full.names = TRUE)

combined_df <- bind_rows(lapply(prox_files, function(f) {
  d <- readRDS(f)
  d <- d[d$celltype != d$Target_CT & !is.na(d$coefficients_2), ]
  d %>% group_by(sample) %>%
    mutate(rank = rank(coefficients_2, ties.method = "first"),
           normalized_rank = if (n() > 1) (rank(rank) - 1) / (n() - 1) else 0.5) %>%
    ungroup()
}))

keep_pairs <- combined_df %>% group_by(Target_CT, celltype) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
  filter(n_samples >= MIN_SAMPLES)
prox <- combined_df %>% semi_join(keep_pairs, by = c("Target_CT", "celltype"))
prox$Response <- factor(RESPONSE_MAP[prox$sample], levels = c("CR", "PD"))
write.csv(prox, file.path(PROX_DIR, "Proximity_Summary_L2.csv"), row.names = FALSE)

# ==============================================================================
# 2. CR vs PD Wilcoxon Rank-Sum Test per Source-Target Pair
# ==============================================================================

prox_stats <- prox %>% group_by(Target_CT, source = celltype) %>%
  summarise(n_CR = sum(Response == "CR"), n_PD = sum(Response == "PD"),
            normrank_median_CR = median(normalized_rank[Response == "CR"]),
            normrank_median_PD = median(normalized_rank[Response == "PD"]),
            normrank_pvalue = if (n_CR >= 2 && n_PD >= 2)
              suppressWarnings(wilcox.test(normalized_rank[Response == "CR"],
                                           normalized_rank[Response == "PD"])$p.value)
              else NA_real_,
            .groups = "drop") %>%
  group_by(Target_CT) %>%
  mutate(normrank_qvalue_within = p.adjust(normrank_pvalue, method = "BH")) %>%
  ungroup()
write.csv(prox_stats, file.path(PROX_DIR, "Proximity_Wilcox_L2.csv"), row.names = FALSE)

# ==============================================================================
# 3. Proximity Boxplots by Response
# ==============================================================================

plot_prox <- function(bdf, stats, title, ord = NULL) {
  if (is.null(ord)) ord <- bdf %>% group_by(celltype) %>%
    summarise(m = median(normalized_rank)) %>% arrange(m) %>% pull(celltype)
  bdf$celltype <- factor(bdf$celltype, levels = ord)
  lab <- stats %>% filter(source %in% ord) %>%
    transmute(celltype = factor(source, levels = ord),
              label = paste0("p=", signif(normrank_pvalue, 2)))
  ggplot(bdf, aes(x = celltype, y = normalized_rank, fill = Response)) +
    geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.85) +
    geom_point(position = position_dodge(0.8), size = 0.8, alpha = 0.6) +
    geom_text(data = lab, aes(x = celltype, y = 1.08, label = label),
              inherit.aes = FALSE, size = 3) +
    scale_fill_manual(values = RESPONSE_COLS) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Source cell type", y = "Normalized proximity rank (0 = closest)",
         title = title) +
    coord_cartesian(ylim = c(0, 1.12))
}

for (target in unique(prox$Target_CT)) {
  p <- plot_prox(prox[prox$Target_CT == target, ],
                 prox_stats[prox_stats$Target_CT == target, ],
                 paste0("Proximity to ", target, " by response"))
  ggsave(p, filename = file.path(PROX_DIR, paste0("ProximityBox_", safe_name(target),
                                                  "_byResponse.pdf")),
         width = 8, height = 5)
}

p <- plot_prox(prox[prox$Target_CT == FIG_TARGET & prox$celltype %in% FIG_SOURCES, ],
               prox_stats[prox_stats$Target_CT == FIG_TARGET, ],
               paste0("Proximity to ", FIG_TARGET), ord = FIG_SOURCES)
ggsave(p, filename = file.path(FIG_DIR, "Proximity_CD8TEM_byResponse.pdf"),
       width = 3.5, height = 4.5)
