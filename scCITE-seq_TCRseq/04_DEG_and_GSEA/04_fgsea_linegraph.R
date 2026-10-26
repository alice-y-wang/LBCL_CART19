################################################################################
## fgsea results: TNF and interferon pathway NES over time
##
## NES score, TNF = Hallmark TNFa only, all interferon sets kept.
## Every pathway is drawn as its own line (version "c").
## Figures: CD8 TEM (left panel) next to CD14 Mono (CAR-, right panel) for
##   - D0 timeline,  CD8 TEM CAR-
##   - TDN timeline, CD8 TEM CAR-
##   - TDN timeline, CD8 TEM CAR+
################################################################################

library(dplyr)
library(ggplot2)


# 1. CONFIG -- edit here ----

fgsea_dir   <- file.path("gsea_analysis_max500", "fgsea_results")
file_suffix <- "_noMT_H_React_fgsea_0.01_pval_FULL.rds"   # FULL = all pathways, not only significant

fig_dir <- file.path("images", "gsea_line")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Timelines (timepoint names as in the .rds lists, in x-axis order)
timelines <- list(
  TDN = c("APH", "TDN", "Peak", "4W"),
  D0  = c("APH", "D0",  "Peak", "4W")
)

# Figures to make: timeline + which CD8 TEM results to use at every timepoint
#   "NonCAR" = CAR- cells, "CAR" = CAR+ cells
figures <- list(
  list(timeline = "D0",  group = "NonCAR"),
  list(timeline = "TDN", group = "NonCAR"),
  list(timeline = "TDN", group = "CAR")
)

# Cell types to load, and the annotation level of each file
cell_levels <- c("CD8 TEM"   = "Level2",
                 "CD14 Mono" = "Level1")

final_cell_type <- "CD8 TEM"     # left panel
final_mono      <- "CD14 Mono"   # right panel, always CAR-

# Set non-significant scores to 0. NULL = keep all scores as they are.
padj_cut <- 0.05

# Pathways in each group. Both Reactome naming styles are listed (plain-text
# from reactome.org, REACTOME_ from MSigDB); section 2 prints which were found.
pathway_sets <- list(
  TNF = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  IFN_alpha = c(
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "Interferon alpha/beta signaling",
    "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"
  ),
  IFN_gamma = c(
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "Interferon gamma signaling",
    "REACTOME_INTERFERON_GAMMA_SIGNALING"
  ),
  IFN_umbrella = c(
    "Interferon Signaling",
    "REACTOME_INTERFERON_SIGNALING"
  )
)

# Short display names for individual pathways
pathway_labels <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"         = "Hallmark TNFa via NF-kB",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"       = "Hallmark IFN-alpha response",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"       = "Hallmark IFN-gamma response",
  "Interferon alpha/beta signaling"          = "Reactome IFN-alpha/beta signaling",
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING" = "Reactome IFN-alpha/beta signaling",
  "Interferon gamma signaling"               = "Reactome IFN-gamma signaling",
  "REACTOME_INTERFERON_GAMMA_SIGNALING"      = "Reactome IFN-gamma signaling",
  "Interferon Signaling"                     = "Reactome IFN signaling",
  "REACTOME_INTERFERON_SIGNALING"            = "Reactome IFN signaling"
)


# 2. Load fgsea results (CD8 TEM and CD14 Mono, CAR- and CAR+, all timepoints) ----

read_fgsea <- function(group, cell_type) {
  f <- file.path(fgsea_dir,
                 paste0(group, "_", cell_levels[[cell_type]], "_", cell_type, file_suffix))
  if (!file.exists(f)) {
    message("No file (skipped): ", basename(f))
    return(NULL)
  }
  x <- readRDS(f)   # named list: one fgsea table per timepoint
  bind_rows(lapply(x, function(d) as.data.frame(d)[, c("pathway", "NES", "pval", "padj")]),
            .id = "timepoint") %>%
    mutate(group = group, cell_type = cell_type)
}

file_grid <- expand.grid(group     = c("NonCAR", "CAR"),
                         cell_type = names(cell_levels),
                         stringsAsFactors = FALSE)

gsea_raw <- bind_rows(mapply(read_fgsea, file_grid$group, file_grid$cell_type, SIMPLIFY = FALSE))

all_timepoints <- unique(unlist(timelines))
unknown_tp <- setdiff(unique(gsea_raw$timepoint), all_timepoints)
if (length(unknown_tp)) warning("Timepoints not in any timeline (dropped): ",
                                paste(unknown_tp, collapse = ", "))

gsea_raw <- gsea_raw %>% filter(timepoint %in% all_timepoints)

# NES score; non-significant scores set to 0
gsea_all <- gsea_raw %>% mutate(score = NES)
if (!is.null(padj_cut)) gsea_all$score[is.na(gsea_all$padj) | gsea_all$padj >= padj_cut] <- 0

# Keep the configured pathways and label them
pathway_lookup <- stack(pathway_sets) %>% dplyr::rename(pathway = values, pathway_group = ind)
print(pathway_lookup %>% mutate(found = pathway %in% gsea_all$pathway))

gsea_sel <- gsea_all %>%
  inner_join(pathway_lookup, by = "pathway") %>%
  mutate(pathway_label = pathway_labels[pathway])

# Rows per cell type x timepoint x analysis
print(with(gsea_sel, table(paste(cell_type, group),
                           factor(timepoint, levels = all_timepoints))))

# Colors: TNF in red, IFN pathways in blues
tnf_label    <- unique(pathway_labels[intersect(pathway_sets$TNF, gsea_sel$pathway)])
ifn_labels   <- setdiff(unique(gsea_sel$pathway_label), tnf_label)
line_colors  <- c(setNames(rep("#e31a1c", length(tnf_label)), tnf_label),
                  setNames(colorRampPalette(c("lightskyblue", "navy"))(length(ifn_labels)), ifn_labels))


# 3. Plot Line -- CD8 TEM and CD14 Mono, every pathway on its own ----

for (fg in figures) {
  tps       <- timelines[[fg$timeline]]
  car_label <- if (fg$group == "CAR") "CAR+" else "CAR-"
  cd8_panel  <- paste0(final_cell_type, " (", car_label, ")")
  mono_panel <- paste0(final_mono, " (CAR-)")
  
  lines_5 <- gsea_sel %>%
    filter(timepoint %in% tps,
           (cell_type == final_cell_type & group == fg$group) |
             (cell_type == final_mono      & group == "NonCAR")) %>%
    mutate(timepoint  = factor(timepoint, levels = tps),
           car_status = ifelse(group == "CAR", "CAR+", "CAR-"),
           panel      = factor(ifelse(cell_type == final_mono, mono_panel, cd8_panel),
                               levels = c(cd8_panel, mono_panel)),
           line       = pathway_label,
           mean_score = score)
  
  p5 <- ggplot(lines_5, aes(timepoint, mean_score, color = line, group = line)) +
    geom_line(linewidth = 1.1) + geom_point(size = 2) +
    facet_wrap(~ panel, nrow = 1) +
    scale_color_manual(values = line_colors) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_classic() +
    labs(x = "Timepoint", y = "NES", color = "Pathway",
         title = paste0("TNF and interferon pathway enrichment over time (", fg$timeline, " timeline)"),
         subtitle = "All pathways shown individually")
  
  file_tag <- paste0(fg$timeline, "_CD8TEM_", if (fg$group == "CAR") "CARpos" else "CARneg", "_AllPathways")
  
  pdf(file.path(fig_dir, paste0("GSEA_Final_", file_tag, "_CD14Mono.pdf")), width = 11, height = 4.5)
  print(p5)
  dev.off()
  
  write.csv(lines_5, file.path(fig_dir, paste0("GSEA_", file_tag, "_scores.csv")), row.names = FALSE)
}