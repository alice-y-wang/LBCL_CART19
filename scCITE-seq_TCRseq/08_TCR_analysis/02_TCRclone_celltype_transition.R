# =============================================================================
# Celltype transition matrices of CAR+ T cell clones along
#
#   Apheresis -> Product -> Peak -> Week 4
#
# -----------------------------------------------------------------------------
# HOW A CLONE GETS ITS CAR LABEL
#
# CAR status is a property of the CELL, but a trajectory needs it to be a
# property of the CLONE, because apheresis predates transduction and every
# apheresis cell is CAR- by construction. So a clone is CAR+ if ANY of its
# cells is CAR.exp TRUE at ANY of Product, Peak or Week 4 (CAR_DEFINING_TPS).
# The label is then carried backwards to apheresis.
#
# Once labelled, a CAR+ clone contributes only its CAR+ cells at post-manufacture
# timepoints (CELL_FILTER_BY_CLONE_CAR). At apheresis all cells of the clone are
# used, since the CAR+ compartment does not exist yet.
#
# -----------------------------------------------------------------------------
# WHAT THE TRANSITION MATRICES MEAN
#
# Cells are not tracked individually -- only clones are. So "movement" from
# celltype i to celltype j is inferred within a clone: if a clone has cells of
# type i at tp1, the destinations are that same clone's celltype composition at
# tp2. E.g. a clone with 3 proliferating cells at tp1 that is 2 CD8 TEM +
# 1 proliferating at tp2 contributes 0.67 to (proliferating -> CD8 TEM) and
# 0.33 to (proliferating -> proliferating).
#
# Each clone contributes in proportion to how many tp1 cells of type i it had
# (sum over clones of n_i(tp1) * p_j(tp2)), so each row is the proportion of
# tp1 cells of type i that end up as each celltype at tp2. Rows sum to 1.
# =============================================================================

set.seed(11)

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

# ---- CONFIG ----------------------------------------------------------------
PATHS_IN <- list(
  full_rds = file.path("seurat_objects", "LBCL_full_dataset_clean_new.RDS"),   # full clean object (metadata included)
  tcr_csv  = file.path("TCR_Data", "LBCL_TCR_metadata_clean_forupload.csv")
)
FIG_DIR  <- "TCR_Figures"
DATA_DIR <- "TCR_Data"
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

CELLTYPE_COL <- "cell.anno"
CLONE_COL    <- "CTaa"          # or CTstrict

# T/NK cell types (values of cell.anno) kept from the full object
TNK_TYPES <- c("CD4 Naive", "CD4 TEM-like", "Treg",
               "CD8 Naive", "CD8 TEM",
               "gdT", "MAIT", "NK", "Proliferating")

PATH_A <- c("Apheresis", "Product", "Peak", "Week 4")

# A clone is CAR+ if it has a CAR+ cell at ANY of these timepoints.
CAR_DEFINING_TPS <- c("Product", "Peak", "Week 4")

# TRUE: a CAR+ clone contributes only its CAR+ cells after manufacture.
# FALSE: all cells of the clone, whatever their call.
CELL_FILTER_BY_CLONE_CAR <- TRUE

# Weighting shown in the figures: "counts" = proportion of cells (each clone
# weighted by its number of tp1 cells); "clone average" = each clone counts equally.
# Both are written to the CSV.
PLOT_WEIGHTING <- "clone average"

# Timepoint labels in the object -> labels used here
TP_LABELS <- c(APH = "Apheresis", TDN = "Product", D0 = "Day 0",
               Peak = "Peak", `4W` = "Week 4", Week4 = "Week 4")

harmonize_tp <- function(x) {
  x   <- as.character(x)
  out <- unname(TP_LABELS[x])
  out[x %in% TP_LABELS] <- x[x %in% TP_LABELS]   # already in the target form
  unmapped <- setdiff(unique(x[is.na(out)]), NA)
  if (length(unmapped)) message("timepoints not mapped (dropped): ", paste(unmapped, collapse = ", "))
  out
}


# =============================================================================
# Load and label
# =============================================================================

stopifnot(file.exists(PATHS_IN$full_rds), file.exists(PATHS_IN$tcr_csv))

# Cell type, CAR.exp, patient and timepoint come from the object's metadata.
# The TCR csv contributes only the clonotype columns, joined on cell.bc.
# Cells with no clonotype drop out.
full_obj <- readRDS(PATHS_IN$full_rds)
base_md  <- full_obj@meta.data
rm(full_obj); gc(verbose = FALSE)

req_cols  <- c("patient.id", "timepoint", "CAR.exp", CELLTYPE_COL)
miss_cols <- setdiff(req_cols, names(base_md))
if (length(miss_cols)) stop("Metadata is missing column(s): ", paste(miss_cols, collapse = ", "))

base_md <- base_md[, !duplicated(names(base_md)), drop = FALSE]
if (!"cell.bc" %in% names(base_md)) base_md$cell.bc <- rownames(base_md)
base_md <- base_md %>%
  filter(.data[[CELLTYPE_COL]] %in% TNK_TYPES) %>%
  mutate(timepoint = harmonize_tp(timepoint))

tcr_cols <- read.csv(PATHS_IN$tcr_csv, row.names = 1) %>%
  dplyr::select(cell.bc, CTaa, CTstrict) %>%
  filter(!duplicated(cell.bc))

# take the clonotype columns from the TCR csv, not the object
base_md <- base_md[, !names(base_md) %in% c("CTaa", "CTstrict"), drop = FALSE]

cells <- base_md %>%
  dplyr::left_join(tcr_cols, by = "cell.bc") %>%
  mutate(patient_id = as.character(patient.id),
         celltype   = .data[[CELLTYPE_COL]],
         clone      = .data[[CLONE_COL]],
         car        = ifelse(CAR.exp, "CAR+", "CAR-")) %>%
  filter(!is.na(clone), clone != "", clone != "NA",
         !is.na(celltype), !is.na(CAR.exp), !is.na(timepoint))

message(nrow(cells), " cells with a celltype, a CAR call and a clonotype")
cat("\ncells per timepoint x CAR call:\n")
print(table(cells$timepoint, cells$car, useNA = "ifany"))

CT_LEVELS <- levels(factor(cells$celltype))
cells$celltype <- factor(cells$celltype, levels = CT_LEVELS)

# =============================================================================
# Clone-level CAR label: CAR+ if any CAR+ cell at Product, Peak OR Week 4.
# =============================================================================

clone_car <- cells %>%
  filter(timepoint %in% CAR_DEFINING_TPS) %>%
  group_by(patient_id, clone) %>%
  summarise(
    clone_car         = if (any(car == "CAR+")) "CAR+" else "CAR-",
    carpos_tps        = paste(sort(unique(as.character(timepoint[car == "CAR+"]))), collapse = "; "),
    carpos_in_product = any(car == "CAR+" & timepoint == "Product"),
    .groups = "drop")

cat("\nclones by CAR label (Product / Peak / Week 4):\n")
print(table(clone_car$clone_car))
cat("\nwhere the CAR+ evidence came from:\n")
print(clone_car %>% filter(clone_car == "CAR+") %>% count(carpos_tps, sort = TRUE))
cat("\nCAR+ clones that would have been MISSED using the product alone: ",
    sum(clone_car$clone_car == "CAR+" & !clone_car$carpos_in_product), "\n")


# =============================================================================
# CAR+ clone cells along the path
# =============================================================================

car_cells <- cells %>%
  filter(timepoint %in% PATH_A) %>%
  dplyr::inner_join(clone_car %>% filter(clone_car == "CAR+") %>%
                      dplyr::select(patient_id, clone, clone_car),
                    by = c("patient_id", "clone"))
if (CELL_FILTER_BY_CLONE_CAR) {
  car_cells <- car_cells %>% filter(timepoint == "Apheresis" | car == clone_car)
}
car_cells <- car_cells %>%
  mutate(timepoint = factor(as.character(timepoint), levels = PATH_A))

cat("\nCAR+ clone cells per timepoint:\n")
print(table(car_cells$timepoint))


# =============================================================================
# CELLTYPE TRANSITION MATRICES
#
# Rows = celltype at tp1, columns = celltype at tp2, within shared clones.
# =============================================================================

transition_matrix <- function(d, tp1, tp2) {
  cnt <- d %>%
    filter(timepoint %in% c(tp1, tp2)) %>%
    count(patient_id, clone, timepoint, celltype, name = "n")
  
  src <- cnt %>% filter(timepoint == tp1) %>%
    dplyr::select(patient_id, clone, from = celltype, n_from = n)
  
  dst <- cnt %>% filter(timepoint == tp2) %>%
    group_by(patient_id, clone) %>%
    mutate(p_to = n / sum(n)) %>%
    ungroup() %>%
    dplyr::select(patient_id, clone, to = celltype, p_to)
  
  j <- inner_join(src, dst, by = c("patient_id", "clone"))
  if (!nrow(j)) return(NULL)
  
  # --- weighting 1: by total counts ---
  m_counts <- j %>%
    group_by(from, to) %>%
    summarise(w = sum(n_from * p_to), .groups = "drop") %>%
    group_by(from) %>%
    mutate(prop = w / sum(w)) %>%
    ungroup() %>%
    mutate(weighting = "counts")
  
  # --- weighting 2: averaged over clones, each clone equal ---
  n_clones <- j %>% distinct(patient_id, clone, from) %>% count(from, name = "n_clones")
  m_avg <- j %>%
    distinct(patient_id, clone, from, to, p_to) %>%
    group_by(from, to) %>%
    summarise(sum_p = sum(p_to), .groups = "drop") %>%
    dplyr::left_join(n_clones, by = "from") %>%
    mutate(prop = sum_p / n_clones, w = NA_real_, weighting = "clone average") %>%
    dplyr::select(from, to, w, prop, weighting)
  
  bind_rows(m_counts, m_avg) %>%
    complete(from = factor(CT_LEVELS, levels = CT_LEVELS),
             to   = factor(CT_LEVELS, levels = CT_LEVELS),
             weighting, fill = list(prop = 0, w = 0)) %>%
    mutate(from = factor(from, levels = CT_LEVELS),
           to   = factor(to,   levels = CT_LEVELS),
           tp1 = tp1, tp2 = tp2) %>%
    dplyr::left_join(n_clones, by = "from")
}

plot_transition <- function(m, title) {
  ggplot(m, aes(to, from, fill = prop)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(prop >= 0.005, sprintf("%.2f", prop), "")),
              size = 2.6,
              colour = ifelse(m$prop > 0.5, "black", "white")) +
    scale_fill_viridis_c(limits = c(0, 1), option = "viridis", name = "proportion") +
    scale_y_discrete(limits = rev(CT_LEVELS)) +
    scale_x_discrete(position = "top") +
    coord_equal() +
    labs(title = title, x = paste("celltype at", unique(m$tp2)),
         y = paste("celltype at", unique(m$tp1))) +
    theme_minimal() +
    theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
          panel.grid = element_blank())
}

all_transitions <- list()
pairs <- lapply(seq_len(length(PATH_A) - 1), function(i) PATH_A[i:(i + 1)])

for (pr in pairs) {
  tp1 <- pr[1]; tp2 <- pr[2]
  m <- transition_matrix(car_cells, tp1, tp2)
  if (is.null(m)) { message("no shared CAR+ clones for ", tp1, " -> ", tp2); next }
  all_transitions[[paste(tp1, tp2)]] <- m
  
  tag <- gsub("[^A-Za-z0-9]", "", paste0(tp1, "_to_", tp2))
  ggsave(file.path(FIG_DIR, paste0("transition_CARpos_", tag, ".pdf")),
         plot_transition(filter(m, weighting == PLOT_WEIGHTING),
                         paste0("CAR+ clones:  ", tp1, "  ->  ", tp2)),
         width = 7, height = 6, limitsize = FALSE)
}

transitions <- bind_rows(all_transitions)
write_csv(transitions, file.path(DATA_DIR, "celltype_transition_matrices_CARpos.csv"))

# how many clones actually support each row? A 1.00 built on two clones is not
# a finding -- read this next to every matrix.
support <- transitions %>%
  distinct(tp1, tp2, from, n_clones) %>%
  arrange(tp1, desc(n_clones))
print(support, n = 50)
write_csv(support, file.path(DATA_DIR, "celltype_transition_support_CARpos.csv"))

saveRDS(list(car_cells   = car_cells,
             clone_car   = clone_car,
             transitions = transitions,
             support     = support),
        file.path(DATA_DIR, "clonal_transition_results_CARpos.rds"))