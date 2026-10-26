library(spatstat)
library(Seurat)

BASE_DIR   <- "/path/to/your/"
SCRIPT_DIR <- "/path/to/your/Scripts"
source(file.path(SCRIPT_DIR, "Utilities", "CODEX_Functions.R"))

SEURAT_DIR   <- file.path(BASE_DIR, "SeuratObj")
DISTANCE_DIR <- file.path(BASE_DIR, "Distance_Analysis", "Distance_Results")
dir.create(DISTANCE_DIR, recursive = TRUE, showWarnings = FALSE)

final_obj_path <- file.path(SEURAT_DIR, "DLBCL_CODEX_Final.rds")
meta_path      <- file.path(DISTANCE_DIR, "metadata_for_distance.rds")

MICRON_PER_PIXEL   <- 0.5069
N_PERMUTATIONS     <- 100
MIN_CELLS          <- 10
MAX_CELLS_PER_TYPE <- 50000

args <- commandArgs(TRUE)
target_cell_name <- args[1]
target_prefix    <- safe_name(target_cell_name)

# ==============================================================================
# 1. Per-Cell Coordinates and Labels
# ==============================================================================

if (!file.exists(meta_path)) {
  obj  <- readRDS(final_obj_path)
  meta <- obj@meta.data
  df <- data.frame(CellID  = rownames(meta),
                   sample  = as.character(meta[[SAMPLE_COL]]),
                   x.coord = as.numeric(meta$x.coord),
                   y.coord = as.numeric(meta$y.coord),
                   label   = as.character(meta[[L2_COL]]),
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$x.coord) & !is.na(df$y.coord), ]
  saveRDS(df, meta_path)
  rm(obj, meta); gc()
}
df <- readRDS(meta_path)
df <- df[!is.na(df$label), ]
sample_names <- sort(unique(df$sample))

# ==============================================================================
# 2. Minimum Distance and Poisson Point-Process Proximity
# ==============================================================================

prox_list <- list()
df_all    <- NULL

for (sample_name in sample_names) {
  sample_data <- df[df$sample == sample_name, ]
  cts <- sample_data[sample_data$label == target_cell_name, ]
  if (nrow(cts) < MIN_CELLS) next
  
  xmin <- min(cts$x.coord); xmax <- max(cts$x.coord)
  ymin <- min(cts$y.coord); ymax <- max(cts$y.coord)
  sample_data_crop <- sample_data[sample_data$x.coord > xmin & sample_data$x.coord < xmax &
                                  sample_data$y.coord > ymin & sample_data$y.coord < ymax, ]
  sample_data_crop_xy <- sample_data_crop[, c("x.coord", "y.coord")]

  win <- owin(c(xmin - 10, xmax + 10), c(ymin - 10, ymax + 10))
  cell_target_cts <- ppp(cts$x.coord, cts$y.coord, window = win, checkdup = FALSE)
  dtarget_fun     <- distfun(cell_target_cts)

  for (cell_type in unique(sample_data$label)) {
    in_crop <- sample_data_crop[sample_data_crop$label == cell_type, ]
    if (nrow(in_crop) <= MIN_CELLS) next

    if (nrow(in_crop) > MAX_CELLS_PER_TYPE) {
      in_crop <- in_crop[sample(nrow(in_crop), MAX_CELLS_PER_TYPE), ]
    }
    x <- in_crop$x.coord
    y <- in_crop$y.coord

    cells <- ppp(x, y, window = win, checkdup = FALSE)
    cm <- list(cell = cells, target = cell_target_cts, dtarget = dtarget_fun)
    f  <- tryCatch(ppm(cell ~ dtarget, data = cm), error = function(e) NULL)

    distance <- dtarget_fun(x, y)
    df_all <- rbind(df_all, data.frame(sample = sample_name, CellID = in_crop$CellID,
                                       source = cell_type, target = target_cell_name,
                                       distance_um = distance * MICRON_PER_PIXEL))

    d <- distance[is.finite(distance)]
    perm_med <- sapply(seq_len(N_PERMUTATIONS), function(k) {
      idx <- sample(nrow(sample_data_crop_xy), length(x))
      pd  <- dtarget_fun(sample_data_crop_xy[idx, 1], sample_data_crop_xy[idx, 2])
      median(pd[is.finite(pd)])
    })

    prox_list[[length(prox_list) + 1]] <- data.frame(
      sample          = sample_name,
      celltype        = cell_type,
      Target_CT       = target_cell_name,
      pop_size        = length(x),
      coefficients_1  = if (is.null(f)) NA else summary(f)$coef$Estimate[1],
      coefficients_2  = if (is.null(f)) NA else summary(f)$coef$Estimate[2],
      distance_median = median(d) * MICRON_PER_PIXEL,
      pvalue          = sum(perm_med < median(d)) / N_PERMUTATIONS)
  }
}

# ==============================================================================
# 3. Save Results
# ==============================================================================

saveRDS(do.call(rbind, prox_list),
        file.path(DISTANCE_DIR, paste0("Proximity_", target_prefix, ".rds")))
saveRDS(df_all, file.path(DISTANCE_DIR, paste0("Distance_ByCell_", target_prefix, ".rds")))
