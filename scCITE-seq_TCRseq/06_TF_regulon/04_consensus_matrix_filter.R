################################################################################
## pySCENIC consensus GRN: count TF|target pairs across 100 GRNBoost runs,
## filter to reproducible pairs, and aggregate their importance scores.
################################################################################

# --- Libraries ---------------------------------------------------------------
library(data.table)
library(dplyr)
library(readr)
library(ggplot2)

# --- CONFIG ------------------------------------------------------------------
pyscenic_dir <- "pyscenic"
grn_dir      <- file.path(pyscenic_dir, "grn_100_run")
images_dir   <- file.path(pyscenic_dir, "images")

dir.create(images_dir, recursive = TRUE, showWarnings = FALSE)

# A TF|target pair is kept if it appears in at least this many of the runs.
min_runs <- 80

# Compartment name -> subdirectory of adjacency matrices, plus a plot label.
compartments <- list(
  CAR    = list(adj_dir = "CAR_adj_mtx",    label = "CAR"),
  NonCAR = list(adj_dir = "NonCAR_adj_mtx", label = "non-CAR"),
  Mono   = list(adj_dir = "mono_adj_mtx",   label = "monocyte")
)

# ============================================================================ #
# 1. Count how many runs each TF|target pair appears in
# ============================================================================ #
count_list <- list()

for (nm in names(compartments)) {
  input_dir <- file.path(grn_dir, compartments[[nm]]$adj_dir)
  files     <- list.files(input_dir, pattern = "\\.csv$")
  message(nm, ": counting TF|target pairs across ", length(files), " runs")
  
  count_dt <- NULL
  
  for (f in files) {
    message("Processing: ", f)
    
    dt <- fread(file.path(input_dir, f), select = c("TF", "target"))
    dt[, TF_target := paste0(TF, "|", target)]
    
    # one vote per run, in case a pair is duplicated within a single file
    dt <- unique(dt, by = "TF_target")
    dt[, count := 1L]
    
    count_dt <- if (is.null(count_dt)) dt else rbind(count_dt, dt)
    
    # collapse incrementally so memory stays bounded across the 100 runs
    count_dt <- count_dt[, .(count = sum(count)), by = .(TF, target, TF_target)]
    print(dim(count_dt))
  }
  
  fwrite(count_dt, file.path(grn_dir, paste0(nm, "_TF_target_counts.csv")))
  count_list[[nm]] <- count_dt
}



# ============================================================================ #
# 2. Keep reproducible pairs and aggregate importance across runs
# ============================================================================ #
keep_list <- list()

for (nm in names(compartments)) {
  keep_list[[nm]] <- count_list[[nm]] %>% filter(count >= min_runs)
  
  input_dir <- file.path(grn_dir, compartments[[nm]]$adj_dir)
  files     <- list.files(input_dir, pattern = "\\.csv$")
  message(nm, ": summing importance for ", nrow(keep_list[[nm]]), " retained pairs")
  
  final_df <- NULL
  
  for (f in files) {
    message("Processing: ", f)
    
    dt <- fread(file.path(input_dir, f), select = c("TF", "target", "importance"))
    dt[, TF_target := paste0(TF, "|", target)]
    dt <- dt %>% filter(TF_target %in% keep_list[[nm]]$TF_target)
    dt$TF <- NULL
    dt$target <- NULL
    print(paste0("temp dataframe dim:", dim(dt)))
    
    final_df <- if (is.null(final_df)) dt else rbind(final_df, dt)
    
    final_df <- final_df %>%
      group_by(TF_target) %>%
      summarise(
        importance = sum(importance, na.rm = TRUE),
        .groups = "drop"
      )
    
    print(dim(final_df))
  }
  
  # mean importance over the runs the pair actually appeared in
  merge_df$importance <- merge_df$importance / merge_df$count
  write_csv(merge_df, file.path(grn_dir, paste0(nm, "_TF_target_importance_avg.csv")))
  
  merge_df$TF_target <- NULL
  merge_df$count <- NULL
  merge_df <- merge_df[, c("TF", "target", "importance")]
  merge_df <- merge_df %>% arrange(desc(importance))
  write_csv(merge_df, file.path(grn_dir, paste0(nm, "_consensus_adj.csv")))
}

