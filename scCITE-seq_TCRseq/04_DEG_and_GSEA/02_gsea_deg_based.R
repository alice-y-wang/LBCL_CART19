################################################################################
## Lymphoma CAR-T single-cell RNA-seq: GSEA (fgsea) and pathway visualization
##
## Run 01_setup_and_degs.R first (loads shared objects).
##
## PATHS:
##   All input/output locations are defined in the CONFIG block below, relative
##   to `analysis_dir`. Edit `analysis_dir` to point at your project root.
################################################################################

library(fgsea)
library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(tidyr)
library(tibble)
library(pheatmap)

# ============================================================================ #
# CONFIG: directories
# ============================================================================ #
# Project / analysis root. Set to an absolute path or leave as "." if the
# working directory is already set by 00_setup.R.
analysis_dir <- "."

# Input data
gmt_dir       <- file.path(analysis_dir, "msigdb_v2024.1.Hs_GMTs")  # MSigDB GMTs
deg_input_dir <- file.path(analysis_dir, "deg_summary")             # summarized DEG CSVs

# Outputs
gsea_dir           <- file.path(analysis_dir, "gsea_analysis")
gsea_table_dir     <- file.path(gsea_dir, "gsea_tables")          # plotGseaTable PDFs
fgsea_out_dir      <- file.path(gsea_dir, "fgsea_results")        # per-celltype fgsea .rds
heatmap_dir        <- file.path(gsea_dir, "heatmaps")             # NES / p-value heatmaps
heatmap_simple_dir <- file.path(gsea_dir, "heatmaps_simplified")  # simplified-pathway heatmaps

# Make sure top-level output directories exist.
for (d in c(gsea_table_dir, fgsea_out_dir, heatmap_dir, heatmap_simple_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Tag appended to output file names.
run_tag <- format(Sys.Date(), "%Y%m%d")

# Single source of truth for fgsea .rds naming, used by both the writer in
# section 1 and the readers in sections 3-5.
fgsea_suffix   <- "_noMT_H_React_fgsea_0.01_pval"
fgsea_rds_name <- function(CAR_status, level, celltype) {
  paste0(CAR_status, "_", level, "_", celltype, fgsea_suffix, ".rds")
}
fgsea_file <- function(...) file.path(fgsea_out_dir, ...)

# ============================================================================ #
# 1. Run GSEA: RPS retained, MT/LINC removed, p < 0.01 DEGs; Hallmark + REACTOME
# ============================================================================ #

# Load MSigDB gene-set collections (downloaded GMT files).
hallmark_pathways <- gmtPathways(file.path(gmt_dir, "h.all.v2024.1.Hs.symbols.gmt"))
reactome_pathways <- gmtPathways(file.path(gmt_dir, "ReactomePathways.gmt"))
combined_pathways <- c(hallmark_pathways, reactome_pathways)

# Directory of summarized DEG tables (one CSV per CAR-status x level).
list_dfs <- list.files(deg_input_dir, pattern = "\\.csv$")

# Cell types already processed at a coarser level; defined in 00_setup.R.
if (!exists("already_seen")) already_seen <- character(0)


#' HELPER FUNCTION FOR Running fgsea per timepoint for a single cell type.
#'
#' @param df          DEG data frame for one cell type (cols: celltype,
#'                     timepoint, geneID, avg_log2FC, ...).
#' @param cell        Cell-type label to filter on.
#' @param CAR_status  "CAR" / "NonCAR", used in output paths.
#' @param level       Annotation-level label, used in output paths.
#' @param image       If TRUE, also writes a plotGseaTable PDF per timepoint.
get_gsea_results <- function(df, cell, CAR_status, level, image = FALSE){
  list_results <- list()
  timepoints <- unique(df$timepoint)
  print(timepoints)
  for (tp in timepoints){
    print(tp)
    
    timepoint.df <- df %>% dplyr::filter(celltype == cell & timepoint == tp)
    print(unique(timepoint.df$celltype))
    
    rank <- timepoint.df$avg_log2FC
    names(rank) <- timepoint.df$geneID
    rank <- sort(rank, decreasing = TRUE)
    
    hallmark.results <- fgsea(pathways = hallmark_pathways,
                              stats = rank,
                              minSize = 5,
                              maxSize = 100)
    
    reactome.results <- fgsea(pathways = reactome_pathways,
                              stats = rank,
                              minSize = 5,
                              maxSize = 100)
    
    full.results <- rbind(hallmark.results, reactome.results)
    full.results <- full.results %>% dplyr::filter(padj < 0.1)
    
    list_results <- append(list_results, list(full.results))
    
    if (image && nrow(full.results) > 0){
      topPathwaysUp   <- full.results[ES > 0][head(order(pval), n = 10), pathway]
      topPathwaysDown <- full.results[ES < 0][head(order(pval), n = 10), pathway]
      topPathways     <- c(topPathwaysUp, rev(topPathwaysDown))
      
      table_subdir <- file.path(gsea_table_dir, paste0(CAR_status, "_", level))
      dir.create(table_subdir, showWarnings = FALSE, recursive = TRUE)
      pdf(file.path(table_subdir, paste0(cell, "_", tp, ".pdf")), height = 12, width = 15)
      print(plotGseaTable(combined_pathways[topPathways], rank, full.results,
                          gseaParam = 0.5, pathwayLabelStyle = list(size = 10),
                          headerLabelStyle = list(color = "red"),
                          valueStyle = list(size = 10),
                          axisLabelStyle = list(size = 6)))
      dev.off()
    }
  }
  names(list_results) <- timepoints
  return(list_results)
}

# run GSEA, per celltype, keeping RPS genes, and saving
for (f in list_dfs){
  
  print(f)
  full_df <- read_csv(file.path(deg_input_dir, f))
  full_df$log10_padj <- -log10(full_df$p_val_adj)
  
  celltypes  <- unique(full_df$celltype)
  timepoints <- unique(full_df$timepoint)
  
  CAR_status <- str_extract(f, "^(NonCAR|CAR)")
  level      <- str_extract(f, "(?i)level\\d") %>% str_to_title()
  
  print(CAR_status)
  print(level)
  
  full_df <- full_df %>%
    mutate(
      log10_padj = pmin(log10_padj, 310)
    )
  
  full_df <- full_df[
    !grepl("^MT", full_df$geneID, ignore.case = TRUE) &
      !grepl("^LINC", full_df$geneID, ignore.case = TRUE),
  ]
  
  full_df <- full_df %>% dplyr::filter(p_val_adj < 0.01)
  
  for (ct in celltypes){
    if ((ct %in% already_seen) & (level == "Level3")){
      next
    }
    ct_df <- full_df %>% dplyr::filter(celltype == ct)
    ct_result <- get_gsea_results(df = ct_df, cell = ct, CAR_status = CAR_status,
                                  level = level, image = TRUE)
    
    print(names(ct_result))
    saveRDS(ct_result, fgsea_file(fgsea_rds_name(CAR_status, level, ct)))
  }
}


# ============================================================================ #
# 2. Revisit GSEA results and collapse redundant (overlapping) pathways
# ============================================================================ #

# Helpers
get_gene_set_collection <- function(pathway,
                                    hallmark_pattern  = "^HALLMARK",
                                    reactome_pattern  = "^REACTOME",
                                    default_collection = "Reactome") {
  pathway <- as.character(pathway)
  out <- rep(default_collection, length(pathway))
  out[str_detect(pathway, regex(reactome_pattern, ignore_case = TRUE))] <- "Reactome"
  out[str_detect(pathway, regex(hallmark_pattern, ignore_case = TRUE))] <- "Hallmark"
  out
}

get_top_fgsea <- function(fgsea_list, top_n = 10) {
  
  top_paths <- imap(
    fgsea_list,
    ~ {
      df <- .x %>% filter(!is.na(padj))
      
      bind_rows(
        df %>%
          filter(NES > 0) %>%
          arrange(padj) %>%
          slice_head(n = top_n),
        
        df %>%
          filter(NES < 0) %>%
          arrange(padj) %>%
          slice_head(n = top_n)
      ) %>%
        pull(pathway)
    }
  ) %>%
    unlist() %>%
    unique()
  
  fgsea_long <- imap_dfr(
    fgsea_list,
    ~ .x %>%
      filter(pathway %in% top_paths) %>%
      mutate(timepoint = .y)
  )
  
  return(fgsea_long)
}

get_NES_matrix <- function(top_fgsea){
  
  nes_mat <- top_fgsea %>%
    select(pathway, timepoint, NES) %>%
    tidyr::pivot_wider(
      names_from = timepoint,
      values_from = NES
    ) %>%
    tibble::column_to_rownames("pathway") %>%
    as.matrix()
  return(nes_mat)
}

get_pval_matrix <- function(top_fgsea){
  
  top_fgsea$log10_padj <- -log10(top_fgsea$padj)
  pval_mat <- top_fgsea %>%
    dplyr::select(pathway, timepoint, log10_padj) %>%
    tidyr::pivot_wider(
      names_from = timepoint,
      values_from = log10_padj
    ) %>%
    tibble::column_to_rownames("pathway") %>%
    as.matrix()
  return(pval_mat)
}

get_signed_pval_matrix <- function(top_fgsea){
  # sign follows NES direction; magnitude is -log10(padj)
  top_fgsea$signed_log10_padj <- sign(top_fgsea$NES) * -log10(top_fgsea$padj)
  
  signed_mat <- top_fgsea %>%
    dplyr::select(pathway, timepoint, signed_log10_padj) %>%
    tidyr::pivot_wider(
      names_from  = timepoint,
      values_from = signed_log10_padj
    ) %>%
    tibble::column_to_rownames("pathway") %>%
    as.matrix()
  return(signed_mat)
}


# ---- Helper: collapse redundant pathways by leading-edge overlap ----------- #

#' Greedily remove pathways whose leading-edge genes heavily overlap a kept one.#-----

simplify_fgsea_best <- function(reactome_df, min_genes = 5, overlap_cutoff = 0.8,
                                use_max_set=TRUE, pval=TRUE, debug_pathway = NULL, desc=TRUE,
                                final_pval_filter = FALSE) {
  results <- reactome_df %>%
    mutate(pathway = as.character(pathway)) %>%
    filter(!is.na(padj))
  
 
  if(pval & desc) {
    results <- results %>%
      mutate(
        # gene_list = strsplit(leadingEdge, ","),
        gene_list = leadingEdge,
        gene_list = map(gene_list, unique),
        num_genes_present = map_int(gene_list, length)
      ) %>%
      filter(num_genes_present >= min_genes) %>%
      arrange(padj, desc(num_genes_present))
    #arrange(desc(num_genes_present),  padj)
  } else if (pval & !desc) {
    results <- results %>%
      mutate(
        # gene_list = strsplit(leadingEdge, ","),
        gene_list = leadingEdge,
        gene_list = map(gene_list, unique),
        num_genes_present = map_int(gene_list, length)
      ) %>%
      filter(num_genes_present >= min_genes) %>%
      arrange(padj, num_genes_present)
    #arrange(num_genes_present, padj)
  } else if (!pval & desc) {
    results <- results %>%
      mutate(
        # gene_list = strsplit(leadingEdge, ","),
        gene_list = leadingEdge,
        gene_list = map(gene_list, unique),
        num_genes_present = map_int(gene_list, length)
      ) %>%
      filter(num_genes_present > min_genes) %>%
      arrange(desc(num_genes_present))
  }
  else {
    results <- results %>%
      mutate(
        # gene_list = strsplit(leadingEdge, ","),
        gene_list = leadingEdge,
        gene_list = map(gene_list, unique),
        num_genes_present = map_int(gene_list, length)
      ) %>%
      filter(num_genes_present > min_genes) %>%
      arrange(num_genes_present)
  }
  
  
  if (nrow(results) == 0) {
    return(list(
      kept_pathways = results,
      removed_pathways = data.frame(),
      summary = "Kept 0/0 pathways (removed 0)"
    ))
  }
  
  is_kept <- rep(TRUE, nrow(results))
  removal_reasons <- rep(NA_character_, nrow(results))
  
  for (i in seq_len(nrow(results) - 1)) {
    if (!is_kept[i]) next
    genes_i <- results$gene_list[[i]]
    pathway_i <- results$pathway[[i]]
    padj_i <- results$padj[[i]]
    for (j in (i + 1):nrow(results)) {
      if (!is_kept[j]) next
      genes_j <- results$gene_list[[j]]
      pathway_j <- results$pathway[[j]]
      padj_j <- results$padj[[j]]
      if(use_max_set) {
        overlap_fraction <- length(intersect(genes_i, genes_j)) / max(length(genes_i), length(genes_j)) # largest set overlap
      } else {
        overlap_fraction <- length(intersect(genes_i, genes_j)) / min(length(genes_i), length(genes_j)) # smallest set overlap
      }
      if (!is.na(overlap_fraction) && overlap_fraction > overlap_cutoff) {
        
        if(final_pval_filter){
          # Keep whichever of the pair has the smaller padj.
          if( padj_i <= padj_j){
            is_kept[j] <- FALSE
            removal_reasons[j] <- sprintf("Removed by %s (%.2f overlap)", pathway_i, overlap_fraction)
          } else {
            is_kept[i] <- FALSE
            removal_reasons[i] <- sprintf("Removed by %s (%.2f p-value)", pathway_j, padj_i - padj_j)
          }
        } else {
          # Default: keep the earlier-ranked pathway, drop the later one.
          is_kept[j] <- FALSE
          removal_reasons[j] <- sprintf("Removed by %s (%.2f overlap)", pathway_i, overlap_fraction)
        }
        
        
        # Debug specific pathway
        if (!is.null(debug_pathway) && results$pathway[[j]] == debug_pathway) {
          cat("DEBUG:", pathway_j, "removed by", pathway_i,
              "- overlap fraction:", overlap_fraction, "\n")
          cat("Genes in", pathway_j, ":", results$gene_list[[j]], "...\n")
          cat("Genes in keeper:", results$gene_list[[i]],  "...\n")
        }
      }
    }
  }
  
  kept_pathways <- results[is_kept, ] %>%
    select(-gene_list, -num_genes_present)
  
  removed_pathways <- results[!is_kept, ] %>%
    select(-gene_list, -num_genes_present)
  
  list(
    kept_pathways = kept_pathways,
    removed_pathways = removed_pathways,
    summary = sprintf("Kept %d/%d pathways (removed %d)",
                      nrow(kept_pathways), nrow(results), nrow(removed_pathways))
  )
  
}

#' Run simplify_fgsea_best() per timepoint and per gene set collection.
simplify_all_fgsea <- function(fgsea_results, colname = NULL, overlap_cutoff = 0.8, min_genes = 5,
                               use_max_set = TRUE, pval = TRUE, desc = TRUE, final_pval_filter = FALSE,
                               debug_pathway = NULL,
                               collection_fun = get_gene_set_collection,
                               collection_levels = c("Hallmark", "Reactome")) {
  list_keep_results   <- list()
  list_remove_results <- list()
  list_summary        <- list()
  timepoints <- names(fgsea_results)
  for (type in timepoints) {
    print(type)
    results_filter <- fgsea_results[[type]]
    
    if (is.null(results_filter) || nrow(results_filter) < 1) {
      list_keep_results[[type]]   <- results_filter
      list_remove_results[[type]] <- results_filter
      next
    }
    
    results_filter <- results_filter %>% filter(padj < 0.1)
    
    if (nrow(results_filter) < 1) {
      list_keep_results[[type]]   <- results_filter
      list_remove_results[[type]] <- results_filter
      next
    }
    
    # Tag collection and prune each collection independently
    results_filter <- results_filter %>%
      mutate(gene_set_collection = collection_fun(pathway))
    
    collections <- unique(results_filter$gene_set_collection)
    collections <- c(intersect(collection_levels, collections),
                     setdiff(collections, collection_levels))
    
    keep_by_collection    <- list()
    remove_by_collection  <- list()
    summary_by_collection <- character(0)
    
    for (coll in collections) {
      sub_res <- results_filter %>% filter(gene_set_collection == coll)
      if (nrow(sub_res) < 1) next
      
      sub_simplify <- simplify_fgsea_best(sub_res, min_genes = min_genes,
                                          overlap_cutoff = overlap_cutoff,
                                          use_max_set = use_max_set,
                                          pval = pval, debug_pathway = debug_pathway,
                                          desc = desc, final_pval_filter = final_pval_filter)
      
      keep_by_collection[[coll]]   <- sub_simplify$kept_pathways
      remove_by_collection[[coll]] <- sub_simplify$removed_pathways
      summary_by_collection[coll]  <- sprintf(
        "%s | %s: %s", type, coll,
        if (is.null(sub_simplify$summary)) "Kept 0/0 pathways (removed 0)" else sub_simplify$summary
      )
      cat("   ", coll, ": ", summary_by_collection[coll], "\n", sep = "")
    }
    
    list_keep_results[[type]]   <- dplyr::bind_rows(keep_by_collection)
    list_remove_results[[type]] <- dplyr::bind_rows(remove_by_collection)
    list_summary[[type]]        <- summary_by_collection
  }
  print(length(list_keep_results))
  
  list(
    full_keep      = list_keep_results,
    full_remove    = list_remove_results,
    full_summaries = unlist(list_summary, use.names = FALSE)
  )
}


# cluster gene sets within databases # -----
cluster_within_groups <- function(mat, groups,
                                  group_order = c("Hallmark", "Reactome"),
                                  distance = "euclidean", method = "complete") {
  groups  <- as.character(groups)
  present <- c(intersect(group_order, unique(groups)),
               setdiff(unique(groups), group_order))
  
  idx   <- integer(0)
  sizes <- integer(0)
  for (g in present) {
    ii <- which(groups == g)
    if (length(ii) > 2) {
      d <- dist(mat[ii, , drop = FALSE], method = distance)
      if (any(is.na(d))) d[is.na(d)] <- max(d, na.rm = TRUE)
      hc <- hclust(d, method = method)
      ii <- ii[hc$order]
    }
    idx   <- c(idx, ii)
    sizes <- c(sizes, length(ii))
  }
  gaps <- cumsum(sizes)
  gaps <- gaps[-length(gaps)]   # no gap after the last block
  
  list(order  = idx,
       gaps   = gaps,
       groups = groups[idx],
       sizes  = setNames(sizes, present))
}

# pretty row labels (kept the HALLMARK/REACTOME prefix so blocks stay readable) # -----
pretty_pathway <- function(x) make.unique(str_to_title(gsub("_", " ", x)))

make_row_ann <- function(mat_ord, ord) {
  data.frame(`Gene set` = factor(ord$groups, levels = unique(ord$groups)),
             row.names = rownames(mat_ord),
             check.names = FALSE)
}
ann_colors <- list(`Gene set` = c(Hallmark = "#8C6BB1", Reactome = "#7FBC79"))

#' Build the NES / pval / signed-pval heatmaps for a set of columns.
plot_grouped_fgsea_heatmaps <- function(combined, celltype, columns_order, column_groups,
                                        out_pdf = NULL, cap = 10,
                                        width = 15, height = 10) {
  
  stopifnot(length(columns_order) == length(column_groups))
  
  df         <- get_top_fgsea(combined)
  nes_mtx    <- get_NES_matrix(df)
  pval_mtx   <- get_pval_matrix(df)
  signed_mtx <- get_signed_pval_matrix(df)
  
  # keep all three matrices on the same row universe / order
  common_rows <- rownames(nes_mtx)
  pval_mtx    <- pval_mtx[common_rows, , drop = FALSE]
  signed_mtx  <- signed_mtx[common_rows, , drop = FALSE]
  
  nes_mtx[is.na(nes_mtx)]       <- 0
  pval_mtx[is.na(pval_mtx)]     <- 0
  signed_mtx[is.na(signed_mtx)] <- 0
  
  # cap magnitude so a near-zero padj (or padj == 0 -> Inf) 
  pval_mtx[pval_mtx > cap]      <-  cap
  signed_mtx[signed_mtx >  cap] <-  cap
  signed_mtx[signed_mtx < -cap] <- -cap
  
  # drop requested columns that don't exist, keep the group labels in sync
  keep_cols     <- columns_order %in% colnames(nes_mtx)
  columns_order <- columns_order[keep_cols]
  column_groups <- column_groups[keep_cols]
  gaps_col      <- head(cumsum(rle(as.character(column_groups))$lengths), -1)
  
  # which collection each heatmap row belongs to
  row_collection <- get_gene_set_collection(rownames(nes_mtx))
  print(table(row_collection))   # sanity check the Hallmark / Reactome split
  
  # cluster rows within Hallmark and within Reactome separately 
  ord_nes    <- cluster_within_groups(nes_mtx[, columns_order, drop = FALSE],    row_collection)
  ord_pval   <- cluster_within_groups(pval_mtx[, columns_order, drop = FALSE],   row_collection)
  ord_signed <- cluster_within_groups(signed_mtx[, columns_order, drop = FALSE], row_collection)
  
  nes_mtx_ord    <- nes_mtx[ord_nes$order, , drop = FALSE]
  pval_mtx_ord   <- pval_mtx[ord_pval$order, , drop = FALSE]
  signed_mtx_ord <- signed_mtx[ord_signed$order, , drop = FALSE]
  
  rownames(nes_mtx_ord)    <- pretty_pathway(rownames(nes_mtx_ord))
  rownames(pval_mtx_ord)   <- pretty_pathway(rownames(pval_mtx_ord))
  rownames(signed_mtx_ord) <- pretty_pathway(rownames(signed_mtx_ord))
  
  if (!is.null(out_pdf)) { pdf(out_pdf, width = width, height = height); on.exit(dev.off()) }
  
  print(pheatmap::pheatmap(
    nes_mtx_ord[, columns_order],
    breaks = seq(-2.5, 2.5, length.out = 99),
    color = colorRampPalette(c("#E7A75E", "white", "#4FB7C5"))(100),
    cluster_rows = FALSE,          # rows already clustered within each gene set
    cluster_cols = FALSE,
    gaps_col = gaps_col,
    gaps_row = ord_nes$gaps,       # gap between Hallmark and Reactome blocks
    annotation_row = make_row_ann(nes_mtx_ord, ord_nes),
    annotation_colors = ann_colors, na_col = "grey90",
    main = paste0(celltype, " fgsea NES of top10 significant at each timepoint")
  ))
  
  grid::grid.newpage()
  print(pheatmap::pheatmap(
    pval_mtx_ord[, columns_order],
    breaks = seq(0, cap, length.out = 99),
    color = colorRampPalette(c("white", "red"))(100),
    cluster_rows = FALSE, cluster_cols = FALSE,
    gaps_col = gaps_col, gaps_row = ord_pval$gaps,
    annotation_row = make_row_ann(pval_mtx_ord, ord_pval),
    annotation_colors = ann_colors, na_col = "grey90",
    main = paste0(celltype, " fgsea top10 at each timepoint pval")
  ))
  
  grid::grid.newpage()
  print(pheatmap::pheatmap(
    signed_mtx_ord[, columns_order],
    breaks = seq(-cap, cap, length.out = 99),
    color = colorRampPalette(c("#E7A75E", "white", "#4FB7C5"))(100),
    cluster_rows = FALSE, cluster_cols = FALSE,
    gaps_col = gaps_col, gaps_row = ord_signed$gaps,
    annotation_row = make_row_ann(signed_mtx_ord, ord_signed),
    annotation_colors = ann_colors, na_col = "grey90",
    main = paste0(celltype, " fgsea signed -log10(padj) [sign(NES)] of top10 at each timepoint")
  ))
  
  invisible(list(nes = nes_mtx_ord, pval = pval_mtx_ord, signed = signed_mtx_ord,
                 columns_order = columns_order, gaps_col = gaps_col,
                 row_collection = row_collection))
}

# ============================================================================ #
# 3. CD8 EM (CAR vs non-CAR)
# ============================================================================ #

simplify <- "SIMPLE"

noncar_cd8em <- readRDS(fgsea_file(fgsea_rds_name("NonCAR", "Level2", "CD8 EM")))
car_cd8em    <- readRDS(fgsea_file(fgsea_rds_name("CAR",    "Level2", "CD8 EM")))

names(noncar_cd8em) <- paste0("NonCAR_", names(noncar_cd8em))
names(car_cd8em)    <- paste0("CAR_",    names(car_cd8em))

combined <- c(
  lapply(car_cd8em,    \(x) x %>% filter(padj < 0.05)),
  lapply(noncar_cd8em, \(x) x %>% filter(padj < 0.05))
)

if (simplify == "SIMPLE") {
  # simplification is per celltype x timepoint x (Hallmark | Reactome)
  combined_simple <- simplify_all_fgsea(combined, overlap_cutoff = 0.8, min_genes = 3)
  print(combined_simple$full_summaries)
  combined <- combined_simple$full_keep
}


cd8em_cols   <- c("CAR_TDN", "CAR_Peak", "CAR_4W",
                  "NonCAR_APH", "NonCAR_TDN", "NonCAR_D0", "NonCAR_Peak")
cd8em_groups <- c(rep("CAR", 3), rep("NonCAR", 4))

plot_grouped_fgsea_heatmaps(
  combined,
  celltype      = "CD8 EM All",
  columns_order = cd8em_cols,
  column_groups = cd8em_groups,
  out_pdf = file.path(heatmap_simple_dir,
                      paste0(simplify, "_CD8_EM_HReact_01pval_heatmap_", run_tag, ".pdf"))
)

# ============================================================================ #
# 4. Myeloid (CD14 / CD16 / DC)
# ============================================================================ #

CD14_fgsea_level1 <- readRDS(fgsea_file(fgsea_rds_name("NonCAR", "Level1", "CD14 Mono")))
CD16_fgsea_Level1 <- readRDS(fgsea_file(fgsea_rds_name("NonCAR", "Level1", "CD16 Mono")))
DC_fgsea_Level1   <- readRDS(fgsea_file(fgsea_rds_name("NonCAR", "Level1", "Dendritic Cell")))

names(CD14_fgsea_level1) <- paste0("CD14_", names(CD14_fgsea_level1))
names(CD16_fgsea_Level1) <- paste0("CD16_", names(CD16_fgsea_Level1))
names(DC_fgsea_Level1)   <- paste0("DC_",   names(DC_fgsea_Level1))

combined <- c(
  lapply(CD14_fgsea_level1, \(x) x %>% filter(padj < 0.05)),
  lapply(CD16_fgsea_Level1, \(x) x %>% filter(padj < 0.05)),
  lapply(DC_fgsea_Level1,   \(x) x %>% filter(padj < 0.05))
)

simplify <- "SIMPLE"
if (simplify == "SIMPLE") {
  # prunes within celltype x timepoint x (Hallmark | Reactome)
  combined_simple <- simplify_all_fgsea(combined, overlap_cutoff = 0.8, min_genes = 3)
  print(combined_simple$full_summaries)
  combined <- combined_simple$full_keep
}

# edit these two to match your actual timepoint names; missing ones are dropped
myeloid_cols   <- c("CD14_APH", "CD14_D0", "CD14_Peak", "CD14_4W",
                    "CD16_APH", "CD16_D0", "CD16_Peak",
                    "DC_APH",   "DC_D0",   "DC_Peak",   "DC_4W")
myeloid_groups <- c(rep("CD14", 4), rep("CD16", 3), rep("DC", 4))

plot_grouped_fgsea_heatmaps(
  combined,
  celltype      = "Myeloid (CD14 / CD16 / DC)",
  columns_order = myeloid_cols,
  column_groups = myeloid_groups,
  out_pdf = file.path(heatmap_simple_dir,
                      paste0(simplify, "_Myeloid_HReact_01pval_heatmap_GroupedGeneSets_", run_tag, ".pdf"))
)

# ============================================================================ #
# 5. CD4 EM-like
# ============================================================================ #

noncar_cd4em <- readRDS(fgsea_file(fgsea_rds_name("NonCAR", "Level2", "CD4 EMCM")))
car_cd4em    <- readRDS(fgsea_file(fgsea_rds_name("CAR",    "Level2", "CD4 EMCM")))

names(noncar_cd4em) <- paste0("NonCAR_", names(noncar_cd4em))
names(car_cd4em)    <- paste0("CAR_",    names(car_cd4em))

combined <- c(
  lapply(car_cd4em,    \(x) x %>% filter(padj < 0.05)),
  lapply(noncar_cd4em, \(x) x %>% filter(padj < 0.05))
)

simplify <- "SIMPLE"
if (simplify == "SIMPLE") {
  combined_simple <- simplify_all_fgsea(combined, overlap_cutoff = 0.8, min_genes = 3)
  print(combined_simple$full_summaries)
  combined <- combined_simple$full_keep
}

cd4_cols   <- c("CAR_TDN", "CAR_Peak", "NonCAR_APH", "NonCAR_TDN", "NonCAR_D0", "NonCAR_Peak")
cd4_groups <- c("CAR", "CAR", "NonCAR", "NonCAR", "NonCAR", "NonCAR")

plot_grouped_fgsea_heatmaps(
  combined,
  celltype      = "CD4 EM-like All",
  columns_order = cd4_cols,
  column_groups = cd4_groups,
  out_pdf = file.path(heatmap_simple_dir,
                      paste0(simplify, "_CD4_EMlike_HReact_01pval_heatmap_GroupedGeneSets_", run_tag, ".pdf"))
)