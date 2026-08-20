
# ---- Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(data.table)
  library(fgsea)
  library(ggvenn)
})

# ---- CONFIG -----------------------------------------------------------------
analysis_dir <- "."
setwd(analysis_dir)

ref_dir  <- "GSE197268_data"
gmt_dir  <- "msigdb_v2024.1.Hs_GMTs"
base_out <- file.path(ref_dir, "gsea_results")

deg_out_dir    <- file.path(base_out, "DEGs")
gsea_out_dir   <- file.path(base_out, "GSEA")
gsea_table_dir <- file.path(base_out, "GSEA_tables")

for (d in c(deg_out_dir, gsea_out_dir, gsea_table_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

f_ref <- "GSE197268_Harmony_Integrated.RDS"

## Global parameters
TEST_USE               <- "wilcox"
MIN_PCT                <- 0.01   # keep low: FindMarkers output feeds GSEA
LOGFC_THRESHOLD        <- 0.01   # near 0 for an unbiased ranking
MIN_CELLS_PER_PATIENT  <- 1      # patients below this are dropped from that arm
MIN_PATIENTS_PER_GROUP <- 1      # need >= this many patients per arm to compare
MIN_CELLS_PER_GROUP    <- 50     # need >= this many cells per arm after balancing
SEED                   <- 1234

GSEA_MIN_SIZE <- 5
GSEA_MAX_SIZE <- 500
GSEA_PADJ     <- 0.1

## Pathways of interest, referenced by name everywhere below
PATHWAY_TNFA <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
PATHWAY_IFNA <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"

hallmark_pathways <- gmtPathways(file.path(gmt_dir, "h.all.v2024.1.Hs.symbols.gmt"))
reactome_pathways <- gmtPathways(file.path(gmt_dir, "ReactomePathways.gmt"))
combined_pathways <- c(hallmark_pathways, reactome_pathways)

# Load and split into compartments # --------

ref.data <- readRDS(file.path(ref_dir, f_ref))

## ---- T cells ----
t.cells <- subset(ref.data, subset = cell_type %in% c("CD4 T", "CD8 T", "Infusion T"))
t.cells <- subset(t.cells, subset = response %in% c("NR", "R"))

split_cols        <- str_split_fixed(t.cells$patient_id, "-", 2)
t.cells$patient   <- split_cols[, 1]
t.cells$timepoint <- split_cols[, 2]

# combine all Day 7 cells into one timepoint
t.cells$timepoint[t.cells$timepoint == "D7"] <- "D7-CART"

## ---- Monocytes / DC ----
monocyte.dataset <- subset(ref.data, subset = cell_type %in% c("mDC", "Monocyte", "pDC"))
monocyte.dataset <- subset(monocyte.dataset, subset = response %in% c("NR", "R"))
monocyte.dataset <- subset(monocyte.dataset,
                           subset = timepoint %in% c("D7-CART", "Infusion"), invert = TRUE)

monocyte.dataset$cell_type[monocyte.dataset$cell_type == "Monocyte"] <- "CD14 Monocyte"
monocyte.dataset$cell_type[monocyte.dataset$integrated_clusters %in% c(16, 26)] <- "CD16 Monocyte"
monocyte.dataset$timepoint_fine[monocyte.dataset$timepoint_fine %in% c("D7")] <- "D7-nonCAR"

# same patient column as the T cells, so both arms use `patient`
split_cols                 <- str_split_fixed(monocyte.dataset$patient_id, "-", 2)
monocyte.dataset$patient   <- split_cols[, 1]

monocyte.dataset$response <- factor(monocyte.dataset$response, levels = c("R", "NR"))

rm(ref.data); gc(verbose = FALSE)


# Patient-balanced downsampling helpers # --------


#' Allocate patient cells as evenly as possible.
#' @param avail named integer vector: cells available per patient
#' @param total integer: number of cells to allocate (<= sum(avail))
#' @return named integer vector of per-patient allocations summing to `total`
water_fill <- function(avail, total) {
  stopifnot(total <= sum(avail), total >= 0)
  alloc  <- setNames(rep(0L, length(avail)), names(avail))
  remain <- as.integer(total)
  
  repeat {
    head_room <- avail - alloc
    active    <- head_room > 0
    if (remain <= 0L || !any(active)) break
    
    share <- remain %/% sum(active)
    
    if (share < 1L) {
      ## fewer slots left than active patients: hand out singles to the
      ## patients with the most headroom
      idx <- which(active)
      idx <- idx[order(head_room[idx], decreasing = TRUE)][seq_len(remain)]
      alloc[idx] <- alloc[idx] + 1L
      remain <- 0L
      break
    }
    
    give <- pmin(share, head_room[active])
    alloc[active] <- alloc[active] + give
    remain <- remain - sum(give)
  }
  alloc
}

#' Select a patient-balanced, group-balanced set of cells for one comparison.
#' @return list(cells, summary, n_per_group) or NULL if not runnable
balanced_downsample <- function(md, group_col, patient_col, group1, group2,
                                min_cells_per_patient  = MIN_CELLS_PER_PATIENT,
                                min_patients_per_group = MIN_PATIENTS_PER_GROUP,
                                seed = SEED) {
  
  md <- md[md[[group_col]] %in% c(group1, group2), , drop = FALSE]
  if (nrow(md) == 0) return(NULL)
  md$.cell <- rownames(md)
  md$.grp  <- as.character(md[[group_col]])
  md$.pat  <- as.character(md[[patient_col]])
  
  cnt <- md %>%
    dplyr::count(.grp, .pat, name = "n_avail") %>%
    dplyr::filter(n_avail >= min_cells_per_patient)
  
  a1 <- dplyr::filter(cnt, .grp == group1)
  a2 <- dplyr::filter(cnt, .grp == group2)
  if (nrow(a1) < min_patients_per_group || nrow(a2) < min_patients_per_group)
    return(NULL)
  
  target <- min(sum(a1$n_avail), sum(a2$n_avail))
  
  al1 <- water_fill(setNames(a1$n_avail, a1$.pat), target)
  al2 <- water_fill(setNames(a2$n_avail, a2$.pat), target)
  
  set.seed(seed)
  pick <- function(g, alloc) {
    unlist(lapply(names(alloc), function(p) {
      pool <- md$.cell[md$.grp == g & md$.pat == p]
      if (alloc[[p]] >= length(pool)) pool else sample(pool, alloc[[p]])
    }), use.names = FALSE)
  }
  
  cells1 <- pick(group1, al1)
  cells2 <- pick(group2, al2)
  
  summ <- rbind(
    data.frame(group = group1, patient = names(al1),
               n_avail = a1$n_avail[match(names(al1), a1$.pat)],
               n_used  = as.integer(al1), stringsAsFactors = FALSE),
    data.frame(group = group2, patient = names(al2),
               n_avail = a2$n_avail[match(names(al2), a2$.pat)],
               n_used  = as.integer(al2), stringsAsFactors = FALSE)
  )
  
  list(cells = c(cells1, cells2), summary = summ, n_per_group = target)
}


# DEG driver: celltype x timepoint, R vs NR # --------

# helper # ------
#' @param obj         Seurat object, already subset to the population of interest
#' @param anno_col    cell-type annotation column
#' @param file_tag    e.g. "CARpos_Tisacel_RvsNR_Tcells"
#' @param level_tag   annotation level label, used in filenames
run_response_degs_balanced <- function(obj,
                                       anno_col,
                                       group_col     = "response",
                                       group1        = "R",
                                       group2        = "NR",
                                       patient_col   = "patient",
                                       timepoint_col = "timepoint",
                                       out_dir       = deg_out_dir,
                                       file_tag      = "RvsNR",
                                       level_tag     = "cell_type",
                                       celltypes     = NULL,
                                       timepoints    = NULL,
                                       assay         = "RNA") {
  
  DefaultAssay(obj) <- assay
  if (inherits(obj[[assay]], "Assay5")) obj <- JoinLayers(obj)
  
  md <- obj@meta.data
  if (is.null(celltypes))  celltypes  <- sort(unique(as.character(md[[anno_col]])))
  if (is.null(timepoints)) timepoints <- sort(unique(as.character(md[[timepoint_col]])))
  
  deg_list <- list()
  log_list <- list()
  
  for (ct in celltypes) {
    for (tp in timepoints) {
      
      tag <- paste0(gsub("[/ ]", "", ct), "_", gsub("[/ ]", "", tp))
      sub_md <- md[as.character(md[[anno_col]])        == ct &
                     as.character(md[[timepoint_col]]) == tp, , drop = FALSE]
      
      if (nrow(sub_md) == 0) {
        message("skip ", tag, ": no cells"); next
      }
      
      bal <- balanced_downsample(sub_md, group_col, patient_col, group1, group2)
      if (is.null(bal)) {
        message("skip ", tag, ": too few patients / cells per arm"); next
      }
      if (bal$n_per_group < MIN_CELLS_PER_GROUP) {
        message("skip ", tag, ": only ", bal$n_per_group, " cells per arm"); next
      }
      
      message("running ", tag, ": ", bal$n_per_group, " cells per arm across ",
              length(unique(bal$summary$patient)), " patients")
      
      bal$summary$celltype  <- ct
      bal$summary$timepoint <- tp
      log_list[[tag]]       <- bal$summary
      
      sub_obj <- subset(obj, cells = bal$cells)
      Idents(sub_obj) <- group_col
      
      degs <- FindMarkers(object          = sub_obj,
                          ident.1         = group1,
                          ident.2         = group2,
                          test.use        = TEST_USE,
                          min.pct         = MIN_PCT,
                          logfc.threshold = LOGFC_THRESHOLD,
                          assay           = assay,
                          verbose         = FALSE)
      
      degs <- degs %>%
        tibble::rownames_to_column("geneID") %>%
        mutate(celltype     = ct,
               timepoint    = tp,
               comparison   = paste0(group1, "_vs_", group2),
               n_per_group  = bal$n_per_group,
               n_patients_1 = sum(bal$summary$group == group1),
               n_patients_2 = sum(bal$summary$group == group2))
      
      write.csv(degs,
                file.path(out_dir, paste0(tag, "_ALLFCs_", file_tag,
                                          "_patientbalanced_", level_tag, ".csv")),
                row.names = FALSE)
      
      deg_list[[tag]] <- degs
    }
  }
  
  all_degs <- if (length(deg_list)) bind_rows(deg_list) else NULL
  all_log  <- if (length(log_list)) bind_rows(log_list)  else NULL
  
  if (!is.null(all_degs))
    write.csv(all_degs,
              file.path(out_dir, paste0("COMBINED_", file_tag, "_", level_tag, ".csv")),
              row.names = FALSE)
  if (!is.null(all_log))
    write.csv(all_log,
              file.path(out_dir, paste0("DOWNSAMPLING_LOG_", file_tag, "_", level_tag, ".csv")),
              row.names = FALSE)
  
  list(degs = all_degs, log = all_log)
}


# HELPER: fgsea per celltype x timepoint # --------


#' Build a ranking vector from a DEG table.
make_rank <- function(df, rank_by = "log2FC") {
  r <- switch(rank_by,
              log2FC   = df$avg_log2FC,
              signed_p = sign(df$avg_log2FC) * -log10(pmax(df$p_val, .Machine$double.xmin)),
              stop("unknown rank_by"))
  names(r) <- df$geneID
  r <- r[is.finite(r)]
  r <- r[!duplicated(names(r))]
  sort(r, decreasing = TRUE)
}

#' Run fgsea for every celltype x timepoint present in `deg_df`.
#' @param sig_col "padj" (within collection) or "padj_combined" (pooled)
#' @return one combined data.table of fgsea results (all pathways, unfiltered)
run_fgsea_by_celltype_timepoint <- function(deg_df, tag, image = FALSE,
                                            rank_by = "log2FC",
                                            sig_col = "padj") {
  
  sig_col <- match.arg(sig_col, c("padj", "padj_combined"))
  
  combos <- deg_df %>% distinct(celltype, timepoint)
  res_list <- list()
  
  for (i in seq_len(nrow(combos))) {
    ct <- combos$celltype[i]; tp <- combos$timepoint[i]
    key <- paste0(gsub("[/ ]", "", ct), "_", gsub("[/ ]", "", tp))
    message("fgsea: ", tag, " | ", key)
    
    d <- deg_df %>%
      filter(celltype == ct, timepoint == tp) %>%
      mutate(log10_padj = pmin(-log10(p_val_adj), 310))
    
    d <- d[
      !grepl("^MT", d$geneID, ignore.case = TRUE) &
        !grepl("^LINC", d$geneID, ignore.case = TRUE),
    ]
    
    d <- d %>% dplyr::filter(p_val_adj < 0.01)
    print(dim(d))
    
    rank <- make_rank(d, rank_by)
    if (length(rank) < 100) { message("  too few genes, skipping"); next }
    
    set.seed(SEED)
    hallmark.results <- fgsea(pathways = hallmark_pathways,
                              stats    = rank,
                              minSize  = GSEA_MIN_SIZE,
                              maxSize  = GSEA_MAX_SIZE,
                              nproc    = 1)
    set.seed(SEED)
    reactome.results <- fgsea(pathways = reactome_pathways,
                              stats    = rank,
                              minSize  = GSEA_MIN_SIZE,
                              maxSize  = GSEA_MAX_SIZE,
                              nproc    = 1)
    
    if (nrow(hallmark.results)) hallmark.results[, collection := "HALLMARK"]
    if (nrow(reactome.results)) reactome.results[, collection := "REACTOME"]
    
    res <- rbind(hallmark.results, reactome.results, fill = TRUE)
    if (!nrow(res)) next
    
    res[, padj_combined := p.adjust(pval, method = "BH")]
    res[, `:=`(celltype = ct, timepoint = tp)]
    setcolorder(res, c("collection", "pathway", "pval", "padj", "padj_combined"))
    
    res_list[[key]] <- res
    
    sig <- res[get(sig_col) < GSEA_PADJ]
    
    ## write a flat CSV (leadingEdge collapsed to a string)
    out <- copy(res)
    out[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
    fwrite(out[order(get(sig_col))],
           file.path(gsea_out_dir, paste0(tag, "_", key, "_fgsea.csv")))
    
    if (image && nrow(sig)) {
      up   <- sig[ES > 0][head(order(pval), 10), pathway]
      down <- sig[ES < 0][head(order(pval), 10), pathway]
      top  <- c(up, rev(down))
      if (length(top)) {
        sub_dir <- file.path(gsea_table_dir, tag)
        dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
        pdf(file.path(sub_dir, paste0(key, ".pdf")), height = 12, width = 15)
        print(plotGseaTable(combined_pathways[top], rank, sig,
                            gseaParam         = 0.5,
                            pathwayLabelStyle = list(size = 10),
                            headerLabelStyle  = list(color = "red"),
                            valueStyle        = list(size = 10),
                            axisLabelStyle    = list(size = 6)))
        dev.off()
      }
    }
  }
  
  if (!length(res_list)) return(NULL)
  all_res <- rbindlist(res_list)
  out <- copy(all_res)
  out[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
  fwrite(out, file.path(gsea_out_dir, paste0("COMBINED_", tag, "_fgsea.csv")))
  all_res
}


# Run on CAR+ Tisa-cel T cells and Tisa-cel monocytes # --------


## CAR+ Tisa-cel only. `CAR` is logical, `generic` is the product.
DefaultAssay(t.cells) <- "RNA"
keep <- Cells(t.cells)[which(as.character(t.cells$CAR)     == "TRUE" &
                               as.character(t.cells$generic) == "Tisa-cel")]
t.car.tisa <- subset(t.cells, cells = keep)

message("T cells kept: ", ncol(t.car.tisa))
print(table(t.car.tisa$cell_type, t.car.tisa$timepoint))
print(table(t.car.tisa$response,  t.car.tisa$timepoint))

t_deg <- run_response_degs_balanced(
  obj         = t.car.tisa,
  anno_col    = "cell_type",
  group_col   = "response",
  group1      = "R", group2 = "NR",
  patient_col = "patient",
  file_tag    = "CARpos_Tisacel_RvsNR_Tcells",
  level_tag   = "cell_type"
)

t_gsea <- run_fgsea_by_celltype_timepoint(
  deg_df = t_deg$degs,
  tag    = "CARpos_Tisacel_Tcells",
  image  = TRUE
)

mono.tisa <- subset(monocyte.dataset, subset = generic == "Tisa-cel")

message("Monocytes kept: ", ncol(mono.tisa))
print(table(mono.tisa$cell_type, mono.tisa$timepoint))

mono_deg <- run_response_degs_balanced(
  obj         = mono.tisa,
  anno_col    = "cell_type",
  group_col   = "response",
  group1      = "R", group2 = "NR",
  patient_col = "patient",
  file_tag    = "Tisacel_RvsNR_Monocytes",
  level_tag   = "cell_type"
)

mono_gsea <- run_fgsea_by_celltype_timepoint(
  deg_df = mono_deg$degs,
  tag    = "Tisacel_Monocytes",
  image  = TRUE
)


#  Overlap of leading-edge genes # --------


## ---- TNFA signaling: External  Day 7 vs This Dataset CAR+ CD8 EM Peak ----
# gene lists obtained from fgsea # ------
external_day7_tnfa <- c("FOS", "SOCS3", "DUSP1", "TUBB2A", "NR4A3", "PPP1R15A", "NR4A2",
                   "BTG2", "TNFSF9", "DUSP2", "KLF6", "CD69", "PHLDA1", "JUNB", "ID2",
                   "FOSB", "SAT1", "JUN", "GADD45B", "TNFAIP3", "ZFP36", "PNRC1",
                   "BTG1", "KLF2", "IER2", "EIF1", "LITAF", "CCNL1", "CD44")

this_dataset_tnfa <- c("PER1", "NR4A2", "PFKFB3", "ZBTB10", "CDKN1A", "DUSP1", "FOSL2",
              "PDE4B", "DUSP4", "IL23A", "CD83", "KDM6B", "CD69", "JUNB", "NR4A3",
              "JUN", "BTG2", "DUSP2", "MAP3K8", "ZFP36", "REL", "TNFAIP3", "ZC3H12A",
              "NINJ1", "BCL3", "FOSB", "IER5", "KLF10", "BTG1", "NFKBIE", "BTG3",
              "RELB", "PTGER4", "NFKBIA", "NFKB1", "PNRC1", "GADD45B", "EIF1",
              "PPP1R15A", "CEBPB", "BHLHE40", "CCNL1", "B4GALT1", "IER2", "DUSP5",
              "KLF6")

shared <- intersect(external_day7_tnfa, this_dataset_tnfa)
p_tnfa_overlap <- phyper(length(shared) - 1, length(external_day7_tnfa),
                         20000 - length(external_day7_tnfa), length(this_dataset_tnfa),
                         lower.tail = FALSE)

pdf(file.path(base_out, "VennDiagram_TNFA_overlap_day7_external_vs_this_dataset.pdf"),
    width = 6, height = 4)
print(
  ggvenn(list(reference_day7 = external_day7_tnfa, this_dataset_cd8em = this_dataset_tnfa),
         fill_color = c("#7F77DD", "#1D9E75")) +
    ggtitle(paste0(PATHWAY_TNFA, "\nhypergeometric overlap p=", signif(p_tnfa_overlap, 3)))
)
dev.off()

## ---- Interferon alpha: External Day 7 vs This Dataset CD14 monocytes ----
external_day7_ifna <- c("RSAD2", "IFI44L", "IFIT3", "IFITM1", "ISG15", "CMPK2", "IFITM3",
                   "IFIT2", "USP18", "IFI44", "EPSTI1", "DHX58", "OASL", "IFI27",
                   "MX1", "SAMD9L", "LY6E", "STAT2", "DDX60", "IRF7", "NCOA7", "OAS1",
                   "HELZ2", "IFI35", "UBE2L6", "PSME2", "BST2", "PSMB9", "GBP2",
                   "TAP1", "PARP14", "HLA-C", "LAP3", "SP110", "CASP1", "ADAR",
                   "PSMB8", "OGFR")

this_dataset_ifna <- c("MX1", "RSAD2", "IFITM1", "IFIT3", "IFITM3", "IFIT2", "ISG20", "EPSTI1",
              "IFI44L", "ISG15", "PSMB9", "SAMD9", "OASL", "CMPK2", "USP18", "OAS1",
              "LY6E", "IFI44", "SAMD9L", "LAP3", "IRF7", "DDX60", "GMPR", "RTP4",
              "PSME2", "IFIH1", "UBE2L6", "IFITM2", "IFI35", "NMI", "PLSCR1", "TRIM5",
              "TRIM21", "SELL", "PARP9", "EIF2AK2", "PNPT1", "DHX58", "IRF2", "PSMB8",
              "BST2", "GBP4", "IRF1", "B2M", "SP110", "TAP1", "TRAFD1", "PSMA3",
              "GBP2", "PSME1", "PARP14", "STAT2", "IFI30", "CASP1", "TRIM14", "HELZ2",
              "LPAR6", "IRF9", "PARP12", "NUB1")

shared <- intersect(external_day7_ifna, this_dataset_ifna)
p_ifna_overlap <- phyper(length(shared) - 1, length(external_day7_ifna),
                         20000 - length(external_day7_ifna), length(this_dataset_ifna),
                         lower.tail = FALSE)

pdf(file.path(base_out, "VennDiagram_IFNA_CD14mono_overlap_day7_external_vs_this_dataset.pdf"),
    width = 6, height = 4)
print(
  ggvenn(list(external_day7 = external_day7_ifna, this_dataset_cd14mono = this_dataset_ifna),
         fill_color = c("#7F77DD", "#1D9E75")) +
    ggtitle(paste0(PATHWAY_IFNA, "\nhypergeometric overlap p=", signif(p_ifna_overlap, 3)))
)
dev.off()


# GSEA Enrichment Plots 

#' Enrichment plot for one pathway in one celltype x timepoint.
#' @param gsea_df optional; used only to annotate NES / padj on the plot
#' @param rank_by must match what was used in the fgsea run
plot_enrichment <- function(deg_df, gsea_df = NULL, pathway_name, cell, tp,
                            sig_col = "padj", rank_by = "log2FC",
                            subtitle_extra = NULL) {
  
  stopifnot(pathway_name %in% names(combined_pathways))
  
  d <- deg_df %>% filter(celltype == cell, timepoint == tp)
  if (!nrow(d)) stop("no DEG rows for ", cell, " / ", tp)
  rank <- make_rank(d, rank_by)
  
  ## stats from the fgsea run, if supplied
  lab <- NULL
  if (!is.null(gsea_df)) {
    s <- as.data.frame(gsea_df) %>%
      filter(pathway == pathway_name, celltype == cell, timepoint == tp)
    if (nrow(s) == 1)
      lab <- sprintf("NES = %.2f    %s = %.3g    size = %d",
                     s$NES, sig_col, s[[sig_col]], s$size)
  }
  
  plotEnrichment(combined_pathways[[pathway_name]], rank) +
    labs(
      title    = gsub("_", " ", pathway_name),
      subtitle = paste(c(paste0(cell, " | ", tp), lab, subtitle_extra),
                       collapse = "\n"),
      x = "Gene rank", y = "Enrichment score"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(size = 10, colour = "grey25"),
      panel.grid.minor = element_blank()
    )
}

p_tnfa <- plot_enrichment(
  deg_df       = t_deg$degs,
  gsea_df      = t_gsea,
  pathway_name = PATHWAY_TNFA,
  cell         = "CD8 T",
  tp           = "D7-CART",
  subtitle_extra = "CAR+ Tisa-cel, R vs NR"
)
ggsave(file.path(gsea_out_dir, "enrichment_CD8T_D7_TNFA_NFKB.pdf"),
       p_tnfa, width = 6, height = 4)

p_ifna <- plot_enrichment(
  deg_df       = mono_deg$degs,
  gsea_df      = mono_gsea,
  pathway_name = PATHWAY_IFNA,
  cell         = "CD14 Monocyte",
  tp           = "D7",
  subtitle_extra = "Tisa-cel, R vs NR"
)
ggsave(file.path(gsea_out_dir, "enrichment_CD14Mono_D7_IFNA.pdf"),
       p_ifna, width = 6, height = 4)

print(p_tnfa)
print(p_ifna)
