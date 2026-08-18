################################################################################
## Lymphoma CAR-T single-cell RNA-seq: differential expression analyses
##
## Run 00_setup.R first. This script assumes the following are already loaded:
##   - Seurat objects: full.dataset, tnk.car, tnk.seurat, mono.object
##   - Helper function: check_downsample()
##
################################################################################

# ---- CONFIG -----------------------------------------------------------------
# Output directory for DEG result CSVs. Edit to match your environment.
results_dir <- file.path("results", "deg")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Tag appended to output filenames (set to "" to omit).
date_tag <- format(Sys.Date(), "%Y%m%d")

# Shared FindMarkers / FindAllMarkers parameters (0.01 thresholds = effectively
# unfiltered, so all fold-changes are retained for downstream filtering).
LOGFC_THRESHOLD <- 0.01
MIN_PCT         <- 0.01
TEST_USE        <- "wilcox"

################################################################################
## Helper functions
################################################################################

#' Pairwise DEGs per cell type, per timepoint, between two grouping levels.
#'
#' Builds an `idents.compare` column of the form
#' "<celltype>_<group>_<timepoint>" and, for each cell type and timepoint,
#' compares `group1` vs `group2` with FindMarkers. Uses check_downsample()
#' to skip underpowered comparisons and to balance large/imbalanced ones by
#' downsampling both groups to the smaller group's size. Writes one CSV per
#' comparison into `out_dir`.
#'
#' @param obj        Seurat object.
#' @param anno_col   Metadata column with cell-type annotations.
#' @param split_col  Metadata column defining the two comparison groups.
#' @param group1     Value of split_col for ident.1 (numerator).
#' @param group2     Value of split_col for ident.2 (denominator).
#' @param out_dir    Directory for output CSVs.
#' @param file_tag   Comparison label embedded in filenames (e.g. "respondervsnon").
#' @param level_tag  Annotation-level label embedded in filenames (e.g. "level2").
#' @param celltypes  Optional vector of cell types; defaults to all in anno_col.
run_pairwise_timepoint_degs <- function(obj, anno_col, split_col,
                                        group1, group2, out_dir, file_tag,
                                        celltypes = NULL) {
  
  obj$idents.compare <- paste0(obj[[anno_col]][, 1], "_",
                               obj[[split_col]][, 1], "_",
                               obj$timepoint)
  Idents(obj) <- "idents.compare"
  
  if (is.null(celltypes)) celltypes <- unique(obj[[anno_col]][, 1])
  timepoints <- unique(obj$timepoint)
  
  for (c in celltypes) {
    for (t in timepoints) {
      ID.1 <- paste0(c, "_", group1, "_", t)
      ID.2 <- paste0(c, "_", group2, "_", t)
      
      downsample.check <- check_downsample(obj, ID.1, ID.2)
      if (identical(downsample.check, "skip")) {
        print("skip this comparison")
        next
      }
      
      subset.obj <- obj
      if (isTRUE(downsample.check)) {
        smallest <- min(sum(obj$idents.compare == ID.1),
                        sum(obj$idents.compare == ID.2))
        subset.obj <- subset(obj, downsample = round(smallest))
        print("subsetting to smallest")
      }
      
      marker.degs <- FindMarkers(
        object          = subset.obj,
        ident.1         = ID.1,
        ident.2         = ID.2,
        test.use        = TEST_USE,
        min.pct         = MIN_PCT,
        logfc.threshold = LOGFC_THRESHOLD,
        verbose         = TRUE,
        assay           = "RNA"
      )
      
      c.sub <- gsub("/", "", x = c)
      write.csv(
        marker.degs,
        file.path(out_dir, paste0(
          c.sub, "_", t, "_ALLFCs_", file_tag,
          "_downsample_", downsample.check,
          "_", date_tag, ".csv")),
        row.names = TRUE
      )
    }
  }
}

#' FindAllMarkers with the project's standard parameters.
#'
#' @param obj       Seurat object.
#' @param ident_col Metadata column to set as identity before testing.
#' @param out_file  Path to write the resulting marker table.
#' @return The marker data frame (also written to disk).
find_all_markers_std <- function(obj, ident_col, out_file) {
  Idents(obj) <- ident_col
  markers <- FindAllMarkers(
    obj,
    assay             = "RNA",
    logfc.threshold   = LOGFC_THRESHOLD,
    test.use          = TEST_USE,
    slot              = "data",
    min.pct           = MIN_PCT,
    min.diff.pct      = 0.01,
    verbose           = TRUE,
    only.pos          = FALSE,
    max.cells.per.ident = Inf,
    random.seed       = 1,
    latent.vars       = NULL,
    min.cells.feature = 3,
    min.cells.group   = 3,
    mean.fxn          = NULL
  )
  write.csv(markers, out_file)
  markers
}

################################################################################
## 1. Responder vs. non-responder DEGs (per timepoint)
################################################################################

# Derive the non-CAR (CAR-negative) subset.
all.cells   <- Cells(full.dataset)
car.only    <- Cells(tnk.car)
noncar.only <- setdiff(all.cells, car.only)
noncar.obj  <- subset(full.dataset, cells = noncar.only)

# Non-CAR, level-2 annotations.
run_pairwise_timepoint_degs(
  obj       = noncar.obj,
  anno_col  = "cell.anno",
  split_col = "response",
  group1    = "CR",
  group2    = "PD",
  out_dir   = file.path(results_dir, "noncar_all"),
  file_tag  = "CRvsPD"
)


# CAR, level-2 annotations.
run_pairwise_timepoint_degs(
  obj       = tnk.car,
  anno_col  = "cell.anno",
  split_col = "response",
  group1    = "CR",
  group2    = "PD",
  out_dir   = file.path(results_dir, "car_all"),
  file_tag  = "CRvsPD"
)



################################################################################
## 2. FindAllMarkers across cell types
################################################################################

# T/NK compartment, level-2 annotations.
tnk.markers <- find_all_markers_std(
  tnk.seurat, "cell.anno",
  file.path(results_dir, "TNK_AllTime_FindAllMarkers.csv")
)

# Monocyte compartment, level-1 and level-2 annotations.
mono.markers <- find_all_markers_std(
  mono.object,  "cell.anno",
  file.path(results_dir, "Mono_AllTime_FindAllMarkers.csv")
)
