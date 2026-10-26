## =============================================================================
## 1) Score T cells with IFN / TNFa modules (AddModuleScore)
## 2) Correlate CD8 TEM and myeloid module scores at each timepoint (patient level)
##
## Myeloid scores come from the myeloid scoring script:
##   <OUT_ROOT>/tables/mono_module_scores.csv
##
## Cell type is always read from `cell.anno`.
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

set.seed(1)
options(stringsAsFactors = FALSE)

## -----------------------------------------------------------------------------
## 0. Config
## -----------------------------------------------------------------------------

PROJ_DIR <- Sys.getenv("PROJ_DIR", unset = ".")
OUT_ROOT <- Sys.getenv("OUT_DIR",  unset = file.path(PROJ_DIR, "signature_scoring"))


PATHS <- list(
  full_rds = file.path(PROJ_DIR, "seurat_objects", "LBCL_full_dataset_clean_new.RDS")  # full clean object (metadata included)
)

# T/NK cell types (values of cell.anno) kept for scoring
TNK_TYPES <- c("CD4 Naive", "CD4 TEM-like", "Treg",
               "CD8 Naive", "CD8 TEM", 
               "gdT", "MAIT", "NK", "Proliferating")

OUT <- list(
  tables = file.path(OUT_ROOT, "tables"),
  corr   = file.path(OUT_ROOT, "tcell_myeloid_correlation")
)
invisible(lapply(OUT, dir.create, showWarnings = FALSE, recursive = TRUE))

PAL <- c(CR = "#4FB7C5", PD = "#E7A75E")

## -----------------------------------------------------------------------------
## Utilities
## -----------------------------------------------------------------------------

log_msg <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

check_inputs <- function(...) {
  f <- unlist(list(...))
  missing <- f[!file.exists(f)]
  if (length(missing)) {
    stop("Missing required input file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}

## =============================================================================
## PART 1. T cell scoring
## =============================================================================

## ---- modules ----------------------------------------------------------------
## Gene lists for the two T cell modules. Part 2 uses both.
car_tnfa <- c(
  "PER1",     "NR4A2",    "PFKFB3",   "ZBTB10",   "CDKN1A",   "DUSP1",    "FOSL2",    "PDE4B",    "DUSP4",    "IL23A",    "CD83",     "KDM6B",   
  "CD69",     "JUNB",     "NR4A3",    "JUN",      "BTG2",     "DUSP2",    "MAP3K8",   "ZFP36",    "REL",      "TNFAIP3",  "ZC3H12A",  "NINJ1",   
  "BCL3",     "FOSB",     "IER5",     "KLF10",    "BTG1",     "NFKBIE",   "BTG3",     "RELB",     "PTGER4",   "NFKBIA",   "NFKB1",    "PNRC1",   
  "GADD45B",  "EIF1",     "PPP1R15A", "CEBPB",    "BHLHE40",  "CCNL1",    "B4GALT1",  "IER2",     "DUSP5",    "KLF6"
)
car_ifn_union <- c(
  "IFIT3",    "MX1",      "EPSTI1",   "RSAD2",    "CMPK2",    "IFITM3",   "USP18",    "PARP9",    "IRF7",     "OAS1",     "TRAFD1",   "TRIM21",  
  "IFI35",    "OASL",     "BST2",     "DDX60",    "PSMB9",    "SAMD9L",   "SAMD9",    "RTP4",     "IFI30",    "ISG15",    "IFITM1",   "TAP1",    
  "UBE2L6",   "PSMB8",    "IFIT1",    "MX2",      "STAT1",    "CD38",     "OAS3",     "OAS2",     "TNFSF10",  "XAF1",     "HLA-DRB1", "GBP1",    
  "FANCF",    "IFI6",     "FKBP5",    "CDK1",     "IRF2",     "GBP5",     "PTPN6",    "HLA-DRA"
)
t_modules <- list(car_ifn_union = car_ifn_union, car_tnfa = car_tnfa)
stopifnot(lengths(t_modules) > 0)   # fill in both gene lists above

check_inputs(PATHS$metadata_csv, PATHS$tnk_rds)


## ---- scoring function -------------------------------------------------------
## Adds a <module>_AMS column to obj@meta.data for each module
score_object <- function(obj, modules, assay = "RNA") {
  
  DefaultAssay(obj) <- assay
  
  log_msg("  AddModuleScore scoring (", length(modules), " modules)")
  for (sig in names(modules)) {
    tmp <- paste0(sig, "_AMS")
    obj <- AddModuleScore(obj, features = modules[sig], name = tmp, seed = 1)
    seurat_col <- paste0(tmp, "1")          # Seurat always appends the index
    obj@meta.data[[tmp]]        <- obj@meta.data[[seurat_col]]
    obj@meta.data[[seurat_col]] <- NULL
  }
  obj
}

## ---- score the T/NK object --------------------------------------------------
log_msg("=== PART 1: T cell scoring ===")

full_obj <- readRDS(PATHS$full_rds)

req_cols <- c("sample.name", "patient.id", "timepoint", "response", "CAR.exp", "cell.anno")
miss_cols <- setdiff(req_cols, colnames(full_obj@meta.data))
if (length(miss_cols)) stop("Metadata is missing column(s): ", paste(miss_cols, collapse = ", "))

# subset T/NK cells using the object's own metadata; drop the QC-failed sample
meta     <- full_obj@meta.data
tnk_keep <- rownames(meta)[meta$cell.anno %in% TNK_TYPES]
log_msg("  keeping ", length(tnk_keep), " T/NK cells of ", ncol(full_obj))
print(table(meta[tnk_keep, "cell.anno"]))

obj <- subset(full_obj, cells = tnk_keep)
rm(full_obj, meta); gc(verbose = FALSE)

obj <- score_object(obj, t_modules)

write.csv(obj@meta.data, file.path(OUT$tables, "tnk_module_scores.csv"))
rm(obj); gc(verbose = FALSE)


## =============================================================================
## PART 2. CD8 TEM vs myeloid score correlation at each timepoint
## =============================================================================
##
## Reuses the saved per-cell scores
##   tables/tnk_module_scores.csv   (part 1)
##   tables/mono_module_scores.csv  (myeloid scoring script)
##
## All scores are AddModuleScore (AMS).
##   pathway  T-cell score          myeloid score
##   IFN      car_ifn_union_AMS     <celltype>_<timepoint>_lg_union_AMS
##   TNFa     car_tnfa_AMS          <celltype>_peak_TNFa_AMS
##
## CD8 TEM by timepoint:  APH, D0 -> CAR- only;  TDN (CAR+ product) -> paired with
## D0 myeloid cells;  Peak, 4W -> CAR- and CAR+ separately.
## Myeloid groups: CD14 Mono, CD16 Mono, DC, and all three pooled.
## Myeloid IFN module per timepoint: APH/D0 -> day0, Peak -> peak,
## 4W -> week4 (falls back to peak when no week4 module).
## Myeloid TNFa: only CD14 and DC peak modules exist; CD16 uses the CD14 one.
##
## Tests: Spearman and Pearson on patient means (BH applied separately per method).
## Primary comparison: CD8 TEM vs All myeloid.
## =============================================================================

log_msg("=== PART 2: CD8 TEM vs myeloid correlation ===")

## ---- config -----------------------------------------------------------------
check_inputs(file.path(OUT$tables, "tnk_module_scores.csv"),
             file.path(OUT$tables, "mono_module_scores.csv"))

T_TYPE_COL <- "cell.anno"                                # T-cell annotation column
T_TYPES    <- c("CD8 TEM")
M_TYPE_COL <- "cell.anno"                                # myeloid annotation column
CAR_COL    <- "CAR.exp"                                  # TRUE / FALSE

M_LABEL <- c("CD14 Mono" = "CD14 Mono", "CD16 Mono" = "CD16 Mono", "Dendritic Cell" = "DC")
M_CODE  <- c("CD14 Mono" = "CD14", "CD16 Mono" = "CD16", "Dendritic Cell" = "DC")
TNF_CT  <- c(CD14 = "CD14", CD16 = "CD14", DC = "DC")    # CD16 has no TNFa module

TP_MAP   <- c(APH = "day0", D0 = "day0", Peak = "peak", `4W` = "week4")
TPS      <- c("APH", "D0", "Peak", "4W")                 # myeloid timepoints
PLOT_TPS <- c("APH", "D0", "TDN", "Peak", "4W")          # plot columns

## which CAR status is used for T cells at each timepoint
## TDN = CAR+ infusion product: its own column, paired with the patient's D0 myeloid cells
T_STATUS <- tibble(timepoint  = c("APH",  "D0",   "TDN",  "Peak", "Peak", "4W",   "4W"),
                   car_status = c("CAR-", "CAR-", "CAR+", "CAR-", "CAR+", "CAR-", "CAR+"))
PRODUCT_TP <- "TDN"

PRIMARY_T <- "CD8 TEM"
PRIMARY_M <- "All myeloid"
MIN_N     <- 4                         # min patients per correlation

## ---- T cells ----------------------------------------------------------------
tc <- read_csv(file.path(OUT$tables, "tnk_module_scores.csv"), show_col_types = FALSE,
               col_select = all_of(c("patient.id", "timepoint", "response", CAR_COL,
                                     T_TYPE_COL, "car_ifn_union_AMS", "car_tnfa_AMS")))

message("TNK timepoints: ", paste(unique(tc$timepoint), collapse = ", "))

tc <- tc %>%
  mutate(car_status = ifelse(as.logical(.data[[CAR_COL]]), "CAR+", "CAR-")) %>%
  filter(.data[[T_TYPE_COL]] %in% T_TYPES, response %in% names(PAL)) %>%
  semi_join(T_STATUS, by = c("timepoint", "car_status")) %>%
  transmute(patient.id, response, timepoint, car_status, t_group = .data[[T_TYPE_COL]],
            IFN = car_ifn_union_AMS, TNFa = car_tnfa_AMS)

t_pt <- tc %>%
  pivot_longer(c(IFN, TNFa), names_to = "pathway", values_to = "score") %>%
  group_by(patient.id, response, timepoint, car_status, pathway, t_group) %>%
  summarise(t_mean = mean(score, na.rm = TRUE), t_cells = n(), .groups = "drop")
rm(tc)

## ---- myeloid ----------------------------------------------------------------
my <- read_csv(file.path(OUT$tables, "mono_module_scores.csv"), show_col_types = FALSE,
               col_select = c(patient.id, timepoint, response, all_of(M_TYPE_COL),
                              matches("_lg_union_AMS$|_TNFa_AMS$"))) %>%
  filter(.data[[M_TYPE_COL]] %in% names(M_CODE), timepoint %in% TPS, response %in% names(PAL))

# pick the matching IFN / TNFa module column for each cell
ct  <- unname(M_CODE[my[[M_TYPE_COL]]])
ifn <- paste0(ct, "_", unname(TP_MAP[my$timepoint]), "_lg_union_AMS")
ifn[!ifn %in% names(my)] <- sub("_week4_", "_peak_", ifn[!ifn %in% names(my)])
tnf <- paste0(unname(TNF_CT[ct]), "_peak_TNFa_AMS")

M    <- as.matrix(select(my, matches("_AMS$")))
getv <- function(cols) M[cbind(seq_len(nrow(M)), match(cols, colnames(M)))]

my <- tibble(patient.id = my$patient.id, response = my$response, timepoint = my$timepoint,
             m_group = unname(M_LABEL[my[[M_TYPE_COL]]]),
             IFN = getv(ifn), TNFa = getv(tnf))
rm(M)
my <- bind_rows(my, mutate(my, m_group = "All myeloid"))

m_pt <- my %>%
  pivot_longer(c(IFN, TNFa), names_to = "pathway", values_to = "score") %>%
  filter(!is.na(score)) %>%
  group_by(patient.id, response, timepoint, pathway, m_group) %>%
  summarise(m_mean = mean(score), m_cells = n(), .groups = "drop")
rm(my)

## ---- helper -----------------------------------------------------------------
## Spearman and Pearson; `rho` holds Spearman rho or Pearson r depending on `method`
cor2 <- function(x, y) {
  ok <- complete.cases(x, y)
  n  <- sum(ok)
  one <- function(m) {
    if (n < MIN_N) return(tibble(method = m, n = n, rho = NA_real_, p = NA_real_))
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = m))
    tibble(method = m, n = n, rho = unname(ct$estimate), p = ct$p.value)
  }
  bind_rows(one("spearman"), one("pearson"))
}

## ---- correlation per timepoint (patient mean) -------------------------------
## T-cell timepoint is kept as the column label; TDN is matched to D0 myeloid
pairs <- t_pt %>%
  mutate(m_timepoint = ifelse(timepoint == PRODUCT_TP, "D0", timepoint)) %>%
  inner_join(rename(m_pt, m_timepoint = timepoint),
             by = c("patient.id", "response", "m_timepoint", "pathway"),
             relationship = "many-to-many")
write_csv(pairs, file.path(OUT$corr, "patient_scores_paired.csv"))

message("Patients with paired CD8 TEM + all-myeloid data:")
print(pairs %>% filter(t_group == PRIMARY_T, m_group == PRIMARY_M, pathway == "IFN") %>%
        count(timepoint, car_status))

cor_tp <- pairs %>%
  group_by(pathway, t_group, car_status, m_group, timepoint) %>%
  group_modify(~ cor2(.x$t_mean, .x$m_mean)) %>%
  ungroup()

## ---- results table ----------------------------------------------------------
cor_all <- cor_tp %>%
  mutate(timepoint = factor(timepoint, levels = PLOT_TPS),
         t_group   = factor(t_group, levels = T_TYPES),
         m_group   = factor(m_group, levels = c(unname(M_LABEL), "All myeloid")),
         primary   = t_group == PRIMARY_T & m_group == PRIMARY_M) %>%
  group_by(method) %>%                                  # BH done separately per method
  mutate(padj_all = p.adjust(p, "BH")) %>%
  group_by(method, primary) %>%
  mutate(padj_within_set = p.adjust(p, "BH")) %>%       # primary tests corrected among themselves
  ungroup() %>%
  arrange(desc(method), desc(primary), pathway, timepoint, car_status, t_group, m_group)
write_csv(cor_all, file.path(OUT$corr, "correlation_CD8TEM_vs_myeloid.csv"))

## ---- plots ------------------------------------------------------------------
theme_set(theme_classic(base_size = 11) +
            theme(plot.title = element_text(face = "bold"),
                  strip.text = element_text(face = "bold"),
                  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)))

fmt_p <- function(p) ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p))

## primary: CD8 TEM vs all myeloid, rows = CAR status, columns = timepoint
prim <- pairs %>%
  filter(t_group == PRIMARY_T, m_group == PRIMARY_M) %>%
  transmute(patient.id, response, pathway, car_status, timepoint, x = m_mean, y = t_mean) %>%
  mutate(timepoint  = factor(timepoint, levels = PLOT_TPS),
         car_status = factor(car_status, levels = c("CAR-", "CAR+")),
         response   = factor(response, levels = names(PAL)))

## panel labels: Spearman and Pearson, each with unadjusted p and
## BH p adjusted within the primary set
prim_lab <- cor_all %>%
  filter(primary) %>%
  mutate(txt = ifelse(is.na(rho),
                      paste0(ifelse(method == "spearman", "Spearman", "Pearson"), ": NA"),
                      sprintf("%s = %.2f, p = %s, padj = %s",
                              ifelse(method == "spearman", "Spearman rho", "Pearson r"),
                              rho, fmt_p(p), fmt_p(padj_within_set)))) %>%
  arrange(desc(method)) %>%                         # Spearman line first
  group_by(pathway, timepoint, car_status) %>%
  summarise(lab = paste(c(paste0("n = ", n[1]), txt), collapse = "\n"), .groups = "drop") %>%
  mutate(car_status = factor(car_status, levels = c("CAR-", "CAR+")))

plot_prim <- function(pw) {
  ggplot(filter(prim, pathway == pw), aes(x, y)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                color = "black", linewidth = 0.6, alpha = 0.15) +
    geom_point(aes(fill = response), shape = 21, size = 2.3, stroke = 0.25) +
    geom_text_repel(aes(label = patient.id), size = 2, max.overlaps = 50,
                    segment.color = "grey60") +
    geom_text(data = filter(prim_lab, pathway == pw), aes(x = -Inf, y = Inf, label = lab),
              hjust = -0.03, vjust = 1.1, size = 2.5, inherit.aes = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.4))) +
    scale_fill_manual(values = PAL) +
    facet_grid(car_status ~ timepoint, scales = "free", drop = FALSE) +
    labs(title = paste0(pw, ": ", PRIMARY_T, " T cells vs ", PRIMARY_M,
                        "  (TDN = CAR+ product vs D0 myeloid)"),
         x = "Myeloid module score", y = "T-cell module score")
}

pdf(file.path(OUT$corr, "correlation_primary_CD8TEM_vs_myeloid.pdf"), width = 15, height = 8)
for (pw in c("IFN", "TNFa")) print(plot_prim(pw))
dev.off()

print(filter(cor_all, primary), n = Inf)
log_msg("Done -> ", OUT$corr)