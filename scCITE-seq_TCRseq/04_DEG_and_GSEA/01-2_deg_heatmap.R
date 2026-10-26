## =============================================================================
## CR vs PD DEG overlaps: CD8 TEM and CD4 TEM-like T cells (CAR+ / CAR-) and myeloid cells
## Input: Seurat FindMarkers (wilcox) output, combined into long CSVs
## Convention: avg_log2FC > 0 = up in CR, avg_log2FC < 0 = up in PD
## Output: gene list / call / category CSVs and the top-N "simple" heatmaps
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(tibble)
  library(ComplexHeatmap); library(circlize)
})
# BiocManager::install("ComplexHeatmap")

## ---- 1. SETTINGS (edit these) -------------------------------------------------
deg_input_dir <- file.path("images", "fig1_supplemental_draft", "files", "new_degs", "summary")

gene_mode <- "all_sig"     # "top_n" = top N per condition | "all_sig" = all significant genes
top_n     <- 50            # used only when gene_mode == "top_n"

run_tag <- if (gene_mode == "top_n") paste0("top", top_n) else "allSig"
out_dir <- file.path("images", "fig1_supplemental_draft", paste0("CR_vs_PD_overlaps_", run_tag, "_01log2FC"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


car_file     <- "CAR_annos_level2_New_unfiltered_fullDEGs.csv"
noncar_file  <- "NonCAR_Annos_Level2_New_unfiltered_fullDEGs.csv"
myeloid_file <- "NonCAR_Annos_Level1_unfiltered_fullDEGs.csv"          # myeloid cells live in the non-CAR file

# Column names in the CSVs
col_gene     <- "geneID"
col_lfc      <- "avg_log2FC"
col_padj     <- "p_val_adj"
col_celltype <- "celltype"
col_time     <- "timepoint"

# Labels exactly as they appear in the celltype column
# T cell types, each run separately; names = tag used in output file names
tcell_types   <- c(CD8TEM = "CD8 TEM", CD4TEMlike = "CD4 TEM-like")
myeloid_types <- c("CD14 Mono", "CD16 Mono", "Dendritic Cell")

# Timepoint order for plots (NULL = order found in file)
time_levels <- c("APH", "D0", "TDN", "Peak", "4W")

padj_cut    <- 0.01
lfc_cut     <- 0.1
pct_cut     <- 0.05
drop_mt     <- TRUE
drop_lncrna <- TRUE
heat_top_n  <- 50      # top N per condition used to build the heatmap gene pools (section 3)
aph_times     <- c("APH")
d0_times      <- c("D0", "Day0")
product_times <- c("TDN", "Product")          # infusion product
post_times    <- c("Peak", "Week4", "4W")
pre_times     <- c(aph_times, d0_times)

# Y-chromosome genes excluded from the heatmap top-N selection
y_genes <- c("RPS4Y1", "RPS4Y2", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "ZFY", "USP9Y",
             "NLGN4Y", "TXLNGY", "PRKY", "TMSB4Y", "SRY", "TBL1Y", "PCDH11Y", "AMELY",
             "TSPY1", "RBMY1A1", "DAZ1", "BCORP1", "MTRNR2L8", "MTRNR2L12", "RPS4X")
y_pattern <- "^TTTY\\d"     # Y-linked non-coding transcripts
heat_exclude_pattern <- "^TRAV|^TRBV"     # TCR variable genes (clonotype-driven)


## ---- 2. Read data --------------------------------------------------------------
# remove lncRNA / non-coding naming patterns (HGNC + GENCODE clone names)
lnc_pattern <- paste(c(
  "^LINC\\d",                 # LINC00861
  "-AS\\d*$",                 # antisense: XXX-AS1
  "-DT$",                     # divergent transcripts: XXX-DT
  "-IT\\d*$", "-OT\\d*$",     # intronic / overlapping transcripts
  "^SNHG\\d",                 # snoRNA host genes
  "^MIR\\d+.*HG",             # miRNA host genes: MIR155HG
  "^[A-Z]{1,2}\\d{5,6}\\.\\d+$"  # clone-based IDs: AC012345.1, AL353135.1, Z93241.1
), collapse = "|")

keep_gene <- function(g) {
  ok <- rep(TRUE, length(g))
  if (drop_mt)     ok <- ok & !grepl("^MT-", g)
  if (drop_lncrna) ok <- ok & !grepl(lnc_pattern, g)
  ok
}

read_deg <- function(f, car_label) {
  z <- read_csv(file.path(deg_input_dir, f), show_col_types = FALSE)
  stopifnot(all(z[[col_lfc]][z$cluster == "CR"] > 0),
            all(z[[col_lfc]][z$cluster == "PD"] < 0))   # sign check
  z %>% transmute(gene     = .data[[col_gene]],
                  lfc      = .data[[col_lfc]],
                  padj     = .data[[col_padj]],
                  pct1     = pct.1,                 # CR
                  pct2     = pct.2,                 # PD
                  celltype = .data[[col_celltype]],
                  time     = as.character(.data[[col_time]]),
                  car      = car_label) %>%
    # CR genes need pct.1 > cutoff, PD genes need pct.2 > cutoff
    mutate(pct_ok = if_else(lfc > 0, pct1 > pct_cut, pct2 > pct_cut)) %>%
    filter(keep_gene(gene)) %>%                       # drop MT- and lncRNA genes
    distinct(gene, celltype, time, car, .keep_all = TRUE)
  
}

tcell_all <- bind_rows(read_deg(car_file, "CAR+"), read_deg(noncar_file, "CAR-")) %>%
  filter(celltype %in% tcell_types)
myel  <- read_deg(myeloid_file, "none") %>% filter(celltype %in% myeloid_types)
stopifnot(all(tcell_types %in% tcell_all$celltype), nrow(myel) > 0)   # fails -> fix labels in SETTINGS

tl <- if (is.null(time_levels)) unique(c(tcell_all$time, myel$time)) else time_levels
# one table per T cell type
tcell_l <- lapply(tcell_types, function(ct)
  tcell_all %>% filter(celltype == ct) %>%
    mutate(time = factor(time, levels = tl)) %>% arrange(car, time))
myel  <- myel %>% mutate(time = factor(time, levels = tl),
                         celltype = factor(celltype, levels = myeloid_types)) %>%
  arrange(celltype, time)

## ---- 3. STEP 1: top-N per condition -> combined pool ----------------------------
top_genes <- function(d, dir) {
  g <- d %>% filter(padj < padj_cut, abs(lfc) > lfc_cut, pct_ok,
                    if (dir == "CR") lfc > 0 else lfc < 0) %>%
    arrange(desc(abs(lfc)), padj)
  if (gene_mode == "top_n") g <- slice_head(g, n = top_n)
  pull(g, gene)
}
top_table <- function(d, keys) {
  d %>% nest(data = -all_of(keys)) %>% expand_grid(dir = c("CR", "PD")) %>%
    mutate(genes = map2(data, dir, top_genes)) %>% select(-data)
}

t_top_l  <- lapply(tcell_l, top_table, keys = c("car", "time"))
m_top    <- top_table(myel,  c("celltype", "time"))
t_pool_l <- lapply(t_top_l, function(x) unique(unlist(x$genes)))
m_pool   <- unique(unlist(m_top$genes))
for (tag in names(tcell_types)) message(tcell_types[[tag]], " pool: ", length(t_pool_l[[tag]]), " genes")
message("Myeloid pool: ", length(m_pool), " genes")

for (tag in names(tcell_types))
  t_top_l[[tag]] %>% unnest(genes) %>%
  write_csv(file.path(out_dir, paste0("Tcell_", tag, "_", run_tag, "_lists.csv")))
m_top %>% unnest(genes) %>% write_csv(file.path(out_dir, paste0("Myeloid_", run_tag, "_lists.csv")))

heat_pool_for <- function(d, keys) {
  d %>% filter(padj < padj_cut, abs(lfc) > lfc_cut, pct_ok,
               !gene %in% y_genes, !grepl(y_pattern, gene),
               !grepl(heat_exclude_pattern, gene)) %>%          # drop TRAV / TRBV
    mutate(dir = if_else(lfc > 0, "CR", "PD")) %>%
    group_by(across(all_of(keys)), dir) %>%
    arrange(desc(abs(lfc)), padj, .by_group = TRUE) %>%
    slice_head(n = heat_top_n) %>% ungroup() %>%
    pull(gene) %>% unique()
}
t_heat_pool_l <- lapply(tcell_l, heat_pool_for, keys = c("car", "time"))
m_heat_pool <- heat_pool_for(myel,  c("celltype", "time"))


## ---- 4. STEP 2: call every pool gene in every condition ------------------------
call_table <- function(d, pool, keys) {
  d %>% filter(gene %in% pool) %>%
    mutate(call = case_when(padj < padj_cut & lfc >  lfc_cut & pct_ok ~ "CR",
                            padj < padj_cut & lfc < -lfc_cut & pct_ok ~ "PD",
                            TRUE ~ "ns")) %>%
    select(gene, all_of(keys), lfc, padj, pct1, pct2, call) %>%
    complete(gene, nesting(!!!rlang::syms(keys)), fill = list(call = "ns"))
}
t_calls_l <- lapply(names(tcell_types), function(tag)
  call_table(tcell_l[[tag]], t_pool_l[[tag]], c("car", "time")))
names(t_calls_l) <- names(tcell_types)
m_calls <- call_table(myel,  m_pool, c("celltype", "time"))
for (tag in names(tcell_types))
  write_csv(t_calls_l[[tag]], file.path(out_dir, paste0("Tcell_", tag, "_pool_calls.csv")))
write_csv(m_calls, file.path(out_dir, "Myeloid_pool_calls.csv"))

## ---- 5. Comparison helpers (all operate on the calls) ---------------------------
all_dir <- function(cl, dir) cl %>% group_by(gene) %>%
  filter(all(call == dir)) %>% pull(gene) %>% unique()

unique_time <- function(cl, dir, tp) cl %>% group_by(gene) %>%
  filter(any(call == dir & time == tp), !any(call == dir & time != tp)) %>%
  pull(gene) %>% unique()

# Typed CR<->PD switches; first match wins (1 -> 2 -> 3 -> 4)
flip_typed <- function(cl) {
  cl %>% mutate(t = as.character(time)) %>%
    group_by(gene) %>%
    summarise(
      cr_pre  = any(call == "CR" & t %in% pre_times),     pd_pre  = any(call == "PD" & t %in% pre_times),
      cr_prod = any(call == "CR" & t %in% product_times), pd_prod = any(call == "PD" & t %in% product_times),
      cr_post = any(call == "CR" & t %in% post_times),    pd_post = any(call == "PD" & t %in% post_times),
      cr_aph  = any(call == "CR" & t %in% aph_times),     pd_aph  = any(call == "PD" & t %in% aph_times),
      cr_d0   = any(call == "CR" & t %in% d0_times),      pd_d0   = any(call == "PD" & t %in% d0_times),
      any_cr  = any(call == "CR"), any_pd = any(call == "PD"),
      up_CR_at = paste(unique(lab[call == "CR"]), collapse = ";"),
      up_PD_at = paste(unique(lab[call == "PD"]), collapse = ";"),
      .groups = "drop") %>%
    mutate(
      aph_x = (cr_aph & (pd_prod | pd_post)) | (pd_aph & (cr_prod | cr_post)),  # APH vs post-infusion
      d0_x  = (cr_d0  & (pd_prod | pd_post)) | (pd_d0  & (cr_prod | cr_post)),  # D0  vs post-infusion
      flip_type = case_when(
        aph_x & !d0_x                             ~ "Flip 1a: APH vs post-infusion",
        d0_x  & !aph_x                            ~ "Flip 1b: D0 vs post-infusion",
        aph_x & d0_x                              ~ "Flip 1c: APH & D0 vs post-infusion",
        (cr_prod & pd_post) | (pd_prod & cr_post) ~ "Flip 2: product vs Peak/Week4",
        (cr_aph & pd_d0) | (pd_aph & cr_d0)       ~ "Flip 3: APH vs D0",
        any_cr & any_pd                           ~ "Flip 4: other")) %>%
    filter(!is.na(flip_type)) %>%
    dplyr::select(gene, flip_type, up_CR_at, up_PD_at)
}
# Exact-combination shared categories (gene's full set of members where call == dir, size >= 2)
add_combos <- function(group, cl, dir, set_col, set_levels, suffix) {
  avail <- intersect(set_levels, unique(as.character(cl[[set_col]])))
  cc <- cl %>% filter(call == dir) %>%
    group_by(gene) %>%
    summarise(members = list(intersect(avail, unique(as.character(.data[[set_col]])))),
              .groups = "drop") %>%
    mutate(n = lengths(members), combo = map_chr(members, paste, collapse = "+")) %>%
    filter(n >= 2)
  if (nrow(cc) == 0) return(invisible(NULL))
  # biggest combos first, then in timepoint / celltype order
  cc <- cc %>%
    mutate(key = map_dbl(members, ~ sum(2^(length(avail) - match(.x, avail))))) %>%
    arrange(desc(n), desc(key))
  for (cb in unique(cc$combo)) {
    n   <- cc$n[cc$combo == cb][1]
    lab <- if (n == length(avail)) paste0("all (", cb, ")") else cb
    add(group, paste0(dir, " | shared ", n, ": ", lab, " | ", suffix), cc$gene[cc$combo == cb])
  }
}

cats <- tibble(group = character(), category = character(), gene = character())
add <- function(group, category, genes)
  if (length(genes)) cats <<- bind_rows(cats, tibble(group, category, gene = genes))

times_in <- function(cl) intersect(tl, unique(as.character(cl$time)))

## ---- 6. T cell comparisons (each T cell type separately) -----------------------
# category group = "Tcell <cell type>", e.g. "Tcell CD8 TEM"
t_flip_l <- list()
for (tag in names(tcell_types)) {
  ct      <- tcell_types[[tag]]
  grp     <- paste("Tcell", ct)
  t_calls <- t_calls_l[[tag]]
  for (d in c("CR", "PD")) {
    add(grp, paste0(d, " | all time x CAR"), all_dir(t_calls, d))
    for (cs in c("CAR+", "CAR-")) {
      cl <- filter(t_calls, car == cs)
      add_combos(grp, cl, d, "time", tl, cs)                     # every timepoint combination
      for (tp in times_in(cl)) add(grp, paste0(d, " | unique ", tp, " | ", cs), unique_time(cl, d, tp))
    }
  }
  # pooled across CAR- and CAR+ so pre-infusion (CAR-) can be compared with CAR+ timepoints
  t_flip <- t_calls %>% mutate(lab = paste(car, time)) %>% flip_typed()
  for (ft in sort(unique(t_flip$flip_type)))
    add(grp, paste0(ft, " | ", ct), t_flip$gene[t_flip$flip_type == ft])
  t_flip_l[[tag]] <- t_flip
}

## ---- 7. Myeloid comparisons ----------------------------------------------------
for (d in c("CR", "PD")) {
  add("Myeloid", paste0(d, " | all time x celltype"), all_dir(m_calls, d))
  for (ct in myeloid_types) {
    cl <- filter(m_calls, celltype == ct)
    add_combos("Myeloid", cl, d, "time", tl, ct)                   # timepoint combos per cell type
    for (tp in times_in(cl)) add("Myeloid", paste0(d, " | unique ", tp, " | ", ct), unique_time(cl, d, tp))
  }
  for (tp in times_in(m_calls))
    add_combos("Myeloid", filter(m_calls, time == tp), d, "celltype", myeloid_types, tp)  # celltype combos per timepoint
}
m_flip <- m_calls %>% mutate(lab = as.character(time)) %>%
  group_by(celltype) %>% group_modify(~ flip_typed(.x)) %>% ungroup()
for (ft in sort(unique(m_flip$flip_type)))
  for (ct in myeloid_types) {
    g <- m_flip$gene[m_flip$flip_type == ft & m_flip$celltype == ct]
    add("Myeloid", paste0(ft, " | ", ct), g)
  }

## ---- 8. Save category tables ---------------------------------------------------
write_csv(cats, file.path(out_dir, "gene_categories_long.csv"))
cats %>% count(group, category, name = "n_genes") %>%
  write_csv(file.path(out_dir, "gene_categories_counts.csv"))
for (tag in names(tcell_types))
  write_csv(t_flip_l[[tag]], file.path(out_dir, paste0("Tcell_", tag, "_flip_genes.csv")))
write_csv(m_flip, file.path(out_dir, "Myeloid_flip_genes.csv"))
print(count(cats, group, category, name = "n_genes"), n = Inf)

## ---- 9. Heatmaps --------------------------------------------------------------
# heatmap gene pools were built in section 3 (top heat_top_n per condition)
fold_unique_post <- FALSE

# Myeloid heatmap: collapse timepoint detail, keep cell type
collapse_cats <- function(category) {
  dir   <- sub(" .*", "", category)                           # CR / PD / Flip
  last  <- sub(".*\\| ", "", category)                        # trailing field: celltype or timepoint
  combo <- sub(".*\\| shared \\d+: (.*) \\|.*", "\\1", category)
  combo <- sub("^all \\((.*)\\)$", "\\1", combo)              # strip "all ( )"
  case_when(
    startsWith(category, "Flip")                                  ~ paste("Flip | CR<->PD switch |", last),
    grepl("all time x", category)                                 ~ paste(dir, "| shared across all"),
    grepl("\\| shared \\d+:", category) & last %in% myeloid_types ~ paste(dir, "| shared timepoints |", last),
    grepl("\\| shared \\d+:", category)                           ~ paste0(dir, " | shared across cell types: ", combo),
    grepl("unique", category)                                     ~ paste(dir, "| unique to one timepoint |", last),
    TRUE                                                          ~ category)
}
# T cell heatmap: collapse timepoint combinations into pre / post-infusion blocks
collapse_tcell <- function(category) {
  all_post <- c(product_times, post_times)
  dir   <- sub(" .*", "", category)
  last  <- sub(".*\\| ", "", category)                               # CAR status
  combo <- sub(".*\\| shared \\d+: (.*) \\|.*", "\\1", category)
  combo <- sub("^all \\((.*)\\)$", "\\1", combo)
  mem   <- strsplit(combo, "+", fixed = TRUE)
  has_aph  <- map_lgl(mem, ~ any(.x %in% aph_times))
  has_d0   <- map_lgl(mem, ~ any(.x %in% d0_times))
  has_post <- map_lgl(mem, ~ any(.x %in% all_post))
  combo_lab <- case_when(
    has_aph & has_d0 & has_post ~ "APH + D0 + post-infusion",
    has_aph & has_post          ~ "APH + post-infusion",
    has_d0  & has_post          ~ "D0 + post-infusion",
    has_aph & has_d0            ~ "APH + D0 (pre-infusion)",
    has_post                    ~ "post-infusion only",
    TRUE                        ~ combo)
  uniq_tp <- sub(".*\\| unique (.*) \\|.*", "\\1", category)
  case_when(
    startsWith(category, "Flip")                 ~ category,
    grepl("all time x", category)                ~ paste(dir, "| shared across all"),
    grepl("\\| shared \\d+:", category)          ~ paste0(dir, " | ", combo_lab, " | ", last),
    grepl("unique", category) & fold_unique_post &
      uniq_tp %in% all_post                      ~ paste0(dir, " | post-infusion only | ", last),
    TRUE                                         ~ category)
}

# Myeloid "simple" categories, priority 1-9
myeloid_simple <- function(cl) {
  cl  <- cl %>% mutate(t = as.character(time), ct = as.character(celltype))
  cts <- intersect(myeloid_types, unique(cl$ct))
  tps <- times_in(cl)
  out <- list()
  
  # 1-2: pre (APH/D0) <-> post (Peak/Week4) switch, per cell type
  out$flip <- cl %>% group_by(gene, ct) %>%
    summarise(flip = (any(call == "CR" & t %in% pre_times) & any(call == "PD" & t %in% post_times)) |
                (any(call == "PD" & t %in% pre_times) & any(call == "CR" & t %in% post_times)),
              .groups = "drop") %>%
    filter(flip) %>% group_by(gene) %>%
    summarise(cts = paste(intersect(cts, ct), collapse = "+"), n = n(), .groups = "drop") %>%
    transmute(gene,
              category = if_else(n >= 2,
                                 paste("Flip | pre vs post-infusion, multiple cell types |", cts),
                                 paste("Flip | pre vs post-infusion |", cts)),
              pri = if_else(n >= 2, 1, 2))
  
  # 3-9: same direction patterns
  for (d in c("CR", "PD")) {
    out[[d]] <- cl %>% group_by(gene) %>%
      summarise(n_hit  = sum(call == d),
                all_ct = list(cts[map_lgl(cts, function(x) all(call[ct == x] == d))]),
                all_tp = list(tps[map_lgl(tps, function(x) all(call[t  == x] == d))]),
                # timepoints where 2+ (but not all) cell types share direction d
                multi  = list(unlist(lapply(tps, function(x) {
                  h <- intersect(cts, ct[call == d & t == x])
                  if (length(h) >= 2 && length(h) < length(cts))
                    paste0(paste(h, collapse = "+"), "|", x)
                }))),
                hit1   = paste0(ct[call == d], " @ ", t[call == d])[1],
                .groups = "drop") %>%
      mutate(n_ct = lengths(all_ct), n_tp = lengths(all_tp), n_multi = lengths(multi),
             ct_lab    = map_chr(all_ct, paste, collapse = "+"),
             tp_lab    = map_chr(all_tp, phase_lab),
             multi_lab = map_chr(multi, function(v) {
               if (!length(v)) return("")
               pr <- sub("\\|.*", "", v); tp <- sub(".*\\|", "", v)
               paste(map_chr(unique(pr), ~ paste0(.x, " @ ", phase_lab(tp[pr == .x]))),
                     collapse = "; ")
             }),
             pri = case_when(n_ct == length(cts) ~ 3,
                             n_ct == 2           ~ 4,
                             n_tp >= 1           ~ 5,
                             n_multi >= 1        ~ 6,
                             n_ct == 1           ~ 7,
                             n_hit == 1          ~ 8,
                             n_hit >= 2          ~ 9),
             category = case_when(
               pri == 3 ~ paste(d, "| all timepoints, all cell types"),
               pri == 4 ~ paste(d, "| all timepoints |", ct_lab),
               pri == 5 ~ paste(d, "| all cell types |", tp_lab),
               pri == 6 ~ paste(d, "| multiple cell types |", multi_lab),
               pri == 7 ~ paste(d, "| all timepoints |", ct_lab),
               pri == 8 ~ paste(d, "| unique |", hit1),
               pri == 9 ~ paste(d, "| other partial pattern"))) %>%
      filter(!is.na(pri)) %>% dplyr::select(gene, category, pri)
  }
  bind_rows(out) %>% arrange(pri, category)
}
cat_priority <- function(category) {
  n_shared <- suppressWarnings(as.numeric(sub(".*\\| shared (\\d+):.*", "\\1", category)))
  case_when(
    startsWith(category, "Flip 1a")    ~ 1.1,
    startsWith(category, "Flip 1b")    ~ 1.2,
    startsWith(category, "Flip 1c")    ~ 1.3,
    startsWith(category, "Flip 2")     ~ 2,
    startsWith(category, "Flip 3")     ~ 3,
    startsWith(category, "Flip 4")     ~ 4,
    grepl("all time x", category)      ~ 5,               # shared across everything
    !is.na(n_shared)                   ~ 6 - n_shared/100, # exact combos, larger first
    grepl("unique", category)          ~ 7,
    startsWith(category, "Flip |")              ~ 1,
    grepl("shared across all", category)        ~ 5,
    grepl("shared across cell types", category) ~ 5.5,
    grepl("shared timepoints", category)        ~ 6,
    grepl("\\| APH \\+ D0 \\+ post-infusion", category) ~ 6.1,
    grepl("\\| APH \\+ post-infusion", category)        ~ 6.2,
    grepl("\\| D0 \\+ post-infusion", category)         ~ 6.3,
    grepl("\\| APH \\+ D0 \\(pre", category)            ~ 6.4,
    grepl("\\| post-infusion only", category)           ~ 6.5,
    grepl("\\| only in ", category)          ~ 7.5,
    TRUE                               ~ 8)
}
phase_lab <- function(tps_hit) {
  ph <- unique(case_when(tps_hit %in% pre_times                      ~ "pre",
                         tps_hit %in% c(product_times, post_times)   ~ "post",
                         TRUE                                        ~ tps_hit))
  if (all(ph == "pre"))  return("pre-infusion (APH/D0)")
  if (all(ph == "post")) return("post-infusion (Peak/Week4)")
  if (setequal(ph, c("pre", "post"))) return("pre + post-infusion")
  paste(ph, collapse = "+")                                         # unmapped timepoint label
}
# names are "<group> <timepoint>" -- must match the data exactly
tcell_col_order <- c("CAR- APH", "CAR- D0",
                     "CAR- TDN", "CAR+ TDN",
                     "CAR- Peak", "CAR+ Peak",
                     "CAR- 4W", "CAR+ 4W")

myel_col_order <- c(paste("Dendritic Cell", c("APH", "D0", "Peak", "4W")),
                    paste("CD14 Mono",      c("APH", "D0", "Peak", "4W")),
                    paste("CD16 Mono",      c("APH", "D0", "Peak", "4W")))
# Rows = genes split by category (gene shown in the first category it falls in),
# clustered within category. Values = log2FC from the full table (0 if not tested).
plot_heat <- function(d, cat_df, col_keys, file, title,
                      col_order = NULL, split_by = NULL,
                      heat_padj = padj_cut, heat_lfc = 0.1, pri_first = FALSE) {
  if (nrow(cat_df) == 0) return(invisible(NULL))
  
  pri_vec <- if ("pri" %in% names(cat_df)) cat_df$pri else cat_priority(cat_df$category)
  cat1 <- cat_df %>%
    mutate(.ord = row_number(), .pri = pri_vec) %>%
    arrange(.pri, .ord) %>% distinct(gene, .keep_all = TRUE) %>%
    mutate(.dir = factor(sub(" .*", "", category), levels = c("CR", "PD", "Flip")))
  cat1 <- if (pri_first) arrange(cat1, .pri, .dir, .ord) else arrange(cat1, .dir, .pri, .ord)
  
  d <- d %>% tidyr::unite("col", all_of(col_keys), sep = " ", remove = FALSE)
  m <- d %>% filter(gene %in% cat1$gene) %>%
    mutate(col = factor(col, levels = unique(col))) %>%
    dplyr::select(gene, col, lfc) %>%
    tidyr::pivot_wider(names_from = col, values_from = lfc, values_fill = 0, names_sort = TRUE) %>%
    tibble::column_to_rownames("gene") %>% as.matrix()
  
  # custom column order
  if (!is.null(col_order)) {
    miss <- setdiff(col_order, colnames(m))
    if (length(miss)) warning("Columns not found (check spelling): ", paste(miss, collapse = ", "),
                              "\nAvailable: ", paste(colnames(m), collapse = ", "))
    m <- m[, intersect(col_order, colnames(m)), drop = FALSE]
  }
  
  # gap at every change in split_by ("time" or "celltype")
  col_split <- NULL
  if (!is.null(split_by)) {
    col_grp   <- as.character(d[[split_by]][match(colnames(m), d$col)])
    col_split <- factor(col_grp, levels = unique(col_grp))
  }
  
  cat1 <- cat_df %>%
    mutate(.ord = row_number(), .pri = cat_priority(category)) %>%
    arrange(.pri, .ord) %>% distinct(gene, .keep_all = TRUE) %>%            # assign by priority
    mutate(.dir = factor(sub(" .*", "", category), levels = c("CR", "PD", "Flip"))) %>%
    arrange(.dir, .pri, .ord)      
  m <- m[cat1$gene, , drop = FALSE]
  
  # direction bar: CR / PD / Flip
  dir_lab <- sub(" .*", "", cat1$category)               # "CR", "PD", "Flip"
  row_ha  <- rowAnnotation(Direction = dir_lab,
                           col = list(Direction = c(CR = "#4FB7C5", PD = "#E7A75E", Flip = "grey50")),
                           show_annotation_name = FALSE)
  
  ht <- Heatmap(m, name = "log2FC\n(+CR / -PD)",
                col = colorRamp2(c(-2, 0, 2), c("#E7A75E", "white", "#4FB7C5")),
                left_annotation = row_ha,
                row_split = factor(cat1$category, levels = unique(cat1$category)),
                column_split = col_split, column_gap = unit(2, "mm"), row_gap = unit(1, "mm"),
                cluster_row_slices = FALSE, cluster_rows = TRUE,
                cluster_columns = FALSE, cluster_column_slices = FALSE,
                show_row_dend = FALSE,
                show_row_names = TRUE, row_names_gp = gpar(fontsize = 6),
                row_title_rot = 0, row_title_gp = gpar(fontsize = 7),
                column_names_rot = 45, column_names_gp = gpar(fontsize = 8),
                column_title = title, border = TRUE)
  
  pdf(file, width = 7 + ncol(m) * 0.4,
      height = min(200, max(6, nrow(m) * 0.11 + 3)))      # ~0.11 in per labelled row
  draw(ht); dev.off()
}


for (tag in names(tcell_types)) {
  ct <- tcell_types[[tag]]
  plot_heat(tcell_l[[tag]],
            cats %>% filter(group == paste("Tcell", ct), gene %in% t_heat_pool_l[[tag]]) %>%
              mutate(category = collapse_tcell(category)),
            c("car", "time"),
            file.path(out_dir, paste0("Heatmap_", tag, "_top50_simple.pdf")),
            paste0(ct, " - CR vs PD, simple (top ", heat_top_n, " per condition)"),
            tcell_col_order, "time")
}

plot_heat(myel,
          cats %>% filter(group == "Myeloid", gene %in% m_heat_pool) %>%
            mutate(category = collapse_cats(category)),
          c("celltype", "time"),
          file.path(out_dir, "Heatmap_Myeloid_all_categories_top50_simple.pdf"),
          paste0("Myeloid - CR vs PD (top ", heat_top_n, " per condition)"), myel_col_order, "celltype")

m_simple <- myeloid_simple(filter(m_calls, gene %in% m_heat_pool))
write_csv(m_simple, file.path(out_dir, "Myeloid_simple_categories.csv"))

