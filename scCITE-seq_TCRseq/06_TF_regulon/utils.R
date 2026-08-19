### Adopted from Seurat Extend Package at https://github.com/huayc09/SeuratExtend
### All functions are modified or adopted from SeuratExtend github
### Citation 
### Hua, Y., Weng, L., Zhao, F., and Rambow, F. (2025). SeuratExtend: streamlining single-cell RNA-seq analysis through an integrated and intuitive framework. Gigascience 14, giaf076. https://doi.org/10.1093/gigascience/giaf076.

#Calculate Statistics 
WaterfallPlot_Calc <- function(
    matr,
    f,
    ident.1,
    ident.2,
    exp.transform,
    order,
    length,
    color,
    len.threshold,
    col.threshold,
    top.n,
    log.base = "e",
    pseudocount = NULL,
    max.cells.per.ident = Inf,
    random.seed = 1
) {
  library(dplyr)
  library(rlist)
  library(tidyr)
  
  # Validate log.base at the beginning of the function
  valid_bases <- c("e", "2", "10")
  use_natural_log <- FALSE
  
  # Convert log.base to numeric if it's a numeric string
  if (is.character(log.base) && !(log.base %in% valid_bases)) {
    # Try to convert to numeric
    numeric_base <- suppressWarnings(as.numeric(log.base))
    if (is.na(numeric_base)) {
      warning("Invalid log.base '", log.base, "' provided, using natural log (base e) instead")
      log.base <- "e"
      use_natural_log <- TRUE
    }
  }
  
  # Validate if it's a non-character, non-numeric value
  if (!is.character(log.base) && !is.numeric(log.base)) {
    warning("Invalid log.base type provided, using natural log (base e)")
    log.base <- "e"
    use_natural_log <- TRUE
  }
  
  f <- factor(f)
  ident.1 <- ident.1 %||% levels(f)[1]
  
  print(ident.1)
  
  # Downsample cells if max.cells.per.ident is set
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    
    # Get indices of cells for each identity
    idx.1 <- which(f == ident.1)
    idx.2 <- if(is.null(ident.2)) which(f != ident.1) else which(f == ident.2)
    
    # Determine ident.2 name for logging
    ident.2.name <- if(is.null(ident.2)) paste0("non-", ident.1) else ident.2
    
    # Track original sizes
    n1_original <- length(idx.1)
    n2_original <- length(idx.2)
    
    # Downsample each group if needed
    if (length(idx.1) > max.cells.per.ident) {
      idx.1 <- sample(x = idx.1, size = max.cells.per.ident)
      message("Downsampling ", ident.1, ": ", n1_original, " -> ", max.cells.per.ident, " cells")
    }
    if (length(idx.2) > max.cells.per.ident) {
      idx.2 <- sample(x = idx.2, size = max.cells.per.ident)
      message("Downsampling ", ident.2.name, ": ", n2_original, " -> ", max.cells.per.ident, " cells")
    }
    
    # Subset matrix and factor to selected cells only
    idx.use <- c(idx.1, idx.2)
    matr <- matr[, idx.use, drop = FALSE]
    f <- f[idx.use]
  }
  
  # Create logical vectors for statistical calculations
  # These are created AFTER potential downsampling, so they always match the current matr
  cell.1 <- (f == ident.1)
  cell.2 <- if(is.null(ident.2)) f != ident.1 else f == ident.2
  
  print(cell.1)
  print(cell.2)
  
  # Check if logFC calculation will be used
  will_use_logfc <- "logFC" %in% c(length, color)
  
  # Prepare matrix for analysis
  if(exp.transform) {
    original_matr <- matr
    matr <- expm1(matr)
  }
  
  # Automatically determine pseudocount if not provided and logFC will be used
  if (is.null(pseudocount) && will_use_logfc) {
    # Check the range on the transformed matrix if applicable
    check_matr <- if(exp.transform) matr else matr
    data_range <- range(check_matr, na.rm = TRUE)
    
    if (data_range[1] >= 0 && data_range[2] <= 1) {
      # 0-1 range data (like AUCell)
      pseudocount <- 0.01
      message("Data range detected as 0-1. Using pseudocount = 0.01 for logFC calculation.")
    } else {
      # Default for count data or other types
      pseudocount <- 1
      message("Using pseudocount = 1 for logFC calculation.")
    }
    
    # Check for negative values and warn if found
    if (data_range[1] < 0) {
      warning("Negative values detected in data. LogFC calculation may be affected.")
    }
  } else if (is.null(pseudocount)) {
    # Set default pseudocount silently if logFC is not used
    pseudocount <- 1
  }
  
  # Helper function to handle t.test with zero variance
  safe_ttest <- function(x, idx1, idx2) {
    # Handle all-zero or constant value cases
    if (all(x == 0) || length(unique(x)) == 1) {
      return(list(statistic = 0, p.value = 1))
    }
    # Perform t.test
    tryCatch({
      t.test(x[idx1], x[idx2])
    }, error = function(e) {
      list(statistic = 0, p.value = 1)
    })
  }
  
  scores <- list()
  if("tscore" %in% c(length, color)){
    scores[["tscore"]] <-
      apply(matr, 1, function(x) safe_ttest(x, cell.1, cell.2)[["statistic"]], simplify = TRUE)
  }
  if("p" %in% c(length, color)){
    print("calculating p")
    scores[["p"]] <-
      apply(matr, 1, function(x){
        p_value <- safe_ttest(x, cell.1, cell.2)[["p.value"]]
        log10p <- pmin(-log10(p_value), 325) * ifelse(mean(x[cell.1]) > mean(x[cell.2]), 1, -1)
        return(log10p)
      }, simplify = TRUE)
  }
  if("logFC" %in% c(length, color)){
    print('calculating logfc')
    # Calculate logFC with the specified base and pseudocount
    scores[["logFC"]] <- apply(matr, 1, function(x) {
      ratio <- mean(x[cell.1] + pseudocount) / mean(x[cell.2] + pseudocount)
      if (log.base == "e" || use_natural_log) {
        return(log(ratio))
      } else if (log.base == "2") {
        return(log2(ratio))
      } else if (log.base == "10") {
        return(log10(ratio))
      } else {
        # Must be a valid numeric at this point
        base <- as.numeric(log.base)
        return(log(ratio, base))
      }
    }, simplify = TRUE)
  }
  
  scores <- as.data.frame(list.cbind(scores))
  scores <- scores %>%
    mutate(rank = rownames(.),
           length = .[[length]],
           color = .[[color]])
  scores <- scores[
    abs(scores$length) > len.threshold &
      abs(scores$color) > col.threshold,
  ]
  
  if(order) {
    scores <- arrange(scores, desc(length))
  }
  if(!is.null(top.n)) {
    if(length(top.n) == 1) top.n <- c(top.n, top.n)
    scores2 <- split(scores, ifelse(scores$length > 0, "pos","neg"))
    scores <- rbind( head(scores2[["pos"]], top.n[1]), tail(scores2[["neg"]], top.n[2]) )
  }
  return(scores)
}

# plot waterfall plot 
WaterfallPlot_Plot_Custom <- function(
    scores,
    color = "p",
    color_theme,
    center_color= TRUE,
    flip = FALSE,
    y.label= NULL,
    angle= NULL,
    hjust= NULL,
    vjust= NULL,
    title= NULL,
    style= c("bar"),
    border = NA
) {
  library(ggplot2)
  library(scales)
  
  # Remove the manual reversal when flip=TRUE
  if(flip){
    scores <- scores[nrow(scores):1, ]
  }
  # Score preprocessing
  scores$rank <- factor(scores$rank, levels = unique(scores$rank))
  
  # Auto-determine angle, hjust, and vjust if not provided
  if (is.null(angle)) {
    max_label_length <- max(nchar(as.character(scores$rank)))
    if (flip) {
      angle <- 0  # Horizontal labels for flipped plot
    } else {
      angle <- if (max_label_length <= 2) 0 else -90  # Vertical for longer labels
    }
  }
  
  if (abs(angle) > 90) {
    warning("Angle should be between -90 and 90 degrees for optimal readability.")
  }
  
  if (is.null(hjust)) {
    if (angle > 0) {
      hjust <- 1  # Right align
    } else if (angle < 0) {
      hjust <- 0  # Left align
    } else {
      hjust <- 0.5  # Center align
    }
  }
  
  if (is.null(vjust)) {
    if (abs(angle) == 90) {
      vjust <- 0.5
    } else {
      vjust <- 1
    }
  }
  
  # Set fill label
  if(color == "p") lab_fill <- "-log10(p)" else lab_fill <- color
  
  # Base plot with common elements
  p <- ggplot(scores, aes(x = rank, y = length, colour = color)) +
    labs(fill = lab_fill, color = lab_fill, x = element_blank(), y = y.label, title = title)
  
  # Add style-specific elements
  if (style == "bar") {
    p <- p + geom_bar(stat = "identity", aes(fill = color), colour = border) +
      theme_classic()
  } else {
    p <- p +
      geom_segment(aes(xend = rank, y = 0, yend = length), size = 0.8) +
      geom_point(aes(y = length), size = 3) +
      theme_minimal()
  }
  p <- p +
    theme(axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  value_range <- range(scores$color)
  if (style == "bar") {
    p <- p + scale_fill_cont_auto(color_theme, center_color = center_color, value_range = value_range)
  } else {
    p <- p + scale_color_cont_auto(color_theme, center_color = center_color, value_range = value_range)
  }
  
  if(flip) p <- p + scale_x_discrete(position = "top") + coord_flip()
  return(p)
}

scale_fill_cont_auto <- function(color_scheme, center_color = NULL, value_range = NULL) {
  names(color_scheme) <- c("low", "mid", "high")
  
  return(ggplot2::scale_fill_gradient2(low = color_scheme["low"],
                              mid = color_scheme["mid"],
                              high = color_scheme["high"]))
         
  }

Seu2Matr <-
  function(
    seu,
    features,
    group.by = NULL,
    split.by = NULL,
    cells = NULL,
    slot = "data",
    assay = NULL,
    priority = c("expr","none"),
    verbose = TRUE
  ) {
    if(!require(SeuratObject)) library(Seurat)
    
    if(!is.null(assay)) DefaultAssay(seu) <- assay
    
    cells <- cells %||% colnames(seu)
    if(is.logical(cells)) {
      if(length(cells) != ncol(seu)) {
        stop("Logical value of 'cells' should be the same length as cells in ",
             "Seurat object")
      } else {
        cells.l <- cells
        cells <- colnames(seu)[cells]
      }
    } else {
      if(all(!cells %in% colnames(seu))) {
        stop("'cells' not found in Seurat object")
      }else if(any(!cells %in% colnames(seu))) {
        cells.out <- setdiff(cells, colnames(seu))
        stop(length(cells.out), " cell(s) not found in Seurat object: '",
             cells.out[1], "'...")
      }
      cells.l <- colnames(seu) %in% cells
    }
    
    if(priority[1] == "expr") {
      check.rep <- intersect(rownames(seu), colnames(seu@meta.data))
      if(any(features %in% check.rep)) {
        feature.rep <- intersect(features, check.rep)
        if(verbose) {
          warning("'Features' found in both expression matrix and 'meta.data' (",
                  feature.rep[1],
                  "...). Ignore the 'meta.data'. If you need to fetch data from ",
                  "'meta.data' please set 'priority' to 'none'")
        }
        seu@meta.data[,feature.rep] <- NULL
      }
    }
    
    if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      matr <- FetchData(object = seu, vars = features, cells = cells, layer = slot, clean = "none")
    } else {
      matr <- FetchData(object = seu, vars = features, cells = cells, slot = slot)
    }
    
    if(is.null(group.by)) {
      f <- factor(Idents(seu)[cells])
    }else if(length(group.by) == 1) {
      if(group.by %in% colnames(seu@meta.data)){
        f <- factor(seu[[group.by]][cells,])
      }else{
        stop("Cannot find '", group.by, "' in meta.data")
      }
    }else if(length(group.by) == ncol(seu)) {
      f <- factor(group.by[cells.l])
    }else if(length(group.by) == length(cells)) {
      f <- factor(group.by)
    }else{
      stop("'group.by' should be variable name in 'meta.data' or ",
           "string with the same length of cells")
    }
    names(f) <- cells
    
    if(is.null(split.by)) {
      f2 <- NULL
    }else if(length(split.by) == 1) {
      if(split.by %in% colnames(seu@meta.data)){
        f2 <- factor(seu[[split.by]][cells,])
        names(f2) <- cells
      }else{
        stop("Cannot find '", split.by, "' in meta.data")
      }
    }else if(length(split.by) == ncol(seu)) {
      f2 <- factor(split.by[cells.l])
      names(f2) <- cells
    }else if(length(split.by) == length(cells)) {
      f2 <- factor(split.by)
      names(f2) <- cells
    }else{
      stop("'split.by' should be variable name in 'meta.data' or ",
           "string with the same length of cells")
    }
    return(list(matr = matr, f = f, f2 = f2))
  }

# Check if ident.1 and ident.2 are in the factor levels and return
# logical vectors
CheckIdent <-
  function(
    f,
    ident.1 = NULL,
    ident.2 = NULL
  ) {
    lv <- levels(f)
    nlv <- 1:nlevels(f)
    
    ident.1 <- ident.1 %||% lv[1]
    if(ident.1 %in% lv) {
      cell.1 <- (f == ident.1)
      name.1 <- ident.1
    } else if(ident.1 %in% nlv) {
      name.1 <- lv[ident.1]
      cell.1 <- (f == name.1)
    } else stop("'ident.1' not found. Possible values: \n",
                paste(lv, collapse = ", "))
    
    if(is.null(ident.2)) {
      name.2 <- paste0("non-", name.1)
      cell.2 <- (f != name.1)
    } else if(ident.2 %in% lv) {
      cell.2 <- (f == ident.2)
      name.2 <- ident.2
    } else if(ident.2 %in% nlv) {
      name.2 <- lv[ident.2]
      cell.2 <- (f == name.2)
    } else stop("'ident.2' not found. Possible values: \n",
                paste(lv, collapse = ", "))
    
    if(identical(name.1, name.2))
      warning("'ident.1' and 'ident.2' are the same??")
    
    return(list(
      cell.1 = cell.1,
      name.1 = name.1,
      cell.2 = cell.2,
      name.2 = name.2
    ))
  }


# imported from rlang

`%||%` <- function (x, y)
{
  if (is.null(x))
    y
  else x
}

is_empty <- function (x) length(x) == 0

# Check species

check_spe <- function(spe) {
  if (is.null(spe)) {
    stop("Species not defined. Please set with options(spe = \"your_species\") or provide directly to the function. Default supported species are 'human' and 'mouse'. For custom species, please refer to the documentation on using custom data.")
  }
  return(TRUE)
}

# Check package version
check_pkg_version <- function(package, version.min) {
  version.min <- package_version(version.min)
  version.pkg <- packageVersion(package)
  if(version.pkg < version.min) {
    stop("Current version of package '", package,"' is ", version.pkg,
         ", but version ", version.min, " is required")
  }
}

# message only once
message_only_once <- function(title, message) {
  option.title <- paste0("seuratextend_warning_",title)
  if(getOption(option.title, TRUE)) {
    message(message, "\nThis message is shown once per session")
    options(setNames(list(FALSE), option.title))
  }
}

# Helper function to handle GetAssayData compatibility between Seurat v4 and v5
.get_assay_data_compat <- function(object, slot = "data", assay = NULL) {
  # Check if we're using SeuratObject v5 or later
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    return(GetAssayData(object, layer = slot, assay = assay))
  } else {
    return(GetAssayData(object, slot = slot, assay = assay))
  }
}

# Internal function to limit values to a specified range
# Similar to Seurat's MinMax function
# Used by DotPlot2, Heatmap, and potentially other visualization functions
# @param data Numeric vector or matrix to be clipped
# @param min Minimum cutoff value. Values below this will be set to min. Default: NULL (no lower limit)
# @param max Maximum cutoff value. Values above this will be set to max. Default: NULL (no upper limit)
# @return Data with values clipped to the specified range
MinMax_internal <- function(data, min = NULL, max = NULL) {
  if (!is.null(min)) {
    data[data < min] <- min
  }
  if (!is.null(max)) {
    data[data > max] <- max
  }
  return(data)
}
