library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(sf)
library(jsonlite)
library(sfheaders)
library(ggrastr)

# ---- Paths (edit, or set CODEX_DIR in the environment) ----------------------
codex_root <- Sys.getenv("CODEX_DIR", unset = "CODEX")
setwd(codex_root)

codex_mapping <- list(
  "PT4"  = "31361C",
  "PT5"  = "SP2223235",
  "PT6"  = "38979",
  "PT7"  = "SP206535",
  "PT8"  = "184419A",
  "PT9"  = "34774",
  "PT10" = "25913A1",
  "PT11" = "CS207490",
  "PT12" = "402641A",
  "PT14" = "10858"
)

# accession -> patient label, used to name every output
accession_to_patient <- setNames(names(codex_mapping), unlist(codex_mapping))

# set data paths #

codex_path <- file.path(codex_root, "CODEX_SeuratAnalysis")
data_path  <- file.path(codex_root, "counts_raw")

# Per-sample inputs. `accession` is the on-disk sample folder / file name; the

samples <- data.frame(
  accession   = c("25913A1", "CS207490", "31361C", "392361A", "SP2223235",
                  "184419A", "38979", "402641A", "SP206535", "10858", "34774"),
  count_dir   = c("62124_slide1_25913A1_CS207490",
                  "62124_slide1_25913A1_CS207490",
                  "62124_slide2_31361C",
                  "71624_slide1_SP2223235_392361A",
                  "71624_slide1_SP2223235_392361A",
                  "71624_slide2_184419A",
                  "71824_slide1_38979",
                  "71824_slide2_SP206535_402641A",
                  "71824_slide2_SP206535_402641A",
                  "101124_slide1_10858",
                  "102224_slide1_34774"),
  geojson     = c("06_21_2024_slide1_259131A.geojson",
                  "06_21_2024_slide1_CS207490.geojson",
                  "06_21_2024_slide2_31361C.geojson",
                  "07_16_2024_slide1_392361A.geojson",
                  "07_16_2024_slide1_SP2223235.geojson",
                  "07_16_2024_slide2_184419A.geojson",
                  "07_18_2024_slide1_38979.geojson",
                  "07_18_2024_slide2_402641A.geojson",
                  "07_18_2024_slide2_SP206535.geojson",
                  "DLBCL_10858.geojson",
                  "DLBCL_34774.geojson"),
  stringsAsFactors = FALSE
)
samples$patient <- unname(accession_to_patient[samples$accession])

count_matrices_list <- samples$count_dir
JSON_list           <- samples$geojson
sample.name         <- samples$patient

for (i in 1:length(JSON_list)){
  i <- 1
  print(sample.name[i])
  print(JSON_list[i])
  print(count_matrices_list[i])
  
  dir.create(file.path(codex_path, "images", sample.name[i]), recursive = TRUE, showWarnings = FALSE)
  
  #Read in channel names
  channel_names <- read_csv(file = "./MarkerList.txt", col_names = FALSE)
  
  #Specify the directory of your data
  # file_path <- paste0(data_path, "count_matrices/", count_matrices_list[i], "/combined_markers.csv")
  file_path <- file.path(data_path, count_matrices_list[i], "combined_markers.csv")
  gc_csd_raw <- read_csv(file.path(file_path))
  print(dim(gc_csd_raw))
  print(colnames(gc_csd_raw))
  
  #renaming stuff, but the segmentation output is a bit different.
  colnames(gc_csd_raw)[5:60] <- channel_names$X1
  colnames(gc_csd_raw)[2] <- "Size"
  colnames(gc_csd_raw)[1] <- "CellID"
  colnames(gc_csd_raw)[3:4] <- c("x.coord", "y.coord")
  gc_csd_raw[,3:60] <- gc_csd_raw[,3:60] / gc_csd_raw$Size # normalize by size.
  
  # check cell size hist
  pdf(file.path(codex_path, "images", sample.name[i], paste0(sample.name[i], "_check_cell_size_hist.pdf")))
  hist(gc_csd_raw$Size[2:length(gc_csd_raw$Size)], breaks = 500)
  dev.off()
  
  #Check QC for DAPI
  topq <- quantile(gc_csd_raw$DAPI, 0.995)
  bottomq <- quantile(gc_csd_raw$DAPI, 0.001) ## change for Patient 5 to 0.1
  checkQC <- ggplot(data = gc_csd_raw, aes(x=`DAPI`)) + geom_histogram(bins = 2000)+
    geom_vline(xintercept = c(topq, bottomq),
               col = c('red', 'blue'),
               lwd = 1,
               linetype = 'dashed')
  
  pdf(file.path(codex_path, "images", sample.name[i], paste0(sample.name[i], "_check_QC_DAPI.pdf")))
  print(checkQC)
  dev.off()
  
  # for 8bit you use >10 and <250
  print(c(bottomq, topq))
  print(nrow(gc_csd_raw))
  
  gc_csd_raw <- gc_csd_raw %>% dplyr::filter(`DAPI`> bottomq, `DAPI`< topq) #quantile values
  gc_csd_raw <- gc_csd_raw %>% dplyr::filter(CellID!=0)
  print(nrow(gc_csd_raw))
  
  
  #View Image
  
  png(file.path(codex_path, "images", sample.name[i], paste0(sample.name[i], "_check_coords_image_2.png")))
  plot <- ggplot(gc_csd_raw, aes(x=`x.coord`, y = -`y.coord`)) + geom_point(alpha = 0.01)
  print(rasterize(plot, layers = 'Point', dpi=300))
  dev.off()
  
  png(file.path(codex_path, "images", sample.name[i], paste0(sample.name[i], "_check_coords_image_2_subset.png")))
  plot <- ggplot(gc_csd_raw, aes(x=`x.coord`, y = -`y.coord`)) + geom_point(alpha = 0.5) +
    coord_cartesian(xlim = c(7500, 10500),
                    ylim = c(-11500, -8000))
  print(rasterize(plot, layers = 'Point', dpi=300))
  dev.off()
  
  
  #Read in QuPath annotations (for polygon) -----
  polygon_list <- c()
  geojson_data <- fromJSON(file.path(data_path, JSON_list[i]))
  features <- geojson_data$features
  num_features <- dim(features)[1]
  for (j in 1:num_features) {
    coordinates <- features$geometry$coordinates[[j]]
    x_coords <- NULL
    y_coords <- NULL
    if (features$geometry$type[j] == "MultiPolygon"){
      print("MultiPolygon")
      num_poly <- length(coordinates)
      for (k in 1:num_poly){
        coordinates_sub <- coordinates[[k]]
        x_coords <- coordinates_sub[,,1]
        y_coords <- coordinates_sub[,,2]
        
        polygon_data <- data.frame(
          x = x_coords,
          y = y_coords
        )
        polygon_list[[length(polygon_list) + 1]] <- polygon_data
      }
      
    }
    else {
      if (length(dim(coordinates)) == 2){
        x_coords <- coordinates[,1]
        y_coords <- coordinates[,2]
      } else {
        x_coords <- coordinates[,,1]
        y_coords <- coordinates[,,2]
      }
      
      polygon_data <- data.frame(
        x = x_coords,
        y = y_coords
      )
      polygon_list[[length(polygon_list) + 1]] <- polygon_data
    }
    
  }
  
  print(length(polygon_list))
  #Function to reorder the polygon, when needed will go in circle
  reorder_polygon <- function(coords_matrix) {
    centroid <- colMeans(coords_matrix)
    #angles <- atan2(coords_matrix[, 2] - centroid[2], coords_matrix[, 1] - centroid[1])
    #sorted_coords <- coords_matrix[order(angles), ]
    sorted_coords <- coords_matrix
    if (!identical(sorted_coords[1, ], sorted_coords[nrow(sorted_coords), ])) {
      sorted_coords <- rbind(sorted_coords, sorted_coords[1, ])
    }
    return(sorted_coords)
  }
  
  #Apply to image
  combined_filtered_data <- data.frame()
  # Loop through each polygon and filter the data frame
  for (k in 1:length(polygon_list)) {
    print(k)
    polygon_coords <- polygon_list[[k]]
    # Check if the polygon is closed (first and last point are the same)
    polygon_1_list <- list(polygon_coords)
    polygon_1_list[[1]] = data.matrix(reorder_polygon(polygon_1_list[[1]]))
    x=polygon_1_list[[1]][,1]
    y=polygon_1_list[[1]][,2]
    plot(c(min(x), max(x)), c(min(-y), max(-y)), type = "n", xlab = "X", ylab = "Y")
    polygon(x, -y, col = "blue")
    polygon <- st_polygon(polygon_1_list)
    gc_csd_raw_sf <- st_as_sf(gc_csd_raw, coords = c("x.coord", "y.coord"))
    filtered_data <- st_intersection(gc_csd_raw_sf, polygon)
    # Append the filtered data to the combined_filtered_data data frame
    combined_filtered_data <- rbind(combined_filtered_data, as.data.frame(filtered_data))
  }
  
  filtered_coordinates <- st_coordinates(combined_filtered_data$geometry)
  combined_filtered_data$x.coord <- filtered_coordinates[, "X"]
  combined_filtered_data$y.coord <- filtered_coordinates[, "Y"]
  combined_filtered_data$geometry <- NULL
  head(combined_filtered_data)
  print(nrow(combined_filtered_data))
  
  #Inspect output
  png(file.path(codex_path, "images", sample.name[i], paste0(sample.name[i], "_check_polygon_image_DAPI_filter.png")))
  plot <- ggplot(combined_filtered_data, aes(x=`x.coord`, y = -`y.coord`)) + geom_point(alpha = 0.01)
  print(rasterize(plot, layers = 'Point', dpi=300))
  dev.off()
  
  #save
  data_path <- file.path(codex_root, "counts_raw")
  dir.create(file.path(data_path, "polygon_filtered", sample.name[i]), recursive = TRUE, showWarnings = FALSE)
  write.csv(combined_filtered_data,
            file.path(data_path, "polygon_filtered", sample.name[i],
                      paste0(sample.name[i], "_polygon_cell_gene.csv")))
  
  #Save output to variable and continue
  gc_csd_raw <- combined_filtered_data
  
  # end of polygon selection #-------
  
  #Remove those cells which are positive for everything and negative for everything
  gc_csd_raw <- tibble(gc_csd_raw) %>% rowwise() %>% dplyr::mutate(RowSum = sum(c_across(cols = ("CD19":"CD22"))))
  plots <- gc_csd_raw %>% ggplot(aes(x=RowSum)) + geom_histogram(bins = 200) + geom_vline(xintercept = quantile(gc_csd_raw$RowSum, probs = (0.99)))
  gc_csd_raw <- as.data.frame(gc_csd_raw)
  rownames(gc_csd_raw) <- gc_csd_raw$CellID
  
  
  png(file.path(codex_path, "images", sample.name[i], paste0(sample.name[i], "_check_rowsums.png")))
  print(plots)
  dev.off()
 
  # Remove the size and x/y coords and other erroneous values
  # (columns named after a file path are leftovers from the quantification step)
  gc_csd_raw_clean <- gc_csd_raw %>% select(-matches("/"))
  rownames(gc_csd_raw_clean) <- rownames(gc_csd_raw_clean)
  
  gc_csd2 <- gc_csd_raw %>% select(-starts_with('x.coord'), -contains("Empty"),
                                   -starts_with('y.coord'),-contains('size'),
                                   -contains("Std.Dev"), -contains("indexes"),
                                   -contains('RowSum'), -contains('Blank'),
                                   -matches("/"), -starts_with("Ch"))
  gc_csd2 <- gc_csd2 %>% select(-contains("Index"),-contains("CellID"), -matches("Python_Index"))
  
  # Create a Seurat object
  t_gc_csd = t(as.matrix(gc_csd2))
  colnames(t_gc_csd) <- rownames(gc_csd2)
  rownames(t_gc_csd) <- colnames(gc_csd2)
  
  
  seurat.obj <- CreateSeuratObject(t_gc_csd, project = paste0("LBCL_", sample.name[i]), assay="CODEX", meta.data = gc_csd_raw_clean)
  print(dim(seurat.obj))
  dir.create(file.path(codex_path, "seurat_objects"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(seurat.obj, file.path(codex_path, "seurat_objects",
                                paste0(sample.name[i], "_CODEX_unfiltered.RDS")))
}