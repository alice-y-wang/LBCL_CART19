# FUNCTIONS FOR TCR ANALYSIS by Martina Bonomi m.bonomi@unsw.edu.au

title_size = 16
labels_size = 14
text_size = 12

#' Plot clonal size distribution for all the samples using samples colors
#' 
#' @param tcr_df dataframes with TCR sequencing information for each patient
#' @return plot
#'
plot_cloneSize_distribution <- function(tcr_df, group){
  
  if(group == 'patient_id'){
    clonal_sizes = tcr_df %>%
      group_by(patient_id, CTaa) %>%
      summarize(clone_size = n()) %>%
      ungroup() %>%
      arrange(patient_id, desc(clone_size))
    clonal_size_freq = clonal_sizes %>%
      group_by(patient_id, clone_size) %>%
      summarize(frequency = n()) %>%
      ungroup() %>%
      mutate(patient_id = factor(patient_id, levels = unique(patient_id)))
  } else {
    clonal_sizes = tcr_df %>%
      group_by(patient_id, CTaa, !!sym(group), timepoint) %>%
      summarize(clone_size = n()) %>%
      ungroup() %>%
      arrange(patient_id, desc(clone_size))
    clonal_size_freq = clonal_sizes %>%
      group_by(patient_id, clone_size, !!sym(group), timepoint) %>%
      summarize(frequency = n()) %>%
      ungroup() %>%
      mutate(patient_id = factor(patient_id, levels = unique(patient_id)))
  }
  
  plot = ggplot(clonal_size_freq, aes(x = clone_size, y = frequency, color = timepoint, group = patient_id)) +
    geom_point(size = 4, alpha = 1) +
    geom_line(linewidth = 0.5, alpha = 0.5) +
    scale_color_manual(values = timepoints_colors_2, name = "Timepoint") +
    scale_x_log10(breaks = c(1, 10, 100, 1000), labels = scales::comma_format()) +
    scale_y_log10(breaks = c(1, 10, 100, 1000), labels = scales::comma_format()) +
    labs(title = "Clonal Size Distribution", x = "Clone Size", y = "Frequency") +
    theme(plot.title = element_text(face = "bold", size = title_size, hjust = 0.5),
          axis.text = element_text(size = text_size),
          axis.title = element_text(face='bold', size=labels_size),
          legend.title = element_text(face='bold',size=labels_size))
  
  return(plot)
}

#' Function to compute Shannon entropy
#' 
#' @param seq TCR repertoire from a sample
#' @return shannon entropy value
#' 
shannon_entropy <- function(seq) {
  freq = table(seq) / length(seq)
  entropy = -sum(freq * log(freq))
  return(entropy)
}

#' Function to compute normalised Shannon entropy
#' 
#' @param seq TCR repertoire from a sample
#' @return shannon entropy value
#' 
norm_shannon_entropy <- function(seq) {
  freq = table(seq) / length(seq)
  entropy = -sum(freq * log(freq))
  max_entropy = log(length(freq))  # log(number of unique clonotypes)
  if (max_entropy == 0) return(0)
  return(entropy / max_entropy)
}

#' Function to compute the Gini-Simpson Index
#' 
#' @param seq TCR repertoire from a sample
#' @return gini-simpson index value
#' 
gini_simpson <- function(seq) {
  freq = table(seq)
  prop = freq / sum(freq)
  gini_simpson_index = 1 - sum(prop^2)
  return(gini_simpson_index)
}

#' Function to compute the Pielou's Evenness
#' 
#' @param seq TCR repertoire from a sample
#' @return pielou's evenness value
#' 
pielous_evenness <- function(seq) {
  entropy = shannon_entropy(seq)
  S = length(unique(seq))
  evenness = entropy / log(S)
  return(evenness)
}

#' Function to compute the richness
#' 
#' @param seq TCR repertoire from a sample
#' @return richness value
#' 
richness <- function(seq) {
  rich = length(unique(seq))/length(seq)
  return(rich)
}

#' Function to compute the Pielou's Evenness
#' 
#' @param seq TCR repertoire from a sample
#' @return pielou's evenness value
#' 
PielouEvenness <- function(seq) {
  entropy = Entropy(table(seq), base = exp(1))
  S = length(unique(seq))
  evenness = entropy / log(S)
  return(evenness)
}

#' Calculate metrics for a sample
#' 
#' @param tcr_df dataframe with TCR repertoires of samples, must contain columns 'sample' and 'CTaa'
#' @param min_size minimum sample size to consider
#' @param n_iterations number of iterations to perform for subsampling big samples
#' @return dataframe with diversity measures (richness, shannon entropy, gini-simpson index, gini coefficient, pielou's evenness)
#' 
calculate_metrics <- function(tcr_df, min_size = 200, n_iterations = 1000){
  
  # Filter to samples with at least min_size
  sample_sizes = tcr_df %>% count(sample2)
  eligible_samples = sample_sizes %>% filter(n >= min_size) %>% pull(sample2)
  
  # Subset TCR data
  tcr_df = tcr_df %>% filter(sample2 %in% eligible_samples)
  
  # Determine the actual min sample size to use
  min_sample_size = tcr_df %>% group_by(sample2) %>% summarise(size = n()) %>% pull(size) %>% min()
  
  metrics_df_list = list()
  
  for (sample_i in unique(tcr_df$sample2)) {
    message("Calculating metrics for sample: ", sample_i)
    
    seq = tcr_df %>% filter(sample2 == sample_i) %>% pull(CTaa)
    patient = str_split(sample_i, "_", simplify = TRUE)[1]
    timepoint = str_split(sample_i, "_", simplify = TRUE)[2]

    vec.richness = numeric(n_iterations)
    vec.shannon = numeric(n_iterations)
    vec.norm_shannon = numeric(n_iterations)
    # vec.ginisimp = numeric(n_iterations)
    # vec.ginicoeff = numeric(n_iterations)
    vec.pielou = numeric(n_iterations)
    vec.inv_simpson = numeric(n_iterations)
    
    for (i in 1:n_iterations) {
      sub_seq = sample(seq, min_sample_size, replace = FALSE)
      tbl = table(sub_seq)
      p = tbl / sum(tbl)
      
      vec.richness[i] = specnumber(tbl)
      vec.shannon[i] = diversity(tbl, index = "shannon")
      vec.norm_shannon[i] = vec.shannon[i] / log(specnumber(tbl))
      # vec.ginisimp[i] = GiniSimpson(factor(sub_seq, levels = unique(sub_seq)))
      # vec.ginicoeff[i] = Gini(table(sub_seq))
      vec.pielou[i] = vec.shannon[i] / log(specnumber(tbl))
      vec.inv_simpson[i] = diversity(tbl, index = "invsimpson")
    }
    
    metrics_df_list[[sample_i]] = data.frame(
      sample = sample_i,
      patient_id = patient,
      timepoint = timepoint,
      richness = round(median(vec.richness), 6),
      shannon_entropy = median(vec.shannon),
      norm_entropy = median(vec.norm_shannon),
      # gini_simpson = median(vec.ginisimp),
      # gini_coefficient = median(vec.ginicoeff),
      evenness = median(vec.pielou),
      inv_simpson = median(vec.inv_simpson)
    )
  }
  
  metrics_df = bind_rows(metrics_df_list)
  return(metrics_df)
}

