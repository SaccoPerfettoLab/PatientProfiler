


#' Extract signaling-driven transcriptomics signature from communities
#'
#' This function analyzes transcriptomic data by performing ANOVA-Tukey on gene expression levels and then generating community-specific gene expression signatures
#'
#' @param base_path Character. Path to the directory containing community folders.
#' @param transcriptomics_file Character. Path to the transcriptomic data file.
#' @param output_dir Character. Path to the directory containing all function's results.
#' @param padj_thres Numeric. Adjusted p-value threshold for filtering significant results. Default: `0.01`.
#' @param mean_exp_clus_thres Numeric. Minimum mean expression threshold for gene filtering. Default: `0`.
#' @param max_val Integer. Maximum number of genes to include in each community-specific signature. Default: `50`.
#' @param cache Boolean. If TRUE uses Anova_result.tsv file in `output_dir` if available. Default: `FALSE`.
#'
#' @return A data frame containing the results of the transcriptomic signature for each community.
#' The function also saves a csv file with the anova results and generates a "Signatures" directory where filtered signature files for each community are stored.
#'
#' @export
#'
#' @examples
#' extract_signatures(base_path = "./vignette/Communities/output_communities/",
#' transcriptomics_file = "./PatientProfiler_processed_input/Transcriptomics_updated.tsv",
#' output_dir = "./Anova_result.tsv",
#' padj_thres = 0.01,
#' mean_exp_clus_thres = 0,
#' max_val = 50,
#' cache = F)
#'
#' extract_signatures(base_path = "./vignette/Communities/output_communities/",
#' transcriptomics_file = "./PatientProfiler_processed_input/Transcriptomics_updated.tsv",
#' output_dir = "./Anova_result.tsv",
#' padj_thres = 0.01,
#' mean_exp_clus_thres = 0,
#' max_val = 10,
#' cache = T)
#'


extract_signatures <- function(base_path,
                               transcriptomics_file,
                               output_dir,
                               padj_thres = 0.01,
                               mean_exp_clus_thres = 0,
                               max_val = 50,
                               cache = FALSE) {
  
  message("Warning: This function may take a long time to execute.")
  
  message("Creating output directories...")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!cache) {
    output_file <- file.path(output_dir, "Anova_result.tsv")
    signature_dir <- file.path(output_dir, "Signatures")
    dir.create(signature_dir, showWarnings = FALSE)
    
    message("Scanning community folders...")
    community_folders <- list.dirs(base_path, recursive = FALSE)
    stratification_table <- data.frame()
    
    for (community_folder in community_folders) {
      community_id <- basename(community_folder)
      patients_file <- file.path(community_folder, "patients.txt")
      
      if (file.exists(patients_file)) {
        message(paste("Processing", patients_file))
        community_data <- read.table(patients_file, header = FALSE, stringsAsFactors = FALSE)
        colnames(community_data)[1] <- "Patient_ID"
        community_data$community <- community_id
        stratification_table <- rbind(stratification_table, community_data)
      } else {
        warning(paste("File not found:", patients_file))
      }
    }
    
    message("Loading transcriptomics data...")
    Transcriptomics_patients <- readr::read_tsv(transcriptomics_file) %>%
      pivot_longer(cols = -1, names_to = "Patient_ID", values_to = "value")
    Transcriptomics_patients <- left_join(Transcriptomics_patients, stratification_table, by = 'Patient_ID')
    
    message("Performing ANOVA analysis...")
    final_results <- perform_anova(Transcriptomics_patients, stratification_table)
    write_tsv(final_results, output_file)
  } else {
    message("Reading previously generated Anova result...")
    if (file.exists(file.path(output_dir, "Anova_result.tsv"))) {
      final_results <- readr::read_tsv(file.path(output_dir, "Anova_result.tsv"), show_col_types = FALSE)
      signature_dir <- file.path(output_dir, "Signatures")
      dir.create(signature_dir, showWarnings = FALSE)
    }
  }
  
  message("Computing global mean of cluster-level diff means...")
  
  cluster_means <- c()
  for (c_i in unique(final_results$cluster)) {
    cluster_data <- final_results %>%
      filter(cluster == paste0("community_", c_i),
             `p adj` < padj_thres,
             mean_exp_clus > mean_exp_clus_thres,
             diff > 0) %>%
      arrange(desc(diff))
    
    top_cluster_data <- cluster_data[1:min(max_val, nrow(cluster_data)), ]
    mean_diff <- mean(top_cluster_data$diff, na.rm = TRUE)
    cluster_means <- c(cluster_means, mean_diff)
    message(paste0("community_", c_i, " mean diff: ", round(mean_diff, 3)))
  }
  
  global_mean_diff <- mean(cluster_means, na.rm = TRUE)
  message(paste("Global mean of cluster means:", round(global_mean_diff, 3)))
  
  message("Extracting signatures with adaptive diff thresholds...")
  for (cluster_i in unique(final_results$cluster)) {
    message(paste("Processing cluster:", cluster_i))
    
    cluster_data <- final_results %>%
      filter(cluster == cluster_i,
             `p adj` < padj_thres,
             mean_exp_clus > mean_exp_clus_thres,
             diff > 0) %>%
      arrange(desc(diff))
    
    if (nrow(cluster_data) == 0) {
      warning(paste("No entries passed initial filter for cluster", cluster_i))
      next
    }
    
    top_entries <- cluster_data[1:min(max_val, nrow(cluster_data)), ]
    mean_diff_cluster <- mean(top_entries$diff, na.rm = TRUE)
    
    diff_threshold <- ifelse(mean_diff_cluster < global_mean_diff, 0.5, 0.7)
    
    message(paste("Cluster", cluster_i,
                  "- mean diff:", round(mean_diff_cluster, 3),
                  "- using diff threshold:", diff_threshold))
    
    final_results_filtered <- cluster_data %>%
      filter(diff >= diff_threshold)
    
    output_signature_file <- file.path(signature_dir, paste0('signature_', cluster_i, '.tsv'))
    write_tsv(final_results_filtered[, 5:6], output_signature_file)
    message(paste("Signature saved to", output_signature_file))
  }
  
  message("Process completed.")
  return(final_results)
}


