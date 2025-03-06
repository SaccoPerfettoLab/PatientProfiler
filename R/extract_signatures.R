

#' Extract signaling-driven transcriptomics signature from communities
#' 
#' This function analyzes transcriptomic data by performing ANOVA-Tukey on gene expression levels and then generating community-specific gene expression signatures
#'
#' @param base_path Character. Path to the directory containing community folders.
#' @param transcriptomics_file Character. Path to the transcriptomic data file. 
#' @param output_dir Character. Path to the directory containing all function's results. 
#' @param padj_thres Numeric. Adjusted p-value threshold for filtering significant results. Default: `0.01`.
#' @param diff_thres Numeric. Minimum difference threshold for gene expression filtering. Default: `0.7`.
#' @param mean_exp_clus_thres Numeric. Minimum mean expression threshold for gene filtering. Default: `0`.
#' @param max_val Integer. Maximum number of genes to include in each community-specific signature. Default: `50`.
#'
#' @return A data frame containing the results of the transcriptomic signature for each community.
#' The function also saves a csv file with the anova results and generates a "Signatures" directory where filtered signature files for each community are stored.
#'
#' @export
#'
#' @examples
#' extract_signatures <- function(base_path = "./vignette/Communities/output_communities/",          
#' transcriptomics_file = "./PatientProfiler_processed_input/Transcriptomics_updated.csv", 
#' output_dir = "./Anova_result.csv",
#' padj_thres = 0.01,
#' diff_thres = 0.5,
#' mean_exp_clus_thres = 0,
#' max_val = 50)



extract_signatures <- function(base_path,          
                               transcriptomics_file, 
                               output_dir,
                               padj_thres = 0.01,
                               diff_thres = 0.7,
                               mean_exp_clus_thres = 0,
                               max_val = 50) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  output_file <- file.path(output_dir, "Anova_result.csv")
  
  signature_dir <- file.path(output_dir, "Signatures")
  dir.create(signature_dir, showWarnings = FALSE)
  
  community_folders <- list.dirs(base_path, recursive = FALSE)
  stratification_table <- data.frame()
  
  for (community_folder in community_folders) {
    community_id <- basename(community_folder)
    patients_file <- file.path(community_folder, "patients.txt") 
    
    if (file.exists(patients_file)) {
      community_data <- read.csv(patients_file, header = TRUE, stringsAsFactors = FALSE)
      colnames(community_data)[1] <- "Patient_ID"  
      community_data$community <- community_id  
      stratification_table <- rbind(stratification_table, community_data)
    } else {
      warning(paste("File not found:", patients_file))
    }
  }
  
  Transcriptomics_patients <- readr::read_csv(transcriptomics_file) %>%
    pivot_longer(cols = -1, names_to = "Patient_ID", values_to = "value")
  Transcriptomics_patients <- left_join(Transcriptomics_patients, stratification_table, by = 'Patient_ID')
  
  final_results <- perform_anova(Transcriptomics_patients, stratification_table)
  write_csv(final_results, output_file)
  
  for (cluster_i in unique(final_results$cluster)) {
    final_results_filtered <- final_results %>%
      dplyr::arrange(`p adj`) %>%
      dplyr::filter(cluster == cluster_i &
                      `p adj` < padj_thres & 
                      diff > diff_thres & 
                      mean_exp_clus > mean_exp_clus_thres)
    
    if (nrow(final_results_filtered) > max_val) {
      final_results_filtered <- final_results_filtered[1:max_val,]
    }
    
    write_csv(final_results_filtered[, 5:6], file.path(signature_dir, paste0('signature_', cluster_i, '.csv')))
  }
  
  return(final_results)
}
