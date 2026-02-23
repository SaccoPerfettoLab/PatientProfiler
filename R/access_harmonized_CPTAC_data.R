#'access_harmonized_CPTAC_data
#'
#'This function provides access to harmonized CPTAC data for 10 tumor types and 3 omics per tumor type. The data has been harmonized using the
#'PatientProfiler function *omics_update*. Users can access the data to use them in the other functions for protein activity inference and
#'for network creation. It downloads the selected datasets locally (by default in the current working directory) and creates variables with the
#'dataframe(s) selected by the user.
#'
#'@param tumors vector of strings with tumor names. The user can choose from these names: "Brca", "Ccrcc", "Coad", "Gbm", "Hnscc", "Lscc",
#'"Luad", "Ov", "Pdac", "Ucec".
#'
#'@param data_types vector of strings containing the data types to extract. The options could be "phospho", "prot", "transc" and "mut"
#'(respectively phosphoproteomics, proteomics, transcriptomics and somatic mutations)
#'
#'@param dest_dir string. Destination directory where files will be downloaded. Default is the current working directory (getwd()).
#'
#'@param overwrite logical. If TRUE, re-download files even if they already exist locally. Default is FALSE.
#'
#'@examples
#' # Extract Phosphoproteomic and Proteomic dataframes for Brca and Ccrcc
#' access_harmonized_CPTAC_data(tumors = c("Brca", "Ccrcc"), data_types = c("phospho", "prot", "mut"))
#'
#' # Six local variables, Brca_phospho, Brca_prot, Brca_mut, Ccrcc_phospho, Ccrcc_prot and Ccrcc_mut will be generated.
#' # The corresponding TSV files will be downloaded to the working directory (unless dest_dir is specified).
#'
#' @export
access_harmonized_CPTAC_data <- function(tumors,
                                         data_types,
                                         dest_dir = getwd(),
                                         overwrite = FALSE) {
  
  base_url <- "https://perfettolab.bio.uniroma1.it/PerfettoLabData/PatientProfiler/CPTAC/"
  
  # Ensure destination directory exists
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  for (tumor in tumors) {
    
    tumor_dir <- file.path(base_url, tumor)
    
    for (data_type in data_types) {
      
      if (data_type == "mut") {
        
        file_name <- paste0(tumor, "_mut_matrix.tsv")
        remote_url <- file.path(tumor_dir, file_name)
        local_path <- file.path(dest_dir, file_name)
        
        # Download to disk (unless already present and overwrite=FALSE)
        if (!file.exists(local_path) || isTRUE(overwrite)) {
          resp <- httr::GET(remote_url, httr::write_disk(local_path, overwrite = TRUE))
          if (httr::status_code(resp) != 200) {
            # Clean up partial file if download failed
            if (file.exists(local_path)) unlink(local_path)
            warning(paste("Unable to download mutation file:", file_name))
            next
          } else {
            cat("Downloaded mutation file:", file_name, "->", local_path, "\n")
          }
        } else {
          cat("Using local mutation file:", file_name, "->", local_path, "\n")
        }
        
        # Read from local disk
        mut_data <- readr::read_tsv(local_path, show_col_types = FALSE)
        
        # Create global variable
        var_name <- paste0(tumor, "_mut")
        assign(var_name, mut_data, envir = .GlobalEnv)
        cat("Read mutation data into:", var_name, "\n")
        
      } else if (data_type %in% c("phospho", "prot", "transc")) {
        
        file_name <- paste0(tumor, "_", data_type, "_updated.tsv")
        remote_url <- file.path(tumor_dir, file_name)
        local_path <- file.path(dest_dir, file_name)
        
        # Download to disk (unless already present and overwrite=FALSE)
        if (!file.exists(local_path) || isTRUE(overwrite)) {
          resp <- httr::GET(remote_url, httr::write_disk(local_path, overwrite = TRUE))
          if (httr::status_code(resp) != 200) {
            # Clean up partial file if download failed
            if (file.exists(local_path)) unlink(local_path)
            warning(paste("Unable to download omic file:", file_name))
            next
          } else {
            cat("Downloaded omic file:", file_name, "->", local_path, "\n")
          }
        } else {
          cat("Using local omic file:", file_name, "->", local_path, "\n")
        }
        
        # Read from local disk
        omic_data <- readr::read_tsv(local_path, show_col_types = FALSE)
        
        # Create global variable
        var_name <- paste0(tumor, "_", data_type)
        assign(var_name, omic_data, envir = .GlobalEnv)
        cat("Read omic data into:", var_name, "\n")
        
      } else {
        warning(paste("Unknown data type:", data_type))
      }
    }
  }
  
  invisible(NULL)
}