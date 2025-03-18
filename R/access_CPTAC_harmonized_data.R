#'access_harmonized_CPTAC_data
#'
#'This function provides access to harmonized CPTAC data for 10 tumor types and 3 omics per tumor type. The data has been harmonized using the
#'PatientProfiler function *omics_update*. Users can access the data to use them in the other functions for protein activity inference and
#'for network creation. It creates local variables with the dataframe(s) selected by the user.

#'@param tumors vector of strings with tumor names. The user can choose from these names: "Brca", "Ccrcc", "Coad", "Gbm", "Hnscc", "Lscc",
#'"Luad", "Ov", "Pdac", "Ucec".
#'
#'@param data_types vector of strings containing the data types to extract. The options could be "phospho", "prot", "transc" and "mut"
#'(respectively phosphoproteomics, proteomics, transcriptomics and somatic mutations)
#'
#'
#'@examples
#' # Extract Phosphoproteomic and Proteomic dataframes for Brca and Ccrcc
#' access_harmonized_CPTAC_data(tumors = c("Brca", "Ccrcc"), data_types = c("phospho", "prot", "mut"))
#'
#' # Six local variables, Brca_phospho, Brca_prot, Brca_mut, Ccrcc_phospho, Ccrcc_prot and Ccrcc_mut will be generated.
#'
#' @export

access_harmonized_CPTAC_data <- function(tumors, data_types) {
  
  base_url <- "http://userver.bio.uniroma1.it/apps/CPTAC"
  
  for (tumor in tumors) {
    
    tumor_dir <- file.path(base_url, tumor)
    
    for (data_type in data_types) {
      
      if (data_type == "mut") {
        mut_file_name <- paste0(tumor, "_mut_matrix.tsv")
        mut_file_path <- file.path(tumor_dir, mut_file_name)
        
        response <- httr::GET(mut_file_path)
        if (status_code(response) == 200) {
          mut_data <- read.tsv(text = content(response, "text"))
          
          var_name <- paste0(tumor, "_mut")
          assign(var_name, mut_data, envir = .GlobalEnv)
          cat("Read mutation data:", mut_file_name, "\n")
        } else {
          warning(paste("Unable to download mutation file:", mut_file_name))
        }
        
      } else if (data_type %in% c("phospho", "prot", "transc")) {
        file_name <- paste0(tumor, "_", data_type, "_updated.tsv")
        file_path <- file.path(tumor_dir, file_name)
        
        response <- httr::GET(file_path)
        if (status_code(response) == 200) {
          temp_file <- tempfile(fileext = ".tsv")
          writeBin(content(response, "raw"), temp_file)
          
          omic_data <- readr::read_tsv(temp_file)
          
          var_name <- paste0(tumor, "_", data_type)
          assign(var_name, omic_data, envir = .GlobalEnv)
          cat("Read omic data:", file_name, "\n")
        } else {
          warning(paste("Unable to download omic file:", file_name))
        }
        
      } else {
        warning(paste("Unknown data type:", data_type))
      }
    }
  }
}
