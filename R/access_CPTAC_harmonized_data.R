#'access_harmonized_CPTAC_data
#'
#'This function provides access to harmonized CPTAC data for 10 tumor types and 3 omics per tumor type. The data has been harmonized using the
#'PatientProfiler function *omics_update*. Users can access the data to use them in the other functions for protein activity inference and
#'for network creation. It creates local variables with the dataframe(s) selected by the user.

#'@param tumors vector of strings with tumor names. The user can choose from these names: "Brca", "Ccrcc", "Coad", "Gbm", "Hnscc", "Lscc",
#'"Luad", "Ov", "Pdac", "Ucec".
#'
#'@param omics vector of strings containing the omics types to extract. The options could be "phospho", "prot" and "transc"
#'(respectively phosphoproteomics, proteomics and transcriptomics)
#'
#'
#'@examples
#' # Extract Phosphoproteomic and Proteomic dataframes for Brca and Ccrcc
#' access_harmonized_CPTAC_data(tumors = c("Brca", "Ccrcc"), omics = c("phospho", "prot"))
#'
#' # Four local variables, Brca_phospho, Brca_prot, Ccrcc_phospho and Ccrcc_prot will be generated.
#'
#' @export

access_harmonized_CPTAC_data <- function(tumors, omics) {
  
  data <- readRDS("http://151.100.141.57/cptac/CPTAC_df.rds") # this link is not valid for now, we'll change it soon

  suffix <- c("phospho", "prot", "transc")

  if (!all(omics %in% suffix)) {
    stop("Omics has to be one or more between 'phospho', 'prot' and 'transc'")
  }

  for (tumor in tumors) {
    if (!tumor %in% names(data)) {
      warning(paste("Tumor not found:", tumor))
      next
    }
    for (omic in omics) {
      if (!omic %in% names(data[[tumor]])) {
        warning(paste("Omic not found for", tumor, ":", omic))
        next
      }
      local_var <- paste0(tumor, "_", omic)

      assign(local_var, data[[tumor]][[omic]], envir = .GlobalEnv)

    }
  }
}
