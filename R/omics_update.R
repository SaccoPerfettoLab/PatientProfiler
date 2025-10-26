#' omics_update
#'
#' This function updates transcriptomics, proteomics, and phosphoproteomics datasets by performing
#' preprocessing steps such as filtering, imputation, and Z-score normalization.
#' Dataframes must contain the following required columns:
#' - df_tr: 'gene_name' and sample columns.
#' - df_pr: 'gene_name' and sample columns.
#' - df_ph: 'gene_name', 'Site' (containing phosphorylation sites in the format: T161), and 'Peptide'
#'   (containing the peptide sequence; the name and the presence of the peptide column is optional).
#'
#' @param df_tr dataframe containing transcriptomics data.
#' @param df_pr dataframe containing proteomics data.
#' @param df_ph dataframe containing phosphoproteomics data.
#' @param threshold numeric, it specifies the threshold percentage of NAs in each row of the dataframe. If there are >= threshold percentage of
#'                  NAs in a row, the row will be removed.
#'@param m number of imputation with mice.
#'@param collapse how to pull dataframes derived from the imputation with mice (default: "median").
#' @param sw_len integer, the sequence window length. It could be 7 or 15 (default: 7).
#' @param uniprot_idx optional integer indicating the column index for UNIPROT IDs to be added at the uniprot column retrieved with AnnotationDbi (default is NULL).
#' @param pep_col_name string, it indicates the name of the peptide sequence column (if it's present) in your phosphoproteomic dataframe (default = "Peptide")
#' @param imp_method string, the method to use for imputation (default: "pmm", but you can choose between "pmm", "norm", "norm.nob",
#'                   "regression", "ri", "logreg", "polyreg", "predictive", "polr", "sample", "cart", "knn", "rf").
#' @param zscore logical, whether to perform Z-score normalization (default: "TRUE").
#' @param zmethod string, specifies whether Z-score normalization is performed by "row" or "column" (default: "column").
#' @param metric string, the centering metric for Z-score normalization. Options are "median" (default) or "mean".
#' @param output_dir a string indicating the updated output folder.
#'
#' @return This function updates the following global variables:
#'
#' transctriptomics_updated -> A data frame containing the updated transcriptomics data.
#' proteomics_updated -> A data frame containing the updated proteomics data.
#' phosphoproteomics_updated -> A data frame containing the updated phosphoproteomics data.
#'
#' @examples
#' transcriptomics_data <- data.frame(matrix(rnorm(100), nrow = 10))
#' proteomics_data <- data.frame(matrix(rnorm(100), nrow = 10))
#' phosphoproteomics_data <- data.frame(matrix(rnorm(100), nrow = 10))
#'
#' omics_update(
#'   df_tr = transcriptomics_data,
#'   df_pr = proteomics_data,
#'   df_ph = phosphoproteomics_data,
#'   sw_len = 7,
#'   uniprot_idx = NULL,
#'   pep_col_name = NULL,
#'   imp_method = "pmm",
#'   zscore = TRUE,
#'   zmethod = "column",
#'   metric = "median"
#' )
#'
#' @export



omics_update <- function(df_tr = NULL,
                         df_pr = NULL,
                         df_ph = NULL,
                         threshold = 80,
                         m= 5,
                         collapse = "median",
                         sw_len = 7,
                         uniprot_idx = NULL,
                         pep_col_name = NULL,
                         imp_method = "pmm",
                         zscore = "TRUE",
                         zmethod = "column",
                         metric = "median",
                         output_dir = "PatientProfiler_processed_input") {
  t0 <- Sys.time()
  
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }

  if (!is.null(sw_len) && !sw_len %in% c(7, 15)) {
    stop("sw_len parameter has to be 7, 15 or NULL (if you do not have the phosphoproteomics dataframe")
  }

  if (!is.null(df_tr)) {
    message("Transcriptomics update started..")

    transcriptomics_updated <<- transcriptomics_update(df_tr,threshold,zscore,zmethod,metric,output_dir)

    readr::write_tsv(transcriptomics_updated, paste0(output_dir,"/","Transcriptomics_updated.tsv"))
    message("Transcriptomics update complete!")
  }

  if(!is.null(df_pr)){
    message("Proteomics update started..")
    proteomics_updated <<- proteomics_update(df_pr,threshold,uniprot_idx,imp_method,m,collapse,zscore,zmethod,metric, output_dir)
    readr::write_tsv(proteomics_updated, paste0(output_dir,"/","Proteomics_updated.tsv"))

    message("Proteomics update complete!")
  }

  if(!is.null(df_ph)){
    message("Phosphoproteomics update started..")
    phosphoproteomics_updated <<- phosphoproteomics_update(df_ph,
                                                           threshold    = threshold,
                                                           sw_len       = sw_len,        
                                                           uniprot_idx  = uniprot_idx,
                                                           pep_col_name = pep_col_name,  
                                                           imp_method   = imp_method,
                                                           m            = m,
                                                           collapse     = collapse,      
                                                           zscore       = zscore,
                                                           zmethod      = zmethod,
                                                           metric       = metric,
                                                           output_dir   = output_dir)

    readr::write_tsv(phosphoproteomics_updated, paste0(output_dir,"/","Phosphoproteomics_updated.tsv"))

    message("Phosphoproteomics update complete!")
  }
  
  t1 <- Sys.time()
  elapsed_secs <- as.numeric(difftime(t1, t0, units = "secs"))
  
  if (elapsed_secs < 60) {
    elapsed_pretty <- paste0(round(elapsed_secs, 2), " seconds")
  } else {
    elapsed_pretty <- paste0(round(elapsed_secs / 60, 2), " minutes")
  }
  
  message(paste0("\n Update completed in ", elapsed_pretty, "."))
  
}



