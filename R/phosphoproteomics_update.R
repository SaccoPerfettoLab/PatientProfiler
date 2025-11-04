#'  phosphoproteomics_update
#'
#' This function processes phosphoproteomics data by performing several steps of data cleaning,
#' peptide modification, and enrichment of sequence information. It is designed to handle datasets
#' with phosphopeptides, modifying peptide sequences, adding multiplicity information, and computing
#' sequence windows around the phosphopeptides.
#'
#'
#' @param df a data frame containing the phosphoproteomics data.
#' @param threshold numeric, it specifies the threshold percentage of NAs in each row of the dataframe. If there are >= threshold percentage of
#'                  NAs in a row, the row will be removed.
#'@param m number of imputation with mice.
#'@param collapse how to pull dataframes derived from the imputation with mice (default: "median").
#' @param seq_len_i an integer specifying the length of the sequence window to be considered around the phosphopeptide.
#' @param uniprot_idx optional integer indicating the column index for UNIPROT IDs to be added at the uniprot column retrieved with AnnotationDbi (default is NULL).
#' @param peptide_col_name (Optional) A string representing the name of the peptide column in the data frame.
#' If not provided, the function assumes no peptide column exists.
#' @param imp_method string, the method to use for imputation (default: "pmm", but you can choose between "pmm", "norm", "norm.nob",
#'                   "regression", "ri", "logreg", "polyreg", "predictive", "polr", "sample", "cart", "knn", "rf").
#' @param output_dir a string indicating the updated output folder.
#'
#' @return A data frame with the processed phosphoproteomics data, including modified peptide sequences,
#' updated amino acid positions, and sequence windows. If no peptide column is provided, the function
#' directly computes amino acid positions.
#'
#' @examples
#' # Example usage:
#' updated_phospho_df <- update_phospho(phospho_data, site_col = 3, gn_idx = 2, seq_len_i = 7, peptide_col_name = "Peptide")
#' head(updated_phospho_df)



phosphoproteomics_update <- function(df_pho,
                                     threshold = 80,
                                     sw_len = 7,
                                     uniprot_idx = NULL,
                                     pep_col_name = NULL,
                                     imp_method = "pmm",
                                     m = 5,
                                     collapse = "median",
                                     group = group,
                                     zscore = TRUE,
                                     zmethod = "column",
                                     metric = "median",
                                     output_dir){

  message("Phosphoproteomics data: removing rows with excessive missing values")
  df_pho <- remove_nas(df_pho, threshold)
  
  message("Phosphoproteomics data: updating phosphorylation site informations")
  df_pho_update <<- update_phospho(df = df_pho,site_col = 2,gn_idx = 1,seq_len_i = 7, uniprot_idx, peptide_col_name = pep_col_name)
  message("Done!")

  message("Phosphoproteomics data: removing duplicates")

  df_pho_update <- df_pho_update %>% distinct(UNIPROT, aminoacid, position,.keep_all = TRUE)
  df_pho_clean <<- remove_duplicates_phosphoproteomics(df_pho_update)

  readr::write_tsv(df_pho_clean, paste0(output_dir,"/","Phosphoproteomics_clean.tsv")) ####
  message("Done!")
  
  # Identify numeric columns, but exclude "position"
  num_cols <- which(sapply(df_pho_clean, is.numeric))
  exclude_cols <- which(colnames(df_pho_clean) %in% c("position"))
  
  valid_num_cols <- setdiff(num_cols, exclude_cols)
  
  if (length(valid_num_cols) == 0) stop("No valid numeric columns found for imputation.")
  
  # Use the first remaining numeric column as start_column
  start_column <- min(valid_num_cols)
  
  
  message("Phosphoproteomics data: missing values imputation")
  df_pho_imputed <<- impute_proteomics(df_pho_clean,
                                       start_column = start_column,
                                       imp_method, 
                                       m,
                                       maxit = 5,
                                       seed = 103,
                                       collapse = collapse,
                                       preserve_observed = TRUE,
                                       clean_patient_names = TRUE, group)
  

  df_pho_imputed <- as.data.frame(df_pho_imputed$imputed_df)
  readr::write_tsv(df_pho_imputed, paste0(output_dir,"/","Phosphoproteomics_imputed.tsv")) ####
  
  message("Done!")

  if (zscore) {
    message("Phosphoproteomics data: computing Z-scores")
    if ("sequence_window" %in% colnames(df_pho_imputed)) {
      df_pho_numeric <- df_pho_imputed %>%
        dplyr::select(-UNIPROT, -sequence_window, -position, -aminoacid, -gene_name)

      metadata_columns <- df_pho_imputed %>%
        dplyr::select(UNIPROT, sequence_window, position, aminoacid, gene_name)
    } else {
      df_pho_numeric <- df_pho_imputed %>%
        dplyr::select(-UNIPROT, -position, -aminoacid, -gene_name)

      metadata_columns <- df_pho_imputed %>%
        dplyr::select(UNIPROT, position, aminoacid, gene_name)

    }

    df_pho_matrix <- as.matrix(df_pho_numeric)
    df_pho_zscore <<- compute_zscore(df_pho_matrix, zmethod, metric)
    df_pho_zscore <<- dplyr::mutate_all(as.data.frame(df_pho_zscore), as.numeric)

    df_pho_zscore <<- cbind(metadata_columns, df_pho_zscore)
    return(df_pho_zscore)
    readr::write_tsv(df_pho_zscore, paste0(output_dir,"/","Phosphoproteomics_updated.tsv"))

    message("Done!")
  } else {
    return(df_pho_imputed)
  }

}
