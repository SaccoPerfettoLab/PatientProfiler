#'  phosphoproteomics_update
#'
#' This function processes phosphoproteomics data by performing several steps of data cleaning,
#' peptide modification, and enrichment of sequence information. It is designed to handle datasets
#' with phosphopeptides, modifying peptide sequences, adding multiplicity information, and computing
#' sequence windows around the phosphopeptides.
#'
#'
#' @param df a data frame containing the phosphoproteomics data.
#' @param seq_len_i an integer specifying the length of the sequence window to be considered around the phosphopeptide.
#' @param peptide_col_name (Optional) A string representing the name of the peptide column in the data frame.
#' If not provided, the function assumes no peptide column exists.
#'
#' @return A data frame with the processed phosphoproteomics data, including modified peptide sequences,
#' updated amino acid positions, and sequence windows. If no peptide column is provided, the function
#' directly computes amino acid positions.
#'
#' @examples
#' # Example usage:
#' updated_phospho_df <- update_phospho(phospho_data, site_col = 3, gn_idx = 2, seq_len_i = 7, peptide_col_name = "Peptide")
#' head(updated_phospho_df)
#'
#' @importFrom dplyr relocate select
#' @importFrom stringr str_remove
#' @importFrom tidyr separate_rows
#' @importFrom AnnotationDbi select
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db


phosphoproteomics_update <- function(df_pho,
                                     sw_len = 7,
                                     pep_col_name = NULL,
                                     impute_method = NULL,
                                     zscore = TRUE,
                                     zmethod = "column",
                                     metric = "median"){

  message("Phosphoproteomics data: updating phosphorylation site informations")
  df_pho_update <<- update_phospho(df = df_pho,site_col = 2,gn_idx = 1,seq_len_i = 7,peptide_col_name = pep_col_name)
  message("Done!")

  message("Phosphoproteomics data: removing duplicates")
  df_pho_clean <<- remove_duplicates_phosphoproteomics(df_pho_update)

  write_xlsx(df_pho_clean, "Phosphoproteomics_clean.xlsx")
  message("Done!")

  message("Phosphoproteomics data: missing values imputation")
  df_pho_imputed <<- impute_proteomics(df_pho_clean,start_column = 6,imputation_method = impute_method)
  write_xlsx(df_pho_imputed, "Phosphoproteomics_imputed.xlsx")
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
    df_pho_zscore <<- mutate_all(as.data.frame(df_pho_zscore), as.numeric)

    df_pho_zscore <<- cbind(metadata_columns, df_pho_zscore)
    return(df_pho_zscore)
    write_xlsx(df_pho_zscore, "Phosphoproteomics_zscore.xlsx")

    message("Done!")
  } else {
    return(df_pho_imputed)
  }

}
