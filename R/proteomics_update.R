#'proteomics_update
#'
#'
#' @param df_pro dataframe, containing a column with protein names.
#' @param imp_method string, the method to use for imputation (default: "pmm", but you can choose between "pmm", "norm", "norm.nob",
#'                   "regression", "ri", "logreg", "polyreg", "predictive", "polr", "sample", "cart", "knn", "rf").
#' @param sequence A boolean (TRUE or FALSE) indicating whether the dataset includes a column named "sequence".
#'                 If TRUE, the function assumes the "sequence" column is present and may process it accordingly.
#'                 If FALSE, the function proceeds without expecting the "sequence" column.
#' @param zscore boolean, indicates whether to calculate z-scores for the data.
#'               If TRUE, z-scores are calculated; if FALSE, a data frame will be returned without z-score calculation.
#' @param zmethod string, it specifies the method for calculating z-scores. It can either be "row" to compute z-scores for each row,
#'                or "column" to compute z-scores for each column of the data.
#' @param metric  string, it indicates the method for centering the data when z-scores are computed.
#'                It could be "mean" or "median".
#' @param output_dir a string indicating the updated output folder.
#' 
#' @return tibble, it contains the filtered and processed proteomics data.
#'
#'
#'
#'@examples
#'prot_updated <- proteomics_update(df, impute_method = "pmm",zscore = TRUE,zmethod = "column",metric = "median")

proteomics_update <- function(df_pro,
                              impute_method = NULL,
                              zscore = TRUE,
                              zmethod = "column",
                              metric = "median",
                              output_dir) {
  
  message("Proteomics data: updating UNIPROT")
  # Seleziona le colonne numeriche da aggregare
  df_pro_agg <- aggregate(df_pro[, 2:ncol(df_pro)],
                          by = list(gene_name = df_pro$gene_name),
                          FUN = mean, na.rm = TRUE)
  
  # Restituisci il dataframe aggregato
  df_pro_agg <- as.data.frame(df_pro_agg)
  
  # Update proteomics data (aggiornamenti ulteriori)
  
  df_pro_update <<- update_proteo(df_pro_agg, 1, sequence = FALSE)
  message("Done!")
  
  # Rimuovi i duplicati
  message("Proteomics data: removing duplicates")
  df_pro_clean <<- remove_duplicates_proteomics(df_pro_update)
  readr::write_tsv(df_pro_clean, paste0(output_dir,"/","Proteomics_clean.tsv"))
  message("Done!")
  
  # Imputazione dei valori mancanti
  message("Proteomics data: missing values imputation")
  df_pro_imputed <<- impute_proteomics(df_pro_clean, 3, impute_method)
  readr::write_tsv(df_pro_imputed, paste0(output_dir,"/","Proteomics_imputed.tsv"))
  message("Done!")
  
  if (zscore) {
    message("Proteomics data: computing Z-scores")
    df_pro_numeric <- df_pro_imputed %>%
      dplyr::select(-UNIPROT, -gene_name)
    
    metadata_columns <- df_pro_imputed %>%
      dplyr::select(UNIPROT,gene_name)
    
    df_pro_matrix <- as.matrix(df_pro_numeric)
    df_pro_zscore <<- compute_zscore(df_pro_matrix, zmethod, metric)
    df_pro_zscore <<- dplyr::mutate_all(as.data.frame(df_pro_zscore), as.numeric)
    
    df_pro_zscore <<- cbind(metadata_columns, df_pro_zscore)
    readr::write_tsv(df_pro_zscore, paste0(output_dir,"/","Proteomics_zscore.tsv"))
    
    return(df_pro_zscore)
    message("Done!")
    
    return(df_pro_zscore)
  } else {
    return(df_pro_imputed)
  }
}
