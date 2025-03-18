#' transcriptomics_update
#'
#' This function processes and filters a transcriptomics dataframe based on user-defined parameters.
#' It can filter rows based on the percentage of missing values (NA), and optionally compute z-scores
#' for the data either by rows or columns.
#'
#' @param df_tra tibble, a data frame containing the transcriptomics data.
#' @param threshold numeric, it specifies the threshold percentage of zeros in each row the transcriptomic dataframe. If there are >= threshold percentage of
#'                  zeros in a row, the row will be removed.
#' @param zscore boolean, indicates whether to calculate z-scores for the data.
#'               If TRUE, z-scores are calculated; if FALSE, a data frame will be returned without z-score calculation.
#' @param zmethod string, it specifies the method for calculating z-scores. It can either be "row" to compute z-scores for each row,
#'                or "column" to compute z-scores for each column of the data.
#' @param metric  string, it indicates the method for centering the data when z-scores are computed.
#'                It could be "mean" or "median".
#' @param output_dir a string indicating the updated output folder.
#'                
#' @return tibble, it contains the filtered and processed transcriptomics data.
#'
#' @examples
#'
#' df_processed <- transcriptomics_update(df_tra, zscore = TRUE, zmethod = "median", metric = "mean")
#'
#' df_cleaned <- transcriptomics_update(df_tra, zscore = FALSE, zmethod = "median", metric = "mean")


transcriptomics_update <- function(df_tra,
                                   threshold = 80,
                                   zscore = TRUE,
                                   zmethod = "column",
                                   metric = "median",
                                   output_dir) {
  
  # df_tra_cod <- retrieve_coding(df_tra)
  
  colnames(df_tra) <- gsub("[.-]", "", colnames(df_tra))
  
  df_tra_clean <<- remove_nas(df_tra, threshold)
  
  message("Transcriptomics data: removing duplicates")
  
  df_tra_agg <<- df_tra_clean %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarize(across(everything(), mean, na.rm = TRUE))
  
  readr::write_tsv(df_tra_agg, paste0(output_dir,"/","Transcriptomics_clean.tsv"))
  
  if (zscore) {
    
    message("Transcriptomics data: computing zscore")
    
    df_tra_numeric <- df_tra_agg %>%
      dplyr::select(-gene_name)
    
    metadata_columns <- df_tra_agg %>%
      dplyr::select(gene_name)
    
    df_tra_matrix <- as.matrix(df_tra_numeric)
    
    df_tra_zscore <<- compute_zscore(df_tra_matrix, zmethod, metric)
    
    df_tra_zscore <<- dplyr::mutate_all(as.data.frame(df_tra_zscore), as.numeric)
    
    df_tra_zscore <<- cbind(metadata_columns, df_tra_zscore)
    
    readr::write_tsv(df_tra_zscore, paste0(output_dir,"/","Transcriptomics_zscore.tsv"))
    return(df_tra_zscore)
  } else {
    return(df_tra_agg)
  }
}

