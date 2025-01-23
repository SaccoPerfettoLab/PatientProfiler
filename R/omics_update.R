#' omics_update
#'
#' This function updates transcriptomics, proteomics, and phosphoproteomics datasets by performing
#' preprocessing steps such as filtering, imputation, and Z-score normalization.
#' Dataframes must contain the following required columns:
#' - df_tr: 'gene_name' and sample columns.
#' - df_pr: 'gene_name' and sample columns.
#' - df_ph: 'gene_name', 'Site' (containing phosphorylation sites in the format: T161), and 'Peptide'
#'   (containing the peptide sequence; the name and the presence of the peptide column is optional).

#' @param df_tr dataframe containing transcriptomics data.
#' @param df_pr dataframe containing proteomics data.
#' @param df_ph dataframe containing phosphoproteomics data.
#' @param threshold numeric, it specifies the threshold percentage of zeros in each row the transcriptomic dataframe. If there are >= threshold percentage of
#'                  zeros in a row, the row will be removed (default: 80).
#' @param sw_len integer, the sequence window length. It could be 7 or 15 (default: 7).
#' @param pep_col_name string, it indicates the name of the peptide sequence column (if it's present) in your phosphoproteomic dataframe (default = "Peptide")
#' @param imp_method string, the method to use for imputation (default: "pmm", but you can choose between "pmm", "norm", "norm.nob",
#'                   "regression", "ri", "logreg", "polyreg", "predictive", "polr", "sample", "cart", "knn", "rf").
#' @param zscore logical, whether to perform Z-score normalization (default: "TRUE").
#' @param zmethod string, specifies whether Z-score normalization is performed by "row" or "column" (default: "column").
#' @param metric string, the centering metric for Z-score normalization. Options are "median" (default) or "mean".

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
#'   imp_method = "pmm",
#'   zscore = TRUE,
#'   zmethod = "column",
#'   metric = "median"
#' )
#'
#' @importFrom mice mice
#' @importFrom dplyr filter
#' @export

omics_update <- function(df_tr = NULL,
                 df_pr = NULL,
                 df_ph = NULL,
                 threshold = 80,
                 sw_len = 7,
                 pep_col_name = NULL,
                 imp_method = "pmm",
                 zscore = "TRUE",
                 zmethod = "column",
                 metric = "median") {


  if (!is.null(sw_len) && !sw_len %in% c(7, 15)) {
    stop("sw_len parameter has to be 7, 15 or NULL (if you do not have the phosphoproteomics dataframe")
  }

  if (!is.null(df_tr)) {
    message("Transcriptomics update started..")

    transcriptomics_updated <<- transcriptomics_update(df_tr,threshold,zscore,zmethod,metric)

    write_xlsx(transcriptomics_updated, "Transcriptomics_updated.xlsx")
    message("Transcriptomics update complete!")
    }

  if(!is.null(df_pr)){
    message("Proteomics update started..")
    proteomics_updated <<- proteomics_update(df_pr,imp_method,zscore,zmethod,metric)
    write_xlsx(proteomics_updated, "Proteomics_updated.xlsx")

    message("Proteomics update complete!")
    }

  if(!is.null(df_ph)){
    message("Phosphoproteomics update started..")
    phosphoproteomics_updated <<- phosphoproteomics_update(df_pho = df_ph,
                                                           sw_len = 7,
                                                           pep_col_name = pep_col_name,
                                                           imp_method,
                                                           zscore,
                                                           zmethod,
                                                           metric)

    write_xlsx(phosphoproteomics_updated, "Phosphoproteomics_updated.xlsx")

    message("Phosphoproteomics update complete!")
    }

}


