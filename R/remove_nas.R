#' remove_nas
#'
#' This function removes rows from a data frame where the percentage of missing values (NAs) exceeds a specified threshold.
#'
#' @param t_df tibble, a dataframe containing transcriptomics data. The first column is assumed to be non-numeric (e.g., identifiers).
#' @param threshold integer, the maximum percentage of missing values allowed per row (e.g., 80 for 80%). Rows with more missing values are removed.
#' @return tibble, a cleaned dataframe with rows containing excessive missing values removed.
#'
#' @examples
#' # Example usage:
#' df_transcriptomics <- data.frame(
#'   ID = c("Gene1", "Gene2", "Gene3"),
#'   Sample1 = c(1, NA, 3),
#'   Sample2 = c(NA, NA, 4),
#'   Sample3 = c(2, 1, NA)
#' )
#'
#' cleaned_data <- remove_nas(df_transcriptomics, threshold = 50)
#' print(cleaned_data)


remove_nas <- function(df, threshold) {
  
  tol = 1e-6
  # Seleziona solo le colonne numeriche
  num_cols <- sapply(df, is.numeric)
  
  # Sostituisci NA, -Inf, Inf e valori prossimi a zero con 0 nelle colonne numeriche
  df[, num_cols] <- lapply(df[, num_cols], function(x) {
    x[is.na(x) | x == -Inf | x == Inf | abs(x) < tol] <- 0
    return(x)
  })
  
  # Calcola il numero massimo di valori "0" per colonna, in base alla soglia
  max_zero <- ncol(df[, num_cols]) * threshold / 100
  
  # Rimuovi le righe che hanno più del limite di "zero" nelle colonne numeriche
  df_clean <- df[rowSums(df[, num_cols] == 0) <= max_zero, ]
  
  return(df_clean)
}




