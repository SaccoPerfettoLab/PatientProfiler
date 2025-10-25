#' remove_nas
#'
#' This function removes rows from a data frame where the percentage of missing values (NAs) 
#' exceeds a specified threshold. It also treats NaN, Inf, and -Inf values as missing (NA).
#'
#' @param df tibble, a dataframe containing omics data. 
#'   The first column is assumed to be non-numeric (e.g., identifiers).
#' @param threshold numeric, the maximum percentage of missing values allowed per row (e.g., 80 for 80%).
#'   Rows with more missing values are removed.
#' @return tibble, a cleaned dataframe with rows containing excessive missing values removed.
#'
#' @examples
#' df_transcriptomics <- data.frame(
#'   ID = c("Gene1", "Gene2", "Gene3"),
#'   Sample1 = c(1, NA, 3),
#'   Sample2 = c(NA, NaN, 4),
#'   Sample3 = c(2, Inf, -Inf)
#' )
#' cleaned_data <- remove_nas(df_transcriptomics, threshold = 50)
#' print(cleaned_data)
#'
remove_nas <- function(df, threshold) {
  
  # Identify numeric columns
  num_cols <- sapply(df, is.numeric)
  
  # Convert NaN, Inf, and -Inf to NA
  df[, num_cols] <- lapply(df[, num_cols], function(x) {
    x[is.nan(x) | is.infinite(x)] <- NA
    return(x)
  })
  
  # Count number of NAs per row (only in numeric columns)
  na_counts <- rowSums(is.na(df[, num_cols]))
  
  # Calculate the maximum number of allowed NAs per row
  max_na <- ncol(df[, num_cols]) * threshold / 100
  
  # Keep only rows with NA count below or equal to the threshold
  df_clean <- df[na_counts <= max_na, ]
  
  return(df_clean)
}
