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

  df[is.na(df) | df == -Inf | df == Inf] <- 0

  max_na <- ncol(df[, -1]) * threshold / 100

  df_clean <- df[rowSums(is.na(df)) <= max_na, ]

  return(df_clean)
}

