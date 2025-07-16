#' impute_proteomics
#'
#' This function imputes missing values in a numeric data frame using the `mice` package and saves the imputed data
#' to a global variable for user access.
#'
#' @param df dataframe with missing data.
#' @param start_column integer, the first numeric column to consider for imputation (default: 1).
#' @param imputation_method string, the imputation method for `mice` (default: "pmm").
#' @return The imputed data frame.
#'
#' @examples
#' result <- impute_data(df, start_column = 1, imputation_method = "pmm")
#' View(imputed_df) # Access the saved variable
#'
#' @import mice

impute_proteomics <- function(df, start_column = 1, imputation_method = NULL) {

  df[sapply(df, is.numeric)] <- lapply(df[sapply(df, is.numeric)], function(x) {x[is.infinite(x) | is.nan(x)] <- NA; x})
  
  column_names <- names(df)[start_column:ncol(df)]
  cleaned_column_names <- lapply(column_names, function(x) gsub("[^A-Za-z0-9]", "", as.character(x)))
  names(df)[start_column:ncol(df)] <- cleaned_column_names

  init = mice(df, maxit = 0)
  meth = init$method
  predM = init$predictorMatrix

  if (!is.null(imputation_method)) {
    meth[colnames(df)[start_column:ncol(df)]] <- imputation_method
  }

  set.seed(103)
  imputed = mice(df, method = meth, predictorMatrix = predM)

  imputed <- complete(imputed)

  return(imputed)
}
