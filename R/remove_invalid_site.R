
#' remove_invalid_site
#'
#' @param phospho_df tibble, the phosphoproteomics dataframe
#' @param site_col integer, the column number of the Site
#'
#' @return phosphoproteomics dataframe with no uncertain site
#'
#' @examples
#'
#' sample_df <- data.frame(
#'   Gene_name = c("UNKL", "ATF7", "SAMD1"),
#'   Site = c("S100T345", "T56", "Y389"),
#'   Value = c(0.1, 0.2, 0.3)
#' )
#'
#' cleaned_df <- remove_invalid_site(sample_df, site_col = 2)

remove_invalid_site <- function(phospho_df, site_col) {

  # Get the name of the Site column using the column number or name
  site_col_name <- if (is.numeric(site_col)) names(phospho_df)[site_col] else site_col

  # Count the number of letters in each string of the Site column
  count_letters <- sapply(gregexpr("[A-Za-z]", phospho_df[[site_col_name]]), function(x) sum(x > 0))

  # Identify invalid rows where there are more than one letter
  invalid_rows <- count_letters > 1

  # Filter out the invalid rows
  phospho_df_filt <- phospho_df[!invalid_rows, ]

  return(phospho_df_filt)
}
