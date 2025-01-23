#' create_aa_position
#'
#' @param phospho_df tibble, the phosphoproteomics dataframe
#' @param site_col integer, the number of the Site column
#'
#' @return phosphoproteomics dataframe with the additional columns 'position' and 'aminoacid' of each phosphosite
#'
#' @examples
#'
#' sample_df <- data.frame(
#'   Gene_name = c("AAAS", "AAGAB", "AAK1"),
#'   Site = c("S462", "S202", "S21"),
#'   Value = c(0.1, 0.2, 0.3)
#' )
#'
#' updated_df <- create_aa_position(sample_df, site_col = 2)

create_aa_position <- function(df, site_col) {

  site <- df[[site_col]]

  aminoacid <- substr(site, 1, 1)

  position <- as.numeric(substring(site, 2))

  df$aminoacid <- aminoacid
  df$position <- position

  return(df)
}
