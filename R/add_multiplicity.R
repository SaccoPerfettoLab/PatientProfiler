#' add_multiplicity
#'
#' @param phospho_df tibble, the phosphoproteomics dataframe
#' @param peptide_col integer, the column number of the Peptide
#'
#' @return phosphoproteomics dataframe with an added 'Multiplicity' column that counts phosphosites in the peptide sequences
#'
#' @examples
#' sample_df <- data.frame(
#'   Gene_name = c("ACIN1", "ZZEF1", "AAK1")
#'   Peptide = c("sRsRsTPVRR", "LLsFRsMEEAR", "SATTTPSG*PR"),
#'   Value = c(0.1, 0.2, 0.3)
#' )
#'
#' modified_df <- add_multiplicity(sample_df, peptide_col = 2)

add_multiplicity <- function(df, peptide_col) {

  pattern <- ifelse(grepl("[*]", df[[peptide_col]]), "[*]", "[a-z]")

  df$Multiplicity <- stringr::str_count(df[[peptide_col]], pattern)

  return(df)
}
