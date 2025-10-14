#' clean_multiple_peptide
#'
#' @param phospho_df tibble, the phosphoproteomics dataframe
#' @param peptide_col name of the column containing the peptide sequence
#'
#' @return phosphoproteomics dataframe with one Peptide, the longest, for each phosphosite
#'
#' @examples
#' sample_phospho_df <- data.frame(
#'   Gene_nname = c("UNKL", "SRGAP3", "NSUN5"),
#'   Peptide = c("AAAAALSGsPPQTEK;AAAAALSGsPPQTEKPTHYR;", "AAACPSsPHK", "AAAGACtPPCT"),
#'   Value = c(0.1, 0.2, 0.3)
#' )
#'
#' cleaned_df <- clean_multiple_peptide(sample_df, peptide_col = "Peptide")
#' @export
clean_multiple_peptide <- function(phospho_df, peptide_col) {

  for (i in 1:nrow(phospho_df)) {

    peptide <- as.character(phospho_df[[peptide_col]][i])
    
    if (!grepl(";", peptide)) {
      next
    }

    sequences <- unlist(strsplit(peptide, ";"))

    longest_sequence <- sequences[which.max(nchar(sequences))]

    phospho_df[[peptide_col]][i] <- longest_sequence
  }

  return(phospho_df)
}
