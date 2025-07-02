#' split_peptide
#'
#' @param phospho_df tibble, the phosphoproteomics dataframe
#' @param mult_col integer, the number of the Multiplicity column
#' @param peptide_col name of the column containing the peptide sequence
#'
#' @return phosphoproteomics dataframe with one Petide for each possible phosphosite
#'
#' @examples
#'
#' sample_df <- data.frame(
#'   Gene_name = c("ACIN1", "ZZEF1", "AAK1"),
#'   Multiplicity = c(3, 2, 1),
#'   Peptide = c("sRsRsTPVRR", "LLsFRsMEEAR", "SATTTPST*PR"),
#'   Value = c(0.1, 0.2, 0.3)
#' )
#'
#' expanded_df <- choose_phosphopeptide(sample_df, mult_col = 2, peptide_col = "Peptide")

split_peptide <- function(phospho_df, mult_col, peptide_col) {

  # Convert column indices to names if they are not already
  mult_col_name <- colnames(phospho_df)[mult_col]
  
  # Initialize an empty dataframe to store the results
  result_df <- data.frame()

  # Loop through each row of the dataframe
  for (i in 1:nrow(phospho_df)) {
    row <- phospho_df[i, ]
    multiplicity <- row[[mult_col_name]]
    peptide <- as.character(row[[peptide_col]])

    # If multiplicity is greater than 1, split the peptide and create new rows where each row has Multiplicity = 1
    if (multiplicity > 1) {

      # Find the positions of the lowercase letters
      lowercase_positions <- which(strsplit(peptide, NULL)[[1]] %in% letters)

      # Loop through each lowercase letter and create a new row with that letter capitalized
      for (pos in lowercase_positions) {
        new_peptide <- peptide
        substr(new_peptide, pos, pos) <- toupper(substr(new_peptide, pos, pos))
        new_row <- row
        new_row[[mult_col_name]] <- 1
        new_row[[peptide_col]] <- new_peptide
        result_df <- rbind(result_df, new_row)
      }
    } else {

      # If multiplicity is 1, just add the row to the result dataframe
      result_df <- rbind(result_df, row)
    }
  }

  return(result_df)
}
