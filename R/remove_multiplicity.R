#' remove_multiplicity
#'
#' @param phospho_df tibble, the phosphoproteomics dataframe
#' @param peptide_col name of the column containing the peptide sequence
#' @param mult_col integer, the number of the Multiplicity column
#' @param gn_idx integer, number of the Gene name column
#'
#' @return phosphoproteomics dataframe with duplicates removed based on multiplicity
#'
#' @examples
#'
#'#' sample_df <- data.frame(
#'   Gene_name = c("ACIN1", "ZZEF1", "AAK1")
#'   Peptide = c("sRsRsTPVRR", "LLsFRsMEEAR", "SATTTPSG*PR"),
#'   Multiplicity = c(3, 2, 1),
#'   Value = c(0.1, 0.2, 0.3)
#' )
#'
#' cleaned_df <- remove_multiplicity(sample_df, peptide_col = "Peptide", mult_col = 3, gn_idx = 1)

remove_multiplicity <- function(phospho_df, peptide_col, mult_col, gn_idx) {

  # Create a new column as a key combining gene name and peptide sequence
  phospho_df$Name_Peptides <- paste0(
    as.character(unlist(phospho_df[, gn_idx])), "_",
    as.character(phospho_df[[peptide_col]])
  )
  
  # Create a vector of unique keys from the combined gene name and peptide sequence
  keys <- as.vector(unlist(unique(paste0(
    as.character(unlist(phospho_df[, gn_idx])), "_", 
    as.character(unlist(phospho_df[, peptide_col]))))))
  
  # Iterate over each unique key
  for (key in keys) {

    # Identify rows corresponding to the current key
    pos <- which(toupper(phospho_df$Name_Peptides) == toupper(key))
    
    # If only one row corresponds to the key, continue to the next iteration
    if (length(pos) == 1) {
      next
    } else {

      # Create a subset of rows with the same key
      subset <- phospho_df[pos, ]

      # Identify rows to remove based on the minimum multiplicity value
      pos_to_remove <- pos[which(subset[[mult_col]] != min(subset[[mult_col]], na.rm = TRUE))]
      
      # Remove the identified rows from the dataframe
      if (length(pos_to_remove) >= 1) {
        phospho_df <- phospho_df[-c(pos_to_remove), ]
      } else {
        next
      }
    }
  }

  # Remove the temporary 'Name_Peptides' column
  phospho_df <- phospho_df[, !colnames(phospho_df) %in% "Name_Peptides"]

  return(phospho_df)
}
