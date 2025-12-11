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
#' @export
remove_multiplicity <- function(phospho_df, peptide_col, mult_col, gn_idx) {

  gene_vec <- toupper(trimws(as.character(phospho_df[[gn_idx]])))
  peptide_vec <- toupper(trimws(as.character(phospho_df[[peptide_col]])))
  multiplicity_vec <- as.numeric(phospho_df[[mult_col]])
  
  # temporary key
  phospho_df$Gene <- gene_vec
  phospho_df$Peptide <- peptide_vec
  phospho_df$Multiplicity <- multiplicity_vec
  
 
  filtered_rows <- c()
  
  for (gene in unique(phospho_df$Gene)) {
    subset <- phospho_df[phospho_df$Gene == gene, ]
    
    # find variants
    while (nrow(subset) > 0) {
      p <- subset$Peptide[1]
      matches <- grep(paste0("^", p), subset$Peptide)
      
      # retain variant with minor multiplicity
      group <- subset[matches, ]
      min_m <- min(group$Multiplicity, na.rm = TRUE)
      keep_row <- group[which(group$Multiplicity == min_m)[1], , drop = FALSE]
      
      filtered_rows <- rbind(filtered_rows, keep_row)
      
      subset <- subset[-matches, ]
    }
  }
  
  # remove temporary columns
  filtered_rows$Gene <- NULL
  filtered_rows$Peptide <- NULL
  filtered_rows$Multiplicity <- NULL
  
  return(filtered_rows)
}
