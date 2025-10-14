#' add_sequence_window
#'
#' @param phospho_df tibble, the phosphoproteomics dataframe
#' @param peptide_col name of the column containing the peptide sequence
#' @param seq_len_i integer, the length of the mers, 7 for 15mers
#' @param sequence_col name of the column containing the protein sequence
#'
#' @return phosphoproteomics dataframe with the phosphosite around the central position of the Peptide
#'
#' @examples
#'sample_df <- data.frame(
#'   Gene_name = c("ACIN1", "ZZEF1", "AAK1"),
#'   Peptide = c("SRsRsTPVR", "LLsFRsMEEAR", "SATTTPST*PR"),
#'   Sequence = c("MDLEGDRNGGAKSRSRSTPVRKKNFFKLN", "MLENLQKHSTLLSFRSMEEARIGEEMSQNSFIKQY", "YKPQKGLSATTTPSTPRNLGMSEQL"),,
#'   Value = c(0.1, 0.2, 0.3)
#' )
#'
#' result_df <- add_sequence_window(sample_df, peptide_col = "Peptide", seq_len_i = 7, sequence_col = "Sequence")
#' @export
add_sequence_window <- function(phospho_df, peptide_col, seq_len_i = 7, sequence_col = "Sequence" ) {
  
  phospho_df$sequence_window <- NA
  
  for (i in 1:nrow(phospho_df)) {
    peptide <- as.character(phospho_df[[peptide_col]][i])
    sequence <- as.character(phospho_df[[sequence_col]][i])
    
    pos_lower <- regexpr("[a-z]", peptide)[1]
    
    pos_star <- regexpr("\\*", peptide)[1]
    
    if (pos_lower != -1) {
      phosphorylated_residue <- substr(peptide, pos_lower, pos_lower)
      central_pos_offset <- pos_lower
    } else if (pos_star != -1) {
      phosphorylated_residue <- substr(peptide, pos_star - 1, pos_star - 1)
      central_pos_offset <- pos_star - 1
    } else {
      phospho_df$sequence_window[i] <- NA
      next
    }
    
    phosphorylated_residue_upper <- toupper(phosphorylated_residue)
    
    start_pos <- regexpr(toupper(peptide), sequence)[1]
  
    if (start_pos == -1) {
      phospho_df$sequence_window[i] <- NA
      next
    }
    
    central_pos <- start_pos + central_pos_offset - 1
    
    start <- max(1, central_pos - seq_len_i)
    end <- min(nchar(sequence), central_pos + seq_len_i)
    window <- substr(sequence, start, end)
    
    if (start == 1) {
      window <- paste0(strrep("_", seq_len_i - (central_pos - 1)), window)
    }
    if (end == nchar(sequence)) {
      window <- paste0(window, strrep("_", seq_len_i - (nchar(sequence) - central_pos)))
    }
    
    required_length <- 2 * seq_len_i + 1
    current_length <- nchar(window)
    
    if (current_length < required_length) {
      total_padding <- required_length - current_length
      padding_left <- total_padding %/% 2
      padding_right <- total_padding - padding_left
      window <- paste0(strrep("_", padding_left), window, strrep("_", padding_right))
    }
    
    phospho_df$sequence_window[i] <- window
  }
  
  return(phospho_df)
}












