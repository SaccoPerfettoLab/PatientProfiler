#' update_position_aa
#'
#' @param phospho_df tibble, the phosphoproteomics dataframe
#' @param sw_col integer, the number of the sequence window column
#' @param sequence_col name of the column containing the protein sequence
#' @param peptide_col name of the column containing the peptide sequence
#'
#' @return phosphoproteomics dataframe with the additional columns 'position' and 'aminoacid' of each phosphosite
#'
#' @examples
#'
#' sample_df <- data.frame(
#'   Gene_name = c("ACIN1", "ZZEF1", "AAK1"),
#'   sequence_window = c("HSRSRSRSTPVRD", "LRLLSFRSMEEARLV", "ATTTPSGSPRTSQQN"),
#'   Sequence = c("MDLEGDRNGGAKSRSRSTPVRKKNFFKLN", "MLENLQKHSTLLSFRSMEEARIGEEMSQNSFIKQY", "YKPQKGLSATTTPSTPRNLGMSEQL"),,
#'   Peptide = c("SRsRsTPVR", "LLsFRsMEEAR", "SATTTPST*PR"),
#'   Value = c(0.1, 0.2, 0.3)
#' )
#'
#' updated_df <- update_position_aa(sample_df, sw_col = 2, sequence_col = 3, peptide_col = 4)
#' @export


update_position_aa <- function(phospho_df, sw_col, sequence_col, peptide_col) {
  phospho_df$position <- NA
  phospho_df$aminoacid <- NA

  for (r in 1:nrow(phospho_df)) {
    sequence_window <- as.character(phospho_df[[sw_col]][r])

    if (!is.na(sequence_window) && nchar(sequence_window) >= 8) {
      aa <- substr(sequence_window, 8, 8)

      sequence_window_clean <- gsub("_", "", sequence_window)

      position_start <- regexpr(sequence_window_clean, phospho_df[[sequence_col]][r])
      position_start <- as.numeric(position_start)

      if (position_start != -1) {
        padding_left <- regexpr("[^_]", sequence_window)[1] - 1

        position <- position_start + (8 - 1) - padding_left

        phospho_df$position[r] <- position
        phospho_df$aminoacid[r] <- aa
      } else {
        phospho_df$position[r] <- NA
        phospho_df$aminoacid[r] <- NA
      }
    } else {
      phospho_df$position[r] <- NA
      phospho_df$aminoacid[r] <- NA
    }
  }

  phospho_df <- phospho_df %>% dplyr::relocate(position, aminoacid)

  return(phospho_df)
}










