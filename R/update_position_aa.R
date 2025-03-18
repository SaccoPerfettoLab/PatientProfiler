#' update_position_aa
#'
#' @param phospho_df tibble, the phosphoproteomics dataframe
#' @param sw_col integer, the number of the sequence window column
#' @param sequence_col integer, the number of the full Sequence column
#' @param peptide_col integer, the number of the Peptide column
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


update_position_aa <- function(phospho_df, sw_col, sequence_col, peptide_col) {
  # Inizializza le colonne 'position' e 'aminoacid'
  phospho_df$position <- NA
  phospho_df$aminoacid <- NA

  for (r in 1:nrow(phospho_df)) {
    # Estrai la finestra di sequenza
    sequence_window <- as.character(phospho_df[[sw_col]][r])

    # Controlla se la finestra di sequenza è NA
    if (!is.na(sequence_window) && nchar(sequence_window) >= 8) {
      # Estrai l'ottavo carattere nella finestra di sequenza
      aa <- substr(sequence_window, 8, 8)

      # Rimuovi i caratteri '_' dalla finestra di sequenza per la ricerca
      sequence_window_clean <- gsub("_", "", sequence_window)

      # Trova la posizione di inizio della finestra di sequenza pulita nella sequenza completa
      position_start <- regexpr(sequence_window_clean, phospho_df[[sequence_col]][r])
      position_start <- as.numeric(position_start)

      # Se la finestra di sequenza pulita è trovata nella sequenza completa
      if (position_start != -1) {
        # Calcola il numero di caratteri '_' a sinistra
        padding_left <- regexpr("[^_]", sequence_window)[1] - 1

        # Calcola la posizione dell'ottavo carattere nella sequenza completa
        position <- position_start + (8 - 1) - padding_left

        # Assegna la posizione e l'aminoacido
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

  # Ricolloca le colonne 'position' e 'aminoacid' all'inizio del dataframe
  phospho_df <- phospho_df %>% dplyr::relocate(position, aminoacid)

  return(phospho_df)
}










