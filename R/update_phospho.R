#' update_phospho
#'
#' This function processes phosphoproteomics data by performing several steps of data cleaning,
#' peptide modification, and enrichment of sequence information. It is designed to handle datasets
#' with phosphopeptides, modifying peptide sequences, adding multiplicity information, and computing
#' sequence windows around the phosphopeptides.
#'
#' @param df dataframe (phosphoproteomics data).
#' @param site_col integer, the column index for phosphorylation site information.
#' @param gn_idx integer, the column index for gene names in the data frame.
#' @param seq_len_i integer, the length of the sequence window to be considered around the phosphopeptide.
#' @param uniprot_idx optional integer indicating the column index for UNIPROT IDs to be added at the uniprot column retrieved with AnnotationDbi (default is NULL).
#' @param peptide_col_name (optional) a string representing the name of the peptide column in the data frame.
#'                         If not provided, the function assumes no peptide column exists.
#'
#' @return A data frame with the processed phosphoproteomics data, including modified peptide sequences,
#' updated amino acid positions, and sequence windows. If no peptide column is provided, the function
#' directly computes amino acid positions.
#' @export
#' @examples
#' # Example usage:
#' updated_phospho_df <- update_phospho(phospho_data, site_col = 3, gn_idx = 2, seq_len_i = 7, peptide_col_name = "Peptide")
#' head(updated_phospho_df)
#'

update_phospho <- function(df, site_col, gn_idx, seq_len_i=7, uniprot_idx = NULL, peptide_col_name = NULL) {

  # Step 1: Clear invalid sites
  df1 <- remove_invalid_site(df, site_col)
    # Step 2: Retrieve primary gene names
  df2 <- retrieve_primary_gene_name(df1, gn_idx, uniprot_idx)

  # Identify the sequence column
  sequence_col <- which(colnames(df2) == "Sequence")

  # Step 3: Check if peptide_col_name is provided and exists in df
  if (!is.null(peptide_col_name) && peptide_col_name %in% colnames(df)) {
    peptide_col <- which(colnames(df) == peptide_col_name)
  } else {
    peptide_col <- NULL
  }

  if (!is.null(peptide_col)) {

    # Step 4: Clean multiple peptides
    df3 <- clean_multiple_peptide(df2, peptide_col)

    # Step 5: Add multiplicity information
    df4 <- add_multiplicity(df3, peptide_col)
    mult_col <- which(colnames(df4) == "Multiplicity")

    # Step 6: Remove duplicate peptides based on multiplicity
    df5 <- remove_multiplicity(df4, peptide_col, mult_col, gn_idx)

    # Re-identify site column after modifications
    site_col <- which(colnames(df5) == "Site")

    # Step 7: Modify peptides to capitalize phosphorylated aminoacids
    df6 <- modify_peptide(df5, mult_col, peptide_col, site_col)

    # Step 8: Choose phosphopeptides
    df7 <- split_peptide(df6, mult_col, peptide_col)

    # Step 9: Add sequence window around the peptide
    df8 <- add_sequence_window(df7, peptide_col, seq_len_i, sequence_col)
    sw_col <- which(colnames(df8) == "sequence_window")

    # Step 10: Update amino acid positions based on the sequence window
    df9 <- update_position_aa(df8, sw_col, sequence_col, peptide_col)

    # Reorder and clean the dataframe columns
    df9 <- df9 %>% relocate(sequence_window, UNIPROT) %>%
      dplyr::select(-Sequence, -Site, -Peptide)
    df9 <- df9[, -ncol(df9)]


  } else {

    # If no peptide column is provided, directly create amino acid position
    df9 <- create_aa_position(df2, site_col)

    # Reorder and clean the dataframe columns
    df9 <- df9 %>%
      dplyr::relocate(UNIPROT, aminoacid, position) %>%
      dplyr::select(-Sequence, -Site)
  }

  return(df9)
}


