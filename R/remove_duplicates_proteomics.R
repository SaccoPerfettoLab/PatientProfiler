#' Remove Duplicate Entries in Proteomics Data
#'
#' This function processes a proteomics data frame to remove duplicate entries
#' based on the `gene_name` column. It ensures that for each gene, a unique entry
#' is retained, and all associated `UNIPROT` identifiers are combined into a single string,
#' separated by semicolons.
#'
#' @param df dataframe (proteomics data). It must include at least the 'gene_name' and 'UNIPROT' columns.
#'
#' @return a cleaned dataframe where each `gene_name` has a unique entry, and
#' associated `UNIPROT` identifiers are aggregated.
#'
#' @examples
#' cleaned_proteomics <- remove_duplicates_proteomics(proteomics_data)
#' head(cleaned_proteomics)
#'


remove_duplicates_proteomics <- function(df) {

  duplicated_genes <- df$gene_name[duplicated(df$gene_name)]

  proteomics_combined <- df %>%
    dplyr::group_by(gene_name) %>%  # Gruppo per 'gene_name'
    dplyr::reframe(UNIPROT = paste(unique(UNIPROT), collapse = ";"))

  proteomics_unique <- df %>%
    dplyr::distinct(gene_name, .keep_all = TRUE)

  proteomics_unique <- proteomics_unique %>%
    dplyr::select(-UNIPROT)

  proteomics_final <- proteomics_unique %>%
    dplyr::left_join(proteomics_combined, by = "gene_name")

  proteomics_final <- proteomics_final %>%
    dplyr::relocate(UNIPROT)

  return(proteomics_final)
}
