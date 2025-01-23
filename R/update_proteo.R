#' update_proteo
#'
#' @param df tibble, the dataframe containing proteomics data
#' @param gn_idx integer, the number of the gene name column
#' @param sequence boolean, TRUE or FALSE
#'
#' @return The proteomics dataframe with updated gene names
#'
#' @examples
#'
#' sample_df <- data.frame(
#'   gene_name = c("AAAS", "AAGAB", "AAK1"),
#'   Value = c(0.1, 0.2, 0.3)
#' )
#'
#' updated_df <- update_proteo(sample_df, gn_idx = 2, sequence = FALSE)


update_proteo <- function(df, gn_idx, sequence = FALSE) {

  df_uniprot <- retrieve_primary_gene_name(df, gn_idx, sequence = FALSE)

  df_uniprot <- dplyr::relocate(df_uniprot,UNIPROT)

  return(df_uniprot)
}
