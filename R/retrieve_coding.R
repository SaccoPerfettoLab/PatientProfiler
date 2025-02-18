#' Filter for coding genes
#'
#' This function filters a data frame to retain only rows corresponding to coding genes based on the `org.Hs.eg.db` database.
#'
#' @param df tibble, a transcriptomics dataframe containing gene names.
#' @param gn_col string, the name of the column in `df` containing gene symbols to match against coding genes (default: "gene_name").
#'
#' @return a filtered dataframe containing only rows with coding genes.
#'
#' @examples
#' # Example usage:
#' filtered_df <- retrieve_coding(my_df, gn_col = "GeneSymbol")
#' head(filtered_df)
#'

retrieve_coding <- function(df) {
  # Check and install BiocManager if not already installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager...")
    install.packages("BiocManager")
  }

  # Check and install AnnotationDbi if not already installed
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    message("Installing AnnotationDbi from Bioconductor...")
    BiocManager::install("AnnotationDbi")
  }

  # Check and install org.Hs.eg.db if not already installed
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message("Installing org.Hs.eg.db from Bioconductor...")
    BiocManager::install("org.Hs.eg.db")
  }

  # Retrieve coding genes
  message("Retrieving coding genes from org.Hs.eg.db...")
  coding_genes <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENSEMBL"),
    columns = c("ENSEMBL", "SYMBOL"),
    keytype = "ENSEMBL"
  )

  # Filter the input dataframe
  df_coding <- dplyr::filter(df, gene_name %in% coding_genes$SYMBOL)
  message("Coding genes obtained!")

  return(df_coding)
}
