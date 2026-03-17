#' Filter for coding genes
#'
#' This function filters a data frame to retain only rows corresponding to
#' protein-coding genes using the Ensembl BioMart database.
#'
#' @param df tibble/data.frame, a transcriptomics dataframe containing a column
#'   named "gene_name" with gene symbols.
#'
#' @return A filtered dataframe containing only rows with protein-coding genes.
#'
#' @examples
#' filtered_df <- retrieve_coding(my_df)
#' head(filtered_df)
#'
retrieve_coding <- function(df) {
  
  # Check that gene_name column exists
  if (!"gene_name" %in% colnames(df)) {
    stop("Column 'gene_name' not found in input dataframe.")
  }
  
  # Check and install BiocManager if needed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Check and install biomaRt if needed
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt")
  }
  
  # Check and install dplyr if needed
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  
  message("Retrieving protein-coding genes from Ensembl BioMart...")
  
  # Connect to Ensembl BioMart
  mart <- biomaRt::useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl"
  )
  
  # Retrieve protein-coding genes
  map_gn <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "transcript_biotype",
    values = "protein_coding",
    mart = mart
  )
  
  # Keep only valid (non-empty) HGNC symbols
  coding_genes <- unique(map_gn$hgnc_symbol)
  coding_genes <- coding_genes[!is.na(coding_genes) & coding_genes != ""]
  
  # Filter dataframe
  df_coding <- df %>%
    dplyr::filter(gene_name %in% coding_genes)
  
  message("Protein-coding genes obtained!")
  
  return(df_coding)
}