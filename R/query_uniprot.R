
#' Query UNIPROT Proteins
#'
#' This function retrieves protein-related information from the UNIPROT database.
#' It performs batch queries to efficiently retrieve data.
#'
#' @param id_input A vector of UNIPROT IDs.
#' @param batch_size An integer specifying the number of IDs per query batch (default is 400).
#' @return A dataframe with UNIPROT information including accession, gene names, protein names, and sequence details.
#' @export
query_uniprot_proteins <- function(id_input, batch_size = 400){
  
  # Filter only valid UNIPROT IDs
  id_input <- unique(id_input)
  id_input <- id_input[grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', id_input)]
  id_input <- id_input[!grepl('CHEBI:|SIGNOR-|NP_|CID:|SID:|PRO_|DB0|_9606', id_input)]
  
  # Define column names for the final dataframe
  header_df <- c("Entry", "Reviewed", "Entry Name", "Protein names", "Gene Names (primary)", "Organism", "Length", "Sequence", "Date_modified")
  df_final <- data.frame(matrix(ncol = length(header_df), nrow = 0))
  colnames(df_final) <- header_df
  
  for (i in seq(1, length(id_input), by = batch_size)) {
    print(paste("Processing batch starting at index:", i))
    batch_ids <- unique(id_input[i:min(i + batch_size - 1, length(id_input))])
    
    query_string <- paste0('accession%3A', paste0(batch_ids, collapse = '%20OR%20accession%3A'))
    url <- paste0('https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_primary%2Corganism_name%2Clength%2Csequence%2Cdate_sequence_modified&format=tsv&query=', query_string)
    
    result <- httr::GET(url)
    file_content <- as.character(httr::content(result, "text"))
    
    df_batch <- readr::read_delim(file_content, delim = '\t', skip_empty_rows = TRUE, show_col_types = FALSE)
    df_final <- dplyr::bind_rows(df_final, df_batch)
  }
  
  return(df_final %>% dplyr::distinct())
}

#' Query UNIPROT
#'
#' Retrieves protein data from UNIPROT, handling isoforms separately.
#'
#' @param id_input A vector of UNIPROT IDs.
#' @param batch_size Number of IDs per batch (default = 400).
#' @return A dataframe with retrieved protein data.
#' @export
query_uniprot <- function(id_input, batch_size = 400) {
  
  id_input <- unique(id_input)
  id_input <- id_input[grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', id_input)]
  id_input <- id_input[!grepl('CHEBI:|SIGNOR-|NP_|CID:|SID:|PRO_|DB0|_9606|VAR|REV', id_input)]
  
  # Separate isoforms
  id_input_iso <- unique(id_input[grepl('-', id_input)])
  
  if (length(id_input_iso) != 0) {
    if (length(id_input_iso) > batch_size) {
      result_iso <- query_uniprot_proteins(id_input = id_input_iso, batch_size = batch_size)
    } else {
      result_iso <- query_uniprot_proteins(id_input = id_input_iso, batch_size = length(id_input_iso))
    }
  }
  
  # Query for non-isoform proteins
  result_proteins <- query_uniprot_proteins(id_input, batch_size)
  
  # Merge results
  if (exists('result_iso')) {
    result_total <- dplyr::bind_rows(result_iso, result_proteins)
  } else {
    result_total <- result_proteins
  }
  
  return(result_total)
}
