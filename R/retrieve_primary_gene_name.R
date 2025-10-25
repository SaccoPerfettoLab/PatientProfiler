#' Retrieve Primary Gene Name from UNIPROT Database
#'
#' This function updates a given omics data frame by retrieving the primary gene name
#' and associated UNIPROT ID for each gene symbol. It queries the UNIPROT database
#' for human or mouse gene symbols and returns the most relevant data based on sequence
#' modifications and review status.
#'
#' @param omic_df dataframe containing omics data. It must have a column with gene symbols.
#' @param gn_idx integer, the column index for gene names in the omics data frame.
#' @param uniprot_idx optional integer indicating the column index for UNIPROT IDs to be added at the uniprot column retrieved with AnnotationDbi (default is NULL).
#' @param organism string, it indicates the organism. Options are 'human' (default) or 'mouse'.
#' @param sequence logical, indicating whether to retain the protein sequence information in the result (default is TRUE).
#'
#' @return A data frame with updated gene names and corresponding UNIPROT IDs. If `sequence = TRUE`, the sequence information is also included.
#'
#' @examples
#' # Example usage:
#' updated_omic_df <- retrieve_primary_gene_name(omic_df, gn_idx = 2, organism = 'human')
#' head(updated_omic_df)
#' @export

retrieve_primary_gene_name <- function(omic_df = br_pr, gn_idx=1, uniprot_idx = NULL, organism = 'human', sequence = TRUE){

  # Rename colnames
  colnames(omic_df)[gn_idx] <- 'gene_name'
  omic_df %>% tidyr::separate_rows(gene_name, sep = ';') -> omic_df

  gene_names <- unlist(unique(omic_df$gene_name))

  if(organism == 'human'){
    uniprot_ids_df <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                            keys = stringr::str_remove(gene_names, '^MT-'), #remove mitocondrial prefix
                                            keytype = "SYMBOL",
                                            column = "UNIPROT",
                                            multiVals = "first") %>%
      dplyr::filter(!is.na(UNIPROT))

  }else{
    uniprot_ids_df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
                                            keys = stringr::str_remove(gene_names, '^MT-'), #remove mitocondrial prefix
                                            keytype = "SYMBOL",
                                            column = "UNIPROT",
                                            multiVals = "first") %>%
      dplyr::filter(!is.na(UNIPROT))
  }

  # if the user provide an uniprot index column it will be added to the uniprot column of AnnotationDbi to do the query in the next step
  if (!is.null(uniprot_idx)) {
    uniprot_col_name <- names(omic_df)[uniprot_idx]
    user_uniprot_df <- omic_df %>%
      dplyr::select(SYMBOL = gene_name, UNIPROT = dplyr::all_of(uniprot_col_name)) %>%
      dplyr::filter(!is.na(UNIPROT), UNIPROT != "", !is.na(SYMBOL), SYMBOL != "") %>%
      dplyr::distinct() %>%
      dplyr::mutate(SYMBOL = as.character(SYMBOL),
                    UNIPROT = as.character(UNIPROT))
    
    uniprot_ids_df <- dplyr::bind_rows(uniprot_ids_df, user_uniprot_df) %>%
      dplyr::mutate(SYMBOL = as.character(SYMBOL),
                    UNIPROT = as.character(UNIPROT)) %>%
      dplyr::distinct()
  }
  
  uniprot_ids <- uniprot_ids_df$UNIPROT

  # Query UNIPROT
  message('  Querying UNIPROT API...')
  res <- SignalingProfiler::query_uniprot(id_input = uniprot_ids)

  # Get a subset of UNIPROT and Primary gene_name and Sequence
  res_sub <- res %>%
    dplyr::select(UNIPROT = Entry, SYMBOL = `Gene Names (primary)`, Sequence, Reviewed, Date_modified) %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::arrange(SYMBOL, dplyr::desc(Reviewed), dplyr::desc(Date_modified)) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice(if (any(Reviewed == "reviewed")) which(Reviewed == "reviewed") else 1) %>%
    dplyr::ungroup()

  # Create a table with UNIPROT, old symbol, primary gene symbol
  old_new_gn <- dplyr::inner_join(res_sub, uniprot_ids_df, by = c('UNIPROT'), suffix = c('.new', '.old'))

  # Update phosphoproteomics dataset
  omic_df_update <- dplyr::inner_join(
    omic_df,
    old_new_gn %>%
      dplyr::filter(!is.na(SYMBOL.new)) %>%
      dplyr::group_by(SYMBOL.old) %>%
      dplyr::reframe(
        UNIPROT = UNIPROT,
        gene_name.new = SYMBOL.new,
        Sequence = Sequence
      ),
    by = c('gene_name' = 'SYMBOL.old'),
    relationship = "many-to-many"
  )

  # Relocate UNIPROT, old and new gene symbol
  omic_df_update <- omic_df_update %>%
    dplyr::relocate(gene_name, gene_name.new, UNIPROT, Sequence) %>%
    tidyr::separate_rows(gene_name)

  # Check proteins having multiple primary new gene symbols
  omic_df_update <- dplyr::left_join(omic_df_update,
                                           omic_df_update %>% dplyr::count(gene_name) %>% dplyr::filter(n>1),
                                           by = 'gene_name')

  # Keep proteins that have an update gene symbol (n != 1) or if they have n > 1 updated genes
  # take the one like the original
  omic_df_update %>% dplyr::filter(is.na(n) | gene_name == gene_name.new) -> omic_df_update1

  omic_df_update1$n <- NULL
  omic_df_update <- omic_df_update1 %>% dplyr::distinct()

  omic_df_update$gene_name.new<- NULL

  # Remove Sequence if not desired
  if(!sequence == TRUE){
    omic_df_update$Sequence <- NULL
  }

  return(omic_df_update)
}

