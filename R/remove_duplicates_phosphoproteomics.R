#' remove_duplicates_phosphoproteomics
#'
#' This function processes a phosphoproteomics data frame by removing duplicates,
#' aggregating metadata, and calculating average values for patient columns.
#'
#' @param df dataframe (phosphoproteomics data). It must include 'UNIPROT', 'gene_name', 'aminoacid
#'
#' @return a cleaned data frame with aggregated metadata and averaged values.
#'
#'
#' @examples
#' phosphoproteomics_data <- data.frame(
#'   gene_name = c("Gene1", "Gene1", "Gene2"),
#'   aminoacid = c("A", "A", "B"),
#'   position = c(1, 1, 2),
#'   sequence_window = c("AAACCC", "AAACCC", "BBBCCC"),
#'   UNIPROT = c("P12345", "P12345", "P67890"),
#'   patient_1 = c(0.5, 0.6, 0.7),
#'   patient_2 = c(0.4, 0.3, 0.8)
#' )
#'
#' cleaned_data <- remove_duplicates_phosphoproteomics(phosphoproteomics_data)
#' head(cleaned_data)
#'
remove_duplicates_phosphoproteomics <- function(df) {

  df_unique <- unique(df)

  if ("sequence_window" %in% colnames(df)) {
    df_unique$key <- paste0(df_unique$gene_name, "_", df_unique$aminoacid, "_",
                            df_unique$position, "_", df_unique$sequence_window)
  } else {
    df_unique$key <- paste0(df_unique$gene_name, "_", df_unique$aminoacid, "_", df_unique$position)
  }

  df_uni <- df_unique %>%
    dplyr::group_by(key) %>%
    dplyr::summarize(UNIPROT = paste(unique(UNIPROT), collapse = ";"), .groups = 'drop')

  df_unique <- df_unique %>% dplyr::select(-UNIPROT)

  df_final <- df_unique %>%
    dplyr::left_join(df_uni, by = "key") %>% dplyr::relocate(UNIPROT)

  aggregate_uniprot <- function(df) {
    df %>%
      dplyr::group_by(gene_name) %>%
      dplyr::summarize(UNIPROT_combined = paste(unique(unlist(strsplit(UNIPROT, ";"))), collapse = ";"), .groups = 'drop')
  }

  df_final <- df_final %>%
    dplyr::left_join(aggregate_uniprot(df_final), by = "gene_name") %>%
    dplyr::mutate(UNIPROT = UNIPROT_combined) %>%
    dplyr::select(-UNIPROT_combined)

  df_final <- unique(df_final)

  if ("sequence_window" %in% colnames(df)) {
    df_final$key <- paste0(df_final$UNIPROT, "_", df_final$aminoacid, "_",
                           df_final$position, "_", df_final$sequence_window)
  } else {
    df_final$key <- paste0(df_final$UNIPROT, "_", df_final$aminoacid, "_", df_final$position)
  }

  df_uni <- df_final %>%
    dplyr::group_by(key) %>%
    dplyr::summarize(gene_name = paste(unique(gene_name), collapse = ";"), .groups = 'drop')

  df_final <- df_final %>% dplyr::select(-gene_name)

  df_final <- df_final %>%
    dplyr::left_join(df_uni, by = "key") %>% dplyr::relocate(gene_name)

  aggregate_gene <- function(df) {
    df %>%
      dplyr::group_by(UNIPROT) %>%
      dplyr::summarize(gene_combined = paste(unique(unlist(strsplit(gene_name, ";"))), collapse = ";"), .groups = 'drop')
  }

  df_final <- df_final %>%
    dplyr::left_join(aggregate_gene(df_final), by = "UNIPROT") %>%
    dplyr::mutate(gene_name = gene_combined) %>%
    dplyr::select(-gene_combined)

  df_final <- unique(df_final)

  if ("sequence_window" %in% colnames(df)) {
    df_final_clean <- df_final %>%
      dplyr::filter(!is.na(sequence_window))
  } else {
    df_final_clean <- df_final
  }

  if ("sequence_window" %in% colnames(df)) {
    df_final_clean$key <- paste0(df_final_clean$gene_name, "_", df_final_clean$aminoacid, "_",
                                 df_final_clean$position, "_", df_final_clean$sequence_window)
  } else {
    df_final_clean$key <- paste0(df_final_clean$gene_name, "_", df_final_clean$aminoacid, "_", df_final_clean$position)
  }

  df_final_clean <- df_final_clean %>% dplyr::relocate(key)


  if ("sequence_window" %in% colnames(df_final_clean)) {

    selected_cols <- colnames(df_final_clean)[!(colnames(df_final_clean) %in% c("key", "gene_name", "UNIPROT", "aminoacid", "position", "sequence_window"))]

  }else{

    selected_cols <- colnames(df_final_clean)[!(colnames(df_final_clean) %in% c("key", "gene_name", "UNIPROT", "aminoacid", "position"))]

  }


  df_final_clean2 <- df_final_clean %>%
    dplyr::group_by(key) %>%
    dplyr::summarize(dplyr::across(selected_cols, mean, na.rm = TRUE), .groups = 'drop')


  if ("sequence_window" %in% colnames(df_final_clean)) {
    metadata_columns <- df_final_clean %>%
      dplyr::select(key, gene_name, UNIPROT, aminoacid, position, sequence_window) %>%
      dplyr::distinct(key, .keep_all = TRUE)

  }else {
    metadata_columns <- df_final_clean %>%
      dplyr::select(key, gene_name, UNIPROT, aminoacid, position) %>%
      dplyr::distinct(key, .keep_all = TRUE)

  }

  if ("sequence_window" %in% colnames(df_final_clean)) {

    df_final_clean2 <- df_final_clean2 %>%
      dplyr::left_join(metadata_columns, by = "key") %>%
      dplyr::select(-key) %>%
      dplyr::relocate(UNIPROT, aminoacid, position, gene_name, sequence_window)

  } else {

    df_final_clean2 <- df_final_clean2 %>%
      dplyr::left_join(metadata_columns, by = "key") %>%
      dplyr::select(-key) %>%
      dplyr::relocate(UNIPROT, aminoacid, position, gene_name)
  }


  if ("sequence_window" %in% colnames(df)) {
    df_final_clean2 <- df_final_clean2 %>% dplyr::relocate(sequence_window, .after = UNIPROT)
  }

  df_final_clean2 <- df_final_clean2 %>%
    dplyr::mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))

  return(df_final_clean2)
}
