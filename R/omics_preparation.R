#' Omics Data Preparation
#'
#' This function prepares transcriptomics, proteomics, and phosphoproteomics data for each patient
#' by creating separate CSV files for individual patient data. It performs various preprocessing steps,
#' including calculation of significance, log-p-values, and reordering of columns. The output files
#' are saved in separate directories for each type of omics data.
#'
#' @param df_tr_updated dataframe containing transcriptomics data. Must include 'gene_ID' and 'gene_name', with patient sample columns.
#' @param df_pr_updated dataframe containing proteomics data. Must include 'UNIPROT', 'gene_name', and patient sample columns.
#' @param transc_dir_name name of the directory for saving transcriptomics data outputs per patient.
#' @param prot_dir_name name of the directory for saving proteomics data outputs per patient.
#' @param phospho_dir_name name of the directory for saving phosphoproteomics data outputs per patient.


#'
#' @examples
#' omics_preparation(df_tr_updated = transcriptomics_data,
#'                   df_pr_updated = proteomics_data,
#'                   df_ph_updated = phosphoproteomics_data,
#'                   transc_dir_name = "Transc_patients",
#'                   prot_dir_name = "Prot_patients",
#'                   phospho_dir_name = "Phospho_patients")
#'
#' @export



omics_preparation <- function(
    df_tr_updated = NULL,
    df_pr_updated = NULL,
    df_ph_updated = NULL,
    transc_dir_name = "Transc_patients",
    prot_dir_name = "Prot_patients",
    phospho_dir_name = "Phospho_patients"
) {

  if (!is.null(df_tr_updated)) {
    message("Transcriptomic data patients preparation...")

    if (dir.exists(transc_dir_name)) {
      unlink(transc_dir_name, recursive = TRUE)  
    }
    dir.create(transc_dir_name)  
    
    for (c in 2:ncol(df_tr_updated)) {
      df_tr_updated_1 = df_tr_updated[, c(1, c)]
      df_tr_updated_2 = cbind(gene_ID = df_tr_updated[, 1], df_tr_updated_1)
      colnames(df_tr_updated_2) = c("gene_ID", "gene_name", "difference")
      df_tr_updated_3 = df_tr_updated_2 %>%
        dplyr::mutate(significant = ifelse(difference >= 1.96 | difference <= -1.96, "+", NA))
      df_tr_updated_4 = df_tr_updated_3 %>%
        dplyr::mutate(logpval = 2 * pnorm(as.numeric(difference), lower.tail = FALSE))
      df_tr_updated_5 = cbind(tstat = df_tr_updated_4[, 3], df_tr_updated_4)
      df_tr_updated_5 = df_tr_updated_5 %>% relocate(gene_ID, gene_name, difference, tstat, logpval, significant)
      readr::write_tsv(df_tr_updated_5, paste0(transc_dir_name, "/Transc_Patient_", as.character(colnames(df_tr_updated)[c]), ".tsv"))
    }
    message("Done!")
  }

  if (!is.null(df_pr_updated)) {
    message("Proteomic data patients preparation...")

    if (dir.exists(prot_dir_name)) {
      unlink(prot_dir_name, recursive = TRUE)
    }
    dir.create(prot_dir_name)  
    
    df_pr_updated <- df_pr_updated %>%
      dplyr::select(UNIPROT, everything())

    for (i in 3:length(df_pr_updated)) {
      pt1 = df_pr_updated[, c(1:2)]
      patients = cbind(pt1, df_pr_updated[, i])
      colnames(patients) = c("UNIPROT", "gene_name", "difference")
      patients = patients %>%
        dplyr::mutate(significant = ifelse(difference >= 1.96 | difference <= -1.96, "+", NA)) %>%
        dplyr::mutate(logpval = 2 * pnorm(as.numeric(difference), lower.tail = FALSE))
      readr::write_tsv(patients, paste0(prot_dir_name, "/Prot_Patient_", as.character(colnames(df_pr_updated)[i]), ".tsv"))
    }
    message("Done!")
  }

  if (!is.null(df_ph_updated)) {
    message("Phosphoproteomics data patients preparation...")

    if (dir.exists(phospho_dir_name)) {
      unlink(phospho_dir_name, recursive = TRUE)  
    }
    dir.create(phospho_dir_name)  
    
    
    if (!"sequence_window" %in% names(df_ph_updated)) {
      df_ph_updated$sequence_window <- NA
    }

    column_order <- c("UNIPROT", "aminoacid", "position", "gene_name", "sequence_window")
    other_columns <- setdiff(names(df_ph_updated), column_order)
    new_order <- c(column_order, other_columns)
    df_ph_updated <- df_ph_updated[, new_order]

    for (i in 6:length(df_ph_updated)) {
      pt1 <- df_ph_updated[, 1:5]
      patients = cbind(pt1, df_ph_updated[, i])
      colnames(patients)[1:6] = c("UNIPROT", "aminoacid", "position", "gene_name", "sequence_window", "difference")
      patients = patients %>%
        dplyr::mutate(logpval = 2 * pnorm(difference, lower.tail = FALSE))
      patients = patients %>%
        dplyr::mutate(significant = ifelse(difference >= 1.96 | difference <= -1.96, "+", NA))
      readr::write_tsv(patients, paste0(phospho_dir_name, "/Phospho_Patient_", as.character(colnames(df_ph_updated)[i]), ".tsv"))
    }
    message("Done!")
  }
}
