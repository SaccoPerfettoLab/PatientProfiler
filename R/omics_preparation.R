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
  suppressPackageStartupMessages({
    library(dplyr); library(stringr); library(readr); library(tidyr)
  })
  
  # helper: trim -> "" -> NA
  .norm <- function(x) {
    if (!is.character(x)) x <- as.character(x)
    x <- str_trim(x)
    na_if(x, "")
  }
  
  # ============= TRANSCRIPTOMICS =============
  if (!is.null(df_tr_updated)) {
    message("Transcriptomic data patients preparation...")
    
    # sanitize IDs and drop rows with NA IDs (gene_ID or gene_name)
    if (!all(c("gene_ID","gene_name") %in% names(df_tr_updated))) {
      stop("df_tr_updated must contain columns: gene_ID, gene_name")
    }
    df_tr_updated <- df_tr_updated %>%
      mutate(
        gene_ID   = .norm(gene_ID),
        gene_name = .norm(gene_name)
      ) %>%
      filter(!is.na(gene_ID), !is.na(gene_name))
    
    if (dir.exists(transc_dir_name)) unlink(transc_dir_name, recursive = TRUE)
    dir.create(transc_dir_name, recursive = TRUE)
    
    for (c in 2:ncol(df_tr_updated)) {
      df_tr_updated_1 <- df_tr_updated[, c(1, c)]
      df_tr_updated_2 <- cbind(gene_ID = df_tr_updated[, 1], df_tr_updated_1)
      colnames(df_tr_updated_2) <- c("gene_ID", "gene_name", "difference")
      df_tr_updated_3 <- df_tr_updated_2 %>%
        mutate(significant = ifelse(difference >= 1.96 | difference <= -1.96, "+", NA))
      df_tr_updated_4 <- df_tr_updated_3 %>%
        mutate(logpval = 2 * pnorm(as.numeric(difference), lower.tail = FALSE))
      df_tr_updated_5 <- cbind(tstat = df_tr_updated_4[, 3], df_tr_updated_4) %>%
        relocate(gene_ID, gene_name, difference, tstat, logpval, significant)
      
      write_tsv(
        df_tr_updated_5,
        file.path(transc_dir_name, paste0("Transc_Patient_", colnames(df_tr_updated)[c], ".tsv"))
      )
    }
    message("Done!")
  }
  
  # ============= PROTEOMICS =============
  if (!is.null(df_pr_updated)) {
    message("Proteomic data patients preparation...")
    
    # sanitize IDs and drop rows with NA IDs (UNIPROT or gene_name)
    if (!all(c("UNIPROT","gene_name") %in% names(df_pr_updated))) {
      stop("df_pr_updated must contain columns: UNIPROT, gene_name")
    }
    df_pr_updated <- df_pr_updated %>%
      mutate(
        UNIPROT   = .norm(UNIPROT),
        gene_name = .norm(gene_name)
      ) %>%
      filter(!is.na(UNIPROT), !is.na(gene_name)) %>%
      select(UNIPROT, everything())  # keep UNIPROT first
    
    if (dir.exists(prot_dir_name)) unlink(prot_dir_name, recursive = TRUE)
    dir.create(prot_dir_name, recursive = TRUE)
    
    # assume columns: 1=UNIPROT, 2=gene_name, 3..=patients
    for (i in 3:ncol(df_pr_updated)) {
      pt1 <- df_pr_updated[, c(1:2)]
      patients <- cbind(pt1, df_pr_updated[, i])
      colnames(patients) <- c("UNIPROT", "gene_name", "difference")
      patients <- patients %>%
        mutate(
          significant = ifelse(difference >= 1.96 | difference <= -1.96, "+", NA),
          logpval     = 2 * pnorm(as.numeric(difference), lower.tail = FALSE)
        )
      
      write_tsv(
        patients,
        file.path(prot_dir_name, paste0("Prot_Patient_", colnames(df_pr_updated)[i], ".tsv"))
      )
    }
    message("Done!")
  }
  
  # ============= PHOSPHOPROTEOMICS =============
  if (!is.null(df_ph_updated)) {
    message("Phosphoproteomics data patients preparation...")
    
    # sanitize IDs and drop rows with NA IDs (UNIPROT or gene_name)
    if (!all(c("UNIPROT","gene_name") %in% names(df_ph_updated))) {
      stop("df_ph_updated must contain columns: UNIPROT, gene_name")
    }
    df_ph_updated <- df_ph_updated %>%
      mutate(
        UNIPROT   = .norm(UNIPROT),
        gene_name = .norm(gene_name)
      ) %>%
      filter(!is.na(UNIPROT), !is.na(gene_name))
    
    if (dir.exists(phospho_dir_name)) unlink(phospho_dir_name, recursive = TRUE)
    dir.create(phospho_dir_name, recursive = TRUE)
    
    # columns to keep first
    if ("sequence_window" %in% names(df_ph_updated)) {
      column_order <- c("UNIPROT", "aminoacid", "position", "gene_name", "sequence_window")
    } else {
      column_order <- c("UNIPROT", "aminoacid", "position", "gene_name")
    }
    other_columns <- setdiff(names(df_ph_updated), column_order)
    start_pos <- ncol(df_ph_updated) - length(other_columns) + 1
    new_order <- c(column_order, other_columns)
    df_ph_updated <- df_ph_updated[, new_order]
    
    for (i in start_pos:ncol(df_ph_updated)) {
      pt1 <- df_ph_updated[, seq_along(column_order)]
      patients <- cbind(pt1, df_ph_updated[, i])
      colnames(patients)[1:(length(column_order) + 1)] <- c(column_order, "difference")
      patients <- patients %>%
        mutate(
          logpval = 2 * pnorm(difference, lower.tail = FALSE),
          significant = ifelse(difference >= 1.96 | difference <= -1.96, "+", NA)
        )
      
      write_tsv(
        patients,
        file.path(phospho_dir_name, paste0("Phospho_Patient_", colnames(df_ph_updated)[i], ".tsv"))
      )
    }
    message("Done!")
  }
}
