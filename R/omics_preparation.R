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
#' @param overwrite logical; if FALSE (default) stops when an output file already exists.
#' @param dir_action what to do if directory already exists: "keep" (default) or "timestamp".
#'   If "timestamp", a new directory with suffix _YYYYmmdd_HHMMSS is created.
#'
#' @export
omics_preparation <- function(
    df_tr_updated = NULL,
    df_pr_updated = NULL,
    df_ph_updated = NULL,
    transc_dir_name = "Transc_patients",
    prot_dir_name = "Prot_patients",
    phospho_dir_name = "Phospho_patients",
    overwrite = FALSE,
    dir_action = c("keep", "timestamp")
) {
  suppressPackageStartupMessages({
    library(dplyr); library(stringr); library(readr); library(tidyr)
  })
  
  dir_action <- match.arg(dir_action)
  
  # helper: trim -> "" -> NA
  .norm <- function(x) {
    if (!is.character(x)) x <- as.character(x)
    x <- str_trim(x)
    na_if(x, "")
  }
  
  # helper: safe dir create (no delete)
  .prepare_dir <- function(dir_name) {
    if (!dir.exists(dir_name)) {
      dir.create(dir_name, recursive = TRUE)
      return(dir_name)
    }
    if (dir_action == "keep") {
      return(dir_name)
    }
    # timestamp mode
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    new_dir <- paste0(dir_name, "_", ts)
    dir.create(new_dir, recursive = TRUE)
    new_dir
  }
  
  # helper: safe write
  .write_tsv_safe <- function(x, path) {
    if (file.exists(path) && !isTRUE(overwrite)) {
      stop("Output file already exists and overwrite=FALSE: ", path, call. = FALSE)
    }
    readr::write_tsv(x, path)
  }
  
  # ============= TRANSCRIPTOMICS =============
  if (!is.null(df_tr_updated)) {
    message("Transcriptomic data patients preparation...")
    
    if (!("gene_name" %in% names(df_tr_updated))) {
      stop("df_tr_updated must contain column: gene_name")
    }
    
    df_tr_updated <- df_tr_updated %>%
      mutate(gene_name = .norm(gene_name)) %>%
      filter(!is.na(gene_name))
    
    out_dir <- .prepare_dir(transc_dir_name)
    
    # Original logic: loop from column 2 to end (assumes col1=gene_ID, col2..=patients OR col2=gene_name + patients)
    # In your original code you bind gene_ID=df_tr_updated[,1] and also take df_tr_updated_1 = cols (1,c),
    # then rename to gene_ID, gene_name, difference. That implies:
    # - col1 = gene_ID
    # - col(c) = patient column (starting at 2)
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
      
      out_path <- file.path(out_dir, paste0("Transc_Patient_", colnames(df_tr_updated)[c], ".tsv"))
      .write_tsv_safe(df_tr_updated_5, out_path)
    }
    
    message("Done!")
  }
  
  # ============= PROTEOMICS =============
  if (!is.null(df_pr_updated)) {
    message("Proteomic data patients preparation...")
    
    if (!all(c("UNIPROT","gene_name") %in% names(df_pr_updated))) {
      stop("df_pr_updated must contain columns: UNIPROT, gene_name")
    }
    
    df_pr_updated <- df_pr_updated %>%
      mutate(
        UNIPROT   = .norm(UNIPROT),
        gene_name = .norm(gene_name)
      ) %>%
      filter(!is.na(UNIPROT), !is.na(gene_name)) %>%
      select(UNIPROT, everything())
    
    out_dir <- .prepare_dir(prot_dir_name)
    
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
      
      out_path <- file.path(out_dir, paste0("Prot_Patient_", colnames(df_pr_updated)[i], ".tsv"))
      .write_tsv_safe(patients, out_path)
    }
    
    message("Done!")
  }
  
  # ============= PHOSPHOPROTEOMICS =============
  if (!is.null(df_ph_updated)) {
    message("Phosphoproteomics data patients preparation...")
    
    if (!all(c("UNIPROT","gene_name") %in% names(df_ph_updated))) {
      stop("df_ph_updated must contain columns: UNIPROT, gene_name")
    }
    
    df_ph_updated <- df_ph_updated %>%
      mutate(
        UNIPROT   = .norm(UNIPROT),
        gene_name = .norm(gene_name)
      ) %>%
      filter(!is.na(UNIPROT), !is.na(gene_name))
    
    out_dir <- .prepare_dir(phospho_dir_name)
    
    # columns to keep first (same logic as your function)
    if ("sequence_window" %in% names(df_ph_updated)) {
      column_order <- c("UNIPROT", "aminoacid", "position", "gene_name", "sequence_window")
    } else {
      column_order <- c("UNIPROT", "aminoacid", "position", "gene_name")
    }
    
    # Keep behavior identical: reorder putting column_order first, then all other columns.
    other_columns <- setdiff(names(df_ph_updated), column_order)
    new_order <- c(column_order, other_columns)
    
    # (Optional safety) if any required column missing, give clearer error
    missing_cols <- setdiff(column_order, names(df_ph_updated))
    if (length(missing_cols) > 0) {
      stop("df_ph_updated is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    df_ph_updated <- df_ph_updated[, new_order]
    
    # patient columns are the "other_columns" portion, i.e. after the metadata block
    start_pos <- length(column_order) + 1
    
    for (i in start_pos:ncol(df_ph_updated)) {
      pt1 <- df_ph_updated[, seq_along(column_order)]
      patients <- cbind(pt1, df_ph_updated[, i])
      
      colnames(patients)[1:(length(column_order) + 1)] <- c(column_order, "difference")
      
      patients <- patients %>%
        mutate(
          logpval = 2 * pnorm(as.numeric(difference), lower.tail = FALSE),
          significant = ifelse(difference >= 1.96 | difference <= -1.96, "+", NA)
        )
      
      out_path <- file.path(out_dir, paste0("Phospho_Patient_", colnames(df_ph_updated)[i], ".tsv"))
      .write_tsv_safe(patients, out_path)
    }
    
    message("Done!")
  }
  
  invisible(TRUE)
}
