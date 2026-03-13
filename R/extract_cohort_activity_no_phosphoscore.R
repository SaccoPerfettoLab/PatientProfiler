#' Extract Cohort Activity Using TFEA and KSEA Only
#'
#' This function loops over each patient in the cohort and extracts protein
#' activity for each patient using `extract_protein_activity_no_phosphoscore()`.
#' It uses transcriptomics and phosphoproteomics data to perform footprint-based
#' analysis only:
#' - TFEA on transcriptomics
#' - KSEA on phosphoproteomics
#'
#' Unlike `extract_cohort_activity()`, this version does **not** compute
#' phosphoscore.
#'
#' *Results output*:
#' Final results for each patient are saved as `.tsv` files inside the output
#' directory using the following format:
#' `Activity_Patient_{patient_id}.tsv`
#'
#' Each file contains these columns:
#' - `UNIPROT`: protein UniProt ID
#' - `gene_name`: gene symbol
#' - `mf`: molecular function
#' - `final_score`: inferred activity score
#' - `method`: method used to compute the score
#'
#' @param prot_dir A string (optional) indicating the proteomics data folder
#'   (default `NULL`).
#' @param trans_dir A string (optional) indicating the transcriptomics data folder
#'   (default `NULL`).
#' @param phospho_dir A string (optional) indicating the phosphoproteomics data
#'   folder (default `NULL`).
#'
#' @param tf_params A list of parameters for TFEA analysis. Available parameters are:
#'   - `omic_data`: dataframe of transcriptomic data for the analysis
#'   - `analysis`: analysis type; default: `'tfea'`
#'   - `organism`: reference species; default: `'human'`
#'   - `reg_minsize`: minimum regulon group size; default: `10`
#'   - `exp_sign`: boolean, whether to consider the sign of expression; default: `FALSE`
#'   - `collectri`: boolean, whether to use CollecTRI regulons; default: `FALSE`
#'   - `hypergeom_corr`: boolean, whether to apply hypergeometric correction; default: `TRUE`
#'   - `GO_annotation`: boolean, whether to include GO annotations; default: `TRUE`
#'   - `correct_proteomics`: boolean, whether to correct using proteomics; default: `FALSE`
#'   - `prot_df`: supporting proteomic dataframe
#'   - `custom`: boolean, whether to use custom regulons; default: `FALSE`
#'   - `custom_path`: path to custom regulons file; default: `NULL`
#'
#' @param kin_params A list of parameters for KSEA analysis. Available parameters are:
#'   - `omic_data`: dataframe of phosphoproteomic data for the analysis
#'   - `analysis`: analysis type; default: `'ksea'`
#'   - `organism`: reference species; default: `'human'`
#'   - `reg_minsize`: minimum regulator group size; default: `5`
#'   - `exp_sign`: boolean, whether to consider the sign of expression; default: `FALSE`
#'   - `integrated_regulons`: boolean, whether to use integrated regulons; default: `TRUE`
#'   - `hypergeom_corr`: boolean, whether to apply hypergeometric correction; default: `TRUE`
#'   - `GO_annotation`: boolean, whether to include GO annotations; default: `TRUE`
#'   - `correct_proteomics`: boolean, whether to correct using proteomics; default: `FALSE`
#'   - `prot_df`: supporting proteomic dataframe
#'   - `custom`: boolean, whether to use custom annotations; default: `FALSE`
#'   - `custom_path`: path to custom annotation file; default: `NULL`
#'
#' @param output_dir A string indicating the output directory that contains result files.
#'
#' @return No return value. Final results are saved as `.tsv` files in the output directory.
#'
#' @examples
#' tf_params <- list(reg_minsize = 5, collectri = TRUE)
#' kin_params <- list(exp_sign = TRUE, GO_annotation = FALSE)
#'
#' extract_cohort_activity_no_phosphoscore(
#'   prot_dir = "Prot_patients/",
#'   trans_dir = "Transc_patients/",
#'   phospho_dir = "Phospho_patients/",
#'   tf_params = tf_params,
#'   kin_params = kin_params,
#'   output_dir = "Activities"
#' )
#'
#' @export
extract_cohort_activity_no_phosphoscore <- function(
    prot_dir = NULL,
    trans_dir = NULL,
    phospho_dir = NULL,
    tf_params = list(),
    kin_params = list(),
    output_dir = "Activities"
) {
  
  message("Starting analysis")
  
  if (is.null(prot_dir) && is.null(trans_dir) && is.null(phospho_dir)) {
    stop("At least one of 'prot_dir', 'trans_dir', or 'phospho_dir' must be provided.")
  }
  
  prot_files <- if (!is.null(prot_dir)) {
    list.files(
      path = prot_dir,
      pattern = "^Prot_Patient_.*\\.tsv$",
      full.names = TRUE
    )
  } else {
    NULL
  }
  
  trans_files <- if (!is.null(trans_dir)) {
    list.files(
      path = trans_dir,
      pattern = "^Transc_Patient_.*\\.tsv$",
      full.names = TRUE
    )
  } else {
    NULL
  }
  
  phospho_files <- if (!is.null(phospho_dir)) {
    list.files(
      path = phospho_dir,
      pattern = "^Phospho_Patient_.*\\.tsv$",
      full.names = TRUE
    )
  } else {
    NULL
  }
  
  prot_ids <- if (!is.null(prot_files) && length(prot_files) > 0) {
    sapply(prot_files, extract_patient_id)
  } else {
    NULL
  }
  
  trans_ids <- if (!is.null(trans_files) && length(trans_files) > 0) {
    sapply(trans_files, extract_patient_id)
  } else {
    NULL
  }
  
  phospho_ids <- if (!is.null(phospho_files) && length(phospho_files) > 0) {
    sapply(phospho_files, extract_patient_id)
  } else {
    NULL
  }
  
  all_patient_ids <- unique(c(prot_ids, trans_ids, phospho_ids))
  total_patients <- length(all_patient_ids)
  
  if (total_patients == 0) {
    warning("No patient files found in the provided directories.")
    return(invisible(NULL))
  }
  
  patient_count <- 0
  
  for (patient_id in all_patient_ids) {
    patient_count <- patient_count + 1
    
    prot_file <- if (!is.null(prot_files) && patient_id %in% prot_ids) {
      prot_files[which(prot_ids == patient_id)[1]]
    } else {
      NULL
    }
    
    trans_file <- if (!is.null(trans_files) && patient_id %in% trans_ids) {
      trans_files[which(trans_ids == patient_id)[1]]
    } else {
      NULL
    }
    
    phospho_file <- if (!is.null(phospho_files) && patient_id %in% phospho_ids) {
      phospho_files[which(phospho_ids == patient_id)[1]]
    } else {
      NULL
    }
    
    message(
      paste("Processing patient", patient_count, "of", total_patients, ":", patient_id)
    )
    
    if (is.null(prot_file) && is.null(trans_file) && is.null(phospho_file)) {
      message(paste("Skipping patient", patient_id, "- no omics data available."))
      next
    }
    
    tryCatch({
      extract_protein_activity_no_phosphoscore(
        prot_file = prot_file,
        trans_file = trans_file,
        phospho_file = phospho_file,
        tf_params = tf_params,
        kin_params = kin_params,
        output_dir = output_dir
      )
      
      message(
        paste("Successfully processed patient", patient_id, "-", patient_count, "of", total_patients)
      )
    }, error = function(e) {
      message(
        paste("Error processing patient", patient_id, "-", patient_count, "of", total_patients, ":", e$message)
      )
    })
  }
  
  invisible(NULL)
}