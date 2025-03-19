#' extract_cohort_activity
#'
#' This function loop over each patient of the cohort and extracts and analyzes the protein activity for each patient (using *extract_protein_activitiìy*
#' function) starting from omic data: transcriptomics, proteomics and phosphoproteomics. It combines footprint based analysis (TFEA and KSEA),
#' with phosphorylation scores (Phosphoscore), generating aggregated results for each patient.
#'
#' *Results output*:
#'    Final results for each patient will be saved as an tsv files inside the 'Activities' directory created by the function with this format:
#'    `Activity_patient_{patient_id}.tsv`.
#'    Each file contains these columns:
#'    - `UNIPROT`: protein Uniprot ID.
#'    - `gene_name`: name of the gene.
#'    - `mf`: molecular functions.
#'    - `final_score`: activity score.
#'    - `method`: how final_score is calculated.
#'
#' @param prot_dir a string (optional) indicating the proteomics data folder (default NULL)
#'
#' @param trans_dir a string (optional) indicating the transcriptomics data folder (default NULL)
#'
#' @param phospho_dir a string (optional) indicating the phosphoproteomics data folder (default NULL)
#'
#' @param tf_params A list of parameters for TFEA analysis. Available parameters are:
#'   - `omic_data`: dataframe of transcriptomic data for the analysis.
#'   - `analysis`: analysis type; default: `'tfea'`.
#'   - `organism`: reference species; default: `'human'`.
#'   - `reg_minsize`: minimum regulons group size; default: `10`.
#'   - `exp_sign`: boolean, indicates whether to consider the sign of the expression; default: `FALSE`.
#'   - `collectri`: boolean, indicates whether to collect regulator-target relationships; default: `FALSE`.
#'   - `hypergeom_corr`: boolean, whether to correct results for hypergeometry; default: `TRUE`.
#'   - `GO_annotation`: boolean, whether to include GO annotations; default: `TRUE`.
#'   - `correct_proteomics`: boolean, whether to correct with proteomic data; default: `FALSE`.
#'   - `prot_df`: supporting proteomic data frame.
#'   **Note:** The `prot_df` parameter will only be used if a `prot_dir` is specified, since the function
#'      automatically loads proteomics data only if `correct_proteomics` is `TRUE` and if it finds files in the specified folder, so you don't need to change this
#'      parameter if you specified the proteomics folder!
#'   - `custom`: boolean, whether to use custom regulons; default: `FALSE`.
#'   - `custom_path`: path to custom regulons file; default: `NULL`.
#'
#' @param kin_params A list of parameters for the KSEA analysis. Available parameters are:
#'   - `omic_data`: dataframe of phosphoproteomic data for analysis.
#'   - `analysis`: analysis type; default: `'ksea'`.
#'   - `organism`: reference species; default: `'human'`.
#'   - `reg_minsize`: minimum size of the regulator group; default: `5`.
#'   - `exp_sign`: boolean, whether to consider the sign of the expression; default: `FALSE`.
#'   - `integrated_regulons`: boolean, whether to consider integrated regulators; default: `TRUE`.
#'   - `hypergeom_corr`: boolean, whether to correct the results for hypergeometry; default: `TRUE`.
#'   - `GO_annotation`: boolean, whether to include GO annotations; default: `TRUE`.
#'   - `correct_proteomics`: boolean, whether to correct with proteomic data; default: `FALSE`.
#'   - `prot_df`: supporting proteomic data frame.
#'      **Note:** The `prot_df` parameter will only be used if a `prot_dir` is specified, since the function
#'      automatically loads proteomics data only if `correct_proteomics` is `TRUE` and if it finds files in the specified folder, so you don't need to change this
#'      parameter if you specified the proteomics folder!
#'   - `custom`: boolean, whether to use custom annotations; default: `FALSE`.
#'   - `custom_path`: path to the custom annotation file; default: `NULL`.
#'
#' @param phosphoscore_params A list of parameters for the calculation of *phosphoscore*. Available parameters are:
#'   - `phosphoproteomic_data`: dataframe of phosphoproteomic data for the calculation of *phosphoscore*.
#'   - `organism`: reference species; default: `'human'`.
#'   - `activatory`: boolean, whether to consider activatory changes; default: `TRUE`.
#'   - `GO_annotation`: boolean, whether to include GO annotations; default: `TRUE`.
#'   - `custom`: boolean, whether to use custom annotations; default: `FALSE`.
#'   - `custom_path`: path to the custom annotation file; default: `NULL`.
#'
#' @param phosphoscore_noseqwin_params A list of parameters for calculating *phosphoscore* without *sequence window*.
#' The available parameters are the same as those in `phosphoscore_params` and are overridden if specified.
#'
#' @param  output_dir a string indicating the output directory that contains results files
#'
#' @return No return value. The final results are saved in .tsv files in the output directory.
#'
#' @examples
#' tf_params <- list(reg_minsize = 5, collectri = TRUE)
#' kin_params <- list(exp_sign = TRUE, GO_annotation = FALSE)
#' extract_cohort_activity(prot_dir = "prot_dir/", trans_dir = "transc_dir/",
#'                            phospho_dir = "phospho_dir/", tf_params, kin_params)
#'
#'
#' @export

extract_cohort_activity <- function(
    prot_dir = NULL,
    trans_dir = "PatientProfiler_test/testing/Transc_patients/",
    phospho_dir = NULL,
    tf_params = list(),
    kin_params = list(),
    phosphoscore_params = list(),
    phosphoscore_noseqwin_params = list(),
    output_dir = "Activities"
) {
  
  print("Starting analysis")
  
  if (is.null(prot_dir) && is.null(trans_dir) && is.null(phospho_dir)) {
    stop("At least one of 'prot_dir', 'trans_dir', or 'phospho_dir' must be provided.")
  }
  
  prot_files <- if (!is.null(prot_dir)) {
    list.files(path = prot_dir, pattern = "^Prot_Patient_.*\\.tsv$", full.names = TRUE)
  } else NULL
  
  trans_files <- if (!is.null(trans_dir)) {
    list.files(path = trans_dir, pattern = "^Transc_Patient_.*\\.tsv$", full.names = TRUE)
  } else NULL
  
  phospho_files <- if (!is.null(phospho_dir)) {
    list.files(path = phospho_dir, pattern = "^Phospho_Patient_.*\\.tsv$", full.names = TRUE)
  } else NULL
  
  prot_ids <- if (!is.null(prot_files)) sapply(prot_files, extract_patient_id) else NULL
  trans_ids <- if (!is.null(trans_files)) sapply(trans_files, extract_patient_id) else NULL
  phospho_ids <- if (!is.null(phospho_files)) sapply(phospho_files, extract_patient_id) else NULL
  
  all_patient_ids <- unique(c(prot_ids, trans_ids, phospho_ids))
  
  patient_count <- 0
  total_patients <- length(all_patient_ids)
  
  for (patient_id in all_patient_ids) {
    patient_count <- patient_count + 1
    
    prot_file <- if (!is.null(prot_files)) prot_files[which(prot_ids == patient_id)[1]] else NULL
    trans_file <- if (!is.null(trans_files)) trans_files[which(trans_ids == patient_id)[1]] else NULL
    phospho_file <- if (!is.null(phospho_files)) phospho_files[which(phospho_ids == patient_id)[1]] else NULL
    
    message(paste("Processing patient", patient_count, "of", total_patients, ":", patient_id))
    
    if (is.null(prot_file) && is.null(trans_file) && is.null(phospho_file)) {
      message(paste("Skipping patient", patient_id, "- no omics data available."))
      next
    }
    
    tryCatch({
      extract_protein_activity(
        prot_file = prot_file,
        trans_file = trans_file,
        phospho_file = phospho_file,
        tf_params = tf_params,
        kin_params = kin_params,
        phosphoscore_params = phosphoscore_params,
        phosphoscore_noseqwin_params = phosphoscore_noseqwin_params,
        output_dir = output_dir
      )
      
      
      message(paste("Successfully processed patient", patient_id, "-", patient_count, "of", total_patients))
    }, error = function(e) {
      message(paste("Error processing patient", patient_id, "-", patient_count, "of", total_patients, ":", e$message))
    })
  }
}
