#' Extract Protein Activity Using TFEA and KSEA Only
#'
#' This wrapper around SignalingProfiler performs protein activity inference for a
#' specific patient starting from omic data: transcriptomics, proteomics, and
#' phosphoproteomics.
#'
#' Unlike `extract_protein_activity()`, this version does **not** compute
#' phosphoscore. Phosphoproteomic data are analyzed only through KSEA.
#'
#' *Results output*:
#' Final results for the patient are saved as `.tsv` files inside the output
#' directory, using the following format:
#' `Activity_patient_{patient_id}.tsv`.
#'
#' Each file contains these columns:
#' - `UNIPROT`: protein UniProt ID
#' - `gene_name`: gene symbol
#' - `mf`: molecular function
#' - `final_score`: inferred activity score
#' - `method`: method used to compute the score
#'
#' @param prot_file A string (optional) indicating the proteomics data file
#'   (default `NULL`).
#' @param trans_file A string (optional) indicating the transcriptomics data file
#'   (default `NULL`).
#' @param phospho_file A string (optional) indicating the phosphoproteomics data
#'   file (default `NULL`).
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
#' @param output_dir A string indicating the output directory containing result files.
#'
#' @return No return value. Final results are saved as `.tsv` files in the output directory.
#'
#' @examples
#' tf_params <- list(reg_minsize = 5, collectri = TRUE)
#' kin_params <- list(exp_sign = TRUE, GO_annotation = FALSE)
#'
#' extract_protein_activity_ksea_only(
#'   prot_file = "Prot_patients/Prot_Patient_CPT000814.tsv",
#'   trans_file = "Transc_patients/Transc_Patient_CPT000814.tsv",
#'   phospho_file = "Phospho_patients/Phospho_Patient_CPT000814.tsv",
#'   tf_params = tf_params,
#'   kin_params = kin_params,
#'   output_dir = "Activities"
#' )
#'
#' @export
extract_protein_activity_no_phosphoscore <- function(
    prot_file = NULL,
    trans_file = NULL,
    phospho_file = NULL,
    tf_params = list(),
    kin_params = list(),
    output_dir = "Activities"
) {
  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Initialize omics data objects
  Prot_P <- NULL
  Trans_P <- NULL
  Phospho_P <- NULL
  
  # Load input files if provided
  if (!is.null(prot_file)) {
    Prot_P <- readr::read_tsv(prot_file, show_col_types = FALSE)
  }
  
  if (!is.null(trans_file)) {
    Trans_P <- readr::read_tsv(trans_file, show_col_types = FALSE)
  }
  
  if (!is.null(phospho_file)) {
    Phospho_P <- readr::read_tsv(phospho_file, show_col_types = FALSE)
  }
  
  # Default TFEA parameters
  default_tf_params <- list(
    omic_data = Trans_P,
    analysis = "tfea",
    organism = "human",
    reg_minsize = 10,
    exp_sign = FALSE,
    collectri = FALSE,
    hypergeom_corr = TRUE,
    GO_annotation = TRUE,
    correct_proteomics = FALSE,
    prot_df = Prot_P,
    custom = FALSE,
    custom_path = NULL
  )
  
  tf_params <- modifyList(default_tf_params, tf_params)
  
  # Default KSEA parameters
  default_kin_params <- list(
    omic_data = Phospho_P,
    analysis = "ksea",
    organism = "human",
    reg_minsize = 5,
    exp_sign = FALSE,
    integrated_regulons = TRUE,
    hypergeom_corr = TRUE,
    GO_annotation = TRUE,
    correct_proteomics = FALSE,
    prot_df = Prot_P,
    custom = FALSE,
    custom_path = NULL
  )
  
  kin_params <- modifyList(default_kin_params, kin_params)
  
  # Run TFEA if transcriptomics data are available
  tf_activity_foot_1 <- if (!is.null(Trans_P)) {
    do.call(SignalingProfiler::run_footprint_based_analysis, tf_params)
  } else {
    data.frame()
  }
  
  # Run KSEA if phosphoproteomics data are available
  kin_phos_activity_foot_1 <- if (!is.null(Phospho_P)) {
    do.call(SignalingProfiler::run_footprint_based_analysis, kin_params)
  } else {
    data.frame()
  }
  
  # Combine only footprint-based results
  toy_activity_df <- dplyr::bind_rows(tf_activity_foot_1, kin_phos_activity_foot_1)
  
  # Standardize output format
  if (!is.null(toy_activity_df) && nrow(toy_activity_df) > 0) {
    toy_activity_df <- toy_activity_df %>%
      dplyr::mutate(method = "VIPER") %>%
      dplyr::rename(final_score = weightedNES) %>%
      dplyr::select(UNIPROT, gene_name, mf, final_score, method)
    
    patient_name <- extract_patient_id(basename(prot_file %||% trans_file %||% phospho_file))
    readr::write_tsv(
      toy_activity_df,
      file.path(output_dir, paste0("Activity_Patient_", patient_name, ".tsv"))
    )
  } else {
    warning("No activity data could be extracted.")
  }
}