#' extract_protein_activity
#'
#' This wrapper to SignalingProfiler preotein activity inference, extracts and analyzes the protein activity for a specific patient starting from omic data: transcriptomics, proteomics and phosphoproteomics.
#' It combines footprint based analysis (TFEA and KSEA), with phosphorylation scores (Phosphoscore), generating aggregated results for the patient.
#'
#' *Results output*:
#'    Final results for the patient will be saved as an tsv files inside the 'Activities' directory created by the function with this format:
#'    `Activity_patient_{patient_id}.tsv`.
#'    Each file contains these columns:
#'    - `UNIPROT`: protein Uniprot ID.
#'    - `gene_name`: name of the gene.
#'    - `mf`: molecular functions.
#'    - `final_score`: activity score.
#'    - `method`: how final_score is calculated.
#'
#' @param prot_file a string (optional) indicating the proteomics data file (default NULL)
#'
#' @param trans_file a string (optional) indicating the transcriptomics data file (default NULL)
#'
#' @param phospho_file a string (optional) indicating the phosphoproteomics data file (default NULL)
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
#'   **Note:** The `prot_df` parameter will only be used if a `prot_file` is specified, since the function
#'      automatically loads proteomics data only if `correct_proteomics` is `TRUE` and if it finds the specified file, so you don't need to change this
#'      parameter if you specified the proteomics file!
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
#'   **Note:** The `prot_df` parameter will only be used if a `prot_file` is specified, since the function
#'      automatically loads proteomics data only if `correct_proteomics` is `TRUE` and if it finds the specified file, so you don't need to change this
#'      parameter if you specified the proteomics file!
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
#' extract_protein_activity(prot_file = "Prot_patients/Prot_Patient_CPT000814.tsv",
#'                          trans_file = "Transc_patients/Transc_Patient_CPT000814.tsv",
#'                          phospho_file = "Phospho_patients/Phospho_Patient_CPT000814.tsv",
#'                          tf_params,
#'                          kin_params,
#'                          output_dir = "Activities")
#'
#' @export


extract_protein_activity <- function(
    prot_file = NULL,
    trans_file = NULL,
    phospho_file = NULL,
    tf_params = list(),
    kin_params = list(),
    phosphoscore_params = list(),
    phosphoscore_noseqwin_params = list(),
    output_dir = "Activities"
) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Initialize variables for omics data
  Prot_P <- NULL
  Trans_P <- NULL
  Phospho_P <- NULL

  # Load files if they exist
  if (!is.null(prot_file)) {
    Prot_P <- readr::read_tsv(prot_file,show_col_types = F)
  }

  if (!is.null(trans_file)) {
    Trans_P <- readr::read_tsv(trans_file,show_col_types = F)
  }

  if (!is.null(phospho_file)) {
    Phospho_P <- readr::read_tsv(phospho_file,show_col_types = F)
  }

  # Initialize activity results
  toy_activity_df <- data.frame()

  # Default TF and KSEA parameters (update dynamically based on available data)
  default_tf_params <- list(
    omic_data = Trans_P,
    analysis = 'tfea',
    organism = 'human',
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

  default_kin_params <- list(
    omic_data = Phospho_P,
    analysis = 'ksea',
    organism = 'human',
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

  default_phosphoscore_params <- list(
    phosphoproteomic_data = Phospho_P,
    organism = 'human',
    activatory = TRUE,
    GO_annotation = TRUE,
    custom = FALSE,
    custom_path = NULL
  )
  phosphoscore_params <- modifyList(default_phosphoscore_params, phosphoscore_params)

  default_phosphoscore_noseqwin_params <- list(
    phosphoproteomic_data = Phospho_P,
    organism = 'human',
    activatory = TRUE,
    GO_annotation = TRUE,
    custom = FALSE,
    custom_path = NULL
  )
  phosphoscore_noseqwin_params <- modifyList(default_phosphoscore_noseqwin_params, phosphoscore_noseqwin_params)

  # Perform TF analysis if transcriptomics data is available
  tf_activity_foot_1 <- if (!is.null(Trans_P)) do.call(SignalingProfiler::run_footprint_based_analysis, tf_params) else data.frame()

  # Perform KSEA analysis if phosphoproteomics data is available
  kin_phos_activity_foot_1 <- if (!is.null(Phospho_P)) do.call(SignalingProfiler::run_footprint_based_analysis, kin_params) else data.frame()

  # Conditional phosphoscore computation
  if (!is.null(Phospho_P)) {
    if (!all(is.na(Phospho_P$sequence_window))) {
      # With `sequence_window`
      phosphoscore_1 <- do.call(SignalingProfiler::phosphoscore_computation, phosphoscore_params)
    } else {
      # Without `sequence_window`
      phosphoscore_1 <- do.call(SignalingProfiler::phosphoscore_computation_aapos, phosphoscore_noseqwin_params)
    }

    combined_tf <- SignalingProfiler::combine_footprint_and_phosphoscore(
      footprint_output = tf_activity_foot_1,
      phosphoscore_df = phosphoscore_1,
      analysis = 'tfea'
    )

    combined_kin_phos <- SignalingProfiler::combine_footprint_and_phosphoscore(
      footprint_output = kin_phos_activity_foot_1,
      phosphoscore_df = phosphoscore_1,
      analysis = 'ksea'
    )

    toy_other <- phosphoscore_1 %>%
      dplyr::filter(mf == 'other') %>%
      dplyr::rename(final_score = phosphoscore) %>%
      dplyr::mutate(method = ifelse(!all(is.na(Phospho_P$sequence_window)), "PhosphoScore", "PhosphoScore_AAPOS"))

    toy_activity_df <- dplyr::bind_rows(combined_tf, combined_kin_phos, toy_other) %>%
      dplyr::select(UNIPROT, gene_name, mf, final_score, method)
  } else {
    # Use VIPER results if phosphoscore is not available
    toy_activity_df <- dplyr::bind_rows(tf_activity_foot_1, kin_phos_activity_foot_1) %>%
      dplyr::mutate(method = "VIPER") %>%
      dplyr::rename(final_score = weightedNES) %>%
      dplyr::select(UNIPROT, gene_name, mf, final_score, method)
  }

  # Save final results
  if (!is.null(toy_activity_df) && nrow(toy_activity_df) > 0) {
    patient_name <- extract_patient_id(basename(prot_file %||% trans_file %||% phospho_file))
    readr::write_tsv(toy_activity_df, paste0(output_dir, "/Activity_Patient_", patient_name, ".tsv"))
  } else {
    warning("No activity data could be extracted.")
  }
}
