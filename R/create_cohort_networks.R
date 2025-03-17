
#' create_cohort_networks
#'
#' This function loops on the folder of omics and activity data of a cohort of patients
#' creating a network of each of the patients from mutations to desired phenotypes.
#' For the files generated for each patient, refer to *create_network* function documentation.
#'
#' @param prot_dir string, path to cohort proteomics data folder, default `NULL`.
#' @param phospho_dir  string, path to cohort proteomics data folder, default `NULL`.
#' @param act_dir  string, path to cohort proteomics data folder, default `NULL`.
#' @param mut_file  string, path to cohort mutations .tsv file, default `NULL`.
#' @param output_dir  string, path to network folder; default `'./Networks_output/'`.
#'
#' @param desired_phenotypes SIGNOR phenotypes vector to infer the activity and include in the model; default:  `NULL`.
#' @param pheno_distances_table character, any character except NULL to internally compute it; default:  `NULL`.
#' @param save_all_files Boolean, if TRUE it will save 13 files per patient, FALSE, just 4; default `FALSE`.
#'
#' @param PKN_options  check *initialize_net_default_params* documentation, in the `$PKN_options` part.
#' @param naive_options check *initialize_net_default_params* documentation, in the `$naive_options` part.
#' @param carnival_options check *initialize_net_default_params* documentation, in the `$carnival_options` part.
#' @param phenoscore_options check *initialize_net_default_params* documentation, in the `$phenoscore_options` part.
#' @param format_options check *initialize_net_default_params* documentation, in the `$format_options` part.
#'
#' @return  No return value. The final results are saved in *RDS* and *sif* files in the output directory.
#' @export
#'
#' @examples
#'
#' PKN_options = list()
#' naive_options = list()
#' carnival_options = list()
#' phenoscore_options = list()
#' format_options = list()
#'
#'create_cohort_networks(trans_dir = './Trans_test/',
#'                       prot_dir = './Prot_test/',
#'                       phospho_dir = './Phos_test/',
#'                       act_dir = './Act_test/',
#'                       mut_file = 'mutations.tsv',
#'                       output_dir = './Networks_output/',
#'                       desired_phenotypes = c('APOPTOSIS', 'PROLIFERATION'),
#'                       pheno_distances_table = TRUE,
#'                       PKN_options = list(direct = FALSE),
#'                       naive_options = list(layers = 2, max_length = c(1,4)),
#'                       carnival_options = list(solver = 'cplex', carnival_type = 'inverse'),
#'                       format_options = list(optimize_on_phenotypes = TRUE, vis_cytoscape = FALSE))
#

create_cohort_networks <- function(
    trans_dir = NULL,
    prot_dir = NULL,
    phospho_dir = NULL,
    act_dir = NULL,
    desired_phenotypes = c(),
    pheno_distances_table = NULL,
    mut_file = NULL, #csv file with mutations...
    output_dir = './Networks_output/',
    save_all_files = FALSE,
    PKN_options = list(),
    naive_options = list(),
    carnival_options = list(),
    phenoscore_options = list(),
    format_options = list()
){

  # CHECK ON INPUTS
  if(length(desired_phenotypes) == 0){
    stop('Please select a vector of phenotypes to connect in your networks')
  }

  # Retrieve omics file lists from specified directories
  prot_files <- if (!is.null(prot_dir)) {
    list.files(path = prot_dir, pattern = "^Prot_Patient_.*\\.tsv$", full.names = TRUE)
  } else character()

  phospho_files <- if (!is.null(phospho_dir)) {
    list.files(path = phospho_dir, pattern = "^Phospho_Patient_.*\\.tsv$", full.names = TRUE)
  } else character()

  trans_files <- if (!is.null(trans_dir)) {
    list.files(path = trans_dir, pattern = "^Transc_Patient_.*\\.tsv$", full.names = TRUE)
  } else character()


  # Retrieve activities file lists from specified directories
  act_files <- if (!is.null(act_dir)) {
    
    constraint_files <- list.files(path = act_dir, pattern = "^Activity_constraints_Patient_.*\\.tsv$", full.names = TRUE)
    
    if (length(constraint_files) > 0) {
      constraint_files 
    } else {
      list.files(path = act_dir, pattern = "^Activity_Patient_.*\\.tsv$", full.names = TRUE) 
    }
  } else character()
  
  # Read mutations file
  sources_df <- if(!is.null(mut_file)){
    readr::read_tsv(mut_file)
  }else character()

  trans_ids <- sapply(trans_files, extract_patient_id)
  prot_ids <- sapply(prot_files, extract_patient_id)
  phospho_ids <- sapply(phospho_files, extract_patient_id)
  act_ids <- sapply(act_files, extract_patient_id)

  all_patient_ids <- unique(c(prot_ids, trans_ids, phospho_ids))

  patient_count <- 0
  total_patients <- length(all_patient_ids)

  successful_pat <- 0

  for (patient_id in all_patient_ids){
    patient_count <- patient_count + 1

    trans_file <- trans_files[which(trans_ids == patient_id)[1]]
    prot_file <- prot_files[which(prot_ids == patient_id)[1]]
    phospho_file <- phospho_files[which(phospho_ids == patient_id)[1]]
    act_file <- act_files[which(act_ids == patient_id)[1]]
    mut_row <- sources_df[sources_df$Patient_ID == patient_id, ]

    message(paste("Processing patient", patient_count, "of", total_patients, ":", patient_id))

    if (is.na(prot_file) && is.na(trans_file) && is.na(phospho_file)) {
      message(paste("Skipping patient", patient_id, "- no omics data available."))
      next
    }else if(is.na(act_file)){
      message(paste("Skipping patient", patient_id, "- no activities data available."))
      next
    }else if(nrow(mut_row) == 0){
      message(paste("Skipping patient", patient_id, "- no mutations data available."))
      next
    }

    successful_pat <- successful_pat + 1

    # Read files
    if(!is.na(trans_file)){
      Trans_P <- readr::read_tsv(trans_file)
    }else{
      Trans_P <- NULL
    }

    if(!is.na(prot_file)){
      Prot_P <- readr::read_tsv(prot_file)
    }else{
      Prot_P <- NULL
    }

    if(!is.na(phospho_file)){
      Phospho_P <- readr::read_tsv(phospho_file)
    }else{
      Phospho_P <- NULL
    }


    if (!is.na(act_file)) {
      act_P <- readr::read_tsv(act_file)
    } else {
      act_P <- NULL
    }


    # For the first successfull patient, we preprocess ProxPath
    if( successful_pat == 1 ){
      if(is.null(pheno_distances_table)){
        warning('No pheno_distances_table provided, using built-in one')
      }else{
        pheno_distances_table <- proxpath_preprocessing(proteomics = Prot_P,
                                                        phosphoproteomics = Phospho_P)

      }
    }

    tryCatch({

      create_network(patient_id = patient_id,
                     sources = mut_row,
                     activities = act_P,
                     transcriptomics = Trans_P,
                     proteomics = Prot_P,
                     phosphoproteomics = Phospho_P,
                     desired_phenotypes = desired_phenotypes,
                     pheno_distances_table = pheno_distances_table,
                     output_dir = output_dir,
                     save_all_files = save_all_files,
                     PKN_options = PKN_options,
                     naive_options = naive_options,
                     carnival_options = carnival_options,
                     phenoscore_options = phenoscore_options,
                     format_options = format_options
                     )

      message(paste("Successfully created network for patient", patient_id, "-", patient_count, "of", total_patients))

    }, error = function(e){
      message(paste("Error processing patient", patient_id, "-", patient_count, "of", total_patients, ":", e$message))
    })
  }
}
