
#' initialize_net_default_params
#'
#' set default PatientProfiler network building parameters
#'
#' @param output_dir user-provided network building directory
#'
#' @return it returns a list of 5 lists for each internal function:
#' - `$PKN_options`: parameters list for SignalingProfiler PKN choice of 7 elements:
#'    - `preprocess`: Boolean, whether to remove interactions involving entities not detected in at least one omic data; default `TRUE`.
#'    - `direct`: Boolean, whether keep only direct interactions; default `TRUE`.
#'    - `organism`: string, available are 'human' or 'mouse'; default `human`.
#'    - `with_atlas`: Boolean, whether integratins Ser/Thr and Tyr atlas inferred interactions; default `FALSE`.
#'    - `custom`: Boolean, whether the user wants to provide a custom PKN;  default `FALSE`.
#'    - `custom_path`: string, path to the custom PKN;  default `FALSE`.
#'
#' - `$naive_options`: parameters list for SignalingProfiler naive network layers of 4 elements:
#'    - `layers`: integer, specifiying the number of layers > 0 and < 4; default  `2`
#'    - `max_length` a vector of int or an int, length == layers; default  `c(1,4)`
#'    - `connect_all` Boolean, whether connecting intermediate nodes of shortest path; default  `TRUE`
#'    - `naive_path` path of the folder where naive network files are saved; default `./Network_output/Naive_`
#'
#' - `$carnival_options`: parameters list for SignalingProfiler CARNIVAL implementations of 4 elements:
#'     - `solver`: string, specifiying the ILP solver for CARNIVAL; default  `cplex`
#'     - `carnival_type`: string, specifying the CARNIVAL implementation among 'inverse', 'vanilla_one_shot', 'vanilla_two_shots', 'vanilla_three_shots';
#'     - `opt_path`: path to the folder were optimized networks files are saved; default `./Network_output/Opt_`
#'     - `carnival_params`: a list of CARNIVAL parameters, set in *SignalingProfiler::default_CARNIVAL_options(solver)*
#'
#' - `$phenoscore_options`: parameters list for SignalingProfiler PhenoScore computation of 3 elements:
#'    - `pheno_path`:  path to the folder were optimized networks files are saved; default `./Network_output/Pheno_`
#'    - `pheno_distance_table`: output of *proxpath_preprocessing* function; default `NULL`.
#'    -  `phenoscore_params`: list of 9 elements for PhenoScore computation:
#'        - `path_length`: max length between proteins and phenotypes; default `3`.
#'        - `stat`: ProxPath z-score statistics for significantly close phenotypes; default `mean`.
#'        - `zscore_threshold`: ProxPath z-score threshold for significantly close phenotypes; default `-1.96`.
#'        - `nrandom`: ProxPath number of randomization of input proteins length;  default `1000`.
#'        - `pvalue_threshold`: ProxPath t-test p-value adjusted randomization threshold;  default `0.05`.
#'        - `remove_cascade`: in linear cascade of proteins regulating a phenotype, consider only the most upstream one in PhenoScore computation; default `TRUE`.
#'        - `node_idx`: whether weighting protein activity contribution to phenotypes' regulation based on number of paths; default  `FALSE`.
#'        - `use_carnival_activity`: whether consider all proteins in the network inferred by carnival (TRUE) or only the experimentally-derived ones (FALSE); default `FALSE`.
#'        - `create_pheno_network`: whether connect model proteins to phenotypes; default `TRUE`.
#'
#' - `$format_options`: parameters list for manipulation of proteins-to-phenotypes models of 3 elements:
#'  - `optimize_on_phenotypes`: Boolean, whether keep only proteins-to-phenotypes edges coherent with the inferred phenotypes activity; default, `FALSE`.
#'  - `circuits_params`: a list of parameters for functional circuits creation:
#'       - `k`: max path length between proteins an phenotypes in the model; if `circuits_params$k=-1` no circuits will be created;
#'       - `start_to_top`: Boolean, whether starting from phenotypes looking for upstream regulators or from sources; default, `TRUE`.
#'  - `vis_cytoscape`: Boolean, whether using *RCy3* package to visualize networks using SignalingProfiler style (if TRUE, ensure to have Cytoscape open); default, `FALSE`.
#'
#' @export
#'

initialize_net_default_params <- function(output_dir) {
  list(
    PKN_options = list(preprocess = TRUE, organism = 'human', direct = TRUE, with_atlas = FALSE, custom = FALSE, custom_path = FALSE),
    naive_options = list(layers = 2, max_length = c(1, 4), connect_all = TRUE, naive_path = paste0(output_dir, '/Naive_')),
    carnival_options = list(solver = 'cplex', carnival_type = 'inverse', opt_path = paste0(output_dir, '/Opt_'), carnival_params = list()),
    phenoscore_options = list(
      pheno_path = paste0(output_dir, '/Pheno_'),
      pheno_distance_table = NULL,
      phenoscore_params = list(
        path_length = 3, stat = 'mean', zscore_threshold = -1.96,
        nrandom = 1000, pvalue_threshold = 0.05,
        remove_cascade = TRUE, node_idx = FALSE, use_carnival_activity = FALSE,
        create_pheno_network = TRUE
      )
    ),
    format_options = list(
      optimize_on_phenotypes = TRUE,
      circuits_params = list(k = 10, start_to_top = TRUE),
      vis_cytoscape = FALSE
    )
  )
}


#' create_network
#'
#' This wrapper to SignalingProfiler network creation step sequentially:
#' 1) Selects a PKN using `get_PKN` function of PatientProfiler
#' 2) Creates a naive network using `create_naive_network` function of PatientProfiler
#' 3) Optimizes the naive network over CARNIVAL activity using `optimize_network_with_carnival` function of PatientProfiler
#' 4) Infers the activity of phenotypes from model proteins and connect them in a proteins-to-phenotypes network using `infer_and_link_phenotypes` function of PatientProfiler
#' 5) Manipulates the proteins-to-phenotypes model for visualization and functional circuits creation using `format_patient_networks` function of PatientProfiler
#'
#' Final results for each patient will be saved as  *RDS* and/or *sif* files in `output_dir` folder, as follows:
#' - Naive network files:
#'      -  `Naive_{patient_id}.RDS` and `Naive_{patient_id}.sif`
#'
#' - Optimized network files:
#'      - `Opt_{patient_id}.RDS` and `Opt_{patient_id}.sif`
#'      - (optional, if `phosphoproteomics` is provided) `Opt_{patient_id}_val.RDS` and `Opt_{patient_id}_val.sif`, phosphoproteomics it is mapped on optimized networks edges;
#'
#' - Phenoscore computation:
#'      - `Pheno_{patient_id}_object.RDS`: output list of *infer_and_link_phenotypes* function.
#'      - (optional, if `phenoscore_options$create_pheno_network == TRUE`): `Pheno_{patient_id}.RDS` and `Pheno_{patient_id}.sif`: proteins-to-phenotypes network.
#'
#' - Formatting patient-specific network:
#'      - (optional, if `format_options$optimize_on_phenotypes == TRUE`):
#'         - `Pheno_{patient_id}_object_opt.RDS`: output list of *infer_and_link_phenotypes* function, overridden in `sp_object_phenotypes` with optimized proteins-to-phenotypes network on phenotypes activities.
#'         - `Pheno_{patient_id}_opt_network.RDS` and `Pheno_{patient_id}_opt_network.sif`: optimized proteins-to-phenotypes network on phenotypes activities.
#'      - (optional, if `format_options$circuits_params$k != -1`): `Pheno_{patient_id}_circuit.RDS`: functional circuit in igraph format.
#'
#' @param patient_id string, patient_ID
#' @param sources dataframe of perturbed nodes (e.g., patient mutations), starting point of the model
#' @param activities dataframe of proteins inferred activity to connect to the starting point
#'
#' @param desired_phenotypes SIGNOR phenotypes vector to infer the activity and include in the model; default:  `NULL`.
#'
#' @param transcriptomics  dataframe of transcriptomics; default: `NULL`.
#' @param proteomics dataframe of proteomics; default: `NULL`.
#' @param phosphoproteomics  dataframe of phosphoproteomics; default: `NULL`.

#' @param pheno_distances_table dataframe of ProxPath protein-to-phenotypes distances; default:  `NULL`.
#' @param output_dir string, path to network folder; default `'./Networks_output/'`.
#' @param save_all_files Boolean, if TRUE it will save 13 files per patient, FALSE, just 4; default `FALSE`.

#' @param PKN_options  check *initialize_net_default_params* documentation, in the `$PKN_options` part.
#' @param naive_options check *initialize_net_default_params* documentation, in the `$naive_options` part.
#' @param carnival_options check *initialize_net_default_params* documentation, in the `$carnival_options` part.
#' @param phenoscore_options check *initialize_net_default_params* documentation, in the `$phenoscore_options` part.
#' @param format_options check *initialize_net_default_params* documentation, in the `$format_options` part.

#' @return No return.
#' @export
#'
#' @examples
#'
#' patient_id <- 'CPT000814'
#'
#' mutations_df <- read_csv('./mutations.csv')
#' sources <- mutations_df[mutations_df$Patient_ID == patient_id, ]
#' activities <- read_csv(paste0('./Activity_Patient_', patient_id))
#'
#' proteomics <- read_csv(paste0('./Prot_Patient_', patient_id))
#' transcriptomics <- read_csv(paste0('./Trans_Patient_', patient_id))
#' phosphoproteomics <- read_csv(paste0('./Phospho_Patient_', patient_id))
#'
#' output_dir <- './Networks_output/'
#' desired_phenotypes <- c('APOPTOSIS', 'PROLIFERATION')
#' pheno_distances_table <- proxpath_preprocessing(proteomics = proteomics, phosphoproteomics = phosphoproteomics)
#'
#' network_params <- initialize_net_default_params(output_dir)

#' create_network(patient_id = 'CPT000814',
#'               sources = sources,
#'               activities = activities,
#'                transcriptomics = transcriptomics,
#'                proteomics = proteomics,
#'                phosphoproteomics = phosphoproteomics,
#'                desired_phenotypes = c('APOPTOSIS', 'PROLIFERATION'),
#'                pheno_distances_table = pheno_distances_table,
#'                output_dir = output_dir,
#'                PKN_options = list(direct = FALSE),
#'                naive_options = list(layers = 3, max_length = c(2,1,3)),
#'                carnival_options = list(carnival_type = 'vanilla_one_shot'),
#'                phenoscore_options = list(),
#'                format_options = list(optimize_on_phenotypes = TRUE,
#'                                      circuits_params = list(k = -1),
#'                                      vis_cytoscape = TRUE))
#'
create_network <- function(patient_id,
                           sources,
                           activities,
                           transcriptomics = NULL,
                           proteomics = NULL,
                           phosphoproteomics = NULL,
                           desired_phenotypes = NULL,
                           pheno_distances_table = NULL,
                           save_all_files = FALSE,
                           output_dir = './Networks_output/',
                           PKN_options = list(),
                           naive_options = list(),
                           carnival_options = list(),
                           phenoscore_options = list(),
                           format_options = list())
{

  # Set output directory
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }

  # Aggregate user_params
  user_params = list(PKN_options = PKN_options,
                     naive_options = naive_options,
                     carnival_options = carnival_options,
                     phenoscore_options = phenoscore_options,
                     format_options = format_options)

  # Initialize and modify default parameters
  default_params <- initialize_net_default_params(output_dir)
  network_params <- modifyList(default_params, user_params)

  # Get Prior Knowledge Network
  omics_list <- list(transcriptomics, proteomics, phosphoproteomics)
  PKN <- get_PKN(PKN_params = network_params$PKN_options, omics_list = omics_list)

  # Generate Naive Network
  network_params$naive_options$naive_path <- paste0(network_params$naive_options$naive_path, patient_id)
  naive_network <- create_naive_network(PKN = PKN,
                                        sources = sources,
                                        activities = activities,
                                        naive_network_parameters = network_params$naive_options)
  message('Naive network building done!')

  # Optimize Network with CARNIVAL
  network_params$carnival_options$opt_path <- paste0(network_params$carnival_options$opt_path, patient_id)
  carnival_output <- optimize_network_with_carnival(sources = sources,
                                                    activities = activities,
                                                    naive_network = naive_network,
                                                    phosphoproteomics = phosphoproteomics,
                                                    save_all_files = save_all_files,
                                                    network_params = network_params)
  message('CARNIVAL optimization done!')

  # PhenoScore analysis
  message('Running phenotypic inference...')
  network_params$phenoscore_options$pheno_path <- paste0(network_params$phenoscore_options$pheno_path, patient_id)
  phenoscore_output <- infer_and_link_phenotypes(carnival_output = carnival_output,
                                                desired_phenotypes = desired_phenotypes,
                                                proteomics = proteomics,
                                                phosphoproteomics = phosphoproteomics,
                                                save_all_files = save_all_files,
                                                network_params = network_params)
  message('Phenotypic inference done!')

  phenoscore_network <- format_patient_network(patient_id,
                                               phenoscore_output = phenoscore_output,
                                               phosphoproteomics = phosphoproteomics,
                                               desired_phenotypes = desired_phenotypes,
                                               sources = sources,
                                               save_all_files = save_all_files,
                                               network_params = network_params)

}



