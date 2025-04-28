#' infer_and_link_phenotypes
#'
#' Wrapper to phenoscore computation algorithm in SignalingProfiler
#'
#' *Return output*
#' - a *list* containing a barplot of phenotypes activity, table of phenotypes regulators and phenotypic activity;
#' - optionally, if `network_params$phenoscore_options$create_pheno_network == TRUE`:
#'    - add `sp_object_phenotypes`, a list containing the proteins to phenotypes network in the form of igraph object, nodes table, and edges table
#'    - creates internally an  *RDS file* and a *SIF file* of the proteins-phenotypes network at `output_dir` directory
#'
#' @param carnival_output output of `optimize_network_with_carnival`, a list of `igraph_network`, `nodes_df`, and `edges_df`
#' @param desired_phenotypes SIGNOR phenotypes vector to infer the activity and include in the model; default:  `NULL`.
#' @param proteomics dataframe of proteomics; default: `NULL`.
#' @param phosphoproteomics  dataframe of phosphoproteomics; default: `NULL`.
#' @param network_params  check *initialize_net_default_params* documentation, in the `$phenoscore_options` part.
#' @param pheno_distances_table dataframe of ProxPath protein-to-phenotypes distances; default:  `NULL`.
#' @param save_all_files Boolean, if TRUE it will save 13 files per patient, FALSE, just 4; default `FALSE`.
#'
#' @return A list of four elements:
#'  - `barplot`: a ggplot2 object representing phenotype modulation inferred from model proteins;
#'  - `table_regulators`: a dataframe reporting for each inferred phenotype the regulators in the model;
#'  - `table_phenotypes`: a dataframe reporting the activity of each inferred phenotype
#'  - `sp_object_phenotypes`: the proteins-to-phenotypes network as a list of:
#'      - `igraph_object`: proteins-to-phenotypes network in igraph format
#'      - `nodes_df`: dataframe of the nodes in the proteins-to-phenotypes network; attributes are described in *optimize_network_with_carnival* documentation
#'      - `edges_df`: dataframe of the edges in the proteins-to-phenotypes network; attributes are described in *optimize_network_with_carnival* documentation
#'
#' @export
#'
#' @examples
#' network_params <- initialize_net_default_params('./Network_output/')
#' network_params$carnival_options <- list(solver = 'cplex', carnival_type = 'inverse')
#' pheno_distance_table <- proxpath_preprocessing(proteomics, phosphoproteomics)
#' network_params$phenoscore_options <- list(pheno_distance_table = pheno_distance_table)
#'
#' omics_list <- list(transcriptomics, proteomics, phosphoproteomics)
#' PKN = get_PKN(omics_list = omics_list)
#' naive_network = create_naive_network(PKN, sources, activities, naive_network_parameters = list(layers = 2, max_length = c(1,4)))
#' carnival_network = optimize_network_with_carnival(sources, activities, phosphoproteomics, network_params = network_params)
#' phenoscore_output = infer_and_link_phenotypes(carnival_output = carnival_network,
#'                                               desired_phenotypes = c('APOPTOSIS', 'PROLIFERATION'),
#'                                               pheno_distance_table = pheno_distance_table,
#'                                               proteomics = proteomics,
#'                                               phosphoproteomics = phosphoproteomics,
#'                                               network_params = network_params)
#'

infer_and_link_phenotypes <- function(carnival_output,
                                      desired_phenotypes = NULL,
                                      pheno_distances_table = NULL,
                                      proteomics = NULL,
                                      phosphoproteomics = NULL,
                                      save_all_files = FALSE,
                                      network_params = list()){

  if(is.null(pheno_distances_table)){
    message('Using SignalingProfiler built-in distance table')
    pheno_distances_table <- SignalingProfiler::access_remote_file(file = 'phenoscore_distances_table.tsv', dir = 'PKN')
  }

  phenoscore_output <- SignalingProfiler::phenoscore_computation(proteins_df = carnival_output$nodes_df,
                                              desired_phenotypes = desired_phenotypes,
                                              pheno_distances_table = pheno_distances_table,
                                              sp_graph = carnival_output$igraph_network,
                                              # closeness of proteins to phenotypes
                                              path_length = network_params$phenoscore_options$phenoscore_params$path_length,
                                              stat = network_params$phenoscore_options$phenoscore_params$stat,
                                              zscore_threshold = network_params$phenoscore_options$phenoscore_params$zscore_threshold,
                                              # exclude random phenotypes
                                              n_random = network_params$phenoscore_options$phenoscore_params$nrandom,
                                              pvalue_threshold =  network_params$phenoscore_options$phenoscore_params$pvalue_threshold,
                                              # optimized network  specificity
                                              remove_cascade =  network_params$phenoscore_options$phenoscore_params$remove_cascade,
                                              node_idx = network_params$phenoscore_options$phenoscore_params$node_idx,
                                              use_carnival_activity = network_params$phenoscore_options$phenoscore_params$use_carnival_activity,
                                              create_pheno_network = network_params$phenoscore_options$phenoscore_params$create_pheno_network
                                              )

  if( network_params$phenoscore_options$phenoscore_params$create_pheno_network){
    igraphToSif(phenoscore_output$sp_object_phenotypes$igraph_network, #to check
                paste0(network_params$phenoscore_options$pheno_path, '.sif'), "sign")

    pp_output <- list(sp_object_phenotypes = phenoscore_output$sp_object_phenotypes)
    saveRDS(pp_output,
            paste0(network_params$phenoscore_options$pheno_path, '.RDS'))
  }

  if(save_all_files){
    saveRDS(phenoscore_output,  paste0(network_params$phenoscore_options$pheno_path, '_object.RDS'))
  }

  return(phenoscore_output)
}

#' proxpath_preprocessing
#'
#' This function calls a Python script in SignalingProfiler that downloads proteins-to-phenotypes interactions in SIGNOR
#' and preprocess it according to (phospho)proteomics data, removing interactions not quantified experimentally
#' If no omics are provided, the function return the built-in object in SignalingProfiler
#'
#' @param proteomics dataframe of proteomics; default: `NULL`.
#' @param phosphoproteomics dataframe of phosphoproteomics; default: `NULL`.
#'
#' @return `pheno_distances_table` a table reporting the distance between each protein-phenotype pair in SIGNOR
#' @export
#'
#' @examples
#' pheno_distances_table <- proxpath_preprocessing(proteomics = NULL, phosphoproteomics = phosphoproteomics)
#'
proxpath_preprocessing <- function(proteomics, phosphoproteomics){
  if(!is.null(proteomics) & !is.null(phosphoproteomics)){
    message('Preprocessing ProxPath distance table')
    pheno_distances_table <- SignalingProfiler::phenoscore_network_preprocessing(proteomics = proteomics,
                                                              phospho = phosphoproteomics)
  }else if(!is.null(proteomics) & is.null(phosphoproteomics)){
    message('Preprocessing ProxPath distance table')
    pheno_distances_table <- SignalingProfiler::phenoscore_network_preprocessing(proteomics = proteomics,
                                                              phospho = proteomics)
  }else if(is.null(proteomics) & !is.null(phosphoproteomics)){
    message('Preprocessing ProxPath distance table')
    pheno_distances_table <- SignalingProfiler::phenoscore_network_preprocessing(proteomics = phosphoproteomics,
                                                              phospho = phosphoproteomics)
  }else{
    pheno_distances_table = get(data("phenoscore_distances_table", package = 'SignalingProfiler'))
  }

  return(pheno_distances_table)
}






