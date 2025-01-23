
#' format_patient_network
#'
#' This function takes as input proteins-to-phenotypes network and can:
#' - (optionally) visualize it in Cytoscape (if `format_options$vis_cytoscape = TRUE`);
#' - (optionally) optimize the network on phenotypes activity (if `format_options$optimize_on_phenotypes = TRUE`);
#' - (optionally) create functional phenotypes (if `format_options$circuits_params$k != -1`);
#'
#' *Result output*
#'  - a *list* containing a barplot of phenotypes activity, table of phenotypes regulators and phenotypic activity, and list of protein-to-phenotypes network (igraph object, nodes, and edges);
#'  - a *igraph object* of the functional circuit connecting user-defined `sources` to `desired_phenotypes`;
#'  - creates internally an *RDS file* containing the functional circuit;
#'  - if `format_options$optimize_on_phenotypes = TRUE`, creates internally an *RDS file* and *SIF file* containing the optimized network on phentoypes.
#'
#' @param patient_id string, patient_ID
#' @param phenoscore_output list of four entities returned by *infer_and_link_phenotypes*
#' @param sources dataframe of perturbed nodes (e.g., patient mutations), starting point of the model
#' @param phosphoproteomics a dataframe of phosphoproteomics
#' @param proteomics a dataframe of proteomics
#' @param desired_phenotypes SIGNOR phenotypes vector to infer the activity and include in the model; default:  `NULL`.
#' @param network_params  check *initialize_net_default_params* documentation, in the `$format_options` part.
#' @return A list of two elements:
#'  - `patient_opt_pheno_network`: a list of four elements with same structure as *infer_and_link_phenotypes* output;  if `format_options$optimize_on_phenotypes = TRUE`, `sp_object_phenotypes` is optimized on phenotypes activity;
#'  - `circuit`: igraph object representing the functional circuit of `k` length from sources to targets
#'
#' @examples
#'
#' output_dir <- './Networks_output/'
#' network_params <- initialize_net_default_params(output_dir)
#' network_params$format_options <- list(optimize_on_phenotypes = TRUE,
#'                                       circuit_params = list(k=-1),
#'                                       vis_cytoscape = TRUE)
#' phenoscore_output = readRDS(paste0(output_dir, './Pheno_', patient_id, '_object.RDS'))
#'
#' format_patient_network(patient_id = 'Patient1',
#'                        phenoscore_output = phenoscore_output,
#'                        sources = sources, # read mutation file
#'                        desired_phenotypes = c('APOPTOSIS', 'PROLIFERATION'),
#'                        proteomics, # read proteomics
#'                        phosphoproteomics, #read phosphoproteomics file
#'                        network_params = network_params)

format_patient_network <- function(patient_id,
                                   phenoscore_output,
                                   phosphoproteomics,
                                   proteomics,
                                   sources,
                                   desired_phenotypes = NULL,
                                   network_params){

  if(is.null(desired_phenotypes)){
    desired_phenotypes <- SignalingProfiler::all_phenotypes
  }

  # Get default parameters
  default_format_options <-  list(optimize_on_phenotypes = TRUE,
                                 circuits_params = list(k = 10,
                                                        start_to_top = TRUE),
                                 vis_cytoscape = FALSE)

  network_params$format_options <-  modifyList(default_format_options,  network_params$format_options)

  # If user wants to optimized over the phenotypes' activity
  if(network_params$format_options$optimize_on_phenotypes){

    message('Optimizing protein-phenotype network on phenotypes activity...')
    default_carnival_params <- SignalingProfiler::default_CARNIVAL_options(solver = network_params$carnival_options$solver)
    network_params$carnival_options$carnival_params <- modifyList(default_carnival_params, network_params$carnival_options$carnival_params)

    phenoscore_output$sp_object_phenotypes <- SignalingProfiler::format_for_visualization(sp_object = phenoscore_output$sp_object_phenotypes)

    phenoscore_output <- SignalingProfiler::optimize_pheno_network(sp_object = phenoscore_output,
                                                                   organism = 'human',
                                                                   phospho_df = phosphoproteomics,
                                                                   carnival_options = network_params$carnival_options$carnival_params,
                                                                   files = TRUE,
                                                                   direct = network_params$PKN_options$direct,
                                                                   with_atlas = network_params$PKN_options$with_atlas,
                                                                   path_sif = paste0(network_params$pheno_options$pheno_path, '_opt_network.sif'),
                                                                   path_rds =  paste0(network_params$pheno_options$pheno_path, '_opt_network.RDS'))

    write_rds(phenoscore_output, paste0(network_params$pheno_options$pheno_path,'_object_opt.RDS'))

  }

  if(network_params$format_options$circuits_params$k == -1){

    if(network_params$format_options$vis_cytoscape){
      # Visualize according to PatientProfiler visualization style
      data_path <- system.file("extdata", "SP_pheno_layout.xml", package = "SignalingProfiler")
      RCy3::importVisualStyles(filename = data_path)

      phenoscore_output$sp_object_phenotypes <- SignalingProfiler::format_for_visualization(sp_object = phenoscore_output$sp_object_phenotypes)
      # Create the network in Cytoscape
      RCy3::createNetworkFromIgraph(igraph=  phenoscore_output$sp_object_phenotypes$igraph_network,
                                    title = patient_id,
                                    collection = 'PatientProfiler_network')
      RCy3::setVisualStyle('SP_pheno_layout')
    }

    return(list(patient_opt_pheno_network = phenoscore_output,
                circuit = NULL))

  } else {

    # Selected souces for circuit distinguishing between inverse and vanilla
    if(network_params$carnival_options$carnival_type == 'inverse' |  is.null(sources)){
      # If inverse CARNIVAL optimization, remove Perturbation node
      graph <- phenoscore_output$sp_object_phenotypes$igraph_network
      pert <- V(graph)[which(V(graph)$name == "Perturbation")]
      graph <- delete_vertices(graph, pert)

      # Select as sources nodes with indegree = 0, removing Perturbation
      in_degree <- degree(graph, mode = "in")
      nodes_with_zero_indegree <- which(in_degree == 0)
      nodes_with_zero_indegree <- nodes_with_zero_indegree[names(nodes_with_zero_indegree) != "Perturbation"]
      start_nodes <- unname(names(nodes_with_zero_indegree))
    }else{
      # Intesect user-selected sources with the connected ones
      start_nodes <- names(from_df_to_list(df_cols = sources))
      start_nodes <- intersect(start_nodes, V(phenoscore_output$sp_object_phenotypes$igraph_network)$name)
    }

    if(length(start_nodes) == 0){
      stop('No connected selected nodes for functional circuits: please change the node list!')
    }

    # Intesect user-selected phenotypes with the connected ones
    phenotypes_clean <- intersect(desired_phenotypes, V(phenoscore_output$sp_object_phenotypes$igraph_network)$name)

    if(length(phenotypes_clean) == 0){
      stop('No connected selected phenotypes for functional circuits: please change the node list!')
    }

    # Create circuit
    message('Building functional circuits...')
    circuit <- SignalingProfiler::pheno_to_start_circuit(SP_object = phenoscore_output$sp_object_phenotype,
                                                         start_nodes = start_nodes,
                                                         phenotypes = phenotypes_clean,
                                                         k = network_params$format_options$circuits_params$k,
                                                         start_to_top = network_params$format_options$circuits_params$start_to_top)

    saveRDS(circuit, file =  paste0(network_params$phenoscore_options$pheno_path, '_circuit.RDS'))

    if(network_params$format_options$vis_cytoscape){
      # Visualize according to PatientProfiler visualization style
      data_path <- system.file("extdata", "SP_pheno_layout.xml", package = "SignalingProfiler")
      RCy3::importVisualStyles(filename = data_path)

      RCy3::createNetworkFromIgraph(igraph= circuit,
                                    title = patient_id,
                                    collection = 'PatientProfiler_circuit')
      RCy3::setVisualStyle('SP_pheno_layout')
    }

    # Return both optimized and circuit
    return(list(patient_opt_pheno_network = phenoscore_output,
                circuit = circuit))

  }
}

