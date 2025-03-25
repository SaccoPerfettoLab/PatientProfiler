#' optimize_network_with_carnival
#'
#' Wrapper function to the four implementations of vanilla and inverse CARNIVAL in SignalingProfiler
#' to retrieve sign-coherent interactions with inferred protein activities
#'
#'  *Results output*
#' An optimized network with sign-coherent edges as:
#' - a *list* of igraph_network, nodes and edges table;
#' - creates internally an *RDS file* and a *SIF file* of the optimized network at `output_dir` directory
#'
#' @param sources dataframe of perturbed nodes (e.g., patient mutations), starting point of the model
#' @param activities dataframe of proteins inferred activity to connect to the starting point
#' @param naive_network igraph object representing the naive network
#' @param network_params check *initialize_net_default_params* documentation, in the `$carnival_options` part.
#' @param phosphoproteomics dataframe of phosphoproteomics; default: `NULL` (no phosphoproteomics mapping on network).
#' @param save_all_files Boolean, if TRUE it will save 13 files per patient, FALSE, just 4; default `FALSE`.
#'
#' @return A list of three elements:
#'  - `igraph_network`: optimized network in igraph format
#'  - `nodes_df`: dataframe of the nodes in the optimized network with the following attributes:
#'    - `gene_name`: HUGO gene symbol of each node
#'    - `UNIPROT`: UNIPROT ID of each node
#'    - `carnival_activity`: node activity predicted by CARNIVAL in the optimization process (usually equal to final_score)
#'    - `final_score`: node activity predicted by extract_protein_activity from experimental data
#'    - `method`: algorithm of the activity; can be CARNIVAL, VIPER (footprint-based analysis), or PhosphoScore
#'    - `discordant`: Boolean, whether `carnival_activity` and `final_score` have opposite signs
#'    - `mf`: molecular function of each node (rec, kin, phos, other, tf)
#'
#'  - `edges_df`: dataframe of the edges in the optimized network with the following attributes:
#'    - `soource`: HUGO gene symbol of starting node
#'    - `target`: HUGO gene symbol of target node
#'    - `sign`: sign of the interactions (1, -1)
#'    - `carnival_weight`: interaction confidence in CARNIVAL optimization (0-100)
#'    - `direct`: Boolean, whether the interaction is direct
#'    - `mechanism`: mechanism associated to the interaction in SIGNOR
#'    - `aminoacid`: if `phosphoproteomics` is provided, represents putative phosphosites modulated in experimental data  modified by the interaction
#'    - `is_quantified`: if `phosphoproteomics` is provided, when TRUE the mapped phosphosite modulation is coherent with interaction sign, but not necessarily significant
#'    - `is_significant`:  if `phosphoproteomics` is provided, when TRUE the mapped phosphosite modulation is coherent with interaction sign and significant
#'    - `FC`: if `phosphoproteomics` is provided, reports the phosphosite's modulation in experimental data
#'
#' @export
#'
#' @examples
#' network_params <- initialize_net_default_params('./Network_output/')
#' network_params$carnival_options <- list(solver = 'cplex', carnival_type = 'inverse')
#'
#' omics_list <- list(transcriptomics, proteomics, phosphoproteomics)
#' PKN = get_PKN(omics_list = omics_list)
#' naive_network = create_naive_network(PKN, sources, activities, naive_network_parameters = list(layers = 2, max_length = c(1,4)))
#' carnival_network = optimize_network_with_carnival(sources, activities, phosphoproteomics, network_params = network_params)

optimize_network_with_carnival <- function(sources,
                                           activities,
                                           naive_network,
                                           phosphoproteomics = NULL,
                                           save_all_files = FALSE,
                                           network_params = list()){

  # Set default technical CARNIVAL parameters
  default_carnival_params <- SignalingProfiler::default_CARNIVAL_options(solver = network_params$carnival_options$solver)
  carnival_params <- modifyList(default_carnival_params, network_params$carnival_options$carnival_params)

  message(paste0('Running CARNIVAL ', network_params$carnival_options$carnival_type, ' optimization...'))

  # No receptor list if inverse CARNIVAL run
  if(network_params$carnival_options$carnival_type == 'inverse'){
    recept_list <- NULL
  }else{
    recept_list <-  from_df_to_list(sources)
  }

  # Prepare CARNIVAL input
  carnival_input <- SignalingProfiler::prepare_carnival_input(naive_network,
                                                              organism = 'human',
                                                              prediction_output = activities,
                                                              recept_list = recept_list)

  # Read naive network
  naive_table <- readr::read_tsv(paste0(network_params$naive_options$naive_path, '.sif'),
                                 col_names = c('source', 'interaction', 'target'))

  # One shot for both vanilla one-shot and inverse
  if(network_params$carnival_options$carnival_type %in% c('vanilla_one_shot', 'inverse')){

    if(network_params$carnival_options$carnival_type == 'inverse'){
      receptors_df <- NULL
    }else{
      receptors_df <- carnival_input %>% dplyr::filter(mf == 'rec')

      if(length(setdiff(receptors_df$gene_name,
                        unique(c(naive_table$source, naive_table$target))))==nrow(receptors_df)){
        stop('No mutations connected to inferred proteins; \nSuggestion: Try expanding the PKN (indirect=T) or increasing layers number')
      }
    }

    carnival_output <- SignalingProfiler::run_carnival_and_create_graph(source_df = receptors_df,
                                                                        target_df = carnival_input %>% dplyr::filter(mf != 'rec'),
                                                                        naive_network = unique(naive_table),
                                                                        proteins_df = carnival_input,
                                                                        organism = 'human',
                                                                        files =TRUE,
                                                                        carnival_options = carnival_params,
                                                                        direct = network_params$PKN_options$direct,
                                                                        with_atlas = network_params$PKN_options$with_atlas,
                                                                        path_sif = paste0(network_params$carnival_options$opt_path, '.sif'),
                                                                        path_rds =  paste0(network_params$carnival_options$opt_path, '.RDS'))


    # Two shots CARNIVAL
  }else if(network_params$carnival_options$carnival_type == 'vanilla_two_shots'){

    # FIRST RUN: RECEPTOR to KIN, PHOS, OTHERS
    receptors_df <- carnival_input %>% dplyr::filter(mf == 'rec')

    target1_df <- carnival_input %>%
      dplyr::filter(mf %in% c('kin', 'phos', 'other'))

    output1 <- SignalingProfiler::run_carnival_and_create_graph(source_df = receptors_df,
                                                                target_df = target1_df,
                                                                naive_network = unique(naive_table),
                                                                proteins_df = carnival_input,
                                                                organism = 'human',
                                                                carnival_options = carnival_params,
                                                                files = FALSE,
                                                                direct = network_params$PKN_options$direct,
                                                                with_atlas = network_params$PKN_options$with_atlas)



    # SECOND RUN: from KIN, PHOS, OTHERS to TFs
    run1_output_nodes <- SignalingProfiler::convert_output_nodes_in_next_input(output1)
    run1_output_nodes$UNIPROT <- ''

    source_df <- run1_output_nodes %>%
      dplyr::filter(mf %in% c('kin', 'phos', 'other'))

    target2_df <- carnival_input %>%
      dplyr::filter(mf == 'tf')

    output2 <- SignalingProfiler::run_carnival_and_create_graph(source_df = source_df,
                                                                target_df = target2_df,
                                                                naive_network = unique(naive_table),
                                                                proteins_df = carnival_input,
                                                                organism = 'human',
                                                                carnival_options = carnival_params,
                                                                files = FALSE,
                                                                direct = network_params$PKN_options$direct,
                                                                with_atlas = network_params$PKN_options$with_atlas)


    # UNION OF RUN1 and RUN2 graphs
    carnival_output <- SignalingProfiler::union_of_graphs(graph_1 = output1$igraph_network,
                                                          graph_2 = output2$igraph_network,
                                                          proteins_df = carnival_input,
                                                          files = TRUE,
                                                          path_sif = paste0(network_params$carnival_options$opt_path, '.sif'),
                                                          path_rds =  paste0(network_params$carnival_options$opt_path, '.RDS'))


    # Three shots CARNIVAL
  }else if(network_params$carnival_options$carnival_type == 'vanilla_three_shots'){

    # FIRST RUN: RECEPTOR to KIN, PHOS, OTHERS
    receptors_df <- carnival_input %>%
      dplyr::filter(mf == 'rec')

    target1_df <- carnival_input %>%
      dplyr::filter(mf %in% c('kin', 'phos', 'other'))

    output1 <- SignalingProfiler::run_carnival_and_create_graph(source_df = receptors_df,
                                                                target_df = target1_df,
                                                                naive_network = unique(naive_table),
                                                                proteins_df = carnival_input,
                                                                organism = 'human',
                                                                carnival_options = carnival_params,
                                                                files = FALSE,
                                                                direct = network_params$PKN_options$direct,
                                                                with_atlas = network_params$PKN_options$with_atlas)


    # SECOND RUN: from KIN, PHOS to OTHER
    run1_output_nodes <- SignalingProfiler::convert_output_nodes_in_next_input(output1)
    run1_output_nodes$UNIPROT <- ''

    source_df <- run1_output_nodes %>%
      dplyr::filter(mf %in% c('kin', 'phos'))

    target2_df <- run1_output_nodes %>%
      dplyr::filter(mf %in% c('other'))

    output2 <- SignalingProfiler::run_carnival_and_create_graph(source_df = source_df,
                                                                target_df = target2_df,
                                                                naive_network = unique(naive_table),
                                                                proteins_df = carnival_input,
                                                                organism = 'human',
                                                                carnival_options = carnival_params,
                                                                files = FALSE,
                                                                direct = network_params$PKN_options$direct,
                                                                with_atlas = network_params$PKN_options$with_atlas)


    #THIRD RUN: from KIN, PHOS to OTHER to TFs
    run2_output_nodes <- SignalingProfiler::convert_output_nodes_in_next_input(output2)
    run2_output_nodes$UNIPROT <- ''

    source2_df <- run2_output_nodes %>%
      dplyr::filter(mf %in% c('kin', 'phos', 'other'))

    target3_df <- carnival_input %>%
      dplyr::filter(mf == 'tf')

    output3 <- SignalingProfiler::run_carnival_and_create_graph(source_df = source2_df,
                                                                target_df = target3_df,
                                                                naive_network = unique(naive_table),
                                                                proteins_df = carnival_input,
                                                                organism = 'human',
                                                                carnival_options = carnival_params,
                                                                files = FALSE,
                                                                direct = network_params$PKN_options$direct,
                                                                with_atlas = network_params$PKN_options$with_atlas)


    # UNION OF RUN1, RUN2 AND RUN3 GRAPHS
    union_run1_run2 <- SignalingProfiler::union_of_graphs(graph_1 = output1$igraph_network,
                                                          graph_2 = output2$igraph_network,
                                                          proteins_df = carnival_input,
                                                          files = FALSE)$igraph_network


    carnival_output <- SignalingProfiler::union_of_graphs(graph_1 = union_run1_run2,
                                                          graph_2 = output3$igraph_network,
                                                          proteins_df = carnival_input,
                                                          files = TRUE,
                                                          path_sif = paste0(network_params$carnival_options$opt_path, '.sif'),
                                                          path_rds =  paste0(network_params$carnival_options$opt_path, '.RDS'))

  }else{
    stop('Please check CARNIVAL_type: allowed CARNIVAL types are \'vanilla_one_shot\', \'vanilla_two_shots\', \'vanilla_three_shots\', \'inverse\\')
  }

  if(is.null(carnival_output)){
    return(carnival_output)
  }

  # Map phosphoproteomics on CARNIVAL output if available
  if(!is.null(phosphoproteomics)){
    carnival_output <- SignalingProfiler::expand_and_map_edges(optimized_object = carnival_output,
                                                               organism = 'human',
                                                               phospho_df = phosphoproteomics,
                                                               files = save_all_files,
                                                               direct = network_params$PKN_options$direct,
                                                               with_atlas = network_params$PKN_options$with_atlas,
                                                               path_sif = paste0(network_params$carnival_options$opt_path, '_val.sif'),
                                                               path_rds =  paste0(network_params$carnival_options$opt_path, '_val.RDS'))
  }
  return(carnival_output)
}
