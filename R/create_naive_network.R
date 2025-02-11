#' create_naive_network
#'
#' Generates hierarchical (multi)layered network from source nodes (perturbed nodes) defining
#' layers by distinguishing inferred proteins by molecular function.
#'
#' *Results output*
#'  The function:
#'  - returns an *igraph object* representing the naive network
#'  - creates internally an *RDS file* and a *SIF file* of the naive network at `output_dir` directory
#'
#' @param PKN dataframe obtained from the function get_PKN representing PKN interactions
#' @param sources dataframe of perturbed nodes (e.g., patient mutations), starting point of the model
#' @param activities dataframe of proteins inferred activity to connect to the starting point
#' @param output_dir string, output naive network directory
#' @param naive_network_parameters A list of parameters for naive network construction. Available parameters are:
#'  - `layers`: integer, specifiying the number of layers > 0 and < 4; default  `2`
#'  - `max_length` a vector of int or an int, length == layers; default  `c(1,4)`
#'  - `connect_all` Boolean, whether connecting intermediate nodes of shortest path; default  `TRUE`
#'  - `naive_path` path of the folder where naive network files are saved; default `./Network_output/Naive_`
#'
#' @return igraph object representing the naive network
#'
#' @examples
#' omics_list <- list(transcriptomics, proteomics)
#' PKN = get_PKN(omics_list = omics_list)
#' create_naive_network(PKN, sources, activities, naive_network_parameters = list(layers = 2, max_length = c(1,4)))


create_naive_network <- function(PKN,
                                 sources,
                                 activities,
                                 output_dir,
                                 naive_network_parameters = list()){

  # Validate input parameters and generate sources vector
  sources_gn <- validate_naive_network_parameters(PKN, sources, activities, output_dir, naive_network_parameters)

  message(paste0('Building ', naive_network_parameters$layers, '-layered naive network...'))

  if(naive_network_parameters$layers == 1){

    naive_network <- SignalingProfiler::one_layer_naive_network(starts_gn = sources_gn,
                                                                targets_gn = activities$gene_name,
                                                                PKN_table = PKN, #or PKN_human
                                                                max_length = naive_network_parameters$max_length[1],
                                                                connect_all = naive_network_parameters$connect_all,
                                                                files = TRUE,
                                                                sif_path = paste0(naive_network_parameters$naive_path, '.sif'),
                                                                rds_path = paste0(naive_network_parameters$naive_path, '.rds')
    )
  }else if(naive_network_parameters$layers == 2){

    kin_phos_other <- activities %>% dplyr::filter(mf %in% c('kin', 'phos', 'other'))
    tfs <- activities %>% dplyr::filter(mf == 'tf')

    naive_network <- SignalingProfiler::two_layer_naive_network(starts_gn = sources_gn,
                                                                intermediate_gn = kin_phos_other$gene_name,
                                                                targets_gn = tfs$gene_name,
                                                                PKN_table = PKN,
                                                                max_length_1 =  naive_network_parameters$max_length[1],
                                                                max_length_2 =  naive_network_parameters$max_length[2],
                                                                connect_all = naive_network_parameters$connect_all,
                                                                files = TRUE,
                                                                sif_path = paste0(naive_network_parameters$naive_path, '.sif'),
                                                                rds_path = paste0(naive_network_parameters$naive_path, '.rds')
    )
  }else if(naive_network_parameters$layers == 3){

    # divide proteins according to the molecular function
    kin_phos <- activities %>% dplyr::filter(mf %in% c('kin', 'phos'))
    other <- activities %>% dplyr::filter(mf == 'other')
    tfs <- activities %>% dplyr::filter(mf == 'tf')

    # create the naïve network
    naive_network <- SignalingProfiler::three_layer_naive_network(starts_gn = sources_gn,
                                                                  intermediate1_gn = kin_phos$gene_name,
                                                                  intermediate2_gn = other$gene_name,
                                                                  targets_gn = tfs$gene_name,
                                                                  PKN_table = PKN,
                                                                  max_length_1 =  naive_network_parameters$max_length[1],
                                                                  max_length_2 =  naive_network_parameters$max_length[2],
                                                                  max_length_3 =  naive_network_parameters$max_length[3],
                                                                  connect_all =  naive_network_parameters$connect_all,
                                                                  files = TRUE,
                                                                  sif_path = paste0(naive_network_parameters$naive_path, '.sif'),
                                                                  rds_path = paste0(naive_network_parameters$naive_path, '.rds')
    )


  }else{stop('Please provide a number of layers > 0 and < 4!')}

  return(naive_network)
}

#' validate_naive_network_parameters
#'
#' Validation function for naive network parameters to plain error handling
#'
#' @param PKN dataframe obtained from the function get_PKN representing PKN interactions
#' @param sources dataframe of perturbed nodes (e.g., patient mutations), starting point of the model
#' @param activities dataframe of proteins inferred activity to connect to the starting point
#' @param output_dir string, output naive network directory
#' @param naive_network_parameters A list of parameters for naive network construction. Available parameters are:
#'  - `layers`: integer, specifiying the number of layers > 0 and < 4; default  `2`
#'  - `max_length` a vector of int or an int, length == layers; default  `c(1,4)`
#'  - `connect_all` Boolean, whether connecting intermediate nodes of shortest path; default  `TRUE`
#'  - `naive_path` path of the folder where naive network files are saved; default `./Network_output/Naive_`
#'
#' @return sources_gn


validate_naive_network_parameters <- function(PKN,
                                              sources,
                                              activities,
                                              output_dir,
                                              naive_network_parameters) {
  if (!is.list(naive_network_parameters)) {
    stop("'naive_network_parameters' must be a list.")
  }

  if (!all(c("layers", "max_length", "connect_all", "naive_path") %in% names(naive_network_parameters))) {
    stop("'naive_network_parameters' must contain 'layers', 'max_length', 'connect_all', and 'naive_path'.")
  }

  if (!is.numeric(naive_network_parameters$layers) || naive_network_parameters$layers %% 1 != 0 ||
      naive_network_parameters$layers < 1 || naive_network_parameters$layers > 3) {
    stop("'layers' must be an integer between 1 and 3.")
  }

  if (!is.numeric(naive_network_parameters$max_length) ||
      length(naive_network_parameters$max_length) != naive_network_parameters$layers) {
    stop(paste0("'max_length' must be a numeric vector with length equal to 'layers' (", naive_network_parameters$layers, ")."))
  }

  if (!is.logical(naive_network_parameters$connect_all)) {
    stop("'connect_all' must be a logical value (TRUE or FALSE).")
  }

  if (!is.character(naive_network_parameters$naive_path)) {
    stop("'naive_path' must be a character string.")
  }

  if (!is.data.frame(sources) || !"Patient_ID" %in% colnames(sources)) {
    stop("'sources' must be a data frame containing a 'Patient_ID' column.")
  }

  if (!is.data.frame(activities) || !all(c("gene_name", "mf") %in% colnames(activities))) {
    stop("'activities' must be a data frame containing 'gene_name' and 'mf' columns.")
  }

  if (!is.data.frame(PKN)) {
    stop("'PKN' must be a data frame.")
  }

  # Generate the a vector of sources for naive network
  sources %>% column_to_rownames('Patient_ID') -> sources_mat
  sources_gn <- colnames(sources_mat)[sources_mat != 0]

  if (length(sources_gn) == 0) {
    stop("No valid sources identified in the 'sources' input.")
  }

  return(sources_gn)
}



