
#' Visualize Community Network in Cytoscape
#'
#' This function loads node and edge data for a specified community,
#' processes the data, creates a graph, and visualizes it in Cytoscape using
#' the PatientProfiler predefined style.
#'
#' @param community_name A string specifying the community name (e.g., "community_7").
#' @param comm_folder (optional) Path to the folder with resulting community from previous step. Default: `"./Communities/output_communities/"`.
#'
#' @return None. The function creates a network visualization in Cytoscape.
#' @export
#'
#' @examples
#'
#' visualize_community("community_7")


visualize_community <- function(community_name, comm_folder = "./Communities/output_communities/") {
  # Extract the community number
  community_number <- stringr::str_extract(community_name, "\\d+$")

  # Construct file paths
  path <- paste0(comm_folder, community_name, "/")
  node_file <- paste0(path, "nodes_", community_number, ".csv")
  edge_file <- paste0(path, "edges_", community_number, ".csv")

  # Load data
  node <- readr::read_csv(node_file)
  edge <- readr::read_csv(edge_file)

  # Modify the gene_name column and assign gene_id
  node$gene_id <- node$gene_name
  node <- node %>%
    dplyr::relocate(gene_id) %>%
    dplyr::mutate(gene_name = stringr::str_remove(gene_name, "_INHIBITED"))

  # Create the graph
  graph <- igraph::graph_from_data_frame(d = edge, vertices = node, directed = TRUE)

  # Import visualization style
  data_path <- system.file("extdata", "pp_style.xml", package = "PatientProfiler")
  RCy3::importVisualStyles(filename = data_path)

  # Create the network in Cytoscape
  RCy3::createNetworkFromIgraph(igraph = graph,
                                title = community_name,
                                collection = "Communities")

  # Set visual style
  RCy3::setVisualStyle("PatientProfiler_communities_style")
}
