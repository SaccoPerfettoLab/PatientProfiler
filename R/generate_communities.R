
#' Generate communities from patient-specif mechanistic models
#'
#' @param dir_path Character. Path to the directory containing the patients_stratification.tsv in Subtype
#' @param network_dir Character. Path to the directory containing patient-specific mechanistic models obtained in Step 3
#' @param output_dir Character. Path to the directory containing the output of this function
#' @param t_lower Numeric. The lower threshold for filtering edges based on their frequency across patients. Default is 4.
#' @param t_upper Numeric. The upper threshold for filtering edges based on their frequency across patients.  If not provided, no upper limit is applied, and all edges above `t_lower` will be included.
#' @param local Boolean. 
#'
#'
#' @return A directory for each community containing nodes table, edges table e patients stratification.
#' @export
#'
#' @examples
#' generate_communities(dir_path = "input_communities", 
#' network_dir = "Networks_output",
#' output_dir = "output_communities")

generate_communities <- function(dir_path, network_dir, output_dir, t_lower, t_upper, local = FALSE) {
  if (local == TRUE) {
    path_package <- './inst/'
  } else {
    path_package <- paste0(.libPaths(), '/PatientProfiler/')
  }
  
  # Import the Python script to find communities
  for (path in path_package) {
    result <<- tryCatch({
      reticulate::source_python(paste0(path, "/python/find_communities.py"))
    }, error = function(e) {
      message("An error occurred: ", e$message, ' with path ', path)
    })
  }
  
  # Take as inputs patient-specific mechanistic models
  file_list <- list.files(network_dir, pattern = "Pheno_(.*)_network\\.RDS", full.names = TRUE) 
  
  # Nodes and Edges creation
  temp_nodes_dir <- file.path(dir_path, "nodes") 
  temp_edges_dir <- file.path(dir_path, "edges")
  dir.create(temp_nodes_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(temp_edges_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (file in file_list) {
    file_name <- basename(file)
    pat_name <- sub("Pheno_(.*)_network\\.RDS", "\\1", file_name)
    
    opt1 <- readRDS(file)
    
    nodes <- opt1$sp_object_phenotypes$nodes_df[, c(1:2,7)]
    nodes_path <- file.path(temp_nodes_dir, paste0("nodes_", pat_name, ".xlsx"))
    writexl::write_xlsx(nodes, nodes_path)
    
    edges <- opt1$sp_object_phenotypes$edges_df[, c(1:4,6)]
    edges_path <- file.path(temp_edges_dir, paste0("edges_", pat_name, ".xlsx"))
    writexl::write_xlsx(edges, edges_path)
  }
  
  # Find communities
  communities <- find_communities(dir_path = dir_path, 
                                  output_dir = output_dir,
                                  t_lower = t_lower, 
                                  t_upper = t_upper) 
  
  # Remove the temporary directory
  unlink(dir_path, recursive = TRUE)
  
  
  return(communities)
}

