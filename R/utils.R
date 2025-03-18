
#' From df to list
#'
#' To use to convert a mutation matrix to a list for vanilla CARNIVAL
#'
#' @param df_cols df with one row and many cols
#'
#' @return list with cols as elements
#' @export
#'
#' @examples
#'
#' patient_id <- 'CPT000814'
#' mutations_df <- read_csv('./mutations.csv')
#' mut_pat <- mutations_df[mutations_df$Patient_ID == patient_id,]
#' from_df_to_list(mut_pat)
#'
from_df_to_list <- function(df_cols){

  mut_list <- as.list(setNames(df_cols, names(df_cols)))[-1] #Removing Patient_ID
  mut_list <- mut_list[ mut_list != 0 ]

  return(mut_list)
}

#' igraphToSif
#'
#' converts an igraph object to SIF file
#'
#' @param inGraph igraph object
#' @param outfile name of output file
#' @param edgeLabel column containing the SIGN information
#'
#' @return no return.
#' igraphToSif(naive_network, outfile='naive.sif', edgeLabel='sign')


igraphToSif <- function(inGraph, outfile="output.sif", edgeLabel="label") {

  if(file.exists(outfile)){
    file.remove(outfile)
  }

  singletons <- as.list(igraph::get.vertex.attribute(inGraph, "name"))
  edgeList <- igraph::get.edgelist(inGraph, names=FALSE)
  nodeNames <- igraph::get.vertex.attribute(inGraph, "name")
  edgeAttribute <- igraph::get.edge.attribute(inGraph, edgeLabel)
  numE <- igraph::ecount(inGraph)

  for (i in 1:numE) {
    node1 <- edgeList[i,1]
    node2 <- edgeList[i,2]
    singletons <- singletons[which(singletons != nodeNames[node1])]
    singletons <- singletons[which(singletons != nodeNames[node2])]
    write(paste0(nodeNames[node1], "\t", edgeAttribute[i], "\t", nodeNames[node2], "\n"),
          outfile, append = TRUE)
  }

  for (single in singletons) {
    write(paste0(single, "\n"),
          outfile, append = TRUE)
  }
}


#' extract_patient_id
#'
#' Extract patients_IDs from PatientProfiler directories structures
#'
#' @param file_path
#'
#' @return it returns a string representing patient id
#' @export
#'
#' @examples
#' extract_patient_id('Prot_Patient_X11.csv')
#'
extract_patient_id <- function(file_path) {
  basename(file_path) %>%
    sub("^(Prot_Patient_|Transc_Patient_|Phospho_Patient_|Activity_Patient_|Activity_constraints_Patient_)(.*)\\.tsv$", "\\2", .)
}
