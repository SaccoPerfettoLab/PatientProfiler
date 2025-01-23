#' get_PKN
#'
#' Selects the desired SIGNOR and PhosphoSitePlus built-in interactome
#'
#' *Results output*
#' A dataframe representing SIGNOR and PhosphoSitePlus built-in interactome of SignalingProfiler package
#'
#' @param PKN_params A list of parameters for PKN selection. Available parameters are:
#'    - `preprocess`: Boolean, whether to remove interactions involving entities not detected in at least one omic data; default `TRUE`.
#'    - `direct`: Boolean, whether keep only direct interactions; default `TRUE`.
#'    - `organism`: string, available are 'human' or 'mouse'; default `human`.
#'    - `with_atlas`: Boolean, whether integratins Ser/Thr and Tyr atlas inferred interactions; default `FALSE`.
#'    - `custom`: Boolean, whether the user wants to provide a custom PKN;  default `FALSE`.
#'    - `custom_path`: string, path to the custom PKN;  default `FALSE`.
#'
#' @param omics_list A list of omics data to use to preprocess PKN.
#'    Each omic file should contain the `gene_name` column reporting the Gene Symbol.
#'
#' @return Returns a dataframe representing the interactions of in the PKN
#' @examples
#' omics_list <- list(transcriptomics, proteomics)
#' PKN_params = list(preprocess = TRUE, direct = TRUE, organism = 'human')


get_PKN <- function(PKN_params = list(),
                    omics_list = NULL){

  message('Selecting PKN...')
  # Separate preprocess from the other parameters
  preprocess <- PKN_params$preprocess
  PKN_params <- PKN_params[  !'preprocess' == names(PKN_params)]

  default_PKN_params = list(organism = 'human', direct = TRUE, with_atlas = FALSE, custom = FALSE, custom_path = FALSE)
  PKN_params = modifyList(default_PKN_params, PKN_params)
  PKN = do.call(choose_PKN, PKN_params)

  if(preprocess){
    message('Preprocessing PKN...')
    PKN =  SignalingProfiler::preprocess_PKN(omics_data = omics_list, PKN_table = PKN)
  }

  return(PKN)
}
