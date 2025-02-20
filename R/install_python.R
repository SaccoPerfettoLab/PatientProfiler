
#' install_pp_py
#'
#' This function must be executed the first time PatientProfiler is installed
#'
#' @param ...
#' @param envname "r-patientprofiler"
#' @param new_env used-specified name for the environment
#'
#' @return no return, r-patientprofiler virtual environment created
#' @export
#'
#' @examples
#'
#' library(PatientProfiler)
#' install_pp_py()
#'
install_pp_py <- function(..., envname = "r-patientprofiler",
                          new_env = identical(envname, "r-patientprofiler")) {

  if(new_env && reticulate::virtualenv_exists(envname)){
    reticulate::virtualenv_remove(envname)
  }

  reticulate::py_install("requests", envname = envname, ...)
  reticulate::py_install("pandas", envname = envname, ...)
  reticulate::py_install("networkx", envname = envname, ...)
  reticulate::py_install("openpyxl", envname = envname, ...)
}
