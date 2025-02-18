
.onLoad <- function(libname, pkgname) {
  reticulate::use_virtualenv("r-patientprofiler", required = FALSE)
}

