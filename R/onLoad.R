.onLoad <- function(libname, pkgname) {
  install_bioconductor_dependencies()
}

install_bioconductor_dependencies <- function() {
  bioc_pkgs <- c("org.Mm.eg.db", "fgsea",
                 "org.Hs.eg.db", "AnnotationDbi")  # List all required Bioconductor packages

  # Identify missing packages
  missing_pkgs <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]

  if (length(missing_pkgs) > 0) {
    message("Installing missing Bioconductor packages: ", paste(missing_pkgs, collapse = ", "))

    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    BiocManager::install(missing_pkgs, ask = FALSE)
  }
}
