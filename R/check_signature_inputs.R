#' Check input format for test_signature_on_dataset
#'
#' @param gene_signatures List of named vectors
#' @param transcriptomic_data Dataframe with gene_name column
#' @param clinical_data Dataframe with patient_id, serial_time, status columns (can be NULL)

check_signature_inputs <- function(gene_signatures, transcriptomic_data, clinical_data = NULL, enrichment_result = NULL) {
  
  if(is.null(enrichment_result)){
    # gene_signatures: list of named vectors
    if (!is.list(gene_signatures) || length(gene_signatures) == 0) {
      stop("gene_signatures must be a non-empty list of named character vectors (gene sets).")
    }
    if (any(sapply(gene_signatures, function(x) !is.character(x)))) {
      stop("Each element in gene_signatures must be a character vector (of gene names).")
    }
    if (is.null(names(gene_signatures)) || any(names(gene_signatures) == "")) {
      stop("gene_signatures must be a *named* list (each signature should have a name).")
    }
    
    # transcriptomic_data: data.frame with 'gene_name'
    if (!is.data.frame(transcriptomic_data)) {
      stop("transcriptomic_data must be a data.frame.")
    }
    if (!"gene_name" %in% colnames(transcriptomic_data)) {
      stop("transcriptomic_data must contain a 'gene_name' column.")
    }
    if (nrow(transcriptomic_data) == 0) {
      stop("transcriptomic_data cannot be empty.")
    }
  }
  
  # clinical_data: data.frame with patient_id, serial_time, status
  if (!is.null(clinical_data)) {
    if (!is.data.frame(clinical_data)) {
      stop("clinical_data must be a data.frame if provided.")
    }
    required_cols <- c("patient_id", "serial_time", "status")
    if (!all(required_cols %in% colnames(clinical_data))) {
      stop(
        sprintf(
          "clinical_data must contain the columns: %s",
          paste(required_cols, collapse = ", ")
        )
      )
    }
    if (!all(clinical_data$status %in% c(0, 1))) {
      stop("clinical_data$status must be coded as 0 (alive) or 1 (deceased).")
    }
  }
  return(invisible(TRUE))
}

