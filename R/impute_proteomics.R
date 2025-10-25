#' impute_proteomics
#'
#' Impute missing values in patient/sample columns using mice.
#' Optionally collapse the completed datasets into a single imputed matrix by taking the
#' cell-wise median (default) or mean.
#'
#' @param df data.frame. Columns before \code{start_column} are identifiers (e.g., \code{gene_name}, \code{UNIPROT}).
#' @param start_column integer. Index of the first patient/sample column to impute (default: \code{1}).
#' @param imp_method character or named vector. \code{mice} imputation method for patient columns
#'   (default: \code{"pmm"}). If a named vector, names must match patient column names (after
#'   optional sanitization).
#' @param m integer. Number of multiple imputations (default: \code{5}).
#' @param maxit integer. Number of \code{mice} iterations (default: \code{5}).
#' @param seed integer. Random seed (default: \code{103}).
#' @param collapse character. How to combine the \eqn{m} completed datasets:
#'   one of \code{c("median","mean","none")}. If \code{"none"}, returns only \code{mids} and the
#'   list of completed datasets; no collapsed \code{imputed_df} is produced. Default: \code{"median"}.
#' @param preserve_observed logical. If \code{TRUE} (default), only original NAs are replaced
#'   by the collapsed values; observed entries remain unchanged. If \code{FALSE}, the collapsed
#'   matrix replaces the entire patient block.
#' @param clean_patient_names logical. If \code{TRUE} (default), temporarily sanitize patient column names
#'   for \code{mice} (alphanumeric + underscores, uniqueness guaranteed) and restore original names
#'   in the returned data.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{imputed_df}: the final imputed data.frame (only if \code{collapse != "none"}), with
#'     identifiers followed by patient columns and \emph{no NAs} in patient columns.
#'   \item \code{mids}: the \code{mids} object returned by \code{mice}.
#'   \item \code{completed_list}: a list of length \code{m} containing the \code{m} completed
#'     patient-only data.frames (\code{complete(mids, "all")}).
#' }
#' @examples
#' # Example: df has columns: gene_name, UNIPROT, then patient columns
#' res <- impute_proteomics(
#'   df                 = df_raw,
#'   start_column       = 3,          # patients start at column 3
#'   imp_method         = "pmm",
#'   m                  = 5,
#'   maxit              = 5,
#'   seed               = 103,
#'   collapse           = "median",   # cell-wise median across m imputations
#'   preserve_observed  = TRUE,
#'   clean_patient_names = TRUE
#' )
#'
#' df_final <- res$imputed_df   # ready for downstream (e.g., z-scores, VIPER)
#'
#' @import mice
impute_proteomics <- function(df,
                              start_column = 1,
                              imp_method = "pmm",
                              m = 5,
                              maxit = 5,
                              seed = 103,
                              collapse = c("median","mean","none"),
                              preserve_observed = TRUE,
                              clean_patient_names = TRUE) {
  
  collapse <- match.arg(collapse)
  stopifnot(is.data.frame(df))
  stopifnot(start_column >= 1, start_column <= ncol(df))
  
  id_cols <- if (start_column > 1) colnames(df)[1:(start_column - 1)] else character(0)
  patient_cols_orig <- colnames(df)[start_column:ncol(df)]
  stopifnot(length(patient_cols_orig) > 0)
  
  # Helper: coerce numeric and turn non-finite into NA (only for patient columns)
  to_numeric_nonfinite_to_na <- function(x) {
    v <- suppressWarnings(as.numeric(x))
    v[!is.finite(v)] <- NA_real_
    v
  }
  
  # Working copy
  df_work <- df
  
  # Coerce patient columns to numeric; replace Inf/-Inf/NaN with NA
  df_work[patient_cols_orig] <- lapply(df_work[patient_cols_orig], to_numeric_nonfinite_to_na)
  
  # Optionally sanitize patient names for mice; keep map to restore
  if (isTRUE(clean_patient_names)) {
    cleaned_names <- vapply(patient_cols_orig, function(x) {
      y <- gsub("[^A-Za-z0-9_]", "_", as.character(x))
      if (!nzchar(y)) y <- "X"
      y
    }, character(1))
    cleaned_names <- make.unique(cleaned_names, sep = "_")
    name_map <- setNames(cleaned_names, patient_cols_orig)
    colnames(df_work)[start_column:ncol(df_work)] <- cleaned_names
    patient_cols_clean <- cleaned_names
  } else {
    name_map <- setNames(patient_cols_orig, patient_cols_orig)
    patient_cols_clean <- patient_cols_orig
  }
  
  # MICE on patient-only block (safety: ensure no non-finites)
  df_pats <- df_work[, patient_cols_clean, drop = FALSE]
  df_pats[] <- lapply(df_pats, function(x){ x[!is.finite(x)] <- NA_real_; x })
  
  set.seed(seed)
  init <- mice::mice(df_pats, maxit = 0, printFlag = FALSE)
  meth <- init$method
  predM <- init$predictorMatrix
  
  if (!is.null(imp_method)) {
    if (length(imp_method) == 1L) {
      meth[colnames(df_pats)] <- imp_method
    } else if (!is.null(names(imp_method))) {
      common <- intersect(names(imp_method), colnames(df_pats))
      meth[common] <- imp_method[common]
    }
  }
  
  mids <- mice::mice(df_pats, m = m, method = meth, predictorMatrix = predM,
                     maxit = maxit, seed = seed, printFlag = TRUE)
  
  completed_list <- mice::complete(mids, action = "all")  # list of length m
  
  if (collapse == "none") {
    return(list(
      imputed_df = NULL,
      mids = mids,
      completed_list = completed_list
    ))
  }
  
  # Collapse cell-wise across m imputations (median or mean)
  # Safety: ensure completed sets have no non-finites
  completed_list <- lapply(completed_list, function(d) {
    d[] <- lapply(d, function(x){ x[!is.finite(x)] <- NA_real_; x })
    d
  })
  arr <- simplify2array(lapply(completed_list, function(d) as.matrix(d)))
  collapsed <- if (collapse == "median") {
    apply(arr, c(1, 2), median, na.rm = TRUE)
  } else {
    apply(arr, c(1, 2), mean, na.rm = TRUE)
  }
  
  # Preserve observed values if requested
  df_obs <- as.matrix(df_work[, patient_cols_clean, drop = FALSE])
  if (isTRUE(preserve_observed)) {
    na_idx <- is.na(df_obs)
    df_mix <- df_obs
    df_mix[na_idx] <- collapsed[na_idx]
  } else {
    df_mix <- collapsed
  }
  
  # Reassemble: ID columns + imputed patients; restore original patient names
  imputed_df <- data.frame(df_work[, id_cols, drop = FALSE],
                           as.data.frame(df_mix), check.names = FALSE)
  colnames(imputed_df)[seq_along(patient_cols_clean) + length(id_cols)] <- names(name_map)
  imputed_df <- imputed_df[, c(id_cols, patient_cols_orig), drop = FALSE]
  
  list(
    imputed_df = imputed_df,
    mids = mids,
    completed_list = completed_list
  )
}
