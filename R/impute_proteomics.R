#' impute_proteomics 
#'
#' Impute missing values in patient/sample columns using MICE.
#' If a `group` vector is provided (e.g., Tumor/Normal), the predictor matrix is
#' "masked" so that samples are imputed only using predictors from the same group.
#'
#' @param df data.frame. Columns before `start_column` are identifiers
#'   (e.g., `gene_name`, `UNIPROT`).
#' @param start_column integer. Index of the first patient/sample column to impute.
#' @param imp_method character or named vector. MICE imputation method
#'   (default: `"pmm"`). If a named vector, names must match patient columns.
#' @param m integer. Number of multiple imputations (default: 5).
#' @param maxit integer. Number of MICE iterations (default: 5).
#' @param seed integer. Random seed (default: 103).
#' @param collapse character. How to combine the m completed datasets:
#'   `"median"`, `"mean"`, or `"none"` (default: `"median"`).
#' @param preserve_observed logical. If TRUE (default), only missing entries
#'   are replaced by imputed values; observed entries remain unchanged.
#' @param clean_patient_names logical. If TRUE (default), temporarily sanitize
#'   patient/sample column names for MICE (alphanumeric + underscores).
#' @param group optional named character vector:
#'   names = ORIGINAL patient column names,
#'   values = group label (e.g., `"Tumor"`, `"Normal"`).  
#'   If provided, the predictor matrix will be masked within each group.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{imputed_df}: final imputed data.frame (if collapse != "none")
#'     \item \code{mids}: the MICE mids object
#'     \item \code{completed_list}: list of m completed datasets
#'   }
#' @import mice



impute_proteomics <- function(df,
                              start_column = 1,
                              imp_method = "pmm",
                              m = 5,
                              maxit = 5,
                              seed = 103,
                              collapse = c("median","mean","none"),
                              preserve_observed = TRUE,
                              clean_patient_names = TRUE,
                              group = NULL) {
  
  collapse <- match.arg(collapse)
  stopifnot(is.data.frame(df))
  stopifnot(start_column >= 1, start_column <= ncol(df))
  
  id_cols <- if (start_column > 1) colnames(df)[1:(start_column - 1)] else character(0)
  patient_cols_orig <- colnames(df)[start_column:ncol(df)]
  stopifnot(length(patient_cols_orig) > 0)
  
  # Helper: coerce to numeric and replace non-finite values with NA
  to_numeric_nonfinite_to_na <- function(x) {
    v <- suppressWarnings(as.numeric(x))
    v[!is.finite(v)] <- NA_real_
    v
  }
  
  # Working copy
  df_work <- df
  df_work[patient_cols_orig] <- lapply(df_work[patient_cols_orig], to_numeric_nonfinite_to_na)
  
  # Optionally sanitize patient/sample column names for MICE
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
  
  # Patient-only block
  df_pats <- df_work[, patient_cols_clean, drop = FALSE]
  df_pats[] <- lapply(df_pats, function(x){ x[!is.finite(x)] <- NA_real_; x })
  
  set.seed(seed)
  init <- mice::mice(df_pats, maxit = 0, printFlag = FALSE)
  meth <- init$method
  predM <- init$predictorMatrix
  
  # Assign imputation method
  if (!is.null(imp_method)) {
    if (length(imp_method) == 1L) {
      meth[colnames(df_pats)] <- imp_method
    } else if (!is.null(names(imp_method))) {
      common <- intersect(names(imp_method), colnames(df_pats))
      meth[common] <- imp_method[common]
    }
  }
  
  # --- MASKING: restrict predictors to within-group only ---
  if (!is.null(group)) {
    stopifnot(all(names(group) %in% patient_cols_orig))
    grp_clean <- setNames(group[names(name_map)], name_map[names(group)])  # map to cleaned names
    for (j in colnames(df_pats)) {
      if (!j %in% names(grp_clean)) next
      same_grp <- names(grp_clean)[grp_clean == grp_clean[j]]
      other_grp <- setdiff(colnames(df_pats), same_grp)
      predM[j, other_grp] <- 0L
      predM[j, j] <- 0L  # never self-predict
    }
  } else {
    diag(predM) <- 0L
  }
  
  # Run MICE
  mids <- mice::mice(df_pats, m = m, method = meth, predictorMatrix = predM,
                     maxit = maxit, seed = seed, printFlag = TRUE)
  
  completed_list <- mice::complete(mids, action = "all")
  
  if (collapse == "none") {
    return(list(imputed_df = NULL, mids = mids, completed_list = completed_list))
  }
  
  # Collapse cell-wise (median/mean) across m imputations
  completed_list <- lapply(completed_list, function(d) {
    d[] <- lapply(d, function(x){ x[!is.finite(x)] <- NA_real_; x }); d
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
  
  # Reassemble final dataframe (ID columns + imputed patients)
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
