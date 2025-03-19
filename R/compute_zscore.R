#' compute_zscore
#'
#' This function computes Z-scores for an omics data matrix either by rows or by columns,
#' using a specified centering metric (mean or median).
#'
#' @param omic_matrix numeric matrix, where rows represent features (e.g., genes, proteins) and columns represent samples.
#' @param by string, it indicates the orientation for Z-score computation. Must be either "row" or "column".
#' @param metric string, it specifies the centering metric to use for Z-score computation. Options are "median" (default) or "mean".
#'
#' @return numeric matrix with the same dimensions as omic_matrix, where values have been scaled to Z-scores.
#'
#' @examples
#' omic_data <- matrix(rnorm(20), nrow = 4, ncol = 5)
#' rownames(omic_data) <- paste0("Gene", 1:4)
#' colnames(omic_data) <- paste0("Sample", 1:5)
#'
#' zscores_row <- compute_zscore(omic_data, by = "row", metric = "median")
#'
#' zscores_col <- compute_zscore(omic_data, by = "column", metric = "mean")
#'
#' @export



compute_zscore <- function(omic_matrix, by, metric = 'median'){


  calculate_zscore_col <- function(mat, metric = 'median') {
    scaled <- matrix(nrow = nrow(mat), ncol = ncol(mat))
    for (i in seq_len(ncol(mat))) {
      column <- mat[, i]

      if (metric == "median") {
        centered <- column - median(column, na.rm = TRUE)
      } else {
        centered <- column - mean(column, na.rm = TRUE)
      }

      scaled[, i] <- centered / sd(centered, na.rm = TRUE)
    }
    colnames(scaled) <- colnames(mat)
    rownames(scaled) <- rownames(mat)
    return(scaled)
  }


  calculate_zscore_row <- function(mat, metric = 'median') {
    centering_func <- ifelse(metric == "median", median, mean)
    centered <- t(apply(mat, 1, function(x) x - centering_func(x, na.rm = TRUE)))
    # centered <- t(apply(mat, 1, function(x) x - median(x, na.rm = TRUE)))
    scaled <- t(apply(centered, 1, function(x) x / sd(x, na.rm = TRUE)))
    return(scaled)

  }

  if(by == 'row'){
    if(metric %in% c('median', 'mean')){
      omic_z_mat <- calculate_zscore_row(omic_matrix, metric)
    }else{
      stop('Please provide a valid metric: mean or median')
    }

  }else if( by == 'column'){
    if(metric %in% c('median', 'mean')){
      omic_z_mat <- calculate_zscore_col(omic_matrix, metric)
    }else{
      stop('Please provide a valid metric: mean or median')
    }
  }else{
    stop('Please provide a valid type of analysis: column or row')
  }

  return(omic_z_mat)
}



