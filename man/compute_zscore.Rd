% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_zscore.R
\name{compute_zscore}
\alias{compute_zscore}
\title{compute_zscore}
\usage{
compute_zscore(omic_matrix, by, metric = "median")
}
\arguments{
\item{omic_matrix}{numeric matrix, where rows represent features (e.g., genes, proteins) and columns represent samples.}

\item{by}{string, it indicates the orientation for Z-score computation. Must be either "row" or "column".}

\item{metric}{string, it specifies the centering metric to use for Z-score computation. Options are "median" (default) or "mean".}
}
\value{
numeric matrix with the same dimensions as omic_matrix, where values have been scaled to Z-scores.
}
\description{
This function computes Z-scores for an omics data matrix either by rows or by columns,
using a specified centering metric (mean or median).
}
\examples{
omic_data <- matrix(rnorm(20), nrow = 4, ncol = 5)
rownames(omic_data) <- paste0("Gene", 1:4)
colnames(omic_data) <- paste0("Sample", 1:5)

zscores_row <- compute_zscore(omic_data, by = "row", metric = "median")

zscores_col <- compute_zscore(omic_data, by = "column", metric = "mean")

}
