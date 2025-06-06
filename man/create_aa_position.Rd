% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_aa_position.R
\name{create_aa_position}
\alias{create_aa_position}
\title{create_aa_position}
\usage{
create_aa_position(df, site_col)
}
\arguments{
\item{site_col}{integer, the number of the Site column}

\item{phospho_df}{tibble, the phosphoproteomics dataframe}
}
\value{
phosphoproteomics dataframe with the additional columns 'position' and 'aminoacid' of each phosphosite
}
\description{
create_aa_position
}
\examples{

sample_df <- data.frame(
  Gene_name = c("AAAS", "AAGAB", "AAK1"),
  Site = c("S462", "S202", "S21"),
  Value = c(0.1, 0.2, 0.3)
)

updated_df <- create_aa_position(sample_df, site_col = 2)
}
