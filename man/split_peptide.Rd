% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_peptide.R
\name{split_peptide}
\alias{split_peptide}
\title{split_peptide}
\usage{
split_peptide(phospho_df, mult_col, peptide_col)
}
\arguments{
\item{phospho_df}{tibble, the phosphoproteomics dataframe}

\item{mult_col}{integer, the number of the Multiplicity column}

\item{peptide_col}{name of the column containing the peptide sequence}
}
\value{
phosphoproteomics dataframe with one Petide for each possible phosphosite
}
\description{
split_peptide
}
\examples{

sample_df <- data.frame(
  Gene_name = c("ACIN1", "ZZEF1", "AAK1"),
  Multiplicity = c(3, 2, 1),
  Peptide = c("sRsRsTPVRR", "LLsFRsMEEAR", "SATTTPST*PR"),
  Value = c(0.1, 0.2, 0.3)
)

expanded_df <- choose_phosphopeptide(sample_df, mult_col = 2, peptide_col = "Peptide")
}
