% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phosphoproteomics_update.R
\name{phosphoproteomics_update}
\alias{phosphoproteomics_update}
\title{phosphoproteomics_update}
\usage{
phosphoproteomics_update(
  df_pho,
  sw_len = 7,
  uniprot_idX = NULL,
  pep_col_name = NULL,
  imp_method = NULL,
  zscore = TRUE,
  zmethod = "column",
  metric = "median",
  output_dir
)
}
\arguments{
\item{imp_method}{string, the method to use for imputation (default: "pmm", but you can choose between "pmm", "norm", "norm.nob",
"regression", "ri", "logreg", "polyreg", "predictive", "polr", "sample", "cart", "knn", "rf").}

\item{output_dir}{a string indicating the updated output folder.}

\item{df}{a data frame containing the phosphoproteomics data.}

\item{seq_len_i}{an integer specifying the length of the sequence window to be considered around the phosphopeptide.}

\item{uniprot_idx}{optional integer indicating the column index for UNIPROT IDs to be added at the uniprot column retrieved with AnnotationDbi (default is NULL).}

\item{peptide_col_name}{(Optional) A string representing the name of the peptide column in the data frame.
If not provided, the function assumes no peptide column exists.}
}
\value{
A data frame with the processed phosphoproteomics data, including modified peptide sequences,
updated amino acid positions, and sequence windows. If no peptide column is provided, the function
directly computes amino acid positions.
}
\description{
This function processes phosphoproteomics data by performing several steps of data cleaning,
peptide modification, and enrichment of sequence information. It is designed to handle datasets
with phosphopeptides, modifying peptide sequences, adding multiplicity information, and computing
sequence windows around the phosphopeptides.
}
\examples{
# Example usage:
updated_phospho_df <- update_phospho(phospho_data, site_col = 3, gn_idx = 2, seq_len_i = 7, peptide_col_name = "Peptide")
head(updated_phospho_df)
}
