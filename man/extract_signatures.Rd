% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_signatures.R
\name{extract_signatures}
\alias{extract_signatures}
\title{Extract signaling-driven transcriptomics signature from communities}
\usage{
extract_signatures(
  base_path,
  transcriptomics_file,
  output_dir,
  padj_thres = 0.01,
  mean_exp_clus_thres = 0,
  max_val = 50,
  cache = FALSE
)
}
\arguments{
\item{base_path}{Character. Path to the directory containing community folders.}

\item{transcriptomics_file}{Character. Path to the transcriptomic data file.}

\item{output_dir}{Character. Path to the directory containing all function's results.}

\item{padj_thres}{Numeric. Adjusted p-value threshold for filtering significant results. Default: \code{0.01}.}

\item{mean_exp_clus_thres}{Numeric. Minimum mean expression threshold for gene filtering. Default: \code{0}.}

\item{max_val}{Integer. Maximum number of genes to include in each community-specific signature. Default: \code{50}.}

\item{cache}{Boolean. If TRUE uses Anova_result.tsv file in \code{output_dir} if available. Default: \code{FALSE}.}
}
\value{
A data frame containing the results of the transcriptomic signature for each community.
The function also saves a csv file with the anova results and generates a "Signatures" directory where filtered signature files for each community are stored.
}
\description{
This function analyzes transcriptomic data by performing ANOVA-Tukey on gene expression levels and then generating community-specific gene expression signatures
}
\examples{
extract_signatures(base_path = "./vignette/Communities/output_communities/",
transcriptomics_file = "./PatientProfiler_processed_input/Transcriptomics_updated.tsv",
output_dir = "./Anova_result.tsv",
padj_thres = 0.01,
mean_exp_clus_thres = 0,
max_val = 50,
cache = F)

extract_signatures(base_path = "./vignette/Communities/output_communities/",
transcriptomics_file = "./PatientProfiler_processed_input/Transcriptomics_updated.tsv",
output_dir = "./Anova_result.tsv",
padj_thres = 0.01,
mean_exp_clus_thres = 0,
max_val = 10,
cache = T)

}
