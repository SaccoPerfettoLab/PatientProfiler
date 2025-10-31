#' Collapse rows by UniProt choosing a representative gene name with clear rules
#'
#' For each \code{UNIPROT}, multiple \code{gene_name} can appear (aliases, readthroughs).
#' This function selects a single representative gene name using the following rules:
#' \enumerate{
#'   \item \strong{Prefer names without hyphens} (readthrough/fusions like \code{A-B} are deprioritized).
#'   \item Among remaining, \strong{prefer names without dots} (versions like \code{MT-ND4.1} are deprioritized).
#'   \item If still multiple, choose the \strong{shortest} symbol.
#'   \item If tie remains, choose \strong{alphabetically}.
#' }
#'
#' If the candidate rows have \emph{identical numeric values} (within \code{tol}), then after applying (1)–(4)
#' if more than one name is still tied, their names are concatenated as \code{"GENE1;GENE2"} to preserve information.
#' If the numeric profiles are \emph{different}, a single representative is always selected.
#'
#' @param df A data.frame/tibble with columns \code{UNIPROT}, \code{gene_name}, and sample/value columns (numeric).
#' @param tol Numeric tolerance to consider two rows equal across all numeric columns. Default \code{1e-8}.
#'
#' @return A tibble collapsed by \code{UNIPROT} with one row per UniProt and a single \code{gene_name}
#'         (or a semicolon-joined name only in the residual tied-identical case). Numeric columns are taken from
#'         the selected representative; if names are concatenated because rows are identical, any of the equal rows is used.
#' @examples
#' # res <- collapse_by_uniprot_choose_gene(df_pro_zscore)
#' # head(res)
#' @export
collapse_by_uniprot_choose_gene <- function(df, tol = 1e-8) {
  stopifnot(all(c("UNIPROT","gene_name") %in% names(df)))
  suppressPackageStartupMessages({
    library(dplyr); library(stringr); library(rlang)
  })
  
  # Identify numeric columns (sample/value columns)
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  
  # Helper: apply (a)-(d) ranking to a character vector of gene names
  rank_genes <- function(v) {
    v <- unique(v)
    # (a) prefer no hyphen
    mask_no_hyphen <- !grepl("-", v)
    if (any(mask_no_hyphen)) v <- v[mask_no_hyphen]
    # (b) prefer no dot
    mask_no_dot <- !grepl("\\.", v)
    if (any(mask_no_dot)) v <- v[mask_no_dot]
    # (c)+(d) shortest, then alphabetical
    v[order(nchar(v), v)]
  }
  
  # Helper: check if all candidate rows have identical numeric values (within tol)
  all_equal_numeric <- function(mat) {
    if (ncol(mat) == 0 || nrow(mat) <= 1) return(TRUE)
    ref <- mat[1, , drop = FALSE]
    # max absolute diff vs first
    max(abs(as.matrix(mat[-1, , drop = FALSE]) - matrix(rep(as.numeric(ref), each = nrow(mat)-1),
                                                        nrow = nrow(mat)-1))) <= tol
  }
  
  df %>%
    mutate(
      UNIPROT   = as.character(UNIPROT),
      gene_name = as.character(gene_name)
    ) %>%
    group_by(UNIPROT) %>%
    group_modify(function(d, key) {
      genes <- d$gene_name
      # Rank by rules (a)-(d)
      ranked <- rank_genes(genes)
      
      # Subset candidates that match the best ranking (may be >1)
      candidates <- d %>% filter(gene_name %in% ranked)
      # Keep only those whose gene_name is among the "best" after ordering
      top_name <- ranked[1]
      top_set  <- ranked[ranked %in% candidates$gene_name]
      
      # Are numeric values identical across remaining candidates?
      if (length(num_cols) > 0) {
        mat <- as.matrix(candidates %>% select(all_of(num_cols)))
        identical_vals <- all_equal_numeric(mat)
      } else {
        identical_vals <- TRUE
      }
      
      # Decide final gene_name
      final_gene <-
        if (length(top_set) == 1) {
          top_set[1]
        } else if (identical_vals) {
          # Only in the residual tied case with identical values: concatenate
          paste(sort(unique(top_set)), collapse = ";")
        } else {
          # Values differ: pick the first by (a)-(d)
          top_set[1]
        }
      
      # Take the row for the chosen/concatenated gene_name:
      # - if concatenated (identical), any row is fine (they’re equal);
      # - otherwise, take the row matching final_gene.
      if (final_gene %in% candidates$gene_name) {
        chosen <- candidates %>% filter(gene_name == final_gene) %>% slice(1)
      } else {
        # concatenated case
        chosen <- candidates %>% slice(1)
        chosen$gene_name <- final_gene
      }
      
      chosen
    }) %>%
    ungroup()
}
