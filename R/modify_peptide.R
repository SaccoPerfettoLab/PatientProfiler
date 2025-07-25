#' modify_peptide
#'
#' @param phospho_df dataframe, phosphoproteomics data
#' @param mult_col integer, number of the Multiplicity column
#' @param peptide_col name of the column containing the peptide sequence
#' @param site_col integer, number  of the Site column
#'
#' @return phosphoproteomics dataframe with just the phosphosite's amminoacid in lowercase
#'
#' @examples
#'
#'#' sample_df <- data.frame(
#'   Gene_name = c("ACIN1", "ZZEF1", "AAK1")
#'   Peptide = c("sRsRsTPVRR", "LLsFRsMEEAR", "SATTTPST*PR"),
#'   Multiplicity = c(3, 2, 1),
#'   Site = c("S245", "S30", "T100"),
#'   Value = c(0.1, 0.2, 0.3)
#' )
#'
#' modified_df <- modify_peptide(sample_df, mult_col = 3, peptide_col = 2, site_col = 4)

modify_peptide <- function(phospho_df, mult_col, peptide_col, site_col) {
  
  
  phospho_df[[peptide_col]] <- ifelse(
    stringr::str_detect(phospho_df[[peptide_col]], "^.+\\..+\\..+$"),  
    stringr::str_extract(phospho_df[[peptide_col]], "(?<=\\.).+?(?=\\.)"), 
    phospho_df[[peptide_col]] 
  )
  
  for (i in 1:nrow(phospho_df)) {

    multiplicity <- phospho_df[i, mult_col]
    peptide <- as.character(phospho_df[[peptide_col]][i])
    site <- as.character(phospho_df[i, site_col])

    if (nchar(gsub("[^a-z]", "", peptide)) > 1) {

      site_letter <- substr(site, 1, 1)
      site_position <- as.numeric(substr(site, 2, nchar(site)))

      modified_peptide <- sapply(1:nchar(peptide), function(j) {
        char <- substr(peptide, j, j)
        if (char %in% letters && toupper(char) != site_letter) {
          return(toupper(char))
        } else {
          return(char)
        }
      })

      phospho_df[[peptide_col]][i] <- paste(modified_peptide, collapse = "")
    }
  }

  return(phospho_df)
}
