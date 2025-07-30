#' test_signature_on_dataset
#'
#' Run fgsea on user-specified set of *gene_signatures* using transcriptomics data.
#' If *enrichment_result* is provided this step is skipped.
#' If *clinical_data* are provided, it draws a Kaplain-Meier curve assigning each patient
#' to the most-enriched signature (lowest padj and higher NES).
#'
#' @param gene_signatures list of named vector specifying the genes for each signature
#' @param transcriptomic_data dataframe, where Gene Symbols are specified in *gene_name* column
#' @param clinical_data dataframe, with *patient_id*, *serial_time* (days_to_death or days_to_last_followup), *status* (0 alive, 1 decesead)
#' @param NES_filter Boolean, if *TRUE* only positive enrichment is considered
#' @param pval_cutoff padj cutoff for fgsea; default `0.05`.
#' @param time_threshold integer, *serial_time* days to filter *clinical_data* based on; default `NULL`.
#' @param exclude_communities character vector, a vector matching the names of gene_signatures, to exclude specific groups from KM; default `NULL`.
#' @param enrichment_result fgsea result; default `NULL`.
#'
#' @return A list of five elements:
#' - *enrichment_results* fgsea enrichment result
#' - *hazard_ratio* hazard ration from survival analysis
#' - *p_value* chi-squared test pvalue
#' - *adjusted_p_value* chi-squared test BH adjusted pvalue
#' - *plot* Kaplan-Meier plot
#'
#' @export
#'
#' @examples
#' #Read signatures from Signatures folder 
#' test_signature_on_dataset(gene_signatures = gene_signatures,
#'                          transcriptomic_data = transc,
#'                          clinical_data = clinical,
#'                          NES_filter = TRUE,
#'                          pval_cutoff = 0.05,
#'                          palette = 'npg',
#'                          time_threshold = NULL,
#'                          exclude_communities = NULL,
#'                          enrichment_result = NULL)
                  
test_signature_on_dataset <- function(gene_signatures,
                                      transcriptomic_data,
                                      clinical_data,
                                      NES_filter = TRUE,
                                      pval_cutoff = 0.05,
                                      palette = 'npg',
                                      time_threshold = NULL,
                                      exclude_communities = NULL,
                                      enrichment_result = NULL) {
  set.seed(1234)
  
  check_signature_inputs(gene_signatures, transcriptomic_data, 
                         clinical_data, enrichment_result)
  
  if(!is.null(enrichment_result)){
    enrichment_sign <- enrichment_result
  }else{
    # Aggregate once
    gene_mat <- stats::aggregate(transcriptomic_data[,-1],
                                 by = list(gene_name = transcriptomic_data$gene_name),
                                 FUN = mean)
    gene_ids <- gene_mat$gene_name
    gene_mat <- gene_mat[,-1]
    rownames(gene_mat) <- gene_ids
    
    # Prepare fgsea results list
    enrichment_list <- vector("list", ncol(gene_mat))
    colnames_mat <- colnames(gene_mat)
    
    message('Running enrichment analysis')
    for (i in seq_along(enrichment_list)) {
      stats <- gene_mat[, i]
      names(stats) <- rownames(gene_mat)
      stats <- stats[!is.na(stats)]
      stats <- sort(stats, decreasing = TRUE)
      
      fgsea_res <- fgsea::fgsea(pathways = gene_signatures,
                                stats = stats,
                                scoreType = "std") # "pos"
      if (nrow(fgsea_res) > 0) {
        fgsea_res$patient_id <- colnames_mat[i]
        enrichment_list[[i]] <- fgsea_res
      }
    }
    
    enrichment_df <- do.call(rbind, enrichment_list)
    enrichment_df <- enrichment_df[, c("pathway", "NES", "padj", "patient_id")]
    
    # Best signature per patient
    enrichment_sign <- enrichment_df %>%
      dplyr::filter(!NES_filter | NES > 0) %>%
      dplyr::group_by(patient_id) %>%
      dplyr::slice_min(order_by = padj, n = 1) %>%
      dplyr::slice_max(order_by = NES, n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::filter(padj < pval_cutoff)
  }
  
  if(length(unique(enrichment_sign$pathway)) == 1){
    stop('Only one signature was enriched, try different parameters for fgsea')
  }
  
  if(!is.null(clinical_data)){
    # Merge enrichment with clinical info
    merged_df <- merge(enrichment_sign, clinical_data, by = 'patient_id') %>%
      dplyr::rename('group' = 'pathway')
    
    if(!is.null(time_threshold)){
      merged_df <- merged_df %>%
        dplyr::filter(serial_time <= time_threshold)
    }
    
    # Exclude unwanted communities
    if (length(exclude_communities) > 0) {
      merged_df <- merged_df[!merged_df$group %in% exclude_communities, ]
    }
    
    message('Running survival analysis...')
    # Survival analysis
    surv_obj <- survival::Surv(time = as.numeric(merged_df$serial_time), event = merged_df$status)
    surv_test <- survival::survdiff(surv_obj ~ group, data = merged_df)
    
    obs <- surv_test$obs
    exp <- surv_test$exp
    HR <- if (length(obs) >= 2) (obs[1] / obs[2]) / (exp[1] / exp[2]) else NA
    
    pval <- stats::pchisq(surv_test$chisq, df = length(surv_test$n) - 1, lower.tail = FALSE)
    pval_adj <- stats::p.adjust(pval, method = "BH")
    
    # Plot
    plot <- survminer::ggsurvplot(
      survival::survfit(survival::Surv(as.numeric(serial_time), status) ~ group, data = merged_df),
      data = merged_df,
      risk.table = TRUE,
      pval = TRUE,
      palette = palette,
      legend = "right",
      ggtheme = ggplot2::theme_minimal()
    )
    
    return(list(
      enrichment_results = enrichment_sign,
      hazard_ratio = HR,
      p_value = pval,
      adjusted_p_value = pval_adj,
      plot = plot
    ))
  }else{
    return(list(
      enrichment_results = enrichment_sign,
      hazard_ratio = NULL,
      p_value = NULL,
      adjusted_p_value = NULL,
      plot = NULL
    ))
  }
}
