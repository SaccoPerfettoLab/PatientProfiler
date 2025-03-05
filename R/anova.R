
#' Perform ANOVA-Tukey Analysis on Transcriptomics Data
#'
#' This function performs an ANOVA-Tukey analysis on transcriptomics data to identify significant gene expression differences across communities.
#'
#' @param Transcriptomics_patients Data frame. Contains transcriptomics data with gene expression values and patient IDs.
#' @param stratification_table Data frame. Contains patient-community assignments.
#'
#' @return A data frame with ANOVA-Tukey results, including adjusted p-values and mean expression levels.
#' @export
#'
#' @examples
#' perform_anova(Transcriptomics_patients, stratification_table)


perform_anova <- function(Transcriptomics_patients, stratification_table) {
  clusters <- unique(stratification_table$community)
  all_results <- list()
  
  for (cluster in clusters) {
    Transcriptomics_patients_cluster <- Transcriptomics_patients
    Transcriptomics_patients_cluster$Cluster_anova <- 'background'
    Transcriptomics_patients_cluster$Cluster_anova[Transcriptomics_patients_cluster$community == cluster] <- cluster
    
    comp_cluster <- Transcriptomics_patients_cluster %>%
      dplyr::group_by(Cluster_anova, gene_name) %>%
      dplyr::summarise(mean = mean(value), .groups = 'drop')
    
    glob_res <- data.frame()
    
    for (gene_i in unique(Transcriptomics_patients$gene_name)) {
      table_summary_gene_i <- Transcriptomics_patients_cluster %>%
        filter(gene_name == gene_i)
      
      clust_value_i <- comp_cluster$mean[comp_cluster$Cluster_anova == cluster & comp_cluster$gene_name == gene_i]  
      
      if (length(unique(table_summary_gene_i$Cluster_anova)) > 1) {
        res.aov.gene.i <- aov(value ~ Cluster_anova, data = table_summary_gene_i)
        res.aov.gene.i.results <- TukeyHSD(res.aov.gene.i)$Cluster_anova
        Comparison <- row.names(res.aov.gene.i.results)
        
        res.aov.gene.i.results <- as_tibble(res.aov.gene.i.results)
        res.aov.gene.i.results$gene_name <- gene_i
        res.aov.gene.i.results$cluster <- cluster
        res.aov.gene.i.results$mean_exp_clus <- clust_value_i
        res.aov.gene.i.results$comparison <- Comparison
        
        glob_res <- rbind(glob_res, res.aov.gene.i.results)
      }
    }
    
    all_results[[cluster]] <- glob_res
  }
  
  final_results <- do.call(rbind, all_results)
  return(final_results)
}

