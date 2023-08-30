### GSEA USING THE FGSEA PACKAGE
# https://rdrr.io/bioc/fgsea/f/vignettes/fgsea-tutorial.Rmd

# Use Reactome pathways

# Use this ranking for fast gsea (fgsea)
#' @param master_df a linear model results data frame with target genes for
#' a given driver ranked by -log(pval) * sign of estimate
#' @param goi the gene name of the driver of interest
#' @param pthres a p-value threshold for significantly enriched pathways
#' @param max_size a maximum size of gene set for GSEA, see fgsea documentation
perform_fgsea <- function(master_df, goi, pthres, max_size = 500) {
  # Subset to a driver of interest
  master_df_driver <- master_df[master_df$R_i.name == goi,]
  
  # Add negative log p-value with directionality
  master_df_driver$negLogPval <- unlist(lapply(1:nrow(master_df_driver), function(i) {
    pval <- master_df_driver$p.value[i]
    estimate <- master_df_driver$estimate[i]
    est_sign <- ifelse(estimate > 0, 1, -1)
    return((-log(pval)) * est_sign)
  }))
  
  # Add entrez IDs
  mapping <- as.data.frame(bitr(master_df_driver$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                                OrgDb=org.Hs.eg.db, drop = TRUE))
  colnames(mapping) <- c("T_k.name", "T_k.entrez")
  master_df_driver <- merge(master_df_driver, mapping, all=TRUE, by="T_k.name")
  
  # Rank targets by negative log p-value
  master_df_driver <- master_df_driver[order(master_df_driver$negLogPval, decreasing = TRUE),]
  
  # Get reactome pathways for the given target genes
  react_pathways <- reactomePathways(master_df_driver$T_k.entrez)
  
  # Get the ranked list
  input <- as.numeric(master_df_driver$negLogPval)
  names(input) <- master_df_driver$T_k.entrez
  input <- input[!is.na(input)]
  print(head(input))
  
  # Perform gsea
  fgseaRes <- fgsea(react_pathways, input, maxSize=max_size)  # Can set eps argument to 0 for a more accurate pval
  fgseaRes <- fgseaRes[order(fgseaRes$padj),]
  
  # Subset to only significant pathways
  fgseaRes_sub <- fgseaRes[fgseaRes$padj < pthres,]
  
  # Collapse pathways
  #collapsed_pathways <- collapsePathways(fgseaRes = fgseaRes_sub, pathways = react_pathways, stats = input)
  foraRes_down <- fora(react_pathways, genes=tail(names(input), 200), universe=names(input))
  foraRes_up <- fora(react_pathways, genes=head(names(input), 200), universe=names(input))
  collapsedPathways_down <- collapsePathwaysORA(foraRes_down[order(pval)][padj < pthres],
                                                react_pathways,
                                                genes=tail(names(input), 200),
                                                universe=names(input), pval.threshold = 0.05)
  collapsedPathways_up <- collapsePathwaysORA(foraRes_up[order(pval)][padj < pthres],
                                              react_pathways,
                                              genes=head(names(input), 200),
                                              universe=names(input), pval.threshold = 0.05)
  mainPathways <- unique(c(foraRes_down[pathway %in% collapsedPathways_down$mainPathways][order(pval), pathway],
                           foraRes_up[pathway %in% collapsedPathways_up$mainPathways][order(pval), pathway]))
  print(mainPathways)
  
  # Plot enrichment
  plotGseaTable(react_pathways[mainPathways], input, fgseaRes_sub, 
                gseaParam = 0.5)
  
  return(fgseaRes)
}

fgseaRes_pik3ca <- perform_fgsea(pc_allGenes, "PIK3CA", 0.05)
