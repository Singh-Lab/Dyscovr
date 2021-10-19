############################################################
### Perform Gene-Set Enrichment Analysis (GSEA)
### Written By: Sara Geraghty, October 2021
############################################################

# This file performs gene set enrichment analysis using the ReactomePA R 
# package (adapted from Josh Wetzel's code)
# ReactomePA vignette: https://yulab-smu.top/biomedical-knowledge-mining-book/

BiocManager::install("ReactomePA")
BiocManager::install("clusterProfiler")

library('ReactomePA')
library('clusterProfiler')
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


#' Takes in a results table from the linear model output and plots GSEA
#' for n upregulated and n downregulated terms using the 
#' ReactomePA package
#' @param results_table a results data frame from my model
#' @param n the number of up/downregulated terms to plot
#' @param all_genes_id_conv a gene ID conversion data frame from BioMart
#' @param output_path a local path to save the GSEA figures to
perform_gsea <- function(results_table, n, all_genes_id_conv, output_path) {
  
  regprot.gsea.rp <- list()
  regprot.gsea.ncg <- list()
  
  # Loop through all the regulatory proteins
  for (regprot in unique(results_table$R_i.name)) {
    
    # Get all of this protein's targets and their associated Betas, sorted in 
    # descending order
    res_table_sub <- results_table[results_table$R_i.name == regprot, c('T_k.name', 'estimate')]
    res_table_sub <- res_table_sub[order(res_table_sub$estimate, decreasing = TRUE),]
    
    res_table_sub$entrez <- bitr(res_table_sub$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                                 OrgDb="org.Hs.eg.db")
    res_table_sub <- res_table_sub[!(is.na(res_table_sub))]
    
    plotdir <- paste0(output_path, "Linear Model/GSEA/")     # TODO: adjust this path
    suppressWarnings(dir.create(plotdir, recursive = TRUE))
    
    regprotBetaScores <- res_table_sub$score; names(regprotBetaScores) <- res_table_sub$entrez
    
    gse.rp <- gsePathway(regprotBetaScores)
    gse.rp <- setReadable(gse.rp, 'org.Hs.eg.db')
    gse.ncg <- gseNCG(regprotBetaScores)
    gse.ncg <- setReadable(gse.ncg, 'org.Hs.eg.db')
    regprot.gsea.rp[[regprot]] <- gse.rp
    regprot.gsea.ncg[[regprot]] <- gse.ncg
    
    # plot GSEA for top 5 upregulated and top 5 downregulated terms
    for (gst in c('rp','ncg')) {
      if (gst == 'rp') res <- gse.rp else res <- gse.ncg
      terms <- res@result$Description[res@result$enrichmentScore > 0][1:n]
      terms <- terms[!is.na(terms)]
      if (length(terms) > 0) {
        termIDs <- sapply(terms, function(x) {res@result$ID[which(res@result$Description == x)]})
        pdf(file = paste0(plotsub,'clust', i,'_',zf,'_gse_',gst,'_gseaCurve_top5terms_up.pdf'),
            width = 15, height = 10)
        plot(gseaplot2(res, geneSetID = termIDs, pvalue_table = TRUE,
                       color = RColorBrewer::brewer.pal(length(termIDs), "Dark2"), ES_geom = "dot"))
        dev.off()
      }
      terms <- res@result$Description[res@result$enrichmentScore < 0][1:n]
      terms <- terms[!is.na(terms)]
      if (length(terms) > 0) {
        termIDs <- sapply(terms, function(x) {res@result$ID[which(res@result$Description == x)]})
        pdf(file = paste0(plotsub,'clust', i,'_',zf,'_gse_',gst,'_gseaCurve_top5terms_down.pdf'),
            width =15, height = 10)
        plot(gseaplot2(res, geneSetID = termIDs, pvalue_table = TRUE,
                       color = RColorBrewer::brewer.pal(length(termIDs), "Dark2"), ES_geom = "dot"))
        dev.off()
      }
    }
  }

if(PERFORM_GSEA) {
  clust.gsea.rp <- list()
  clust.gsea.ncg <- list()
  for (i in sort(unique(exemplars$clustNum))) {
    if (i == 0) next
    zf <- exemplars[clustNum == i]$zf
    x <- betaMat[,zf]
    x <- sort(x, decreasing = TRUE)
    zfBetas <- data.table(symbol = names(x), score = x)
    zfBetas$entrez <- mapIds(org.Hs.eg.db, keys = zfBetas$symbol,
                             keytype = "SYMBOL", column="ENTREZID")
    zfBetas <- zfBetas[!(is.na(zfBetas$entrez))]
    plotsub <- paste0(PLOT_DIR, 'reactome_exemplars/')
    suppressWarnings(dir.create(plotsub, recursive = TRUE))
    zfBetaScores <- zfBetas$score; names(zfBetaScores) <- zfBetas$entrez
    gse.rp <- gsePathway(zfBetaScores)
    gse.rp <- setReadable(gse.rp, 'org.Hs.eg.db')
    gse.ncg <- gseNCG(zfBetaScores)
    gse.ncg <- setReadable(gse.ncg, 'org.Hs.eg.db')
    clust.gsea.rp[[zf]] <- gse.rp
    clust.gsea.ncg[[zf]] <- gse.ncg
    
    # plot GSEA for top 5 upregulated and top 5 downregulated terms
    for (gst in c('rp','ncg')) {
      if (gst == 'rp') res <- gse.rp else res <- gse.ncg
      terms <- res@result$Description[res@result$enrichmentScore > 0][1:5]
      terms <- terms[!is.na(terms)]
      if (length(terms) > 0) {
        termIDs <- sapply(terms, function(x) {res@result$ID[which(res@result$Description == x)]})
        pdf(file = paste0(plotsub,'clust', i,'_',zf,'_gse_',gst,'_gseaCurve_top5terms_up.pdf'),
            width = 15, height = 10)
        plot(gseaplot2(res, geneSetID = termIDs, pvalue_table = TRUE,
                       color = RColorBrewer::brewer.pal(length(termIDs), "Dark2"), ES_geom = "dot"))
        dev.off()
      }
      terms <- res@result$Description[res@result$enrichmentScore < 0][1:5]
      terms <- terms[!is.na(terms)]
      if (length(terms) > 0) {
        termIDs <- sapply(terms, function(x) {res@result$ID[which(res@result$Description == x)]})
        pdf(file = paste0(plotsub,'clust', i,'_',zf,'_gse_',gst,'_gseaCurve_top5terms_down.pdf'),
            width =15, height = 10)
        plot(gseaplot2(res, geneSetID = termIDs, pvalue_table = TRUE,
                       color = RColorBrewer::brewer.pal(length(termIDs), "Dark2"), ES_geom = "dot"))
        dev.off()
      }
    }
  }
}



