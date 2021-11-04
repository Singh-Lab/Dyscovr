############################################################
### Perform Gene-Set Enrichment Analysis (GSEA)
### Written By: Sara Geraghty, October 2021
############################################################

# This file performs gene set enrichment analysis using the ReactomePA R 
# package (adapted from Josh Wetzel's code)
# ReactomePA vignette: https://yulab-smu.top/biomedical-knowledge-mining-book/

BiocManager::install("ReactomePA")
BiocManager::install("clusterProfiler")
BiocManager::install("DOSE")
BiocManager::install("enrichplot")

library(data.table)
library('ReactomePA')
library('clusterProfiler')
library("DOSE")
library("enrichplot")
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)


output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/BRCA/GSEA/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/Pan-Cancer/GSEA/"

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
  #regprot.gsea.go <- list()
  #regprot.gsea.kegg <- list()
  
  # Loop through all the regulatory proteins
  for (regprot in unique(results_table$R_i.name)) {
    
    # Get all of this protein's targets and their associated Betas, sorted in 
    # descending order
    res_table_sub <- results_table[results_table$R_i.name == regprot, c('T_k.name', 'estimate')]
    
    mapping <- as.data.frame(bitr(res_table_sub$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb="org.Hs.eg.db", drop = FALSE))
    colnames(mapping) <- c("T_k.name", "T_k.entrez")
    res_table_sub <- merge(res_table_sub, mapping, all=TRUE, by="T_k.name")
    res_table_sub <- res_table_sub[order(res_table_sub$estimate, decreasing = TRUE),]
    print(res_table_sub)
    #res_table_sub <- res_table_sub[!(is.na(res_table_sub))]
    #print(res_table_sub)
    
    plotdir <- paste0(output_path, "Linear Model/GSEA/")    
    suppressWarnings(dir.create(plotdir, recursive = TRUE))

    regprotBetaScores <- res_table_sub$estimate
    names(regprotBetaScores) <- res_table_sub$T_k.entrez

    gse.rp <- gsePathway(regprotBetaScores)
    gse.rp <- setReadable(gse.rp, 'org.Hs.eg.db')
    gse.ncg <- gseNCG(regprotBetaScores)
    gse.ncg <- setReadable(gse.ncg, 'org.Hs.eg.db')
    #gse.go <- gseGO(regprotBetaScores, OrgDb = 'org.Hs.eg.db')
    #gse.go <- setReadable(gse.go, 'org.Hs.eg.db')
    #gse.kegg <- gseKEGG(regprotBetaScores)
    #gse.kegg <- setReadable(gse.kegg, 'org.Hs.eg.db')
    
    regprot.gsea.rp[[regprot]] <- gse.rp
    regprot.gsea.ncg[[regprot]] <- gse.ncg
    #regprot.gsea.go[[regprot]] <- gse.go
    #regprot.gsea.kegg[[regprot]] <- gse.kegg
    
    # plot GSEA for top N upregulated and top N downregulated terms
    #for (gst in c('rp','ncg', 'go', 'kegg')) {
    for (gst in c('rp','ncg')) {
      if (gst == 'rp') res <- gse.rp
      #else if (gst == 'kegg') res <- gse.kegg
      #else if (gst == 'go') res <- gse.go
      else res <- gse.ncg
        
      terms <- res@result$Description[res@result$enrichmentScore > 0][1:n]
      terms <- terms[!is.na(terms)]
      if (length(terms) > 0) {
        termIDs <- sapply(terms, function(x) {res@result$ID[which(res@result$Description == x)]})
        pdf(file = paste0(plotdir, regprot,'_gse_',gst,'_gseaCurve_top', n, 'terms_up.pdf'),
            width = 15, height = 10)
        plot(gseaplot2(res, geneSetID = termIDs, pvalue_table = TRUE,
                       color = RColorBrewer::brewer.pal(length(termIDs), "Dark2"), ES_geom = "dot"))
        dev.off()
      }
      terms <- res@result$Description[res@result$enrichmentScore < 0][1:n]
      terms <- terms[!is.na(terms)]
      if (length(terms) > 0) {
        termIDs <- sapply(terms, function(x) {res@result$ID[which(res@result$Description == x)]})
        pdf(file = paste0(plotdir, regprot,'_gse_',gst,'_gseaCurve_top', n, 'terms_down.pdf'),
            width =15, height = 10)
        plot(gseaplot2(res, geneSetID = termIDs, pvalue_table = TRUE,
                       color = RColorBrewer::brewer.pal(length(termIDs), "Dark2"), ES_geom = "dot"))
        dev.off()
      }
    }
  }
}

perform_gsea(master_df_mut_corrected, n = 5, all_genes_id_conv, output_path)
perform_gsea(master_df_cna_corrected, n = 5, all_genes_id_conv, output_path)





#####################################
# ORIGINAL CODE

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



