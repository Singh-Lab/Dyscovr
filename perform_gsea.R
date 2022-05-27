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
BiocManager::install(c("gage", "gageData"))

library(data.table)
library("gage")
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
#' @param results_table a results data frame from my model with Beta estimates + pvalues
#' @param n the number of up/downregulated terms to plot
#' @param sort_by either "estimate" or "p.value" to indicate whether we want to 
#' sort our top hits by Beta estimate or p-value when performing GSEA
#' @param all_genes_id_conv a gene ID conversion data frame from BioMart
#' @param output_path a local path to save the GSEA figures to
perform_gsea <- function(results_table, n, sort_by, output_path) {
  
  regprot.gsea.rp <- list()
  regprot.gsea.ncg <- list()
  regprot.gsea.go <- list()
  regprot.gsea.kegg <- list()
  
  # Loop through all the regulatory proteins
  for (regprot in unique(results_table$R_i.name)) {
    
    # Get all of this protein's targets and their associated Betas, sorted in 
    # descending order
    res_table_sub <- results_table[results_table$R_i.name == regprot, c('T_k.name', 'estimate', 'p.value')]
    print(head(bitr(res_table_sub$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb=org.Hs.eg.db, drop = TRUE)))
    mapping <- as.data.frame(bitr(res_table_sub$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb=org.Hs.eg.db, drop = TRUE))
    colnames(mapping) <- c("T_k.name", "T_k.entrez")
    res_table_sub <- merge(res_table_sub, mapping, all=TRUE, by="T_k.name")
    if(sort_by == "estimate") {
      res_table_sub <- res_table_sub[order(res_table_sub$estimate, decreasing = TRUE),]
    } else {
      res_table_sub$negLogPval <- unlist(lapply(1:nrow(res_table_sub), function(i) {
        pval <- res_table_sub$p.value[i]
        estimate <- res_table_sub$estimate[i]
        est_sign <- ifelse(estimate > 0, 1, -1)
        return(-log(pval) * est_sign)
      }))
      res_table_sub <- res_table_sub[order(res_table_sub$negLogPval, decreasing = TRUE),]
    }
    #res_table_sub <- res_table_sub[!(is.na(res_table_sub))]
    #print(res_table_sub)
    
    plotdir <- paste0(output_path, "Linear Model/GSEA/")    
    suppressWarnings(dir.create(plotdir, recursive = TRUE))
    
    regprotBetaScores <- res_table_sub$estimate
    if(sort_by == "p.value") regprotBetaScores <- res_table_sub$negLogPval
    names(regprotBetaScores) <- res_table_sub$T_k.entrez

    gse.rp <- gsePathway(regprotBetaScores, pvalueCutoff = 1, pAdjustMethod = "BH")
    gse.rp <- setReadable(gse.rp, org.Hs.eg.db)
    gse.ncg <- gseNCG(regprotBetaScores, pvalueCutoff = 1, pAdjustMethod = "BH")
    gse.ncg <- setReadable(gse.ncg, org.Hs.eg.db)
    gse.go <- gseGO(regprotBetaScores, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, 
                    pAdjustMethod = "BH", keyType = "ENTREZID")
    gse.go <- setReadable(gse.go, org.Hs.eg.db)
    #gse.kegg <- gseKEGG(regprotBetaScores, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, 
                        #pAdjustMethod = "BH", keyType = "ENTREZID")
    #gse.kegg <- setReadable(gse.kegg, org.Hs.eg.db)
    
    regprot.gsea.rp[[regprot]] <- gse.rp
    regprot.gsea.ncg[[regprot]] <- gse.ncg
    regprot.gsea.go[[regprot]] <- gse.go
    #regprot.gsea.kegg[[regprot]] <- gse.kegg
    res_tab <- data.table()
    
    
    # plot GSEA for top N upregulated and top N downregulated terms
    for (gst in c('rp','ncg', 'go')) {
    #for (gst in c('rp','ncg')) {
      if (gst == 'rp') res <- gse.rp
      #else if (gst == 'kegg') res <- gse.kegg
      else if (gst == 'go') res <- gse.go
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
            width = 15, height = 10)
        plot(gseaplot2(res, geneSetID = termIDs, pvalue_table = TRUE,
                       color = RColorBrewer::brewer.pal(length(termIDs), "Dark2"), ES_geom = "dot"))
        dev.off()
      }
      if(gst == "go") {res_tab <- res}
    }
  }
  return(res_tab)
}

enriched_terms <- perform_gsea(master_df_mut_corrected, n = 5, "estimate", output_path)
enriched_terms <- perform_gsea(master_df_cna_corrected, n = 5, "estimate", output_path)


enriched_terms <- enriched_terms@result


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


#' Takes in a results table from the linear model output and plots gene set enrichment 
#' using the 'gage' method and package
#' @param results_table a results data frame from my model with Beta estimates + pvalues
#' @param sort_by either "estimate" or "p.value" to indicate whether we want to 
#' sort our top hits by Beta estimate or p-value when performing GSEA
#' @param output_path a local path to save the GSEA figures to
perform_gage <- function(results_table, sort_by, all_genes_id_conv, output_path) {
  
  # Create an input matrix that matches gage's input type (rows are genes, columns are samples)
  # In our case, target genes will be rows, and there will be one column for Beta estimate or q-value
  if(sort_by == "estimate") {
    input_table <- results_table[, c("T_k.name", "estimate")]
    input_table <- input_table[order(input_table$estimate, decreasing = TRUE),]
  } else {
    input_table <- results_table[, c("T_k.name", "estimate", "p.value")]
    input_table$negLogPval <- unlist(lapply(1:nrow(input_table), function(i) {
      pval <- input_table$p.value[i]
      estimate <- input_table$estimate[i]
      est_sign <- ifelse(estimate > 0, 1, -1)
      return(-log(pval) * est_sign)
    }))
    input_table <- input_table[order(input_table$negLogPval, decreasing = TRUE),]
  }
  input_table <- input_table[!(duplicated(input_table$T_k.name)) & !(is.na(input_table$T_k.name)) &
                               !(input_table$T_k.name == ""),]
  print(head(input_table))
  
  # Convert gene ids to entrez IDs
  mapping <- as.data.frame(bitr(input_table$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                                OrgDb=org.Hs.eg.db, drop = FALSE))
  colnames(mapping) <- c("T_k.name", "T_k.entrez")
  input_table <- merge(input_table, mapping, all=TRUE, by="T_k.name")
  input_table <- input_table[!(duplicated(input_table$T_k.entrez)) & !(input_table$T_k.entrez == "") &
                               !(is.na(input_table$T_k.entrez)),]
  print(head(input_table))
  rownames(input_table) <- input_table$T_k.entrez
  if(sort_by == "estimate") {input_table <- input_table[, "estimate", drop = FALSE]}
  else {input_table <- input_table[, "negLogPval", drop = FALSE]}
  print(head(input_table))
  
  # Option: Replace the actual values with ranks
  #input_ranks <- rank(input_table[,1])
  
  # Import KEGG and GO data
  data("kegg.gs")
  data("go.gs")
  kg.hsa=kegg.gsets("hsa")
  kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx]
  
  # Call gage
  kegg_res <- gage(input_table, gsets = kegg.gs, ref = NULL, samp = NULL, rank.test = TRUE)
  go_res <- gage(input_table, gsets = go.gs, ref = NULL, samp = NULL, rank.test = TRUE)
  
  # Print the top results
  print("KEGG, Upregulated:")
  print(head(kegg_res$greater[,1:4]), 4)
  print("KEGG, Downregulated:")
  print(head(kegg_res$less[,1:4]), 4)
  
  print("GO, Upregulated:")
  print(head(go_res$greater[,1:4]), 4)
  print("GO, Downregulated:")
  print(head(go_res$less[,1:4]), 4)
  
  # Return both sets of results
  return(list("KEGG" = kegg_res, "GO" = go_res))
}

gage_results <- perform_gage(allgenes_p53, "estimate", output_path)
gage_results_kegg <- gage_results[[1]]
gage_results_go <- gage_results[[2]]

gage_results_pval <- perform_gage(allgenes_p53, "p.value", output_path)
gage_results_pval_kegg <- gage_results_pval[[1]]
gage_results_pval_go <- gage_results_pval[[2]]
