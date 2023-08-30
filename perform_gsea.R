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
BiocManager::install('fgsea')

# If needed:
# Check libPaths using .libPaths()
# Manually set library using "lib = "C:/Users/sarae/AppData/Local/R/win-library/4.3"

library(data.table)
library(ggnewscale)
library("gage")
library('ReactomePA')
library('clusterProfiler')
library("DOSE")
library("enrichplot")
library('fgsea')
library(BiocParallel)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

# Run when "Error in summary.connection(connection) : invalid connection"
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

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
#' @param output_path a local path to save the GSEA figures to
perform_gsea <- function(results_table, n, sort_by, output_path) {
  
  regprot.gsea.rp <- list()
  regprot.gsea.ncg <- list()
  regprot.gsea.go <- list()
  #regprot.gsea.kegg <- list()
  
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
        return((-log(pval)) * est_sign)
      }))
      res_table_sub <- res_table_sub[order(res_table_sub$negLogPval, decreasing = TRUE),]
    }
    #res_table_sub <- res_table_sub[!(is.na(res_table_sub))]
    #print(res_table_sub)
    
    plotdir <- paste0(output_path, "Linear Model/GSEA/")    
    suppressWarnings(dir.create(plotdir, recursive = TRUE))
    
    #na_ind <- c()
    regprotBetaScores <- res_table_sub$estimate
    if(sort_by == "p.value") {
      regprotBetaScores <- res_table_sub$negLogPval
      #na_ind <- which(!is.finite(regprotBetaScores))
      #regprotBetaScores[na_ind] <- NA
    }
    #regprotBetaScores <- regprotBetaScores[-na_ind]
    #names(regprotBetaScores) <- res_table_sub$T_k.entrez[-na_ind]
    print(head(regprotBetaScores))

    gse.rp <- gsePathway(regprotBetaScores, pvalueCutoff = 1, pAdjustMethod = "BH")
    gse.rp <- setReadable(gse.rp, org.Hs.eg.db)
    print("hi")
    gse.ncg <- gseNCG(regprotBetaScores, pvalueCutoff = 1, pAdjustMethod = "BH")
    gse.ncg <- setReadable(gse.ncg, org.Hs.eg.db)
    print("howdy")
    gse.go <- gseGO(regprotBetaScores, pvalueCutoff = 1, 
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

enriched_terms <- perform_gsea(master_df_mut_corrected, n = 20, "p.value", output_path)
enriched_terms <- perform_gsea(master_df_cna_corrected, n = 20, "p.value", output_path)


enriched_terms <- enriched_terms@result

# Separate into up and down regulated pathways and calculate Jaccard coefficient
tp53_enriched_terms_res <- list()
for (d in c('up','down')) {
  if (d == 'up') {
    flag <-  tp53_enriched_terms@result$enrichmentScore > 0
  } else {
    flag <-  tp53_enriched_terms@result$enrichmentScore < 0
  }
  if (any(flag)) {
    tp53_enriched_terms_res[[d]] <- tp53_enriched_terms
    tp53_enriched_terms_res[[d]]@result <- tp53_enriched_terms_res[[d]]@result[which(flag),]
    tp53_enriched_terms_res[[d]] <- pairwise_termsim(tp53_enriched_terms_res[[d]]) # Computed matrix of pathway similarity (default is Jaccard coef)
  }
}
MIN_JACCARD <- 0.1

# Can look at other parameters in emapplot to make the plot more readable
p <- emapplot(tp53_enriched_terms_res[["up"]], edge.params = list(min = MIN_JACCARD),
              cex.params = list(category_label = 0.6), color = "p.adjust")
p <- emapplot(tp53_enriched_terms_res[["down"]], edge.params = list(min = MIN_JACCARD),
              cex.params = list(category_label = 0.6), color = "p.adjust")

# Look further down the list
tp53_enriched_terms_res[["up2"]] <- tp53_enriched_terms_res[["up"]]
tp53_enriched_terms_res[["up2"]]@result <- tp53_enriched_terms_res[["up2"]]@result[30:nrow(tp53_enriched_terms_res[["up"]]@result),]
emapplot(tp53_enriched_terms_res[["up2"]], edge.params = list(min = MIN_JACCARD),
         cex.params = list(category_label = 0.6), color = "p.adjust")

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
  sig_hits <- unlist(master_df_driver[master_df_driver$q.value < 0.1,])
  foraRes <- fora(react_pathways, genes=sig_hits, universe=names(input))
  collapsedPathways <- collapsePathwaysORA(foraRes[order(pval)][padj < pthres],
                                           react_pathways,
                                           genes=sig_hits,
                                           universe=names(input), pval.threshold = 0.05)

  mainPathways <- foraRes[pathway %in% collapsedPathways$mainPathways][order(pval), pathway]
  print(mainPathways)
  
  # Plot enrichment
  plotGseaTable(react_pathways[mainPathways], input, fgseaRes_sub)
  
  return(list(fgseaRes, mainPathways))
}

fgseaRes_pik3ca <- perform_fgsea(pc_allGenes, "PIK3CA", 0.05)
fgseaRes_pik3ca_df <- fgseaRes_pik3ca[[1]]
fgseaRes_pik3ca_mainPathways <- fgseaRes_pik3ca[[2]]


# Create a heat map with the results for significant, overlapping pathways 
#' @param list_of_results_gsea a named list of GSEA results from fgsea
#' @param sig_thres an adjusted p-value significance threshold
create_sig_pathway_overlap_grid_fgsea <- function(list_of_results_gsea, sig_thres) {
  # Get all the pathways significant in at least one driver
  all_sig_pws <- unique(unlist(lapply(list_of_results_gsea, function(df) df[df$padj < sig_thres, 'pathway'])))
  
  # Set up the driver x pathway matrix that we will fill and use to plot
  input_df <- data.frame(matrix(nrow = length(all_sig_pws), ncol = length(list_of_results_gsea)))
  rownames(input_df) <- all_sig_pws
  colnames(input_df) <- names(list_of_results_gsea)
  
  for(i in 1:length(list_of_results_gsea)) {
    driver <- names(list_of_results_gsea)[i]
    gsea_res <- list_of_results_gsea[[i]]
    sig_hits <- gsea_res[gsea_res$padj < sig_thres, 'pathway']
    input_df[,driver] <- unlist(lapply(all_sig_pws, function(pw) 
      ifelse(pw %fin% sig_hits, 1, 0)))  # could change this to return the actual directional pval
  }
  
  # Limit to just pathways with a 1 in more than one driver
  print(input_df)
  input_df <- input_df[which(rowSums(input_df) > 1),]
  #input_df <- apply(input_df, MARGIN = 2, function(y) as.integer(y))
  print(input_df)
  
  # Make heatmap
  pheatmap(input_df, angle_col = "45", main = "",
           color=c("white", "#0072B5FF"), fontsize_row = 9, fontsize_col = 14,
           fontface_col = "bold", legend = FALSE, cellwidth=20, cellheight = 10) 
  
}

create_sig_pathway_overlap_grid(list("TP53" = results_gsea_tp53, "PIK3CA" = results_gsea_pik3ca, 
                                     "KRAS" = results_gsea_kras, "IDH1" = results_gsea_idh1), 0.2)

# Alternatively, make a heatmap just from the main pathways
# Make heatmap
unique_main_pathways <- unique(c(fgseaRes_tp53_mainPathways, fgseaRes_pik3ca_mainPathways,
                                 fgseaRes_kras_mainPathways, fgseaRes_idh1_mainPathways))
driver_pathways <- list("TP53" = fgseaRes_tp53_mainPathways, "PIK3CA" = fgseaRes_pik3ca_mainPathways, 
                        "KRAS" = fgseaRes_kras_mainPathways, "IDH1" = fgseaRes_idh1_mainPathways)
input_df <- data.frame(matrix(nrow = length(unique_main_pathways), ncol = length(driver_pathways)))
rownames(input_df) <- unique_main_pathways
colnames(input_df) <- names(driver_pathways)
for (i in 1:ncol(input_df)) {
  driver <- colnames(input_df)[i]
  binary_vals <- unlist(lapply(rownames(input_df), function(pw) ifelse(pw %fin% driver_pathways[[i]], 1, 0)))
  input_df[,i] <- binary_vals
}
print(input_df)
pheatmap(input_df, angle_col = "45", main = "",
         color=c("white", "#0072B5FF"), fontsize_row = 9, fontsize_col = 14,
         fontface_col = "bold", legend = FALSE, cellwidth=20, cellheight = 10) 

#####################################
# JOSH ORIGINAL CODE

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
