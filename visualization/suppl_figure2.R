############################################################
# Code to create Suppl. Figure 2 Visualizations
# Written by Sara Geraghty
# PUBLICATION INFORMATION
############################################################

library(data.table)
library(ggplot2)
library(ReactomePA)
library("clusterProfiler")
library("DOSE")
library(GOSemSim)
library("enrichplot")

# Local PATH to directory containing Dyscovr output files
PATH <-  paste0(getwd(), "Output/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
### IMPORT PAN-CANCER OUTPUT FILE(S)
############################################################
outfn <- "res_top_0.05_allGenes_quantile_rawCNA_methMRaw_3PCs_Nonsyn.Drivers.Vogel.elim.vif5.sp0.7_corrected_MUT.csv"
pc_allGenes <- read.csv(paste0(PATH, outfn))

outfn_wttn <- "res_top_0.05_allGenes_quantile_rawCNA_methMRaw_3PCs_Nonsyn.All.Vogel.elim.vif5.sp0.7_corrected_MUT.csv"
pc_allGenes_wtnn <- read.csv(paste0(PATH, outfn_wtnn))


############################################################
### PART A: GENE SET ENRICHMENT 
############################################################
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
# Call if needed

# DOCUMENTATION FOR CLUSTER PROFILER: https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html

#' Takes in a results table from the Dyscovr and runs GSEA using the ReactomePA 
#' package, with options for sorting
#' Defaults to using all gene sets (Reactome, NCG, GO, KEGG, MKEGG, and WikiPathways),
#' but will exhibit speedups if undesired gene sets are commented out
#' @param results_table a results data frame from Dyscovr with mutation coefficient 
#' estimates + q-values
#' @param goi a driver gene name of interest
#' @param sort_by either "estimate" or "p.value" to indicate whether we want to 
#' sort our top hits by Beta estimate or p-value when performing GSEA
perform_gsea <- function(results_table, goi, sort_by) {

  # Get all of this driver protein's targets and their associated mutation 
  # coefficients, sorted descending order
  res_table_sub <- results_table[results_table$R_i.name == goi, 
                                 c('T_k', 'T_k.name', 'estimate', 'p.value')]
  mapping <- as.data.frame(bitr(res_table_sub$T_k.name, fromType = "SYMBOL", 
                                toType = "ENTREZID", OrgDb=org.Hs.eg.db, drop = T))
  mapping_kegg <- as.data.frame(bitr_kegg(res_table_sub$T_k, fromType = "uniprot", 
                                          toType = "kegg", drop = T, 
                                          organism = "hsa"))
  colnames(mapping) <- c("T_k.name", "T_k.entrez")
  colnames(mapping_kegg) <- c("T_k", "T_k.kegg")
  
  res_table_sub <- merge(res_table_sub, mapping, all=T, by="T_k.name")
  res_table_sub_kegg <- merge(res_table_sub, mapping_kegg, all=T, by="T_k")
  
  if(sort_by == "estimate") {
    res_table_sub <- res_table_sub[order(res_table_sub$estimate, decreasing = T),]
    res_table_sub_kegg <- res_table_sub_kegg[order(res_table_sub_kegg$estimate, 
                                                   decreasing = T),]
  } else {
    res_table_sub$negLogPval <- unlist(lapply(1:nrow(res_table_sub), function(i) {
      pval <- res_table_sub$p.value[i]
      estimate <- res_table_sub$estimate[i]
      est_sign <- ifelse(estimate > 0, 1, -1)
      return((-log(pval)) * est_sign)
    }))
    res_table_sub <- res_table_sub[order(res_table_sub$negLogPval, 
                                         decreasing = T),]
    
    res_table_sub_kegg$negLogPval <- unlist(lapply(1:nrow(res_table_sub_kegg), function(i) {
      pval <- res_table_sub_kegg$p.value[i]
      estimate <- res_table_sub_kegg$estimate[i]
      est_sign <- ifelse(estimate > 0, 1, -1)
      return((-log(pval)) * est_sign)
    }))
    res_table_sub_kegg <- res_table_sub_kegg[order(res_table_sub_kegg$negLogPval, 
                                                   decreasing = T),]
    res_table_sub_kegg <- na.omit(res_table_sub_kegg)
  }
  
  driverBetaScores <- res_table_sub$estimate
  driverBetaScores_kegg <- res_table_sub_kegg$estimate
  
  if(sort_by == "p.value") {
    driverBetaScores <- res_table_sub$negLogPval
    driverBetaScores_kegg <- res_table_sub_kegg$negLogPval
  }
  names(driverBetaScores) <- res_table_sub$T_k.entrez
  names(driverBetaScores_kegg) <- res_table_sub_kegg$T_k.kegg
  
  gse.rp <- gsePathway(driverBetaScores, pvalueCutoff = 1, pAdjustMethod = "BH")
  gse.rp <- setReadable(gse.rp, org.Hs.eg.db)
  
  gse.ncg <- gseNCG(driverBetaScores, pvalueCutoff = 1, pAdjustMethod = "BH")
  gse.ncg <- setReadable(gse.ncg, org.Hs.eg.db)
  
  gse.go <- gseGO(driverBetaScores, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, 
                  pAdjustMethod = "BH", keyType = "ENTREZID", ont = "BP")  # BP = biological process
  gse.go <- setReadable(gse.go, org.Hs.eg.db)
  
  gse.kegg <- gseKEGG(driverBetaScores_kegg, pvalueCutoff = 1, 
                      pAdjustMethod = "BH", keyType = "kegg")
  
  gse.mkegg <- gseMKEGG(driverBetaScores_kegg, pvalueCutoff = 1, 
                        pAdjustMethod = "BH", keyType = "kegg")
  
  gse.wp <- gseWP(driverBetaScores, organism = "Homo sapiens", pvalueCutoff = 1, 
                  pAdjustMethod = "BH")
  
  return(list("rp" = gse.rp, "ncg" = gse.ncg, "go" = gse.go, "kegg" = gse.kegg,
              "mkegg" = gse.mkegg, "wp" = gse.wp))
}

# Call function
results_gsea_tp53 <- perform_gsea(pc_allGenes[pc_allGenes$R_i.name == "TP53",], 
                             "TP53", "p.value")
results_gsea_tp53 <- results_gsea_tp53@result

results_gsea_pik3ca <- perform_gsea(pc_allGenes[pc_allGenes$R_i.name == "PIK3CA",], 
                                  "PIK3CA", "p.value")
results_gsea_pik3ca <- results_gsea_pik3ca@result

results_gsea_kras <- perform_gsea(pc_allGenes[pc_allGenes$R_i.name == "KRAS",], 
                                    "KRAS", "p.value")
results_gsea_kras <- results_gsea_kras@result

results_gsea_idh1 <- perform_gsea(pc_allGenes[pc_allGenes$R_i.name == "IDH1",], 
                                  "IDH1", "p.value")
results_gsea_idh1 <- results_gsea_kras@result


#' Merge functionally similar pathways using ReactomePA's built-in semantic 
#' similarity analysis
#' Link to documentation: https://yulab-smu.top/biomedical-knowledge-mining-book/GOSemSim.html
#' @param gsea_res a GSEA file from ReactomePA
#' @param sim_thres a threshold for similarity metric
#' @param hsGO a hsGO file for goSim
#' @param measure a method of semantic similarity measurement, either 
#' 'Wang' or 'Jiang' 
#' @param order_by either 'q.value' or 'enrichmentScore', what to order by
merge_pathways_semantic_sim <- function(gsea_res, sim_thres, hsGO, 
                                        measure, order_by) {
  
  # Make sure we are ordered by correct quantity
  if(order_by == "q.value") {
    gsea_res <- gsea_res[order(gsea_res$qvalue, decreasing=F),]
  } else if (order_by == "enrichmentScore") {
    gsea_res <- gsea_res[order(abs(gsea_res$enrichmentScore), decreasing=T),]
  } else {
    print("Only implemented for 'qvalue' and '-log(qval)*dir', and 'enrichmentScore', 
          please try again with one of these order_by quantities.")
  }
  
  rows_to_keep <- c(1)
  go_term_list <- gsea_res[1, 'ID']
  go_term_list_sizes <- gsea_res[1, 'setSize']
  
  for (i in 2:nrow(gsea_res)) {
    curr_go <- gsea_res[i, 'ID']
    size_go <- gsea_res[i, 'setSize']
    similarity_bool <- F
    for (j in 1:length(go_term_list)) {
      sim <- goSim(go_term_list[[j]], curr_go, semData=hsGO, measure=measure) 
      if(is.na(sim)) {sim <- 0}
      if(sim > sim_thres) {
        similarity_bool <- T
        break
      }
    }
    if(!similarity_bool) {
      rows_to_keep <- c(rows_to_keep, i)
      go_term_list <- append(go_term_list, curr_go)
      go_term_list_sizes <- append(go_term_list_sizes, size_go)
    }
  }
  
  return(gsea_res[rows_to_keep, ])
}

# Import hsGO data
hsGO <- godata('org.Hs.eg.db', ont="BP")  # BP = biological process

# Call merging function
results_gsea_tp53[["go"]] <- merge_pathways_semantic_sim(results_gsea_tp53[["go"]],
                                                    0.7, hsGO, "Wang", 'q.value')
results_gsea_pik3ca[["go"]] <- merge_pathways_semantic_sim(results_gsea_pik3ca[["go"]],
                                                    0.7, hsGO, "Wang", 'q.value')
results_gsea_kras[["go"]] <- merge_pathways_semantic_sim(results_gsea_kras[["go"]],
                                                    0.7, hsGO, "Wang", 'q.value')
results_gsea_idh1[["go"]] <- merge_pathways_semantic_sim(results_gsea_idh1[["go"]],
                                                    0.7, hsGO, "Wang", 'q.value')


#' Takes the output from the ReactomePA gene set enrichment analysis (above) and
#' converts it into a bar chart format for visualization, with -log(p-value) on
#' the x-axis (for multiple drivers)
#' @param list_of_results_gsea a list of results tables from GSEA using the 
#' ReactomePA package, with each entry in the list being the results for a given
#' driver gene. List names are driver gene names.
#' @param n the number of top pathways to display
#' @param qval_thres a q-value threshold for a significantly enriched pathway
#' @param sort_by a method by which to sort the pathways; defaults to q-value,
#' but can also specify enrichment score ("ES"), leading edge percentage 
#' ("leading_edge"), or -log10(q-value) times direction ("neglog10(qval)*dir")
#' @param db the name of the database from which the pathways originated 
#' (e.g. "GO")
create_gsea_barchart_multGenes <- function(list_of_results_gsea, n, qval_thres, 
                                           sort_by, db) {
  
  # Subset the results table based on the number of pathways
  input_dfs <- lapply(1:length(list_of_results_gsea), function(i) {
    results_gsea <- list_of_results_gsea[[i]]
    driver <- names(list_of_results_gsea)[i]
    
    if(!is.na(qval_thres)) {results_gsea <- results_gsea[results_gsea$qvalue < 
                                                           qval_thres, ]}
    
    if(sort_by == "ES") {results_gsea <- results_gsea[order(-abs(
      results_gsea$enrichmentScore)),]}
    
    else if (sort_by == "neglog10(qval)*dir") {
      dir <- unlist(lapply(results_gsea$enrichmentScore, function(x) 
        ifelse(x>0, 1, -1)))
      results_gsea$neglog10qvaltimesdir <- unlist(lapply(1:nrow(results_gsea), function(i) {
        qval <- results_gsea$qvalue[i]
        d <- dir[i]
        return((-log10(qval)) * d)
      }))
    }
    else if (sort_by == "leading_edge") {
      results_gsea$leading_edge_signal <- unlist(lapply(results_gsea$leading_edge, function(x) {
        spl <- unlist(strsplit(x, ", ", fixed = T))[3]
        spl <- unlist(strsplit(unlist(strsplit(x, "=", fixed = T))[2], 
                               "%", fixed = T))[1]
        return(as.numeric(spl))
      }))
      results_gsea <- results_gsea[order(-results_gsea$leading_edge_signal),]
    }
    else {
      results_gsea <- results_gsea[order(results_gsea$qvalue),]
      
      # If desired, uncomment to sort pathways within each q-value bucket by 
      # leading edge signal
      #buckets <- lapply(unique(results_gsea$qvalue), function(q) {
      #  sub <- results_gsea[results_gsea$qvalue == q,]
      #  sub$leading_edge_signal <- unlist(lapply(sub$leading_edge, function(x) {
      #    spl <- unlist(strsplit(x, ", ", fixed = T))[3]
      #    spl <- unlist(strsplit(unlist(strsplit(x, "=", fixed = T))[2], 
      #                           "%", fixed = T))[1]
      #    return(as.numeric(spl))
      #  }))
      #  return(sub[order(-sub$leading_edge_signal),])
      #})
      #results_gsea <- do.call(rbind, buckets)
      #print(head(results_gsea))
    }
    
    results_gsea_sub <- results_gsea[1:min(n, nrow(results_gsea)), ]

    # Convert the p-values to -log10(pvalue)
    negLog10_pvalues <- -log10(results_gsea_sub$qvalue)
    
    # Create an input data frame for ggplot
    input_df <- data.frame("Enriched.Pathway" = results_gsea_sub$Description,
                           "Neg.Log10.Pval" = negLog10_pvalues, 
                           "ES" = results_gsea_sub$enrichmentScore,
                           "Driver" = driver)
    return(input_df)
  })
  input_df <- do.call(rbind, input_dfs)
  
  roles <- function(x) sub("[^_]*_","",x) 
  input_df$Enriched.Pathway.Driver <- unlist(lapply(1:nrow(input_df), function(i) 
    paste(input_df[i, 'Driver'], input_df[i,'Enriched.Pathway'], sep = "_")))
  
  p <- ggplot(input_df, aes(x =  Enriched.Pathway.Driver, #reorder(Enriched.Pathway.Driver, abs(ES)), 
                            y = ES, fill = Driver)) +
    geom_col(width = 0.7, color = "black") + coord_flip() + theme_minimal() + 
    theme(legend.position = "none", 
          axis.title = element_text(face = "bold", size = 14), 
          axis.text = element_text(face = "bold", size = 10), 
          strip.text.x = element_text(face="bold", size=14, 
                                      margin = margin(.2, 0, .2, 0, "cm"))) + 
    xlab(paste("Enriched", paste(db, "Pathways"))) + ylab("Enrichment Score") + 
    scale_fill_manual(values = c("#FFDC91FF","#20854EFF","#BC3C29FF","#0072B5FF")) + 
    scale_x_discrete(labels=roles) +
    facet_wrap(~factor(Driver, levels = c("TP53", "PIK3CA", "KRAS", "IDH1")), 
               ncol = 1, scales = "free_y") 

  p
}

list_of_results_gsea <- list("TP53" = results_gsea_tp53[["go"]], 
                             "PIK3CA" = results_gsea_pik3ca[["go"]], 
                             "KRAS" = results_gsea_kras[["go"]], 
                             "IDH1" = results_gsea_idh1[["go"]])

create_gsea_barchart_multGenes(list_of_results_gsea, 8, 0.01, "q.value", "GO")


############################################################
### PART D: BARPLOT OF HITS AT Q < 0.01
############################################################
qval_thres <- 0.2
pc_allGenes_sig <- pc_allGenes_wtnn[pc_allGenes_wtnn$q.value < qval_thres,]
pc_allGenes_sig_freq <- melt(table(pc_allGenes_sig$R_i.name))
colnames(pc_allGenes_sig_freq) <- c("Driver", "Num.Hits")
pc_allGenes_sig_freq <- pc_allGenes_sig_freq[order(pc_allGenes_sig_freq$Num.Hits),]
ggplot(pc_allGenes_sig_freq, aes(y = Num.Hits, fill = Driver, 
                                 x = reorder(Driver, -Num.Hits, mean))) + 
  geom_bar(position = "dodge", width = 0.95, stat = "identity", 
           show.legend = F, color = "black") + 
  scale_fill_manual(values = c("#FFDC91FF", "#20854EFF", "#BC3C29FF", 
                               "#0072B5FF", "gray")) + 
  xlab("Driver") + 
  ylab(paste0("\n", paste0("Number of hits (q < ", paste0(qval_thres, ")")))) +
  theme(axis.text = element_text(face="bold", size = 16), 
        axis.title=element_text(size=18, face="bold"), 
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = 'white'))
