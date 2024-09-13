############################################################
# Code to Create Figure 3 Visualizations
# Written by Sara Geraghty
# PUBLICATION INFORMATION
############################################################

library(data.table)
library(ggplot2)
library(pheatmap)
library("RColorBrewer")
library(ReactomePA)
library("clusterProfiler")
library("DOSE")
library(broom)
library(dplyr)
library(metap)
library(VennDiagram)
library(ggsignif)
library(qvalue)
library(ggforce)
library(graphlayouts)
library(oaqc)
library(concaveman)
library(ggrepel)
library(ggsci)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(gage)

# Local PATH to directory containing Dyscovr output files
PATH <-  paste0(getwd(), "Output/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
### IMPORT PAN- and PER-CANCER OUTPUT FILE(S)
############################################################
# Pan-Cancer
outfn <- "res_top_0.05_allGenes_quantile_rawCNA_methMRaw_3PCs_Nonsyn.Drivers.Vogel.elim.vif5.sp0.7_corrected_MUT.csv"
pc_allGenes <- read.csv(paste0(PATH, outfn), header = T, check.names = F, 
                        row.names = 1)

# Per-Cancer
perCancer_fns <- intersect(list.files(path = PATH, recursive = T,
                                      pattern = "_corrected_MUT"), 
                           intersect(list.files(path = PATH, recursive = T,
                                                pattern = "allGenes_"),
                                     list.files(path = PATH, recursive = T,
                                                pattern = "Nonsyn.Drivers.Vogel.elim.vif.5")))
perCancer <- lapply(perCancer_fns, function(f) 
  fread(paste0(PATH, f), header = T))
names(perCancer) <- unlist(lapply(perCancer_fns, function(x)
  unlist(strsplit(x, "/", fixed = T))[1]))


############################################################
### PART A: PAN-CANCER METABOLIC GSEA HEATMAP
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

# Call function and write KEGG results to file. For metabolic analysis, we'll 
# use KEGG due to the greater prevalance of metabolic pathways.
tp53_pc_gsea <- perform_gsea(pc_allGenes[pc_allGenes$R_i.name == "TP53",], 
                                  "TP53", "p.value")
tp53_pc_gsea_kegg <- tp53_pc_gsea[["kegg"]]@result
write.csv(tp53_pc_gsea_kegg, paste0(PATH, "GSEA_Results/TP53/PC_gsea_kegg_results.csv"))

pik3ca_pc_gsea <- perform_gsea(pc_allGenes[pc_allGenes$R_i.name == "PIK3CA",], 
                                    "PIK3CA", "p.value")
pik3ca_pc_gsea_kegg <- pik3ca_pc_gsea[["kegg"]]@result
write.csv(pik3ca_pc_gsea_kegg, paste0(PATH, "GSEA_Results/PIK3CA/PC_gsea_kegg_results.csv"))

kras_pc_gsea <- perform_gsea(pc_allGenes[pc_allGenes$R_i.name == "KRAS",], 
                                  "KRAS", "p.value")
kras_pc_gsea_kegg <- kras_pc_gsea[["kegg"]]@result
write.csv(kras_pc_gsea_kegg, paste0(PATH, "GSEA_Results/KRAS/PC_gsea_kegg_results.csv"))

idh1_pc_gsea <- perform_gsea(pc_allGenes[pc_allGenes$R_i.name == "IDH1",], 
                                  "IDH1", "p.value")
idh1_pc_gsea_kegg <- idh1_pc_gsea[["kegg"]]@result
write.csv(idh1_pc_gsea_kegg, paste0(PATH, "GSEA_Results/IDH1/PC_gsea_kegg_results.csv"))


# Import as needed
tp53_pc_gsea_kegg <- read.csv(paste0(PATH, "GSEA_Results/TP53/PC_gsea_kegg_results.csv"), 
                              header = T, check.names = F)
pik3ca_pc_gsea_kegg <- read.csv(paste0(PATH, "GSEA_Results/PIK3CA/PC_gsea_kegg_results.csv"), 
                                header = T, check.names = F)
#pik3ca_pc_gsea_kegg[pik3ca_pc_gsea_kegg$Description == "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate", 'Description'] <- "Glycosaminoglycan biosynthesis"
kras_pc_gsea_kegg <- read.csv(paste0(PATH, "GSEA_Results/KRAS/PC_gsea_kegg_results.csv"), 
                              header = T, check.names = F)
idh1_pc_gsea_kegg <- read.csv(paste0(PATH, "GSEA_Results/IDH1/PC_gsea_kegg_results.csv"), 
                              header = T, check.names = F)

#' Compute the Jaccard correlation coefficient of a and b
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#' Merge KEGG pathways in a GSEA result by Jaccard coefficient. Any pathways 
#' exceeding the Jaccard threshold for similarity will be merged, with the 
#' pathway with the more significant value providing the name of the merged pathways
#' @param gsea_res a GSEA results file from ReactomePA
#' @param jaccard_thres a threshold for the Jaccard coefficient, beyond which
#' two pathways should be considered functionally similar and merged
#' @param order_by either 'qvalue', '-log(qval)*dir', or 'enrichmentScore' to
#' denote how GSEA file should be ordered prior to merging
merge_pathways <- function(gsea_res, jaccard_thres, order_by) {
  
  # Make sure we are ordered by correct quantity
  if(order_by == "qvalue") {
    gsea_res <- gsea_res[order(gsea_res$qvalue, decreasing=F),]
  } else if (order_by == "enrichmentScore") {
    gsea_res <- gsea_res[order(abs(gsea_res$enrichmentScore), decreasing=T),]
  } else {
    print("Only implemented for 'qvalue' and '-log(qval)*dir', and 'enrichmentScore', 
          please try again with one of these order_by quantities.")
  }
  
  rows_to_keep <- c(1)
  pathway_member_list <- list(strsplit(gsea_res[1, 'core_enrichment'], 
                                       "/", fixed = T))
  
  for (i in 2:nrow(gsea_res)) {
    curr_members <- strsplit(gsea_res[i, 'core_enrichment'], "/", fixed = T)
    similarity_bool <- F
    for (j in 1:length(pathway_member_list)) {
      jacc <- jaccard(pathway_member_list[[j]], curr_members[[1]]) 
      if(jacc > jaccard_thres) {
        similarity_bool <- T
        break
      }
    }
    if(!similarity_bool) {
      rows_to_keep <- c(rows_to_keep, i)
      pathway_member_list <- append(pathway_member_list, curr_members)
    }
  }
  
  return(gsea_res[rows_to_keep, ])
}

tp53_pc_gsea_kegg_merged_0.1 <- merge_pathways(tp53_pc_gsea_kegg, 0.1, 'qvalue')  # from 337 to 110 pathways
pik3ca_pc_gsea_kegg_merged_0.1 <- merge_pathways(pik3ca_pc_gsea_kegg, 0.1, 'qvalue')  # from 337 to 108 pathways
kras_pc_gsea_kegg_merged_0.1 <- merge_pathways(kras_pc_gsea_kegg, 0.1, 'qvalue') # from 337 to 107 pathways
idh1_pc_gsea_kegg_merged_0.1 <- merge_pathways(idh1_pc_gsea_kegg, 0.1, 'qvalue') # from 337 to 107 pathways


# Create a heat map with the results for significant, overlapping pathways 
#' @param list_of_results_gsea a named list of GSEA results from Reactome PA
#' @param sig_thres a q-value significance threshold
create_sig_pathway_overlap_grid <- function(list_of_results_gsea, sig_thres) {
  
  # Get all the pathways significant in at least one driver
  all_sig_pws <- unique(unlist(lapply(list_of_results_gsea, function(df) 
    df[df$qvalue < sig_thres, 'Description'])))
  
  # Set up the driver x pathway matrix that we will fill and use to plot
  input_df <- data.frame(matrix(nrow = length(all_sig_pws), 
                                ncol = length(list_of_results_gsea)))
  rownames(input_df) <- all_sig_pws
  colnames(input_df) <- names(list_of_results_gsea)
  
  for(i in 1:length(list_of_results_gsea)) {
    driver <- names(list_of_results_gsea)[i]
    gsea_res <- list_of_results_gsea[[i]]
    sig_hits <- gsea_res[gsea_res$qvalue < sig_thres, 'Description']
    
    input_df[,driver] <- unlist(lapply(all_sig_pws, function(pw) {
      if(pw %fin% sig_hits) {return(i)} 
      else{return(0)}
    })) 
  }
  # Limit to just pathways with a 1 in more than one driver
  #input_df <- input_df[which(rowSums(input_df) > 1),]
  #input_df <- input_df[rowSums(input_df > 0) > 1,]
  #input_df <- apply(input_df, MARGIN = 2, function(y) as.integer(y))
  print(input_df)
  
  # Make heatmap
  pheatmap(t(input_df), angle_col = "45", main = "",
           #color=c("#E18727FF", "white", "#0072B5FF"),
           c("white", "#0072B5FF", "#BC3C29FF", "#20854EFF", "#E18727FF"),
           #color=c("-1" = "#0072B5FF", "-2" = "#BC3C29FF", "-3" = "#20854EFF", 
           #"-4" = "#E18727FF", "1" = "#0072B599", "2" = "#BC3C2999", 
           #"3" = "#20854E99", "4" = "#E1872799", "0" = "white"), 
           fontsize_row = 14, fontsize_col = 9, fontface_row = "bold", 
           fontface_col = "bold", legend = F, cellwidth=18, cellheight = 23,
           cluster_rows = F, cluster_cols = F) 
  
  return(input_df)
}

pathway_matrix <- create_sig_pathway_overlap_grid(list("TP53" = tp53_pc_gsea_kegg_merged_0.1, 
                                                       "PIK3CA" = pik3ca_pc_gsea_kegg_merged_0.1, 
                                                       "KRAS" = kras_pc_gsea_kegg_merged_0.1, 
                                                       "IDH1" = idh1_pc_gsea_kegg_merged_0.1), 
                                                  0.05)

# Import the KEGG set of metabolic pathways, downloaded from the KEGG database
kegg_metabolic_pws <- read.table(paste0(PATH, "kegg_pathways_metabolism.keg.txt"), 
                                 header = F, sep = "\t")

# Adjust this file for proper parsing
curr_pw_broad <- ""   # the broad, B level pathway
curr_pw_narrow <- ""  # the narrow, C level pathway
kegg_metabolic_pws_list <- list()
tmp_narrow_pw_list <- list()
tmp_entry_list <- c()

for(i in 1:nrow(kegg_metabolic_pws)) {
  entry <- kegg_metabolic_pws[i, 1]
  entry_spl <- unlist(strsplit(entry, " ", fixed = T))
  entry_spl <- entry_spl[!(entry_spl == "")]
  # Broad pathway label
  if(entry_spl[1] == "B") {
    if(!(curr_pw_broad == "")) {
      kegg_metabolic_pws_list <- append(kegg_metabolic_pws_list, tmp_narrow_pw_list)
      tmp_narrow_pw_list <- list()
    }
    curr_pw_broad <- paste(entry_spl[3:length(entry_spl)], collapse = " ")
    #names(kegg_metabolic_pws_list)[length(kegg_metabolic_pws_list)] <- curr_pw_broad
    print(curr_pw_broad)
    print(kegg_metabolic_pws_list)
    # Narrow pathway label
  } else if (entry_spl[1] == "C") {
    if(!(curr_pw_narrow == "")) {
      print(head(tmp_entry_list))
      tmp_narrow_pw_list <- append(tmp_narrow_pw_list, paste(unique(tmp_entry_list), collapse = ";"))
      tmp_entry_list <- c()
    }
    curr_pw_narrow <- paste(entry_spl[3:length(entry_spl)], collapse = " ")
    print(curr_pw_narrow)
    names(tmp_narrow_pw_list)[length(tmp_narrow_pw_list)] <- curr_pw_narrow
    print(tmp_narrow_pw_list)
    # Gene entry
  } else {
    gene <- unlist(strsplit(unlist(strsplit(entry, ";", fixed = T))[1], " ", fixed = T))
    gene <- gene[!(gene == "")]
    gene <- toupper(gene[length(gene)])
    if(!grepl(".", gene, fixed = T)) {
      if(grepl("_", gene, fixed = T)) {
        gene <- unlist(strsplit(gene, "_", fixed = T))
        gene[2:length(gene)] <- unlist(lapply(gene[2:length(gene)], function(x) {
          core <- gsub("^\\d+|\\d+$", "", gene[1])  
          return(paste0(core, x))
        }))
      }
      tmp_entry_list <- c(tmp_entry_list, gene)
      print(gene)
    }
  }
}

kegg_genes <- unique(unlist(strsplit(paste(as.character(unlist(kegg_metabolic_pws_list)), 
                                           collapse = ";"), ";", fixed = T)))  # 2904
kegg_genes <- kegg_genes[kegg_genes != "METABOLISM"]

kegg_metabolic_pw_names <- unlist(lapply(kegg_metabolic_pws$V1, function(x) {
  x <- unlist(strsplit(x, " [PATH", fixed = T))[1]
  spl <- unlist(strsplit(x, " ", fixed = T))
  spl <- spl[spl != ""]
  if(!(spl[1] %in% c("A", "D"))) {
    spl <- spl[spl != "NA"]
    new <- paste(spl[3:length(spl)], collapse = " ")
    return(new)
  }
}))
kegg_metabolic_pw_names <- kegg_metabolic_pw_names[kegg_metabolic_pw_names != 
                                                     "NA Metabolism"]
write.table(kegg_metabolic_pw_names, paste0(PATH, "kegg_pathways_metabolism.txt"),
            quote = F, row.names = F)
kegg_metabolic_pw_names <- read.table(paste0(PATH, "kegg_pathways_metabolism.txt"), 
                                      header = T, sep = "\n")[,1]

# Call function again with only metabolic KEGG pathways
pathway_matrix_metabol <- create_sig_pathway_overlap_grid(
  list("TP53" = tp53_pc_gsea_kegg[tp53_pc_gsea_kegg$Description %fin% 
                                    kegg_metabolic_pw_names,], 
       "PIK3CA" = pik3ca_pc_gsea_kegg[pik3ca_pc_gsea_kegg$Description %fin%
                                        kegg_metabolic_pw_names,], 
       "KRAS" = kras_pc_gsea_kegg[kras_pc_gsea_kegg$Description %fin%
                                    kegg_metabolic_pw_names,], 
       "IDH1" = idh1_pc_gsea_kegg[idh1_pc_gsea_kegg$Description %fin% 
                                    kegg_metabolic_pw_names,]), 
  0.1)


#' Make a more detailed version containing information about the enrichment score
#' in heatmap form, per driver or for one driver per-cancer type
#' @param list_of_gsea_files a list of GSEA files for each driver or for each 
#' cancer type, named accordingly
#' @param qval_thres a q-value threshold for a significantly enriched pathway
#' @param terms an optional subset of pathways of interest, e.g. KEGG metabolic
#' pathways
create_gsea_heatmap <- function(list_of_gsea_files, qval_thres, terms) {
  
  list.gsea <- list()
  list.gsea.signif <- list()
  
  signif_pws <- c()
  
  for(i in 1:length(list_of_gsea_files)) {
    entry <- names(list_of_gsea_files)[i]
    f <- list_of_gsea_files[[i]]
    
    if(length(terms) > 1) {
      f <- f[f$Description %fin% terms, ]
    }
    
    signif_pws <- c(signif_pws, as.character(unlist(f[f$qvalue < qval_thres, 
                                                      'Description'])))
    
    matrix <- data.frame("description" = f$Description, 
                         "enrichment.score" = f$enrichmentScore)
    colnames(matrix)[2] <- paste0("enrichment.score.", entry)
    dir <- unlist(lapply(matrix$enrichment.score, function(x) 
      ifelse(x > 0, 1, -1)))
    signif_matrix <- data.frame("description" = f$Description, 
                                "qval" = -log10(f$qvalue))
    signif_matrix$negLogqvaltimesDir <- unlist(lapply(1:length(signif_matrix$qval), 
                                                      function(i) {
      q <- signif_matrix$qval[i]
      d <- dir[i]
      return(q * d)
    }))
    colnames(signif_matrix)[2] <- paste0("negLogqvaltimesDir.", entry)
    
    list.gsea[[i]] <- matrix
    list.gsea.signif[[i]] <- signif_matrix
  }
  names(list.gsea) <- names(list_of_gsea_files)
  
  matrix <- Reduce(function(x, y) merge(x, y, by="description"), list.gsea)
  signif_matrix <- Reduce(function(x, y) merge(x, y, by="description"), 
                          list.gsea.signif)
  print(length(signif_pws))
  
  matrix <- matrix[matrix$description %fin% signif_pws,]
  signif_matrix <- signif_matrix[signif_matrix$description %fin% signif_pws,]
  
  rownames(matrix) <- matrix$description
  rownames(signif_matrix) <- signif_matrix$description
  matrix <- matrix[,colnames(matrix) != "description"]
  signif_matrix <- signif_matrix[,colnames(signif_matrix) != "description"]
  #print(head(matrix))
  
  matrix <- t(matrix)
  signif_matrix <- t(signif_matrix)
  #print(head(matrix))
  
  rownames(matrix) <- names(list_of_gsea_files)
  #rownames(signif_matrix) <- names(list_of_gsea_files)
  
  # Fix any excessively long pathway names
  colnames(matrix)[which(grepl("Glycosaminoglycan biosynthesis", colnames(matrix)))] <- 
    "Glycosaminoglycan biosynthesis*"
  
  print(head(signif_matrix))
  
  # Optionally keep only those significant in at least 2 cancer types/ driver genes
  #cols_to_keep <- colSums(is.na(matrix))<nrow(matrix)
  #cols_to_keep <- cols_to_keep & (colSums(signif_matrix < padj_thres) >= 2)
  #matrix <- matrix[,cols_to_keep]
  #signif_matrix <- signif_matrix[,cols_to_keep]
  
  #matrix[is.na(matrix)] <- 0
  
  # Remove rows that are not entirely 0
  rows_to_keep <- !(rowSums(matrix) == 0)
  matrix <- matrix[rows_to_keep,]
  signif_matrix <- signif_matrix[rows_to_keep,]
  #print(head(matrix))
  
  rc1 <- colorRampPalette(c("#0072B5FF", "white"))(50)
  rc2 <- colorRampPalette(c("white", "#BC3C29FF"))(50)           
  rampcols <- c(rc1, rc2)
  
  rb1 <- seq(min(matrix, na.rm = T), 0, length.out = 51)
  rb2 <- seq(0, max(matrix, na.rm = T), length.out = 51)[-1]
  rampbreaks <- c(rb1, rb2)
  hm <- pheatmap(matrix, angle_col = "45", main = "", cellheight=25, 
                 cellwidth=22, color= rampcols, breaks = rampbreaks, bias = 1, 
                 fontsize_row = 12, fontsize_col = 10, na_col = "gray90") 
  
  #hm <- pheatmap(matrix, angle_col = "45", main = "", cellheight=15, cellwidth=12,
  #color=colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(50),
  #fontsize_row = 10, fontsize_col = 8) 
  #ylab = "\n\nCancer Type", xlab = "KEGG Metabolic Pathway") #, cellheight = 25, cellwidth = 25)
  
  print(hm)
}

list_of_gsea_files_metabol <- list(
  "TP53" = tp53_pc_gsea_kegg[tp53_pc_gsea_kegg$Description %fin% 
                               kegg_metabolic_pw_names,], 
  "PIK3CA" = pik3ca_pc_gsea_kegg[pik3ca_pc_gsea_kegg$Description %fin%
                                   kegg_metabolic_pw_names,], 
  "KRAS" = kras_pc_gsea_kegg[kras_pc_gsea_kegg$Description %fin%
                               kegg_metabolic_pw_names,], 
  "IDH1" = idh1_pc_gsea_kegg[idh1_pc_gsea_kegg$Description %fin% 
                               kegg_metabolic_pw_names,])

create_gsea_heatmap(list_of_gsea_files_metabol, 0.05, kegg_metabolic_pw_names)


############################################################
### PART B: PER-CANCER PYRIMIDINE METABOLIC HEATMAP, KRAS, AND INSULIN SIGNALING
### METABOLIC HEATMAP, PIK3CA
############################################################
#' Create a heatmap for a particular genetic pathway and particular driver that
#' hones in on the -log10(qval)*dir for each target gene in that pathway from Dyscovr
#' @param perCancer list of results files from Dyscovr, named according to cancer
#' type
#' @param gene_set a vector of genes in the particular pathway of interest
#' @param driver the name of the driver of interest
#' @param metric_to_display either '-log10(qval)*dir' or 'tstat' to denote
#' which metric will be displayed in the heatmap
create_specific_pathway_hm <- function(perCancer, gene_set, driver, metric_to_display) {
  # Create input list for heatmap
  matrix_list <- lapply(perCancer, function(x) {
    if(!(driver %fin% x$R_i.name)) {return(NA)}
    
    vals <- unlist(lapply(gene_set, function(g) {
      if(!(g %fin% x$T_k.name)){return(NA)}
      if(metric_to_display == "tstat") {
        tstat <- as.numeric(x[(x$T_k.name == g) & (x$R_i.name == driver), 
                              'statistic'])
        return(tstat)
      } else {
        est <- as.numeric(x[(x$T_k.name == g) & (x$R_i.name == driver), 
                            'estimate'])
        negLogQval <- -log10(as.numeric(x[(x$T_k.name == g) & 
                                            (x$R_i.name == driver), 'q.value']))
        if(!is.na(est)) {
          if(est > 0) {return(1*negLogQval)}
          else {return(-1*negLogQval)}
        } else {return(NA)}
      }
    }))
    return(vals) 
  })
  matrix_list <- matrix_list[!is.na(matrix_list)]

  # Bind these together across the different cancer types
  matrix <- do.call(rbind, matrix_list)
  colnames(matrix) <- gene_set

  # Remove columns that are at least 50% non-zero
  matrix <- matrix[,colSums(is.na(matrix)) < (nrow(matrix)*0.5)]
  print(head(matrix))
  
  # Optionally restrict to those rows and columns that have sufficiently high
  # values (most significant)
  #matrix_sub <- matrix[, colSums(abs(matrix) > 2, na.rm = T) >= 1]
  #matrix_sub <- matrix_sub[rowSums(abs(matrix_sub) > 2, na.rm = T) >= 1, ]
  #print(head(matrix_sub))
  
  # Get the colors to be centered on 0
  rc1 <- colorRampPalette(c("#0072B5FF", "white"))(50)
  rc2 <- colorRampPalette(c("white", "#BC3C29FF"))(50)           
  rampcols <- c(rc1, rc2)
  
  rb1 <- seq(-max(matrix, na.rm = T), 0, length.out = 51)
  rb2 <- seq(0, max(matrix, na.rm = T), length.out = 51)[-1]
  rampbreaks <- c(rb1, rb2)
  p <- pheatmap(matrix, angle_col = "45", main = "", color= rampcols, 
           breaks = rampbreaks, bias = 1, na_col = "grey90") 
  print(p)
}

# Import pyrimidine metabolism genes from various sources
## CURATED ##
pyrimidine_pw_genes <- c("AICDA", "APOBEC2", "CDA", "CDADC1", "DCK", "DPYD", 
                         "DTYMK", "PYCRL", "TK1", "TK2", "TYMP", "UPB1", "UPP1", 
                         "UPP2")
## KEGG ##
#pyrimidine_pw_genes_kegg <- c("NT5E", "CANT1", "TK2", "ENTPD4", "CTPS2", "NME4", 
#                              "NME2", "NME1", "DCTD", "ENPP3", "NT5C", "NME6", 
#                              "NUDT2")
## MKEGG ##
#pyrimidine_pw_genes_mkegg <- c("CTPS2", "NME4", "NME2", "NME6", "NME1", "NME7", 
#                               "RRM1", "RRM2B", "TYMS", "RRM2", "DUT", "DTYMK",
#                               "CMPK2", "CTPS1")

#unique_pyrimidine_set <- unique(c(pyrimidine_pw_genes, pyrimidine_pw_genes_kegg, 
#                                  pyrimidine_pw_genes_mkegg))


# Using gage to access KEGG gene sets
kegg_gsets <- kegg.gsets(species = "hsa", id.type = 'entrez')
pyrimidine_pw_genes_kegg <- kegg_gsets$kg.sets$`hsa00240 Pyrimidine metabolism`
pyrimidine_pw_genes_kegg <- as.data.frame(bitr(pyrimidine_pw_genes_kegg, 
                                               fromType = "ENTREZID", 
                                               toType = "SYMBOL", 
                                               OrgDb=org.Hs.eg.db, drop = T))[,"SYMBOL"] # 58 genes

# Using MSigDB to access KEGG gene sets
kegg_pyrimidine_metabolism <- read.csv(paste0(PATH, "KEGG/KEGG_PYRIMIDINE_METABOLISM.v2023.1.Hs.csv"))
kegg_pyrimidine_metabolism_genes <- unlist(strsplit(kegg_pyrimidine_metabolism[17,2], 
                                                    ",", fixed = T))
kegg_pyrimidine_metabolism_genes <- kegg_pyrimidine_metabolism_genes[
  kegg_pyrimidine_metabolism_genes != ""]  #98 genes

# Identify which of these is significant in at least one cancer type
driver <- "KRAS"
unique_pyrimidine_set_sig <- unlist(lapply(unique_pyrimidine_set, function(g) {
  sig_v <- c()
  for (i in 1:length(perCancer)) {
    x <- perCancer[[i]]
    if(driver %fin% x$R_i.name) {
      if(g %fin% x$T_k.name) {
        if(x[(x$T_k.name == g) & (x$R_i.name == driver), 'q.value'] < 0.2) {
          sig_v <- c(sig_v, T)
        }
        else {sig_v <- c(sig_v, F)}
      }
    }
  }
  if(length(sig_v[sig_v == T]) < 1) {return(NA)}
  else {return(g)}
}))
unique_pyrimidine_set_sig <- unique_pyrimidine_set_sig[!is.na(
  unique_pyrimidine_set_sig)]

# Call function
create_specific_pathway_hm(perCancer, unique_pyrimidine_set_sig, driver, "tstat")


# Import insulin signaling genes from KEGG/MKEGG
# Downloaded from: https://www.gsea-msigdb.org/gsea/msigdb/cards/KEGG_INSULIN_SIGNALING_PATHWAY
insulin_signaling_genes <- read.csv(paste0(PATH, "KEGG/KEGG_INSULIN_SIGNALING_PATHWAY.v2023.2.Hs.csv"), 
                                    header = T, check.names = F)[17,2]
insulin_signaling_genes <- unlist(strsplit(insulin_signaling_genes, ",", fixed = T)) 
insulin_signaling_genes <- insulin_signaling_genes[insulin_signaling_genes != "PIK3CA"] # 136 genes
insulin_signaling_genes <- c(insulin_signaling_genes, "KBTBD2")

# Import insulin secretion genes (hsa04911), which is significantly downregulated
# Downloaded from: https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+path:hsa04911
insulin_secretion_table <- read.table(paste0(PATH, "KEGG/KEGG_INSULIN_SECRETION_hsa04911.txt"),
                                    sep = "\t", header = F)
insulin_secretion_ids_kegg <- unlist(lapply(insulin_secretion_table[,1], function(entry) {
  entry_spl <- unlist(strsplit(entry, " ", fixed = T))
  entry_spl <- entry_spl[!(entry_spl == "")]
  entry <- entry_spl[grepl("hsa:", entry_spl)]
  entry_id <- unlist(strsplit(entry, ":", fixed = T))[2]
  return(entry_id)
}))
mapping_kegg <- as.data.frame(bitr_kegg(insulin_secretion_ids_kegg, 
                                        fromType = "kegg", toType = "uniprot", 
                                        drop = T, organism = "hsa"))
mapping <- as.data.frame(bitr(mapping_kegg$uniprot, fromType = "UNIPROT", 
                              toType = "SYMBOL", drop = T, OrgDb = "org.Hs.eg.db"))
insulin_secretion_genes <- unique(mapping$SYMBOL)

driver <- "PIK3CA"
unique_insulin_secretion_sig <- unlist(lapply(insulin_secretion_genes, function(g) {
  sig_v <- c()
  for (i in 1:length(perCancer)) {
    x <- perCancer[[i]]
    if(driver %fin% x$R_i.name) {
      if(g %fin% x$T_k.name) {
        if(x[(x$T_k.name == g) & (x$R_i.name == driver), 'q.value'] < 0.2) {
          sig_v <- c(sig_v, T)
        }
        else {sig_v <- c(sig_v, F)}
      }
    }
  }
  if(length(sig_v[sig_v == T]) < 1) {return(NA)}
  else {return(g)}
}))
unique_insulin_secretion_sig <- unique_insulin_secretion_sig[!is.na(
  unique_insulin_secretion_sig)]
#unique_insulin_secretion_sig <- c(unique_insulin_secretion_sig, "KBTBD2")

# Call function
create_specific_pathway_hm(perCancer[!(names(perCancer) %fin% c("GBM"))], 
                           unique_insulin_secretion_sig, driver, "tstat")


############################################################
### PART C: DEPMAP PER-GENE DIFFERENTIAL DEPENDENCY BY DRIVER MUTATION STATUS
############################################################
# DepMap files can be downloaded directly from their website using custom 
# (https://depmap.org/portal/download/custom/) or full (https://depmap.org/portal/download/all/)
# download webpages

# Import CRISPRi, expression, mutation, and CNA data
crispr <- read.csv(paste0(PATH, "DepMap/CRISPR_(DepMap_Public_23Q2+Score,_Chronos).csv"), 
                                    header = T, check.names = F, row.names = 1)
#expression <- read.csv(paste0(PATH, "DepMap/Expression_Public_23Q2.csv"), 
#                                        header = T, check.names = F)
# Alternatively, quantile-normalized expression
expression_qn <- read.csv(paste0(PATH, "DepMap/expression_quantile_normalized_sklearn.csv"), 
                       header = T, check.names = F, row.names = 1)

mutations <- read.csv(paste0(PATH, "DepMap/OmicsSomaticMutations.csv"),
                      header = T, check.names = F)
cnas <- read.csv(paste0(PATH, "DepMap/Copy_Number_(Absolute).csv"), 
                 header = T, check.names = F, row.names = 1)

# Get mutations specifically for driver genes of interest in the analysis
types <- c("MISSENSE", "NONSENSE", "SPLICE_SITE", "NONSTOP")
mutations_idh1 <- mutations[(mutations$HugoSymbol == "IDH1") & 
                              (mutations$VariantInfo %fin% types),]
mutations_kras <- mutations[(mutations$HugoSymbol == "KRAS") & 
                              (mutations$VariantInfo %fin% types),]
mutations_tp53 <- mutations[(mutations$HugoSymbol == "TP53") & 
                              (mutations$VariantInfo %fin% types),]
mutations_pik3ca <- mutations[(mutations$HugoSymbol == "PIK3CA") & 
                                (mutations$VariantInfo %fin% types),]
mutations_ctnnb1 <- mutations[(mutations$HugoSymbol == "CTNNB1") & 
                                (mutations$VariantInfo %fin% types),]
mutations_fbxw7 <- mutations[(mutations$HugoSymbol == "FBXW7") & 
                               (mutations$VariantInfo %fin% types),]

# Optionally limit to only hotspot mutations
hotspot_mut_pik3ca <- c(179234297, 179218303, 179218294, 179203765)
mutations_pik3ca_hotspot <- mutations_pik3ca[mutations_pik3ca$Pos %fin% hotspot_mut_pik3ca,]

# Get CNAs specifically for driver genes of interest in the analysis
drivers <- c("IDH1", "KRAS", "TP53", "PIK3CA", "CTNNB1", "FBXW7")
cnas_drivers <- cnas[, colnames(cnas) %fin% drivers]

# Cell line metadata information
cell_line_sample_info <- read.csv(paste0(PATH, "DepMap/cell_line_sample_info.csv"),
                                 header = T, check.names = F)
cell_line_cancer_type_mapping <- list(
  "ACC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Adrenal Cancer", 'DepMap_ID'],
  "BLCA" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Bladder Cancer", 'DepMap_ID'],
  "BRCA" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Breast Cancer", 'DepMap_ID'],
  "CESC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Cervical Cancer", 'DepMap_ID'],
  "COAD" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Colon/Colorectal Cancer", 'DepMap_ID'],
  "ESCA" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Esophageal Cancer", 'DepMap_ID'],
  "HNSC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Head and Neck Cancer", 'DepMap_ID'],
  "KICH" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Kidney Cancer", 'DepMap_ID'],
  "KIRC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Kidney Cancer", 'DepMap_ID'],
  "KIRP" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Kidney Cancer", 'DepMap_ID'],
  "LGG" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Brain Cancer", 'DepMap_ID'],
  "LIHC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Liver Cancer", 'DepMap_ID'],
  "LUAD" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Lung Cancer", 'DepMap_ID'],
  "LUSC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Lung Cancer", 'DepMap_ID'],
  "MESO" = NA,
  "PAAD" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Pancreatic Cancer", 'DepMap_ID'],
  "PRAD" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Prostate Cancer", 'DepMap_ID'],
  "SARC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Sarcoma", 'DepMap_ID'],
  "STAD" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Gastric Cancer", 'DepMap_ID'],
  "THCA" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Thyroid Cancer", 'DepMap_ID'],
  "UCEC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Endometrial/Uterine Cancer", 'DepMap_ID']
)


#' Adjust the knockout and expression data frames to merge them with the mutation
#' data and "melt" them so that we can more easily make plots
#' @param depmap_df a DepMap knockout, expression, or CNA data frame
#' @param mutation_df a DepMap mutation data frame
#' @param genes_of_interest a vector of the Hugo IDs of genes whose mutation
#' status is of interest
adjust_depmap_df <- function(depmap_df, mutation_df, genes_of_interest, 
                             cna_df, del_or_amp_vect) {
  
  # Melt data frame to get it ready for plots
  depmap_df_melt <- melt(as.matrix(depmap_df))
  print(head(depmap_df_melt))
  
  # Adjust column names
  depmap_id_col <- ifelse(grepl("ACH", depmap_df_melt[1,1]), 1, 2)
  colnames(depmap_df_melt)[depmap_id_col] <- "depmap_id"
  gene_id_col <- ifelse(depmap_id_col == 1, 2, 1)
  colnames(depmap_df_melt)[gene_id_col] <- "gene"
  print(head(depmap_df_melt))
  
  # Make the "value" column numeric, the "gene" column a character
  depmap_df_melt$value <- as.numeric(unlist(depmap_df_melt$value))
  depmap_df_melt$gene <- as.character(unlist(depmap_df_melt$gene))
  
  # Add mutation status of genes of interest 
  unique_cls <- unique(intersect(mutation_df$ModelID, depmap_df_melt$depmap_id))
  print(length(unique_cls))
  
  # Get the cell lines that have a mutation
  unique_drivers <- unique(mutation_df$HugoSymbol)
  mutation_df_driver_cols <- as.data.frame(lapply(unique_drivers, function(driver) {
    cls <- mutation_df[mutation_df$HugoSymbol == driver, 'ModelID']
    return(as.factor(unlist(lapply(unique_cls, function(cl) 
      ifelse(cl %fin% cls, 1, 0)))))
  }))
  colnames(mutation_df_driver_cols) <- unique_drivers
  mutation_df_new <- cbind(data.frame("depmap_id" = unique_cls), 
                           mutation_df_driver_cols)
  print(head(mutation_df_new))

  depmap_df_melt <- merge(depmap_df_melt, mutation_df_new, by = "depmap_id")
  #depmap_df_melt[is.na(depmap_df_melt)] <- 0
  
  # Get the cell lines that have an amplification or deletion
  if(length(cna_df) > 1) {
    cna_df_goi <- cna_df[, genes_of_interest]
    colnames(cna_df_goi) <- unlist(lapply(colnames(cna_df_goi), function(x)
      paste0(x, ".CNA")))
    genes_of_interest_cna <- colnames(cna_df_goi)
    cna_df_goi$depmap_id <- rownames(cna_df_goi)
    
    depmap_df_melt <- merge(depmap_df_melt, cna_df_goi, by = "depmap_id")
    cna_cols <- which(colnames(depmap_df_melt) %fin% genes_of_interest_cna)
    depmap_df_melt[, cna_cols] <- lapply(1:length(cna_cols), function(i) {
      col <- depmap_df_melt[, cna_cols[i]]
      bucketed_col <- bucket_cna(col, del_or_amp_vect[i])
      bucketed_col <- as.factor(bucketed_col)
      return(bucketed_col)
    })
  }
  
  return(depmap_df_melt)
}

#' Bucket CNA values, given whether we want to prioritize deletions or 
#' amplifications
#' @param cna_vals vector of CNA values 
#' @param deletionOrAmp either "deletion" or "amplification" to indicate which 
#' we are prioritizing
#' Define thresholds based on this thread: https://forum.depmap.org/t/defining-deep-deletions-and-amplifications/710/3,
#' based on definitions from the TCGA: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/#copy-number-estimation
bucket_cna <- function(cna_vals, deletionOrAmp) {
  bucketed_cna_vals <- unlist(lapply(cna_vals, function(x) {
    # Reported as log2(CN + 1). We want just the CN.
    x_mod <- 2^(x) - 1
    if(deletionOrAmp == "deletion") {
      # we use 3 here because there is a pseudocount of 1; if no pseudocount, use log2(2)
      #if (x < log2(3)) {return(1)}
      if(x < 2^(-1.2)) {return(1)}
      else {return(0)}
    } else {
      #if (x > log2(3)) {return(1)}
      if(x > (2^0.75)) {return(1)}
      
      else {return(0)}
    }
  }))
  return(bucketed_cna_vals)
}

mutations_drivers <- do.call(rbind, list(#mutations_fbxw7, mutations_ctnnb1, 
                                         mutations_idh1, mutations_kras, 
                                         mutations_pik3ca, mutations_tp53))
onc_or_ts <- list("TP53" = "deletion", "PIK3CA" = "amplification", 
                  "KRAS" = "amplification", "IDH1" = "amplification",
                  "CTNNB1" = "amplification", "FBXW7" = "deletion")
crispr <- adjust_depmap_df(crispr, mutations_drivers,  
                           c("TP53", "PIK3CA", "KRAS", "IDH1"),
                           cnas_drivers, onc_or_ts)
expression_qn <- adjust_depmap_df(expression_qn, mutations_drivers, 
                               c("TP53", "PIK3CA", "KRAS", "IDH1"),
                               cnas_drivers, onc_or_ts)
cnas <- adjust_depmap_df(cnas, mutations_drivers, 
                         c("TP53", "PIK3CA", "KRAS", "IDH1"),
                         NA, NA)
# Convert CNA values to log2(CN + 1)
cnas$value <- unlist(lapply(cnas$value, function(cn) log2(cn + 1)))


#' Plot boxplot of gene(s) dependency across cell lines, by driver mutation status
#' @param depmap_df a DepMap knockout or expression data frame
#' @param gene_of_interest a Hugo ID for a given driver gene of interest
#' @param targ the name of the target gene(s) of interest
#' @param yaxis_lab y-axis label; e.g. cell viability
make_dependency_boxplot <- function(depmap_df, gene_of_interest, targ, yaxis_lab) {
  
  # Subset to just the target gene(s)
  depmap_df <- depmap_df[depmap_df$gene %fin% targ,]
  
  # Uncomment for simpler version
  #bxp_pub <- ggplot(depmap_df, aes_string(x = "gene", y = "value", 
                                            #fill = gene_of_interest)) + 
  #  geom_boxplot() + ylab(yaxis_lab) +
  #  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #  stat_compare_means(method = "wilcox.test")
    
  col_of_interest <- which(colnames(depmap_df) == gene_of_interest)
  depmap_df[,col_of_interest] <- unlist(lapply(
    as.numeric(unlist(depmap_df[,col_of_interest])), function(x) {
      return(ifelse(x == 1, paste0(gene_of_interest,".WT"), 
                    paste0(gene_of_interest, ".Mut")))
  }))
  #print(head(depmap_df))
  
  bxp_pub <- ggplot(depmap_df, aes_string(x = gene_of_interest, y = "value", 
                                                 fill = gene_of_interest)) + 
      geom_boxplot() + ylab(yaxis_lab) + theme_minimal() +
      theme(axis.text.x = element_text(
        face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold", size = 12),
        axis.ticks = element_blank(), axis.line = element_blank(),
        legend.position = "none", 
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5)) + 
      ggtitle(paste0("\u0394", targ)) + coord_flip() +
      scale_fill_manual(values=c("#BC3C29FF","#0072B5FF")) +  
      #scale_fill_manual(values=c("#6F99ADFF","#E18727FF")) + 
      scale_x_discrete(labels = c(paste0(gene_of_interest, "\nMutated"), 
                                  paste0(gene_of_interest, "\nUnmutated"))) +
      geom_signif(comparisons = list(c(paste0(gene_of_interest, ".WT"),
                                       paste0(gene_of_interest, ".Mut"))), 
                  na.rm = T, test = "wilcox.test", map_signif_level = T) 
    #stat_compare_means(method = "wilcox.test") #comparisons = list(c(0,1)) #, map_signif_level =c("***"=0.001, "**"=0.01, "*"=0.05)
    

  print(bxp_pub)
}

# Call function for various drivers, such as KRAS with one or more target
# genes, such as NT5E 
make_dependency_boxplot(crispr, "KRAS", c("NT5E"), "Cell Viability")
make_dependency_boxplot(crispr, "PIK3CA", c("KBTBD2"), "Cell Viability")

# Repeat just in a particular cell line, like colon or breast
make_dependency_boxplot(crispr[crispr$depmap_id %fin% cell_line_cancer_type_mapping[["COAD"]],], 
                        "KRAS", c("NT5E"), "Cell Viability")
make_dependency_boxplot(crispr[crispr$depmap_id %fin% cell_line_cancer_type_mapping[["BRCA"]],], 
                        "PIK3CA", c("KBTBD2"), "Cell Viability")

# Just driver-relevant cancer types
pik3ca_relevant_cancer_types <- c("BLCA", "BRCA", "COAD", "ESCA", "GBM", "HNSC", 
                                  "LGG", "LUSC", "UCEC")
pik3ca_relevant_depmap_ids <- unlist(lapply(pik3ca_relevant_cancer_types, function(ct)
  cell_line_cancer_type_mapping[[ct]]))
make_dependency_boxplot(crispr[crispr$depmap_id %fin% pik3ca_relevant_depmap_ids,], 
                        "PIK3CA", c("KBTBD2"), "Cell Viability")


############################################################
### PART D: DEPMAP PER-GENE EXPRESSION OF DRIVER V. VIABILITY 
### OF CELLS UPON TARGET KNOCKOUT
############################################################
#' Run synthetic lethality analysis for a driven gene across all of its cancer 
#' types according to a given q-value threshold and multiple ways of trimming
#' the number of tests (e.g. limiting to only targets with drugs from DrugBank,
#' and limiting to those that have significantly different dependency by driver
#' mutation status)
#' @param driver_name the gene name of a given driver gene of interest
#' @param perCancer the named list of Dyscovr output files, limited to
#' cancer types containing the provided driver gene
#' @param qval_thres a q-value threshold for per-cancer significance
#' @param use_pc_hits a T/F value indicating whether we are limiting to targets
#' that are statistically significant pan-cancer
#' @param pc_hits if restricting to PC hits, a vector of PC hits
#' @param use_drugbank a T/F value indicating whether we are limiting to targets
#' with a drug in DrugBank'
#' @param drugbank_gns if restricting by DrugBank, a vector of gene names with 
#' known associated drugs from the database
#' @param use_dependency_check a T/F value indicating whether we are checking for
#' differential dependency by mutation status prior to looking for a correlation
#' between unmutated expression of driver and target dependency
#' @param cell_line_cancer_type_mapping a mapping between depmap IDs and cancer type
#' @param method either "Spearman", "Pearson", "Regression", or "Wilcoxon" to 
#' denote the method we are using to determine significant synthetic lethality
#' @param ts_or_oncogene either "ts" or "oncogene" to denote if the driver is a 
#' tumor suppressor or oncogene, for Wilcoxon analysis
#' @param thres_cls an integer denoting the minimum number of cell lines in a 
#' given cancer type that do not have a mutation in the given driver
#' @param pval_thres a p-value threshold for significance
get_synthetic_lethals <- function(driver_name, perCancer, qval_thres, use_pc_hits, 
                                  pc_hits, use_drugbank, drugbank_gns, 
                                  use_dependency_check, crispr, expression,
                                  cell_line_cancer_type_mapping, method, 
                                  ts_or_oncogene, thres_cls = 15, pval_thres = 1) {
  
  per_cancer_results <- lapply(1:length(perCancer), function(i) {
    ct <- names(perCancer)[i]
    master_df <- perCancer[[i]]
    
    tophits_ct <- as.character(unlist(master_df[(master_df$q.value < qval_thres) & 
                                                  #(master_df$estimate > 0) & 
                                                  (master_df$R_i.name == driver_name), 
                                                'T_k.name']))
    #print(head(tophits_ct))
    
    # Limit CRISPR and expression DepMap data to just this cancer type
    crispr_sub <- crispr[crispr$depmap_id %fin% 
                           c(cell_line_cancer_type_mapping[[ct]]),]
    expression_sub <- expression[expression$depmap_id %fin% 
                                   c(cell_line_cancer_type_mapping[[ct]]),]
    
    # Ensure there are at least a certain number (e.g. 10-15) lines that do not have
    # a mutation in the given driver
    if(nrow(expression_sub[which(expression_sub[,driver_name] == 0),]) < thres_cls) {
      print(paste("There are fewer than", 
                  paste(thres_cls, 
                        paste("cell lines that lack a driver mutation in",
                              paste0(driver_name, 
                                     paste(" in", paste0(ct, ". Returning NA.")))))))
      return(NA)
    }
    
    
    # If desired, limit to DrugBank genes
    if(use_drugbank) {
      tophits_ct <- tophits_ct[tophits_ct %fin% drugbank_gns]
    }
    
    # If desired, limit to pan-cancer hits
    if(use_pc_hits) {
      tophits_ct <- tophits_ct[tophits_ct %fin% pc_hits]
    }
    print(head(tophits_ct))
    
    # If desired, limit be differential dependency; particularly, significantly  
    # lower viability in the unmutated case
    if(use_dependency_check) {
      wilcoxon_res <- data.table("gene" = tophits_ct)
      wilcoxon_res$pval <- unlist(lapply(wilcoxon_res$gene, function(g) {
        dep_mut <- crispr_sub[(crispr_sub$gene == g) & 
                                (crispr_sub[,driver_name] == 1), 'value']
        dep_noMut <- crispr_sub[(crispr_sub$gene == g) &
                                  (crispr_sub[,driver_name] == 0), 'value']
        if((length(dep_noMut) < 1) | (length(dep_mut) < 1)) {return(NA)}
        if(is.na(dep_noMut) | is.na(dep_mut)) {return(NA)}
        #res <- wilcox.test(dep_mut, dep_noMut, alternative = "greater", exact = F)
        res <- wilcox.test(dep_mut, dep_noMut, exact = F)
        return(res$p.value)
      }))
      wilcoxon_res$qval <- qvalue(wilcoxon_res$pval)$qvalues
      wilcoxon_res <- wilcoxon_res[order(wilcoxon_res$pval),]
      signif_hits <- wilcoxon_res[wilcoxon_res$qval < qval_thres, 'gene']  
      tophits_ct <- tophits_ct[tophits_ct %fin% signif_hits]
    }
    
    tophits_depmap <- NA
    
    # Calculate correlation results (Spearman or Pearson)
    if(length(tophits_ct) > 0) {
      if((method == "Spearman") | (method == "Pearson")) {
        corr_res_sig <- dependency_by_expression_analysis(crispr_sub,
                                                          expression_sub,
                                                          tophits_ct, method,
                                                          driver_name, 0)
        corr_res_sig <- na.omit(corr_res_sig)
        if(length(corr_res_sig) == 0) {return(NA)}
        if(nrow(corr_res_sig) == 0) {return(NA)}
        #corr_res_sig$qval <- qvalue(corr_res_sig$pval)$qvalues
        corr_res_sig$padj <- p.adjust(corr_res_sig$pval, "BH")
        
        #corr_res_sig <- corr_res_sig[order(corr_res_sig$qval),]
        corr_res_sig <- corr_res_sig[order(corr_res_sig$padj),]
        
        tophits_depmap <- corr_res_sig[(corr_res_sig$pval < pval_thres),] #& 
                                         #(corr_res_sig$stat > 0),]
        
      } else if (method == "Wilcoxon") {
        wilcox_res_sig <- perform_wilcoxon_analysis(crispr_sub, expression_sub,
                                            tophits_ct, driver_name, NA, 
                                            ts_or_oncogene)
        print(wilcox_res_sig)
        wilcox_res_sig <- na.omit(wilcox_res_sig)
        if(length(wilcox_res_sig) == 0) {return(NA)}
        if(nrow(wilcox_res_sig) == 0) {return(NA)}
        wilcox_res_sig$padj <- p.adjust(wilcox_res_sig$pval, "BH")
        wilcox_res_sig <- wilcox_res_sig[order(wilcox_res_sig$padj),]
        
        #tophits_depmap <- wilcox_res_sig
        tophits_depmap <- wilcox_res_sig[(wilcox_res_sig$pval < pval_thres),]
        
      } else if (method == "Regression") {
        lm_res_sig <- perform_lm_analysis(crispr_sub, expression_sub,
                                          tophits_ct, driver_name, 0)
        
        lm_res_sig <- na.omit(lm_res_sig)
        if(length(lm_res_sig) == 0) {return(NA)}
        if(nrow(lm_res_sig) == 0) {return(NA)}
        unique_terms <- unique(lm_res_sig$term)
        
        adj_dfs <- lapply(unique_terms, function(term) {
          lm_res_sig_sub <- lm_res_sig[lm_res_sig$term == term,]
          lm_res_sig_sub$padj <- p.adjust(lm_res_sig_sub$pval, "BH")
          return(lm_res_sig_sub)
        })
        lm_res_sig <- do.call(rbind, adj_dfs)
        #lm_res_sig$padj <- p.adjust(lm_res_sig$pval, "BH")
        lm_res_sig <- lm_res_sig[order(lm_res_sig$padj),]
        
        tophits_depmap <- lm_res_sig[(lm_res_sig$pval < pval_thres),] 
                                      #& (lm_res_sig$stat > 0),]
      } else {
        print("Only implemented methods are Spearman, Pearson, Wilcoxon, and 
              Regression. Please try again with one of these methods specified.")
        return(NA)
      }
      
      if(length(tophits_depmap) == 0) {return(NA)}
      if(nrow(tophits_depmap) == 0) {return(NA)}
      #dependency_by_expression_analysis_plot(crispr_sub, expression_sub, 
                                             #tophits_depmap$gene, driver_name, 0)
      
      tophits_depmap$cancer.type <- rep(ct, times = nrow(tophits_depmap))
      return(tophits_depmap)
    } else {return(NA)}
  })
  per_cancer_results <- per_cancer_results[!is.na(per_cancer_results)]
  full_results <- do.call(rbind, per_cancer_results)
  
  return(full_results)
}

#' Explore the relationship between the expression of a driver gene (mutated or
#' unmutated) and cell viability when each of a given set of targets is knocked
#' out via CRISPRi using Spearman or Pearson correlation
#' @param crispr a DepMap data frame with CRISPRi knockout data
#' @param expression a DepMap data frame with gene expression data
#' @param list_of_genes a vector of target genes for plotting
#' @param method either "Spearman" or "Pearson"
#' @param driver_mut the name of the driver gene of interest
#' @param mut_status the mutation status of the driver gene, 0 or 1
dependency_by_expression_analysis <- function(crispr, expression, list_of_genes, 
                                              method, driver_mut, mut_status) {
  
  expression <- expression[expression$gene %fin% c(driver_mut),]
  
  if(!is.na(driver_mut)) {
    expression <- expression[expression[,which(
      colnames(expression) == driver_mut)] == mut_status,]
  }
  expression_driver <- expression[expression$gene == driver_mut, 
                                  c("depmap_id", "gene", "value")]
  
  corr_dfs <- lapply(list_of_genes, function(gene) {
    print(gene)
    if(!(gene %fin% crispr$gene)) {return(NA)}
    crispr <- crispr[crispr$gene %fin% c(gene),]
    crispr_gene <- crispr[crispr$gene == gene, c("depmap_id", "gene", "value")]
    
    data_input <- merge(crispr_gene, expression_driver, by = "depmap_id")
    
    if(nrow(data_input) < 5) {return(NA)}
    
    if(method == "Spearman") {
      spearman <- cor.test(data_input$value.x, data_input$value.y, 
                           method = "spearman")
      spearman_stat <- as.numeric(spearman$estimate)
      spearman_pval <- spearman$p.value
      
      print(spearman_pval)
      return(data.frame("stat" = spearman_stat, "pval" = spearman_pval, 
                        "gene" = gene))
    } else {
      pearson <- cor.test(data_input$value.x, data_input$value.y, 
                           method = "pearson")
      pearson_stat <- as.numeric(pearson$estimate)
      pearson_pval <- pearson$p.value
      
      print(pearson_pval)
      return(data.frame("stat" = pearson_stat, "pval" = pearson_pval, 
                        "gene" = gene))
    }
  })
  
  corr_dfs <- corr_dfs[!is.na(corr_dfs)]
  corr_df <- do.call(rbind, corr_dfs)
  
  return(corr_df)
}

#' Explore the relationship between the expression of a driver gene (mutated or
#' unmutated) and cell viability when each of a given set of targets is knocked
#' out via CRISPRi using Wilcoxon rank sum test
#' @param crispr a DepMap data frame with CRISPRi knockout data
#' @param expression a DepMap data frame with gene expression data
#' @param list_of_genes a vector of target genes for plotting
#' @param driver_mut the name of the driver gene of interest
#' @param mut_status the mutation status of the driver gene, 0 or 1
#' @param ts_or_oncogene either "ts" or "oncogene" to indicate if driver is a
#' tumor suppressor or oncogene
perform_wilcoxon_analysis <- function(crispr, expression, list_of_genes, driver_mut, 
                              mut_status, ts_or_oncogene) {
  
  expression <- expression[expression$gene %fin% c(driver_mut),]
  
  if(!is.na(mut_status)) {
    expression <- expression[expression[,which(
      colnames(expression) == driver_mut)] == mut_status,]
  }
  expression_driver <- expression[expression$gene == driver_mut, 
                                  c("depmap_id", "gene", "value")]
  
  # Split into "over-" or "under-expression" groups based on criteria given in 
  # https://www.cell.com/cell/pdf/S0092-8674(14)00977-5.pdf; "under-expressed" if
  # below the 10th percentile of its expression levels across all unmutated 
  # samples in the dataset (for tumor suppressors), and similarly, as 
  # "over-expressed" if its expression level is above its 90th percentile (for
  # oncogenes)
  if(ts_or_oncogene == "oncogene") {
    percentile_90 <- as.numeric(quantile(expression_driver$value, c(.90)))
    expression_driver_overexpr <- expression_driver[expression_driver$value > 
                                                           percentile_90,]
    expression_driver_reg <- expression_driver[(expression_driver$value <= 
                                                  percentile_90) ,]
    
    wilcox_dfs <- lapply(list_of_genes, function(gene) {
      print(gene)
      if(!(gene %fin% crispr$gene)) {return(NA)}
      crispr <- crispr[crispr$gene %fin% c(gene),]
      crispr_gene <- crispr[crispr$gene == gene, c("depmap_id", "gene", "value")]
      
      data_input_reg <- merge(crispr_gene, expression_driver_reg, 
                              by = "depmap_id", all.x = F)
      data_input_overexpr <- merge(crispr_gene, expression_driver_overexpr, 
                                   by = "depmap_id", all.x = F)
      print(head(data_input_reg))
      print(head(data_input_overexpr))

      #if((nrow(data_input_overexpr) < 5) | (nrow(data_input_reg) < 5)) {
      #  return(NA)
      #}
      
      wilcox.res <- wilcox.test(data_input_overexpr$value.x, data_input_reg$value.x,
                                alternative = "greater")
      wilcox_stat <- as.numeric(wilcox.res$statistic)
      wilcox_pval <- wilcox.res$p.value
      
      print(wilcox_pval)
      return(data.frame("stat" = wilcox_stat, "pval" = wilcox_pval, 
                        "gene" = gene))
    })
    
    wilcox_dfs <- wilcox_dfs[!is.na(wilcox_dfs)]
    wilcox_df <- do.call(rbind, wilcox_dfs)
    
    return(wilcox_df)
    
  } else if (ts_or_oncogene == "ts") {
    percentile_10 <- as.numeric(quantile(expression_driver$value, c(.10)))
    expression_driver_underexpr <- expression_driver[expression_driver$value < 
                                                           percentile_10,]
    expression_driver_reg <- expression_driver[expression_driver$value >= 
                                                 percentile_10,]
    
    wilcox_dfs <- lapply(list_of_genes, function(gene) {
      print(gene)
      if(!(gene %fin% crispr$gene)) {return(NA)}
      crispr <- crispr[crispr$gene %fin% c(gene),]
      crispr_gene <- crispr[crispr$gene == gene, c("depmap_id", "gene", "value")]
      
      data_input_reg <- merge(crispr_gene, expression_driver_reg, 
                              by = "depmap_id", all.x = F)
      data_input_underexpr <- merge(crispr_gene, expression_driver_underexpr, 
                                   by = "depmap_id", all.x = F)
      print(head(data_input_reg))
      print(head(data_input_underexpr))
      
      #if((nrow(data_input_underexpr) < 5) | (nrow(data_input_reg) < 5)) {
      #  return(NA)
      #}
      
      wilcox.res <- wilcox.test(data_input_underexpr$value.x, data_input_reg$value.x,
                                alternative = "greater")
      wilcox_stat <- as.numeric(wilcox.res$statistic)
      wilcox_pval <- wilcox.res$p.value
      
      print(wilcox_pval)
      return(data.frame("stat" = wilcox_stat, "pval" = wilcox_pval, 
                        "gene" = gene))
    })
    
    wilcox_dfs <- wilcox_dfs[!is.na(wilcox_dfs)]
    wilcox_df <- do.call(rbind, wilcox_dfs)
    
    return(wilcox_df)
    
  } else {
    print("Only implemented for 'oncogene' or 'ts'. Please pick one.")
    return(NA)
  }
}

#' Explore the relationship between the expression of a driver gene (mutated or
#' unmutated) and cell viability when each of a given set of targets is knocked
#' out via CRISPRi using a simple linear regression model
#' @param crispr a DepMap data frame with CRISPRi knockout data
#' @param expression a DepMap data frame with gene expression data
#' @param list_of_genes a vector of target genes for plotting
#' @param driver_mut the name of the driver gene of interest
#' @param mut_status the mutation status of the driver gene, 0 or 1
#' @param restrict_or_incorp either "restrict" to limit to samples with the 
#' particular mutation status in the given driver gene, or "incorp" to include
#' a term for mutation status in the model itself (the default)
perform_lm_analysis <- function(crispr, expression, list_of_genes, driver_mut, 
                                mut_status, restrict_or_incorp = "incorp") {
  
  expression <- expression[expression$gene %fin% c(driver_mut),]
  
  if(!(is.na(driver_mut)) & (restrict_or_incorp == "restrict")) {
    expression <- expression[expression[,which(
      colnames(expression) == driver_mut)] == mut_status,]
  } 
  expression_driver <- expression[expression$gene == driver_mut, ]
  
  driver_col <- as.integer(which(colnames(expression_driver) == driver_mut))
  expression_driver <- cbind(expression_driver[, c(driver_col)],
                             expression_driver[, c("depmap_id", "gene", "value")])
  colnames(expression_driver)[1] <- "driver_mut_status"
  colnames(expression_driver)[
    which(colnames(expression_driver) == "value")] <- "expression_val"

  lm_dfs <- lapply(list_of_genes, function(gene) {
    print(gene)
    if(!(gene %fin% crispr$gene)) {return(NA)}
    crispr <- crispr[crispr$gene %fin% c(gene),]
    crispr_gene <- crispr[crispr$gene == gene, c("depmap_id", "gene", "value")]
    colnames(crispr_gene)[
      which(colnames(crispr_gene) == "value")] <- "crispr_val"
    
    data_input <- merge(crispr_gene, expression_driver, by = "depmap_id")
    data_input$driver_mut_status <- as.factor(data_input$driver_mut_status)
    print(head(data_input))
    
    if(nrow(data_input) < 5) {return(NA)}
    
    lm_res <- NA
    tryCatch({
      if((restrict_or_incorp == "incorp") & 
         (length(unique(data_input$driver_mut_status)) > 1)) {
        lm_res <- tidy(lm(crispr_val ~ expression_val + driver_mut_status, 
                          data = data_input))
      } else {
        lm_res <- tidy(lm(crispr_val ~ expression_val, data = data_input))
      }
    }, error = function(cond) {
      print(cond)
    })
    
    print(lm_res)
    if(length(lm_res) > 1) {
      lm_res <- lm_res[lm_res$term != "(Intercept)",]
      lm_res_stats <- as.numeric(lm_res$estimate)
      lm_res_pvals <- as.numeric(lm_res$p.value)
      lm_res_terms <- lm_res$term
      return(data.frame("stat" = lm_res_stats, "pval" = lm_res_pvals, 
                        "term" = lm_res_terms,
                        "gene" = rep(gene, times = length(lm_res_stats))))
    }
  })
  
  lm_dfs <- lm_dfs[!is.na(lm_dfs)]
  lm_df <- do.call(rbind, lm_dfs)
  print(head(lm_df))
  print(dim(lm_df))
  
  return(lm_df)
}


#' Explore the relationship between the expression of a driver gene (mutated or
#' unmutated) and cell viability when each of a given set of targets is knocked
#' out via CRISPRi using a simple linear regression model across all cancer types
#' with a one-hot encoded term for cancer type
#' @param crispr a DepMap data frame with CRISPRi knockout data
#' @param expression a DepMap data frame with gene expression data
#' @param cna a DepMap data frame with CNA data (optional)
#' @param driver_mut the name of the driver gene of interest
#' @param perCancer the named list of Dyscovr output files, limited to
#' cancer types containing the provided driver gene
#' @param qval_thres a q-value threshold
#' @param use_pc_hits a T/F value indicating whether we are limiting to targets
#' that are statistically significant pan-cancer
#' @param pc_hits if restricting to PC hits, a vector of PC hits
#' @param use_drugbank a T/F value indicating whether we are limiting to targets
#' with a drug in DrugBank'
#' @param drugbank_gns if restricting by DrugBank, a vector of gene names with 
#' known associated drugs from the database
#' @param use_dependency_check a T/F value indicating whether we are checking for
#' differential dependency by mutation status prior to looking for a correlation
#' between unmutated expression of driver and target dependency
#' @param cell_line_cancer_type_mapping a mapping between depmap IDs and cancer type
#' @param terms_of_interest a vector of terms we'd like to include as variables
#' (options include 'mutation', 'cna', and 'expression')
#' @param thres_cls an integer denoting the minimum number of cell lines in a 
#' given cancer type that do not have a mutation in the given driver
perform_lm_analysis_crosscancers <- function(crispr, expression, cna, driver_mut, 
                                             perCancer, qval_thres, use_pc_hits, 
                                             pc_hits, use_drugbank, drugbank_gns, 
                                             use_dependency_check, 
                                             cell_line_cancer_type_mapping,
                                             terms_of_interest, thres_cls = 15) {
  
  tophits <- unique(unlist(lapply(perCancer, function(master_df) {
    hits <- as.character(unlist(master_df[(master_df$q.value < qval_thres) & 
                                            (master_df$R_i.name == driver_name),
                                          'T_k.name']))
    return(hits)
  })))

  # Limit CRISPR, CNA, and expression DepMap data to just cancer types in which 
  # driver is mutated
  cts <- names(perCancer)
  depmap_ids <- unlist(lapply(1:length(cell_line_cancer_type_mapping), function(i) {
    if(names(cell_line_cancer_type_mapping)[i] %fin% cts) {
      return(unlist(cell_line_cancer_type_mapping[[i]]))
    } else {return(NA)}
  }))
  depmap_ids <- depmap_ids[!is.na(depmap_ids)]

  crispr_sub <- crispr[crispr$depmap_id %fin% depmap_ids,]
  expression_sub <- expression[expression$depmap_id %fin% depmap_ids,]
  cna_sub <- NA
  if(length(cna) > 1) {
    cna_sub <- cna[cna$depmap_id %fin% depmap_ids, ]
  }
  
  # Add cancer type
  new_cols <- lapply(1:length(cell_line_cancer_type_mapping), function(i) {
    ct <- names(cell_line_cancer_type_mapping)[i]
    print(ct)
    depmaps_ct <- unlist(cell_line_cancer_type_mapping[[i]])
    intersect_depmaps <- intersect(expression_sub$depmap_id, 
                                   intersect(crispr_sub$depmap_id, depmaps_ct))
    if(length(cna_sub) > 1) {
      intersect_depmaps <- intersect(intersect_depmaps, cna_sub$depmap_id)
    }
    return(data.table("depmap_id" = intersect_depmaps, 
                      "cancer_type" = rep(ct, times = length(intersect_depmaps))))
  })
  ct_df <- distinct(do.call(rbind, new_cols))
  print(head(ct_df))

  #crispr_sub <- merge(crispr_sub, ct_df, by = "depmap_id", all = T)
  crispr_sub <- left_join(crispr_sub, ct_df, by = "depmap_id", 
                          relationship = "many-to-many")
  expression_sub <- left_join(expression_sub, ct_df, by = "depmap_id", 
                              relationship = "many-to-many")
  if(length(cna_sub) > 1) {
    cna_sub <- left_join(cna_sub, ct_df, by = "depmap_id", 
                         relationship = "many-to-many")
  }

  # Ensure there are at least a certain number (e.g. 10) lines that do not have
  # a mutation in the given driver AND that do have a mutation in the given driver
  if(nrow(expression_sub[which(expression_sub[,driver_name] == 0),]) < thres_cls) {
    print(paste("There are fewer than", 
                paste(thres_cls, 
                      paste("cell lines that lack a driver mutation in",
                            paste0(driver_name, 
                                   paste(" in", paste0(ct, ". Returning NA.")))))))
    return(NA)
  }
  if(nrow(expression_sub[which(expression_sub[,driver_name] == 1),]) < thres_cls) {
    print(paste("There are fewer than", 
                paste(thres_cls, 
                      paste("cell lines that have a driver mutation in",
                            paste0(driver_name, 
                                   paste(" in", paste0(ct, ". Returning NA.")))))))
    return(NA)
  }
  
  # If desired, limit to DrugBank genes
  if(use_drugbank) {
    tophits <- tophits[tophits %fin% drugbank_gns]
  }
  
  # If desired, limit to pan-cancer hits
  if(use_pc_hits) {
    tophits <- tophits[tophits %fin% pc_hits]
  }
  print(head(tophits))
  
  # If desired, limit be differential dependency; particularly, significantly  
  # lower viability in the unmutated case
  if(use_dependency_check) {
    wilcoxon_res <- data.table("gene" = tophits)
    wilcoxon_res$pval <- unlist(lapply(wilcoxon_res$gene, function(g) {
      dep_mut <- crispr_sub[(crispr_sub$gene == g) & 
                              (crispr_sub[,driver_name] == 1), 'value']
      dep_noMut <- crispr_sub[(crispr_sub$gene == g) &
                                (crispr_sub[,driver_name] == 0), 'value']
      if((length(dep_noMut) < 1) | (length(dep_mut) < 1)) {return(NA)}
      if(is.na(dep_noMut) | is.na(dep_mut)) {return(NA)}
      res <- wilcox.test(dep_mut, dep_noMut, alternative = "greater", exact = F)
      return(res$p.value)
    }))
    wilcoxon_res$qval <- qvalue(wilcoxon_res$pval)$qvalues
    wilcoxon_res <- wilcoxon_res[order(wilcoxon_res$pval),]
    signif_hits <- wilcoxon_res[wilcoxon_res$qval < qval_thres, 'gene']  
    tophits <- tophits[tophits %fin% signif_hits]
  }
  
  tophits_depmap <- NA
  
  # Perform regression
  if(length(tophits) > 0) {
    expression_sub <- expression_sub[expression_sub$gene %fin% c(driver_mut),]
    expression_driver <- expression_sub[expression_sub$gene == driver_mut, ]
    
    driver_col <- as.integer(which(colnames(expression_driver) == driver_mut))
    driver_cna_col <- as.integer(which(colnames(expression_driver) == paste0(driver_mut, ".CNA")))
    expression_driver <- cbind(expression_driver[, c(driver_col, driver_cna_col)],
                               expression_driver[, c("depmap_id", "gene", 
                                                     "value", "cancer_type")])
    colnames(expression_driver)[1] <- "driver_mut_status"
    colnames(expression_driver)[2] <- "driver_cna_status"
    colnames(expression_driver)[
      which(colnames(expression_driver) == "value")] <- "expression_val"
    
    lm_dfs <- lapply(tophits, function(gene) {
      print(gene)
      if(!(gene %fin% crispr_sub$gene)) {return(NA)}
      crispr_sub <- crispr_sub[crispr_sub$gene %fin% c(gene),]
      crispr_gene <- crispr_sub[crispr_sub$gene == gene, c("depmap_id", "gene", 
                                                           "value", "cancer_type")]
      colnames(crispr_gene)[
        which(colnames(crispr_gene) == "value")] <- "crispr_val"
      df_list <- list(crispr_gene, expression_driver)
      
      data_input <- merge(crispr_gene, expression_driver, 
                          by = c("depmap_id", "cancer_type"), all = T)
      data_input$driver_mut_status <- as.factor(data_input$driver_mut_status)
      
      if(length(cna_sub) > 1) {
        if(!(gene %fin% cna_sub$gene)) {return(NA)}
        cna_sub <- cna_sub[cna_sub$gene %fin% c(gene), ]
        cna_gene <- cna_sub[cna_sub$gene == gene, c("depmap_id", "gene",
                                                    "value", "cancer_type")]
        colnames(cna_gene)[
          which(colnames(cna_gene) == "value")] <- "cna_val"
        df_list <- list(crispr_gene, cna_gene, expression_driver)
        
        data_input <- merge(crispr_gene, cna_gene, by = c("depmap_id", "gene", 
                                                          "cancer_type"), all = T)
        data_input <- merge(data_input, expression_driver, 
                            by = c("depmap_id", "cancer_type"), all = T)
        data_input$driver_mut_status <- as.factor(data_input$driver_mut_status)
        data_input$driver_cna_status <- as.factor(data_input$driver_cna_status)
      }
      print(head(data_input))
      
      if(nrow(data_input) < 5) {return(NA)}
      
      lm_res <- NA
      
      # TODO: create a separate function to build formulas as they get long + more complicated
      tryCatch({
        formula <- create_formula(terms_of_interest, data_input)
        #lm_res <- tidy(lm(formula, data = data_input))
        lm_res <- speedglm::speedlm(formula = formula, data = data_input)
        
        summary_table <- tidy(lm_res)
        summary_table <- as.data.table(
          summary_table[, colSums(is.na(summary_table)) < nrow(summary_table)])
        
      }, error = function(cond) {
        print(cond)
      })
      
      #print(summary_table)
      
      if(length(summary_table) > 1) {
        summary_table <- summary_table[summary_table$term != "(Intercept)",]
        lm_res_stats <- as.numeric(summary_table$estimate)
        lm_res_pvals <- summary_table$p.value
        lm_res_terms <- summary_table$term
        print(lm_res_pvals)
        return(data.frame("stat" = lm_res_stats, "pval" = lm_res_pvals, 
                          "term" = lm_res_terms,
                          "gene" = rep(gene, times = length(lm_res_stats))))
      } else {return(NA)}
    })
    
    lm_dfs <- lm_dfs[!is.na(lm_dfs)]
    lm_res <- do.call(rbind, lm_dfs)
    
    lm_res <- na.omit(lm_res)
    if(length(lm_res) == 0) {return(NA)}
    if(nrow(lm_res) == 0) {return(NA)}
    
    unique_terms <- unique(lm_res$term)
    adj_dfs <- lapply(unique_terms, function(term) {
      lm_res_sub <- lm_res[lm_res$term == term,]
      #lm_res_sub$padj <- p.adjust(lm_res_sub$pval, "BH")
      lm_res_sub$qvalue <- qvalue(lm_res_sub$pval)$qvalues
      return(lm_res_sub)
    })
    lm_res_sig <- do.call(rbind, adj_dfs)
    #lm_res_sig <- lm_res_sig[order(lm_res_sig$padj),]
    lm_res_sig <- lm_res_sig[order(lm_res_sig$qvalue),]
    
    tophits_depmap <- lm_res_sig[(lm_res_sig$pval < 1), ] #& 
                                  # (lm_res_sig$stat > 0),]
    if(length(tophits_depmap) == 0) {return(NA)}
    if(nrow(tophits_depmap) == 0) {return(NA)}
    #dependency_by_expression_analysis_plot(crispr_sub, expression_sub, 
    #                                       tophits_depmap[!(grepl("cancer_type", 
    #                                                              tophits_depmap$term)), 
    #                                                        'gene'], 
    #                                       driver_name, 0)
    
    return(tophits_depmap)
  } else {return(NA)}
}


#' Construct linear regression formula using a set of variables indicating which
#' terms we are including in a given model
#' @param terms_of_interest a vector of driver terms of interest to include in the
#' model, possible options are 'mutation', 'cna', and 'expression'
#' @param data_input a data frame with input to the regression, each column
#' containing values for a particular term
create_formula <- function(terms_of_interest, data_input) {
  formula <- "crispr_val ~ "
  
  # Driver mutation status
  if(('mutation' %fin% terms_of_interest) & 
     (length(unique(data_input$driver_mut_status)) > 1)) {
    formula <- paste0(formula, 'driver_mut_status ')
  }
  
  # Driver CNA status
  if(('cna' %fin% terms_of_interest) & 
     (length(unique(data_input$driver_cna_status)) > 1)) {
    formula <- paste(formula, 'driver_cna_status ', sep = "+ ")
  }
  
  # Driver expression
  if('expression' %fin% terms_of_interest) {
    formula <- paste(formula, 'expression_val ', sep = "+ ")
  }
  
  # Target CNA
  if('cna_val' %fin% colnames(data_input)) {
    formula <- paste(formula, 'cna_val ', sep = "+ ")
  }
  
  # Cancer type
  # if(length(unique(data_input$cancer_type)) > 1) {
  if('cancer_type' %fin% colnames(data_input)) {
    formula <- paste(formula, 'cancer_type', sep = "+ ")
  }
  
  print(formula)
  if(formula == "crispr_val ~ ") { return(NA) }
  else {return(formula)}
}

#' Create per-driver visualization of the relationship between driver expression
#' (mutated or unmutated) and cell viability upon putative target gene knockout.
#' @param crispr a DepMap data frame with CRISPRi knockout data
#' @param expression a DepMap data frame with gene expression data
#' @param list_of_genes a vector of target genes for plotting
#' @param driver_mut the name of the driver gene of interest
#' @param mut_status the mutation status of the driver gene, 0 or 1
dependency_by_expression_analysis_plot <- function(crispr,expression, 
                                                   list_of_genes, driver_mut, 
                                                   mut_status) {
  
  expression <- expression[expression$gene %fin% c(driver_mut),]
  
  if(!is.na(driver_mut)) {
    expression <- expression[
      expression[,which(colnames(expression) == driver_mut)] == mut_status,]
  }
  
  expression_driver <- expression[expression$gene == driver_mut, 
                                                   c("depmap_id", "gene", "value")]
  colnames(expression_driver) <- c("depmap_id", "driver", "expression")
  
  crispr_dfs <- lapply(list_of_genes, function(gene) {
    print(gene)
    if(!(gene %fin% crispr$gene)) {return(NA)}
    crispr <- crispr[crispr$gene %fin% c(gene),]
    
    crispr_gene <- crispr[crispr$gene == gene, c("depmap_id", "gene", "value")]
    colnames(crispr_gene) <- c("depmap_id", "gene", "dependency")
    
    data_input <- merge(crispr_gene, expression_driver, by = "depmap_id")
  })
  
  crispr_dfs <- crispr_dfs[!is.na(crispr_dfs)]
  data_input <- do.call(rbind, crispr_dfs)
  #print(head(data_input))
  
  mut_lab <- "Unmutated"
  if(mut_status == 1) {mut_lab <- "Mutated"}
  
  # Create a plot to visualize the correlations
  g <- ggplot(data_input, aes(x = expression, y = dependency)) + geom_point() +
    theme_minimal() + ylab("Cell Viability") + 
    xlab(paste(driver_mut, paste0("(", paste0(mut_lab, paste(")", "Expression"))))) + 
    facet_wrap (~gene, scales = "free_y", labeller = label_bquote(Delta~.(gene))) + 
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(face = "bold", size = 12),
          strip.text = element_text(face = "bold", size = 12)) + 
    geom_smooth(method = "lm", col = "#6F99ADFF") + 
    scale_y_continuous(labels = formatter(nsmall = 1), 
                       breaks = scales::pretty_breaks(n = 5)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 

  print(g)
  
  return(data_input)
}

#' Helper function to round values for graph labels
formatter <- function(...){
  function(x) format(round(x, 1), ...)
}


#' Using the metap package, recombine the p-values from across 
#' cancer types of subtypes to explore global effects
#' @param results a data frame of results from get_synthetic_lethals
#' without subsetting to significant hits- change p-value threshold to 1), with 
#' a 'cancer_type' column denoting the cancer type
#' @param by_term if we are recombining on a regression output with multiple 
#' terms, recombine for each term separately
recombine_pvals <- function(results, by_term = T) {
  # Get all the unique target genes 
  list_of_results <- NA
  if(by_term) {
    unique_terms <- unique(results$term)
    list_of_results <- lapply(unique_terms, function(t) results[results$term == t,])
    names(list_of_results) <- unique_terms
  } else {list_of_results <- list(results)}
  
  unique_targets <- unique(unlist(lapply(list_of_results, function(results) 
    results$gene)))
  
  # For each of these, get the p-values from across the various cancer types 
  new_pvals <- lapply(unique_targets, function(targ) {
    print(targ) 
    pvals <- lapply(1:length(list_of_results), function(i) {
      results <- list_of_results[[i]]
      feature <- names(list_of_results)[i]
      if(targ %fin% results$gene) {
        results_targ <- results[results$gene == targ, ]
        if(nrow(results_targ) == 1) {return(data.frame("pval" = results_targ$pval,
                                                       "feature" = feature,
                                                       "gene" = targ))}
        pvals <- results_targ$pval
        # Do the correction using these p-values
        pval_new <- metap::sumlog(pvals, log.p = F)$p
        return(data.frame("pval" = pval_new, "feature" = feature, 
                          "gene" = targ))
      } else {return(NA)}
    })

    if(length(pvals) == 1) {
      if(is.na(pvals)) {return(NA)}
    }
    pvals <- pvals[!is.na(pvals)]
    print(pvals)
    pval_df <- do.call(rbind, pvals)
    return(pval_df)
  })
  new_pvals <- new_pvals[!is.na(new_pvals)]
  #print(new_pvals)
  new_df <- do.call(rbind, new_pvals)
  
  # Create a new data frame with all targets and their p-values
  #new_df <- data.frame("gene" = unique_targets, "pval" = new_pvals)
  
  # Multiple hypothesis testing correction
  new_df$qval <- qvalue(new_df$pval)$qvalues
  #new_df$p.adj <- p.adjust(new_df$pval, "BH")
  
  return(new_df)
}


# Import DrugBank data (to limit analysis to genes that have at least one 
# associated entry in DrugBank). File can be accessed via the DrugBank website 
# with a request using a University email (drugbank.com)
drugbank_data <- read.csv("D:/DrugBank/all_drugbank_uniprot_links.csv", 
                          header = T, check.names = F)
targets_with_drugs <- unique(drugbank_data$`UniProt ID`)
targets_with_drugs_gns <- unlist(mclapply(targets_with_drugs, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 
                                        'external_gene_name'])), collapse = ";")))
targets_with_drugs_gns <- unique(targets_with_drugs_gns[targets_with_drugs_gns != ""])  # 2814

# Alternatively, limit to just Recon3D metabolic genes
recon3d_genes_df <- read.csv(paste0(PATH, "metabolic_targets.csv"), 
                             header = T, check.names = F)
recon3d_genes <- unique(unlist(lapply(recon3d_genes_df$ensg, function(x) 
  paste(unique(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == x, 
                                 'external_gene_name']), collapse = ";"))))

# Define a given driver and subset Dyscovr output to those cancer types with 
# that driver
driver_name <- "KRAS"
perCancer_driver <- lapply(perCancer, function(x) {
  if(driver_name %fin% x$R_i.name) return(x)
  else {return(NA)}
})
perCancer_driver <- perCancer_driver[!is.na(perCancer_driver)]

# Get pan-cancer hits at a given threshold
pc_qval_thres <- 0.2
pc_hits_driver <- pc_allGenes[(pc_allGenes$R_i.name == driver_name) &
                                (pc_allGenes$q.value < pc_qval_thres), 
                              'T_k.name']

# Call function
synthetic_lethals <- get_synthetic_lethals(driver_name, perCancer_driver, 
                                           0.2, T, pc_hits_driver, T, 
                                           targets_with_drugs_gns, F, 
                                           crispr, expression,
                                           cell_line_cancer_type_mapping, 
                                           "Regression",
                                           "oncogene")
synthetic_lethals <- synthetic_lethals[order(synthetic_lethals$padj),]


# Run on cancer types individually and recombine their p-values
synthetic_lethals_comb <- recombine_pvals(synthetic_lethals)
synthetic_lethals_comb <- synthetic_lethals_comb[order(synthetic_lethals_comb$p.adj),]


# Perform a cross-cancer LM analysis that includes covariates for cancer type
lm_results <- perform_lm_analysis_crosscancers(crispr, expression_qn, cnas, driver_name, 
                                               perCancer_driver, 0.2, T, 
                                               pc_hits_driver, F, NA, 
                                               F, cell_line_cancer_type_mapping,
                                               c('mutation', 'expression'))
lm_results <- lm_results[!grepl("cancer_type", lm_results$term),]
lm_results <- lm_results[!grepl("cna_val", lm_results$term),]


# Visualize overlap between hits found from various sources
sig_hits <- list("PAN_LM" = unique(lm_results[lm_results$padj < 0.2, 'gene']),
                 "PER_LM" = unique(synthetic_lethals[synthetic_lethals$padj < 0.2, 'gene']),
                 "PER_RECOMB" = unique(synthetic_lethals_comb[synthetic_lethals_comb$p.adj < 0.2, 'gene']))
plt <- venn.diagram(sig_hits, category.names = names(sig_hits), 
                    filename = NULL, output = TRUE, lwd = 2, lty = 'blank', 
                    fill = c("#0072B5FF", "#BC3C29FF", "#20854EFF"), #, "#E18727FF", "#7876B1FF""#FFDC91FF", "#EE4C97FF"), 
                    cex = 2, fontface = "bold", fontfamily = "sans", cat.cex = 2, 
                    cat.fontface = "bold", cat.fontfamily = "sans")
grid::grid.draw(plt)

# Display statistical significance of enrichment in STRING targets (See figure1.R)
# NOTE: Prior to doing enrichment analysis, adjust the 'universe' of tests (e.g. 
# intersect target set with a) pan- or per-cancer hits that were tested, b) only
# genes tested in DepMap CRISPR and/or expression analyses, c) only genes present
# in any other outside sets used, such as Recon3D or DrugBank)
# e.g. kras_string_nw_targs_pc_depmap <- intersect(kras_string_nw_targs, intersect(crispr$gene, pc_hits_driver)) 
# per_hits <- unique(unlist(lapply(perCancer_driver, function(x) unique(x[x$q.value < 0.2, 'T_k.name']))))
# kras_string_nw_targs_per_depmap <- intersect(kras_string_nw_targs, intersect(crispr$gene, per_hits)) 
compute_statistical_enrichment(lm_results, kras_string_nw_targs, "both", 0.01)
# ...

input <- melt(as.data.frame(list("TP53" = -log10(1.09E-11), "PIK3CA" = -log10(9.53E-14), 
                                 "KRAS" = -log10(1E-99)))) #"IDH1" = -log10())))
colnames(input) <- c("Driver", "negLog10pval")

ggplot(input, aes(y = negLog10pval, fill = Driver, x = reorder(Driver, negLog10pval, mean))) + 
  geom_bar(position = "dodge", width = 0.95, stat = "identity", show.legend = FALSE, color = "black") + 
  scale_fill_manual(values = c("#0072B5FF", "#BC3C29FF", "#20854EFF", "#FFDC91FF")) + #scale_color_nejm() +
  xlab("Driver") + ylab(paste0("-log10(pval) of K-S Enrichment,\n STRING Neighbors")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black") + coord_flip() + #theme_minimal() +
  scale_y_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text(face="bold", size = 14), 
        axis.title=element_text(size=16, face="bold"), panel.grid.major = element_blank(),
        panel.background = element_rect(fill = 'white'))


# For enrichment analyses, optionally limit to only the first (most significant)
# occurrence of each target gene
seen_genes <- c()
rows_to_keep <- c()
#synthetic_lethals <- synthetic_lethals[order(synthetic_lethals$pvalue),] #option to order by original pval, rather than adjusted
for (i in 1:nrow(synthetic_lethals)) {
  gene <- synthetic_lethals[i, 'gene']
  if(!(gene %fin% seen_genes)) {
    seen_genes <- c(seen_genes, gene)
    rows_to_keep <- c(rows_to_keep, i)
  }
}
synthetic_lethals_distinct <- synthetic_lethals[rows_to_keep,]

#' For enrichment analyses, optionally re-rank genes by the average of their
#' rankings for each feature
#' @param synthetic_lethals an output data frame from synthetic lethal analysis
rerank_by_average_of_features <- function(synthetic_lethals) {
  unique_genes <- unique(synthetic_lethals$gene)
  gene_avgs <- unlist(lapply(unique_genes, function(g) {
    row_nums <- which(synthetic_lethals$gene == g)
    return(mean(row_nums))
  }))
  gene_avgs_df <- data.frame("index" = gene_avgs, "gene" = unique_genes)
  synthetic_lethals <- merge(synthetic_lethals, gene_avgs_df, by = "gene")
  synthetic_lethals <- synthetic_lethals[order(synthetic_lethals$index),]
  synthetic_lethals <- synthetic_lethals[match(unique(synthetic_lethals$gene), 
                                               synthetic_lethals$gene),]
  return(synthetic_lethals)
}

# Alternatively, recombine p-values for targets across modalities
#' Using the metap package, recombine the p-values for target gene results for
#' both driver expression and driver mutation status
#' @param results a data frame of results from get_synthetic_lethals
#' without subsetting to significant hits- change p-value threshold to 1)
recombine_pvals_byterm <- function(results) {
  unique_terms <- unique(results$term)
  unique_targets <- unique(results$gene)
  
  # For each of these, get the p-values from across the various cancer types 
  new_pvals <- lapply(unique_targets, function(targ) {
    print(targ) 
    results_targ <- results[results$gene == targ, ]
    if(nrow(results_targ) == 1) {
      return(data.frame("pval" = results_targ$pval, "gene" = targ))}
    pvals <- results_targ$pval
    
    # Do the correction using these p-values
    pval_new <- metap::sumlog(pvals, log.p = F)$p
    return(data.frame("pval" = pval_new, "gene" = targ))
  })
  
  new_df <- do.call(rbind, new_pvals)

  # Multiple hypothesis testing correction
  new_df$qval <- qvalue(new_df$pval)$qvalues
  #new_df$p.adj <- p.adjust(new_df$pval, "BH")
  
  return(new_df)
}

lm_results_comb <- recombine_pvals_byterm(lm_results)
lm_results_comb <- lm_results_comb[order(lm_results_comb$qval),]

#' Limit recombined version to only those that are in the correct direction
#' (both expression and mutation, or just mutation)
#' @param lm_res a results table from DepMap linear regression analysis
#' @param lm_res_fisher a results table from DepMap LM analysis that has been
#' recombined using Fisher's method (one p-value per target gene)
#' @param onco_or_ts either "oncogene" or "ts" to indicate whether the driver
#' gene is an oncogene or a tumor suppressor gene
#' @param type either "both", to signify we want to limit by directionality of
#' both the mutation and the expression terms, or "mut" to signify we want to 
#' limit by the directionality of the mutation term only
limit_by_directionality <- function(lm_res, lm_res_fisher, onco_or_ts, type) {
  rows_to_keep <- unlist(lapply(1:nrow(lm_res_fisher), function(i) {
    t <- lm_res_fisher[i, 'gene']
    
    # Get the directionality of target t
    dir_mut <- NA
    if(onco_or_ts == "onco") {
      dir_mut <- ifelse(lm_res[(lm_res$gene == t) & (lm_res$term == 'driver_mut_status1'), 
                               'stat'] > 0, 1, -1)
    } else if(onco_or_ts == "ts") {
      dir_mut <- ifelse(lm_res[(lm_res$gene == t) & (lm_res$term == 'driver_mut_status1'), 
                               'stat'] < 0, 1, -1)
    } else {print("Must supply either onco or ts status.")}

    if(type == "both") {
      dir_exp <- ifelse(lm_res[lm_res$gene == t & (lm_res$term == 'expression_val'), 
                               'stat'] > 0, 1, -1)
      if((dir_mut == 1) & (dir_exp == 1)) {return(i)}
      else{return(NA)}
    } else if (type == "mut") {
      if(dir_mut == 1) {return(i)}
      else{return(NA)}
    } else {
      print(paste("Type not implemented:", type))
      return(NA)
    }
  }))
  rows_to_keep <- rows_to_keep[!is.na(rows_to_keep)]
  
  lm_res_fisher_sub <- lm_res_fisher[rows_to_keep,]
  return(lm_res_fisher_sub)
}

# Call function
lm_res_kras_panq0.2_perq0.2_fishersp_subB <- limit_by_directionality(lm_res_kras_panq0.2_perq0.2,
                                                                    lm_res_kras_panq0.2_perq0.2_fishersp,
                                                                    "onco", "both")
lm_res_tp53_panq0.2_perq0.2_fishersp_subB <- limit_by_directionality(lm_res_tp53_panq0.2_perq0.2,
                                                                     lm_res_tp53_panq0.2_perq0.2_fishersp,
                                                                     "ts", "both")

# Validate on experimentally verified, gold-standard SLs from
# https://www.nature.com/articles/s41467-018-04647-1#citeas; Suppl. Information 1
gold_standard_sls <- read.csv(paste0(PATH, "Validation_Files/GoldStandard_SLs_10.1038s41467-018-04647-1.csv"),
                              header = T, check.names = F)
gold_standard_sls <- gold_standard_sls[gold_standard_sls$SL == 1,]

# Keep only gene pairs where both are present in our analysis
gold_standard_sls <- gold_standard_sls[(gold_standard_sls$gene1 %fin% crispr$gene) & 
                                         (gold_standard_sls$gene2 %fin% crispr$gene),]

# Keep only gene pairs where at least one of the genes in the pair is a driver
# gene of interest
perCancer_drivers <- unique(unlist(lapply(perCancer, function(x) unique(x$R_i.name))))
gold_standard_sls_drivers <- gold_standard_sls[(gold_standard_sls$gene1 %fin% perCancer_drivers) |
                                                 (gold_standard_sls$gene2 %fin% perCancer_drivers),]

driver <- "KRAS"
gold_standard_sls_driver <- unique(unlist(
  gold_standard_sls_drivers[(gold_standard_sls_drivers$gene1 == driver) | 
                              (gold_standard_sls_drivers$gene2 == driver), 
                            c('gene1', 'gene2')]))
gold_standard_sls_driver <- gold_standard_sls_driver[gold_standard_sls_driver != driver]
compute_statistical_enrichment(synthetic_lethals, intersect(gold_standard_sls_driver, synthetic_lethals$gene), 
                               "both", 0.01, NA)

# Try this for pairings from SynLethDB as well (human): https://synlethdb.sist.shanghaitech.edu.cn/#/download
gold_standard_sls_synlethdb <- read.csv(paste0(PATH, "Validation_Files/Human_SL_SynLethDB.csv"),
                              header = T, check.names = F)
# Optional thresholding for confidence 
#gold_standard_sls_synlethdb <- gold_standard_sls_synlethdb[gold_standard_sls_synlethdb$r.statistic_score > 0.5,]
# Optional thresholding to keep only experimental sources and text mining
experimental_sources <- c("High Throughput", "GenomeRNAi", "Text Mining", "RNAi Screen",
                          "Synlethality;GenomeRNAi", "CRISPR/CRISPRi", "Drug Screen",
                          "Synlethality;Text Mining", "Low Throughput", "Text Mining;Synlethality",
                          "Decipher;Text Mining", "GenomeRNAi;Text Mining", "Text Mining;Daisy",
                          "High Throughput|Low Throughput", "GenomeRNAi;Decipher")
experimental_sources_restrict <- c("High Throughput", "GenomeRNAi", "RNAi Screen", 
                                   "Synlethality;GenomeRNAi", "CRISPR/CRISPRi", 
                                   "Drug Screen", "Low Throughput", "GenomeRNAi;Text Mining",
                                   "High Throughput|Low Throughput", "GenomeRNAi;Decipher")
gold_standard_sls_synlethdb_exp <- gold_standard_sls_synlethdb[gold_standard_sls_synlethdb$r.source %fin%
                                                                 experimental_sources,]  # eliminates ~10K entries
gold_standard_sls_synlethdb_drivers <- gold_standard_sls_synlethdb[
  (gold_standard_sls_synlethdb$n1.name %fin% perCancer_drivers) |
    (gold_standard_sls_synlethdb$n2.name %fin% perCancer_drivers),]
gold_standard_sls_synlethdb_driver <- unique(unlist(
  gold_standard_sls_synlethdb_drivers[
    (gold_standard_sls_synlethdb_drivers$n1.name == driver) | 
      (gold_standard_sls_synlethdb_drivers$n2.name == driver), 
                            c('n1.name', 'n2.name')]))
gold_standard_sls_synlethdb_driver <- gold_standard_sls_synlethdb_driver[
  gold_standard_sls_synlethdb_driver != driver]  # 2446 for KRAS, 1357 at 0.5 thres., 263 at 0.7 thres, 5 at 0.8 thres.

gene_universe <- unique(synthetic_lethals[, 'gene'])
compute_statistical_enrichment(synthetic_lethals, intersect(gold_standard_sls_synlethdb_driver, 
                                                            gene_universe), 
                               "both", 0.2, NA)
# Other options - SLOAD is cancer-specific, but you need to email them to get their results

#' Creates a Venn diagram with the overlapping significant hits from our SL pipeline
#' and the SynLethDB database
#' @param synleth_df a data frame containing synthetic lethals for a given driver
#'from our pipeline 
#' @param gene_universe the universe of tested genes from Dyscovr, prior to 
#' directionality filtering
#' @param gold_std_sls_driver a list of gold standard SLs for driver from SynLethDB
#' @param goi the name of the gene-of-interest, for labeling purposes
#' @param qval_thres a q-value threshold for significance
plot_hit_overlap_synleth <- function(synleth_df, gene_universe,
                                     gold_std_sls_driver, goi, qval_thres) {
  
  # Get significant SLs from our pipeline
  sig_hits <- synleth_df[synleth_df$qvalue < qval_thres, 'gene']
  print(length(sig_hits))
  
  # Keep only gold standard hits that were tested in our pipeline
  gold_std_sls_driver <- intersect(gene_universe, gold_std_sls_driver)
  
  # Calculate whether overlap is significant using an approx. of an exact hypergeometric test
  len_hits1 <- length(sig_hits)
  len_hits2 <- length(gold_std_sls_driver)
  len_ol <- length(intersect(sig_hits, gold_std_sls_driver))
  n <- length(gene_universe)
  
  print(setdiff(gold_std_sls_driver, sig_hits))
  
  # n is size of urn; hits1 is the number of genes from set 1. Say we were to draw
  # the number of balls the size of hits2 (# of genes from set 2), w/o replacement. 
  # M of these balls are also found in hits2. What are the odds that M >= len_ol?
  phyper_res <- phyper(len_ol-1, len_hits1, n-len_hits1, len_hits2, 
                       lower.tail = FALSE, log.p = FALSE)
  print(paste("p-value:", phyper_res))
  
  # Plot Venn Diagram
  #group.colors <- c("#0072B5FF", "#FFDC91FF")
  #names(group.colors) <- names(c("Dyscovr", "SynLethDB"))
  
  #ggVennDiagram(sig_hits_list, edge_size = 5, edge_lty = "solid", category.names = names(sig_hits_list),
  #                  set_size = 8, label_size = 8, label_alpha = 0,
  #                 label = "count", set_color = c("#E18727FF", "#0072B5FF")) + #rep("black", times = length(sig_hits_list))) +
  #theme(legend.position = "none") + #scale_fill_manual(values = group.colors) +
  #scale_fill_gradient(low="white", high = "white") + #labs(fill = paste("# Hits, q<", qval_thres)) +
  #scale_fill_manual(values = c("#FFDC91FF", "#0072B5FF")) + # "#E18727FF"
  #scale_color_manual(values = rep("black", times = length(sig_hits_list))) +
  #scale_color_manual(values = c("#E18727FF", "#0072B5FF")) +
  #labs(title = paste("q <", paste0(qval_thres, paste(", HG p =", paste0(format(phyper_res, digits = 4, scientific = TRUE), "\n")))))
  #text(labels = paste("HG p=", format(phyper_res, digits = 4, scientific = TRUE)),
  #x = 0.1, y = 0.1)
  
  sig_hits_list <- list("SynLethDB" = gold_std_sls_driver, "Dyscovr" = sig_hits)
  VennDiag <- euler(sig_hits_list)
  plot(VennDiag, counts = TRUE, font=2, cex=6, cex.main=1, alpha=0.9, quantities=T,
       fill=c("#20854EFF", "#FFDC91FF"), col="darkgray", lwd=3,
       main = paste("q <", paste0(qval_thres, paste0(", HG p=", 
                                                     format(phyper_res, digits = 4, 
                                                            scientific = T)))))
}

# Call function
gene_universe <- unique(kras_synleth[, 'gene'])
plot_hit_overlap_synleth(kras_synleth_dirsub, gene_universe,
                         gold_standard_sls_synlethdb_0.7_driver,
                         "KRAS", 0.2)


# Add rankings for mutation and expression terms 
rankings <- c()
curr_mut <- 1
curr_expr <- 1
for (i in 1:nrow(pik3ca_lm_res_panq0.2_perq1_wtargcna)) {
  term <- pik3ca_lm_res_panq0.2_perq1_wtargcna[i, 'term']
  if(term == 'expression_val') {
    rankings <- c(rankings, curr_expr)
    curr_expr <- curr_expr + 1
  } else if (term == 'driver_mut_status1') {
    rankings <- c(rankings, curr_mut)
    curr_mut <- curr_mut + 1
  } else {rankings <- c(rankings, NA)}
}
pik3ca_lm_res_panq0.2_perq1_wtargcna$ranking <- rankings

# Create a star-like network visualization for a given driver gene based on 
# functional annotations (e.g. STRING, HumanBase) and cluster based on pathway

# Useful tutorial: https://mr.schochastics.net/material/netVizR/

#' Generate results graph display for a driver of interest 
#' This graphical display takes confidence interaction relationships from 
#' a network like STRING or HumanBase, and displays these relationships
#' as a graph. The graph will only include hits with a q-value below the
#' given threshold. The edges and genes are colored by whether activation or 
#' repression is predicted by Dyscovr when the driver is mutated. 
#' @param synthetic_lethals an output DF from synthetic lethal analysis
#' @param goi the name of the driver gene of interest
#' @param qval_thres a q-value threshold, below which a target-driver
#' pairing is considered significant
#' @param network a table with the network information from either 
#' HumanBase or STRING
#' @param network_label either "STRING" or "HumanBase" to denote 
#' what type of network this is
#' @param conf_thres a confidence threshold for a connection between two 
#' genes in the given network
#' @param top_n if not NA, gives a value of the top N most significant targets 
#' to use
#' @param label_genes if TRUE, labels individual genes and colors by cluster; if 
#' FALSE, only colors by cluster and wraps clusters
#' @param num_clust an option to specify the number of clusters, default is 15
#' @param disconn_comp_size remove disconnected components that are smaller than
#' the given size
#' @param pathway_label_list a list of all pathways of interest, with each pathway
#' entry containing a vector of gene names for all genes in that pathway. All genes
#' in the given 'synthetic_lethals' data frame should be contained in one and 
#' only one of the pathways (optional)
create_graphical_representation <- function(synthetic_lethals, goi, qval_thres, 
                                            network, network_label, conf_thres, 
                                            top_n, label_genes, pathway_label_list,
                                            num_clust = 15, disconn_comp_size = 5) {
  
  # Subset the synthetic lethal DF to the given q-value threshold/ ranking
  synthetic_lethals_sub <- synthetic_lethals[synthetic_lethals$qval < qval_thres,]
  # synthetic_lethals_sub <- synthetic_lethals[synthetic_lethals$padj < qval_thres,]
  if(!is.na(top_n)) {synthetic_lethals_sub <- synthetic_lethals_sub[1:top_n, ]}
  nodes <-  c(goi, unique(synthetic_lethals_sub$gene))
  
  network_model_targs <- NA
  # Subset the network to the significant targets, plus the GOI
  if(network_label == "STRING") {
    network_model_targs <- network[(network$node1_name %fin% nodes) & 
                                     (network$node2_name %fin% nodes),]
  } else if (network_label == "HumanBase") {
    network_model_targs <- network[(network$node1_name %fin% nodes) & 
                                     (network$node2_name %fin% nodes),]
    
  } else {
    print("Error in given network label. Only implemented for STRING and 
          HumanBase networks.")
  }
  # Remove cases where nodes 1 and 2 are the same, or rows where the values in 
  # the row are duplicated in a different order. Divide all confidence scores by 
  # 1000 for STRING
  network_model_targs <- subset(network_model_targs, node1_name != node2_name & 
                                  !duplicated(cbind(pmin(node1_name, node2_name), 
                                                    pmax(node1_name, node2_name))))
  if(network_label == "STRING") {
    network_model_targs$combined_score <- network_model_targs$combined_score / 1000
  }
  print(head(network_model_targs))
  
  # For each of the targets, get its a) pathway that contains it, 
  # b) confidence score in the given network, in relation to the driver 
  edge_rows <- lapply(1:nrow(network_model_targs), function(i) {
    conf <- network_model_targs[i, 'combined_score']   
    if(conf > conf_thres) {
      g1 <- network_model_targs[i, 'node1_name']
      g2 <- network_model_targs[i, 'node2_name']
      return(data.frame("src" = g1, "target" = g2, 
                        "conf" = as.numeric(conf)))
    } else {return(NA)}
  })
  edge_table <- do.call(rbind, edge_rows)
  colnames(edge_table) <- c("src", "target", "conf")
  edge_table <- na.omit(edge_table)

  unique_nodes <- unique(c(edge_table$src, edge_table$target))
  node_table <- as.data.frame(unique_nodes)
  colnames(node_table) <- "nodes"
  
  # Create a mapping from pathways to integers, if we have a priori pathways
  # Otherwise, we are just clustering based on connectivity and calling clusters later
  if(length(pathway_label_list) > 1) {
    pathway_label_list_int <- data.frame("pw" = names(pathway_label_list),
                                         "marker" = 1:length(pathway_label_list))
    # Add the pathway to the node table
    node_table$pw <- as.factor(unlist(lapply(unique_nodes, function(n) {
      pw_int <- NA
      if(n == goi) {pw_int <- nrow(pathway_label_list_int)+1}  # Driver gene is excluded from pathway clusters
      else {
        pw <- names(which(sapply(pathway_label_list, function(X) n %fin% X)))
        pw <- pw[!is.na(pw)][1]  #current hack to just take the first pw
        if(length(pw) == 0) {return(NA)} 
        if(is.na(pw)) {return(NA)}
        else {
          pw_int <- as.integer(unlist(pathway_label_list_int[
            pathway_label_list_int$pw == pw, "marker"]))
        }
      }
      return(pw_int)
    })))
  }
  node_table <- na.omit(node_table)
  
  node_table$node_size <- unlist(lapply(node_table$nodes, function(n) {
    if(n == goi) {return(2)}
    else {return(1)}
  }))
  print(node_table)
  
  # Ensure that node and edge tables contain the same genes
  intersect_genes <- intersect(node_table$nodes, 
                               unique(c(edge_table$src, edge_table$target)))
  node_table <- node_table[node_table$nodes %fin% intersect_genes,]
  edge_table <- edge_table[(edge_table$src %fin% intersect_genes) & 
                             (edge_table$target %fin% intersect_genes),]
  
  #return(list(edge_table, node_table))

  # Convert to ggraph format
  network_table <- tbl_graph(nodes = node_table, edges = edge_table, directed = F)
  #network_table <- igraph::graph_from_data_frame(edge_table, vertices = node_table) %>% as_tbl_graph()

  # If we do not have predefined pathways, group nodes based on community structure
  # Options for algorithms to use: https://tidygraph.data-imaginist.com/reference/group_graph.html
  if(length(pathway_label_list) < 2) {
    #network_table <- network_table %>% mutate(group = group_components("strong"))  # community structure algorithm
    #network_table <- network_table %>% mutate(group = group_edge_betweenness(
    #  directed = F, weights = (1-conf), n_groups = num_clust))  # community structure algorithm
    # Note that weights are interpreted as edge lengths during betweenness calculations 
    # (high weight = weak connection, low weight = strong connection)
    network_table <- network_table %>% mutate(group = group_fast_greedy(weights = (1-conf), 
        n_groups = num_clust))
    
    # Modify the table to make the driver gene its own group with a unique color
    network_table_df <- network_table %>% activate(nodes) %>% as_data_frame()
    network_table_df[network_table_df$nodes == goi, 'group'] <- max(
      as.integer(network_table_df$group)) + 1
    network_table_df$group <- as.factor(as.integer(network_table_df$group))
    print(head(network_table_df$group))

    # Put this back into tbl_graph form
    network_table <- network_table %>% activate(nodes) %>% 
      mutate(group = network_table_df$group)
  }
  
  # Create ggraph
  if(length(pathway_label_list) > 1) {
    driver_targ_network <- ggraph(network_table, layout = "fr") + 
                                  #layout = "centrality", cent = graph.strength(network_table)) +   # other options: stress, fr
                                  #layout = "manual", x = bb$xy[,1], y = bb$xy[,2]) + 
      geom_edge_link0(edge_colour = "darkgray", edge_width = 1) +  
      geom_node_point(aes(fill = pw, size = node_size), shape = 21) +    
      scale_size(range = c(5,10), guide = 'none') + 
      geom_node_label(aes(label = nodes), repel = T, show.legend = F, 
                      label.size = NA, label.padding = 0.05) +
      #geom_mark_hull(mapping = aes(x, y, group = pw, fill = pw, label = pw), 
      #               concavity = 4, expand = unit(2, 'mm'), alpha = 0.25) +
      theme_graph() +
      scale_fill_nejm(labels = names(pathway_label_list), name = "Pathway") +
      scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
      theme(legend.title = element_text(face = "bold", size = 14),
            legend.text = element_text(size = 12), legend.position="bottom") 
    
    # Render the network
    show(driver_targ_network)
    
  } else {
    # List of 20 distinct colors:
    # https://sashamaps.net/docs/resources/20-colors/
    distinct_colors <- c('#4363d8', '#e6194B', '#3cb44b', '#ffe119', '#f58231', 
                         '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', 
                         "#ffffff",
                         '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', 
                         '#aaffc3', '#808000', '#ffd8b1', "slategray1", '#000075',
                         "cyan2", "#7876B1FF", "#BC3C29FF", '#a9a9a9')

    # Option to add pathway annotations using GSEA (right now just taking the 
    # most statistically significantly enriched KEGG pathway)
    if(!label_genes) {
      network_table <- add_pathway_annotations(network_table, edge_table,
                                               synthetic_lethals)
    }
    print(network_table %>% activate(nodes) %>% as_data_frame())
    
    # Remove small disconnected components
    network_table_sub <- induced_subgraph(
      network_table, V(network_table)[ave(1:vcount(network_table), 
                                          membership(components(network_table)), 
                                          FUN = length) > disconn_comp_size])
    print(igraph::as_data_frame(network_table_sub, what = "vertices"))
    
    # Use 'backbone structure' to separate out what might otherwise be a hairball
    bb <- layout_as_backbone(network_table_sub, keep = conf_thres)
    E(network_table_sub)$col <- F
    E(network_table_sub)$col[bb$backbone] <- T
    print(igraph::as_data_frame(network_table_sub, what = "vertices"))

    if(label_genes) {
      driver_targ_network <- ggraph(network_table_sub, #, layout = "fr") + 
        #layout = "centrality", cent = graph.strength(network_table)) +   # other options: stress, fr
        layout = "manual", x = bb$xy[,1], y = bb$xy[,2]) + 
        geom_edge_link0(edge_colour = "darkgray", edge_width = 1) +  
        geom_node_point(aes(fill = group, size = node_size), shape = 21) +    
        scale_size(range = c(5,10), guide = 'none') + 
        geom_node_label(aes(label = nodes), repel = T, show.legend = F, 
                        label.size = NA, label.padding = 0.05) +
        #geom_mark_hull(aes(x, y, group = group, fill = group, label = group), 
        #               concavity = 4, expand = unit(2, 'mm'), alpha = 0.25) +
        theme_graph() + #scale_fill_nejm() + 
        scale_fill_manual(values = distinct_colors) + 
        scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
        theme(legend.position = "none")
      
    } else {
      driver_targ_network <- ggraph(network_table_sub, #, layout = "fr") + 
        #layout = "centrality", cent = graph.strength(network_table)) +   # other options: stress, fr
        layout = "manual", x = bb$xy[,1], y = bb$xy[,2]) + 
        geom_edge_link0(aes(col = col), edge_colour = "darkgray", edge_width = 1) +  
        geom_node_point(aes(fill = group, size = node_size), shape = 21) +    
        scale_size(range = c(5,10), guide = 'none') + 
        #geom_node_label(aes(label = pw), repel = T, show.legend = F, 
        #                label.size = NA, label.padding = 0.05) +
        geom_mark_hull(aes(x, y, group = group, fill = group, label = pw), 
                       expand = unit(4, 'mm'), concavity = 4, alpha = 0.25,
                       label.fontsize = 10) +
                       #alpha = 0.25) + 
                       #con.type = "none") +
        theme_graph() + #scale_fill_nejm() + 
        scale_fill_manual(values = distinct_colors) + 
        scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
        theme(legend.position = "none")
    }

    # Render the network
    show(driver_targ_network)
  }
  return(network_table %>% activate(nodes) %>% as_data_frame())
}


#' Helper function to perform GSEA on the genes in each cluster separately to get
#' the single most statistically significantly enriched term, which we will use
#' to label each cluster. Currently only implemented for KEGG pathways. Already
#' subsetted to only connected components.
#' @param network_table_sub a tbl_graph object with node and edge components
#' @param edge_table an edge table for creating an updated tbl_graph object
#' @param synthetic_lethals an output DF from synthetic lethal analysis
add_pathway_annotations <- function(network_table, edge_table, synthetic_lethals) {
  # Get the nodes as a data frame
  network_table_df <- network_table %>% activate(nodes) %>% as_data_frame()
  
  # Perform GSEA separately for each group
  max_pw_grp <- max(as.integer(network_table_df$group))
  group_dfs <- lapply(1:(max_pw_grp-1), function(grp) {
    # Extract the nodes for this group
    network_table_grp <- network_table_df[network_table_df$group == grp, ]
    genes <- network_table_grp$nodes
    genes_df <- data.frame("gene" = genes)
    genes_df$negLog10pval <- unlist(lapply(genes_df$gene, function(g)
      -log10(synthetic_lethals[synthetic_lethals$gene == g, 'pval'])))
    print(head(genes_df))
    
    mapping <- as.data.frame(bitr(genes, fromType = "SYMBOL", 
                                  toType = "UNIPROT", OrgDb=org.Hs.eg.db, drop = T))
    if(length(mapping) < 2) {
      network_table_grp$pw <- rep("None", times = nrow(network_table_grp))
      return(network_table_grp)
    }
    mapping_kegg <- as.data.frame(bitr_kegg(mapping$UNIPROT, fromType = "uniprot", 
                                            toType = "kegg", drop = T, 
                                            organism = "hsa"))
    colnames(mapping) <- c("gene", "uniprot")
    colnames(mapping_kegg) <- c("uniprot", "kegg")
    
    mapping <- merge(mapping, mapping_kegg, by = "uniprot")
    print(head(mapping))
    
    results_table <- merge(genes_df, mapping, all=T, by="gene")
    results_table <- distinct(results_table[, c("gene", "kegg", "negLog10pval")])
    print(head(results_table))
    
    results_table <- results_table[order(results_table$negLog10pval, decreasing = T),]
    betas <- results_table$negLog10pval
    names(betas) <- results_table$kegg
    print(betas)
    
    enrichment_res <- character(0)
    tryCatch({
      enrichment_res <- gseKEGG(betas, pvalueCutoff = 1, pAdjustMethod = "BH", 
                                keyType = "kegg")
    }, error = function(cond) {
      print(cond)
    })
    
    # Run enrichment analysis using overenrichment analysis
    #mapping <- as.data.frame(bitr(genes, fromType = "SYMBOL", toType = "UNIPROT", 
    #                              OrgDb=org.Hs.eg.db, drop = T))
    #mapping_kegg <- as.data.frame(bitr_kegg(mapping$UNIPROT, fromType = "uniprot", 
    #                                        toType = "kegg", drop = T, 
    #                                        organism = "hsa"))

    #enrichment_res <- enrichKEGG(gene = mapping$UNIPROT, organism = 'hsa', 
    #                             pAdjustMethod = "BH", pvalueCutoff = 0.2,
    #                             keyType = "uniprot")

    if(length(enrichment_res) < 2) {
      network_table_grp$pw <- rep("None", times = nrow(network_table_grp))
      return(network_table_grp)
    }
    #print(head(enrichment_res@result))

    top_pw <- ""
    if(nrow(enrichment_res) > 0) {top_pw <- enrichment_res@result[1, 'Description']}
    print(top_pw)
    
    network_table_grp$pw <- rep(top_pw, times = nrow(network_table_grp))
    return(as.data.frame(network_table_grp))
  })

  # Bind data tables together
  network_table_df_pw <- do.call(rbind, group_dfs)
  network_table_df_driver <- network_table_df[network_table_df$group == max_pw_grp,]
  network_table_df_driver$pw <- network_table_df_driver[1, 'nodes']
  network_table_df_pw <- rbind(network_table_df_pw, network_table_df_driver)

  # Factorize the pathway annotations
  network_table_df_pw$pw <- as.character(unlist(network_table_df_pw$pw))
  print(network_table_df_pw)

  # Put this back into tbl_graph form
  network_table <- tbl_graph(nodes = network_table_df_pw, edges = edge_table, 
                             directed = F)
  #network_table <- network_table %>% activate(nodes) %>% 
    #mutate(pw = network_table_df_pw$pw)
  
  return(network_table)
}


# Read back, once completed
string_nw_full <- fread(paste0(PATH, "Network_Data/STRING/string.9606.protein.links.v11.5.namesAdded.txt"))

# Put together a pathway label list; get KEGG gene sets using gage
kegg_gsets <- kegg.gsets(species = "hsa", id.type = 'entrez')
# Convert to gene names
kegg_gsets_gns <- kegg_gsets
kegg_gsets_gns$kg.sets <-  lapply(kegg_gsets$kg.sets, function(pw) {
  mapping <- as.data.frame(bitr(pw, fromType = "ENTREZID", toType = "SYMBOL", 
                                OrgDb=org.Hs.eg.db, drop = T))
  return(mapping$SYMBOL)
})
kegg_gsets_gns_list <- kegg_gsets_gns$kg.sets
# Reformat pathway names to remove hsa
names(kegg_gsets_gns_list) <- unlist(lapply(names(kegg_gsets_gns_list), function(n) {
  n_spl <- unlist(strsplit(n, " ", fixed = T))
  return(paste(n_spl[2:length(n_spl)], collapse = " "))
}))

# Call function
set.seed(5)
create_graphical_representation(pik3ca_lm_res_panq0.2_perq1_wtargcna_comb, 
                                "PIK3CA", 0.2, string_nw_full, "STRING", 0.4, NA,
                                F, kegg_gsets_gns_list)
create_graphical_representation(pik3ca_lm_res_panq0.2_perq0.2_wtargcna_comb, 
                                "PIK3CA", 0.2, string_nw_full, "STRING", 0.4, NA,
                                T, NA)

create_graphical_representation(pik3ca_lm_res_panq0.2_perq0.2_wtargcna_comb, 
                                "PIK3CA", 0.2, string_nw_full, "STRING", 0.4, NA,
                                T, NA, num_clust = 10)


#' Perform GSEA on the whole synthetic lethal table together, sorted by p-value
#' @param results_table a result data frame from synthetic lethality pipeline
#' @param minggssize the minimium gene set size for gseGO
perform_gsea_synleth <- function(results_table, mingssize = 50) {
  
  # Get all of this driver protein's targets and their associated mutation 
  # coefficients, sorted descending order
  #mapping <- as.data.frame(bitr(results_table$gene, fromType = "SYMBOL", 
  #                              toType = "UNIPROT", OrgDb=org.Hs.eg.db, drop = T))
  #mapping_kegg <- as.data.frame(bitr_kegg(mapping$UNIPROT, fromType = "uniprot", 
  #                                        toType = "kegg", drop = T, 
  #                                        organism = "hsa"))
  #colnames(mapping_kegg) <- c("uniprot", "kegg")
  
  mapping <- as.data.frame(bitr(results_table$gene, fromType = "SYMBOL", 
                                toType = "ENTREZID", OrgDb=org.Hs.eg.db, drop = T))
  colnames(mapping) <- c("gene", "entrez")
  

  #mapping <- merge(mapping, mapping_kegg, by = "uniprot")
  print(head(mapping))
  
  results_table <- merge(results_table, mapping, all=T, by="gene")
  #results_table <- distinct(results_table[, c("gene", "kegg", "pval")])
  results_table <- distinct(results_table[, c("gene", "entrez", "pval", "qvalue")])
  results_table$negLog10pval <- unlist(lapply(results_table$pval, function(p) -log10(p)))

  results_table <- results_table[order(results_table$negLog10pval, decreasing = T),]
  #results_table <- results_table[order(results_table$stat, decreasing = T),]
  results_table <- na.omit(results_table)
  
  #betas <- results_table$stat
  betas <- results_table$negLog10pval
  #names(betas) <- results_table$kegg
  names(betas) <- results_table$entrez
  
  #gse.kegg <- gseKEGG(betas, pvalueCutoff = 1, pAdjustMethod = "BH", 
                      #keyType = "kegg")
  #gse.kegg <- enrichKEGG(gene = results_table[results_table$qval < 0.2, 'uniprot'], 
  #                       organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.2,
  #                       keyType = "uniprot")
  #gse.kegg <- gseMKEGG(betas, pvalueCutoff = 1, 
                        #pAdjustMethod = "BH", keyType = "kegg")
  #return(gse.kegg)
  
  gse.go <- gseGO(betas, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, 
                  pAdjustMethod = "BH", keyType = "ENTREZID", ont = "BP",
                  minGSSize = mingssize)# BP = biological process
  gse.go <- setReadable(gse.go, org.Hs.eg.db)
  
  #gse.go <- enrichGO(gene = unlist(results_table[results_table$qvalue < 0.2, 'entrez']),
  #                    universe = results_table$entrez,
  #                    OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
  #                    pAdjustMethod = "BH", pvalueCutoff = 1, readable = T)
  #gse.go <- setReadable(gse.go, org.Hs.eg.db)
  
  return(gse.go)
}

gsea_res_tp53_panq0.2_perq0.2_fisherp_subB <- perform_gsea_synleth(lm_res_tp53_panq0.2_perq0.2_fishersp_subB)
gsea_res_tp53_panq0.2_perq0.2_fisherp_subB <- gsea_res_tp53_panq0.2_perq0.2_fisherp_subB@result
gsea_res_pik3ca_panq0.2_perq0.2_fisherp_subB <- perform_gsea_synleth(lm_res_pik3ca_panq0.2_perq0.2_fishersp_subB)
gsea_res_pik3ca_panq0.2_perq0.2_fisherp_subB <- gsea_res_pik3ca_panq0.2_perq0.2_fisherp_subB@result
gsea_res_kras_panq0.2_perq0.2_fisherp_subB <- perform_gsea_synleth(lm_res_kras_panq0.2_perq0.2_fishersp_subB)
gsea_res_kras_panq0.2_perq0.2_fisherp_subB <- gsea_res_kras_panq0.2_perq0.2_fisherp_subB@result

#' Add the average statistic for expression and mutation and use for GSEA
add_avg_stat <- function(lm_res, lm_res_fisher) {
  lm_res_fisher$stat <- unlist(lapply(lm_res_fisher$gene, function(g) {
    stat_mut <- lm_res[(lm_res$gene == g) & (lm_res$term == 'driver_mut_status1'), 'stat']
    stat_exp <- lm_res[(lm_res$gene == g) & (lm_res$term == 'expression_val'), 'stat']
    avg_stat <- mean(stat_mut, stat_exp)
    return(avg_stat)
  }))
  return(lm_res_fisher)
}
lm_res_tp53_panq0.2_perq0.2_fishersp_subB <- add_avg_stat(lm_res_tp53_panq0.2_perq0.2, 
                                                          lm_res_tp53_panq0.2_perq0.2_fishersp_subB)
lm_res_kras_panq0.2_perq0.2_fishersp_subB <- add_avg_stat(lm_res_kras_panq0.2_perq0.2, 
                                                          lm_res_kras_panq0.2_perq0.2_fishersp_subB)
lm_res_pik3ca_panq0.2_perq0.2_fishersp_subB <- add_avg_stat(lm_res_pik3ca_panq0.2_perq0.2, 
                                                          lm_res_pik3ca_panq0.2_perq0.2_fishersp_subB)


# Correlate drug sensitivity to target gene dependency 
drug_sensitivity <- read.csv(paste0(PATH, "DepMap/Drug_sensitivity_(PRISM_Re
                                    purposing_Primary_Screen).csv"),
                             header = T, check.names = F)
# AUC(CTD^2) -- better quality and more cell lines
drug_sensitivity <- read.csv(paste0(PATH, "DepMap/Drug_sensitivity_AUC_(CTD^2).csv"),
                             header = T, check.names = F)
colnames(drug_sensitivity)[1] <- "depmap_id"

# Annotate file with GOI dependency for each cell line
goi <- "KBTBD2"
crispr_sub <- crispr[crispr$gene == goi, c("depmap_id", "value")]
colnames(crispr_sub) <- c("depmap_id", paste0(goi, "_dependency"))
drug_sensitivity <- merge(crispr_sub, drug_sensitivity, by = "depmap_id")
drug_sensitivity <- drug_sensitivity[,c(1,2,10:ncol(drug_sensitivity))]

#' Calculate correlation between GOI dependency and sensitivity to each drug
#' @param drug_sensitivity a merged DF with drug sensitivity and GOI dependency data
#' @param corr_type either 'spearman' or 'pearson' to indicate the type of correlation
#' @param qval_thres a q-value threshold for significance
#' @param target_label the name of the target gene of interest
calc_corr_dep_sens <- function(drug_sensitivity, corr_type, qval_thres, target_label) {
  dep_vals <- drug_sensitivity[,2]
  print(head(dep_vals))
  
  corr_dfs <- lapply(3:ncol(drug_sensitivity), function(i) {
    drug <- colnames(drug_sensitivity)[i]
    print(drug)
    
    sens_vals <- drug_sensitivity[,i]
    corr_res <- NA
    
    if(tolower(corr_type) == "spearman") {
      corr_res <- cor.test(dep_vals, sens_vals, method = "spearman")
    } else if(tolower(corr_type) == "pearson") {
      corr_res <- cor.test(dep_vals, sens_vals, method = "pearson")
    } else {
      print("Only implemented for Spearman and Pearson. Please choose one.")
    }
    outdf <- data.frame("drug" = drug, "corr" = as.numeric(unlist(corr_res$estimate)), 
                        "pval" = as.numeric(unlist(corr_res$p.value)))
    return(outdf)
  })
  corr_df <- do.call(rbind, corr_dfs)
  
  # MHT correction
  corr_df$qval <- qvalue(corr_df$pval)$qvalues
  corr_df$neglog10qval <- unlist(lapply(corr_df$qval, function(q) -log10(q)))
  #corr_df$neglog10pval <- unlist(lapply(corr_df$pval, function(p) -log10(p)))
  
  corr_df$up_or_down <- unlist(lapply(1:nrow(corr_df), function(i) {
    if(corr_df[i, 'qval'] < qval_thres) {
      beta_sign <- ifelse(corr_df[i, 'corr'] > 0, 1, 0) 
      if(beta_sign == 1) {return("pos")}
      else {return("neg")}
    } else {return("ns")}
  }))
  
  corr_df$drug_label <- unlist(lapply(1:nrow(corr_df), function(i) 
    ifelse(corr_df$qval[i] < (qval_thres*0.3), 
           unlist(strsplit(corr_df$drug[i], " (", fixed = T))[1], NA)))
  print(head(corr_df))
  
  # Create visualization
  cols <- c("pos" = "#BC3C29FF", "neg" = "#0072B5FF", "ns" = "grey") 
  alphas <- c("pos" = 0.95, "neg" = 0.95, "ns" = 0.5)
  
  p <- ggplot(corr_df, aes(x = corr, y = neglog10qval, fill = up_or_down, 
                           alpha = up_or_down, label = drug_label)) + 
    geom_point(shape = 21, size = 3) + 
    xlab(paste(corr_type, "Correlation,", 
               paste(target_label, "CRISPRi Dependency +\n Drug Sensitivity (AUC)"))) + 
    ylab("-log10(q-value)") + theme_minimal() + 
    geom_hline(yintercept = -log10(qval_thres), linetype = "dashed") +
    scale_fill_manual(values = cols) + scale_alpha_manual(values = alphas) + 
    theme(legend.position = "none", 
          axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 12)) + 
    geom_text_repel(min.segment.length = unit(0.1, "lines"), size = 5, 
                    force = 3.5, box.padding = 0.4) 
  print(p)
  
  return(corr_df[order(corr_df$qval),])
}

# Call function
corr_df <- calc_corr_dep_sens(drug_sensitivity, "Spearman", 0.2, goi)


# Visualize the "mutation" and "expression" term q-values for a given target
# as a simple barplot
target <- "KBTBD2"
target_comb_qval <- pik3ca_lm_res_fishers_recombB[pik3ca_lm_res_fishers_recombB$gene == target, 
                                                  'qval']
target_precomb_beta_mut <- pik3ca_lm_res[(pik3ca_lm_res$gene == target) & 
                                           (pik3ca_lm_res$term == 'driver_mut_status1'), 
                                         'stat']
target_precomb_beta_expr <- pik3ca_lm_res[(pik3ca_lm_res$gene == target) & 
                                           (pik3ca_lm_res$term == 'expression_val'), 
                                         'stat']
target_precomb_qval_mut <- pik3ca_lm_res[(pik3ca_lm_res$gene == target) & 
                                           (pik3ca_lm_res$term == 'driver_mut_status1'), 
                                         'qvalue']
target_precomb_qval_expr <- pik3ca_lm_res[(pik3ca_lm_res$gene == target) & 
                                            (pik3ca_lm_res$term == 'expression_val'), 
                                          'qvalue']
input <- data.frame("coeff" = c(target_precomb_beta_mut, target_precomb_beta_expr),
                    "qval" = c(target_precomb_qval_mut, target_precomb_qval_expr),
                    "label" = c("Mutation", "Expression"))
ggplot(input, aes(y = coeff, fill = label, x = reorder(label, coeff, mean))) + 
  geom_bar(position = "dodge", width = 0.95, stat = "identity", show.legend = F, 
           color = "black") + 
  scale_fill_manual(values = c("#0072B5FF", "#BC3C29FF")) + #scale_color_nejm() +
  xlab("Driver Term") + ylab("Estimated Coefficient") + coord_flip() + #theme_minimal() +
  scale_y_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text(face="bold", size = 14), 
        axis.title=element_text(size=16, face="bold"), panel.grid.major = element_blank(),
        panel.background = element_rect(fill = 'white'))


#' Generate results graph display for a driver of interest 
#' This graphical display takes confidence interaction relationships from 
#' a network like STRING or HumanBase, and displays these relationships
#' as a graph. The graph will only include hits with a q-value below the
#' given threshold. The edges and genes are colored by whether activation or 
#' repression is predicted by Dyscovr when the driver is mutated. 
#' @param synthetic_lethals an output DF from synthetic lethal analysis
#' @param dyscovr_out the raw output from Dyscovr, which has directionality predictions
#' @param goi the name of the driver gene of interest
#' @param qval_thres a q-value threshold, below which a target-driver
#' pairing is considered significant
#' @param network a table with the network information from either 
#' HumanBase or STRING
#' @param network_label either "STRING" or "HumanBase" to denote 
#' what type of network this is
#' @param conf_thres a confidence threshold for a connection between two 
#' genes in the given network
#' @param top_n if not NA, gives a value of the top N most significant targets 
#' to use
#' @param disconn_comp_size the size of disconnected components to remove (default = 5)
create_graphical_representation_dir <- function(synthetic_lethals, dyscover_out,
                                                goi, qval_thres, network, 
                                                network_label, conf_thres, top_n,
                                                disconn_comp_size = 5) {
  
  # Subset the master DF to the given q-value threshold
  synthetic_lethals_sub <- synthetic_lethals[(synthetic_lethals$qval < qval_thres),]
  if(!is.na(top_n)) {synthetic_lethals_sub <- synthetic_lethals_sub[1:top_n, ]}
  nodes <-  c(goi, unique(synthetic_lethals_sub$gene))
  
  network_model_targs <- NA
  
  # Subset the network to the significant targets, plus the GOI
  if(network_label == "STRING") {
    network_model_targs <- network[(network$node1_name %fin% nodes) & 
                                     (network$node2_name %fin% nodes),]
  } else if (network_label == "HumanBase") {
    network$node1_name <- rep(goi, times = nrow(network))
    colnames(network)[which(colnames(network) == "SYMBOL")] <- "node2_name"
    colnames(network)[which(colnames(network) == "SCORE")] <- "combined_score"
    print(head(network))
    network_model_targs <- network[(network$node1_name %fin% nodes) & 
                                     (network$node2_name %fin% nodes),]
    print(head(network_model_targs))
    
  } else {
    print("Error in given network label. Only implemented for STRING and 
          HumanBase networks.")
  }
  
  # Remove cases where nodes 1 and 2 are the same, or rows where the values in 
  # the row are duplicated in a different order. Divide all confidence scores by 
  # 1000 for STRING
  if(network_label == "STRING") {
    network_model_targs <- subset(network_model_targs, node1_name != node2_name & 
                                    !duplicated(cbind(pmin(node1_name, node2_name), 
                                                      pmax(node1_name, node2_name))))
    network_model_targs$combined_score <- network_model_targs$combined_score / 1000
  }
  print(head(network_model_targs))
  
  # For each of the targets, get its a) directionality (up- or down-regulated), 
  # b) confidence score in the given network, in relation to the driver 
  edge_rows <- lapply(1:nrow(network_model_targs), function(i) {
    conf <- network_model_targs[i, 'combined_score']   
    if(conf > conf_thres) {
      g1 <- network_model_targs[i, 'node1_name']
      g2 <- network_model_targs[i, 'node2_name']
      return(data.frame("src" = g1, "target" = g2, 
                        "conf" = as.numeric(conf)))
    } else {return (NA)}
  })
  edge_rows <- edge_rows[!is.na(edge_rows)]
  edge_table <- do.call(rbind, edge_rows)
  colnames(edge_table) <- c("src", "target", "conf")
  
  unique_nodes <- unique(c(edge_table$src, edge_table$target))
  node_table <- as.data.frame(unique_nodes)
  colnames(node_table) <- "nodes"
  node_table$dir <- as.factor(unlist(lapply(unique_nodes, function(n) {
    dir <- NA
    if(n == goi) {dir <- 0}
    else {
      dir <- ifelse(as.numeric(dyscover_out[dyscover_out$T_k.name == n,
                                            'estimate'])[1] > 0, 1, -1)
    }
    return(dir)
  })))
  node_table$node_size <- unlist(lapply(node_table$dir, function(d) {
    if(d == 0) {return(2)}
    else {return(1)}
  }))
  print(node_table)
  
  # Convert to ggraph format
  network_table <- tbl_graph(nodes = node_table, edges = edge_table, directed = F)
  
  # Remove small disconnected components
  network_table_sub <- induced_subgraph(
    network_table, V(network_table)[ave(1:vcount(network_table), 
                                        membership(components(network_table)), 
                                        FUN = length) > disconn_comp_size])
  print(igraph::as_data_frame(network_table_sub, what = "vertices"))
  
  # Create ggraph
  driver_targ_network <- ggraph(network_table_sub, layout = "fr") +   # other options: stress
    geom_edge_link0(edge_colour = "darkgray", edge_width = 1) +  
    geom_node_point(aes(fill = dir, size = node_size), shape = 21) +    
    scale_size(range = c(5,10), guide = 'none') + 
    geom_node_label(aes(label = nodes), repel = TRUE, show.legend = F, 
                    label.size = NA, label.padding = 0.05) +
    theme_graph() +
    scale_fill_manual(values = c("#0072B5FF", "lightgreen", "#E18727FF"), 
                      labels = c("Pred. Downregulation", "Driver", 
                                 "Pred. Upregulation"),
                      name = "Direction of Regulation") +
    theme(legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12), legend.position="bottom") 
  
  # Render the network
  show(driver_targ_network)
  
  return(node_table)
}

pik3ca_nodes <- create_graphical_representation_dir(pik3ca_lm_res_fishers_recombB, 
                                                    pc_allGenes[pc_allGenes$R_i.name == "PIK3CA",],
                                                    "PIK3CA", 0.2, string_nw_full, 
                                                    "STRING", 0.4, NA)
kras_nodes <- create_graphical_representation_dir(kras_lm_res_fishers_recombB, 
                                                  pc_allGenes[pc_allGenes$R_i.name == "KRAS",],
                                                  "KRAS", 0.2, string_nw_full, 
                                                  "STRING", 0.2, NA)


# Test for differential response to driver inhibitors when SL targets are KO'd vs. when
# Dyscovr targets are KO'd vs. non-Dyscovr targets are KO'd
# AUC(CTD^2) -- better quality and more cell lines
drug_sensitivity <- read.csv(paste0(PATH, "DepMap/Drug_sensitivity_AUC_(CTD^2).csv"),
                             header = T, check.names = F)
colnames(drug_sensitivity)[1] <- "depmap_id"

# Annotate file with GOI dependency for each cell line
drug_sensitivity_crispr <- merge(crispr, drug_sensitivity, by = "depmap_id")
#drug_sensitivity <- drug_sensitivity[,c(1,2,10:ncol(drug_sensitivity))] #TODO: potentially use lineage info

# Import driver-specific inhibitors
kras_inhibitors <- unlist(read.table(paste0(PATH, "DepMap/kras_inhibitors_noSOS1interactionblockers.txt"), 
                                     header = F)[,1])
pik3ca_inhibitors <- unlist(read.table(paste0(PATH, "DepMap/pik3ca_inhibitors.txt"), 
                                       header = F)[,1])

driver <- "PIK3CA"
qval_thres <- 0.01
driver_hits <- pc_allGenes[(pc_allGenes$R_i.name == driver) &
                             (pc_allGenes$q.value < qval_thres), 'T_k.name']
not_hits <- setdiff(unique(pc_allGenes$T_k.name), driver_hits)

driver_sl_hits <- kras_synleth[kras_synleth$qvalue < 0.2, 'gene']
driver_hits_not_sl <- setdiff(driver_hits, driver_sl_hits)

# TODO: Limit to only cases where driver is mutated?

# Get drug response for each of the hit categories across all the driver-targeting drugs
colnames_drugs_noctrp <- unlist(lapply(colnames(drug_sensitivity_crispr)[15:ncol(drug_sensitivity_crispr)], function(x)
  unlist(strsplit(x, " (", fixed = T))[1]))
colnames(drug_sensitivity_crispr)[15:ncol(drug_sensitivity_crispr)] <- colnames_drugs_noctrp
drug_sensitivity_crispr_kras_cols <- c(1:3, 5, 8, which(colnames(drug_sensitivity_crispr) %fin% 
                                                          kras_inhibitors))  # none
drug_sensitivity_crispr_pik3ca_cols <- c(1:3, 6, 8, which(colnames(drug_sensitivity_crispr) %fin% 
                                                          pik3ca_inhibitors))  # 4 columns
drug_sensitivity_crispr_pik3ca <- drug_sensitivity_crispr[,drug_sensitivity_crispr_pik3ca_cols]

drug_sensitivity_crispr_pik3ca_mut <- drug_sensitivity_crispr_pik3ca[
  drug_sensitivity_crispr_pik3ca$PIK3CA == 1,]

crispr_driver_sl_hits <- drug_sensitivity_crispr_pik3ca_mut[drug_sensitivity_crispr_pik3ca_mut$gene %fin% 
                                                              driver_sl_hits, 'value']
crispr_driver_hits_not_sl <- drug_sensitivity_crispr_pik3ca_mut[drug_sensitivity_crispr_pik3ca_mut$gene %fin% 
                                                                  driver_hits_not_sl, 'value']
crispr_driver_not_hits <- drug_sensitivity_crispr_pik3ca_mut[drug_sensitivity_crispr_pik3ca_mut$gene %fin% 
                                                                  not_hits, 'value']
input <- list("SL.hits" = crispr_driver_sl_hits, "Dyscovr.Non.SL.hits" = crispr_driver_hits_not_sl,
              "Not.hits" = crispr_driver_not_hits)
input_m <- melt(input)
colnames(input_m) <- c("Viability.Upon.KO", "Category")
ggplot(input_m, aes(Viability.Upon.KO, fill = Category)) + 
  geom_density(alpha = 0.2, linewidth = 1) +
  theme_minimal() + ylab("Density") + xlab("Cell Line Viability Upon CRISPR KO") +
  theme(axis.text.x = element_text(
    face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 12),
    axis.ticks = element_blank(), axis.line = element_blank(),
    legend.text = element_text(size = 12), 
    legend.title = element_text(face = "bold", size = 12)) +
  scale_fill_manual(values = c("#0072B5FF", "#BC3C29FF", "#20854EFF"))
  
# Is there a correlation between sensitivity to drugs and sensitivity to KO of
# these targets, but not the other targets?
crispr_driver_sl_hits_df <- drug_sensitivity_crispr_pik3ca_mut[drug_sensitivity_crispr_pik3ca_mut$gene %fin% 
                                                              driver_sl_hits, c(3, 6:9)]
crispr_driver_sl_hits_df$Category <- "SL.hit"
crispr_driver_hits_not_sl_df <- drug_sensitivity_crispr_pik3ca_mut[drug_sensitivity_crispr_pik3ca_mut$gene %fin% 
                                                                  driver_hits_not_sl, c(3, 6:9)]
crispr_driver_hits_not_sl_df$Category <- "Dyscovr.Non.SL.hits"
crispr_driver_not_hits_df <- drug_sensitivity_crispr_pik3ca_mut[drug_sensitivity_crispr_pik3ca_mut$gene %fin% 
                                                               not_hits, c(3, 6:9)]
crispr_driver_not_hits_df$Category <- "Non.hit"
input_df <- do.call(rbind, list(crispr_driver_sl_hits_df, crispr_driver_hits_not_sl_df,
                                crispr_driver_not_hits_df))
input_df$avg.sensitivity <- unlist(lapply(1:nrow(input_df), function(i) {
  vals <- input_df[i, 2:5]
  vals <- vals[!is.na(vals)]
  if(length(vals) > 0) {return(mean(vals))}
  else {return(NA)}
}))
input_df$Category <- as.factor(input_df$Category)
input_df$value <- as.numeric(input_df$value)
input_df$avg.sensitivity <- as.factor(input_df$avg.sensitivity)

# Plot Scatterplot
ggplot(input_df, aes(x = avg.sensitivity, y = value, color = Category)) + 
  geom_line() +
  #geom_bar(position = "stack", stat = "identity") +
  #geom_point() + geom_labelsmooth(aes(label = Category), fill="white", method="lm",
   #                               size = 3, linewidth = 1, boxlinewidth = 0.4) +  # se=F, fullrange=T, 
  scale_color_nejm() + theme_minimal() + 
  xlab("PIK3CA-Inhibitor Avg. Sensitivity") + ylab("Cell Viability Upon CRISPR KO") +
  labs(col = "Category of Target") + 
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12))


