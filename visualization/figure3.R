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
pc_allGenes <- read.csv(paste0(PATH, outfn))

# Per-Cancer
perCancer_fns <- intersect(list.files(path = PATH, recursive = T,
                                      pattern = "_corrected_MUT"), 
                           intersect(list.files(path = PATH, recursive = T,
                                                pattern = "allGenes_"),
                                     list.files(path = PATH, recursive = T,
                                                pattern = "Nonsyn.Drivers.Vogel.elim.vif.5")))
perCancer <- lapply(perCancer_fns, function(f) 
  fread(paste0(PATH, f), header = T))
names(perCancer) <- unlist(lapply(pc_perCancer, function(x)
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
### PART B: PER-CANCER PYRIMIDINE METABOLIC HEATMAP, KRAS 
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
pyrimidine_pw_genes_kegg <- c("NT5E", "CANT1", "TK2", "ENTPD4", "CTPS2", "NME4", 
                              "NME2", "NME1", "DCTD", "ENPP3", "NT5C", "NME6", 
                              "NUDT2")
## MKEGG ##
pyrimidine_pw_genes_mkegg <- c("CTPS2", "NME4", "NME2", "NME6", "NME1", "NME7", 
                               "RRM1", "RRM2B", "TYMS", "RRM2", "DUT", "DTYMK",
                               "CMPK2", "CTPS1")

unique_pyrimidine_set <- unique(c(pyrimidine_pw_genes, pyrimidine_pw_genes_kegg, 
                                  pyrimidine_pw_genes_mkegg))

kegg_pyrimidine_metabolism <- read.csv(paste0(PATH, "KEGG_PYRIMIDINE_METABOLISM.v2023.1.Hs.csv"))
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


############################################################
### PART C: DEPMAP PER-GENE DIFFERENTIAL DEPENDENCY BY DRIVER MUTATION STATUS
############################################################
# DepMap files can be downloaded directly from their website using custom 
# (https://depmap.org/portal/download/custom/) or full (https://depmap.org/portal/download/all/)
# download webpages

# Files from this analysis include:
# 1. CRISPR (DepMap Public 23Q4+Score, Chronos)
# 2. Expression Public 23Q4
# 3. Omics_Somatic_Mutations
# 4. Cell line metadata

# Import CRISPRi, expression, and mutation data
crispr <- read.csv(paste0(PATH, "DepMap/CRISPR_(DepMap_Public_23Q2+Score,_Chronos).csv"), 
                                    header = T, check.names = F)
expression <- read.csv(paste0(PATH, "DepMap/Expression_Public_23Q2.csv"), 
                                        header = T, check.names = F)
mutations <- read.csv(paste0(PATH, "DepMap/OmicsSomaticMutations.csv"),
                      header = T, check.names = F)

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
#' @param depmap_df a DepMap knockout or expression data frame
#' @param mutation_df a DepMap mutation data frame
#' @param genes_of_interest a vector of the Hugo IDs of genes whose mutation
#' status is of interest
adjust_depmap_df <- function(depmap_df, mutation_df, genes_of_interest) {
  
  # Melt data frame to get it ready for plots
  depmap_df_melt <- melt(depmap_df)
  
  # Adjust column names
  colnames(depmap_df_melt)[which(colnames(depmap_df_melt) == "variable")] <- "gene"
  colnames(depmap_df_melt)[1] <- "depmap_id"
  
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

  depmap_df_melt <- merge(depmap_df_melt, mutation_df_new, by = "depmap_id")
  #depmap_df_melt[is.na(depmap_df_melt)] <- 0
  
  return(depmap_df_melt)
}

mutations_drivers <- do.call(rbind, list(mutations_fbxw7, mutations_ctnnb1, 
                                         mutations_idh1, mutations_kras, 
                                         mutations_pik3ca, mutations_tp53))
crispr <- adjust_depmap_df(crispr, mutations_drivers,  
                           c("TP53", "PIK3CA", "KRAS", "IDH1", "CTNNB1", "FBXW7"))
expression <- adjust_depmap_df(expression, mutations_drivers, 
                               c("TP53", "PIK3CA", "KRAS", "IDH1", "CTNNB1", "FBXW7"))


#' Plot boxplot of gene(s) dependency across cell lines, by driver mutation status
#' @param depmap_df a DepMap knockout or expression data frame
#' @param gene_of_interest a Hugo ID for a given driver gene of interest
#' @param targ the name of the target gene(s) of interest
#' @param yaxis_lab y-axis label; either expression, CRISPRi, or RNAi
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
      ggtitle(paste0("\u0394", "NT5E")) + 
      scale_fill_manual(values=c("#6F99ADFF","#E18727FF")) + 
      scale_x_discrete(labels = c(paste0(gene_of_interest, "\nMutated"), 
                                  paste0(gene_of_interest, "\nUnmutated"))) +
      geom_signif(comparisons = list(c("KRAS.WT","KRAS.Mut")), na.rm = T, 
                  test = "wilcox.test", map_signif_level = T) 
    #stat_compare_means(method = "wilcox.test") #comparisons = list(c(0,1)) #, map_signif_level =c("***"=0.001, "**"=0.01, "*"=0.05)
    

  print(bxp_pub)
}

# Call function for various drivers, such as KRAS with one or more target
# genes, such as NT5E 
make_dependency_boxplot(crispr, "KRAS", c("NT5E"), "Cell Viability")

# Repeat just in a particular cell line, like colon (shown)
make_dependency_boxplot(crispr[crispr$depmap_id %fin% cell_line_cancer_type_mapping[["COAD"]],], 
                        "KRAS", c("NT5E"), "Cell Viability")


############################################################
### PART D: DEPMAP PER-GENE EXPRESSION OF DRIVER V. VIABILITY 
### OF CELLS UPON TARGET KNOCKOUT
############################################################
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

#' Run synthetic lethality analysis for a driven gene across all of its cancer 
#' types according to a given q-value threshold and multiple ways of trimming
#' the number of tests (e.g. limiting to only targets with drugs from DrugBank,
#' and limiting to those that have significantly different dependency by driver
#' mutation status)
#' @param driver_name the gene name of a given driver gene of interest
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
#' @param method either "Spearman", "Pearson", "Regression", or "Wilcoxon" to 
#' denote the method we are using to determine significant synthetic lethality
#' @param ts_or_oncogene either "ts" or "oncogene" to denote if the driver is a 
#' tumor suppressor or oncogene, for Wilcoxon analysis
#' @param thres_cls an integer denoting the minimum number of cell lines in a 
#' given cancer type that do not have a mutation in the given driver
get_synthetic_lethals <- function(driver_name, perCancer, qval_thres, use_pc_hits, 
                                  pc_hits, use_drugbank, drugbank_gns, 
                                  use_dependency_check, crispr, expression,
                                  cell_line_cancer_type_mapping, method, 
                                  ts_or_oncogene, thres_cls = 15) {
  
  per_cancer_results <- lapply(1:length(perCancer), function(i) {
    ct <- names(perCancer)[i]
    master_df <- perCancer[[i]]
    
    tophits_ct <- as.character(unlist(master_df[(master_df$q.value < qval_thres) & 
                                                  (master_df$estimate > 0) & 
                                                  (master_df$R_i.name == driver_name), 
                                                'T_k.name']))
    #print(head(tophits_ct))
    
    # Limit CRISPR and expression DepMap data to just this cancer type
    crispr_sub <- crispr[crispr$depmap_id %fin% 
                           c(cell_line_cancer_type_mapping[[ct]]),]
    expression_sub <- expression[expression$depmap_id %fin% 
                                   c(cell_line_cancer_type_mapping[[ct]]),]
    
    # Ensure there are at least a certain number (e.g. 10) lines that do not have
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
        res <- wilcox.test(dep_mut, dep_noMut, alternative = "greater", exact = F)
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
        
        tophits_depmap <- corr_res_sig[(corr_res_sig$pval < 0.05) & 
                                         (corr_res_sig$stat > 0),]
      } else if (method == "Wilcoxon") {
        wilcox_res_sig <- perform_wilcoxon_analysis(crispr_sub, expression_sub,
                                            tophits_ct, driver_name, 0, 
                                            ts_or_oncogene)
        print(wilcox_res_sig)
        wilcox_res_sig <- na.omit(wilcox_res_sig)
        if(length(wilcox_res_sig) == 0) {return(NA)}
        if(nrow(wilcox_res_sig) == 0) {return(NA)}
        wilcox_res_sig$padj <- p.adjust(wilcox_res_sig$pval, "BH")
        wilcox_res_sig <- wilcox_res_sig[order(wilcox_res_sig$padj),]
        
        #tophits_depmap <- wilcox_res_sig
        tophits_depmap <- wilcox_res_sig[(wilcox_res_sig$pval < 0.05),]
      } else if (method == "Regression") {
        lm_res_sig <- perform_lm_analysis(crispr_sub, expression_sub,
                                          tophits_ct, driver_name, 0)
        
        lm_res_sig <- na.omit(lm_res_sig)
        if(length(lm_res_sig) == 0) {return(NA)}
        if(nrow(lm_res_sig) == 0) {return(NA)}
        lm_res_sig$padj <- p.adjust(lm_res_sig$pval, "BH")
        lm_res_sig <- lm_res_sig[order(lm_res_sig$padj),]
        
        tophits_depmap <- lm_res_sig[(lm_res_sig$pval < 0.05) & 
                                       (lm_res_sig$stat > 0),]
      } else {
        print("Only implemented methods are Spearman, Pearson, Wilcoxon, and 
              Regression. Please try again with one of these methods specified.")
        return(NA)
      }
      
      if(length(tophits_depmap) == 0) {return(NA)}
      if(nrow(tophits_depmap) == 0) {return(NA)}
      dependency_by_expression_analysis_plot(crispr_sub, expression_sub, 
                                             tophits_depmap$gene, driver_name, 0)
      
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
    expression_driver_reg <- expression_driver[expression_driver$value <= 
                                                 percentile_90,]
    
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
#' @param ts_or_oncogene either "ts" or "oncogene" to indicate if driver is a
#' tumor suppressor or oncogene
perform_lm_analysis <- function(crispr, expression, list_of_genes, driver_mut, 
                                mut_status) {
  
  expression <- expression[expression$gene %fin% c(driver_mut),]
  
  if(!is.na(driver_mut)) {
    expression <- expression[expression[,which(
      colnames(expression) == driver_mut)] == mut_status,]
  }
  expression_driver <- expression[expression$gene == driver_mut, 
                                  c("depmap_id", "gene", "value")]
  
  lm_dfs <- lapply(list_of_genes, function(gene) {
    print(gene)
    if(!(gene %fin% crispr$gene)) {return(NA)}
    crispr <- crispr[crispr$gene %fin% c(gene),]
    crispr_gene <- crispr[crispr$gene == gene, c("depmap_id", "gene", "value")]
    
    data_input <- merge(crispr_gene, expression_driver, by = "depmap_id")
    print(head(data_input))
    
    if(nrow(data_input) < 5) {return(NA)}
    
    lm_res <- tidy(lm(value.x ~ value.y, data = data_input))
    print(lm_res)
    lm_res <- lm_res[lm_res$term == 'value.y',]
    lm_res_stat <- as.numeric(lm_res$estimate[1])
    lm_res_pval <- lm_res$p.value[1]
    
    print(lm_res_pval)
    return(data.frame("stat" = lm_res_stat, "pval" = lm_res_pval, 
                      "gene" = gene))
    
  })
  
  lm_dfs <- lm_dfs[!is.na(lm_dfs)]
  lm_df <- do.call(rbind, lm_dfs)
  
  return(lm_df)
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


# Define a given driver and subset Dyscovr output to those cancer types with 
# that driver
driver_name <- "KRAS"
perCancer_driver <- lapply(perCancer, function(x) {
  if(driver_name %fin% x$R_i.name) return(x)
  else {return(NA)}
})
perCancer_driver <- perCancer_driver[!is.na(perCancer_driver)]

# Get pan-cancer hits at a given threshold
pc_qval_thres <- 0.01
pc_hits_driver <- pc_allGenes[(pc_allGenes$R_i.name == driver_name) &
                                (pc_allGenes$q.value < pc_qval_thres), 
                              'T_k.name']

# Call function
synthetic_lethals <- get_synthetic_lethals(driver_name, perCancer_driver, 
                                           0.2, T, pc_hits_driver, T, 
                                           targets_with_drugs_gns, F, 
                                           crispr, expression,
                                           cell_line_cancer_type_mapping, 
                                           "Spearman",
                                           "oncogene")

