############################################################
# Code to Create Figure 4 Visualizations
# Written by Sara Geraghty
# PUBLICATION INFORMATION
############################################################

library(data.table)
library(ggplot2)
library("RColorBrewer")
library(ReactomePA)
library("clusterProfiler")
library("DOSE")
library(broom)
library(dplyr)
library(qvalue)

# Local PATH to directory containing Dyscovr output files
PATH <-  paste0(getwd(), "Output/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
### IMPORT RESULTS FROM SYNTHETIC LETHALITY ANALYSIS
############################################################
lm_res_tp53_panq0.2_perq0.2_fishersp_subB <- fread(paste0(PATH, "Synthetic_Lethality/TP53/synthetic_lethals_pancancer_regression_qn_panq0.2_perq0.2_browns_dirsubB.csv"),
                                                   header = T)
lm_res_pik3ca_panq0.2_perq0.2_fishersp_subB <- fread(paste0(PATH, "Synthetic_Lethality/PIK3CA/synthetic_lethals_pancancer_regression_qn_panq0.2_perq0.2_browns_dirsubB.csv"),
                                                   header = T)
lm_res_kras_panq0.2_perq0.2_fishersp_subB <- fread(paste0(PATH, "Synthetic_Lethality/KRAS/synthetic_lethals_pancancer_regression_qn_panq0.2_perq0.2_browns_dirsubB.csv"),
                                                     header = T)

############################################################
### PART B: GSEA ON SYNTHETIC LETHALS
############################################################
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


############################################################
### CHECK FOR ENRICHMENT IN GOLD STANDARD SL SETS (IN-TEXT)
############################################################
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

# Try this for pairings from SynLethDB as well (Human): https://synlethdb.sist.shanghaitech.edu.cn/#/download
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

# Use 'computate_statistical_enrichment' function from figure2.R
compute_statistical_enrichment(synthetic_lethals, intersect(gold_standard_sls_synlethdb_driver, 
                                                            gene_universe), 
                               "both", 0.2, NA)
