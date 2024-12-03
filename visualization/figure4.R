############################################################
# Code to Create Figure 4 Visualizations
# Written by Sara Geraghty
# https://www.biorxiv.org/content/10.1101/2024.11.20.624509v1
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
#' @param kegg_or_go either 'kegg' or 'go' to indicate which database we are using
#' @param minggssize the minimium gene set size for gseGO (default is 10)
perform_gsea_synleth <- function(results_table, kegg_or_go, mingssize = 10) {
  
  # Get all of this driver protein's targets and their associated mutation 
  # coefficients, sorted descending order
  if(kegg_or_go == "kegg") {

    mapping_kegg <- as.data.frame(bitr_kegg(mapping$UNIPROT, fromType = "uniprot", 
                                            toType = "kegg", drop = T, 
                                            organism = "hsa"))
    mapping <- as.data.frame(bitr(results_table$gene, fromType = "SYMBOL", 
                                  toType = "UNIPROT", OrgDb=org.Hs.eg.db, drop = T))
    mapping <- merge(mapping, mapping_kegg, by = "uniprot")
    
    colnames(mapping_kegg) <- c("uniprot", "kegg")
    print(head(mapping))
    
    results_table <- merge(results_table, mapping, all=T, by="gene")
    results_table <- distinct(results_table[, c("gene", "kegg", "pval")])
    results_table$negLog10pval <- unlist(lapply(results_table$pval, function(p) -log10(p)))
    
    results_table <- results_table[order(results_table$negLog10pval, decreasing = T),]
    #results_table <- results_table[order(results_table$stat, decreasing = T),]
    results_table <- na.omit(results_table)
    
    #betas <- results_table$stat
    betas <- results_table$negLog10pval
    names(betas) <- results_table$kegg
    
    gse.kegg <- gseKEGG(betas, pvalueCutoff = 1, pAdjustMethod = "BH", 
                        keyType = "kegg")
    gse.kegg <- enrichKEGG(gene = results_table[results_table$qval < 0.2, 'uniprot'], 
                           organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.2,
                           keyType = "uniprot")
    #gse.kegg <- gseMKEGG(betas, pvalueCutoff = 1, 
    #pAdjustMethod = "BH", keyType = "kegg")
    return(gse.kegg)
    
  } else if (kegg_or_go == "go") {

    mapping <- as.data.frame(bitr(results_table$gene, fromType = "SYMBOL", 
                                  toType = "ENTREZID", OrgDb=org.Hs.eg.db, drop = T))
    colnames(mapping) <- c("gene", "entrez")
    print(head(mapping))
    
    results_table <- merge(results_table, mapping, all=T, by="gene")
    results_table <- distinct(results_table[, c("gene", "entrez", "pval", "qvalue")])
    results_table$negLog10pval <- unlist(lapply(results_table$pval, function(p) -log10(p)))
    
    results_table <- results_table[order(results_table$negLog10pval, decreasing = T),]
    #results_table <- results_table[order(results_table$stat, decreasing = T),]
    results_table <- na.omit(results_table)
    
    #betas <- results_table$stat
    betas <- results_table$negLog10pval
    names(betas) <- results_table$entrez
    
    gse.go <- gseGO(betas, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, 
                    pAdjustMethod = "BH", keyType = "ENTREZID", ont = "BP", # BP = biological process
                    minGSSize = mingssize)
    gse.go <- setReadable(gse.go, org.Hs.eg.db)
    
    return(gse.go)
    
  } else {
    print("Only implemented for 'kegg' or 'go' database inputs. Please try again.")
    return(NA)
  }
}

gsea_res_tp53_panq0.2_perq0.2_fisherp_subB <- perform_gsea_synleth(lm_res_tp53_panq0.2_perq0.2_fishersp_subB, "go")
gsea_res_tp53_panq0.2_perq0.2_fisherp_subB <- gsea_res_tp53_panq0.2_perq0.2_fisherp_subB@result
gsea_res_pik3ca_panq0.2_perq0.2_fisherp_subB <- perform_gsea_synleth(lm_res_pik3ca_panq0.2_perq0.2_fishersp_subB, "go")
gsea_res_pik3ca_panq0.2_perq0.2_fisherp_subB <- gsea_res_pik3ca_panq0.2_perq0.2_fisherp_subB@result
gsea_res_kras_panq0.2_perq0.2_fisherp_subB <- perform_gsea_synleth(lm_res_kras_panq0.2_perq0.2_fishersp_subB, "go")
gsea_res_kras_panq0.2_perq0.2_fisherp_subB <- gsea_res_kras_panq0.2_perq0.2_fisherp_subB@result

list_of_results_gsea <- list("TP53" = gsea_res_tp53_panq0.2_perq0.2_fisherp_subB, 
                             "PIK3CA" = gsea_res_pik3ca_panq0.2_perq0.2_fisherp_subB, 
                             "KRAS" = gsea_res_kras_panq0.2_perq0.2_fisherp_subB)

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
      buckets <- lapply(unique(results_gsea$qvalue), function(q) {
        sub <- results_gsea[results_gsea$qvalue == q,]
        sub$leading_edge_signal <- unlist(lapply(sub$leading_edge, function(x) {
          spl <- unlist(strsplit(x, ", ", fixed = T))[3]
          spl <- unlist(strsplit(unlist(strsplit(x, "=", fixed = T))[2], 
                                 "%", fixed = T))[1]
          return(as.numeric(spl))
        }))
        return(sub[order(-sub$leading_edge_signal),])
      })
      results_gsea <- do.call(rbind, buckets)
      print(head(results_gsea))
    }
    
    results_gsea_sub <- results_gsea[1:min(n, nrow(results_gsea)), ]
    
    # Convert the p-values to -log10(pvalue)
    negLog10_pvalues <- -log10(results_gsea_sub$pvalue + 1E10^(-10))
    print(head(negLog10_pvalues))
    
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
  
  # Split up pathways longer than N characters with a newline
  input_df$Enriched.Pathway.Driver <- unlist(lapply(input_df$Enriched.Pathway.Driver, function(pw) {
    nchar_pw <- nchar(pw)
    if(nchar_pw > 50) {
      pw_spl <- unlist(strsplit(pw, " ", fixed = T))
      half <- ceiling(length(pw_spl) / 2)
      pw_spl[half] <- paste0(pw_spl[half], "\n")
      return(paste(pw_spl, collapse = " "))
    }
    else{return(pw)}
  }))
  
  #p <- ggplot(input_df, aes(x = reorder(Enriched.Pathway.Driver, abs(ES)), # Enriched.Pathway.Driver,
  p <- ggplot(input_df, aes(x = reorder(Enriched.Pathway.Driver, Neg.Log10.Pval),
                            #y = ES, fill = Driver)) +
                            y = Neg.Log10.Pval, fill = Driver)) +
    geom_col(width = 0.7, color = "black") + coord_flip() + theme_minimal() + 
    theme(legend.position = "none", 
          axis.title = element_text(face = "bold", size = 14), 
          axis.text = element_text(face = "bold", size = 10), 
          strip.text.x = element_text(face="bold", size=14, 
                                      margin = margin(.2, 0, .2, 0, "cm"))) + 
    xlab(paste("Enriched", paste(db, "Pathways"))) + 
    #ylab("Enrichment Score") + 
    ylab("-log10(pval)") + 
    #geom_hline(yintercept=-log10(qval_thres), linetype="dashed", color = "black") +
    scale_fill_manual(values = c("#20854EFF","#BC3C29FF","#0072B5FF")) + 
    scale_x_discrete(labels=roles) +
    facet_wrap(~factor(Driver, levels = c("TP53", "PIK3CA", "KRAS")), 
               ncol = 1, scales = "free_y") 
  
  p
}

create_gsea_barchart_multGenes(list_of_results_gsea, 3, 0.2, "q.value", "GO")


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
