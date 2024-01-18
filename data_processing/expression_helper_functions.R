############################################################
### Expression (RNA-Seq) Helper Functions
### ASSOCIATED PUBLICATION INFORMATION
############################################################

library(stringr)
library(edgeR)
library(dplyr)
library(tibble)
library(parallel)
library(gdata)
library(matrixStats)
library(data.table)
library("gplots")

# Includes functions for additional visualization and processing of TCGA
# RNA-Seq data, not included in publication

# Local path to directory containing gene expression data files
PATH <- getwd()
OUTPUT_PATH <- paste0(getwd(), "Output/Expression/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)



############################################################
### ADJUST EDGER GROUPS TO GROUP BY MUTATION STATUS
### OF A GIVEN PROTEIN OF INTEREST 
############################################################
#' Takes a subsetted edgeR targets DF and changes the grouping
#' to be by the mutation status of a given protein of interest
#' (e.g., TP53) so that we can look at DE genes between the 
#' mutated and non-mutated sample groups.
#' @param targets_sub a subsetted edgeR targets DF
#' @param prot a protein whose mutation status we are 
#' interested in (Hugo Symbol)
#' @param mutation_df a mutation count DF that correlates patient
#' ID to the mutation status of the given protein of interest
group_by_mutation_status <- function(targets_sub, prot, mutation_df) {
  
  # Get the patients with/ without a mutation in the given protein
  mutation_df_sub <- mutation_df[mutation_df$Gene_Symbol == prot, 
                                 2:ncol(mutation_df)]
  patients_mut <- colnames(mutation_df_sub)[which(mutation_df_sub[1,] > 0)]
  #print(colnames(mutation_df_sub)[which(mutation_df_sub[1,] > 0)])
  
  # Add these new groups to the edgeR targets DF under the 'group' label 
  # (patients with and without mutation)
  targets_sub_patients <- unlist(lapply(targets_sub$description, function(x)
    paste(unlist(strsplit(x, "-", fixed = T))[3:4], collapse = "-")))
  targets_sub$group <- unlist(lapply(targets_sub_patients, function(x) 
    ifelse((x %fin% patients_mut), 1, 0)))
  
  return(targets_sub)
}

# Define the protein of interest
prot_of_interest <- "TP53"  # TP53

# Import the mutation count DF (see process_mutation_data.R)
mut_count_matrix <- read.csv(paste0(OUTPUT_PATH, 
                                    "Mutation/mut_count_matrix_nonsynonymous.csv"),
                             header = T, check.names = F, row.names = 1)

# Call function
targets_groupedByTP53Mut <- group_by_mutation_status(targets, prot_of_interest, 
                                                     mut_count_matrix)
targets_groupedByTP53Mut <- distinct(targets_groupedByTP53Mut)


############################################################
### ESTIMATE & VISUALIZE NAIVE DIFF. GENE EXPRESSION ###
############################################################
#' Estimate and visualize top differentially expressed genes,
#' without using a linear modeling approach. Examine the 
#' differences between the DE genes by mutation status of a 
#' given protein (e.g. TP53).  
#' @param d edgeR output data structure to visualize
#' @param target_group an optional vector of gene ensg IDs of interest 
#' (e.g. metabolic genes)
#' @param method BH or Qvalue, to denote our MHT correction method
#' @param thres an adjusted p/ qvalue threshold for significance
#' @param all_genes_id_conv an ID conversion table from bioMart
visualize_de_by_mut_status <- function(d, target_group, method, 
                                       thres, all_genes_id_conv) {
  
  # Subset to the target group
  d_rownames <- unlist(lapply(unlist(rownames(d$counts)), function(x) 
    unlist(strsplit(x, ".", fixed = T))[1]))
  d <- d[which(d_rownames %fin% target_group),]
  
  d <- visualize_dispersion(d)
  
  # Get the DE genes across groups (1 is mutated, 0 is not mutated)
  count_norm <- cpm(d)
  count_norm <- as.data.frame(count_norm)
  
  conditions <- d$samples[, c("description", "group")]
  
  # Run the Wilcoxon rank-sum test for each gene
  pvalues <- sapply(1:nrow(count_norm), function(i){
    data <- cbind.data.frame(gene.exp = as.numeric(t(count_norm[i,])), 
                             mut.stat = conditions$group)
    p <- wilcox.test(gene.exp~mut.stat, data)$p.value
    return(p)
  })
  pval_df <- data.frame("gene" = unlist(lapply(rownames(count_norm), function(x)
    unlist(strsplit(x, ".", fixed = T))[1])), "Pvalue" = pvalues)  
  
  if(method == "BH") {
    topTags <- topTags(et, n = 10000, adjust.method = "BH", p.value = thres)
    de1 <- decideTestsDGE(et, adjust.method = "BH", p.value = thres)
    print(summary(de1))
    
    # Adjust the topTags table
    rownames(topTags$table) <- unlist(lapply(rownames(topTags$table), function(x) 
      unlist(strsplit(as.character(x), ".", fixed = T))[1]))
    print(head(topTags$table))
    topTags$table$gene.name <- unlist(lapply(
      rownames(topTags$table), function(ensg)
        paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == 
                                                ensg, 'external_gene_name'])), 
              collapse = ";")))
    
    return(topTags)
  } 
  
  else if(method == "Qvalue") {
    pval_df <- na.omit(pval_df)
    pval_df <- pval_df[!(is.nan(pval_df$Pvalue) | is.infinite(pval_df$Pvalue)),]
    
    qobj <- qvalue(as.numeric(pval_df$Pvalue))
    pval_df$Qvalue <- qobj$qvalues
    
    pval_df_sub <- pval_df[pval_df$Qvalue < thres, ]
    
    # Add gene names
    pval_df_sub$gene.name <- unlist(lapply(pval_df_sub$gene, function(ensg) {
      return(paste(unique(unlist(
        all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == ensg, 
                          'external_gene_name'])), collapse = ";"))
    }))
    
    return(pval_df_sub)
  }
  
  else {
    print("Only implemented for BH and Qvalue")
    return(NA)
  }
}

# Read in a gene set of interest
metabolism_gene_list_ensg <- read.table(paste0(PATH, 
                                               "metabolism_gene_list_ensg.txt"),
                                        header = T)[,1]
# Convert to external gene names, if needed
metabolism_gene_list <- unlist(lapply(metabolism_gene_list_ensg, function(x) 
  paste(unique(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == 
                                   unlist(strsplit(x, ";", fixed = T))[1], 
                                 'external_gene_name']), collapse = ";")))

# Call function for an individual cancer type
topHits_tp53 <- visualize_de_by_mut_status(d, metabolism_gene_list_ensg, 
                                           "Qvalue", 0.2, all_genes_id_conv)

# Or, if we are looking pan-cancer, apply to each item in list
lapply(d_list, visualize_de_by_mut_status)

# Adjust to perform this analysis by subtypes (shown for BRCA)
d$samples$description <- unlist(lapply(d$samples$description, function(x) 
  unlist(strsplit(x, "-", fixed = T))[3]))
d_lumA <- d[,which(unlist(d$samples$description) %fin% patient_set_lumA)]
d_lumB <- d[,which(unlist(d$samples$description) %fin% patient_set_lumB)]
d_basal <- d[,which(unlist(d$samples$description) %fin% patient_set_basal)]
d_her2 <- d[,which(unlist(d$samples$description) %fin% patient_set_her2)]

topHits_tp53_lumA <- visualize_de_by_mut_status(d_lumA, metabolism_gene_list_ensg, 
                                                "Qvalue", 0.2, all_genes_id_conv)
topHits_tp53_lumB <- visualize_de_by_mut_status(d_lumB, metabolism_gene_list_ensg, 
                                                "Qvalue", 0.2, all_genes_id_conv)
topHits_tp53_basal <- visualize_de_by_mut_status(d_basal, metabolism_gene_list_ensg, 
                                                 "Qvalue", 0.2, all_genes_id_conv)
topHits_tp53_her2 <- visualize_de_by_mut_status(d_her2, metabolism_gene_list_ensg, 
                                                "Qvalue", 0.2, all_genes_id_conv)


#' Visualize the top DE genes by mutation status of multiple genes of interest 
#' using a heat map
#' @param topHits_gene1 a formal class TopTags object with the top hits from the 
#' first gene of interest
#' @param topHits_gene2 a formal class TopTags object with the top hits from the 
#' second gene of interest
#' @param gene1_name the external gene name of the first gene of interest
#' @param gene2_name the external gene name of the second gene of interest
#' @param subset_genes an optional subset of genes that we are interested in 
#' looking at expression effects of (e.g. metabolic genes)
create_de_by_mut_status_heatmap <- function(topHits_gene1, topHits_gene2, 
                                            gene1_name, gene2_name, subset_genes) {
  
  # Get the tables from the TopTags object
  topHits_gene1_tab <- topHits_gene1$table
  topHits_gene2_tab <- topHits_gene2$table
  
  # Limit them to just the gene.name and the logFC
  topHits_gene1_tab <- topHits_gene1_tab[,c("logFC", "gene.name")]
  topHits_gene2_tab <- topHits_gene2_tab[,c("logFC", "gene.name")]
  
  # Adjust the logFC column name
  colnames(topHits_gene1_tab)[which(colnames(topHits_gene1_tab) == "logFC")] <- 
    paste(gene1_name, "logFC", sep = ".")
  colnames(topHits_gene2_tab)[which(colnames(topHits_gene2_tab) == "logFC")] <- 
    paste(gene2_name, "logFC", sep = ".")
  
  # Remove any gene.name entries that are blank
  topHits_gene1_tab <- topHits_gene1_tab[-which(topHits_gene1_tab$gene.name == ""),]
  topHits_gene2_tab <- topHits_gene2_tab[-which(topHits_gene2_tab$gene.name == ""),]
  
  # Merge them together by gene name
  topHits_tab <- merge(topHits_gene1_tab, topHits_gene2_tab, 
                       by = "gene.name", all = T)
  
  # Deal with duplicates by adding .X to them (where X is a number)
  duplicates <- unique(topHits_tab$gene.name[which(duplicated(
    topHits_tab$gene.name))])
  
  for (d in duplicates) {
    indices_of_d <- which(topHits_tab$gene.name == d)
    for (i in 1:length(indices_of_d)) {
      index <- indices_of_d[i]
      topHits_tab$gene.name[index] <- paste(d, i, sep = ".")
    }
  }
  
  # Make the gene.name a rowname rather than a column for plotting
  rownames(topHits_tab) <- topHits_tab$gene.name
  topHits_tab <- topHits_tab[, which(colnames(topHits_tab) != "gene.name")]
  
  # Turn into a matrix
  topHits_matrix <- data.matrix(topHits_tab)
  
  # For the purposes of the heatmap, convert all NAs to 0
  topHits_matrix[is.na(topHits_matrix)] <- 0
  
  # If "subset_genes" is given, limit the matrix to just these genes 
  try(topHits_matrix <- topHits_matrix[which(rownames(topHits_matrix) %fin% 
                                               subset_genes),], silent = T)
  
  # Plot an enhanced heatmap
  heatmap.2(topHits_matrix, scale = "none", col = bluered(100), trace = "none",
            density.info = "none", labCol = "", dendrogram = c("row"), 
            add.expr = text(x = seq_along(colnames(topHits_matrix)), 
                            y = -2, srt = 0, labels = colnames(topHits_matrix), 
                            xpd = NA, cex = 2, pos = 1))
  
  return(topHits_matrix)
}

# Call function for genes of interest (shown for TP53, PIK3CA)
logFC_matrix <- create_de_by_mut_status_heatmap(topHits_tp53, topHits_pik3ca, 
                                                "TP53", "PIK3CA")
logFC_matrix_metabolism <- create_de_by_mut_status_heatmap(topHits_tp53, 
                                                           topHits_pik3ca, 
                                                           "TP53", "PIK3CA", 
                                                           metabolic_genes_gn)


############################################################
#### ESTIMATE SPEARMAN CORRELATION OF LOGFC ###
############################################################
#' Compute and print the Spearman correlation of the coefficients and/or the 
#' T-statistic, given two groups of interest
#' @param results_table the output matrix from "create_de_by_mut_status_heatmap" 
#' function
#' @param ri_1 the external gene name of the first protein of interest
#' @param ri_2 the external gene name of the second protein of interest
compute_and_print_spearman <- function(results_table, ri_1, ri_2) {
  
  target_genes <- unique(rownames(results_table))
  
  grp1_logFC <- results_table[,1]
  grp2_logFC <- results_table[,2]
  
  # Get logFC spearman
  logFC_spearman <- cor.test(grp1_logFC, grp2_logFC, method = "spearman")
  logFC_spearman_stat <- as.numeric(logFC_spearman$estimate)
  logFC_spearman_pval <- logFC_spearman$p.value
  
  # Print the results
  print(paste("Spearman results for", paste(ri_1, paste("and", ri_2))))
  print(paste("logFC correlation of", 
              paste(logFC_spearman_stat, 
                    paste(", p-value of", logFC_spearman_pval))))
  
  # Create a plot to visualize the correlations
  plot(grp1_logFC, grp2_logFC, pch = 19, col = "lightblue", 
       xlab = paste(ri_1, "logFC"), ylab = paste(ri_2, "logFC"))
  abline(lm(grp2_logFC ~ grp1_logFC), col = "red", lwd = 3)
  text(labels = paste("Correlation:", 
                      paste(round(logFC_spearman_stat, 6), 
                            paste(", p-value:", round(logFC_spearman_pval, 6)))), 
       x = max(grp1_logFC, na.rm = T)-sd(grp1_logFC, na.rm = T)*5, 
       y = max(grp2_logFC, na.rm = T)-sd(grp2_logFC, na.rm = T), col = "black")
}



# Call function for genes of interest (shown for TP53 and PIK3CA)
compute_and_print_spearman(logFC_matrix, "TP53", "PIK3CA")
compute_and_print_spearman(logFC_matrix_metabolism, "TP53", "PIK3CA")


############################################################
#### VISUALIZE OVERLAP BETWEEN NAIVE DE GENES AND THOSE 
#### PREDICTED BY DYSCOVR 
############################################################
# Import results from Dyscovr (sample file name: "results_file_mutation.csv")
dyscovr_df <- read.csv(paste0(OUTPATH, "Dyscovr/results_file_mutation.csv"), 
                       header = T, check.names = F)

# Use processing functions from process_dyscovr_output.R to get 
# top hits from model results
dyscovr_df_sig <- get_signif_correl(dyscovr_df, qval_thres = 0.1)
tophits_tp53 <- dyscovr_df[dyscovr_df$R_i.name == "TP53", 'T_k.name']

# Make Venn diagram
vd_tp53 <- venn.diagram(list(tophits_tp53, 
                             logFC_matrix_metabolism$TP53.logFC[
                               logFC_matrix_metabolism$TP53.logFC > 0]),
                        category.names = c("Model Hits", "Naive DE Hits"), 
                        filename = NULL, output = T, lwd = 2, lty = 'blank', 
                        fill = c("red", "blue"), cex = 2, fontface = "bold",
                        fontfamily = "sans", cat.cex = 2, cat.fontface = "bold",
                        cat.fontfamily = "sans")
grid::grid.draw(vd_tp53)


############################################################
### GET A TUMOR-ONLY VERSION OF AN EXPRESSION DF ###
############################################################
#' In case a tumor-only version is needed (e.g., for PEER),
#' this function takes an expression DF and limits it to just cancer
#' samples
#' @param expression_df
limit_to_tumor_only <- function(expression_df) {
  cols_to_keep <- which(grepl("-0", colnames(expression_df)))
  expression_df <- expression_df[, cols_to_keep]
  return(expression_df)
}

# Call function on an imported expression DF
expression_df_tumor_only <- limit_to_tumor_only(expression_df)

write.csv(expression_df_tumor_only, paste0(OUTPUT_PATH, "expression_DF_TO.csv"))


############################################################
### IMPORT AND PROCESS NORMAL GTEx DATA FOR COMPARISON ###
############################################################
# Import GTEx expression data for a given cancer type or pan-cancer
gtex_reads <- read.table(paste0(PATH, "GTEx/gene_reads.gct"), 
                                header = T, check.names = F)

# Remove the version from the ENSG IDs
gtex_reads$Name <- unlist(lapply(gtex_reads$Name, function(x)
  unlist(strsplit(x, ".", fixed = T))[1]))

# Import results from Dyscovr (sample file name: "results_file_mutation.csv")
dyscovr_df <- read.csv(paste0(OUTPATH, "Dyscovr/results_file_mutation.csv"), 
                       header = T, check.names = F)

#' Process GTEx expression data using the same pipeline as TCGA expression data,
#' for uniform comparison
#' @param gtex_reads a data frame of RNA-seq reads from GTEx
#' @param dyscovr_df a results data frame from Dyscovr for the same tissue type
process_gtex_data <- function(gtex_reads, dyscovr_df) {
  # Limit to only genes found in Dyscovr's output
  gtex_reads_sub <- gtex_reads[gtex_reads$Description %fin% dyscovr_df$T_k.name, ]
  gtex_reads_foredgeR <- gtex_reads_sub
  rownames(gtex_reads_foredgeR) <- gtex_reads_foredgeR$Name
  gtex_reads_foredgeR <- gtex_reads_foredgeR[,4:ncol(gtex_reads_foredgeR)]
  
  d <- DGEList(counts = gtex_reads_foredgeR)
  
  # Use edgeR's filterByExpr() function
  # Function documentation: https://rdrr.io/bioc/edgeR/man/filterByExpr.html
  keep <- filterByExpr(cpm(d)) 
  d <- d[keep,]
  
  # Reset the library sizes
  d$samples$lib.size <- colSums(d$counts)
  
  # Reset the colnames to the sample names
  colnames(d$counts) <- d$samples$description
  
  # Perform TMM normalization
  d <- calcNormFactors(d, method = "TMM")
  tmm <- cpm(d)
  colnames(tmm) <- rownames(d$samples)
  
  counts_filtByExpr <- d$counts
  colnames(counts_filtByExpr) <- rownames(d$samples)
  
  return(list(counts_filtByExpr, tmm))
}

# Call function
res <- process_gtex_data(gtex_reads, dyscover_df)


# Write to file
write.csv(res[[1]], paste0(OUTPUT_PATH, "GTEx/gene_reads_filtByExpr.csv"))
write.csv(res[[2]], paste0(OUTPUT_PATH, "GTEx/gene_reads_filtByExpr_tmm.csv"))

# Read back, as needed
counts_filtByExpr <- read.csv(paste0(OUTPUT_PATH, 
                                     "GTEx/gene_reads_filtByExpr.csv"),
                              header = T, check.names = F, row.names = 1)
tmm <- read.csv(paste0(OUTPUT_PATH, "GTEx/gene_reads_filtByExpr_tmm.csv"),
                header = T, check.names = F, row.names = 1)

# Get the average expression per gene in normal tissue
avgGeneExpr <- data.frame("avg.expr" = as.numeric(rowMeans(tmm)), 
                          "ensg" = rownames(tmm))
avgGeneExpr$gene_name <- unlist(lapply(avgGeneExpr$ensg, function(ensg) 
  gtex_reads[gtex_reads$Name == ensg, 'Description']))

# Get the median expression per gene in normal tissue
medianGeneExpr <- data.frame("median.expr" = as.numeric(rowMedians(tmm)), 
                             "ensg" = rownames(tmm))
medianGeneExpr$gene_name <- unlist(lapply(medianGeneExpr$ensg, function(ensg) 
  gtex_reads[gtex_reads$Name == ensg, 'Description']))


# Multiply these values by the coefficients fit by Dyscovr, plus the original 
# expression value (for how things change when a driver is mutated)
#' @param dyscovr_df a results data frame from Dyscovr for the same tissue type
#' @param goi the gene-of-interest external name (e.g. a given driver)
#' @param geneExpr either the average or median expression per gene in normal
#' tissue
mult_expr_by_coeff <- function(dyscovr_df, goi, geneExpr) {
  dyscovr_df_goi <- dyscovr_df[dyscovr_df$R_i.name == goi,]
  colname <- paste("avg.expr.mult.by.betas", goi, sep = ".")
  geneExpr[,colname] <- unlist(lapply(1:nrow(geneExpr), function(i) {
    name <- geneExpr$gene_name[i]
    expr_val <- geneExpr$avg.expr[i]
    beta <- dyscovr_df_goi[dyscovr_df_goi$T_k.name == name, 'estimate']
    val <- ((beta * expr_val)[1]) + expr_val[1]
    return(val)
  }))
  return(geneExpr)
}

# Call function for driver genes of interest
avgGeneExpr <- mult_expr_by_coeff(master_df, "TP53", avgGeneExpr)
avgGeneExpr <- mult_expr_by_coeff(master_df, "PIK3CA", avgGeneExpr)

# Compare this to the TMM-normalized TCGA cancer expression averages in the TCGA;
# import processed TMM expression DF from process_expression_data.R
expression_df_tmm <- read.csv(paste0(OUTPUT_PATH, 
                                     "tmm_normalized_expression_counts.csv"),
                              header = T, check.names = F, row.names = 1)

# Limit to just the same set of genes
expression_df_tmm_sub <- expression_df_tmm[rownames(expression_df_tmm) %fin% 
                                             avgGeneExpr$ensg,]
colnames(expression_df_tmm_sub) <- unlist(lapply(
  colnames(expression_df_tmm_sub), function(x) 
    paste(unlist(strsplit(x, "-", fixed = T))[3:4], collapse = "-")))

# Import the mutation count DF (see process_mutation_data.R)
mut_count_matrix <- read.csv(paste0(OUTPUT_PATH, 
                                    "Mutation/mut_count_matrix_nonsynonymous.csv"),
                             header = T, check.names = F, row.names = 1)

#' Limit to just samples with a nonsynonymous mutation in the driver gene
#' @param mutation_count_matrix a sample (column) by gene (row) nonsynonymous
#' mutation count matrix from process_mutation_data.R
#' @param goi the name of a gene-of-interest whose mutations we are interested in
get_mutant_patients <- function(mutation_count_matrix, goi) {
  regprot_row <- mutation_count_matrix[mutation_count_matrix$Gene_Symbol == goi, ]
  ind <- which(as.integer(regprot_row) > 0)
  return(colnames(regprot_row)[ind])
}

# Call for a given driver of interest
tp53_mutant_patients <- get_mutant_patients(mut_count_matrix, "TP53")
tp53_nonmutant_patients <- setdiff(colnames(expression_df_tmm), 
                                   tp53_mutant_patients)

expression_df_tmm_p53Mut <- expression_df_tmm[, colnames(expression_df_tmm) %fin% 
                                                tp53_mutant_patients]
expression_df_tmm_p53Nonmut <- expression_df_tmm[, colnames(expression_df_tmm) %fin% 
                                                   tp53_nonmutant_patients]

avgGeneExpr_Cancer <- data.frame("avg.expr.cancer" = as.numeric(rowMeans(
  expression_df_tmm_p53Mut)), "ensg" = rownames(expression_df_tmm_p53Mut))
medianGeneExpr_Cancer <- data.frame("median.expr.cancer" = as.numeric(rowMedians(
  as.matrix(expression_df_tmm_p53Mut))), 
  "ensg" = rownames(expression_df_tmm_p53Mut))

avgGeneExpr_merged <- merge(avgGeneExpr, avgGeneExpr_Cancer, by = "ensg")
medianGeneExpr_merged <- merge(medianGeneExpr, medianGeneExpr_Cancer, 
                               by = "ensg")


#' Calculate spearman correlation between GTEx and TCGA expression of a given
#' gene of interest 
#' @param log_gtex the log2 of the average gene expression for GTEx
#' @param log_tcga the log2 of the average gene expression for TCGA
#' @param goi the external name of a gene of interest
calc_and_plot_spearman <- function(log_gtex, log_tcga, goi) {
  logFC_spearman_avg <- cor.test(log_gtex, log_tcga, method = "spearman")
  logFC_spearman_avg_stat <- as.numeric(logFC_spearman_avg$estimate)
  logFC_spearman_avg_pval <- logFC_spearman_avg$p.value
  
  # Print the results
  print(paste("Correlation of", 
              paste(logFC_spearman_avg_stat, 
                    paste(", p-value of", logFC_spearman_avg_pval))))
  
  # Create a plot to visualize the correlations
  plot(log_gtex, log_tcga, pch = 19, col = "lightblue", 
       xlab = paste("GTEx TMM-Norm Expr in Breast *", paste(goi, "Betas + 1")), 
       ylab = paste("TCGA TMM-Norm Expr in Breast,", paste(goi, "Mutant Samp")))
  abline(lm(log_tcga ~ log_gtex), col = "red", lwd = 3)
  text(labels = paste("Correlation:", 
                      paste(round(logFC_spearman_avg_stat, 6), 
                            paste(", p-value:", 
                                  round(logFC_spearman_avg_pval, 6)))), 
       x = max(log_gtex, na.rm = T) - sd(log_gtex, na.rm = T)*5, 
       y = max(log_tcga, na.rm = T)-sd(log_tcga, na.rm = T), col = "black")
  lm_res <- summary(lm(log_tcga ~ log_gtex))
  r2 <- lm_res$r.squared
  pval <- as.numeric(lm_res$coefficients[,4])[2]
  print(paste("R2:", paste(r2, paste(", p-value:", pval))))
}

# Calculate Spearman correlation and visualize
log_gtex <- log2(avgGeneExpr_merged$avg.expr.mult.by.tp53.betas+1)
log_tcga <- log2(avgGeneExpr_merged$avg.expr.cancer+1)
gtex <- avgGeneExpr_merged$avg.expr.mult.by.tp53.betas
tcga <- avgGeneExpr_merged$avg.expr.cancer

# Call function
calc_and_plot_spearman(log_gtex, log_tcga)
calc_and_plot_spearman(gtex, tcga)

# Limit to just significant hits (at some q-value threshold)
avgGeneExpr_merged_sig <- avgGeneExpr_merged[avgGeneExpr_merged$gene_name %fin% 
                                               dyscovr_df[dyscovr_df$q.value < 
                                                            0.2, 'T_k.name'],]

# Calculate Spearman correlation (subsetted to significant hits) and visualize
log_gtex_sig <- log2(avgGeneExpr_merged_sig$avg.expr.mult.by.tp53.betas+1)
log_tcga_sig <- log2(avgGeneExpr_merged_sig$avg.expr.cancer+1)
gtex_sig <- avgGeneExpr_merged_sig$avg.expr.mult.by.betas.TP53
tcga_sig <- avgGeneExpr_merged_sig$avg.expr.cancer

# Call function
calc_and_plot_spearman(log_gtex_sig, log_tcga_sig)
calc_and_plot_spearman(gtex_sig, tcga_sig)


# Try with patients without a driver mutation (ideally should 
# see no correlation)
avgGeneExpr_cancer_WT <- data.frame("avg.expr.cancer" = as.numeric(rowMeans(
  expression_df_tmm_p53Nonmut)), "ensg" = rownames(expression_df_tmm_p53Nonmut))
avgGeneExpr_merged_WT <- merge(avgGeneExpr, avgGeneExpr_cancer_WT, by = "ensg")
avgGeneExpr_merged_WT_sig <- avgGeneExpr_merged_WT[
  avgGeneExpr_merged_WT$gene_name %fin% dyscovr_df[dyscovr_df$q.value < 0.2, 
                                                   'T_k.name'],]

# Calculate Spearman correlation and visualize
log_gtex_p53WT <- log2(avgGeneExpr_merged_p53WT$avg.expr.mult.by.tp53.betas+1)
log_tcga_p53WT <- log2(avgGeneExpr_merged_p53WT$avg.expr.cancer+1)
gtex_WT <- avgGeneExpr_merged_WT$avg.expr.mult.by.betas.TP53
tcga_WT <- avgGeneExpr_merged_WT$avg.expr.cancer

calc_and_plot_spearman(log_gtex_p53WT, log_tcga_p53WT)
calc_and_plot_spearman(gtex_WT, tcga_WT)

log_gtex_p53WT_sig <- log2(avgGeneExpr_merged_p53WT_sig$avg.expr.mult.by.tp53.betas+1)
log_tcga_p53WT_sig <- log2(avgGeneExpr_merged_p53WT_sig$avg.expr.cancer+1)
gtex_WT_sig <- avgGeneExpr_merged_WT_sig$avg.expr.mult.by.tp53.betas
tcga_WT_sig <- avgGeneExpr_merged_WT_sig$avg.expr.cancer

calc_and_plot_spearman(log_gtex_p53WT_sig, log_tcga_p53WT_sig)
calc_and_plot_spearman(gtex_WT_sig, tcga_WT_sig)


############################################################
### VISUALIZE EXPRESSION OF TARGET GENE
############################################################
#' Visualize the expression of a particular target gene using a scatterplot
#' across all patients
#' @param expression_df a data frame with expression values (columns are patient 
#' IDs-sample IDs, rows are ENSG IDs)
#' @param ensg the ENSG ID of the target gene of interest
visualize_tg_expression <- function(expression_df, ensg) {
  plot(as.numeric(expression_df[expression_df$ensg_id == ensg,]), 
       xlab = "Patient", ylab = "Expression")
  abline(h = mean(as.numeric(expression_df[expression_df$ensg_id == ensg,])), 
         col = "red")
}

############################################################
### VISUALIZE EXPRESSION OF TARGET GENE, BY CANCER V. NORMAL
############################################################
#' Get expression profiles for patient cancer samples vs. normal samples
#' and visualize using a box plot
#' @param expression_df a data frame with expression values (columns are patient 
#' IDs-sample IDs, rows are ENSG IDs)
#' @param ensg the ENSG ID of the target gene of interest
visualize_tg_expression_cvn <- function(expression_df, ensg) {
  expression_cancer <- as.numeric(
    expression_df[rownames(expression_df) == ensg, 
                  grepl("-0", colnames(expression_df))])
  expression_normal <- as.numeric(
    expression_df[rownames(expression_df) == ensg, 
                  grepl("-11", colnames(expression_df))])
  plot(expression_cancer, xlab = "Patient", ylab = "Expression", col = "red")
  points(expression_normal, col = "blue")
  legend(500, 3000, legend = c("Cancer", "Normal"), col = c("red", "blue"), 
         lty = 1)  # Adjust placement of legend 
  
  # Make box plot versions
  dataList <- list("Cancer" = expression_cancer, "Normal" = expression_normal)
  boxplot(dataList, ylab = "Expression")
  
  #simple t-test
  ttest_res <- t.test(expression_normal, expression_cancer, alternative = "less")
  print(ttest_res$p.value)
}


############################################################
### VISUALIZE EXPRESSION OF TARGET GENE, BY MUTATION/ CNA/ 
### METHYLATION STATUS OF REGULATORY GENE OF INTEREST
############################################################
#' Limit to just samples with a nonsynonymous mutation in the driver gene
#' @param mutation_count_matrix a sample (column) by gene (row) nonsynonymous
#' mutation count matrix from process_mutation_data.R
#' @param goi the name of a gene-of-interest whose mutations we are 
#' interested in
get_mutant_patients <- function(mutation_count_matrix, goi) {
  regprot_row <- mutation_count_matrix[mutation_count_matrix$Gene_Symbol == goi, ]
  ind <- which(as.integer(regprot_row) > 0)
  return(colnames(regprot_row)[ind])
}

#' Get a vector of patients that have a deletion in this given 
#' gene-of-interest
#' @param cna_df a data frame with information about CNAs in given genes
#' @param goi_ensg the name of a gene-of-interest whose CNAs we are 
#' interested in
get_deleted_patients <- function(cna_df, goi_ensg) {
  cnas <- unlist(cna_df[cna_df$ensg_id == goi_ensg, ])
  to_keep <- which(cnas < 2)
  samps <- colnames(cna_df)[to_keep]
  return(samps)
}

#' Get a vector of patients that have an amplification in this given 
#' gene-of-interest
#' @param cna_df a data frame with information about CNAs in given genes
#' @param goi_ensg the name of a gene-of-interest whose CNAs we are 
#' interested in
get_amplified_patients <- function(cna_df, goi_ensg) {
  cnas <- unlist(cna_df[cna_df$ensg_id == goi_ensg, ])
  to_keep <- which(cnas > 2)
  samps <- colnames(cna_df)[to_keep]
  return(samps)
}

#' Get a vector of patients that have a methylated version of this given 
#' gene-of-interest
#' @param methylation_df a data frame with information about methylation 
#' in given genes
#' @param goi the name of a gene-of-interest whose methylation we are 
#' interested in
get_methylated_patients <- function(methylation_df, goi) {
  methylation_df_goi <- methylation_df[methylation_df$Gene_Symbol == goi, 
                                           3:ncol(methylation_df), with = F]
  samps <- unlist(lapply(3:ncol(methylation_df_goi), function(i) {
    x <- as.numeric(unlist(methylation_df_goi[, i, with = F]))
    # 0 as threshold -- see paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
    if(x > 0) {return(colnames(methylation_df)[i])}  
    else {return(NA)}
  }))
  samps <- samps[!is.na(samps)]
  return(samps)
}


# IMPORT NEEDED MATRICES

# Call for a given driver of interest
tp53_mutant_patients <- get_mutant_patients(mut_count_matrix, "TP53")
tp53_deleted_patients <- get_deleted_patients(cna_df, tp53_ensg)
pik3ca_amplified_patients <- get_amplified_patients(cna_df, pik3ca_ensg) 
kras_methylated_patients <- get_methylated_patients(methylation_df_M, "KRAS")


#' Get the expression in each of the mutational groups (mutated/ unmutated)
#' @param expression_df_cancer a data frame with expression values (columns are 
#' patient IDs-sample IDs, rows are ENSG IDs), subsetted to only cancer samples
#' @param mutant_patients a vector of all the patients that have a mutation in 
#' the given gene-of-interest
#' @param ensg the ENSG ID of the target gene of interest
get_expression_by_mut_group <- function(expression_df_cancer, mutant_patients, 
                                        ensg) {
  expression_mutants <- as.numeric(expression_df_cancer[
    rownames(expression_df_cancer) == ensg, 
    colnames(expression_df_cancer) %fin% mutant_patients])
  expression_normal <- as.numeric(expression_df_cancer[
    rownames(expression_df_cancer) == ensg, 
    !(colnames(expression_df_cancer) %fin% mutant_patients)])
  return(list("mutants" = expression_mutants, "normal" = expression_normal))
}

#' Get the expression in each of the mutational groups (mutated/ unmutated), 
#' including normal samples as a separate category
#' @param expression_df_cancer a data frame with expression values (columns are 
#' patient IDs-sample IDs, rows are ENSG IDs), subsetted to only cancer samples
#' @param mutant_patients a vector of all the patients that have a mutation in 
#' the given gene-of-interest
#' @param ensg the ENSG ID of the target gene of interest
get_expression_by_mut_group_incl_norm <- function(expression_df, 
                                                  mutant_patients, ensg) {
  expression_normal <- as.numeric(expression_df[
    rownames(expression_df) == ensg, grepl("-11", colnames(expression_df))])
  expression_df_cancer <- expression_df[,!(grepl("-11", 
                                                 colnames(expression_df)))]
  print(head(expression_df_cancer))
  
  expression_mutants <- as.numeric(expression_df_cancer[
    rownames(expression_df_cancer) == ensg, 
    colnames(expression_df_cancer) %fin% mutant_patients])
  expression_wt <- as.numeric(expression_df_cancer[
    rownames(expression_df_cancer) == ensg, 
    !(colnames(expression_df_cancer) %fin% mutant_patients)])
  return(list("normal" = expression_normal, "cancer_wt" = expression_wt, 
              "cancer_mutants" = expression_mutants))
}

#' Visualize expression in mutation vs. normal patient groups using a boxplot
#' @param expression_df an expression data frame from process_expression_data.R
#' @param mutant_patients a vector of TCGA IDs of patients with mutations in the 
#' given gene-of-interest
#' @param goi the external name of the gene-of-interest whose mutation status
#' we are using as a comparator
#' @param ensg_target the ENSG ID of the target gene-of-interest
#' @param target the external name of the target gene-of-interest
display_expression_by_mut_group <- function(expression_df, mutant_patients, goi,
                                                       ensg_target, target) {
  expression_by_mut_group <- get_expression_by_mut_group(expression_df,
                                                         mutant_patients, 
                                                         ensg_target)
  expression_mutants <- expression_by_mut_group[[1]]
  expression_normal <- expression_by_mut_group[[2]]
  dataList_expression <- list("No Mutation" = expression_normal, 
                              "Mutation" = expression_mutants)
  bxp <- boxplot(dataList_expression, ylab = "Expression (Quantile-Normalized)", 
          main = paste(target, 
                       paste("Expression By", 
                             paste(goi, "Mutation Status"))))
  print(bxp)
}

# NOTE: option to take the intersection of "mutant_patients" with patients 
# that have a particular CNA or methylation status in the target gene or another
# gene; see example below
ensg_target <- "ENSG00000105281"
expression_mutants_krasmeth <- as.numeric(expression_df[
  rownames(expression_df) == ensg_target, 
  colnames(expression_df_quantile_norm) %fin% 
    intersect(idh1_mutant_patients, kras_methylated_patients)])
expression_mutants_notkrasmeth <- as.numeric(expression_df[
  rownames(expression_df) == ensg_target, colnames(expression_df) %fin% 
    setdiff(idh1_mutant_patients, kras_methylated_patients)])
expression_normal_krasmeth <- as.numeric(expression_df[
  rownames(expression_df) == ensg_target, 
  !((colnames(expression_df) %fin% idh1_mutant_patients) | 
      (colnames(expression_df) %fin% kras_methylated_patients))])
expression_normal_notkrasmeth <- as.numeric(expression_df[
  rownames(expression_df) == ensg_target, 
  !(colnames(expression_df) %fin% idh1_mutant_patients) & 
    (colnames(expression_df) %fin% kras_methylated_patients)])  
boxplot(list("No Mutation, KRAS Methylated" = expression_normal_krasmeth, 
             "No Mutation, KRAS Not Methylated" = expression_normal_notkrasmeth, 
             "Mutation, KRAS Methylated" = expression_mutants_krasmeth, 
             "Mutation, KRAS Not Methylated" = expression_mutants_notkrasmeth), 
        ylab = "Expression", 
        main = paste("KRAS", paste("Expression By", paste("IDH1", 
                                                          "Mutation Status"))))

############################################################
### VISUALIZE EXPRESSION OF TARGET GENE, BY FUNCTIONAL CNA NUMBER
### OF REGULATORY GENE OF INTEREST
############################################################
#' Get the expression in each of the functional copy number groups 
#' @param expression_df_cancer a data frame with expression values (columns are 
#' patient IDs-sample IDs, rows are ENSG IDs), subsetted to only cancer samples
#' @param funct_copy_dict a two-column data frame with sample ID and # of 
#' functional copies, from CNA_helper_functions.R
#' @param ensg the ENSG ID of the target gene of interest
#' @param ts_or_onco a character value indicating whether or not this is a 
#' "tumor_suppressor" or an "oncogene"
get_expression_by_fnc <- function(expression_df_cancer, funct_copy_dict, ensg, 
                                  ts_or_onco) {
  funct_copy_list <- list()
  num_funct_copy_options <- length(unique(funct_copy_dict$Num.Functional.Copies))
  
  for (i in 1:num_funct_copy_options) {
    samp_ids <- funct_copy_dict[funct_copy_dict$Num.Functional.Copies == (i-1), 
                                'sample.id']
    funct_copy_list[[i]] <- as.numeric(express_df_cancer[
      rownames(express_df_cancer) == ensg, 
      colnames(express_df_cancer) %fin% samp_ids])
  }
  
  # Tumor suppressor
  if(ts_or_onco == "tumor_suppressor") {
    grouped_list <- list("Full KO" = funct_copy_list[[1]], 
                         "Partial KO" = funct_copy_list[[2]],
                         "Normal/Amp" = unlist(funct_copy_list[
                           3:length(funct_copy_list)]))
  # Oncogene
  } else if (ts_or_onco == "oncogene") {
    grouped_list <- list("Full or Partial KO" = c(funct_copy_list[[1]], 
                                                  funct_copy_list[[2]]), 
                         "Normal" = funct_copy_list[[3]],
                         "Amplification" = unlist(funct_copy_list[
                           4:length(funct_copy_list)]))
    
  } else {print("Must supply either 'tumor_suppressor' or 'oncogene'.")}
  return(grouped_list)
}


############################################################
### VISUALIZE EXPRESSION OF TARGET GENE, BY MUTATION
### OF REGULATORY GENE OF INTEREST, FOR PAIRED CANCER-NORMAL
### EXPRESSION DATA
############################################################
#' Visualization of a paired analysis of expression in cancer / expression in 
#' normal for each driver mutational category. Also calculates and displays 
#' results of a paired Wilcoxon test between these mutational categories.
#' @param mutant_patients vector of patients that have a mutation in the 
#' regulatory gene of interest 
#' @param ensg the ENSG ID of the target gene of interest
#' @param targ_name the name of the target gene of interest
#' @param goi_name the gene name of the regulatory gene of interest
#' @param expression_df the expression data frame, limited to only patients in 
#' a given cancer type or subtype of interest
make_paired_expr_by_mut_group_vis <- function(mutant_patients, ensg, targ_name, 
                                              goi_name, expression_df) {
  
  expression_df <- expression_df[rownames(expression_df) == ensg, ]
  
  # Get normal samples
  samples_normal <- colnames(expression_df)[grepl("-11", 
                                                  colnames(expression_df))]
  patients_normal <- unlist(lapply(samples_normal, function(x) 
    unlist(strsplit(x, "-", fixed = T))[1]))
  
  # Get cancer samples
  samples_cancer <- colnames(expression_df)[grepl("-0", 
                                                  colnames(expression_df))]
  patients_cancer <- unlist(lapply(samples_cancer, function(x) 
    unlist(strsplit(x, "-", fixed = T))[1]))
  
  # Get the patients that have both tumor and normal samples
  intersecting_pats <- intersect(patients_normal, patients_cancer)
  print(length(intersecting_pats))
  
  # Limit to only those patients with both tumor and normal samples
  samples_normal_sub <- unlist(lapply(samples_normal, function(x) 
    ifelse(unlist(strsplit(x, "-", fixed = T))[1] %fin% 
             intersecting_pats, x, NA)))
  samples_normal_sub <- samples_normal_sub[!is.na(samples_normal_sub)]
  
  samples_cancer_sub <- unlist(lapply(samples_cancer, function(x) 
    ifelse(unlist(strsplit(x, "-", fixed = T))[1] %fin% 
             intersecting_pats, x, NA)))
  samples_cancer_sub <- samples_cancer_sub[!is.na(samples_cancer_sub)]
  
  cancer_norm_ratios_mut <- c()
  cancer_norm_ratios_wt <- c()
  paired_vals_df <- data.frame(matrix(nrow = length(intersecting_pats), 
                                      ncol = 3))
  colnames(paired_vals_df) <- c("cancer", "normal", "status")
  
  # For each pair of samples, get the ratio of cancer to normal expression and
  # store this with information on whether that patient has a mutation in the goi
  for(i in 1:length(samples_normal_sub)) {
    samp_norm <- samples_normal_sub[i]
    pat <- unlist(strsplit(samp_norm, "-", fixed = T))[1]
    samp_cancer <- unique(samples_cancer_sub[grepl(pat, samples_cancer_sub)])
    if(length(samp_cancer) > 1) {
      samp_cancer <- samp_cancer[grepl("-01", samp_cancer)][1]
    }
    
    expression_normal <- as.numeric(expression_df[, colnames(expression_df) == 
                                                    samp_norm])
    expression_cancer <- as.numeric(expression_df[, colnames(expression_df) == 
                                                    samp_cancer])
    
    if(samp_cancer %fin% mutant_patients) {
      cancer_norm_ratios_mut <- c(cancer_norm_ratios_mut, 
                                  (expression_cancer / expression_normal))
      paired_vals_df[i,] <- c(expression_cancer, expression_normal, "Mut")
    } else {
      cancer_norm_ratios_wt <- c(cancer_norm_ratios_wt, 
                                 (expression_cancer / expression_normal))
      paired_vals_df[i,] <- c(expression_cancer, expression_normal, "WT")
    }
  }
  
  dataList <- list("TP53 WT (Cancer / Normal)" = cancer_norm_ratios_wt,
                   "TP53 Mut. (Cancer / Normal)" = cancer_norm_ratios_mut)
  
  # Visualize these ratios using a boxplot
  boxplot(dataList, ylab = "Cancer to Normal TMM Expression Ratio", 
          xlab = paste(goi_name, "Mutation Status"),
          main = paste(targ_name, paste("Tumor-to-Normal Expression Ratio By", 
                                        paste(goi_name, "Mutation Status"))))
  abline(h = 1, col = "red")
  
  wilcox_res <- wilcox.test(x = cancer_norm_ratios_wt, 
                            y = cancer_norm_ratios_mut, paired = FALSE, 
                            alternative = "less")
  print(wilcox_res)
  
  # Add results from Wilcoxon test to plot
  text(1, max(c(cancer_norm_ratios_wt, cancer_norm_ratios_mut)), 
       paste("Wilcoxon p-value: ", round(wilcox_res$p.value, digits = 4)))
  
  paired_wilcox_res <- wilcox.test(as.numeric(paired_vals_df$cancer), 
                                   as.numeric(paired_vals_df$normal), 
                                   paired = T, alternative = "greater")
  print(paired_wilcox_res)
  
  paired_wilcox_res_wt <- wilcox.test(
    as.numeric(paired_vals_df[paired_vals_df$status == "WT", 'cancer']), 
    as.numeric(paired_vals_df[paired_vals_df$status == "WT", 'normal']), 
    paired = T, alternative = "greater")
  print(paired_wilcox_res_wt)
  
  paired_wilcox_res_mut <- wilcox.test(
    as.numeric(paired_vals_df[paired_vals_df$status == "Mut", 'cancer']), 
    as.numeric(paired_vals_df[paired_vals_df$status == "Mut", 'normal']),
    paired = T, alternative = "greater")
  print(paired_wilcox_res_mut)
}

# Adjust headers of expression DF to be only patient-sample ID (XXXX-XX), 
# if needed
colnames(expression_df) <- unlist(lapply(colnames(expression), function(x) 
  paste(unlist(strsplit(x, "-", fixed = T))[3:4], collapse = "-")))

# Call function
make_paired_expr_by_mut_group_vis(mutant_patients, ensg, "SHMT2", "TP53",
                                  expression_df)

# Can also perform this analysis on particular subtypes; an example for breast 
# cancer shown below (subtype IDs stored in separate files)
# Import the list of unique patient IDs that have all data types 
patient_set <- read.table(paste0(PATH, "UNIQUE_PATIENT_IDS.txt"), 
                          header = T)[,1]
patient_set_lumA <- intersect(
  read.table(paste0(PATH, "Patient_Subsets/Luminal.A_patient_ids.txt"), 
             header = T)[,1], patient_set)
patient_set_lumB <- intersect(
  read.table(paste0(PATH, "Patient_Subsets/Luminal.B_patient_ids.txt"), 
             header = T)[,1], patient_set)
patient_set_lumAB <- intersect(
  read.table(paste0(PATH, "Patient_Subsets/Luminal.A.B_patient_ids.txt"), 
             header = T)[,1], patient_set)
patient_set_basal <- intersect(
  read.table(paste0(PATH, "Patient_Subsets/Basal_patient_ids.txt"), 
             header = T)[,1], patient_set)
patient_set_her2 <- intersect(
  read.table(paste0(PATH, "Patient_Subsets/HER2_patient_ids.txt"), 
             header = T)[,1], patient_set)


expression_df_pats <- unlist(lapply(colnames(expression_df), function(x)
  unlist(strsplit(x, "-", fixed = T))[1]))
expression_df_lumA <- expression_df[, which(expression_df_pats %fin% 
                                              patient_set_lumA)]
expression_df_lumB <- expression_df[, which(expression_df_pats %fin% 
                                              patient_set_lumB)]
expression_df_basal <- expression_df[, which(expression_df_pats %fin% 
                                               patient_set_basal)]
expression_df_her2 <- expression_df[, which(expression_df_pats %fin% 
                                              patient_set_her2)]

