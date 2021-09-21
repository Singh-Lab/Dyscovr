############################################################
### Process Linear Model Output
### Written By: Sara Camilli, May 2021
############################################################

# Given an output 'master file' from a linear model run, this
# file performs the following functions:
# 1. Basic Output Visualization (Beta Distribution, P-Value Distribution,
# Q-Q Plot, and SE Distribution)
# 2. Heat map of t-statistics (clustered) for cases with multiple regulatory proteins
# 3. Perform multiple hypothesis testing correction
# 4. Obtain significant correlations (q-value thresholding)
# 5. Visualize enrichment of top hit genes in CGC/ Vogelstein,
# and potentially other cancer gene lists
# 6. Annotate output with relevance to cancer using TFcancer database

# install.packages("gplots")
# install.packages("pheatmap")
# BiocManager::install("ComplexHeatmap")
library("gplots")
library("pheatmap")
library(ComplexHeatmap)
library(broom)
library(qvalue)
library(ggplot2)
library(stringr)
library(dplyr)
library(VennDiagram)
library("RColorBrewer")

# Path to output files
main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

# Path to where output figures should be saved
output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/BRCA/Linear Model/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/Pan-Cancer/Linear Model/"


# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


############################################################
############################################################
#### BASIC VISUALIZTION OF OUTPUT
############################################################
############################################################

############################################################
#### VISUALIZE BETA VALUE DISTRIBUTION
############################################################
#' Function plots a histogram to visualize the distribution
#' of regulatory protein r_i beta values across all linear models.
#' @param results_table a master DF produced from run_linear_model()
visualize_beta_distrib <- function(results_table) {
  betas <- results_table$estimate
  hist(betas, main = "Histogram of Beta Coefficient Values Across all Reg. Proteins",
       xlab = "Beta Coefficient Value", ylab = "Frequency")
}

# Call this function & save output
fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/Mutation, eQTL/Beta Distribution (I-Protein, TMM, bucketCNA, iciTotFrac).png", sep = "")
png(fn, width = 450, height = 350)
visualize_beta_distrib(master_df_mut)
dev.off()

fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/CNA, eQTL/Beta Distribution (I-Protein, TMM, bucketCNA, iciTotFrac).png", sep = "")
png(fn, width = 450, height = 350)
visualize_beta_distrib(master_df_cna)
dev.off()


############################################################
#### VISUALIZE P-VALUE DISTRIBUTION
############################################################
#' Function plots a histogram to visualize the distribution
#' of regulatory protein r_i beta values across all linear models.
#' @param results_table a master DF produced from run_linear_model()
visualize_pval_distrib <- function(results_table) {
  pvals <- results_table$p.value[!is.na(results_table$p.value) & 
                                   !is.infinite(results_table$p.value)]
  hist(pvals, main = "Histogram of p-Values Across all Reg. Proteins",
       xlab = "p-value", ylab = "Frequency", col = "blueviolet")
}

# Call this function & save output
fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/Mutation, eQTL/P-Value Distribution (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 450, height = 350)
visualize_pval_distrib(master_df_mut)
dev.off()

fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/CNA, eQTL/P-Value Distribution (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 450, height = 350)
visualize_pval_distrib(master_df_cna)
dev.off()

############################################################
#### VISUALIZE Q-Q PLOT
############################################################
#' Function plots a Q-Q plot to visualize the distribution of p-values
#' and assess whether they come from a uniform distribution
#' @param results_table a master DF produced from run_linear_model()
qqplot_pvals <- function(results_table) {
  qqnorm(results_table$p.value, pch = 1, frame = FALSE)
  qqline(results_table$p.value, col = "steelblue", lwd = 2)
}

# Call this function
#qqplot_pvals(master_df)


############################################################
#### VISUALIZE ERROR DISTRIBUTION
############################################################
#' Function plots a histogram of the standard errors produced from
#' all LM runs to assess whether the errors derive from a normal distribution
#' @param results_table a master DF produced from run_linear_model()
visualize_error_distrib <- function(results_table) {
  hist(results_table$std.error, main = "Standard Error Distribution Across All Tests",
       xlab = "Standard Error (SE)", ylab = "Frequency")
}

# Call this function
#visualize_error_distrib(master_df)



############################################################
############################################################
#### PEFORM MULTIPLE HYPOTHESIS TESTING CORRECTION
############################################################
############################################################
#' Function takes in an output results table and applies multiple
#' hypothesis testing correction (Storey's q-value correction) to 
#' all p-values in order to add a column of q-values. Returns 
#' the results table with a column for q-values.
#' @param results_table a master DF produced from run_linear_model()
mh_correct <- function(results_table) {
  
  # Get the qvalue object
  qobj <- qvalue(p = results_table$p.value)
  
  # OPT: plot some useful plots & print some useful information
  plot(qobj)
  print(summary(qobj))
  #print(paste("Pi0 (Propr. of true null hypotheses):", qobj$pi0))
  
  qvals <- qobj$qvalues # extract qvalues
  
  results_table$q.value <- qvals # add the qvalues back to the data frame
  
  # Return the tidied linear model fit with q-values
  return(results_table)
}

# Call this function
master_df_mut_corrected <- mh_correct(master_df_mut)
master_df_cna_corrected <- mh_correct(master_df_cna)

#Q-Value Visualization (I-Protein, FPKM, rawCNA, methBeta, iciTotFrac)

############################################################
############################################################
#### ADD GENE NAMES TO FILE
############################################################
############################################################
#' Given a master data frame result from model, add a column for target gene name
#' and regulatory protein gene name. Return the updated data frame.
#' @param master_df_sig a data.table object produced from the linear_model.R function
#' @param all_genes_id_conv a bioMart file with conversions between different gene ID types
add_targ_regprot_gns <- function(master_df_sig, all_genes_id_conv) {
  master_df_sig$T_k.name <- unlist(lapply(master_df_sig$T_k, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == unlist(strsplit(x, ";", fixed = TRUE))[1], 
                                   'external_gene_name']), collapse = ";")))
  
  # Add a column for the regulatory protein name
  master_df_sig$R_i.name <- unlist(lapply(master_df_sig$R_i, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == unlist(strsplit(x, ";", fixed = TRUE))[1], 
                                   'external_gene_name']), collapse = ";")))
  return(master_df_sig)
}

# Call this function
master_df_mut_corrected <- add_targ_regprot_gns(master_df_mut_corrected, all_genes_id_conv)
master_df_cna_corrected <- add_targ_regprot_gns(master_df_cna_corrected, all_genes_id_conv)

# Write this to a new file
fwrite(master_df_mut_corrected, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/eQTL/output_results_P53_chipeat_iprotein_TMM_rawCNA_cibersortTotFrac_corrected_MUT.csv", sep = ""))
fwrite(master_df_cna_corrected, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/eQTL/output_results_P53_chipeat_iprotein_TMM_rawCNA_cibersortTotFrac_corrected_CNA.csv", sep = ""))


############################################################
############################################################
#### HEAT MAP VISUALIZATION
############################################################
############################################################
#' Function plots a regulatory protein vs. gene target
#' clustered heat map, with entries being the t-statistic
#' of the hypothesis test
#' @param results_table a master DF produced from linear_model.R
#' @param outpath a path to a directory to write the t-statistic matrix to
create_heat_map <- function(results_table, outpath) {
  # Create the regulatory protein vs. gene target table
  matrix <- data.frame(matrix(ncol = length(unique(results_table$R_i.name)),
                              nrow = length(unique(results_table$T_k.name))))
  colnames(matrix) <- unique(results_table$R_i.name)
  rownames(matrix) <- unique(results_table$T_k.name)
  
  # Fill in this table with t-statistics
  for (i in 1:nrow(results_table)) {
    tstat <- results_table$statistic[i]
    regprot <- results_table$R_i.name[i]
    targ <- results_table$T_k.name[i]
    
    matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- tstat
  }
  # NOTE: leave all the unfilled pairings as NA
  
  # Write to file
  fwrite(matrix, paste(outpath, "Linear Model/Highly Deleted TFs/Non-Tumor-Normal Matched/tstatistic_matrix_delTFs_metabolicTargs_TMM_bucketCNA_iprot_iciTotFrac_uncorrected_MUT_RANDOMIZED.csv", sep = ""))
  
  # Convert from data frame to matrix
  matrix <- data.matrix(matrix)
  print(head(matrix))
  
  # Plot a default heatmap
  #col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
  #heatmap(matrix, scale = "none", col = col)
  
  # Plot an enhanced heatmap
  heatmap.2(matrix, scale = "none", col = bluered(100), trace = "none",
            density.info = "none") # clusters by default using hclust, but can specify others using param 'hclustfun'
  
  # Plot a pretty heatmap
  # pheatmap(matrix, cutree_rows = 4) # options are available for changing clustering metric & method
  # defaults to euclidean and complete
  
  # Plot a complex heatmap using Bioconductor
  #Heatmap(matrix, name = "Results Heatmap", column_title = "Regulatory Proteins",
          #row_title = "Gene Targets", row_names_gp = gpar(fontsize = 6))
  # Additional arguments: show_row_names, show_column_names, show_row_hclust, 
  # clustering_distance_rows, clustering_distance_columns (the metric for clustering,
  # e.g. "euclidean" or "pearson"),
  # clustering_method_rows, clustering_method_columns (the method for clustering, 
  # e.g. "complete" or "average")
}

# Call this function
fn <- paste(output_path, "ALL/Non-Tumor-Normal-Matched/Mutation, eQTL/Heatmap (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 550, height = 450)
create_heat_map(master_df_mut_corrected, main_path)
dev.off()

fn <- paste(output_path, "ALL/Non-Tumor-Normal-Matched/CNA, eQTL/Heatmap (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 550, height = 450)
create_heat_map(master_df_cna_corrected, main_path)
dev.off()


############################################################
############################################################
#### OBTAIN SIGNIFICANT CORRELATIONS
############################################################
############################################################
#' Function takes in an output results table with q-values and 
#' restricts it to only models that exceed the q-value threshold
#' (are statistically significant correlations). Then ranks the 
#' remaining by q-values and returns the ranked top hits list.
#' @param results_table a master DF produced from run_linear_model() that has q-values added
#' from the mh_correct() function
#' @param qval_thres a threshold for significance for q-values
get_signif_correl <- function(results_table, qval_thres) {
  
  # Limit to only entries that exceed the given qvalue threshold
  results_table_sig <- results_table %>% filter(q.value < qval_thres)
  
  # Sort the table by qvalue
  results_table_sig_ordered <- results_table_sig[order(results_table_sig$q.value, 
                                                       decreasing = FALSE),]
  
  return(results_table_sig_ordered)
}

# Call this function
qval_thres <- 0.1
master_df_mut_sig <- get_signif_correl(master_df_mut_corrected, qval_thres)
master_df_cna_sig <- get_signif_correl(master_df_cna_corrected, qval_thres)

# Write these results to a new file
fwrite(master_df_mut_sig, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/eQTL/significant_output_P53_chipeat_iprotein_tmm_rawCNA_methBetaRaw_cibersortTotFrac_MUT.csv", sep = ""))                                                                                                                     # 118 (including CNAi and Methi as covar.)
fwrite(master_df_cna_sig, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/eQTL/significant_output_P53_chipeat_iprotein_tmm_rawCNA_methBetaRaw_cibersortTotFrac_CNA.csv", sep = ""))                                                                                                                     # 118 (including CNAi and Methi as covar.)




############################################################
############################################################
#### PLOT OVERLAP BETWEEN MUTATION & CNA TOP HITS
############################################################
############################################################
#' Plots the overlap in significant top target gene hits for mutation 
#' and CNA results. 
#' @param fn an output path/filename for the visualization to be saved to
#' @param master_df_mut_sig mutation master DF, with gene names added and 
#' thresholded for significance
#'@param master_df_cna_sig CNA master DF, with gene names added and 
#' thresholded for significance
plot_tophit_overlap <- function(fn, master_df_mut_sig, master_df_cna_sig) {
  overlap_genes <- intersect(master_df_mut_sig$T_k.name, master_df_cna_sig$T_k.name)
  print(overlap_genes)
  
  # Plot a Venn Diagram
  #myCol <- brewer.pal(2, "Pastel2")
  venn.diagram(list(master_df_mut_sig$T_k.name, master_df_cna_sig$T_k.name),
               category.names = c("Mutation", "CNA"), filename = fn, output = TRUE,
               lwd = 2, lty = 'blank', fill = c("red", "blue"), cex = 0.6, fontface = "bold",
               fontfamily = "sans")#, cat.cex = 0.6, cat.fontface = "bold",
               #cat.default.pos = "outer", cat.fontfamily = "sans", rotation = 1)
}

# Call this function and write to PNG file
fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/Top Gene Hit Overlap (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 450, height = 350)
plot_tophit_overlap(fn, master_df_mut_sig, master_df_cna_sig)
dev.off()

############################################################
############################################################
#### VISUALIZE ENRICHMENT OF TOP HITS IN CGC/ VOGELSTEIN/ etc.
############################################################
############################################################

############################################################
#### IMPORT CANCER GENE LISTS 
############################################################
# Import a table containing a compiled list of known cancer genes
known_cancer_genes_table <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/GRCh38_driver_gene_list.tsv", sep = "\t",
                                       header = TRUE, check.names = FALSE, comment.char = "#")

# Import TF cancer data table(s)
tfcancer_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/"

# FOR BRCA
tfcancer_df <- read.csv(paste(tfcancer_path, "BRCA Data/Validation_Files/BRCA_TFcancer.csv", sep = ""),
                        header = TRUE, check.names = FALSE)

# FOR PAN-CANCER (must combine files)
#' Creates a TFcancer file that is the combination of all the individual cancer
#' type TFcancer files (33 cancer types). Adds a column for the cancer type.
#' @param tfcancer_path local path to the TF cancer files
merge_tfcancer_dfs <- function(tfcancer_path) {
  # List all the files in the path
  tfcancer_files <- list.files(paste(tfcancer_path, "TCGA Data (ALL)/Validation_Files/TFcancer/"))
  
  # Read in the first DF 
  tfcancer_df <- read.csv(paste(tfcancer_path, tfcancer_files[1], sep = ""), header = TRUE,
                          check.names = FALSE)
  tfcancer_df$cancer_type <- rep(unlist(strsplit(tfcancer_files[1], ".", fixed = TRUE))[1], 
                                 nrow(tfcancer_df))
  # Merge the rest into this DF with a column for cancer type
  for (i in 2:length(tfcancer_files)) {
    new_tfcancer_df <- read.csv(paste(tfcancer_path, tfcancer_files[i], sep = ""), header = TRUE,
                                check.names = FALSE)
    new_tfcancer_df$cancer_type <- rep(unlist(strsplit(tfcancer_files[i], ".", fixed = TRUE))[1], 
                                       nrow(new_tfcancer_df))
    tfcancer_df <- rbind(tfcancer_df, new_tfcancer_df)
  }
  return(tfcancer_df)
}

# Call function for pan-cancer
#tfcancer_df <- merge_tfcancer_dfs(tfcancer_path)

# Write this back to a CSV
#fwrite(tfcancer_df, paste(path, "Validation_Files/TFcancer/ALL_CANCERS.csv", sep = ""))


############################################################
#### ADJUST TFCANCER DF
############################################################
#' Function adjusts the TFcancer data frame to 1. make all the
#' gene names uppercase; 2. separate out the semicolon-separated
#' gene target names into individual rows
#' @param tfcancer_df the imported data frame with the information
#' from the TFcancer database
adjust_tfcancer_df <- function(tfcancer_df) {
  # First, break up all the comma-separated targets into new rows
  new_dfs <- lapply(1:nrow(tfcancer_df), function(i) {
    # Get the target genes into a vector with whitespace stripped
    targ_genes <- unlist(strsplit(stringr::str_replace_all(tfcancer_df[i, 'gene'], fixed(" "), ""), 
                                  ";", fixed = TRUE))
    targ_genes <- targ_genes[!is.na(targ_genes)]
    print(targ_genes)
    
    # For each, create a data frame row with everything duplicated from the original
    # DF except for the gene column. Return a new data frame with all these rows.
    samp_df <- tfcancer_df[i,]
    samp_df$gene <- targ_genes[1]
    for (j in 2:length(targ_genes)) {
      samp_df <- rbind(samp_df, tfcancer_df[i,])
      samp_df[j, 'gene'] <- targ_genes[j]
    }
    return(samp_df)
  })
  
  # Re-combine all these partial DFs back into one
  new_tfcancer_df <- do.call(rbind, new_dfs)
  
  return(new_tfcancer_df)
}

# Call this function
#new_tfcancer_df <- adjust_tfcancer_df(tfcancer_df)

# Write to a new file

# For BRCA
#fwrite(new_tfcancer_df, paste(tfcancer_path, "BRCA Data/Validation_Files/BRCA_TFcancer_adj.csv", sep = ""))
# For pan-cancer
#fwrite(new_tfcancer_df, paste(path, "Validation_Files/TFcancer/ALL_CANCERS_adj.csv", sep = ""))


############################################################
#### CREATE ENRICHMENT PLOT
############################################################
#' Function to plot the enrichment of cancer genes among the significant target genes
#' @param master_df_sig a master DF produced from run_linear_model() that has q-values
#' and has been thresholded to only those pairings that exceed a significance threshold
#' @param known_cancer_genes_table a data frame containing a compiled list of known 
#' cancer genes (CGC, Vogelstein, etc.) with various types of IDs
#' @param tfcancer_df a data frame from TFcancer, which has lists of transcription
#' factors shown to be associated with cancer from curated literature
plot_cancer_enrichment <- function(master_df_sig, known_cancer_genes_table, tfcancer_df) {
  # Get all the significant target genes 
  significant_target_hits <- unlist(master_df_sig$T_k.name)
  print(significant_target_hits)
  
  # Fill in a vector of 0 and 1 for each significant target hit to indicate if it 
  # is a known cancer gene
  cancer_vect <- unlist(lapply(significant_target_hits, function(x) {
    val1 <- ifelse(x %fin% known_cancer_genes_table$primary_gene_names, 1, 0)
    val2 <- ifelse(x %fin% tfcancer_df$gene, 1, 0)
  }))
  
  # Plot the enrichment (fraction of CGC genes at/ above given rank)
  ranks <- 1:length(significant_target_hits)
  frac_cancer_vect <- c()
  count_of_cancer_genes <- 0
  for (i in 1:length(significant_target_hits)) {
    val_of_curr_tg <- cancer_vect[i]
    count_of_cancer_genes <- count_of_cancer_genes + val_of_curr_tg
    frac <- count_of_cancer_genes / i
    frac_cancer_vect <- c(frac_cancer_vect, frac)
  }
  
  plot(ranks, frac_cancer_vect, main = "Enrichment of Significant Target Genes in Known Cancer Genes",
       xlab = "Rank", ylab = "Fraction of Target Genes that are Known Cancer Genes")
}

# Run function
fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/Mutation, eQTL/Enrichment in Known Cancer Genes (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 650, height = 350)
plot_cancer_enrichment(master_df_mut_sig, known_cancer_genes_table, tfcancer_df)
dev.off()

fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/CNA, eQTL/Enrichment in Known Cancer Genes (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 650, height = 350)
plot_cancer_enrichment(master_df_cna_sig, known_cancer_genes_table, tfcancer_df)
dev.off()



############################################################
#### ANNOTATE WHETHER TF REGULATION PREDICTIONS ALIGN WITH
#### TFCANCER FINGDINGS
############################################################
#' Add columns to the master results data frame to signify 
#' whether the prediction aligns with something found in the
#' TFcancer database. If the TF in question is found in the
#' database, and the gene in question is found to be 'regulated
#' by' that gene, then the regulation type is given, along with
#' all the corresponding information. Similarly, if a gene in 
#' question is found to be in the database, and the TF in 
#' question is found to be 'targeted by' that gene, the 
#' information will also be given. Adjusted master DF is returned.
#' @param master_df_sig the master results DF to be annotated
#' @param tfcancer_df the TFcancer data frame with the information
#' to be used for annotation
annotate_master_w_tfcancer <- function(master_df_sig, tfcancer_df) {
  # Look at each top hit individually, creating a new DF to cbind
  # to the master DF
  tfcancerdf_list <- lapply(1:nrow(master_df_sig), function(i) {
    #tf <- master_df_sig[i, 'R_i.name']
    tf <- "p53"
    gene <- master_df_sig[i, 'T_k.name']
    
    # Subset the TFcancer database to see if the TF or gene are in there
    tf_sub <- tfcancer_df[tfcancer_df$tf == tf,]
    gene_sub <- tfcancer_df[tfcancer_df$gene == gene,]
    
    # Row entries to fill in
    df_to_return <- data.frame(tf_present = FALSE, gene_present= FALSE, any_relationship = FALSE,
                               characteristics = "", regulation = "", hallmark = "", 
                               original_text = "")
    
    if (!(nrow(tf_sub) == 0)) {df_to_return$tf_present <- TRUE}
    if (!(nrow(gene_sub) == 0)) {df_to_return$gene_present <- TRUE}
    
    both_sub <- tf_sub[grepl(gene, tf_sub$gene),]
    if (!nrow(both_sub) == 0) {
      df_to_return$any_relationship <- TRUE
      df_to_return$characteristics <- both_sub$characteristics
      df_to_return$regulation <- both_sub$regulation
      df_to_return$hallmark <- both_sub$hallmark
      df_to_return$original_text <- both_sub$original_text
    }
    
    return(df_to_return)
  })
  
  # Rbind all these DFs together
  tfcancer_res_df <- do.call(rbind, tfcancerdf_list)
  
  # Add this DF to the original master results DF and return
  master_df_sig <- cbind(master_df_sig, tfcancer_res_df)
  
  return(master_df_sig)
}

# Run function
#master_df_sig <- annotate_master_w_tfcancer(master_df_sig, tfcancer_df)

# Write to a file
#fwrite(master_df_sig, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/iprotein_significant_output_TP53_TMM_rawCNA_iprot_iciTotFrac.csv", sep = ""))


#### NOTE: ANOTHER OPTION IS TO COPY THE TOP HIT GENES INTO GENEONTOLOGY: http://geneontology.org/


############################################################
#### ASIDE: DETERMINE WHAT COVARIATES ARE IMPORTANT
############################################################

# Read back a master DF
master_df <- fread(paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/iprotein_output_results_TP53_alltargs_TMM_bucketCNA_iprot_iciTotFrac_uncorrected.csv", sep = ""), 
                   header = TRUE)

# Remove the (Intercept) terms
master_df <- master_df[master_df$term != "(Intercept)",]

# Of the top X remaining tests, get what the terms are and plot the proportions
x <- 500
master_df_topx <- master_df[1:x,]
terms <- unique(master_df_topx$term)
terms_counts <- unlist(lapply(terms, function(x) nrow(master_df_topx[master_df_topx$term == x,])))
terms_counts_df <- data.frame('term' = terms, 'freq' = terms_counts)
pie(terms_counts_df$freq, labels = terms_counts_df$term, main = paste("Most Significant Covariates (Top", paste(x, "from All Tests)")))


# Get the proportions of all the tests with p-value <0.05
master_df_sig <- master_df[master_df$p.value < 0.05,]
terms <- unique(master_df_sig$term)
terms_counts <- unlist(lapply(terms, function(x) nrow(master_df_sig[master_df_sig$term == x,])))
terms_counts_df <- data.frame('term' = terms, 'freq' = terms_counts)
pie(terms_counts_df$freq, labels = terms_counts_df$term, main = "Categories of Significant Covariates (All Tests)")

