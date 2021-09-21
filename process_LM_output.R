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

#!/usr/bin/env Rscript

# install.packages("gplots")
# install.packages("pheatmap")
# BiocManager::install("ComplexHeatmap")
library("gplots")
#library("pheatmap")
#library(ComplexHeatmap)
library(broom)
library(qvalue)
library(ggplot2)
library(stringr)
library(dplyr)
library("RColorBrewer")


############################################################
# SET PATH
############################################################
# Path to where output figures should be saved
output_vis_path <- paste(source_path, paste("output_visualizations/", args$cancerType, sep = ""), sep = "")
output_vis_path <- paste(output_vis_path, args$tester_name, sep = "/")
if(tumNormMatched) {
  output_vis_path <- paste(output_vis_path, "tumor_normal_matched", sep = "/")
} else {output_vis_path <- paste(output_vis_path, "tumor_only", sep = "/")}
output_vis_path <- paste(output_vis_path, args$QTLtype, sep = "/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/all_genes_id_conv.csv", header = TRUE)






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
       xlab = "Beta Coefficient Value", ylab = "Frequency", col = "darkseagreen2")
}


# Create specific filenames
fn_mut <- paste(output_vis_path, paste("Beta Distrib. (", paste(paste(outfn, "_MUT", sep = ""), ").png", sep = ""), sep = ""), sep = "/")
fn_cna <- paste(output_vis_path, paste("Beta Distrib. (", paste(paste(outfn, "_CNA", sep = ""), ").png", sep = ""), sep = ""), sep = "/")

# Call this function & save output
png(fn_mut, width = 450, height = "350")
visualize_beta_distrib(master_df_mut)
dev.off()

png(fn_cna, width = 450, height = "350")
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

# Create specific filenames
fn_mut <- paste(output_vis_path, paste("P-Value Distrib. (", paste(paste(outfn, "_MUT", sep = ""), ").png", sep = ""), sep = ""), sep = "/")
fn_cna <- paste(output_vis_path, paste("P-Value Distrib. (", paste(paste(outfn, "_CNA", sep = ""), ").png", sep = ""), sep = ""), sep = "/")

# Call this function & save output
png(fn_mut, width = 450, height = "350")
visualize_pval_distrib(master_df_mut)
dev.off()

png(fn_cna, width = 450, height = "350")
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
#' @param fn_qvalvis a filename with path for a q-value visualization
#' produced from the qvalue package
#' @param fn_qvalsum a filename with path for a q-value summary 
#' with useful information produced from the qvalue package
mh_correct <- function(results_table, fn_qvalvis, fn_qvalsum) {
  
  # Get the qvalue object
  qobj <- qvalue(p = results_table$p.value)
  
  # Plot some useful plots & print some useful information
  png(fn_qvalvis, width = 450, height = "350")
  plot(qobj)
  dev.off()
  
  png(fn_qvalsum, width = 450, height = "350")
  print(summary(qobj))
  dev.off()
  
  #print(paste("Pi0 (Propr. of true null hypotheses):", qobj$pi0))
  
  qvals <- qobj$qvalues # extract qvalues
  
  results_table$q.value <- qvals # add the qvalues back to the data frame
  
  # Return the tidied linear model fit with q-values
  return(results_table)
}

# Call this function

fn_mut_qvalvis <- paste(output_vis_path, paste("Q-Value Visualization (", paste(paste(outfn, "_MUT", sep = ""), ").png", sep = ""), sep = ""), sep = "/")
fn_cna_qvalvis <- paste(output_vis_path, paste("Q-Value Visualization (", paste(paste(outfn, "_CNA", sep = ""), ").png", sep = ""), sep = ""), sep = "/")
fn_mut_qvalsum <- paste(output_vis_path, paste("Q-Value Summary (", paste(paste(outfn, "_MUT", sep = ""), ").png", sep = ""), sep = ""), sep = "/")
fn_cna_qvalsum <- paste(output_vis_path, paste("Q-Value Summary (", paste(paste(outfn, "_CNA", sep = ""), ").png", sep = ""), sep = ""), sep = "/")

master_df_mut_corrected <- mh_correct(master_df_mut, fn_mut_qvalvis, fn_mut_qvalsum)
master_df_cna_corrected <- mh_correct(master_df_cna, fn_cna_qvalvis, fn_cna_qvalsum)

# Write this to a new file
outfn <- str_replace(outfn, "uncorrected", "corrected") 
fwrite(master_df_mut_corrected, paste(outpath, paste(outfn, paste("_MUT", ".csv", sep = ""), sep = ""), sep = "/"))
fwrite(master_df_cna_corrected, paste(outpath, paste(outfn, paste("_CNA", ".csv", sep = ""), sep = ""), sep = "/"))

#fwrite(master_df_mut_corrected, paste(main_path, "TP53/Non-Tumor-Normal Matched/iprotein_output_results_TP53_TMM_bucketCNA_iprot_iciTotFrac_MUT_corrected.csv", sep = ""))
#fwrite(master_df_cna_corrected, paste(main_path, "TP53/Non-Tumor-Normal Matched/iprotein_output_results_TP53_TMM_bucketCNA_iprot_iciTotFrac_CNA_corrected.csv", sep = ""))


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


############################################################
############################################################
#### ADD GENE NAMES TO FINAL FILE VERSION
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
master_df_mut_sig <- add_targ_regprot_gns(master_df_mut_sig, all_genes_id_conv)
master_df_cna_sig <- add_targ_regprot_gns(master_df_cna_sig, all_genes_id_conv)


# Write these results to a new file
outfn <- str_replace(outfn, "corrected", "") 
outfn <- str_replace(outfn, "output_results", "significant_output")
if(length(master_df_mut_sig) > 0) {
  fwrite(master_df_mut_sig, paste(outpath, paste(outfn, paste("_MUT", ".csv", sep = ""), sep = ""), sep = "/"))
}
if(length(master_df_cna_sig) > 0) {
  fwrite(master_df_cna_sig, paste(outpath, paste(outfn, paste("_CNA", ".csv", sep = ""), sep = ""), sep = "/"))
}

#fwrite(master_df_mut_sig, paste(main_path, "TP53/Non-Tumor-Normal Matched/iprotein_significant_output_TP53_TMM_bucketCNA_iprot_iciTotFrac_MUT.csv", sep = ""))                                                                                                                     # 118 (including CNAi and Methi as covar.)
#fwrite(master_df_cna_sig, paste(main_path, "TP53/Non-Tumor-Normal Matched/iprotein_significant_output_TP53_TMM_bucketCNA_iprot_iciTotFrac_CNA.csv", sep = ""))                                                                                                                     # 118 (including CNAi and Methi as covar.)

############################################################
############################################################
#### HEAT MAP VISUALIZATION
############################################################
############################################################
#' Function plots a regulatory protein vs. gene target
#' clustered heat map, with entries being the t-statistic
#' of the hypothesis test
#' @param results_table
create_heat_map <- function(results_table) {
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
    
    tryCatch({
      matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- tstat
    }, error=function(cond){
      print(cond)
      matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- 0
    })
  }
  # NOTE: leave all the unfilled pairings as NA
  
  # Convert from data frame to matrix
  matrix <- data.matrix(matrix)
  
  # Replace NA/NaN with 0
  matrix[is.na(matrix)] <- 0
  matrix[is.nan(matrix)] <- 0
  
  print(head(matrix))
  
  # Plot a default heatmap
  #col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
  #heatmap(matrix, scale = "none", col = col)
  
  # Plot an enhanced heatmap
  heatmap.2(matrix, scale = "none", col = bluered(100), trace = "none",
            density.info = "none") # clusters by default using hclust, but can specify others using param 'hclustfun'
  
  # Plot a pretty heatmap
  #pheatmap(matrix, cutree_rows = 4) # options are available for changing clustering metric & method
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

if(!test) {
  fn_mut <- paste(output_vis_path, paste("Heatmap (", paste(paste(outfn, "_MUT", sep = ""), ").png", sep = ""), sep = ""), sep = "/")
  fn_cna <- paste(output_vis_path, paste("Heatmap (", paste(paste(outfn, "_CNA", sep = ""), ").png", sep = ""), sep = ""), sep = "/")
  
  # Call this function
  png(fn_mut, width = 450, height = "350")
  create_heat_map(master_df_mut_corrected)
  dev.off()
  
  png(fn_cna, width = 450, height = "350")
  create_heat_map(master_df_cna_corrected)
  dev.off()
}


############################################################
############################################################
#### VISUALIZE ENRICHMENT OF TOP HITS IN CGC/ VOGELSTEIN/ etc.
############################################################
############################################################

############################################################
#### IMPORT CANCER GENE LISTS 
############################################################
# Import a table containing a compiled list of known cancer genes
known_cancer_genes_table <- read.table("/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files//GRCh38_driver_gene_list.tsv", sep = "\t",
                                       header = TRUE, check.names = FALSE, comment.char = "#")

# Import TF cancer data table(s), if using
#tfcancer_path <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/"
#tfcancer_path <- paste(tfcancer_path, args$cancerType, sep = "")

#if(args$cancerType == "BRCA") {
#tfcancer_df <- read.csv(paste(tfcancer_path, "BRCA_TFcancer.csv", sep = ""),
#header = TRUE, check.names = FALSE)
#} else {
#tfcancer_df <- read.csv(paste(tfcancer_path, "Pan-Cancer_TFcancer.csv", sep = ""),
#header = TRUE, check.names = FALSE)
#}

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
#fwrite(tfcancer_df, paste(main_path, "Validation_Files/TFcancer/ALL_CANCERS.csv", sep = ""))


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
  print(head(significant_target_hits))
  
  # Fill in a vector of 0 and 1 for each significant target hit to indicate if it 
  # is a known cancer gene
  cancer_vect <- unlist(lapply(significant_target_hits, function(x) {
    val1 <- ifelse(x %fin% known_cancer_genes_table$primary_gene_names, 1, 0)
    if(!missing(tfcancer_df)) {val2 <- ifelse(x %fin% tfcancer_df$gene, 1, 0)}
    if(val1 | val2) {return(1)}
    else {return(0)}
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


fn_mut <- paste(output_vis_path, paste("Known Cancer Gene Enrichment (", 
                                       paste(paste(outfn, "_MUT", sep = ""), ").png", sep = ""), sep = ""), sep = "/")
fn_cna <- paste(output_vis_path, paste("Known Cancer Gene Enrichment (", 
                                       paste(paste(outfn, "_CNA", sep = ""), ").png", sep = ""), sep = ""), sep = "/")


# Run function
png(fn, width = 450, height = "350")
plot_cancer_enrichment(master_df_mut_sig, known_cancer_genes_table)
#plot_cancer_enrichment(master_df_mut_sig, known_cancer_genes_table, tfcancer_df)
dev.off()

png(fn, width = 450, height = "350")
plot_cancer_enrichment(master_df_cna_sig, known_cancer_genes_table)
#plot_cancer_enrichment(master_df_cna_sig, known_cancer_genes_table, tfcancer_df)
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


############################################################
#### ASIDE: INVESTIGATE WHICH RUNS HAVE THE MOST EXTREME 
#### VALUES
############################################################

# The case: we are seeing outliers in the Beta values that are skewing the values 
# non-normal. We want to see if this problem is widespread and under what conditions
# the Beta values are the most well-behaved.

# Path to location out output files (we'll be looking at the uncorrected files).
output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Linear Model/TP53/Non-Tumor-Normal Matched/"

# Import the uncorrected full, MUT and CNA results of a few sample files
master_df <- read.csv(paste(output_path, "eQTL/output_results_P53_chipeat_iprotein_tmm_CNAbucket_inclAmp_methBetaBucketed_cibersortTotalFrac_uncorrected.csv", sep = ""), header = TRUE, check.names = FALSE)
master_df_mut <- read.csv(paste(output_path, "eQTL/output_results_P53_chipeat_iprotein_tmm_CNAbucket_inclAmp_methBetaBucketed_cibersortTotalFrac_uncorrected_MUT.csv", sep = ""), header = TRUE, check.names = FALSE)
master_df_cna <- read.csv(paste(output_path, "eQTL/output_results_P53_chipeat_iprotein_tmm_CNAbucket_inclAmp_methBetaBucketed_cibersortTotalFrac_uncorrected_CNA.csv", sep = ""), header = TRUE, check.names = FALSE)

# Test for normality of Betas -- shapiro.test only valid up to 5000 samples
ks.test(master_df_mut$estimate, y = 'pnorm', alternative = 'two.sided')
ks.test(master_df_cna$estimate, y = 'pnorm', alternative = 'two.sided')
ks.test(master_df$estimate, y = 'pnorm', alternative = 'two.sided')

# All results files
results_files_mut_eQTL <- list.files(paste(output_path, "eQTL/", sep = ""), pattern = "uncorrected_MUT.csv")
results_files_cna_eQTL <- list.files(paste(output_path, "eQTL/", sep = ""), pattern = "uncorrected_CNA.csv")

#' Get the top p-values and return a ranked list of these files with the 
#' corresponding p-value and the top gene
#' @param results_files a vector of the names of the results files (uncorrected)
#' @param path a local path to the results files directory
get_files_with_outliers <- function(results_files, path) {
  output_tab <- data.frame(matrix(nrow = length(results_files), ncol = 8))
  
  switch <- 0
  for (i in 1:length(results_files)) {
    # import the file
    filename <- results_files[i]
    file <- fread(paste0(path, filename), header = TRUE)
    
    # if the first file, get the column names
    if(switch == 0) {
      colnames(output_tab) <- c("file.name", colnames(file)[2:ncol(file)])
      switch <- 1
    }
    
    # Get the first entry and add to the output table with the file name
    fn_sub <- unlist(strsplit(filename, "results_", fixed = TRUE))[2]
    
    output_tab[i,] <- c(fn_sub, file[1, 2:ncol(file)])
  }
  return(output_tab)
}

mut_files_with_outliers <- get_files_with_outliers(results_files_mut_eQTL, 
                                                   paste0(output_path, "eQTL/"))
cna_files_with_outliers <- get_files_with_outliers(results_files_cna_eQTL, 
                                                   paste0(output_path, "eQTL/"))

# To look at PEER factor trends
mut_peer_vals <- mut_files_with_outliers[grepl("PEER", mut_files_with_outliers$file.name), 'p.value']
mut_peer_vals <- c(mut_peer_vals[length(mut_peer_vals)], mut_peer_vals[1:(length(mut_peer_vals)-1)])
plot(mut_peer_vals)
plot(log(mut_peer_vals))

# Check which have Beta distributions that are non-normal
get_nonnormal_betaDistr <- function(results_files, path) {
  list_nonnorm_betas <- c()
  
  for (i in 1:length(results_files)) {
    # import the file
    filename <- results_files[i]
    file <- fread(paste0(path, filename), header = TRUE)
    print(head(file))
    
    # test for normality of the Betas
    ks_res <- ks.test(file$estimate, y = 'pnorm', alternative = 'two.sided')
    print(ks_res$p.value)
    if(ks_res$p.value < 1*10^(-20)) {
      fn_sub <- unlist(strsplit(filename, "results_", fixed = TRUE))[2]
      list_nonnorm_betas <- c(list_nonnorm_betas, fn_sub)
    }
  }
  return(list_nonnorm_betas)
}

mut_files_nonnormalBetas <- get_nonnormal_betaDistr(results_files_mut_eQTL, 
                                                   paste0(output_path, "eQTL/"))
cna_files_nonnormalBetas <- get_nonnormal_betaDistr(results_files_cna_eQTL, 
                                                   paste0(output_path, "eQTL/"))
