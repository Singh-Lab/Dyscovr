############################################################
### Process Linear Model Output
### Written By: Sara Geraghty, May 2021
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
# install.packages("meta", "metasens")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("pathview")
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
library("pathview")
library("gage")
library(dendextend)
library(gclus)
library(meta)
library(metasens)
library(metafor)

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
  hist(betas, main = "Histogram of Beta Coefficient Values",
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
  hist(pvals, main = "Histogram of p-Values",
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
  qobj <- NA
  if(length(results_table$p.value) < 100) {
    # With a small number of pvalues we may not be able to accurately estimate pi0,
    # so we set to 1 (the equivalent of B-H correction)
    qobj <- qvalue(p = results_table$p.value, pi0 = 1)
  } else {
    qobj <- qvalue(p = results_table$p.value)
  }
  
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
create_heat_map <- function(results_table) {
  
  matrix <- create_hm_input_matrix(results_table)
  
  # Plot a default heatmap
  #col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
  #heatmap(matrix, scale = "none", col = col)
  
  # Plot an enhanced heatmap
  # Clusters by default using hclust, but can specify others using param 'hclustfun'
  hm <- heatmap.2(matrix, scale = "none", col = bluered(100), trace = "none",
                  density.info = "none", dendrogram = "row", Colv = FALSE, 
                  Rowv = TRUE, key = TRUE, key.title = NA, key.xlab = "Beta") 
  plot(hm)
  #heatmap.2(matrix, scale = "none", col = bluered(100), trace = "none",
  #density.info = "none", labCol = "", dendrogram = c("row"), 
  #add.expr = text(x = seq_along(colnames(matrix)), 
  #y = -2, srt = 0, labels = colnames(matrix), 
  #xpd = NA, cex = 2, pos = 1))
  
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
  
  # Opt: Get and return the gene order after clustering
  # return(data.frame(gene = rownames(matrix)[hm$rowInd]))

  return(matrix)
}


#' Helper function to create the heat map input matrix
#' @param results_table a master DF produced from linear_model.R
create_hm_input_matrix <- function(results_table) {
  # Create the regulatory protein vs. gene target table
  matrix <- data.frame(matrix(ncol = length(unique(results_table$R_i.name)),
                              nrow = length(unique(results_table$T_k.name))))
  colnames(matrix) <- unique(results_table$R_i.name)
  rownames(matrix) <- unique(results_table$T_k.name)
  
  # Fill in this table with t-statistics/ Betas
  for (i in 1:nrow(results_table)) {
    #tstat <- results_table$statistic[i]
    beta <- results_table$estimate[i]
    regprot <- results_table$R_i.name[i]
    targ <- results_table$T_k.name[i]
    
    tryCatch({
      #matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- tstat
      matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- beta
    }, error=function(cond){
      print(cond)
      matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- 0
    })
  }
  # NOTE: leave all the unfilled pairings as NA
  
  # Replace NA/NaN with 0, or alternatively use na.omit
  #matrix[is.na(matrix)] <- 0
  #matrix[is.nan(matrix)] <- 0
  matrix <- na.omit(matrix)
  
  # Convert from data frame to matrix
  matrix <- data.matrix(matrix)
  print(head(matrix))
  
  return(matrix)
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


# Additional playing with clustering
matrix <- create_heat_map(master_df_mut_corrected, main_path)
hm <- heatmap.2(matrix, scale = "none", col = bluered(100), trace = "none",
                density.info = "none")

row_clust <- hclust(dist(matrix, method = "euclidean"), method = 'ward.D2')
plot(row_clust)
h2_sort <- sort(cutree(row_clust, h = 2))
plot(row_clust)
abline(h = 2, col = "red2", lty = 2, lwd = 2)

h2_sort_df <- as.data.frame(h22_sort)
colnames(h2_sort_df)[1] <- "cluster"
h2_sort_df$targ_gene <- rownames(h2_sort_df)

dendrogram <- as.dendrogram(row_clust)
cols_branches <- c("darkred", "forestgreen", "orange", "firebrick1", "yellow", 
                   "deeppink1", "cyan2", "darkslategray", "chartreuse")
# Set the colors of 9 branches
dendrogram <- color_branches(dendrogram, k = 9, col = cols_branches)
col_labels <- get_leaves_branches_col(dendrogram)
#col_labels <- list("1" = "Serine glycine biosynthesis", "2" ="GABA-B receptor II signaling", 
                  # "3" = "Fructose galactose metabolism/ Glycolysis", "4" = "N/A", 
                  # "5" = "S-adenosylmethionine biosynthesis", "6" = "Sulfate assimilation", 
                  # "7" = "Fructose galactose metabolism/ Androgen/etrogene/progesterone biosynthesis", 
                  # "8" = "2-arachidonoylglycerol biosynthesis", "9" = "O-antigen biosynthesis")
col_labels <- col_labels[order(order.dendrogram(dendrogram))]


h2_sort_df$labels <- unlist(lapply(1:length(unique(h2_sort_df$cluster)), function(i) {
  num_entries <- nrow(h2_sort_df[h2_sort_df$cluster == i,])
  return(rep(col_labels[[i]], times = num_entries))
}))
h22_sort_df$colors <- unlist(lapply(1:length(unique(h2_sort_df$cluster)), function(i) {
  num_entries <- nrow(h2_sort_df[h2_sort_df$cluster == i,])
  return(rep(cols_branches[[i]], times = num_entries))
}))


heatmap.2(matrix, scale = "none", col = bluered(100), trace = "none", density.info = "none",
          Rowv = dendrogram, key.xlab = "Beta", RowSideColors = col_labels,
          colRow = col_labels)



# Cluster when sorting the x-axis by order from another heat map (e.g. all BRCA)
matrix <- create_heat_map(master_df_mut_corrected, main_path)
hm <- heatmap.2(matrix, scale = "none", col = bluered(100), trace = "none",
                density.info = "none")
gene_order <- data.frame(gene = rownames(matrix)[hm$rowInd])

#Opt: write this to a file
gene_order_r <- t(gene_order)
gene_order_r <- rev(gene_order_r)
write.csv(gene_order_r, paste0(output_path, "gene_order_allBRCA_recon3d_mut.csv"))

# Read back, if needed
gene_order_r <- read.csv(paste0(output_path, "gene_order_allBRCA_recon3d_mut.csv"))[,2]

# Replace "master_df_mut_corrected_subtype" with the subtype DF we want to reorder
master_df_mut_corrected_subtype <- left_join(data.frame(T_k.name = gene_order_r), 
                                             master_df_mut_corrected_subtype, by = "T_k.name")
create_heat_map(master_df_mut_corrected_subtype)


# Cluster when sorting the x-axis by known pathway genes (cluster within each known group)
clusters <- list(
  "Serine.synthesis" = c("PHGDH", "PSAT1", "PSPH"),
  "Folate.cycle" = c("SHMT1", "SHMT2", "MTHFD1", "MTHFD1L", "MTHFD2", "MTHFD2L", 
                     "MTHFR", "TYMS", "DHFR", "FTCD", "MTFMT", "ALDH1L1"),
  "Methionine.SAM.SAH" = c("MAT1A", "MAT2A", "MAT2B", "GNMT", "AHCY", "MTR"),
  "Transsulferation" = c("CBS", "CTH"),
  "Taurine.synthesis" = c("CDO1", "CSAD"),
  "Glutathione.prod" = c("GCLC", "GSS", "GCLM", "GPX1", "GSTP1"),
  "Glycine.cleavage" = c("GCSH", "GLDC", "AMT", "DLD", "GLYAT"),
  "Choline" = c("CHDH", "ALDH7A1", "BHMT", "BHMT2", "DMGDH", "SARDH"),
  "Lipid.metabolism" = c("CHKA", "CHKB", "PCYT1A", "PCYT1B", "CHPT1", "PEMT"),
  "Glycolysis" = c("HK1", "HK2", "PGI", "PFK1", "PFK2", "ALDOA", "TPI", "GAPDH", 
                   "PGK1", "PGK2", "PGAM1", "PGAM2", "ENO1", "PKM1", "PKM2"),
  "TCA.cycle" = c("CS", "ACO1", "ACO2", "IDH1", "IDH2", "OGDH", "SUCLA2", "SDHA", 
                  "SDHB", "SDHC", "FH", "MDH1", "MDH2", "PC"),
  "Nucleot.Biosyn" = c("PRPS", "PPAT", "GART", "PFAS", "PAICS", "ADSL", "ATIC", 
                       "IMPDH1", "GMPS", "ADSS", "ADSL", "AK1", "NME"),
  "Pentose.phos" = c("ALDOA", "ALDOB", "ALDOC", "DERA", "FBP1", "FBP2", "G6PD", 
                     "GPI", "H6PD", "PFKL", "PFKM", "PFKP", "PGD", "PGLS", "PGM1", 
                     "PGM3", "PRPS1", "PRPS1L1", "PRPS2", "RBKS")
)
clusters$Other <- setdiff(res_recon3d_mut$T_k.name, as.character(unlist(clusters)))


#' Create the regulatory protein vs. gene target ordered matrix for clustering,
#' using functions from the gclus package
#' @param results_table master DF results from linear model
#' @param clusters a list with each known cluster and the corresponding gene names
create_clust_matrix <- function(results_table, clusters) {
  matrix <- create_hm_input_matrix(results_table)
    
  #dist_vect <- matrix[,1]
  #names(dist_vect) <- rownames(matrix)
  
  # Create a symmetrical distance matrix from this vector
  #dist_mat <- dist(dist_vect)
  
  # Get numeric vector corresponding to which cluster each gene is in
  #print(dist_vect)
  clusters_numeric <- as.numeric(unlist(lapply(rownames(matrix), function(x) {
    cluster_num <- which(sapply(clusters, function(y) x %in% y))
    return(cluster_num)
  })))
  
  # Create an ordered matrix
  #ordered_mat1 <- order.clusters(-matrix[,1], clusters_numeric)
  
  # For each given cluster, for each regprot, do kmeans clustering within this cluster
  ri_dfs <- lapply(unique(results_table$R_i.name), function(regprot) {
    cluster_dfs <- lapply(1:length(clusters), function(i) {
      # Get the genes from this cluster
      genes <- clusters[[i]]
      tab_genes <- res_recon3d_mut[which(res_recon3d_mut$T_k.name %fin% genes) & 
                                     which(res_recon3d_mut$R_i.name == regprot),]
      tab_genes <- tab_genes[, c("estimate", "T_k.name", "R_i.name")]
      
      # Cluster this to get "inner clusters" within this cluster
      if(nrow(tab_genes) < 2) {return(tab_genes)}
      tab_clust <- hclust(dist(tab_genes, method = "euclidean"), method = 'ward.D2')
      tab_genes_sort <- as.data.frame(sort(cutree(tab_clust, h = 2)))
      colnames(tab_genes_sort)[1] <- "inner.cluster"
      tab_genes_sort$targ_gene <- rownames(tab_genes_sort)
      
      # Reorder the genes based on these inner clusters

      
      # Return this clustered DF
      print(tab_genes_sort)
      return(tab_genes_sort)
    })
    # Recombine these into one
    print(cluster_dfs)
    combined_cluster_df <- do.call(rbind, cluster_dfs)
    print(combined_cluster_df)
    return(combined_cluster_df)
  })
  
  # Recombine this by column
  ordered_df <- do.call(cbind, ri_dfs)
  
  # Return this ordered matrix
  return(ordered_mat)
}


tp53_dist_mat <- create_dist_matrix(res_recon3d_mut, "TP53", clusters)
pik3ca_dist_mat <- create_dist_matrix(res_recon3d_mut, "PIK3CA", clusters)


#' Plot the ordered matrix as a heat map, using k-means clustering within each known cluster
#' @param dist_mat1 distance matrix produced from the prior function for goi 1
#' @param dist_mat2 distance matrix produced from the prior function for goi 2
plot_grouped_heatmap <- function(grouped_dist_mat, grouped_ordered_mat) {
  # Get the k-means clustering
  k_means <- kmeans(grouped_dist_mat, length(clusters))$cluster
  
  # Get order obtained from the k-means clustering
  #grouped_ordered_mat <- unlist(memship2clus(k_means))
  
  # Plot
  #layout(matrix(1:3,nrow=1,ncol=3),widths=c(0.1,1,1))
  #op = par(mar = 1,1,1,0.2)
  #colors <- cbind(k_means, k_means)
  #plotcolors(colors[grouped_ordered_mat,])
  
  #par(mar=c(1,6,1,1))
  #rlabels <- names(grouped_dist_mat)[grouped_ordered_mat]
  
  #plotcolors(cmat[grouped_ordered_mat, grouped_ordered_mat], rlabels = rlabels)
  
  #par(op)
  #layout(matrix(1,2))
  
  # Alternative
  grouped_ordered_mat <- order.clusters(-grouped_dist_mat, k_means)
  
  par(mar = c(1,1,1,1))
  layout(matrix(1:3,nrow=1,ncol=3),widths=c(0.1,1,1))

  colors <- cbind(k_means, k_means)
  plotcolors(colors[grouped_ordered_mat,])
  
  par(mar=c(1,6,1,3))
  rlabels <- names(grouped_dist_mat)[grouped_ordered_mat]
  cmat <- dmat.color(as.matrix(grouped_dist_mat), rev(cm.colors(length(clusters))))
  plotcolors(cmat[grouped_ordered_mat, grouped_ordered_mat], rlabels = rlabels)
  
  par(op)
  layout(matrix(1,1))
}


############################################################
############################################################
#### HEAT MAP VISUALIZATION FOR TOP MUTATED GENES
############################################################
############################################################

#' Given a master DF for a particular set of patients (e.g. full set or subtype)
#' that was run on multiple candidate or query genes on a set of target genes,
#' create an output figure that visualizes the Betas where query genes are labeled
#' on the left and names of target genes (on x-axis) are ignored 
#' @param results_table results_table a master DF produced from linear_model.R
#' @param outpath a path to a directory to write the t-statistic matrix to
create_hm_for_top_mutated <- function(results_table) {
  
  # Hack- limit to just the top 4 mutated
  results_table <- results_table[results_table$R_i.name %in% c("TP53", "PIK3CA", "TTN", "RYR2", "KMT2C"),]
  
  # Create the regulatory protein vs. gene target table
  matrix <- t(create_hm_input_matrix(results_table))
  matrix <- data.matrix(matrix)
  #print(matrix[1,])
  
  # Plot an enhanced heatmap
  # Clusters by default using hclust, but can specify others using param 'hclustfun'
  heatmap.2(matrix, scale = "none", col = bluered(100), trace = "none",
                  density.info = "none", Colv = TRUE, dendrogram = "none",
                  Rowv = FALSE, key = TRUE, key.title = NA, key.xlab = "Beta", 
                  rowsep = 1:nrow(matrix), sepwidth = c(0.1, 0.1), 
                  labCol = rep("", ncol(matrix)), offsetRow = -37)
  
  return(matrix)
}


# Call this function
create_hm_for_top_mutated(master_df_mut_corrected)
create_hm_for_top_mutated(master_df_cna_corrected)

# Get the absolute value of the sum of Betas for given R_i
print(rowSums(apply(master_df_mut_corrected, 2, abs)))

#' Perform q-value correction on each R_i individually 
#' @param master_df
get_qval_dict <- function(master_df) {
  
  master_df <- master_df[!(is.na(master_df$p.value)),]
  
  qval_dict <- lapply(1:length(unique(master_df$R_i.name)), function(i) {
    r_i <- unique(master_df$R_i.name)[i]
    #print(r_i)
    pvals <- master_df[master_df$R_i.name == r_i, 'p.value']
    #print(pvals)
    qobj <- NA
    if(length(pvals) < 100) {
      # With a small number of pvalues we may not be able to accurately estimate pi0,
      # so we set to 1 (the equivalent of B-H correction)
      qobj <- qvalue(p = pvals, pi0 = 1)
    } else {
      qobj <- qvalue(p = pvals)
    }
    
    # OPT: plot some useful plots & print some useful information
    #plot(qobj)
    #print(summary(qobj))
    #print(paste("Pi0 (Propr. of true null hypotheses):", qobj$pi0))
    
    qvals <- qobj$qvalues
    return(list("Names" = master_df[master_df$R_i.name == r_i, 'T_k.name'],
                "Qvals" = qvals))
  })
  names(qval_dict) <- unique(master_df$R_i.name)
  return(qval_dict)
}


#' Use this qval to print the number of significant hits for each R_i
#' @param qval_dict a q-value dictionary as output above (Names and Qvals
#' for each query protein)
print_num_sig_hits <- function(qval_dict, thres) {
  for(i in 1:length(qval_dict)) {
    print(paste("Protein name:", names(qval_dict)[i]))
    print(paste("Number of significant hits, q <", thres))
    curr_entry <- qval_dict[[i]]
    print(length(curr_entry$Names[which(curr_entry$Qvals < thres)]))
  }
}


############################################################
############################################################
#### GROUPED BAR CHART TO SHOW NUMBERS OF SIGNIFICANT HITS
#### FOR EACH OF QUERY GENES
############################################################
############################################################
#' Creates a stacked bar plot to show the number of significant hits
#' (at given q-value threshold) found for each query gene, for 
#' each particular subtype of interest
#' @param qval_dict_list a list of qvalue dictionaries (as output above),
#' one for each subtype of interest. Each entry should have the same 
#' query genes.
#' @param thres a qvalue threshold for significance
create_stacked_sig_hits_bar_plot <- function(qval_dict_list, thres) {
  sub_dfs <- lapply(1:length(qval_dict_list), function(i) {
    subtype_entry <- qval_dict_list[[i]]
    prot_dfs <- lapply(1:length(subtype_entry), function(j) {
      prot_entry <- subtype_entry[[j]]
      # Turn this list into a DF
      prot_df <- data.frame("Query.Prot" = names(subtype_entry)[j],
                            "Num.Sig.Hits" = length(prot_entry$Names[which(prot_entry$Qvals < thres)]))
      #prot_df$Query.Prot <- rep(names(prot_entry)[j], times = ncol(prot_df))
      return(prot_df)
    })
    combined_prot_df <- do.call(rbind, prot_dfs)
    print(head(combined_prot_df))
    # Add the name of the subtype
    print(names(qval_dict_list)[i])
    combined_prot_df$Subtype <- rep(names(qval_dict_list)[i], 
                                    times = nrow(combined_prot_df))
    return(combined_prot_df)
  }) 
  # Combined these into one 
  full_df <- do.call(rbind, sub_dfs)
  #print(head(full_df))
  
  # Create a stacked bar plot from this
  ggplot(full_df, aes(fill=Query.Prot, y=Num.Sig.Hits, x=Subtype)) + 
    geom_bar(position="stack", stat="identity")
  
  return(full_df)
}

# Create a list of qval dictionaries for each subtype
qval_dict_lumA <- get_qval_dict(master_df_lumA)
qval_dict_lumB <- get_qval_dict(master_df_lumB)
qval_dict_basal <- get_qval_dict(master_df_basal)
qval_dict_her2 <- get_qval_dict(master_df_her2)

qval_dict_list <- list("LumA" = qval_dict_lumA, "LumB" = qval_dict_lumB,
                       "Basal" = qval_dict_basal, "HER2" = qval_dict_her2)

combined_df <- create_stacked_sig_hits_bar_plot(qval_dict_list, 0.1)

driver_genes <- c("TP53", "PIK3CA", "KMT2C")
combined_df$Status <- unlist(lapply(combined_df$Query.Prot, function(g) 
  ifelse(g %in% driver_genes, "Driver", "Non-Driver")))

nb = length(unique(combined_df$Query.Prot))
nm = length(unique(combined_df$Status))
colors = apply(expand.grid(seq(70,40,length=nm), 100, seq(15,375,length=nb+1)[1:nb]), 1, 
               function(x) hcl(x[3],x[2],x[1]))
colors[3] <- "#E8909C"

ggplot(combined_df, aes(fill=interaction(Query.Prot, Status), y=Num.Sig.Hits, x=Subtype)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=colors) + theme_classic()


# Examine significant hit overlap
p53_top_allgenes <- list("LumA" = qval_dict_lumA_allgenes$TP53$Names[which(qval_dict_lumA_allgenes$TP53$Qvals < 0.2)],
                         "LumB" = qval_dict_lumB_allgenes$TP53$Names[which(qval_dict_lumB_allgenes$TP53$Qvals < 0.2)],
                         "Basal" = qval_dict_basal_allgenes$TP53$Names[which(qval_dict_basal_allgenes$TP53$Qvals < 0.2)],
                         "HER2" = qval_dict_her2_allgenes$TP53$Names[which(qval_dict_her2_allgenes$TP53$Qvals < 0.2)])
plt <- venn.diagram(p53_top_allgenes, category.names = c("LumA", "LumB", "Basal", "HER2"), 
                    filename = NULL, output = TRUE, lwd = 2, lty = 'blank', 
                    fill = c("red", "blue", "green", "yellow"), cex = 2, fontface = "bold",
                    fontfamily = "sans", cat.cex = 2, cat.fontface = "bold",cat.fontfamily = "sans")
                    #cat.default.pos = "outer", cat.pos = c(180, 180, 180, 180)) #, cat.fontfamily = "sans", rotation = 1)
grid::grid.draw(plt)


############################################################
############################################################
#### SPEARMAN CORRELATION OF BETAS AND T-STATISTICS
############################################################
############################################################
#' Compute and print the Spearman correlation of the Betas and of the T-statistic,
#' given two groups of interest
#' @param results_table the output master DF from the linear model
#' @param ri_1 the external gene name of the first protein of interest
#' @param ri_2 the external gene name of the second protein of interest
compute_and_print_spearman <- function(results_table, ri_1, ri_2) {
  if((ri_1 %fin% results_table$R_i.name) & (ri_2 %fin% results_table$R_i.name)) {
    
    target_genes <- unique(results_table$T_k.name)
    
    # Mini functions to get Betas/t-statistics for each target gene
    #' @param results_table the master DF 
    #' @param target_genes the list of target genes
    #' @param ri the given regulatory protein
    #' @param type "Betas" or "t-statistics" to indicate what value we are returning
    get_values <- function(results_table, target_genes, ri, type) {
      vals <- unlist(lapply(target_genes, function(tg) {
        if(type == "Betas") {
          est <- results_table[(results_table$R_i.name == ri) & (results_table$T_k.name == tg), "estimate"]
        } else {
          est <- results_table[(results_table$R_i.name == ri) & (results_table$T_k.name == tg), "statistic"]
        }
        if (length(est) == 0) {est <- 0}  # To ensure the lengths of the Beta vectors are the same
        return(est)
      }))
    }
    
    grp1_Betas <- get_values(results_table, target_genes, ri_1, "Betas")
    grp2_Betas <- get_values(results_table, target_genes, ri_2, "Betas")
    
    # Get Betas spearman
    betas_spearman <- cor.test(grp1_Betas, grp2_Betas, method = "spearman")
    betas_spearman_stat <- as.numeric(betas_spearman$estimate)
    betas_spearman_pval <- betas_spearman$p.value
    
    # Print the results
    print(paste("Spearman results for", paste(ri_1, paste("and", ri_2))))
    print(paste("Beta correlation of", paste(betas_spearman_stat, paste(", p-value of", betas_spearman_pval))))
    
    # Create a plot to visualize the correlations
    plot(grp1_Betas, grp2_Betas, pch = 19, col = "lightblue", #main = "Betas Spearman Correlation",
         xlab = paste(ri_1, "Betas"), ylab = paste(ri_2, "Betas"))
    abline(lm(grp2_Betas ~ grp1_Betas), col = "red", lwd = 3)
    text(labels = paste("Correlation:", paste(round(betas_spearman_stat, 6), 
                                              paste(", p-value:", round(betas_spearman_pval, 6)))), 
         x = max(grp1_Betas, na.rm = TRUE)-sd(grp1_Betas, na.rm = TRUE)*3, 
         y = max(grp2_Betas, na.rm = TRUE)-sd(grp2_Betas, na.rm = TRUE), col = "black")
    
    
  } else {print("Error. Provided regprots are not in the given master DF.")}
}

ri_1 <- "TP53"
ri_2 <- "PIK3CA"

# Call function
compute_and_print_spearman(master_df, ri_1, ri_2)



#' Compute and print Spearman, given two separate results tables
#' @param results_table1 first output master DF from the linear model
#' @param results_table2 second output master DF from the linear model
#' @param ri_1 the external gene name of the first protein of interest
#' @param ri_2 the external gene name of the second protein of interest
compute_and_print_spearman_multDF <- function(results_table1, results_table2, ri_1, ri_2) {
  if((ri_1 %fin% results_table1$R_i.name) & (ri_2 %fin% results_table2$R_i.name)) {
    
    target_genes <- intersect(unique(results_table1$T_k.name), unique(results_table2$T_k.name))
    
    # Mini functions to get Betas/t-statistics for each target gene
    #' @param results_table the master DF 
    #' @param target_genes the list of target genes
    #' @param ri the given regulatory protein
    #' @param type "Betas" or "t-statistics" to indicate what value we are returning
    get_values <- function(results_table, target_genes, ri, type) {
      vals <- unlist(lapply(target_genes, function(tg) {
        if(type == "Betas") {
          est <- results_table[(results_table$R_i.name == ri) & (results_table$T_k.name == tg), "estimate"]
        } else {
          est <- results_table[(results_table$R_i.name == ri) & (results_table$T_k.name == tg), "statistic"]
        }
        if (length(est) == 0) {est <- 0}  # To ensure the lengths of the Beta vectors are the same
        return(est)
      }))
    }
    
    grp1_Betas <- get_values(results_table1, target_genes, ri_1, "Betas")
    grp2_Betas <- get_values(results_table2, target_genes, ri_2, "Betas")
    
    # Get Betas spearman
    betas_spearman <- cor.test(grp1_Betas, grp2_Betas, method = "spearman")
    betas_spearman_stat <- as.numeric(betas_spearman$estimate)
    betas_spearman_pval <- betas_spearman$p.value
    
    # Print the results
    print(paste("Spearman results for", paste(ri_1, paste("and", ri_2))))
    print(paste("Beta correlation of", paste(betas_spearman_stat, paste(", p-value of", betas_spearman_pval))))
    
    # Create a plot to visualize the correlations
    plot(grp1_Betas, grp2_Betas, pch = 19, col = "lightblue", #main = "Betas Spearman Correlation",
         xlab = paste(ri_1, "Betas"), ylab = paste(ri_2, "Betas"))
    abline(lm(grp2_Betas ~ grp1_Betas), col = "red", lwd = 3)
    text(labels = paste("Correlation:", paste(round(betas_spearman_stat, 6), 
                                              paste(", p-value:", round(betas_spearman_pval, 6)))), 
         x = max(grp1_Betas, na.rm = TRUE)-sd(grp1_Betas, na.rm = TRUE)*3, 
         y = max(grp2_Betas, na.rm = TRUE)-sd(grp2_Betas, na.rm = TRUE), col = "black")
    
    
  } else {print("Error. Provided regprots are not in the given master DF.")}
}

ri_1 <- "TP53"
ri_2 <- "PIK3CA"

# Call function
compute_and_print_spearman_multDF(master_df1, master_df2, ri_1, ri_2)


#' Compute and print the Spearman correlation of the T-statistics for two separate
#' subpopulations
#' @param results_table_subpop1 the output master DF from the linear model for 
#' the first subgroup
#' #' @param results_table_subpop2 the output master DF from the linear model for 
#' the second subgroup
#' @param ri the external gene name of the protein of interest
#' @param subpop1 the name of the first subgroup
#' @param subpop2 the name of the second subgroup
#' @param sub_to_top if given, subsets to the genes with smallest X and largest X
#' t-statistics from the two subpopulations
compute_and_print_spearman_subpops <- function(results_table_subpop1, results_table_subpop2, 
                                               ri, subpop1, subpop2, sub_to_top) {
  if((ri %fin% results_table_subpop1$R_i.name) & (ri %fin% results_table_subpop2$R_i.name)) {
    
    if(is.infinite(sub_to_top)) {
      target_genes <- intersect(unique(results_table_subpop1$T_k.name), 
                                unique(results_table_subpop2$T_k.name))
    } else {
      # Use a helper function to subset the target genes
      target_genes <- subset_targ_genes(results_table_subpop1, results_table_subpop2, 
                                        ri, sub_to_top)
    }

    # Mini functions to get t-statistics for each target gene
    #' @param results_table the master DF 
    #' @param target_genes the list of target genes
    #' @param ri the given regulatory protein
    get_values <- function(results_table, target_genes, ri) {
      vals <- unlist(lapply(target_genes, function(tg) {
        est <- results_table[(results_table$R_i.name == ri) & (results_table$T_k.name == tg), "statistic"]
        if (length(est) == 0) {est <- 0}  # To ensure the lengths of the Beta vectors are the same
        return(est)
      }))
    }
    
    subpop1_tstats <- get_values(results_table_subpop1, target_genes, ri)
    subpop2_tstats <- get_values(results_table_subpop2, target_genes, ri)
    
    # Get Spearman
    spearman <- cor.test(subpop1_tstats, subpop2_tstats, method = "spearman")
    spearman_stat <- as.numeric(spearman$estimate)
    spearman_pval <- spearman$p.value
    
    # Print the results
    print(paste("Spearman results for", paste(subpop1, paste("and", paste(subpop2, paste("in", ri))))))
    print(paste("T-statistic correlation of", paste(spearman_stat, paste(", p-value of", spearman_pval))))
    
    # Create a plot to visualize the correlations
    plot(subpop1_tstats, subpop2_tstats, pch = 19, col = "lightblue", #main = "Betas Spearman Correlation",
         xlab = paste(subpop1, paste(ri, "T-statistics")), ylab = paste(subpop2, paste(ri, "T-statistics")))
    abline(lm(subpop1_tstats ~ subpop2_tstats), col = "red", lwd = 3)
    text(labels = paste("Correlation:", paste(round(spearman_stat, 6), 
                                              paste(", p-value:", round(spearman_pval, 6)))), 
         x = max(subpop1_tstats, na.rm = TRUE)-sd(subpop1_tstats, na.rm = TRUE)*3, 
         y = max(subpop2_tstats, na.rm = TRUE)-sd(subpop2_tstats, na.rm = TRUE), col = "black")
    
    
  } else {print("Error. Provided regprot not in the given master DF.")}
}

#' Helper function for subsetting the target genes to the top X/ bottom X t-statistics 
#' in the two subpopulations
#' @param results_table_subpop1 the output master DF from the linear model for 
#' the first subgroup
#' @param results_table_subpop2 the output master DF from the linear model for 
#' the second subgroup
#' @param ri the external gene name of the protein of interest
#' @param sub_to_top if given, subsets to the genes with smallest X and largest X
#' t-statistics from the two subpopulations
subset_targ_genes <- function(results_table_subpop1, results_table_subpop2, ri, sub_to_pop) {
  
  # First, subset to the regprot of interest
  results_table_subpop1_sub <- results_table_subpop1[results_table_subpop1$R_i.name == ri,]
  results_table_subpop2_sub <- results_table_subpop2[results_table_subpop2$R_i.name == ri,]
  
  # Next, get the top X and bottom X t-statistics 
  results_table_subpop1_sub <- arrange(results_table_subpop1_sub, desc(estimate))
  results_table_subpop2_sub <- arrange(results_table_subpop2_sub, desc(estimate))
  
  top100_subpop1 <- results_table_subpop1_sub[1:sub_to_pop, 'T_k.name']
  bottom100_subpop1 <- results_table_subpop1_sub[(nrow(results_table_subpop1_sub)-sub_to_pop):nrow(results_table_subpop1_sub), 
                                                 'T_k.name']
  top100_subpop2 <- results_table_subpop2_sub[1:sub_to_pop, 'T_k.name']
  bottom100_subpop2 <- results_table_subpop2_sub[(nrow(results_table_subpop2_sub)-sub_to_pop):nrow(results_table_subpop2_sub), 
                                                 'T_k.name']
  
  top_genes <- unique(c(top100_subpop1, bottom100_subpop1, top100_subpop2, bottom100_subpop2))
  
  return(top_genes)
}

ri <- "TP53"
ri <- "PIK3CA"

# Call function
compute_and_print_spearman_subpops(lumA_allgenes_mut, lumB_allgenes_mut, ri, 
                                   "LumA", "LumB", 100)


############################################################
############################################################
#### VISUALIZE ON RELEVANT PATHWAYS
############################################################
############################################################
#' Visualize our results on pathways using the Pathview package
#' Link to vignette: https://pathview.r-forge.r-project.org/pathview.pdf
#' 
visualize_pathway <- function(master_df_corrected, pathway_id, label) {
  gene_data <- master_df_corrected$estimate
  names(gene_data) <- master_df_corrected$T_k.name
  
  # Run using pathview
  pv.out <- pathview(gene.data = gene_data, pathway.id = pathway_id,
                     species = "hsa", out.suffix = label, kegg.native = T,
                     gene.idtype = "SYMBOL")
}

pathway_tp53 <- "04115"
pathway_breastcancer <- "05224"
pathway_cancerCarbonMetabolism <- "05230"
pathway_cancerCholineMetabolism <- "05231"

visualize_pathway(master_df_corrected, pathway_tp53, "cancerRelated_metabolicTargs")
visualize_pathway(master_df_corrected, pathway_breastcancer, "cancerRelated_metabolicTargs")
visualize_pathway(master_df_corrected, pathway_cancerCarbonMetabolism, "cancerRelated_metabolicTargs")
visualize_pathway(master_df_corrected, pathway_cancerCholineMetabolism, "cancerRelated_metabolicTargs")



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
#' @param master_df_mut_sig mutation master DF, with gene names added and 
#' thresholded for significance
#'@param master_df_cna_sig CNA master DF, with gene names added and 
#' thresholded for significance
plot_tophit_overlap <- function(master_df_mut_sig, master_df_cna_sig) {
  #overlap_genes <- intersect(master_df_mut_sig$T_k.name, master_df_cna_sig$T_k.name)
  #print(overlap_genes)
  
  # Alternatively, look at only combinations of overlap
  pairs_mut_sig <- unlist(lapply(1:nrow(master_df_mut_sig), function(i)
    paste(master_df_mut_sig[i, 'R_i.name'], master_df_mut_sig[i, 'T_k.name'], sep = ":")))
  pairs_cna_sig <- unlist(lapply(1:nrow(master_df_cna_sig), function(i)
    paste(master_df_cna_sig[i, 'R_i.name'], master_df_cna_sig[i, 'T_k.name'], sep = ":")))
  print(paste("Length intersect:", length(intersect(pairs_mut_sig, pairs_cna_sig))))
  
  pairs_list <- list("Mutation" = pairs_mut_sig, "CNA" = pairs_cna_sig)
  
  # Plot a Venn Diagram
  #myCol <- brewer.pal(2, "Pastel2")
  plt <- venn.diagram(pairs_list, category.names = c("Mutation", "CNA"), filename = NULL, output = TRUE,
               lwd = 2, lty = 'blank', fill = c("red", "blue"), cex = 2, fontface = "bold",
               fontfamily = "sans", cat.cex = 2, cat.fontface = "bold",cat.fontfamily = "sans",
               cat.default.pos = "outer", cat.pos = c(180, 180)) #, cat.fontfamily = "sans", rotation = 1)
  grid::grid.draw(plt)
}

# Call this function and write to PNG file
fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/Top Gene Hit Overlap (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 450, height = 350)
plot_tophit_overlap(master_df_mut_sig, master_df_cna_sig)
dev.off()


############################################################
############################################################
#### VISUALIZE THE RELATIVE RANKS OF TOP HITS FOR TWO REGPROTS
############################################################
############################################################
#' Creates a line graph with two lines - one for each regulatory protein
#' of interest. Plots the rank (of the top hit from the master DF) vs.
#' the cumulative sum to that rank of hits that were for the given regprot.
#' @param master_df a DF produced as output from the linear model, MHT corrected
#' @param regprot1 the gene name of the first regulatory protein of interest
#' @param regprot2 the gene name of the second regulatory protein of interest
visualize_tophit_relative_ranks <- function(master_df, regprot1, regprot2) {
  
  regprot1_hits <- unlist(lapply(1:nrow(master_df), function(i) 
    ifelse(master_df$R_i.name[i] == regprot1, 1, 0)))
  regprot2_hits <- unlist(lapply(1:nrow(master_df), function(i) 
    ifelse(master_df$R_i.name[i] == regprot2, 1, 0)))
  
  regprot1_cumul <- unlist(lapply(1:length(regprot1_hits), function(i) {
    return(sum(as.numeric(regprot1_hits[1:i])))
  }))
  regprot2_cumul <- unlist(lapply(1:length(regprot2_hits), function(i) {
    return(sum(as.numeric(regprot2_hits[1:i])))
  }))

  df_tmp <- data.frame(regprot1 = regprot1_cumul, regprot2 = regprot2_cumul)
  df_tmp$rank <- 1:nrow(df_tmp)
  
  plot(df_tmp$rank, df_tmp[,1], col = "orange", lty = 1, xlab = "rank", 
       ylab = "cumulative # of top hits")
  points(df_tmp$rank, df_tmp[,2], col = "darkgreen", lty = 1)
  
  legend(x = "bottomright", y=NULL, legend = c(regprot1, regprot2), fill = c("orange", "darkgreen"))
}

visualize_tophit_relative_ranks(master_df_mut_corrected, "TP53", "PIK3CA")


#' Creates a line graph with two lines - one for each regulatory protein
#' of interest. Plots the rank (of the top hit from the master DF) vs.
#' the Beta to that rank of hits that were for the given regprot.
#' @param master_df a DF produced as output from the linear model, MHT corrected
#' @param regprot1 the gene name of the first regulatory protein of interest
#' @param regprot2 the gene name of the second regulatory protein of interest
visualize_tophit_relative_Betas <- function(master_df, regprot1, regprot2) {
  
  regprot1_hits <- unlist(lapply(1:nrow(master_df), function(i) 
    ifelse(master_df$R_i.name[i] == regprot1, 1, 0)))
  regprot2_hits <- unlist(lapply(1:nrow(master_df), function(i) 
    ifelse(master_df$R_i.name[i] == regprot2, 1, 0)))
  
  get_beta_by_rank <- function(regprot_hits, master_df) {
    prev_Beta <- 0
    regprot_cumul <- c()
    for (i in 1:length(regprot_hits)) {
      if(regprot_hits[i] == 0) {
        regprot_cumul <- c(regprot_cumul, prev_Beta)
      } else {
        beta <- master_df$estimate[i]
        regprot_cumul <- c(regprot_cumul, beta)
        prev_Beta <- beta
      }
    }
    return(regprot_cumul)
  }
 
  regprot1_cumul <- get_beta_by_rank(regprot1_hits, master_df)
  regprot2_cumul <- get_beta_by_rank(regprot2_hits, master_df)
  
  df_tmp <- data.frame(regprot1 = regprot1_cumul, regprot2 = regprot2_cumul)
  df_tmp$rank <- 1:nrow(df_tmp)
  
  plot(df_tmp$rank, df_tmp[,1], col = "orange", lty = 1, xlab = "rank", 
       ylab = "Relative Betas")
  points(df_tmp$rank, df_tmp[,2], col = "darkgreen", lty = 1)
  abline(a = 0, b = 0)
  
  legend(x = "bottomright", y=NULL, legend = c(regprot1, regprot2), fill = c("orange", "darkgreen"))
}

visualize_tophit_relative_Betas(master_df_mut_corrected, "TP53", "PIK3CA")


############################################################
############################################################
#### PERFORM A META-ANALYSIS TO AGGREGATE SIGNAL
############################################################
############################################################
# Use the package "meta" for perform a meta-analysis, treating the various subtype
# results as results from separate studies
# Link to documentation: https://cran.r-project.org/web/packages/meta/meta.pdf

#' Create an input data frame of the type expected by meta
#' @param list_of_master_dfs a list of the master DFs with results we are interested
#' in combining in our meta-analysis. These are already subsetted to only have results
#' for one regulatory protein of interest. List names are the names of the subgroups.
#' @param vect_of_n a vector the same length as the above list, with the number of 
#' samples for each subgroup (n)
#' @param targ the external gene name of the target gene of interest
create_meta_input_tab <- function(list_of_master_dfs, vect_of_n, targ) {
  meta_input_tab <- data.frame("study" = names(list_of_master_dfs), 
                               "n_patients" = vect_of_n)
  # Get the Beta estimate and t-statistic value for each of the subgroups, for
  # the given target
  meta_input_tab$beta <- lapply(list_of_master_dfs, function(x) 
    x[x$T_k.name == targ, 'estimate']) 
  meta_input_tab$tstat <- lapply(list_of_master_dfs, function(x) 
    x[x$T_k.name == targ, 'statistic']) 
  meta_input_tab$std.error <- lapply(list_of_master_dfs, function(x) 
    x[x$T_k.name == targ, 'std.error']) 
  
  return(meta_input_tab)
}

master_df_list_p53 <- list("LumA"= lumA_allgenes_mut_p53, "LumB" = lumB_allgenes_mut_p53,
                       "Basal" = basal_allgenes_mut_p53, "HER2" = her2_allgenes_mut_p53)

patient_set <- read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt"), header = TRUE)[,1]
patient_set_lt20pamp <- intersect(read.table(paste0(main_path, "Patient Subsets/LessThan20PercAmp_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumA <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.A_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumB <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.B_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_basal <- intersect(read.table(paste0(main_path, "Patient Subsets/Basal_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_her2 <- intersect(read.table(paste0(main_path, "Patient Subsets/HER2_patient_ids.txt"), header = TRUE)[,1], patient_set)

n_of_subtypes <- c(length(patient_set_lumA), length(patient_set_lumB), length(patient_set_basal), length(patient_set_her2))
  
meta_input_tab_p53 <- create_meta_input_tab(master_df_list_p53, n_of_subtypes, "DDB2")
meta_input_tab_p53 <- create_meta_input_tab(master_df_list_p53, n_of_subtypes, "SRD5A1")

#' Perform the meta-analysis using a random-effects model
#' @param meta_input_tab an input table for meta-analysis as created by helper function
#' @param targ the external gene name of the target gene of interest
perform_meta_analysis <- function(meta_input_tab, targ) {
  ma_model <- rma(yi = as.numeric(meta_input_tab$tstat), 
                  sei = as.numeric(meta_input_tab$std.error))
  print(summary(ma_model))
  
  # Make a forest plot 
  forest(ma_model, slab = meta_input_tab$study)
  
  # Add details
  text(-10, -1.5, pos=4, cex=0.75, bquote(paste("RE Model (Q = ", 
                                              .(formatC(ma_model$QE, digits=2, format="f")), ", df = ", .(ma_model$k - ma_model$p),
                                              ", p = ", .(formatC(ma_model$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                              .(formatC(ma_model$I2, digits=1, format="f")), "%)")))
  # Add title
  grid.text(paste(targ, "Meta-Analysis"), .5, .8, gp=gpar(cex=1.5))
}

perform_meta_analysis(meta_input_tab_p53, "DDB2")
perform_meta_analysis(meta_input_tab_p53, "SRD5A1")

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

# Import TP53-specific ChIP-eat targets
tp53_chipeat_targets <- read.table(paste0(path, "Saved Output Data Files/BRCA/ChIP-eat/tp53_chipeat_targets.txt"), header = FALSE)[,1]

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
#### CREATE ENRICHMENT PLOTS
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


#' Function to plot the enrichment of cancer genes among the significant target genes
#' @param master_df a master DF produced from run_linear_model() that has q-values
#' @param chipeat_targs a vector of known ChIP-eat targets for a gene of interest (master
#' DF should already be subsetted to just this gene)
#' @param goi a string denoting the gene-of-interest, for labeling graph
plot_chipseq_enrichment <- function(master_df, chipeat_targs, goi) {
  # Get all the target genes 
  target_hits <- unlist(master_df$T_k.name)
  print(length(target_hits))
  print(length(chipeat_targs))
  
  # Fill in a vector of 0 and 1 for each significant target hit to indicate if it 
  # is a ChIP-eat target
  chipeat_vect <- unlist(lapply(target_hits, function(x) {
    return(ifelse(x %fin% chipeat_targs, 1, 0))
  }))
  
  # Plot the enrichment (fraction of ChIP-eat genes at/ above given rank)
  ranks <- 1:length(target_hits)
  frac_chipeat_vect <- c()
  count_of_chipeat_genes <- 0
  for (i in 1:length(target_hits)) {
    val_of_curr_tg <- chipeat_vect[i]
    count_of_chipeat_genes <- count_of_chipeat_genes + val_of_curr_tg
    frac <- count_of_chipeat_genes / i
    frac_chipeat_vect <- c(frac_chipeat_vect, frac)
  }
  
  plot(ranks[1:100], frac_chipeat_vect[1:100], pch = 16, 
       main = paste("Enrichment of Significant Target Genes in", paste(goi, "ChIP-eat Targets")),
       xlab = "Rank", ylab = "Fraction of Target Genes that are ChIP-eat Targets")
}


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
#### ASIDE: EXAMINE THE STATISTICAL OVERLAP OF MULTIPLE 
#### MASTER DF TOP HITS
############################################################
#' Looks at only top hits for each given master DF (only one regprot), below 
#' a given q-value threshold, and makes a Venn diagram to show overlap, along with 
#' reporting whether there is more overlap than might be expected given
#' the size of the total target pool
#' @param master_df1 the first master DF
#' @param master_df2 the second master DF
#' @param qval_thres the threshold below which we will consider the pairing significant
#' @param size_targgene_pool the number of target genes we are looking at
#' @param master_df1_label a label for what the first master DF is
#' @param master_df2_label a label for what the second master DF is
get_overlap_of_top_hits <- function(master_df1, master_df2, qval_thres, size_targgene_pool,
                                    master_df1_label, master_df2_label) {
  # Limit the master DFs to only significant hits
  master_df1_sig <- master_df1[master_df1$q.value < qval_thres,]
  master_df2_sig <- master_df1[master_df2$q.value < qval_thres,]
  
  # Plot a Venn diagram of the intersection
  myCol <- brewer.pal(3, "Accent")
  grid.newpage()
  v <- venn.diagram(x = list(master_df1_sig$T_k.name, master_df2_sig$T_k.name), 
               filename = NULL, #paste(master_df1_label, paste(master_df2_label, "BRCA.LumAB.png", sep = "_"), sep = "_"),
               #lwd = 2, lty = 'blank', col = c("#440154ff", '#21908dff'),  
               cex = 3, fill = myCol[1:2],
               fontface = "bold", fontfamily = "sans",
               cat.cex = 1.5, cat.fontface = "bold", #cat.default.pos = "outer", 
               cat.fontfamily = "sans", #rotation = 1,
               category.names = c(master_df1_label, master_df2_label),
               hyper.test = TRUE, lower.tail = FALSE)
  grid.draw(v)
  
  # Get the intersection between the top genes
  intersecting_genes <- intersect(master_df1_sig$T_k.name, master_df2_sig$T_k.name)
  num_intersecting_genes <- length(intersecting_genes)
  print(paste("Number of Significant Mutation Hits:", nrow(master_df1_sig)))
  print(paste("Number of Significant CNA Hits:", nrow(master_df2_sig)))
  print(paste("Length Overlap:", num_intersecting_genes))
  
  # Use the hypergeometric distribution to determine if the overlap is more
  # than we would expect by chance
  print(phyper(q = num_intersecting_genes - 1, m = nrow(master_df1_sig),
               n = size_targgene_pool - nrow(master_df1_sig), k = nrow(master_df2_sig),
               lower.tail = FALSE))
}

master_df_tp53_mut <- read.csv(paste0(main_path, "Linear Model/TP53/Tumor_Only/eQTL/output_results_LumAB_P53_metabolicTargs_iprotein_tmm_CNAbucket_justDel_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_uncorrected_MUT.csv"), header = TRUE, check.names = FALSE)
master_df_tp53_cna <- read.csv(paste0(main_path, "Linear Model/TP53/Tumor_Only/eQTL/output_results_LumAB_P53_metabolicTargs_iprotein_tmm_CNAbucket_justDel_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_uncorrected_CNA.csv"), header = TRUE, check.names = FALSE)

# Call function
# Num. metabolic targets in Recon3D: 1837
get_overlap_of_top_hits(master_df_tp53_mut, master_df_tp53_cna, 0.1, 1837, "TP53 Mutation", "TP53 Deletion")