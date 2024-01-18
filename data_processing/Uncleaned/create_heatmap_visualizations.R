############################################################
### CREATE HEAT MAP VISUALIZATIONS
### Written By: Sara Geraghty, July 2020
############################################################

# Visualize output of linear_model.R using a variety of heat map methods,
# including a basic heat map, a clustered and colored heat map, and a heat
# map of only the significant hit target genes

# install.packages("pheatmap")
# install.packages("gplots")
# BiocManager::install("gage")
# BiocManager::install("ComplexHeatmap")
library("gplots")
library("pheatmap")
library(ComplexHeatmap)
library(ggplot2)
library("gage")
library(dendextend)
library(gclus)

# NEJM color palatte: https://nanx.me/ggsci/reference/pal_nejm.html

# Path to output files
main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

# Path to where output figures should be saved
output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/BRCA/Linear Model/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/Pan-Cancer/Linear Model/"


# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


############################################################
#### BASIC HEAT MAP VISUALIZATION
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
  #plot(hm)
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


############################################################
#### CLUSTERED HEAT MAP WITH GENE TARGET CLUSTERS COLORED
############################################################
#' Create a clustered heatmap with gene target clusters colored
#' @param master_df_mut_corrected a master DF output from model
#' @param height_thres a height threshold for defining clusters
create_colored_heatmap <- function(master_df_mut_corrected, height_thres) {
  matrix <- create_heat_map(master_df_mut_corrected)
  hm <- heatmap.2(matrix, scale = "none", col = bluered(100), trace = "none",
                  density.info = "none")
  
  row_clust <- hclust(dist(matrix, method = "euclidean"), method = 'ward.D2')
  plot(row_clust)
  h_sort <- sort(cutree(row_clust, h = height_thres))
  plot(row_clust)
  abline(h = height_thres, col = "red2", lty = 2, lwd = 2)
  
  h_sort_df <- as.data.frame(h_sort)
  colnames(h_sort_df)[1] <- "cluster"
  h_sort_df$targ_gene <- rownames(h_sort_df)
  
  dendrogram <- as.dendrogram(row_clust)
  cols_branches <- c("darkred", "forestgreen", "orange", "firebrick1", "yellow", 
                     "deeppink1", "cyan2", "darkslategray", "chartreuse")
  # Set the colors of 9 branches
  dendrogram <- color_branches(dendrogram, k = length(row_clust), col = cols_branches)
  col_labels <- get_leaves_branches_col(dendrogram)
  #col_labels <- list("1" = "Serine glycine biosynthesis", "2" ="GABA-B receptor II signaling", 
  # "3" = "Fructose galactose metabolism/ Glycolysis", "4" = "N/A", 
  # "5" = "S-adenosylmethionine biosynthesis", "6" = "Sulfate assimilation", 
  # "7" = "Fructose galactose metabolism/ Androgen/etrogene/progesterone biosynthesis", 
  # "8" = "2-arachidonoylglycerol biosynthesis", "9" = "O-antigen biosynthesis")
  col_labels <- col_labels[order(order.dendrogram(dendrogram))]
  
  leaf_cols <- as.data.frame(leaf_Colors(dendrogram))
  leaf_cols$targ_gene <- rownames(leaf_cols)
  colnames(leaf_cols)[1] <- "colors"
  h_sort_df <- merge(h_sort_df, leaf_cols, by = "targ_gene")
  h_sort_df <- h_sort_df[order(h_sort_df$cluster),]
  
  hm_col <- heatmap.2(matrix, scale = "none", col = bluered(100), trace = "none", density.info = "none",
                      Rowv = dendrogram, key.xlab = "Beta", RowSideColors = col_labels,
                      colRow = col_labels)
  plot(hm_col)
}


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
#### HEAT MAP VISUALIZATION FOR TOP MUTATED GENES
############################################################
#' Given a master DF for a particular set of patients (e.g. full set or subtype)
#' that was run on multiple candidate or query genes on a set of target genes,
#' create an output figure that visualizes the Betas where query genes are labeled
#' on the left and names of target genes (on x-axis) are ignored 
#' @param results_table results_table a master DF produced from linear_model.R
#' @param outpath a path to a directory to write the t-statistic matrix to
create_hm_for_top_mutated <- function(results_table) {
  
  # Hack- limit to just the top 4 mutated
  results_table <- results_table[results_table$R_i.name %in% 
                                   c("TP53", "PIK3CA", "TTN", "RYR2", "KMT2C"),]
  
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


