##############################################################################
# CREATE FIGURES FOR TP53/ PIK3CA PAPER
# Created By: Sara Geraghty, Dec. 2021
##############################################################################

# Compiles the code for making these figures into one file
# All other plots were made in BioRender.com

library(ggplot2)
library(gplots)
library(RColorBrewer)
library(ggVennDiagram)
library("gage")
library(dendextend)
library("pheatmap")
library(ComplexHeatmap)
library(data.table)
library('ReactomePA')
library('clusterProfiler')
library("DOSE")
library("enrichplot")
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"


# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


##############################################################################
### VENN DIAGRAMS TP53 & PIK3CA MUTATIONS
##############################################################################
#' Creates a Venn Diagram to show how many patients have missense mutations in the
#' given patient set in the given genes 1 and 2
#' @param patient_set a vector of TCGA IDs of patients of interest (XXXX)
#' @param mut_count_matrix a mutation count matrix (patients are columns, rows are genes)
#' @param gene1 the external name of gene 1
#' @param gene2 the external name of gene 2
plot_mutational_overlap <- function(patient_set, mut_count_matrix, gene1, gene2) {
  # Limit the mutation count matrix to the given patient set
  colnames(mut_count_matrix) <- unlist(lapply(colnames(mut_count_matrix), function(x) {
    return(unlist(strsplit(x, "-", fixed = TRUE))[3])
  }))
  mut_count_matrix_sub <- mut_count_matrix[,colnames(mut_count_matrix) %in% patient_set]

  # Get the number of tumor samples with a mutation in each gene (count >= 1)
  mut_gene1_counts <- as.numeric(unlist(mut_count_matrix_sub[rownames(mut_count_matrix_sub) == gene1, ]))
  gene1_mut_samps <- which(mut_gene1_counts >= 1)
  
  mut_gene2_counts <- as.numeric(unlist(mut_count_matrix_sub[rownames(mut_count_matrix_sub) == gene2, ]))
  gene2_mut_samps <- which(mut_gene2_counts >= 1)
  
  # Get the number of tumor samples without a mutation in either gene (count = 0)
  #mut_neither_counts <- mut_count_matrix_sub[(rownames(mut_count_matrix_sub) == gene1) | 
                                                              #(rownames(mut_count_matrix_sub) == gene2), ]
  #neither_mut_samps <- as.numeric(unlist(which(colSums(mut_neither_counts) == 0)))
  
  mut_samps <- list("gene1" = gene1_mut_samps, "gene2" = gene2_mut_samps, "all" = 1:ncol(mut_count_matrix_sub))

  # Make the Venn Diagram
  ggVennDiagram(mut_samps, label_alpha = 0, category.names = c(gene1, gene2, "All"), set_color = "black",
                set_size = 10, label_size = 8, edge_size = 0) +
    ggplot2::scale_fill_gradient(low="cornsilk1", high = "cadetblue3")
}

patient_set <- read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt"), header = TRUE)[,1]
patient_set_lt20pamp <- read.table(paste0(main_path, "Patient Subsets/LessThan20PercAmp_patient_ids.txt"), header = TRUE)[,1]
patient_set_lumA <- read.table(paste0(main_path, "Patient Subsets/Luminal.A_patient_ids.txt"), header = TRUE)[,1]
patient_set_lumB <- read.table(paste0(main_path, "Patient Subsets/Luminal.B_patient_ids.txt"), header = TRUE)[,1]
patient_set_lumAB <- read.table(paste0(main_path, "Patient Subsets/Luminal.A.B_patient_ids.txt"), header = TRUE)[,1]
patient_set_basal <- read.table(paste0(main_path, "Patient Subsets/Basal_patient_ids.txt"), header = TRUE)[,1]
patient_set_her2 <- read.table(paste0(main_path, "Patient Subsets/HER2_patient_ids.txt"), header = TRUE)[,1]
patient_set_normLike <- read.table(paste0(main_path, "Patient Subsets/Normal-like_patient_ids.txt"), header = TRUE)[,1]

patient_set_pc <- read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt"), header = TRUE)[,1]

blca_subtype <- TCGAquery_subtype(tumor = "BLCA")
hnsc_subtype <- TCGAquery_subtype(tumor = "HNSC")
blca_subtype$justPat <- unlist(lapply(blca_subtype$patient, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))
hnsc_subtype$patient <- as.character(hnsc_subtype$patient)
hnsc_subtype$justPat <- unlist(lapply(hnsc_subtype$patient, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))

patient_set_blca <- intersect(patient_set_pc, blca_subtype$justPat)
patient_set_hnsc <- intersect(patient_set_pc, hnsc_subtype$justPat)

mutation_count_matrix <- read.csv(paste0(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense.csv"), 
                                  header = TRUE, row.names = 1, check.names = FALSE)
#mutation_count_matrix <- read.csv(paste0(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense_ALL.csv"), 
                                  #header = TRUE, row.names = 1, check.names = FALSE)

gene1 <- "TP53"
gene2 <- "PIK3CA"

# Call function
plot_mutational_overlap(patient_set, mutation_count_matrix, gene1, gene2)


##############################################################################
### MODEL RESULT HEAT MAPS
##############################################################################
#' Function plots a regulatory protein vs. gene target
#' clustered heat map, with entries being the t-statistic
#' of the hypothesis test
#' @param results_table a master DF produced from linear_model.R
#' @param outpath a path to a directory to write the t-statistic matrix to
create_heat_map <- function(results_table) {
  # Create the regulatory protein vs. gene target table
  matrix <- data.frame(matrix(ncol = length(unique(results_table$R_i.name)),
                              nrow = length(unique(results_table$T_k.name))))
  colnames(matrix) <- unique(results_table$R_i.name)
  rownames(matrix) <- unique(results_table$T_k.name)
  
  # Fill in this table with t-statistics
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

# Recon3D results for TP53/ PIK3CA (No overlap)
master_df_mut <- read.csv(paste0(main_path, "Linear Model/PIK3CA_TP53/eQTL/output_results_NoMutOverlap_P53_PIK3CA_metabolicTargs_iprotein_tmmSDGr1_rawCNA_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_corrected_MUT.csv"),
                          header = TRUE, check.names = FALSE)
master_df_mut_lumA <- read.csv(paste0(main_path, "Linear Model/PIK3CA_TP53/eQTL/output_results_LumA.NoMutOverlap_P53_PIK3CA_metabolicTargs_iprotein_tmmSDGr1_rawCNA_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_corrected_MUT.csv"),
                               header = TRUE, check.names = FALSE)
master_df_mut_lumB <- read.csv(paste0(main_path, "Linear Model/PIK3CA_TP53/eQTL/output_results_LumB.NoMutOverlap_P53_PIK3CA_metabolicTargs_iprotein_tmmSDGr1_rawCNA_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_corrected_MUT.csv"),
                               header = TRUE, check.names = FALSE)
master_df_mut_lumAB <- read.csv(paste0(main_path, "Linear Model/PIK3CA_TP53/eQTL/output_results_LumAB.NoMutOverlap_P53_PIK3CA_metabolicTargs_iprotein_tmmSDGr1_rawCNA_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_corrected_MUT.csv"),
                                header = TRUE, check.names = FALSE)
master_df_mut_lt20pamp <- read.csv(paste0(main_path, "Linear Model/PIK3CA_TP53/eQTL/output_results_LT20PAmp.NoMutOverlap_P53_PIK3CA_metabolicTargs_iprotein_tmmSDGr1_rawCNA_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_corrected_MUT.csv"),
                                   header = TRUE, check.names = FALSE)

# Call this function
create_heat_map(master_df_mut, main_path)


#' Function plots a regulatory protein vs. gene target
#' clustered heat map, with entries being the t-statistic
#' of the hypothesis test. Adds an additional dimension of coloring
#' the clusters, eliminating x-axis clustering
#' @param master_df a master DF produced from linear_model.R
#' @param height a height at which to cut the dendrogram to create clusters
create_heatmap_with_clustering <- function(master_df, height) {
  matrix <- create_heat_map(master_df)
  print(head(matrix))
  hm <- heatmap.2(matrix, scale = "none", col = bluered(1000), trace = "none",
                  density.info = "none", key = TRUE, key.title = NA, key.xlab = "Beta")
  rownames(matrix)[hm$rowInd]
  row_clust <- hclust(dist(matrix, method = "euclidean"), method = 'ward.D2')
  plot(row_clust)
  
  # Cut the dendrogram using the given height
  h_sort <- sort(cutree(row_clust, h = height))
  #plot(row_clust)
  #abline(h = height, col = "red2", lty = 2, lwd = 2)
  
  # Add the new clusters to the DF
  h_sort_df <- as.data.frame(h_sort)
  colnames(h_sort_df)[1] <- "cluster"
  h_sort_df$targ_gene <- rownames(h_sort_df)
  num_clust <- length(unique(h_sort_df$cluster))
  dendrogram <- as.dendrogram(row_clust)
  
  # Use preset colors, or color according to the number of clusters
  #cols_branches <- c("darkred", "forestgreen", "orange", "firebrick1", "yellow", 
                     #"deeppink1", "cyan2", "darkslategray", "chartreuse")
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  print(col_vector)
  cols_branches <- sample(col_vector, num_clust)
  
  # Add the colors to the dendrogram
  dendrogram <- color_branches(dendrogram, k = num_clust, col = cols_branches)
  col_labels <- get_leaves_branches_col(dendrogram)
  col_labels <- col_labels[order(order.dendrogram(dendrogram))]
  
  # Add the labels and colors to the DF so we can ID them and search them later
  h_sort_df$labels <- unlist(lapply(1:num_clust, function(i) {
    num_entries <- nrow(h_sort_df[h_sort_df$cluster == i,])
    return(rep(col_labels[[i]], times = num_entries))
  }))
  h_sort_df$colors <- unlist(lapply(1:num_clust, function(i) {
    num_entries <- nrow(h_sort_df[h_sort_df$cluster == i,])
    return(rep(cols_branches[[i]], times = num_entries))
  }))
  
  # Plot the new and improved heatmap
  new_hm <- heatmap.2(matrix, scale = "none", col = bluered(1000), trace = "none", 
                      density.info = "none", dendrogram = "row", Rowv = dendrogram, 
                      Colv = "none", key.title = NA, key.xlab = "Beta", 
                      RowSideColors = col_labels, colRow = col_labels) 
  # col = brewer.pal(n = 100, name = "YlGnBu")
  print(new_hm)
  
  # Return the labeled DF
  return(h_sort_df)
}

# Recon3D results for TP53/ PIK3CA (No overlap)
master_df_mut <- read.csv(paste0(main_path, "Linear Model/PIK3CA_TP53/eQTL/output_results_NoMutOverlap_P53_PIK3CA_metabolicTargs_iprotein_tmmSDGr1_rawCNA_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_corrected_MUT.csv"),
                          header = TRUE, check.names = FALSE)
master_df_mut_lumA <- read.csv(paste0(main_path, "Linear Model/PIK3CA_TP53/eQTL/output_results_LumA.NoMutOverlap_P53_PIK3CA_metabolicTargs_iprotein_tmmSDGr1_rawCNA_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_corrected_MUT.csv"),
                          header = TRUE, check.names = FALSE)
master_df_mut_lumB <- read.csv(paste0(main_path, "Linear Model/PIK3CA_TP53/eQTL/output_results_LumB.NoMutOverlap_P53_PIK3CA_metabolicTargs_iprotein_tmmSDGr1_rawCNA_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_corrected_MUT.csv"),
                          header = TRUE, check.names = FALSE)
master_df_mut_lumAB <- read.csv(paste0(main_path, "Linear Model/PIK3CA_TP53/eQTL/output_results_LumAB.NoMutOverlap_P53_PIK3CA_metabolicTargs_iprotein_tmmSDGr1_rawCNA_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_corrected_MUT.csv"),
                          header = TRUE, check.names = FALSE)
master_df_mut_lt20pamp <- read.csv(paste0(main_path, "Linear Model/PIK3CA_TP53/eQTL/output_results_LT20PAmp.NoMutOverlap_P53_PIK3CA_metabolicTargs_iprotein_tmmSDGr1_rawCNA_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_corrected_MUT.csv"),
                          header = TRUE, check.names = FALSE)


# Call this function
create_heatmap_with_clustering(master_df_mut, height = 2)



##############################################################################
### SPEARMAN CORRELATION PLOTS
##############################################################################
#' Compute and print the Spearman correlation of the Betas and of the T-statistic,
#' given two groups of interest
#' @param results_table the output master DF from the linear model
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
compute_and_print_spearman(master_df_mut, ri_1, ri_2)

##############################################################################
### COMPARATIVE BOXPLOTS FOR INDIVIDUAL CASES
##############################################################################

##############################################################################
### VENN DIAGRAMS TOP HIT OVERLAP WITH HYPERGEOMETRIC TEST STATISTICS
##############################################################################
#' Plots the overlap in significant top target gene hits for mutation 
#' and CNA results. 
#' @param master_df_mut_sig mutation master DF, with gene names added and 
#' thresholded for significance
#'@param master_df_cna_sig CNA master DF, with gene names added and 
#' thresholded for significance
plot_tophit_overlap <- function(master_df_mut_sig, master_df_cna_sig) {
  overlap_genes <- intersect(master_df_mut_sig$T_k.name, master_df_cna_sig$T_k.name)
  print(overlap_genes)
  
  # Plot a Venn Diagram
  #myCol <- brewer.pal(2, "Pastel2")
  plt <- venn.diagram(list(master_df_mut_sig$T_k.name, master_df_cna_sig$T_k.name),
                      category.names = c("Mutation", "CNA"), filename = NULL, output = TRUE,
                      lwd = 2, lty = 'blank', fill = c("red", "blue"), cex = 0.6, fontface = "bold",
                      fontfamily = "sans") #, cat.cex = 0.6, cat.fontface = "bold",
                      #cat.default.pos = "outer", cat.fontfamily = "sans", rotation = 1)
  grid::grid.draw(plt)
}

# Call this function and write to PNG file
plot_tophit_overlap(fn, master_df_mut_sig, master_df_cna_sig)



##############################################################################
### GENE SET ENRICHMENT ANALYSIS CURVES
##############################################################################
#' Takes in a results table from the linear model output and plots GSEA
#' for n upregulated and n downregulated terms using the 
#' ReactomePA package
#' @param results_table a results data frame from my model
#' @param n the number of up/downregulated terms to plot
#' @param all_genes_id_conv a gene ID conversion data frame from BioMart
#' @param output_path a local path to save the GSEA figures to
perform_gsea <- function(results_table, n, all_genes_id_conv, output_path) {
  
  regprot.gsea.rp <- list()
  regprot.gsea.ncg <- list()
  #regprot.gsea.go <- list()
  #regprot.gsea.kegg <- list()
  
  # Loop through all the regulatory proteins
  for (regprot in unique(results_table$R_i.name)) {
    
    # Get all of this protein's targets and their associated Betas, sorted in 
    # descending order
    res_table_sub <- results_table[results_table$R_i.name == regprot, c('T_k.name', 'estimate')]
    
    mapping <- as.data.frame(bitr(res_table_sub$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                                  OrgDb="org.Hs.eg.db", drop = FALSE))
    colnames(mapping) <- c("T_k.name", "T_k.entrez")
    res_table_sub <- merge(res_table_sub, mapping, all=TRUE, by="T_k.name")
    res_table_sub <- res_table_sub[order(res_table_sub$estimate, decreasing = TRUE),]
    print(res_table_sub)
    #res_table_sub <- res_table_sub[!(is.na(res_table_sub))]
    #print(res_table_sub)
    
    plotdir <- paste0(output_path, "Linear Model/GSEA/")    
    suppressWarnings(dir.create(plotdir, recursive = TRUE))
    
    regprotBetaScores <- res_table_sub$estimate
    names(regprotBetaScores) <- res_table_sub$T_k.entrez
    
    gse.rp <- gsePathway(regprotBetaScores)
    gse.rp <- setReadable(gse.rp, 'org.Hs.eg.db')
    gse.ncg <- gseNCG(regprotBetaScores)
    gse.ncg <- setReadable(gse.ncg, 'org.Hs.eg.db')
    #gse.go <- gseGO(regprotBetaScores, OrgDb = 'org.Hs.eg.db')
    #gse.go <- setReadable(gse.go, 'org.Hs.eg.db')
    #gse.kegg <- gseKEGG(regprotBetaScores)
    #gse.kegg <- setReadable(gse.kegg, 'org.Hs.eg.db')
    
    regprot.gsea.rp[[regprot]] <- gse.rp
    regprot.gsea.ncg[[regprot]] <- gse.ncg
    #regprot.gsea.go[[regprot]] <- gse.go
    #regprot.gsea.kegg[[regprot]] <- gse.kegg
    
    # plot GSEA for top N upregulated and top N downregulated terms
    #for (gst in c('rp','ncg', 'go', 'kegg')) {
    for (gst in c('rp','ncg')) {
      if (gst == 'rp') res <- gse.rp
      #else if (gst == 'kegg') res <- gse.kegg
      #else if (gst == 'go') res <- gse.go
      else res <- gse.ncg
      
      terms <- res@result$Description[res@result$enrichmentScore > 0][1:n]
      terms <- terms[!is.na(terms)]
      if (length(terms) > 0) {
        termIDs <- sapply(terms, function(x) {res@result$ID[which(res@result$Description == x)]})
        pdf(file = paste0(plotdir, regprot,'_gse_',gst,'_gseaCurve_top', n, 'terms_up.pdf'),
            width = 15, height = 10)
        plot(gseaplot2(res, geneSetID = termIDs, pvalue_table = TRUE,
                       color = RColorBrewer::brewer.pal(length(termIDs), "Dark2"), ES_geom = "dot"))
        dev.off()
      }
      terms <- res@result$Description[res@result$enrichmentScore < 0][1:n]
      terms <- terms[!is.na(terms)]
      if (length(terms) > 0) {
        termIDs <- sapply(terms, function(x) {res@result$ID[which(res@result$Description == x)]})
        pdf(file = paste0(plotdir, regprot,'_gse_',gst,'_gseaCurve_top', n, 'terms_down.pdf'),
            width = 15, height = 10)
        plot(gseaplot2(res, geneSetID = termIDs, pvalue_table = TRUE,
                       color = RColorBrewer::brewer.pal(length(termIDs), "Dark2"), ES_geom = "dot"))
        dev.off()
      }
    }
  }
}

perform_gsea(master_df_mut_corrected, n = 5, all_genes_id_conv, output_path)
perform_gsea(master_df_cna_corrected, n = 5, all_genes_id_conv, output_path)


