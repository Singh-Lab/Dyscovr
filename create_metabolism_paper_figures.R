##############################################################################
# CREATE FIGURES FOR TP53/ PIK3CA PAPER
# Created By: Sara Geraghty, Dec. 2021
##############################################################################

# Compiles the code for making these figures into one file
# All other plots were made in BioRender.com

library(ggplot2)
library(gplots)
library(maftools)
library(stats)
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
library(ggrepel)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)


# NEJM color palatte: https://nanx.me/ggsci/reference/pal_nejm.html


main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/BRCA/GSEA/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/Pan-Cancer/GSEA/"


# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)



##############################################################################
##############################################################################
### FIGURE 1
##############################################################################
##############################################################################

# PART A: METHODOLOGICAL OVERVIEW: CREATED IN BIORENDER

# PART B: TOP MUTATED DRIVERS BY CANCER TYPE/SUBTYPE

#' Takes a MAF file from the TCGA, along with a list of known driver genes,
#' and creates a barplot of the mutation frequency of the top n drivers for
#' each cancer type or subtype. If a synonymous MAF file is given, also plots
#' the synonymous mutation frequency for each of these genes
#' @param maf a nonsynonymous MAF file from the TCGA 
#' @param driver_gene_df a data frame with known driver genes
#' @param n_drivers the number of top mutated drivers to include per cancer type or subtype
#' @param clinical_df a clinical DF from the TCGA to relate sample IDs to their corresponding
#' cancer type or subtype
#' @param ct_or_subtype either "cancer_type" or "subtype" to denote whether we want to look 
#' across cancer types, or across subtypes within a given cancer type
#' @param synony_maf_filename OPT: another synonymous MAF file; if given, will create a stacked barplot 
#' with nonsynonymous vs. synonymous mutation frequency
create_top_mutated_drivers_by_cancer_type_barplot <- function(maf, driver_gene_df, n_drivers, clinical_df, 
                                                              ct_or_subtype, synony_maf_filename) {
  
  # Get all the IDs in the MAF file
  maf_file_barcodes <- unlist(lapply(as.character(unlist(maf@data$Tumor_Sample_Barcode)), function(x)
    paste(unlist(strsplit(x, "-", fixed = TRUE))[1:3], collapse = "-")))
  
  # Split MAF file according to cancer type or subtypes
  maf_list <- list()
  if (ct_or_subtype == "cancer_type") {
    
    # Get all the different cancer types
    cts <- unique(clinical_df$project_id)
    
    # Split MAFs into each cancer type
    maf_list <- lapply(cts, function(ct) {
      cancer_type_patients <- as.character(unlist(clinical_df[clinical_df$project_id == ct, 
                                                              'tcga_barcode']))
      maf_sub <- maf
      maf_sub@data <- maf_sub@data[which(maf_file_barcodes %fin% cancer_type_patients),]
      return(maf_sub)
    })
    names(maf_list) <- cts
  
    
  } else if (ct_or_subtype == "subtype") {
    
    subtype_file <- TCGAquery_subtype(tumor = unlist(strsplit(clinical_df$project_id[1], "-", fixed = TRUE))[2])
    subtypes <- unique(subtype_file$BRCA_Subtype_PAM50)    #TODO: make this generalizable to other cancer types if needed
    subtypes <- subtypes[subtypes != "NA"]
    maf_list <- lapply(subtypes, function(st) {
      subtype_patients <- as.character(unlist(subtype_file[subtype_file$BRCA_Subtype_PAM50 == st, 'patient']))
      maf_sub <- maf
      maf_sub@data <- maf@data[which(maf_file_barcodes %fin% subtype_patients),]
      return(maf_sub)
    })
    names(maf_list) <- subtypes
  
  } else {print("Error: must be either 'cancer_type' or 'subtype'.")}
  
  #print(head(maf_list))
  
  # Create a stacked plot
  line = 7
  cex = 1.5
  side = 2
  las = 3
  par(mfrow = c(length(maf_list), 1), oma=c(1,6,1,1))
  
  # Order alphabetically
  maf_list <- maf_list[order(names(maf_list))]
  maf_list_copy <- maf_list
  
  # Create plots for each
  for (i in 1:length(maf_list)) {
    # Subset to only driver genes
    m_sub <- subsetMaf(maf_list[[i]], genes = driver_gene_df$primary_gene_names)
    #m_sub_data <- m_sub@data[as.character(unlist(m_sub@data$Gene)) %fin% driver_gene_df$ensembl_gene_id,]

    # Create barplot
    colors <- c("#0072B5FF", "#BC3C29FF", "#20854EFF", "#E18727FF")
    names(colors) <- c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site")
    mafbarplot(m_sub, n = n_drivers, fontSize = 1.5, legendfontSize = 1.25, borderCol = NA, 
               showPct = FALSE, color = colors)
    mtext(names(maf_list)[i], side=side, line=line, cex=cex, las=las)
  }
  
  return(maf_list_copy)
}

# Get maf
maf_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/Somatic_Mut_Data/"
maf_filename <- paste0(maf_path, "TCGA.BRCA.muse.b8ca5856-9819-459c-87c5-94e91aca4032.DR-10.0.somatic.maf")
# Import the MAF file using maftools
maf <- read.maf(maf_filename)
# Remove Translation_Start_Site mutations (not using)
maf <- subsetMaf(maf, query = "Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation',
                                                              'Nonstop_Mutation', 'Splice_Site')")

# Get clinical DF
clinical_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/clinical_data_subset.csv", 
                        header = TRUE, check.names = FALSE)

# Call function
maf_list_by_subtype <- create_top_mutated_drivers_by_cancer_type_barplot(maf, driver_gene_df, 5, 
                                                                         clinical_df, "subtype", NA)


dev.off()


# PART C: SIGNIFICANT CORRELATIONS BY CANCER TYPE/ SUBTYPE SIZE AND # OF SAMPLE MUTATIONS

#' Takes a list of result master data frames under consideration, along with a dictionary
#' with the size of each corresponding cancer type/ subtype, and a MAF file with nonsynonymous
#' mutations, and produces a scatter plot relating signif. count to sample size and mutation
#' frequency for a given gene.
#' @param list_of_master_dfs a list of master DFs, with the name being the cancer type or subtype
#' @param sample_size_dict a list relating each cancer type/ subtype to its sample size
#' @param list_of_maf_files a list of MAF file with mutation frequencies for all genes, one per cancer type/subtype
#' @param goi a driver gene-of-interest
#' @param qval_thres a q-value threshold for significance
relate_sample_size_to_sig_eQTL <- function(list_of_master_dfs, sample_size_dict, list_of_maf_files, goi, qval_thres) {
  
  # Create a data frame that will link these factors
  input_df <- data.frame("type" = names(sample_size_dict), "sample.size" = as.numeric(sample_size_dict))
  
  # Get the mutation frequency of the given GOI for each cancer type/subtype
  # Import get_mut_count_matrix from maf_helper_functions.R
  total.mut.count <- unlist(lapply(1:length(list_of_maf_files), function(i) {
    m <- list_of_maf_files[[i]]
    mut_count_matrix <- get_mut_count_matrix(m)
    mut_count <- rowSums(mut_count_matrix[rownames(mut_count_matrix) == goi, ])
    return(mut_count)
  }))
  names(total.mut.count) <- names(list_of_maf_files)
  total.mut.count <- as.data.frame(total.mut.count)
  total.mut.count$type <- rownames(total.mut.count)
  print(head(total.mut.count))

  # Merge this with our input DF by type
  input_df <- merge(input_df, total.mut.count, by = "type")
  
  # Get the number of significant eQTLs for this goi
  num.sig.hits <- unlist(lapply(1:length(list_of_master_dfs), function(i) {
    m <- list_of_master_dfs[[i]]
    m_sig <- m[(m$q.value < qval_thres) & (m$R_i.name == goi),]
    return(nrow(m_sig))
  }))
  names(num.sig.hits) <- names(list_of_master_dfs)
  num.sig.hits <- as.data.frame(num.sig.hits)
  num.sig.hits$type <- rownames(num.sig.hits)
  print(head(num.sig.hits))

  # Merge this with our input DF by type
  input_df <- merge(input_df, num.sig.hits, by = "type")
  print(head(input_df))
  
  # Now, create the plot with this input DF
  p <- ggplot(input_df, aes(x = sample.size, y = num.sig.hits, color = type)) + 
    geom_point(aes(size = total.mut.count)) + scale_color_nejm() +
    theme_minimal() + xlab("Sample Size") + ylab(paste("Num. Significant Correlations (q <", paste0(qval_thres, ")"))) +
    labs(size = "Total Nonsyn. Mut. Count", color = "BRCA Subtype") + 
    theme(axis.title = element_text(face = "bold", size = 14), 
          axis.text = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12))
  p
  
  # Alternative version with total mut count vs. num sig hits
  #p <- ggplot(input_df, aes(x = total.mut.count, y = num.sig.hits, color = type)) + 
  #  geom_point(aes(size = sample.size)) + scale_color_nejm() +
  #  theme_minimal() + xlab("Total Nonsyn. Mut. Count") + ylab(paste("Num. Significant Correlations (q <", paste0(qval_thres, ")"))) +
  #  labs(size = "Sample Size", color = "BRCA Subtype")
  #p
}

#TODO: IMPORT RELEVANT MASTER DFs

allgenes_p53 <- read.csv("C:/Users/sarae/Documents/res_P53_allGenes_iprotein_tmmSDGr1filtByExpr_rawCNA_methMRaw_cibersortTotalFrac_2PCs_rmCis_rmMetast_PIK3CA_covs.inclCT.inclGenoPCs.MN_corrected_MUT.csv",
                         header = TRUE, check.names = FALSE)
allgenes_pik3ca <- read.csv("C:/Users/sarae/Documents/res_PIK3CA_allGenes_iprotein_tmmSDGr1filtByExpr_rawCNA_methMRaw_cibersortTotalFrac_2PCs_rmCis_rmMetast_TP53_Covs.incl2GenoPCs.inclCT.MN_corrected_MUT.csv",
                            header = TRUE, check.names = FALSE)
allgenes_sf3b1 <- read.csv("C:/Users/sarae/Documents/res_SF3B1_allGenes_iprotein_tmmSDGr1filtByExpr_rawCNA_methMRaw_cibersortTotalFrac_2PCs_rmCis_rmMetast_TP53_PIK3CA_Covs.inclCT.MN_corrected_MUT.csv",
                           header = TRUE, check.names = FALSE)


metabol_p53_basal <- read.csv("C:/Users/sarae/Documents/res_Basal_P53_metabolicTargs_iprotein_tmmSDGr1filtByExpr_rawCNA_methMRaw_cibersortTotalFrac_rmCis_rmMetast_PIK3CA_covs.inclCT.MN_corrected_MUT.csv",
                               header = TRUE, check.names = FALSE)
metabol_p53_her2 <- read.csv("C:/Users/sarae/Documents/res_HER2_P53_metabolicTargs_iprotein_tmmSDGr1filtByExpr_rawCNA_methMRaw_cibersortTotalFrac_rmCis_rmMetast_PIK3CA_covs.inclCT.MN_corrected_MUT.csv",
                              header = TRUE, check.names = FALSE)
metabol_p53_lumA <- read.csv("C:/Users/sarae/Documents/res_LumA_P53_metabolicTargs_iprotein_tmmSDGr1filtByExpr_rawCNA_methMRaw_cibersortTotalFrac_rmCis_rmMetast_PIK3CA_covs.inclCT.MN_corrected_MUT.csv",
                             header = TRUE, check.names = FALSE)
metabol_p53_lumB <- read.csv("C:/Users/sarae/Documents/res_LumB_P53_metabolicTargs_iprotein_tmmSDGr1filtByExpr_rawCNA_methMRaw_cibersortTotalFrac_rmCis_rmMetast_PIK3CA_covs.inclCT.MN_corrected_MUT.csv",
                             header = TRUE, check.names = FALSE)

list_of_master_dfs <- list("LumA" = metabol_p53_lumA, "LumB" = metabol_p53_lumB, 
                           "Her2" = metabol_p53_her2, "Basal" = metabol_p53_basal)

# Import intersecting patients
intersecting_patients <- read.table(paste(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt", sep = ""))[,1]

# Import a subtype DF to relate patients to subtype & create a dictionary of sample size per subtype
brca_subtype_file <- TCGAquery_subtype(tumor = "BRCA")

subtype_file_patient_ids <- unlist(lapply(brca_subtype_file$patient, function(x)
  unlist(strsplit(x, "-", fixed = TRUE))[3]))
brca_subtype_file_sub <- brca_subtype_file[which(subtype_file_patient_ids %in% intersecting_patients),]

subtypes_of_interest <- c("LumA", "LumB", "Her2", "Basal")
list_of_sample_sizes <- lapply(subtypes_of_interest, function(st) 
  nrow(brca_subtype_file_sub[brca_subtype_file_sub$BRCA_Subtype_PAM50 == st,]))
names(list_of_sample_sizes) <- subtypes_of_interest

# Use the generated list of MAFs subsetted by subtype from above (part B)

# Call function
relate_sample_size_to_sig_eQTL(list_of_master_dfs, list_of_sample_sizes, maf_list_by_subtype, "TP53", 0.1)



# PART D: OVERLAP IN PARTICULAR CORRELATIONS BETWEEN SUBGROUPS

#' Plot overlap in significant hits for a GOI between different subgroups using a Venn
#' diagram format
#' @param list_of_master_dfs a list of master DFs, with the name being the cancer type or subtype
#' @param goi a driver gene-of-interest
#' @param qval_thres a q-value threshold for significance
plot_overlap_between_subgroups <- function(list_of_master_dfs, goi, qval_thres) {
  sig_hits <- lapply(list_of_master_dfs, function(m) 
    m[(m$R_i.name == goi) & (m$q.value < qval_thres), 'T_k.name'])
  
  plt <- venn.diagram(sig_hits, category.names = names(list_of_master_dfs), 
                      filename = NULL, output = TRUE, lwd = 2, lty = 'blank', 
                      fill = c("#0072B5FF", "#BC3C29FF", "#20854EFF", "#E18727FF"), cex = 2, fontface = "bold",
                      fontfamily = "sans", cat.cex = 2, cat.fontface = "bold", cat.fontfamily = "sans")
  #cat.default.pos = "outer", cat.pos = c(180, 180, 180, 180)) #, cat.fontfamily = "sans", rotation = 1)
  grid::grid.draw(plt)
}


plot_overlap_between_subgroups(list_of_master_dfs, "TP53", 0.1)


# PART E: BETA HISTOGRAM WITH NUMBER OF SIGNIFICANT CORRELATIONS LABELED

#' Plot a multi-layer histogram (per cancer type or subtype) showing all the
#' Beta values for a particular GOI with a significant hits (defined as being 
#' below a given q-value threshold) highlighted
#' @param list_of_master_dfs a list of master DFs, with the name being the cancer type or subtype
#' @param goi a driver gene-of-interest
#' @param qval_thres a q-value threshold for significance
plot_beta_histograms_by_subtype <- function(list_of_master_dfs, goi, qval_thres) {
  betas <- lapply(list_of_master_dfs, function(m) 
    m[m$R_i.name == goi, 'estimate'])
  names(betas) <- names(list_of_master_dfs)
  betas <- as.data.frame(betas)
  betas_m <- melt(betas)
  print(head(betas_m))
  
  means_of_betas <- unlist(lapply(1:ncol(betas), function(i) mean(betas[,i])))
  print(head(means_of_betas))
  
  p <- ggplot(betas_m, aes(x=value, fill = variable)) + geom_density(alpha = 0.3) +
    scale_color_nejm() + theme_minimal() + xlab("Beta estimate") + ylab("Density") + 
    labs(fill = "BRCA Subtype") + 
    theme(axis.title = element_text(face = "bold", size = 14), 
          axis.text = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12))
    #+ geom_vline(xintercept=means_of_betas, size=0.25, color="black", linetype = "dotted", alpha = 0.2) # add the mean of the betas for each
  
  p
}


plot_beta_histograms_by_subtype(list_of_master_dfs, "TP53", 0.1)


##############################################################################
##############################################################################
### FIGURE 2
##############################################################################
##############################################################################

# PART A: HEAT MAP OF DRIVER VS. METABOLIC TARGET BETAS, CLUSTERED BY GENE

#' Function plots a regulatory protein vs. gene target
#' clustered heat map, with entries being the t-statistic
#' of the hypothesis test
#' @param results_table a master DF produced from linear_model.R
#' @param outpath a path to a directory to write the t-statistic matrix to
create_heat_map <- function(results_table) {
  
  # Use helper function to create input matrix from results table
  matrix <- create_regprot_v_genetarg_matrix(results_table)
  
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


#' Helper function to create regprot vs. gene target matrix and fill with Beta values 
#' @param results_table
create_regprot_v_genetarg_matrix <- function(results_table) {
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
  
  return(matrix)
}



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


# Call this function
example_df <- read.csv("C:/Users/sarae/Documents/output_results_cancerRelated_mattMetabolicTargs_iprotein_tmmSDGr1_rawCNA_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_corrected_MUT.csv",
                       header = TRUE, check.names = FALSE)
create_heatmap_with_clustering(example_df, height = 2)



# PART B: TP53 UP- AND DOWN-REGULATED METABOLIC HITS
#' Create a volcano plot for a given GOI in a particular cancer type or subtype,
#' with significant hits highlighted and labeled
#' @param master_df a data frame produced by the model
#' @param goi a driver gene-of-interest
#' @param qval_thres a q-value threshold for significance
create_volcano_plot <- function(master_df, goi, qval_thres) {
  
  # Subset to the GOI 
  master_df <- master_df[master_df$R_i.name == goi,]
  master_df <- na.omit(master_df)

  # Get the log2(Beta), which we are considering here to be like log2(fold change)
  log2_beta <- unlist(lapply(master_df$estimate, function(e) { 
    if((is.nan(e)) | (length(e) == 0)) {return(0)}
    #else {return(log2(abs(e)))}
    else {return(e)}
  }))
  #print(head(log2_beta))
  
  # Get the -log10(qvalue)
  neg_log10_qval <- unlist(lapply(master_df$q.value, function(q) { 
    if((is.nan(q)) | (length(q) == 0)) {return(-log10(1))}
    else {return(-log10(q))}
  }))
  #print(head(neg_log10_qval))
  
  up_or_down <- unlist(lapply(1:nrow(master_df), function(i) {
    if(master_df[i, 'q.value'] < qval_thres) {
      beta_sign <- ifelse(master_df[i, 'estimate'] > 0, 1, 0) 
      if(beta_sign == 1) {return("up")}
      else {return("down")}
    } else {return("ns")}
  }))
  
  gene_names <- unlist(lapply(1:nrow(master_df), function(i) 
    ifelse(-log10(master_df$q.value[i]) > 1.75, master_df$T_k.name[i], NA)))
  print(head(gene_names))
  
  # Put this into a DF for plotting
  plot_df <- data.frame("log2_beta" = log2_beta, "neg_log10_qval" = neg_log10_qval,
                        "up_or_down" = up_or_down, "gene" = gene_names)
  

  # Set color and alpha params
  cols <- c("up" = "#0072B5FF", "down" = "#E18727FF", "ns" = "grey") 
  alphas <- c("up" = 0.95, "down" = 0.95, "ns" = 0.5)

  p <- ggplot(plot_df, aes(x = log2_beta, y = neg_log10_qval, fill = up_or_down, 
                           alpha = up_or_down, label = gene)) + 
    geom_point(shape = 21, size = 2.75) + xlab("Beta Coefficient") + 
    ylab("-log10(q-value)") + theme_minimal() + 
    geom_hline(yintercept = -log10(qval_thres), linetype = "dashed") +
    scale_fill_manual(values = cols) + scale_alpha_manual(values = alphas) + 
    theme(legend.position = "none") + geom_text_repel()
  
  p
}

create_volcano_plot(allgenes_p53, "TP53", 0.1)


# PART C: INDIVIDUAL BOX OR VIOLIN PLOTS FOR A SPECIFIC PATHWAY OF INTEREST

#' Given a driver gene of interest and a particular target gene, creates a boxplot
#' or violin plot of the gene expression differences between mutant and non-mutant
#' versions of the given driver gene
#' @param expression_df a normalized and filtered expression DF
#' @param mutation_regprot_df a mutation DF produced from process_mutation_data.R
#' @param driver_ensg the driver gene ENSG ID
#' @param driver_uniprot the driver gene Uniprot ID
#' @param target the target gene ENSG ID
#' @param box_or_violin either "box" or "violin" to denote whether we want a 
#' boxplot or a violin version of the plot
#' @param driver_name the gene name of the driver gene, for labeling
#' @param target_name the gene name of the target gene, for labeling
create_driver_target_boxplot <- function(expression_df, mutation_regprot_df, 
                                         driver_ensg, driver_uniprot, target, 
                                         box_or_violin, driver_name, target_name) {
  
  # Get mutant/ non-mutant patients
  mutant_patients <- get_mutant_patients(mutation_regprot_df, driver_uniprot)
  
  # Get expression of target gene within each mutational group
  expression_by_group <- get_expression_by_mut_group(expression_df, mutant_patients, target)
  expression_mutants <- expression_by_group[[1]]
  expression_normal <- expression_by_group[[2]]
  
  # Pad them with NAs as needed, in order to put them in a data frame
  if(length(expression_mutants) > length(expression_normal)) {
    nas <- rep(NA, times = length(expression_mutants) - length(expression_normal))
    expression_normal <- c(expression_normal, nas)
  } else if (length(expression_normal) > length(expression_mutants)) {
    nas <- rep(NA, times = length(expression_normal) - length(expression_mutants))
    expression_mutants <- c(expression_mutants, nas)
  } else {print("Same length. No NAs added.")}
    
  dataframe_expr <- data.frame("No.Mutation" = expression_normal, "Mutation" = expression_mutants)
  dataframe_expr_m <- na.omit(melt(dataframe_expr))
  print(head(dataframe_expr_m))
  
  # Create a boxplot or violin plot
  if(box_or_violin == "box") {
    #boxplot(dataList_expr, ylab = "Expression", main = paste(target, paste("Expression By", paste(driver_ensg, "Mutation Status"))))
    
    p <- ggplot(dataframe_expr_m, aes(x = variable, y = value, fill = variable)) + 
      geom_boxplot(size = 1) + geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.2) + 
      theme_minimal() + scale_fill_manual(values=c("#E18727FF", "#6F99ADFF")) +
      xlab(paste("\n", paste(driver_name, "Mutation Status"))) + ylab(paste(target_name, "Expression\n")) +
      theme(legend.position = "none", axis.title = element_text(face = "bold")) + 
      scale_x_discrete(labels=c(paste(driver_name, "Wild-Type"), paste(driver_name, "Mutant")))
  
  } else {
    p <- ggplot(dataframe_expr_m, aes(x = variable, y = value, fill = variable)) + 
      geom_violin(trim = FALSE) + geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.2) + 
      theme_minimal() + scale_fill_manual(values=c("#E18727FF", "#6F99ADFF")) +
      xlab(paste("\n", paste(driver_name, "Mutation Status"))) + ylab(paste(target_name, "Expression\n")) +
      theme(legend.position = "none", axis.title = element_text(face = "bold")) + 
      scale_x_discrete(labels=c(paste(driver_name, "Wild-Type"), paste(driver_name, "Mutant"))) +
      stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black", size = 1)  # FOR POINT RANGE
      #stat_summary(fun.data=mean_sdl, geom="crossbar", width=0.08)  # FOR CROSSBAR
  }
  p
}


#' Helper function that gets patients with mutation in the given driver
#' @param mutation_regprot_df a mutation DF produced from process_mutation_data.R
#' @param driver_uniprot the driver gene Uniprot ID
get_mutant_patients <- function(mutation_regprot_df, driver_uniprot) {
  driver_rows <- unlist(lapply(mutation_regprot_df$Query, function(x) {
    ifelse(unlist(strsplit(x, "|", fixed = TRUE))[2] == driver_uniprot, TRUE, FALSE)
  }))
  mutant_patients <- unique(unlist(strsplit(mutation_regprot_df[driver_rows, "Patient"], ";", fixed = TRUE)))
  return(mutant_patients)
}

#' Get the expression in each of the groups (mutated/ unmutated)
#' @param expression_df a data frame with expression values (columns are patient IDs-sample IDs, 
#' rows are ENSG IDs)
#' @param mutant_patients a vector of all the patients that have a mutation in the given
#' regulatory protein
#' @param ensg the ENSG ID of the target gene of interest
get_expression_by_mut_group <- function(expression_df, mutant_patients, ensg) {
  #colnames(expression_df) <- unlist(lapply(colnames(expression_df), function(x) unlist(strsplit(x, "-", fixed = TRUE))[1]))
  expression_mutants <- as.numeric(expression_df[rownames(expression_df) == ensg, 
                                                 colnames(expression_df) %in% mutant_patients])
  print(head(expression_mutants))
  expression_normal <- as.numeric(expression_df[rownames(expression_df) == ensg, 
                                                !(colnames(expression_df) %in% mutant_patients)])
  print(head(expression_normal))
  return(list("mutants" = expression_mutants, "normal" = expression_normal))
}



# Uniprot and ENSG IDs for driver gene of interest
tp53_uniprot <- "P04637" # TP53
tp53_ensg <- "ENSG00000141510"

pik3ca_uniprot <- "P42336"
pik3ca_ensg <- "ENSG00000121879"

# ENSG ID for target gene of interest
pfkfb3 <- "ENSG00000170525"

# Import expression DF
expression_df_tmm <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_tmmSDGr1filtByExpr_CancerOnly_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
expression_df_quantile_norm <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_quantile_norm_edgeRfilt_SDGr1_inverseNtransform_CancerOnly_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
rownames(expression_df_quantile_norm) <- expression_df_quantile_norm$ensg_id
expression_df_quantile_norm <- expression_df_quantile_norm[,3:ncol(expression_df_quantile_norm)]

# Import mutation DF
mutation_regprot_df <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Regprot_Mutation/iprotein_results_missense_nonsense_IntersectPatients.csv", sep = ""), 
                                header = TRUE, row.names = 1, check.names = FALSE)
mutation_regprot_df <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Regprot_Mutation/iprotein_results_nonsynonymous_IntersectPatients.csv", sep = ""), 
                                header = TRUE, row.names = 1, check.names = FALSE)

# Call function
create_driver_target_boxplot(expression_df_quantile_norm, mutation_regprot_df, tp53_ensg, tp53_uniprot, pfkfb3, "box", "TP53", "PFKFB3")
create_driver_target_boxplot(expression_df_quantile_norm, mutation_regprot_df,tp53_ensg, tp53_uniprot, pfkfb3, "violin", "TP53", "PFKFB3")


##############################################################################
##############################################################################
### FIGURE 3
##############################################################################
##############################################################################

# PART A: ENRICHMENT PLOT

#' Function to plot the enrichment of various target groups among the significant gene hits,
#' together on one plot, where the enrichment in the given source is denoted by color and line form
#' @param master_df a master DF produced from run_linear_model() that has q-values
#' @param target_sets_list a list of each of the target sets we are interested in doing enrichment on,
#' with each set's name/ label as the name for the list item
#' @param goi a string denoting the gene-of-interest, for labeling graph
#' @param thres a threshold for the number of hit genes we show in the full plot; in NA then we plot
#' enrichment across all gene targets
#' @param cutout_thres a threshold for the inset/cutout enrichment plot ('zoomed' in on)
#' @param silent_df OPT: an additional, corresponding DF with silent mutation results
plot_combined_enrichment <- function(master_df, target_sets_list, goi, thres, cutout_thres, silent_df) {
  
  # For each source, create a 0 and 1 vector for each  target hit to indicate if it 
  # is a a hit in that source
  target_set_vector_list <- lapply(target_sets_list, function(source) {
    vect <- unlist(lapply(1:nrow(master_df), function(i) {
      x <- master_df[i, 'T_k.name']
      return(ifelse(x %fin% source, 1, 0))
    }))
    return(vect)
  })
  
  target_set_vector_list_silent <- NA
  if(!is.na(silent_df)) {
    names(target_set_vector_list) <- paste0("MN.", names(target_sets_list))
    silent_df$T_k.name <- as.factor(silent_df$T_k.name)
    target_set_vector_list_silent <- lapply(target_sets_list, function(source) {
      vect <- unlist(lapply(1:nrow(silent_df), function(i) {
        x <- silent_df[i, 'T_k.name']
        return(ifelse(x %fin% source, 1, 0))
      }))
      return(vect)
    })
    names(target_set_vector_list_silent) <- paste0("S.", names(target_sets_list))
  } else {
    names(target_set_vector_list) <- names(target_sets_list)
  }
  
  #print(target_set_vector_list)
  #print(target_set_vector_list_silent)
  
  # Use helper function to get fraction of genes at given rank
  target_set_fraction_df <- get_frac_of_genes_at_rank(target_set_vector_list, master_df)
  target_set_fraction_df_silent <- NA
  if(!is.na(silent_df)) {
    target_set_fraction_df_silent <- get_frac_of_genes_at_rank(target_set_vector_list_silent,
                                                               silent_df)
    target_set_fraction_df <- merge(target_set_fraction_df, target_set_fraction_df_silent, by = "Rank")
  }
  
  
  target_sets_fraction_df_m <- melt(target_set_fraction_df, "Rank")
  colnames(target_sets_fraction_df_m) <- c("Rank", "Source", "Frac")
  target_sets_fraction_df_m$Mutation.Type <- unlist(lapply(target_sets_fraction_df_m$Source, function(x) {
    val <- "Nonsynonymous"
    if(grepl("S.", x, fixed = TRUE)) {val <- "Synonymous"}
    return(val)
  }))
  
  print(head(target_sets_fraction_df_m))
  print(unique(target_sets_fraction_df_m$Source))
  
  # Create inset plot
  p2 <- ggplot(target_sets_fraction_df_m[target_sets_fraction_df_m$Rank %in% 1:cutout_thres,], 
               mapping = aes(x = Rank, y = Frac, color = Source)) + 
    geom_line(aes(linetype=Mutation.Type), size = 1.25) +
    #geom_point(size = 1) + 
    scale_x_continuous(limits = c(1,cutout_thres)) + 
    #scale_color_nejm() + 
    #scale_alpha_manual(values=c(0.9,0.9,0.4,0.4)) +
    #scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#BC3C29FF", "#0072B5FF")) + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.minor = element_blank(), axis.text.x = element_text(size = 14, face = "bold"), 
          axis.text.y =  element_text(size = 14, face = "bold"))
  
  p <- ggplot(target_sets_fraction_df_m[target_sets_fraction_df_m$Rank %in% 1:thres,], 
              mapping = aes(x = Rank, y = Frac, color = Source)) + 
    geom_line(aes(linetype=Mutation.Type), size = 1.25) +
    #theme_minimal() +
    theme(axis.text = element_text(size = 16, face = "bold"), 
          axis.title = element_text(size = 18, face = "bold"), #panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), #panel.border = element_blank(),
          panel.background = element_rect(fill = 'white'), legend.text=element_text(size=14), 
          legend.title = element_text(size = 16)) +
    #geom_point(size = 1) + 
    #ggtitle(paste("Fraction of Model-Prioritized Target Genes in", paste(goi, "Sourced Targets"))) +
    xlab("Gene Rank") + ylab("Fraction in Set of Known Targets") + 
    scale_x_continuous(limits = c(1,thres)) +
    #scale_color_nejm() + 
    #scale_alpha_manual(values=c(0.9,0.9,0.4,0.4)) +
    #scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#BC3C29FF", "#0072B5FF"),
    #labels = c(names(target_sets_list), rep("", times = length(names(target_sets_list))))) + 
    annotation_custom(ggplotGrob(p2), xmin =(thres-round(thres/1.5)), xmax = thres, ymin = 0.30, ymax = 1.0) +
    geom_rect(data=target_sets_fraction_df_m[1:cutout_thres,], 
              aes(xmin = 1, xmax = cutout_thres, ymin = 0, ymax = 1, fill = "gray"),
              color = NA, fill = alpha("gray", .01))
  
  #vp <- viewport(width = 0.4, height = 0.4, x = 0.8, y = 0.2)
  #print(p)
  #print(p, vp = vp)
  
  return(p)
}


#' Get the fraction of genes at or above the given rank, for each set
#' @param target_set_vector_list a named list of each source with hit/not binary vector
#' @param master_df master df with ordered top hits
get_frac_of_genes_at_rank <- function(target_set_vector_list, master_df) {
  target_set_fraction_list <- lapply(target_set_vector_list, function(vect) {
    frac_nw_vect <- c()
    count_of_nw_genes <- 0
    for (i in 1:nrow(master_df)) {
      val_of_curr_tg <- vect[i]
      count_of_nw_genes <- count_of_nw_genes + val_of_curr_tg
      frac <- count_of_nw_genes / i
      frac_nw_vect <- c(frac_nw_vect, frac)
    }
    return(frac_nw_vect)
  })
  
  names(target_set_fraction_list) <- names(target_set_vector_list)
  target_sets_fraction_df <- as.data.frame(target_set_fraction_list)
  
  target_sets_fraction_df$Rank <- 1:nrow(master_df)
  
  return(target_sets_fraction_df)
}

# See create_enrichment_visualization.R for import lines for all validation datasets

# Create a list of all the gene sets we want to look for enrichment in
tp53_target_set_list <- list("ChIP-eat" = tp53_chipeat_targets, "Compiled.Cancer" = known_cancer_genes,
                             "STRING" = tp53_string_nw_targs, "HumanBase" = tp53_nw_targs, 
                             "Curated.Fischer.2017" = tp53_curated_targets)
tp53_target_set_list <- list("Compiled.Cancer" = known_cancer_genes,
                             "STRING" = tp53_string_nw_targs, "HumanBase" = tp53_nw_targs, 
                             "Curated.Fischer.2017" = tp53_curated_targets, 
                             "Curated.TRRUST" = tp53_trrust_targets, "Curated.TF-Target" = tf_target_tp53_full)

pik3ca_target_set_list <- list("Compiled.Cancer" = known_cancer_genes,
                               "STRING" = pik3ca_string_nw_targs, "HumanBase" = pik3ca_nw_targs,
                               "Curated.Cizkova.2017" = pik3ca_curated_targs)

# For grant
tp53_select_target_set_list <- list("Curated.Fischer.2017" = tp53_curated_targets, "KEGG.hsa04115" = tp53_kegg_pathway_genes,
                                    "TRRUST" = tp53_trrust_targets, "DoRothEA" = tp53_dorothea_targets)
tp53_select_target_set_list <- list("Curated.Fischer.2017" = tp53_curated_targets,
                                    "DoRothEA" = tp53_dorothea_targets)
tp53_select_target_set_list_cna <- list("STRING" = tp53_string_nw_targs,
                                        "TRRUST" = tp53_trrust_targets_upstr)
pik3ca_select_target_set_list <- list("STRING" = pik3ca_string_nw_targs_top500, "KEGG.hsa04151" = pik3ca_kegg_pathway_genes)


# Call function
plot_combined_enrichment(allgenes_p53, tp53_target_set_list, "TP53", 500, 50, NA)
plot_combined_enrichment(allgenes_pik3ca, pik3ca_target_set_list, "PIK3CA", 500, 50, NA)


# PART B: GENE SET ENRICHMENT ANALYSIS

# OPTION 1: GENE SET ENRICHMENT CURVES USING REACTOMEPA
#' Takes in a results table from the linear model output and plots GSEA
#' for n upregulated and n downregulated terms using the 
#' ReactomePA package
#' @param results_table a results data frame from my model with Beta estimates + pvalues
#' @param n the number of up/downregulated terms to plot
#' @param sort_by either "estimate" or "p.value" to indicate whether we want to 
#' sort our top hits by Beta estimate or p-value when performing GSEA
#' @param output_path a local path to save the GSEA figures to
perform_gsea <- function(results_table, n, sort_by, output_path) {
  
  regprot.gsea.rp <- list()
  regprot.gsea.ncg <- list()
  regprot.gsea.go <- list()
  regprot.gsea.kegg <- list()
  
  # Loop through all the regulatory proteins
  for (regprot in unique(results_table$R_i.name)) {
    
    # Get all of this protein's targets and their associated Betas, sorted in 
    # descending order
    res_table_sub <- results_table[results_table$R_i.name == regprot, c('T_k.name', 'estimate', 'p.value')]
    print(head(bitr(res_table_sub$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb=org.Hs.eg.db, drop = TRUE)))
    mapping <- as.data.frame(bitr(res_table_sub$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                                  OrgDb=org.Hs.eg.db, drop = TRUE))
    colnames(mapping) <- c("T_k.name", "T_k.entrez")
    res_table_sub <- merge(res_table_sub, mapping, all=TRUE, by="T_k.name")
    if(sort_by == "estimate") {
      res_table_sub <- res_table_sub[order(res_table_sub$estimate, decreasing = TRUE),]
    } else {
      res_table_sub$negLogPval <- unlist(lapply(1:nrow(res_table_sub), function(i) {
        pval <- res_table_sub$p.value[i]
        estimate <- res_table_sub$estimate[i]
        est_sign <- ifelse(estimate > 0, 1, -1)
        return((-log(pval)) * est_sign)
      }))
      res_table_sub <- res_table_sub[order(res_table_sub$negLogPval, decreasing = TRUE),]
    }
    #res_table_sub <- res_table_sub[!(is.na(res_table_sub))]
    #print(res_table_sub)
    
    plotdir <- paste0(output_path, "Linear Model/GSEA/")    
    suppressWarnings(dir.create(plotdir, recursive = TRUE))
    
    regprotBetaScores <- res_table_sub$estimate
    if(sort_by == "p.value") regprotBetaScores <- res_table_sub$negLogPval
    names(regprotBetaScores) <- res_table_sub$T_k.entrez
    
    gse.rp <- gsePathway(regprotBetaScores, pvalueCutoff = 1, pAdjustMethod = "BH")
    gse.rp <- setReadable(gse.rp, org.Hs.eg.db)
    gse.ncg <- gseNCG(regprotBetaScores, pvalueCutoff = 1, pAdjustMethod = "BH")
    gse.ncg <- setReadable(gse.ncg, org.Hs.eg.db)
    gse.go <- gseGO(regprotBetaScores, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, 
                    pAdjustMethod = "BH", keyType = "ENTREZID")
    gse.go <- setReadable(gse.go, org.Hs.eg.db)
    #gse.kegg <- gseKEGG(regprotBetaScores, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, 
    #pAdjustMethod = "BH", keyType = "ENTREZID")
    #gse.kegg <- setReadable(gse.kegg, org.Hs.eg.db)
    
    regprot.gsea.rp[[regprot]] <- gse.rp
    regprot.gsea.ncg[[regprot]] <- gse.ncg
    regprot.gsea.go[[regprot]] <- gse.go
    #regprot.gsea.kegg[[regprot]] <- gse.kegg
    res_tab <- data.table()
    
    
    # plot GSEA for top N upregulated and top N downregulated terms
    for (gst in c('rp','ncg', 'go')) {
      #for (gst in c('rp','ncg')) {
      if (gst == 'rp') res <- gse.rp
      #else if (gst == 'kegg') res <- gse.kegg
      else if (gst == 'go') res <- gse.go
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
      if(gst == "go") {res_tab <- res}
    }
  }
  return(res_tab)
}

results_gsea <- perform_gsea(master_df_mut_corrected, n = 5, "p.value", output_path)
results_gsea <-perform_gsea(master_df_cna_corrected, n = 5, "p.value", output_path)

results_gsea <- results_gsea@result


# OPTION 2: BAR CHART REPORTING THE ENRICHMENT FOR TOP PATHWAYS
#' Takes the output from the ReactomePA gene set enrichment analysis (above) and
#' converts it into a bar chart format for visualization, with -log(p-value) on
#' the x-axis
#' @param results_gsea a results table from GSEA using the ReactomePA package
#' @param n the number of top pathways to display
#' @param db the name of the database from which the pathways originated (e.g. "GO")
create_gsea_barchart <- function(results_gsea, n, db) {
  # Subset the results table based on the number of pathways
  results_gsea_sub <- results_gsea[1:n, ]
  
  # Convert the p-values to -log10(pvalue)
  negLog10_pvalues <- -log10(results_gsea_sub$p.adjust)
  
  # Create an input data frame for ggplot
  input_df <- data.frame("Enriched.Pathway" = results_gsea_sub$Description,
                         "Neg.Log10.Pval" = negLog10_pvalues)
  
  #color_palette <- rep(pal_nejm("default")(8), times = ceiling(n/8))
  
  p <- ggplot(input_df, aes(x = reorder(Enriched.Pathway, negLog10_pvalues), 
                            y = Neg.Log10.Pval)) +
    geom_col(width = 0.7, fill = "#0072B5FF") + coord_flip() + theme_minimal() + 
    theme(legend.position = "none") +  ylab("-log10(adj.pvalue)") + 
    xlab(paste("Enriched", paste(db, "Pathway")))
    #+ scale_fill_manual(values = color_palette)
    #+ scale_fill_nejm() 
  
  p
}

# Call function
create_gsea_barchart(results_gsea, 25, "GO")


# OPTION 3: GENE SET ENRICHMENT PATHWAY MAP (TODO)


# PART C: COMPARING REAL TO RANDOMIZED VERSIONS OF MODEL

# OPTION 1: DENSITY PLOT OF BETA DISTRIBUTIONS (H0 vs. HA)
#' Plot a multi-layer histogram showing all the distribution of Beta values for 
#' a particular GOI, both "real" and "randomized" versions (H0 vs. HA) with
#' significant hits (defined as being below a given q-value threshold) highlighted
#' @param real_master_df a master DF for the "real" run (HA)
#' @param random_master_df a master DF for the "randomized" run (H0)
#' @param goi a driver gene-of-interest
#' @param qval_thres a q-value threshold for significance
plot_beta_densities_real_and_random <- function(real_master_df, random_master_df, goi, qval_thres) {
  
  betas_real <- real_master_df[real_master_df$R_i.name == goi, c('estimate', 'T_k.name')]
  betas_random <- random_master_df[random_master_df$R_i.name == goi, c('estimate', 'T_k.name')]

  betas <- merge(betas_real, betas_random, by = "T_k.name")
  colnames(betas)[colnames(betas) == "estimate.x"] <- "Real"
  colnames(betas)[colnames(betas) == "estimate.y"] <- "Random"
  betas_m <- melt(betas)

  
  # Use a two-sided Smirnov test to test if there are differences between these two distributions 
  ks <- ks.test(betas_real$estimate, betas_random$estimate, alternative = "two.sided")
  pval <- ks$p.value
  print(paste("K-S Test P-Value:", pval))
  
  my_comparisons <- c("Real", "Random")
  
  p <- ggplot(betas_m, aes(x=value, fill = variable)) + geom_density(alpha = 0.2) +
    scale_color_nejm() + theme_minimal() + xlab("Beta estimate") + ylab("Density") + 
    labs(fill = "Real or Randomized") + annotate("text", x = 0.5, y = 4, label = paste("K-S Test, P-Value:", pval))
    #stat_compare_means(comparisons = my_comparisons)
    #geom_vline(xintercept=means_of_betas, size=0.25, color="black", linetype = "dotted", alpha = 0.2) # add the mean of the betas for each
  
  p
}

# Call function
plot_beta_densities_real_and_random(allgenes_p53, allgenes_p53_random, "TP53", 0.1)


# OPTION 2: OVERLAID P-VALUE HISTOGRAMS OF REAL AND RANDOM
#' Plot a multi-layer histogram showing the p-value distributions for 
#' a particular GOI, both "real" and "randomized" versions (H0 vs. HA) with
#' significant hits (defined as being below a given q-value threshold) highlighted
#' @param real_master_df a master DF for the "real" run (HA)
#' @param random_master_df a master DF for the "randomized" run (H0)
#' @param goi a driver gene-of-interest
#' @param qval_thres a q-value threshold for significance
plot_pvalue_histograms_real_and_random <- function(real_master_df, random_master_df, 
                                                   goi, qval_thres) {
  
  pvals_real <- real_master_df[real_master_df$R_i.name == goi, c('p.value', 'T_k.name')]
  pvals_random <- random_master_df[random_master_df$R_i.name == goi, c('p.value', 'T_k.name')]
  
  pvals <- merge(pvals_real, pvals_random, by = "T_k.name")
  colnames(pvals)[colnames(pvals) == "p.value.x"] <- "Real"
  colnames(pvals)[colnames(pvals) == "p.value.y"] <- "Random"
  pvals_m <- melt(pvals)
  
  print(head(pvals_m))
  
  # Get the p-value corresponding to q-value threshold
  max_qval_under_thres <- max(real_master_df[real_master_df$q.value < qval_thres, 'q.value'])
  print(max_qval_under_thres)
  pval_thres <- real_master_df[real_master_df$q.value == max_qval_under_thres, 'p.value']
  
  p <- ggplot(pvals_m, aes(x = value, fill = variable)) + scale_fill_nejm() +
    geom_histogram(alpha = 0.3, position = "identity", bins = 50) + 
    theme_minimal() + xlab("P-Value") + ylab("Frequency") + 
    labs(fill = "Real or Randomized") + theme(axis.title = element_text(face = "bold", size = 14),
                                              axis.text = element_text(face = "bold", size = 12),
                                              legend.text = element_text(size=12),
                                              legend.title = element_text(face = "bold", size = 14)) +
    geom_vline(xintercept = pval_thres, linetype='dotted')+
    annotate("text", y = 1000, x = pval_thres, 
             label = paste("q < 0.1; p <", round(pval_thres, digits = 4)), 
             hjust = -0.5, size = 5, fontface = "bold")
  p
}

# Call function
plot_pvalue_histograms_real_and_random(allgenes_p53, allgenes_p53_random, "TP53", 0.1)


# PART D: PLOTTING OVERLAP WITH METABRIC RESULTS
#' Creates a Venn diagram with the overlapping significant hits when model is 
#' run on the same GOI using two independent data sources
#' @param master_df_list a list of the master DFs to compare, with names being
#' labels for the data source
#' @param goi the name of the gene-of-interest, for labeling purposes
#' @param qval_thres a q-value threshold for significance
plot_hit_overlap_between_data_sources <- function(master_df_list, goi, qval_thres) {
  
  # Get the significant hits for each
  sig_hits_list <- lapply(master_df_list, function(df) df[df$q.value < qval_thres, 'T_k.name'])
  names(sig_hits_list) <- names(master_df_list)
  #print(head(sig_hits_list))
  
  # Plot Venn Diagram
  ggVennDiagram(sig_hits_list, label_alpha = 0, category.names = names(sig_hits_list),
                set_size = 9, label_size = 9, edge_size = 0, 
                label = "count", set_color = rep("black", times = length(sig_hits_list))) +
    scale_fill_gradient(low="#FFDC91FF", high = "#0072B5FF") + labs(fill = paste("# Hits, q<", qval_thres)) +
    scale_color_manual(values = rep("black", times = length(sig_hits_list)))
}

# Import relevant results
allgenes_p53_metabric <- read.csv("C:/Users/sarae/Documents/res_P53_allGenes_iprotein_log2intensity_CNAbucket_justDel_methlog2Raw_cibersortTotalFrac_rmCis_PIK3CA_covs.inclCT.MN_uncorrected_MUT.csv",
                                  header = TRUE, check.names = FALSE)

master_df_list <- list("TCGA" = allgenes_p53, "METABRIC" = allgenes_p53_metabric)

# Call function
plot_hit_overlap_between_data_sources(master_df_list, "TP53", 0.1)























##############################################################################
##############################################################################
### ARCHIVES
##############################################################################
##############################################################################

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
  colnames(mut_count_matrix)[2:ncol(mut_count_matrix)] <- unlist(lapply(colnames(mut_count_matrix)[2:ncol(mut_count_matrix)], function(x) {
    return(unlist(strsplit(x, "-", fixed = TRUE))[1])
  }))
  mut_count_matrix_sub <- mut_count_matrix[,c(1, (which(colnames(mut_count_matrix)[2:ncol(mut_count_matrix)] %in% patient_set) + 1))]

  # Get the number of tumor samples with a mutation in each gene (count >= 1)
  mut_gene1_counts <- as.numeric(unlist(mut_count_matrix_sub[mut_count_matrix_sub$Gene_Symbol == gene1, 2:ncol(mut_count_matrix_sub)]))
  gene1_mut_samps <- which(mut_gene1_counts >= 1)
  
  mut_gene2_counts <- as.numeric(unlist(mut_count_matrix_sub[mut_count_matrix_sub$Gene_Symbol == gene2, 2:ncol(mut_count_matrix_sub)]))
  gene2_mut_samps <- which(mut_gene2_counts >= 1)
  
  mut_samps <- list("gene1" = gene1_mut_samps, "gene2" = gene2_mut_samps, "all" = 1:ncol(mut_count_matrix_sub))

  # Make the Venn Diagram
  ggVennDiagram(mut_samps, label_alpha = 0, category.names = c(gene1, gene2, "All"), set_color = "black",
                set_size = 10, label_size = 8, edge_size = 0) +
    ggplot2::scale_fill_gradient(low="cornsilk1", high = "cadetblue3")
}

patient_set <- read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt"), header = TRUE)[,1]
patient_set_lt20pamp <- intersect(read.table(paste0(main_path, "Patient Subsets/LessThan20PercAmp_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumA <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.A_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumB <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.B_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumAB <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.A.B_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_basal <- intersect(read.table(paste0(main_path, "Patient Subsets/Basal_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_her2 <- intersect(read.table(paste0(main_path, "Patient Subsets/HER2_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_normLike <- intersect(read.table(paste0(main_path, "Patient Subsets/Normal-like_patient_ids.txt"), header = TRUE)[,1], patient_set)

patient_set_pc <- read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt"), header = TRUE)[,1]

blca_subtype <- TCGAquery_subtype(tumor = "BLCA")
hnsc_subtype <- TCGAquery_subtype(tumor = "HNSC")
blca_subtype$justPat <- unlist(lapply(blca_subtype$patient, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))
hnsc_subtype$patient <- as.character(hnsc_subtype$patient)
hnsc_subtype$justPat <- unlist(lapply(hnsc_subtype$patient, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))

patient_set_blca <- intersect(patient_set_pc, blca_subtype$justPat)
patient_set_hnsc <- intersect(patient_set_pc, hnsc_subtype$justPat)

mutation_count_matrix <- read.csv(paste0(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_CancerOnly_IntersectPatients.csv"), 
                                  header = TRUE, row.names = 1, check.names = FALSE)
mutation_count_matrix <- read.csv(paste0(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_nonsense_CancerOnly_IntersectPatients.csv"), 
                                  header = TRUE, row.names = 1, check.names = FALSE)
#mutation_count_matrix <- read.csv(paste0(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_ALL.csv"), 
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
  
  # Use helper function to create input matrix from results table
  matrix <- create_regprot_v_genetarg_matrix(results_table)
  
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


#' Helper function to create regprot vs. gene target matrix and fill with Beta values 
#' @param results_table
create_regprot_v_genetarg_matrix <- function(results_table) {
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
### COMBINED SUBTYPE MODEL RESULT HEAT MAPS
##############################################################################
#' Create a partitioned heat map with the results from all breast cancer subtypes,
#' partitioned by subtype and set to a uniform scale
#' @param master_df_lumA luminal A results DF
#' @param master_df_lumB luminal B results DF
#' @param master_df_basal basal (TN) results DF
#' @param master_df_her2 HER2 results DF
#' @param signif_only a TRUE/ FALSE value indicating whether or not to limit
#' heatmap to only significant hits (q < 0.2)
create_subtype_partitioned_heatmap <- function(master_df_lumA, master_df_lumB,
                                               master_df_basal, master_df_her2, 
                                               signif_only) {
  
  # Limit to significant hits only, if desired
  if(signif_only) {
    master_df_lumA <- master_df_lumA[master_df_lumA$q.value < 0.2,]
    master_df_lumB <- master_df_lumB[master_df_lumB$q.value < 0.2,]
    master_df_basal <- master_df_basal[master_df_basal$q.value < 0.2,]
    master_df_her2 <- master_df_her2[master_df_her2$q.value < 0.2,]
  }
  
  # User helper function to create an input matrix for each master DF
  lumA_incl <- FALSE
  if((nrow(master_df_lumA) > 0) & (length(unique(master_df_lumA$R_i.name)) >= 2)) {
    matrix_lumA <- create_regprot_v_genetarg_matrix(master_df_lumA)
    matrix_lumA <- cbind(matrix_lumA, rep(1, times = nrow(matrix_lumA)))
    colnames(matrix_lumA)[ncol(matrix_lumA)] <- "subtype"
    lumA_incl <- TRUE
  }
  
  lumB_incl <- FALSE
  if((nrow(master_df_lumB) > 0) & (length(unique(master_df_lumB$R_i.name)) >= 2)) {
    matrix_lumB <- create_regprot_v_genetarg_matrix(master_df_lumB)
    matrix_lumB <- cbind(matrix_lumB, rep(2, times = nrow(matrix_lumB)))
    colnames(matrix_lumB)[ncol(matrix_lumB)] <- "subtype"
    lumB_incl <- TRUE
  }
  
  basal_incl <- FALSE
  if((nrow(master_df_basal) > 0) & (length(unique(master_df_basal$R_i.name)) >= 2)) {
    matrix_basal <- create_regprot_v_genetarg_matrix(master_df_basal)
    matrix_basal <- cbind(matrix_basal, rep(3, times = nrow(matrix_basal)))
    colnames(matrix_basal)[ncol(matrix_basal)] <- "subtype"
    basal_incl <- TRUE
  }
  
  her2_incl <- FALSE
  if((nrow(master_df_her2) > 0) & (length(unique(master_df_her2$R_i.name)) >= 2)) {
    matrix_her2 <- create_regprot_v_genetarg_matrix(master_df_her2)
    matrix_her2 <- cbind(matrix_her2, rep(4, times = nrow(matrix_her2)))
    colnames(matrix_her2)[ncol(matrix_her2)] <- "subtype"
    her2_incl <- TRUE
  }
  
  
  # Combine all these into one
  matrices <- list(matrix_lumA, matrix_lumB, matrix_basal, matrix_her2)
  incl_vect <- c(lumA_incl, lumB_incl, basal_incl, her2_incl)
  incl_vect_num <- unlist(lapply(1:length(incl_vect), function(i) ifelse(incl_vect[i] == TRUE, i, NA)))
  incl_vect_num <- incl_vect_num[!is.na(incl_vect_num)]
  matrices_to_keep <- matrices[incl_vect_num]
  print(matrices_to_keep)
  master_matrix <- do.call(rbind, matrices_to_keep)
  print(head(master_matrix))
  print(unique(master_matrix[,'subtype']))
  #master_matrix <- rbind(rbind(rbind(matrix_lumA, matrix_lumB), matrix_basal), matrix_her2)
  
  # Get row clustering
  #row_dend <- hclust(dist(master_matrix[,1:2]))
  
  hm_annot <- HeatmapAnnotation(df = data.frame(subtype = c(rep("Luminal A", nrow(matrix_lumA)), 
                                                            rep("Luminal B", nrow(matrix_lumB)),
                                                            rep("Basal", nrow(matrix_basal)),
                                                            rep("HER2", nrow(matrix_her2)))),
                                col = list(subtype = c("Luminal A" = "pink", "Luminal B" = "lightblue",
                                                       "Basal" = "orange", "HER2" = "lightgreen")),
                                show_legend = TRUE, annotation_name_side = "right")
  draw(hm_annot)
  
  # Create partitioned heatmap
  Heatmap(master_matrix[, 1:2], name = "Beta", column_title = "Driver Gene", column_title_side = "bottom",
          row_title = "Target Gene", split = master_matrix[, 'subtype'], row_gap = unit(2, "mm"), 
          row_names_gp = gpar(fontsize = 6), border = c("black"), cluster_columns = FALSE, 
          column_names_rot = 0, show_row_names = FALSE)
          #rowAnnotation = hm_annot, show_annotation_legend = TRUE)
  
}

create_subtype_partitioned_heatmap(res_metabol_lumA_mut, res_metabol_lumB_mut, res_metabol_basal_mut,
                                   res_metabol_her2_mut, TRUE)

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


