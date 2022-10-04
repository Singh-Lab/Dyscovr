################################################################
### ANALYZE LI ET AL CELL LINE METABOLOMIC DATA ###
### Written by: Sara Geraghty, Feb. 2022
### Link to Paper: https://www.nature.com/articles/s41591-019-0404-8#Sec31
################################################################

# Analyzes the data from Li et al. (2019), Nature Medicine, from
# Supplementary tables 1-4.

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/CL Metabolic Data (Li et al. 2019)/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


# Import the supplementary data table information

# Cell line annotations
cell_line_annot <- read.csv(paste0(path, "cell_line_annotations.csv"), 
                            header = TRUE, check.names = FALSE)
cell_line_annot_brca <- cell_line_annot[(grepl("breast", cell_line_annot$Classifications, ignore.case = TRUE)) & 
                                           (grepl("female", cell_line_annot$Gender)),]
brca_cell_lines <- cell_line_annot_brca$Name

# Clean metabolite data from LCMS
metabolite_df <- read.csv(paste0(path, "clean_metabolite_data.csv"), 
                          header = TRUE, check.names = FALSE, row.names = 1)
metabolite_df_cl_names <- unlist(lapply(rownames(metabolite_df), function(x) 
  unlist(strsplit(x, "_", fixed = TRUE))[1]))
metabolite_df_brca <- metabolite_df[which(metabolite_df_cl_names %in% brca_cell_lines),]

# Mutation, CNA, & methylation correlations to metabolite levels
mut_metabol_corr <- read.csv(paste0(path, "mutation_metabolite_correlations.csv"), 
                             header = TRUE, check.names = FALSE)
cna_metabol_corr <- read.csv(paste0(path, "cna_metabolite_correlations.csv"), 
                             header = TRUE, check.names = FALSE)
methylation_metabol_corr <- read.csv(paste0(path, "methylation_metabolite_correlations.csv"), 
                             header = TRUE, check.names = FALSE)


################################################################
# INVESTIGATE MUTATION CORRELATIONS FOR PARTICULAR GOIS
################################################################
#' For a particular given gene of interest, find metabolic targets
#' that are predicted by their model to have mutation correlation
#' t-statistics above or below some given threshold
#' @param mutation_corr a mutation-metabolite correlation DF
#' @param goi a gene whose mutation status is of interest
#' @param t_thres a t-statistic threshold for significance
#' @param outpath a path to write output files to 
get_mutation_metabol_corr_for_goi <- function(mutation_corr, goi, t_thres, outpath) {
  # Subset the mutation-correlation data frame to this GOI
  mutation_corr_goi <-  mutation_corr[grepl(goi, mutation_corr$genes),]
  
  # Get the names of the metabolites that exceed the given t-statistic
  # threshold in either direction (look both at all mutations and deleterious
  # mutations and report both)
  mutation_corr_goi_allMut <- mutation_corr_goi[grepl("all_mutations", mutation_corr_goi$genes),]
  mutation_corr_goi_delMut <- mutation_corr_goi[grepl("deleterious_mutations", mutation_corr_goi$genes),]
  
  allMut_tvals <- as.numeric(unlist(mutation_corr_goi_allMut[,3:ncol(mutation_corr_goi_allMut)]))
  delMut_tvals <- as.numeric(unlist(mutation_corr_goi_delMut[,3:ncol(mutation_corr_goi_delMut)]))
  
  print(head(allMut_tvals))
  
  allMut_sig_tvals <- as.numeric(unlist(which(abs(allMut_tvals) > t_thres))) + 2
  delMut_sig_tvals <- as.numeric(unlist(which(abs(delMut_tvals) > t_thres))) + 2
  
  print(head(allMut_sig_tvals))
  
  # Report the names of the significant metabolites
  allMut_sig_metabolites <- colnames(mutation_corr_goi)[allMut_sig_tvals]
  delMut_sig_metabolites <- colnames(mutation_corr_goi)[delMut_sig_tvals]
  
  # Print these results to console
  print(paste(goi, paste("mutation (all mutation) correlated metabolites with t-statistic threshold of", 
                         paste(t_thres, ":"))))
  print(allMut_sig_metabolites)
  print(paste(goi, paste("mutation (deleterious mutation) correlated metabolites with t-statistic threshold of", 
                         paste(t_thres, ":"))))
  print(delMut_sig_metabolites)
  
  # Write these results to a file
  allMut_df <- data.frame("Metabolite" = allMut_sig_metabolites, "t-statistic" = unlist(mutation_corr_goi_allMut[,allMut_sig_tvals]))
  write.csv(allMut_df, paste(outpath, paste(goi, "all_mutation_top_metabolic_correlations.csv", sep = "_")))
  delMut_df <- data.frame("Metabolite" = delMut_sig_metabolites, "t-statistic" = unlist(mutation_corr_goi_delMut[,delMut_sig_tvals]))
  write.csv(delMut_df, paste(outpath, paste(goi, "deleterious_mutation_top_metabolic_correlations.csv", sep = "_")))
  
  return(allMut_df)
}

tp53_mut_df <- get_mutation_metabol_corr_for_goi(mut_metabol_corr, "TP53", 2.0, path)
pik3ca_mut_df <- get_mutation_metabol_corr_for_goi(mut_metabol_corr, "PIK3CA", 2.0, path)

################################################################
# INVESTIGATE CNA CORRELATIONS FOR PARTICULAR GOIS
################################################################
#' For a particular given gene of interest, find metabolic targets
#' that are predicted by their model to have CNA correlation
#' t-statistics above or below some given threshold
#' @param cna_corr a CNA-metabolite correlation DF
#' @param goi a gene whose mutation status is of interest
#' @param t_thres a t-statistic threshold for significance
#' @param outpath a path to write output files to 
get_cna_metabol_corr_for_goi <- function(cna_corr, goi, t_thres, outpath) {
  # Subset the CNA-correlation data frame to this GOI
  cna_corr_goi <-  cna_corr[grepl(goi, cna_corr$genes),]
  
  # Get the names of the metabolites that exceed the given t-statistic
  # threshold in either direction 
  cna_tvals <- as.numeric(unlist(cna_corr_goi[,3:ncol(cna_corr_goi)]))
  print(head(cna_tvals))
  cna_sig_tvals <- as.numeric(unlist(which(abs(cna_tvals) > t_thres))) + 2
  print(head(cna_sig_tvals))
  
  # Report the names of the significant metabolites
  cna_sig_metabolites <- colnames(cna_corr_goi)[cna_sig_tvals]

  # Print these results
  print(paste(goi, paste("CNA correlated metabolites with t-statistic threshold of", 
                         paste(t_thres, ":"))))
  print(cna_sig_metabolites)
  
  # Write these results to a file
  cna_df <- data.frame("Metabolite" = cna_sig_metabolites, "t-statistic" = unlist(cna_corr_goi[,cna_sig_tvals]))
  write.csv(cna_df, paste(outpath, paste(goi, "cna_top_metabolic_correlations.csv", sep = "_")))
}

get_cna_metabol_corr_for_goi(cna_metabol_corr, "TP53", 1.0, path)  # not part of these results
get_cna_metabol_corr_for_goi(cna_metabol_corr, "PIK3CA", 1.0, path) # not part of these results


################################################################
# RELATE METABOLITES TO METABOLIC ENZYMES
################################################################
# Use output from Recon3D metabolic network to relate top metabolites to genes of interest
# Output file generated by Antonio, see "gene-to-metabolite-mapping" folder

gene_to_metabol_mapping <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/gene-to-metabolite-mapping/output.tsv", 
                                    header = TRUE, check.names = FALSE, sep = "\t")

# Also get a mapping between BIGG metabolite IDs and common metabolite names, downloaded from: http://bigg.ucsd.edu/data_access
metabolite_id_dict <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/gene-to-metabolite-mapping/bigg_models_metabolites.tsv", 
                               header = TRUE, check.names = FALSE, sep = "\t")
# Add the common names to the mapping
gene_to_metabol_mapping$metabolite_list_names <- unlist(lapply(gene_to_metabol_mapping$metabolite_list, function(metabols){
  spl_metabols <- unique(unlist(strsplit(metabols, split = ",", fixed = TRUE)))
  names <- unlist(lapply(spl_metabols, function(x) metabolite_id_dict[metabolite_id_dict$bigg_id == x, 'name']))
  return(paste(names, collapse = ","))
}))

#' Gets and returns a vector of ENSG IDs associated with a particular given
#' metabolite of interest
#' @param gene_to_metabol_mapping a two-column mapping between ENSG IDs and associated
#' metabolites produced from Recon3D and Antonio's code
#' @param metabolite the name of a metabolite of interest
get_prots_assoc_w_metabol <- function(gene_to_metabol_mapping, metabolite) {
  return(paste(unlist(gene_to_metabol_mapping[grepl(metabolite, gene_to_metabol_mapping$metabolite_list_names, ignore.case = TRUE), 'ensembl_gene_id']),
               collapse = ","))
}

# Use this function to add the assoc. proteins of interest for each significant metabolite
# to the table(s) above
tp53_mut_df$Assoc.Genes <- lapply(tp53_mut_df$Metabolite, function(metabol) 
  return(get_prots_assoc_w_metabol(gene_to_metabol_mapping, metabol)))
pik3ca_mut_df$Assoc.Genes <- lapply(pik3ca_mut_df$Metabolite, function(metabol) 
  return(get_prots_assoc_w_metabol(gene_to_metabol_mapping, metabol)))

# Add another column with associated gene names
tp53_mut_df$Assoc.Gene.Names <- unlist(lapply(tp53_mut_df$Assoc.Genes, function(x) {
  ids <- unlist(strsplit(x, ",", fixed = TRUE))
  names <- lapply(ids, function(id) {
    return(paste(unique(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == id, 
                                   'external_gene_name']), collapse = ";"))
  })
  names_char <- paste(names, collapse = ",")
  return(names_char)
}))
pik3ca_mut_df$Assoc.Gene.Names <- unlist(lapply(pik3ca_mut_df$Assoc.Genes, function(x) {
  ids <- unlist(strsplit(x, ",", fixed = TRUE))
  names <- lapply(ids, function(id) {
    return(paste(unique(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == id, 
                                          'external_gene_name']), collapse = ";"))
  })
  names_char <- paste(names, collapse = ",")
  return(names_char)
}))  

# Add a column that, for each protein name, annotates whether it's a significant hit
# from the model
tp53_signif_hits <- unlist(master_df[master_df$q.value < 0.2, 'T_k.name'])
tp53_mut_df$Signif.Hit <- unlist(lapply(tp53_mut_df$Assoc.Gene.Names, function(x) {
  names <- unlist(strsplit(x, ",", fixed = TRUE))
  in_signif_hits <- unlist(lapply(names, function(n) ifelse(n %in% tp53_signif_hits, TRUE, FALSE)))
  return(paste(in_signif_hits, collapse = ","))
}))
pik3ca_signif_hits <- unlist(master_df[master_df$q.value < 0.2, 'T_k.name'])
pik3ca_mut_df$Signif.Hit <- unlist(lapply(pik3ca_mut_df$Assoc.Gene.Names, function(x) {
  names <- unlist(strsplit(x, ",", fixed = TRUE))
  in_signif_hits <- unlist(lapply(names, function(n) ifelse(n %in% pik3ca_signif_hits, TRUE, FALSE)))
  return(paste(in_signif_hits, collapse = ","))
}))

print(intersect(tp53_signif_hits, unlist(strsplit(as.character(unlist(tp53_mut_df$Assoc.Gene.Names)), ",", fixed = TRUE))))
print(intersect(pik3ca_signif_hits, unlist(strsplit(as.character(unlist(pik3ca_mut_df$Assoc.Gene.Names)), ",", fixed = TRUE))))

#' Identify the significant correlations and print them
#' @param mut_df a metabolism DF for a particular protein of interest 
print_signif_corr <- function(mut_df) {
  for(i in 1:nrow(mut_df)) {
    if(grepl("TRUE", mut_df$Signif.Hit[i])) {
      true_false_vect <- as.integer(as.logical(unlist(strsplit(mut_df$Signif.Hit[i], ",", fixed = TRUE))))
      str <- "Mutation correlated with"
      estimate <- mut_df$t.statistic[i]
      if(estimate > 0) {str <- paste(str, "upregulation of")}
      else {str <- paste(str, "downregulation of")}
      str <- paste(str, mut_df$Metabolite[i])
      str <- paste0(str, paste(", with t-statistic of", estimate))
      str <- paste0(str, ". This metabolite is associated with expression of")
      genes <- unlist(strsplit(mut_df$Assoc.Gene.Names[i], ",", fixed = TRUE))[which(true_false_vect == 1)]
      genes <- paste(genes, collapse = ",")
      str <- paste(str, genes)
      str <- paste0(str, ", which is(are) significant hits from our model.")
      print(str)
    }
  }
}

print_signif_corr(tp53_mut_df)
print_signif_corr(pik3ca_mut_df)


# Check length of overlap between metabolic gene list inputs and genes from the paper
metabolism_gene_list_ensg <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/metabolism_gene_list_ensg.txt", header = TRUE)[,1]
length(intersect(metabolism_gene_list_ensg, unlist(strsplit(as.character(unlist(tp53_mut_df$Assoc.Genes)), ",", fixed = TRUE))))  
# For TP53; 141 genes of 166 genes from paper; 1837 metabolic genes from Recon3D (there are 25 genes here that were not tested in the metabolism model)


################################################################
# PLOT ENRICHMENT IN GENES RELATED TO TOP HIT METABOLITES
################################################################
#' For the top hits of a given protein-of-interest, plot an 
#' enrichment curve for whether these hits are associated with 
#' a significantly dysregulated metabolite from this paper
#' @param master_df a master DF produced from our model, subsetted
#' to the given protein of interest
#' @param mut_df genes from a metabolism DF for a particular protein of interest 
#' @param thres a qvalue threshold for significance
plot_target_assoc_w_metabol_enrichment <- function(master_df, mut_df_genes, thres) {
  master_df_sig <- master_df[master_df$q.value < thres,]
  
  sig_hits <- unique(unlist(master_df_sig$T_k.name))

  rank <- 1:length(sig_hits)
  
  cumulat_val <- 0
  enrich_vals <- c()
  
  for(i in 1:length(sig_hits)) {
    sig_hit <- sig_hits[i]
    if(sig_hit %in% mut_df_genes) {
        cumulat_val <- cumulat_val + 1
    }
    frac <- cumulat_val / i
    print(frac)
    enrich_vals <- c(enrich_vals, frac)
  }
  
  plot(rank, enrich_vals, main = "Enrichment of Significant Target Genes in P53-Assoc. Metabolic Regulators from Li. et al. (2019)",
       xlab = "Rank", ylab = "Fraction of Target Genes that are Signif. Metabolic Reg.", pch = 19)
  
}

tp53_mut_df_genes <- unique(unlist(strsplit(as.character(unlist(tp53_mut_df$Assoc.Genes)), ",", fixed = TRUE)))
tp53_mut_df_genes_overlap <- tp53_mut_df_genes[tp53_mut_df_genes %in% metabolism_gene_list_ensg]
tp53_mut_df_genes_overlap_names <- unique(unlist(lapply(tp53_mut_df_genes_overlap, function(x) 
  all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == x, 'external_gene_name'])))


plot_target_assoc_w_metabol_enrichment(master_df_mut, tp53_mut_df_genes_overlap_names, 0.2)


