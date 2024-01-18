##############################################################################
# CREATE FIGURES FOR METABOLISM PAPER
# Created By: Sara Geraghty, Dec. 2021
##############################################################################

# Compiles the code for making these figures into one file
# All other plots were made in BioRender.com

library(ggplot2)
library(TCGAbiolinks)
library(gplots)
library(stats)
library(RColorBrewer)
library(ggVennDiagram)
#library("gage")
library(dendextend)
library("pheatmap")
#library(ComplexHeatmap)
library(data.table)
library('ReactomePA')
library('clusterProfiler')
library("DOSE")
library(GOSemSim)
library("enrichplot")
library(ggrepel)
library(ggsci)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(UpSetR)
library(reshape2)
library(ggbeeswarm)
library(survival)
library(survminer)
library(eulerr)
library(DescTools)
library(lsa)
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
               showPct = TRUE, color = colors)
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

# Import intersecting patients
intersecting_patients <- read.table(paste(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt", sep = ""))[,1]

# Subset to only patients of interest
maf_file_barcodes <- as.character(unlist(maf@data$Tumor_Sample_Barcode))
patient_ids <- unlist(lapply(maf_file_barcodes, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))
intersecting_patient_barcodes <- unique(maf_file_barcodes[which(patient_ids %in% intersecting_patients)])
maf_sub <- subsetMaf(maf, tsb = intersecting_patient_barcodes)

# Get clinical DF
clinical_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/clinical_data_subset.csv", 
                        header = TRUE, check.names = FALSE)

# Call function
maf_list_by_subtype <- create_top_mutated_drivers_by_cancer_type_barplot(maf_sub, driver_gene_df, 5, 
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
relate_sample_size_to_sig_eQTL <- function(list_of_master_dfs, sample_size_dict, 
                                           list_of_maf_files, goi, qval_thres) {
  
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
  
  # OPTION: use sample size to convert this to a fraction (% mutated)
  total.mut.count$frac.mutated <- unlist(lapply(1:nrow(total.mut.count), function(i) {
    type <- total.mut.count[i, 'type']
    samp.size <- input_df[input_df$type == type, 'sample.size']
    frac.mut <- total.mut.count[i, 'total.mut.count'] / samp.size
    return(frac.mut)
  }))
  
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
    geom_point(aes(size = frac.mutated)) + scale_color_nejm() +
    theme_minimal() + xlab("Sample Size") + ylab(paste("Num. Significant Correlations (q <", paste0(qval_thres, ")"))) +
    labs(size = paste("Frac. with Nonsyn. Mut in", goi), color = "BRCA Subtype") + 
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

subtypes_of_interest <- c("LumA", "LumB", "Her2", "Basal", "Normal")
list_of_sample_sizes <- lapply(subtypes_of_interest, function(st) 
  nrow(brca_subtype_file_sub[brca_subtype_file_sub$BRCA_Subtype_PAM50 == st,]))
names(list_of_sample_sizes) <- subtypes_of_interest

# Use the generated list of MAFs subsetted by subtype from above (part B)

# Call function
relate_sample_size_to_sig_eQTL(list_of_master_dfs, list_of_sample_sizes, maf_list_by_subtype, "TP53", 0.1)



#' Takes a list of result master data frames under consideration, along with a dictionary
#' with the size of each corresponding cancer type/ subtype, and a MAF file with nonsynonymous
#' mutations, and produces a scatter plot relating signif. count to sample size and mutation
#' frequency for a given gene.
#' @param list_of_master_dfs a list of master DFs, with the name being the cancer type or subtype
#' @param sample_size_dict a list relating each cancer type/ subtype to its sample size
#' @param patient_cancer_mapping a mapping of cancer types to the patient IDs of that cancer type
#' @param mutation_count_df a mutation count matrix from maftools
#' @param goi a driver gene-of-interest
#' @param qval_thres a q-value threshold for significance
relate_sample_size_to_sig_eQTL2 <- function(list_of_master_dfs, sample_size_dict, 
                                            patient_cancer_mapping, mutation_count_df,
                                            goi, qval_thres) {
  
  # Create a data frame that will link these factors
  input_df <- data.frame("cancer.type" = names(sample_size_dict), "sample.size" = as.numeric(sample_size_dict))
  
  mut_df_pats <- unlist(lapply(colnames(mutation_count_df)[2:(ncol(mutation_count_df)-1)], function(x)
    unlist(strsplit(x, "-", fixed = T))[1]))
  print(mut_df_pats)
  
  # Get the mutation frequency of the given GOI for each cancer type/subtype
  # Import get_mut_count_matrix from maf_helper_functions.R
  total.mut.count <- unlist(lapply(1:length(patient_cancer_mapping), function(i) {
    pats <- patient_cancer_mapping[[i]]
    ct <- names(patient_cancer_mapping)[i]
    mut_stat <- as.integer(mutation_count_df[mutation_count_df$Gene_Symbol == goi, (which(mut_df_pats %fin% pats)+1)])
    mut_count <- sum(unlist(lapply(mut_stat, function(m) ifelse(m > 0, 1, 0))))
    return(mut_count)
  }))
  names(total.mut.count) <- names(patient_cancer_mapping)
  total.mut.count <- as.data.frame(total.mut.count)
  total.mut.count$cancer.type <- names(patient_cancer_mapping)
  print(head(total.mut.count))
  
  # OPTION: use sample size to convert this to a fraction (% mutated)
  total.mut.count$frac.mutated <- unlist(lapply(1:nrow(total.mut.count), function(i) {
    type <- total.mut.count[i, 'cancer.type']
    print(type)
    samp.size <- as.numeric(input_df[input_df$cancer.type == type, 'sample.size'])
    print(samp.size)
    frac.mut <- total.mut.count[i, 'total.mut.count'] / samp.size
    return(frac.mut)
  }))
  
  # Merge this with our input DF by type
  input_df <- merge(input_df, total.mut.count, by = "cancer.type")
  
  # Get the number of significant eQTLs, for this goi if desired
  num.sig.hits <- unlist(lapply(1:length(list_of_master_dfs), function(i) {
    m <- list_of_master_dfs[[i]]
    m_sig <- m[(m$q.value < qval_thres) & (m$R_i.name == goi),]
    return(nrow(m_sig))
  }))
  names(num.sig.hits) <- names(list_of_master_dfs)
  num.sig.hits <- as.data.frame(num.sig.hits)
  num.sig.hits$cancer.type <- rownames(num.sig.hits)
  print(head(num.sig.hits))
  
  # Merge this with our input DF by type
  input_df <- merge(input_df, num.sig.hits, by = "cancer.type")
  print(head(input_df))
  
  # Now, create the plot with this input DF
  p <- ggplot(input_df, aes(x = frac.mutated, y = num.sig.hits, color = cancer.type)) + 
    geom_point(aes(size = sample.size)) + #scale_color_nejm() +
    scale_color_manual(values = c("#BC3C29FF", "#E18727FF", "#FFDC91FF", "khaki1", "seagreen2", "#20854EFF", "mediumaquamarine",
                                  "cyan2", "skyblue1", "#0072B5FF",  "darkblue", "#7876B1FF", "mediumorchid", "plum1",
                                  "lightpink", "#EE4C97FF", "palevioletred4", "beige", "wheat3", "saddlebrown", "gray", "black")) +
    theme_minimal() + xlab(paste("Frac. with Nonsyn. Mut in", goi)) + ylab(paste("Num. Significant Correlations (q <", paste0(qval_thres, ")"))) +
    labs(size = "Sample Size", color = "Cancer Type") + 
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
# For cancer types
# Use patient cancer mapping, see generate_protein_input_files.R
patient_cancer_mapping <- get_patient_cancer_mapping(clinical_df)

list_of_sample_sizes <- lapply(patient_cancer_mapping, function(x) length(x))
names(list_of_sample_sizes) <- unlist(lapply(names(patient_cancer_mapping), function(x) unlist(strsplit(x, "-", fixed = T))[2]))
list_of_sample_sizes <- list_of_sample_sizes[names(list_of_sample_sizes) %fin% names(list_of_master_dfs)]

# Call function
relate_sample_size_to_sig_eQTL2(list_of_master_dfs, list_of_sample_sizes, 
                                            patient_cancer_mapping, mutation_count_df,
                                            goi, qval_thres)

# PART D: OVERLAP IN PARTICULAR CORRELATIONS BETWEEN SUBGROUPS

#' Plot overlap in significant hits for a GOI between different subgroups using a Venn
#' diagram format
#' @param list_of_master_dfs a list of master DFs, with the name being the cancer type or subtype
#' @param goi a driver gene-of-interest's uniprot ID
#' @param qval_thres a q-value threshold for significance
#' @param top_n a threshold for the number of top hits to look at
plot_overlap_between_subgroups <- function(list_of_master_dfs, goi, qval_thres, top_n) {
  sig_hits <- lapply(1:length(list_of_master_dfs), function(i) {
    m <- list_of_master_dfs[[i]]
    #m[(m$R_i.name == goi) & (m$q.value < qval_thres), 'T_k.name'])
    sig_hits_tmp <- as.character(unlist(m[((grepl(goi, m$term)) & (m$q.value < qval_thres)), 'T_k.name']))
    print(paste(names(list_of_master_dfs)[i], length(sig_hits_tmp)))
    if(!is.na(top_n)) {
      if(length(sig_hits_tmp) > top_n) {sig_hits_tmp <- sig_hits_tmp[1:top_n]}
    }
    if(length(sig_hits_tmp) > 0) {return(sig_hits_tmp)}
    else {return(NA)}
  })
  names(sig_hits) <- names(list_of_master_dfs)
  sig_hits <- sig_hits[!is.na(sig_hits)]

  plt <- venn.diagram(sig_hits, category.names = names(sig_hits), 
                      filename = NULL, output = TRUE, lwd = 2, lty = 'blank', 
                      fill = c("#0072B5FF", "#BC3C29FF", "#20854EFF", "#E18727FF", "#7876B1FF"),#, "#FFDC91FF"), "#EE4C97FF"), 
                      cex = 2, fontface = "bold",
                      fontfamily = "sans", cat.cex = 2, cat.fontface = "bold", cat.fontfamily = "sans")
  #cat.default.pos = "outer", cat.pos = c(180, 180, 180, 180)) #, cat.fontfamily = "sans", rotation = 1)
  grid::grid.draw(plt)
}


#plot_overlap_between_subgroups(list_of_master_dfs, "TP53", 0.1, 50)
plot_overlap_between_subgroups(list_of_master_dfs, "P04637", 0.2, 50)

list_of_master_dfs <- list("BLCA" = blca_metabolic_res, "BRCA" = brca_metabolic_res, 
                           "CESC" = cesc_metabolic_res, "COAD" = coad_metabolic_res, 
                           "GBM" = gbm_metabolic_res, "HNSC" = hnsc_metabolic_res, 
                           "KICH" = kich_metabolic_res, "LGG" = lgg_metabolic_res, 
                           "LIHC" = lihc_metabolic_res, "LUAD" = luad_metabolic_res, 
                           "LUSC" = lusc_metabolic_res, "THCA" = thca_metabolic_res, 
                           "UCEC" = ucec_metabolic_res)


# Alternatively, display this using an UpSet plot
plot_overlap_between_subgroups_upset <- function(list_of_master_dfs, goi, qval_thres, 
                                                 top_n, target_set) {
  sig_hits <- lapply(1:length(list_of_master_dfs), function(i) {
    m <- list_of_master_dfs[[i]]
    #m[(m$R_i.name == goi) & (m$q.value < qval_thres), 'T_k.name'])
    sig_hits_tmp <- as.character(unlist(m[(m$R_i == goi) & (m$q.value < qval_thres), 
                             'T_k.name']))
    print(paste(names(list_of_master_dfs)[i], length(sig_hits_tmp)))
    if(!is.na(top_n)) {
      if(length(sig_hits_tmp) > top_n) {sig_hits_tmp <- sig_hits_tmp[1:top_n]}
    }
    return(sig_hits_tmp)
  })
  names(sig_hits) <- names(list_of_master_dfs)
  sig_hits <- sig_hits[lengths(sig_hits) != 0]
  print(sig_hits)
  print.table(intersect(sig_hits[["LIHC"]], intersect(sig_hits[["BRCA"]], intersect(sig_hits[["MESO"]],
                                         intersect(sig_hits[["LGG"]], sig_hits[["LUAD"]])))), 
              quote = F, row.names = F)
  
  # Optionally, use a Fisher's exact test to check whether or not the intersecting hits are
  # enriched in a particular target group
  #intersecting_genes <- lapply(1:length(sig_hits), function(i) {
    #if(!(i == length(sig_hits))) {
      #return(intersect(sig_hits[[i]], sig_hits[[i+1]]))
    #}
  #})
  #intersecting_genes <- unique(unlist(intersecting_genes))
  #all_targets <- unique(unlist(lapply(list_of_master_dfs, function(x) unlist(x$T_k.name))))
  #print(length(all_targets))
  #calculate_intersect_hg(intersecting_genes, unique(unlist(sig_hits)), target_set, all_targets)
  
  ylab <- paste("TP53 Hit Intersections, q <", qval_thres)
  
  sig_hits_top5 <- sig_hits[order(lengths(sig_hits), decreasing = T)][1:5]
  
  main_bar_col <- c("#0072B5FF")
  sets_bar_col <- c("#FFDC91FF")
  matrix_col <- c("#E18727FF")
  shade_col <- c("wheat4")
  if(!is.na(top_n)) {ylab <- paste0(ylab, paste(", top", paste(top_n, "hits per cancer type")))}
  plt <- upset(fromList(sig_hits_top5),
               order.by = 'freq', point.size = 3.5, line.size = 2, 
               mainbar.y.label = paste0("Intersection Size, q<", qval_thres), 
               sets.x.label = "Number of Hits", 
               text.scale = c(2, 2, 2, 2, 2, 2),
               main.bar.color = main_bar_col,
               sets.bar.color = sets_bar_col,
               matrix.color = matrix_col,
               shade.color = shade_col)
  print(plt)
}
plot_overlap_between_subgroups_upset(list_of_master_dfs, "P04637", 0.2, NA, tp53_curated_targets)

# Aside, for coloring hits:
tp53_intersect_hits <- "MDM2      SPATA18   SMS       G6PD      CDKN1A    EDA2R     CKS2      PGS1      PLK1      CEP95     HSD17B6   ERCC6L    RBM39     TMBIM4   
MAU2      TTLL3     BTG2      PNISR     POLA1     F5        FBXO38    SESN1     SUPT3H    CSTF2     ACOT9     ZNF404    PRDM4     SLC10A3  
TLK2      CGNL1     WHAMM     HHAT      RBM6      BHLHE40   NPTX2     CDKN3     RPL15     PGGT1B    STAT5B    TUBA4A    MISP      PAN3     
WDR6      TNFRSF10B BBS7      TRIM3     KANSL1    RHOB      PIP5K1B   USP36     CHD2"
tp53_intersect_hits <- unlist(strsplit(tp53_intersect_hits, "\n", fixed = T))
tp53_intersect_hits <- unlist(strsplit(tp53_intersect_hits, " ", fixed = T))
tp53_intersect_hits <- tp53_intersect_hits[tp53_intersect_hits != ""]
print(paste(tp53_intersect_hits[order(tp53_intersect_hits)], collapse = ", "), 
      quote = F, row.names = F)
for(h in tp53_intersect_hits) {
  print(h)
  e_meso <- top_drivers_0.05[["MESO"]][(top_drivers_0.05[["MESO"]]$R_i.name == "TP53") & (top_drivers_0.05[["MESO"]]$T_k.name == h), 'estimate']
  e_brca <- top_drivers_0.05[["BRCA"]][(top_drivers_0.05[["BRCA"]]$R_i.name == "TP53") & (top_drivers_0.05[["BRCA"]]$T_k.name == h), 'estimate']
  e_luad <- top_drivers_0.05[["LUAD"]][(top_drivers_0.05[["LUAD"]]$R_i.name == "TP53") & (top_drivers_0.05[["LUAD"]]$T_k.name == h), 'estimate']
  e_lihc <- top_drivers_0.05[["LIHC"]][(top_drivers_0.05[["LIHC"]]$R_i.name == "TP53") & (top_drivers_0.05[["LIHC"]]$T_k.name == h), 'estimate']
  e_lgg <- top_drivers_0.05[["LGG"]][(top_drivers_0.05[["LGG"]]$R_i.name == "TP53") & (top_drivers_0.05[["LGG"]]$T_k.name == h), 'estimate']
  if((e_meso > 0) & (e_brca > 0) & (e_luad > 0) & (e_lihc > 0) & (e_lgg > 0)) {print("all up")}
  else if((e_meso < 0) & (e_brca < 0) & (e_luad < 0) & (e_lihc < 0) & (e_lgg < 0)) {print("all down")}
  else {print("mixed")}
}

#' Stacked bar chart to show the number of significant genes in the same direction 
#' in at least two cancer types, in opposite directions in at least two cancer types,
#' and unique to only one cancer type
#' @param list_of_listed_master_dfs a nested list, outer list of driver genes, and
#' for each driver gene a list of master DFs by cancer type
#' @param qval_thres q-value threshold for significance
#' @param top_n an optional specification of the top n ranked targets, rather 
#' than using q-value
create_stacked_bar_directionality <- function(list_of_listed_master_dfs, 
                                              qval_thres, top_n) {
  drivers <- names(list_of_listed_master_dfs)
  
  dfs_list <- lapply(1:length(drivers), function(i) {
    d <- drivers[i]
    print(d)
    
    # For each cancer type, identify the significantly up- and down-regulated
    # target genes for the given driver
    list_of_master_dfs <- list_of_listed_master_dfs[[i]]
    sig_hits <- lapply(1:length(list_of_master_dfs), function(i) {
      m <- list_of_master_dfs[[i]]
      m_sub <- m[(m$R_i.name == d), ]
      if(!is.na(qval_thres)) {m_sub <- m_sub[m_sub$q.value < qval_thres, ]}
      if(!is.na(top_n)) {m_sub <- m_sub[1:min(nrow(m_sub), top_n), ]}
      
      sig_hits_up <- as.character(unlist(m_sub[m_sub$estimate > 0, 'T_k.name']))
      sig_hits_down <- as.character(unlist(m_sub[m_sub$estimate < 0, 'T_k.name']))
      
      sig_hits_list <- list("up.sig" = sig_hits_up, "down.sig" = sig_hits_down)
      
      return(sig_hits_list)
    })
    names(sig_hits) <- names(list_of_master_dfs)
    
    # Find hits that are in up.sig/ down.sig in at least 2 cases
    up.sig.full <- unlist(lapply(sig_hits, function(x) unlist(x$up.sig)))
    down.sig.full <- unlist(lapply(sig_hits, function(x) unlist(x$down.sig)))
    if((length(up.sig.full) == 0) & (length(down.sig.full) == 0)) {return(NA)}
    
    up.sig.full.dup <- unique(up.sig.full[duplicated(up.sig.full)])
    down.sig.full.dup <- unique(down.sig.full[duplicated(down.sig.full)])
    
    # Find hits that are up.sig/ down.sig in only 1 case
    up.sig.full.unique <- unique(up.sig.full[!duplicated(up.sig.full)])
    down.sig.full.unique <- unique(down.sig.full[!duplicated(down.sig.full)])
    
    # Get percentages
    total_num_hits <- length(c(up.sig.full.dup, down.sig.full.dup, 
                               up.sig.full.unique, down.sig.full.unique))

    # Add to data frame and return
    df <- data.frame('Category' = c("Up.Multi-Cancer", "Down.Multi-Cancer", 
                                    "Up.Single-Cancer", "Down.Single-Cancer"),
                     'Num.Sig.Hits' = c(length(up.sig.full.dup), 
                                        length(down.sig.full.dup),
                                        length(up.sig.full.unique), 
                                        length(down.sig.full.unique)),
                     'Perc.of.Hits' = c((length(up.sig.full.dup)/
                                           total_num_hits * 100),
                                        (length(down.sig.full.dup)/
                                           total_num_hits * 100),
                                        (length(up.sig.full.unique)/
                                           total_num_hits * 100),
                                        (length(down.sig.full.unique)/
                                           total_num_hits * 100)),
                     'Driver' = rep(d, times = 4))
    return(df)
  })
  dfs_list <- dfs_list[!is.na(dfs_list)]
  full_df <- do.call(rbind, dfs_list)
  full_df$Category <- factor(full_df$Category)
  
  full_df$Num.Sig.Hits.Label <- unlist(lapply(1:nrow(full_df), function(i) {
    perc <- full_df$Perc.of.Hits[i]
    num <- full_df$Num.Sig.Hits[i]
    return(ifelse(perc > 5, as.character(num), ""))
  }))
  
  cols <- c("Up.Multi-Cancer" = "#E18727FF", "Up.Single-Cancer" = "#E1872799", 
            "Down.Multi-Cancer" = "#0072B5FF", "Down.Single-Cancer" = "#0072B599")
  ylab <- "Percentage of Significant Hits"
  if(!is.na(qval_thres)) {ylab <- paste0(ylab, paste0(", q < ", qval_thres))}
  if(!is.na(top_n)) {ylab <- paste0(ylab, paste0(", Top", top_n))}
  
  # Create a stacked bar plot from this
  g <- ggplot(full_df, aes(fill=Category, y=Perc.of.Hits, 
                           x=reorder(Driver, -Num.Sig.Hits), 
                           label = Num.Sig.Hits.Label)) + # y=Num.Sig.Hits
    geom_bar(position="stack", stat="identity", color = "black") + 
    xlab("Driver Gene") + ylab(ylab) +
    theme_minimal() + #scale_fill_nejm() +
    scale_fill_manual(values = cols, 
                      name = "Direction of Regulation,\nMulti- or Single-Cancer") +
    geom_text(size = 4, position = position_stack(vjust = 0.5)) +
    theme(axis.title = element_text(face = "bold", size = 14), 
          axis.text.x = element_text(face = "bold", size = 12, angle = 45, 
                                     vjust = 1, hjust=1),
          axis.text.y = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12), legend.position = "bottom") #, legend.position = c(0.75,0.75))
  print(g)
  
  return(full_df)
}

top_drivers_0.05_tp53 <- lapply(top_drivers_0.05, function(x) x[x$R_i.name == "TP53",])

list_of_listed_master_dfs <- list("TP53" = top_drivers_0.05_tp53, "ARID1A" = top_drivers_0.05_arid1a,
                                  "CDKN2A" = top_drivers_0.05_cdkn2a, "CTNNB1" = top_drivers_0.05_ctnnb1,
                                  "EP300" = top_drivers_0.05_ep300, "FBXW7" = top_drivers_0.05_fbxw7,
                                  "KRAS" = top_drivers_0.05_kras, "NFE2L2" = top_drivers_0.05_nfe2l2,
                                  "PIK3CA" = top_drivers_0.05_pik3ca, "PTEN" = top_drivers_0.05_pten,
                                  "SMAD4" = top_drivers_0.05_smad4)

intersection_df <- create_stacked_bar_directionality(list_of_listed_master_dfs, 
                                                     0.2, NA)


#' Helper function to calculate Fisher's (HG) enrichment of intersecting genes in a 
#' target set
#' @param intersecting_genes the genes that intersect between any two groups
#' @param all_hits the set of total hits across all the groups
#' @param target_set a vector of known targets of the given GOI
#' @param all_targets the total set of gene targets tested
calculate_intersect_hg <- function(intersecting_genes, all_hits, target_set, all_targets) {

  num_intersect <- length(intersecting_genes)
  print(num_intersect)
  num_known_targets <- length(target_set)

  num_known_intersect_targets <- length(intersect(intersecting_genes, target_set))
  num_known_nonintersect_targets <- length(intersect(all_hits, target_set))
  num_notknown_intersect_targets <- length(intersect(intersecting_genes, all_targets))
  num_notknown_nonintersect_targets <- length(intersect(all_hits, all_targets))
  
  contigency_table <- matrix(c(num_known_intersect_targets, num_known_nonintersect_targets,
                               num_notknown_intersect_targets, num_notknown_nonintersect_targets),
                             nrow = 2, ncol = 2)
  print(contigency_table)
  
  fisher_res <- fisher.test(contigency_table, alternative = "greater")
  print(fisher_res$p.value)
}

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


# FOR PER-CANCER RESULTS: CREATE A BAR CHART SHOWING THE NUMBER OF SIGNIFICANT HITS 
# FOR EACH CANCER TYPE AT A GIVEN Q-VALUE THRESHOLD
#' @param list_of_master_dfs a list of master DFs, with the name being the cancer type or subtype
#' @param goi optional: a driver gene-of-interest's uniprot ID (count only hits for this GOI)
#' @param qval_thres a q-value threshold for significance
create_per_ct_sig_hits_barplot <- function(list_of_master_dfs, goi, qval_thres) {
  sig_hits_dts <- lapply(1:length(list_of_master_dfs), function(i) {
    m <- list_of_master_dfs[[i]]
    #m[(m$R_i.name == goi) & (m$q.value < qval_thres), 'T_k.name'])
    sig_hits <- length(unlist(m[(m$q.value < qval_thres), 'T_k.name']))
    if(!is.na(goi)) {
      sig_hits <- length(unlist(m[((m$R_i.name == goi) & (m$q.value < qval_thres)), 'T_k.name']))
    }
    return(data.table("Cancer.Type" = names(list_of_master_dfs)[i], "Num.Sig.Hits" = sig_hits))
  })
  sig_hits_dt_full <- do.call(rbind, sig_hits_dts)
  
  print(sig_hits_dt_full)
  
  p <- ggplot(sig_hits_dt_full, aes(x = reorder(Cancer.Type, Num.Sig.Hits), 
                            y = Num.Sig.Hits)) +
    geom_col(width = 0.7, fill = "#0072B5FF") + coord_flip() + theme_minimal() + 
    theme(legend.position = "none", axis.title = element_text(face = "bold", size = 14), 
          axis.text = element_text(face = "bold", size = 12)) + 
    ylab(paste("Number of hits at q <", qval_thres)) + xlab("Cancer Type") 
  #+ scale_fill_manual(values = color_palette)
  #+ scale_fill_nejm() 
  
  p
}
create_per_ct_sig_hits_barplot(top_drivers_0.05, NA, 0.1)

# ANOTHER OPTION: NETWORK REPRESENTATION
# Useful tutorial: https://mr.schochastics.net/material/netVizR/

#' Generate results graph display for a driver of interest 
#' This graphical display takes confidence interaction relationships from 
#' a network like STRING or HumanBase, and displays these relationships
#' as a graph. The graph will only include hits with a q-value below the
#' given threshold (and any intermediate genes?). The edges and genes
#' are colored by whether activation or repression is predicted when the 
#' driver is mutated. 
#' @param master_df an output DF from linear_model.R
#' @param goi the name(s) of the driver gene(s) of interest
#' @param qval_thres a q-value threshold, below which a target-driver
#' pairing is considered significant
#' @param network a table with the network information from either 
#' HumanBase or STRING
#' @param network_label either "STRING" or "HumanBase" to denote 
#' what type of network this is
#' @param conf_thres a confidence threshold for a connection between two 
#' genes in the given network
#' @param top_n if not NA, gives a value of the top N most significant targets to use
create_graphical_representation <- function(master_df, goi, qval_thres, 
                                            network, network_label, conf_thres, top_n) {
  
  # Subset the master DF to the given goi and q-value threshold
  master_df_sub <- master_df[(master_df$R_i.name == goi) & (master_df$q.value < qval_thres),]
  if(!is.na(top_n)) {master_df_sub <- master_df_sub[1:top_n, ]}
  nodes <-  c(goi, unique(master_df_sub$T_k.name))
  
  network_model_targs <- NA

  # Subset the network to the significant targets, plus the GOI
  if(network_label == "STRING") {
    network_model_targs <- network[(network$node1_name %fin% nodes) & 
                                     (network$node2_name %fin% nodes),]
  } else if (network_label == "HumanBase") {
    #network_model_targs <- network[(network$GENE1 %fin% nodes) & 
    #                                 (network$GENE2 %fin% nodes),]
    #colnames(network_model_targs) <- c("node1_name", "node2_name", "combined_score")
    network_model_targs <- network[(network$node1_name %fin% nodes) & 
                                     (network$node2_name %fin% nodes),]
    
  } else {
    print("Error in given network label. Only implemented for STRING and HumanBase networks.")
  }
  print(head(network_model_targs))

  # Remove cases where nodes 1 and 2 are the same, or rows where the values in the row are
  # duplicated in a different order. Divide all confidence scores by 1000 for STRING
  network_model_targs <- subset(network_model_targs, node1_name != node2_name & 
                                  !duplicated(cbind(pmin(node1_name, node2_name), 
                                                    pmax(node1_name, node2_name))))
  if(network_label == "STRING") {
    network_model_targs$combined_score <- network_model_targs$combined_score / 1000
  }

  #nodes_sub <- unique(c(network_model_targs$node1_name, network_model_targs$node2_name))

  # For each of the targets, get its a) directionality (up- or down-regulated), 
  # b) confidence score in the given network, in relation to the driver (and/or to the other targets?)
  edge_rows <- lapply(1:nrow(network_model_targs), function(i) {
    conf <- network_model_targs[i, 'combined_score']   
    if(conf > conf_thres) {
      g1 <- network_model_targs[i, 'node1_name']
      g2 <- network_model_targs[i, 'node2_name']
      return(data.frame("src" = g1, "target" = g2, 
                        "conf" = as.numeric(conf)))
    } else {return (NA)}
  })
  edge_rows <- edge_rows[!is.na(edge_rows)]
  edge_table <- do.call(rbind, edge_rows)
  colnames(edge_table) <- c("src", "target", "conf")
  print(edge_table)
  
  unique_nodes <- unique(c(edge_table$src, edge_table$target))
  node_table <- as.data.frame(unique_nodes)
  colnames(node_table) <- "nodes"
  node_table$dir <- as.factor(unlist(lapply(unique_nodes, function(n) {
    dir <- NA
    if(n == goi) {dir <- 0}
    else {
      dir <- ifelse(as.numeric(master_df_sub[master_df_sub$T_k.name == n, 'estimate'])[1] > 0, 1, -1)
      #dir <- ifelse(master_df_sub[(master_df_sub$R_i.name %in% c(g1, g2)) & 
       #                             (master_df_sub$T_k.name %in% c(g1, g2)), 'estimate'] > 0, 1, 2)
    }
    return(dir)
  })))
  node_table$node_size <- unlist(lapply(node_table$dir, function(d) {
    if(d == 0) {return(2)}
    else {return(1)}
  }))
  print(node_table)

  #links <- unique(edge_table[c("src", "target")])
  #links <- mutate(links, dir = edge_table$dir, weight = edge_table$conf)
  
  # Convert to ggraph format
  network_table <- tbl_graph(nodes = node_table, 
                             edges = edge_table, 
                             directed = FALSE)
  
  #print(network_table)
  
  # Create ggraph
  driver_targ_network <- ggraph(network_table, layout = "fr") +   # other options: stress
    geom_edge_link0(edge_colour = "darkgray", edge_width = 1) +  # aes(edge_width = conf), 
    geom_node_point(aes(fill = dir, size = node_size), shape = 21) +    #size = 5    
    scale_size(range = c(5,10), guide = 'none') + 
    #scale_edge_width(range = c(0.7, 1.25)) +
    geom_node_label(aes(label = nodes), repel = TRUE, show.legend = F, label.size = NA, label.padding = 0.05) +
    #geom_node_text(aes(label = unique_nodes), nudge_y = 0.2, nudge_x = 0.05, repel = T) + 
    theme_graph() +
    #scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#20854EFF"), 
                       #labels = c("Pred.Downreg", "Driver", "Pred.Upreg")) +
    scale_fill_manual(values = c("cornflowerblue", "lightgreen", "#E18727FF"), 
                        labels = c("Pred. Downregulation", "Driver", "Pred. Upregulation"),
                      name = "Direction of Regulation") +
    theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12), legend.position="bottom") 
  print(driver_targ_network)
  # Render the network
  show(driver_targ_network)
  
  return(node_table)
}


adjust_string_nw_files <- function(all_genes_id_conv) {
  # Import the String network file and adjust it
  string_nw <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/STRING/string.9606.protein.links.v11.5.txt",
                     header = TRUE)
  colnames(string_nw) <- c("node1", "node2", "combined_score")
  string_nw$index <- 1:nrow(string_nw)
  
  string_nw_info <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/STRING/string.9606.protein.info.v11.5.txt", 
                               header = T, sep = "\t")
  # Add info for any missing IDs
  missing_ids <- setdiff(string_nw_info$string_protein_id, unique(c(string_nw$node1, string_nw$node2)))
  missing_info <- lapply(missing_ids, function(id) {
    new_id <- unlist(strsplit(id, ".", fixed = T))[2]
    print(new_id)
    if(new_id %fin% all_genes_id_conv$ensembl_peptide_id) {
      id_row <- all_genes_id_conv[all_genes_id_conv$ensembl_peptide_id == new_id,]
      print(id_row)
      name <- unique(unlist(id_row[,'external_gene_name']))
      #length <- unlist(id_row[, 'end_position']) - unlist(id_row[, 'start_position'])
      return(data.frame("string_protein_id" = id, "preferred_name" = name, 
                        "protein_size" = NA, "annotation" = NA))
    } else {return(NA)}
  })
  missing_info <- missing_info[!is.na(missing_info)]
  missing_info_df <- do.call(rbind, missing_info)
  string_nw_info <- rbind(string_nw_info, missing_info_df)
  
  # Convert any ENSG IDs to external gene names
  string_nw_info$preferred_name <- unlist(lapply(string_nw_info$preferred_name, function(n) {
    new_n <- n
    if(grepl("ENSG", n, fixed = T)) {
      new_n <- paste(unique(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == n, 
                                              'external_gene_name']), collapse = ";")
    }
    return(new_n)
  }))
  
  missing_ids_from_run <- setdiff(pc_allGenes$T_k.name, string_nw_info$preferred_name)
  missing_info_from_run <- lapply(missing_ids_from_run, function(name) {
    ensp <- unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == name, 'ensembl_peptide_id'])
    ensp <- paste0("9606.", ensp)
    return(data.frame("string_protein_id" = ensp, "preferred_name" = name, 
                      "protein_size" = NA, "annotation" = NA))
  })
  missing_info_from_run_df <- do.call(rbind, missing_info_from_run)
  missing_info_from_run_df <- missing_info_from_run_df[!(missing_info_from_run_df$string_protein_id == "9606."),]
  string_nw_info <- rbind(string_nw_info, missing_info_from_run_df)
  
  string_nw_pt1 <- merge(string_nw[, c('node1', 'combined_score', 'index')], string_nw_info[, c('string_protein_id', 'preferred_name')], 
                         by.x = 'node1', by.y = 'string_protein_id')
  colnames(string_nw_pt1) <- c("node1", "combined_score", "index", "node1_name")
  string_nw_pt1 <- string_nw_pt1[order(string_nw_pt1$index),]
  
  string_nw_pt2 <- merge(string_nw[, c('node2', 'combined_score', 'index')], string_nw_info[, c('string_protein_id', 'preferred_name')], 
                         by.x = 'node2', by.y = 'string_protein_id')
  colnames(string_nw_pt2) <- c("node2", "combined_score", "index", "node2_name")
  string_nw_pt2 <- string_nw_pt2[order(string_nw_pt2$index),]
  
  string_nw_full <- merge(string_nw_pt1, string_nw_pt2, by = c('combined_score', 'index'))
  
  string_nw_full$node1 <- unlist(lapply(string_nw_full$node1, function(x) unlist(strsplit(x, ".", fixed = T))[2]))
  string_nw_full$node2 <- unlist(lapply(string_nw_full$node2, function(x) unlist(strsplit(x, ".", fixed = T))[2]))
  
  return(string_nw_full)
}

string_nw_full <- adjust_string_nw_files(all_genes_id_conv)
fwrite(string_nw_full, "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/STRING/string.9606.protein.links.v11.5.namesAdded.txt")

# Read back
string_nw_full <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/STRING/string.9606.protein.links.v11.5.namesAdded.txt")

# Import the HumanBase tissue-specific networks and adjust

# Import the network interaction files for the genes of interest
nw_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/Network_Data/"
nw_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/"


# BREAST
goi <- "TP53"
goi <- "PIK3CA"

tp53_mamm_epi_nw <- read.table(paste0(nw_path, "hb_mammary_epithelium_tp53_network.txt"), 
                               header = TRUE, sep = ",")
tp53_mamm_gl_nw <- read.table(paste0(nw_path, "hb_mammary_gland_tp53_network.txt"), 
                              header = TRUE, sep = ",")
pik3ca_mamm_epi_nw <- read.table(paste0(nw_path, "hb_mammary_epithelium_pik3ca_network.txt"), 
                                 header = TRUE, sep = ",")
pik3ca_mamm_gl_nw <- read.table(paste0(nw_path, "hb_mammary_gland_pik3ca_network.txt"), 
                                header = TRUE, sep = ",")

tp53_mamm_epi_nw <- as.data.frame(lapply(tp53_mamm_epi_nw, trimws))
tp53_mamm_gl_nw <- as.data.frame(lapply(tp53_mamm_gl_nw, trimws))
pik3ca_mamm_epi_nw <- as.data.frame(lapply(pik3ca_mamm_epi_nw, trimws))
pik3ca_mamm_gl_nw <- as.data.frame(lapply(pik3ca_mamm_gl_nw, trimws))

conf_thres <- 0.4
tp53_mamm_epi_nw_sig <- tp53_mamm_epi_nw[tp53_mamm_epi_nw$WEIGHT > conf_thres,]
tp53_mamm_gl_nw_sig <- tp53_mamm_gl_nw[tp53_mamm_gl_nw$WEIGHT > conf_thres,]
pik3ca_mamm_epi_nw_sig <- pik3ca_mamm_epi_nw[pik3ca_mamm_epi_nw$WEIGHT > conf_thres,]
pik3ca_mamm_gl_nw_sig <- pik3ca_mamm_gl_nw[pik3ca_mamm_gl_nw$WEIGHT > conf_thres,]

tp53_nw_sig <- rbind(tp53_mamm_epi_nw_sig, tp53_mamm_gl_nw_sig)
pik3ca_nw_sig <- rbind(pik3ca_mamm_epi_nw_sig, pik3ca_mamm_gl_nw_sig)

# Full network
hb_mammary_epithelium <- read.table(paste0(nw_path, "HumanBase/mammary_epithelium.txt"))
colnames(hb_mammary_epithelium) <- c("node1_entrez", "node2_entrez", "combined_score")
hb_mammary_epithelium_conf <- hb_mammary_epithelium[hb_mammary_epithelium$combined_score >= 0.4, ]

print(length(unique(c(hb_mammary_epithelium_conf$node1_entrez, hb_mammary_epithelium_conf$node2_entrez))))
unique_ids <- unique(c(hb_mammary_epithelium_conf$node1_entrez, hb_mammary_epithelium_conf$node2_entrez))
mapping <- data.frame("entrez" = unique_ids)
mapping$gene_name <- unlist(lapply(mapping$entrez, function(e)
  tryCatch({
    return(bitr(e, fromType = "ENTREZID", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)[['SYMBOL']])
  }, error = function(cond) {
    print(cond)
    return(NA)
  })))

hb_mammary_epithelium_conf$node1_name <- unlist(lapply(hb_mammary_epithelium_conf$node1_entrez, function(e)
  return(mapping[mapping$entrez == e, 'gene_name'])))
hb_mammary_epithelium_conf$node2_name <- unlist(lapply(hb_mammary_epithelium_conf$node2_entrez, function(e)
  return(mapping[mapping$entrez == e, 'gene_name'])))
write.csv(hb_mammary_epithelium_conf, paste0(nw_path, "HumanBase/mammary_epithelium_0.4conf.csv"))

hb_mammary_epithelium_conf <- read.csv(paste0(nw_path, "HumanBase/mammary_epithelium_0.4conf.csv"))


# Call function
create_graphical_representation(top_drivers_0.05_brca, "TP53", 0.2, string_nw_full, "STRING", 0.4, 100)
create_graphical_representation(top_drivers_0.05_brca, "TP53", 0.2, hb_mammary_epithelium_conf, "HumanBase", 0.4, 100)

# TESTER VERSION
nodes <- c("ATV", "MTV", "ETX", "MMH", "ELOV", "READ")
network_model_targs <- data.frame("node1" = sample(nodes, size = 24, replace = TRUE), 
                                  "node2" = sample(nodes, size = 24, replace = TRUE), 
                                  "combined_score" = runif(24), "estimate" = runif(24, min = -1, max = 1))


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
                  Rowv = TRUE, key = TRUE, key.title = NA, key.xlab = "T-statistic") 
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
  #matrix <- data.frame(matrix(ncol = length(unique(results_table$R_i.name)),
                              #nrow = length(unique(results_table$T_k.name))))
  matrix <- data.frame(matrix(ncol = length(unique(results_table$term)),
                              nrow = length(unique(results_table$T_k.name))))
  unique_ri_terms <- unique(results_table$term)
  r_i.names <- unlist(lapply(unique_ri_terms, function(x) {
    regprot <- unlist(strsplit(x, "_", fixed = TRUE))[1]
    name <- unique(unlist(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == regprot, 
                                            'external_gene_name']))
    return(name)
  }))
  mapping <- data.frame("swissprot" = unlist(lapply(unique_ri_terms, function(x)
    unlist(strsplit(x, "_", fixed = TRUE))[1])), "name" = r_i.names)
  colnames(matrix) <- r_i.names
  rownames(matrix) <- unique(results_table$T_k.name)
  
  # Fill in this table with t-statistics or -log10(qval)*dir
  for (i in 1:nrow(results_table)) {
    tstat <- results_table$statistic[i]
    #beta <- results_table$estimate[i]
    dir <- ifelse(results_table$estimate[i] > 0, 1, -1)
    neglog10qvaltimesdir <- (-log10(results_table$q.value[i]))*dir
    regprot <- results_table$R_i.name[i]
    #regprot_uniprot <- unlist(strsplit(results_table$term[i], "_", fixed = TRUE))[1]
    #regprot <- mapping[mapping$swissprot == regprot_uniprot, 'name']
    print(regprot)
    targ <- results_table$T_k.name[i]
    
    tryCatch({
      #matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- tstat
      #matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- beta
      matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- neglog10qvaltimesdir
    }, error=function(cond){
      print(cond)
      matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- 0
    })
  }
  # NOTE: leave all the unfilled pairings as NA
  print(matrix)
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
                  density.info = "none", key = TRUE, key.title = NA, key.xlab = "-log10(qval)*direction") #key.xlab = "T-statistic")
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
  cols_branches <- sample(col_vector, num_clust, replace = T)
  
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
  #col_pal = brewer.pal(n = 100, name = "YlGnBu")
  new_hm <- heatmap.2(matrix, scale = "none", col = colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100), 
                      trace = "none", density.info = "none", dendrogram = "row", Rowv = dendrogram, 
                      Colv = "none", key.title = NA,  key.xlab = "-log10(qval)*direction", #key.xlab = "T-statistic", 
                      RowSideColors = col_labels, colRow = col_labels) 
  #new_hm <- pheatmap(matrix, angle_col = 45)
   
  print(new_hm)
  
  # Return the labeled DF
  return(h_sort_df)
}


# Call this function
example_df <- read.csv("C:/Users/sarae/Documents/output_results_cancerRelated_mattMetabolicTargs_iprotein_tmmSDGr1_rawCNA_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_corrected_MUT.csv",
                       header = TRUE, check.names = FALSE)
create_heatmap_with_clustering(example_df, height = 2)


one_carbon_metabolism_genes <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/1c_metabolic_targets.csv")[, 'swissprot']
one_carbon_metabolism_genes <- unlist(lapply(one_carbon_metabolism_genes, function(g) unlist(strsplit(g, ";", fixed = T))))
glycolysis_genes <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/glycolysis_metabolic_targets.csv")[, 'swissprot']
tca_cycle_genes <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/tca_cycle_metabolic_targets.csv")[, 'swissprot']
aa_metabolism_genes <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/aa_metabolism_targets.csv")[, 'Identifier']
# oxidative phosphorylation
kegg_oxidative_phosphorylation <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/KEGG/KEGG_OXIDATIVE_PHOSPHORYLATION.v2023.1.Hs.csv")
kegg_oxidative_phosphorylation_genes <- unlist(strsplit(kegg_oxidative_phosphorylation[17,2], ",", fixed = T))
kegg_oxidative_phosphorylation_genes <- kegg_oxidative_phosphorylation_genes[kegg_oxidative_phosphorylation_genes != ""]  #132 genes

# valine, leucine, and isoleucine degradation
kegg_vli_deg <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/KEGG/KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION.v2023.1.Hs.tsv.csv")
kegg_vli_deg_genes <- unlist(strsplit(kegg_vli_deg[17,2], ",", fixed = T))
kegg_vli_deg_genes <- kegg_vli_deg_genes[kegg_vli_deg_genes != ""]  #44 genes

# Pyrimidine metabolism
kegg_pyrimidine_metabolism <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/KEGG/KEGG_PYRIMIDINE_METABOLISM.v2023.1.Hs.csv")
kegg_pyrimidine_metabolism_genes <- unlist(strsplit(kegg_pyrimidine_metabolism[17,2], ",", fixed = T))
kegg_pyrimidine_metabolism_genes <- kegg_pyrimidine_metabolism_genes[kegg_pyrimidine_metabolism_genes != ""]  #98 genes


create_heatmap_with_clustering(example_df[example_df$T_k %fin% one_carbon_metabolism_genes,], height = 2)


# PART B: UP- AND DOWN-REGULATED METABOLIC HITS
#' Create a volcano plot for a given GOI in a particular cancer type or subtype,
#' with significant hits highlighted and labeled
#' @param master_df a data frame produced by the model
#' @param goi a driver gene-of-interest uniprot ID
#' @param qval_thres a q-value threshold for significance
create_volcano_plot <- function(master_df, goi, qval_thres) {
  
  # Subset to the GOI 
  master_df <- master_df[master_df$R_i == goi,]
  #master_df <- master_df[grepl(goi, master_df$term),]
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
    ifelse(-log10(master_df$q.value[i]) > 5, master_df$T_k.name[i], NA)))
  print(head(gene_names))
  
  # Put this into a DF for plotting
  plot_df <- data.frame("log2_beta" = log2_beta, "neg_log10_qval" = neg_log10_qval,
                        "up_or_down" = up_or_down, "gene" = gene_names)
  print(max(plot_df$neg_log10_qval[plot_df$neg_log10_qval != Inf]) + 5)
  plot_df[which(plot_df$neg_log10_qval == Inf), 'neg_log10_qval'] <- max(
    plot_df$neg_log10_qval[plot_df$neg_log10_qval != Inf]) + 5
  
  print(head(plot_df))
  
  # Optionally squish each point outside the given x-axis limits
  xlim_min <- -1.5
  xlim_max <- 1.5
  plot_df[which(plot_df$log2_beta > xlim_max), 'log2_beta'] <- xlim_max
  plot_df[which(plot_df$log2_beta < xlim_min), 'log2_beta'] <- xlim_min

  # Set color and alpha params
  cols <- c("up" = "#E18727FF", "down" = "#0072B5FF", "ns" = "grey") 
  alphas <- c("up" = 0.95, "down" = 0.95, "ns" = 0.5)

  p <- ggplot(plot_df, aes(x = log2_beta, y = neg_log10_qval, fill = up_or_down, 
                           alpha = up_or_down, label = gene)) + 
    geom_point(shape = 21, size = 2.5) + xlab("Mutation Coefficient") + 
    ylab("-log10(q-value)") + theme_minimal() + 
    ylim(0, max(plot_df$neg_log10_qval[plot_df$neg_log10_qval != Inf]) + 1) +
    geom_hline(yintercept = -log10(qval_thres), linetype = "dashed") +
    scale_fill_manual(values = cols) + scale_alpha_manual(values = alphas) + 
    theme(legend.position = "none", 
          axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 12)) + 
    geom_text_repel(min.segment.length = unit(0.1, "lines"), size = 4, 
                    force = 3.5, box.padding = 0.4) + xlim(xlim_min, xlim_max)
                    #segment.curvature = -0.1, segment.ncp = 3, segment.angle = 20) 
  
  p
}

create_volcano_plot(allgenes_p53, "P04637", 0.1)


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
  mutant_patients <- get_mutant_patients2(mutation_regprot_df, driver_name)
  print(head(mutant_patients))
  
  # Get expression of target gene within each mutational group
  expression_by_group <- get_expression_by_mut_group(expression_df, mutant_patients, target)
  expression_mutants <- expression_by_group$mutants
  expression_normal <- expression_by_group$normal

  # Pad them with NAs as needed, in order to put them in a data frame
  if(length(expression_mutants) > length(expression_normal)) {
    nas <- rep(NA, times = length(expression_mutants) - length(expression_normal))
    expression_normal <- c(expression_normal, nas)
  } else if (length(expression_normal) > length(expression_mutants)) {
    nas <- rep(NA, times = length(expression_normal) - length(expression_mutants))
    expression_mutants <- c(expression_mutants, nas)
  } else {print("Same length. No NAs added.")}
    
  dataframe_expr <- data.frame("No.Mutation" = expression_normal, "Mutation" = expression_mutants)
  print(mean(na.omit(dataframe_expr$No.Mutation)))
  print(mean(na.omit(dataframe_expr$Mutation)))
  dataframe_expr_m <- na.omit(melt(dataframe_expr))
  print(head(dataframe_expr_m))
  
  # Create a boxplot or violin plot
  if(box_or_violin == "box") {
    #boxplot(dataList_expr, ylab = "Expression", main = paste(target, paste("Expression By", paste(driver_ensg, "Mutation Status"))))
    
    p <- ggplot(dataframe_expr_m, aes(x = variable, y = value, fill = variable)) + 
      geom_boxplot(size = 1) + #geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.2) + 
      theme_minimal() + scale_fill_manual(values=c("#6F99ADFF","#E18727FF")) +
      xlab(paste(driver_name, "Mutation Status")) + ylab(paste(target_name, "Expression")) +
      theme(legend.position = "none", axis.title = element_text(face = "bold", size = 14), 
            axis.text = element_text(face = "bold", size = 12)) +
      stat_compare_means(method = "wilcox.test")
      #scale_x_discrete(labels=c(paste(driver_name, "Wild-Type"), paste(driver_name, "Mutant")))
  
  } else {
    p <- ggplot(dataframe_expr_m, aes(x = variable, y = value, fill = variable)) + 
      geom_violin(trim = FALSE) + #geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.2) + 
      theme_minimal() + scale_fill_manual(values=c("#6F99ADFF","#E18727FF")) +
      xlab(paste(driver_name, "Mutation Status")) + ylab(paste(target_name, "Expression")) +
      theme(legend.position = "none", axis.title = element_text(face = "bold", size = 14), 
            axis.text = element_text(face = "bold", size = 12)) + 
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

get_mutant_patients2 <- function(mutation_count_matrix, regprot_name) {
  regprot_row <- mutation_count_matrix[mutation_count_matrix$Gene_Symbol == regprot_name, ]
  ind <- which(as.integer(regprot_row) > 0)
  return(colnames(regprot_row)[ind])
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


create_driver_target_boxplot(expression_df_qn, mutation_count_matrix, "ENSG00000141510", "P04637", 
                             "ENSG00000131747", "box", "TP53", "TOP2A")

clinical_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/clinical_data_subset.csv",
                        header = T, check.names = F)

# Generate a sample/patient to cancer type mapping
get_patient_cancer_mapping <- function(clinical_df) {
  specific_types <- unique(clinical_df$project_id)
  patient_cancer_mapping <- lapply(specific_types, function(ct) {
    pats <- clinical_df[grepl(ct, clinical_df$project_id),'case_submitter_id']
    pats_ids <- unlist(lapply(pats, function(pat) 
      unlist(strsplit(pat, "-", fixed = TRUE))[3]))
    return(unique(pats_ids))
  })
  names(patient_cancer_mapping) <- specific_types
  return(patient_cancer_mapping)
}

# Create mapping, if applicable
patient_cancer_mapping <- get_patient_cancer_mapping(clinical_df)

mutation_count_matrix_pats <- unlist(lapply(colnames(mutation_count_matrix)[2:ncol(mutation_count_matrix)], function(x)
  unlist(strsplit(x, "-", fixed = T))[1]))
expression_df_quantile_norm_pats <- unlist(lapply(colnames(expression_df_quantile_norm), function(x)
  unlist(strsplit(x, "-", fixed = T))[1]))
create_driver_target_boxplot(expression_df_quantile_norm[, which(expression_df_quantile_norm_pats %fin% patient_cancer_mapping$`TCGA-COAD`)], 
                             mutation_count_matrix[, c(1, (which(mutation_count_matrix_pats %fin% patient_cancer_mapping$`TCGA-COAD`)+1))], 
                             "ENSG00000133703", "P01116", "ENSG00000135318", "box", "KRAS", "NT5E")

#' Create a cross-cancer type boxplot to highlight driver-target relationships that are unique to
#' a particular cancer type
#' @param expression_df a normalized and filtered expression DF
#' @param patient_cancer_mapping a mapping of patients to their cancer type
#' @param mutation_count_df a mutation count matrix
#' @param driver_ensg the driver gene ENSG ID
#' @param target the target gene ENSG ID
#' @param box_or_violin either "box" or "violin" to denote whether we want a 
#' boxplot or a violin version of the plot
#' @param driver_name the gene name of the driver gene, for labeling
#' @param target_name the gene name of the target gene, for labeling
#' @param ct_of_interest optionally, a cancer type of interest, so that we can color it differently
create_cross_cancer_expr_boxplot <- function(expression_df, patient_cancer_mapping, mutation_count_df,
                                             driver_ensg, target, box_or_violin, driver_name, 
                                             target_name, ct_of_interest) {
  
  expr_df_pats <- unlist(lapply(colnames(expression_df), function(x) 
    unlist(strsplit(x, "-", fixed = T))[1]))
  mut_count_df_pats <- unlist(lapply(colnames(mutation_count_df)[2:ncol(mutation_count_df)], function(x) 
    unlist(strsplit(x, "-", fixed = T))[1]))
  
  per_cancer_logfc <- lapply(1:length(patient_cancer_mapping), function(i) {
    ct <- names(patient_cancer_mapping)[i]
    pats <- patient_cancer_mapping[[i]]
    
    expr_df_ct <- expression_df[, which(expr_df_pats %fin% pats)]
    mut_df_ct <- mutation_count_df[, c(1, (which(mut_count_df_pats %fin% pats)+1), ncol(mutation_count_df))]
    
    # Get mutant/ non-mutant patients
    mutant_patients <- get_mutant_patients2(mut_df_ct, driver_name)
    #print(head(mutant_patients))
    
    # Get logfc of expression of target gene between mutational groups
    expression_by_group <- get_expression_by_mut_group(expr_df_ct, mutant_patients, target)
    expression_mutants <- expression_by_group[[1]]
    expression_normal <- expression_by_group[[2]]
    #expression_logfc <- log2(as.numeric(expression_mutants) / as.numeric(expression_normal))
    expression_logfc <- log2(as.numeric(median(expression_mutants)) / as.numeric(median(expression_normal)))
    
    print(paste(ct, expression_logfc))
    return(expression_logfc)
    
    #return(list("mutant" = expression_mutants, "normal" = expression_normal))
  })
  names(per_cancer_logfc) <- names(patient_cancer_mapping)
  
  # Pad them with NAs as needed, in order to put them in a data frame
  #max_length <- max(lengths(per_cancer_logfc))
  #print(max_length)
  #per_cancer_logfc <- lapply(per_cancer_logfc, function(x) if(length(x) < max_length) {
  #  nas <- rep(NA, times = max_length - length(x))
  #  return(c(x, nas))
  #})
  #print(per_cancer_logfc)
  
  #per_cancer_logfc_df <- as.data.frame(do.call(cbind, per_cancer_logfc))
  per_cancer_logfc_df <- as.data.frame(do.call(rbind, per_cancer_logfc))
  colnames(per_cancer_logfc_df) <- 'logfc.expression'
  per_cancer_logfc_df$Cancer.Type <- names(per_cancer_logfc)
  per_cancer_logfc_df$Cancer.Type.of.Interest <- unlist(lapply(per_cancer_logfc_df$Cancer.Type, function(ct)
    if(ct == ct_of_interest) {return(TRUE)} else{return(FALSE)}))
  per_cancer_logfc_df <- per_cancer_logfc_df[order(per_cancer_logfc_df$logfc.expression, decreasing = T),]
  print(per_cancer_logfc_df)
  #per_cancer_logfc_df_m <- melt(per_cancer_logfc_df)
  #print(head(per_cancer_logfc_df_m))
  #per_cancer_logfc_df_m$Cancer.Type.of.Interest <- unlist(lapply(per_cancer_logfc_df_m$variable, function(v) 
    #if(v == ct_of_interest) {return(T)} else{return(F)}))
  #per_cancer_logfc_df_m <- na.omit(per_cancer_logfc_df_m)
  #print(per_cancer_logfc_df_m)
  
  # Create a boxplot or violin plot
  if(box_or_violin == "box") {
    #boxplot(dataList_expr, ylab = "Expression", main = paste(target, paste("Expression By", paste(driver_ensg, "Mutation Status"))))
    
    p <- ggplot(per_cancer_logfc_df, aes(x = reorder(Cancer.Type, -logfc.expression), y = logfc.expression, fill = Cancer.Type.of.Interest)) + 
      geom_bar(stat = "identity") + #geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.2) + 
      theme_minimal() + #scale_fill_manual(values=c("#6F99ADFF", "#E18727FF")) +
      xlab("Cancer Type") + ylab(paste(target_name, paste("LogFC Expression, by", paste(driver_name, "Mutation Status\n")))) +
      theme(legend.position = "none", axis.title = element_text(face = "bold")) #+ 
      #scale_x_discrete(labels=c(paste(driver_name, "Wild-Type"), paste(driver_name, "Mutant")))
    
  } else {
    p <- ggplot(per_cancer_logfc_df_m, aes(x = variable, y = value, fill = Cancer.Type.of.Interest)) + 
      geom_violin(trim = FALSE) + #geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.2) + 
      theme_minimal() + scale_fill_manual(values=c("#6F99ADFF", "#E18727FF")) +
      xlab("Cancer Type") + ylab(paste(target_name, paste("LogFC Expression, by", paste(driver_name, "Mutation Status\n")))) +
      theme(legend.position = "none", axis.title = element_text(face = "bold")) + 
      #scale_x_discrete(labels=c(paste(driver_name, "Wild-Type"), paste(driver_name, "Mutant"))) +
      stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black", size = 1)  # FOR POINT RANGE
    #stat_summary(fun.data=mean_sdl, geom="crossbar", width=0.08)  # FOR CROSSBAR
  }
  p
  
}


# Uniprot and ENSG IDs for driver gene of interest
tp53_uniprot <- "P04637" # TP53
tp53_ensg <- "ENSG00000141510"

pik3ca_uniprot <- "P42336"
pik3ca_ensg <- "ENSG00000121879"

idh1_uniprot <- "O75874"
idh1_ensg <- "ENSG00000138413"
  

# ENSG ID for target gene of interest
pfkfb3 <- "ENSG00000170525"

# Import expression DF
expression_df_tmm <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_tmmSDGr1filtByExpr_CancerOnly_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
expression_df_quantile_norm <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
rownames(expression_df_quantile_norm) <- expression_df_quantile_norm$ensg_id
expression_df_quantile_norm <- expression_df_quantile_norm[,3:ncol(expression_df_quantile_norm)]

# Import mutation DF
mutation_count_df <- read.csv(paste(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv", sep = ""), 
                                header = TRUE, check.names = FALSE)
# Call function
create_driver_target_boxplot(expression_df_quantile_norm, mutation_count_df, tp53_ensg, tp53_uniprot, pfkfb3, "box", "TP53", "PFKFB3")
create_driver_target_boxplot(expression_df_quantile_norm, mutation_count_df,tp53_ensg, tp53_uniprot, pfkfb3, "violin", "TP53", "PFKFB3")


# Find candidates and call function
other_ct_hits <- unique(unlist(lapply(top_drivers_0.05_metabol[!names(top_drivers_0.05_metabol) == "BRCA"], function(x) 
  if("PIK3CA" %in% x$R_i.name) {unlist(x[(x$R_i.name == "PIK3CA") & (x$q.value < 0.1), 'T_k.name'])} else{return(NA)})))
other_ct_hits <- other_ct_hits[!is.na(other_ct_hits)]
setdiff(unlist(top_drivers_0.05_metabol$BRCA[(top_drivers_0.05_metabol$BRCA$R_i.name == "PIK3CA") & 
                                               (top_drivers_0.05_metabol$BRCA$q.value < 0.1), 'T_k.name']), other_ct_hits)

patient_cancer_mapping_sub_p53mut <- patient_cancer_mapping_sub[names(patient_cancer_mapping_sub) %in% c("LIHC", "LUAD", "BLCA", "PRAD", "ACC", "LUSC", "KICH", 
                                                                                                         "SARC", "UCEC", "STAD", "ESCA", "LGG", "BRCA", "HNSC", 
                                                                                                         "PAAD", "COAD", "MESO")]
create_cross_cancer_expr_boxplot(expression_df_quantile_norm, patient_cancer_mapping, mutation_count_df, 
                                 tp53_ensg, pfkfb3, "box", "TP53", "PFKFB3", "COAD")

patient_cancer_mapping_sub_pik3camut <- patient_cancer_mapping_sub[names(patient_cancer_mapping_sub) %in% c("BRCA", "COAD", "BLCA", 
                                                                                                            "CESC", "UCEC")]
create_cross_cancer_expr_boxplot(expression_df_quantile_norm, patient_cancer_mapping_sub_pik3camut, mutation_count_df, 
                                 pik3ca_ensg, "ENSG00000100596", "box", "PIK3CA", "SPTLC2", "BRCA")

# PART D: META-ANALYSIS FOREST PLOT <SEE PERFORM_META_ANALYSIS.R> 



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
  
  master_df <- master_df[order(master_df$q.value),]
  if(length(silent_df) > 1) {
    silent_df <- silent_df[order(silent_df$q.value),]
  }
  
  # For each source, create a 0 and 1 vector for each  target hit to indicate if it 
  # is a a hit in that source
  target_set_vector_list <- lapply(target_sets_list, function(source) {
    source <- setdiff(source, setdiff(source, master_df$T_k.name))
    vect <- unlist(lapply(1:nrow(master_df), function(i) {
      x <- master_df[i, 'T_k.name']
      return(ifelse(x %fin% source, 1, 0))
    }))
    return(vect)
  })
  
  target_set_vector_list_silent <- NA
  if(length(silent_df) > 1) {
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
  if(length(silent_df) > 1) {
    target_set_fraction_df_silent <- get_frac_of_genes_at_rank(target_set_vector_list_silent,
                                                               silent_df)
    target_set_fraction_df <- merge(target_set_fraction_df, target_set_fraction_df_silent, 
                                    by = "Rank")
  }
  
  
  target_sets_fraction_df_m <- melt(target_set_fraction_df, "Rank")
  colnames(target_sets_fraction_df_m) <- c("Rank", "Source", "Frac")
  target_sets_fraction_df_m$Mutation.Type <- unlist(lapply(target_sets_fraction_df_m$Source, function(x) {
    val <- "Nonsynonymous"
    if(grepl("S.", x, fixed = TRUE)) {val <- "Synonymous"}
    return(val)
  }))
  #target_sets_fraction_df_m$Source.Broad <- unlist(lapply(as.character(target_sets_fraction_df_m$Source), function(x)
  #  unlist(strsplit(x, ".", fixed = T))[2]))
  
  print(head(target_sets_fraction_df_m))
  print(unique(target_sets_fraction_df_m$Source))
  #print(hist(target_sets_fraction_df_m[target_sets_fraction_df_m$Mutation.Type == "Synonymous", 'Frac']))
  print(ggplot(target_sets_fraction_df_m[(target_sets_fraction_df_m$Mutation.Type == "Synonymous") & 
                                           (target_sets_fraction_df_m$Source == "S.Curated_Fischer_2017"),],
               aes(x = Rank, y = Frac, group = 1)) + geom_line())
  print(ggplot(target_sets_fraction_df_m[(target_sets_fraction_df_m$Mutation.Type == "Synonymous") & 
                                           (target_sets_fraction_df_m$Source == "S.DoRothEA"),],
               aes(x = Rank, y = Frac, group = 1)) + geom_line())
  return(target_sets_fraction_df_m)
  
  #color_list <- c("#BC3C29FF", "#0072B5FF")
  #unique_sources <- unique(unlist(lapply(as.character(target_sets_fraction_df_m$Source), function(x)
  #  unlist(strsplit(x, ".", fixed = T))[2])))
  #names(color_list) <- unique_sources
  
  # Create inset plot
  p2 <- ggplot(target_sets_fraction_df_m[target_sets_fraction_df_m$Rank %in% 1:cutout_thres,], 
               mapping = aes(x = Rank, y = Frac, color = Source)) + 
    geom_line(aes(linetype=Mutation.Type, alpha = Mutation.Type), size = 1.25) +
    #geom_point(size = 1) + 
    scale_x_continuous(limits = c(1,cutout_thres)) + 
    scale_color_nejm() + 
    scale_alpha_manual(values=c(0.9,0.9,0.4,0.4)) +
    scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#BC3C29FF", "#0072B5FF")) + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.minor = element_blank(), axis.text.x = element_text(size = 14, face = "bold"), 
          axis.text.y =  element_text(size = 14, face = "bold"))
  
  p <- ggplot(target_sets_fraction_df_m[target_sets_fraction_df_m$Rank %in% 1:thres,], 
              mapping = aes(x = Rank, y = Frac, color = Source)) +  #, color = Source)) + 
    geom_line(aes(linetype=Mutation.Type, alpha = Mutation.Type), size = 1.25) +
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
    scale_color_nejm() + 
    scale_alpha_manual(values=c(0.9,0.9,0.4,0.4)) +
    scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#BC3C29FF", "#0072B5FF"),
                       labels = c(names(target_sets_list), rep("", times = length(names(target_sets_list))))) + 
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
tp53_select_target_set_list <- list("Curated_Fischer_2017" = tp53_curated_targets, "KEGG.hsa04115" = tp53_kegg_pathway_genes,
                                    "TRRUST" = tp53_trrust_targets, "DoRothEA" = tp53_dorothea_targets)
tp53_select_target_set_list <- list("Curated_Fischer_2017" = tp53_curated_targets,
                                    "DoRothEA" = tp53_dorothea_targets)
tp53_select_target_set_list_cna <- list("STRING" = tp53_string_nw_targs,
                                        "TRRUST" = tp53_trrust_targets_upstr)
pik3ca_select_target_set_list <- list("STRING" = pik3ca_string_nw_targs_top500, "KEGG.hsa04151" = pik3ca_kegg_pathway_genes)


# Call function
plot_combined_enrichment(allgenes_p53, tp53_target_set_list, "TP53", 500, 50, NA)
plot_combined_enrichment(allgenes_pik3ca, pik3ca_target_set_list, "PIK3CA", 500, 50, NA)


# Same as above, but one target set and multiple cancer types
plot_combined_enrichment_percancer <- function(master_df_list, target_set, goi, thres, cutout_thres) {
  
  master_df_list <- lapply(master_df_list, function(master_df) 
    master_df[order(master_df$q.value),])
  
  # For each source, create a 0 and 1 vector for each  target hit to indicate if it 
  # is a a hit in that source, and use that to create a fraction of genes at rank DF
  target_set_fraction_dfs <- lapply(1:length(master_df_list), function(i) {
    master_df <- master_df_list[[i]]
    target_set <- setdiff(target_set, setdiff(target_set, master_df$T_k.name))
    vect <- unlist(lapply(1:nrow(master_df), function(i) {
      x <- master_df[i, 'T_k.name']
      v <- 0
      if(x %fin% target_set) {v <- 1}
      return(v)
      #return(ifelse(x %fin% target_set, 1, 0))
    }))
    vect <- list(vect)
    names(vect) <- names(master_df_list)[i]
    target_set_frac_df <- get_frac_of_genes_at_rank(vect, master_df)
    target_set_frac_df_m <- melt(target_set_frac_df, "Rank")
    colnames(target_set_frac_df_m) <- c("Rank", "Source", "Frac")
    target_set_frac_df_m$Mutation.Type <- rep("Nonsynonymous", times = nrow(target_set_frac_df_m))
    print(head(target_set_frac_df_m))
    return(target_set_frac_df_m)
  })
  target_set_fraction_df <- do.call(rbind, target_set_fraction_dfs)

  print(head(target_set_fraction_df))
  print(unique(target_set_fraction_df$Source))
  
  # Create inset plot
  p2 <- ggplot(target_set_fraction_df[target_set_fraction_df$Rank %in% 1:cutout_thres,], 
               mapping = aes(x = Rank, y = Frac, color = Source)) + 
    geom_line(aes(linetype=Mutation.Type), size = 1.25) +
    #geom_point(size = 1) + 
    scale_x_continuous(limits = c(1,cutout_thres)) + 
    scale_color_nejm() + 
    #scale_alpha_manual(values=c(0.9,0.9,0.4,0.4)) +
    #scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#BC3C29FF", "#0072B5FF")) + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.minor = element_blank(), axis.text.x = element_text(size = 14, face = "bold"), 
          axis.text.y =  element_text(size = 14, face = "bold"))
  
  p <- ggplot(target_set_fraction_df[target_set_fraction_df$Rank %in% 1:thres,], 
              mapping = aes(x = Rank, y = Frac, color = Source)) +  #, color = Source)) + 
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
    scale_color_nejm() + 
    #scale_alpha_manual(values=c(0.9,0.9,0.4,0.4)) +
    #scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#BC3C29FF", "#0072B5FF"),
                       #labels = c(names(target_sets_list), rep("", times = length(names(target_sets_list))))) + 
    annotation_custom(ggplotGrob(p2), xmin =(thres-round(thres/1.5)), xmax = thres, ymin = 0.30, ymax = 1.0) +
    geom_rect(data=target_set_fraction_df[1:cutout_thres,], 
              aes(xmin = 1, xmax = cutout_thres, ymin = 0, ymax = 1, fill = "gray"),
              color = NA, fill = alpha("gray", .01))
  
  #vp <- viewport(width = 0.4, height = 0.4, x = 0.8, y = 0.2)
  #print(p)
  #print(p, vp = vp)
  
  return(p)
}

tp53_top_cancer_types <- c("LGG", "LIHC", "LUAD", "BRCA", "MESO", "CESC", "PRAD", 
                           "HNSC", "COAD", "BLCA", "UCEC", "KICH", "ACC")
top_drivers_0.05_tp53 <- lapply(top_drivers_0.05, function(x) x[x$R_i.name == "TP53",])
plot_combined_enrichment_percancer(top_drivers_0.05_tp53[names(top_drivers_0.05_tp53) %in% tp53_top_cancer_types],
                                   tp53_curated_targets, "TP53", 500, 50)

# PART B: GENE SET ENRICHMENT ANALYSIS

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
# Call if needed

# DOCUMENTATION FOR THIS PACKAGE: https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html

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
  regprot.gsea.mkegg <- list()
  regprot.gsea.wp <- list()
  
  # Loop through all the regulatory proteins
  for (regprot in unique(results_table$R_i.name)) {
    
    # Get all of this protein's targets and their associated Betas, sorted in 
    # descending order
    res_table_sub <- results_table[results_table$R_i.name == regprot, c('T_k', 'T_k.name', 'estimate', 'p.value')]
    print(head(bitr(res_table_sub$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb=org.Hs.eg.db, drop = TRUE)))
    mapping <- as.data.frame(bitr(res_table_sub$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                                  OrgDb=org.Hs.eg.db, drop = TRUE))
    mapping_kegg <- as.data.frame(bitr_kegg(res_table_sub$T_k, fromType = "uniprot", toType = "kegg",
                                            drop = TRUE, organism = "hsa"))
    colnames(mapping) <- c("T_k.name", "T_k.entrez")
    colnames(mapping_kegg) <- c("T_k", "T_k.kegg")
    res_table_sub <- merge(res_table_sub, mapping, all=TRUE, by="T_k.name")
    res_table_sub_kegg <- merge(res_table_sub, mapping_kegg, all=TRUE, by="T_k")
    if(sort_by == "estimate") {
      res_table_sub <- res_table_sub[order(res_table_sub$estimate, decreasing = TRUE),]
      res_table_sub_kegg <- res_table_sub_kegg[order(res_table_sub_kegg$estimate, decreasing = TRUE),]
    } else {
      res_table_sub$negLogPval <- unlist(lapply(1:nrow(res_table_sub), function(i) {
        pval <- res_table_sub$p.value[i]
        estimate <- res_table_sub$estimate[i]
        est_sign <- ifelse(estimate > 0, 1, -1)
        return((-log(pval)) * est_sign)
      }))
      res_table_sub <- res_table_sub[order(res_table_sub$negLogPval, decreasing = TRUE),]
      
      res_table_sub_kegg$negLogPval <- unlist(lapply(1:nrow(res_table_sub_kegg), function(i) {
        pval <- res_table_sub_kegg$p.value[i]
        estimate <- res_table_sub_kegg$estimate[i]
        est_sign <- ifelse(estimate > 0, 1, -1)
        return((-log(pval)) * est_sign)
      }))
      res_table_sub_kegg <- res_table_sub_kegg[order(res_table_sub_kegg$negLogPval, decreasing = TRUE),]
      res_table_sub_kegg <- na.omit(res_table_sub_kegg)
    }
    
    plotdir <- paste0(output_path, "Linear Model/GSEA/")    
    suppressWarnings(dir.create(plotdir, recursive = TRUE))
    
    regprotBetaScores <- res_table_sub$estimate
    regprotBetaScores_kegg <- res_table_sub_kegg$estimate
    
    if(sort_by == "p.value") {
      regprotBetaScores <- res_table_sub$negLogPval
      regprotBetaScores_kegg <- res_table_sub_kegg$negLogPval
    }
    names(regprotBetaScores) <- res_table_sub$T_k.entrez
    names(regprotBetaScores_kegg) <- res_table_sub_kegg$T_k.kegg
    
    gse.rp <- gsePathway(regprotBetaScores, pvalueCutoff = 1, pAdjustMethod = "BH")
    gse.rp <- setReadable(gse.rp, org.Hs.eg.db)
    gse.ncg <- gseNCG(regprotBetaScores, pvalueCutoff = 1, pAdjustMethod = "BH")
    gse.ncg <- setReadable(gse.ncg, org.Hs.eg.db)
    gse.go <- gseGO(regprotBetaScores, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, 
                    pAdjustMethod = "BH", keyType = "ENTREZID", ont = "BP")  # BP = biological process
    gse.go <- setReadable(gse.go, org.Hs.eg.db)
    gse.kegg <- gseKEGG(regprotBetaScores_kegg, pvalueCutoff = 1, 
                        pAdjustMethod = "BH", keyType = "kegg")
    #gse.kegg <- setReadable(gse.kegg, org.Hs.eg.db, keyType = "kegg")
    gse.mkegg <- gseMKEGG(regprotBetaScores_kegg, pvalueCutoff = 1, 
                          pAdjustMethod = "BH", keyType = "kegg")
    gse.wp <- gseWP(regprotBetaScores, organism = "Homo sapiens", pvalueCutoff = 1, 
                    pAdjustMethod = "BH")
    
    regprot.gsea.rp[[regprot]] <- gse.rp
    regprot.gsea.ncg[[regprot]] <- gse.ncg
    regprot.gsea.go[[regprot]] <- gse.go
    regprot.gsea.kegg[[regprot]] <- gse.kegg
    regprot.gsea.mkegg[[regprot]] <- gse.mkegg
    regprot.gsea.wp[[regprot]] <- gse.wp
    res_tab <- NA
    
    
    # plot GSEA for top N upregulated and top N downregulated terms
    for (gst in c('rp','ncg', 'go', 'kegg', 'mkegg', 'wp')) {
      #for (gst in c('rp','ncg')) {
      if (gst == 'rp') res <- gse.rp
      else if (gst == 'kegg') res <- gse.kegg
      else if (gst == 'mkegg') res <- gse.mkegg
      else if (gst == 'go') res <- gse.go
      else if (gst == 'wp') res <- gse.wp
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
      if(gst == "kegg") {res_tab <- res}
    }
  }
  return(res_tab)
}

results_gsea <- perform_gsea(master_df_mut_corrected, n = 5, "p.value", output_path)
results_gsea <-perform_gsea(master_df_cna_corrected, n = 5, "p.value", output_path)

results_gsea <- results_gsea@result

kras_top_drivers_0.05 <- lapply(top_drivers_0.05, function(x) {
  if("KRAS" %fin% x$R_i.name) {
     return(x[x$R_i.name == "KRAS",])
  } else {return(NA)}
})
kras_top_drivers_0.05 <- kras_top_drivers_0.05[!is.na(kras_top_drivers_0.05)]

kras_gsea_results <- lapply(kras_top_drivers_0.05, function(x) 
  perform_gsea(x, n = 5, "p.value", output_path))
names(kras_gsea_results) <- names(kras_top_drivers_0.05)

lapply(1:length(kras_gsea_results), function(i) {
  name <- names(kras_gsea_results)[i]
  write.csv(kras_gsea_results[[i]], paste0("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/GSEA_Results/KRAS/",
                                         paste0(name, "_gsea_mkegg.csv")))
})


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
    theme(legend.position = "none", axis.title = element_text(face = "bold", size = 14), 
          axis.text = element_text(face = "bold", size = 12)) + 
    ylab("-log10(adj.pvalue)") + xlab(paste("Enriched", paste(db, "Pathways"))) 
    #+ scale_fill_manual(values = color_palette)
    #+ scale_fill_nejm() 
  
  p
}

# Call function
create_gsea_barchart(results_gsea, 25, "GO")

# Same as above, but for multiple drivers
create_gsea_barchart_multGenes <- function(list_of_results_gsea, n, qval_thres, 
                                           sort_by, db) {
  # Subset the results table based on the number of pathways
  input_dfs <- lapply(1:length(list_of_results_gsea), function(i) {
    results_gsea <- list_of_results_gsea[[i]]
    driver <- names(list_of_results_gsea)[i]
    
    if(!is.na(qval_thres)) {results_gsea <- results_gsea[results_gsea$qvalue < 
                                                           qval_thres, ]}
    
    if(sort_by == "ES") {results_gsea <- results_gsea[order(-abs(results_gsea$enrichmentScore)),]}
    else if (sort_by == "neglog10(qval)*dir") {
      dir <- unlist(lapply(results_gsea$enrichmentScore, function(x) ifelse(x>0, 1, -1)))
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
      
      # If desired, sort pathways within each q-value bucket by leading edge signal
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
    print(head(results_gsea_sub))
    
    # Convert the p-values to -log10(pvalue)
    negLog10_pvalues <- -log10(results_gsea_sub$qvalue)
    
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
  
  #color_palette <- rep(pal_nejm("default")(8), times = ceiling(n/8))
  p <- ggplot(input_df, aes(x =  Enriched.Pathway.Driver, #reorder(Enriched.Pathway.Driver, abs(ES)), 
                            y = ES, fill = Driver)) +
    geom_col(width = 0.7, color = "black") + coord_flip() + theme_minimal() + 
    theme(legend.position = "none", axis.title = element_text(face = "bold", size = 14), 
          axis.text = element_text(face = "bold", size = 10), 
          strip.text.x = element_text(face="bold", size=14, margin = margin(.2, 0, .2, 0, "cm"))) + 
    xlab(paste("Enriched", paste(db, "Pathways"))) + ylab("Enrichment Score") + #ylab("-log10(qval)") +
    scale_fill_manual(values = c("#FFDC91FF","#20854EFF","#BC3C29FF","#0072B5FF")) + 
    scale_x_discrete(labels=roles) +
    facet_wrap(~factor(Driver, levels = c("TP53", "PIK3CA", "KRAS", "IDH1")), 
               ncol = 1, scales = "free_y") 
  #+ scale_fill_nejm() 
  
  p
}

tp53_pc_gsea_kegg <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/GSEA_Results/TP53/PC_gsea_kegg_results.csv", 
                              header = T, check.names = F)
pik3ca_pc_gsea_kegg <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/GSEA_Results/PIK3CA/PC_gsea_kegg_results.csv", 
                                header = T, check.names = F)
pik3ca_pc_gsea_kegg[pik3ca_pc_gsea_kegg$Description == "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate", 'Description'] <- "Glycosaminoglycan biosynthesis"
kras_pc_gsea_kegg <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/GSEA_Results/KRAS/PC_gsea_kegg_results.csv", 
                              header = T, check.names = F)
idh1_pc_gsea_kegg <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/GSEA_Results/IDH1/PC_gsea_kegg_results.csv", 
                              header = T, check.names = F)

list_of_results_gsea <- list("TP53" = tp53_pc_gsea_kegg, "PIK3CA" = pik3ca_pc_gsea_kegg, 
                             "KRAS" = kras_pc_gsea_kegg, "IDH1" = idh1_pc_gsea_kegg)

create_gsea_barchart_multGenes(list_of_results_gsea, 6, 0.01, "ES", "KEGG")


# Option to merge together pathways that share a Jaccard index of > threshold
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

merge_pathways <- function(gsea_res, jaccard_thres, order_by) {
  
  # Make sure we are ordered by correct quantity
  if(order_by == "qvalue") {
    gsea_res <- gsea_res[order(gsea_res$qvalue, decreasing=F),]
  } else if (order_by == "enrichmentScore") {
    gsea_res <- gsea_res[order(abs(gsea_res$enrichmentScore), decreasing=T),]
  } else {
    print("Only implemented for 'qvalue' and '-log(qval)*dir', and 'enrichmentScore', 
          please try again with one of these order_by quantities.")
  }
  
  rows_to_keep <- c(1)
  pathway_member_list <- list(strsplit(gsea_res[1, 'core_enrichment'], "/", fixed = T))
  
  for (i in 2:nrow(gsea_res)) {
    curr_members <- strsplit(gsea_res[i, 'core_enrichment'], "/", fixed = T)
    similarity_bool <- F
    for (j in 1:length(pathway_member_list)) {
      jacc <- jaccard(pathway_member_list[[j]], curr_members[[1]]) 
      print(jacc)
      if(jacc > jaccard_thres) {
        similarity_bool <- T
        break
      }
    }
    if(!similarity_bool) {
      rows_to_keep <- c(rows_to_keep, i)
      pathway_member_list <- append(pathway_member_list, curr_members)
    }
  }
  
  return(gsea_res[rows_to_keep, ])
}

tp53_pc_gsea_kegg_merged_0.1 <- merge_pathways(tp53_pc_gsea_kegg, 0.1)  # from 337 to 110 pathways
pik3ca_pc_gsea_kegg_merged_0.1 <- merge_pathways(pik3ca_pc_gsea_kegg, 0.1)  # from 337 to 108 pathways
kras_pc_gsea_kegg_merged_0.1 <- merge_pathways(kras_pc_gsea_kegg, 0.1) # from 337 to 107 pathways
idh1_pc_gsea_kegg_merged_0.1 <- merge_pathways(idh1_pc_gsea_kegg, 0.1) # from 337 to 107 pathways

list_of_results_gsea_merged <- list("TP53" = tp53_pc_gsea_kegg_merged_0.1, "PIK3CA" = pik3ca_pc_gsea_kegg_merged_0.1, 
                             "KRAS" = kras_pc_gsea_kegg_merged_0.1, "IDH1" = idh1_pc_gsea_kegg_merged_0.1)
create_gsea_barchart_multGenes(list_of_results_gsea_merged, 6, 0.01, "ES", "KEGG")

#' Alternate version of this using ReactomePA's built-in semantic similarity analysis
#' Link to documentation: https://yulab-smu.top/biomedical-knowledge-mining-book/GOSemSim.html
#' @param gsea_res a GSEA file from ReactomePA
#' @param sim_thres a threshold for similarity metric
#' @param hsGO a hsGO file for goSim
#' @param measure a method of semantic similarity measurement, either 
#' 'Wang' or 'Jiang' 
#' @param order_by either 'q.value' or 'enrichmentScore', what to order by
merge_pathways_semantic_sim <- function(gsea_res, sim_thres, hsGO, 
                                        measure, order_by) {
  
  # Make sure we are ordered by correct quantity
  if(order_by == "q.value") {
    gsea_res <- gsea_res[order(gsea_res$qvalue, decreasing=F),]
  } else if (order_by == "enrichmentScore") {
    gsea_res <- gsea_res[order(abs(gsea_res$enrichmentScore), decreasing=T),]
  } else {
    print("Only implemented for 'qvalue' and '-log(qval)*dir', and 'enrichmentScore', 
          please try again with one of these order_by quantities.")
  }
  
  rows_to_keep <- c(1)
  go_term_list <- gsea_res[1, 'ID']
  go_term_list_sizes <- gsea_res[1, 'setSize']

  for (i in 2:nrow(gsea_res)) {
    curr_go <- gsea_res[i, 'ID']
    size_go <- gsea_res[i, 'setSize']
    similarity_bool <- F
    for (j in 1:length(go_term_list)) {
      sim <- goSim(go_term_list[[j]], curr_go, semData=hsGO, measure=measure) 
      print(sim)
      if(is.na(sim)) {sim <- 0}
      if(sim > sim_thres) {
        similarity_bool <- T
        break
      }
    }
    if(!similarity_bool) {
      rows_to_keep <- c(rows_to_keep, i)
      go_term_list <- append(go_term_list, curr_go)
      go_term_list_sizes <- append(go_term_list_sizes, size_go)
    }
  }
  
  return(gsea_res[rows_to_keep, ])
}

#' Generate pairwise similarity matrix
#' @param gsea_res a GSEA file from ReactomePA
#' @param hsGO a hsGO file for goSim
#' @param measure a method of semantic similarity measurement, either 
#' 'Wang' or 'Jiang' 
#' @param order_by either 'q.value' or 'enrichmentScore', what to order by
#' @param thres q-value threshold for significance of given pathway
generate_pw_go_sim_mat <- function(gsea_res, hsGO, measure, order_by, thres) {
  
  gsea_res <- gsea_res[gsea_res$qvalue < thres,]
  
  # Make sure we are ordered by correct quantity
  if(order_by == "q.value") {
    gsea_res <- gsea_res[order(gsea_res$qvalue, decreasing=F),]
  } else if (order_by == "enrichmentScore") {
    gsea_res <- gsea_res[order(abs(gsea_res$enrichmentScore), decreasing=T),]
  } else {
    print("Only implemented for 'qvalue' and '-log(qval)*dir', and 'enrichmentScore', 
          please try again with one of these order_by quantities.")
  }
  
  output_mat <- matrix(nrow=nrow(gsea_res), ncol=nrow(gsea_res))
  print(dim(output_mat))
  
  for(i in 1:nrow(gsea_res)) {
    for (j in 1:nrow(gsea_res)) {
      if (i==j) {output_mat[i,j] <- 1}
      else {
        go_pw1 <- gsea_res[i, 'ID']
        go_pw2 <- gsea_res[j, 'ID']
        sim <- goSim(go_pw1, go_pw2, semData=hsGO, measure=measure) 
        print(sim)
        output_mat[i,j] <- sim
      }
    }
  }
  
  output_mat <- as.data.frame(output_mat)
  rownames(output_mat) <- gsea_res$ID
  colnames(output_mat) <- gsea_res$ID
  return(output_mat)
}

hsGO <- godata('org.Hs.eg.db', ont="BP")  # BP = biological process

gsea_go_pc_tp53_merged <- merge_pathways_semantic_sim(gsea_go_pc_tp53, 0.4, hsGO,
                                                      "Wang", 'q.value')

idh1_go_sim_mat <- generate_pw_go_sim_mat(gsea_go_pc_idh1, hsGO, "Wang", 
                                          "q.value", 0.1)
# Make the top half of the triangle (which is identical) NA
for (i in 1:nrow(go_sim_mat)) {
  go_sim_mat[i, i:ncol(go_sim_mat)] <- NA
}


#' Merge pathways based on this similarity matrix, from most to least similar,
#' not dependent on q-value ordering
#' @param go_sim_mat a pathway similarity matrix between all significant GO PWs,
#' with the top half of the identity NA'd
#' @param thres a similarity threshold, beyond which we will no longer merge
#' two pathways (similarity ranges from 0-1) #TODO: to we want to boost the rank
#' of pathways that have been merged with several other pathways?
#' @param gsea_df a GO GSEA data frame with the new, merged PWs
#' @param metric_to_keep either "setSize" to indicate we are choosing the 
#' pathway with the smaller set size, "qvalue" to indicate we are choosing the
#' pathway with the smaller q-value, or "enrichmentScore", to indicate we are
#' choosing the pathway with the larger enrichment score
merge_pathways_by_sim_mat <- function(go_sim_mat, thres, gsea_df, 
                                      metric_to_keep) {
  # Sort the values from most to least similarity
  go_sim_mat_values <- unlist(lapply(1:ncol(go_sim_mat), function(i) {
    vals <- as.numeric(unlist(go_sim_mat[,i]))
    return(vals[!is.na(vals)])
  }))
  go_sim_mat_values <- go_sim_mat_values[order(go_sim_mat_values, 
                                               decreasing = T)]
  # Limit to those with similarity above the given threshold
  go_sim_mat_values <- go_sim_mat_values[go_sim_mat_values >= thres]
  
  # Create a version of GSEA DF we will subset
  gsea_df_sub <- gsea_df
  
  # For each of the remaining similarity values, identify the GO pathways they
  # are associated with
  for (v in go_sim_mat_values) {
    print(v)
    mat <- which(go_sim_mat == v, arr.ind = T)
    rown <- rownames(go_sim_mat)[mat[,1]]
    coln <- colnames(go_sim_mat)[mat[,2]]
    
    for(i in 1:length(rown)) {
      go1 <- rown[i]
      go2 <- coln[i]
      
      if((go1 %fin% gsea_df_sub$ID) & (go2 %fin% gsea_df_sub$ID)) {
        # Get which GO term has the larger set size OR a larger
        # enrichment score OR a smaller q-value
        if(metric_to_keep == "setSize") {
          setsize1 <- gsea_df_sub[gsea_df_sub$ID == go1, 'setSize']
          setsize2 <- gsea_df_sub[gsea_df_sub$ID == go2, 'setSize']
          if(setsize1 < setsize2) {
            gsea_df_sub <- gsea_df_sub[!(gsea_df_sub$ID == go1),]
          } else {
            gsea_df_sub <- gsea_df_sub[!(gsea_df_sub$ID == go2),]
          }
        } else if (metric_to_keep == "qvalue") {
          qval1 <- gsea_df[gsea_df$ID == go1, 'qvalue']
          qval2 <- gsea_df[gsea_df$ID == go2, 'qvalue']
          if(qval1 > qval2) {
            gsea_df_sub <- gsea_df_sub[!(gsea_df_sub$ID == go1),]
          } else if (qval1 < qval2) {
            gsea_df_sub <- gsea_df_sub[!(gsea_df_sub$ID == go2),]
          } else {
            # Q-values are the same, use ES
            es1 <- abs(gsea_df[gsea_df$ID == go1, 'enrichmentScore'])
            es2 <- abs(gsea_df[gsea_df$ID == go2, 'enrichmentScore'])
            if(es1 > es2) {
              gsea_df_sub <- gsea_df_sub[!(gsea_df_sub$ID == go2),]
            }
            else {
              gsea_df_sub <- gsea_df_sub[!(gsea_df_sub$ID == go1),]
            }
          }
        } else if (metric_to_keep == "enrichmentScore") {
          es1 <- abs(gsea_df[gsea_df$ID == go1, 'enrichmentScore'])
          es2 <- abs(gsea_df[gsea_df$ID == go2, 'enrichmentScore'])
          if(es1 > es2) {
            gsea_df_sub <- gsea_df_sub[!(gsea_df_sub$ID == go2),]
          }
          else {
            gsea_df_sub <- gsea_df_sub[!(gsea_df_sub$ID == go1),]
          }
        } else {
          print("Only implemented for 'setSize', 'qvalue', and 'enrichmentScore'. 
              Using q-value.")
        }
      }
    }
  }
  print(head(gsea_df_sub))
  return(gsea_df_sub)
}

gsea_go_pc_kras_sub_setSize0.7 <- merge_pathways_by_sim_mat(kras_go_sim_mat2, 0.7, 
                                                  gsea_go_pc_kras, "setSize")

## USE THE REACTOMEPA compareCluster() FUNCTION TO COMPARE GENE SETS
run_compare_cluster <- function(results_table_list, qval_thres, sort_by) {
  
  list_of_scores <- lapply(1:length(results_table_list), function(i) {
    regprot <- names(results_table_list)[i]
    results_table <- results_table_list[[i]]
    results_table <- results_table[results_table$q.value < qval_thres,]
    
    res_table_sub <- results_table[results_table$R_i.name == regprot, c('T_k', 'T_k.name', 'estimate', 'p.value')]
    print(head(bitr(res_table_sub$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb=org.Hs.eg.db, drop = TRUE)))
    mapping <- as.data.frame(bitr(res_table_sub$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                                  OrgDb=org.Hs.eg.db, drop = TRUE))
    mapping_kegg <- as.data.frame(bitr_kegg(res_table_sub$T_k, fromType = "uniprot", toType = "kegg",
                                            drop = TRUE, organism = "hsa"))
    colnames(mapping) <- c("T_k.name", "T_k.entrez")
    colnames(mapping_kegg) <- c("T_k", "T_k.kegg")
    
    res_table_sub <- merge(res_table_sub, mapping, all=TRUE, by="T_k.name")
    res_table_sub_kegg <- merge(res_table_sub, mapping_kegg, all=TRUE, by="T_k")
    
    if(sort_by == "estimate") {
      res_table_sub <- res_table_sub[order(res_table_sub$estimate, decreasing = TRUE),]
      res_table_sub_kegg <- res_table_sub_kegg[order(res_table_sub_kegg$estimate, decreasing = TRUE),]
      
    } else {
      res_table_sub <- res_table_sub[order(res_table_sub$p.value, decreasing = TRUE),]
      res_table_sub_kegg <- res_table_sub_kegg[order(res_table_sub_kegg$p.value, decreasing = TRUE),]
    }

    return(res_table_sub$T_k.entrez)
    return(res_table_sub_kegg$T_k.kegg)
  })
  names(list_of_scores) <- names(results_table_list)
  
  ck <- compareCluster(geneCluster = list_of_scores, fun = enrichGO, OrgDb = org.Hs.eg.db, ont = "BP")
  ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  ck <- compareCluster(geneCluster = list_of_scores, fun = enrichPathway)
  ck <- setReadable(ck, org.Hs.eg.db)
  
  ck <- compareCluster(geneCluster = list_of_scores, fun = enrichKEGG)

  print(dotplot(ck))

}

results_table_list <- list("TP53" = pc_allgenes[pc_allgenes$R_i.name == "TP53",],
                           "PIK3CA" = pc_allgenes[pc_allgenes$R_i.name == "PIK3CA",],
                           "KRAS" = pc_allgenes[pc_allgenes$R_i.name == "KRAS",],
                           "IDH1" = pc_allgenes[pc_allgenes$R_i.name == "IDH1",])

results_table_list_up <- lapply(results_table_list, function(x) x[x$estimate > 0,])
run_compare_cluster(results_table_list, 0.01, "p.value")


# Version for STRINGdB results:
create_gsea_barchart <- function(results_gsea, n, db) {
  # Subset the results table based on the number of pathways
  results_gsea_sub <- results_gsea[1:n, ]
  results_gsea_sub <- results_gsea_sub[!is.na(results_gsea_sub$description),]
  #results_gsea_sub <- results_gsea_sub[results_gsea_sub$category == db,]
  
  # Convert the p-values to -log10(pvalue)
  negLog10_pvalues <- -log10(results_gsea_sub$fdr)
  
  # Create an input data frame for ggplot
  input_df <- data.frame("Enriched.Pathway" = results_gsea_sub$description,
                         "Neg.Log10.Pval" = negLog10_pvalues)
  
  #color_palette <- rep(pal_nejm("default")(8), times = ceiling(n/8))
  
  p <- ggplot(input_df, aes(x = reorder(Enriched.Pathway, negLog10_pvalues), 
                            y = Neg.Log10.Pval)) +
    geom_col(width = 0.7, fill = "#0072B5FF") + coord_flip() + theme_minimal() + 
    theme(legend.position = "none", axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 12)) + 
    ylab("-log10(adj.pvalue)") + xlab(paste("Enriched", paste(db, "Pathways"))) 
  #+ scale_fill_manual(values = color_palette)
  #+ scale_fill_nejm() 
  
  p
}

# Instantiate the string DB object 
#string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=0, 
                          #input_directory="")

# Create a heat map with the results for significant, overlapping pathways 
#' @param list_of_results_gsea_up a named list of GSEA results from Reactome PA
#' @param list_of_results_gsea_down an optional second named list of GSEA results 
#' from Reactome PA, typically downregulated (same length/ # of drivers)
#' @param sig_thres an adjusted p-value significance threshold
create_sig_pathway_overlap_grid <- function(list_of_results_gsea_up, 
                                            list_of_results_gsea_down, sig_thres) {
  
  # Get all the pathways significant in at least one driver
  all_sig_pws_up <- unique(unlist(lapply(list_of_results_gsea_up, function(df) 
    df[df$fdr < sig_thres, 'description'])))
  all_sig_pws_down <- c()
  if(!(length(list_of_results_gsea_down) == 0)) {
    all_sig_pws_down <- unique(unlist(lapply(list_of_results_gsea_down, function(df) 
      df[df$fdr < sig_thres, 'description'])))
  }
  
  # Set up the driver x pathway matrix that we will fill and use to plot
  all_pws <- unique(c(all_sig_pws_up, all_sig_pws_down))
  input_df <- data.frame(matrix(nrow = length(all_pws), 
                                ncol = length(list_of_results_gsea_up)))
  rownames(input_df) <- all_pws
  colnames(input_df) <- names(list_of_results_gsea_up)
  
  for(i in 1:length(list_of_results_gsea_up)) {
    driver <- names(list_of_results_gsea_up)[i]
    gsea_res_up <- list_of_results_gsea_up[[i]]
    sig_hits_up <- gsea_res_up[gsea_res_up$fdr < sig_thres, 'description']
    
    sig_hits_down <- NA
    if(!(length(all_sig_pws_down) == 0)) {
      gsea_res_down <- list_of_results_gsea_down[[i]]
      sig_hits_down <- gsea_res_down[gsea_res_down$fdr < sig_thres, 'description']
    }
    
    input_df[,driver] <- unlist(lapply(all_pws, function(pw) {
      if(pw %fin% sig_hits_up) {return(i)} 
      if(pw %fin% sig_hits_down) {return(-i)} 
      else{return(0)}
    })) 
  }

  # Limit to just pathways with a 1 in more than one driver
  #input_df <- input_df[which(rowSums(input_df) > 1),]
  #input_df <- input_df[rowSums(input_df > 0) > 1,]
  #input_df <- apply(input_df, MARGIN = 2, function(y) as.integer(y))
  print(input_df)
  
  # Make heatmap
  pheatmap(input_df, angle_col = "45", main = "",
           #color=c("#E18727FF", "white", "#0072B5FF"),
           c("white", "#0072B5FF", "#BC3C29FF", "#20854EFF", "#E18727FF"),
           #color=c("-1" = "#0072B5FF", "-2" = "#BC3C29FF", "-3" = "#20854EFF", 
                   #"-4" = "#E18727FF", "1" = "#0072B599", "2" = "#BC3C2999", 
                   #"3" = "#20854E99", "4" = "#E1872799", "0" = "white"), 
           fontsize_row = 9, fontsize_col = 14, fontface_col = "bold", 
           legend = FALSE, cellwidth=20, cellheight = 9,
           cluster_rows = F, cluster_cols = F) 
  
  return(input_df)
}

pathway_matrix <- create_sig_pathway_overlap_grid(list("TP53" = results_gsea_tp53, "PIK3CA" = results_gsea_pik3ca, 
                                     "KRAS" = results_gsea_kras, "IDH1" = results_gsea_idh1), 0.2)

pathway_matrix <- create_sig_pathway_overlap_grid(list("TP53" = stringdb_enrichment_tp53_0.1_up[stringdb_enrichment_tp53_0.1_up$category == "RCTM",], "PIK3CA" = stringdb_enrichment_pik3ca_0.1_up[stringdb_enrichment_pik3ca_0.1_up$category == "RCTM",], 
                                                       "KRAS" = stringdb_enrichment_kras_0.1_up[stringdb_enrichment_kras_0.1_up$category == "RCTM",], "IDH1" = stringdb_enrichment_idh1_0.1_up[stringdb_enrichment_idh1_0.1_up$category == "RCTM",]),
                                                  list("TP53" = stringdb_enrichment_tp53_0.1_down[stringdb_enrichment_tp53_0.1_down$category == "RCTM",], "PIK3CA" = stringdb_enrichment_pik3ca_0.1_down[stringdb_enrichment_pik3ca_0.1_down$category == "RCTM",], 
                                                       "KRAS" = stringdb_enrichment_kras_0.1_down[stringdb_enrichment_kras_0.1_down$category == "RCTM",], "IDH1" = stringdb_enrichment_idh1_0.1_down[stringdb_enrichment_idh1_0.1_down$category == "RCTM",]), 
                                                  0.05)

create_sig_pathway_overlap_grid(list("TP53" = stringdb_enrichment_tp53_up_new[stringdb_enrichment_tp53_up_new$category == "KEGG",], "PIK3CA" = stringdb_enrichment_pik3ca_up_new[stringdb_enrichment_pik3ca_up_new$category == "KEGG",], 
                                     "KRAS" = stringdb_enrichment_kras_up_new[stringdb_enrichment_kras_up_new$category == "KEGG",], "IDH1" = stringdb_enrichment_idh1_up_new[stringdb_enrichment_idh1_up_new$category == "KEGG",]), list(), 
                                0.2)


drivers <- c("TP53", "PIK3CA", "KRAS", "IDH1")
all_sig_pws <- unique(c(tp53_jaccard0.1_pathways, pik3ca_jaccard0.1_pathways, 
                        kras_jaccard0.1_pathways, idh1_jaccard0.1_pathways))
driver_enrich_pws <- list(tp53_jaccard0.1_pathways, pik3ca_jaccard0.1_pathways, 
                            kras_jaccard0.1_pathways, idh1_jaccard0.1_pathways)
names(driver_enrich_pws) <- drivers

input_df <- data.frame(matrix(nrow = length(all_sig_pws), ncol = 4))
rownames(input_df) <- all_sig_pws
colnames(input_df) <- drivers
for(i in 1:length(drivers)) {
  driver <- drivers[i]
  gsea_res <- driver_enrich_pws[[which(names(driver_enrich_pws) == driver)]]
  input_df[,driver] <- unlist(lapply(all_sig_pws, function(pw) 
    ifelse(pw %fin% gsea_res, i, 0)))  # could change this to return the actual directional pval
}
#input_df <- input_df[which(rowSums(input_df) > 1),]
print(input_df)

# Make heatmap
pheatmap(input_df, angle_col = "45", main = "",
         color=c("white", "#0072B5FF", "#BC3C29FF", "#20854EFF", "#E18727FF"), fontsize_row = 9, fontsize_col = 14,
         fontface_col = "bold", legend = FALSE, cellwidth=20, cellheight = 10) 


# OPTION 3: INDIVIDUAL BARPLOTS
# Create a combined bar plot with the results for significant pathways 
#' @param list_of_results_gsea a named list of GSEA results from STRINGdb or ReactomePA
#' @param n number of pathways to include per driver
#' @param source either "STRINGdb" or "ReactomePA", to indicate which column names we are using
#' @param label the name of the pathway source, for labeling barplot
create_sig_pathway_barplot <- function(list_of_results_gsea, n, source, label) {

  input_df_full <- data.frame()
  
  # Create plots for each
  for (i in 1:length(list_of_results_gsea)) {
    
    results_gsea <- list_of_results_gsea[[i]]

    # Subset the results table based on the number of pathways
    results_gsea_sub <- results_gsea[1:n, ]
    
    if(source == "STRINGdb") {
      results_gsea_sub <- results_gsea_sub[!is.na(results_gsea_sub$description),]
      # Convert the p-values to -log10(pvalue)
      negLog10_pvalues <- -log10(results_gsea_sub$fdr)
      # Create an input data frame for ggplot
      input_df <- data.frame("Enriched.Pathway" = results_gsea_sub$description,
                             "Neg.Log10.Pval" = as.numeric(negLog10_pvalues), 
                             "Driver" = names(list_of_results_gsea)[i])
      input_df_full <- rbind(input_df_full, input_df)
      
    }
    else if (source == "ReactomePA") {
      results_gsea_sub <- results_gsea_sub[!is.na(results_gsea_sub$Description),]
      # Convert the p-values to -log10(pvalue)
      negLog10_pvalues <- -log10(results_gsea_sub$qvalue)
      # Create an input data frame for ggplot
      input_df <- data.frame("Enriched.Pathway" = results_gsea_sub$Description,
                             "Neg.Log10.Pval" = as.numeric(negLog10_pvalues), 
                             "Driver" = names(list_of_results_gsea)[i])
      input_df_full <- rbind(input_df_full, input_df)
    }
    else {
      print("Only implemented for STRINGdb or ReactomePA. Please retry with a valid source.")
      return(NA)
    }
    #results_gsea_sub <- results_gsea[!is.na(results_gsea$description),]

  }
  ggplot(input_df_full, aes(x = reorder(Enriched.Pathway, Neg.Log10.Pval), 
                       y = Neg.Log10.Pval, fill = Driver)) +
    geom_col() + coord_flip() + theme_minimal() + 
    scale_fill_manual(values = c("#E18727FF", "#20854EFF", "#BC3C29FF", "#0072B5FF")) + 
    theme(legend.position = "none", axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 11), strip.text.x = element_text(size = 14, face = "bold")) + 
    facet_grid(~fct_rev(Driver), scales = "free_y", space = "free_y") + #ncol = 1,
    ylab("-log10(qvalue)") + xlab(paste("Enriched", label)) 
}

create_sig_pathway_barplot(list("TP53" = stringdb_enrichment_tp53_up_new[stringdb_enrichment_tp53_up_new$category == "WikiPathways",],
                                "PIK3CA" = stringdb_enrichment_pik3ca_up_new[stringdb_enrichment_pik3ca_up_new$category == "WikiPathways",],
                                "KRAS" = stringdb_enrichment_kras_up_new[stringdb_enrichment_kras_up_new$category == "WikiPathways",],
                                "IDH1" = stringdb_enrichment_idh1_up_new[stringdb_enrichment_idh1_up_new$category == "WikiPathways",]),
                           n = 10, "STRINGdb", "WikiPathways")

create_sig_pathway_barplot(list("TP53" = tp53_pc_gsea_kegg, "PIK3CA" = pik3ca_pc_gsea_kegg,
                                "KRAS" = kras_pc_gsea_kegg, "IDH1" = idh1_pc_gsea_kegg),
                           n = 10, "ReactomePA", "KEGG Pathways")


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
  
  pvals_real <- real_master_df[, c('p.value', 'T_k.name')]
  pvals_random <- random_master_df[, c('p.value', 'T_k.name')]
  if(!is.na(goi)) {
    pvals_real <- real_master_df[real_master_df$R_i.name == goi, c('p.value', 'T_k.name')]
    pvals_random <- random_master_df[random_master_df$R_i.name == goi, c('p.value', 'T_k.name')]
  } 
  
  pvals <- merge(pvals_real, pvals_random, by = "T_k.name")
  colnames(pvals)[colnames(pvals) == "p.value.x"] <- "Real"
  colnames(pvals)[colnames(pvals) == "p.value.y"] <- "Random"
  pvals_m <- melt(pvals)
  
  print(head(pvals_m))
  
  # Get the p-value corresponding to q-value threshold
  max_qval_under_thres <- max(real_master_df[real_master_df$q.value < qval_thres, 'q.value'])
  print(max_qval_under_thres)
  pval_thres <- real_master_df[real_master_df$q.value == max_qval_under_thres, 'p.value']
  
  p <- ggplot(pvals_m, aes(x = value, fill = variable)) + 
    geom_histogram(alpha = 0.4, position = "identity", bins = 50) + 
    theme_minimal() + xlab("P-Value") + ylab("Frequency") + 
    labs(fill = "") + theme(axis.title = element_text(face = "bold", size = 14),
                                              axis.text = element_text(face = "bold", size = 12),
                                              legend.text = element_text(size=12), legend.position = c(0.8, 0.8),
                                              legend.title = element_text(face = "bold", size = 14)) +
    #geom_vline(xintercept = pval_thres, linetype='dotted', linewidth = 1) + 
    scale_fill_nejm() 
    #annotate("text", y = 2000, x = pval_thres, 
             #label = paste("q <", paste(qval_thres, paste0("; p <", round(pval_thres, digits = 4)))), 
             #hjust = -0.05, size = 5, fontface = "bold")
  p
}

# Call function
plot_pvalue_histograms_real_and_random(allgenes_p53, allgenes_p53_random, "TP53", 0.1)


# A non-overlaid version
ggplot(allgenes_p53, aes(x = p.value)) + 
       geom_histogram(alpha = 1, position = "identity", bins = 50, fill = "#0072B5FF") + 
       theme_minimal() + xlab("P-Value") + ylab("Frequency") + 
       labs(fill = "") + theme(axis.title = element_text(face = "bold", size = 14),
                                                             axis.text = element_text(face = "bold", size = 12),
                                                             legend.text = element_text(size=12),
                                                             legend.title = element_text(face = "bold", size = 14)) +
       geom_vline(xintercept = 0.05, linetype='dashed', size = 1)+
       annotate("text", y = 1000, x = 0.01, label = "p < 0.05", hjust = -0.5, size = 5, fontface = "bold")

# PART D: PLOTTING OVERLAP WITH METABRIC RESULTS
#' Creates a Venn diagram with the overlapping significant hits when model is 
#' run on the same GOI using two independent data sources
#' @param master_df_list a list of the master DFs to compare, with names being
#' labels for the data source
#' @param goi the name of the gene-of-interest, for labeling purposes
#' @param qval_thres a q-value threshold for significance
plot_hit_overlap_between_data_sources <- function(master_df_list, goi, qval_thres) {
  
  target_sets <- lapply(master_df_list, function(df) unique(df$T_k.name))
  target_overlap <- intersect(target_sets[[1]], target_sets[[2]])
  print(length(target_overlap))
  
  # Limit just to genes found in both sets, then get the significant hits for each
  sig_hits_list <- lapply(master_df_list, function(df) {
    df <- df[df$R_i.name == goi,]
    df <- df[df$T_k.name %fin% target_overlap, ]
    df[df$q.value < qval_thres, 'T_k.name']
  })
  names(sig_hits_list) <- names(master_df_list)
  print(head(sig_hits_list))
  
  # Calculate whether overlap is significant using an approx. of an exact hypergeometric test
  len_hits1 <- length(sig_hits_list[[1]])
  len_hits2 <- length(sig_hits_list[[2]])
  len_ol <- length(intersect(sig_hits_list[[1]], sig_hits_list[[2]]))
  n <- length(target_overlap)
  
  # n is size of urn; hits1 is the number of genes from set 1. Say we were to draw
  # the number of balls the size of hits2 (# of genes from set 2), w/o replacement. 
  # M of these balls are also found in hits2. What are the odds that M >= len_ol?
  phyper_res <- phyper(len_ol-1, len_hits1, n-len_hits1, len_hits2, 
                       lower.tail = FALSE, log.p = FALSE)
  print(paste("p-value:", phyper_res))
  
  # Plot Venn Diagram
  #group.colors <- c("#0072B5FF", "#FFDC91FF")
  #names(group.colors) <- names(sig_hits_list)
  
  #ggVennDiagram(sig_hits_list, edge_size = 5, edge_lty = "solid", category.names = names(sig_hits_list),
   #                  set_size = 8, label_size = 8, label_alpha = 0,
    #                 label = "count", set_color = c("#E18727FF", "#0072B5FF")) + #rep("black", times = length(sig_hits_list))) +
    #theme(legend.position = "none") + #scale_fill_manual(values = group.colors) +
    #scale_fill_gradient(low="white", high = "white") + #labs(fill = paste("# Hits, q<", qval_thres)) +
    #scale_fill_manual(values = c("#FFDC91FF", "#0072B5FF")) + # "#E18727FF"
    #scale_color_manual(values = rep("black", times = length(sig_hits_list))) +
    #scale_color_manual(values = c("#E18727FF", "#0072B5FF")) +
    #labs(title = paste("q <", paste0(qval_thres, paste(", HG p =", paste0(format(phyper_res, digits = 4, scientific = TRUE), "\n")))))
  #text(labels = paste("HG p=", format(phyper_res, digits = 4, scientific = TRUE)),
       #x = 0.1, y = 0.1)
  
  VennDiag <- euler(sig_hits_list)
  plot(VennDiag, counts = TRUE, font=2, cex=6, cex.main=1, alpha=0.9, quantities=T,
       fill=c("#20854EFF", "#FFDC91FF"), col="darkgray", lwd=3,
       main = paste("q <", paste0(qval_thres, paste0(", HG p=", format(phyper_res, digits = 4, scientific = TRUE)))))
}

# Import relevant results
allgenes_p53_metabric <- read.csv("C:/Users/sarae/Documents/res_top_0.05_allGenes_df_rawCNA_methMRaw_Nonsyn.Driver.Vogel.sklearn.elim.vif.5_corrected_MUT.csv",
                                  header = TRUE, check.names = FALSE)

master_df_list <- list("TCGA-BRCA" = allgenes_p53, "METABRIC" = allgenes_p53_metabric)

# Call function
plot_hit_overlap_between_data_sources(master_df_list, "TP53", 0.1)

# Hack to look at the top 5% of hits in both
allgenes_p53_top5perc <- allgenes_p53[1:as.integer(nrow(allgenes_p53) * 0.05),]
allgenes_metabric_top5perc <- allgenes_p53_metabric[1:as.integer(nrow(allgenes_p53_metabric) * 0.05),]
master_df_list_top5perc <- list("TCGA" = allgenes_p53_top5perc, "METABRIC" = allgenes_p53_metabric_top5perc)
plot_hit_overlap_between_data_sources(master_df_list_top5perc, "TP53", 1)

# Reuse this function to look at overlap between male & female
plot_hit_overlap_between_data_sources(list("Male" = allgenes_pc_male, 
                                           "Female" = allgenes_pc_female), 
                                      "TP53", 0.1)


# ASIDE: PLOTTING HIT OVERLAP FOR EACH DRIVER MUTATED AT >= 5% IN AT LEAST 3 
# CANCER TYPES BETWEEN THOSE CANCER TYPES 
#' Take the spearman correlation of the coefficients fit in each cancer types 
#' for hits at or above some threshold, or a HG test of the hits that are 
#' shared at that threshold. Plot barplot per driver
#' @param top_drivers_0.05 list of output data frames from Dyscover per-cancer
#' @param qval_thres a q-value threshold for significance
#' @param method either 'Spearman', 'DP', 'CS', or 'HG' to indicate what kind of
#'  test we are performing
#' @param drivers a vector of the driver gene names that are highly mutated in
#' at least 3 cancer types
create_perdriver_percancer_barplot <- function(top_drivers_0.05, qval_thres, method,
                                           drivers) {

  if((method == "Spearman") | (method == "DP") | (method == "CS") ) {
    vals_per_driver <- lapply(drivers, function(d) {
      top_drivers_0.05_sub <- subset_output_list(top_drivers_0.05, d, qval_thres)
      
      # Merge estimates by gene
      if(length(top_drivers_0.05_sub) > 1) {
        top_drivers_0.05_merged <- Reduce(function(...) 
          merge(..., by = 'T_k.name', all = T), top_drivers_0.05_sub)
        #top_drivers_0.05_merged <- na.omit(top_drivers_0.05_merged)
        #print(head(top_drivers_0.05_merged))
        
        # Keep only rows that have at least two non-NA values
        #top_drivers_0.05_merged <- top_drivers_0.05_merged[rowSums(!is.na(
        #  top_drivers_0.05_merged)) > 1,]

        if(nrow(top_drivers_0.05_merged) > 9) {
          coeff_vects <- lapply(2:ncol(top_drivers_0.05_merged), function(i) {
            vals <- unlist(top_drivers_0.05_merged[,i,with=F])
            #vals[is.na(vals)] <- 0
            return(as.numeric(vals))
          })
          # Take the pairwise spearman/ DP of each vector
          vals <- c()
          for (i in 1:(length(coeff_vects)-1)) {
            for (j in (i+1):length(coeff_vects)) {
              v1 <- coeff_vects[[i]]
              v2 <- coeff_vects[[j]]
              if(method == "Spearman") {
                spearman <- cor.test(v1, v2, method = "spearman")
                spearman_stat <- as.numeric(spearman$estimate)
                vals <- c(vals, spearman_stat)
              } else if (method == "DP") {
                dp <- as.numeric(sum(v1*v2, na.rm = T))
                vals <- c(vals, dp)
              } else {
                vals_to_rm <- unique(c(which(is.na(v1)), which(is.na(v2))))
                if(length(vals_to_rm) > 0) {
                  v1 <- v1[-vals_to_rm]
                  v2 <- v2[-vals_to_rm]
                }
                cos <- cosine(v1, v2)
                vals <- c(vals, cos)
              }
            }
          }
          
          # Alternatively... 
          if(method == "Spearman") {
            corr_mat <- as.data.frame(rcorr(x = as.matrix(
              top_drivers_0.05_merged[,2:ncol(top_drivers_0.05_merged)]), 
              type = "spearman")$r)
            print(head(corr_mat))
            vals <- unlist(lapply(1:(ncol(corr_mat)-1), function(i) {
              return(corr_mat[(i+1):nrow(corr_mat),i])
            }))
          }
          
          return(vals)
        } else {return(NA)}
      }
      else{return(NA)}
    })
    names(vals_per_driver) <- drivers
    vals_per_driver <- vals_per_driver[!is.na(vals_per_driver)]
    print(head(vals_per_driver))
    
    input_df <- melt(vals_per_driver)
    if(method == "Spearman") {
      colnames(input_df) <- c("Spearman", "Driver")
      print(head(input_df))
      g <- ggplot(input_df, aes(x = Driver, y = Spearman, fill = Driver)) + 
        geom_boxplot() + theme_minimal() +
        scale_fill_manual(values = c("TP53" = "#0072B5FF", "PIK3CA" = "#BC3C29FF", 
                                     "KRAS" = "#20854EFF", "CTNNB1" = "#E18727FF", 
                                     "FBXW7" = "palevioletred3")) +
        geom_hline(yintercept = 0, linetype = "dashed") + 
        ylab("Pairwise Spearman Correlation") + 
        scale_x_discrete(limits = rev) + 
      theme(axis.text = element_text(face = "bold", size = 12),
            axis.title = element_text(face = "bold", size = 14),
            axis.ticks = element_blank(), axis.line = element_blank(),
            panel.border = element_blank(), legend.position = "none",
            panel.grid.major = element_line(color='#eeeeee'))
      print(g)
    }
    if(method == "DP") {
      colnames(input_df) <- c("Dot.Product", "Driver")
      print(head(input_df))
      g <- ggplot(input_df, aes(x = Driver, y = Dot.Product, fill = Driver)) + 
        geom_boxplot() + theme_minimal() +
        scale_fill_manual(values = c("TP53" = "#0072B5FF", "PIK3CA" = "#BC3C29FF", 
                                     "KRAS" = "#20854EFF", "CTNNB1" = "#E18727FF", 
                                     "FBXW7" = "palevioletred3")) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        ylab("Pairwise Dot Product") + 
        scale_x_discrete(limits = rev) + 
      theme(axis.text = element_text(face = "bold", size = 12),
            axis.title = element_text(face = "bold", size = 14),
            axis.ticks = element_blank(), axis.line = element_blank(),
            panel.border = element_blank(), legend.position = "none",
            panel.grid.major = element_line(color='#eeeeee'))
      print(g)
    }
    if (method == "CS") {
      colnames(input_df) <- c("Cosine.Similarity", "Driver")
      print(head(input_df))
      g <- ggplot(input_df, aes(x = Driver, y = Cosine.Similarity, fill = Driver)) + 
        geom_boxplot() + theme_minimal() +
        scale_fill_manual(values = c("TP53" = "#0072B5FF", "PIK3CA" = "#BC3C29FF", 
                                     "KRAS" = "#20854EFF", "CTNNB1" = "#E18727FF", 
                                     "FBXW7" = "palevioletred3")) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        ylab("Pairwise Cosine Similarity") + 
        scale_x_discrete(limits = rev) + 
        theme(axis.text = element_text(face = "bold", size = 12),
              axis.title = element_text(face = "bold", size = 14),
              axis.ticks = element_blank(), axis.line = element_blank(),
              panel.border = element_blank(), legend.position = "none",
              panel.grid.major = element_line(color='#eeeeee'))
      print(g)
    }

  } else if (method == "HG") {
    hg_per_driver <- lapply(drivers, function(d) {
      targets_by_ct <- lapply(1:length(top_drivers_0.05), function(i) {
        x <- top_drivers_0.05[[i]]
        if(!is.na(qval_thres)) {
          if(nrow(x[(x$q.value < qval_thres) & (x$R_i.name == d),]) > 0) {
            return(as.character(unlist(x[(x$q.value < qval_thres) & (x$R_i.name == d), 
                     'T_k.name', with=F])))
          } else {return(NA)}
        } else {return(as.character(unlist(x[x$R_i.name == d, 'T_k.name', with=F])))}
      })
      names(targets_by_ct) <- names(top_drivers_0.05)
      targets_by_ct <- targets_by_ct[!is.na(targets_by_ct)]
      targets_by_ct <- targets_by_ct[lengths(targets_by_ct) != 0]
      print(targets_by_ct)
      
      if(length(targets_by_ct) > 1) {
          # Take the pairwise HG of each set of cancer types
          hg <- c()
          for (i in 1:(length(targets_by_ct)-1)) {
            len_hits1 <- length(targets_by_ct[[i]])
            ct1 <- names(targets_by_ct)[i]
            dfi <- top_drivers_0.05[[ct1]]
            
            for (j in (i+1):length(targets_by_ct)) {
              # Calculate whether overlap is significant using an approx. of an 
              # exact hypergeometric test

              len_hits2 <- length(targets_by_ct[[j]])
              target_overlap <- intersect(targets_by_ct[[i]], 
                                          targets_by_ct[[j]])
              len_ol <- length(target_overlap)
              print(len_ol)
              
              ct2 <- names(targets_by_ct)[j]
              dfj <- top_drivers_0.05[[ct2]]
              n <- length(intersect(as.character(unlist(dfi[dfi$R_i.name == d, 'T_k.name'])), 
                                    as.character(unlist(dfj[dfj$R_i.name == d, 'T_k.name']))))
              print(n)
              
              # n is size of urn; hits1 is the number of genes from set 1. Say we were to draw
              # the number of balls the size of hits2 (# of genes from set 2), w/o replacement. 
              # M of these balls are also found in hits2. What are the odds that M >= len_ol?
              phyper_res <- phyper(len_ol-1, len_hits1, n-len_hits1, len_hits2, 
                                   lower.tail = FALSE, log.p = FALSE)
              print(paste("p-value:", phyper_res))
              hg <- c(hg, phyper_res)
            }
          }
          return(hg)
        } else {return(NA)}
    })
    names(hg_per_driver) <- drivers
    hg_per_driver <- hg_per_driver[!is.na(hg_per_driver)]
    print(head(hg_per_driver))
    
    input_df <- melt(hg_per_driver)
    colnames(input_df) <- c("HG.pvalue", "Driver")
    print(head(input_df))
    ggplot(input_df, aes(x = Driver, y = HG.pvalue)) + geom_boxplot()
  } else {
    print("Only implemented for methods 'Spearman', 'DP', and 'HG'. Please try 
    again with one of these inputs.")
  }
}
    

#' Helper function to subset each output data frame in the list by driver d and
#' by q-value; keep only the necessary columns (estimate and target name)
#' @param top_drivers_0.05 list of output data frames from Dyscover per-cancer
#' @param d the name of the given driver
#' @param qval_thres a q-value threshold for significance
subset_output_list <- function(top_drivers_0.05, d, qval_thres) {
  top_drivers_0.05_sub <- lapply(1:length(top_drivers_0.05), function(i) {
    x <- top_drivers_0.05[[i]]
    ct <- names(top_drivers_0.05)[i]
    
    if(d %fin% x$R_i.name) {
      # Keep only if there are at least 15 hits for this driver
      if(nrow(x[(x$R_i.name == d) & (x$q.value < 0.2),]) < 15) {return(NA)}
      #est_newname <- paste0('estimate.', ct)
      #colnames(x)[which(colnames(x) == 'estimate')] <- est_newname
      est_newname <- paste0('statistic.', ct)
      colnames(x)[which(colnames(x) == 'statistic')] <- est_newname
      cols_to_keep <- c(est_newname, 'T_k.name')
      
      if(!is.na(qval_thres)) {
        if(nrow(x[(x$q.value < qval_thres) & (x$R_i.name == d),]) > 0) {
          return(x[(x$q.value < qval_thres) & (x$R_i.name == d), 
                   which(colnames(x) %fin% cols_to_keep), with=F])
        } else {return(NA)}
      } else {return(x[x$R_i.name == d, which(colnames(x) %fin% cols_to_keep), 
                       with=F])}
    } else {return(NA)}
  })
  top_drivers_0.05_sub <- top_drivers_0.05_sub[!is.na(top_drivers_0.05_sub)]
  #print(head(top_drivers_0.05_sub))
  
  return(top_drivers_0.05_sub)
}

# Get the drivers that are highly mutated in at least 3 cancer types
drivers_all <- unlist(lapply(top_drivers_0.05, function(x) unique(x$R_i.name)))
drivers_all <- drivers_all[!is.na(drivers_all)]
drivers_tab <- table(drivers_all)
drivers_freq_atleast3 <- names(drivers_tab)[drivers_tab >= 3]

# Call function
create_perdriver_percancer_barplot(top_drivers_0.05, 0.2, 'Spearman', 
                               drivers_freq_atleast3)
create_perdriver_percancer_barplot(top_drivers_0.05, NA, 'HG', drivers_freq_atleast3)


#' Related function that does this for one driver gene across all cancer types,
#' in order to determine which cancer types are the most similar
#' @param top_drivers_0.05 list of output DFs from Dyscovr for each cancer type
#' @param qval_thres q-value threshold for significance
#' @param driver the name of a driver gene of interest
create_similarity_hm_percancer <- function(top_drivers_0.05, qval_thres, driver) {

  # Create a list for each cancer type that, for each target, extracts its beta
  up_and_down_list <- lapply(1:length(top_drivers_0.05), function(i) {
    
    x <- top_drivers_0.05[[i]]
    ct <- names(top_drivers_0.05)[i]
    
    # Subset list of outputs to the driver gene of interest, as well as the appr. 
    # q-value threshold if desired
    x <- x[x$R_i.name == driver,]
    if(!is.na(qval_thres)) {
      x <- x[x$q.value < qval_thres,]
    }
    if(nrow(x) < 10) {return(NA)}  # remove cancer types with <10 hits for this target
    
    # Determine what is up/down
    #dirs <- unlist(lapply(x$estimate, function(e) ifelse(e > 0, 1, -1)))
    dirs <- x$estimate
    df <- data.frame("direction" = dirs, "targets" = x$T_k.name)
    colnames(df)[which(colnames(df) == "direction")] <- ct
    return(df)
  })
  names(up_and_down_list) <- names(top_drivers_0.05)
  up_and_down_list <- up_and_down_list[!is.na(up_and_down_list)]
  
  up_and_down_df <- Reduce(function(...) 
    merge(..., by = 'targets', all=T), up_and_down_list)

  rownames(up_and_down_df) <- up_and_down_df$targets
  up_and_down_df <- up_and_down_df[,2:ncol(up_and_down_df)]
  print(head(as.matrix(up_and_down_df)))
  
  # Note that Hmisc removes NAs according to the pairwise comparisons
  sim_mat <- rcorr(x = as.matrix(up_and_down_df), type = "spearman")
  print(head(sim_mat))
  
  df_long <- inner_join(
    melt(sim_mat$r, value.name = "r"),
    melt(sim_mat$P, value.name = "p"),
    by = c("Var1", "Var2")
  ) %>%
    dplyr::mutate(r_lab = abbreviateSTR(r, "r"), p_lab = abbreviateSTR(p, "p")) %>% 
    dplyr::mutate(label = paste(r_lab, p_lab, sep = "\n")) %>%
    rowwise() %>%
    dplyr::mutate(pair = sort(c(Var1, Var2)) %>% paste(collapse = ",")) %>%
    group_by(pair) %>%
    distinct(pair, .keep_all = T)
  df_long <- as.data.frame(df_long)
  print(head(df_long))
  df_long[which(df_long$Var1 == df_long$Var2), 'r'] <- NA 
  df_long$label <- unlist(lapply(df_long$label, function(x) 
    ifelse(x == "r = 1\n", NA, x)))
  #df_long <- df_long[-which(df_long$Var1 == df_long$Var2),] 
  print(head(df_long))

  # Create the visualization
  ggplot(df_long, aes(Var1, Var2, fill = r)) + xlab("Cancer Type") + 
    ylab("Cancer Type") + geom_raster() +
    #geom_tile(aes(fill = Avg.Spearman), color = 'white') +
    geom_text(aes(label = label)) + theme_minimal() +
    scale_fill_gradient(low = '#eeeeee', high = '#7876B1FF', space = 'Lab',
                        na.value = "transparent") +
    theme(axis.text.x = element_text(angle = 90, face = "bold", size = 12),
          axis.text.y = element_text(face = "bold", size = 12),
          axis.title = element_text(face = "bold", size = 14),
          legend.title = element_text(face = "bold", size = 14),
          axis.ticks = element_blank(), axis.line = element_blank(),
          #panel.border = element_blank(), 
          panel.grid.major = element_line(color='#eeeeee'))

}

top_drivers_0.05_tp53 <- lapply(top_drivers_0.05, function(x) {
  if("TP53" %in% x$R_i.name) return(x) else{return(NA)} 
})
top_drivers_0.05_tp53 <- top_drivers_0.05_tp53[!is.na(top_drivers_0.05_tp53)]
create_similarity_hm_percancer(top_drivers_0.05_tp53, 0.2, "TP53")


abbreviateSTR <- function(value, prefix){  # format string more concisely
  lst = c()
  for (item in value) {
    if (is.nan(item) || is.na(item)) { # if item is NaN return empty string
      lst <- c(lst, '')
      next
    }
    item <- round(item, 2) # round to two digits
    item_char <- ""
    if (item == 0) { # if rounding results in 0 clarify
      item_char = '< 0.01'
    } else {item_char <- as.character(item)}
    #item_char <- sub("(^[0])+", "", item_char)    # remove leading 0: 0.05 -> .05
    #item_char <- sub("(^-[0])+", "-", item_char)  # remove leading -0: -0.05 -> -.05
    if(!(item == 0)) {lst <- c(lst, paste(prefix, item_char, sep = " = "))}
    else {lst <- c(lst, paste(prefix, item_char, sep = " "))}
  }
  return(lst)
}


#' Create a similar heatmap using the dot product or cosine similarity
#' @param top_drivers_0.05 list of output DFs from Dyscovr for each cancer type
#' @param driver the name of a driver gene of interest
#' @param method either 'dot' or 'cosine' to indicate dot product or cosine similarity
create_dotproduct_cosinesim_hm_percancer <- function(top_drivers_0.05, driver, method) {
  
  # Create a list for each cancer type that, for each target, determines if it is
  # being up- or down-regulated by the driver
  ct_list <- lapply(1:length(top_drivers_0.05), function(i) {
    
    x <- top_drivers_0.05[[i]]
    ct <- names(top_drivers_0.05)[i]
    
    # Subset list of outputs to the driver gene of interest, as well as the appr. 
    # q-value threshold if desired
    if(driver %fin% x$R_i.name) {
      x <- x[x$R_i.name == driver,]
      
      if(nrow(x[x$q.value < 0.2,]) < 10) {return(NA)}  # remove cancer types with <20 hits for this target
      
      betas <- x$estimate
      df <- data.frame("estimate" = betas, "targets" = x$T_k.name)
      colnames(df)[which(colnames(df) == "estimate")] <- ct
      return(df)
    } else {return(NA)}
  })
  names(ct_list) <- names(top_drivers_0.05)
  ct_list <- ct_list[!is.na(ct_list)]
  
  ct_df <- Reduce(function(...) merge(..., by = 'targets', all=T), ct_list)
  rownames(ct_df) <- ct_df$targets
  ct_df <- ct_df[,2:ncol(ct_df)]
  print(head(ct_df))
  
  pairs_df <- as.data.frame(expand.grid.unique(names(ct_list), names(ct_list)))
  colnames(pairs_df) <- c("Var1", "Var2")
  print(head(pairs_df))
  
  # Calculate pairwise dot product or cosine similarity
    pairs_df$res <- unlist(lapply(1:nrow(pairs_df), function(i) {
      ct1 <- as.character(pairs_df[i, 'Var1'])
      ct2 <- as.character(pairs_df[i, 'Var2'])
      if(ct1 == ct2) {return(NA)}
      else {
        vals_ct1 <- as.numeric(ct_df[, names(ct_df) == ct1])
        vals_ct2 <- as.numeric(ct_df[, names(ct_df) == ct2])
        if(method == "dot") {
          dp <- sum(vals_ct1*vals_ct2, na.rm = T)
          print(dp)
          return(dp)
        } else if(method == "cosine") {
          vals_to_rm <- unique(c(which(is.na(vals_ct1)), which(is.na(vals_ct2))))
          vals_ct1 <- vals_ct1[-vals_to_rm]
          vals_ct2 <- vals_ct2[-vals_to_rm]
          cos <- cosine(vals_ct1, vals_ct2)
          print(cos)
          return(cos)
        } else {
          print("Only implemented for 'dot' or 'cosine' methods. Please try again.")
          return(NA)
        }
      }
    }))

  print(head(pairs_df))
  
  # Create the visualization
  ggplot(pairs_df, aes(Var1, Var2, fill = res)) + xlab("Cancer Type") + 
    ylab("Cancer Type") + geom_raster() +
    #geom_tile(aes(fill = Avg.Spearman), color = 'white') +
    geom_text(aes(label = round(res, 2))) + theme_minimal() +
    scale_fill_gradient2(mid = '#eeeeee', low = "#FFDC91FF", high = '#7876B1FF', 
                         space = 'Lab', na.value = "transparent",
                         name = "DP", position = "bottom") +
    theme(axis.text.x = element_text(angle = 90, face = "bold", size = 12),
          axis.text.y = element_text(face = "bold", size = 12),
          axis.title = element_text(face = "bold", size = 14),
          legend.title = element_text(face = "bold", size = 14),
          axis.ticks = element_blank(), axis.line = element_blank(),
          #panel.border = element_blank(), 
          panel.grid.major = element_line(color='#eeeeee'))
  
}

expand.grid.unique <- function(x, y, include.equals=TRUE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i) {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

create_dotproduct_cosinesim_hm_percancer(top_drivers_0.05, "TP53", "dot")


# PART E: PLOTTING ENRICHMENT IN LI. ET AL. TARGETS

# Analyzes the data from Li et al. (2019), Nature Medicine, from
# Supplementary tables 1-4.
# Link to Paper: https://www.nature.com/articles/s41591-019-0404-8#Sec31

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/CL Metabolic Data (Li et al. 2019)/"

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


# Call helper function to get the correlations
tp53_mut_df <- get_mutation_metabol_corr_for_goi(mut_metabol_corr, "TP53", 2.0, path)
pik3ca_mut_df <- get_mutation_metabol_corr_for_goi(mut_metabol_corr, "PIK3CA", 2.0, path)


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


#' For the top hits of a given protein-of-interest, plot an 
#' enrichment curve for whether these hits are associated with 
#' a significantly dysregulated metabolite from this paper
#' @param master_df a master DF produced from our model, subsetted
#' to the given protein of interest
#' @param mut_df genes from a metabolism DF for a particular protein of interest 
#' @param qval_thres a qvalue threshold for significance
#' @param tval_thres a t-value threshold that was used for inclusion of Li et al. 
#' metabolic genes
plot_target_assoc_w_metabol_enrichment <- function(master_df, mut_df_genes, 
                                                   qval_thres, tval_thres) {
  master_df_sig <- master_df[master_df$q.value < qval_thres,]
  
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
    enrich_vals <- c(enrich_vals, frac)
  }
  
  input_df <- data.frame("Rank" = rank, "Enrich.Vals" = enrich_vals)
  
  p <- ggplot(data = input_df, aes(x = Rank, y = Enrich.Vals)) + geom_line(linewidth = 2, color = "#0072B5FF") + 
    xlab("Rank") + theme_minimal() + 
    theme(axis.title = element_text(face = "bold", size = 14), 
          axis.text = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12)) +
    ylab(paste("Frac. of Hits (q <", paste0(qval_thres, paste(") that are Metabolic Reg (t >", paste0(tval_thres, ")"))))) 
    
  p
}

metabolism_gene_list_ensg <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/metabolism_gene_list_ensg.txt",
                                        header = T)[,1]

tp53_mut_df_genes <- unique(unlist(strsplit(as.character(unlist(tp53_mut_df$Assoc.Genes)), ",", fixed = TRUE)))
tp53_mut_df_genes_overlap <- tp53_mut_df_genes[tp53_mut_df_genes %in% metabolism_gene_list_ensg]
tp53_mut_df_genes_overlap_names <- unique(unlist(lapply(tp53_mut_df_genes_overlap, function(x) 
  all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == x, 'external_gene_name'])))


plot_target_assoc_w_metabol_enrichment(master_df_mut, tp53_mut_df_genes_overlap_names, 0.1, 2.0)


# PART F: SHAP PLOTS FOR BETAS, INSTEAD OF A VOLCANO PLOT; FOR EACH TOP TARGET FOR
# A GIVEN DRIVER, WE CAN SHOW WHAT THE BETAS ARE IN EACH CANCER TYPE

#' Function to create a SHAP-style plot, adapted from Antonio's code
#' @param target_names vector of targets of a given driver we'd like to plot
#' @param top_driver_list a list of output master DFs from model, names correspond to cancer types
#' @param goi the gene-of-interest, or a list
create_shap_plot <- function(target_names, top_driver_list, goi) {
  
  # Subset the top driver list to this goi, and to the given targets
  top_driver_dfs_list <- lapply(1:length(top_driver_list), function(i) {
    x <- top_driver_list[[i]]
    if(length(intersect(x$R_i.name, goi)) > 0) {
      x_full <- data.frame()
      for(g in goi) {
        x_sub <- x[(x$R_i.name == g) & (x$T_k.name %fin% target_names), 
                   c('T_k.name', 'estimate', 'p.value')]
        x_sub$negLogPval <- as.numeric(unlist(lapply(x_sub$p.value, function(x) 
          ifelse(is.numeric(x), -log2(x), NA))))
        x_sub$dir <- as.integer(unlist(lapply(x_sub$estimate, function(x) 
          ifelse(x > 0, 1, -1))))
        x_sub$Dir.negLogPval <- as.numeric(unlist(lapply(1:nrow(x_sub), function(i) {
          dir <- x_sub[i, 'dir'] 
          return(x_sub[i, 'negLogPval'] * dir)
        })))
        x_sub$ct <- as.factor(rep(names(top_driver_list)[i], times = nrow(x_sub)))
        x_sub$driver <- as.factor(rep(g, times = nrow(x_sub)))
        colnames(x_sub) <- c("Gene", "Beta", "P.value", "Neg.Log.P.value", "Dir", 
                             "Dir.Neg.Log.P.value", "Cancer.Type", "Driver")
        
        x_sub$Gene <- as.factor(x_sub$Gene)
        x_sub$Beta <- as.numeric(x_sub$Beta)
        
        x_full <- rbind(x_full, x_sub)
      }
      print(head(x_full))
      return(x_full)
    } else {return(NA)}
  })
  top_driver_dfs_list <- top_driver_dfs_list[!is.na(top_driver_dfs_list)]
  
  top_driver_df <- do.call(rbind, top_driver_dfs_list)
  print(top_driver_df)
  #top_driver_df_m <- melt(top_driver_df)
  #print(top_driver_df_m)
  
  g <- ggplot(top_driver_df, aes(x=Dir.Neg.Log.P.value, y=Gene, color=Dir.Neg.Log.P.value)) + 
  #g <- ggplot(top_driver_df, aes(x=Dir.Neg.Log.P.value, y=Gene, color=Driver)) + 
    geom_beeswarm(method = "compactswarm", corral = "wrap", corral.width=0.9, size = 3.5) + #alpha=0.9, 
    #scale_color_manual(values = list("TP53" = "#0072B5FF", "PIK3CA" = "#BC3C29FF", 
                              #"KRAS" = "#20854EFF", "IDH1" = "#E18727FF")) + #, "#7876B1FF"), "#FFDC91FF"), "#EE4C97FF"),))
    scale_color_gradient2(low="#0072B5FF", mid="darkgrey", high="#BC3C29FF") + 
    theme_minimal() + scale_y_discrete(limits=rev) + 
    labs(x="-log2(p.value) * direction", y="Gene", color="Driver") + 
    theme(axis.text.x=element_text(size=14), axis.title.x=element_text(size=14, face = "bold"), 
          axis.title.y=element_text(size=14, face = "bold"), axis.text.y=element_text(size=14),
          legend.title = element_text(size=14, face="bold"), 
          legend.text = element_text(size=14), legend.position = "bottom") 
          #legend.position = "none")
  
  print(g)
  
  return(top_driver_df)
}

top_metabolic_targets_tp53 <- top_drivers_metabol[top_drivers_metabol$R_i.name == "TP53", 'T_k.name'][1:10]
top_allgene_targets_tp53 <- top_drivers_allgenes[top_drivers_allgenes$R_i.name == "TP53", 'T_k.name'][1:10]

#' Function to get the top targets for a particular driver of interest across different cancer types.
#' Two methods are 'any', which chooses the most significant hits for the given driver that 
#' originate from any cancer type (do not need to be shared) and 'across', which chooses the 
#' significant hits that are found in the most cancer types.
#' @param top_drivers_0.05 list of master DFs, named by cancer type or subtype
#' @param goi the name of the driver gene of interest
#' @param qval_thres the minimum q-value to be considered a significant hit
#' @param any_or_across either 'any' or 'across', to indicate the method used (as described above)
get_top_targets_per_cancer <- function(top_drivers_0.05, goi, qval_thres, any_or_across) {
  top_allgene_targets <- c() 
  
  for(df in top_drivers_0.05) {
    top_allgene_targets <- c(top_allgene_targets, 
                             df[((df$R_i.name == goi) & (df$q.value < qval_thres)), 'T_k.name'])
  }
  
  # Choose the most significant hits from any cancer type
  if(any_or_across == "any") {
    top_allgene_targets <- unique(unlist(top_allgene_targets))  
    
  # Alternatively, choose the targets that are most common among the cancer types
  } else if (any_or_across == "across") {
    top_allgene_targets_t <- melt(table(unlist(top_allgene_targets)))  
    top_allgene_targets_t <- top_allgene_targets_t[order(top_allgene_targets_t$value, decreasing = T),]
    top_allgene_targets <- as.character(top_allgene_targets_t[top_allgene_targets_t$value == max(top_allgene_targets_t$value), "Var1"])
    while(length(top_allgene_targets) < 5) {
      top_allgene_targets_t <- top_allgene_targets_t[top_allgene_targets_t$value != max(top_allgene_targets_t$value),]
      top_allgene_targets <- c(top_allgene_targets, as.character(top_allgene_targets_t[top_allgene_targets_t$value == 
                                                                                         max(top_allgene_targets_t$value), "Var1"]))
    }
    
  } else {
    print("Error, only options are 'any' or 'across'.")
  }
  
  return(top_allgene_targets)
}

# Call helper 
top_allgene_targets_ctnnb1 <- get_top_targets_per_cancer(top_drivers_0.05, "CTNNB1", 0.13, "across")
top_allgene_targets_pik3ca <- get_top_targets_per_cancer(top_drivers_0.05, "PIK3CA", 0.06, "across")
top_allgene_targets_tp53 <- get_top_targets_per_cancer(top_drivers_0.05, "TP53", 0.05, "across")


shap_res <- create_shap_plot(top_metabolic_targets_tp53, top_drivers_0.05_metabol, "TP53")
shap_res <- create_shap_plot(top_allgene_targets_ctnnb1, top_drivers_0.05, "CTNNB1")

pc_allGenes_sig <- pc_allGenes[pc_allGenes$q.value < 0.08,]
sig_genes_alldrivers <- names(which(table(pc_allGenes_sig$T_k.name) == 4))[1:10]

shap_res <- create_shap_plot(intersecting_genes, pc_allGenes, c("TP53", "IDH1", "PIK3CA", "KRAS"))

# Limit to cancer types with the given driver
tp53_cancer_types <- c("ACC", "BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KICH", "LGG", 
                       "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "PRAD", "SARC", "STAD", "UCEC")
create_shap_plot(top_metabolic_targets_tp53, top_drivers_0.05_metabol[names(top_drivers_0.05_metabol) %in% tp53_cancer_types], "TP53")


top_metabolic_targets_pik3ca <- pc_metabol[pc_metabol$R_i.name == "PIK3CA", 'T_k.name'][1:10]
pik3ca_cancer_types <- c("BLCA", "BRCA", "COAD", "CESC", "UCEC")
create_shap_plot(top_metabolic_targets_pik3ca, top_drivers_0.05_metabol[names(top_drivers_0.05_metabol) %in% pik3ca_cancer_types], "PIK3CA")

kras_cancer_types <- c("LUAD", "COAD", "PAAD")

# Do this for a specific GOI and specific target set (e.g. amino acid metabolism, one-carbon metabolism)
top_aa_metabolism_targets_tp53 <- top_drivers_metabol[(top_drivers_metabol$R_i.name == "TP53") & 
                                                        (top_drivers_metabol$T_k %fin% aa_metabolism_genes), 
                                                      'T_k.name'][1:10]
create_shap_plot(top_aa_metabolism_targets_tp53, top_drivers_0.05_metabol[names(top_drivers_0.05_metabol) %in% tp53_cancer_types], "TP53")

top_1c_metabolism_targets_tp53 <- top_drivers_metabol[(top_drivers_metabol$R_i.name == "TP53") & 
                                                        (top_drivers_metabol$T_k %fin% one_carbon_metabolism_genes), 
                                                      'T_k.name'][1:10]
create_shap_plot(top_1c_metabolism_targets_tp53, top_drivers_0.05_metabol[names(top_drivers_0.05_metabol) %in% tp53_cancer_types], "TP53")


########################################################################
### SURVIVAL CURVES AND DRUG RESPONSE
########################################################################
#' Create a survival curve for a given driver gene using TCGA clinical data
#' @param clinical_df a clinical DF from the TCGA for a certain patient set
#' @param driver_name the gene name of the driver of interest
#' @param mutation_count_matrix a mutation count matrix from maftools for TCGA data
create_survival_curves <- function(clinical_df, driver_name, mutation_count_matrix) {
  clinical_df_survival <- clinical_df[,colnames(clinical_df) %fin% c("case_submitter_id", "vital_status","days_to_death","days_to_last_follow_up")]
  clinical_df_survival <- distinct(clinical_df_survival[which(!is.na(clinical_df_survival$vital_status)),])
  #table(clinical_df_survival$vital_status)   #11532 alive, 4182 dead, 14 not reported
  clinical_df_survival$status <- ifelse(clinical_df_survival$vital_status == 'Alive', 0,1)
  clinical_df_survival$patient <- unlist(lapply(clinical_df_survival$case_submitter_id, function(x)
    unlist(strsplit(x, "-", fixed = T))[3]))
  
  clinical_df_survival$days <- ifelse(clinical_df_survival$status == 0,
                                      clinical_df_survival$days_to_last_follow_up,
                                      clinical_df_survival$days_to_death)
  
  #print(clinical_df_survival)
  
  # Merge the driver mutation status information with the clinical information, if needed
  if(!("patient" %in% colnames(mutation_count_matrix))) {
    mutation_count_matrix_driver <- melt(mutation_count_matrix[mutation_count_matrix$Gene_Symbol == driver_name,])[,c("variable", "value")]
    colnames(mutation_count_matrix_driver) <- c('patient', 'mutation_status')
    #print(head(mutation_count_matrix_driver))
    mutation_count_matrix_driver$patient <- unlist(lapply(mutation_count_matrix_driver$patient, function(x)
      unlist(strsplit(as.character(x), "-", fixed = T))[1]))
    mutation_count_matrix_driver$mutation_status <- unlist(lapply(mutation_count_matrix_driver$mutation_status, function(x)
      ifelse(x > 0, 1, 0)))
  } else {
    mutation_count_matrix_driver <- mutation_count_matrix
  }
  
  
  survival_df <- merge(mutation_count_matrix_driver, clinical_df_survival[,c("patient","days","status")],
                       by = "patient")
  survival_df <- survival_df[survival_df$days != "'--",]
  survival_df$days <- as.numeric(survival_df$days)
  survival_df$treated <- unlist(lapply(survival_df$response, function(x) 
    ifelse(is.na(x), 0, 1))) 
  print(head(survival_df))
  
  # Compute survival
  #s_mut <- Surv(time = survival_df[survival_df$mutation_status == 1, "days"], 
  #              event = survival_df[survival_df$mutation_status == 1, "status"])
  #s_nonmut <- Surv(time = survival_df[survival_df$mutation_status == 0, "days"], 
  #              event = survival_df[survival_df$mutation_status == 0, "status"])
  
  # Use the survminer package to create the curves
  fit <- survfit(Surv(time = days, event = status) ~ mutation_status + treated, 
                 data = survival_df)
  p <- ggsurvplot(fit, data = survival_df, palette = c("#0072B5FF", "#0072B599", "#BC3C29FF", "#BC3C2999"),
                  pval = T, ggtheme = theme_minimal(), 
                  #legend.title = paste(driver_name, "Mutation Status"), legend = "bottom",
                  #legend.labs = c("No Mutation", "Mutation"), font.main = c(16, "bold", "black"),
                  legend.title = paste(driver_name, "Mutation Status + Treatment Status"), legend = "bottom",
                  legend.labs = c("No Mutation, No Treatment", "No Mutation, Treatment", 
                                  "Mutation, No Treatment", "Mutation, Treatment"), font.main = c(16, "bold", "black"),
                  font.x = c(14, "bold", "black"), font.y = c(14, "bold", "black"), 
                  font.tickslab = c(12, "bold", "black"), font.legend = c(12, "bold", "black")) + 
    xlab("Time (Days)") + ylab("Survival Probability")
  #ggpar(p, font.x = c(12, "bold", "black"), font.y = c(12, "bold", "black"))
  print(p)
  
  # Print s
  #print(s)
  return(survival_df)
}

# Import pan-cancer clinical DF
clinical_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/clinical_data_subset.csv", 
                        header = TRUE, check.names = FALSE)

# Import the pan-cancer mutation count matrix
mutation_count_matrix <- read.csv(paste(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv", sep = ""), 
                                  header = TRUE, check.names = FALSE)

# Call function
tp53_survival_res <- create_survival_curves(clinical_df, "TP53", mutation_count_matrix)


### DRUG RESPONSE ###
# Data source: http://lifeome.net/supp/drug_response/
drug_response <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Drug_Response/drug_response.txt", 
                            header = T, sep = "\t")

# How many patients in each cancer type were treated with drugs?
patients <- distinct(drug_response[,c(1,2)])
cancer_freq_table <- melt(table(patients$cancers))
colnames(cancer_freq_table) <- c("Cancer.Type", "Num.Unique.Patients")
cancer_freq_table <- cancer_freq_table[order(cancer_freq_table$Num.Unique.Patients, decreasing = T),]
cancer_freq_table$Num.Unique.Patients <- as.numeric(cancer_freq_table$Num.Unique.Patients)
ggplot(cancer_freq_table, aes(x = reorder(Cancer.Type, -Num.Unique.Patients), 
                              y = Num.Unique.Patients)) + 
  geom_bar(position="stack", stat="identity", fill = "#0072B5FF") + theme_minimal() +
  xlab("TCGA Cancer Type") + ylab("Number of Unique Patients") +
  theme(axis.text = element_text(face="bold", size = 16), axis.text.x = element_text(angle = 45, vjust =1, hjust=1), 
        axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14, face="bold"),
        legend.text = element_text(size=12))

# In how many cancer types was each drug used?
drugs <- distinct(drug_response[,c(1,3)])
drug_freq_table <- melt(table(drugs$drug.name))
colnames(drug_freq_table) <- c("Drug.Name", "Num.Cancer.Types")
drug_freq_table <- drug_freq_table[order(drug_freq_table$Num.Cancer.Types, decreasing = T),]
drug_freq_table$Num.Cancer.Types <- as.numeric(drug_freq_table$Num.Cancer.Types)
drug_freq_table$Drug.Name <- as.character(gsub("\\\\", ".", gsub(">", ".", gsub("<", ".", drug_freq_table$Drug.Name))))
ggplot(drug_freq_table, aes(x = reorder(Drug.Name, -Num.Cancer.Types), 
                              y = Num.Cancer.Types)) + 
  geom_bar(position="stack", stat="identity", fill = "#BC3C29FF") + theme_minimal() +
  xlab("Drug Name") + ylab("Number of Cancer Types Treated") +
  theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 45, vjust =1, hjust=1), 
        axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14, face="bold"),
        legend.text = element_text(size=12))

# Histogram of the frequency values
ggplot(drug_freq_table, aes(x = Num.Cancer.Types)) + geom_histogram(position = "identity", bins = 5, fill = "#BC3C29FF") + 
  theme_minimal() + xlab("Number of Cancer Types Treated") + ylab("Frequency") +
  theme(axis.text = element_text(size = 16, face="bold"), axis.text.x = element_text(angle = 45, vjust =1, hjust=1), 
        axis.title=element_text(size=18,face="bold"))

# Get a table of the top drugs with their cancer types
drugs_w_cancer_types <- data.table("drug.name" = unique(drugs$drug.name),
                                   "cancer.types" = unlist(lapply(unique(drugs$drug.name), function(drug)
                                     paste(unique(drugs[drugs$drug.name == drug, 'cancers']), collapse = ", "))))

# Call function
tp53_survival_res <- create_survival_curves(clinical_df, "TP53", mutation_count_matrix, drug_response)

# Plot drug response by mutation status 
mutation_count_matrix_driver <- melt(mutation_count_matrix[mutation_count_matrix$Gene_Symbol == driver_name,])[,c("variable", "value")]
colnames(mutation_count_matrix_driver) <- c('patient', 'mutation_status')
#print(head(mutation_count_matrix_driver))
mutation_count_matrix_driver$patient <- unlist(lapply(mutation_count_matrix_driver$patient, function(x)
  unlist(strsplit(as.character(x), "-", fixed = T))[1]))
mutation_count_matrix_driver$mutation_status <- unlist(lapply(mutation_count_matrix_driver$mutation_status, function(x)
  ifelse(x > 0, 1, 0)))

# Merge with drug treatment and response information
drug_response_cancer <- drug_response[drug_response$cancers %in% c("LGG", "LUAD"),]
drug_response_drug <- drug_response[drug_response$drug.name == "Pemetrexed",]

drug_response_drug <- drug_response[drug_response$drug.name %fin% c("Vinblastine", "Vincristine", "Vinorelbine"),]

mutation_count_matrix_driver$vin.response <- unlist(lapply(mutation_count_matrix_driver$patient, function(x) {
  print(x)
  val <- unique(drug_response_drug[grepl(x, drug_response_drug$patient.arr), 'response'])
  if(length(val) == 0) {return(NA)}
  if(length(val) > 1) {return(val[1])}
  else{return(val)}
}))

#' Create a bar plot for patient drug response by mutation status and conduct a 
#' Chi-squared goodness of fit test and Fisher's exact test
#' @param mutation_count_matrix_driver driver mutation status for each sample
#' @param col_name the name of the mutation count matrix that has the drug response information
#' @param lab a label with the name of the drug
#' @param goi the name of the driver gene-of-interest
create_drug_response_barplot <- function(mutation_count_matrix_driver, col_name, lab, goi) {
  input_df_nomut <- melt(table(mutation_count_matrix_driver[mutation_count_matrix_driver$mutation_status == 0,
                                                            col_name]))
  colnames(input_df_nomut) <- c("Response", "Frequency")
  input_df_nomut$Mutation.Status <- rep(0, times = nrow(input_df_nomut))
  input_df_mut <- melt(table(mutation_count_matrix_driver[mutation_count_matrix_driver$mutation_status == 1,
                                                          col_name]))
  colnames(input_df_mut) <- c("Response", "Frequency")
  input_df_mut$Mutation.Status <- rep(1, times = nrow(input_df_mut))
  input_df <- rbind(input_df_nomut, input_df_mut)
  input_df$Response <- factor(input_df$Response, levels = c("Complete Response", "Partial Response",
                                                            "Stable Disease", "Clinical Progressive Disease"))
  input_df$Response.Aggr <- unlist(lapply(input_df$Response, function(x) ifelse(grepl("Response", x), "Response", "Disease")))
  input_df$Response.Aggr <- factor(input_df$Response.Aggr, levels = c("Response", "Disease"))
  input_df$Response.Aggr <- as.character(input_df$Response.Aggr)
  
  g <- ggplot(input_df, aes(x = Response.Aggr, y = Frequency, fill = as.factor(Mutation.Status))) +
    geom_bar(stat = "identity", position = position_dodge()) +
    xlab(paste(lab, "Drug Response")) + ylab("Frequency") + labs(fill = paste(goi, "Mutation Status")) +
    theme_minimal() + scale_fill_nejm() +
    theme(axis.text.x=element_text(size=14), axis.title.x=element_text(size=14, face="bold"), 
          axis.title.y=element_text(size=14, face="bold"), axis.text.y=element_text(size=14),
          legend.title = element_text(size=14, face="bold"), 
          legend.text = element_text(size=14), legend.position = "bottom") 
  print(g)
  
  
  # Chi-squared googness-of-fit test
  
  # Difference in overall response vs. disease?
  print(chisq.test(c(sum(input_df[(input_df$Response.Aggr == "Response") & (input_df$Mutation.Status == 0), 'Frequency']),
               sum(input_df[(input_df$Response.Aggr == "Response") & (input_df$Mutation.Status == 1), 'Frequency'])), 
             p = c(0.5, 0.5))$p.value)
  print(chisq.test(c(sum(input_df[(input_df$Response.Aggr == "Disease") & (input_df$Mutation.Status == 0), 'Frequency']),
               sum(input_df[(input_df$Response.Aggr == "Disease") & (input_df$Mutation.Status == 1), 'Frequency'])), 
             p = c(0.5, 0.5))$p.value)
  print(chisq.test(c(sum(input_df[(input_df$Response.Aggr == "Response") & (input_df$Mutation.Status == 1), 'Frequency']),
               sum(input_df[(input_df$Response.Aggr == "Disease") & (input_df$Mutation.Status == 1), 'Frequency'])), 
             p = c(0.5, 0.5))$p.value)
  
  # Fisher's exact test
  conting_table_response <- data.frame("No.Mut" = sum(input_df[(input_df$Response.Aggr == "Response") & (input_df$Mutation.Status == 0), 'Frequency']),
                                       "Mut" = sum(input_df[(input_df$Response.Aggr == "Response") & (input_df$Mutation.Status == 1), 'Frequency']))
  conting_table_disease <- data.frame("No.Mut" = sum(input_df[(input_df$Response.Aggr == "Disease") & (input_df$Mutation.Status == 0), 'Frequency']),
                                      "Mut" = sum(input_df[(input_df$Response.Aggr == "Disease") & (input_df$Mutation.Status == 1), 'Frequency']))
  conting_table <- rbind(conting_table_response, conting_table_disease)
  rownames(conting_table) <- c("Response", "Disease")
  
  test <- fisher.test(conting_table)
  print("Fisher p-value:")
  print(test$p.value)   #0.0536
  
  return(input_df)
}

vin_input_df <- create_drug_response_barplot(mutation_count_matrix_tp53, "vin.response", "\nTUBB-Targeting", "TP53")
vin_input_df <- create_drug_response_barplot(mutation_count_matrix_tp53, "vin.response", "\nTUBB-Targeting", "TP53")
etop_input_df <- create_drug_response_barplot(mutation_count_matrix_tp53, "etop.response", "\nEtoposide", "TP53")
irino_input_df <- create_drug_response_barplot(mutation_count_matrix_idh1, "irino.response", "\nIrinotexan", "IDH1")

##############################################################################
### VENN DIAGRAMS TP53 & PIK3CA MUTATIONS
##############################################################################
#' Creates a Venn Diagram to show how many patients have missense mutations in the
#' given patient set in the given genes 1 and 2
#' @param patient_set a vector of TCGA IDs of patients of interest (XXXX)
#' @param mut_count_matrix a mutation count matrix (patients are columns, rows are genes)
#' @param genes the external names of driver genes
plot_mutational_overlap <- function(patient_set, mut_count_matrix, genes) {
  # Limit the mutation count matrix to the given patient set, if given
  if(length(patient_set) > 2) {
    colnames(mut_count_matrix)[2:ncol(mut_count_matrix)] <- unlist(lapply(colnames(mut_count_matrix)[2:ncol(mut_count_matrix)], function(x) {
      return(unlist(strsplit(x, "-", fixed = TRUE))[1])
    }))
    mut_count_matrix <- mut_count_matrix[,c(1, (which(colnames(mut_count_matrix)[2:ncol(mut_count_matrix)] %fin% patient_set) + 1))]
  }

  # Get the number of tumor samples with a mutation in each gene (count >= 1)
  #mut_samps <- list("all" = 1:ncol(mut_count_matrix))
  mut_samps <- list()
  for(g in genes) {
    mut_gene_counts <- as.numeric(unlist(mut_count_matrix[mut_count_matrix$Gene_Symbol == g, 
                                                              2:ncol(mut_count_matrix)]))
    mut_gene_vect <- list(which(mut_gene_counts >= 1))
    mut_samps <- append(mut_samps, mut_gene_vect)
  }
  print(head(mut_samps))
  names(mut_samps)[2:length(mut_samps)] <- genes
  
  # Make the Venn Diagram
  ggVennDiagram(mut_samps, label_alpha = 0, category.names = genes, set_color = "black",
                set_size = 9, label_size = 9, edge_size = 0, label = "count") + 
    theme(legend.position = "none") +
    ggplot2::scale_fill_gradient(low = "#FFDC91FF", high = "#6F99ADFF")
}

patient_set <- read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids_washu.txt"), header = TRUE)[,1]
patient_set_lt20pamp <- intersect(read.table(paste0(main_path, "Patient Subsets/LessThan20PercAmp_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumA <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.A_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumB <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.B_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumAB <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.A.B_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_basal <- intersect(read.table(paste0(main_path, "Patient Subsets/Basal_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_her2 <- intersect(read.table(paste0(main_path, "Patient Subsets/HER2_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_normLike <- intersect(read.table(paste0(main_path, "Patient Subsets/Normal-like_patient_ids.txt"), header = TRUE)[,1], patient_set)

patient_set_pc <- read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt"), header = TRUE)[,1]

patient_set_male <- intersect(patient_set_pc, read.table(paste0(main_path, "Patient Subsets/male_patients.txt"), header = TRUE)[,1])
patient_set_female <- intersect(patient_set_pc, read.table(paste0(main_path, "Patient Subsets/female_patients.txt"), header = TRUE)[,1])


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
mutation_count_matrix <- read.csv(paste0(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv"), 
                                  header = TRUE, check.names = FALSE)
#mutation_count_matrix <- read.csv(paste0(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_ALL.csv"), 
                                  #header = TRUE, row.names = 1, check.names = FALSE)

gene1 <- "TP53"
gene2 <- "PIK3CA"

# Call function
plot_mutational_overlap(patient_set, mutation_count_matrix, c(gene1, gene2))


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
                  Rowv = TRUE, key = TRUE, key.title = NA, key.xlab = "T-statistic") 
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
  if(!("R_i.name" %in% colnames(results_table))) {
    mapping <- data.frame("uniprot" = unlist(lapply(unique(results_table$term), function(x) 
      unlist(strsplit(x, "_", fixed = TRUE))[1])))
    mapping$name <- unlist(lapply(mapping$uniprot, function(x) 
      unique(unlist(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'external_gene_name']))))
    results_table$R_i.name <- unlist(lapply(results_table$term, function(term) {
      u <- unlist(strsplit(term, "_", fixed = TRUE))[1]
      mapping[mapping$uniprot == u, 'name']
    }))
  }
  print(head(results_table))
  
  # Create the regulatory protein vs. gene target table
  matrix <- data.frame(matrix(ncol = length(unique(results_table$R_i.name)),
                              nrow = length(unique(results_table$T_k.name))))
  colnames(matrix) <- unique(results_table$R_i.name)
  rownames(matrix) <- unique(results_table$T_k.name)
  
  # Fill in this table with t-statistics
  for (i in 1:nrow(results_table)) {
    tstat <- results_table$statistic[i]
    #beta <- results_table$estimate[i]
    regprot <- results_table$R_i.name[i]
    targ <- results_table$T_k.name[i]
    
    tryCatch({
      matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- tstat
      #matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- beta
    }, error=function(cond){
      print(cond)
      matrix[rownames(matrix) == targ, colnames(matrix) == regprot] <- 0
    })
  }
  # NOTE: leave all the unfilled pairings as NA
  print(head(matrix))
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
                  density.info = "none", key = TRUE, key.title = NA, key.xlab = "T-statistic")
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
  new_hm <- heatmap.2(matrix, scale = "none", col = colorRampPalette(c())(50), trace = "none", 
                      density.info = "none", dendrogram = "row", Rowv = dendrogram, 
                      Colv = "none", key.title = NA, key.xlab = "T-statistic", 
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

master_df_mut_sigTargs <- master_df_mut[master_df_mut$q.value < 0.2, 'T_k.name']
master_df_mut_onlySigTargs <- master_df_mut[master_df_mut$T_k.name %in% master_df_mut_sigTargs,]

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
    #master_df_her2 <- master_df_her2[master_df_her2$q.value < 0.2,]
  }
  
  matrix_lumA <- NA
  matrix_lumB <- NA
  matrix_basal <- NA
  
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
    matrix_lumB <- cbind(matrix_lumB, as.matrix(data.table("CDH1" = rep(0, times = nrow(matrix_lumB)))))
    lumB_incl <- TRUE
  }
  
  basal_incl <- FALSE
  if((nrow(master_df_basal) > 0)) { #& (length(unique(master_df_basal$R_i.name)) >= 2)) {
    matrix_basal <- create_regprot_v_genetarg_matrix(master_df_basal)
    matrix_basal <- cbind(matrix_basal, rep(3, times = nrow(matrix_basal)))
    colnames(matrix_basal)[ncol(matrix_basal)] <- "subtype"
    matrix_basal<- cbind(matrix_basal, as.matrix(data.table("PIK3CA" = rep(0, times = nrow(matrix_basal)))))
    matrix_basal <- cbind(matrix_basal, as.matrix(data.table("CDH1" = rep(0, times = nrow(matrix_basal)))))
    matrix_basal <- cbind(matrix_basal, as.matrix(data.table("KMT2C" = rep(0, times = nrow(matrix_basal)))))
    basal_incl <- TRUE
  }
  
  # her2_incl <- FALSE
  # if((nrow(master_df_her2) > 0) & (length(unique(master_df_her2$R_i.name)) >= 2)) {
  #   matrix_her2 <- create_regprot_v_genetarg_matrix(master_df_her2)
  #   matrix_her2 <- cbind(matrix_her2, rep(4, times = nrow(matrix_her2)))
  #   colnames(matrix_her2)[ncol(matrix_her2)] <- "subtype"
  #   her2_incl <- TRUE
  # }
  
  
  # Combine all these into one
  matrices <- list(matrix_lumA, matrix_lumB, matrix_basal) #, matrix_her2)
  incl_vect <- c(lumA_incl, lumB_incl, basal_incl) #, her2_incl)
  incl_vect_num <- unlist(lapply(1:length(incl_vect), function(i) 
    ifelse(incl_vect[i] == TRUE, i, NA)))
  incl_vect_num <- incl_vect_num[!is.na(incl_vect_num)]
  matrices_to_keep <- matrices[incl_vect_num]
  
  # Ensure all columns are the same
  #unique_cols <- unique(c(unlist(lapply(matrices, function(m) return(colnames(m))))))
  #print(unique_cols)
  #matrices <- lapply(matrices, function(m) {
  #  if(length(m) > 1) {
  #    if(length(unique_cols) != ncol(m)) {
  #      missing_cols <- setdiff(unique_cols, colnames(m))
  #      print(missing_cols)
  #      missing_df <- data.frame(matrix(nrow = nrow(m), ncol = length(missing_cols)))
  #      missing_df[is.na(missing_df)] <- 0
  #      colnames(missing_df) <- missing_cols
  #      m <- cbind(m, missing_df)
  #      print(head(m))
  #    }
  #    return(m)
  #  }
  #})
  
  
  print(head(matrices_to_keep))
  master_matrix <- do.call(rbind, matrices_to_keep)
  print(head(master_matrix))
  print(unique(master_matrix[,'subtype']))
  #master_matrix <- rbind(rbind(rbind(matrix_lumA, matrix_lumB), matrix_basal), matrix_her2)
  
  # Get row clustering
  #row_dend <- hclust(dist(master_matrix[,1:2]))
  hm_annot <- HeatmapAnnotation(df = data.frame(subtype = c(rep("Luminal A", nrow(matrix_lumA)), 
                                                            rep("Luminal B", nrow(matrix_lumB)),
                                                            rep("Basal", nrow(matrix_basal)))),
                                                            #rep("HER2", nrow(matrix_her2)))),
                                col = list(subtype = c("Luminal A" = "#BC3C29FF", "Luminal B" = "#0072B5FF",
                                                       "Basal" = "#20854EFF")), #, "HER2" = "#E18727FF")),
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
                                   res_metabol_her2_mut, FALSE)

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
### GENE SET ENRICHMENT ANALYSIS HEAT MAPS, METABOLISM
##############################################################################
#' Conduct a K-S enrichment test on each of these metabolic pathways 
#' Color each box by the enrichment ratio, and outline the box if p is less than
#' some threshold.
#' @param kegg_pws_list a list of KEGG pathways with their entries
#' @param list_of_master_dfs a list of master DFs from various cancer types, with the 
#' names corresponding to the cancer type name
#' @param goi the gene name of the driver of interest
#' @param one_or_two_sided either "one" or "two" to indicate if we are doing a 
#' one- or two-sided K-S test
create_gsea_heatmap <- function(kegg_pws_list, list_of_master_dfs, goi, one_or_two_sided) {
  
  matrix <- data.frame(matrix(nrow = length(list_of_master_dfs), ncol = length(kegg_pws_list)))
  signif_matrix <- data.frame(matrix(nrow = length(list_of_master_dfs), ncol = length(kegg_pws_list)))
  
  for(i in 1:length(list_of_master_dfs)) {
    df <- list_of_master_dfs[[i]]
    # Restrict to just the GOI
    df_goi <- df[df$R_i.name == goi, ]
    df_goi <- df_goi[order(df_goi$q.value), ]
    
    for (j in 1:length(kegg_pws_list)) {
      genes <- unlist(strsplit(kegg_pws_list[[j]], ";", fixed = T))
      genes <- intersect(genes, df_goi$T_k.name)   # Restrict only to genes we tested
      
      if(length(genes) > 0) {
        # Perform the K-S test
        # Create a uniform distribution the same length as the number of genes
        uniform <- 1:length(df_goi$T_k.name)
        
        # Get the ranks of the known targets within our ranked list
        known_targs <- setdiff(genes, setdiff(genes, df_goi$T_k.name))
        ranks_of_targs <- unlist(lapply(genes, function(x) {
          if(x %fin% df_goi$T_k.name) {return(min(which(df_goi$T_k.name == x,)))}
          else {return(NA)}
        }))
        ranks_of_targs <- ranks_of_targs[!(is.infinite(ranks_of_targs) | (is.na(ranks_of_targs)))]
        
        # Perform the K-S test, with the uniform distribution
        tryCatch({
          if(one_or_two_sided == "one") {
            ks_res <- ks.test(ranks_of_targs, uniform, alternative = "greater")
          } else {
            ks_res <- ks.test(ranks_of_targs, uniform, alternative = "two.sided")
          }
          print(ks_res)
          print(ks_res$p.value)
          
          # Add the statistic to the matrix
          #matrix[i,j] <- ks_res$statistic
          # Get the direction of regulation of most of the genes in the pathway; just the
          # significant hits?
          #gene_betas <- as.numeric(unlist(df_goi[(df_goi$T_k.name %fin% known_targs) & (df_goi$q.value < 0.1), 
                                                 #'estimate']))
          gene_betas <- as.numeric(unlist(df_goi[(df_goi$T_k.name %fin% known_targs), 'estimate']))
          print(head(gene_betas))
          dir <- 1
          #if(length(gene_betas[gene_betas > 0]) < length(gene_betas[gene_betas < 0])) {dir <- -1} 
          if(sum(gene_betas) < 0) {dir <- -1}
          matrix[i,j] <- (-log2(ks_res$p.value)) * dir
          #matrix[i,j] <- ks_res$p.value
          signif_matrix[i,j] <- ks_res$p.value
        }, error = function(cond) {
          print(cond)
          matrix[i,j] <- NA
          signif_matrix[i,j] <- NA
        })
      } else {
        matrix[i,j] <- NA
        signif_matrix[i,j] <- NA
      }
    }
  }
  rownames(matrix) <- names(list_of_master_dfs)
  kegg_labels <- unlist(lapply(names(kegg_pws_list), function(x) 
    unlist(strsplit(x, " [", fixed = T))[1]))
  colnames(matrix) <- kegg_labels
  
  rownames(signif_matrix) <- names(list_of_master_dfs)
  colnames(signif_matrix) <- kegg_labels
  
  cols_to_keep <- colSums(is.na(matrix))<nrow(matrix)
  matrix <- matrix[,cols_to_keep]
  signif_matrix <- signif_matrix[,cols_to_keep]
  
  matrix[is.na(matrix)] <- 0
  
  rows_to_keep <- !(rowSums(matrix) == 0)
  matrix <- matrix[rows_to_keep,]
  signif_matrix <- signif_matrix[rows_to_keep,]
  print(head(matrix))
  
  # Optionally subset the matrix only to those pathways that are p < 0.05 in at least one ct
  #cols_w_signif_val <- as.logical(colSums(signif_matrix < 0.05) >= 1)
  #cols_w_signif_val[is.na(cols_w_signif_val)] <- FALSE
  #matrix_plot <- matrix[, cols_w_signif_val]
  #print(matrix_plot)
  # otherwise, 
  matrix_plot <- matrix
  
  # Create the heatmap itself
  # Clusters by default using hclust, but can specify others using param 'hclustfun'
  #hm <- heatmap.2(as.matrix(matrix), scale = "none", col = bluered(100), trace = "none",
  #                density.info = "none", dendrogram = "row", Colv = TRUE, 
  #                Rowv = TRUE, key = TRUE, key.title = NA, key.xlab = "Enrichment Ratio")
  hm <- pheatmap(matrix_plot, angle_col = "45", main = "", cellheight=15, cellwidth=12,
                 #legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, max(matrix)), 
                 #legend_labels = c("0", "0.2", "0.4", "0.6", "0.8", "1", "-log2(pval)\n\n"),
                 color=colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(50),
                 fontsize_row = 10, fontsize_col = 8) #,filename = paste0(getwd(), "/top_0.05_drivers_kegg_gsea_hm.png"), width = 400, height = 100) #legend_labels = "-log(pval) * direction")
  print(hm)
  #grid.text("KEGG Metabolic Pathway", y=-0.1, gp=gpar(fontsize=16))
  #grid.text("Cancer Type", x=-0.1, rot=90, gp=gpar(fontsize=16))
  
  return(list(matrix, signif_matrix))
}

# Call function
#The 'mar' argument of 'par' sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'
#Default is par(mar=c(5,4,4,1)+.1)
#par(mar = c(6, 6, 4, 2) + 0.1)
output_matrices <- create_gsea_heatmap(kegg_metabolic_pws_list, top_drivers_0.1, "TP53", "two")

# Q-value correction
signif_matrix <- output_matrices[[2]]
pvals <- as.numeric(unlist(as.data.frame(signif_matrix)))
qvals <- qvalue(pvals)$qvalues
padj_vals <- p.adjust(pvals)
signif_matrix_corr <- as.data.frame(matrix(qvals, nrow = nrow(signif_matrix), ncol = ncol(signif_matrix)))
rownames(signif_matrix_corr) <- rownames(signif_matrix)
colnames(signif_matrix_corr) <- colnames(signif_matrix)
colnames(signif_matrix_corr)[which(signif_matrix_corr < 0.1, arr.ind = T)[2]]

colnames(signif_matrix)[which(signif_matrix < 0.05, arr.ind = T)[,'col']]

matrix <- output_matrices[[1]]
print(matrix[6,19])

# Import the KEGG set of metabolic pathways
kegg_metabolic_pws <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/KEGG/kegg_pathways_metabolism.keg.txt", header = T, sep = "\t")

# Adjust this file for proper parsing
curr_pw_broad <- ""   # the broad, B level pathway
curr_pw_narrow <- ""  # the narrow, C level pathway
kegg_metabolic_pws_list <- list()
tmp_narrow_pw_list <- list()
tmp_entry_list <- c()

for(i in 1:nrow(kegg_metabolic_pws)) {
  entry <- kegg_metabolic_pws[i, 1]
  entry_spl <- unlist(strsplit(entry, " ", fixed = T))
  entry_spl <- entry_spl[!(entry_spl == "")]
  # Broad pathway label
  if(entry_spl[1] == "B") {
    if(!(curr_pw_broad == "")) {
      kegg_metabolic_pws_list <- append(kegg_metabolic_pws_list, tmp_narrow_pw_list)
      tmp_narrow_pw_list <- list()
    }
    curr_pw_broad <- paste(entry_spl[3:length(entry_spl)], collapse = " ")
    #names(kegg_metabolic_pws_list)[length(kegg_metabolic_pws_list)] <- curr_pw_broad
    print(curr_pw_broad)
    print(kegg_metabolic_pws_list)
  # Narrow pathway label
  } else if (entry_spl[1] == "C") {
    if(!(curr_pw_narrow == "")) {
      print(head(tmp_entry_list))
      tmp_narrow_pw_list <- append(tmp_narrow_pw_list, paste(unique(tmp_entry_list), collapse = ";"))
      tmp_entry_list <- c()
    }
    curr_pw_narrow <- paste(entry_spl[3:length(entry_spl)], collapse = " ")
    print(curr_pw_narrow)
    names(tmp_narrow_pw_list)[length(tmp_narrow_pw_list)] <- curr_pw_narrow
    print(tmp_narrow_pw_list)
  # Gene entry
  } else {
    gene <- unlist(strsplit(unlist(strsplit(entry, ";", fixed = T))[1], " ", fixed = T))
    gene <- gene[!(gene == "")]
    gene <- toupper(gene[length(gene)])
    if(!grepl(".", gene, fixed = T)) {
      if(grepl("_", gene, fixed = T)) {
        gene <- unlist(strsplit(gene, "_", fixed = T))
        gene[2:length(gene)] <- unlist(lapply(gene[2:length(gene)], function(x) {
          core <- gsub("^\\d+|\\d+$", "", gene[1])  
          return(paste0(core, x))
        }))
      }
      tmp_entry_list <- c(tmp_entry_list, gene)
      print(gene)
    }
  }
}

# Get the intersection of the genes in all these pathways with our metabolic set
# from Recon3D
kegg_genes <- unique(unlist(strsplit(paste(as.character(unlist(kegg_metabolic_pws_list)), 
                                           collapse = ";"), ";", fixed = T)))  # 2904
recon3d_genes_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/metabolic_targets.csv", 
                          header = T, check.names = F)
recon3d_genes <- unique(unlist(lapply(recon3d_genes_df$ensg, function(x) 
  paste(unique(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == x, 'external_gene_name']), collapse = ";"))))  #1753 genes
length(intersect(kegg_genes, recon3d_genes))  #371 oof

# Same as above, but uses Josh's method for ranking p-values by negative log * sign
# of Beta estimate prior to performing the enrichment test

create_gsea_files <- function(list_of_master_dfs, goi) {
  
  for(i in 1:length(list_of_master_dfs)) {
    df <- list_of_master_dfs[[i]]
    
    if(goi %fin% df$R_i.name) {
      # Restrict to just the GOI
      df_goi <- df[df$R_i.name == goi, ]
      
      # Set up a symbol to Entrez ID mapping and add
      #mapping <- as.data.frame(bitr(df_goi$T_k.name, fromType = "SYMBOL", toType = "ENTREZID",
                                    #OrgDb=org.Hs.eg.db, drop = TRUE))
      #colnames(mapping) <- c("T_k.name", "T_k.entrez")
      mapping <- as.data.frame(bitr_kegg(df_goi$T_k, fromType = "uniprot", toType = "kegg",
                                              drop = TRUE, organism = "hsa"))
      colnames(mapping) <- c("T_k", "T_k.kegg")

      #df_goi <- merge(df_goi, mapping, all=TRUE, by="T_k.name")
      df_goi <- merge(df_goi, mapping, all=TRUE, by="T_k")
      
      # Add the negative log p-value to the DF
      df_goi$negLogPval <- unlist(lapply(1:nrow(df_goi), function(i) {
        pval <- df_goi$p.value[i]
        estimate <- df_goi$estimate[i]
        est_sign <- ifelse(estimate > 0, 1, -1)
        return((-log2(pval)) * est_sign)
      }))
      df_goi <- df_goi[order(df_goi$negLogPval, decreasing = T), ]
      df_goi <- na.omit(df_goi)
      
      scores <- df_goi$negLogPval
      #names(scores) <- df_goi$T_k.entrez
      names(scores) <- df_goi$T_k.kegg
      
      #gse.go <- gseGO(scores, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, 
                      #pAdjustMethod = "BH", keyType = "ENTREZID", ont = "BP")
      #gse.go <- setReadable(gse.go, org.Hs.eg.db)
      #print(head(gse.go@result))

      #write.csv(gse.go@result, paste0("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/GSEA_Results/", 
               #paste(goi, paste0(names(list_of_master_dfs)[i], "_gsea_go_results.csv"), sep = "/")))
      
      gse.kegg <- gseKEGG(scores, pvalueCutoff = 1, 
                          pAdjustMethod = "BH", keyType = "kegg")
      print(head(gse.kegg@result))
      write.csv(gse.kegg@result, paste0("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/GSEA_Results/", 
                                      paste(goi, paste0(names(list_of_master_dfs)[i], "_gsea_kegg_results.csv"), sep = "/")))
    }
  }
}

create_gsea_files(top_drivers_0.05, "CTNNB1")


# Can get a strange error about connection trouble; https://stackoverflow.com/questions/64519640/error-in-summary-connectionconnection-invalid-connection
# Use this function to eliminate this error

# Alternatively, make the heat map from pre-constructed gsea files
create_gsea_heatmap_from_file <- function(list_of_gsea_files, padj_thres, terms) {
  
  list.gsea <- list()
  list.gsea.signif <- list()
  
  signif_pws <- c()
  
  for(i in 1:length(list_of_gsea_files)) {
    ct <- names(list_of_gsea_files)[i]
    f <- list_of_gsea_files[[i]]
    
    f <- f[f$Description %fin% terms, ]
    signif_pws <- c(signif_pws, as.character(unlist(f[f$p.adjust < padj_thres, 
                                                      'Description'])))
    
    matrix <- data.frame("description" = f$Description, 
                         "enrichment.score" = f$enrichmentScore)
    colnames(matrix)[2] <- paste0("enrichment.score.", ct)
    signif_matrix <- data.frame("description" = f$Description, 
                                "pval" = f$p.adjust)
    colnames(signif_matrix)[2] <- paste0("pval.", ct)
    
    print(head(matrix))
    list.gsea[[i]] <- matrix
    list.gsea.signif[[i]] <- signif_matrix
  }
  names(list.gsea) <- names(list_of_gsea_files)

  matrix <- Reduce(function(x, y) merge(x, y, by="description"), list.gsea)
  signif_matrix <- Reduce(function(x, y) merge(x, y, by="description"), list.gsea.signif)
  print(length(signif_pws))
  
  matrix <- matrix[matrix$description %fin% signif_pws,]
  signif_matrix <- signif_matrix[signif_matrix$description %fin% signif_pws,]
  
  rownames(matrix) <- matrix$description
  rownames(signif_matrix) <- signif_matrix$description
  matrix <- matrix[,colnames(matrix) != "description"]
  signif_matrix <- signif_matrix[,colnames(signif_matrix) != "description"]
  #print(head(matrix))
  
  matrix <- t(matrix)
  signif_matrix <- t(signif_matrix)
  
  rownames(matrix) <- names(list_of_gsea_files)
  rownames(signif_matrix) <- names(list_of_gsea_files)
  
  print(head(signif_matrix))
  
  # Alternatively, keep only those significant in at least 2 cancer types
  cols_to_keep <- colSums(is.na(matrix))<nrow(matrix)
  cols_to_keep <- cols_to_keep & (colSums(signif_matrix < padj_thres) >= 2)
  matrix <- matrix[,cols_to_keep]
  signif_matrix <- signif_matrix[,cols_to_keep]

  #matrix[is.na(matrix)] <- 0
  head(matrix)
  print(rowSums(matrix))
  
  rows_to_keep <- !(rowSums(matrix) == 0)
  matrix <- matrix[rows_to_keep,]
  signif_matrix <- signif_matrix[rows_to_keep,]
  #print(head(matrix))
  
  
  rc1 <- colorRampPalette(c("#0072B5FF", "white"))(50)
  rc2 <- colorRampPalette(c("white", "#BC3C29FF"))(50)           
  rampcols <- c(rc1, rc2)
  
  rb1 <- seq(min(matrix, na.rm = T), 0, length.out = 51)
  rb2 <- seq(0, max(matrix, na.rm = T), length.out = 51)[-1]
  rampbreaks <- c(rb1, rb2)
  hm <- pheatmap(matrix, angle_col = "45", main = "", cellheight=25, cellwidth=18,
           color= rampcols, breaks = rampbreaks, bias = 1, fontsize_row = 12, 
           fontsize_col = 10, na_col = "gray90") 
  
  #hm <- pheatmap(matrix, angle_col = "45", main = "", cellheight=15, cellwidth=12,
                 #color=colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(50),
                 #fontsize_row = 10, fontsize_col = 8) 
  #ylab = "\n\nCancer Type", xlab = "KEGG Metabolic Pathway") #, cellheight = 25, cellwidth = 25)
  
  print(hm)
}


# Import the GSEA files
gsea_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/GSEA_Results/"
kras_gsea_fns <- list.files(paste0(gsea_path, "KRAS/"), include.dirs = F, recursive = F, pattern = "gsea_go_results.csv")
kras_gsea_fns <- list.files(paste0(gsea_path, "KRAS/"), include.dirs = F, recursive = F, pattern = "gsea_kegg_results.csv")

kras_gsea_files <- lapply(kras_gsea_fns, function(x) fread(paste0(gsea_path, paste0("KRAS/", x)), header = T))
# NOTE: add cancer type names
names(kras_gsea_files) <- unlist(lapply(kras_gsea_fns, function(x) unlist(strsplit(x, "_", fixed = T))[1]))
names(tp53_gsea_files) <- unlist(lapply(tp53_gsea_fns, function(x) unlist(strsplit(x, "_", fixed = T))[1]))

# For GO, import a list of metabolic GO terms (retrieved from https://www.informatics.jax.org/go/term/GO:0008152)
go_metabolic_terms <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/GO/GO_term_summary_20230506_145527.csv",
                               header = T, check.names = F)$`Annotated Term`

create_gsea_heatmap_from_file(kras_gsea_files, 0.01, go_metabolic_terms)

kegg_metabolic_pws <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/KEGG/kegg_pathways_metabolism.keg.txt", 
                                 header = F, sep = "\t")
kegg_metabolic_pw_names <- unlist(lapply(kegg_metabolic_pws$V1, function(x) {
  x <- unlist(strsplit(x, " [", fixed = T))[1]
  spl <- unlist(strsplit(x, " ", fixed = T))
  spl <- spl[spl != ""]
  if(!(spl[1] %in% c("A", "D"))) {
    #spl <- paste(spl, collapse = "")
    #new <- paste(unlist(strsplit(spl, ""))[-c(1,2,3,4,5,6)], collapse = "")
    print(spl)
    spl <- spl[spl != "NA"]
    new <- paste(spl[3:length(spl)], collapse = " ")
    print(new)
    return(new)
  }
}))
kegg_metabolic_pw_names <- kegg_metabolic_pw_names[kegg_metabolic_pw_names != "NA Metabolism"]
write.table(kegg_metabolic_pw_names, "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/KEGG/kegg_pathways_metabolism.txt",
            quote = F, row.names = F)

kegg_metabolic_pw_names <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/KEGG/kegg_pathways_metabolism.txt", 
                                      header = T, sep = "\n")[,1]

create_gsea_heatmap_from_file(kras_gsea_files, 0.1, kegg_metabolic_pw_names)

# Alternatively:
create_gsea_heatmap_from_file(kras_gsea_files, 0.1, unique(kras_pc_gsea_kegg$Description))


#Specific terms of interest
terms_of_interest_kras <- c("pyrimidine nucleoside catabolic process", "ceramide metabolic process",
                            "sphingolipid metabolic process", "phospholipid catabolic process", "collagen metabolic process")


# Validate the genes in a particular pathway
pyrimidine_pw_genes <- c("AICDA", "APOBEC2", "CDA", "CDADC1", "DCK", "DPYD", 
                         "DTYMK", "PYCRL", "TK1", "TK2", "TYMP", "UPB1", "UPP1", "UPP2")
pyrimidine_pw_genes_kegg <- c("NT5E", "CANT1", "TK2", "ENTPD4", "CTPS2", "NME4", 
                              "NME2", "NME1", "DCTD", "ENPP3", "NT5C", "NME6", "NUDT2")
pyrimidine_pw_genes_mkegg <- c("CTPS2", "NME4", "NME2", "NME6", "NME1", "NME7", 
                               "RRM1", "RRM2B", "TYMS", "RRM2", "DUT", "DTYMK",
                               "CMPK2", "CTPS1")

unique_pyrimidine_set <- unique(c(pyrimidine_pw_genes, pyrimidine_pw_genes_kegg, 
                                  pyrimidine_pw_genes_mkegg))
unique_pyrimidine_set_sig <- unlist(lapply(unique_pyrimidine_set, function(g) {
  sig_v <- c()
  for (i in 1:length(top_drivers_0.05)) {
    x <- top_drivers_0.05[[i]]
    if("KRAS" %fin% x$R_i.name) {
      if(g %fin% x$T_k.name) {
        if(x[(x$T_k.name == g) & (x$R_i.name == "KRAS"), 'estimate'] < 0.2) {
          sig_v <- c(sig_v, T)
        }
        else {sig_v <- c(sig_v, F)}
      }
    }
  }
  if(length(sig_v[sig_v == T]) < 1) {return(NA)}
  else {return(g)}
}))
unique_pyrimidine_set_sig <- unique_pyrimidine_set_sig[!is.na(unique_pyrimidine_set_sig)]

nucleotide_metabolism_kegg <- c("NT5E", "CANT1", "TK2", "ENTPD4", "CTPS2", "GMPR",
                                "AK1", "ENTPD2", "NME4", "NME2", "DCTD", "GUK1", 
                                "PNP", "ENPP3", "NT5C", "NME6")
matrix_pyrimidine_list <- lapply(top_drivers_0.05[names(top_drivers_0.05) %in% 
                                                    c("LUAD", "COAD", "PAAD", "UCEC")], function(x) {
  vals <- unlist(lapply(pyrimidine_pw_genes, function(g) {
    if(!(g %fin% x$T_k.name)){return(NA)}
    est <- as.numeric(x[(x$T_k.name == g) & (x$R_i.name == "KRAS"), 'estimate'])
    negLogQval <- -log10(as.numeric(x[(x$T_k.name == g) & (x$R_i.name == "KRAS"), 'q.value']))
    if(!is.na(est)) {
      if(est > 0) {return(1*negLogQval)}
      else {return(-1*negLogQval)}
    } else {return(NA)}
  }))
  return(vals)
})
matrix_pyrimidine <- do.call(rbind, matrix_pyrimidine_list)
colnames(matrix_pyrimidine) <- pyrimidine_pw_genes
matrix_pyrimidine <- matrix_pyrimidine[,colSums(is.na(matrix_pyrimidine)) < 
                                         (nrow(matrix_pyrimidine)*0.5)]

matrix_pyrimidine_sub <- matrix_pyrimidine[, colSums(abs(matrix_pyrimidine) > 2, na.rm = T) >= 1]
matrix_pyrimidine_sub <- matrix_pyrimidine_sub[rowSums(abs(matrix_pyrimidine_sub) > 2, na.rm = T) >= 1, ]

# Get the colors to be centered on 0
rc1 <- colorRampPalette(c("#0072B5FF", "white"))(50)
rc2 <- colorRampPalette(c("white", "#BC3C29FF"))(50)           
rampcols <- c(rc1, rc2)

#rb1 <- seq(min(matrix_pyrimidine_sub, na.rm = T), 0, length.out = 51)
rb1 <- seq(-max(matrix_pyrimidine_sub, na.rm = T), 0, length.out = 51)
rb2 <- seq(0, max(matrix_pyrimidine_sub, na.rm = T), length.out = 51)[-1]
rampbreaks <- c(rb1, rb2)
pheatmap(matrix_pyrimidine_sub, angle_col = "45", main = "",
         color= rampcols, breaks = rampbreaks, bias = 1, na_col = "grey90") 


branched_chain_aa_catabolic_process_genes <- c("ACAD8", "AUH", "BCAT1", "BCAT2", "BCKDHA", "BCKDHD", "BCKDK",
                                               "DBT", "DLD", "HIBADH", "HIBCH", "IVD", "SLC25A44", "ACADSB",
                                               "ACAT1", "ALDH6A1", "HMGCL", "HMGCLL1", "HSD17B10", "MCCC2")
matrix_aa_list <- lapply(top_drivers_0.05[!(names(top_drivers_0.05) %in% c("CESC", "KIRC", "KIRP", "PCPG", "PRAD", "SARC", "THCA"))], function(x) {
  vals <- unlist(lapply(branched_chain_aa_catabolic_process_genes, function(g) {
    est <- as.numeric(x[(x$T_k.name == g) & (x$R_i.name == "TP53"), 'estimate'])
    negLogPval <- -log2(as.numeric(x[(x$T_k.name == g) & (x$R_i.name == "TP53"), 'p.value']))
    print(est)
    if(!is.na(est)) {
      if(est > 0) {return(1*negLogPval)}
      else {return(-1*negLogPval)}
    } else {return(NA)}
  }))
  return(vals)
})
matrix_aa <- do.call(rbind, matrix_aa_list)
colnames(matrix_aa) <- branched_chain_aa_catabolic_process_genes
matrix_aa <- matrix_aa[,colSums(is.na(matrix_aa))<nrow(matrix_aa)]

matrix_aa_sub <- matrix_aa[, colSums(abs(matrix_aa) > 8, na.rm = T) >= 1]
matrix_aa_sub <- matrix_aa_sub[rowSums(abs(matrix_aa_sub) > 8, na.rm = T) >= 1, ]

rb1 <- seq(min(matrix_aa_sub, na.rm = T), 0, length.out = 51)
rb2 <- seq(0, max(matrix_aa_sub, na.rm = T), length.out = 51)[-1]
rampbreaks <- c(rb1, rb2)
pheatmap(matrix_aa_sub, angle_col = "45", main = "",
         color= rampcols, breaks = rampbreaks, bias = 1) 


# Same as above, but for STRINGdb-based enrichment
create_gsea_files_stringdb <- function(list_of_master_dfs, goi, background, qval_thres) {
  
  # Instantiate the string DB object 
  string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=0, 
                            input_directory="")
  
  # Set the background, if needed
  if(!is.na(background)) {
    background <- as.data.frame(background)
    colnames(background) <- "gene"
    background <- string_db$map(background, "gene", removeUnmappedRows = TRUE)
    
    string_db <- string_db$set_background(background$STRING_id)
  }
  
  for(i in 1:length(list_of_master_dfs)) {
    df <- list_of_master_dfs[[i]]
    
    if(goi %fin% df$R_i.name) {
      # Restrict to just the GOI
      df_goi <- df[df$R_i.name == goi, ]
      df_goi_sub <- df_goi[df_goi$q.value < qval_thres, ]
      
      # Calculate the enrichment for all sources  
      if(nrow(df_goi_sub) > 1) {
        stringdb_enrichment <- string_db$get_enrichment(df_goi_sub$T_k.name)
        #print(head(stringdb_enrichment, n=20))
        
        write.csv(stringdb_enrichment, paste0("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/GSEA_Results/", 
                                              paste(goi, paste0(names(list_of_master_dfs)[i], "_gsea_stringdb_results.csv"), sep = "/")))
      } 
    }
  }
}

create_gsea_files_stringdb(top_drivers_0.05, "TP53", NA, 0.1)

create_gsea_heatmap_from_file_stringdb <- function(list_of_gsea_files, fdr_thres, db) {
  
  list.gsea <- list()

  signif_pws <- c()
  cts_added <- c()
  
  for(i in 1:length(list_of_gsea_files)) {
    ct <- names(list_of_gsea_files)[i]
    f <- list_of_gsea_files[[i]]
    
    f <- f[f$category == db, ]
    if(nrow(f > 0)) {
      signif_pws <- c(signif_pws, as.character(f$description))
      if(!is.na(fdr_thres)) {
        signif_pws <- c(signif_pws, as.character(unlist(f[f$fdr < fdr_thres, 'description'])))
      }
      #signif_matrix <- data.frame("description" = f$description, 
                                  #"fdr" = f$fdr)
      signif_matrix <- data.frame("description" = f$description, 
                                  "fdr" = 1)
      colnames(signif_matrix)[2] <- paste0("fdr.", ct)
      
      print(head(signif_matrix))
      list.gsea[[i]] <- signif_matrix
      cts_added <- c(cts_added, names(list_of_gsea_files)[i])
    } else {list.gsea[[i]] <- NA}
  }
  list.gsea <- list.gsea[!is.na(list.gsea)]
  names(list.gsea) <- cts_added
  print(list.gsea)
  
  matrix <- Reduce(function(x, y) merge(x, y, by="description", all = T), list.gsea)
  #print(length(signif_pws))
  
  matrix <- matrix[matrix$description %fin% signif_pws,]
  rownames(matrix) <- matrix$description
  matrix <- matrix[,colnames(matrix) != "description"]
  print(head(matrix))
  
  matrix <- t(matrix)
  print(head(matrix))
  
  rownames(matrix) <- cts_added
  
  print(head(matrix))

  cols_to_keep <- colSums(is.na(matrix))<nrow(matrix)
  matrix <- matrix[,cols_to_keep]

  matrix[is.na(matrix)] <- 0
  
  # Keep only pathways where at least 2 cancer types have nonzero values
  matrix <- matrix[, colSums(matrix > 0) >= 2]
  
  rows_to_keep <- !(rowSums(matrix) == 0)
  matrix <- matrix[rows_to_keep,]
  
  print(head(matrix))
  
  hm <- pheatmap(matrix, angle_col = "45", main = "", cellheight=20, cellwidth=17,
                 color=colorRampPalette(c("white", "#BC3C29FF"))(50),
                 fontsize_row = 10, fontsize_col = 8, legend = F) 
  #ylab = "\n\nCancer Type", xlab = "KEGG Metabolic Pathway") #, cellheight = 25, cellwidth = 25)
  
  print(hm)
  
  return(matrix)
}


# Import the GSEA files
gsea_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/GSEA_Results/"
kras_gsea_fns <- list.files(paste0(gsea_path, "KRAS/"), include.dirs = F, recursive = F, pattern = "gsea_stringdb_results.csv")
kras_gsea_files <- lapply(kras_gsea_fns, function(x) fread(paste0(gsea_path, paste0("KRAS/", x)), header = T))
# NOTE: add cancer type names
names(kras_gsea_files) <- unlist(lapply(kras_gsea_fns, function(x) unlist(strsplit(x, "_", fixed = T))[1]))
names(tp53_gsea_files) <- unlist(lapply(tp53_gsea_fns, function(x) unlist(strsplit(x, "_", fixed = T))[1]))

create_gsea_heatmap_from_file_stringdb(tp53_gsea_files, NA, "WikiPathways")


