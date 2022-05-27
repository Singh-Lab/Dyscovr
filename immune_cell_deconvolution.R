############################################################
### Immune Cell Deconvolution
### Written By: Sara Camilli, April 2021
############################################################

#BiocManager::install("EPIC")
library(EPIC)
library(dplyr)
library(parallel)
library(immunedeconv)
library(knitr)
library(RColorBrewer)
library(ggplot2)
library(TCGAbiolinks)

# This file takes in varying types of expression files and uses them to run through
# an immune cell deconvolution pipeline. This is primarily done using the R package
# immunedeconv (https://icbi-lab.github.io/immunedeconv/index.html)

# At the bottom of the file, also presents an alternative method using pre-calculated
# estimates from TIMER on TCGA data
# http://timer.cistrome.org/infiltration_estimation_for_tcga.csv.gz

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/"
#path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/"

output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Expression/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Expression/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


################################################################################
### IMPORT EXPRESSION RAW COUNT MATRIX
################################################################################
expression_dataframe_counts <- read.csv(paste0(output_path, "expression_counts_DF.csv"), header = TRUE, 
                                        row.names = 1, check.names = FALSE)
# For METABRIC, use the unlogged intensity matrix
#expression_dataframe_counts <- read.table(paste0(path, "cBioPortal/brca_metabric/data_mrna_agilent_microarray_unlogged.txt"), header = TRUE, 
                                        #row.names = 1, check.names = FALSE)
#expression_dataframe_counts <- distinct(expression_dataframe_counts)
#rownames(expression_dataframe_counts) <- make.names(expression_dataframe_counts$Hugo_Symbol, unique = TRUE)
#expression_dataframe_counts <- expression_dataframe_counts[, -c(1,2)]
#expression_dataframe_counts <- na.omit(expression_dataframe_counts)

# If desired, subset this expression count matrix to a desired cancer type
gdc_sample_sheet <- read.csv(paste(path, "Expression_Data/Counts/gdc_sample_sheet_total.csv", sep = ""),
                             header = TRUE, check.names = FALSE)

gbm_patients <- gdc_sample_sheet[grepl("GBM", gdc_sample_sheet$Project.ID), 'Sample.ID']
expression_dataframe_gbm <- expression_dataframe_counts[,colnames(expression_dataframe_counts) %fin% 
                                                          gbm_patients]  #52 patients

luad_patients <- gdc_sample_sheet[grepl("LUAD", gdc_sample_sheet$Project.ID), 'Sample.ID']
expression_dataframe_luad <- expression_dataframe_counts[,colnames(expression_dataframe_counts) %fin% 
                                                          luad_patients]  #470 patients

brca_patients <- gdc_sample_sheet[grepl("BRCA", gdc_sample_sheet$Project.ID), 'Sample.ID']
expression_dataframe_brca <- expression_dataframe_counts[,colnames(expression_dataframe_counts) %fin% 
                                                           brca_patients]  

################################################################################
### CREATE TPM EXPRESSION MATRIX FROM COUNT MATRIX
################################################################################
# For several methods, counts must first be converted to TPM (transcripts per million)
# Counts to TPM:
  # Divide counts by the length of each gene in KB (reads per kilobase, RPK)
  # Sum RPK values across whole sample and divide by 1,000,000 to get "per million" scaling factor
  # Divide RPK values by "per million" scaling factor
  # Returns a TPM-normalized dataframe

#' Function to convert counts to TPM (transcripts per million)
#' @param count_df a data frame containing expression counts (columns are patient ID-sample ID, 
#' rows are ENSG IDs)
#' @param gene_lengths a data frame containing the lengths of all the genes (in bp)
convert_counts_to_tpm <- function(count_df, gene_lengths) {
  # Get the lengths (in order) for all the genes in the count dataframe
  ensg_lengths <- unlist(lapply(rownames(expression_dataframe_counts), function(x) {
    if(x %fin% gene_lengths$ensg) {
      return(gene_lengths[gene_lengths$ensg == x, 'length'])
    } else {return(NA)}
  }))
  
  tpm_df_cols <- mclapply(1:ncol(count_df), function(i) {
    print(paste(i, ncol(count_df), sep = "/"))
    # Get this sample's counts
    samp_i_counts <- count_df[,i]
    # For every count in this column, divide by the length of the gene
    rpk_counts <- samp_i_counts / ensg_lengths
    # Divide all values by the per million scaling factor
    scaling_fact <- sum(rpk_counts, na.rm = TRUE) / 1000000
    tpm_counts <- rpk_counts / scaling_fact
    return(tpm_counts)
  })
  tpm_df <- as.data.frame(do.call("cbind", tpm_df_cols))
  rownames(tpm_df) <- rownames(count_df)
  colnames(tpm_df) <- colnames(count_df)
  return(tpm_df)
}

#' Using the BioMart Gene ID conversion data frame, creates and returns a data frame
#' with two columns that have the ensg ID ("ensg") and the length of the gene in BP ("length")
#' @param all_genes_id_conv a master data frame with all genes and their various IDs, along
#' with their start and end positions and other data (from BioMart)
get_gene_length_df <- function(all_genes_id_conv) {
  sub_conv_df <- all_genes_id_conv[,c('ensembl_gene_id', 'start_position', 'end_position')]
  lengths <- unlist(lapply(1:nrow(sub_conv_df), function(i) sub_conv_df[i, 'end_position'] - sub_conv_df[i, 'start_position']))
  gene_length_df <- data.frame("ensg" = sub_conv_df$ensembl_gene_id, "length" = lengths)
  return(distinct(gene_length_df))
}

gene_lengths <- get_gene_length_df(all_genes_id_conv)
tpm_dataframe <- convert_counts_to_tpm(expression_dataframe_counts, gene_lengths)

# Write this TPM dataframe to a TSV/CSV file
write.table(tpm_dataframe, paste(output_path, "expression_tpm_DF.tsv", sep = ""), 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.csv(tpm_dataframe, paste(output_path, "expression_tpm_DF.csv", sep = ""))

# NOTE: This file can then used by quanTIseq on the command line 

# Read back
tpm_dataframe <- read.table(paste(output_path, "expression_tpm_DF.tsv", sep = ""), 
                            sep = "\t", header = FALSE, row.names = 1)


################################################################################
### UPDATE TPM MATRIX WITH GENE NAMES, RATHER THAN ENSG IDs
################################################################################
#' Gene names need to be names, not ENSG IDs; convert all ENSG IDs to gene names 
#' in a given data frame
#' @param df a data frame with gene names as the row names
#' @param all_genes_id_conv a master data frame with all genes and their various IDs, along
#' with their start and end positions and other data (from BioMart)
convert_ensg_to_gn <- function(df, all_genes_id_conv) {
  gene_names <- unlist(lapply(rownames(df), function(x) {
    names <- all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == x, 'external_gene_name']
    if(length(names) > 1) {names <- names[1]}
    if(length(names) == 0) {return(NA)}
    return(names)
  }))
  
  # Remove NA rows
  ind_to_rm <- which(is.na(gene_names))
  ind_to_rm <- c(ind_to_rm, which(duplicated(gene_names)))
  gene_names <- gene_names[-ind_to_rm]
  dataframe_gn <- df[-ind_to_rm,]
  
  # Add gene names as row names
  rownames(dataframe_gn) <- gene_names
  
  return(dataframe_gn)
}

tpm_dataframe_gn <- convert_ensg_to_gn(tpm_dataframe, all_genes_id_conv)
fpkm_dataframe_gn <- convert_ensg_to_gn(expression_dataframe_fpkm, all_genes_id_conv)

# Write these files to DF
write.csv(tpm_dataframe_gn, paste(output_path, "expression_tpm_gn_DF.csv", sep = ""))
write.csv(fpkm_dataframe_gn, paste(output_path, "expression_fpkm_gn_DF.csv", sep = ""))
#
#
# Read back if needed
tpm_dataframe_gn <- read.csv(paste(output_path, "expression_tpm_gn_DF.csv", sep = ""), 
                             header = TRUE, check.names = FALSE, row.names = 1)


################################################################################
### FIX DUPLICATE COLUMN NAMES IN TPM FILE
################################################################################
#' Remove duplicate columns in the expression data frame
#' @param dataframe_gn an expression data frame where the row names are gene names
#' (rather than ENSG IDs)
eliminate_dupl_columns <- function(dataframe_gn) {
  new_colnames <- unlist(lapply(1:ncol(dataframe_gn), function(i) {
    name <- colnames(dataframe_gn)[i]
    if(i %in% which(duplicated(colnames(dataframe_gn)))) {
      name_new <- paste(name, "1", sep = "-")
      return(name_new)
    } else {return(name)}
  }))
  colnames(dataframe_gn) <- new_colnames
  return(dataframe_gn)
}
tpm_dataframe_gn <- eliminate_dupl_columns(tpm_dataframe_gn)
fpkm_dataframe_gn <- eliminate_dupl_columns(fpkm_dataframe_gn)


################################################################################
### OPTIONAL: STRATIFY BY CANCER SUBTYPE
################################################################################

# Import the files with information about patient subtype using TCGAbiolinks
general_cancer_types <- c("ACC", "BRCA", "BLCA", "CESC", "CHOL", "COAD", "ESCA",
                          "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC",
                          "LUAD", "LUSC", "PAAD", "PCPG", "PRAD", "READ", "SKCM",
                          "SARC", "STAD", "THCA", "UCEC", "UCS", "UVM")

import_subtype_information <- function(general_cancer_types) {
  subtype_files <- lapply(general_cancer_types, function(type) {
    subtype <- TCGAquery_subtype(tumor = type)
    return(subtype)
  })
  names(subtype_files) <- general_cancer_types
  
  # Note: To combine into one file, we'll have to do some work to ID matching 
  # columns with differing names
  
  return(subtype_files)
}

# Or, if just using BRCA:
brca_subtype <- TCGAquery_subtype(tumor = "BRCA")

#' Given a TPM data frame and a file with information about patient subtype,
#' stratify the TPM data frame into one data frame per subtype. Return this 
#' as a list of data frames, with the name corresponding to the subtype.
#' @param tpm_dataframe_gn a patient (column) by gene (row) data frame that 
#' contains TPM-normalized expression values
#' @param subtype_file a clinical file from the TCGA that has patient IDs and
#' their corresponding cancer subtypes
stratify_by_subtype <- function(tpm_dataframe, subtype_file) {
  # NOTE: the name of the column will vary by cancer type
  subtypes <- unique(subtype_file$BRCA_Subtype_PAM50)
  
  tpm_dataframe_list <- lapply(subtypes, function(st) {
    # Identify all the patients that have this subtype
    patients <- unlist(subtype_file[subtype_file$BRCA_Subtype_PAM50 == st, 'patient'])
    
    # Subset the TPM data frame to include only these patients
    tpm_adjcols <- unlist(lapply(colnames(tpm_dataframe), function(cn) paste(unlist(strsplit(cn, "-", fixed = TRUE))[1:3], 
                                                                 collapse = "-")))
    tpm_df_adjcols <- tpm_dataframe
    colnames(tpm_df_adjcols) <- tpm_adjcols
    
    cols_to_keep <- which(colnames(tpm_df_adjcols) %fin% patients)
    tpm_sub <- tpm_dataframe[,cols_to_keep]
    
    # Return this subsetted DF 
    return(tpm_sub)
  })
  
  names(tpm_dataframe_list) <- subtypes
  
  return(tpm_dataframe_list)
}

# Call this function (for BRCA)
tpm_dfs_by_st <- stratify_by_subtype(tpm_dataframe_gn, brca_subtype)



################################################################################
### SEPARATE TUMOR AND NORMAL SAMPLES
################################################################################
tpm_dataframe_gn_tumor_only <- tpm_dataframe_gn[,grepl("-0", colnames(tpm_dataframe_gn))]
tpm_dataframe_gn_normal_only <- tpm_dataframe_gn[,grepl("-1", colnames(tpm_dataframe_gn))]

# OPTIONAL: Further segregate the normal samples into 'blood derived normal' and 
# 'solid tissue normal'

# Import the 'sample' file
sample_file <- read.csv(paste(path, "sample.csv"))

################################################################################
### RUN IMMUNEDECONV METHODS
################################################################################
# The immunedeconv R package integrates multiple immune cell deconv. packages, including:
# quantiseq, TIMER, CIBERSORT, MCPCounter, xCell, EPIC


# quanTIseq
results_quantiseq <- immunedeconv::deconvolute(tpm_dataframe_gn_tumor_only, 'quantiseq', tumor = TRUE)
results_quantiseq_tab <- as.data.frame(results_quantiseq)
write.csv(results_quantiseq_tab, paste(output_path, "quanTIseq Cell Deconvolution/quantiseq_cellFrac.csv", sep = ""))

# EPIC
results_epic <- immunedeconv::deconvolute(tpm_dataframe_gn_tumor_only, "epic", tumor = TRUE)
results_epic_tab <- as.data.frame(results_epic)
write.csv(results_epic_tab, paste(output_path, "EPIC Cell Deconvolution/epic_cellFrac_from_immunedeconv.csv", sep = ""))

# xCell
results_xCell <- immunedeconv::deconvolute(tpm_dataframe_gn_tumor_only, "xcell")
write.csv(results_xCell, paste(output_path, "xCell Cell Deconvolution/xcell_cellFrac.csv", sep = ""))

# mcp_counter
results_mcpcounter <- immunedeconv::deconvolute(tpm_dataframe_gn_tumor_only, "mcp_counter")
write.csv(results_xCell, paste(output_path, "MCP-Counter Cell Deconvolution/mcpcounter_cellFrac.csv", sep = ""))

# TIMER
indic <- rep("brca", times = ncol(tpm_dataframe_gn_tumor_only))
results_timer <- immunedeconv::deconvolute(tpm_dataframe_gn_tumor_only, "timer", indications = indic)
write.csv(results_timer, paste(output_path, "TIMER Cell Deconvolution/timer_cellFrac.csv", sep = ""))

# CIBERSORT
# Set location of necessary source code files
cibersort_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/CIBERSORT Source Code/"
set_cibersort_binary(paste(cibersort_path, "CIBERSORT.R", sep = ""))
set_cibersort_mat(paste(cibersort_path, "CIBERSORT_package/LM22.txt", sep = ""))

results_cibersort <- deconvolute(tpm_dataframe_gn_tumor_only, "cibersort")
write.csv(results_cibersort, paste(output_path, "CIBERSORT Cell Deconvolution/cibersort_cellFrac.csv", sep = ""))
results_cibersort_abs <- as.data.frame(deconvolute(tpm_dataframe_gn_tumor_only, "cibersort_abs"))
write.csv(results_cibersort_abs, paste(output_path, "CIBERSORT Cell Deconvolution/cibersort_abs_cellFrac.csv", sep = ""))

# Run a version optimized for microarray data (METABRIC)
output_path_metabric <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/cBioPortal/METABRIC/Expression"
results_cibersort <- deconvolute(expression_dataframe_counts, method = "cibersort", arrays = TRUE)
write.csv(results_cibersort, paste(output_path_metabric, "CIBERSORT Cell Deconvolution/cibersort_cellFrac_metabric.csv", sep = ""))
results_cibersort_abs <- as.data.frame(deconvolute(expression_dataframe_counts, method = "cibersort_abs", arrays = TRUE))
write.csv(results_cibersort_abs, paste(output_path_metabric, "CIBERSORT Cell Deconvolution/cibersort_abs_cellFrac_metabric.csv", sep = ""))



#
#
#
# Read back files if needed
results_cibersort_abs_tumor <- read.csv(paste(output_path, "CIBERSORT Cell Deconvolution/cibersort_abs_cellFrac.csv", sep = ""),
                                        header = TRUE, check.names = FALSE, row.names = 1)
results_cibersort_abs_normal <- read.csv(paste(output_path, "CIBERSORT Cell Deconvolution/cibersort_abs_cellFrac_normalSamp.csv", sep = ""),
                                        header = TRUE, check.names = FALSE, row.names = 1)
results_epic_tumor <- read.csv(paste(output_path, "EPIC Cell Deconvolution/epic_cellFrac_from_immunedeconv.csv", sep = ""),
                                        header = TRUE, check.names = FALSE, row.names = 1)
results_epic_normal <- read.csv(paste(output_path, "EPIC Cell Deconvolution/epic_cellFrac_from_immunedeconv_normalSamp.csv", sep = ""),
                               header = TRUE, check.names = FALSE, row.names = 1)
results_quantiseq_tumor <- read.csv(paste(output_path, "quanTIseq Cell Deconvolution/quantiseq_cellFrac.csv", sep = ""),
                                header = TRUE, check.names = FALSE, row.names = 1)
results_quantiseq_normal <- read.csv(paste(output_path, "quanTIseq Cell Deconvolution/quantiseq_cellFrac_normalSamp.csv", sep = ""),
                                    header = TRUE, check.names = FALSE, row.names = 1)


################################################################################
### CONVERT RESULTS TABLES INTO A FORMAT ACCEPTABLE TO GGPLOT2
################################################################################
# Version 1 
convert_res_for_vis <- function(results_tab) {
  reorg_cols <- lapply(1:nrow(results_tab), function(i) as.numeric(results_tab[i, 2:ncol(results_tab)]))
  reorg_df <- do.call("cbind", reorg_cols)
  colnames(reorg_df) <- results_tab$cell_type
  rownames(reorg_df) <- colnames(results_tab)[2:length(colnames(results_tab))]
  reorg_df <- as.data.frame(reorg_df)
  return(reorg_df)
}

# Version 2 (Most useful for boxplots, etc.)
convert_res_for_boxplot <- function(results_tab) {
  reorg_tabs <- lapply(2:ncol(results_tab), function(i){
    samp_name <- colnames(results_tab)[i]
    samp_vals <- as.numeric(results_tab[,i])
    samp_types <- results_tab$cell_type
    return(data.frame("Sample" = rep(samp_name, times = length(samp_vals)),
                      "Absolute.Cell.Fraction" = samp_vals,
                      "Cell.Type" = samp_types))
  })
  reorg_df <- do.call("rbind", reorg_tabs)
  return(reorg_df)
}

epic_reorg_df <- convert_res_for_boxplot(results_epic_tab)
quantiseq_reorg_df <- convert_res_for_boxplot(results_quantiseq_tab)
cibersort_abs_reorg_df <- convert_res_for_boxplot(as.data.frame(results_cibersort_abs))
xcell_reorg_df <- convert_res_for_boxplot(as.data.frame(results_xCell))
mcpcounter_reorg_df <- convert_res_for_boxplot(as.data.frame(results_mcpcounter))
timer_reorg_df <- convert_res_for_boxplot(as.data.frame(results_timer))
cibersort_reorg_df <- convert_res_for_boxplot(as.data.frame(results_cibersort))

# OR
epic_reorg_df_tumor <- convert_res_for_boxplot(results_epic_tumor)
epic_reorg_df_normal <- convert_res_for_boxplot(results_epic_normal)
quantiseq_reorg_df_tumor <- convert_res_for_boxplot(results_quantiseq_tumor)
quantiseq_reorg_df_normal <- convert_res_for_boxplot(results_quantiseq_normal)
cibersort_abs_reorg_df_tumor <- convert_res_for_boxplot(results_cibersort_abs_tumor)
cibersort_abs_reorg_df_normal <- convert_res_for_boxplot(results_cibersort_abs_normal)


################################################################################
### MERGE PATIENTS WITH NUMEROUS TUMOR OR NORMAL SAMPLES
################################################################################
#' Takes patients that have multiple samples and averages the absolute cell fraction
#' for each of the cell types together
#' @param results a cancer or normal results file from one of the programs
average_mult_samples <- function(results) {
  unique_patients <- unique(unlist(lapply(results$Sample, function(x) 
    unlist(strsplit(x, "-", fixed = TRUE))[3])))
  
  new_subtables <- lapply(unique_patients, function(x) {
    # Get the various samples for that patient
    patient_subset_res <- results[grepl(x, results$Sample),]
    unique_samples <- unique(patient_subset_res$Sample)
    
    # If there is more than one unique sample name, we need to do the averaging
    if (length(unique_samples) > 1) {
      # Take the information for the first sample
      merged_subset_res <- patient_subset_res[patient_subset_res$Sample == unique_samples[1],]
      cell_types <- merged_subset_res$Cell.Type
      # Incorporate the other samples into that first sample
      for (i in 2:length(unique_samples)) {
        samplei_subset_res <- patient_subset_res[patient_subset_res$Sample == unique_samples[i],]
        average_fract_per_ct <- unlist(lapply(cell_types, function(ct) {
          prev_mean <- merged_subset_res[merged_subset_res$Cell.Type == ct, 'Absolute.Cell.Fraction']
          new_val <- samplei_subset_res[samplei_subset_res$Cell.Type == ct, 'Absolute.Cell.Fraction']
          return(mean(c(prev_mean, new_val)))
        }))
        merged_subset_res$Absolute.Cell.Fraction <- average_fract_per_ct
      }
      return(merged_subset_res)
      
      # Otherwise, just return the table with that one sample
    } else {
      return(patient_subset_res)
    }
  })
  
  # Bind together the new subtables
  new_results_df <- do.call(rbind, new_subtables)
  return(new_results_df)
}

# Call this function
epic_reorg_df_tumor <- average_mult_samples(epic_reorg_df_tumor)
epic_reorg_df_normal <- average_mult_samples(epic_reorg_df_normal)
quantiseq_reorg_df_tumor <- average_mult_samples(quantiseq_reorg_df_tumor)
quantiseq_reorg_df_normal <- average_mult_samples(quantiseq_reorg_df_normal)
cibersort_abs_reorg_df_tumor <- average_mult_samples(cibersort_abs_reorg_df_tumor)
cibersort_abs_reorg_df_normal <- average_mult_samples(cibersort_abs_reorg_df_normal)


################################################################################
### VISUALIZE ABSOLUTE CELL FRACTION RESULTS
################################################################################
#' Plot a boxplot of the absolute cell fraction output
#' Methods that produce absolute cell fractions: EPIC, quanTIseq, and CIBERSORT abs
#' @param reorg_df a results data frame from an immune cell deconvolution method
#' that produces absolute cell fractions, reorganized to be suitable for a box plot
#' @param name_method the name of the immune cell deconvolution method that produced
#' the given results data frame
plot_boxplot <- function(reorg_df, name_method) {
  ggplot(reorg_df, aes(x=Cell.Type, y=Absolute.Cell.Fraction)) + geom_boxplot(outlier.colour="black", outlier.shape=16,
                                                                              outlier.size=2, notch=FALSE, fill="cornflowerblue") +
    xlab("Cell Type") + ylab("Absolute Cell Fraction, Across Samples") + ggtitle(paste(name_method, "Immune Cell Deconvolution Results", sep = " ")) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
}
plot_boxplot(epic_reorg_df_tumor, "EPIC Tumor")
plot_boxplot(epic_reorg_df_normal, "EPIC Normal")
plot_boxplot(quantiseq_reorg_df_tumor, "quanTIseq Tumor")
plot_boxplot(quantiseq_reorg_df_normal, "quanTIseq Normal")
plot_boxplot(cibersort_abs_reorg_df_tumor, "cibersort_abs Tumor")
plot_boxplot(cibersort_abs_reorg_df_normal, "cibersort_abs Normal")


# OPT: Make a stacked barplot to show variation between samples
ggplot(cibersort_abs_reorg_df, aes(fill = Cell.Type, y = Absolute.Cell.Fraction, x = Sample)) + 
  geom_bar(position = 'stack', stat = 'identity')


################################################################################
### VISUALIZE RELATIVE CELL FRACTION RESULTS
################################################################################
#' Plot a boxplot of the relative cell fraction output
#' Methods that produce relative fractions: xCell, MCP-Counter, TIMER, and CIBERSORT
#' @param reorg_df a results data frame from an immune cell deconvolution method
#' that produces relative cell fractions, reorganized to be suitable for a box plot
#' @param name_method the name of the immune cell deconvolution method that produced
#' the given results data frame
plot_boxplot_relative <- function(reorg_df, name_method) {
  ggplot(reorg_df, aes(x=Cell.Type, y=Absolute.Cell.Fraction)) + geom_boxplot(outlier.colour="black", outlier.shape=16,
                                                                              outlier.size=2, notch=FALSE, fill="cornflowerblue") +
    xlab("Cell Type") + ylab("Relative Cell Fraction, Across BRCA Samples") + ggtitle(paste(name_method, "Immune Cell Deconvolution Results, BRCA", sep = " ")) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
}
plot_boxplot_relative(xcell_reorg_df[xcell_reorg_df$Cell.Type %in% c("microenvironment score", "stroma score", "immune score"),], "xCell")
plot_boxplot_relative(mcpcounter_reorg_df, "MCP-Counter")
plot_boxplot_relative(timer_reorg_df, "TIMER")
plot_boxplot_relative(cibersort_reorg_df, "CIBERSORT")

#' Plot a dotplot using a results data frame from an immune cell deconvolution method
#' @param results the raw results file from the immune cell deconvolution method,
#' not reorganized in any way
plot_dotplots <- function(results) {
  results %>% gather(sample, score, -cell_type) %>%
    ggplot(aes(x=sample, y=score, color=cell_type)) +
    geom_point(size=1) +
    facet_wrap(~cell_type, scales="free_x", ncol=3) +
    #scale_color_brewer(palette="Paired", guide=FALSE) +
    coord_flip() + theme_bw() + theme(axis.text.y = element_blank())
}
plot_dotplots(results_mcpcounter)
plot_dotplots(results_xCell)
plot_dotplots(results_timer)
plot_dotplots(results_cibersort)


################################################################################
### VISUALIZE PREIDICTED TUMOR VS. NORMAL LEVEL OF IMMUNE CELL INFILTRATION
################################################################################
#' Make a dot plot that shows the estimated fraction of particular immune cell 
#' types in tumor vs. normal samples for a given absolute immune cell 
#' deconvolution method (EPIC, CIBERSORT abs, OR quanTIseq)
#' @param results_cancer a cancer results table from the given method
#' @param results_normal a normal results table from the given method
#' @param method a string with the name of the method for labeling ('EPIC', 
#' 'quanTIseq', or 'CIBERSORT')
plot_tumor_vs_norm_infilt_estimates <- function(results_cancer, results_normal, method) {

  # Adjust the sample names so they are only the patient IDs (and we can pair them)
  results_cancer$Sample <- unlist(lapply(results_cancer$Sample, function(x) 
    unlist(strsplit(x, "-", fixed = TRUE))[3]))
  results_normal$Sample <- unlist(lapply(results_normal$Sample, function(x) 
    unlist(strsplit(x, "-", fixed = TRUE))[3]))
  
  # Add a tumor and normal column (optional) to make the table more clear
  #results_cancer$Tum.Or.Norm <- rep("Tumor", nrow(results_cancer))
  #results_normal$Tum.Or.Norm <- rep("Normal", nrow(results_normal))
  
  # Merge the cancer and normal results files
  results_merged <- merge(results_cancer, results_normal, by = c('Sample', 'Cell.Type')) 
  
  print(results_merged)

  # Plot the tumor vs. normal cell fraction for each cell type
  ggplot(results_merged, aes(x=Absolute.Cell.Fraction.x, y = Absolute.Cell.Fraction.y)) + 
    geom_point() + labs(x = "Cell Fraction (Tumor)", y = "Cell Fraction (Normal)") + 
    geom_abline(slope=1, intercept=0) + facet_wrap(~Cell.Type) + ggtitle(method)
  
  # Plot the tumor vs. normal fraction for immune cells as a whole
  #cell_fractions <- lapply(unique(results_merged$Sample), function(x) {
    #if(method == "CIBERSORT") {
      # CIBERSORT doesn't have an 'uncharacterized' cell type category (just sum everything)
      #samp_data <- results_merged[results_merged$Sample == x, c('Absolute.Cell.Fraction.x', 
                                                              #'Absolute.Cell.Fraction.y')]
    #} else {
      #samp_data <- results_merged[((results_merged$Sample == x) & 
                                   #!(results_merged$Cell.Type == "uncharacterized cell")), 
                                #c('Absolute.Cell.Fraction.x', 'Absolute.Cell.Fraction.y')]
    #}
    #return(data.frame('Sample' = x, 'Absolute.Cell.Fraction.x' = sum(samp_data$Absolute.Cell.Fraction.x, na.rm = TRUE),
                      #'Absolute.Cell.Fraction.y' = sum(samp_data$Absolute.Cell.Fraction.y, na.rm = TRUE)))
  #})
  #cell_frac_table <- do.call(rbind, cell_fractions)
  #print(cell_frac_table)
  
  #ggplot(cell_frac_table, aes(x=Absolute.Cell.Fraction.x, y = Absolute.Cell.Fraction.y)) + 
    #geom_point() + labs(x = "Cell Fraction (Tumor)", y = "Cell Fraction (Normal)") + 
    #geom_abline(slope=1, intercept=0) + ggtitle(method)
}



plot_tumor_vs_norm_infilt_estimates(epic_reorg_df_tumor, epic_reorg_df_normal, 
                                    "EPIC")
plot_tumor_vs_norm_infilt_estimates(quantiseq_reorg_df_tumor, quantiseq_reorg_df_normal, 
                                    "quanTIseq")
plot_tumor_vs_norm_infilt_estimates(cibersort_abs_reorg_df_tumor, cibersort_abs_reorg_df_normal, 
                                    "CIBERSORT Abs")

################################################################################
### OPTIONAL CODE FOR RUNNING EPIC USING ITS OWN R PACKAGE
################################################################################
# Use EPIC to get tumor cell deconvolution results
epic_results_tpm <- EPIC(bulk = tpm_dataframe_gn)
epic_results_fpkm <- EPIC(bulk = fpkm_dataframe_gn)

# Take the pieces from these results
epic_mRNAprop_tpm <- as.data.frame(epic_results_tpm[[1]])
epic_cellFrac_tpm <- as.data.frame(epic_results_tpm[[2]])
epic_gof_tpm <- epic_results_tpm[[3]]

epic_mRNAprop_fpkm <- as.data.frame(epic_results_fpkm[[1]])
epic_cellFrac_fpkm <- as.data.frame(epic_results_fpkm[[2]])
epic_gof_fpkm <- epic_results_fpkm[[3]]

# Write these results to files
write.csv(epic_mRNAprop_tpm, paste(output_path, "EPIC Cell Deconvolution/epic_mRNAprop_tpm.csv", sep = ""))
write.csv(epic_cellFrac_tpm, paste(output_path, "EPIC Cell Deconvolution/epic_cellFrac_tpm.csv", sep = ""))
write.csv(epic_gof_tpm, paste(output_path, "EPIC Cell Deconvolution/epic_gof_tpm.csv", sep = ""))

write.csv(epic_mRNAprop_fpkm, paste(output_path, "EPIC Cell Deconvolution/epic_mRNAprop_fpkm.csv", sep = ""))
write.csv(epic_cellFrac_fpkm, paste(output_path, "EPIC Cell Deconvolution/epic_cellFrac_fpkm.csv", sep = ""))
write.csv(epic_gof_fpkm, paste(output_path, "EPIC Cell Deconvolution/epic_gof_fpkm.csv", sep = ""))







