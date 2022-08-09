############################################################
### Additional Pre-model Inputs & Processing 
### Written By: Sara Geraghty, April 2021
############################################################

library(dplyr)
library(VennDiagram)
library(data.table)
library(TCGAbiolinks)


# This file contains additional pre-processing that is applied to all the data file
# types that will be input to the linear model function.
# This includes:
  # 1. Combine patient and sample DFs into one
  # 2. Limit methylation (non-tumor-normal matched) to only cancer files
  # 3. For non-tumor-normal matched data, limiting to only overlapping samples
  # 4. Creating an regulatory protein data frame for each level of specificity (I-Protein, I-Domain,
      # and I-Binding Position) in the form of matched Uniprot and ENSG ID inputs
  # 5. Construct a data frame of targets for a given regulatory protein of interest (for model
      # testing purposes)
  # 6. Find samples that have both an amplification and mutation, or deletion and mutation, in the same gene

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/cBioPortal/METABRIC/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)

############################################################
############################################################
# IMPORT NECESSARY FILES FOR MODEL (EXPRESSION, METHYLATION,
# TARGET GENE MUTATION, REGULATORY PROTEIN MUTATION, CNA, &
# PATIENT INFORMATION)
############################################################
############################################################

############################################################
# IMPORT EXPRESSION FILES
############################################################
### TUMOR-NORMAL MATCHED ###

# 1. Filtered FPKM

# 2. Filtered TMM

# 3. Quantile-Normalized

# 4. Rank-Normalized

# TODO: Run through the TN-Matched pipeline with all of these expression datasets
# and then import each variation of the file
# (see "limit_to_tumor_normal_matched.R")


### NON-TUMOR-NORMAL MATCHED ###
# 1. Filtered FPKM 
expression_df_fpkm <- fread(paste(main_path, "Expression/expression_fpkm_filt_gn_DF.csv", sep = ""), 
                          header = TRUE)  # fpkm-normalized and filtered
colnames(expression_df_fpkm)[1] <- 'ensg_id'

# 2. Filtered TMM 
expression_df_tmm <- fread(paste(main_path, "Expression/tmm_normalized_expression_counts.csv", sep = ""), 
                          header = TRUE)
# expression_df_tmm <- fread(paste(main_path, "Expression/ALL_CT_tmm_normalized_expression_filtByExpr_counts.csv", sep = ""), 
                        #header = TRUE)
colnames(expression_df_tmm)[1] <- 'ensg_id'

# 3. Quantile-Normalized
expression_df_qn <- fread(paste(main_path, "Expression/expression_quantile_norm_DF.csv", sep = ""),
                          header = TRUE)
colnames(expression_df_qn)[1] <- 'ensg_id'

# 4. Rank-Normalized
expression_df_rn <- fread(paste(main_path, "Expression/expression_rank_norm_DF.csv", sep = ""),
                          header = TRUE)
colnames(expression_df_rn)[1] <- 'ensg_id'


############################################################
# IMPORT METHYLATION FILES
############################################################
### TUMOR-NORMAL MATCHED ###

#TODO: run through TN-Matched pipeline for methylation files and import each variation
# (see "limit_to_tumor_normal_matched.R")


### NON-TUMOR-NORMAL MATCHED ###
# Raw
methylation_df_Beta <- fread(paste(main_path, "Methylation/methylation_DF_Beta_cancer_only.csv", sep = ""), 
                        header = TRUE)
methylation_df_M <- fread(paste(main_path, "Methylation/methylation_DF_M_cancer_only.csv", sep = ""), 
                        header = TRUE)

# Thresholded
methylation_df_Beta_0.8 <- fread(paste(main_path, "Methylation/methylation_DF_Beta_0.8.csv", sep = ""), 
                           header = TRUE)
colnames(methylation_df_Beta_0.8)[1] <- 'Gene_Symbol'
methylation_df_Beta_0.5 <- fread(paste(main_path, "Methylation/methylation_DF_Beta_0.5.csv", sep = ""), 
                        header = TRUE)
colnames(methylation_df_Beta_0.5)[1] <- 'Gene_Symbol'


# Bucketed
methylation_df_bucketed_Beta <- fread(paste(main_path, "Methylation/methylation_DF_bucketed_Beta.csv", sep = ""), 
                        header = TRUE)
methylation_df_bucketed_Beta$Gene_Symbol <- methylation_df_Beta$Gene_Symbol # If needed

methylation_df_bucketed_M <- fread(paste(main_path, "Methylation/methylation_DF_bucketed_M.csv", sep = ""), 
                        header = TRUE)
methylation_df_bucketed_M$Gene_Symbol <- methylation_df_M$Gene_Symbol # If needed


############################################################
# IMPORT GENE TARGET MUTATION FILES
############################################################
### TUMOR-NORMAL MATCHED ###
#TODO: run through TN-Matched pipeline for mutation files and import each variation
# (see "limit_to_tumor_normal_matched.R")

# If Swissprot IDs are not already added:
#mutation_targ_df$Swissprot <- unlist(lapply(rownames(mutation_targ_df), function(x) 
# paste(unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'uniprot_gn_id']), collapse = ";")))

### NON-TUMOR-NORMAL MATCHED ###
mutation_targ_df <- fread(paste(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense.csv", sep = ""), 
                             header = TRUE)
#mutation_targ_df <- fread(paste(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense_ALL.csv", sep = ""), 
                          #header = TRUE)
colnames(mutation_targ_df)[1] <- 'Gene_Symbol'
colnames(mutation_targ_df)[2:ncol(mutation_targ_df)] <- unlist(lapply(colnames(mutation_targ_df)[2:ncol(mutation_targ_df)], function(x) 
  paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-"))) # Fix the sample IDs

# Add the additional un-included genes (those with no mutations in any patient)
allgene_targets <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/allgene_targets.csv",
                            header = TRUE, check.names = FALSE, row.names = 1)
all_gene_names <- unlist(lapply(allgene_targets$ensg, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == x,'external_gene_name'])), collapse = ";")))
all_gene_names <- unique(all_gene_names)
unincluded_genes <- setdiff(all_gene_names, mutation_targ_df$Gene_Symbol)
unincluded_gene_df <- data.frame(matrix(nrow = length(unincluded_genes), ncol = ncol(mutation_targ_df)))
unincluded_gene_df[is.na(unincluded_gene_df)] <- 0
colnames(unincluded_gene_df) <- colnames(mutation_targ_df)
unincluded_gene_df$Gene_Symbol <- unincluded_genes

mutation_targ_df <- rbind(mutation_targ_df, unincluded_gene_df)

# If Swissprot IDs are not already added:
#mutation_targ_df$Swissprot <- unlist(lapply(mutation_targ_df$Gene_Symbol, function(x) 
# paste(unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'uniprot_gn_id']), collapse = ";")))



############################################################
# IMPORT REGULATORY PROTEIN MUTATION FILES
############################################################
### TUMOR-NORMAL MATCHED ###
#TODO: run through TN-Matched pipeline for mutation files and import each variation
# (see "limit_to_tumor_normal_matched.R")


### NON-TUMOR-NORMAL MATCHED ###
mutation_regprot_df <- fread(paste(main_path, "Mutation/iprotein_results_nonsynonymous.csv", sep = ""),
                             header = TRUE)   # I-Protein, all ligands
mutation_regprot_df <- fread(paste(main_path, "Mutation/iprotein_results_missense_nonsense.csv", sep = ""),
                                header = TRUE)   # I-Protein, all ligands
mutation_regprot_df <- fread(paste(main_path, "Mutation/iprotein_results_missense_nucacids.csv", sep = ""),
                                header = TRUE)   # I-Protein, nucleic acids only

# For I-Protein only, add a Swissprot ID column
mutation_regprot_df$Swissprot <- unlist(lapply(mutation_regprot_df$Query, function(x) 
  unlist(strsplit(x, "|", fixed = TRUE))[2]))


mutation_regprot_df <- fread(paste(main_path, "Mutation/idomain_results_missense.csv", sep = ""),
                                header = TRUE)   # I-Domain, all ligands
mutation_regprot_df <- fread(paste(main_path, "Mutation/idomain_results_missense_nucacids.csv", sep = ""),
                                header = TRUE)   # I-Domain, nucleic acids only
mutation_regprot_df <- fread(paste(main_path, "Mutation/ibindingpos_results_missense.csv", sep = ""),
                                header = TRUE)   # I-Binding position, all ligands
mutation_regprot_df <- fread(paste(main_path, "Mutation/ibindingpos_results_missense_nucacids.csv", sep = ""),
                                header = TRUE)   # I-Binding position, nucleic acids only


############################################################
# IMPORT CNA FILES
############################################################
### TUMOR-NORMAL MATCHED ###
#TODO: run through TN-Matched pipeline for CNA files and import each variation
# (see "limit_to_tumor_normal_matched.R")


### NON-TUMOR-NORMAL MATCHED ###
# Gene-Level Score (GISTIC2 Data)
# cna_df <- fread(paste(main_path, "CNV/Gene-level Score/CNV_DF_AllGenes.csv", sep = ""), 
      # header = TRUE, row.names = 1)

# Raw Copy Number 
cna_df_raw <- read.csv(paste(main_path, "CNV/Gene-level Raw/CNV_DF_AllGenes.csv", sep = ""), 
                       header = TRUE, row.names = 1, check.names = FALSE)
# Bucketed Copy Number (Incl. Amp)
cna_df_bucket_inclAmp <- read.table(paste(main_path, "CNV/Gene-level Raw/CNV_DF_bucketed_inclAmp_AllGenes.tsv", sep = ""), 
                header = TRUE, row.names = 1, check.names = FALSE)
# Bucketed Copy Number (Excl. Amp)
cna_df_bucket_exclAmp <- read.table(paste(main_path, "CNV/Gene-level Raw/CNV_DF_bucketed_exclAmp_AllGenes.tsv", sep = ""), 
                header = TRUE, row.names = 1, check.names = FALSE)

#' Adjust the gene names to remove the part after the '.'
#' @param cna_df one of the above CNA DFs
adj_rownams <- function(cna_df) {
  rownames(cna_df) <- unlist(lapply(rownames(cna_df), function(x) 
    unlist(strsplit(x, ".", fixed = TRUE))[1]))
  return(cna_df)
}
cna_df_raw <- adj_rownams(cna_df_raw)
cna_df_bucket_inclAmp <- adj_rownams(cna_df_bucket_inclAmp)
cna_df_bucket_exclAmp <- adj_rownams(cna_df_bucket_exclAmp)


############################################################
# IMPORT PATIENT FILES
############################################################
### TUMOR-NORMAL MATCHED ###
#TODO: run through TN-Matched pipeline for patient files and import each variation
# (see "limit_to_tumor_normal_matched.R")
patient_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/patient_dataframe_tnm.csv", sep = ""), 
                       header = TRUE, row.names = 1)

### NON-TUMOR-NORMAL MATCHED ###
#patient_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/patient_dataframe_ntnm.csv", sep = ""), 
                       #header = TRUE, row.names = 1)
patient_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/patient_dataframe_ntnm_MN_normAge.csv", sep = ""), 
                       header = TRUE, row.names = 1)

patient_df_brca_blca_hnsc <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/patient_dataframe_ntnm_MN_normAge_inclBRCAsubtypes.csv", sep = ""), 
                                      header = TRUE, row.names = 1)

# If BRCA, eliminate the gender column
is_brca <- TRUE
#is_brca <- FALSE
if(is_brca) {patient_df <- patient_df[, -which(colnames(patient_df) == "Gender")]}


############################################################
# IMPORT SAMPLE FILES
############################################################
### TUMOR-NORMAL MATCHED ###
#TODO: run through TN-Matched pipeline for sample files and import each variation
# (see "limit_to_tumor_normal_matched.R")

### NON-TUMOR-NORMAL MATCHED ###

# Cibersort Abs, Total Fraction, Bucketed
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_total_frac_tmm.csv", sep = ""), 
                       header = TRUE, row.names = 1)  # Total fraction
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_total_frac_fpkm.csv", sep = ""), 
                      header = TRUE, row.names = 1)  # Total fraction
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_total_frac_fpkm_top10k.csv", sep = ""), 
                      header = TRUE, row.names = 1)  # Total fraction
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_total_frac_qn_top10k.csv", sep = ""), 
                      header = TRUE, row.names = 1)  # Total fraction
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_total_frac_rn_mean_MedGr10.csv", sep = ""), 
                      header = TRUE, row.names = 1)  # Total fraction
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_total_frac_rn_mean_MedGr10_top10k.csv", sep = ""), 
                      header = TRUE, row.names = 1)  # Total fraction

# Cibersort Abs
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_abs_tmm.csv", sep = ""), 
                      header = TRUE, row.names = 1)   # Per-Cell-Type, CIBERSORT Abs
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_abs_fpkm.csv", sep = ""), 
                      header = TRUE, row.names = 1)   # Per-Cell-Type, CIBERSORT Abs
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_abs_fpkm_top10k.csv", sep = ""), 
                      header = TRUE, row.names = 1)   # Per-Cell-Type, CIBERSORT Abs
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_abs_qn_top10k.csv", sep = ""), 
                      header = TRUE, row.names = 1)   # Per-Cell-Type, CIBERSORT Abs
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_abs_rn_mean_MedGr10.csv", sep = ""), 
                      header = TRUE, row.names = 1)   # Per-Cell-Type, CIBERSORT Abs
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_abs_rn_mean_MedGr10_top10k.csv", sep = ""), 
                      header = TRUE, row.names = 1)   # Per-Cell-Type, CIBERSORT Abs

sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_timer_tmm.csv", sep = ""), 
                      header = TRUE, row.names = 1)   # Per-Cell-Type, TIMER
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_timer_fpkm.csv", sep = ""), 
                      header = TRUE, row.names = 1)   # Per-Cell-Type, TIMER
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_timer_fpkm_top10k.csv", sep = ""), 
                      header = TRUE, row.names = 1)   # Per-Cell-Type, TIMER
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_timer_qn_top10k.csv", sep = ""), 
                      header = TRUE, row.names = 1)   # Per-Cell-Type, TIMER
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_timer_rn_mean_MedGr10.csv", sep = ""), 
                      header = TRUE, row.names = 1)   # Per-Cell-Type, TIMER
sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_timer_rn_mean_MedGr10_top10k.csv", sep = ""), 
                      header = TRUE, row.names = 1)   # Per-Cell-Type, TIMER


############################################################
############################################################
# 1. COMBINE PATIENT AND SAMPLE FILES
############################################################
############################################################
#' Combines the patient and sample data frames such that there
#' is one row per patient sample. Multiple samples from the 
#' same patient will have duplicate patient information.
#' @param patient_df the patient data frame (patients are rows,
#' patient characteristics are columns)
#' @param sample_df the sample data frame (patient samples are
#' rows, sample characteristics are columns)
#' @param dataset either 'tcga', 'metabric', 'icgc', or 'cptac3'
combine_patient_and_samp_dfs <- function(patient_df, sample_df, dataset) {
  # For each sample, retrieve the corresponding patient info
  sample_dfs <- lapply(1:nrow(sample_df), function(i) {
    # Get the sample barcode and its row
    sample_barcode <- rownames(sample_df)[i]
    sample_row <- sample_df[i,]
    print(sample_barcode)
    #print(sample_row)
    
    # Get the corresponding row from the patient DF, if tcga
    patient <- sample_barcode
    if(dataset == "tcga") {
      patient <- unlist(strsplit(sample_barcode, "-", fixed = TRUE))[3]
    }
    print(patient)
    
    if(patient %fin% rownames(patient_df)) {
      patient_row <- patient_df[rownames(patient_df) == patient,]
      #print(patient_row)
      
      # Combine these together, so that the sample ID is maintained as the row name
      if(!((length(sample_row) == 0) | (length(patient_row) == 0))) {
        return(cbind(sample_row, patient_row))
      } else {
        print("Something went wrong.")
        return(NA)
      }
    } else {
      print("Patient from sample DF is not in patient DF.")
      return(NA)
    }
  })
  
  # Recombine all these new, combined rows into a full DF
  #sample_dfs <- sample_dfs[!is.na(sample_dfs)]
  combined_df <- do.call(rbind, sample_dfs)
  rownames(combined_df) <- rownames(sample_df)
  
  return(combined_df)
}

# Call this function
combined_pat_samp_df <- combine_patient_and_samp_dfs(patient_df, sample_df, 'tcga')

# Remove rows that are entirely NA
combined_pat_samp_df <- combined_pat_samp_df[rowSums(is.na(combined_pat_samp_df)) != ncol(combined_pat_samp_df),]
  

# Write the results to a file
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_tmm_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_fpkm_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_fpkm_top10k_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_qn_top10k_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_rn_mean_MedGr10_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_rn_mean_MedGr10_top10k_ntnm.csv", sep = ""))

write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_abs_tmm_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_abs_fpkm_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_abs_fpkm_top10k_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_abs_qn_top10k_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_abs_rn_mean_MedGr10_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_abs_rn_mean_MedGr10_top10k_ntnm.csv", sep = ""))

write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_timer_tmm_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_timer_fpkm_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_timer_fpkm_top10k_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_timer_qn_top10k_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_timer_rn_mean_MedGr10_ntnm.csv", sep = ""))
write.csv(combined_pat_samp_df, paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_timer_rn_mean_MedGr10_top10k_ntnm.csv", sep = ""))


#
#
#

# Read this back (generic)
patient_sample_df <- fread(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF.csv", sep = ""),
                              header = TRUE)

############################################################
############################################################
# 2. LIMIT METHYLATION FILE TO ONLY CANCER SAMPLES
############################################################
############################################################
# Limit to only cancer samples
methylation_df <- methylation_df[, c(TRUE, grepl("-0", colnames(methylation_df)[2:ncol(methylation_df)])),
                                 with = FALSE]
#colnames(methylation_df) <- unlist(lapply(colnames(methylation_df)[2:ncol(methylation_df)], 
                                          #function(x) unlist(strsplit(x, "-", fixed = TRUE))[1]))

# Save this file again for future use
write.csv(methylation_df, paste(main_path, "Methylation/methylation_DF_Beta_cancer_only.csv", sep = ""))
write.csv(methylation_df, paste(main_path, "Methylation/methylation_DF_M_cancer_only.csv", sep = ""))
write.csv(methylation_df, paste(main_path, "Methylation/methylation_DF_Beta_0.8_cancer_only.csv", sep = ""))
write.csv(methylation_df, paste(main_path, "Methylation/methylation_DF_bucketed_Beta_cancer_only.csv", sep = ""))
write.csv(methylation_df, paste(main_path, "Methylation/methylation_DF_bucketed_M_cancer_only.csv", sep = ""))


############################################################
############################################################
# 3. FOR NON-TUMOR-NORMAL MATCHED, LIMIT TO OVERLAPPING PATIENTS
############################################################
############################################################
# Get intersecting patients

#' Given a TCGA sample ID (XXXX-XXX, e.g. A0WY-01A), splits it 
#' apart and returns a list of unique patient IDs (XXXX, e.g. A0WY)
#' that have at least one cancer sample
#' @param sample_ids a vector of TCGA sample IDs (XXXX-XXX, e.g. A0WY-01A)
get_unique_patients <- function(sample_ids) {
  unique_patients <- unique(unlist(lapply(sample_ids, function(x) {
    spl_name <- unlist(strsplit(x, "-", fixed = TRUE))
    if(startsWith(spl_name[2], "0")) {return(spl_name[1])}
  })))
  return(unique_patients)
}

# Call this function for all representative data frames (check that other CNA, expression,
# etc. data frames have the same patients)
cna_patients <- get_unique_patients(colnames(cna_df_raw))    # 732; 6123
methylation_patients <- get_unique_patients(colnames(methylation_df)[2:ncol(methylation_df)])  #-1]) # 734; 6411
exp_samp_ids <- unlist(lapply(colnames(expression_df)[2:ncol(expression_df)], function(x)
  paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))
expression_patients <- get_unique_patients(exp_samp_ids)  # 738; 7424
mutation_targ_patients <- unlist(get_unique_patients(colnames(mutation_targ_df)[2:ncol(mutation_targ_df)])) # 978; 10002
mutation_regprot_patients <- unique(unlist(lapply(mutation_regprot_df$Patient, function(x) 
  unlist(strsplit(x, ";", fixed = TRUE)))))   # 619; 6732 -- ignore this (6875 with nonsense)

# If needed
#patient_sample_df <- as.data.frame(patient_sample_df)
#rownames(patient_sample_df) <- unlist(patient_sample_df[,1]) 
#patient_sample_df <- patient_sample_df[,2:ncol(patient_sample_df)]

clin_samp_ids <- unlist(lapply(rownames(patient_sample_df), function(x)
  paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))
clinical_patients <- get_unique_patients(clin_samp_ids)  # 739; 7745

intersecting_patients <- intersect(cna_patients, 
                                   intersect(methylation_patients,
                                             intersect(expression_patients,
                                                       intersect(mutation_targ_patients, clinical_patients))))
print(length(intersecting_patients)) # 664 in BRCA; 4545 PC

#' Given a list of intersecting patient IDs (XXXX, e.g. A0WY), it subsets a 
#' given data frame with column names that contain sample IDs (XXXX-XXX, e.g. A0WY-01A)
#' to include only those given patients
#' @param intersecting_patients a vector of intersecting patient IDs
#' @param df the data frame to subset using the patient ids
#' @param colNames a TRUE/FALSE value indicating whether the patient/sample IDs are
#' found on the rows or the columns of the DF (TRUE is columns, FALSE is rows)
subset_by_intersecting_ids <- function(intersecting_patients, df, colNames) {
  if(colNames == TRUE) {
    # Keep only intersecting patients
    just_patients <- unlist(lapply(colnames(df)[2:ncol(df)], function(x)
      unlist(strsplit(x, "-", fixed = TRUE))[1]))
    cols_to_keep <- as.numeric(c(1,(which(just_patients %fin% intersecting_patients)+1)))
    print(length(cols_to_keep))
    df <- as.data.frame(df)
    df_adj <- df[, cols_to_keep]
    

    # Keep only cancer samples
    just_samples <- unlist(lapply(colnames(df_adj)[2:ncol(df_adj)], function(x)
      unlist(strsplit(x, "-", fixed = TRUE))[2]))
    df_adj <- df_adj[,c(1,which(unlist(lapply(just_samples, function(x) 
      startsWith(x, "0")))))]
    # Remove any duplicate samples
    df_adj <- df_adj[, !(grepl("-1", colnames(df_adj), fixed = TRUE))]

  } else {
    if(length(unlist(strsplit(rownames(df)[1], "-", fixed = TRUE))) == 4) {
      rownames(df) <- unlist(lapply(rownames(df), function(x) paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], 
                                                                    collapse = "-")))
    }
    # Keep only intersecting patients
    just_patients <- unlist(lapply(rownames(df), function(x)
      unlist(strsplit(x, "-", fixed = TRUE))[1]))
    df_adj <- df[which(just_patients %fin% intersecting_patients),]
    print(head(df_adj))
    # Keep only cancer samples
    just_samples <- unlist(lapply(rownames(df_adj), function(x)
      unlist(strsplit(x, "-", fixed = TRUE))[2]))
    df_adj <- df_adj[which(unlist(lapply(just_samples, function(x) 
      startsWith(x, "0")))),]
    # Remove any duplicate samples
    df_adj <- df_adj[!(grepl('-1', rownames(df_adj), fixed = TRUE)),]
  }
  return(df_adj)
}

# Call this function
cna_df_sub <- subset_by_intersecting_ids(intersecting_patients, cna_df, TRUE)
methylation_df_sub <- subset_by_intersecting_ids(intersecting_patients, methylation_df, TRUE)
mutation_targ_df_sub <- subset_by_intersecting_ids(intersecting_patients, mutation_targ_df, TRUE)
colnames(expression_df)[2:ncol(expression_df)] <- unlist(lapply(colnames(expression_df)[2:ncol(expression_df)], function(x) 
  paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse= "-")))
expression_df_sub <- subset_by_intersecting_ids(intersecting_patients, expression_df, TRUE)
rownames(patient_sample_df) <- unlist(lapply(rownames(patient_sample_df), function(x) 
  paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse= "-")))
patient_sample_df_sub <- subset_by_intersecting_ids(intersecting_patients, patient_sample_df, FALSE)


# Handle the regulatory protein DF separately
regprot_df_new_patient_labels <- lapply(mutation_regprot_df$Patient, function(x) {
  # Split apart the semicolon separated sample IDs
  spl_patients <- unlist(strsplit(x, ";", fixed = TRUE))
  # Extract just the 4-digit patient ID
  spl_patients_nosamp <- unlist(lapply(spl_patients, function(x) 
    unlist(strsplit(x, "-", fixed = TRUE))[1]))
  # Check if this patient is in the intersecting patients list and return TRUE if so, F o.w.
  matching_patient_ind <- unlist(lapply(spl_patients_nosamp, function(y) {
    if(y %fin% intersecting_patients) {
      return(TRUE)
    } else {return(FALSE)}
  }))
  # If there are overlapping patients in this row, then we want to keep this row
  if(TRUE %in% matching_patient_ind) {
    samps_to_keep <- spl_patients[matching_patient_ind]
    
    # Keep only cancer samples
    samp_ids <- unlist(lapply(samps_to_keep, function(x) unlist(strsplit(x, "-", fixed = TRUE))[2]))
    samps_to_keep <- samps_to_keep[unlist(lapply(samp_ids, function(x) startsWith(x, "0")))]

    # Keep only non-duplicate samples
    samps_to_keep <- samps_to_keep[!grepl(".1", samps_to_keep, fixed = TRUE)]

    # If we still have samples, return them in this row 
    if(length(samps_to_keep) > 0) {
      return(paste(samps_to_keep, collapse = ";"))
    } else {return(NA)}
  } else {return(NA)}
})

# Remove NA
mutation_regprot_df_sub <- mutation_regprot_df[unlist(lapply(regprot_df_new_patient_labels, function(x) 
  ifelse(is.na(x), FALSE, TRUE))),]
# Update the DF with the new intersecting patient labels
mutation_regprot_df_sub$Patient <- regprot_df_new_patient_labels[!is.na(regprot_df_new_patient_labels)]
mutation_regprot_df_sub$Patient <- unlist(lapply(mutation_regprot_df_sub$Patient, function(x) return(unlist(x))))

# Add ENSG IDs to methylation DF
methylation_df_sub$ensg_ids <- unlist(lapply(methylation_df_sub$Gene_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'ensembl_gene_id'])), collapse = ";")))

# Re-write these all to to files
write.csv(cna_df_sub, paste(main_path, "Linear Model/Tumor_Only/CNV/CNA_AllGenes_CancerOnly_IntersectPatients.csv", sep = ""))
write.csv(methylation_df_sub, paste(main_path, "Linear Model/Tumor_Only/Methylation/methylation_M_CancerOnly_IntersectPatients.csv", sep = ""))
write.csv(expression_df_sub, paste(main_path, "Linear Model/Tumor_Only/Expression/expression_tmm_IntersectPatients.csv", sep = ""))
write.csv(mutation_targ_df_sub, paste(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_IntersectPatients.csv", sep = ""))
write.csv(mutation_regprot_df_sub, paste(main_path, "Linear Model/Tumor_Only/Regprot_Mutation/iprotein_results_missense_IntersectPatients.csv", sep = ""))
write.csv(patient_sample_df_sub, paste(main_path, "Linear Model/Tumor_Only/Patient/combined_patient_sample_IntersectPatients.csv", sep = ""))

#OPT:
write.table(intersecting_patients, paste(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt", sep = ""))
intersecting_patients <- read.table(paste(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt", sep = ""))[,1]

############################################################
############################################################
# 4. MAKE REGULATORY PROTEIN INPUT DATAFRAMES
############################################################
############################################################
prot_path <- paste(main_path, "Mutation/", sep = "")

iprotein_prots <- as.character(unlist(fread(paste(prot_path, "swissprot_ids_missense_iprotein.csv", sep = ""), header = TRUE)[,2]))
iprotein_prots_nucacids <- as.character(unlist(fread(paste(prot_path, "swissprot_ids_missense_iprotein_nucacids.csv", sep = ""), header = TRUE)[,2]))

idomain_prots <- as.character(unlist(fread(paste(prot_path, "swissprot_ids_missense_idomain.csv", sep = ""), header = TRUE)[,2]))  # all ligands
idomain_prots_nucacids <- as.character(unlist(fread(paste(prot_path, "swissprot_ids_missense_idomain_nucacids.csv", sep = ""), header = TRUE)[,2]))  # DNA/RNA-binding only

ibindingpos_prots <- as.character(unlist(fread(paste(prot_path, "swissprot_ids_missense_ibindingpos.csv", sep = ""), header = TRUE)[,2]))  # all ligands
ibindingpos_prots_nucacids <- as.character(unlist(fread(paste(prot_path, "swissprot_ids_missense_ibindingpos_nucacids.csv", sep = ""), header = TRUE)[,2]))  # DNA/RNA-binding only

#' Converts protein uniprot IDs to ensembl IDs
#' @param protein_ids a vector of uniprot IDs to be converted
convert_to_ensembl <- function(protein_ids) {
  return(unlist(lapply(protein_ids, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'ensembl_gene_id']), collapse = ";"))))
}
iprotein_prots_ensg <- convert_to_ensembl(iprotein_prots)
iprotein_prots_nucacids_ensg <- convert_to_ensembl(iprotein_prots_nucacids)
idomain_prots_ensg <- convert_to_ensembl(idomain_prots)
idomain_prots_nucacids_ensg <- convert_to_ensembl(idomain_prots_nucacids)
ibindingpos_prots_ensg <- convert_to_ensembl(ibindingpos_prots)
ibindingpos_prots_nucacids_ensg <- convert_to_ensembl(ibindingpos_prots_nucacids)


#' Combines uniprot and ENSG IDs into a protein IDs data frame, which it writes to a CSV
#' @param protein_ids a vector of uniprot IDs for regulatory proteins
#' @param protein_ids_ensembl a vector of ENSG IDs for regulatory proteins
#' @param label the level of specificity in question ("iprotein", "idomain", etc.)
#' @param prot_path a path for where to write the output data frame
recombine_into_df_and_write <- function(protein_ids, protein_ids_ensembl, label, prot_path) {
  protein_ids_df <- data.frame("swissprot_ids" = protein_ids, "ensg_ids" = protein_ids_ensembl)
  # Remove any that do not have ENSG IDs
  protein_ids_df <- protein_ids_df[!protein_ids_df$ensg_ids == "",]
  # Write this as a CSV file for later use
  write.csv(protein_ids_df, paste(prot_path, paste("Files for Linear Model/", paste(label, "protein_ids_df.csv", sep = "_"), sep = ""), sep = ""))
}

recombine_into_df_and_write(iprotein_prots, iprotein_prots_ensg, "iprotein", prot_path)
recombine_into_df_and_write(iprotein_prots_nucacids, iprotein_prots_nucacids_ensg, "iprotein_nucacids", prot_path)
recombine_into_df_and_write(idomain_prots, idomain_prots_ensg, "idomain", prot_path)
recombine_into_df_and_write(idomain_prots_nucacids, idomain_prots_nucacids_ensg, "idomain_nucacids", prot_path)
recombine_into_df_and_write(ibindingpos_prots, ibindingpos_prots_ensg, "ibindingpos", prot_path)
recombine_into_df_and_write(ibindingpos_prots_nucacids, ibindingpos_prots_nucacids_ensg, "ibindingpos_nucacids", prot_path)



#' Make regulatory input data frame using a mutational threshold (e.g. >=5%), with an option
#' to limit to known driver genes
#' @param regprot_df a regulatory mutation data frame for a given specificity/ mutation type
#' @param mut_freq_thres a mutation frequency threshold (e.g. 5%) for which to include the 
#' given gene
#' @param patient_ids a set of unique patient/ sample IDs for the given cohort
#' @param driver_df if not NA, use known driver DF to limit to only proteins with known drivers
#' @param all_genes_id_conv a gene ID conversion file from bioMart
create_regulatory_prot_input_df <- function(regprot_df, mut_freq_thres, patient_ids, 
                                            driver_df, all_genes_id_conv) {
  
  # Calculate the minimum number of mutated samples needed to keep a given gene
  mut_count_thres <- mut_freq_thres * length(patient_ids)
  print(mut_count_thres)
  
  # Subset the regprot data frame using this
  regprot_df_rows <- lapply(1:nrow(regprot_df), function(i) {
    patients <- unlist(strsplit(as.character(unlist(regprot_df[i, 'Patient'])), ";", fixed = TRUE))
    # Limit to just patients in the intersecting set
    patients_justID <- unlist(lapply(patients, function(x) unlist(strsplit(x, "-", fixed = TRUE))[1]))
    patients <- patients[which(patients_justID %in% patient_ids)]
    if (length(patients) > mut_count_thres) {return(regprot_df[i,])}
    else {return(NA)}
  })
  regprot_df_rows <- regprot_df_rows[!is.na(regprot_df_rows)]
  
  # Extract the ENSG and Swissprot IDs
  uniprot_ids <- unique(unlist(lapply(regprot_df_rows, function(x) x$Swissprot)))
  ensg_ids <- unlist(lapply(uniprot_ids, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 
                                   'ensembl_gene_id']), collapse = ";")))
  
  output_df <- data.frame("swissprot_ids" = uniprot_ids, "ensg_ids" = ensg_ids)
  
  # If driver_df is not NA, then we want to limit to only driver genes
  if(!is.na(driver_df)) {
    output_df <- output_df[which(ensg_ids %fin% driver_df$ensembl_gene_id),]
  }
  
  return(output_df)
}

# Call function
regprot_input_df_all_5perc <- create_regulatory_prot_input_df(mutation_regprot_df, 0.05, 
                                                              intersecting_patients, NA,
                                                              all_genes_id_conv)
regprot_input_df_driver_5perc <- create_regulatory_prot_input_df(mutation_regprot_df, 0.05, 
                                                              intersecting_patients, driver_gene_df,
                                                              all_genes_id_conv)

# Limit just to particular subtypes
patient_set_lumA <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.A_patient_ids.txt"), header = TRUE)[,1], intersecting_patients)
patient_set_lumB <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.B_patient_ids.txt"), header = TRUE)[,1], intersecting_patients)
patient_set_basal <- intersect(read.table(paste0(main_path, "Patient Subsets/Basal_patient_ids.txt"), header = TRUE)[,1], intersecting_patients)
patient_set_her2 <- intersect(read.table(paste0(main_path, "Patient Subsets/HER2_patient_ids.txt"), header = TRUE)[,1], intersecting_patients)


regprot_input_df_all_5perc_missense_lumA <- create_regulatory_prot_input_df(mutation_regprot_df_missense, 0.05, 
                                                                            patient_set_lumA, NA,
                                                                            all_genes_id_conv)

# Write to files
write.csv(regprot_input_df_all_5perc, paste0(prot_path, "Files for Linear Model/iprotein_protein_ids_df_gr0.05Freq_drivers_missense.csv"))


############################################################
############################################################
# 5. MAKE TESTER TARGETS DATAFRAME FOR GIVEN REG. PROTEIN
############################################################
############################################################
sample_protein_uniprot <- "P04637" 
sample_protein_ensg <- "ENSG00000141510"
sample_protein_name <- "P53"


### OPTION 1: CURATED TARGETS ###
# Use the I-Protein curated targets data frame, since it is the most comprehensive
curated_targets_df <- fread(paste(main_path, "Curated TF Data/iprotein_curated_targets_df_full.csv", sep = ""), 
                               header = TRUE)
sample_targets <- curated_targets_df[,colnames(curated_targets_df) == sample_protein_name] 
sample_targets <- sample_targets[!is.na(sample_targets)] 

# For TP53, get targets from publication (Fischer, 2017, https://doi.org/10.1038/onc.2016.502)
tp53_targets_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Validation_Files/TP53_targets.csv", 
                            header = TRUE)
tp53_targets_df <- tp53_targets_df[tp53_targets_df$Gene.Symbol != "",]
sample_targets <- tp53_targets_df$Gene.Symbol


### OPTION 2: ChIP-eat file targets ###
sample_targets <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/ChIP-eat/tp53_chipeat_targets.txt")[,1]
sample_targets <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/ChIP-eat/foxa1_chipeat_targets.txt")[,1]

# For options 1 and 2: get target swissprot IDs and ENSG IDs from gene names
sample_targets_swissprot <- unlist(lapply(sample_targets, function(x) 
  paste(unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'uniprot_gn_id']), collapse = ";")))
sample_targets_swissprot <- unique(sample_targets_swissprot[sample_targets_swissprot != ""])
# Get target ENSG IDs
sample_targets_ensg <- unlist(lapply(sample_targets_swissprot, function(x) 
  paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == unlist(strsplit(x, ";", fixed = TRUE))[1], 
                                 'ensembl_gene_id']), collapse = ";")))

### OPTION 3: ALL POSSIBLE GENE TARGETS ###
all_genes_id_conv_sub <- all_genes_id_conv[,c('ensembl_gene_id', 'uniprot_gn_id')]
all_genes_id_conv_sub <- distinct(all_genes_id_conv_sub[all_genes_id_conv_sub$uniprot_gn_id != "",])
write.csv(all_genes_id_conv_sub, "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv_subset.csv")
sample_targets_ensg <- unique(all_genes_id_conv_sub$ensembl_gene_id)   # Get the unique ENSG ids
sample_targets_swissprot <- unlist(lapply(sample_targets_ensg, function(x) 
  paste(unlist(all_genes_id_conv_sub[all_genes_id_conv_sub$ensembl_gene_id == x, 'uniprot_gn_id']), collapse = ";")))

# Options 1-3: Make into a dataframe
sample_targets_DF <- data.frame("swissprot" = sample_targets_swissprot, "ensg" = sample_targets_ensg)

# Write to CSV
# Option 1
write.csv(sample_targets_DF, paste(main_path, "Linear Model/TP53/tp53_curated_targets.csv", sep = ""))
write.csv(sample_targets_DF, paste(main_path, "Linear Model/TP53/foxa1_curated_targets.csv", sep = ""))

# Option 2
write.csv(sample_targets_DF, paste(main_path, "Linear Model/TP53/tp53_chipeat_targets.csv", sep = ""))
write.csv(sample_targets_DF, paste(main_path, "Linear Model/FOXA1/foxa1_chipeat_targets.csv", sep = ""))

# Option 3
write.csv(sample_targets_DF, paste(main_path, "Linear Model/allgene_targets.csv", sep = ""))

############################################################
############################################################
# 6. FIND SAMPLES WITH MUTATION AND CNA IN SAME GENE
############################################################
############################################################
#' Given a sample set of interest, and whether we are interested in amplifications or
#' deletions, identifies the overlap of samples that have both a mutation event 
#' and an amplification or deletion in the given gene. Eliminate overlap in the 
#' given gene and return this new patient list
#' @param mut_count_matrix a mutation count matrix with sample IDs as column names (XXXX-XX)
#' and a column called Gene_Symbol with the external gene name
#' @param cna_df a CNA data frame with the raw copy number of each gene (rownames are gene
#' ENSG IDs) for each sample (colnames are sample IDs, XXXX-XX)
#' @param gn the external gene name of the gene of interest
#' @param ensg the ENSG ID of the gene of interest
#' @param samples a vector of patient samples (XXXX-XX) in the cohort of interest
#' @param ampOrDel either "Amplification" to indicate we are interested in overlap with 
#' amplifications, or "Deletion" to indicate we are interested in overlap with deletions
find_samps_w_muts_and_cnas_in_same_gene <- function(mut_count_matrix, cna_df, gn, 
                                                    ensg, samples, ampOrDel) {
  
  colnames(mut_count_matrix) <- unlist(lapply(colnames(mut_count_matrix), function(x) {
    return(unlist(strsplit(x, "-", fixed = TRUE))[1])
  }))
  colnames(cna_df) <- unlist(lapply(colnames(cna_df), function(x) {
    return(unlist(strsplit(x, "-", fixed = TRUE))[1])
  }))
  
  samples_mut <- colnames(mut_count_matrix)[intersect(which(as.numeric(unlist(
    mut_count_matrix[mut_count_matrix$Gene_Symbol == gn,])) == 1), 
    which(colnames(mut_count_matrix) %in% samples))]
  
  if(ampOrDel == "Amplification") {
    samples_cna <- colnames(cna_df)[intersect(which(as.numeric(unlist(cna_df[rownames(cna_df) == ensg,])) > 2), 
                                              which(colnames(cna_df) %in% samples))]
  } else {
    samples_cna <- colnames(cna_df)[intersect(which(as.numeric(unlist(cna_df[rownames(cna_df) == ensg,])) < 2), 
                                              which(colnames(cna_df) %in% samples))]
  }
  
  # Plot a Venn Diagram to visualize the overlap between these two groups
  plt <- venn.diagram(list(samples_mut, samples_cna), category.names = c("Mutation", ampOrDel), 
               filename = NULL, output = TRUE, lwd = 2, lty = 'blank', fill = c("red", "blue"), 
               cex = 2, fontface = "bold", fontfamily = "sans", cat.pos = 180,
               cat.fontfamily = "sans", cat.cex = 1)
              #, cat.cex = 0.6, cat.fontface = "bold",
              #, rotation = 1)
  grid::grid.draw(plt)
  
  # Get & return the samples without overlap
  samples_no_overlap <- c(setdiff(samples_mut, samples_cna), setdiff(samples_cna, samples_mut))
  return(samples_no_overlap)
}

gn <- "TP53"
ensg <- "ENSG00000141510"
  
gn <- "PIK3CA"
ensg <- "ENSG00000121879"


############################################################
############################################################
# 7. CREATE SAMPLE ID FILES WITH NO OVERLAP IN MUTATION 
# BETWEEN TWO GIVEN GENES
############################################################
############################################################
mutation_count_matrix <- read.csv(paste0(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_IntersectPatients.csv"), 
                                  header = TRUE, row.names = 1, check.names = FALSE)
mutation_count_matrix_MN <- read.csv(paste0(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_nonsense_IntersectPatients.csv"), 
                                  header = TRUE, row.names = 1, check.names = FALSE)

gene1 <- "TP53"
gene2 <- "PIK3CA"

#' Get patients without a mutation the given gene1 AND gene2 (eliminate overlap)
#' @param mutation_count_matrix a mutation count matrix of all samples and genes
#' @param gene1 the gene name of the first gene 
#' @param gene2 the gene name of the second gene 
eliminate_mutation_overlap <- function(mutation_count_matrix, gene1, gene2) {
  mut_count_mat_sub <- mutation_count_matrix[(mutation_count_matrix$Gene_Symbol == gene1) | 
                                               (mutation_count_matrix$Gene_Symbol == gene2),]
  g1_mut_indices <- unlist(lapply(1:ncol(mut_count_mat_sub), function(i) 
    ifelse(as.numeric(mut_count_mat_sub[1,i]) > 0, i, NA)))
  g2_mut_indices <- unlist(lapply(1:ncol(mut_count_mat_sub), function(i) 
    ifelse(as.numeric(mut_count_mat_sub[2,i]) > 0, i, NA)))
  
  g1_mut_indices <- g1_mut_indices[!is.na(g1_mut_indices)]
  g2_mut_indices <- g2_mut_indices[!is.na(g2_mut_indices)]
  
  both_mut <- intersect(g1_mut_indices, g2_mut_indices)
  
  mut_count_mat_sub_excl_doub_mut <- mut_count_mat_sub[, -both_mut]
  patients_no_mut_overlap <- unlist(lapply(colnames(mut_count_mat_sub_excl_doub_mut), function(x)
    unlist(strsplit(x, "-", fixed = TRUE))[1]))
  patients_no_mut_overlap <- patients_no_mut_overlap[!(patients_no_mut_overlap == "Gene_Symbol")]
  
  return(patients_no_mut_overlap)
}


# Call function
patients_no_mut_overlap <- eliminate_mutation_overlap(mutation_count_matrix, gene1, gene2)
patients_no_mut_overlap_MN <- eliminate_mutation_overlap(mutation_count_matrix_MN, gene1, gene2)

# write to files
write(patients_no_mut_overlap, paste0(main_path, "Patient Subsets/NotTP53andPIK3CAmut_patient_ids.txt"))
write(patients_no_mut_overlap_MN, paste0(main_path, "Patient Subsets/NotTP53andPIK3CAMNmut_patient_ids.txt"))


# Get overlap with other patient subsets
patient_set_lt20pamp <- read.table(paste0(main_path, "Patient Subsets/LessThan20PercAmp_patient_ids.txt"), header = TRUE)[,1]
patient_set_lumA <- read.table(paste0(main_path, "Patient Subsets/Luminal.A_patient_ids.txt"), header = TRUE)[,1]
patient_set_lumB <- read.table(paste0(main_path, "Patient Subsets/Luminal.B_patient_ids.txt"), header = TRUE)[,1]
patient_set_lumAB <- read.table(paste0(main_path, "Patient Subsets/Luminal.A.B_patient_ids.txt"), header = TRUE)[,1]
patient_set_basal <- read.table(paste0(main_path, "Patient Subsets/Basal_patient_ids.txt"), header = TRUE)[,1]
patient_set_her2 <- read.table(paste0(main_path, "Patient Subsets/HER2_patient_ids.txt"), header = TRUE)[,1]
patient_set_normLike <- read.table(paste0(main_path, "Patient Subsets/Normal-like_patient_ids.txt"), header = TRUE)[,1]

lumA_noMutOL <- intersect(patient_set_lumA, patients_no_mut_overlap_MN)
lumB_noMutOL <- intersect(patient_set_lumB, patients_no_mut_overlap_MN)
basal_noMutOL <- intersect(patient_set_basal, patients_no_mut_overlap_MN)
her2_noMutOL <- intersect(patient_set_her2, patients_no_mut_overlap_MN)
normLike_noMutOL <- intersect(patient_set_normLike, patients_no_mut_overlap_MN)
lt20pamp_noMutOL <- intersect(patient_set_lt20pamp, patients_no_mut_overlap_MN)

lt20pamp_lumA_noMutOL <- intersect(patient_set_lt20pamp, lumA_noMutOL)
lt20pamp_lumB_noMutOL <- intersect(patient_set_lt20pamp, lumB_noMutOL)
lt20pamp_basal_noMutOL <- intersect(patient_set_lt20pamp, basal_noMutOL)
lt20pamp_her2_noMutOL <- intersect(patient_set_lt20pamp, her2_noMutOL)
lt20pamp_normLike_noMutOL <- intersect(patient_set_lt20pamp, normLike_noMutOL)

# Write these subsets to files
write(lumA_noMutOL, paste0(main_path, "Patient Subsets/Luminal.A.NotTP53andPIK3CAMNmut_patient_ids.txt"))
write(lumB_noMutOL, paste0(main_path, "Patient Subsets/Luminal.B.NotTP53andPIK3CAMNmut_patient_ids.txt"))
write(basal_noMutOL, paste0(main_path, "Patient Subsets/Basal.NotTP53andPIK3CAMNmut_patient_ids.txt"))
write(her2_noMutOL, paste0(main_path, "Patient Subsets/HER2.NotTP53andPIK3CAMNmut_patient_ids.txt"))
write(normLike_noMutOL, paste0(main_path, "Patient Subsets/Normal-like.NotTP53andPIK3CAMNmut_patient_ids.txt"))

write(lt20pamp_lumA_noMutOL, paste0(main_path, "Patient Subsets/LT20PAmp.NotTP53andPIK3CAMNmut.LumA_patient_ids.txt"))
write(lt20pamp_lumB_noMutOL, paste0(main_path, "Patient Subsets/LT20PAmp.NotTP53andPIK3CAMNmut.LumB_patient_ids.txt"))
write(lt20pamp_basal_noMutOL, paste0(main_path, "Patient Subsets/LT20PAmp.NotTP53andPIK3CAMNmut.Basal_patient_ids.txt"))
write(lt20pamp_her2_noMutOL, paste0(main_path, "Patient Subsets/LT20PAmp.NotTP53andPIK3CAMNmut.HER2_patient_ids.txt"))
write(lt20pamp_normLike_noMutOL, paste0(main_path, "Patient Subsets/LT20PAmp.NotTP53andPIK3CAMNmut.NormLike_patient_ids.txt"))


############################################################
############################################################
# 8. FIND SAMPLES WITH MUTATION IN ONE GENE & NO MUTATION
# IN ANOTHER (INDEPENDENT POSITIVE AND NEGATIVE SETS)
############################################################
############################################################

mutation_count_matrix <- read.csv(paste0(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense.csv"), 
                                  header = TRUE, row.names = 1, check.names = FALSE)

gene1 <- "TP53"
gene2 <- "PIK3CA"

#' Get patients with a mutation just in given gene1 and not in given gene2
#' @param mutation_count_matrix a mutation count matrix of all samples and genes
#' @param gene1 the gene name of the first gene (with mutation)
#' @param gene2 the gene name of the second gene (no mutation)
get_independent_mut_patient_sets <- function(mutation_count_matrix, gene1, gene2) {
  # Subset the mutation count matrix to just these two genes
  mutation_count_matrix_sub <- mutation_count_matrix[rownames(mutation_count_matrix) %in% c(gene1, gene2),]
  
  # Get the patients with a mutation in the first gene but no mutation in the second
  g1_mut_g2_noMut <- mutation_count_matrix_sub[, intersect(which(mutation_count_matrix_sub[rownames(mutation_count_matrix_sub) == gene1,] == 1),
                                                           (which(mutation_count_matrix_sub[rownames(mutation_count_matrix_sub) == gene2,] == 0)))]
  g1_mut_g2_noMut_pats <- unique(unlist(lapply(colnames(g1_mut_g2_noMut), function(x) 
    unlist(strsplit(x, "-", fixed = TRUE))[3])))
  
  return(g1_mut_g2_noMut_pats)
}


g1_mut_g2_noMut_pats <- get_independent_mut_patient_sets(mutation_count_matrix, gene1, gene2)
g2_mut_g1_noMut_pats <- get_independent_mut_patient_sets(mutation_count_matrix, gene2, gene1)


#' Get and return patients with no mutation in the given gene
#' @param mutation_count_matrix a mutation count matrix of all samples and genes
#' @param gene the gene name of the gene whose patients (with a mutation in this gene)
#' we want to exclude
get_patients_with_no_mut_in_gene <- function(mutation_count_matrix, gene) {
  # Subset the mutation count matrix to just this gene
  mutation_count_matrix_sub <- mutation_count_matrix[rownames(mutation_count_matrix) == gene,]
  
  # Get the patients without a mutation in the gene
  gene_noMut <- mutation_count_matrix_sub[, which(mutation_count_matrix_sub[rownames(mutation_count_matrix_sub) == gene,] == 0)]
  gene_noMut_pats <- unique(unlist(lapply(colnames(gene_noMut), function(x) 
    unlist(strsplit(x, "-", fixed = TRUE))[3])))
  
  return(gene_noMut_pats)
}

g1_noMut_pats <- get_patients_with_no_mut_in_gene(mutation_count_matrix, gene1)
g2_noMut_pats <- get_patients_with_no_mut_in_gene(mutation_count_matrix, gene2)

# Get the intersection with our intersecting ids
intersecting_pats <- read.table(paste(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt", sep = ""))[,1]

g1_noMut_pats <- intersect(g1_noMut_pats, intersecting_pats)
g2_noMut_pats <- intersect(g2_noMut_pats, intersecting_pats)

noMut_pats_list <- list(gene1 = g1_noMut_pats, gene2 = g2_noMut_pats)

# Plot a Venn diagram of the overlap
ggVennDiagram(noMut_pats_list, label_alpha = 0, category.names = c(paste("No Mut,", gene1), paste("No Mut,", gene2)), set_color = "black",
              set_size = 10, label_size = 8, edge_size = 0) +
  ggplot2::scale_fill_gradient(low="cornsilk1", high = "cadetblue3")

# Write these to files
write(g1_noMut_pats, paste0(main_path, paste0("Patient Subsets/", paste0(gene1, ".NoMut_patient_ids.txt"))))
write(g2_noMut_pats, paste0(main_path, paste0("Patient Subsets/", paste0(gene2, ".NoMut_patient_ids.txt"))))

# Get the intersection with LT20PAmp patients
lt20pamp_pats <- read.table(paste0(main_path, "Patient Subsets/LessThan20PercAmp_patient_ids.txt"))[,1]

g1_noMut_pats_lt20pamp <- intersect(g1_noMut_pats, lt20pamp_pats)
g2_noMut_pats_lt20pamp <- intersect(g2_noMut_pats, lt20pamp_pats)

noMut_pats_lt20pamp_list <- list(gene1 = g1_noMut_pats_lt20pamp, gene2 = g2_noMut_pats_lt20pamp)

# Plot a Venn diagram of the overlap
ggVennDiagram(noMut_pats_lt20pamp_list, label_alpha = 0, category.names = c(paste("No Mut,", gene1), paste("No Mut,", gene2)), set_color = "black",
              set_size = 10, label_size = 8, edge_size = 0) +
  ggplot2::scale_fill_gradient(low="cornsilk1", high = "cadetblue3")

# Write these to files
write(g1_noMut_pats_lt20pamp, paste0(main_path, paste0("Patient Subsets/", paste0("LT20PAmp.", paste0(gene1, ".NoMut_patient_ids.txt")))))
write(g2_noMut_pats_lt20pamp, paste0(main_path, paste0("Patient Subsets/", paste0("LT20PAmp.", paste0(gene2, ".NoMut_patient_ids.txt")))))



