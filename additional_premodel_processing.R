############################################################
### Additional Pre-model Inputs & Processing 
### Written By: Sara Geraghty, April 2021
############################################################

options("install.lock"=FALSE)
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
  # 4. Construct a data frame of targets for a given regulatory protein of interest (for model
      # testing purposes)
  # 5. Find samples that have both an amplification and mutation, or deletion and mutation, in the same gene

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
expression_df_qn <- fread(paste(main_path, "Expression/expression_quantile_normalized_sklearn.csv", sep = ""),
                          header = TRUE)
expression_df_qn <- fread(paste(main_path, "Expression/ALL_CT_expression_quantile_normalized_sklearn.csv", sep = ""),
                          header = TRUE)
colnames(expression_df_qn)[1] <- 'ensg_id'

# 4. Rank-Normalized
expression_df_rn <- fread(paste(main_path, "Expression/expression_rank_norm_DF.csv", sep = ""),
                          header = TRUE)
colnames(expression_df_rn)[1] <- 'ensg_id'


# For pan-cancer, limit to just cancer TCGA samples
expression_df_qn_sub <- expression_df_qn[, c(1, (which(grepl("TCGA", colnames(expression_df_qn)[2:ncol(expression_df_qn)]))+1)), with = F]
colnames(expression_df_qn_sub)[2:ncol(expression_df_qn_sub)] <- unlist(lapply(colnames(expression_df_qn_sub)[2:ncol(expression_df_qn_sub)], function(x)
  paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))



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


# For pan-cancer, make sure we are not missing any BRCA data
methylation_df_M_brca <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Methylation/methylation_DF_M_cancer_only.csv", header = T)
missing_brca_cols <- setdiff(colnames(methylation_df_M_brca), colnames(methylation_df_M))
methylation_df_M_full <- merge(methylation_df_M, methylation_df_M_brca[,which(colnames(methylation_df_M_brca) %in% c("Gene_Symbol", missing_brca_cols)), with = F], by = "Gene_Symbol")
write.csv(methylation_df_M_full, paste(main_path, "Methylation/methylation_DF_M_cancer_only.csv", sep = ""))


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
mutation_targ_df <- fread(paste(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_nonsynonymous_ALL_inclNonmut.csv", sep = ""), 
                             header = TRUE)
#colnames(mutation_targ_df) <- unlist(lapply(colnames(mutation_targ_df), function(x) str_replace_all(string = x, pattern = " ", repl = "")))
mutation_targ_df <- mutation_targ_df[,2:ncol(mutation_targ_df)]
#mutation_targ_df <- fread(paste(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_nonsynonymous_ALL_inclNonmut.csv", sep = ""), 
                          #header = TRUE)
colnames(mutation_targ_df)[1] <- 'Gene_Symbol'
#colnames(mutation_targ_df)[2:ncol(mutation_targ_df)] <- unlist(lapply(colnames(mutation_targ_df)[2:ncol(mutation_targ_df)], function(x) 
#  paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-"))) # Fix the sample IDs

# If needed, add the additional un-included genes (those with no mutations in any patient)
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

# If needed, add the additional un-included patients (those with no mutations in any gene)
add_missing_patients_to_mut_count_mat <- function(mut_count_matrix, intersecting_patients) {
  existing_patients <- unlist(lapply(colnames(mut_count_matrix)[2:ncol(mut_count_matrix)], function(x) 
    unlist(strsplit(x, "-", fixed = TRUE))[1]))
  print(head(existing_patients))
  missing_patients <- setdiff(intersecting_patients, existing_patients)
  print(length(missing_patients))
  new_pat_df <- data.frame(matrix(nrow = nrow(mut_count_matrix), ncol = length(missing_patients)))
  new_pat_df[is.na(new_pat_df)] <- 0
  colnames(new_pat_df) <- paste0(missing_patients, "-01A")
  
  new_mut_count_mat <- cbind(mut_count_matrix, new_pat_df)
  return(new_mut_count_mat)
}
intersecting_patients_brca <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/unique_brca_patient_ids_2.txt", 
                                         header = TRUE, row.names = 1, check.names = FALSE)[,1]
intersecting_patients_pc <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/unique_patients_ids_2.txt", 
                                       header = FALSE, check.names = FALSE)[,1]

mutation_targ_df <- add_missing_patients_to_mut_count_mat(mutation_targ_df, intersecting_patients)

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
cna_df_raw <- read.csv(paste(main_path, "CNV/Gene-level Raw/CNV_DF_AllGenes_CancerOnly_MergedIsoforms.csv", sep = ""), 
                       header = TRUE, row.names = 1, check.names = FALSE)
# Bucketed Copy Number (Incl. Amp)
cna_df_bucket_inclAmp <- read.table(paste(main_path, "CNV/Gene-level Raw/CNV_DF_bucketed_inclAmp_AllGenes.tsv", sep = ""), 
                header = TRUE, row.names = 1, check.names = FALSE)
# Bucketed Copy Number (Excl. Amp)
cna_df_bucket_exclAmp <- read.table(paste(main_path, "CNV/Gene-level Raw/CNV_DF_bucketed_exclAmp_AllGenes.tsv", sep = ""), 
                header = TRUE, row.names = 1, check.names = FALSE)


# If needed, combine CNA isoforms for pan-cancer
cna_df_unique_ensg_ids <- unique(unlist(lapply(rownames(cna_df_raw), function(e)
  unlist(strsplit(e, ".", fixed = T))[1])))
ensg_id_row_indices <- lapply(cna_df_unique_ensg_ids, function(e) 
  which(grepl(e, rownames(cna_df_raw))))
names(ensg_id_row_indices) <- cna_df_unique_ensg_ids

cna_df_new_rows <- mclapply(1:length(ensg_id_row_indices), function(i) {
  ensg <- names(ensg_id_row_indices)[i]
  ind <- ensg_id_row_indices[[i]]
  print(ensg)
  print(paste(i, paste("/", length(ensg_id_row_indices))))
  if(length(ind) > 1) {
    counter <- 1
    vals <- NA
    while(counter < length(ind)) {
      if(is.na(vals)) {
        row1 <- unlist(cna_df_raw[ind[counter], ])
        row2 <- unlist(cna_df_raw[ind[counter+1], ])
        vals <- dplyr::coalesce(row1, row2)
        counter <- counter + 2
      } else {
        row <- unlist(cna_df_raw[ind[counter], ])
        vals <- dplyr::coalesce(vals, row2)
      }
    }
    return(vals)
  } else {
    return(unlist(cna_df_raw[ind,]))
  }
})
cna_df_new <- as.data.frame(do.call(rbind, cna_df_new_rows))
rownames(cna_df_new) <- cna_df_unique_ensg_ids

write.csv(cna_df_new, paste(main_path, "CNV/Gene-level Raw/CNV_DF_AllGenes_CancerOnly_MergedIsoforms.csv", sep = ""))

cna_df_new <- read.csv(paste(main_path, "CNV/Gene-level Raw/CNV_DF_AllGenes_CancerOnly_MergedIsoforms.csv", sep = ""),
                       header = T, check.names = F, row.names = 1)

#' If there are no isoforms, just adjust the gene names to remove the part after the '.'
#' @param cna_df one of the above CNA DFs
adj_rownams <- function(cna_df) {
  rownames(cna_df) <- make.names(unlist(lapply(rownames(cna_df), function(x) 
    unlist(strsplit(x, ".", fixed = TRUE))[1])), unique = T)
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



patient_df_brca_blca_hnsc <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/patient_dataframe_ntnm_normAge_inclBRCA.BLCA.HNSCsubtypes.csv", sep = ""), 
                                      header = TRUE, row.names = 1)

patient_df_coad_read <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/patient_dataframe_ntnm_normAge_COAD_READ.csv", sep = ""), 
                                      header = TRUE, row.names = 1)
# With subtypes
patient_df_coad_read_subtypes <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/patient_dataframe_ntnm_MN_normAge_inclCOADREADsubtypes.csv", sep = ""), 
                                 header = TRUE, row.names = 1)
patient_df_coad_read_cmssubtypes <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/patient_dataframe_ntnm_MN_normAge_inclCOADREADcmssubtypes.csv", sep = ""), 
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
#' @param dataset either 'tcga', 'metabric', 'chinese_tn', icgc', or 'cptac3'
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
      patient <- unlist(strsplit(sample_barcode, "-", fixed = TRUE))[1]
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

# Remove rows/ columns that are entirely NA
combined_pat_samp_df <- combined_pat_samp_df[rowSums(is.na(combined_pat_samp_df)) != ncol(combined_pat_samp_df),]
combined_pat_samp_df <- combined_pat_samp_df[,colSums(is.na(combined_pat_samp_df)) != nrow(combined_pat_samp_df)]

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


# For per-cancer files (patient files stored in list 'pan_cancer_patient_dfs')
combined_pat_samp_dfs <- lapply(pan_cancer_patient_dfs, function(pat_df) {
  comb_df <- combine_patient_and_samp_dfs(pat_df, sample_df_to, 'tcga')
  comb_df <- comb_df[rowSums(is.na(comb_df)) != ncol(comb_df),]
  comb_df <- comb_df[,colSums(is.na(comb_df)) != nrow(comb_df)]
  return(comb_df)
})
names(combined_pat_samp_dfs) <- names(pan_cancer_patient_dfs)

# Write to file
lapply(1:length(combined_pat_samp_dfs), function(i) {
  cancer_name <- names(combined_pat_samp_dfs)[i]
  fn <- paste0("Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_normAge_ntnm_", 
               paste0(cancer_name, "_inclSubtypes.csv"))
  write.csv(combined_pat_samp_dfs[[i]], paste0(main_path, fn))
})

#
#
#

# Read this back (generic)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF.csv", sep = ""),
                              header = TRUE, row.names = 1)

patient_sample_files <- list.files(paste0(main_path, "Linear Model/Patient and Sample DFs/"), pattern = "_inclSubtypes.csv")
patient_sample_files <- patient_sample_files[grepl("combined", patient_sample_files)]
patient_sample_dfs <- lapply(patient_sample_files, function(f)
  read.csv(paste0(main_path, paste0("Linear Model/Patient and Sample DFs/", f)), header = T, row.names = 1))
names(patient_sample_dfs) <- unlist(lapply(patient_sample_files, function(f) 
  unlist(strsplit(f, "_", fixed = T))[10]))

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
cna_patients <- get_unique_patients(colnames(cna_df_new))    # 732; 7500
methylation_patients <- get_unique_patients(colnames(methylation_df)[2:ncol(methylation_df)])  #-1]) # 734; 6491
#exp_samp_ids <- unlist(lapply(colnames(expression_df)[2:ncol(expression_df)], function(x)
  #paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))
#expression_patients <- get_unique_patients(exp_samp_ids)  # 738; 7375
expression_patients <- get_unique_patients(colnames(expression_df)[2:ncol(expression_df)])  # 732; 7375
mutation_targ_patients <- unlist(get_unique_patients(colnames(mutation_targ_df)[!(colnames(mutation_targ_df) %in% c('Swissprot', 'Gene_Symbol'))])) # 978; 7258; #737
mutation_regprot_patients <- unique(unlist(lapply(mutation_regprot_df$Patient, function(x) 
  unlist(strsplit(x, ";", fixed = TRUE)))))   # 615; 6732 -- ignore this (6875 with nonsense, 6889 with nonsynonymous)

# If needed
#patient_sample_df <- as.data.frame(patient_sample_df)
#rownames(patient_sample_df) <- unlist(patient_sample_df[,1]) 
#patient_sample_df <- patient_sample_df[,2:ncol(patient_sample_df)]

#clin_samp_ids <- unlist(lapply(rownames(patient_sample_df), function(x)
#  paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))
clinical_patients <- get_unique_patients(rownames(patient_sample_df))  # 732; 6960; 313 for COAD/READ
clinical_patients_perCt <- unique(unlist(lapply(patient_sample_dfs, function(df) get_unique_patients(rownames(df)))))  # 6960

intersecting_patients <- intersect(cna_patients, 
                                   intersect(methylation_patients,
                                             intersect(expression_patients,
                                                       intersect(mutation_targ_patients, clinical_patients))))
print(length(intersecting_patients)) # 727 in BRCA (605 using UCSF genotype PCs, 609 if using WashU genotype PCs for BRCA; 5722 if using WashU PCs P-C)

#182 for COAD/READ

# Remove hypermutator patients
intersecting_patients_hypermutFilt <- intersecting_patients[!(intersecting_patients %in% total_excluded_pats)]

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

# Do this for the per-cancer patient/sample files
patient_sample_dfs_sub <- lapply(patient_sample_dfs, function(df) 
  subset_by_intersecting_ids(intersecting_patients, df, F))
names(patient_sample_dfs_sub) <- names(patient_sample_dfs)

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

# Add Swissprot IDs to mutation DF
mutation_targ_df_sub$Swissprot <- unlist(lapply(mutation_targ_df_sub$Gene_Symbol, function(x) 
  paste(unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'uniprot_gn_id']), collapse = ";")))

# Re-write these all to to files
write.csv(cna_df_sub, paste(main_path, "Linear Model/Tumor_Only/CNV/CNA_AllGenes_CancerOnly_IntersectPatients.csv", sep = ""))
write.csv(methylation_df_sub, paste(main_path, "Linear Model/Tumor_Only/Methylation/methylation_M_CancerOnly_IntersectPatients.csv", sep = ""))
write.csv(expression_df_sub, paste(main_path, "Linear Model/Tumor_Only/Expression/expression_tmm_IntersectPatients.csv", sep = ""))
write.csv(mutation_targ_df_sub, paste(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_IntersectPatients.csv", sep = ""))
write.csv(mutation_regprot_df_sub, paste(main_path, "Linear Model/Tumor_Only/Regprot_Mutation/iprotein_results_missense_IntersectPatients.csv", sep = ""))
write.csv(patient_sample_df_sub, paste(main_path, "Linear Model/Tumor_Only/Patient/combined_patient_sample_IntersectPatients.csv", sep = ""))

lapply(1:length(patient_sample_dfs_sub), function(i) {
  ct <- names(patient_sample_dfs_sub)[i]
  fn <- paste0("Linear Model/Tumor_Only/Patient/combined_patient_sample_cibersort_total_frac_", 
               paste0(ct, "_inclSubtypes_IntersectPatientsWashU.csv"))
  write.csv(patient_sample_dfs_sub[[i]], paste0(main_path, fn))
})

#OPT:
write.table(intersecting_patients, paste(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt", sep = ""))
intersecting_patients <- read.table(paste(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt", sep = ""))[,1]


############################################################
############################################################
# 4. MAKE TESTER TARGETS DATAFRAME FOR GIVEN REG. PROTEIN
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


############################################################
############################################################
# 9. MERGE PAN-CANCER EXPRESSION DFs INTO ONE
############################################################
############################################################

expression_dfs <- list.files(paste0(main_path, "Expression/"), 
                             pattern = "_expression_quantile_normalized_sklearn.csv")

merged_expression_df <- NA

for(i in 1:length(expression_dfs)) {
  expression_df <- distinct(fread(paste0(main_path, paste0("Expression/", expression_dfs[i])), 
                         header = TRUE))
  colnames(expression_df)[1] <- "ensg_id"
  if(i > 1) {
    merged_expression_df <- merge(merged_expression_df, expression_df, 
                                  by = "ensg_id", all = TRUE)
  } else {merged_expression_df <- expression_df}
}

write.csv(merged_expression_df, paste0(main_path, "Expression/ALL_CT_expression_quantile_normalized_sklearn.csv"))