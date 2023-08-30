############################################################
### Process METABRIC Data
### Written By: Sara Geraghty, April 2022
############################################################

library(GenomicRanges)
library(TRONCO)

# METABRIC Data Retrieved from cBioPortal and downloaded as a .tar file, which was unzipped

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/cBioPortal/brca_metabric/"

output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/cBioPortal/METABRIC/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


# The goal of this file is to take in METABRIC data files in order to run our model on another
# dataset aside from the TCGA (as a validation of our findings). This file imports each data type of 
# interest and formats them to match the TCGA files, such that they will be able to flow through the 
# various existing pipelines we have for each data type.


############################################################
### IMPORT MUTATION DATA
############################################################
# Import mutation data file
metabric_mutation_df <- read.table(paste0(main_path, "data_mutations.txt"), header = TRUE, 
                                   sep = "\t", check.names = FALSE)

# Generate the mutation count matrices, using helper functions from maf_helper_functions.R
clin_filename <- paste(main_path, "data_clinical_sample_sub.txt", sep = "")
clinical_df <- read.table(clin_filename, header = TRUE, sep = ",")
metabric_mutation_df <- metabric_mutation_df[metabric_mutation_df$Tumor_Sample_Barcode %fin% clinical_df$SAMPLE_ID,]
metabric_mutation_df_nonsynonymous <- metabric_mutation_df[(metabric_mutation_df$Variant_Classification == "Missense_Mutation") | 
                                           (metabric_mutation_df$Variant_Classification == "Nonsense_Mutation") |
                                           (metabric_mutation_df$Variant_Classification == "Splice_Site") |
                                           (metabric_mutation_df$Variant_Classification == "Nonstop_Mutation"),]
metabric_mutation_df_silent <- metabric_mutation_df[metabric_mutation_df$Variant_Classification == "Silent",]


maf_nonsynon <- read.maf(metabric_mutation_df_nonsynonymous)
maf_silent <- read.maf(metabric_mutation_df_silent, vc_nonSyn = "Silent")
  
mut_count_matrix_nonsyn <- get_mut_count_matrix(maf_nonsynon)
mut_count_matrix_silent <- get_mut_count_matrix(maf_silent)

# Now, append rows of '0s' for those genes that are not mutated in any patient
allgene_targets <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/allgene_targets.csv",
                            header = TRUE, check.names = FALSE, row.names = 1)
all_gene_names <- unlist(lapply(allgene_targets$ensg, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == x,'external_gene_name'])), collapse = ";")))

get_unincluded_gene_df <- function(all_gene_names, mut_count_matrix) {
  unincluded_genes <- setdiff(all_gene_names, rownames(mut_count_matrix))
  unincluded_gene_df <- data.frame(matrix(nrow = length(unincluded_genes), 
                                          ncol = ncol(mut_count_matrix)))
  unincluded_gene_df[is.na(unincluded_gene_df)] <- 0
  rownames(unincluded_gene_df) <- unincluded_genes
  colnames(unincluded_gene_df) <- colnames(mut_count_matrix)
  return(unincluded_gene_df)
}

mut_count_matrix_nonsyn_full <- rbind(mut_count_matrix_nonsyn, 
                                 get_unincluded_gene_df(all_gene_names, mut_count_matrix_nonsyn))
mut_count_matrix_silent_full <- rbind(mut_count_matrix_silent, 
                                 get_unincluded_gene_df(all_gene_names, mut_count_matrix_silent))

write.csv(mut_count_matrix_nonsyn_full, paste0(output_path, "Mutation/Mutation Count Matrices/mut_count_matrix_nonsynonymous.csv"))
write.csv(mut_count_matrix_silent_full, paste0(output_path, "Mutation/Mutation Count Matrices/mut_count_matrix_silent.csv"))

# Read back if needed
mut_count_matrix_nonsyn_full <- read.csv(paste0(output_path, "Mutation/Mutation Count Matrices/mut_count_matrix_nonsynonymous.csv"),
                                         header = T, check.names = F, row.names = 1)
mut_count_matrix_silent_full <- read.csv(paste0(output_path, "Mutation/Mutation Count Matrices/mut_count_matrix_silent.csv"),
                                         header = T, check.names = F, row.names = 1)

############################################################
### FILTER PATIENTS THAT EXCEED A FIXED HYPERMUTATOR THRESHOLD, IF NEEDED
############################################################
# Using uniform threshold of 368, as from paper (see TCGA analysis)
thres <- 368

# Limit to those samples below this threshold
mut_count_matrix_nonsyn_colsums <- unlist(lapply(1:ncol(mut_count_matrix_nonsyn_full), function(i) {
  vals <- as.integer(unlist(mut_count_matrix_nonsyn_full[,i]))
  return(length(vals[vals > 0]))
}))

# Visualize this 
hist(mut_count_matrix_nonsyn_colsums, main = "Total Mutation Counts per Patient, METABRIC", xlab = "Total # of Mutations")

#mut_count_matrix_nonsyn_sub <- mut_count_matrix_nonsyn_full[, which(mut_count_matrix_nonsyn_colsums < thres)]

############################################################
### IMPORT CNA DATA
############################################################
# Import CNA data file (0 is normal, -1 is deletion of 1 copy, -2 is deletion of 2 copies, 2+ is amplification)
# Will have to use bucketing method for this CNA data
metabric_cna_df <- read.table(paste0(main_path, "data_cna.txt"), header = TRUE, 
                              sep = "\t", check.names = FALSE)

# Remove the Entrez ID & HUGO ID columns, use ENSG IDs instead
ensg_ids <- unlist(lapply(metabric_cna_df$Hugo_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'ensembl_gene_id'])), collapse = ";")))
metabric_cna_df <- metabric_cna_df[ensg_ids != "",]
ensg_ids <- ensg_ids[ensg_ids != ""]
metabric_cna_df <- metabric_cna_df[,3:ncol(metabric_cna_df)]
rownames(metabric_cna_df) <- make.names(ensg_ids, unique = TRUE)

# Write to file
write.csv(metabric_cna_df, paste0(output_path, "CNV/CNA_DF_AllGenes_CancerOnly.csv"))


############################################################
### IMPORT EXPRESSION DATA
############################################################
# Note that that expression data is microarray data, which needs to be handled differently -
# it appears to already be normalized, filtered, and log2 transformed (Potentially an RMA normalization?)
metabric_expression_df <- read.table(paste0(main_path, "data_mrna_agilent_microarray.txt"), 
                                     header = TRUE, sep = "\t", check.names = FALSE)

# Reverse the log2-normalization in order to feed into immunedeconv ("immune_cell_deconvolution.R" pipeline)
metabric_expression_df_unlogged <- apply(metabric_expression_df[,3:ncol(metabric_expression_df)], c(1,2), function(x) 2^x)
metabric_expression_df_unlogged <- cbind(metabric_expression_df[,1:2], metabric_expression_df_unlogged)
write.table(metabric_expression_df_unlogged, paste0(main_path, "data_mrna_agilent_microarray_unlogged.txt"))

# Filter out genes with more than 50% missing values
rows_with_more_than_50p <- unlist(lapply(1:nrow(metabric_expression_df), function(i) {
    r <- as.numeric(metabric_expression_df[i, 3:ncol(metabric_expression_df)])
    print(head(r))
    if (length(r[(r == "") | (is.na(r))]) < (0.5 * length(r))) {return(i)}
    else {return(NA)}
  }))
rows_with_more_than_50p <- rows_with_more_than_50p[!is.na(rows_with_more_than_50p)]
metabric_expression_df <- metabric_expression_df[c(1:2, (rows_with_more_than_50p + 2)),]

# Convert the Hugo and Entrez symbols to ENSG IDs
ensg_ids <- unlist(lapply(metabric_expression_df$Hugo_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'ensembl_gene_id'])), collapse = ";")))
metabric_expression_df <- metabric_expression_df[ensg_ids != "",]
ensg_ids <- ensg_ids[ensg_ids != ""]
metabric_expression_df <- metabric_expression_df[,3:ncol(metabric_expression_df)]
rownames(metabric_expression_df) <- make.names(ensg_ids, unique = TRUE)

# Filter genes with mean expression < 5 or standard deviation smaller than 0.3
# Method taken from 10.1186/s13058-015-0618-8, Supp. Methods p. 3
metabric_row_means <- rowMeans(metabric_expression_df, na.rm = TRUE)
metabric_expression_df <- metabric_expression_df[which(metabric_row_means > 5),] # doesn't actually filter any in our case; range is 5.08 to 13.99
metabric_sds <- rowSds(as.matrix(metabric_expression_df))  # ranges from 0.106 to 3.0669
metabric_expression_df <- metabric_expression_df[which(metabric_sds > 0.3),]

# Write to file
write.csv(metabric_expression_df, paste0(output_path, "Expression/expression_log2intensity_filtBySD0.3_DF.csv"))


############################################################
### IMPORT METHYLATION DATA
############################################################
# This is promoter methylation (RRBS), bisulfite sequencing ONLY. This is different from the 
# TCGA, which uses Illumina450 sequencing arrays combined with SeSAMe methylation Beta estimation.
# (Here, a promoter region is defined as 500 bp upstream and 50 bp downstream from TSS).
# Quality control and filtering of regions with few reads has already been performed: https://www.nature.com/articles/s41467-021-25661-w#Sec11
# RRBS processing was done using the gpatterns package (https://github.com/tanaylab/gpatterns). File contains individual CpG methylation for each read.
metabric_methylation_df <- read.table(paste0(main_path, "data_methylation_promoters_rrbs.txt"), 
                                      header = TRUE, sep = "\t", check.names = FALSE)

rownames(metabric_methylation_df) <- metabric_methylation_df$Hugo_Symbol
metabric_methylation_df <- metabric_methylation_df[,2:ncol(metabric_methylation_df)]

write.csv(metabric_methylation_df, paste0(output_path, "Methylation/methylation_DF_RRBS.csv"))

# The range of these values is [0,1]. Can take the logit of this to get a normal distribution
metabric_methylation_df_M <- apply(metabric_methylation_df, c(1,2), function(x) {
  if((is.na(x)) | (length(x) == 0)) {return(NA)}
  if(x == 1) {x <- 0.9999999999}    # Logit of 1 is error
  if(x == 0) {x <- 0.0000000001}    # Logit of 0 is -Inf
  m <- log2(x / (1-x))
  return(m)
})

write.csv(metabric_methylation_df_M, paste0(output_path, "Methylation/methylation_DF_RRBS_logit.csv"))


############################################################
### IMPORT CLINICAL/ SAMPLE DATA
############################################################
# NOTE: all patients are female (no filtering needed)
metabric_clinical_patient_df <- read.table(paste0(main_path, "data_clinical_patient.txt"), 
                                           header = TRUE, sep = "\t", check.names = FALSE)
metabric_clinical_sample_df <- read.table(paste0(main_path, "data_clinical_sample.txt"), 
                                          header = TRUE, sep = "\t", check.names = FALSE)

# NOTE: Only primary tumor samples are included here (no normal or metastatic samples)
# Put these through the "create_sample_df" and "create_patient_df" pipelines

############################################################
### GET THE SAMPLES THAT OVERLAP BETWEEN ALL DATA TYPES
############################################################
case_info <- read.table(paste0(main_path, "case_lists/cases_all.txt"), header = TRUE,
                        sep = "\t", check.names = FALSE)

patient_df_patient_ids <- metabric_clinical_patient_df$PATIENT_ID
sample_df_patient_ids <- metabric_clinical_sample_df$PATIENT_ID
mutation_df_patient_ids <- metabric_mutation_df_nonsynonymous$Tumor_Sample_Barcode
expression_df_patient_ids <- colnames(metabric_expression_df)[3:ncol(metabric_expression_df)]
cna_df_patient_ids <- colnames(metabric_cna_df)[3:ncol(metabric_cna_df)]
methylation_df_patient_ids <- colnames(metabric_methylation_df)[2:ncol(metabric_methylation_df)]

intersecting_patient_ids <- intersect(patient_df_patient_ids, 
                                      intersect(sample_df_patient_ids,
                                                intersect(mutation_df_patient_ids, 
                                                          intersect(expression_df_patient_ids, 
                                                                    intersect(cna_df_patient_ids, methylation_df_patient_ids)))))
# There are 886 intersecting BRCA patient IDs
write.table(intersecting_patient_ids, "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/cBioPortal/brca_metabric/brca_intersecting_ids.txt",
            quote = FALSE, row.names = FALSE)

# Read back
intersecting_patient_ids <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/cBioPortal/brca_metabric/brca_intersecting_ids.txt")[,1]


############################################################
### LIMIT THESE DATA FILES TO THOSE SAMPLES WITH ALL DATA TYPES
############################################################
metabric_mutation_df_sub <- metabric_mutation_df[metabric_mutation_df$Tumor_Sample_Barcode %fin% 
                                                         intersecting_patient_ids,]
metabric_cna_df_sub <- metabric_cna_df[, c(1,2, (which(colnames(metabric_cna_df)[3:ncol(metabric_cna_df)] %fin% 
                                                         intersecting_patient_ids) + 2))]
metabric_expression_df_sub <- metabric_expression_df[, c(1,2, (which(colnames(metabric_expression_df)[3:ncol(metabric_expression_df)] %fin% 
                                                         intersecting_patient_ids) + 2))]
metabric_methylation_df_sub <- metabric_methylation_df[, c(1,2, (which(colnames(metabric_methylation_df)[2:ncol(metabric_methylation_df)] %fin% 
                                                                       intersecting_patient_ids) + 1))]
metabric_clinical_patient_df_sub <- metabric_clinical_patient_df[metabric_clinical_patient_df$PATIENT_ID %fin% 
                                                                   intersecting_patient_ids,]
metabric_clinical_sample_df_sub <- metabric_clinical_sample_df[metabric_clinical_sample_df$PATIENT_ID %fin% 
                                                                   intersecting_patient_ids,]


############################################################
### WRITE BACK
############################################################
write.table(metabric_clinical_patient_df_sub, paste0(main_path, "data_clinical_patient_sub.txt"))
write.table(metabric_clinical_sample_df_sub, paste0(main_path, "data_clinical_sample_sub.txt"))


############################################################
### COMBINE CLINICAL AND SAMPLE DFs INTO ONE
############################################################
patient_df <- read.csv(paste(output_path, "Patient and Sample DFs/patient_dataframe.csv", sep = ""), 
                       header = TRUE, row.names = 1)
# If BRCA, eliminate the gender column
is_brca <- TRUE
#is_brca <- FALSE
if(is_brca) {patient_df <- patient_df[, -which(colnames(patient_df) == "Gender")]}

sample_df <- read.csv(paste(output_path, "Patient and Sample DFs/sample_dataframe_cibersort_total_frac.csv", sep = ""), 
                      header = TRUE, row.names = 1)  # Total fraction

# Use the combine_patient_and_samp_dfs() function from additional_premodel_processing.R to combine
combined_pat_samp_df <- combine_patient_and_samp_dfs(patient_df, sample_df, "metabric")

# Remove rows that are entirely NA
combined_pat_samp_df <- combined_pat_samp_df[rowSums(is.na(combined_pat_samp_df)) != ncol(combined_pat_samp_df),]

# Write to CSV
write.csv(combined_pat_samp_df, paste0(output_path, "Patient and Sample DFs/combined_patient_sample_dataframe_cibersort_total_frac.csv"))


############################################################
### IMPORT ALL FILES BACK AND ADJUST HEADERS
############################################################

## MUTATION REGPROT DF ##
metabric_mutation_df_nonsynonymous <- read.csv(paste0(output_path, "Mutation/iprotein_results_nonsynonymous.csv"), 
                                           header = TRUE, check.names = FALSE)
metabric_mutation_df_silent <- read.csv(paste0(output_path, "Mutation/iprotein_results_silent.csv"), 
                                        header = TRUE, check.names = FALSE)

metabric_mutation_df_nonsynonymous$Swissprot <- unlist(lapply(metabric_mutation_df_nonsynonymous$Query, function(x) 
  unlist(strsplit(x, "|", fixed = TRUE))[2]))
metabric_mutation_df_silent$Swissprot <- unlist(lapply(metabric_mutation_df_silent$Query, function(x) 
  unlist(strsplit(x, "|", fixed = TRUE))[2]))

## MUTATION GENE TARG DF ##
metabric_mutation_count_df_nonsynonymous <- read.csv(paste0(output_path, "Mutation/Mutation Count Matrices/mut_count_matrix_nonsynonymous.csv"), 
                                                header = TRUE, check.names = FALSE)
metabric_mutation_count_df_silent <- read.csv(paste0(output_path, "Mutation/Mutation Count Matrices/mut_count_matrix_silent.csv"), 
                                                 header = TRUE, check.names = FALSE)
colnames(metabric_mutation_count_df_nonsynonymous)[1] <- 'Gene_Symbol'
colnames(metabric_mutation_count_df_silent)[1] <- 'Gene_Symbol'

#' Function to add patients with no mutation of the the given specificity in any protein
#' @param mut_count_matrix a maftools-produced mutation count matrix
#' @param intersecting_patients a list of patients that have all data types within the given
#' cancer type, before filtering is done
add_missing_patients_to_mut_count_mat <- function(mut_count_matrix, intersecting_patients) {
  missing_patients <- setdiff(intersecting_patients, colnames(mut_count_matrix))
  new_pat_df <- data.frame(matrix(nrow = nrow(mut_count_matrix), ncol = length(missing_patients)))
  new_pat_df[is.na(new_pat_df)] <- 0
  colnames(new_pat_df) <- missing_patients
  
  new_mut_count_mat <- cbind(mut_count_matrix, new_pat_df)
  return(new_mut_count_mat)
}

metabric_mutation_count_df_nonsynonymous <- add_missing_patients_to_mut_count_mat(metabric_mutation_count_df_nonsynonymous, intersecting_patient_ids)
metabric_mutation_count_df_silent <- add_missing_patients_to_mut_count_mat(metabric_mutation_count_df_silent, intersecting_patient_ids)

# If the gene symbols are currently rownames and not the first column
metabric_mutation_count_df_nonsynonymous$Swissprot <- unlist(lapply(metabric_mutation_count_df_nonsynonymous$Gene_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'uniprot_gn_id'])), collapse = ";")))
metabric_mutation_count_df_silent$Swissprot <- unlist(lapply(metabric_mutation_count_df_silent$Gene_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'uniprot_gn_id'])), collapse = ";")))


## CNA DF ##
metabric_cna_df <- read.csv(paste0(output_path, "CNV/CNA_DF_AllGenes.csv"), 
                                          header = TRUE, check.names = FALSE, row.names = 1)

## EXPRESSION DF ##
metabric_expression_df <- read.csv(paste0(output_path, "Expression/expression_log2intensity_filtBySD0.2_DF.csv"), 
                                          header = TRUE, check.names = FALSE)
colnames(metabric_expression_df)[1] <- 'ensg_id'


## METHYLATION DF ##
metabric_methylation_df_log2 <- read.csv(paste0(output_path, "Methylation/methylation_DF_RRBS.csv"), 
                                         header = TRUE, check.names = FALSE)
metabric_methylation_df_logit <- read.csv(paste0(output_path, "Methylation/methylation_DF_RRBS_logit.csv"), 
                                         header = TRUE, check.names = FALSE)
colnames(metabric_methylation_df_log2)[1] <- "Gene_Symbol"
colnames(metabric_methylation_df_logit)[1] <- "Gene_Symbol"

metabric_methylation_df_log2$ensg_ids <- unlist(lapply(metabric_methylation_df_log2$Gene_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'ensembl_gene_id'])), collapse = ";")))
metabric_methylation_df_logit$ensg_ids <- unlist(lapply(metabric_methylation_df_logit$Gene_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'ensembl_gene_id'])), collapse = ";")))

## PATIENT SAMPLE DF ## 
metabric_patient_sample_df <- read.csv(paste0(output_path, "Patient and Sample DFs/combined_patient_sample_dataframe_cibersort_total_frac.csv"), 
                                       header = TRUE, check.names = FALSE)
colnames(metabric_patient_sample_df)[1] <- "sample_id"

############################################################
### ENSURE ALL IDS STILL OVERLAP
############################################################
mutation_ids_nonsyn <- unique(unlist(lapply(metabric_mutation_df_nonsynonymous$Patient, function(x) 
  unlist(strsplit(x, ";", fixed = TRUE)))))  #819
mutation_ids_silent <- unique(unlist(lapply(metabric_mutation_df_silent$Patient, function(x) 
  unlist(strsplit(x, ";", fixed = TRUE)))))  #512

mutation_count_mat_ids_nonsyn <- colnames(metabric_mutation_count_df_nonsynonymous)[2:ncol(metabric_mutation_count_df_nonsynonymous)]
mutation_count_mat_ids_silent <- colnames(metabric_mutation_count_df_silent)[2:ncol(metabric_mutation_count_df_silent)]

cna_ids <- colnames(metabric_cna_df)

expression_ids <- colnames(metabric_expression_df)[2:ncol(metabric_expression_df)]

methylation_ids_log2 <- colnames(metabric_methylation_df_log2)[2:ncol(metabric_methylation_df_log2)]
methylation_ids_logit <- colnames(metabric_methylation_df_logit)[2:ncol(metabric_methylation_df_logit)]

patient_sample_ids <- metabric_patient_sample_df$sample_id


intersecting_patients <- intersect(cna_ids, 
                                   intersect(methylation_ids_log2,
                                             intersect(expression_ids,
                                                       intersect(mutation_count_mat_ids_missense, patient_sample_ids))))
print(length(intersecting_patients)) # 854 patients 

write.table(intersecting_patients, "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/cBioPortal/METABRIC/intersecting_ids.txt",
            quote = FALSE, row.names = FALSE)

# Read back
intersecting_patients <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/cBioPortal/METABRIC/intersecting_ids.txt", 
                                    header = T, check.names = F)[,1]



############################################################
### IF NECESSARY, SUBSET FINALIZED DATASETS BY INTERSECTING IDs
############################################################
metabric_mutation_count_df_nonsynonymous_sub <- metabric_mutation_count_df_nonsynonymous[, c(1, which(colnames(metabric_mutation_count_df_nonsynonymous) %in% 
                                                                                              intersecting_patients))]
metabric_mutation_count_df_silent_sub <- metabric_mutation_count_df_silent[, c(1, which(colnames(metabric_mutation_count_df_silent) %in% 
                                                                                                intersecting_patients))]

metabric_cna_df_sub <- metabric_cna_df[, c(1, which(colnames(metabric_cna_df) %in% intersecting_patients))]


metabric_expression_df_sub <- metabric_expression_df[, c(1, which(colnames(metabric_expression_df) %in% intersecting_patients))]

metabric_methylation_df_log2_sub <- metabric_methylation_df_log2[, c(1, which(colnames(metabric_methylation_df_log2) %in% intersecting_patients))]
metabric_methylation_df_logit_sub <- metabric_methylation_df_logit[, c(1, which(colnames(metabric_methylation_df_logit) %in% intersecting_patients))]

metabric_patient_sample_df_sub <- metabric_patient_sample_df[metabric_patient_sample_df$sample_id %in% intersecting_patients,]
  

# Do this for the labels on the regprot data frame as well
get_mutation_df_new_patient_labels <- function(metabric_mutation_df, intersecting_patients) {
  mutation_df_new_patient_labels <- lapply(metabric_mutation_df$Patient, function(x) {
    # Split apart the semicolon separated sample IDs
    spl_patients <- unlist(strsplit(x, ";", fixed = TRUE))
    # Check if this patient is in the intersecting patients list and return TRUE if so, F o.w.
    matching_patient_ind <- unlist(lapply(spl_patients, function(y) {
      if(y %fin% intersecting_patients) {
        return(TRUE)
      } else {return(FALSE)}
    }))
    # If there are overlapping patients in this row, then we want to keep this row
    if(TRUE %in% matching_patient_ind) {
      samps_to_keep <- spl_patients[matching_patient_ind]

      # If we still have samples, return them in this row 
      if(length(samps_to_keep) > 0) {
        return(paste(samps_to_keep, collapse = ";"))
      } else {return(NA)}
    } else {return(NA)}
  })
  return(mutation_df_new_patient_labels)
}

metabric_mutation_df_nonsynonymous_labels <- get_mutation_df_new_patient_labels(metabric_mutation_df_nonsynonymous, 
                                                                                 intersecting_patients)
metabric_mutation_df_silent_labels <- get_mutation_df_new_patient_labels(metabric_mutation_df_silent, 
                                                                                 intersecting_patients)

# Remove NA
metabric_mutation_df_nonsynonymous_sub <- metabric_mutation_df_nonsynonymous[unlist(lapply(metabric_mutation_df_nonsynonymous_labels, function(x) 
  ifelse(is.na(x), FALSE, TRUE))),]
metabric_mutation_df_silent_sub <- metabric_mutation_df_missense[unlist(lapply(metabric_mutation_df_silent_labels, function(x) 
  ifelse(is.na(x), FALSE, TRUE))),]

# Update the DF with the new intersecting patient labels
metabric_mutation_df_nonsynonymous_sub$Patient <- metabric_mutation_df_nonsynonymous_labels[!is.na(metabric_mutation_df_nonsynonymous_labels)]
metabric_mutation_df_nonsynonymous_sub$Patient <- unlist(lapply(metabric_mutation_df_nonsynonymous_sub$Patient, function(x) return(unlist(x))))

metabric_mutation_df_silent_sub$Patient <- metabric_mutation_df_silent_labels[!is.na(metabric_mutation_df_silent_labels)]
metabric_mutation_df_silent_sub$Patient <- unlist(lapply(metabric_mutation_df_silent_sub$Patient, function(x) return(unlist(x))))

############################################################
### WRITE THESE NEWLY SUBSETTED FILES BACK
############################################################
write.csv(metabric_mutation_count_df_nonsynonymous_sub, paste0(output_path, "Linear Model/mutation_count_df_nonsynonymous_IntersectPatients.csv"))
write.csv(metabric_mutation_count_df_silent_sub, paste0(output_path, "Linear Model/mutation_count_df_silent_IntersectPatients.csv"))

write.csv(metabric_cna_df_sub, paste0(output_path, "Linear Model/cna_df_allgenes_IntersectPatients.csv"))

write.csv(metabric_expression_df_sub, paste0(output_path, "Linear Model/expression_df_log2intensity_filtBySD0.2_IntersectPatients.csv"))

write.csv(metabric_methylation_df_log2_sub, paste0(output_path, "Linear Model/methylation_df_RRBS_log2_IntersectPatients.csv"))
write.csv(metabric_methylation_df_logit_sub, paste0(output_path, "Linear Model/methylation_df_RRBS_logit_IntersectPatients.csv"))

write.csv(metabric_patient_sample_df_sub, paste0(output_path, "Linear Model/patient_sample_df_cibersort_total_frac_IntersectPatients.csv"))

write.csv(metabric_mutation_df_missense_sub, paste0(output_path, "Linear Model/iprotein_mutation_df_missense.csv"))
write.csv(metabric_mutation_df_misAndNon_sub, paste0(output_path, "Linear Model/iprotein_mutation_df_missense_nonsense.csv"))
write.csv(metabric_mutation_df_nonsynonymous_sub, paste0(output_path, "Linear Model/iprotein_mutation_df_nonsynonymous.csv"))
write.csv(metabric_mutation_df_silent_sub, paste0(output_path, "Linear Model/iprotein_mutation_df_silent.csv"))


# Read back
metabric_mutation_df_nonsynonymous_sub <- read.csv(paste0(output_path, "Linear Model/mutation_count_df_nonsynonymous_IntersectPatients.csv"),
                                                   header = T, check.names = F)

#############################################################
### GENERATE PROTEIN INPUT TABLES (DRIVERS)
#############################################################
# Import DF with known driver genes
driver_gene_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/GRCh38_driver_gene_list.tsv", 
                           sep = "\t", header = TRUE, comment.char = "#", skip = 10)

iprotein_prots_nonsyn <- unique(metabric_mutation_df_nonsynonymous_sub$Swissprot)  # 122 unique mutated proteins
iprotein_prots_silent <- unique(metabric_mutation_df_silent_sub$Swissprot)  # 122 unique mutated proteins

# Proceed with nonsynonymous mutations; use functions from 'generate_protein_input_tables.R'
#iprotein_prots_ensg <- convert_to_ensembl(iprotein_prots_nonsyn)
# Create initial DF with all proteins
#recombine_into_df_and_write(iprotein_prots, iprotein_prots_ensg, "iprotein", prot_path)

# Call filter to include only proteins mutated in at least 5 patients, and at least 5% of patients
# Nonsynonymous
regprot_input_df_all_5perc <- create_regulatory_prot_input_df(as.data.table(metabric_mutation_df_nonsynonymous_sub), 
                                                              0.05, 5, NA, intersecting_patients, "i-protein", intersecting_patients, 
                                                              NA, all_genes_id_conv, NA, NA, "METABRIC")
regprot_input_df_driver_5perc <- create_regulatory_prot_input_df(as.data.table(metabric_mutation_df_nonsynonymous_sub), 
                                                                 0.05, 5, NA, intersecting_patients, "i-protein", intersecting_patients, 
                                                                 driver_gene_df, all_genes_id_conv, NA, NA, "METABRIC")

# Print out the names of the drivers that were written to a given file
print(paste(unique(unlist(lapply(regprot_input_df_driver_5perc$swissprot_ids, function(id) 
  all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == id, 'external_gene_name']))), collapse = ", "))


write.csv(regprot_input_df_all_5perc, paste0(output_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.05Freq_all_nonsynonymous.csv"))
write.csv(regprot_input_df_driver_5perc, paste0(output_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.05Freq_drivers_nonsynonymous.csv"))


# Repeat this for individual subtypes
patient_set_lumA <- intersect(intersecting_patients, metabric_clinical_patient_df[metabric_clinical_patient_df$CLAUDIN_SUBTYPE == "LumA", 'PATIENT_ID'])
patient_set_lumB <- intersect(intersecting_patients, metabric_clinical_patient_df[metabric_clinical_patient_df$CLAUDIN_SUBTYPE == "LumB", 'PATIENT_ID'])
patient_set_basal <- intersect(intersecting_patients, metabric_clinical_patient_df[metabric_clinical_patient_df$CLAUDIN_SUBTYPE %in% c('claudin-low', 'Basal'), 'PATIENT_ID'])
patient_set_her2 <- intersect(intersecting_patients, metabric_clinical_patient_df[metabric_clinical_patient_df$CLAUDIN_SUBTYPE == 'Her2', 'PATIENT_ID'])


#' Function to keep only patients in the regprot mutation DF that are a member of the
#' given patient set
#' @param metabric_mutation_df a regprot mutation DF with METABRIC IDs
#' @param patient_set a vector of METABRIC patient IDs to keep
get_regprot_df_subsetted_by_subtype <- function(metabric_mutation_df, patient_set) {
  new_rows <- lapply(1:nrow(metabric_mutation_df), function(i) {
    pats <- unlist(strsplit(metabric_mutation_df[i, 'Patient'], ";", fixed = TRUE))
    pats_in_set <- paste(intersect(pats, patient_set), collapse = ";")
    
    row <- as.data.frame(metabric_mutation_df[i,])
    colnames(row) <- colnames(metabric_mutation_df)
    row[,'Patient'] <- pats_in_set
    return(row)
  })
  new_df <- do.call("rbind", new_rows)
  return(new_df)
}

metabric_mutation_df_nonsynonymous_lumA <- get_regprot_df_subsetted_by_subtype(metabric_mutation_df_nonsynonymous,
                                                                               patient_set_lumA)
metabric_mutation_df_nonsynonymous_lumB <- get_regprot_df_subsetted_by_subtype(metabric_mutation_df_nonsynonymous,
                                                                               patient_set_lumB)
metabric_mutation_df_nonsynonymous_basal <- get_regprot_df_subsetted_by_subtype(metabric_mutation_df_nonsynonymous,
                                                                               patient_set_basal)
metabric_mutation_df_nonsynonymous_her2 <- get_regprot_df_subsetted_by_subtype(metabric_mutation_df_nonsynonymous,
                                                                               patient_set_her2)

regprot_input_df_all_5perc_lumA <- create_regulatory_prot_input_df(metabric_mutation_df_nonsynonymous_lumA, 
                                                              0.05, 5, "i-protein", patient_set_lumA, 
                                                              NA, all_genes_id_conv, NA, NA, "METABRIC")
regprot_input_df_driver_5perc_lumA <- create_regulatory_prot_input_df(metabric_mutation_df_nonsynonymous_lumA, 
                                                                 0.05, 5, "i-protein", patient_set_lumA, 
                                                                 driver_gene_df, all_genes_id_conv, NA, NA, "METABRIC")

write.csv(regprot_input_df_all_5perc, paste0(output_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.05Freq_all_nonsynonymous.csv"))
write.csv(regprot_input_df_driver_5perc, paste0(output_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.05Freq_drivers_nonsynonymous.csv"))
