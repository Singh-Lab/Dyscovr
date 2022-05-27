############################################################
### Process METABRIC Data
### Written By: Sara Geraghty, April 2022
############################################################

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
#metabric_MN_mutation_df <- metabric_mutation_df[(metabric_mutation_df$Variant_Classification == "Missense_Mutation") |
                                                  #(metabric_mutation_df$Variant_Classification == "Nonsense_Mutation"),]
# Use this file in the "process_mutation_data.R" pipeline


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
mutation_df_patient_ids <- metabric_MN_mutation_df$Tumor_Sample_Barcode
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
metabric_MN_mutation_df_sub <- metabric_MN_mutation_df[metabric_MN_mutation_df$Tumor_Sample_Barcode %fin% 
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

sample_df <- read.csv(paste(main_path, "Patient and Sample DFs/sample_dataframe_cibersort_total_frac.csv", sep = ""), 
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
metabric_mutation_df_missense <- read.csv(paste0(output_path, "Mutation/iprotein_results_missense.csv"), 
                                          header = TRUE, check.names = FALSE)
metabric_mutation_df_misAndNon <- read.csv(paste0(output_path, "Mutation/iprotein_results_missense_nonsense.csv"), 
                                          header = TRUE, check.names = FALSE)
metabric_mutation_df_missense$Swissprot <- unlist(lapply(metabric_mutation_df_missense$Query, function(x) 
  unlist(strsplit(x, "|", fixed = TRUE))[2]))
metabric_mutation_df_misAndNon$Swissprot <- unlist(lapply(metabric_mutation_df_misAndNon$Query, function(x) 
  unlist(strsplit(x, "|", fixed = TRUE))[2]))


metabric_mutation_count_df_missense <- read.csv(paste0(output_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense.csv"), 
                                       header = TRUE, check.names = FALSE)
metabric_mutation_count_df_misAndNon <- read.csv(paste0(output_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense_nonsense.csv"), 
                                                header = TRUE, check.names = FALSE)
colnames(metabric_mutation_count_df_missense)[1] <- 'Gene_Symbol'
colnames(metabric_mutation_count_df_misAndNon)[1] <- 'Gene_Symbol'


metabric_cna_df <- read.csv(paste0(output_path, "CNV/CNA_DF_AllGenes.csv"), 
                                          header = TRUE, check.names = FALSE, row.names = 1)

metabric_expression_df <- read.csv(paste0(output_path, "Expression/expression_log2intensity_filtBySD0.2_DF.csv"), 
                                          header = TRUE, check.names = FALSE)
colnames(metabric_expression_df)[1] <- 'ensg_id'


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


metabric_patient_sample_df <- read.csv(paste0(output_path, "Patient and Sample DFs/combined_patient_sample_dataframe_cibersort_total_frac.csv"), 
                                       header = TRUE, check.names = FALSE)
colnames(metabric_patient_sample_df)[1] <- "sample_id"

############################################################
### ENSURE ALL IDS STILL OVERLAP
############################################################
mutation_ids_missense <- unique(unlist(lapply(metabric_mutation_df_missense$Patient, function(x) 
  unlist(strsplit(x, ";", fixed = TRUE)))))  #826
mutation_ids_misAndNon <- unique(unlist(lapply(metabric_mutation_df_misAndNon$Patient, function(x) 
  unlist(strsplit(x, ";", fixed = TRUE)))))  #834

mutation_count_mat_ids_missense <- colnames(metabric_mutation_count_df_missense)[2:ncol(metabric_mutation_count_df_missense)]
mutation_count_mat_ids_misAndNon <- colnames(metabric_mutation_count_df_misAndNon)[2:ncol(metabric_mutation_count_df_misAndNon)]

cna_ids <- colnames(metabric_cna_df)

expression_ids <- colnames(metabric_expression_df)[2:ncol(metabric_expression_df)]

methylation_ids_log2 <- colnames(metabric_methylation_df_log2)[2:ncol(metabric_methylation_df_log2)]
methylation_ids_logit <- colnames(metabric_methylation_df_logit)[2:ncol(metabric_methylation_df_logit)]

patient_sample_ids <- metabric_patient_sample_df$sample_id


intersecting_patients <- intersect(cna_ids, 
                                   intersect(methylation_ids_log2,
                                             intersect(expression_ids,
                                                       intersect(mutation_count_mat_ids_missense, patient_sample_ids))))
print(length(intersecting_patients)) # 826 patients

write.table(intersecting_patients, "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/cBioPortal/METABRIC/intersecting_ids.txt",
            quote = FALSE, row.names = FALSE)


############################################################
### IF NECESSARY, SUBSET FINALIZED DATASETS BY INTERSECTING IDs
############################################################

metabric_mutation_count_df_missense_sub <- metabric_mutation_count_df_missense[, c(1, which(colnames(metabric_mutation_count_df_missense) %in% 
                                                                                         intersecting_patients))]
metabric_mutation_count_df_misAndNon_sub <- metabric_mutation_count_df_misAndNon[, c(1, which(colnames(metabric_mutation_count_df_misAndNon) %in% 
                                                                                           intersecting_patients))]

metabric_cna_df_sub <- metabric_cna_df[, c(1, which(colnames(metabric_cna_df) %in% intersecting_patients))]


metabric_expression_df_sub <- metabric_expression_df[, c(1, which(colnames(metabric_expression_df) %in% intersecting_patients))]

metabric_methylation_df_log2_sub <- metabric_methylation_df_log2[, c(1, which(colnames(metabric_methylation_df_log2) %in% intersecting_patients))]
metabric_methylation_df_logit_sub <- metabric_methylation_df_logit[, c(1, which(colnames(metabric_methylation_df_logit) %in% intersecting_patients))]

metabric_patient_sample_df_sub <- metabric_patient_sample_df[metabric_patient_sample_df$sample_id %in% intersecting_patients,]
  

# Do this for the labels on the regprot data frame as well
mutation_df_new_patient_labels <- lapply(metabric_mutation_df_missense$Patient, function(x) {
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

# Remove NA
metabric_mutation_df_missense_sub <- metabric_mutation_df_missense[unlist(lapply(mutation_df_new_patient_labels, function(x) 
  ifelse(is.na(x), FALSE, TRUE))),]
# Update the DF with the new intersecting patient labels
metabric_mutation_df_missense_sub$Patient <- mutation_df_new_patient_labels[!is.na(mutation_df_new_patient_labels)]
metabric_mutation_df_missense_sub$Patient <- unlist(lapply(metabric_mutation_df_missense_sub$Patient, function(x) return(unlist(x))))


############################################################
### WRITE THESE NEWLY SUBSETTED FILES BACK
############################################################
write.csv(metabric_mutation_count_df_missense_sub, paste0(output_path, "Linear Model/mutation_count_df_missense_IntersectPatients.csv"))
write.csv(metabric_mutation_count_df_misAndNon_sub, paste0(output_path, "Linear Model/mutation_count_df_missense_nonsense_IntersectPatients.csv"))
write.csv(metabric_cna_df_sub, paste0(output_path, "Linear Model/cna_df_allgenes_IntersectPatients.csv"))
write.csv(metabric_expression_df_sub, paste0(output_path, "Linear Model/expression_df_log2intensity_filtBySD0.2_IntersectPatients.csv"))
write.csv(metabric_methylation_df_log2_sub, paste0(output_path, "Linear Model/methylation_df_RRBS_log2_IntersectPatients.csv"))
write.csv(metabric_methylation_df_logit_sub, paste0(output_path, "Linear Model/methylation_df_RRBS_logit_IntersectPatients.csv"))
write.csv(metabric_patient_sample_df_sub, paste0(output_path, "Linear Model/patient_sample_df_cibersort_total_frac_IntersectPatients.csv"))
write.csv(metabric_mutation_df_missense_sub, paste0(output_path, "Linear Model/iprotein_mutation_df_missense.csv"))
write.csv(metabric_mutation_df_misAndNon_sub, paste0(output_path, "Linear Model/iprotein_mutation_df_missense_nonsense.csv"))




