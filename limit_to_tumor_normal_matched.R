############################################################
### Filter to Tumor-Normal Matched
### Written By: Sara Geraghty, March 2021
############################################################

# This file takes all linear model input files and limits them
# to only patients with tumor-normal matched data.

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


############################################################
# IMPORT EXPRESSION FILES
############################################################
# 1. Filtered FPKM
expression_df_fpkm <- fread(paste(main_path, "Expression/expression_fpkm_filt_DF.csv", sep = ""), 
                       header = TRUE)  # fpkm-normalized and filtered

# 2. Filtered TMM 
expression_df_tmm <- fread(paste(main_path, "Expression/tmm_normalized_expression_counts.csv", sep = ""), 
                       header = TRUE)

# 3. Quantile-Normalized
expression_df_qn <- fread(paste(main_path, "Expression/expression_quantile_norm_DF.csv", sep = ""),
                       header = TRUE)

# 4. Rank-Normalized
expression_df_rn <- fread(paste(main_path, "Expression/expression_rank_norm_DF.csv", sep = ""),
                       header = TRUE)

# Are there differences in the samples of each of these? Do we need to do them 
# separately, or can we use the tumor-normal matched patients from one and apply to the others?
diff1 <- setdiff(colnames(expression_df_fpkm), colnames(expression_df_tmm))
diff2 <- setdiff(colnames(expression_df_fpkm), colnames(expression_df_qn))
diff3 <- setdiff(colnames(expression_df_fpkm), colnames(expression_df_rn))

print(length(diff1) > 0) # FALSE
print(length(diff2) > 0) # FALSE 
print(length(diff3) > 0) # FALSE


############################################################
# GET PATIENTS WITH TUMOR-NORMAL MATCHED EXPRESSION
############################################################
# Get all patients that have a normal expression file (typically the limiting factor)
patients_with_matched <- unlist(lapply(colnames(expression_df)[grepl("-1", colnames(expression_df))], 
                                       function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))

# Make sure no patients *just* have normal sequenced data and no cancer; eliminate if so
exp_colnames_patients_only <- unlist(lapply(colnames(expression_df), function(x) 
  unlist(strsplit(x, "-", fixed = TRUE))[3]))

# Get patients with only a normal column
#patients_w_only_norm <- exp_colnames_patients_only[ave(exp_colnames_patients_only, 
                                                       #exp_colnames_patients_only, FUN = length) == 1]
patients_w_only_norm <- patients_with_matched[which(unlist(lapply(patients_with_matched, function(pat) {
  ifelse(exp_colnames_patients_only[exp_colnames_patients_only == pat] == 1, TRUE, FALSE)})))]

if(length(patients_w_only_norm) > 0) {
  # Remove these patients from the list
  patients_with_matched <- setdiff(patients_with_matched, patients_w_only_norm)
}
# 87 in BRCA

############################################################
# IMPORT METHYLATION FILES
############################################################
# 1. Methylation status in cancer as a covariate
# Raw
methylation_df_beta <- fread(paste(main_path, "Methylation/ methylation_DF_Beta.csv", sep = ""), 
                        header = TRUE)
methylation_df_m <- fread(paste(main_path, "Methylation/ methylation_DF_M.csv", sep = ""), 
                        header = TRUE)

# Thresholded
methylation_df_thres_0.8 <- fread(paste(main_path, "Methylation/methylation_DF_Beta_0.8.csv", sep = ""), 
                        header = TRUE)
methylation_df_thres_0.5 <- fread(paste(main_path, "Methylation/methylation_DF_Beta_0.5.csv", sep = ""), 
                        header = TRUE)

# Bucketed
methylation_df_bucket_Beta <- fread(paste(main_path, "Methylation/methylation_DF_bucketed_Beta.csv", sep = ""), 
                        header = TRUE)
methylation_df_bucket_M <- fread(paste(main_path, "Methylation/methylation_DF_bucketed_M.csv", sep = ""), 
                        header = TRUE)

# Are there differences in the samples of each of these? Do we need to do them 
# separately, or can we use the tumor-normal matched patients from one and apply to the others?
diff1 <- setdiff(colnames(methylation_df_beta), colnames(methylation_df_m))
diff2 <- setdiff(colnames(methylation_df_beta), colnames(methylation_df_thres_0.8))
diff3 <- setdiff(colnames(methylation_df_beta), colnames(methylation_df_thres_0.5))
diff4 <- setdiff(colnames(methylation_df_beta), colnames(methylation_df_bucket_Beta))
diff5 <- setdiff(colnames(methylation_df_beta), colnames(methylation_df_bucket_M))

print(length(diff1) > 0) # FALSE
print(length(diff2) > 0) # FALSE 
print(length(diff3) > 0) # FALSE
print(length(diff4) > 0) # FALSE 
print(length(diff5) > 0) # FALSE

############################################################
# GET PATIENTS WITH TUMOR-NORMAL MATCHED METHYLATION
############################################################
# Get non-matches from the methylation file and eliminate them as well
methyl_patients <- unlist(lapply(colnames(methylation_df), function(x) {
  if (grepl("-0", x)) {
    return(unlist(strsplit(x, "-", fixed = TRUE))[1])
  } else {return(NA)}
}))
methyl_patients <- methyl_patients[!is.na(methyl_patients)]

# Get patients that have at least a cancer sample (for eQTL)
patients_with_matched <- intersect(patients_with_matched, methyl_patients)  # 82


# Get patients that have a paired cancer and normal sample (for meQTL)
# Make sure no patients *just* have normal sequenced data and no cancer; eliminate if so
meth_colnames_patients_only <- unlist(lapply(colnames(methylation_df), function(x) 
  unlist(strsplit(x, "-", fixed = TRUE))[1]))

patients_w_only_norm <- patients_with_matched[which(unlist(lapply(patients_with_matched, function(pat) {
  ifelse(meth_colnames_patients_only[meth_colnames_patients_only == pat] == 1, TRUE, FALSE)})))]

if(length(patients_w_only_norm) > 0) {
  # Remove these patients from the list
  patients_with_matched_meQTL <- setdiff(patients_with_matched, patients_w_only_norm)
}


############################################################
# IMPORT CNA FILES
############################################################
# cna_df <- read.csv(paste(main_path, "CNV/Gene-level Score/CNV_DF_AllGenes.csv", sep = ""), 
# header = TRUE, row.names = 1)
# Raw CNA Values
cna_df_raw <- read.csv(paste(main_path, "CNV/Gene-level Raw/CNV_DF_AllGenes.csv", sep = ""), 
                   header = TRUE, row.names = 1, check.names = FALSE)
# Bucketed CNA Values (Excl. Amp)
cna_df_bucket_exclAmp <- read.csv(paste(main_path, "CNV/Gene-level Raw/CNV_DF_bucketed_exclAmp_AllGenes.csv", sep = ""), 
                   header = TRUE, check.names = FALSE)

# Bucketed CNA Values (Incl. Amp)
cna_df_bucket_inclAmp <- read.csv(paste(main_path, "CNV/Gene-level Raw/CNV_DF_bucketed_inclAmp_AllGenes.csv", sep = ""), 
                   header = TRUE, check.names = FALSE)


# Are there differences in the samples of each of these? Do we need to do them 
# separately, or can we use the tumor-normal matched patients from one and apply to the others?
diff1 <- setdiff(colnames(cna_df_raw), colnames(cna_df_bucket_exclAmp))
diff2 <- setdiff(colnames(cna_df_raw), colnames(cna_df_bucket_inclAmp))

print(length(diff1) > 0) # FALSE
print(length(diff2) > 0) # FALSE


# Get the patients that have CNA files & expression/ methylation
cna_patients <- unlist(lapply(colnames(cna_df), function(x) 
  unlist(strsplit(x, "-", fixed = TRUE))[1]))
patients_with_matched <- intersect(patients_with_matched, cna_patients)  # 81
if(exists(patients_with_matched_meQTL)) {
  patients_with_matched_meQTL <- intersect(patients_with_matched_meQTL, cna_patients)  
}


############################################################
# MUTATION TARGET DATA FRAME
############################################################
# Mutation count matrix, for checking if the downstream target t_k has a mutation in a given patient
mutation_targ_df <- read.csv(paste(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense.csv", sep = ""), 
                             header = TRUE, row.names = 1, check.names = FALSE)

mut_patients <- unlist(lapply(colnames(mutation_targ_df), function(x) 
  unlist(strsplit(x, "-", fixed = TRUE))[3]))

# Add a column of Swissprot IDs
mutation_targ_df$swissprot <- unlist(lapply(rownames(mutation_targ_df), function(x) 
  paste(unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'uniprot_gn_id']), collapse = ";")))

patients_with_matched <- intersect(patients_with_matched, mut_patients) # 79
if(exists(patients_with_matched_meQTL)) {
  patients_with_matched_meQTL <- intersect(patients_with_matched_meQTL, mut_patients) 
}

############################################################
# REGULATORY PROTEIN MUTATION DF
############################################################
# The output results file, for checking which patients have a mutation at the given level of 
# specificity for a particular regulatory protein r_i
mutation_regprot_df <- read.csv(paste(main_path, "Mutation/iprotein_results_missense.csv", sep = ""))   # I-Protein, all ligands
mutation_regprot_df <- read.csv(paste(main_path, "Mutation/iprotein_results_missense_nucacids.csv", sep = ""))   # I-Protein, DNA/RNA-binding only
# For I-Protein, add a Swissprot column
mutation_regprot_df$Swissprot <- unlist(lapply(mutation_regprot_df$Query, function(x) 
  unlist(strsplit(x, "|", fixed = TRUE))[2]))
mutation_regprot_df <- read.csv(paste(main_path, "Mutation/idomain_results_missense.csv", sep = ""))   # I-Domain, all ligands
mutation_regprot_df <- read.csv(paste(main_path, "Mutation/idomain_results_missense_nucacids.csv", sep = ""))   # I-Domain, DNA/RNA-binding only
mutation_regprot_df <- read.csv(paste(main_path, "Mutation/ibindingpos_results_missense.csv", sep = ""))   # I-Binding Position, all ligands
mutation_regprot_df <- read.csv(paste(main_path, "Mutation/ibindingpos_results_missense_nucacids.csv", sep = ""))   # I-Binding Position, DNA/RNA-binding only


#' Function that gets the patient overlap from a mutation regulatory protein DF 
#' and the current list of patients that have matched tumor-normal samples in
#' all other data types so far
#' @param mutation_regprot_df a mutation regulatory protein DF from the 
#' process_mutation_data.R file
#' @param patients_with_matched a vector of current patients that have tumor-
#' normal matched data for the other cancer types
get_regprot_patients_w_matched <- function(mutation_regprot_df, patients_with_matched) {
  patients_with_matched <- unlist(lapply(mutation_regprot_df$Patient, function(x) {
    patients <- unlist(strsplit(x, ";", fixed = TRUE))
    patients <- unlist(lapply(patients, function(x) unlist(strsplit(x, "-", fixed = TRUE))[1]))
    return(intersect(patients, patients_with_matched))
  }))
}

patients_with_matched <- unique(get_regprot_patients_w_matched(mutation_regprot_df, 
                                                               patients_with_matched))  # 73
patients_with_matched_meQTL <- unique(get_regprot_patients_w_matched(mutation_regprot_df, 
                                                                     patients_with_matched_meQTL))  


############################################################
# IMPORT PATIENT/SAMPLE FILES
############################################################
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_abs_fpkm_tnm.csv", sep = ""), 
                       header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_abs_fpkm_top10k_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_abs_tmm_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_abs_qn_top10k_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_abs_rn_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_abs_rn_top10k_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)

patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_tmm_tnm.csv", sep = ""), 
                       header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_fpkm_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_fpkm_top10k_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_qn_top10k_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_rn_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_cibersort_total_frac_rn_top10k_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)


patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_timer_tmm_tnm.csv", sep = ""), 
                       header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_timer_fpkm_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_timer_fpkm_top10k_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_timer_qn_top10k_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_timer_rn_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)
patient_sample_df <- read.csv(paste(main_path, "Linear Model/Patient and Sample DFs/combined_patient_sample_DF_timer_rn_top10k_tnm.csv", sep = ""), 
                              header = TRUE, row.names = 1, check.names = FALSE)

patients_clin <- unlist(lapply(rownames(patient_sample_df), function(x) {
  unlist(strsplit(x, "-", fixed = TRUE))[3]}))
patients_with_matched <- intersect(patients_with_matched, patients_clin)  # 73
patients_with_matched_meQTL <- intersect(patients_with_matched_meQTL, patients_clin) 


############################################################
# LIMIT ALL FILES TO JUST PATIENTS WITH TUMOR-NORMAL MATCHED 
# EXPRESSION & METHYLATION DATA: eQTL
############################################################
#' Limit all files to tumor-normal matched, given a set of patients that have matched
#' (for either the eQTL or meQTL case)
#' @param patients_with_matched a vector of patients that have T-N matched data for
#' all data types
#' @param expression_df the expression data frame to limit to only these patients
#' @param methylation_df the methylation data frame to limit to only these patients
#' @param cna_df the CNA data frame to limit to only these patients
#' @param mutation_targ_df the mutation gene target DF to limit to only these patients
#' @param mutation_regprot_df the mutation regulatory protein DF to limit to only 
#' these patients
#' @param patient_sample_df the DF with patient and sample clinical information to
#' limit to only these patients
limit_files_to_TN_matched <- function(patients_with_matched, expression_df, methylation_df,
                                      cna_df, mutation_targ_df, mutation_regprot_df, 
                                      patient_sample_df) {
  # Expression DF
  expression_df_cols_to_keep <- unlist(lapply(colnames(expression_df), function(x) 
    ifelse(unlist(strsplit(x, "-", fixed = TRUE))[3] %fin% patients_with_matched, TRUE, FALSE)))
  expression_df <- expression_df[, expression_df_cols_to_keep]
  
  # Methylation DF
  methylation_df_cols_to_keep <- unlist(lapply(colnames(methylation_df), function(x) 
    ifelse(unlist(strsplit(x, "-", fixed = TRUE))[1] %fin% patients_with_matched, TRUE, FALSE)))
  methylation_df <- methylation_df[, methylation_df_cols_to_keep]
  ensg_ids <- unlist(lapply(rownames(methylation_df), function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'ensembl_gene_id']), collapse = ";")))
  methylation_df <- cbind(methylation_df, ensg_ids)
  #differential_methylation_df <- differential_methylation_df[,colnames(differential_methylation_df) %fin% patients_with_matched]
  #ensg_ids <- unlist(lapply(rownames(differential_methylation_df), function(x) paste(unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'ensembl_gene_id']), collapse = ";")))
  #differential_methylation_df <- cbind(differential_methylation_df, ensg_ids)
  
  # CNA DF
  cna_df_cols_to_keep <- unlist(lapply(colnames(cna_df), function(x) 
    ifelse(unlist(strsplit(x, "-", fixed = TRUE))[1] %fin% patients_with_matched, TRUE, FALSE)))
  cna_df <- cna_df[, cna_df_cols_to_keep]
  
  # Mutation DFs
  mutation_targ_cols_to_keep <- unlist(lapply(colnames(mutation_targ_df), function(x) 
    ifelse(unlist(strsplit(x, "-", fixed = TRUE))[3] %fin% patients_with_matched, TRUE, FALSE)))
  mutation_targ_df <- mutation_targ_df[, mutation_targ_cols_to_keep]  
  
  mutation_regprot_df <- limit_mutation_to_tumNormMatched(mutation_regprot_df, 
                                                          patients_with_matched, 
                                                          "I-Protein")
  # Patient DF
  patient_sample_rows_to_keep <- unlist(lapply(rownames(patient_sample_df), function(x)
    ifelse(unlist(strsplit(x, "-", fixed = TRUE))[3] %fin% patients_with_matched, TRUE, FALSE)))
  patient_sample_df <- patient_sample_df[patient_sample_rows_to_keep,]
  
  # Return all of these as a list
  return(list("Expression_DF" = expression_df, "Methylation_DF" = methylation_df,
              "CNA_DF" = cna_df, "Mutation_Targ_DF" = mutation_targ_df,
              "Mutation_Regprot_DF" = mutation_regprot_df, "Patient_Samp_DF" = patient_sample_df))
}

# Call this function
eQTL_tables <- limit_files_to_TN_matched(patients_with_matched, expression_df,
                                         methylation_df, cna_df, mutation_targ_df,
                                         mutation_regprot_df, patient_sample_df)
expression_df <- eQTL_tables[['Expression_DF']]
methylation_df <- eQTL_tables[['Methylation_DF']]
cna_df <- eQTL_tables[['CNA_DF']]
mutation_targ_df <- eQTL_tables[['Mutation_Targ_DF']]
mutation_regprot_df <- eQTL_tables[['Mutation_Regprot_DF']]
patient_sample_df <- eQTL_tables[['Patient_Samp_DF']]


meQTL_tables <- limit_files_to_TN_matched(patients_with_matched_meQTL, expression_df,
                                          methylation_df, cna_df, mutation_targ_df,
                                          mutation_regprot_df, patient_sample_df)
expression_df <- meQTL_tables[['Expression_DF']]
methylation_df <- meQTL_tables[['Methylation_DF']]
cna_df <- meQTL_tables[['CNA_DF']]
mutation_targ_df <- meQTL_tables[['Mutation_Targ_DF']]
mutation_regprot_df <- meQTL_tables[['Mutation_Regprot_DF']]
patient_sample_df <- meQTL_tables[['Patient_Samp_DF']]


#' Helper function for filtering regulatory protein mutation file (based on specificity) 
#' to only tumor-normal matched patients
#' @param mutation_regprot_df the regulatory protein mutation file to be filtered
#' @param specificity the specificity in question (either "I-Protein", "I-Domain", or "I-Binding Position")
limit_mutation_to_tumNormMatched <- function(mutation_regprot_df, patients_with_matched, specificity) {
  if (specificity == "I-Domain" | specificity == "I-Binding Position") {
    mutation_regprot_df <- mutation_regprot_df[mutation_regprot_df$Patient %fin% patients_with_matched,]
  } else if (specificity == "I-Protein") {
    mut_reg_df_new_rows <- lapply(1:nrow(mutation_regprot_df), function(i) {
      patients <- unlist(strsplit(mutation_regprot_df$Patient[i], ";", fixed = TRUE))
      patients <- unlist(lapply(patients, function(x) 
        unlist(strsplit(x, "-", fixed = TRUE))[1]))
      overlap <- intersect(patients, patients_with_matched)
      
      if (length(overlap) > 0) {
        # If there is overlap, keep only the intersecting patients in the label
        row <- mutation_regprot_df[i,]
        row$Patient <- paste(overlap, collapse = ";")
        return(row)
      
      } else {
        # If there is no overlap, we don't want to keep this row (no T-N matched patients)
        return(NA)
      }
    })
    mut_reg_df_new_rows <- mut_reg_df_new_rows[!is.na(mut_reg_df_new_rows)]
    mutation_regprot_df_new <- do.call(rbind, mut_reg_df_new_rows)
    return(mutation_regprot_df_new)
  } else {
    print(paste("Unknown specificity:", specificity))
  }
}


############################################################
# WRITE & READ TUMOR-NORMAL MATCHED FILES
############################################################
# Write these to files 

# EXPRESSION
write.csv(expression_df, paste(main_path, "Expression/TNM/expression_fpkm_DF_tumNormMatched.csv", sep = ""))
# write.csv(expression_df, paste(main_path, "Expression/TNM/expression_fpkm_filt_DF_tumNormMatched.csv", sep = ""))
# write.csv(expression_df, paste(main_path, "Expression/TNM/expression_tmm_DF_tumNormMatched.csv", sep = ""))
write.csv(expression_df, paste(main_path, "Expression/TNM/expression_fpkm_DF_tumNormMatched_meQTL.csv", sep = ""))
# write.csv(expression_df, paste(main_path, "Expression/TNM/expression_fpkm_filt_DF_tumNormMatched_meQTL.csv", sep = ""))
# write.csv(expression_df, paste(main_path, "Expression/TNM/expression_tmm_DF_tumNormMatched_meQTL.csv", sep = ""))

# METHYLATION
write.csv(methylation_df, paste(main_path, "Methylation/TNM/methylation_DF_0.8_tumNormMatched.csv", sep = ""))
write.csv(differential_methylation_df, paste(main_path, "Methylation/TNM/methylation_DF_0.8_tumNormMatched.csv", sep = ""))
write.csv(methylation_df_dependent, paste(main_path, "Methylation/TNM/methylation_DF_0.8_tumNormMatched_dependent.csv", sep = ""))

# CNA
write.csv(cna_df, paste(main_path, "CNV/Gene-level Raw/TNM/CNV_DF_AllGenes_CancerOnly_tumNormMatched.csv", sep = ""))
write.csv(cna_df, paste(main_path, "CNV/Gene-level Raw/TNM/CNV_DF_AllGenes_CancerOnly_tumNormMatched_meQTL.csv", sep = ""))

# MUTATION
write.csv(mutation_targ_df, paste(main_path, "Mutation/Mutation Count Matrices/TNM/mut_count_matrix_missense_tumNormMatched.csv", sep = ""))
write.csv(mutation_targ_df, paste(main_path, "Mutation/Mutation Count Matrices/TNM/mut_count_matrix_missense_tumNormMatched_meQTL.csv", sep = ""))
write.csv(mutation_regprot_df, paste(main_path, "Mutation/TNM/iprotein_results_missense_tumNormMatched.csv", sep = ""))
write.csv(mutation_regprot_df, paste(main_path, "Mutation/TNM/iprotein_results_missense_tumNormMatched_meQTL.csv", sep = ""))

# PATIENT/ SAMPLE
write.csv(patient_sample_df, paste(main_path, "Linear Model/Patient and Sample DFs/TNM/combined_patient_sample_tmm_totalFracICD_rawCNA_tumNormMatched.csv", sep = ""))
write.csv(patient_sample_df, paste(main_path, "Linear Model/Patient and Sample DFs/TNM/combined_patient_sample_tmm_totalFracICD_rawCNA_tumNormMatched_meQTL.csv", sep = ""))

# LIST OF PATIENTS WITH TUMOR-NORMAL MATCHED EXPRESSION & METHYLATION DATA
write(patients_with_matched, paste(main_path, "patients_with_tumNormMatched_eQTL.txt", sep = ""))
write(patients_with_matched_meQTL, paste(main_path, "patients_with_tumNormMatched_meQTL.txt", sep = ""))


#
#
#


# Read them back
# Unfiltered FPKM
#expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_DF_tumNormMatched.csv", sep = ""), header = TRUE)
#expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_DF_tumNormMatched_meQTL.csv", sep = ""), header = TRUE)
#methylation_df <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_tumNormMatched.csv", sep = ""))
#methylation_df_dependent <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_tumNormMatched_dependent.csv", sep = ""))

# Filtered FPKM
expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_filt_DF_tumNormMatched.csv", sep = ""), header = TRUE)
expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_filt_DF_tumNormMatched_meQTL.csv", sep = ""), header = TRUE)
methylation_df <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_filtFPKM_tumNormMatched.csv", sep = ""))
methylation_df_dependent <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_filtFPKM_tumNormMatched_dependent.csv", sep = ""))

patients_with_matched <- read.csv(paste(main_path, "patients_with_tumNormMatched_eQTL.csv", sep = ""), header = FALSE)[,2]
patients_with_matched_meQTL <- read.csv(paste(main_path, "patients_with_tumNormMatched_meQTL.csv", sep = ""), header = FALSE)[,2]


# DOWNSTREAM TARGETS (Will be integrated with full gene set for proteins that have no studied targets)
#downstream_targs_df <- read.csv(paste(main_path, "Linear Model/downstream_targets.csv", sep = ""), header = TRUE)   #TODO: GET THIS
