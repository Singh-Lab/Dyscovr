############################################################
### Process Chinese Triple Negative Breast Cancer Cohort Data
### Written By: Sara Geraghty, March 2023
############################################################

# Data retrieved from the following source: https://www.cell.com/cancer-cell/fulltext/S1535-6108(19)30096-0
# "Genomic and Transcriptomic Landscape of Triple-Negative Breast Cancers: Subtypes and Treatment Strategies" (Jiang et al., 2019)

# Data is accessible via NODE at: http://www.biosino.org/node/project/detail/OEP000155

library(GenomicRanges)
library(TRONCO)
library(maftools)
library(matrixStats)

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/TripleNeg_Chinese/"
output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/TripleNeg_Chinese/"


# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


# This file imports each data type of interest and formats them to match the TCGA files, 
# such that they will be able to flow through the various existing pipelines we have for each data type.

############################################################
### IMPORT CLINICAL PATIENT/ SAMPLE DATA
############################################################
tn_clinical_df <- read.csv(paste0(main_path, "clinical.csv"), header = T, check.names = F)

# Limit to female patients
tn_clinical_df <- tn_clinical_df[tn_clinical_df$Sex %in% c("Female", "Female "),]

# Run through create_sample_dataframe.R and create_patient_dataframe.R pipelines.
# Merge these files using combine_patient_and_sample_dfs function in additional_premodel_processing.R

############################################################
### IMPORT EXPRESSION DATA
############################################################
# Note that that expression data is microarray data, which needs to be handled differently -
# it appears to already be normalized, filtered, and log2 transformed (Potentially an RMA normalization?)
tn_expression_df <- read.table(paste0(main_path, "xiao_et_al_full_gene_expression_matrix.tsv"), 
                               header = T, sep = "\t", check.names = F, row.names = 1)
tn_expression_df <- t(tn_expression_df)

# TODO: revert FPKM-normalized expression data back to raw count data, and have a uniform pipeline?

# Filter genes with mean expression < 2 or standard deviation smaller than 1
# Method taken from 10.1186/s13058-015-0618-8, Supp. Methods p. 3
tn_row_means <- rowMeans(tn_expression_df, na.rm = TRUE)
tn_sds <- rowSds(as.matrix(tn_expression_df))
to_keep <- unique(c(which(tn_row_means > 2), which(tn_sds > 1)))
tn_expression_df <- as.data.frame(tn_expression_df[to_keep,]) # keeps 18,359 genes
  # ranges from 0.106 to 3.0669

# Convert the Hugo Symbols to ENSG IDs
ensg_ids <- unlist(lapply(rownames(tn_expression_df), function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'ensembl_gene_id'])), collapse = ";")))
tn_expression_df <- tn_expression_df[which(ensg_ids != ""),]
ensg_ids <- ensg_ids[ensg_ids != ""]
tn_expression_df <- cbind(data.frame("ensg_ids" = ensg_ids), tn_expression_df)

# Write to file
write.csv(tn_expression_df, paste0(output_path, "Expression/expression_FPKM_filtBySD1_DF.csv"))

############################################################
### IMPORT MUTATION DATA
############################################################
# Import mutation data file
tn_mutation_df <- read.table(paste0(main_path, "FUSCCTNBC_Mutations_V15.1.txt"), header = TRUE, 
                                   sep = "\t", check.names = FALSE)
# Use this file in the "process_mutation_data.R" pipeline to get the iprotein-level DF

# Generate the mutation count matrices, using helper functions from maf_helper_functions.R
tn_mutation_df <- tn_mutation_df[tn_mutation_df$Tumor_Sample_Barcode %fin% tn_clinical_df$Project_ID,]
tn_mutation_df_nonsynonymous <- tn_mutation_df[(tn_mutation_df$Variant_Classification == "Missense_Mutation") | 
                                                             (tn_mutation_df$Variant_Classification == "Nonsense_Mutation") |
                                                             (tn_mutation_df$Variant_Classification == "Splice_Site") |
                                                             (tn_mutation_df$Variant_Classification == "Nonstop_Mutation"),]
tn_mutation_df_silent <- tn_mutation_df[tn_mutation_df$Variant_Classification == "Silent",]

maf_nonsynon <- read.maf(tn_mutation_df_nonsynonymous)
maf_silent <- read.maf(tn_mutation_df_silent, vc_nonSyn = "Silent")

mut_count_matrix_nonsyn <- get_mut_count_matrix(maf_nonsynon)
mut_count_matrix_silent <- get_mut_count_matrix(maf_silent)

# Now, append rows of '0s' for those genes that are not mutated in any patient
all_gene_names <- rownames(tn_expression_df)

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

# Add total nonsynonymous mutational burden to clinical DF, using get_num_mutations_per_patient function
# from maf_helper_functions.R
tmb <- get_num_mutations_per_patient(mut_count_matrix_nonsyn_full)
tmb$Project_ID <- rownames(tmb)
tn_clinical_df <- merge(tn_clinical_df, tmb, by = "Project_ID")
colnames(tn_clinical_df)[which(colnames( tn_clinical_df) == "Total.Num.Mut")] <- "Total.Num.Muts"


############################################################
### IMPORT CNA DATA
############################################################
# Import CNA data file (0 is normal, -1 is deletion of 1 copy, -2 is deletion of 2 copies, 2+ is amplification)
# Will have to use bucketing method for this CNA data
tn_cna_df <- read.table(paste0(main_path, "FUSCCTNBC_MergedASCATSegmentations.txt"), header = TRUE, 
                              sep = "\t", check.names = FALSE)

# Remove the Entrez ID & HUGO ID columns, use ENSG IDs instead
ensg_ids <- unlist(lapply(tn_cna_df$Hugo_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'ensembl_gene_id'])), collapse = ";")))
tn_cna_df <- tn_cna_df[ensg_ids != "",]
ensg_ids <- ensg_ids[ensg_ids != ""]
tn_cna_df <- tn_cna_df[,3:ncol(tn_cna_df)]
rownames(tn_cna_df) <- make.names(ensg_ids, unique = TRUE)

# Write to file
write.csv(tn_cna_df, paste0(output_path, "CNV/CNA_DF_AllGenes_CancerOnly.csv"))


############################################################
### READ PROCESSED FILES BACK AND ADJUST HEADERS
############################################################

## MUTATION GENE TARG DF ##
tn_mutation_count_df_nonsynonymous <- read.csv(paste0(output_path, "Mutation/Mutation Count Matrices/mut_count_matrix_nonsynonymous.csv"), 
                                                     header = TRUE, check.names = FALSE)
tn_mutation_count_df_silent <- read.csv(paste0(output_path, "Mutation/Mutation Count Matrices/mut_count_matrix_silent.csv"), 
                                              header = TRUE, check.names = FALSE)
colnames(tn_mutation_count_df_nonsynonymous)[1] <- 'Gene_Symbol'
colnames(tn_mutation_count_df_silent)[1] <- 'Gene_Symbol'

tn_mutation_count_df_nonsynonymous$Swissprot <- unlist(lapply(tn_mutation_count_df_nonsynonymous$Gene_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'uniprot_gn_id'])), collapse = ";")))
tn_mutation_count_df_silent$Swissprot <- unlist(lapply(tn_mutation_count_df_silent$Gene_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,'uniprot_gn_id'])), collapse = ";")))


## CNA DF ##
tn_cna_df <- read.csv(paste0(output_path, "CNV/CNA_DF_AllGenes.csv"), 
                            header = TRUE, check.names = FALSE, row.names = 1)

## EXPRESSION DF ##
tn_expression_df <- read.csv(paste0(output_path, "Expression/expression_FPKM_filtBySD1_DF.csv"), 
                                   header = TRUE, check.names = FALSE)
colnames(tn_expression_df)[1] <- 'ensg_id'


## PATIENT SAMPLE DF ## 
tn_patient_sample_df <- read.csv(paste0(main_path, "Patient and Sample DFs/combined_patient_sample_DF.csv"), 
                                 header = T, check.names = F)
colnames(tn_patient_sample_df)[1] <- "sample_id"


############################################################
### ENSURE ALL IDS STILL OVERLAP
############################################################
mutation_ids_nonsyn <- unique(unlist(lapply(tn_mutation_df_nonsynonymous$Tumor_Sample_Barcode, function(x) 
  unlist(strsplit(x, ";", fixed = TRUE)))))  #271
mutation_ids_silent <- unique(unlist(lapply(tn_mutation_df_silent$Tumor_Sample_Barcode, function(x) 
  unlist(strsplit(x, ";", fixed = TRUE)))))  #271

cna_ids <- colnames(tn_cna_df) #401

expression_ids <- colnames(tn_expression_df)[2:ncol(tn_expression_df)]  #258

patient_sample_ids <- tn_patient_sample_df$sample_id #270

intersecting_patients <- intersect(cna_ids, intersect(expression_ids, intersect(mutation_ids_nonsyn, patient_sample_ids)))
print(length(intersecting_patients)) # 168 patients 

write.table(intersecting_patients, "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/TripleNeg_Chinese/intersecting_ids.txt",
            quote = FALSE, row.names = FALSE)

# Read back
intersecting_patients <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/TripleNeg_Chinese/intersecting_ids.txt", 
                                    header = FALSE)[,1]

############################################################
### IF NECESSARY, SUBSET FINALIZED DATASETS BY INTERSECTING IDs
############################################################
tn_mutation_count_df_missense_sub <- tn_mutation_count_df_missense[, c(1, which(colnames(tn_mutation_count_df_missense) %in% 
                                                                                              intersecting_patients))]
tn_mutation_count_df_misAndNon_sub <- tn_mutation_count_df_misAndNon[, c(1, which(colnames(tn_mutation_count_df_misAndNon) %in% 
                                                                                                intersecting_patients))]
tn_mutation_count_df_nonsynonymous_sub <- tn_mutation_count_df_nonsynonymous[, c(1, which(colnames(tn_mutation_count_df_nonsynonymous) %in% 
                                                                                                        intersecting_patients))]
tn_mutation_count_df_silent_sub <- tn_mutation_count_df_silent[, c(1, which(colnames(tn_mutation_count_df_silent) %in% 
                                                                                          intersecting_patients))]

tn_cna_df_sub <- tn_cna_df[, c(1, which(colnames(tn_cna_df) %in% intersecting_patients))]


tn_expression_df_sub <- tn_expression_df[, c(1, which(colnames(tn_expression_df) %in% intersecting_patients))]

tn_methylation_df_log2_sub <- tn_methylation_df_log2[, c(1, which(colnames(tn_methylation_df_log2) %in% intersecting_patients))]
tn_methylation_df_logit_sub <- tn_methylation_df_logit[, c(1, which(colnames(tn_methylation_df_logit) %in% intersecting_patients))]

tn_patient_sample_df_sub <- tn_patient_sample_df[tn_patient_sample_df$sample_id %in% intersecting_patients,]


# Do this for the labels on the regprot data frame as well
get_mutation_df_new_patient_labels <- function(tn_mutation_df, intersecting_patients) {
  mutation_df_new_patient_labels <- lapply(tn_mutation_df$Patient, function(x) {
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

tn_mutation_df_missense_labels <- get_mutation_df_new_patient_labels(tn_mutation_df_missense, 
                                                                           intersecting_patients)
tn_mutation_df_misAndNon_labels <- get_mutation_df_new_patient_labels(tn_mutation_df_misAndNon, 
                                                                            intersecting_patients)
tn_mutation_df_nonsynonymous_labels <- get_mutation_df_new_patient_labels(tn_mutation_df_nonsynonymous, 
                                                                                intersecting_patients)
tn_mutation_df_silent_labels <- get_mutation_df_new_patient_labels(tn_mutation_df_silent, 
                                                                         intersecting_patients)

# Remove NA
tn_mutation_df_missense_sub <- tn_mutation_df_missense[unlist(lapply(tn_mutation_df_missense_labels, function(x) 
  ifelse(is.na(x), FALSE, TRUE))),]
tn_mutation_df_misAndNon_sub <- tn_mutation_df_misAndNon[unlist(lapply(tn_mutation_df_misAndNon_labels, function(x) 
  ifelse(is.na(x), FALSE, TRUE))),]
tn_mutation_df_nonsynonymous_sub <- tn_mutation_df_nonsynonymous[unlist(lapply(tn_mutation_df_nonsynonymous_labels, function(x) 
  ifelse(is.na(x), FALSE, TRUE))),]
tn_mutation_df_silent_sub <- tn_mutation_df_missense[unlist(lapply(tn_mutation_df_silent_labels, function(x) 
  ifelse(is.na(x), FALSE, TRUE))),]

# Update the DF with the new intersecting patient labels
tn_mutation_df_missense_sub$Patient <- tn_mutation_df_missense_labels[!is.na(tn_mutation_df_missense_labels)]
tn_mutation_df_missense_sub$Patient <- unlist(lapply(tn_mutation_df_missense_sub$Patient, function(x) return(unlist(x))))

tn_mutation_df_misAndNon_sub$Patient <- tn_mutation_df_misAndNon_labels[!is.na(tn_mutation_df_misAndNon_labels)]
tn_mutation_df_misAndNon_sub$Patient <- unlist(lapply(tn_mutation_df_misAndNon_sub$Patient, function(x) return(unlist(x))))

tn_mutation_df_nonsynonymous_sub$Patient <- tn_mutation_df_nonsynonymous_labels[!is.na(tn_mutation_df_nonsynonymous_labels)]
tn_mutation_df_nonsynonymous_sub$Patient <- unlist(lapply(tn_mutation_df_nonsynonymous_sub$Patient, function(x) return(unlist(x))))

tn_mutation_df_silent_sub$Patient <- tn_mutation_df_silent_labels[!is.na(tn_mutation_df_silent_labels)]
tn_mutation_df_silent_sub$Patient <- unlist(lapply(tn_mutation_df_silent_sub$Patient, function(x) return(unlist(x))))

############################################################
### WRITE THESE NEWLY SUBSETTED FILES BACK
############################################################
write.csv(tn_mutation_count_df_missense_sub, paste0(output_path, "Linear Model/mutation_count_df_missense_IntersectPatients.csv"))
write.csv(tn_mutation_count_df_misAndNon_sub, paste0(output_path, "Linear Model/mutation_count_df_missense_nonsense_IntersectPatients.csv"))
write.csv(tn_mutation_count_df_nonsynonymous_sub, paste0(output_path, "Linear Model/mutation_count_df_nonsynonymous_IntersectPatients.csv"))
write.csv(tn_mutation_count_df_silent_sub, paste0(output_path, "Linear Model/mutation_count_df_silent_IntersectPatients.csv"))

write.csv(tn_cna_df_sub, paste0(output_path, "Linear Model/cna_df_allgenes_IntersectPatients.csv"))

write.csv(tn_expression_df_sub, paste0(output_path, "Linear Model/expression_df_log2intensity_filtBySD0.2_IntersectPatients.csv"))

write.csv(tn_methylation_df_log2_sub, paste0(output_path, "Linear Model/methylation_df_RRBS_log2_IntersectPatients.csv"))
write.csv(tn_methylation_df_logit_sub, paste0(output_path, "Linear Model/methylation_df_RRBS_logit_IntersectPatients.csv"))

write.csv(tn_patient_sample_df_sub, paste0(output_path, "Linear Model/patient_sample_df_cibersort_total_frac_IntersectPatients.csv"))

write.csv(tn_mutation_df_missense_sub, paste0(output_path, "Linear Model/iprotein_mutation_df_missense.csv"))
write.csv(tn_mutation_df_misAndNon_sub, paste0(output_path, "Linear Model/iprotein_mutation_df_missense_nonsense.csv"))
write.csv(tn_mutation_df_nonsynonymous_sub, paste0(output_path, "Linear Model/iprotein_mutation_df_nonsynonymous.csv"))
write.csv(tn_mutation_df_silent_sub, paste0(output_path, "Linear Model/iprotein_mutation_df_silent.csv"))


#############################################################
### GENERATE PROTEIN INPUT TABLES (DRIVERS)
#############################################################
# Import DF with known driver genes
driver_gene_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/GRCh38_driver_gene_list.tsv", 
                           sep = "\t", header = TRUE, comment.char = "#", skip = 10)

iprotein_prots_missense <- unique(tn_mutation_df_missense_sub$Swissprot)  # 122 unique mutated proteins
iprotein_prots_misAndNon <- unique(tn_mutation_df_misAndNon_sub$Swissprot)  # 122 unique mutated proteins
iprotein_prots_nonsyn <- unique(tn_mutation_df_nonsynonymous_sub$Swissprot)  # 122 unique mutated proteins
iprotein_prots_silent <- unique(tn_mutation_df_silent_sub$Swissprot)  # 122 unique mutated proteins

# Proceed with nonsynonymous mutations; use functions from 'generate_protein_input_tables.R'
iprotein_prots_ensg <- convert_to_ensembl(iprotein_prots_nonsyn)

# Create initial DF with all proteins
recombine_into_df_and_write(iprotein_prots, iprotein_prots_ensg, "iprotein", prot_path)

# Call filter to include only proteins mutated in at least 5 patients, and at least 2.5% or 5% of patients

# Nonsynonymous
regprot_input_df_all_5perc <- create_regulatory_prot_input_df(tn_mutation_df_nonsynonymous, 
                                                              0.05, 5, "i-protein", intersecting_patients, 
                                                              NA, all_genes_id_conv, NA, NA, "tn")
regprot_input_df_driver_5perc <- create_regulatory_prot_input_df(tn_mutation_df_nonsynonymous, 
                                                                 0.05, 5, "i-protein", intersecting_patients, 
                                                                 driver_gene_df, all_genes_id_conv, NA, NA, "tn")

regprot_input_df_all_2.5perc <- create_regulatory_prot_input_df(tn_mutation_df_nonsynonymous, 
                                                                0.025, 5, "i-protein", intersecting_patients, 
                                                                NA, all_genes_id_conv, NA, NA, "tn")
regprot_input_df_driver_2.5perc <- create_regulatory_prot_input_df(tn_mutation_df_nonsynonymous, 
                                                                   0.025, 5, "i-protein", intersecting_patients, 
                                                                   driver_gene_df, all_genes_id_conv, NA, NA, "tn")


write.csv(regprot_input_df_all_5perc, paste0(output_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.05Freq_all_nonsynonymous.csv"))
write.csv(regprot_input_df_driver_5perc, paste0(output_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.05Freq_drivers_nonsynonymous.csv"))
write.csv(regprot_input_df_all_2.5perc, paste0(output_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.025Freq_all_nonsynonymous.csv"))
write.csv(regprot_input_df_driver_2.5perc, paste0(output_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.025Freq_drivers_nonsynonymous.csv"))


# Repeat this for individual subtypes
patient_set_lumA <- intersect(intersecting_patients, tn_clinical_patient_df[tn_clinical_patient_df$CLAUDIN_SUBTYPE == "LumA", 'PATIENT_ID'])
patient_set_lumB <- intersect(intersecting_patients, tn_clinical_patient_df[tn_clinical_patient_df$CLAUDIN_SUBTYPE == "LumB", 'PATIENT_ID'])
patient_set_basal <- intersect(intersecting_patients, tn_clinical_patient_df[tn_clinical_patient_df$CLAUDIN_SUBTYPE %in% c('claudin-low', 'Basal'), 'PATIENT_ID'])
patient_set_her2 <- intersect(intersecting_patients, tn_clinical_patient_df[tn_clinical_patient_df$CLAUDIN_SUBTYPE == 'Her2', 'PATIENT_ID'])


#' Function to keep only patients in the regprot mutation DF that are a member of the
#' given patient set
#' @param tn_mutation_df a regprot mutation DF with tn IDs
#' @param patient_set a vector of tn patient IDs to keep
get_regprot_df_subsetted_by_subtype <- function(tn_mutation_df, patient_set) {
  new_rows <- lapply(1:nrow(tn_mutation_df), function(i) {
    pats <- unlist(strsplit(tn_mutation_df[i, 'Patient'], ";", fixed = TRUE))
    pats_in_set <- paste(intersect(pats, patient_set), collapse = ";")
    
    row <- as.data.frame(tn_mutation_df[i,])
    colnames(row) <- colnames(tn_mutation_df)
    row[,'Patient'] <- pats_in_set
    return(row)
  })
  new_df <- do.call("rbind", new_rows)
  return(new_df)
}

tn_mutation_df_nonsynonymous_lumA <- get_regprot_df_subsetted_by_subtype(tn_mutation_df_nonsynonymous,
                                                                               patient_set_lumA)
tn_mutation_df_nonsynonymous_lumB <- get_regprot_df_subsetted_by_subtype(tn_mutation_df_nonsynonymous,
                                                                               patient_set_lumB)
tn_mutation_df_nonsynonymous_basal <- get_regprot_df_subsetted_by_subtype(tn_mutation_df_nonsynonymous,
                                                                                patient_set_basal)
tn_mutation_df_nonsynonymous_her2 <- get_regprot_df_subsetted_by_subtype(tn_mutation_df_nonsynonymous,
                                                                               patient_set_her2)

regprot_input_df_all_5perc_lumA <- create_regulatory_prot_input_df(tn_mutation_df_nonsynonymous_lumA, 
                                                                   0.05, 5, "i-protein", patient_set_lumA, 
                                                                   NA, all_genes_id_conv, NA, NA, "tn")
regprot_input_df_driver_5perc_lumA <- create_regulatory_prot_input_df(tn_mutation_df_nonsynonymous_lumA, 
                                                                      0.05, 5, "i-protein", patient_set_lumA, 
                                                                      driver_gene_df, all_genes_id_conv, NA, NA, "tn")

regprot_input_df_all_2.5perc_lumA <- create_regulatory_prot_input_df(tn_mutation_df_nonsynonymous_lumA, 
                                                                     0.025, 5, "i-protein", patient_set_lumA, 
                                                                     NA, all_genes_id_conv, NA, NA, "tn")
regprot_input_df_driver_2.5perc_lumA <- create_regulatory_prot_input_df(tn_mutation_df_nonsynonymous_lumA, 
                                                                        0.025, 5, "i-protein", patient_set_lumA, 
                                                                        driver_gene_df, all_genes_id_conv, NA, NA, "tn")


write.csv(regprot_input_df_all_5perc, paste0(output_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.05Freq_all_nonsynonymous.csv"))
write.csv(regprot_input_df_driver_5perc, paste0(output_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.05Freq_drivers_nonsynonymous.csv"))
write.csv(regprot_input_df_all_2.5perc, paste0(output_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.025Freq_all_nonsynonymous.csv"))
write.csv(regprot_input_df_driver_2.5perc, paste0(output_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.025Freq_drivers_nonsynonymous.csv"))

