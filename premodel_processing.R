############################################################
### Additional Pre-Dyscovr Processing 
### Written by Sara Geraghty, Princeton University
### https://www.biorxiv.org/content/10.1101/2024.11.20.624509v1
############################################################

options("install.lock"=F)
library(dplyr)
library(ggplot2)
library(grid)
library(parallel)
library(VennDiagram)
library(data.table)
library(TCGAbiolinks)

# Contains final pre-Dyscovr processing of data files to use as input in Dyscovr
# framework; is run post-individual data type processing, prior to 
# create_lm_input_tables.R 

# This includes:
# 1. Combine patient and sample DFs into one
# 2. Limiting to only overlapping samples (final check for intersection)
# 3. Construct a data frame of targets for a given protein of interest (often 
# for testing purposes)
# 4. Helper functions to find specific patient sets for particular purposes: 
  # Samples that have both an amplification and mutation or deletion and mutation 
  # in the same gene; samples with mutation in one gene and no mutation in another;
  # male or female patients 

# Local PATH to directory containing processed data files
PATH <- paste0(getwd(), "Output/")
OUTPUT_PATH <- paste0(getwd(), "Output/Dyscovr_Input/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
############################################################
# IMPORT NECESSARY FILES FOR DYSCOVR (EXPRESSION, METHYLATION,
# GENE MUTATION CNA, & PATIENT/ SAMPLE INFORMATION)
############################################################
############################################################

############################################################
# IMPORT EXPRESSION FILES
############################################################
expression_df <- fread(paste0(
  PATH, "Expression/expression_quantile_normalized_sklearn.csv"), header = T)
colnames(expression_df)[1] <- 'ensg_id'

# Adjust TCGA sample IDs to be just patient and sample (e.g. XXXX-XX) 
expression_df <- expression_df[, c(1, (which(grepl(
  "TCGA", colnames(expression_df)[2:ncol(expression_df)]))+1)), with = F]
colnames(expression_df)[2:ncol(expression_df)] <- unlist(lapply(
  colnames(expression_df)[2:ncol(expression_df)], function(x)
    paste(unlist(strsplit(x, "-", fixed = T))[3:4], collapse = "-")))

############################################################
# IMPORT METHYLATION FILES
############################################################
methylation_df_M <- fread(paste0(
  PATH, "Methylation/methylation_DF_M_cancer_only.csv"), header = T)
methylation_df_Beta <- fread(paste0(
  PATH, "Methylation/methylation_DF_Beta_cancer_only.csv"), header = T)

# Set whichever methylation DF type is desired
methylation_df <- methylation_df_M
methylation_df <- methylation_df[,2:ncol(methylation_df)]

############################################################
# IMPORT MUTATION FILES
############################################################
mutation_df <- fread(paste0(
  PATH, "Mutation/mut_count_matrix_nonsynonymous.csv"), header = T)
#mutation_df <- mutation_df[,2:ncol(mutation_df)]
colnames(mutation_df)[1] <- 'Gene_Symbol'

# If needed, fix the sample IDs (e.g. XXXX-XX-XXXX-XX, instead of XXXX-XX)
colnames(mutation_df)[2:ncol(mutation_df)] <- unlist(lapply(
  colnames(mutation_df)[2:ncol(mutation_df)], function(x) 
    paste(unlist(strsplit(x, "-", fixed = T))[3:4], collapse = "-"))) 

############################################################
# IMPORT CNA FILES
############################################################
cna_df <- read.csv(paste0(
  PATH, "CNA/CNA_DF_AllGenes_CancerOnly_MergedIsoforms.csv"), 
  header = T, row.names = 1, check.names = F)

############################################################
# IMPORT PATIENT FILES
############################################################
patient_df <- read.csv(paste0(
  OUTPUT_PATH, "Patient/patient_dataframe.csv"), 
  header = T, row.names = 1)
patient_df <- read.csv(paste0(
  OUTPUT_PATH, "Patient/patient_dataframe_inclSubtypes.csv"), 
  header = T, row.names = 1)

# If BRCA or another gender-specific cancer type, eliminate the gender column
is_brca <- F
if(is_brca) {patient_df <- patient_df[, -which(colnames(patient_df) == "Gender")]}

# Import list of per-cancer patient DFs
per_cancer_patient_df_fns <- intersect(list.files(paste0(OUTPUT_PATH, "Patient/"), 
                                        pattern = "_inclSubtypes.csv"),
                                       list.files(paste0(OUTPUT_PATH, "Patient/"), 
                                                  pattern = "patient_dataframe"))
per_cancer_patient_df_fns <- per_cancer_patient_df_fns[
  !(per_cancer_patient_df_fns == "patient_dataframe_inclSubtypes.csv")]
per_cancer_patient_dfs <- lapply(per_cancer_patient_df_fns, function(fn)
  read.csv(paste0(OUTPUT_PATH, paste0("Patient/", fn)), header = T, row.names = 1))

############################################################
# IMPORT SAMPLE FILES
############################################################
sample_df <- read.csv(paste0(
  OUTPUT_PATH, "Sample/sample_dataframe_cibersort_total_frac_washu.csv"), 
  header = T, row.names = 1)  

############################################################
############################################################
# COMBINE PATIENT AND SAMPLE FILES
############################################################
############################################################
#' Combines the patient and sample data frames such that there is one row per 
#' patient sample. Multiple samples from the same patient will have duplicate 
#' patient information.
#' @param patient_df the patient data frame (patients are rows, patient 
#' characteristics are columns)
#' @param sample_df the sample data frame (patient samples are rows, sample 
#' characteristics are columns)
#' @param dataset either 'tcga' or 'metabric'
combine_patient_and_samp_dfs <- function(patient_df, sample_df, dataset) {
  # For each sample, retrieve the corresponding patient info
  sample_dfs <- lapply(1:nrow(sample_df), function(i) {
    # Get the sample barcode and its row
    sample_barcode <- rownames(sample_df)[i]
    sample_row <- sample_df[i,]

    # Get the corresponding row from the patient DF, if TCGA
    if(length(sample_barcode) > 0) {
      patient <- sample_barcode
      if(dataset == "tcga") {
        patient <- unlist(strsplit(sample_barcode, "-", fixed = T))[1]
      }
      print(patient)

      if(patient %fin% rownames(patient_df)) {
        patient_row <- patient_df[rownames(patient_df) == patient,]

        # Combine these together, so that the sample ID is maintained as rowname
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
    } else {return(NA)}
  })
  is_na <- which(is.na(sample_dfs))
  sample_dfs <- sample_dfs[!is_na]
  
  # Recombine all these new, combined rows into a full DF
  if(length(sample_dfs) > 0) {
    combined_df <- do.call(rbind, sample_dfs)
    rownames(combined_df) <- rownames(sample_df)[-is.na]
  } else {combined_df <- NA}
  
  return(combined_df)
}

# Call this function for one cancer type, or pan-cancer
combined_pat_samp_df <- combine_patient_and_samp_dfs(patient_df, sample_df, 
                                                     'tcga')
# Remove rows/ columns that are entirely NA
combined_pat_samp_df <- combined_pat_samp_df[
  rowSums(is.na(combined_pat_samp_df)) != ncol(combined_pat_samp_df),]
combined_pat_samp_df <- combined_pat_samp_df[
  ,colSums(is.na(combined_pat_samp_df)) != nrow(combined_pat_samp_df)]

# Write the results to file
write.csv(combined_pat_samp_df, paste0(
  OUTPUT_PATH, "Sample/combined_patient_sample_DF_cibersort_total_frac_washu.csv"))

#
#
#

# Call this function per-cancer across many cancer types
sample_df_pats <- unlist(lapply(rownames(sample_df), function(x) 
  unlist(strsplit(x, "-", fixed = T))[1]))

combined_pat_samp_dfs <- lapply(pan_cancer_patient_dfs, function(pat_df) {
  sample_df_sub <- sample_df[which(sample_df_pats %fin% rownames(pat_df)),]
  if((length(sample_df_to_sub) > 0) & (length(pat_df) > 0)) {
    comb_df <- combine_patient_and_samp_dfs(pat_df, sample_df, 'tcga')
    if(length(comb_df) == 1) {
      if(!is.na(comb_df)) {
        comb_df <- comb_df[rowSums(is.na(comb_df)) != ncol(comb_df),]
        comb_df <- comb_df[,colSums(is.na(comb_df)) != nrow(comb_df)]
      }
    } else {
      comb_df <- comb_df[rowSums(is.na(comb_df)) != ncol(comb_df),]
      comb_df <- comb_df[,colSums(is.na(comb_df)) != nrow(comb_df)]
    }
    return(comb_df)
  } else{return(NA)}
})
names(combined_pat_samp_dfs) <- names(pan_cancer_patient_dfs)

# Write to files
lapply(1:length(combined_pat_samp_dfs), function(i) {
  cancer_name <- names(combined_pat_samp_dfs)[i]
  fn <- paste0("Sample/combined_patient_sample_DF_cibersort_total_frac_washu_", 
               paste0(cancer_name, "_inclSubtypes.csv"))
  write.csv(combined_pat_samp_dfs[[i]], paste0(OUTPUT_PATH, fn))
})

#
#
#

# Read these files back, as needed
patient_sample_df <- read.csv(paste0(
  OUTPUT_PATH, "Sample/combined_patient_sample_DF_cibersort_total_frac_washu.csv"),
  header = T, row.names = 1)

patient_sample_files <- intersect(list.files(paste0(OUTPUT_PATH, "Sample/"), 
                                   pattern = "_inclSubtypes.csv"), 
                                  list.files(paste0(OUTPUT_PATH, "Sample/"), 
                                             pattern = "combined_patient_sample"))
patient_sample_dfs <- lapply(patient_sample_files, function(f)
  read.csv(paste0(OUTPUT_PATH, paste0("Sample/", f)), header = T, row.names = 1))
names(patient_sample_dfs) <- unlist(lapply(patient_sample_files, function(f) 
  unlist(strsplit(f, "_", fixed = T))[10]))


############################################################
############################################################
# LIMIT TO OVERLAPPING/ INTERSECTING PATIENTS, POST-PROCESSING
############################################################
############################################################
#' Given a TCGA sample ID (XXXX-XXX, e.g. A0WY-01A), splits it 
#' apart and returns a list of unique patient IDs (XXXX, e.g. A0WY)
#' that have at least one cancer sample
#' @param sample_ids a vector of TCGA sample IDs (XXXX-XXX, e.g. A0WY-01A)
get_unique_patients <- function(sample_ids) {
  unique_patients <- unique(unlist(lapply(sample_ids, function(x) {
    spl_name <- unlist(strsplit(x, "-", fixed = T))
    if(startsWith(spl_name[2], "0")) {return(spl_name[1])}
  })))
  return(unique_patients)
}

# Call this function for all data type (check that other CNA, expression,
# etc. data frames have the same patients)
cna_patients <- get_unique_patients(colnames(cna_df))    
methylation_patients <- get_unique_patients(
  colnames(methylation_df)[2:ncol(methylation_df)])  
expression_patients <- get_unique_patients(
  colnames(expression_df)[2:ncol(expression_df)])  
mutation_patients <- unlist(get_unique_patients(
  colnames(mutation_df)[!(colnames(mutation_df) %fin% 
                                 c('Swissprot', 'Gene_Symbol'))])) 
clinical_patients <- get_unique_patients(rownames(patient_sample_df))  
#clinical_patients_perCt <- unique(unlist(lapply(patient_sample_dfs, function(df) 
 # get_unique_patients(rownames(df)))))  

intersecting_patients <- intersect(cna_patients, 
                                   intersect(methylation_patients,
                                             intersect(expression_patients,
                                                       intersect(mutation_patients, 
                                                                 clinical_patients))))
print(length(intersecting_patients)) 

#
#
#

#' Given a list of intersecting patient IDs (XXXX, e.g. A0WY), subsets a 
#' given data frame with column names that contain sample IDs (XXXX-XXX, 
#' e.g. A0WY-01A) to include only those given patients
#' @param intersecting_patients a vector of intersecting patient IDs
#' @param df the data frame to subset using the patient ids
#' @param colNames a T/F value indicating whether the patient/sample IDs are
#' found on the rows or the columns of the DF (T is columns, F is rows)
subset_by_intersecting_ids <- function(intersecting_patients, df, colNames) {
  if(colNames == T) {
    # Keep only intersecting patients
    just_patients <- unlist(lapply(colnames(df)[2:ncol(df)], function(x)
      unlist(strsplit(x, "-", fixed = T))[1]))
    cols_to_keep <- as.numeric(c(1,(which(just_patients %fin% 
                                            intersecting_patients)+1)))
    df <- as.data.frame(df)
    df_adj <- df[, cols_to_keep]
    
    # Keep only cancer samples (final check)
    just_samples <- unlist(lapply(colnames(df_adj)[2:ncol(df_adj)], function(x)
      unlist(strsplit(x, "-", fixed = T))[2]))
    df_adj <- df_adj[,c(1,which(unlist(lapply(just_samples, function(x) 
      startsWith(x, "0"))))+1)]
    
    # Remove any duplicate samples
    df_adj <- df_adj[, !(grepl("-1", colnames(df_adj), fixed = T))]
    
  } else {
    if(length(unlist(strsplit(rownames(df)[1], "-", fixed = T))) == 4) {
      rownames(df) <- unlist(lapply(rownames(df), function(x) 
        paste(unlist(strsplit(x, "-", fixed = T))[3:4], collapse = "-")))
    }
    # Keep only intersecting patients
    just_patients <- unlist(lapply(rownames(df), function(x)
      unlist(strsplit(x, "-", fixed = T))[1]))
    df_adj <- df[which(just_patients %fin% intersecting_patients),]
    if(length(df_adj) == 0) {return(NA)}

    # Keep only cancer samples (final check)
    just_samples <- unlist(lapply(rownames(df_adj), function(x)
      unlist(strsplit(x, "-", fixed = T))[2]))
    df_adj <- df_adj[which(unlist(lapply(just_samples, function(x) 
      startsWith(x, "0")))),]
    
    # Remove any duplicate samples
    df_adj <- df_adj[!(grepl('-1', rownames(df_adj), fixed = T)),]
  }
  return(df_adj)
}

# Call this function
cna_df_sub <- subset_by_intersecting_ids(intersecting_patients, cna_df, T)  
methylation_df_sub <- subset_by_intersecting_ids(intersecting_patients, 
                                                 methylation_df, T) 
mutation_df_sub <- subset_by_intersecting_ids(intersecting_patients, 
                                              mutation_df, T)  
expression_df_sub <- subset_by_intersecting_ids(intersecting_patients, 
                                                expression_df, T) 
rownames(patient_sample_df) <- unlist(lapply(
  rownames(patient_sample_df), function(x) 
    paste(unlist(strsplit(x, "-", fixed = T))[3:4], collapse = "-")))
patient_sample_df_sub <- subset_by_intersecting_ids(
  intersecting_patients, patient_sample_df, F) 

# Apply this to the per-cancer patient/sample files
patient_sample_dfs_sub <- lapply(patient_sample_dfs, function(df) 
  subset_by_intersecting_ids(intersecting_patients, df, F))
names(patient_sample_dfs_sub) <- names(patient_sample_dfs)


############################################################
# ADD ADDITIONAL GENE IDENTIFIERS
############################################################
# Add ENSG IDs to methylation DF
methylation_df_sub$ensg_ids <- unlist(lapply(methylation_df_sub$Gene_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x,
                                        'ensembl_gene_id'])), collapse = ";")))

# Add ENSG IDs to expression DF
expression_df_sub$ensg_id <- unlist(lapply(expression_df_sub$ensg_id, function(x) 
    unlist(strsplit(x, ".", fixed = T))[1]))

# Add Swissprot IDs to mutation DF
mutation_df_sub$Swissprot <- unlist(lapply(mutation_df_sub$Gene_Symbol, function(x) 
  paste(unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 
                                 'uniprot_gn_id']), collapse = ";")))


############################################################
# WRITE TO FILES FOR USE IN DYSCOVR
############################################################
write.csv(cna_df_sub, paste0(
  OUTPUT_PATH, "CNA/CNA_AllGenes_IntersectPatients.csv"))
write.csv(methylation_df_sub, paste0(
  OUTPUT_PATH, "Methylation/methylation_M_CancerOnly_IntersectPatients.csv"))
write.csv(expression_df_sub, paste0(
  OUTPUT_PATH, "Expression/expression_quantile_norm_IntersectPatients.csv"))
write.csv(mutation_df_sub, paste0(
  OUTPUT_PATH, "Mutation/mut_count_matrix_nonsynonymous_IntersectPatients.csv"))
write.csv(patient_sample_df_sub, paste(
  OUTPUT_PATH, "Sample/combined_patient_sample_IntersectPatients.csv"))

lapply(1:length(patient_sample_dfs_sub), function(i) {
  ct <- names(patient_sample_dfs_sub)[i]
  fn <- paste0("Sample/combined_patient_sample_cibersort_total_frac_washu", 
               paste0(ct, "_inclSubtypes_IntersectPatients.csv"))
  write.csv(patient_sample_dfs_sub[[i]], paste0(OUTPUT_PATH, fn))
})

# Write these intersecting IDs
write.table(intersecting_patients, paste0(OUTPUT_PATH, "intersecting_ids.txt"))

# Read back, as needed
intersecting_patients <- read.table(paste0(
  OUTPUT_PATH, "intersecting_ids.txt"))[,1]


############################################################
############################################################
# MAKE TESTER TARGETS DATAFRAME FOR GIVEN PROTEIN-OF-INTEREST
############################################################
############################################################
sample_protein_uniprot <- "P04637" 
sample_protein_ensg <- "ENSG00000141510"
sample_protein_name <- "P53"

### OPTION 1: CURATED TARGETS ###
# Download a table of curated targets and import here (column names are 
# protein names)
curated_targets_df <- fread(paste0(
  PATH, "Curated_Target_Data/curated_targets_df_full.csv"), header = T)
sample_targets <- curated_targets_df[,colnames(curated_targets_df) == 
                                       sample_protein_name] 
sample_targets <- sample_targets[!is.na(sample_targets)] 

# For TP53, get targets from publication 
# (Fischer, 2017, https://doi.org/10.1038/onc.2016.502)
tp53_targets_df <- fread(paste0(PATH, "Validation_Files/TP53_targets.csv"), 
                         header = T)
tp53_targets_df <- tp53_targets_df[tp53_targets_df$Gene.Symbol != "",]
sample_targets <- tp53_targets_df$Gene.Symbol

### OPTION 2: ChIP-eat file targets ###
sample_targets <- read.table(paste0(PATH, "ChIP-eat/tp53_chipeat_targets.txt"))[,1]

#
#
#

# Get target Swissprot IDs and ENSG IDs from gene names
sample_targets_swissprot <- unlist(lapply(sample_targets, function(x) 
  paste(unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 
                                 'uniprot_gn_id']), collapse = ";")))
sample_targets_swissprot <- unique(
  sample_targets_swissprot[sample_targets_swissprot != ""])

sample_targets_ensg <- unlist(lapply(sample_targets_swissprot, function(x) 
  paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == 
                                   unlist(strsplit(x, ";", fixed = T))[1], 
                                 'ensembl_gene_id']), collapse = ";")))


### OPTION 3: ALL POSSIBLE GENE TARGETS ###
all_genes_id_conv_sub <- all_genes_id_conv[,c('ensembl_gene_id', 'uniprot_gn_id')]
all_genes_id_conv_sub <- distinct(all_genes_id_conv_sub[
  all_genes_id_conv_sub$uniprot_gn_id != "",])
write.csv(all_genes_id_conv_sub, paste0(PATH, "all_genes_id_conv_subset.csv"))

# Get the unique ENSG and Swissprot ids from gene names
sample_targets_ensg <- unique(all_genes_id_conv_sub$ensembl_gene_id)   
sample_targets_swissprot <- unlist(lapply(sample_targets_ensg, function(x) 
  paste(unlist(all_genes_id_conv_sub[all_genes_id_conv_sub$ensembl_gene_id == x, 
                                     'uniprot_gn_id']), collapse = ";")))


# For all options, make ENSG and Swissprot IDs into a data frame
sample_targets_DF <- data.frame("swissprot" = sample_targets_swissprot, 
                                "ensg" = sample_targets_ensg)

# Write to CSV
# Option 1
write.csv(sample_targets_DF, paste0(PATH, "TP53/tp53_curated_targets.csv"))
# Option 2
write.csv(sample_targets_DF, paste0(PATH, "TP53/tp53_chipeat_targets.csv"))
# Option 3
write.csv(sample_targets_DF, paste(PATH, "allgene_targets.csv"))


############################################################
############################################################
# FIND SAMPLES WITH MUTATION AND CNA IN SAME GENE
############################################################
############################################################
#' Given a sample set of interest, and whether we are interested in 
#' amplifications or deletions, identifies the overlap of samples that have both 
#' a mutation event and an amplification or deletion in the given gene. 
#' Eliminate overlap in the given gene and return this new patient list
#' @param mut_count_matrix a mutation count matrix with sample IDs as column 
#' names (XXXX-XX) and a column called Gene_Symbol with the external gene name
#' @param cna_df a CNA data frame with the raw copy number of each gene (row 
#' names are gene ENSG IDs) for each sample (colnames are sample IDs, XXXX-XX)
#' @param gn the external gene name of the gene of interest
#' @param ensg the ENSG ID of the gene of interest
#' @param samples a vector of patient samples (XXXX-XX) in the cohort of 
#' interest
#' @param ampOrDel either "Amplification" to indicate we are interested in 
#' overlap with amplifications, or "Deletion" to indicate we are interested in
#'overlap with deletions
find_samps_w_muts_and_cnas_in_same_gene <- function(mut_count_matrix, cna_df, gn, 
                                                    ensg, samples, ampOrDel) {
  
  colnames(mut_count_matrix) <- unlist(lapply(
    colnames(mut_count_matrix), function(x) {
      return(unlist(strsplit(x, "-", fixed = T))[1])
  }))
  colnames(cna_df) <- unlist(lapply(colnames(cna_df), function(x) {
    return(unlist(strsplit(x, "-", fixed = T))[1])
  }))
  
  samples_mut <- colnames(mut_count_matrix)[intersect(which(as.numeric(unlist(
    mut_count_matrix[mut_count_matrix$Gene_Symbol == gn,])) == 1), 
    which(colnames(mut_count_matrix) %fin% samples))]
  
  if(ampOrDel == "Amplification") {
    samples_cna <- colnames(cna_df)[intersect(
      which(as.numeric(unlist(cna_df[rownames(cna_df) == ensg,])) > 2), 
      which(colnames(cna_df) %fin% samples))]
  } else {
    samples_cna <- colnames(cna_df)[intersect(
      which(as.numeric(unlist(cna_df[rownames(cna_df) == ensg,])) < 2), 
      which(colnames(cna_df) %fin% samples))]
  }
  
  # Plot a Venn Diagram to visualize the overlap between these two groups
  plt <- venn.diagram(list(samples_mut, samples_cna), 
                      category.names = c("Mutation", ampOrDel), 
                      filename = NULL, output = T, lwd = 2, lty = 'blank', 
                      fill = c("red", "blue"), cex = 2, fontface = "bold", 
                      fontfamily = "sans", cat.pos = 180, 
                      cat.fontfamily = "sans", cat.cex = 1)
  grid::grid.draw(plt)
  
  # Get & return the samples without overlap
  samples_no_overlap <- c(setdiff(samples_mut, samples_cna), 
                          setdiff(samples_cna, samples_mut))
  return(samples_no_overlap)
}


gn <- "TP53"
ensg <- "ENSG00000141510"
gn <- "PIK3CA"
ensg <- "ENSG00000121879"


############################################################
############################################################
# CREATE SAMPLE ID FILES WITH NO OVERLAP IN MUTATION 
# BETWEEN TWO GIVEN GENES
############################################################
############################################################
gene1 <- "TP53"
gene2 <- "PIK3CA"

#' Get patients without a mutation the given gene1 AND gene2 (eliminate overlap)
#' @param mutation_count_matrix a mutation count matrix of all samples and genes
#' @param gene1 the gene name of the first gene 
#' @param gene2 the gene name of the second gene 
eliminate_mutation_overlap <- function(mutation_count_matrix, gene1, gene2) {
  mut_count_mat_sub <- mutation_count_matrix[
    (mutation_count_matrix$Gene_Symbol == gene1) | 
      (mutation_count_matrix$Gene_Symbol == gene2),]
  g1_mut_indices <- unlist(lapply(1:ncol(mut_count_mat_sub), function(i) 
    ifelse(as.numeric(mut_count_mat_sub[1,i]) > 0, i, NA)))
  g2_mut_indices <- unlist(lapply(1:ncol(mut_count_mat_sub), function(i) 
    ifelse(as.numeric(mut_count_mat_sub[2,i]) > 0, i, NA)))
  
  g1_mut_indices <- g1_mut_indices[!is.na(g1_mut_indices)]
  g2_mut_indices <- g2_mut_indices[!is.na(g2_mut_indices)]
  
  both_mut <- intersect(g1_mut_indices, g2_mut_indices)
  
  mut_count_mat_sub_excl_doub_mut <- mut_count_mat_sub[, -both_mut]
  patients_no_mut_overlap <- unlist(lapply(
    colnames(mut_count_mat_sub_excl_doub_mut), function(x)
      unlist(strsplit(x, "-", fixed = T))[1]))
  patients_no_mut_overlap <- patients_no_mut_overlap[!(patients_no_mut_overlap == 
                                                         "Gene_Symbol")]
  
  return(patients_no_mut_overlap)
}


# Call function
patients_no_mut_overlap <- eliminate_mutation_overlap(mutation_df, gene1, gene2)

# write to files
write(patients_no_mut_overlap, paste0(
  PATH, "Patient_Subsets/NotTP53andPIK3CAmut_patient_ids.txt"))

# Get overlap with other patient subsets
patient_set_lumA <- read.table(paste0(
  PATH, "Patient_Subsets/Luminal.A_patient_ids.txt"), header = T)[,1]
patient_set_lumB <- read.table(paste0(
  PATH, "Patient_Subsets/Luminal.B_patient_ids.txt"), header = T)[,1]
patient_set_lumAB <- read.table(paste0(
  PATH, "Patient_Subsets/Luminal.A.B_patient_ids.txt"), header = T)[,1]
patient_set_basal <- read.table(paste0(
  PATH, "Patient_Subsets/Basal_patient_ids.txt"), header = T)[,1]
patient_set_her2 <- read.table(paste0(
  PATH, "Patient_Subsets/HER2_patient_ids.txt"), header = T)[,1]
patient_set_normLike <- read.table(paste0(
  PATH, "Patient_Subsets/Normal-like_patient_ids.txt"), header = T)[,1]

lumA_noMutOL <- intersect(patient_set_lumA, patients_no_mut_overlap)
lumB_noMutOL <- intersect(patient_set_lumB, patients_no_mut_overlap)
basal_noMutOL <- intersect(patient_set_basal, patients_no_mut_overlap)
her2_noMutOL <- intersect(patient_set_her2, patients_no_mut_overlap)
normLike_noMutOL <- intersect(patient_set_normLike, patients_no_mut_overlap)
lt20pamp_noMutOL <- intersect(patient_set_lt20pamp, patients_no_mut_overlap)

lt20pamp_lumA_noMutOL <- intersect(patient_set_lt20pamp, lumA_noMutOL)
lt20pamp_lumB_noMutOL <- intersect(patient_set_lt20pamp, lumB_noMutOL)
lt20pamp_basal_noMutOL <- intersect(patient_set_lt20pamp, basal_noMutOL)
lt20pamp_her2_noMutOL <- intersect(patient_set_lt20pamp, her2_noMutOL)
lt20pamp_normLike_noMutOL <- intersect(patient_set_lt20pamp, normLike_noMutOL)

# Write these subsets to files
write(lumA_noMutOL, paste0(
  PATH, "Patient_Subsets/Luminal.A.NotTP53andPIK3CAMNmut_patient_ids.txt"))
write(lumB_noMutOL, paste0(
  PATH, "Patient_Subsets/Luminal.B.NotTP53andPIK3CAMNmut_patient_ids.txt"))
write(basal_noMutOL, paste0(
  PATH, "Patient_Subsets/Basal.NotTP53andPIK3CAMNmut_patient_ids.txt"))
write(her2_noMutOL, paste0(
  PATH, "Patient_Subsets/HER2.NotTP53andPIK3CAMNmut_patient_ids.txt"))
write(normLike_noMutOL, paste0(
  PATH, "Patient_Subsets/Normal-like.NotTP53andPIK3CAMNmut_patient_ids.txt"))


############################################################
############################################################
# FIND SAMPLES WITH MUTATION IN ONE GENE & NO MUTATION
# IN ANOTHER (INDEPENDENT POSITIVE AND NEGATIVE SETS)
############################################################
############################################################
gene1 <- "TP53"
gene2 <- "PIK3CA"

#' Get patients with a mutation just in given gene1 and not in given gene2
#' @param mutation_count_matrix a mutation count matrix of all samples and genes
#' @param gene1 the gene name of the first gene (with mutation)
#' @param gene2 the gene name of the second gene (no mutation)
get_independent_mut_patient_sets <- function(mutation_count_matrix, gene1, gene2) {
  # Subset the mutation count matrix to just these two genes
  mutation_count_matrix_sub <- mutation_count_matrix[
    rownames(mutation_count_matrix) %fin% c(gene1, gene2),]
  
  # Get the patients with a mutation in the first gene but no mutation in the second
  g1_mut_g2_noMut <- mutation_count_matrix_sub[, intersect(
    which(mutation_count_matrix_sub[rownames(mutation_count_matrix_sub) == gene1,] == 1),
    (which(mutation_count_matrix_sub[rownames(mutation_count_matrix_sub) == gene2,] == 0)))]
  g1_mut_g2_noMut_pats <- unique(unlist(lapply(colnames(g1_mut_g2_noMut), function(x) 
    unlist(strsplit(x, "-", fixed = T))[3])))
  
  return(g1_mut_g2_noMut_pats)
}


g1_mut_g2_noMut_pats <- get_independent_mut_patient_sets(mutation_count_matrix, 
                                                         gene1, gene2)
g2_mut_g1_noMut_pats <- get_independent_mut_patient_sets(mutation_count_matrix, 
                                                         gene2, gene1)


#' Get and return patients with no mutation in the given gene
#' @param mutation_count_matrix a mutation count matrix of all samples and genes
#' @param gene the gene name of the gene whose patients (with a mutation in this gene)
#' we want to exclude
get_patients_with_no_mut_in_gene <- function(mutation_count_matrix, gene) {
  
  # Subset the mutation count matrix to just this gene
  mutation_count_matrix_sub <- mutation_count_matrix[
    rownames(mutation_count_matrix) == gene,]
  
  # Get the patients without a mutation in the gene
  gene_noMut <- mutation_count_matrix_sub[, which(
    mutation_count_matrix_sub[rownames(mutation_count_matrix_sub) == gene,] == 0)]
  gene_noMut_pats <- unique(unlist(lapply(colnames(gene_noMut), function(x) 
    unlist(strsplit(x, "-", fixed = T))[3])))
  
  return(gene_noMut_pats)
}

g1_noMut_pats <- get_patients_with_no_mut_in_gene(mutation_count_matrix, gene1)
g2_noMut_pats <- get_patients_with_no_mut_in_gene(mutation_count_matrix, gene2)

# Get the intersection with our intersecting ids
intersecting_pats <- read.table(paste0(OUTPUT_PATH, "intersecting_ids.txt"))[,1]

g1_noMut_pats <- intersect(g1_noMut_pats, intersecting_pats)
g2_noMut_pats <- intersect(g2_noMut_pats, intersecting_pats)

noMut_pats_list <- list(gene1 = g1_noMut_pats, gene2 = g2_noMut_pats)


# Plot a Venn diagram of the overlap
ggVennDiagram(noMut_pats_list, label_alpha = 0, 
              category.names = c(paste("No Mut,", gene1), paste("No Mut,", gene2)), 
              set_color = "black", set_size = 10, label_size = 8, edge_size = 0) +
  ggplot2::scale_fill_gradient(low="cornsilk1", high = "cadetblue3")

# Write these to files
write(g1_noMut_pats, paste0(PATH, paste0("Patient_Subsets/", 
                                         paste0(gene1, ".NoMut_patient_ids.txt"))))
write(g2_noMut_pats, paste0(PATH, paste0("Patient_Subsets/", 
                                         paste0(gene2, ".NoMut_patient_ids.txt"))))


############################################################
############################################################
# IDENTIFY MALE OR FEMALE PATIENTS 
############################################################
############################################################
# In the pan-cancer case, male is 0 and female is 1
patient_sample_df <- read.csv(paste0(
  OUTPUT_PATH, "Sample/
  combined_patient_sample_washu_inclSubtypes_IntersectPatients.csv"), 
  header = T, check.names = F, row.names = 1)

male_pats <- unique(unlist(lapply(
  rownames(patient_sample_df[patient_sample_df$Gender == 0,]), function(x)
    unlist(strsplit(x, "-", fixed = T))[1])))   
female_pats <- unique(unlist(lapply(
  rownames(patient_sample_df[patient_sample_df$Gender == 1,]), function(x)
    unlist(strsplit(x, "-", fixed = T))[1]))) 

# If needed, use a patient-cancer mapping to eliminate patients from a 
# sex-specific cancer. These include BRCA (breast), OV (ovarian), CESC (cervical), 
# PRAD (prostate), TCGT (testicular), UCS and UCEC (uterine)
sex_specific_cts <- c("BRCA", "OV", "CESC", "PRAD", "TCGT", "UCS", "UCEC")
sex_specific_pats <- unlist(patient_cancer_mapping[
  names(patient_cancer_mapping) %fin% sex_specific_cts])

male_pats <- setdiff(male_pats, sex_specific_pats)  
female_pats <- setdiff(female_pats, sex_specific_pats)  

# Write these to file
write(male_pats, paste0(PATH, "Patient_Subsets/male_patients.txt"))
write(female_pats, paste0(PATH, "Patient_Subsets/female_patients.txt"))
