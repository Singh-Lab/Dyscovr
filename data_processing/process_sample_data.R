############################################################
### PROCESS SAMPLE DATA
### PUBLICATION INFORMATION
############################################################

# This file uses a sample-level clinical data frame to create an output "sample 
# data frame" that will be used in the Dyscovr model. This data frame has the 
# following format:
  # Columns: Linear model input variables (Library size (TCGA only), tumor purity 
    # estimate, immune cell infiltration estimates)
  # Rows: Patients that have all necessary data types 
  # Entries: 0 or 1 values

library(TCGAbiolinks)
library(dplyr)


# Local PATH to directory containing sample data file
PATH <- getwd()
OUTPUT_PATH <- paste0(getwd(), "Output/Sample/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "Input Data Files/all_genes_id_conv.csv"), 
                              header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
### IMPORT ALIQUOT/ CLINICAL AND PATIENT FILES
############################################################
### Aliquot data, for TCGA ###
aliquot_df <- read.csv(paste0(PATH, "aliquot.csv"), header = T, check.names = F)

# Subset for a particular cancer type of interest, like BRCA
aliquot_df <- aliquot_df[grepl("BRCA", aliquot_df$project_id),]

# Import the patients of interest and subset the aliquot DF accordingly
patient_id_list <- read.table(paste0(PATH, "intersecting_ids.txt"), 
                              header = T)[,1]
aliquot_df$patient_id <- unlist(lapply(aliquot_df$sample_submitter_id, function(x) {
  unlist(strsplit(x, "-", fixed = T))[3]}))
aliquot_df <- aliquot_df[which(aliquot_df$patient_id %fin% patients),]

### Clinical sample DF, for METABRIC ###
clin_samp_df <- read.table(paste0(PATH, "data_clinical_sample_sub.txt"), 
                           header = T, check.names = F, sep = ",", row.names = 1)


############################################################
### IMPORT TUMOR PURITY/ CELLULARITY FILES
############################################################
# Import tumor purity data frame for TCGA
tumor_purity_df <- read.csv(paste0(PATH, "tumor_purity_estimates.csv"),
                            header = T, check.names = F)
colnames(tumor_purity_df)[1] <- 'Sample.ID'

# Import tumor cellularity for METABRIC, found in the patient clinical DF
tumor_purity_df <- read.table(paste0(PATH, "data_clinical_patient_sub.txt"), 
                              header = T, check.names = F, sep = ",")
tumor_purity_df <- tumor_purity_df[, c("PATIENT_ID", "CELLULARITY")]


############################################################
### IMPORT IMMUNE CELL INFILTRATION FILES
############################################################
# Import the pre-calculated immune cell infiltration estimates from the provided
# CSV file from the TIMER website ("Infiltration Estimation for all TCGA tumors':
# http://timer.cistrome.org/)
icd_precomp_est <- read.csv(
  paste0(PATH, "infiltration_estimation_for_tcga_TIMER.csv"), 
  header = T, check.names = F)

# Select the columns of interest from the tools of interest
TIMER_cols <- c("B cell_TIMER", "T cell CD4+_TIMER", "T cell CD8+_TIMER", 
                "Neutrophil_TIMER", "Macrophage_TIMER", 
                "Myeloid dendritic cell_TIMER")

CIBERSORT_cols <- c("B cell naive_CIBERSORT", "B cell memory_CIBERSORT", 
                    "B cell plasma_CIBERSORT", "T cell CD8+_CIBERSORT", 
                    "T cell CD4+ naive_CIBERSORT", 
                    "T cell CD4+ memory resting_CIBERSORT",
                    "T cell CD4+ memory activated_CIBERSORT", 
                    "T cell follicular helper_CIBERSORT",
                    "T cell regulatory (Tregs)_CIBERSORT", 
                    "T cell gamma delta_CIBERSORT", "NK cell resting_CIBERSORT", 
                    "NK cell activated_CIBERSORT", "Monocyte_CIBERSORT", 
                    "Macrophage M0_CIBERSORT", "Macrophage M1_CIBERSORT", 
                    "Macrophage M2_CIBERSORT", 
                    "Myeloid dendritic cell resting_CIBERSORT", 
                    "Myeloid dendritic cell activated_CIBERSORT", 
                    "MAST cell activated_CIBERSORT", "MAST cell resting_CIBERSORT",
                    "Eosinophil_CIBERSORT", "Neutrophil_CIBERSORT")

CIBERSORT_abs_cols <- c("B cell naive_CIBERSORT-ABS", "B cell memory_CIBERSORT-ABS", 
                        "B cell plasma_CIBERSORT-ABS", "T cell CD8+_CIBERSORT-ABS", 
                        "T cell CD4+ naive_CIBERSORT-ABS", 
                        "T cell CD4+ memory resting_CIBERSORT-ABS",
                        "T cell CD4+ memory activated_CIBERSORT-ABS", 
                        "T cell follicular helper_CIBERSORT-ABS",
                        "T cell regulatory (Tregs)_CIBERSORT-ABS", 
                        "T cell gamma delta_CIBERSORT-ABS",
                        "NK cell resting_CIBERSORT-ABS", 
                        "NK cell activated_CIBERSORT-ABS", 
                        "Monocyte_CIBERSORT-ABS", "Macrophage M0_CIBERSORT-ABS", 
                        "Macrophage M1_CIBERSORT-ABS", "Macrophage M2_CIBERSORT-ABS", 
                        "Myeloid dendritic cell resting_CIBERSORT-ABS", 
                        "Myeloid dendritic cell activated_CIBERSORT-ABS", 
                        "MAST cell activated_CIBERSORT-ABS", 
                        "MAST cell resting_CIBERSORT-ABS", 
                        "Eosinophil_CIBERSORT-ABS", "Neutrophil_CIBERSORT-ABS")

# For METABRIC, import one of the files we created (see 'immune_cell_infiltration.R')
icd_precomp_est <- read.csv(
  paste0(OUTPUT_PATH, "Expression/CIBERSORT/cibersort_abs_cellFrac_metabric.csv"),
  header = T, check.names = F, row.names = 1)
icd_precomp_est <- read.csv(
  paste0(OUTPUT_PATH, "Expression/CIBERSORT/cibersort_cellFrac_metabric.csv"),
  header = T, check.names = F, row.names = 1)
icd_precomp_est <- read.csv(
  paste0(OUTPUT_PATH, "Expression/quanTIseq/quantiseq_cellFrac.csv"),
  header = T, check.names = F, row.names = 1)


############################################################
# IMPORT PRE-CALLED GENOTYPE PCs (TCGA ONLY)
############################################################
# Link to file downloads: https://gdc.cancer.gov/about-data/publications/CCG-AIM-2020
# Link to associated Carrot-Zhang Paper: https://www.cell.com/cancer-cell/pdf/S1535-6108(20)30211-7.pdf

# 1. Washington University
washu_pcs <- read.table(paste0(PATH, "Genotype/WashU_PCA_ethnicity_assigned.tsv"), 
                        header = T, check.names = F, sep = "\t")
# 2. UCSF 
ucsf_pcs <- read.csv(paste0(PATH, "UCSF_Ancestry_Calls.csv"), 
                     header = T, check.names = F)

# 3. The Broad Institute
broad_pcs <- read.csv(paste0(PATH, "Broad_ancestry_PCA.txt"), 
                      header = T, check.names = F, sep = "\t")


############################################################
### CALL CREATE SAMPLE DATA FRAME FUNCTION
############################################################
# Sample Call (TCGA)
sample_dataframe <- create_sample_dataframe(aliquot_df, tumor_purity_df, 
                                            icd_precomp_est, CIBERSORT_abs_cols, 
                                            washu_pcs, "washu", "tcga")

# Sample Call (METABRIC)
sample_dataframe <- create_sample_dataframe(clin_samp_df, tumor_purity_df, 
                                            icd_precomp_est, CIBERSORT_abs_cols, 
                                            NA, NA, "metabric")

# Write results to file
write.csv(sample_dataframe, paste0(
  OUTPUT_PATH, "Dyscovr_Input/Sample/sample_dataframe_cibersort_total_frac_washu.csv"))
write.csv(sample_dataframe, paste0(
  OUTPUT_PATH, "Dyscovr_Input/Sample/sample_dataframe_cibersort_total_frac.csv"))


############################################################
### CREATE SAMPLE DATA FRAME
############################################################
#' Creates a skeleton "sample data frame" that will contain characteristics of 
#' interest for each patient sample. Calls helper functions to fill patient data 
#' frame and then returns it.
#' @param aliquot_df an aliquot data file from the TCGA for the given patient 
#' cohort (or, alternatively, a clinical sample data file from METABRIC)
#' @param tumor_purity_df an estimate of tumor purity per sample
#' @param immune_cell_infl_df an estimate of immune cell infiltration, per sample
#' @param ici_columns a vector of strings indicating the column names of the 
#' tool of interest for the immune cell infiltration DF
#' @param genotype_pc_file a file with pre-called genotype PCs from an external source
#' (either the Broad, UCSF, or WashU). NA if not using.
#' @param genotype_pc_file_label either "broad", "ucsf", or "washu" to indicate the 
#' source of the file and dictate how we parse it. NA if not using.
#' @param num_geno_pcs the number of genotype PCs we wish to include (default is 3)
#' @param dataset either 'tcga' or 'metabric'to indicate what column names to 
#' reference/ data types to include
create_sample_dataframe <- function(aliquot_df, tumor_purity_df, 
                                    immune_cell_infl_df, ici_columns,  
                                    genotype_pc_file, genotype_pc_file_label, 
                                    dataset, num_geno_pcs = 3) {
  
  # Create relational data frame to hold patient-dependent inputs to Dyscovr
  # Columns: Dyscovr input variables ; Rows: Samples that have all necessary 
  # data types (Unique 4-digit ID)
  unique_sample_ids <- c()
  if(dataset == 'tcga') {
    unique_sample_ids <- unique(aliquot_df$sample_submitter_id)
  } else if (dataset == 'metabric') {
    unique_sample_ids <- unique(aliquot_df$SAMPLE_ID)
  } else if (dataset == "chinese_tn") {
    unique_sample_ids <- unique(aliquot_df$Project_ID)
    tumor_purity_df <- aliquot_df
  } else {print("Only implemented for TCGA and METABRIC.")}
  
  # Create starter DF with tumor purity data
  sample_characteristics <- c("Tumor_purity_b1", "Tumor_purity_b2", 
                              "Tumor_purity_b3")
  sample_dataframe <- data.frame(matrix(ncol = length(sample_characteristics), 
                                        nrow = length(unique_sample_ids)))
  colnames(sample_dataframe) <- sample_characteristics
  rownames(sample_dataframe) <- unique_sample_ids
  
  # Input all the patient-specific information 
  sample_dataframe <- input_sample_specific_info(sample_dataframe, aliquot_df, 
                                                 tumor_purity_df, 
                                                 immune_cell_infl_df, 
                                                 ici_columns, genotype_pc_file, 
                                                 genotype_pc_file_label, 
                                                 num_geno_pcs, dataset)
  
  # Remove any rows that are entirely NA
  sample_dataframe <- sample_dataframe[rowSums(is.na(sample_dataframe)) != 
                                         ncol(sample_dataframe),]
  
  return(sample_dataframe)
}

############################################################

#' Takes the empty sample data frame, along with the aliquot data frame, and 
#' updates & returns sample DF with all the sample-specific information, e.g.
#' tumor purity and immune cell infiltration. 
#' @param input_df a starter data frame with the dimensions of the number
#' of samples in the cohort (rows)
#' @param aliquot_df aliquot data file from the TCGA for the given patient cohort
#' @param tumor_purity_df an estimate of tumor purity per sample
#' @param immune_cell_infl_df an estimate of immune cell infiltration per sample
#' @param ici_columns a vector of strings indicating the column names of the 
#' tool of interest for the immune cell infiltration DF
#' @param genotype_pc_file a file with pre-called genotype PCs from an external 
#' source (either the Broad, UCSF, or WashU). NA if not using.
#' @param genotype_pc_file_label either "broad", "ucsf", or "washu" to indicate 
#' the source of the file and dictate how we parse it. NA if not using.
#' @param num_geno_pcs the # of genotype PCs we wish to include (default is 3)
#' @param dataset either 'tcga' or 'metabric' to indicate what column
#' names to reference/ data types to include
input_sample_specific_info <- function(input_df, aliquot_df, tumor_purity_df, 
                                       immune_cell_infl_df, ici_columns, 
                                       genotype_pc_file, genotype_pc_file_label, 
                                       num_geno_pcs, dataset) {
  input_dataframe_updated <- input_df
  sample_ids <- rownames(input_df)
  
  # TUMOR PURITY
  if(dataset == "tcga") {
    tumor_purity_df <- convert_tumor_purity_to_cpe(
      rownames(input_dataframe_updated), tumor_purity_df)
  } 
  tumor_purity_df_conv_bucketed <- bucket_tumor_purity(tumor_purity_df, dataset)
  input_dataframe_updated[,1:3] <- tumor_purity_df_conv_bucketed
  
  # IMMUNE CELL INFILTRATION
  # OPTION 1: INCLUDE TOTAL PERCENTAGE OF IMMUNE CELLS AS A COVARIATE 
  # (need absolute fractions from Cibersort Abs) -- default method
  input_dataframe_updated <- add_total_immune_cell_percentage(
    input_dataframe_updated, immune_cell_infl_df, ici_columns, dataset)
  
  # OPTION 2: INCLUDE THE FRACTION OF EACH IMMUNE CELL TYPE AS INDIVIDUAL 
  # COVARIATES -- uncomment to replace option 1
  #input_dataframe_updated <- add_immune_cell_fractions(input_dataframe_updated, 
  #immune_cell_infl_df, ici_columns, dataset)
  
  # PRE-CALLED GENOTYPE PCs
  if(length(genotype_pc_file) > 1) {
    genotype_pc_file_processed <- process_genotype_pc_file(genotype_pc_file, 
                                                           genotype_pc_file_label, 
                                                           num_geno_pcs, sample_ids)
    input_dataframe_updated <- merge(input_dataframe_updated, 
                                     genotype_pc_file_processed, 
                                     by = 'row.names', all = T) 
    rownames(input_dataframe_updated) <- input_dataframe_updated$Row.names
    input_dataframe_updated <- input_dataframe_updated[,setdiff(colnames(
      input_dataframe_updated), "Row.names")]
  }
  
  return(input_dataframe_updated)
}

############################################################ 

#' Takes in the tumor purity estimates file from Aran et al. 2015 and 
#' for each sample, extracts either the CPE value (if available) or
#' the median of the other estimates (if CPE not available). 
#' @param sample_ids a vector of all the sample IDs of interest
#' @param tumor_purity_df a data frame from Arun et al., with columns
#' for the sample ID, the cancer ID, four estimates from various tumor
#' purity algorithms, and a CPE column with a combined purity estimate.
convert_tumor_purity_to_cpe <- function(sample_ids, tumor_purity_df) {
  
  # For each sample, get its tumor purity 
  cpe_values <- unlist(lapply(sample_ids, function(id) {
    if (id %fin% tumor_purity_df$Sample.ID) {
      
      sample_sub <- tumor_purity_df[tumor_purity_df$Sample.ID == id,]
      
      # If CPE is there, use this value
      if (!(is.nan(sample_sub$CPE))) {return(sample_sub$CPE)}
      
      # If CPE is NaN, use the median of the other estimates
      else {
        vals <- unlist(sample_sub[,3:6])
        vals <- vals[!is.nan(vals)]
        return(median(vals))
      }
    } else {return(NA)}
  }))
  
  cpe_df <- data.frame(cpe_values, row.names = sample_ids)
  colnames(cpe_df) <- 'Tumor.purity'
  
  return(cpe_df)
}

############################################################ 

#' Takes in the tumor purity CPE estimates data frame, created
#' from the file from Aran et al. 2015. For each sample, it
#' converts the CPE value to a bucketed value: low purity (0.0-0.5), 
#' medium purity (0.5-0.75), and high purity (0.75-1.0). METABRIC creates
#' 3 buckets based on "low", "medium", and "high" estimates of tumor purity.
#' @param tumor_purity_df a data frame produced from the
#' convert_tumor_purity_to_cpe() function (row names are sample TCGA IDs,
#' Tumor.purity is a column of the CPE purity value) for TCGA, or a 
#' sample clinical data frame for METABRIC
#' @param dataset either 'tcga' or 'metabric' to indicate what column names to 
#' reference/ data types to include
bucket_tumor_purity <- function(tumor_purity_df, dataset) {
  
  # For each TCGA sample, get its tumor purity and classify it into a bucket
  tumor_purity_df_bucket <- data.frame()
  if(dataset == "tcga") {
    bucketed_rows <- lapply(tumor_purity_df$Tumor.purity, function(tp) {
      if(is.na(tp)) {return(c(NA,NA,NA))}
      else {
        if (tp <= 0.5) {return(c(1,0,0))} 
        else if ((tp > 0.5) & (tp <= 0.75)) {return(c(0,1,0))} 
        else if (tp > 0.75) {return(c(0,0,1))} 
        else {
          print(paste("Inappropriate purity value:", tp))
          return(c(NA,NA,NA))
        }
      }
    })
    
  } else if (dataset == "metabric") {
    bucketed_rows <- lapply(tumor_purity_df$CELLULARITY, function(x) {
      if(x == "Low") {return(c(1,0,0))}
      else if (x == "Moderate") {return(c(0,1,0))}
      else if (x == "High") {return(c(0,0,1))}
      else {print(paste("Inappropriate purity value:", x)); return(c(NA,NA,NA))}
    })
  } else {print("Only implemented for TCGA and METABRIC.")}
  
  tumor_purity_df_bucket <- do.call(rbind, bucketed_rows)
  rownames(tumor_purity_df_bucket) <- rownames(tumor_purity_df)
  colnames(tumor_purity_df_bucket) <- c('Tumor_purity_b1', 'Tumor_purity_b2',
                                        'Tumor_purity_b3')
  
  #print(tumor_purity_df_bucket)
  return(tumor_purity_df_bucket)
}


############################################################ 
### OPTION 1: GET THE TOTAL IMMUNE CELL PERCENTAGE FOR EACH SAMPLE
#' For each sample, gets the absolute percentage of immune cells within that 
#' sample. Return input data frame with, for each row (sample), a column with the 
#' absolute fraction of immune cells from CIBERSORT Abs, bucketed. Groups it
#' into one of three buckets: low immune cell infiltration (0-0.3), medium immune
#' cell infiltration (0.3-0.7), high immune cell infiltration (0.7-1.0).
#' @param sample_df a sample data frame for which to add immune cell 
#' infiltration estimates 
#' @param immune_cell_infl_df a data frame of immune cell infiltration estimates
#' @param ici_columns a vector of the columns of the tool of interest (e.g. 
#' CIBERSORT Abs)
#' @param dataset either 'tcga', or 'metabric' to indicate what column names to 
#' reference/ data types to include
add_total_immune_cell_percentage <- function(sample_df, immune_cell_infl_df, 
                                             ici_columns, dataset) {
  
  immune_cell_est <- lapply(rownames(sample_df), function(sample) {
    
    # Remove the final letter to match, if TCGA
    if (dataset == "tcga") {sample <- substr(sample, 1, nchar(sample)-1)}
    
    # For this sample, get a subset of the IC infiltration data frame with the 
    # columns of interest from the tool of interest
    ici_df_sub <- c()
    total_ic_frac <- NA
    
    if(dataset == 'tcga') {
      if (sample %fin% immune_cell_infl_df$cell_type) {
        ici_df_sub <- immune_cell_infl_df[
          grepl(sample, immune_cell_infl_df$cell_type),
          colnames(immune_cell_infl_df) %fin% ici_columns]
      } else {return (c(NA, NA, NA))}
      # Sum the fraction of all these non-tumor cell types together
      total_ic_frac <- sum(as.numeric(ici_df_sub), na.rm = T)
    }
    else if (dataset == 'metabric') {
      if (sample %fin% colnames(immune_cell_infl_df)) {
        ici_df_sub <- immune_cell_infl_df[, grepl(
          sample, colnames(immune_cell_infl_df))]
      } else {return (c(NA, NA, NA))}
      # Sum the fraction of all these non-tumor cell types together
      total_ic_frac <- sum(as.numeric(ici_df_sub), na.rm = T)
    }
    else {print("Only implemented for TCGA and METABRIC.")}
    
    print(total_ic_frac)
    thresh <- c(0.3, 0.7)
    if(dataset == "metabric") {thresh <- c(1.1, 1.3)}
    
    # Return this value, bucketed 
    if(total_ic_frac <= thresh[1]) {return(c(1,0,0))}
    else if ((total_ic_frac > thresh[1]) & (total_ic_frac <= thresh[2])) {
      return(c(0,1,0))
    }
    else if (total_ic_frac > thresh[2]) {return(c(0,0,1))}
    else {return(c(NA, NA, NA))}
  })
  
  immune_cell_est_df <- do.call(rbind, immune_cell_est)
  colnames(immune_cell_est_df) <- c("Tot_IC_Frac_b1", "Tot_IC_Frac_b2", 
                                    "Tot_IC_Frac_b3")
  
  sample_df <- cbind(sample_df, immune_cell_est_df)
  
  return(sample_df)
}

############################################################ 
### OPTION 2: GET THE FRACTION OF EACH CELL TYPE IN EACH SAMPLE, WITH A COVARIATE 
### FOR EACH
#' For each sample, get the immune cell fractions for all cell types of interest.
#' Return an input data frame with, for each row (sample), has additional columns
#' with the fraction of each cell type in that sample.
#' @param sample_df a sample data frame for which to add immune cell 
#' infiltration estimates 
#' @param immune_cell_infl_df a data frame of immune cell infiltration estimates
#' @param ici_columns a vector of the columns of the tool of interest (e.g. 
#' CIBERSORT Abs)
add_immune_cell_fractions <- function(sample_df, immune_cell_infl_df, 
                                      ici_columns) {
  
  immune_cell_est <- lapply(rownames(sample_df), function(sample) {
    
    # Remove the final letter to match
    sample <- substr(sample, 1, nchar(sample)-1)
    
    # For this sample, get a subset of the IC infiltration data frame with 
    # the columns of interest from the tool of interest
    if(sample %fin% immune_cell_infl_df$cell_type) {
      ici_df_sub <- immune_cell_infl_df[
        grepl(sample, immune_cell_infl_df$cell_type),
        colnames(immune_cell_infl_df) %fin% ici_columns]
    } else { return(NA)}
    
    # Return just this subset (1 row x length(ici_columns) dimensional DF)
    return(ici_df_sub)
  })
  
  # Bind all these rows together
  ici_dataframe <- do.call(rbind, immune_cell_est)
  sample_df <- cbind(sample_df, ici_dataframe)
  
  return(sample_df)
}


############################################################ 
#' Takes a pre-processed genotype PC file, from either the Broad, UCSF, or
#' WashU, and reformats it to have the desired number of PCs in the correct 
#' format to add to the patient data frame
#' @param genotype_pc_file the file with patient IDs and associated genotype PCs
#' @param genotype_pc_file_label either "broad", "ucsf", or "washu" to denote
#' what the file format is
#' @param num_pcs the number of pcs we would like to record in the patient file
#' (default is 3)
#' @param clinical_df a clinical DF with sample names from the TCGA
process_genotype_pc_file <- function(genotype_pc_file, genotype_pc_file_label, 
                                     num_pcs, sample_ids) {
  
  outfile <- NA
  
  # Process each file type differently based on its format
  if(genotype_pc_file_label == "washu") {
    genotype_pc_file$Sample <- unlist(lapply(genotype_pc_file$Sample, function(x)
      paste(unlist(strsplit(x, "-", fixed = T))[1:4], collapse = "-")))
    genotype_pc_file_sub <- distinct(genotype_pc_file[
      genotype_pc_file$Sample %fin% sample_ids, 
      c("Sample", paste0("PC", 1:num_pcs))])
    outfile <- genotype_pc_file_sub
    
    # Handle any duplicates (e.g., two plates for the same sample) by 
    # averaging the PCs
    duplicates <- unlist(outfile$Sample[which(duplicated(outfile$Sample))])
    for (d in duplicates) {
      outfile_d <- outfile[outfile$Sample == d,]
      outfile_d_colMeans <- colMeans(outfile_d[grepl("PC", colnames(outfile_d))])
      new_row <- outfile_d[1, , drop = F]
      new_row[1,grepl("PC", colnames(new_row))] <- as.numeric(
        unlist(outfile_d_colMeans))
      outfile <- outfile[outfile$Sample != d, ]
      outfile <- rbind(outfile, new_row)
    }
    
    rownames(outfile) <- outfile$Sample
    outfile <- outfile[,setdiff(colnames(outfile), "Sample")]
    
    outfile <- add_missing_ids(outfile, sample_ids)
    
  } else if (genotype_pc_file_label == "ucsf") {
    genotype_pc_file$Sample <- unlist(lapply(genotype_pc_file$Aliquot_ID, function(x)
      paste(unlist(strsplit(x, "-", fixed = T))[1:4], collapse = "-")))
    genotype_pc_file_sub <- genotype_pc_file[genotype_pc_file$Sample %fin% sample_ids, 
                                             c("Sample", paste0("PC", 1:num_pcs))]
    outfile <- genotype_pc_file_sub
    rownames(outfile) <- outfile$Sample
    outfile <- outfile[,setdiff(colnames(outfile), "Sample")]
    
    # Use helper function to add missing IDs
    outfile <- add_missing_ids(outfile, sample_ids)
    
  } else if (genotype_pc_file_label == "broad") {
    # Split apart the PCs first
    genotype_pc_file$PC1 <- unlist(lapply(genotype_pc_file$PC1:PC2:PC3, function(x) 
      unlist(strsplit(x, ":", fixed = T))[1]))
    genotype_pc_file$PC1 <- unlist(lapply(genotype_pc_file$PC1:PC2:PC3, function(x) 
      unlist(strsplit(x, ":", fixed = T))[2]))
    genotype_pc_file$PC1 <- unlist(lapply(genotype_pc_file$PC1:PC2:PC3, function(x) 
      unlist(strsplit(x, ":", fixed = T))[3]))
    
    patient_ids <- unlist(lapply(sample_ids, function(x) 
      paste(unlist(strsplit(x, "-", fixed = T))[1:3], collapse = "-")))
    
    genotype_pc_file_sub <- NA
    if(num_pcs > 3) {
      print("Error - only 3 PCs provided for Broad.")
      genotype_pc_file_sub <- genotype_pc_file[
        genotype_pc_file$SampleID %fin% patient_ids, 
        c("SampleID", paste0("PC", 1:3))]
      
    } else {  #UCSF
      genotype_pc_file_sub <- genotype_pc_file[
        genotype_pc_file$SampleID %fin% patient_ids, 
        c("SampleID", paste0("PC", 1:num_pcs))]
    }
    outfile <- genotype_pc_file_sub
    rownames(outfile) <- outfile$SampleID
    outfile <- outfile[,setdiff(colnames(outfile), "SampleID")]
    
    outfile <- add_missing_ids(outfile, sample_ids)
    
  } else {
    print("Error. Only implemented for WashU, UCSF, and the Broad input files. 
          Returning NA.")
    return(data.frame("PC1" = rep(NA, times = )))
  }
  
  return(outfile)
}


#' Copies over the PC values for normal samples to the corresponding tumor 
#' sample, if missing
#' @param outfile a file with row names as samples with data, and columns for PCs
#' @param sample_ids a vector of all the sample IDs
add_missing_ids <- function(outfile, sample_ids) {
  missing_ids <- setdiff(sample_ids, rownames(outfile))
  #print(paste("Missing IDs:", head(missing_ids)))
  
  new_rows <- lapply(missing_ids, function(id) {
    id_trim <- unlist(strsplit(id, "-", fixed = T))[3]
    # Check if there is another entry for this patient with data
    if(T %fin% grepl(id_trim, rownames(outfile))) {
      row <- outfile[grepl(id_trim, rownames(outfile)), , drop = F]
      if(nrow(row) > 1) {row <- row[1,]}
      return(row)
    } else {return(NA)}
  })
  new_df <- do.call(rbind, new_rows)
  rownames(new_df) <- missing_ids
  
  outfile <- rbind(outfile, new_df)
  
  return(outfile)
}


############################################################ 
#' Takes a main data frame for all sample characteristics, along with a sub-data 
#' frame, and merges the two by their column names
#' @param main_df the patient data frame we are constructing
#' @param partial_df a partial data frame with matching column names we'd like 
#' to merge with our main patient data frame
merge_dataframes_by_colname <- function(main_df, partial_df) {
  
  # Ensure that the patients are in the same order
  partial_df <- partial_df[order(match(rownames(partial_df), 
                                       rownames(main_df))),]
  
  # Replace the empty columns with the new ones
  for (i in 1:ncol(partial_df)) {
    pos <- as.numeric(which(colnames(partial_df)[i] == colnames(main_df)))
    main_df[,pos] <- partial_df[,i]
  }
  
  return(main_df)
}
