############################################################
### Create Sample Data Frame
### Written By: Sara Camilli, June 2021
############################################################

# This file uses a sample data frame, and a variable declaring whether the data 
  # of interest is or is not breast cancer, to create an output "sample data frame" 
  # that will be used in the linear  model. This data frame has the following format:
    # Columns: Linear model input variables (Library size, tumor purity estimate,
    # immune cell infiltration estimates, and PEER factors)
    # Rows: Patients that have all necessary data types 1...j...L (Unique 4-digit ID)
    # Entries: 0 or 1 values

library(TCGAbiolinks)

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

input_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/"
#input_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/"


is_brca <- TRUE
# is_brca <- FALSE


############################################################
### IMPORT ALIQUOT AND PATIENT FILES
############################################################
# Aliquot data
aliquot_df <- read.csv(paste(input_path, "aliquot.csv", sep = ""), header = TRUE, check.names = FALSE)
if(is_brca) {
  aliquot_df <- aliquot_df[grepl("BRCA", aliquot_df$project_id),]
}

# Import the patients of interest and subset the aliquot DF to only these patients
patients <- read.table(paste(input_path, "unique_brca_patient_ids_2.txt", sep = ""), header = TRUE,
                     row.names = 1)[,1]
# patients <- read.table(paste(input_path, "unique_patient_ids_2.txt", sep = ""), header = TRUE,
  # row.names = 1)[,1]

aliquot_df$patient_id <- unlist(lapply(aliquot_df$sample_submitter_id, function(x) {
  unlist(strsplit(x, "-", fixed = TRUE))[3]}))
aliquot_df <- aliquot_df[which(aliquot_df$patient_id %in% patients),]

############################################################
### IMPORT EXPRESSION FILES
############################################################
# Import expression counts data frame, or samples data frame from edgeR normalization
expression_df <- read.csv(paste(main_path, "Expression/expression_counts_DF.csv", sep = ""), 
                          header = TRUE, check.names = FALSE)
# FOR BRCA
tmm_norm_samples_df <- read.csv(paste(main_path, "Expression/tmm_normalized_expression_samples.csv", sep = ""),
                                header = TRUE, check.names = FALSE)
# FOR PAN-CANCER
tmm_norm_samples_df <- read.csv(paste(main_path, "Expression/TMM/ALL_CT_tmm_normalized_expression_samples.csv", sep = ""),
                                header = TRUE, check.names = FALSE)

# Import tumor purity data frame
tumor_purity_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/tumor_purity_estimates.csv",
                            header = TRUE, check.names = FALSE)
# if needed:
colnames(tumor_purity_df)[1] <- 'Sample.ID'

############################################################
### IMPORT PEER FILES
############################################################
# TMM
peer_factors_tmm <- read.csv(paste(main_path, "PEER/TMM/peer_factors (TMM, noCov, noMean, allGenes).csv", sep = ""),
                             header = TRUE, row.names = 1, check.names = FALSE)

# FPKM
peer_factors_fpkm <- read.csv(paste(main_path, "PEER/FPKM/peer_factors (FPKM, noCov, noMean, allGenes).csv", sep = ""),
                              header = TRUE, row.names = 1, check.names = FALSE)
peer_factors_fpkm_top10k <- read.csv(paste(main_path, "PEER/FPKM/peer_factors (FPKM, noCov, noMean, top10k).csv", sep = ""),
                              header = TRUE, row.names = 1, check.names = FALSE)

# Quantile-Normalized
peer_factors_qn_top10k <- read.csv(paste(main_path, "PEER/Q-N/TiebreakerMean/peer_factors (Q-N, noCov, noMean, top10k).csv", sep = ""),
                              header = TRUE, row.names = 1, check.names = FALSE)

# Rank-Normalized
peer_factors_rn <- read.csv(paste(main_path, "PEER/R-N/TiebreakerMean/peer_factors (R-N_mean_MedGr10_TO, noCov, noMean, allGenes).csv", sep = ""),
                              header = TRUE, row.names = 1, check.names = FALSE)
peer_factors_rn_top10k <- read.csv(paste(main_path, "PEER/R-N/TiebreakerMean/peer_factors (R-N_mean_MedGr10_TO, noCov, noMean, top10k).csv", sep = ""),
                            header = TRUE, row.names = 1, check.names = FALSE)

############################################################
### IMPORT IMMUNE CELL INFILTRATION FILES
############################################################
# Import the pre-calculated immune cell infiltration estimates from the provided
# CSV file from the TIMER website
icd_precomp_est <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/infiltration_estimation_for_tcga_TIMER.csv",
                              header = TRUE, check.names = FALSE)

# Select the columns of interest from the tools of interest
TIMER_cols <- c("B cell_TIMER", "T cell CD4+_TIMER", "T cell CD8+_TIMER", 
                "Neutrophil_TIMER", "Macrophage_TIMER", "Myeloid dendritic cell_TIMER")

CIBERSORT_cols <- c("B cell naive_CIBERSORT", "B cell memory_CIBERSORT", 
                    "B cell plasma_CIBERSORT", "T cell CD8+_CIBERSORT", 
                    "T cell CD4+ naive_CIBERSORT", "T cell CD4+ memory resting_CIBERSORT",
                    "T cell CD4+ memory activated_CIBERSORT", "T cell follicular helper_CIBERSORT",
                    "T cell regulatory (Tregs)_CIBERSORT", "T cell gamma delta_CIBERSORT",
                    "NK cell resting_CIBERSORT", "NK cell activated_CIBERSORT", 
                    "Monocyte_CIBERSORT", "Macrophage M0_CIBERSORT", 
                    "Macrophage M1_CIBERSORT", "Macrophage M2_CIBERSORT", 
                    "Myeloid dendritic cell resting_CIBERSORT", 
                    "Myeloid dendritic cell activated_CIBERSORT", 
                    "MAST cell activated_CIBERSORT", "MAST cell resting_CIBERSORT",
                    "Eosinophil_CIBERSORT", "Neutrophil_CIBERSORT")

CIBERSORT_abs_cols <- c("B cell naive_CIBERSORT-ABS", "B cell memory_CIBERSORT-ABS", 
                        "B cell plasma_CIBERSORT-ABS", "T cell CD8+_CIBERSORT-ABS", 
                        "T cell CD4+ naive_CIBERSORT-ABS", "T cell CD4+ memory resting_CIBERSORT-ABS",
                        "T cell CD4+ memory activated_CIBERSORT-ABS", "T cell follicular helper_CIBERSORT-ABS",
                        "T cell regulatory (Tregs)_CIBERSORT-ABS", "T cell gamma delta_CIBERSORT-ABS",
                        "NK cell resting_CIBERSORT-ABS", "NK cell activated_CIBERSORT-ABS", 
                        "Monocyte_CIBERSORT-ABS", "Macrophage M0_CIBERSORT-ABS", 
                        "Macrophage M1_CIBERSORT-ABS", "Macrophage M2_CIBERSORT-ABS", 
                        "Myeloid dendritic cell resting_CIBERSORT-ABS", 
                        "Myeloid dendritic cell activated_CIBERSORT-ABS", 
                        "MAST cell activated_CIBERSORT-ABS", "MAST cell resting_CIBERSORT-ABS",
                        "Eosinophil_CIBERSORT-ABS", "Neutrophil_CIBERSORT-ABS")


############################################################
### CALL CREATE SAMPLE DATA FRAME FUNCTION
############################################################
# CIBERSORT abs
sample_dataframe <- create_sample_dataframe(aliquot_df, expression_df, tumor_purity_df, 
                                            icd_precomp_est, CIBERSORT_abs_cols, peer_factors_tmm)
sample_dataframe <- create_sample_dataframe(aliquot_df, expression_df, tumor_purity_df, 
                                            icd_precomp_est, CIBERSORT_abs_cols, peer_factors_fpkm)
sample_dataframe <- create_sample_dataframe(aliquot_df, expression_df, tumor_purity_df, 
                                            icd_precomp_est, CIBERSORT_abs_cols, peer_factors_fpkm_top10k)
sample_dataframe <- create_sample_dataframe(aliquot_df, expression_df, tumor_purity_df, 
                                            icd_precomp_est, CIBERSORT_abs_cols, peer_factors_qn_top10k)
sample_dataframe <- create_sample_dataframe(aliquot_df, expression_df, tumor_purity_df, 
                                            icd_precomp_est, CIBERSORT_abs_cols, peer_factors_rn)
sample_dataframe <- create_sample_dataframe(aliquot_df, expression_df, tumor_purity_df, 
                                            icd_precomp_est, CIBERSORT_abs_cols, peer_factors_rn_top10k)

# TIMER
sample_dataframe <- create_sample_dataframe(aliquot_df, expression_df, tumor_purity_df, 
                                            icd_precomp_est, TIMER_cols, peer_factors_tmm)
sample_dataframe <- create_sample_dataframe(aliquot_df, expression_df, tumor_purity_df, 
                                            icd_precomp_est, TIMER_cols, peer_factors_fpkm)
sample_dataframe <- create_sample_dataframe(aliquot_df, expression_df, tumor_purity_df, 
                                            icd_precomp_est, TIMER_cols, peer_factors_fpkm_top10k)
sample_dataframe <- create_sample_dataframe(aliquot_df, expression_df, tumor_purity_df, 
                                            icd_precomp_est, TIMER_cols, peer_factors_qn_top10k)
sample_dataframe <- create_sample_dataframe(aliquot_df, expression_df, tumor_purity_df, 
                                            icd_precomp_est, TIMER_cols, peer_factors_rn)
sample_dataframe <- create_sample_dataframe(aliquot_df, expression_df, tumor_purity_df, 
                                            icd_precomp_est, TIMER_cols, peer_factors_rn_top10k)

write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_timer_tmm.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_timer_fpkm.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_timer_fpkm_top10k.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_timer_qn_top10k.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_timer_rn.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_timer_rn_top10k.csv", sep = ""))

write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_tmm.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_fpkm.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_fpkm_top10k.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_qn_top10k.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_rn.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_rn_top10k.csv", sep = ""))

write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_abs_tmm.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_abs_fpkm.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_abs_fpkm_top10k.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_abs_qn_top10k.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_abs_rn.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_abs_rn_top10k.csv", sep = ""))

write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_total_frac_tmm.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_total_frac_fpkm.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_total_frac_fpkm_top10k.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_total_frac_qn_top10k.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_total_frac_rn.csv", sep = ""))
write.csv(sample_dataframe, paste(main_path, "Linear Model/Patient and Sample DFs/sample_dataframe_cibersort_total_frac_rn_top10k.csv", sep = ""))


############################################################
### CREATE SAMPLE DATA FRAME
############################################################
#' Creates a skeleton "sample data frame" that will contain
#' all the characteristics of interest for each patient sample. Calls helper 
#' function to fill patient data frame and then returns it.
#' @param aliquot_df an aliquot data file from the TCGA for the given patient cohort
#' @param expression_df a counts expression data frame, used to add patient library sizes 
#' @param tumor_purity_df an estimate of tumor purity per sample
#' @param immune_cell_infl_df an estimate of immune cell infiltration for the TIMER 
#' tool, per sample
#' @param ici_columns a vector of strings indicating the column names of the tool of
#' interest for the immune cell infiltration DF
create_sample_dataframe <- function(aliquot_df, expression_df, tumor_purity_df, 
                                    immune_cell_infl_df, ici_columns, peer_df) {
  
  # Create relational data frame to hold patient-dependent inputs to linear model
  # Columns: Linear model input variables 
  # Rows: Samples that have all necessary data types 1...j...L (Unique 4-digit ID)
  sample_characteristics <- c("Tumor_purity_b1", "Tumor_purity_b2", "Tumor_purity_b3")
  unique_sample_ids <- unique(aliquot_df$sample_submitter_id)
  sample_dataframe <- data.frame(matrix(ncol = length(sample_characteristics), 
                                         nrow = length(unique_sample_ids)))
  colnames(sample_dataframe) <- sample_characteristics
  rownames(sample_dataframe) <- unique_sample_ids
  
  # Input all the patient-specific information 
  sample_dataframe <- input_sample_specific_info(sample_dataframe, aliquot_df, 
                                                 expression_df, tumor_purity_df,
                                                 immune_cell_infl_df, ici_columns,
                                                 peer_df)
  
  # Remove any rows that are entirely NA
  sample_dataframe <- sample_dataframe[rowSums(is.na(sample_dataframe)) != ncol(sample_dataframe),]
  
  return(sample_dataframe)
}

############################################################

#' Takes the empty sample data frame, along with the aliquot data frame and 
#' whether or not we are looking at BRCA, and updates sample DF with all the 
#' sample-specific information: total number of mutations, tumor purity,  
#' immune cell infiltration, and library sizes. 
#' Returns data frame with these fields filled.
#' @param input_df a blank data frame with the dimensions of the the number
#' of sample characteristics (columns) and number of samples in the 
#' cohort (rows)
#' @param aliquot_df aliquot data file from the TCGA for the given patient cohort
#' @param mut_count_matrix a mutation count data frame created from refiguring of 
#' maftools results
#' @param expression_df a counts expression data frame, used to add patient library sizes 
#' @param tumor_purity_df an estimate of tumor purity per sample
#' @param immune_cell_infl_df an estimate of immune cell infiltration for the TIMER 
#' tool, per sample
#' @param ici_columns a vector of strings indicating the column names of the tool of
#' interest for the immune cell infiltration DF
input_sample_specific_info <- function(input_df, aliquot_df, expression_df, 
                                       tumor_purity_df, immune_cell_infl_df,
                                       ici_columns, peer_df) {
  input_dataframe_updated <- input_df
  sample_ids <- rownames(input_df)
  
  # TUMOR PURITY
  tumor_purity_df_conv <- convert_tumor_purity_to_cpe(rownames(input_dataframe_updated), 
                                                      tumor_purity_df)
  tumor_purity_df_conv_bucketed <- bucket_tumor_purity(tumor_purity_df_conv)
  input_dataframe_updated[,1:3] <- tumor_purity_df_conv_bucketed
  
  
  # IMMUNE CELL INFILTRATION
  # OPTION 1: INCLUDE TOTAL PERCENTAGE OF IMMUNE CELLS AS A COVARIATE (need absolute fractions)
  #input_dataframe_updated <- add_total_immune_cell_percentage(input_dataframe_updated, 
                                                              #immune_cell_infl_df,
                                                              #ici_columns)
  
  # OPTION 2: INCLUDE THE FRACTION OF EACH IMMUNE CELL TYPE AS INDIVIDUAL COVARIATES
  input_dataframe_updated <- add_immune_cell_fractions(input_dataframe_updated, 
                                                      immune_cell_infl_df, ici_columns)
  
  # PEER FACTORS
  input_dataframe_updated <- add_peer(input_dataframe_updated, peer_df)
  
  # LIBRARY SIZES
  input_dataframe_updated <- add_library_sizes(input_dataframe_updated, expression_df)
  
  
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
    if (id %in% tumor_purity_df$Sample.ID) {
      
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
  
  #print(cpe_df)
  
  return(cpe_df)
}

############################################################ 

#' Takes in the tumor purity CPE estimates data frame, created
#' from the file from Aran et al. 2015. For each sample, it
#' converts the CPE value to a bucketed value: low purity (0.0-0.5), 
#' medium purity (0.5-0.75), and high purity (0.75-1.0)
#' @param tumor_purity_df a data frame produced from the
#' convert_tumor_purity_to_cpe function (row names are sample TCGA IDs,
#' Tumor.purity is a column of the CPE purity value)
bucket_tumor_purity <- function(tumor_purity_df) {
  
  # For each sample, get its tumor purity and classify it into a bucket
  bucketed_rows <- lapply(tumor_purity_df$Tumor.purity, function(tp) {
    if(is.na(tp)) {return(c(NA,NA,NA))}
    else {
      if (tp <= 0.5) {return(c(1,0,0))} 
      else if ((tp > 0.5) & (tp <= 0.75)) {return(c(0,1,0))} 
      else if (tp > 0.75) {return(c(0,0,1))} 
      else {print(paste("Inappropriate purity value:", tp))}
    }
    
  })
  
  # Recombine these rows into a data frame
  tumor_purity_df_bucket <- do.call(rbind, bucketed_rows)
  colnames(tumor_purity_df_bucket) <- c('Tumor_purity_b1', 'Tumor_purity_b2',
                                        'Tumor_purity_b3')
  rownames(tumor_purity_df_bucket) <- rownames(tumor_purity_df)
  
  print(tumor_purity_df_bucket)
  
  return(tumor_purity_df_bucket)
}


############################################################ 
### OPTION 1: GET THE TOTAL IMMUNE CELL PERCENTAGE FOR EACH SAMPLE
# TODO: do we include all cell types, or just infiltrating ones (like macrophages
# and CD8+ T-Cells?) 
#' For each sample, gets the absolute percentage of immune cells within that sample.
#' Return input data frame with, for each row (sample), a column with the 
#' absolute fraction of immune cells from CIBERSORT Abs, bucketed. Groups it
#' into one of three buckets: low immune cell infiltration (0-0.3), medium immune
#' cell infiltration (0.3-0.7), high immune cell infiltration (0.7-1.0).
#' @param sample_df a sample data frame to add immune cell infiltration estimates to
#' @param immune_cell_infl_df a data frame from the TIMER website (see link above)
#' that has pre-computed immune cell infiltration estimates for the TCGA
#' @param ici_columns a vector of the columns of CIBERSORT Abs
add_total_immune_cell_percentage <- function(sample_df, immune_cell_infl_df, 
                                             ici_columns) {
  immune_cell_est <- lapply(rownames(sample_df), function(sample) {
    
    # Remove the final letter to match
    sample <- substr(sample, 1, nchar(sample)-1)
    
    # For this sample, get a subset of the IC infiltration data frame with the columns 
    # of interest from the tool of interest
    if (sample %in% immune_cell_infl_df$cell_type) {
      ici_df_sub <- immune_cell_infl_df[grepl(sample, immune_cell_infl_df$cell_type),
                                        colnames(immune_cell_infl_df) %in% ici_columns]
      print(head(ici_df_sub))
      
      # Sum the fraction of all these non-tumor cell types together
      total_ic_frac <- sum(as.numeric(ici_df_sub), na.rm = TRUE)
      
      #return(total_ic_frac)
      
      # Return this value, bucketed 
      if(total_ic_frac <= 0.3) {return(c(1,0,0))}
      else if ((total_ic_frac > 0.3) & (total_ic_frac <= 0.7)) {return(c(0,1,0))}
      else if (total_ic_frac > 0.7) {return(c(0,0,1))}
      else {return(c(NA, NA, NA))}
      
    } else {return (c(NA, NA, NA))}
  })
  immune_cell_est_df <- do.call(rbind, immune_cell_est)
  colnames(immune_cell_est_df) <- c("Tot_IC_Frac_b1", "Tot_IC_Frac_b2", "Tot_IC_Frac_b3")
  
  sample_df <- cbind(sample_df, immune_cell_est_df)

  return(sample_df)
}

############################################################ 
### OPTION 2: GET THE FRACTION OF EACH CELL TYPE IN EACH SAMPLE, WITH A COVARIATE 
### FOR EACH
#' For each sample, get the immune cell fractions for all cell types of interest.
#' Return an input data frame with, for each row (sample), has additional columns
#' with the fraction of each cell type in that sample.
#' @param sample_df a sample data frame to add immune cell infiltration estimates to
#' @param immune_cell_infl_df a data frame from the TIMER website (see link above)
#' that has pre-computed immune cell infiltration estimates for the TCGA
#' @param ici_columns a vector of the columns of the immune cell deconvolution
#' tool of interest
add_immune_cell_fractions <- function(sample_df, immune_cell_infl_df, ici_columns) {
  
  immune_cell_est <- lapply(rownames(sample_df), function(sample) {
    
    # Remove the final letter to match
    sample <- substr(sample, 1, nchar(sample)-1)

    # For this sample, get a subset of the IC infiltration data frame with the columns 
    # of interest from the tool of interest
    if(sample %fin% immune_cell_infl_df$cell_type) {
      ici_df_sub <- immune_cell_infl_df[grepl(sample, immune_cell_infl_df$cell_type),
                                        colnames(immune_cell_infl_df) %in% ici_columns]
    } else { return(NA)}
    
    # Return just this subset (1 row by length(ici_columns) dimensional data frame)
    
    return(ici_df_sub)
  })
  # Bind all these rows together
  ici_dataframe <- do.call(rbind, immune_cell_est)
  print(ici_dataframe)
  #colnames(ici_dataframe) <- ici_columns
  sample_df <- cbind(sample_df, ici_dataframe)
  
  return(sample_df)
}

############################################################ 
#' For each sample, get the first 10 inferred covariate values from PEER.
#' Return an input data frame with, for each row (sample), has additional columns
#' with the PEER factors for that sample.
#' @param sample_df a sample data frame to add PEER estimates to
#' @param peer_df a data frame from the results of PEER (factors vs. samples)
add_peer <- function(sample_df, peer_df) {
  
  peer_rows <- lapply(rownames(sample_df), function(sample) {
    
    # For this sample, get & return the first 10 inferred covariates
    if(sample %fin% colnames(peer_df)) {
      peer_factors <- as.numeric(unlist(peer_df[1:10, colnames(peer_df) == sample]))
      print(peer_factors)
      return(peer_factors)
    } else {return(rep(NA, times=10))}
  })
  
  names(peer_rows) <- rownames(sample_df)
  print(head(peer_rows))

  # Bind all these rows together
  peer_dataframe <- do.call(rbind, peer_rows)[,1:10]
  print(head(peer_dataframe))
  print(dim(peer_dataframe))
  colnames(peer_dataframe) <- c("PEER1", "PEER2", "PEER3", "PEER4", "PEER5", 
                                "PEER6", "PEER7", "PEER8", "PEER9", "PEER10")
  print(peer_dataframe)
  
  sample_df <- cbind(sample_df, peer_dataframe)
  
  return(sample_df)
}

############################################################ 

### OPTION 1: GET THE LIBRARY SIZES DIRECTLY FROM THE RAW COUNT EXPRESSION DATA
#' Add the effective library size for each sample to this sample data frame as 
#' another column  Effective library size (N_j) is the total number of reads for 
#' a given sample after cross-sample count normalization. 
#' @param sample_df a sample data frame to add library sizes to
#' @param expression_df a cross-sample normalized expression data frame with read
#' counts to sum
add_library_sizes <- function(sample_df, expression_df) {
  library_sizes <- lapply(rownames(sample_df), function(sample) {
    # Get the column(s) of the sample of interest in the expression data frame
    if (sample %in% colnames(expression_df)) {
      expression_df_sub <- expression_df[,grepl(sample, colnames(expression_df))]

      if(is.data.frame(expression_df_sub)) {
        expression_df_sub <- expression_df_sub[,1]
      }
      
      # Sum all the reads in this sample column
      sum <- sum(as.numeric(expression_df_sub))
      
      return(sum)
    } else {return(NA)}
  })
  
  # Combine all results into new DF and cbind this to the patient DF
  lib_size_df <- do.call(rbind, library_sizes)
  sample_df <- cbind(sample_df, lib_size_df)
  colnames(sample_df)[ncol(sample_df)] <- "Lib_Size"
  
  return(sample_df)
}


### OPTION 2: GET THE LIBRARY SIZES FOR TMM NORMALIZED SAMPLES RIGHT FROM THE 
### SAMPLE FILES PRODUCED BY EDGER
#' Add the effective library size for each sample to this sample data frame as 
#' another column. Effective library size (N_j) is the total number of reads for
#' a given sample after cross-sample count normalization. 
#' @param sample_df a tissue sample data frame to add library sizes to
#' @param samples_norm_df a samples DF with library size column produced from edgeR
#' during TMM-normalization
add_library_sizes <- function(sample_df, samples_norm_df) {
  library_sizes <- unlist(lapply(rownames(sample_df), function(sample) {
    #print(sample)
    # Get the row of the sample of interest in the expression data frame
    if (grepl(sample, samples_norm_df$description)) {
      lib_size <- samples_norm_df[grepl(sample, samples_norm_df$description), 'lib.size']
      return(lib_size)
    } else {return(NA)}
  }))
  
  # Combine all results into new DF and cbind this to the patient DF
  #lib_size_df <- do.call(rbind, library_sizes)
  #print(library_sizes)
  sample_df <- cbind(sample_df, library_sizes)
  colnames(sample_df)[ncol(sample_df)] <- "Lib_Size"
  
  return(sample_df)
}

############################################################ 

#' Takes a main data frame for all sample characteristics,
#' along with a sub-data frame, and merges the two by their
#' column names
#' @param main_df the patient data frame we are constructing
#' @param partial_df a partial data frame with matching column
#' names we'd like to merge with our main patient data frame
merge_dataframes_by_colname <- function(main_df, partial_df) {
  for (i in 1:ncol(partial_df)) {
    pos <- as.numeric(which(colnames(partial_df)[i] == colnames(main_df)))
    main_df[,pos] <- partial_df[,i]
  }
  return(main_df)
}
