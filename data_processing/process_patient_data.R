############################################################
### PROCESS PATIENT DATA
### PUBLICATION INFORMATION
############################################################

# This file uses a clinical data frame to create an output "patient data frame" 
# that will be used in the Dyscovr model. This DF has the following format:
  # Columns: Linear model input variables (Gender, Age (buckets 1-10, or 
    # normalized), Prior malignancies, Treatment-radiation, 
    # Treatment-pharmacological, Cancer Type or Subtype)
  # Rows: Patients that have all necessary data types 
  # Entries: 0 or 1 values

library(TCGAbiolinks)

# Local PATH to directory containing patient data file
PATH <- getwd()
OUTPUT_PATH <- paste0(getwd(), "Output/Patient/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "Input Data Files/all_genes_id_conv.csv"), 
                              header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
# IMPORT CLINICAL DATA
############################################################
### TCGA ###
clinical_df <- read.csv(paste0(PATH, "clinical_data.csv"), 
                             header = T, check.names = F)

### METABRIC ###
clinical_df <- read.table(paste0(PATH, "data_clinical_patient_sub.txt"), 
                          header = T, check.names = F)
# Add TMB data for METABRIC only
sample_df <- read.table(paste0(PATH, "data_clinical_sample_sub.txt"), 
                        header = T, check.names = F, sep = ",")
clinical_df <- merge(clinical_df, sample_df[,c('PATIENT_ID', 
                                               "TMB_NONSYNONYMOUS")], 
                     by = 'PATIENT_ID')

############################################################
# IMPORT SUBTYPE INFORMATION (TCGA)
############################################################
# Useful resource for compilation of TCGA molecular subtypes: 
# https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/subtypes.html
cancer_types <- c("ACC", "BRCA", "BLCA", "CESC", "CHOL", "COAD", "ESCA",
                          "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC",
                          "LUAD", "LUSC", "PAAD", "PCPG", "PRAD", "READ", "SKCM",
                          "SARC", "STAD", "THCA", "UCEC", "UCS", "UVM")
cancer_types_no_subtype_info <- c("DLBC", "MESO", "THYM", "OV", "TGCT", "LAML")

import_subtype_information <- function(cancer_types) {
  subtype_files <- lapply(cancer_types, function(type) {
    subtype <- TCGAquery_subtype(tumor = type)
    return(subtype)
  })
  names(subtype_files) <- general_cancer_types
  
  # Note: To combine into one file, we'll have to do some work to ID matching 
  # columns with differing names
  return(subtype_files)
}

# Pan-cancer subtype info
pc_subtype <- import_subtype_information(cancer_types)

# Or, if just using one cancer type, e.g. BRCA:
brca_subtype <- TCGAquery_subtype(tumor = "BRCA")


############################################################
# CREATE SAMPLE DATA FRAME FOR DYSCOVR
############################################################
### Sample Call (TCGA) ###
# For a single cancer type, e.g. BRCA
brca_patient_df <- create_patient_dataframe(clinical_df, is_pc, brca_subtype, 
                                            NA, F, 'tcga')

# Pan-cancer
pc_patient_df <- create_patient_dataframe(clinical_df, NA, is_pc, pc_subtype, 
                                          NA, F, 'tcga')

# Write to CSV
write.csv(brca_patient_df, paste0(
  OUTPUT_PATH, "Dyscovr_Input/Patient/patient_dataframe.csv"))
write.csv(pc_patient_df, paste0(
  OUTPUT_PATH, "Dyscovr_Input/Patient/patient_dataframe.csv"))

#
#
#

# Run individually for each cancer type, across all cancer types, and get back a 
# list of patient DFs with subtype information
pan_cancer_patient_dfs <- lapply(1:length(pc_subtype), function(i) {
  st_df <- list(pc_subtype[[i]])
  cancer_name <- names(pc_subtype)[i]
  names(st_df) <- cancer_name
  
  clin_df_sub <- clinical_df[grepl(cancer_name, clinical_df$project_id),]
  
  if(!((length(clin_df_sub) == 0) | (nrow(clin_df_sub) == 0))) {
    pat_df <- create_patient_dataframe(clin_df_sub, NA, is_pc, st_df, 
                                       NA, F, 'tcga')
    return(pat_df)
  }
  return(NA)
})
names(pan_cancer_patient_dfs) <- names(pc_subtype)


# Run for cancer types without subtype information
pan_cancer_patient_dfs_pt2 <- lapply(
  cancer_types_no_subtype_info, function(cancer_name){
    clin_df_sub <- clinical_df[grepl(cancer_name, clinical_df$project_id),]
    if(!((length(clin_df_sub) == 0) | (nrow(clin_df_sub) == 0))) {
      pat_df <- create_patient_dataframe(clin_df_sub, NA, is_pc, NA, 
                                         NA, F, 'tcga')
      return(pat_df)
    }
})
names(pan_cancer_patient_dfs_pt2) <- cancer_types_no_subtype_info

# Join together for patients with and without subtype information
pan_cancer_patient_dfs_full <- c(pan_cancer_patient_dfs, pan_cancer_patient_dfs_pt2)

# Write to files
lapply(1:length(pan_cancer_patient_dfs_full), function(i) {
  fn <- paste0("patient_dataframe_", 
               paste0(names(pan_cancer_patient_dfs_full)[i], 
                      "_inclSubtypes.csv"))
  write.csv(pan_cancer_patient_dfs_full[[i]], 
            paste0(OUTPUT_PATH, paste0("Dyscovr_Input/Patient/", fn)))
})

# Add patients without subtype info to the pan-cancer subtype DF
last_ct <- colnames(pc_patient_df)[ncol(pc_patient_df)]  #143
curr_b <- 143+1
new_rows <- list(pc_patient_df)
new_rows_index <- 2

for(i in 1:length(pan_cancer_patient_dfs_pt2)) {
  df <- pan_cancer_patient_dfs_pt2[[i]]
  
  # Add in the missing cancer type info
  df_blank <- data.frame(matrix(nrow = nrow(df), 
                                ncol = ncol(pc_patient_df[,grepl(
                                  "Cancer_type", colnames(pc_patient_df))]) + 
                                  length(pan_cancer_patient_dfs_pt2)))
  df_blank[is.na(df_blank)] <- 0
  df_full <- cbind(df, df_blank)
  colnames(df_full) <- c(
    colnames(df), paste0("Cancer_type_b", 
                         1:(ncol(pc_patient_df[,grepl("Cancer_type", 
                                                      colnames(pc_patient_df))]) + 
                              length(pan_cancer_patient_dfs_pt2))))
  
  # Add in a new cancer type variable for this cancer type
  df_full[, curr_b] <- rep(1, times = nrow(df_full))
  
  new_rows[[new_rows_index]] <- as.data.frame(df_full)
  curr_b <- curr_b + 1
  new_rows_index <- new_rows_index + 1
}

df_blank <- data.frame(matrix(nrow = nrow(pc_patient_df), 
                              ncol = length(pan_cancer_patient_dfs_pt2)))
df_blank[is.na(df_blank)] <- 0

# Reset curr_b
curr_b <- 143+1
colnames(df_blank) <- paste0("Cancer_type_b", 
                             curr_b:(curr_b+(length(pan_cancer_patient_dfs_pt2)-1)))
pc_patient_df <- cbind(pc_patient_df, df_blank)
pc_patient_df_new <- do.call(rbind, new_rows)

# Write to CSV
write.csv(pc_patient_df_new, pc_patient_df, paste0(
  OUTPUT_PATH, "Dyscovr_Input/Patient/patient_dataframe_inclSubtypes.csv"))

#
#
#

### METABRIC ###
brca_patient_df <- create_patient_dataframe(clinical_df, NA, is_pc, NA, NA, 
                                            F, 'metabric')
# Write to CSV
write.csv(brca_patient_df, paste0(
  OUTPUT_PATH, "Dyscovr_Input/Patient/patient_dataframe.csv"))


############################################################

#' Creates a skeleton "patient data frame" that will contain all the 
#' characteristics of interest for each patient. Calls helper function to fill 
#' patient data frame and then returns it.
#' @param clinical_df clinical data file from the TCGA for the given patient 
#' cohort
#' @param is_pc a T/F value indicating whether the cohort is per-cancer or not
#' @param subtype_file a file from TCGAbiolinks that has information about 
#' subtypes or, alternatively (if pan-cancer) a file that identifies the cancer 
#' type of each sample
#' @param additional_subtype_info if is_pc is F, a file from TCGAbiolinks that 
#' has information about subtypes or, alternatively (if pan-cancer) a file that 
#' identifies the cancer type of each sample. If is_pc is T, is a list of 
#' separate individual cancer subtype types (the names of list items are the 
#' cancer type that subtype file is for)
#' @param age_bucketing a T/F value to indicate whether or not we are bucketing
#' age (if not, we are normalizing it to be between 0 and 1)
#' @param dataset either 'tcga' or 'metabric', to refer to proper column names
create_patient_dataframe <- function(clinical_df, is_pc, additional_subtype_info, 
                                     age_bucketing, dataset) {
  
  # Create relational data frame to hold patient-dependent inputs to linear model
  # Columns: Dyscovr input variables ; Rows: Samples that have all necessary 
  # data types (Unique 4-digit ID)
  if(age_bucketing) {
    patient_characteristics <- c("Gender", "Age_b1", "Age_b2", "Age_b3", 
                                 "Age_b4", "Age_b5", "Age_b6", "Age_b7", 
                                 "Age_b8", "Age_b9", "Age_b10", "Prior_malig", 
                                 "Treatment_rad", "Treatment_pharm", "Tot_Mut_b1", 
                                 "Tot_Mut_b2", "Tot_Mut_b3", "PC1", "PC2", "PC3")
  } else {
    patient_characteristics <- c("Gender", "Age", "Prior_malig", "Treatment_rad", 
                                 "Treatment_pharm", "Tot_Mut_b1", "Tot_Mut_b2", 
                                 "Tot_Mut_b3", "PC1", "PC2", "PC3")
  }
  if(is_pc) {
    patient_characteristics <- add_cancer_type(patient_characteristics, is_pc, 
                                               additional_subtype_info, 
                                               clinical_df, dataset)
  }
  print(patient_characteristics)
  
  unique_patient_ids <- c()
  if(dataset == 'tcga') {
    unique_patient_ids <- unique(unlist(lapply(
      clinical_df$case_submitter_id, function(x) 
        unlist(strsplit(x, split = "-", fixed = T))[3])))
  } else if (dataset == "metabric") {
    unique_patient_ids <- unique(clinical_df$PATIENT_ID)
  } else {print("Only implemented for TCGA and METABRIC datasets.")}
  
  # Create skeleton patient data frame of the correct dimensions
  patient_dataframe <- data.frame(matrix(ncol = length(patient_characteristics), 
                                         nrow = length(unique_patient_ids)))
  colnames(patient_dataframe) <- patient_characteristics
  rownames(patient_dataframe) <- unique_patient_ids
  
  # Input all the patient-specific information 
  patient_dataframe <- input_patient_specific_info(patient_dataframe, 
                                                   clinical_df, is_pc,
                                                   additional_subtype_info, 
                                                   age_bucketing, dataset)
  
  return(patient_dataframe)
}


############################################################

#' Takes the empty patient data frame, along with the clinical data frame and 
#' updates patient DF with all the patient-specific information: gender,
#' age, prior malignancies, treatment, and total number of mutations. 
#' @param input_df a blank data frame with the dimensions of the the number
#' of patient characteristics (columns) and # of patients in the cohort (rows)
#' @param clinical_df clinical data file from the TCGA for the given patient 
#' cohort
#' @param is_pc a T/F value indicating whether the cohort is per-cancer or not
#' @param additional_subtype_info if is_pc is F, a file from TCGAbiolinks that 
#' has information about subtypes or, alternatively (if pan-cancer) a file that 
#' identifies the cancer type of each sample. If is_pc is T, is a list of 
#' separate individual cancer subtype types (the names of list items are the 
#' cancer type that subtype file is for)
#' @param age_bucketing a T/F value to indicate whether or not we are bucketing
#' age (if not, we are normalizing it to be between 0 and 1)
#' @param dataset either 'tcga' or 'metabric' to refer to proper column names
input_patient_specific_info <- function(input_df, clinical_df, is_pc, 
                                        additional_subtype_info, age_bucketing, 
                                        dataset) {
  input_dataframe_updated <- input_df
  
  # GENDER
  gender_vect <- c()
  if(dataset == 'tcga') {
    gender_vect <- clinical_df$gender[seq(1, length(clinical_df$gender), 2)] 
  } else if (dataset == "metabric") {
    gender_vect <- clinical_df$SEX
  } else {print("Only implemented for TCGA and METABRIC.")}
  gender_vect_binary <- ifelse((gender_vect == "male") | 
                                 (gender_vect == "Male"), 0, 1)
  input_dataframe_updated[,'Gender'] <- gender_vect_binary
  
  # AGE
  if(age_bucketing) {
    age_bin_df <- data.frame()
    if(dataset == 'tcga') {
      age_bin_df <- convert_age_to_bins(clinical_df$age_at_index, dataset)
    } else if (dataset == 'metabric') {
      age_bin_df <- convert_age_to_bins(clinical_df$AGE_AT_DIAGNOSIS, dataset)
    } else {print("Only implemented for TCGA and METABRIC.")}
    input_dataframe_updated <- merge_dataframes_by_colname(
      input_dataframe_updated, age_bin_df) 
  } else {
    age_vect <- c()
    if(dataset == 'tcga') {
      age_vect <- as.numeric(clinical_df$age_at_index)[
        seq(1, length(clinical_df$age_at_index), 2)] / 100
    } else if (dataset == 'metabric') {
      age_vect <- clinical_df$AGE_AT_DIAGNOSIS / 100
    } else {print("Only implemented for TCGA and METABRIC.")}
    input_dataframe_updated[,'Age'] <- age_vect
  }

  # PRIOR_MALIG
  if(dataset == 'tcga') {
    prior_malig_vect <- clinical_df$prior_malignancy[
      seq(1, length(clinical_df$prior_malignancy), 2)]
    input_dataframe_updated[,'Prior_malig'] <- ifelse(
      prior_malig_vect == "no" | prior_malig_vect == "not reported", 0, 1)
  } else if (dataset == 'metabric') {
    input_dataframe_updated[,'Prior_malig'] <- unlist(lapply(
      clinical_df$RFS_STATUS, function(x) unlist(strsplit(x, ":", fixed = T))[1]))
  } else {print("Only implemented for TCGA and METABRIC.")}
  
  # TREATMENT
  treatment_bin_df <- data.frame()
  if(dataset == 'tcga') {
    treatment_or_ther_vect <- clinical_df$treatment_or_therapy
    treatment_type_vect <- clinical_df$treatment_type
    treatment_bin_df <- convert_treatments_to_bins(treatment_or_ther_vect, 
                                                   treatment_type_vect)
  } else if (dataset == 'metabric') {
    treatment_bin_df <- convert_treatments_to_bins_metabric(clinical_df)
  } else {print("Only implemented for TCGA and METABRIC.")}
  input_dataframe_updated <- merge_dataframes_by_colname(input_dataframe_updated, 
                                                         treatment_bin_df)
  
  # TOTAL NUM MUTATIONS
  if(dataset == 'tcga') {
    tot_num_mut_df <- convert_tot_num_mut_to_bins(clinical_df$Total.Num.Muts, 
                                                  dataset)
  } else if (dataset == 'metabric') {
    tot_num_mut_df <- convert_tot_num_mut_to_bins(
      log2(clinical_df$TMB_NONSYNONYMOUS + 1), dataset)
  } else {print("Only implemented for TCGA or METABRIC.")}
  input_dataframe_updated <- merge_dataframes_by_colname(input_dataframe_updated, 
                                                         tot_num_mut_df)
  
  # CANCER TYPE OR SUBTYPE
  if(length(additional_subtype_info) >= 1) {
    cancer_type_bin_df <- convert_ct_to_bins(clinical_df, 
                                             input_dataframe_updated, 
                                             is_pc, additional_subtype_info, 
                                             dataset)
    input_dataframe_updated <- merge_dataframes_by_colname( 
      input_dataframe_updated, cancer_type_bin_df)
  }
  
  return(input_dataframe_updated)
}


############################################################
#' Takes a main data frame for all patient characteristics,
#' along with a sub-data frame, and merges the two by their
#' column names
#' @param main_df the patient data frame we are constructing
#' @param partial_df a partial data frame with matching column
#' names we'd like to merge with our main patient data frame
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


############################################################ 
#' Takes a vector of patient ages and creates a DF with separate columns 
#' (one for each of the bins) consisting of 0 and 1 values for each patient 
#' (0 if they are not in that age bin, 1 if they are)
#' @param age_vector a vector of all patient ages in the cohort, two entries 
#' per patient
#' @param dataset either 'tcga' or 'metabric'; to refer to proper column names
convert_age_to_bins <- function(age_vector, dataset) {
  
  # For TCGA, get only every second item in the vector (patients all have two 
  # almost identical entries; all age entries are exactly duplicated)
  if(dataset == "tcga") {age_vector <- age_vector[seq(1, length(age_vector), 2)]} 
  
  # Create data frame (columns are bins, rows are patients)
  num_bins <- 10
  age_df <- data.frame(matrix(nrow = length(age_vector), ncol = num_bins))
  colnames(age_df) <- c("Age_b1", "Age_b2", "Age_b3", "Age_b4", "Age_b5", 
                        "Age_b6", "Age_b7", "Age_b8", "Age_b9", "Age_b10")
  
  for (i in 1:length(age_vector)) {
    age <- age_vector[i]
    age_df[i,] <- rep(0, num_bins)        # set all cols for that row equal to 0
    if (age <= 9) {
      age_df[i,"Age_b1"] <- 1      # set this bucket for that row equal to 1
    } else if (age > 9 & age <= 19) {
      age_df[i,"Age_b2"] <- 1      # set this bucket for that row equal to 1
    } else if (age > 19 & age <= 29) {
      age_df[i,"Age_b3"] <- 1      # set this bucket for that row equal to 1
    } else if (age > 29 & age <= 39) {
      age_df[i,"Age_b4"] <- 1      # set this bucket for that row equal to 1
    } else if (age > 39 & age <= 49) {
      age_df[i,"Age_b5"] <- 1      # set this bucket for that row equal to 1
    } else if (age > 49 & age <= 59) {
      age_df[i,"Age_b6"] <- 1      # set this bucket for that row equal to 1
    } else if (age > 59 & age <= 69) {
      age_df[i,"Age_b7"] <- 1      # set this bucket for that row equal to 1
    } else if (age > 69 & age <= 79) {
      age_df[i,"Age_b8"] <- 1      # set this bucket for that row equal to 1
    } else if (age > 79 & age <= 89) {
      age_df[i,"Age_b9"] <- 1      # set this bucket for that row equal to 1
    } else if (age > 89) {
      age_df[i,"Age_b10"] <- 1      # set this bucket for that row equal to 1
    } else {
      print("This is an invalid age. Please address this issue & re-try.")
    }
  }
  return(age_df)
}


############################################################
#' Takes vectors of whether the patient has been treated (yes/ no) 
#' for a given therapy (radiation/ pharmaceutical). Combines these 
#' vectors into one DF (rows are patients, columns are bins : treatment_rad 
#' and treatment_pharm) which have a 0 if patient has not been treated with 
#' that particular kind of therapy and 1 if they have
#' @param treatment_or_ther_vect a vector of whether or not a patient has been 
#' treated with a given therapy ('yes' or 'no' strings)
#' @param treatment_type_vect a vector of the treatment type corresponding
#' to that 'yes' or 'no' (either 'Pharmaceutical Therapy, NOS' or 'Radiation
#' Therapy, NOS) 
convert_treatments_to_bins <- function(treatment_or_ther_vect, 
                                       treatment_type_vect) {
  
  # Create data frame (columns are bins, rows are patients)
  num_bins <- 2
  num_pat <- length(treatment_or_ther_vect) / 2
  treatment_df <- data.frame(matrix(nrow = num_pat, ncol = num_bins))
  colnames(treatment_df) <- c("Treatment_rad", "Treatment_pharm")
  
  treatment_vect_ind <- 1         # index into our treatment vector
  for (i in 1:nrow(treatment_df)) {
    # Extract the two sets of treatment and rad/pharm combos for patient i
    treatment_1 <- treatment_or_ther_vect[treatment_vect_ind]
    rad_or_pharm_1 <- treatment_type_vect[treatment_vect_ind]
    treatment_2 <- treatment_or_ther_vect[treatment_vect_ind+1]
    rad_or_pharm_2 <- treatment_type_vect[treatment_vect_ind+1]
    
    # Populate patient i's row with 0's to start
    treatment_df[i,] <- rep(0, num_bins) 
    
    # Use helper function to update the data frame with the treatment & rad/pharm
    # information for patient i
    treatment_df <- rad_or_pharm_helper(rad_or_pharm_1, treatment_1, 
                                        treatment_df, i)
    treatment_df <- rad_or_pharm_helper(rad_or_pharm_2, treatment_2, 
                                        treatment_df, i)
    
    # Update indices
    treatment_vect_ind <- treatment_vect_ind + 1
  }
  return(treatment_df)
}

############################################################ 
#' Helper function for converting treatment to bins. Takes a treatment category 
#' variable (radiation or pharma) and a yes or no treatment variable and 
#' updates/returns given treatment data frame
#' @param rad_or_pharm a string that is either 'Radiation Therapy,
#' NOS' or 'Pharmaceutical Therapy, NOS')
rad_or_pharm_helper <- function(rad_or_pharm, treatment, treatment_df, i) {
  if (rad_or_pharm == "Radiation Therapy, NOS") {
    if (treatment == "yes") {
      treatment_df[i, 'Treatment_rad'] <- 1
    }
  } else if (rad_or_pharm == "Pharmaceutical Therapy, NOS") {
    if (treatment == "yes") {
      treatment_df[i, 'Treatment_pharm'] <- 1
    }
  } else {
    print("Unseen treatment type. Look into this and retry.")
  }
  return(treatment_df)
}


############################################################ 
#' METABRIC-specific function for converting treatment to bins. Takes a
#' clinical data frame and converts to a binned treatment DF. Groups 
#' hormone therapy and chemotherapy together as "pharmacological" therapy.
#' @param clinical_df a patient clinical DF from METABRIC
convert_treatments_to_bins_metabric <- function(clinical_df) {
  treatment_df <- data.frame(matrix(nrow = nrow(clinical_df), ncol = 2))
  rownames(treatment_df) <- clinical_df$PATIENT_ID
  colnames(treatment_df) <- c("Treatment_rad", "Treatment_pharm")
  
  for(i in 1:nrow(clinical_df)) {
    # Get radiation status
    treatment_df[i, 'Treatment_rad'] <- ifelse(
      clinical_df[i, 'RADIO_THERAPY'] == "YES", 1, 0)
    # Get pharmacological status
    treatment_df[i, 'Treatment_pharm'] <- ifelse(
      (clinical_df[i, 'CHEMOTHERAPY'] == "YES") | 
        (clinical_df[i, 'HORMONE_THERAPY'] == "YES"), 1, 0)
  }
  #print(head(treatment_df))
  
  return(treatment_df)
}

############################################################
#' Takes the vector of patient characteristics and adds buckets for the 
#' appropriate number of subtypes (for a specific cancer type) or types 
#' (for pan-cancer analysis). Returns the updated vector.
#' @param patient_characterisitics a vector of existing patient characteristics
#' to be updated with cancer type/subtype
#' @param is_pc a T/F value indicating whether the cohort is per-cancer or not
#' @param additional_subtype_info if is_pc is F, a file from TCGAbiolinks that 
#' has information about subtypes or, alternatively (if pan-cancer) a file that 
#' identifies the cancer type of each sample. If is_pc is T, is a list of 
#' separate individual cancer subtype types (the names of list items are the 
#' cancer type that subtype file is for)
#' @param clinical_df clinical data file from the TCGA for the given patient 
#' cohort
#' @param dataset either 'tcga' or 'metabric'; to refer to proper column names
add_cancer_type <- function(patient_characteristics, is_pc, 
                            additional_subtype_info, clinical_df, dataset) {
  
  # Looking pan-cancer broadly, or at one individual cancer type
  if (!is_pc) {
    
    # If pan-cancer, find how many unique cancer types there are in total
    unique_cts <- unique(clinical_df$project_id)
    unique_cts <- unlist(lapply(unique_cts, function(ct) 
      unlist(strsplit(ct, "-", fixed = T))[2]))
    
    # Add subtype information for either pan-cancer or each cancer type in the 
    # pan-cancer file
    if(length(additional_subtype_info) >= 1) {
      # Remove the broader "cancer type" label for the cancer types that have 
      # subtype information we want to include
      unique_cts <- unique_cts[!grepl(paste(names(additional_subtype_info), 
                                            collapse = "|"), unique_cts)]
      unique_subtypes <- unique(unlist(lapply(
        1:length(additional_subtype_info), function(i) {
          subtype_df <- additional_subtype_info[[i]]
          cancer_name <- names(additional_subtype_info)[i]
          subtypes <- get_subtype(NA, subtype_df, cancer_name)
          subtypes_full <- paste(trimws(cancer_name), trimws(subtypes), 
                                 sep = ".")
        return(subtypes_full)
      })))
      unique_cts <- c(unique_cts, unique_subtypes)
    }
    length_cts <- length(unique_cts)
  }
  
  # Not pan-cancer: look at cancer subtypes for each cancer type individually
  # within the set of cancer types under consideration
  else {
    # Find the subtype column and check how many types are given
    unique_subtypes <- NA
    
    if (dataset == 'tcga') {
      if(T %fin% grepl("subtype", colnames(additional_subtype_info), 
                       ignore.case = T)) {
        vals <- grepl("subtype", colnames(additional_subtype_info), 
                      ignore.case = T)
        if(length(vals[vals == T]) == 1) {
          unique_subtypes <- unique(as.character(unlist(
            additional_subtype_info[ ,grepl("subtype", 
                                            colnames(additional_subtype_info), 
                                            ignore.case = T)])))
        } 
      } else {
        unique_subtypes <- unique(as.character(unlist(
          additional_subtype_info[ ,grepl("mRNA cluster", 
                                          colnames(additional_subtype_info), 
                                          ignore.case = T)])))
      }
    } else if (dataset == 'metabric') {
      unique_subtypes <- unique(unlist(clinical_df[, 'CLAUDIN_SUBTYPE']))
      unique_subtypes <- unique_subtypes[unique_subtypes != "claudin-low"]
      
    } else {print("Only implemented for TCGA and METABRIC.")}
    
    unique_subtypes <- unique_subtypes[!(is.na(unique_subtypes) | 
                                           unique_subtypes == "NA")]
    length_cts <- length(unique_subtypes)
    
  }
  
  # Create buckets equivalent to the number of types or subtypes
  new_buckets <- unlist(lapply(1:length_cts, function(x) 
    paste("Cancer_type_b", x, sep = "")))
  
  # Add new buckets to the patient characteristics list and return
  patient_characteristics <- c(patient_characteristics, new_buckets)
  
  return(patient_characteristics)
}


############################################################ 
#' Takes a file with information about cancer type or subtype per patient and 
#' creates a DF with separate columns (one per bin) consisting of 0 and 1 values 
#' for each patient (0 if they are not in that bin, 1 if they are)
#' @param clinical_df clinical data file from the TCGA for the given patient 
#' cohort
#' @param input_dataframe_updated the work-in-progress input data frame with all
#' the patient characteristics for the linear model
#' @param is_pc a T/F value indicating whether the cohort is per-cancer or not
#' @param additional_subtype_info if is_pc is F, a file from TCGAbiolinks that 
#' has information about subtypes or, alternatively (if pan-cancer) a file that 
#' identifies the cancer type of each sample. If is_pc is T, is a list of 
#' separate individual cancer subtype types (the names of list items are the 
#' cancer type that subtype file is for)
#' @param dataset either 'tcga' or 'metabric'; to refer to proper column names
convert_ct_to_bins <- function(clinical_df, input_dataframe_updated, 
                               is_pc, additional_subtype_info, dataset) { 
  
  # Create a blank data frame with all the same patient IDs as the input data frame
  cancer_type_bin_df <- input_dataframe_updated[, grepl(
    "Cancer_type", colnames(input_dataframe_updated))]
  
  # Keep a growing map of which type will correspond to which bucket
  ct_to_bucket_map <- data.frame(ct = rep(NA, ncol(cancer_type_bin_df)), 
                                 bucket_num = 1:ncol(cancer_type_bin_df))
  
  patients_w_subtype_info <- c()
  cancer_types_no_subtype_info <- NA
  

  
  if(is_pc) {
    patients_w_subtype_info <- unlist(lapply(additional_subtype_info, function(x) {
      pats <- unlist(lapply(as.character(x$patient), function(p) 
        unlist(strsplit(p, "-", fixed = T))[3]))
      return(pats)                       
    }))
    cancer_types_no_subtype_info <- setdiff(unique(unlist(lapply(
      clinical_df$project_id, function(x) 
        unlist(strsplit(x, "-", fixed = T))[2]))), 
      names(additional_subtype_info))
  } else {
    patients_w_subtype_info <- as.character(additional_subtype_info$patient)
  }
  
  # Go through each patient and find their subtype/type info, then add to the 
  # appropriate bucket
  for (i in 1:nrow(cancer_type_bin_df)) {
    patient_id <- rownames(cancer_type_bin_df)[i]
    
    ct <- as.character(unlist(clinical_df[grepl(patient_id, 
                                                clinical_df$case_submitter_id), 
                                          'project_id']))
    ct <- unlist(strsplit(ct, "-", fixed = T))[2]
    
    # Per-cancer: get the cancer type first, then the subtype
    if (is_pc) {
      if(length(patients_w_subtype_info) != 0) {
        pat_found <- 0
        for(v in 1:length(additional_subtype_info)) {
          subtype_df <- additional_subtype_info[[v]]
          cancer_name <- names(additional_subtype_info)[v]
          if((T %fin% grepl(patient_id, subtype_df$patient)) | 
             (T %fin% grepl(patient_id, subtype_df$sample))) {
            st <- get_subtype(patient_id, subtype_df, cancer_name)[1]
            # Create a compound label with cancer type.subtype
            ct <- paste(trimws(cancer_name), trimws(st), sep = ".")
            pat_found <- 1
          }
        }
        if(pat_found == 0) {
          print("Patient not found in additional subtype info")
          if(!(ct %fin% cancer_types_no_subtype_info)) {
            ct <- paste(ct, "NA", sep = ".")
          }
        } 
      } 
      if(length(ct) == 0) {ct <- "Not Evaluable"}
      print(paste("CT:", ct))
    } 
    
    # Either pan-cancer or in an individual cancer type: get the cancer subtype
    else {
      if(dataset == "tcga") {
        # Automatically detect which column is the subtype variable based on 
        # common subtype column names; can alternatively adjust this code 
        # directly to search for the name of the exact subtype column desired
        if(T %fin% grepl("subtype", colnames(subtype_file), ignore.case = T)) {
          vals <- grepl("subtype", colnames(subtype_file), ignore.case = T)
          if(length(vals[vals == T]) == 1) {
            ct <- as.character(unlist(subtype_file[grepl(patient_id, 
                                                         subtype_file$patient), 
                                                   grepl("subtype", 
                                                         colnames(subtype_file), 
                                                         ignore.case = T)]))
          } 
        } else {
          ct <- as.character(unlist(subtype_file[grepl(patient_id, 
                                                       subtype_file$patient), 
                                                 grepl("mRNA cluster", 
                                                       colnames(subtype_file), 
                                                       ignore.case = T)]))
        }
        
      } else if (dataset == "metabric") {
        # Use the CLAUDIN_SUBTYPE variable for subtype, and merge basal and 
        # claudin-low
        ct <- as.character(clinical_df[grepl(patient_id, clinical_df$PATIENT_ID), 
                                       'CLAUDIN_SUBTYPE'])
        if (ct == 'claudin-low') {ct <- "Basal"}
      } else {print("Only implemented for TCGA and METABRIC.")}
      
    }

    # If this cancer type isn't already in the cancer type-to-bucket mapping, 
    # add it
    if(!(ct %fin% ct_to_bucket_map$ct)) {
      # Add a new entry to the map (replace the next NA)
      pos_of_first_NA <- min(which(is.na(ct_to_bucket_map$ct)))
      ct_to_bucket_map$ct[pos_of_first_NA] <- ct
    }
    
    # Get the bucket that it's in
    bucket_num <- ct_to_bucket_map[ct_to_bucket_map$ct == ct, 'bucket_num']
    bucket_num <- bucket_num[!is.na(bucket_num)]
    
    # Create a new row where there is a '1' for this bucket and a '0' for 
    # everything else
    new_row <- rep(0, times = ncol(cancer_type_bin_df))
    new_row[bucket_num] <- 1
    
    # Update the blank DF for this patient
    cancer_type_bin_df[i,] <- new_row
  }
  
  # Return the newly binned DF (with patient IDs that match the order of the 
  # growing input DF)
  return(cancer_type_bin_df)
}


#' Helper function that gets the subtype for a given patient of interest
#' @param patient_id the patient ID; if patient ID is NA, return all the unique 
#' subtypes for the given CT
#' @param subtype_df the subtype DF for the particular subtype in which the 
#' patient is found
#' @param name the name of the cancer type
get_subtype <- function(patient_id, subtype_df, name) {
  if(!is.na(patient_id)) {
    st <- switch(
      name,
      "COAD" = as.character(unlist(subtype_df[grepl(patient_id, 
                                                    subtype_df$patient), 
                                              grepl("MSI", colnames(subtype_df), 
                                                    ignore.case = T)])),
      "READ" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("MSI", colnames(subtype_df), 
                                                    ignore.case = T)])),
      "ESCA" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("MSI", colnames(subtype_df), 
                                                    ignore.case = T)])),
      "ACC" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                             grepl("COC", colnames(subtype_df), 
                                                   ignore.case = T)])),
      "BRCA" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("PAM50", colnames(subtype_df),
                                                    ignore.case = T)])),
      "BLCA" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("mRNA cluster", 
                                                    colnames(subtype_df), 
                                                    ignore.case = T)])),
      "CESC" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient),
                                              grepl("SAMP:CIMP_call", 
                                                    colnames(subtype_df), 
                                                    ignore.case = T)])),
      "CHOL" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("mRNA-seq", colnames(subtype_df), 
                                                    ignore.case = T)])),
      "PRAD" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("SUBTYPE", colnames(subtype_df), 
                                                    ignore.case = T)])),
      "SKCM" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("RNASEQ.CLUSTER", 
                                                    colnames(subtype_df), 
                                                    ignore.case = T)])),
      "GBM" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                             grepl("Original.Subtype", 
                                                   colnames(subtype_df), 
                                                   ignore.case = T)])),
      "LGG" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                             grepl("Original.Subtype", 
                                                   colnames(subtype_df), 
                                                   ignore.case = T)])),
      "HNSC" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("RNA", colnames(subtype_df), 
                                                    ignore.case = T)])),
      "KICH" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("Histological.Subtype", 
                                                    colnames(subtype_df), 
                                                    ignore.case = T)])),
      "KIRC" = "", # no real subtypes, is already a subtype of RCC
      "KIRP" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("tumor_type.KIRP.path.", 
                                                    colnames(subtype_df), 
                                                    ignore.case = T)])),
      "KICH" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("Histological.Subtype", 
                                                    colnames(subtype_df), 
                                                    ignore.case = T)])),
      "LIHC" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient),
                                              grepl("iCluster", 
                                                    colnames(subtype_df), 
                                                    ignore.case = T)])),
      "LUAD" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("iCluster", 
                                                    colnames(subtype_df), 
                                                    ignore.case = T)])),
      "LUSC" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("Expression.Subtype", 
                                                    colnames(subtype_df), 
                                                    ignore.case = T)])),
      "PAAD" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient),
                                              grepl("Histological type by RHH", 
                                                    colnames(subtype_df),
                                                    ignore.case = T)])),
      "PCPG" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("mRNA Subtype Clusters", 
                                                    colnames(subtype_df), 
                                                    ignore.case = T)])),
      "SARC" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("short histo", 
                                                    colnames(subtype_df), 
                                                    ignore.case = T)])),
      "STAD" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("Molecular.Subtype", 
                                                    colnames(subtype_df),
                                                    ignore.case = T)])),
      "THCA" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$sample),
                                              (colnames(subtype_df) == 
                                                 "mRNA_Cluster_number")])),
      "UCEC" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                              grepl("msi", colnames(subtype_df), 
                                                    ignore.case = T)])),
      "UCS" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                             colnames(subtype_df) == 
                                               "histologic subtype"])),
      "UVM" = as.character(unlist(subtype_df[grepl(patient_id, subtype_df$patient), 
                                             grepl("mRNA Cluster", 
                                                   colnames(subtype_df), 
                                                   ignore.case = T)])))
  } else {
    st <- switch(name,
                 "COAD" = unique(as.character(unlist(subtype_df[, grepl(
                   "MSI", colnames(subtype_df), ignore.case = T)]))),
                 "READ" = unique(as.character(unlist(subtype_df[, grepl(
                   "MSI", colnames(subtype_df), ignore.case = T)]))),
                 "ESCA" = unique(as.character(unlist(subtype_df[, grepl(
                   "MSI", colnames(subtype_df), ignore.case = T)]))),
                 "ACC" = unique(as.character(unlist(subtype_df[, grepl(
                   "COC", colnames(subtype_df), ignore.case = T)]))),
                 "BRCA" = unique(as.character(unlist(subtype_df[, grepl(
                   "PAM50", colnames(subtype_df), ignore.case = T)]))),
                 "BLCA" = unique(as.character(unlist(subtype_df[, grepl(
                   "mRNA cluster", colnames(subtype_df), ignore.case = T)]))),
                 "CESC" = unique(as.character(unlist(subtype_df[, grepl(
                   "SAMP:CIMP_call", colnames(subtype_df), ignore.case = T)]))),
                 "CHOL" = unique(as.character(unlist(subtype_df[, grepl(
                   "mRNA-seq", colnames(subtype_df), ignore.case = T)]))),
                 "PRAD" = unique(as.character(unlist(subtype_df[, grepl(
                   "SUBTYPE", colnames(subtype_df), ignore.case = T)]))),
                 "SKCM" = unique(as.character(unlist(subtype_df[, grepl(
                   "RNASEQ.CLUSTER", colnames(subtype_df), ignore.case = T)]))),
                 "GBM" = unique(as.character(unlist(subtype_df[, grepl(
                   "Original.Subtype", colnames(subtype_df), ignore.case = T)]))),
                 "LGG" = unique(as.character(unlist(subtype_df[, grepl(
                   "Original.Subtype", colnames(subtype_df), ignore.case = T)]))),
                 "HNSC" = unique(as.character(unlist(subtype_df[, grepl(
                   "RNA", colnames(subtype_df), ignore.case = T)]))),
                 "KICH" = unique(as.character(unlist(subtype_df[, grepl(
                   "Histological.Subtype", colnames(subtype_df), ignore.case = T)]))),
                 "KIRC" = "", # no real subtypes, is already a subtype of RCC
                 "KIRP" = unique(as.character(unlist(subtype_df[, grepl(
                   "tumor_type.KIRP.path.", colnames(subtype_df), ignore.case = T)]))),
                 "KICH" = unique(as.character(unlist(subtype_df[, grepl(
                   "Histological.Subtype", colnames(subtype_df), ignore.case = T)]))),
                 "LIHC" = unique(as.character(unlist(subtype_df[, grepl(
                   "iCluster", colnames(subtype_df), ignore.case = T)]))),
                 "LUAD" = unique(as.character(unlist(subtype_df[, grepl(
                   "iCluster", colnames(subtype_df), ignore.case = T)]))),
                 "LUSC" = unique(as.character(unlist(subtype_df[, grepl(
                   "Expression.Subtype", colnames(subtype_df), ignore.case = T)]))),
                 "PAAD" = unique(as.character(unlist(subtype_df[, grepl(
                   "Histological type by RHH", colnames(subtype_df), ignore.case = T)]))),
                 "PCPG" = unique(as.character(unlist(subtype_df[, grepl(
                   "mRNA Subtype Clusters", colnames(subtype_df), ignore.case = T)]))),
                 "SARC" = unique(as.character(unlist(subtype_df[, grepl(
                   "short histo", colnames(subtype_df), ignore.case = T)]))),
                 "STAD" = unique(as.character(unlist(subtype_df[, grepl(
                   "Molecular.Subtype", colnames(subtype_df), ignore.case = T)]))),
                 "THCA" = unique(as.character(unlist(subtype_df[, colnames(subtype_df) == 
                                                                  "mRNA_Cluster_number"]))),
                 "UCEC" = unique(as.character(unlist(subtype_df[, grepl(
                   "msi", colnames(subtype_df), ignore.case = T)]))),
                 "UCS" = unique(as.character(unlist(subtype_df[, colnames(subtype_df) == 
                                                                 "histologic subtype"]))),
                 "UVM" = unique(as.character(unlist(subtype_df[, grepl(
                   "mRNA Cluster", colnames(subtype_df), ignore.case = T)]))))
    st <- unique(c(st, "NA"))
  }
  return(st)
}


############################################################ 
#' Takes a vector of patient total number of mutations and creates a DF 
#' with separate columns (one for each of the bins) consisting of 0 and 1 
#' values for each patient (0 if they are not in that bin, 1 if they are)
#' @param total_num_muts a vector of the total number of mutations for each
#' patients' primary tumor sample
#' @param dataset either 'tcga' or 'metabric'
convert_tot_num_mut_to_bins <- function(total_num_muts, dataset) {
  
  # For TCGA, get only every second item in the vector (patients all have two 
  # almost identical entries)
  if(dataset == "tcga") {
    total_num_muts <- total_num_muts[seq(1, length(total_num_muts), 2)]
  }
  
  # Create data frame (columns are bins, rows are samples)
  num_bins <- 3
  tot_mut_df <- data.frame(matrix(nrow = length(total_num_muts), ncol = num_bins))
  colnames(tot_mut_df) <- c("Tot_Mut_b1", "Tot_Mut_b2", "Tot_Mut_b3")
  
  # For each item in the vector, determine which bin the mutation number 
  # goes in and add it
  thresh <- c(30, 60)
  if (dataset == "metabric") {
    #thresh <- c(4, 8)
    # If using log2(TMB+1), the thresholds are different
    thresh <- c(2.5, 3.5)
  }
  for (i in 1:length(total_num_muts)) {
    tnm <- total_num_muts[i]
    
    tot_mut_df[i,] <- rep(0, num_bins)    # set all cols for that row equal to 0
    
    if (tnm <= thresh[1]) {
      tot_mut_df[i,"Tot_Mut_b1"] <- 1     # set this bucket for that row equal to 1
    } else if (tnm > thresh[1] & tnm <= thresh[2]) {
      tot_mut_df[i,"Tot_Mut_b2"] <- 1     # set this bucket for that row equal to 1
    } else if (tnm > thresh[2]) {
      tot_mut_df[i,"Tot_Mut_b3"] <- 1     # set this bucket for that row equal to 1
    } else {
      print("This is an invalid tot_num_mut. Please address this issue & re-try.")
    }
    
  }
  return(tot_mut_df)
}

############################################################ +-
#' Takes the number of rows in the patient data frame (# of patients), 
#' as well as the downstream target DF, and calculates how many rows 
#' the overall linear model input data frame will need
#' @param len_patient_df the number of rows in the patient data frame
#' @param downstream_targ_df the data frame of all regulatory proteins
#' at the given level of specificity (columns) and their known 
#' regulatory gene targets (entries)
get_total_num_rows <- function(len_patient_df, downstream_targ_df) {
  # Extract all the gene targets from the downstream target data frame
  all_targets <- c()
  for (i in 1:ncol(downstream_targ_df)) {
    all_targets <- c(targets, downstream_targ_df[,i])
  }
  all_targets <- all_targets[!is.na(all_targets)]
  
  # Get the number of these targets
  num_targets <- length(all_targets)
  
  # Multiply this by the number of patients, as we will examine
  # all of these targets in each patient
  total_num_rows <- num_targets * len_patient_df
  print(paste("TOTAL NUM. ROWS:", total_num_rows))
  
  return(total_num_rows)
}
