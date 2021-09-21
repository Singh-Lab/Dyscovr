############################################################
### Create Patient Data Frame
### Written By: Sara Camilli, July 2020
############################################################

# This file uses a clinical data frame, and a variable declaring whether the data of interest
# is or is not breast cancer, to create an output "patient data frame" that will be used in the linear
# model. This data frame has the following format:
  # Columns: Linear model input variables (Gender, Age (buckets 1-10), Race (buckets 1-5), 
    # Prior malignancies, Treatment-radiation, Treatment-pharmacological, Cancer Type or Subtype,
    # and first 2 PCs generated from genotype matrices)
  # Rows: Patients that have all necessary data types 1...j...L (Unique 4-digit ID)
  # Entries: 0 or 1 values

library(TCGAbiolinks)

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

input_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/"
#input_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/"


# Import clinical data
# Tumor-Normal Matched
clinical_df_tnm <- read.csv(paste(input_path, "clinical_data_subset.csv", sep = ""), 
                        header = TRUE, row.names = 1, check.names = FALSE)
# Non-Tumor-Normal Matched
clinical_df_ntnm <- read.csv(paste(input_path, "clinical_wMutCounts.csv", sep = ""), 
                        header = TRUE, check.names = FALSE)

is_brca <- TRUE
# is_brca <- FALSE

is_pc <- FALSE
# is_pc <- TRUE

general_cancer_types <- c("ACC", "BRCA", "BLCA", "CESC", "CHOL", "COAD", "ESCA",
                          "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC",
                          "LUAD", "LUSC", "PAAD", "PCPG", "PRAD", "READ", "SKCM",
                          "SARC", "STAD", "THCA", "UCEC", "UCS", "UVM")

import_subtype_information <- function(general_cancer_types) {
  subtype_files <- lapply(general_cancer_types, function(type) {
    subtype <- TCGAquery_subtype(tumor = type)
    return(subtype)
  })
  names(subtype_files) <- general_cancer_types
  
  # Note: To combine into one file, we'll have to do some work to ID matching 
  # columns with differing names
  
  return(subtype_files)
}

# Pan-Cancer subtype info
pc_subtype <- import_subtype_information(general_cancer_types)

# Or, if just using BRCA:
brca_subtype <- TCGAquery_subtype(tumor = "BRCA")


# Import PCA information
### BRCA ###
# Tumor-Normal Matched
brca_smartpca_output_tnm <- read.table(paste(main_path, "smartpca/Tumor-Normal Matched/BRCA_output_eigenvect.pca", sep = ""),
                                   header = FALSE, check.names = FALSE, row.names = 1)
# Non-Tumor-Normal Matched
brca_smartpca_output_ntnm <- read.table(paste(main_path, "smartpca/Non-Tumor-Normal Matched/BRCA_output_eigenvect.pca", sep = ""),
                                   header = FALSE, check.names = FALSE, row.names = 1)
### PAN-CANCER ###
# Tumor-Normal Matched
pc_smartpca_output_tnm <- read.table(paste(main_path, "smartpca/Tumor-Normal Matched/pc_output_eigenvect.pca", sep = ""),
                                   header = FALSE, check.names = FALSE, row.names = 1)
# Non-Tumor-Normal Matched
pc_smartpca_output_ntnm <- read.table(paste(main_path, "smartpca/Non-Tumor-Normal Matched/pc_output_eigenvect.pca", sep = ""),
                                   header = FALSE, check.names = FALSE, row.names = 1)


# Create data frame with all patient characteristics
brca_patient_df_tnm <- create_patient_dataframe(clinical_df_tnm, is_brca, 
                                                  brca_subtype, is_pc, brca_smartpca_output_tnm)
brca_patient_df_ntnm <- create_patient_dataframe(clinical_df_ntnm, is_brca, 
                                                  brca_subtype, is_pc, brca_smartpca_output_ntnm)
pc_patient_df_tnm <- create_patient_dataframe(clinical_df_tnm, is_brca, 
                                                pc_subtype, is_pc, pc_smartpca_output_tnm)
pc_patient_df_ntnm <- create_patient_dataframe(clinical_df_ntnm, is_brca, 
                                                 pc_subtype, is_pc, pc_smartpca_output_ntnm)


# Write this to a CSV
write.csv(brca_patient_df_tnm, paste(main_path, "Linear Model/Patient and Sample DFs/patient_dataframe_tnm.csv", sep = ""))
write.csv(brca_patient_df_ntnm, paste(main_path, "Linear Model/Patient and Sample DFs/patient_dataframe_ntnm.csv", sep = ""))
write.csv(pc_patient_df_tnm, paste(main_path, "Linear Model/Patient and Sample DFs/patient_dataframe_tnm.csv", sep = ""))
write.csv(pc_patient_df_tnm, paste(main_path, "Linear Model/Patient and Sample DFs/patient_dataframe_ntnm.csv", sep = ""))


############################################################

#' Creates a skeleton "patient data frame" that will contain
#' all the characteristics of interest for each patient. Calls helper 
#' function to fill patient data frame and then returns it.
#' @param clinical_df clinical data file from the TCGA for the given patient cohort
#' @param is_brca a TRUE/FALSE value indicating whether the cohort is breast
#' cancer or not
#' @param is_pc a TRUE/FALSE value indicating whether the cohort is pan-cancer
#' or not
#' @param subtype_file a file from TCGAbiolinks that has information about subtypes
#' or, alternatively (if pan-cancer) a file that identifies the cancer type of 
#' each sample
#' @param smartpca_output a file with output values for the top 10 PCs for each 
#' sample, created from running smartpca
create_patient_dataframe <- function(clinical_df, is_brca, subtype_file, is_pc,
                                     smartpca_output) {
  
  # Create relational data frame to hold patient-dependent inputs to linear model
  # Columns: Linear model input variables 
  # Rows: Samples that have all necessary data types 1...j...L (Unique 4-digit ID)
  patient_characteristics <- c("Gender", "Age_b1", "Age_b2", "Age_b3", "Age_b4", "Age_b5", "Age_b6",
                               "Age_b7", "Age_b8", "Age_b9", "Age_b10", "Race_b1", "Race_b2", "Race_b3",
                               "Race_b4", "Race_b5", "Prior_malig", "Treatment_rad", "Treatment_pharm",
                               "Tot_Mut_b1", "Tot_Mut_b2", "Tot_Mut_b3", "PC1", "PC2")
  patient_characteristics <- add_cancer_type(patient_characteristics, is_pc, subtype_file)
  unique_patient_ids <- unique(unlist(lapply(clinical_df$case_submitter_id, function(x) 
    unlist(strsplit(x, split = "-", fixed = TRUE))[3])))
  patient_dataframe <- data.frame(matrix(ncol = length(patient_characteristics), 
                                         nrow = length(unique_patient_ids)))
  colnames(patient_dataframe) <- patient_characteristics
  rownames(patient_dataframe) <- unique_patient_ids
  
  # Input all the patient-specific information 
  patient_dataframe <- input_patient_specific_info(patient_dataframe, clinical_df, 
                                                   is_brca, subtype_file, is_pc,
                                                   smartpca_output)
  
  return(patient_dataframe)
}

############################################################

#' Takes the vector of patient characteristics and adds buckets for the 
#' appropriate number of subtypes (for a specific cancer type) or types 
#' (for pan-cancer analysis). Returns the updated vector.
#' @param patient_characterisitics a vector of existing patient characteristics
#' to be updated with cancer type/subtype
#' @param is_pc a TRUE/FALSE value indicating whether the cohort is pan-cancer or not
#' @param subtype_file a file from TCGAbiolinks that has information about subtypes
#' or, alternatively (if pan-cancer) a file that identifies the cancer type of 
#' each sample (TODO: make this file!)
add_cancer_type <- function(patient_characteristics, is_pc, subtype_file) {
  # Pan-Cancer: look at cancer type
  if (is_pc) {
    # Find how many unique cancer types there are in total
    unique_cts <- unique(subtype_file$cancer.type)
    length_cts <- length(unique_cts)
  }
  # Not pan-cancer: look at cancer subtype
  else {
    # Find the subtype column and check how many types are given
    unique_subtypes <- unique(unlist(subtype_file[ ,grepl("subtype", colnames(subtype_file),  
                                                      ignore.case = TRUE)]))
    unique_subtypes <- unique_subtypes[!(is.na(unique_subtypes) | unique_subtypes == "NA")]
    length_cts <- length(unique_subtypes)
    
  }
  # Create buckets equivalent to the number of types or subtypes
  new_buckets <- unlist(lapply(1:length_cts, function(x) paste("Cancer_type_b", x, sep = "")))
  
  # Add new buckets to the patient characteristics list and return
  patient_characteristics <- c(patient_characteristics, new_buckets)
  
  return(patient_characteristics)
}

############################################################

############################################################

#' Takes the empty patient data frame, along with the clinical data frame and 
#' whether or not we are looking at BRCA, and updates patient DF with all the 
#' patient-specific information: gender (if not BRCA), age, race, prior 
#' malignancies, treatment, and total number of mutations. 
#' Returns data frame with these fields filled.
#' @param input_df a blank data frame with the dimensions of the the number
#' of patient characteristics (columns) and number of patients in the 
#' cohort (rows)
#' @param clinical_df clinical data file from the TCGA for the given patient cohort
#' @param is_brca a TRUE/FALSE value indicating whether the cohort is breast
#' cancer or not
#' @param subtype_file a file from TCGAbiolinks that has information about subtypes
#' or, alternatively (if pan-cancer) a file that identifies the cancer type of 
#' each sample (TODO: make this file!)
#' @param is_pc a TRUE/FALSE value indicating whether the cohort is pan-cancer or not
#' @param smartpca_output a file with output values for the top 10 PCs for each 
#' sample, created from running smartpca
input_patient_specific_info <- function(input_df, clinical_df, is_brca, subtype_file, 
                                        is_pc, smartpca_output) {
  input_dataframe_updated <- input_df
  
  # GENDER
  if (!is_brca) {
    gender_vect <- clinical_df$gender[seq(1, length(clinical_df$gender), 2)]   
    gender_vect_binary <- ifelse(gender_vect == "male", 0, 1)
    input_dataframe_updated[,'Gender'] <- gender_vect_binary
  } 
  # AGE
  age_bin_df <- convert_age_to_bins(clinical_df$age_at_index)
  input_dataframe_updated <- merge_dataframes_by_colname(input_dataframe_updated, 
                                                         age_bin_df)   
  # RACE
  race_bin_df <- convert_race_to_bins(clinical_df$race, clinical_df$ethnicity)
  input_dataframe_updated <- merge_dataframes_by_colname(input_dataframe_updated, 
                                                         race_bin_df) 
  # PRIOR_MALIG
  prior_malig_vect <- clinical_df$prior_malignancy[seq(1, length(clinical_df$prior_malignancy), 2)]
  input_dataframe_updated[,'Prior_malig'] <- ifelse(prior_malig_vect == "no" | 
                                                      prior_malig_vect == "not reported", 0, 1)
  # TREATMENT
  treatment_or_ther_vect <- clinical_df$treatment_or_therapy
  treatment_type_vect <- clinical_df$treatment_type
  treatment_bin_df <- convert_treatments_to_bins(treatment_or_ther_vect, 
                                                 treatment_type_vect)
  input_dataframe_updated <- merge_dataframes_by_colname(input_dataframe_updated, 
                                                         treatment_bin_df)
  # CANCER TYPE OR SUBTYPE
  cancer_type_bin_df <- convert_ct_to_bins(subtype_file, input_dataframe_updated, is_pc)
  input_dataframe_updated <- merge_dataframes_by_colname(input_dataframe_updated,
                                                         cancer_type_bin_df)
  # TOTAL NUM MUTATIONS
  tot_num_mut_df <- convert_tot_num_mut_to_bins(clinical_df$Total.Num.Muts)
  input_dataframe_updated <- merge_dataframes_by_colname(input_dataframe_updated, 
                                                         tot_num_mut_df)
  # SMARTPCA PCs
  smartpca_output <- reorg_smartpca_output(smartpca_output)
  smartpca_output <- pad_nas(smartpca_output, rownames(input_dataframe_updated))
  print(dim(smartpca_output))
  input_dataframe_updated <- merge_dataframes_by_colname(input_dataframe_updated,
                                                         smartpca_output)

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
convert_age_to_bins <- function(age_vector) {
  
  # Get only every second item in the vector (patients all have two almost identical entries)
  age_vector <- age_vector[seq(1, length(age_vector), 2)]
  
  # Create data frame (columns are bins, rows are patients)
  num_bins <- 10
  age_df <- data.frame(matrix(nrow = length(age_vector), ncol = num_bins))
  colnames(age_df) <- c("Age_b1", "Age_b2", "Age_b3", "Age_b4", "Age_b5", "Age_b6",
                        "Age_b7", "Age_b8", "Age_b9", "Age_b10")
  for (i in 1:length(age_vector)) {
    age <- age_vector[i]
    age_df[i,] <- rep(0, num_bins)              # set all cols for that row equal to 0
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

#' Takes vectors of patient races and ethnicities and creates separate columns 
#' (one for each of the bins) consisting of 0 and 1 values for 
#' each patient (0 if they are not in that race/ ethnicity bin, 1 if they are)
#' @param race_vector a vector of patient races (either 'white', 'black or african
#' american', 'asian', or other)
#' @param ethnicity_vector a vector of patient ethnicities (either 'hispanic or
#' latino' or 'not hispanic or latino')
convert_race_to_bins <- function(race_vector, ethnicity_vector) {
  
  # Get only every second item in the vector (patients all have two almost identical entries)
  race_vector <- race_vector[seq(1, length(race_vector), 2)]
  ethnicity_vector <- ethnicity_vector[seq(1, length(ethnicity_vector), 2)]
  
  # Create data frame (columns are bins, rows are patients)
  num_bins <- 5
  race_df <- data.frame(matrix(nrow = length(race_vector), ncol = num_bins))
  colnames(race_df) <- c("Race_b1", "Race_b2", "Race_b3", "Race_b4", "Race_b5")
  
  for (i in 1:length(race_vector)) {
    race <- race_vector[i]
    ethnicity <- ethnicity_vector[i]
    race_df[i,] <- rep(0, num_bins)               # set all cols for that row equal to 0
    if (race == "white") {
      if (ethnicity == "not hispanic or latino") {
        race_df[i,"Race_b1"] <- 1      # set this bucket for that row equal to 1
      }
    } else if (race == "white") {
      if (ethnicity == "hispanic or latino") {
        race_df[i,"Race_b2"] <- 1      # set this bucket for that row equal to 1
      }
    } else if (race == "black or african american") {
      race_df[i, "Race_b3"] <- 1      # set this bucket for that row equal to 1
    } else if (race == "asian") {
      race_df[i, "Race_b4"] <- 1      # set this bucket for that row equal to 1
    } else {
      race_df[i, "Race_b5"] <- 1      # set this bucket for that row equal to 1
    }
  }
  return(race_df)
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
convert_treatments_to_bins <- function(treatment_or_ther_vect, treatment_type_vect) {
  
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
    treatment_df <- rad_or_pharm_helper(rad_or_pharm_1, treatment_1, treatment_df, i)
    treatment_df <- rad_or_pharm_helper(rad_or_pharm_2, treatment_2, treatment_df, i)
    
    # Update indices
    treatment_vect_ind <- treatment_vect_ind + 1
  }
  return(treatment_df)
}

############################################################ 

#' Helper function for converting treatment to bins. Takes a
#' treatment category variable (radiation or pharma) and a 
#' yes or no treatment variable and updates/returns given 
#' treatment data frame
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

#' Takes a file with information about cancer type or subtype
#' per patient and creates a DF with separate columns (one per bin)
#' consisting of 0 and 1 values for each patient (0 if they are not 
#' in that bin, 1 if they are)
#' @param subtype_file a file from TCGAbiolinks that has information about subtypes
#' or, alternatively (if pan-cancer) a file that identifies the cancer type of 
#' each sample (TODO: make this file!)
#' @param input_dataframe_updated the work-in-progress input data frame with all
#' the patient characteristics for the linear model
#' @param is_pc a TRUE/FALSE value indicating whether the cohort is pan-cancer or not
convert_ct_to_bins <- function(subtype_file, input_dataframe_updated, is_pc) {
  
  # Create a blank data frame with all the same patient IDs as the input data frame
  cancer_type_bin_df <- input_dataframe_updated[, grepl("Cancer_type", colnames(input_dataframe_updated))]
  
  # Keep a growing map of which type will correspond to which bucket
  ct_to_bucket_map <- data.frame(ct = rep(NA, ncol(cancer_type_bin_df)), 
                                 bucket_num = 1:ncol(cancer_type_bin_df))
  
  # Go through each patient and find their subtype/type info, then add to the appropriate bucket
  for (i in 1:nrow(cancer_type_bin_df)) {
    patient_id <- rownames(cancer_type_bin_df)[i]
    # Pan-cancer: get the cancer type
    if (is_pc) {
      ct <- as.character(unlist(subtype_file[grepl(patient_id, subtype_file$patient), 
                                             'cancer.type']))
    } 
    # Not pan-cancer: get the cancer subtype
    else {
      ct <- as.character(unlist(subtype_file[grepl(patient_id, subtype_file$patient), 
                                             grepl("subtype", colnames(subtype_file), 
                                                   ignore.case = TRUE)]))
    }

    # If this cancer type isn't already in the cancer type-to-bucket mapping, add it
    if(!(ct %in% ct_to_bucket_map$ct)) {
      # Add a new entry to the map (replace the next NA)
      pos_of_first_NA <- min(which(is.na(ct_to_bucket_map$ct)))
      ct_to_bucket_map$ct[pos_of_first_NA] <- ct
    }
    # Get the bucket that it's in
    bucket_num <- ct_to_bucket_map[ct_to_bucket_map$ct == ct, 'bucket_num']
    # Create a new row where there is a '1' for this bucket and a '0' for 
    # everything else
    new_row <- rep(0, times = ncol(cancer_type_bin_df))
    new_row[bucket_num] <- 1
    # Update the blank DF for this patient
    cancer_type_bin_df[i,] <- new_row
  }
  
  # Return the newly binned DF (with patient IDs that match the order of the growing
  # input DF)
  return(cancer_type_bin_df)
}


############################################################ 

#' Takes a vector of patient total number of mutations and creates a DF 
#' with separate columns (one for each of the bins) consisting of 0 and 1 
#' values for each patient (0 if they are not in that bin, 1 if they are)
#' @param total_num_muts a vector of the total number of mutations for each
#' patients' primary tumor sample
convert_tot_num_mut_to_bins <- function(total_num_muts) {
  
  # Get only every second item in the vector (patients all have two almost identical entries)
  total_num_muts <- total_num_muts[seq(1, length(total_num_muts), 2)]
  
  # Create data frame (columns are bins, rows are samples)
  num_bins <- 3
  tot_mut_df <- data.frame(matrix(nrow = length(total_num_muts), ncol = num_bins))
  colnames(tot_mut_df) <- c("Tot_Mut_b1", "Tot_Mut_b2", "Tot_Mut_b3")
  
  # For each item in the vector, determine which bin the mutation number goes in and add it
  for (i in 1:length(total_num_muts)) {
    tnm <- total_num_muts[i]
    
    tot_mut_df[i,] <- rep(0, num_bins)              # set all cols for that row equal to 0
    
    if (tnm <= 30) {
      tot_mut_df[i,"Tot_Mut_b1"] <- 1      # set this bucket for that row equal to 1
    } else if (tnm > 30 & tnm <= 60) {
      tot_mut_df[i,"Tot_Mut_b2"] <- 1      # set this bucket for that row equal to 1
    } else if (tnm > 60) {
      tot_mut_df[i,"Tot_Mut_b3"] <- 1      # set this bucket for that row equal to 1
    } else {
      print("This is an invalid tot_num_mut. Please address this issue & re-try.")
    }

  }
  return(tot_mut_df)
}

############################################################ 

#' Takes a smartpca output data frame (rows are patient IDs and columns are
#' PCs) and averages the PC values across patient samples. Also reformats the 
#' row and column names so it can be added to the patient DF
#' @param smartpca_output a file with output values for the top 10 PCs for each 
#' sample, created from running smartpca
reorg_smartpca_output <- function(smartpca_output) {
  
  smartpca_output <- smartpca_output[,1:2]   # take only the first 2 PCs
  
  # First, average together the PCs from patients with multiple samples
  indic_of_dupl <- unlist(lapply(1:nrow(smartpca_output), function(i) {
    rownam <- rownames(smartpca_output)[i]
    length_spl <- length(unlist(strsplit(rownam, "-", fixed = TRUE)))
    if(length_spl > 3) {return(i)}
    else {return(NA)}
  }))
  indic_of_dupl <- indic_of_dupl[!is.na(indic_of_dupl)]
  
  seen_pats <- c()
  new_smartpca_output <- smartpca_output[-indic_of_dupl,]
  for (index in indic_of_dupl) {
    # For each duplicate, get the corresponding patient ID
    pat <- unlist(strsplit(rownames(smartpca_output)[index], "-", fixed = TRUE))[3]
    
    # If we haven't already handled this patient, average their rows
    if(!pat %in% seen_pats) {
      # Get all the rows with data for this patient
      pat_rows <- smartpca_output[grepl(pat, rownames(smartpca_output)),]
      pc_avgs <- unlist(colMeans(pat_rows))
      # Update the patient in the new output DF with the averages
      new_smartpca_output[grepl(pat, rownames(new_smartpca_output)),] <- pc_avgs
    }
    seen_pats <- c(seen_pats, pat)
  }
  
  # Then, fix the row and column names to match the input data frame
  rownames(new_smartpca_output) <- unlist(lapply(rownames(new_smartpca_output), function(x)
    unlist(strsplit(x, "-", fixed = TRUE))[3]))
  colnames(new_smartpca_output) <- c("PC1", "PC2")
  
  return(new_smartpca_output)
}


############################################################ 

#' Takes a smartpca output data frame (rows are patient IDs and columns are
#' PCs) creates false entries with "NA" values for missing patients
#' @param smartpca_output a file with output values for the top 10 PCs for each 
#' sample, created from running smartpca
#' @param patient_ids all the patient IDs in the patient DF
pad_nas <- function(smartpca_output, patient_ids) {
  
  # Remove patients in the smartpca output but not the patient list
  smartpca_output <- smartpca_output[rownames(smartpca_output) %fin% patient_ids,]
  
  # Which patients are in the patient ID list but not the smartpca output?
  nonoverlap_patients <- setdiff(patient_ids, rownames(smartpca_output))
  print(nonoverlap_patients)
  
  # Create a new entry for each of these with 'NA' values and add to smartpca df
  new_rows <- lapply(nonoverlap_patients, function(pat) {
    row <- data.frame("PC1" = NA, "PC2" = NA)
    rownames(row) <- pat
    return(row)
  })
  new_df <- rbindlist(new_rows)
  smartpca_output <- rbind(smartpca_output, new_df)
  
  return(smartpca_output)
}

############################################################ 

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
  
  return(total_num_rows)
}


