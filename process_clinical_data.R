############################################################
### Process Clinical Data 
### Written By: Sara Camilli, August 2020
############################################################

# This file processes the various clinical data files for each data type in order to:
# 1. Filter them to include only patients that have all types of data
# 2. Filter them to include only data columns that all or most patients possess

# NOTE: the clinical files have already been imported and converted
# to dataframes in "get_TCGA_barcodes.R", where they were saved with a TCGA
 # barcode column attached. 

library(TCGAbiolinks)
library(ggplot2)

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/"

input_path <- paste(path, "Input Data Files/BRCA Data/", sep = "")
#input_path <- paste(path, "Input Data Files/TCGA Data (ALL)/", sep = "")

############################################################
# RETAIN ONLY PATIENTS WITH ALL DATATYPES IN GIVEN CLINICAL FILE
############################################################

# Read in clinical files with TCGA barcode column (add this column by going through get_TCGA_barcodes.R)
mut_clin_file <- read.csv(paste(input_path, "Somatic_Mut_Data/mut_clinical_data.csv", sep = ""), 
                          header = TRUE, na.strings = "\'--")
exp_clin_file <- read.csv(paste(input_path, "Expression_Data/exp_clinical_data.csv", sep = ""),
                          header = TRUE, na.strings = "\'--")
# opt 1. for CNV (CNV scores)
cnv_clin_file <- read.csv(paste(input_path, "CNV_Data/Gene-Level Copy Number Scores/cnv_clinical_data.csv", sep = ""),
                          header = TRUE, na.strings = "\'--")
# opt 2. for CNV (raw CNV values)
cnv_clin_file <- read.csv(paste(input_path, "CNV_Data/Gene-Level Copy Number Raw/cnv_clinical_data.csv", sep = ""),
                          header = TRUE, na.strings = "\'--")
meth_clin_file <- read.csv(paste(input_path, "Methylation_Data/meth_clinical_data.csv", sep = ""),
                                 header = TRUE, na.strings = "\'--")

# Read in all the unique patient IDs from saved file, if using a specific TCGA project
# NOTE: no way to use TCGA biolinks to get pan-cancer information. Have to just get all the clinical files and 
# find the intersecting IDs. Save this to file. Otherwise, file originated from analyses in TCGAbiolinks_assay.R
unique_patient_ids <- get_common_patients(input_path, mut_clin_file, exp_clin_file, cnv_clin_file, meth_clin_file) # Length: 7,864; OR
unique_patient_ids <- read.table(paste(input_path, "unique_brca_patient_ids.txt", sep = ""), header = FALSE)[,1] # Length: 747 (OPT. 1) or 732 (OPT. 2)

#' Takes clinical files (with a column for tcga_barcode) and identifies 
#' a list of unique TCGA patient identifiers common across all of them. 
#' Writes this list to a file and returns it.
#' @param path the main path to all the clinical files
#' @param mut_clin_file mutation clinical file from the TCGA
#' @param exp_clin_file expression clinical file from the TCGA
#' @param cnv_clin_file copy number variation clinical file from the TCGA
#' @param meth_clin_file methylation clinical file from the TCGA
get_common_patients <- function(path, mut_clin_file, exp_clin_file, cnv_clin_file, meth_clin_file) {
  # Extract all the IDs
  mut_ids <- unlist(lapply(mut_clin_file$case_submitter_id, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))
  exp_ids <- unlist(lapply(exp_clin_file$case_submitter_id, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))
  cnv_ids <- unlist(lapply(cnv_clin_file$case_submitter_id, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))
  meth_ids <- unlist(lapply(meth_clin_file$case_submitter_id, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))
  
  # Get just the overlapping IDs
  intersecting_ids <- intersect(intersect(intersect(mut_ids, exp_ids), cnv_ids), meth_ids)
  
  # Write these IDs to a file
  #write.table(intersecting_ids, paste(input_path, "unique_brca_patient_ids.txt", sep = ""))
  #write.table(intersecting_ids, paste(input_path, "unique_patient_ids_2.txt", sep = ""))
  
  # Return the IDs
  return(intersecting_ids)
}

#' Creates a subsetted clinical file from a given clinical file, containing 
#' only the given patients
#' @param clin_file a clinical file (with a column for tcga_barcode)
#' @param unique_patient_ids a list of patient IDs for patients that have all data types of interest
subset_clin_file_by_patient_ids <- function(clin_file, unique_patient_ids) {
  clin_df_sub <- data.frame(matrix(nrow = 0, ncol = ncol(clin_file)))
  colnames(clin_df_sub) <- colnames(clin_file)
  
  for (i in 1:nrow(clin_file)) {
    barcode <- clin_file$case_submitter_id[i]
    patient_id <- unlist(strsplit(barcode, split = "-", fixed = TRUE))[3]
    if (patient_id %fin% unique_patient_ids) {
      clin_df_sub <- rbind(clin_df_sub, clin_file[i,])
    }
  }
  return(clin_df_sub)
}

# Use helper function to subset the clinical files by patient IDs for 
# patients that have all data types of interest
mut_clin_df_sub <- subset_clin_file_by_patient_ids(mut_clin_file, unique_patient_ids)
  # Original BRCA file: 2089, Subset BRCA file: 1494 (all patients have 2 samples)
  # Original Pan-Cancer file: 19730, Subset Pan-Cancer file: 15728
exp_clin_df_sub <- subset_clin_file_by_patient_ids(exp_clin_file, unique_patient_ids)
  # Original BRCA file: 2182, Subset BRCA file: 1494 (all patients have 2 samples)
  # Original Pan-Cancer file: 21642, Subset Pan-Cancer file: 17222
cnv_clin_df_sub <- subset_clin_file_by_patient_ids(cnv_clin_file, unique_patient_ids)
  # OPT 1: Original BRCA file: 2176, Subset BRCA file: 1492 (all but 2 patients have 2 samples)
  # OPT 2: Original BRCA file: 2110, Subset BRCA file: 1464 (all patients have 2 samples)
  # Original Pan-Cancer file: 21968, Subset Pan-Cancer file: 15728
meth_clin_df_sub <- subset_clin_file_by_patient_ids(meth_clin_file, unique_patient_ids)
  # Original BRCA file: 1581, Subset BRCA file: 1490 (all but 4 patients have 2 samples)
  # Original Pan-Cancer file: 16779, Subset Pan-Cancer file: 15731


############################################################
# RETAIN ONLY FEMALE PATIENTS (BRCA DATA ONLY)
############################################################
mut_clin_df_sub <- mut_clin_df_sub[mut_clin_df_sub$gender == "female",]
  # Number of rows: 1477
exp_clin_df_sub <- exp_clin_df_sub[exp_clin_df_sub$gender == "female",]
  # Number of rows: 1478 (OPT.1); 1464 (OPT.2)
cnv_clin_df_sub <- cnv_clin_df_sub[cnv_clin_df_sub$gender == "female",]
  # Number of rows: 1478
meth_clin_df_sub <- meth_clin_df_sub[meth_clin_df_sub$gender == "female",]
  # Number of rows: 1479

############################################################
# RETAIN ONLY DATA COLUMNS WITH >N% NON-NULL ENTRIES
############################################################
#' Takes in a clinical data frame and returns it subsetted
#' to include only columns that have >N% non-null entries
#' @param clin_df the clinical data frame to be subsetted
#' @param N the percentage of non-null entries needed to 
#' retain a given column
retain_nonnull_cols <- function(clin_df, N) {
  # Get N% of the number of rows
  n_perc <- N * nrow(clin_df)  # ex. 747
  # Subset the dataframe using this
  clin_df_subset <- clin_df[,colSums(is.na(clin_df)) < n_perc]
  return(clin_df_subset)
}

N <- 0.5 # At N <- 0.5, all went from 155 to 35 columns

mut_clin_df_sub <- retain_nonnull_cols(mut_clin_df_sub, N)
exp_clin_df_sub <- retain_nonnull_cols(exp_clin_df_sub, N)
cnv_clin_df_sub <- retain_nonnull_cols(cnv_clin_df_sub, N)
meth_clin_df_sub <- retain_nonnull_cols(meth_clin_df_sub, N)

# Write these subsetted clinical files to CSV
write.csv(mut_clin_df_sub, paste(input_path, "Somatic_Mut_Data/mut_clinical_data_subset.csv", sep = ""))
write.csv(exp_clin_df_sub, paste(input_path, "Expression_Data/exp_clinical_data_subset.csv", sep = ""))
write.csv(cnv_clin_df_sub, paste(input_path, "CNV_Data/Gene-Level Copy Number Score/cnv_clinical_data_subset.csv", sep = "")) # OPT. 1 
write.csv(cnv_clin_df_sub, paste(input_path, "CNV_Data/Gene-Level Copy Number Raw/cnv_clinical_data_subset.csv", sep = "")) # OPT. 2 
write.csv(meth_clin_df_sub, paste(input_path, "Methylation_Data/meth_clinical_data_subset.csv", sep = ""))

# Print the remaining columns with greater than N% non-null entries
colnames(mut_clin_df_sub)

# We will only use one clinical file for the linear model, since these should all
# now be mostly identical (excepting a few patients lacking a matched tumor-normal for
# a particular datatype) -- we will use the most conservative file
clinical_df <- mut_clin_df_sub 
write.csv(clinical_df, paste(input_path, "clinical_data_subset.csv", sep = ""))

# When retrieving this information, use:
clinical_df <- read.csv(paste(input_path, "clinical_data_subset.csv", sep = ""), header = TRUE)


############################################################
# VISUALIZE THE NUMBER OF PATIENT CANCER SAMPLES PER CANCER TYPE
############################################################
#' Creates a bar plot to visualize the number of patients within each cancer type
#' category, to decide which ones we should run individually
#' @param clinical_df a clinical data frame from the TCGA
visualize_num_samps_per_ct <- function(clinical_df) {
  
  # Get the number of unique cancer types
  unique_cts <- unique(unlist(lapply(unlist(clinical_df$project), function(x) 
    unlist(strsplit(x, "-", fixed = TRUE))[2])))
  
  # Count up the number of patients that fall into each CT
  count_list <- lapply(unique_cts, function(ct) {
    num_ct <- length(unique(clinical_df[grepl(ct, clinical_df$project), 'case_submitter_id']))
    return(num_ct)
  })
  names(count_list) <- unique_cts
  count_mat <- data.frame("cancer_type" = names(count_list), "patient_count" = as.integer(unlist(count_list)))
  print(head(count_mat))
  
  # Create an ordered bar chart for visualization
  count_mat$cancer_type <- factor(count_mat$cancer_type, levels = count_mat$cancer_type[order(count_mat$patient_count, 
                                                                                              decreasing = FALSE)])
  ggplot(count_mat, aes(x=cancer_type, y=patient_count)) + geom_bar(stat = "identity") + coord_flip()
}


# Import the clinical DF, if needed
clinical_df <- read.csv(paste(input_path, "clinical_wMutCounts.csv", sep = ""), 
                        header = TRUE, check.names = FALSE)

visualize_num_samps_per_ct(clinical_df)


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

# ADDITIONAL OPTION USING TCGA BIOLINKS:
brca_clin <- GDCquery_clinic(project = "TCGA-BRCA", type = "Clinical")
brca_query_clin <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical")
GDCdownload(brca_query_clin)
brca_clin_prep <- GDCprepare(brca_query_clin)

# Bind all clinical files together
tcga_projects <- c("LAML", "ACC", "BCLA", "LGG", "BRCA", "CESC", "CHOL", "LCML", 
                   "COAD", "CNTL", "ESCA", "FPPP", "GBM", "HNSC", "KICH", "KIRC",
                   "KIRP", "KIHC", "KIRC", "KIRP", "KIHC", "LUAD", "LUSC", "DLBC",
                   "MESO", "MISC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC",
                   "SKCM", "STAD", "TGCT", "THYM", "THCA", "UCS", "UCEC", "UVM")



