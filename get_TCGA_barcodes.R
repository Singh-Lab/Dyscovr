############################################################
### GET TCGA BARCODES FROM UUIDS
### Written By: Sara Camilli, July 2020
############################################################

# BiocManager::install("TCGAutils")
library(TCGAutils)

# This file reads in all clinical files from GDC for a particular cancer type
# and extracts the UUIDs. It then uses TCGAUtils package to convert all UUIDs to
# TCGA barcodes, adding a column for TCGA barcodes to the clinical files and 
# creating an output CSV containing matching TCGA barcodes and UUIDs.


path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/"
#path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/"


############################################################
### READ IN CLINICAL DATA FILES 
############################################################
# Somatic Mutation Clinical File
mut_clin_file <- read.csv(paste(path, "Somatic_Mut_Data/mut_clinical_data.csv", sep = ""), 
                          header = TRUE, check.names = FALSE) 

# RNA-Seq Clinical File
exp_clin_file <- read.csv(paste(path, "Expression_Data/Counts/exp_clinical_data.csv", sep = ""), 
                          header = TRUE, check.names = FALSE)
#colnames(exp_clin_file)[1] <- "case_id"

# CNV Clinical File
cnv_clin_file <- read.csv(paste(path, "CNV_Data/Gene-Level Copy Number Raw/cnv_clinical_data.csv", sep = ""), 
                          header = TRUE, check.names = FALSE)
#colnames(cnv_clin_file)[1] <- "case_id"

# Methylation Clinical File
meth_clin_file <- read.csv(paste(path, "Methylation_Data/meth_clinical_data.csv", sep = ""), 
                           header = TRUE, check.names = FALSE)
#colnames(meth_clin_file)[1] <- "case_id"


############################################################
### CREATE TCGA BARCODE-UUID CONVERSION FILE
############################################################
#' This function takes in a clinical file data frame, along with an output file
#' path and filename, and converts each UUID to a TCGA ID. It then adds a 
#' column in the data frame for this TCGA ID in the original clinical DF 
#' AND outputs a two-column UUID-TCGA conversion CSV to the desired path
#' @param clin_file clinical file produced from the TCGA for the given cohort
#' @param output_path a path to write the barcode-UUID conversion document to
create_barcode_uuid_conv <- function(clin_file, output_path) {
  
  # Extract all the file's UUIDs from the clinical file and 
  # convert them to a list of characters
  file_uuids <- unique(as.character(clin_file$case_id))   # may also be manifest$entity_id
  
  # Create a data frame to store all the UUIDs and corresponding TCGA Barcodes
  uuids_and_barcodes = data.frame(matrix(nrow = length(file_uuids), ncol = 2))
  colnames(uuids_and_barcodes) <- c("uuid", "tcga_barcode")
  uuids_and_barcodes$uuid <- file_uuids
  
  uuids_and_barcodes$tcga_barcode <- unlist(lapply(file_uuids, function(uuid) {
    print(uuid)
    tryCatch(
      expr = {
        # Use function from TCGAutils package
        tcga_barcode <- UUIDtoBarcode(uuid, from_type = "case_id")
        return(tcga_barcode[1,2])
      },
      error = function(e) {
        print(paste("Unable to convert UUID ", uuid))
        return(NA)
      }
    )
  }))

  
  # Save this as a separate file
  write.csv(uuids_and_barcodes, output_path)
  return(uuids_and_barcodes)
}

mut_uuids_and_barcodes <- create_barcode_uuid_conv(mut_clin_file, 
                                                   output_path = paste(path, "Somatic_Mut_Data/mut_uuid_barcode_conversion.csv", sep = ""))
exp_uuids_and_barcodes <- create_barcode_uuid_conv(exp_clin_file, 
                                                   output_path = paste(path, "Expression_Data/exp_uuid_barcode_conversion.csv", sep = ""))
cnv_uuids_and_barcodes <- create_barcode_uuid_conv(cnv_clin_file, 
                                                   output_path = paste(path, "CNV_Data/cnv_uuid_barcode_conversion.csv", sep = ""))
meth_uuids_and_barcodes <- create_barcode_uuid_conv(meth_clin_file, 
                                                    output_path = paste(path, "Methylation_Data/meth_uuid_barcode_conversion.csv", sep = ""))

# Add list of TCGA IDs as a column in the original clinical file
mut_clin_file$tcga_barcode <- mut_uuids_and_barcodes$tcga_barcode
exp_clin_file$tcga_barcode <- exp_uuids_and_barcodes$tcga_barcode
cnv_clin_file$tcga_barcode <- cnv_uuids_and_barcodes$tcga_barcode
meth_clin_file$tcga_barcode <- meth_uuids_and_barcodes$tcga_barcode

# OPT: Save over the original clinical file with this barcode column added
write.csv(mut_clin_file, paste(path, "Somatic_Mut_Data/mut_clinical_data.csv", sep = ""))
write.csv(exp_clin_file, paste(path, "Expression_Data/exp_clinical_data.csv", sep = ""))
write.csv(cnv_clin_file, paste(path, "CNV_Data/cnv_clinical_data.csv", sep = ""))
write.csv(meth_clin_file, paste(path, "Methylation_Data/meth_clinical_data.csv", sep = ""))

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
### NOTE: THIS FUNCTION IS INACCURATE AND SHOULD NO LONGER BE USED ###
#
#
############################################################
### INFER NUMBER OF PATIENTS WITH TUMOR-NORMAL PAIRINGS 
############################################################
#' Now that we have all the file TCGA barcodes, we want to match TCGA 
#' barcodes in order to infer tumor-normal pairings. This function
#' extracts the "participant" portion of the barcode
#' @param uuids_mapping a data frame that matches UUIDs and TCGA barcodes
get_num_matches <- function(uuids_mapping) {
  seen_indiv <- c()     # A list of unique individuals we've already seen
  num_matches <- 0
  for (i in 1:nrow(uuids_mapping)) {
    barcode <- uuids_mapping[i,2]
    indiv_id <- unlist(strsplit(barcode, "-"))
    if (!TRUE %in% is.na(indiv_id)) {
      if (indiv_id[1] == "TCGA") {
        indiv_id <- indiv_id[3]
        if (indiv_id %in% seen_indiv) {
          num_matches <- num_matches + 1
        } else {
          seen_indiv <- c(seen_indiv, indiv_id)
        }
      }
    }
  }
  return(num_matches)
}

# Mutation
num_matches_mut <- get_num_matches(mut_uuids_and_barcodes)      # Will hold the number of healthy-tumor matched samples we have
print(paste("Number of Healthy-Tumor Matches, Mutation Data:", num_matches_mut))

# Expression
num_matches_exp <- get_num_matches(exp_uuids_and_barcodes)      # Will hold the number of healthy-tumor matched samples we have
print(paste("Number of Healthy-Tumor Matches, Expression Data:", num_matches_exp))

# CNV
num_matches_cnv <- get_num_matches(cnv_uuids_and_barcodes)      # Will hold the number of healthy-tumor matched samples we have
print(paste("Number of Healthy-Tumor Matches, CNV Data:", num_matches_cnv))

num_matches_meth <- get_num_matches(meth_uuids_and_barcodes)      # Will hold the number of healthy-tumor matched samples we have
print(paste("Number of Healthy-Tumor Matches, Methylation Data:", num_matches_meth))
