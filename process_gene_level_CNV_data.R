############################################################
### Process Gene-Level CNA Data 
### Written By: Sara Camilli, August 2020
############################################################

# This file processes the gene-level CNA (Score or Raw Value) data file(s) in order to:
  # 1. Filter to include only patients that have all types of data
  # 2. Filter to include only the transcription factors of interest

# Format of output dataframe: 
  # Rows : Transcription Factors (ENSG ID) OR all genes in genome
  # Columns : Patient TCGA ID (4-digit) OPT: (-Sample ID (01A is tumor, 11A is normal))
  # Entries : CNA value (positive represents an amplification event, negative represents deletion event,
    # 0 represents no amplification or deletion) OR the actual CNA value

library(TCGAbiolinks)
library(data.table)

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/"
#path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/"

output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/CNV/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/CNV/"


############################################################
############################################################
# OPTION 1: USE GISTIC2 SCORING FILES (AMPLIFIED, DELETED,
# OR UNALTERED)
############################################################
############################################################

############################################################
# IMPORT GISTIC2 AND RELATED FILE(S)
############################################################
# Import the CNV GISTIC2 file for particular project
#GISTIC_file <- read.table(paste(path, "CNV_Data/Gene-Level Copy Number Scores/BRCA.focal_score_by_genes.txt", sep = ""), header = TRUE, sep = "\t", 
                          #colClasses = "character")
  #  19729 x 1109 table (rows: ENSG ID gene names; columns: sample UUIDs)

# For pan-cancer, get the names of all the CNV files across cancers and merge them into one
GISTIC_filenames <- list.files(paste(path, "CNV_Data/", sep = ""), pattern = ".focal_score_by_genes.txt")

#' For pan-cancer analyses, this function aggregates all the GISTIC2 CNV
#' files across projects into one file
#' @param path a local path to the GISTIC filenames
#' @param GISTIC_filenames a vector a filenames for all the GISTIC2 files of interest
create_aggregate_cnv_file <- function(path, GISTIC_filenames) {
  # Create template from the first file
  GISTIC_file <- read.table(paste(path, paste("CNV_Data/", GISTIC_filenames[1], sep = ""), sep = ""), header = TRUE, sep = "\t", 
                            colClasses = "character", check.names = FALSE)

  # Merge each of the preceding CNV files to this first file
  for (i in 2:length(GISTIC_filenames)) {
    # Get the next file in the list 
    file <- read.table(paste(path, paste("CNV_Data/", GISTIC_filenames[i], sep = ""), sep = ""), header = TRUE, sep = "\t", 
                       colClasses = "character", check.names = FALSE)
    GISTIC_file <- merge(GISTIC_file, file, by = c("Gene Symbol", "Gene ID", "Cytoband"))
  }
  # Return aggregate file
  return(GISTIC_file)
}

GISTIC_file <- create_aggregate_cnv_file(GISTIC_filenames)

# Write this back to a CSV for future use
fwrite(GISTIC_file, paste(path, "CNV_Data/Aggregated.focal_score_by_genes.csv", sep = ""))

# Read back if this is already done
GISTIC_file <- read.csv(paste(path, "CNV_Data/Aggregated.focal_score_by_genes.csv", sep = ""))

# Also import the annotations file, if needed
# GISTIC_annotations <- read.table(paste(path, "CNV_Data/cnv_annotations.txt", sep = ""), header = TRUE, sep = "\t")

# From process CNV data:
# Import the list of patient IDs that have all four data types
#patient_id_list <- read.table(paste(path, "unique_brca_patient_ids.txt", sep = ""), header = TRUE)[,1]
patient_id_list <- read.table(paste(path, "unique_patient_ids.txt", sep = ""), header = FALSE)[,2]

# Import Biospecimen information to do the sample-patient ID conversions
#biosp_df <- read.csv(paste(path, "CNV_Data/Gene-Level Copy Number Scores/aliquot.csv", sep = ""), header = TRUE, na.strings = "\'--")
biosp_df <- read.csv(paste(path, "CNV_Data/aliquot.csv", sep = ""), header = TRUE, na.strings = "\'--")

# Import UUID-TCGA Barcode Conversion document
#uuid_barcode_conv <- read.csv(paste(path, "CNV_Data/Gene-Level Copy Number Scores/cnv_uuid_barcode_conversion.csv", sep = ""), header = TRUE)
uuid_barcode_conv <- read.csv(paste(path, "CNV_Data/cnv_uuid_barcode_conversion.csv", sep = ""), header = TRUE)


############################################################
# CONVERT GENE TARGETS LIST TO ENSG IDs (MATCH GISTIC2)
############################################################
# USE BIOMART R TOOL
# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)

protein_ids_iprotein <- read.csv(paste(output_path, "Mutation/swissprot_ids_missense_iprotein.csv", sep = ""), header = TRUE)[,2]
protein_ids_idomain <- read.csv(paste(output_path, "Mutation/swissprot_ids_missense_idomain.csv", sep = ""), header = TRUE)[,2]
protein_ids_ibindingpos <- read.csv(paste(output_path, "Mutation/swissprot_ids_missense_ibindingpos.csv", sep = ""), header = TRUE)[,2]

ensg_ids_iprotein <- unlist(lapply(protein_ids_iprotein, function(x) all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'ensembl_gene_id']))
ensg_ids_idomain <- unlist(lapply(protein_ids_idomain, function(x) all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'ensembl_gene_id']))
ensg_ids_ibindingpos <- unlist(lapply(protein_ids_ibindingpos, function(x) all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'ensembl_gene_id']))


############################################################
# FILTERING FUNCTION (RETAIN ONLY PATIENTS WITH ALL DATATYPES)
############################################################
#' Function takes a GISTIC2-style CNV filename, patient ID list, either
#' a regulatory protein list or a list of all genes in the human genome, 
#' a biospecimen file, and a UUID-TCGA barcode conversion table. Creates a 
#' data frame for the file that retains the CNV data for only 
#' those patients that also have all other data types of interest.
#' @param GISTIC_file filename of the given GISTIC2 file
#' @param patient_id_list vector of all patient IDs of interest
#' @param gene_list a vector of regulatory proteins (or all genes in the
#' human genome)
#' @param biosp_file an associated biospecimen file from the GDC
#' @param uuid_barcode_conv a conversion document linking UUIDs to TCGA IDs
filter_CNV <- function(GISTIC_file, patient_id_list, gene_list, biosp_file, uuid_barcode_conv) {
  GISTIC_df_sub <- data.frame(matrix(ncol = 0, nrow = nrow(GISTIC_file)))
  
  # Strip the versions off the end of the gene ids
  stripped_gene_ids <- strip_gene_ids(GISTIC_file$Gene.Symbol)

  # Subset columns by patient ID
  for (i in 4:ncol(GISTIC_file)) {
    # Extract the tumor sample ID and convert the delimiter from "." to "-", 
    # remove leading "X" if applicable using helper function
    tcga_id <- colnames(GISTIC_file)[i]
    new_id <- fix_samp_id(tcga_id)
    
    # Extract the corresponding case (patient) ID for this tumor sample
    if (TRUE %fin% grepl(new_id, biosp_file$aliquot_id)) {
      # Get the first index at which the tumor sample ID matches
      index <- grep(new_id, biosp_file$aliquot_id)[1]
    }
    else {
      print(paste("Unable to locate sample UUID ", new_id))
      next
    }
    # Get the corresponding case (patient) ID for this match
    case_id <- biosp_file$case_id[index]
    print(paste("case_id: ", case_id))
    
    # Get the corresponding TCGA barcode for this case UUID
    barcode <- unique(uuid_barcode_conv[uuid_barcode_conv$uuid == case_id, 'tcga_barcode'])
    print(paste("barcode: ", barcode))
    
    # Get the sample ID for this case
    samp <- biosp_file[biosp_file$case_id == case_id, 'portion_submitter_id']
    samp_id <- unlist(strsplit(samp, "-", fixed = TRUE))[4]
    print(paste("samp id: ", samp_id))
    
    # Extract just the patient identifier and sample identifier
    patient_id <- unlist(strsplit(barcode, split = "-", fixed = TRUE))[3]
    print(paste("patient_id:", patient_id))
    
    # If this ID is in the patient ID list, add this column to the new DF
    if (!length(patient_id) == 0) {
      if (patient_id %fin% patient_id_list) {
        GISTIC_df_sub <- cbind(GISTIC_df_sub, GISTIC_file[,i])
        print(head(GISTIC_df_sub))
        colnames(GISTIC_df_sub)[ncol(GISTIC_df_sub)] <- paste(patient_id, samp_id, sep = "-")
      }
    }
  }
  GISTIC_df_sub <- GISTIC_df_sub[rownames(GISTIC_df_sub) %fin% gene_list,]  # Keep only the rows with gene names in the TF list
  rownames(GISTIC_df_sub) <- stripped_gene_ids  # Add the gene symbol as row names
  return(GISTIC_df_sub)
}

#' Simple helper function that adjusts the formatting of the sample UUID and 
#' returns the UUID with appropriate formatting
#' @param samp_id a sample UUID to be adjusted
fix_samp_id <- function(samp_id) {
  if (startsWith(samp_id, "X")) {
    samp_id <- substr(samp_id, 2, nchar(samp_id)-1)
  }
  split_id <- unlist(strsplit(samp_id, split = ".", fixed = TRUE))
  new_id <- paste(split_id, collapse = "-")
  print(paste("new_id: ", new_id))
  return(new_id)
}

#' Takes a vector of ENSG IDs and strips the version from them
#' (the period and the bit after)
#' @param gene_ids vector of ENSG IDs to be stripped
strip_gene_ids <- function(gene_ids) {
  stripped_ids <- unlist(lapply(gene_ids, FUN = function(id) unlist(strsplit(id, split = ".", fixed = TRUE))[1]))
  return(stripped_ids)
}

############################################################
# RUN FILTERING FUNCTION ON CNV FILES
############################################################
# Run filtering function (TF targets)
cnv_df_iprotein <- filter_CNV(GISTIC_file, patient_id_list, ensg_ids_iprotein, biosp_df, uuid_barcode_conv) 
cnv_df_idomain <- filter_CNV(GISTIC_file, patient_id_list, ensg_ids_idomain, biosp_df, uuid_barcode_conv)
cnv_df_ibindingpos <- filter_CNV(GISTIC_file, patient_id_list, ensg_ids_ibindingpos, biosp_df, uuid_barcode_conv)

# Run filtering function (all genes)
genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg38") 
gene_ids <- genes$ensembl_gene_id
cnv_df_allgenes <- filter_CNV(GISTIC_file, patient_id_list, gene_ids, biosp_df, uuid_barcode_conv)

# The values in the dataframe themselves don't need any further processing, so
  # we're all ready to go to input this into the linear model!

# Write to CSV files
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/CNV/"
output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/CNV/"
fwrite(cnv_df_iprotein, paste(output_path, "CNV_DF_iprotein.csv", sep = ""), row.names = TRUE)
fwrite(cnv_df_idomain, paste(output_path, "CNV_DF_idomain.csv", sep = ""), row.names = TRUE)
fwrite(cnv_df_ibindingpos, paste(output_path, "CNV_DF_ibindingpos.csv", sep = ""), row.names = TRUE)
fwrite(cnv_df_allgenes, paste(output_path, "CNV_DF_AllGenes.csv", sep = ""), row.names = TRUE)

# Read back from CSV files
cnv_df_iprotein <- read.csv(paste(output_path, "CNV_DF_iprotein.csv", sep = ""), header = TRUE, row.names = 1)
cnv_df_idomain <- read.csv(paste(output_path, "CNV_DF_idomain.csv", sep = ""), header = TRUE, row.names = 1)
cnv_df_ibindingpos <- read.csv(paste(output_path, "CNV_DF_ibindingpos.csv", sep = ""), header = TRUE, row.names = 1)
cnv_df_allgenes <- read.csv(paste(output_path, "CNV_DF_AllGenes.csv", sep = ""), header = TRUE, row.names = 1)

#
#
#
#
#
#


############################################################
############################################################
# OPTION 2: USE RAW GENE-LEVEL COPY NUMBER VALUE (NUMBER
# OF AMPLIFICATIONS) FROM ASCAT2
############################################################
############################################################

############################################################
# IMPORT COPY NUMBER VALUE FILENAMES & RELATED FILES
############################################################
# Import the per-patient CNV files for particular project
ascat2_filenames <- list.files(paste(path, "CNV_Data/Gene-Level Copy Number Raw/", sep = ""), pattern = ".gene_level_copy_number")

# Import the list of patient IDs that have all four data types (the '2' indicates that they were created using the raw CNV
  # files rather than the GISTIC2 score CNV files)
patient_id_list <- read.table(paste(path, "unique_brca_patient_ids_2.txt", sep = ""), header = TRUE)[,1]
#patient_id_list <- read.table(paste(path, "unique_patient_ids_2.txt", sep = ""), header = FALSE)[,2]

# Import Biospecimen information to do the sample-patient ID conversions
biosp_df <- read.csv(paste(path, "CNV_Data/Gene-Level Copy Number Raw/aliquot.csv", sep = ""), header = TRUE, na.strings = "\'--")
# biosp_df <- read.csv(paste(path, "aliquot.csv", sep = ""), header = TRUE, na.strings = "\'--")


############################################################
# MAKE PER-GENE OUTPUT DATAFRAME
############################################################
#' Takes all the patient ASCAT2, raw CNV value filenames and, for each patient
#' that is in the patient ID list, adds the copy number value per gene, per patient 
#' to an output data frame (see top for format) which it then returns
#' @param path the local path to the directory of ascat2 files
#' @param ascat2_filenames vector of filenames in the path's directory
#' @param patient_id_list a vector of patient 4-letter/number TCGA IDs with all 
#' four data types of interest
#' @param biosp_df the biospecimen information for sample-patient ID conversions
make_output_cnv_dataframe <- function(path, ascat2_filenames, patient_id_list, 
                                      biosp_df) {

  # Loop through all patients, to get a column of CNV values for each gene for each patient
  cnv_cols <- lapply(1:length(ascat2_filenames), function(i) {
    filename <- ascat2_filenames[i]
    print(paste(i, paste("/", length(ascat2_filenames))))
    
    # Is this patient in the patient ID list?
    uuid_aliquot <- unlist(strsplit(filename, ".", fixed = TRUE))[2]
    #print(uuid_aliquot)
    samp_submitter_id <- biosp_df[biosp_df$aliquot_id == uuid_aliquot, 'sample_submitter_id']
    #tcga_id <- unique(uuid_barcode_conv[uuid_barcode_conv$uuid he== uuid_case, 'tcga_barcode'])
    patient_id <- unlist(strsplit(samp_submitter_id, "-", fixed = TRUE))[3]
    #print(patient_id)
    
    if(length(patient_id) == 0) {return(NA)}
    else if (patient_id %fin% patient_id_list) {
      # Upload this CNV file
      cnv_file <- data.table::fread(paste(path, paste("CNV_Data/Gene-Level Copy Number Raw/", filename, sep = ""),
                                          sep = ""), sep = "\t", header = TRUE)
      
      # Get the copy number info for this
      copy_num_col <- cnv_file$copy_number
      
      # Prepend the patient ID/sample id & return
      full_patient_id <- paste(unlist(strsplit(samp_submitter_id, "-", fixed = TRUE))[3:4], collapse = "-")
      
      #return(c(patient_id, copy_num_col))
      return(c(full_patient_id, copy_num_col))
    } 
    else {return(NA)}
  })
  cnv_cols <- cnv_cols[!is.na(cnv_cols)]
  
  # Combine these columns into a dataframe
    # Rows : Transcription Factors (ENSG ID) OR all genes in genome
    # Columns : Patient TCGA ID (4-digit) OPT: (-Sample ID (01A is tumor, 11A is normal)); double check that the CNA values for normal samples exist
    # Entries : CNA value 
  output_cnv_df <- as.data.frame(do.call("cbind", cnv_cols))

  # Convert the top row to column names (patient IDs)
  colnames(output_cnv_df) <- as.character(unlist(output_cnv_df[1,]))
  output_cnv_df <- output_cnv_df[-1,]
  
  # Get the gene IDs from a sample to add as rownames
  cnv_file <- data.table::fread(paste(path, paste("CNV_Data/Gene-Level Copy Number Raw/", ascat2_filenames[1], sep = ""),
                                      sep = ""), sep = "\t", header = TRUE)
  rownames(output_cnv_df) <- cnv_file$gene_id
  
  return(output_cnv_df)
}

# Run the function
output_cnv_df <- make_output_cnv_dataframe(path, ascat2_filenames, 
                                           patient_id_list, biosp_df)

# Write the output data frame to a file
fwrite(output_cnv_df, paste(output_path, "Gene-level Raw/CNV_DF_AllGenes.csv", sep = ""))


############################################################
# CREATE A BUCKETED CNA FILE
############################################################
#' Given an input CNA file with the raw copy number value for all patient samples, 
#' convert these numbers to a corresponding bucket:
#' 1. Including Amplifications
#' Bucket 1 : Both copies deleted (0 copies)
#' Bucket 2 : One copy deleted (1 copy)
#' Bucket 3 : Normal number of copies, or more (2+ copies)
#' 2. Excluding amplifications
#' Bucket 1 : Both copies deleted (0 copies)
#' Bucket 2 : One copy deleted (1 copy)
#' Bucket 3 : Normal number of copies, (2 copies); >2 copies results in NA
#' @param cna_df the CNA data frame with the raw copy numbers for all patient samples
#' @param inclAmp a TRUE/FALSE value to indicate whether or not we are including 
#' amplification events (>2 copies)
create_bucketed_cna_df <- function(cna_df, inclAmp) {
  new_df <- apply(cna_df[,1:ncol(cna_df)], MARGIN = c(1,2), function(x) {
    
    # Establish the bucket for this value and return it 
    if(is.na(x)) {return(NA)}
    if(x == 0) {return(c(1,0,0))}
    else if (x == 1) {return(c(0,1,0))}
    else if (x == 2) {return(c(0,1,0))}
    else {
      if (inclAmp) {return(c(0,0,1))}
      else {return(NA)}
    }
  })
  print(head(new_df))
  return(new_df)
}


# Call this function
bucketed_cna_df <- create_bucketed_cna_df(output_cnv_df, inclAmp = TRUE)
#bucketed_cna_df <- create_bucketed_cna_df(output_cnv_df, inclAmp = FALSE)

rownames(bucketed_cna_df) <- rownames(output_cnv_df)

# Write to file (tab-separated so they are easier to read back)
write.table(bucketed_cna_df, paste(output_path, "Gene-level Raw/CNV_DF_bucketed_inclAmp_AllGenes.tsv", sep = ""),
          sep = "\t")
#write.table(bucketed_cna_df, paste(output_path, "Gene-level Raw/CNV_DF_bucketed_exclAmp_AllGenes.tsv", sep = ""),
          #sep = "\t")


############################################################
### VISUALIZE THE FREQUENCY OF THE GENE-LEVEL CNA DATA
############################################################
#TODO: fill this part in!


