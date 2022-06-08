############################################################
### Process Gene-Level CNA Data 
### Written By: Sara Geraghty, August 2020
############################################################

# This file processes the gene-level CNA (Score or Raw Value) data file(s) in order to:
  # 1. Filter to include only patients that have all types of data
  # 2. Filter to include only the transcription factors of interest

# Format of output data frame: 
  # Rows : Transcription Factors (ENSG ID) OR all genes in genome
  # Columns : Patient ID
      # TCGA ID (4-digit) OPT: (-Sample ID (01A is tumor, 11A is normal))
      # METABRIC ID (XX-XXXX), all are primary tumor
  # Entries : CNA value (positive represents an amplification event, negative represents deletion event,
    # 0 represents no amplification or deletion) OR the actual CNA value

library(TCGAbiolinks)
library(data.table)
library(EnsDb.Hsapiens.v86)

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
write.csv(output_cnv_df, paste(output_path, "Gene-level Raw/CNV_DF_AllGenes.csv", sep = ""))


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
### CREATE A DATA FRAME WITH THE AMPLIFICATION/ DELETION 
### STATUS OF NEIGHBORING TUMOR DRIVER GENES
############################################################
#' Using the ASCAT2 CNA files, along with an expression file and 
#' a table of known oncogenes/ tumor suppressors, creates a new table that,
#' for all genes and samples, gives the CNA
#' status of the two most proximal cancer drivers genes on the same chromosome,
#' for which the CNA status results in a significant expression change 
#' across the population as a whole.
#' @param path_to_raw_cna_files a path to the raw, patient-level ASCAT2 files, 
#' which have positions of the CNAs in addition to the gene
#' @param expression_df an expression data frame with the TMM-normalized 
#' expression value for each sample (y-axis) and each gene (x-axis)
#' @param driver_gene_df a table of known driver genes with their tumor
#' suppressor/ oncogene status
#' @param sample_sheet a GDC sample sheet to link patient UUID to TCGA ID
#' @param cna_df a compiled CNA data frame from the previous section
create_flanking_cna_status_df <- function(path_to_raw_cna_files, expression_df, 
                                          driver_gene_df, sample_sheet, cna_df) {
  # Get the names of the CNA files for each patient
  cna_files <- list.files(path_to_raw_cna_files, pattern = ".gene_level_copy_number")
  
  # For each of these files, import the DF and extract the necessary info
  # to create a mini-DF for each patient
  patient_dfs <- lapply(1:length(cna_files), function(i) {
    fn <- cna_files[i]
    print(paste("File", paste(i, paste("of", length(cna_files)))))
    curr_file <- fread(paste0(path_to_raw_cna_files, fn), header = TRUE, 
                          check.names = FALSE, sep = "\t")
    curr_file$gene_id <- unlist(lapply(curr_file$gene_id, function(x)
      unlist(strsplit(x, ".", fixed = TRUE))[1]))
    curr_file <- curr_file[curr_file$gene_id %fin% rownames(expression_df),]
    
    sample_id <- unlist(strsplit(sample_sheet[sample_sheet$`File Name` == fn, 'Sample ID'], 
                                 ",", fixed = TRUE))[1]
    sample_id_short <- paste(unlist(strsplit(sample_id, "-", fixed = TRUE))[3:4], collapse = "-")
    print(sample_id_short)
    
    outdf <- data.frame("gene_id" = curr_file$gene_id, "gene_name" = curr_file$gene_name,
                        "copy_number" = curr_file$copy_number)
    #outdf$gene_id <- unlist(lapply(outdf$gene_id, function(x) unlist(strsplit(x, ".", fixed = TRUE))[1]))
    outdf$sample_id <- rep(sample_id_short, times = nrow(outdf))
    outdf <- outdf[!is.na(outdf$copy_number),]
    print(head(outdf))

    outdf <- merge(outdf, get_flanking_gene_cnas(outdf, curr_file, driver_gene_df, 
                                                 expression_df, cna_df), by = "gene_name")
    return(outdf)
  })
  
  full_df <- do.call(rbind, patient_dfs)
  return(full_df)
}


#' Helper function that, given the output DF in-progress for a given sample,
#' along with the driver gene information and an expression DF, returns a 
#' DF with: the names of each closest flanking driver gene on the same chromosome
#' with a matching, non-normal CNA status that has differential expression effects in the population
#' @param outdf the outfile for this patient, in progress
#' @param curr_file the CNA file for this patient
#' @param driver_gene_df a table of known driver genes with their tumor
#' suppressor/ oncogene status
#' @param expression_df an expression data frame with the TMM-normalized 
#' expression value for each sample (y-axis) and each gene (x-axis)
#' @param cna_df a compiled CNA data frame from the previous section
get_flanking_gene_cnas <- function(outdf, curr_file, driver_gene_df, 
                                   expression_df, cna_df) {
  # Loop through each gene in the outdf
  gene_level_dfs <- lapply(1:nrow(outdf), function(i) {
    curr_gene <- outdf$gene_name[i]
    curr_cna_stat <- outdf$copy_number[i]
    print(paste(curr_gene, curr_cna_stat))
    
    # Get the two most flanking genes on the same chr
    curr_chrom <- curr_file[curr_file$gene_name == curr_gene, 'chromosome']
    curr_start <- curr_file[curr_file$gene_name == curr_gene, 'start']
    curr_end <- curr_file[curr_file$gene_name == curr_gene, 'end']
    
    # Limit the possibilities to those genes on the same chromosome, that are
    # not the current gene
    curr_file_sub <- curr_file[(curr_file$chromosome == curr_chrom) &
                               (curr_file$gene_name != curr_gene) &
                                 (!is.na(curr_file$copy_number)),]
    
    gene_df <- data.frame("gene_name" = curr_gene, "CN" = curr_cna_stat, 
                    "gene_start" = curr_start, "gene_end" = curr_end)
    print("gene DF")
    print(gene_df)
    print(nrow(gene_df))

    # Use helper function to get the upstream and downstream neighbors
    if(!is.na(curr_cna_stat)) {
      upstr_neighbor <- get_neighbor(curr_file_sub, curr_start, curr_end, curr_cna_stat, 
                                     driver_gene_df, cna_df, expression_df, "upstream")
      print("Upstr Neighbor")
      print(upstr_neighbor)
      print("Downstr Neighbor")
      downstr_neighbor <- get_neighbor(curr_file_sub, curr_start, curr_end, curr_cna_stat,
                                       driver_gene_df, cna_df, expression_df, "downstream")
      print(downstr_neighbor)
      gene_df <- cbind(gene_df, upstr_neighbor)
      gene_df <- cbind(gene_df, downstr_neighbor)
    }
    print("gene DF new")
    print(head(gene_df))
    return(gene_df)
  })
  
  # Bind together the gene level DFs 
  full_df <- do.call(rbind, gene_level_dfs)
  print("FULL DF")
  print(head(full_df))
  return(full_df)
}


#' Helper function that, given a CNA file for a given patient that has been
#' subsetted to the relevant chromosome, looks either up or downstream to find
#' the closest driver gene neighbor. Returns the name of the neighbor.
#' @param curr_file_sub a CNA file for given patient subsetted to chromosome
#' of interest
#' @param curr_start the start position of the CNA for this gene
#' @param curr_end the end position of the CNA for this gene
#' @param curr_cna_stat the raw CNA value for this gene
#' @param driver_gene_df a table of known driver gene with their tumor
#' suppressor/ oncogene status
#' @param expression_df an expression data frame with the TMM-normalized 
#' expression value for each sample (y-axis) and each gene (x-axis)
#' @param cna_df a compiled CNA data frame from the previous section
#' @param up_or_down whether we are looking "upstream" or "downstream" 
get_neighbor <- function(curr_file_sub, curr_start, curr_end, curr_cna_stat,
                         driver_gene_df, cna_df, expression_df, up_or_down) {
  
  driver_stat <- FALSE  # loop until we have found a neighboring driver gene or
  matching_cna <- TRUE  # loop until we no longer have a matching CNA value (i.e. the end of the amplif. or del.)
  neighbor_df <- data.frame(matrix(nrow = 1, ncol = 3))
  
  if(length(curr_cna_stat) == 0) {return(neighbor_df)}
  
  if(up_or_down == "upstream") {
    neighbor_df <- data.frame("upstr_neighbor" = NA, "upstr_neighbor_start" = NA, 
                              "upstr_neighbor_end" = NA)
    if(curr_cna_stat == 2) {return(neighbor_df)}
    tryCatch({
      curr_file_sub <- curr_file_sub[(curr_file_sub$start < curr_start) &
                                       (curr_file_sub$end <= curr_end)]
    }, error = function(cond) {return(neighbor_df)})
  }
  else if (up_or_down == "downstream") {
    neighbor_df <- data.frame("downstr_neighbor" = NA, "downstr_neighbor_start" = NA, 
                              "downstr_neighbor_end" = NA)
    if(curr_cna_stat == 2) {return(neighbor_df)}
    tryCatch({
      curr_file_sub <- curr_file_sub[(curr_file_sub$start >= curr_start) &
                                       (curr_file_sub$end > curr_end)]
    }, error = function(cond) {return(neighbor_df)})
  } else {print("Must be either 'upstream' or 'downstream'.")}
  
  # Search up the chromosome until we find a driver gene with a matching CNA (if it exists)
  while ((driver_stat == FALSE) & (matching_cna == TRUE)) {
    nearest <- NA
    ind_nearest <- NA
    
    if(up_or_down == "upstream") {
      ind_nearest <- which.min(curr_file_sub$end - curr_start)
      nearest <- curr_file_sub[ind_nearest,]
      
    } else if (up_or_down == "downstream") {
      ind_nearest <- which.min(curr_end - curr_file_sub$start)
      nearest <- curr_file_sub[ind_nearest,]
      
    } else {print("Must input either 'upstream' or 'downstream'.")}
    
    #print("nearest")
    #print(nearest)
    #print(paste("ind nearest", ind_nearest))
    #print(paste("curr cna stat", curr_cna_stat))
    
    if (!(length(nearest) == 0)) {
      if ((!(nrow(nearest) == 0)) & (!(ncol(nearest) < 8))) {
        neighbor_df[,1] <- nearest$gene_name
        neighbor_df[,2] <- nearest$start
        neighbor_df[,3] <- nearest$end
      }
    }
    print("Neighbor DF")
    print(neighbor_df)
    
    # As long as they have the same copy number, check if this is a driver gene
    if (!(length(nearest) == 0)) {
      if ((!(nrow(nearest) == 0)) & (!(ncol(nearest) < 8)) & 
          ("copy_number" %fin% colnames(nearest))) {
        if (as.integer(nearest$copy_number) == curr_cna_stat) {
          if (nearest$gene_name %fin% driver_gene_df$primary_gene_names) {
            # Make sure that, for an amplification, this is an oncogene and, for a 
            # deletion, this is a tumor suppressor
            status <- driver_gene_df[driver_gene_df$primary_gene_names == nearest$gene_name, 
                                     'cancer_driver_status']
            status <- unlist(strsplit(status, ",", fixed = TRUE))
            
            if (((curr_cna_stat > 2) & ("O" %fin% status)) | ((curr_cna_stat < 2) & ("T" %fin% status)) |
                (length(intersect(c("O", "T"), status)) == 0)) {
              
              amp_or_del <- NA
              if (curr_cna_stat > 2) {amp_or_del <- "amp"}
              else if (curr_cna_stat < 2) {amp_or_del <- "del"}
              else {amp_or_del <- "norm"}
              
              # If we are looking at an amplification or deletion event, then make sure there
              # is an expression difference in the right direction
              if(!(amp_or_del == "norm")) {
                exp_change <- get_expr_change(nearest, expression_df, cna_df, amp_or_del)
                if (exp_change) {
                  driver_stat <- TRUE
                } else {
                  curr_file_sub <- curr_file_sub[-ind_nearest,]
                }
              } else {driver_stat <- TRUE}
            } else {
              curr_file_sub <- curr_file_sub[-ind_nearest,]
            } 
          } else {
            curr_file_sub <- curr_file_sub[-ind_nearest,]
          } 
        } else {
          matching_cna <- FALSE
          neighbor_df[1,] <- rep(NA, times = 3)
        }
      } else {
        matching_cna <- FALSE
        neighbor_df[1,] <- rep(NA, times = 3)
      }
    } else {
      matching_cna <- FALSE
      neighbor_df[1,] <- rep(NA, times = 3)
    }
  }
  
  #print(paste("Returning Neighbor"))
  #print(neighbor_df)
  #print(nrow(neighbor_df))
  return(neighbor_df)
}

#' Given a neighbor gene and an expression DF, checks to see if there is a 
#' significant difference in expression in the proper direction when this gene
#' is amplified/ deleted. Returns TRUE if so, FALSE otherwise.
#' @param nearest a single-row table with the name, CNA, etc. of the neighbor gene
#' @param expression_df a TMM-normalized expression DF
#' @param cna_df a CNA DF produced from the previous section
#' @param amp_or_del whether we expect "amp", "del"
get_expr_change <- function(nearest, expression_df, cna_df, amp_or_del) {
  gene_id <- nearest$gene_id
  cna_df_sub <- cna_df[rownames(cna_df) == gene_id,]
  expression_df_sub <- expression_df[rownames(expression_df) == gene_id,]
  
  # For amplification, split the patient population into those with amplifications 
  # and those without
  if(amp_or_del == "amp") {
    ampl_patients <- colnames(cna_df_sub)[which(unlist(cna_df_sub[1,]) > 2)]
    nonampl_patients <- setdiff(colnames(cna_df_sub), ampl_patients)
    
    express_ampl_patients <- as.numeric(expression_df_sub[,colnames(expression_df_sub) %fin% ampl_patients])
    express_nonampl_patients <- as.numeric(expression_df_sub[,colnames(expression_df_sub) %fin% nonampl_patients])
    
    wilcox_res <- wilcox.test(express_ampl_patients, express_nonampl_patients, alternative = "greater")
  } else {
    del_patients <- colnames(cna_df_sub)[which(unlist(cna_df_sub[1,]) < 2)]
    nondel_patients <- setdiff(colnames(cna_df_sub), del_patients)
    
    express_del_patients <- as.numeric(expression_df_sub[,colnames(expression_df_sub) %fin% del_patients])
    express_nondel_patients <- as.numeric(expression_df_sub[,colnames(expression_df_sub) %fin% nondel_patients])
    
    wilcox_res <- wilcox.test(express_del_patients, express_nondel_patients, alternative = "less")
  }
  
  print(wilcox_res)
  if(wilcox_res$p.value < 0.05) {return(TRUE)}
  else(return(FALSE))
}


# Import necessary files
path_to_cna_files <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/CNV_Data/Gene-Level Copy Number Raw/"
sample_sheet <- read.csv(paste0(path_to_cna_files, "gdc_sample_sheet.2021-02-17.tsv"), header = TRUE, 
                         check.names = FALSE, sep = "\t")

local_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
cna_df <- read.csv(paste0(local_path, "CNV/Gene-level Raw/CNV_DF_AllGenes_CancerOnly.csv"), 
                   header = TRUE, check.names = FALSE, row.names = 1)
expression_df <- read.csv(paste0(local_path, "Expression/expression_tmm_DF_filtByExpr_2groups_SDGr1_TO.csv"), 
                          header = TRUE, check.names = FALSE, row.names = 1)
colnames(expression_df) <- unlist(lapply(colnames(expression_df), function(x) 
  paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))

driver_gene_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/GRCh38_driver_gene_list.tsv",
                           header = TRUE, comment.char = "#", sep = "\t")

# Call function
flanking_cna_status_df <- create_flanking_cna_status_df(path_to_cna_files, expression_df, driver_gene_df,
                                                        sample_sheet, cna_df)

# If desired, add a column for the ENSG ID as well
flanking_cna_status_df$ensg_id <- unlist(lapply(flanking_cna_status_df$gene_name, function(x) 
  paste(unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'ensembl_gene_id']), collapse = ";")))

# Write this to a file
write.csv(flanking_cna_status_df, paste0(local_path, "CNV/Gene-level Raw/Neighboring_Driver_CNA_DF_AllGenes_CancerOnly.csv"))


############################################################
### CREATE A DATA FRAME WITH BINARY AMPLIFICATION/ DELETION 
### STATUS OF NEIGHBORING TUMOR DRIVER GENES
############################################################
#' From this, create a binary sample x gene matrix that, as the value,
#' has two binary values separated by a semicolon (e.g. 0;0, 0;1, 1;0,
#' or 1;1) to indicate whether, if the gene in question is amplified or
#' deleted in the given sample, there is an upstream;downstream cancer
#' driver gene on the same chromosome that is likewise amplified or 
#' deleted and has a resultant expression change
#' @param flanking_cna_status_df a CNA data frame created from
#' the above functions that has details of neighboring cancer driver
#' gene amplification and deletion events
create_binary_flanking_cna_status_df <- function(flanking_cna_status_df) {
  samples <- unique(flanking_cna_status_df$sample_id)
  genes <- unique(flanking_cna_status_df$ensg_id)
  
  # Get all the non-NA entries from the flanking CNA status DF and use them
  # to create binary entries for new DF
  new_dfs <- lapply(1:length(samples), function(i) {
    sample <- samples[i]
    
    flanking_df_sub <- flanking_cna_status_df[flanking_cna_status_df$sample_id == sample,]
    
    vals <- unlist(lapply(1:nrow(flanking_df_sub), function(j) {
      row <- flanking_df_sub[j,]
      if (!is.na(row$upstr_neighbor)) {
        (!is.na(row$downstr_neighbor)) {
          val <- "1;1"
        } else {
          val <- "1;0"
        }
      } else if (!is.na(row$downstr_neighbor)) {val <- "0;1"}
      else {val <- "0;0"}
      return(val)
    }))
    tmp_df <- data.frame("gene_id" = flanking_df_sub$gene_id, "value" = vals)
    return(tmp_df)
  })
  
  combined_df <- Reduce(function(x,y) merge(x = x, y = y, by = "gene_id"), new_dfs)
  colnames(combined_df) <- samples
  rownames(combined_df) <- combined_df$gene_id
  combined_df <- combined_df[, -c("gene_id")]

  return(combined_df)
}

# Call function
binary_flanking_cna_status_df <- create_binary_flanking_cna_status_df(flanking_cna_status_df)

# Write to a file
write.csv(binary_flanking_cna_status_df, paste0(local_path, "CNV/Gene-level Raw/Binary_Neighboring_Driver_CNA_DF_AllGenes_CancerOnly.csv"))


################################


#' As an alternative, creates a data frame containing all the CGC oncogenes and tumor 
#' suppressors with their chromosome, start and end position, whether they are differently 
#' expressed between amplified/ deleted patient cohorts within the larger group, and, for 
#' each patient, their copy number. Then, we can reference this mapping for each gene
#' of interest in order to obtain if there are matching CGC neighbors that are co-amplified or
#' deleted.
#' @param expression_df an expression data frame with the TMM-normalized 
#' expression value for each sample (y-axis) and each gene (x-axis)
#' @param driver_gene_df a table of known driver genes with their tumor
#' suppressor/ oncogene status
#' @param cna_df a compiled CNA data frame from the previous section
create_cgc_cna_mapping <- function(expression_df, driver_gene_df, cna_df) {
  driver_gene_mapping <- data.frame("driver.gene" = driver_gene_df$ensembl_gene_id, 
                                    "onco.or.ts" = driver_gene_df$cancer_driver_status)
  driver_gene_mapping$onco.or.ts <- unlist(lapply(driver_gene_mapping$onco.or.ts, function(x) {
    x_spl <- unlist(strsplit(x, ",", fixed = TRUE))
    x_spl <- x_spl[x_spl %in% c("O", "T")]  # keep only T or O labels
  }))
  
  # Get the chromosome, start position, and strand for this gene
  edb <- EnsDb.Hsapiens.v86
  driver_gene_symbols <- driver_gene_df$primary_gene_name
  driver_gene_mapping$chrom <- mapIds(edb, keys = driver_gene_symbols, 
                             keytype = "SYMBOL", column = "SEQNAME")
  driver_gene_mapping$txStart <- mapIds(edb, keys = driver_gene_symbols, 
                               keytype = "SYMBOL", column = "TXSEQSTART")
  driver_gene_mapping$txStrand <- mapIds(edb, keys = driver_gene_symbols, 
                                keytype = "SYMBOL", column = "SEQSTRAND")
  driver_gene_mapping <- driver_gene_mapping[!(is.na(driver_gene_mapping$chrom))]
  
  # Get the copy number for each of these genes in each patient
  cna_df_drivers <- cna_df[rownames(cna_df) %in% driver_gene_mapping$driver.gene,]
  cna_df_drivers$driver.gene <- rownames(cna_df_drivers)
  rownames(cna_df_drivers) <- 1:nrow(cna_df_drivers)
  
  driver_gene_mapping_full <- merge(driver_gene_mapping, cna_df_drivers, by = "driver.gene")
  
  return(driver_gene_mapping_full)
}




