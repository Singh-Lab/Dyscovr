############################################################
### Process Gene-Level CNA Data 
### Written By: Sara Camilli, August 2020
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
library(Rfast)
library(parallel)
library(dplyr)

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/"
#path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/"

output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/CNV/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/CNV/"


# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


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
#patient_id_list <- read.table(paste(path, "unique_patient_ids_2.txt", sep = ""), header = T)[,1]

# Import Biospecimen information to do the sample-patient ID conversions
biosp_df <- read.csv(paste(path, "CNV_Data/Gene-Level Copy Number Raw/aliquot.csv", sep = ""), header = TRUE, na.strings = "\'--")
# biosp_df <- read.csv(paste(path, "aliquot.csv", sep = ""), header = TRUE, na.strings = "\'--")

# All genes we are interested in
allgene_targets <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/allgene_targets.csv", 
                            header = TRUE, row.names = 1, check.names = FALSE)

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
#' @param allgene_targets the genes we are specifically interested in profiling
make_output_cnv_dataframe <- function(path, ascat2_filenames, patient_id_list, 
                                      biosp_df, allgene_targets) {

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

    if(length(patient_id) == 0) {return(NA)}
    else if (patient_id %fin% patient_id_list) {
      # Upload this CNV file
      cnv_file <- data.table::fread(paste(path, paste("CNV_Data/Gene-Level Copy Number Raw/", filename, sep = ""),
                                          sep = ""), sep = "\t", header = TRUE)
      gene_ids <- unlist(lapply(cnv_file$gene_id, function(x) unlist(strsplit(x, ".", fixed = T))[1]))
      cnv_file <- cnv_file[which(gene_ids %fin% allgene_targets$ensg), ]
      
      # Get the copy number info for this
      copy_num_col <- cnv_file$copy_number
      
      # Prepend the patient ID/sample id & return
      full_patient_id <- paste(unlist(strsplit(samp_submitter_id, "-", fixed = TRUE))[3:4], collapse = "-")
      
      #return(c(patient_id, copy_num_col))
      copy_num_df <- as.data.frame(copy_num_col)
      colnames(copy_num_df) <- full_patient_id
      rownames(copy_num_df) <- cnv_file$gene_id

      return(copy_num_df)
    } 
    else {return(NA)}
  })
  cnv_cols <- cnv_cols[!is.na(cnv_cols)]
  
  # Combine these columns into a dataframe
    # Rows : Transcription Factors (ENSG ID) OR all genes in genome
    # Columns : Patient TCGA ID (4-digit) OPT: (-Sample ID (01A is tumor, 11A is normal)); double check that the CNA values for normal samples exist
    # Entries : CNA value 
  #output_cnv_df <- as.data.frame(do.call("cbind", cnv_cols))
  
  #output_cnv_df <- as.data.frame(Reduce(function(x, y) dplyr::full_join(x, y, by = "ensg_id"), cnv_cols))

  # Convert the top row to column names (patient IDs)
  #colnames(output_cnv_df) <- as.character(unlist(output_cnv_df[1,]))
  #output_cnv_df <- output_cnv_df[-1,]

  #return(output_cnv_df)
  return(cnv_cols)
}

# Run the function
output_cnv_df_cols <- make_output_cnv_dataframe(path, ascat2_filenames, 
                                           patient_id_list, biosp_df,
                                           allgene_targets)

if(length(unique(unlist(lapply(output_cnv_df_cols, function(df) nrow(df))))) > 1) {
  lengths <- unique(unlist(lapply(output_cnv_df_cols, function(df) nrow(df))))
  output_cnv_dfs_perLength <- lapply(lengths, function(l) {
    to_keep_l <- unlist(lapply(1:length(output_cnv_df_cols), function(i) 
      if(nrow(output_cnv_df_cols[[i]]) == l) {return(i)}))
    output_cnv_df_l <- output_cnv_df_cols[to_keep_l]
    output_cnv_df_l_full <- do.call(cbind, output_cnv_df_l)
    output_cnv_df_l_full$ensg_id <- rownames(output_cnv_df_l_full)
    return(output_cnv_df_l_full)
  })
  output_cnv_df_full <- merge(output_cnv_dfs_perLength[[1]],  output_cnv_dfs_perLength[[2]], 
                              by = 'ensg_id', all = T)   #TODO: this does not combine multiple gene isoforms, which results in many NAs
  rownames(output_cnv_df_full) <- output_cnv_df_full$ensg_id
  output_cnv_df_full <- output_cnv_df_full[,!(colnames(output_cnv_df_full) %in% 'ensg_id')]
  
} else {
  output_cnv_df_full <- do.call(cbind, output_cnv_df_cols)
  rownames(output_cnv_df_full) <- output_cnv_df_full$ensg_id
  output_cnv_df_full <- output_cnv_df_full[,!(colnames(output_cnv_df_full) %in% 'ensg_id')]
}
# Get the gene IDs from a sample to add as rownames
#cnv_file <- data.table::fread(paste(path, paste("CNV_Data/Gene-Level Copy Number Raw/", ascat2_filenames[1], sep = ""),
 #                                   sep = ""), sep = "\t", header = TRUE)
#rownames(output_cnv_df) <- cnv_file$gene_id

# Write the output data frame to a file
write.csv(output_cnv_df_full, paste(output_path, "Gene-level Raw/CNV_DF_AllGenes.csv", sep = ""))


# Check the number of patients that have an ASCAT file
ascat2_filename_uuids <- unlist(lapply(ascat2_filenames, function(f) unlist(strsplit(f, ".", fixed = T))[2]))
ascat2_filename_case_ids <- unlist(lapply(ascat2_filename_uuids, function(uuid) 
  biosp_df[biosp_df$aliquot_id == uuid, 'case_submitter_id']))
ascat2_filename_patient_ids <- unlist(lapply(ascat2_filename_case_ids, function(id) unlist(strsplit(id, "-", fixed = TRUE))[3]))
length(unique(ascat2_filename_patient_ids))  # original 8297 files; now 9924 files, pan-cancer; 1054 files BRCA

# what is the intersection with our patient ID list?  # 6123 PC; 732 BRCA
length(unique(intersect(ascat2_filename_patient_ids, patient_id_list)))

# Are there CNA files we should have that appear to be missing? If so, use these IDs to make sure all files are downloaded
length(setdiff(patient_id_list, ascat2_filename_patient_ids))
case_ids_missing_cnv <- unlist(lapply(setdiff(patient_id_list, ascat2_filename_patient_ids), function(x) 
  unique(cnv_clin_df_sub[grepl(x, cnv_clin_df_sub$case_submitter_id), 'case_submitter_id'])))
print(paste(case_ids_missing_cnv, collapse = ","))

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
#' For every human gene, get the two most proximal cancer driver genes on the same
#' chromosome (return a table with rownames as the human genes, column one as the most
#' left-flanking cancer-driver gene, and column two as the most right-flanking cancer 
#' driver gene)
#' @param all_genes_id_conv a bioMart file with positional info for all genes of interest
#' @param expression_df an expression DF with all genes of interest
#' @param driver_gene_df a DF with driver gene information
create_proximal_driver_table <- function(all_genes_id_conv, expression_df, driver_gene_df, check_for_expr) {
  genes <- rownames(expression_df)
  
  neighbor_dfs <- lapply(1:length(genes), function(i) {
    g <- genes[i]
    print(paste("Gene", paste(i, paste("/", length(genes)))))
    
    all_genes_id_conv_sub <- all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == g, ]
    
    g_start <- all_genes_id_conv_sub$start_position
    g_end <- all_genes_id_conv_sub$end_position
    g_chrom <- all_genes_id_conv_sub$chromosome_name
    g_strand <- all_genes_id_conv_sub$strand

    all_genes_id_conv_chrom <- all_genes_id_conv[(all_genes_id_conv$chromosome_name == g_chrom) & 
                                                   (all_genes_id_conv$strand == g_strand),]
    upstr_neighbors <- all_genes_id_conv_chrom[all_genes_id_conv_chrom$end_position <= g_end, ]
    downstr_neighbors <- all_genes_id_conv_chrom[all_genes_id_conv_chrom$start_position >= g_start,]
    
    #' Helper function that gets the closest neighboring up- or downstream
    #' driver gene to the given gene
    #' @param neighbors a table with only the neighbors of the given target
    #' gene in a particular direction
    #' @param driver_gene_df a driver gene DF
    #' @param direction "up" or "down" to indicate the direction of search
    get_neighbor <- function(neighbors, driver_gene_df, direction) {
      
      if(nrow(neighbors) == 0) {return(NA)}
      
      if(direction == "up") {neighbors <- neighbors[order(neighbors$start_position, decreasing = TRUE),]}
      else {neighbors <- neighbors[order(neighbors$end_position, decreasing = FALSE),]}
      
      # Get the nearest driving neighbor
      for(i in 1:nrow(neighbors)) {
        nearest_neighbor <- as.character(neighbors[i, 'ensembl_gene_id'])
        if(nearest_neighbor %fin% driver_gene_df$ensembl_gene_id) {
          return(nearest_neighbor)
        }
      }
      return(NA)
    }
    
    upstr_neighbor <- get_neighbor(upstr_neighbors, driver_gene_df, "up")
    downstr_neighbor <-  get_neighbor(downstr_neighbors, driver_gene_df, "down")
    
    neighbor_df <- data.frame("Upstream_Neighbor" = upstr_neighbor, 
                              "Downstream_Neighbor" = downstr_neighbor)
    rownames(neighbor_df) <- g
    return(neighbor_df)
  })
  
  full_neighbor_df <- do.call(rbind, neighbor_dfs)
  full_neighbor_df <- na.omit(full_neighbor_df)
  return(full_neighbor_df)
}



#dummy_cna_file <- read.csv(paste0(path_to_cna_files, "TCGA-BRCA.0c86eab9-cfc8-4a03-823b-4331c43d2038.gene_level_copy_number.tsv"),
                          # header = TRUE, check.names = FALSE, row.names = 1, sep = "\t")

all_genes_id_conv_sub <- as.data.table(all_genes_id_conv[,c("ensembl_gene_id", "chromosome_name", 
                                                            "start_position", "end_position", "strand")])
all_genes_id_conv_sub <- distinct(all_genes_id_conv_sub)
all_genes_id_conv_sub <- all_genes_id_conv_sub[all_genes_id_conv_sub$ensembl_gene_id %fin% rownames(expression_df),]

proximal_driver_gene_table <- create_proximal_driver_table(all_genes_id_conv_sub, expression_df, driver_gene_df)

proximal_driver_gene_table <- fread(paste0(output_path, "upstream_downstream_driver_neighbors.csv"), 
                                       header = TRUE, check.names = FALSE)
colnames(proximal_driver_gene_table)[1] <- "Gene_ID"



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
#' @param proximal_driver_gene_df a data frame from the create_proximal_driver_table
#' function (neighboring upstream and downstream drivers)
#' @param sample_sheet a GDC sample sheet to link patient UUID to TCGA ID
#' @param cna_df a compiled CNA data frame from the previous section
create_flanking_cna_status_df <- function(path_to_raw_cna_files, expression_df, 
                                          proximal_driver_gene_df, sample_sheet, cna_df) {
  # Get the names of the CNA files for each patient
  cna_files <- list.files(path_to_raw_cna_files, pattern = ".gene_level_copy_number")
  outpath <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/CNV/Patient-Specific Flanking Driver CNA DFs/"
  
  
  # For each of these files, import the DF and extract the necessary info
  # to create a mini-DF for each patient
  patient_dfs <- mclapply(1:length(cna_files), function(i) {
    fn <- cna_files[i]
    curr_file <- fread(paste0(path_to_raw_cna_files, fn), header = TRUE, sep = "\t")
    curr_file$gene_id <- unlist(lapply(as.character(curr_file$gene_id), function(x)
      unlist(strsplit(x, ".", fixed = TRUE))[1]))
    curr_file <- curr_file[curr_file$gene_id %fin% expression_df$Gene_ID,]
    curr_file <- na.omit(curr_file)

    sample_id <- unlist(strsplit(as.character(sample_sheet[sample_sheet$`File Name` == fn, 'Sample ID']), 
                                 ",", fixed = TRUE))[1]
    sample_id_short <- unlist(strsplit(sample_id, "-", fixed = TRUE))[3]
    print(sample_id_short)
    
    outdf <- data.frame("gene_id" = curr_file$gene_id, "gene_name" = curr_file$gene_name,
                        "copy_number" = curr_file$copy_number)
    outdf$gene_id <- unlist(lapply(curr_file$gene_id, function(x) 
      unlist(strsplit(x, ".", fixed = TRUE))[1]))
    outdf$sample_id <- rep(sample_id_short, times = nrow(outdf))
    outdf <- outdf[!is.na(outdf$copy_number),]
    
    outdf <- merge(outdf, get_flanking_gene_cnas(outdf, curr_file, proximal_driver_gene_df, 
                                                 expression_df, cna_df), by = c("gene_id", "copy_number"))
    
    print(head(outdf))
    
    fwrite(outdf, paste0(outpath, paste0(sample_id_short, "_flanking_driver_df.csv")))
    return(outdf)
  })
  
  full_df <- do.call(rbind, patient_dfs)
  return(full_df)
}

#' Helper function that, given the output DF in-progress for a given sample,
#' along with the driver gene information and an expression DF, returns a 
#' DF with: the names of each closest flanking driver gene on the same chromosome
#' with a matching CNA status that has differential expression effects in the population
#' @param outdf the outfile for this patient, in progress
#' @param curr_file the CNA file for this patient
#' @param proximal_driver_gene_df a data frame from the create_proximal_driver_table
#' function (neighboring upstream and downstream drivers)
#' @param expression_df an expression data frame with the TMM-normalized 
#' expression value for each sample (y-axis) and each gene (x-axis)
#' @param cna_df a compiled CNA data frame from the previous section
get_flanking_gene_cnas <- function(outdf, curr_file, proximal_driver_gene_df, 
                                   expression_df, cna_df) {
  
  outdf <- outdf[outdf$gene_id %fin% proximal_driver_gene_df$Gene_ID,]
  
  median_chrom_copy_number_dict <- sapply(unique(curr_file$chromosome), FUN = function(c) {
    median <- Rfast::med(curr_file[curr_file$chromosome == c, 'copy_number'][[1]])
    if(median < 2) {median <- ceiling(median)}
    if(median > 2) {median <- floor(median)}
    return(median)
  }, simplify = FALSE, USE.NAMES = TRUE)

  # Loop through each gene in the outdf
  gene_level_dfs <- mclapply(1:nrow(outdf), function(i) {
    curr_id <- outdf$gene_id[i]
    curr_cna_stat <- outdf$copy_number[i]

    #print(paste("gene", paste(i, paste("/", nrow(outdf)))))
    
    # Get the copy number information for this gene in this sample
    gene_df <- curr_file[curr_file$gene_id == curr_id, c('gene_id', 'chromosome', 'start', 'end')]
    gene_df$copy_number <- curr_cna_stat
    
    # For the current chromosome, get the median copy number (so we know what 
    # "normal" is and can judge if there was an amplification or deletion event)
    median_chrom_copy_number <- NA
    if(!is.na(gene_df$chromosome)) {
      median_chrom_copy_number <- median_chrom_copy_number_dict[[gene_df$chromosome]]
    }
    
    #print("gene DF")
    #print(gene_df)
    #print(nrow(gene_df))
    
    # Use helper function to get the upstream and downstream neighbors, if there is 
    # an amplification or deletion affecting this gene
    if((!is.na(curr_cna_stat)) & (median_chrom_copy_number != curr_cna_stat)) {
      
      proximal_driver_gene_df_sub <- proximal_driver_gene_df[proximal_driver_gene_df$Gene_ID == curr_id, ]
      upstr_neighbor <- as.character(proximal_driver_gene_df_sub[,'Upstream_Neighbor'])
      downstr_neighbor <- as.character(proximal_driver_gene_df_sub[,'Downstream_Neighbor'])
      if(length(upstr_neighbor) == 0) {upstr_neighbor <- NA}
      if(length(downstr_neighbor) == 0) {downstr_neighbor <- NA}
      
      # Get the copy numbers of these neighbors
      upstr_copy_number <- NA
      downstr_copy_number <- NA

      # Get the copy number status of the neighbors
      exp_diff_res <- FALSE
      if(!is.na(upstr_neighbor)) {
        upstr_copy_number <- as.integer(curr_file[curr_file$gene_id == upstr_neighbor, 
                                                  'copy_number'][[1]])
        if(length(upstr_copy_number) == 0) {upstr_copy_number <- NA}
        else {
          if(upstr_copy_number == curr_cna_stat) {
            amp_or_del <- "del"
            if(curr_cna_stat > median_chrom_copy_number) {amp_or_del <- "amp"}
            exp_diff_res <- get_expr_change(upstr_neighbor, expression_df, 
                                            cna_df, amp_or_del, median_chrom_copy_number)
          }
        }
      }
      if(!is.na(downstr_neighbor)) {
        downstr_copy_number <- as.integer(curr_file[curr_file$gene_id == downstr_neighbor, 
                                                    'copy_number'][[1]])
        if(length(downstr_copy_number) == 0) {downstr_copy_number <- NA}
        else {
          if(downstr_copy_number == curr_cna_stat) {
            amp_or_del <- "del"
            if(curr_cna_stat > median_chrom_copy_number) {amp_or_del <- "amp"}
            exp_diff_res <- get_expr_change(downstr_neighbor, expression_df, 
                                            cna_df, amp_or_del, median_chrom_copy_number)
          }
        }
      }
      #print(upstr_neighbor)
      #print(downstr_neighbor)
      #print(upstr_copy_number)
      #print(downstr_copy_number)
      if (exp_diff_res) {
        gene_df_new <- bind_cols(gene_df, list("upstr_neighbor" = upstr_neighbor, 
                                               "upstr_copy_number" = upstr_copy_number, 
                                               "downstr_neighbor" = downstr_neighbor, 
                                               "downstr_copy_number" = downstr_copy_number))
        #colnames(gene_df_new)[ncol(gene_df)+1:ncol(gene_df_new)] <- c("upstr_neighbor", "upstr_copy_number", 
        #"downstr_neighbor", "downstr_copy_number")
      } else {
        gene_df_new <- bind_cols(gene_df, list("upstr_neighbor" = "No.Exp.Diff", "upstr_copy_number" = NA,
                                               "downstr_neighbor" = "No.Exp.Diff", "downstr_copy_number" = NA))
      }
      
    } else {
      gene_df_new <- bind_cols(gene_df, list("upstr_neighbor" = "No.CNA", "upstr_copy_number" = NA,
                                             "downstr_neighbor" = "No.CNA", "downstr_copy_number" = NA))
    }
    #print("gene DF new")
    #print(head(gene_df_new))
    #return(gene_df_new)
  })
  # Bind together the gene level DFs 
  full_df <- do.call(rbind, gene_level_dfs)
  print("FULL DF")
  print(head(full_df))
  return(full_df)
}


#' Given a neighbor gene and an expression DF, checks to see if there is a 
#' significant difference in expression in the proper direction when this gene
#' is amplified/ deleted. Returns TRUE if so, FALSE otherwise.
#' @param nearest a gene ID of a neighboring driver gene
#' @param expression_df a TMM-normalized expression DF
#' @param cna_df a CNA DF produced from the previous section
#' @param amp_or_del whether we expect "amp", "del"
#' @param med_cna the median CN for the given chromosome
get_expr_change <- function(gene_id, expression_df, cna_df, amp_or_del, med_cna) {
  #gene_id <- nearest$gene_id
  cna_df_sub <- cna_df[cna_df$Gene_ID == gene_id,]
  expression_df_sub <- expression_df[expression_df$Gene_ID == gene_id,]
  
  # For amplification, split the patient population into those with amplifications 
  # and those without
  if(amp_or_del == "amp") {
    ampl_patients <- colnames(cna_df_sub)[which(unlist(cna_df_sub[1,]) > med_cna)]
    ampl_patients <- ampl_patients[ampl_patients != "Gene_ID"]
    nonampl_patients <- setdiff(colnames(cna_df_sub), ampl_patients)

    express_ampl_patients <- mean(as.numeric(expression_df_sub[,colnames(expression_df_sub) %fin% ampl_patients]))
    express_nonampl_patients <- mean(as.numeric(expression_df_sub[,colnames(expression_df_sub) %fin% nonampl_patients]))
    
    if(express_ampl_patients >= express_nonampl_patients) {return(TRUE)}
    else{return(FALSE)}
    #wilcox_res <- wilcox.test(express_ampl_patients, express_nonampl_patients, alternative = "greater")
    
  } else {
    del_patients <- colnames(cna_df_sub)[which(unlist(cna_df_sub[1,]) < med_cna)]
    del_patients <- del_patients[del_patients != "Gene_ID"]
    nondel_patients <- setdiff(colnames(cna_df_sub), del_patients)
    
    express_del_patients <- mean(as.numeric(expression_df_sub[,colnames(expression_df_sub) %fin% del_patients]))
    express_nondel_patients <- mean(as.numeric(expression_df_sub[,colnames(expression_df_sub) %fin% nondel_patients]))
    
    if(express_del_patients <= express_nondel_patients) {return(TRUE)}
    else{return(FALSE)}
    #wilcox_res <- wilcox.test(express_del_patients, express_nondel_patients, alternative = "less")
  }
  
  #print(wilcox_res)
  #if(wilcox_res$p.value < 0.05) {return(TRUE)}
  #else(return(FALSE))
  return(NA)
}




# Import necessary files
path_to_cna_files <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/CNV_Data/Gene-Level Copy Number Raw/"
sample_sheet <- fread(paste0(path_to_cna_files, "gdc_sample_sheet.2021-02-17.tsv"), 
                      header = TRUE, sep = "\t")

local_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
driver_gene_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/GRCh38_driver_gene_list.tsv",
                           header = TRUE, comment.char = "#", sep = "\t")

expression_df <- fread(paste0(local_path, "Expression/expression_tmm_DF_filtByExpr_2groups_SDGr1_TO.csv"), 
                       header = TRUE)
colnames(expression_df) <- c("Gene_ID", unlist(lapply(colnames(expression_df)[2:ncol(expression_df)], function(x) 
  paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-"))))

cna_df <- fread(paste0(local_path, "CNV/Gene-level Raw/CNV_DF_AllGenes_CancerOnly.csv"), 
                header = TRUE)
colnames(cna_df)[1] <- "Gene_ID"
cna_df <- cna_df[cna_df$Gene_ID %fin% expression_df$Gene_ID,]


# Call function
cna_neighbor_df <- create_flanking_cna_status_df(path_to_cna_files, expression_df, 
                                                 proximal_driver_gene_table, sample_sheet, cna_df)


# Write this to a file
fwrite(cna_neighbor_df, paste0(local_path, "CNV/neighboring_driver_cna_df_allgenes.csv"))

# Read back 
cna_neighbor_df <- fread(paste0(local_path, "CNV/neighboring_driver_cna_df_allgenes.csv"), header = TRUE, check.names = FALSE)


#' Given a flanking status CNA data frame (produced from functions above),
#' create a 0/1 bucketed data frame with genes as rows and patients as 
#' columns, denoting whether or not there is a neighboring driver with
#' a) shared amplification or deletion status (defined as a CN value that
#' is above or below the median CN for the chromosome), and b) an expression
#' difference in the affected driver gene that matches the direction of the
#' effect (upregulated if amplified, downregulated if deleted).
#' @param cna_neighbor_df a flanking status CNA data frame as described above
bucket_flanking_cna_status_df <- function(cna_neighbor_df) {
  
  cna_neighbor_df_sub <- cna_neighbor_df[!(is.na(cna_neighbor_df$upstr_copy_number) & 
                                             is.na(cna_neighbor_df$downstr_copy_number)),]
  
  # Initialize a new data frame that we will fill
  bucketed_cna_neighbor_df <- data.frame(matrix(ncol = length(unique(cna_neighbor_df$sample_id)),
                                                nrow = length(unique(cna_neighbor_df$gene_id))))
  colnames(bucketed_cna_neighbor_df) <- unique(cna_neighbor_df$sample_id)
  rownames(bucketed_cna_neighbor_df) <- unique(cna_neighbor_df$gene_id)
  
  # Loop through all genes, and for each, get a 0/1 value for each patient indicating whether or
  # not that patient has a driver gene co-amplified or co-deleted on the same chromosome (e.g. 
  # whether the patient is in the DF after NA filtering)
  for (i in 1:nrow(bucketed_cna_neighbor_df)) {
    print(paste("Gene", paste(i, paste("/", nrow(bucketed_cna_neighbor_df)))))
    
    gene <- rownames(bucketed_cna_neighbor_df)[i]
    neighbor_df_sub <- cna_neighbor_df_sub[cna_neighbor_df_sub$gene_id == gene,]
    binary_vals <- unlist(lapply(colnames(bucketed_cna_neighbor_df), function(samp) {
      if(samp %fin% neighbor_df_sub$sample_id) {return(1)}
      else {return(0)}
    }))
    bucketed_cna_neighbor_df[i,] <- binary_vals
  }
  
  return(bucketed_cna_neighbor_df)
}

# First, just eliminate all rows in the data frame that have two or more NAs (because they do not 
# have a neighboring co-amplified or co-deleted driver)
cna_neighbor_df_sub <- cna_neighbor_df[!(is.na(cna_neighbor_df$upstr_copy_number) & 
                                           is.na(cna_neighbor_df$downstr_copy_number)),]
fwrite(cna_neighbor_df_sub, paste0(local_path, "CNV/neighboring_driver_cna_df_onlyAmpAndDel.csv"))

# Call function
bucketed_cna_neighbor_df <- bucket_flanking_cna_status_df(cna_neighbor_df)

# Write to a file
write.csv(bucketed_cna_neighbor_df, paste0(local_path, "CNV/neighboring_driver_cna_df_binary.csv"))
