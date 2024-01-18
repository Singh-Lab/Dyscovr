############################################################
### Process Methylation Data
### ASSOCIATED PUBLICATION INFORMATION
############################################################

library(rlang)
library(dplyr)
library(rlist)
library(parallel)
library(tictoc)
library(data.table)
library(abind)
library(doParallel)

# This file takes Level 3 methylation data from the GDC, processed via the 
# LiftOver pipeline

# See this paper for information about the differences between methylation data 
# processed by the LiftOver pipeline and TCGA Legacy Data: 
# https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30201-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471219302017%3Fshowall%3DT

# This is the link to the NIG GDC Legacy Archive, which has level 1 .idat files
# https://portal.gdc.cancer.gov/legacy-archive/search/f?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_format%22,%22value%22:%5B%22idat%22%5D%7D%7D%5D%7D


# Local path to directory containing Level 3 Methylation files
PATH <- getwd()
OUTPUT_PATH <- paste0(getwd(), "Output/Methylation/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)

# INPUT: Each patient has 1-2 individual Illumina methylation files, which can 
# be read in as large tables.

# OUTPUT: A large compiled table with the following format:
# Rows : Human genes (ENSG ID) 
# Columns : Patient TCGA ID (4-Digit)-Sample ID 
# Entries : Methylation M-value value 


############################################################
# PRE-DOWNLOAD FILTERING
############################################################
# Import the list of unique patient IDs that have all data types 
patient_id_list <- read.table(paste(PATH, "UNIQUE_PATIENT_IDS.txt", sep = ""), 
                              header = T)[,1]

# For an individual cancer type, import that manifest
manifest_filenames <- unlist(fread(paste0(PATH, "Methylation/MANIFEST.txt"), 
                                   header = T, sep = "\t")[,2])

# For pan-cancer analyses, combine all individual cancer type manifests into one
manifest_names <- list.files(paste0(PATH, "Methylation/", pattern = "MANIFEST"))
manifest <- fread(paste0(paste0(PATH, "Methylation/"), manifest_names[1]), 
                  header = T)
for (i in 2:length(manifest_names)) {
  mani <- fread(paste0(paste0(PATH, "Methylation/"), manifest_names[i]), 
                header = T)
  manifest <- rbind(manifest, mani)
}

# Write this to a file
fwrite(manifest, paste0(OUTPUT_PATH, "pan-cancer_methylation_manifest.csv"))

# Read back, as needed
manifest <- fread(paste0(OUTPUT_PATH, "pan-cancer_methylation_manifest.csv"), 
                  header = T)
manifest <- manifest[!(manifest$id == "\\N"),]
manifest_filenames <- manifest$filename

# Import the clinical DF
clinical_df <- read.csv(paste0(PATH, "Methylation/clinical_data.csv"), 
                             header = T, check.names = F)


#' Since methylation files are so large, we want to avoid downloading patient 
#' methylation files that we won't end up needing. So, we want to find the
#' intersection of all the case IDs we'll be using and all the possible
#' methylation case IDs (after we apply filters) and only download this set.
#' @param manifest the manifest file with all the filenames of individual 
#' patient methylation files
#' @param patient_ids a vector of patient TCGA IDs of interest
#' @param meth_clin_df_sub a clinical DF for methylation
get_patient_id_overlap <- function(manifest, patient_ids, meth_clin_df_sub) {
  manifest_ids <- unlist(lapply(1:length(manifest), function (i) {
    manifest_barcode <- unlist(strsplit(manifest[i], split = ".", fixed = T))[6]
    manifest_id <- unlist(strsplit(manifest_barcode, split = "-", fixed = T))[3]
    return(manifest_id)
  }))
  
  inter_ids <- intersect(manifest_ids, patient_ids)
  print(length(inter_ids))
  noninter_ids <- setdiff(manifest_ids, patient_ids)
  print(length(noninter_ids))
  
  missing_files <- setdiff(patient_ids, manifest_ids)
  print(length(missing_files))  # make sure this is 0
  
  # if not, print the IDs so that we can download them and try again
  if(length(missing_files) != 0) {
    missing_file_barcodes <- unlist(lapply(missing_files, function(f) {
      unique(meth_clin_df_sub[grepl(f, meth_clin_df_sub$case_submitter_id), 
                              'case_submitter_id'])
    }))
    
    print(paste(missing_file_barcodes, collapse = ","))
  }
  return(inter_ids)
}

# Get the IDs that intersect between the patient ID vector and the 
# manifest filenames
intersecting_ids <- get_patient_id_overlap(manifest_filenames, patient_ids)
# 7116 of 7864 intersect pan-cancer

# Write these IDs to a file
write(intersecting_ids, paste0(OUTPUT_PATH, "methylation_intersecting_ids.txt"))

# Read back if needed
intersecting_ids <- as.character(unlist(
  fread(paste(OUTPUT_PATH, "methylation_intersecting_ids.txt", sep = ""), 
                                              header = F)[,1]))

# Import the full set of human genes, or a subset of genes of interest
genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg38") 
gene_ids <- genes$ensembl_gene_id
gene_names <- unique(genes$external_gene_name)


############################################################
# IMPORT THE PROCESSED (LIFTOVER) METHYLATION TEXT FILES
############################################################
# Pan-Cancer
filenames <- list.files(paste0(PATH, "Methylation/"), pattern = ".gdc_hg38")  

# Limit to only intersecting IDs, if we haven't already filtered
filename_patient_ids <- unlist(lapply(filenames, function(x) 
  unlist(strsplit(unlist(strsplit(x, ".", fixed = T))[6], "-", fixed = T))[3]))
filenames <- filenames[which(filename_patient_ids %fin% intersecting_ids)] 


############################################################
# CLEAN METHYLATION DATA FRAMES
############################################################
#' Currently, the data has multiple genes listed per given probe in
#' the same row of the data frame. This function separates them out
#' into new rows (one gene per row, each with that probe's information)
#' @param filenames list of character filenames of methylation files, no 
#' attached path
#' @param path the local path to the methylation files
#' @param output_path a local path to write the cleaned methylation files to
clean_all_meth_dfs <- function(filenames, path, output_path)   {
  # Run through each patient to clean their files individually
  meth_vals_all_patients <- mclapply(filenames, function(x) {
    
    # Get the patient-sample ID for this file
    patient_sample_id <- paste(unlist(strsplit(
      unlist(strsplit(x, ".", fixed = T))[6], "-", fixed = T))[3:4], 
      collapse = "-")
    print(patient_sample_id)
    
    # Import the methylation table 
    meth_table <- fread(paste0(path, x), sep = "\t", header = T, 
                        select = c(Beta_value='numeric', 
                                   Gene_Symbol='character'))
    
    # Clean the table using helper function
    meth_table_clean <- clean_tab(meth_table)
    
    # Write the cleaned methylation table to a file
    fwrite(meth_table_clean, paste0(output_path, 
                                    paste0(patient_sample_id, 
                                           "_clean_methylation.csv")))
  })
} 


#' Helper function for cleaning a methylation table (separates gene names out
#' into new rows (one gene per row, each with that probe's information)
#' @param meth_table a Liftover methylation data frame
clean_tab <- function(meth_table) {
  # Filter out NA or empty values 
  meth_table <- meth_table %>% filter(!(is.na(Beta_value)) & 
                                        !(Gene_Symbol == "."))
  
  # For each row that may have more than one gene symbol, split apart
  # into multiple rows
  res <- mclapply(1:nrow(meth_table), function(i) {
    gene_symb <- meth_table$Gene_Symbol[i]
    spl <- unique(unlist(strsplit(gene_symb, ";", fixed = T)))
    if(length(spl) > 1) {
      new_tab <- data.frame(matrix(nrow=length(spl), ncol=2))
      colnames(new_tab) <- c("Beta_value", "Gene_Symbol")
      new_tab$Gene_Symbol <- spl
      new_tab$Beta_value <- rep(meth_table$Beta_value[i], nrow(new_tab))
      return(new_tab)
    } else {
      row <- meth_table[i,]
      row$Gene_Symbol <- spl
      return(row)
    }
  })
  tab <- as.data.frame(rbindlist(res, use.names = T, fill = T))
  
  # If there are duplicated genes and values, average them
  tab <- as.data.table(tab)
  tab <- tab[, list(Beta_value = mean(Beta_value)), Gene_Symbol]
  
  return(tab)
}

# Call function 
LOCAL_PATH <- paste0(PATH, "Methylation/")

# Potentially only look at a subset of new files we haven't already cleaned
patient_ids_curr <- unlist(lapply(list.files(OUTPUT_PATH, function(x) 
  unlist(strsplit(x, "-", fixed = T))[1])))
filename_patient_ids <- unlist(lapply(filenames, function(x) 
  unlist(strsplit(unlist(strsplit(x, ".", fixed = T))[6], "-", fixed = T))[3]))
filenames_curr <- filenames[which(filename_patient_ids %fin% patient_ids_curr)]
filenames_needed <- setdiff(filenames, filenames_curr)

clean_all_meth_dfs(filenames = filenames_needed, PATH = LOCAL_PATH, 
                   OUTPUT_PATH = OUTPUT_PATH)


############################################################
# CREATE COMPILED METHYLATION DF FROM THE SEPARATE FILES
############################################################
#' Takes in a list of cleaned methylation file names (each corresponding to 
#' one patient), optional lists of genes of interest or all genes (as well as a 
#' label to distinguish whether we're given a "full" set or a "subset"), and a 
#' list of intersecting patient IDs whose files we want to include. Creates the 
#' final compiled data frame and returns it.
#' @param filenames list of character file  names of cleaned methylation files, 
#' no attached path
#' @param path the local path to the methylation files
#' @param label "full" for using all genes, "subset" for using only certain 
#' genes of interest
#' @param returnBeta a T/F value; if T, return average Beta values, if F
#' return average M-values
#' @param gene_subset OPTIONAL: the ENSG IDs of a select subset of genes of 
#' interest; only when label is "subset"
#' @param beta_threshold OPTIONAL: beta threshold (above this value, we will 
#' consider this probe to be fully methylated)
create_methylation_df <- function(filenames, path, label, returnBeta, 
                                  gene_subset, beta_threshold) {
  
  methylation_dfs <- lapply(1:length(filenames), function(i) {
    filename <- filenames[i]
    print(paste(i, paste("/", length(filenames))))
    
    patient_sample_id <- unlist(strsplit(filename, "_", fixed = T))[1]
    print(patient_sample_id)
    
    # Read in their cleaned methylation table
    tic("upload table")
    meth_table <- fread(paste0(path, filename), header = T, 
                        select = c(Beta_value='numeric', 
                                   Gene_Symbol='character'))
    toc()
    
    # OPTIONAL: Subset the table to only beta values > threshold
    try (meth_table <- meth_table %>% dplyr::filter(Beta_value > beta_threshold),
         silent = T)
    toc()
    
    tic("make methylation state array")
    if (label == "subset") {
      meth_table <- meth_table[grep(paste(gene_subset, collapse = "|"), 
                                    meth_table$Gene_Symbol),]
    } 
    toc()
    
    # Convert Beta values to M-values, if desired
    if(!returnBeta) {
      # Logit of 1 is error
      meth_table$Beta_value[meth_table$Beta_value == 1] <- 0.9999999999
      meth_table$Beta_value <- unlist(lapply(meth_table$Beta_value, function(x) 
        log2(x / (1-x))))
      #print(head(meth_table))
    }
    
    # Replace the 'Beta_value' label with the patient sample ID
    tic("add patient sample id")
    colnames(meth_table)[which(colnames(meth_table) == 
                                 "Beta_value")] <- patient_sample_id
    toc()
    
    return(meth_table)
  })
  
  # Add the methylation arrays for this patient sample to the DF
  tic("combine into one DF")
  methylation_df <- Reduce(function(x, y) merge(x, y, by = "Gene_Symbol", 
                                                all = T), methylation_dfs)
  toc()
  
  print(head(methylation_df))
  
  # Return finished data frame
  return(methylation_df)
}


# Remove some big files and gc() to clear memory
rm(genes)
rm(gene_ids)
gc()


# Call function
cleaned_meth_filenames <- list.files(OUTPUT_PATH, pattern = "_clean_methylation")

#beta_thres <- 0.8
methylation_df <- create_methylation_df(filenames = cleaned_meth_filenames, 
                                        #beta_threshold = beta_thres,
                                        PATH = OUTPUT_PATH,
                                        label = "full", returnBeta = F)

# Replace NAs with 0's
methylation_df[is.na(methylation_df)] <- 0

# Limit to only cancer samples
methylation_df <- methylation_df[, c(T, grepl(
  "-0", colnames(methylation_df)[2:ncol(methylation_df)])), with = F]

# Write this table to a CSV 
#output_filename <- paste0(OUTPUT_PATH, "methylation_DF_Beta.csv")
output_filename <- paste0(OUTPUT_PATH, "methylation_DF_M.csv")
fwrite(methylation_df, output_filename)

# Read back this CSV, as needed
methylation_df <- fread(output_filename, header = T)

# OPTIONAL: check for missing values
apply(methylation_df_Beta, 2, function(y) 
  {if(length(y[y == 0]) > length(y)*0.5) print("missing > 50% of values")})


############################################################
# CREATE A BUCKETED METHYLATION FILE
############################################################
#' Given an input methylation file with the average Beta or
#' M-values for all patient samples, convert these numbers 
#' to a corresponding bucket:
#' Bucket 1 : low methylation (Beta < 0.3 or M < -0.8)
#' Bucket 2 : medium methylation (0.3 <= Beta < 0.7 or -0.8 <= M < 0.8)
#' Bucket 3 : high methylation (Beta >= 0.7 or M >= 0.8)
#' Note: Betas range from 0-1, Ms range from -0.85 to approx. 3 
#' (then goes up exponentially)
#' Rows are gene symbols and the columns are patient sample IDs (XXXX-XX),
#' while entries are a comma-separated list of 3 integers, with two 0's and 
#' one 1 corresponding to the methylation bucket, e.g. 0,1,0 for bucket 2
#' @param methylation_df the methylation data frame with the average Beta 
#' or M-values for all patient samples
#' @param isBeta a T/F value to indicate whether or not the input file
#' consists of Beta or M-values
create_bucketed_methylation_df <- function(methylation_df, isBeta) {
  new_df <- as.data.frame(apply(methylation_df[,2:ncol(methylation_df)], 
                                MARGIN = c(1,2), function(x) {
                                  
    # Establish the lower and upper bounds
    # Beta values
    lower_bound <- 0.3
    upper_bound <- 0.7
    # M-values
    if (!isBeta) {
      lower_bound <- -0.8
      upper_bound <- 0.8
    }
    
    # Establish the bucket for this value and return it 
    if(x < lower_bound) {return("c(1,0,0)")}
    else if ((x >= lower_bound) & (x < upper_bound)) {return("c(0,1,0)")}
    else {return("c(0,0,1)")}
    
  }))
  new_df$Gene_Symbol <- methylation_df$Gene_Symbol
  
  #print(head(new_df))
  return(new_df)
}


# Call this function
bucketed_methylation_df <- create_bucketed_methylation_df(methylation_df, 
                                                          isBeta = T)
rownames(bucketed_methylation_df) <- methylation_df$Gene_Symbol

# Write to file
fwrite(bucketed_methylation_df, paste0(OUTPUT_PATH, 
                                      "methylation_DF_bucketed_Beta.csv"))

