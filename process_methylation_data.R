############################################################
### Process Methylation Data
### Written By: Sara Geraghty, August 2020
############################################################

#library(methylumi)
#library(lumi)
#library(wateRmelon)
#library(minfi)
library(TCGAbiolinks)
library(rlang)
library(dplyr)
library(rlist)
library(parallel)
library(tictoc)
library(data.table)
library(abind)
library(doParallel)
library(sva)

# This file takes either:
# 1. Raw legacy methylation data from the GDC (level 1) OR
# 2. Level 3 methylation data from the GDC, processed via the LiftOver pipeline

# See this paper for information about the differences between methylation data 
# processed by the LiftOver pipeline and TCGA Legacy Data: 
# https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30201-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471219302017%3Fshowall%3Dtrue

# This is the link to the NIG GDC Legacy Archive, which has level 1 .idat files
# https://portal.gdc.cancer.gov/legacy-archive/search/f?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_format%22,%22value%22:%5B%22idat%22%5D%7D%7D%5D%7D


path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/"
# path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/"

output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Methylation/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Methylation/"

outpath_clean <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Methylation/Clean Methylation Files/"
#outpath_clean <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Methylation/Clean Methylation Files/"


############################################################
############################################################
### PART ONE: METHYLATION STATUS FILE
### Create a file that, for each gene and each patient, 
### lists a "0" or "1" value to indicate methylated or non-
### methylated
############################################################
############################################################
# A NOTE OF CAUTION: Beta values tend to be the most unreliable at extremes. The
# files produced from this part used 0.8 as a threshold for "fully methylated". 
# This may be an unreliable measure; differential methylation or bucketed 
# methylation is preferred if possible.

# INPUT: Each patient has 1-2 individual Illumina methylation files, which can be 
# read in as large tables.

# OUTPUT: A large compiled table with the following format:
  # Rows : Transcription Factors (ENSG ID) OR all genes in the human genome (ENSG ID) 
    # OR all genes in the human genome (gene name)
  # Columns : Patient TCGA ID (4-Digit)-Sample ID (01A is primary tumor, 11A is normal)
  # Entries : Methylation value 
    # ALT 1: 1 represents a likely methylation event on that gene or its promoter/ 
    # enhancer regions in that patient, 0 represents no methylation event
    # ALT 2: each entry is a vector of 3, with two 0s and one 1. Where the one is in
    # the vector indicates the bucket of the methylation 


############################################################
# PRE-DOWNLOAD FILTERING
############################################################
patient_ids <- fread(paste(path, "unique_brca_patient_ids.txt", sep = ""), 
                          header = TRUE)[,'x']
# If needed: patient_ids <- unlist(lapply(patient_ids, function(x) 
        # unlist(strsplit(x, split = " ", fixed = TRUE))[2]))
# patient_ids <- fread(paste(path, "unique_patient_ids.txt", sep = ""), 
                           # header = FALSE)[,1]
manifest_filenames <- unlist(fread(paste(path, "Methylation_Data/Level3/MANIFEST.txt", sep = ""), 
                                 header = TRUE, sep = "\t")[,2])

# For pan-cancer, combine all manifests into one
input_path_pc_meth <- "D:/"
manifest_names <- list.files(paste("Pan-Cancer Methylation Data/Level3/", 
                                   input_path_pc_meth, sep = ""), pattern = "MANIFEST")
manifest <- fread(paste(paste("Pan-Cancer Methylation Data/Level3/", 
                                   input_path_pc_meth, sep = ""), manifest_names[1], sep = ""), 
                       header = TRUE)
for (i in 2:length(manifest_names)) {
  mani <- fread(paste(paste("Pan-Cancer Methylation Data/Level3/", input_path_pc_meth, sep = ""), 
                           manifest_names[i], sep = ""), header = TRUE)
  manifest <- rbind(manifest, mani)
}

# Write this to a file
fwrite(manifest, paste(output_path, "pan-cancer_methylation_manifest.csv", sep = ""))

# Read back
manifest <- fread(paste(output_path, "pan-cancer_methylation_manifest.csv", sep = ""), 
                     header=TRUE)
manifest_filenames <- manifest$filename


#' Since methylation files are so large, we want to NOT download any patient 
#' methylation files that we won't end up needing. So, we want to find the
#' intersection of all the case IDs we'll be using and all the possible
#' methylation case IDs (after we apply filters) and only download this set.
#' @param manifest the manifest file with all the filenames of individual 
#' patient methylation files
#' @param patient_ids a vector of patient TCGA IDs of interest
get_patient_id_overlap <- function(manifest, patient_ids) {
  manifest_ids <- unlist(lapply(1:length(manifest), function (i) {
    manifest_barcode <- unlist(strsplit(manifest[i], split = ".", fixed = TRUE))[6]
    manifest_id <- unlist(strsplit(manifest_barcode, split = "-", fixed = TRUE))[3]
    return(manifest_id)
  }))

  inter_ids <- intersect(manifest_ids, patient_ids)
  noninter_ids <- setdiff(manifest_ids, patient_ids)
  return(inter_ids)
  #return(noninter_ids)
}

# Get the IDs that intersect between the patient ID vector and the manifest filenames
intersecting_ids <- get_patient_id_overlap(manifest_filenames, patient_ids)
  # 739 of 747 intersect in BRCA
  # 7116 of 7864 intersect pan-cancer

# Write these IDs to a file
write(intersecting_ids, paste(output_path, "methylation_intersecting_ids.txt", sep = ""))

# Read back if needed
intersecting_ids <- as.character(unlist(fread(paste(output_path, "methylation_intersecting_ids.txt", sep = ""), 
                          header = FALSE)[,1]))


############################################################
# OPTIONAL: IMPORT THE LIST OF TFs WE ARE INTERESTED IN
############################################################
idomain_prots <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_idomain.csv", sep = ",")[,2]
ibindingpos_prots <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_ibindingpos.csv", sep = ",")[,2]

# Otherwise, import the full set of human genes
genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg38") 
gene_ids <- genes$ensembl_gene_id
gene_names <- unique(genes$external_gene_name)


############################################################
# IMPORT THE PROCESSED (LIFTOVER) METHYLATION TEXT FILES
############################################################
# BRCA 
filenames <- list.files(paste(path, "Methylation_Data/Level3/", sep = ""), pattern = ".gdc_hg38")

# Pan-Cancer
# input_path_pc_meth <- "D:/"
# filenames <- list.files(paste(input_path_pc_meth, "Pan-Cancer Methylation Data/", sep = ""), 
  # pattern = ".gdc_hg38")  # 5121 files, excluding BRCA

# Limit to only intersecting IDs, if we haven't already filtered
filename_patient_ids <- unlist(lapply(filenames, function(x) 
  unlist(strsplit(unlist(strsplit(x, ".", fixed = TRUE))[6], "-", fixed = TRUE))[3]))
filenames <- filenames[which(filename_patient_ids %fin% intersecting_ids)]


############################################################
# CLEAN METHYLATION DATA FRAMES
############################################################
#' Currently, the data has multiple genes listed per given probe in
#' the same row of the data frame. This function separates them out
#' into new rows (one gene per row, each with that probe's information)
#' @param filenames list of character filenames of methylation files, no attached path
#' @param path the local path to the methylation files
#' @param output_path a local path to write the cleaned methylation files to
clean_all_meth_dfs <- function(filenames, path, output_path)   {
  # Run through each patient to clean their files individually
  meth_vals_all_patients <- mclapply(filenames, function(x) {
    
    # Get the patient-sample ID for this file
    patient_sample_id <- paste(unlist(strsplit(unlist(strsplit(x, ".", fixed = TRUE))[6], "-", fixed = TRUE))[3:4],
                               collapse = "-")
    print(patient_sample_id)

    # Extract the methylation table from this file
    meth_table <- fread(paste(path, x, sep = ""), sep = "\t", header = TRUE, 
            select = c(Beta_value='numeric', Gene_Symbol='character'))
    
    # Clean the table using helper function
    meth_table_clean <- clean_tab(meth_table)
    
    # Write the cleaned methylation table to a file
    fwrite(meth_table_clean, paste(output_path, paste(patient_sample_id, "_clean_methylation.csv", sep = ""), sep = ""))
  })
} 

#' Helper function for cleaning a methylation table (separates gene names out
#' into new rows (one gene per row, each with that probe's information)
#' @param meth_table a methylation data frame
clean_tab <- function(meth_table) {
  # Filter out NA or empty values 
  meth_table <- meth_table %>% filter(!(is.na(Beta_value)) & !(Gene_Symbol == "."))
  # For each row that may have more than one gene symbol, split apart
  # into multiple rows
  res <- mclapply(1:nrow(meth_table), function(i) {
    gene_symb <- meth_table$Gene_Symbol[i]
    spl <- unique(unlist(strsplit(gene_symb, ";", fixed = TRUE)))
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
  tab <- as.data.frame(rbindlist(res, use.names = TRUE, fill = TRUE))
  
  # If there are duplicated genes and values, average them
  #tab <- aggregate(tab$Beta_value, by = c(tab$Gene_Symbol), data = tab, FUN = mean)
  tab <- as.data.table(tab)
  tab <- tab[, list(Beta_value = mean(Beta_value)), Gene_Symbol]
  
  print(dim(tab))
  
  return(tab)
}

# Call function 
local_path <- paste(path, "Methylation_Data/Level3/", sep = "")
# local_path <- paste(input_path_pc_meth, "Pan-Cancer Methylation Data/", sep = "")
clean_all_meth_dfs(filenames = filenames, path = local_path, output_path = outpath_clean)


############################################################
# CREATE COMPILED METHYLATION DF FROM THE SEPARATE FILES
############################################################
#' Takes in a list of cleaned methylation filenames (each corresponding to one patient), 
#' optional lists of proteins of interest or all genes (as well as a label to 
#' distinguish whether we're given a "full" set or a "subset"), and a list of 
#' intersecting patient IDs whose files we want to include. Creates the final 
#' data frame described at top of file, fills it in, and returns it.
#' @param filenames list of character filenames of cleaned methylation files, 
#' no attached path
#' @param path the local path to the methylation files
#' @param label "full" for using all genes, "subset" for using only certain genes of interest
#' @param returnBeta a TRUE/FALSE value; if TRUE, return average Beta values, if FALSE
#' return average M-values
#' @param binding_prots OPTIONAL: the ENSG IDs of a select subset of genes of interest; 
#' only provided when label is "subset"
#' @param beta_threshold OPTIONAL: beta threshold (above this value, we will consider this
#' probe to be fully methylated)
create_methylation_df <- function(filenames, path, label, returnBeta, 
                                  binding_prots, beta_threshold) {

  methylation_dfs <- lapply(1:length(filenames), function(i) {
    filename <- filenames[i]
    print(paste(i, paste("/", length(filenames))))
    
    patient_sample_id <- unlist(strsplit(filename, "_", fixed = TRUE))[1]
    print(patient_sample_id)
    
    # Read in their cleaned methylation table
    tic("upload table")
    meth_table <- fread(paste(path, filename, sep = ""), header = TRUE, 
                        select = c(Beta_value='numeric', Gene_Symbol='character'))
    toc()
    
    # OPTIONAL: Subset the table to only beta values > threshold
    try (meth_table <- meth_table %>% dplyr::filter(Beta_value > beta_threshold),
         silent = TRUE)
    toc()
    
    tic("make methylation state array")
    if (label == "subset") {
      meth_table <- meth_table[grep(paste(binding_prots, collapse = "|"), 
                                            meth_table$Gene_Symbol),]
      # Alternative 1: If there are some methylation markings for each gene in this patient, 
      # add 1 (otherwise add 0)
      #meth_state_array <- as.character(unlist(lapply(binding_prots, function(id) 
        #if(id %fin% meth_table$Gene_Symbol) {return(1)} else{return(0)})))
    } 
    #else if (label == "full") {
      # Alternative 1: If there are some methylation markings for each gene in this patient, 
      # add 1 (otherwise add 0)
      # meth_state_array <- as.character(unlist(lapply(gene_names, function(id) 
        #if (id %fin% meth_table$Gene_Symbol) {return(1)} else {return(0)})))
      
      #return(meth_table)
    #} else {
      #print(paste("Invalid label:", label))
      #return(0)
    #}
    toc()
    
    # Convert Beta values to M-values, if desired
    if(!returnBeta) {
      meth_table$Beta_value[meth_table$Beta_value == 1] <- 0.9999999999    # Logit of 1 is error
      meth_table$Beta_value <- unlist(lapply(meth_table$Beta_value, function(x) log2(x / (1-x))))
      print(head(meth_table))
    }
    
    # Replace the 'Beta_value' label with the patient sample ID
    tic("add patient sample id")
    #meth_state_array <- prepend(meth_state_array, patient_id_full, before = 1) 
    colnames(meth_table)[which(colnames(meth_table) == "Beta_value")] <- patient_sample_id
    toc()
    
    return(meth_table)
  })
  
  print(head(methylation_dfs))
  
  # Add the methylation arrays for this patient sample to the DF
  tic("combine into one DF")
  #methylation_df <- do.call("cbind", methylation_columns)
  methylation_df <- Reduce(function(x, y) merge(x, y, by = "Gene_Symbol", all = T), methylation_dfs)
  toc()
  
  print(head(methylation_df))
  
  # Return finished data frame
  return(methylation_df)
}


# Run the above function to get compiled data frame
#beta_thres <- 0.8

# Remove some big files and gc() to clear memory
rm(genes)
rm(gene_ids)
gc()

cleaned_meth_filenames <- list.files(outpath_clean, pattern = "_clean_methylation")


methylation_df <- create_methylation_df(filenames = cleaned_meth_filenames, 
                                        #beta_threshold = beta_thres,
                                        path = outpath_clean,
                                        label = "full", returnBeta = TRUE)

# Replace NAs with 0's
methylation_df[is.na(methylation_df)] <- 0

# Write this table to a CSV 
# output_filename <- paste(output_path, paste("methylation_DF_Beta_", paste(beta_thres, ".csv", sep = ""), 
                              #sep = ""), sep = "")
output_filename <- paste0(output_path, "methylation_DF_Beta.csv")
output_filename <- paste0(output_path, "methylation_DF_M.csv")

fwrite(methylation_df, output_filename)

# Read back this CSV
methylation_df <- fread(output_filename, header = TRUE)

# OPT: check for missing values
apply(methylation_df_Beta, 2, function(y) {if(length(y[y == 0]) > length(y)*0.5) print("missing > 50% of values")})

############################################################
# CREATE A BUCKETED METHYLATION FILE
############################################################
#' Given an input methylation file with the average Beta or
#' M-values for all patient samples, convert these numbers 
#' to a corresponding bucket:
#' Bucket 1 : low methylation (Beta < 0.3 or M < -0.8)
#' Bucket 2 : medium methylation (0.3 <= Beta < 0.7 or -0.8 <= M < 0.8)
#' Bucket 3 : high methylation (Beta >= 0.7 or M >= 0.8)
#' Note: Betas range from 0-1, Ms range from -0.85 to approx. 3 (then goes up exponentially)
#' Rows will be gene symbol and the columns will be patient sample IDs (XXXX-XX),
#' while entries will be a comma-separated list of 3 integers, with two 0's and 
#' one 1 corresponding to the bucket the methylation is in, e.g. 0,1,0 for bucket 2
#' @param methylation_df the methylation data frame with the average Beta or M-values
#' for all patient samples
#' @param isBeta a TRUE/FALSE value to indicate whether or not the input file
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
  
  print(head(new_df))
  return(new_df)
}


# Call this function
bucketed_methylation_df <- create_bucketed_methylation_df(methylation_df, isBeta = TRUE)
rownames(bucketed_methylation_df) <- methylation_df$Gene_Symbol

# Write to file
fwrite(bucketed_methylation_df, paste(output_path, "methylation_DF_bucketed_Beta.csv", sep = ""))


















# NOTE: I NO LONGER NEED THIS SECTION, SINCE WE ARE RUNNING MODEL ON A 
# PER-SAMPLE RATHER THAN PER-PATIENT BASIS. THE SECTION REMAINS HERE FOR REFERENCE.


############################################################
############################################################
### PART TWO: DIFFERENTIAL METHYLATION FILE
### Create files that, for each gene and each patient, 
### calculates differential methylation between patient 
### tumor and normal samples 
############################################################
############################################################


############################################################
# GET TUMOR-NORMAL MATCHED PATIENT IDs
############################################################
#' Get the number of patients with tumor-normal matched files and their IDs
#' (return the IDs of these patients)
#' @param intersecting_ids a vector of the patient IDs of interest that have
#' corresponding methylation filenames in the manifest
#' @param filenames a vector of the methylation filenames
get_matched_patient_IDs <- function(intersecting_ids, filenames) {
  matched_patient_IDs <- unlist(lapply(intersecting_ids, function(id) {
    filenames_patient <- filenames[grepl(id, filenames)]
    
    # Remove any "character(0)" elements 
    filenames_patient <- unlist(lapply(filenames_patient, function(x) 
      if(length(x) > 0) {return(x)}))
    
    # Check how many filenames are here; if it's more than 1, it's possible
    # this patient might have tumor-normal matched files
    if(length(filenames_patient) > 1) {
      
      # Get just the TCGA ID for each filename (TCGA-XX-XXXX-XXX)
      tcga_ids <- unlist(lapply(filenames_patient, function(nam) {
        return(unlist(strsplit(unlist(strsplit(nam, ".", fixed = TRUE))[6], "-", 
                               fixed = TRUE))[1:4])
      }))

      # Ensure that these IDs are in the right place in the filename (in the patient ID slot); 
      # will be in the subset
      tcga_ids_matching <- unlist(lapply(tcga_ids, function(x) 
        ifelse(id %fin% x, TRUE, FALSE)))
      if (length(tcga_ids_matching[tcga_ids_matching == TRUE]) > 1) {
        # Check if there is at least one of each file type (cancer & normal)
        if ((length(filenames_patient[grepl("11", tcga_ids)]) >= 1) &    # normal
            ((length(filenames_patient[grepl("01", tcga_ids)]) >= 1) |   # tumor
             (length(filenames_patient[grepl("06", tcga_ids)]) >= 1))) {  # metastatic tumor
          return(id)
        }
      } 
    }
  }))
  return(matched_patient_IDs)
}

matched_patient_IDs <- get_matched_patient_IDs(intersecting_ids, filenames)
# BRCA: 87 patients
# Pan-Cancer: 352 patients, excluding BRCA (of 7116)


############################################################
# GET SUMMARIZING METHYLATION DF
############################################################
# INPUT: Patients with multiple, matched tumor-normal methylation files
# OUTPUT: A large compiled table with the following format:
# Rows : Transcription Factors (ENSG ID) OR all genes in the human genome (ENSG ID)
# Columns : Patient TCGA ID (4-Digit)-Sample ID plus 3-digit sample ID (-01A, -11A, etc.) 
# Entries: Average beta value across the length of the gene; average M-value across the 
# length of the gene (M = log2(Beta/(1-Beta)))

#' Looks at all patients that have matched tumor-normal data and calculates the 
#' average Beta value and M-value, on a per-gene basis, for both tumor and normal. 
#' Writes these files on a patient-level basis. Also outputs a data frame of the 
#' form described above.
#' @param filenames list of character filenames of cleaned methylation files, 
#' no attached path
#' @param path the local path to the methylation files
#' @param matched_patient_IDs the patient IDs (4-digit TCGA ID) that have all four 
#' data types & tumor-normal matched files
#' @param thres OPTIONAL: the Beta threshold for a fully methylated probe
get_summarizing_methylation_df <- function(filenames, path, matched_patient_IDs, thres) {
    
    # Loop through all patient files
    avg_tabs <- mclapply(1:length(matched_patient_IDs), function(i) {
      
      patient_id <- matched_patient_IDs[i]
      
      # Get the filenames corresponding to this patient
      filenames_patient <- filenames[grepl(patient_id, filenames)]
      meth_table_normal_clean <- filenames_patient[grepl("-11")]   # or, grepl("-norm")
      meth_table_tumor_clean <- filenames_patient[grepl("-0")]     # or, grepl("-tum")
      
      # Import the cleaned tumor & normal methylation files for this patient
      meth_table_normal_clean <- fread(paste(path, paste("Methylation/Clean Methylation Files/", 
                                                            meth_table_normal_clean, sep = ""), sep = ""))
      meth_table_tumor_clean <- fread(paste(path, paste("Methylation/Clean Methylation Files/", 
                                                           meth_table_tumor_clean, sep = ""), sep = ""))
  
      # Get the average beta and M-values for each of these and write them to files
      # using helper function
      tic("get & write files")
      avg_tab_normal <- merge_and_write_files(meth_table_normal_clean, output_path, 
                                                   thres, "11", patient_id) 
      avg_tab_tumor <- merge_and_write_files(meth_table_tumor_clean, output_path, 
                                                  thres, "0", patient_id) 
      toc()
      
      # OPT: merge the two tables based on Gene_Symbol, for tumor and for normal separately
      #full_normal_tab <- merge(avg_tab_normal, num_probes_tab_normal, by = "Gene_Symbol")
      #full_tumor_tab <- merge(avg_tab_tumor, num_probes_tab_tumor, by = "Gene_Symbol")
      
      # Merge the tumor and normal data
      if(!(length(avg_tab_normal) == 0 & length(avg_tab_tumor) == 0)) {
        colnames(avg_tab_normal)[colnames(avg_tab_normal) == "Beta_value"] <- paste(patient_id, "-11", sep = "")
        colnames(avg_tab_tumor)[colnames(avg_tab_tumor) == "Beta_value"] <- paste(patient_id, "-0", sep = "")
        
        tic("merge tumor and normal")
        full_beta_tab <- merge(avg_tab_normal, avg_tab_tumor, by = "Gene_Symbol")
        print(head(full_beta_tab))
        toc()
        
        return(full_beta_tab)
      }
      else {return(NA)}
      # colnames(num_probes_tab_tumor)[colnames(num_probes_tab_tumor) == "Num.Fully.Methylated.Probes"] <- paste(patient_id, tumor_marker, sep = "")
      # colnames(num_probes_tab_normal)[colnames(num_probes_tab_normal) == "Num.Fully.Methylated.Probes"] <- paste(patient_id, normal_marker, sep = "")
      # full_num_probes_tab <- merge(num_probes_tab_tumor, num_probes_tab_normal, by = "Gene_Symbol")
      
      
      # return(full_num_probes_tab)
    
  }) 
  meth_vals_all_patients <- meth_vals_all_patients[!is.na(meth_vals_all_patients)]
  
  # Add the methylation arrays for this patient to the DF
  methylation_df <- Reduce(merge.all, meth_vals_all_patients)
  print(methylation_df)
  
  # Return finished data frame
  return(methylation_df)
}

thres <- 0.8

# Call this function to get the differential methylation table
input_path <- paste(input_path_pc_meth, "Pan-Cancer Methylation Data/", sep = "")
summarizing_methylation_df <- get_summarizing_methylation_df(filenames = filenames, 
                                                             path = input_path, 
                                                             matched_patient_IDs = matched_patient_IDs, 
                                                             thres = thres)
matched_patient_IDs_w_labs <- unlist(lapply(matched_patient_IDs, function(x) 
  c(paste(x, "-norm", sep = ""), paste(x, "-tumor", sep = ""))))
colnames(summarizing_methylation_df) <- c("Gene_Symbol", matched_patient_IDs_w_labs)

# Write to a file
output_path <- 'C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Methylation/'
fwrite(summarizing_methylation_df, paste(output_path, "Differential Methylation/average_beta_per_gene.csv", sep = ""))



### HELPER FUNCTIONS ###


#' Given a cleaned methylation table(s), either tumor or normal, it calculates
#' the average Beta table (where Beta values are averaged from all probes across
#' the length of the gene), average M-value, and the number of fully methylated 
#' probes table. Writes the cleaned tables and calculated tables to files at the 
#' given output path and returns the average Beta table.
#' @param meth_table_clean a single or list of cleaned methylation tables
#' @param output_path a path to write the cleaned methylation files, average Beta
#' files, and number of probes files
#' @param thres a Beta threshold for a fully methylated probe
#' @param tum_or_norm a value indicating if the methylation files are for tumor ("0")
#' or normal ("11") samples
#' @param patient_id the 4-digit patient identifier for the given patient
merge_and_write_files <- function(meth_table_clean, output_path, thres, 
                                  tum_or_norm, patient_id) {
  # Create an empty object that will hold the merged average beta table
  avg_beta_table <- NA

  # If there are multiple data frames, average them together before writing to a file
  if(!is.data.frame(meth_table_clean) & !length(meth_table_clean) == 0) {
    for (i in 1:length(meth_table_clean)) {
      meth_tab <- meth_table_clean[[i]]
      #print(meth_tab)
      avg_beta_tab <- get_and_write_comp_files(meth_tab, tum_or_norm, patient_id, 
                                               output_path, thres)
      # If this is the first table we've looked at, we don't have to do any merging yet
      if (is.na(avg_beta_table)) {avg_beta_table <- avg_beta_tab}
      # Otherwise, merge this table into the main table
      else {
        avg_beta_table <- merge(avg_beta_table, avg_beta_tab, by.x = "Gene_Symbol", 
                                by.y = "Gene_Symbol")
        avg_beta_table <- data.frame(Gene_Symbol = avg_beta_tab$Gene_Symbol, 
                                     Beta_value = rowMeans(avg_beta_tab[,-1]))
      }
    }
  # If there is only one data frame, write this to a file 
  } else {
    avg_beta_table <- get_and_write_comp_files(meth_tab, tum_or_norm, patient_id, output_path, thres)
  }
  return(avg_beta_table)
}

#' A helper function for the above to write the files to the given path and to 
#' actively call the helper functions that will calculate the average Betas 
#' and the number of fully methylated probes
#' @param meth_tab a cleaned methylation table
#' @param tum_or_norm a value indicating if the methylation files are for tumor ("0")
#' or normal ("11") samples
#' @param patient_id the 4-digit patient identifier for the given patient
#' @param output_path a path to write the cleaned methylation files, average Beta
#' files, and number of probes files
#' @param thres a Beta threshold for a fully methylated probe
get_and_write_comp_files <- function(meth_tab, tum_or_norm, patient_id, output_path, thres) {
  # Get the average Beta value table using helper function 
  # (average of all probes across gene length)
  avg_beta_tab <- get_avg_betas(meth_tab)
  
  # OPT: Convert these to M-values and add that as a column
  # avg_beta_tab <- get_avg_mvals(avg_beta_tab)
  
  # Write to a file
  fwrite(avg_beta_tab, paste(output_path, paste(patient_id, paste(tum_or_norm, paste(i, "_avg_betas_per_gene.csv", sep = ""), sep = "_"), sep = "-"), sep = ""))
  
  # Get the number of fully methylated probes for each gene using helper function
  # and write this to a file
  num_probes_tab <- get_abs_num(meth_tab, thres)
  fwrite(num_probes_tab, paste(output_path, paste(patient_id, paste(tum_or_norm, paste(i, "_num_methylated_probes_per_gene.csv", sep = ""), sep = "_"), sep = "-"), sep = ""))
  
  return(avg_beta_tab)
}

#' Takes a cleaned methylation table and averages the beta values across genes 
#' that have multiple probes
#' @param meth_table a cleaned methylation table
get_avg_betas <- function(meth_table) {
  # Aggregate, grouping by gene name and applying the mean.default() function
  new_tab <- aggregate(meth_table[,'Beta_value'], by = list(meth_table$Gene_Symbol), 
                       FUN = mean.default)
  colnames(new_tab)[1] <- "Gene_Symbol"
  return(new_tab)
}

#' Takes in an average Beta table, converts all the Beta values to M-values
#' (logit) and returns the table with a column for the mean M-value
#' @param avg_beta_tab
get_avg_mvals <- function(avg_beta_tab) {
  mvals <- unlist(lapply(meth_table$Beta_value, function(x) log(x / (1-x))))
  avg_beta_tab$M_value <- mvals
  return(avg_beta_tab)
}

#' Takes a cleaned methylation table and counts the number of significantly 
#' methylated probes per gene
#' @param meth_table a cleaned methylation table
#' @param thres a Beta threshold for a fully methylated probe
get_abs_num <- function(meth_table, thres) {
  filt_table <- meth_table %>% filter(Beta_value > thres)
  new_tab <- data.frame('Gene_Symbol' = unique(filt_table$Gene_Symbol))
  abs_num_col <- unlist(mclapply(1:nrow(new_tab), function(i) {
    gene <- new_tab$Gene_Symbol[i]
    return(nrow(meth_table %>% filter(Gene_Symbol == gene)))
  }))
  new_tab$Num.Fully.Methylated.Probes <- abs_num_col
  return(new_tab)
}

#' Helper function for merging methylation tables by gene symbol
#' @param x methylation table 1
#' @param y methylation table 2
merge.all <- function(x, y) {
  merge(x, y, all=TRUE, by="Gene_Symbol")
}


############################################################
# ALT: ADD M-VALUES TO THE AVG BETA PER-GENE FILES
############################################################
#' Takes the average beta per gene files for each patient 
#' sample and adds another column for the M-value
#' M = log2(Beta / 1-Beta)
#' @param path a local path to the folder that has the 
#' avg_betas_per_gene files
add_mval_column <- function(path) {
  filenames <- list.files(path, pattern = "avg_betas_per_gene")
  for (i in 1:length(filenames)) {
    file <- data.table::fread(paste(path, filenames[i], sep = ""), sep = ",", header = TRUE)
    colnames(file)[2] <- "Avg.Beta"
    mvals <- unlist(lapply(file$Avg.Beta, function(x) log(x / (1-x))))
    file$Avg.M <- mvals
    fwrite(file, paste(path, filenames[i], sep = ""))
  }
}


############################################################
# GET DIFFERENTIAL AVG. BETA PER-GENE FILES 
############################################################
#' Takes an output file of the form given above and, for each patient, calculates 
#' the absolute difference in the average Beta value between tumor and normal
#' for each gene. 
#' @param summarizing_methylation_df a data frame like that produced from the above function
get_diff_beta_file <- function(summarizing_methylation_df) {
  # Loop over all patient columns
  new_cols <- lapply(seq(2, ncol(summarizing_methylation_df), by=2), function(i) {
    # Get the corresponding tumor and normal columns
    col_norm <- summarizing_methylation_df[,i]
    col_tumor <- summarizing_methylation_df[,i+1]
    # Calculate the difference in the average Beta between them
    col_diff <- col_tumor - col_norm
    # Add the patient ID 
    patient <- unlist(strsplit(colnames(summarizing_methylation_df)[i], "-", fixed = TRUE))[1]
    return(c(patient, col_diff))
  })
  
  df <- as.data.frame(do.call("cbind", new_cols))
  
  # Make the top row (patient IDs) the column names
  colnames(df) <- df[1,]
  df <- df[-1,]
  
  # Add the gene symbol column to this new data frame
  df <- cbind(summarizing_methylation_df[,1], df)
  colnames(df)[1] <- "Gene.Symbol"
  
  return(df)
}

diff_beta_file <- get_diff_beta_file(summarizing_methylation_df)

# Write to a CSV
fwrite(diff_beta_file, paste(output_path, "Differential Methylation/differential_avg_beta_per_gene.csv", sep = ""))


############################################################
# GET NUM. FULLY METHYLATED PROBES PER-GENE FILES
############################################################
# INPUT: Patients with multiple, matched tumor-normal methylation files
# OUTPUT: A large compiled table with the following format:
# Rows : Transcription Factors (ENSG ID) OR all genes in the human genome (ENSG ID)
# Columns : Patient TCGA ID (4-Digit)-Sample ID plus sample (-norm, -tumor), which contain 
  # absolute number of fully methylated probes, plus one column with
  # the differential number of probes column per patient

#' Looks at all patients that have tumor-normal matched data and calculates the 
#' number of confidently methylated probes in the sample. Outputs a data frame 
#' of the form described above.
#' @param path the local path to the methylation files
#' @param num_diff_probes_files a data frame like that produced from the above function
get_diff_num_meth_probes_file <- function(path, num_diff_probes_files) {
  # Get all the patient IDs
  patients <- unlist(unique(lapply(num_diff_probes_files, 
                                   function(x) unlist(strsplit(x, "-", fixed = TRUE))[1])))
  
  # Loop through each and get files
  cols <- lapply(1:length(patients), function(i) {
    patient_filenames <- num_diff_probes_files[grepl(patients[i], num_diff_probes_files)]
    patient_filenames_tum <- patient_filenames[grepl("-0", patient_filenames)]
    patient_filenames_norm <- patient_filenames[grepl("-11", patient_filenames)]
    
    tum_df <- data.frame(matrix(nrow = 0, ncol = 0))
    norm_df <- data.frame(matrix(nrow = 0, ncol = 0))
    
    if(!(length(patient_filenames_norm) == 0 | length(patient_filenames_tum) == 0)) {
      if (length(patient_filenames_tum) > 1) {
        tum_df <- get_average_within_cat(patient_filenames_tum)
      } else {
        file <- data.table::fread(paste(path, patient_filenames_tum[1], sep = "") , sep = ",", header = TRUE)
        tum_df <- file
      }
      if(length(patient_filenames_norm) > 1) {
        norm_df <- get_average_within_cat(patient_filenames_norm)
      } else {
        file <- data.table::fread(paste(path, patient_filenames_norm[1], sep = "") , sep = ",", header = TRUE)
        norm_df <- file
      }
      
      # Merge the tumor and normal data frames by Gene_Symbol
      merged_df <- as.data.frame(merge(tum_df, norm_df, by.x = "Gene_Symbol", 
                                       by.y = "Gene_Symbol", all.x = FALSE, all.y = FALSE)) 

      # Get the difference (tumor - normal) of these columns
      print(merged_df)
      merged_df$Diff_Tumor_Norm <- merged_df[,2] - merged_df[,3]
      
      # Add the patient ID
      tum_lab <- paste(patients[i], "-tumor", sep = "")
      norm_lab <- paste(patients[i], "-norm", sep = "")
      diff_lab <- paste(patients[i], "-diff_Tumor_Minus_Norm", sep = "")
      colnames(merged_df)[2:4] <- c(tum_lab, norm_lab, diff_lab)
      
      print(merged_df)
      
      # Write to file
      filename <- paste(patients[i], "_indivAndDiffMethylation.csv", sep = "")
      fwrite(merged_df, paste(path, filename, sep=""))
      
      # Return
      return(merged_df)
    }
  })
  
  # Bind these columns together into a data.frame
  avg_df <- cols[1]
  for (i in 2:length(cols)) {
    avg_df <- merge(avg_df, cols[i], by.x = "Gene_Symbol", by.y = "Gene_Symbol", all.x = TRUE, all.y = TRUE)
  }
  
  # Get the gene symbols
  #samp_file <- data.table::fread(paste(path, patient_filenames[1], sep = "") , sep = ",", header = TRUE)
  #gene_symb <- samp_file$Gene_Symbol
  #avg_df <- cbind(gene_symb, avg_df)
  #colnames(avg_df)[1] <- "Gene_Symbol"
  
  return(avg_df)
}


#' If there is more than one tumor or normal file for a given patient, this function
#' is called to merge them into one file (averaging the columns together for each gene)
#' @param filenames the vector of methylation filenames for a given patient and 
#' category (tumor or normal)
get_average_within_cat <- function(filenames) {
  # Get the first file to start with
  avg_df <- fread(paste(path, filenames[1]) , sep = ",", header = TRUE)
  
  # Merge the rest into the first one by Gene_Symbol
  for (filename in filenames[2:length(filenames)]) {
    # Upload file
    file <- fread(paste(path, filename) , sep = ",", header = TRUE)
    # Merge into avg DF
    avg_df <- merge(avg_df, file, by.x = "Gene_Symbol", by.y = "Gene_Symbol", 
                    all.x = TRUE, all.y = TRUE)
  }
  
  # Average all the columns together for each gene
  if (ncol(avg_df) > 2) {
    avg_df$Means <- rowMeans(avg_df[,-1])
  } else {avg_df$Means <- avg_df[,2]}
  
  avg_df <- as.data.frame(avg_df)
  print(avg_df)
  return(avg_df)
}
#
#
#

# Main 
num_diff_probes_files <- list.files(paste(output_path, "Clean Methylation Files/", sep = ""), 
                          pattern = "num_methylated_probes")
samp_path <- paste(output_path, "Clean Methylation Files/", sep = "")

# Call function
num_diff_probes_file <- get_diff_num_meth_probes_file(samp_path, num_diff_probes_files)

# Write result to CSV
fwrite(num_diff_probes_file, paste(output_path, "Differential Methylation/differential_num_probes_per_gene.csv", sep = ""))
































































#
#
#
#
############################################################
############################################################
# ARCHIVES: GET BETAS AND DIFFERENTIAL METHYLATION USING PACKAGES
############################################################
############################################################
# NOTE: THESE ARE OUT OF DATE AND DO NOT APPLY TO INDIVIDUAL TUMOR-NORMAL
# DIFFERENTIAL METHYLATION (ONLY DIFFERENTIAL METHYLATION BETWEEN GROUPS). 
# DO NOT USE.



############################################################
############################################################
# OPT 1: IMPORT THE RAW (LEGACY) METHYLATION TEXT FILES
############################################################
############################################################

suppressPackageStartupMessages(library(methylumi, quietly=TRUE))

barcodes <- manifest # FIX THIS (get all the barcodes of the patients we want)
methylumi_set <- methylumIDAT(barcodes = barcodes, pdat = NULL, parallel = TRUE, 
                              idatPath = paste(path, "Methylation_Data/Level1/", sep = ""))
# input_path_pc_meth <- "D:/"
#methylumi_set <- methylumIDAT(barcodes = barcodes, pdat = NULL, parallel = TRUE, 
#idatPath = paste(input_path_pc_meth, 
# "Pan-Cancer Methylation Data/Level1/", sep = ""))

############################################################
# PREPROCESS THE LEGACY METHYLATION FILES USING METHYLUMI
############################################################

#TODO: Fill this section in if desired. Here are tutorials to follow:
# https://bioconductor.riken.jp/packages/3.3/bioc/vignettes/lumi/inst/doc/methylationAnalysis.pdf
# https://master.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#differential-methylation-analysis


############################################################
# NORMALIZE RAW VALUES USING BMIQ & CREATE COMPILED DF
############################################################
#' Apply BMIQ normalization to these files of raw Beta values
#' @param path a local path to the methylation files
#' @param filenames a vector of methylation filenames of interest
#' @param intersecting_ids all the patient IDs of interest that
#' also have an associated methylation file listed in the manifest
normalize_betas <- function(path, filenames, intersecting_ids) {
  methylation_objects <- lapply(1:length(filenames), function(i) {
    print(paste(i, paste("/", length(filenames))))
    filename <- filenames[i]
    
    patient_id <- unlist(strsplit(filename, "-", fixed = TRUE))[5]
    sample_id <- unlist(strsplit(filename, "-", fixed = TRUE))[6]
    
    # Ensure that this patient is one with all types of data
    if (patient_id %fin% intersecting_ids) {
      # We will use the methylumi R package to read Illumina methylation data
      # into a MethylumiSet object
      full_filename <- paste(path, paste("Methylation_Data/", filename, sep = ""), sep = "")
      print(full_filename)
      methyl_obj <- methylumiR(system.file(full_filename, package = "methylumi"), 
                               qcfile = system.file(full_filename),
                               sep = "\t", header = TRUE)
      
      # Apply the BMIQ function
      bmiq_results <- BMIQ(methyl_obj)
      print(bmiq_results)
      
      # Write the results to a file
      output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Methylation/BMIQ-Normalized/"
      write(bmiq_results, paste(output_path, paste(patient_id, 
                                                   paste(sample_id, "_normalizationSummary.txt"), sep = "_"), sep = ""))
      
      # Return this object
      return(bmiq_results)
      
      # Extract just the normalized beta values for the sample
      #bmiq_betas <- bmiq_results[[1]]
      #print(head(bmiq_betas))
      
      # Prepend the patient_id-sample_id
      #patient_id_full <- paste(patient_id, sample_id, sep = "-")
      #beta_array <- prepend(bmiq_betas, patient_id_full, before = 1)
      
      # Return this column
      #return(beta_array)
    }
  })
  #methylation_df <- do.call("cbind", methylation_columns)
  
  # Add patient IDs as column headers
  #colnames(methylation_df) <- methylation_df[1,]
  #methylation_df <- methylation_df[-1,]  
  
  # Return this list of normalized methylation objects
  return(methylation_objects)
  
}

methylation_objects <- normalize_betas(path, filenames, intersecting_ids)


# Protocol using R package 'minfi'
# https://bioconductor.org/packages/release/bioc/manuals/minfi/man/minfi.pdf

#' Import Level 3methylation files using minfi function 'readTCGA'. These files are 
#' then returned as a list
#' @param path the local path to the folder containing all the level 3 
#' methylation files
#' @param sample_sheet a data frame with all the sample information for this cohort
import_methfiles_minfi <- function(path, sample_sheet) {
  
  filenames <- list.files(path, pattern = ".gdc_hg38")
  filenames_full <- unlist(lapply(filenames, function(x) paste(path, x, sep = "")))
  
  # Reads in the file as a GenomicRatioSet object & extracts important components
  dfs <- lapply(filenames_full, function(fn) {
    
    # Subset the sample sheet to just this sample
    sample_id <- paste(unlist(strsplit(unlist(strsplit(fn, ".", fixed = TRUE))[6], "-", fixed = TRUE))[1:4], 
                        collapse = "-")
    #sample_sheet_sub <- sample_sheet[sample_sheet$Sample.ID == sample_id,]
    group <- ifelse(grepl(sample_id, "-0"), "TUMOR", "NORMAL")
    
    # Get the relevant pData information (sample ID and status) to add
    pData_df <- data.frame("Sample_Group" = rep(group, times = 10), 
                           "Sample_ID" =rep(sample_id, times = 10))

    GRset <- minfi::readTCGA(fn, sep = "\t", keyName = "Composite Element REF",
                             Betaname = "Beta_value", pData = pData_df,
                             array = "IlluminaHumanMethylation450k", 
                             mergeManifest = FALSE, showProgress = TRUE)
    

    # Remove probes with SNPs
    GRset <- remove_probes_wSnps(GRset)
    
    # Optional function to output some potentially useful information
    print_useful_information(GRset)
    
    # Opt: Batch effects corrections
    # sva_res <- batch_effects_corr(GRset, status, pheno)
    
    # Identifying DMRs (differentially methylated regions)/ DMPs (differentially methylated positions)
    registerDoParallel(cores = 3)
    dmps <- get_DMPs(GRset)
    dmrs <- getDMRs(GRset)
    
  })
}


#' Another option: import level 1 methylation files and go through the preprocessing steps
#' 
#' 


# The path to the methylation files of interest
methfile_path <- paste(path, "Methylation_Data/Level3/", sep = "")

# Import the sample sheet (pData information)
sample_sheet <- fread(paste(methfile_path, "gdc_sample_sheet.level3_methylation.csv", sep = ""))


#' Prints some useful methylation GRanges information that we can 
#' extract using built-in minfi functions
#' @param GRset a GenomicRatioSet object
print_useful_information <- function(GRset) {
  ann <- getAnnotation(GRset)  # full annotation
  beta <- getBeta(GRset)
  #M <- getM(GRset)
  #CN <- getCN(GRSet)
  sampleNames <- sampleNames(GRset)
  probeNames <- featureNames(GRset)
  pheno < pData(GRset)
  
  print("Annotation")
  print(head(ann))
  print("Beta")
  print(head(beta))
  print("Sample Names")
  print(sampleNames)
  print("Probe Names")
  print(probeNames)
  print("Pheno")
  print(head(pheno))
}


#' Batch effects correction using built-in minfi functions
#' @param GRset a GenomicRatioSet object
#' @param status
#' @param pheno the tumor/normal phenotype
batch_effects_corr <- function(GRset, status, pheno) {
  mval <- getM(GRset)[1:5000,]
  mod <- model.matrix(~as.factor(status), data=pheno)
  mod0 <- model.matrix(~1, data=pheno)
  sva.results <- sva(mval, mod, mod0)
  return(sva.results)
}


#' Removes probes with a SNP at the CpG interrogation using minfi built-ins
#' @param GRset a GenomicRatioSet object
remove_probes_wSnps <- function(GRset) {
  snps <- getSnpInfo(GRset)
  GRset <- addSnpInfo(GRset)
  GRset <- dropLociWithSnps(GRset, snps = c("SBE", "CpG"), maf = 0)  #any minor allele frequency
  return(GRset)
}


#' Get differentially methylated positions (DMPs) using minfi
#' @param GRset_funnorm a normalized GenomicRatioSet object
get_DMPs <- function(GRset) {
  beta <- getBeta(GRset)
  cancer_or_norm <- pData(GRset)$Sample_Group
  dmp <- dmpFinder(beta, pheno = cancer_or_norm, type = "categorical",
                   qCutoff = 0.05)
  head(dmp)
  return(dmp)
}


#' Get differentially methylated regions (DMRs) using minfi
#' @param GRset_funnorm a normalized GenomicRatioSet object
get_DMRs <- function(GRset) {
  cancer_or_norm <- pData(GRset)$Sample_Group
  designMatrix <- model.matrix(~ cancer_or_norm)
  dmrs <- bumphunter(GRset, design = designMatrix, cutoff = 0.2, B = 0, type = "Beta")
  # If # of candidate bumps > 30000, increase cutoff to, say, 0.3
  # Once cutoff is found, run with permutations
  dmrs <- bumphunter(GRset, design = designMatrix, cutoff = 0.2, B = 1000, type = "Beta")
  print(names(dmrs))
  print(head(dmrs$table, n=3))
  return(dmrs)
}


#' Get the beta value of specific SNPs using minfi
#' @param GRset a GenomicRatioSet object
get_snp_beta <- function(GRset) {
  snps <- getSnpBeta(GRSet)
  print(head(snps))
  return(snps)
}

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
# ALTERNATIVE USING TCGAbiolinks:
library(ComplexHeatmap)

meth_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "DNA Methylation",
                       data.type = "Methylation Beta Value",
                       experimental.strategy = "Methylation Array",
                       workflow.type = "Liftover",
                       data.format = "txt",
                       platform = c("Illumina Human Methylation 450"),
                       access = "open")
GDCdownload(meth_query)
meth_brca <- GDCprepare(meth_query, save = FALSE)

# Remove probes with NA
meth_brca <- subset(meth_brca, subset = (rowSums(is.na(assay(meth_brca))) == 0))

# Remove probes in chromosomes X, Y, and NA
meth_brca <- subset(meth_brca, subset = !as.character(seqnames(meth_brca)) %in% 
                      c("chrNA", "chrX", "chrY"))

# Analyze differentially methylated regions (DMRs)
meth_brca_dmr <- TCGAanalyze_DMR(meth_brca.chr1,
                            groupCol = "disease_type", # a column in the colData matrix; will likely have to switch
                            group1 = "Tumor", # a type of the the disease type column
                            group2= "Normal", # a type of the the disease column
                            p.cut = 10^(-2),
                            diffmean.cut = 0.25,
                            legend = "State",
                            plot.filename = "BRCA_metvolcano.png",
                            cores = 1) # if set to 1 there will be a progress bar


# Import the clinical data
brca_clin <- GDCquery_clinic("TCGA-BRCA", "Clinical")

# Get hyper- or hypomethylated probes
signif_meth_brca <- meth_brca_dmr[values(meth_brca_dmr)[,"status.Breast.carcinoma"] %in%
                      c("Hypermethylated","Hypomethylated"),]

# Add clinical data to annotate the samples, and sort the clinical data to have the 
# same order as the DNA methylation matrix
clinical.order <- clinical[match(substr(colnames(signif_meth_brca),1,12),
                                 clinical$bcr_patient_barcode),]
ta = HeatmapAnnotation(df = clinical.order[,c("disease","gender","vital_status","race")],
                       col = list(disease = c("Tumor" = "grey", "Normal" = "black")))

                            

