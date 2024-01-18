############################################################
### Process Gene-Level Copy Number Alteration (CNA) Data 
### ASSOCIATED PUBLICATION INFORMATION
############################################################

library(data.table)
library(dplyr)

# Local path to directory containing ASCAT2-generated CNA data files
PATH <- getwd()
OUTPUT_PATH <- paste0(getwd(), "Output/CNA/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
# IMPORT COPY NUMBER VALUE FILES & RELATED FILES
############################################################
# Import the per-patient CNA files for particular project
ascat2_filenames <- list.files(paste0(PATH, "Input/CNA/"), 
                               pattern = ".gene_level_copy_number")

# Import the list of unique patient IDs that have all data types 
patient_id_list <- read.table(paste0(PATH, "intersecting_ids.txt"), 
                              header = T)[,1]

# Import Biospecimen information to do the sample-patient ID conversions
biosp_df <- read.csv(paste0(PATH, "aliquot.csv"), header = T, na.strings = "\'--")

# Import a list of all genes we are interested in 
allgene_targets <- read.csv(paste0(PATH, "allgene_targets.csv"), header = T, 
                                   row.names = 1, check.names = F)


############################################################
# MAKE PER-GENE OUTPUT DATA FRAME
############################################################
#' Takes all the patient ASCAT2, raw CNA value file names and, for each patient
#' that is in the patient ID list, adds the copy number per gene (row), per 
#' patient (column) to an output data frame that is returned
#' @param PATH the local path to the directory of ascat2 files
#' @param ascat2_filenames vector of file names in the PATH's directory
#' @param patient_id_list a vector of patient 4-letter/number TCGA IDs with all 
#' data types of interest
#' @param biosp_df the biospecimen information for sample-patient ID conversions
#' @param allgene_targets the genes we are specifically interested in profiling
make_output_cna_df <- function(PATH, ascat2_filenames, patient_id_list, 
                                      biosp_df, allgene_targets) {

  # Loop through all patients, to get a column of CNA values for each gene for 
  # each patient
  cna_cols <- lapply(1:length(ascat2_filenames), function(i) {
    filename <- ascat2_filenames[i]
    print(paste(i, paste("/", length(ascat2_filenames))))
    
    # Is this patient in the patient ID list?
    uuid_aliquot <- unlist(strsplit(filename, ".", fixed = T))[2]
    samp_submitter_id <- biosp_df[biosp_df$aliquot_id == uuid_aliquot, 
                                  'sample_submitter_id']
    patient_id <- unlist(strsplit(samp_submitter_id, "-", fixed = T))[3]

    if(length(patient_id) == 0) {return(NA)}
    else if (patient_id %fin% patient_id_list) {
      # Upload this CNA file
      cna_file <- fread(paste0(PATH, paste0("Input/CNA/", filename)), 
                        sep = "\t", header = T)
      gene_ids <- unlist(lapply(cna_file$gene_id, function(x) 
        unlist(strsplit(x, ".", fixed = T))[1]))
      cna_file <- cna_file[which(gene_ids %fin% allgene_targets$ensg), ]
      
      # Get the copy number info for this
      copy_num_col <- cna_file$copy_number
      
      # Prepend the patient ID/sample id & return
      full_patient_id <- paste(unlist(strsplit(samp_submitter_id, "-", 
                                               fixed = T))[3:4], collapse = "-")
      
      copy_num_df <- as.data.frame(copy_num_col)
      colnames(copy_num_df) <- full_patient_id
      rownames(copy_num_df) <- cna_file$gene_id

      return(copy_num_df)
    } 
    else {return(NA)}
  })
  cna_cols <- cna_cols[!is.na(cna_cols)]
  
  # Use helper function to condense multiple isoforms
  output_cna_df <- address_isoforms(cna_cols)

  return(output_cna_df)
}

# Run the function
output_cna_df <- make_output_cna_df(PATH, ascat2_filenames, 
                                           patient_id_list, biosp_df,
                                           allgene_targets)

# Write the resulting data frame to a CSV file
write.csv(output_cna_df, paste0(OUTPUT_PATH, "CNA/CNA_DF_AllGenes.csv"))


############################################################
# ADDRESS CNA ISOFORMS
############################################################
#' Helper function that combines CNA per-patient columns into a data frame. 
#' In the case where there are multiple isoforms of a gene, this function
#' takes... CHECK THIS
#' @param cna_cols a list the length of the number of tumor samples, with each
#' entry containing a single-column data frame of per-gene CNA values
address_isoforms <- function(cna_cols) {
  num_rows_per_entry <- unique(unlist(lapply(cna_cols, function(df) nrow(df))))
  output_cna_df_full <- NA
  
  # If the number of rows per entry is not equivalent across all entries
  if(length(num_rows_per_entry) > 1) {
    
    output_cna_dfs_perLength <- lapply(num_rows_per_entry, function(l) {
      to_keep_l <- unlist(lapply(1:length(cna_cols), function(i) 
        if(nrow(cna_cols[[i]]) == l) {return(i)}))
      output_cna_df_l <- cna_cols[to_keep_l]
      output_cna_df_l_full <- do.call(cbind, cna_cols)
      output_cna_df_l_full$ensg_id <- rownames(output_cna_df_l_full)
      return(output_cna_df_l_full)
    })
    output_cna_df_full <- Reduce(function(df1, df2) 
      merge(df1, df2, by = 'ensg_id', all = T), output_cna_dfs_perLength)
    rownames(output_cna_df_full) <- output_cna_df_full$ensg_id
    output_cna_df_full <- output_cna_df_full[,!(colnames(output_cna_df_full) 
                                                %fin% 'ensg_id')]
  
  # Otherwise, simply combine together into a single data frame
  } else {
    output_cna_df_full <- do.call(cbind, cna_cols)
    rownames(output_cna_df_full) <- output_cna_df_full$ensg_id
    output_cna_df_full <- output_cna_df_full[,!(colnames(output_cna_df_full) %fin% 'ensg_id')]
  }
  
  return(output_cna_df_full)
}
