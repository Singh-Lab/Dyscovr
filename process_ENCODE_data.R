############################################################
### Process ENCODE Data Files
### Written By: Sara Camilli, July 2020
############################################################

# This files explores 2 options for using ENCODE data:

# 1. Use downstream targets for TFs in our list that are found in the manually curated lists
# 2. Use raw ChIP-seq data from ENCODE to infer targets for all the transcription
# factors we are interested in that have an experiment in ENCODE2/3

############################################################
############################################################
# OPTION ONE: USE CURATED TF DATASETS

# ENCODE Transcription Factor Targets Dataset (from ChIP-seq data): https://amp.pharm.mssm.edu/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets
# TRANSFAC Curated Transcription Factor Targets Dataset (Manually curated): https://amp.pharm.mssm.edu/Harmonizome/dataset/TRANSFAC+Curated+Transcription+Factor+Targets
# TRUUST Curated Transcription Factor Targets Dataset (Using sentence-based text mining): https://www.grnpedia.org/trrust/
############################################################
############################################################

library(Matrix)
options(max.print=130000)

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)

############################################################
# IMPORT DATASETS AS MATRICES
############################################################
# ENCODE: 22819 x 182 matrix, rows are targets, columns are TFs
encode_gene_attrib_mat <- read.csv(paste(path, "Input Data Files/BRCA Data/ENCODE Gene-Attrib Matrices/gene_attrib_matrix_encode_clean.csv", sep = ""),
                                   row.names = 1)
# TRANSFAC: 13216 x 201 matrix, rows are targets, columns are TFs
transfac_gene_attrib_mat <- read.csv(paste(path, "Input Data Files/BRCA Data/ENCODE Gene-Attrib Matrices/gene_attrib_matrix_transfac_clean.csv", sep = ""),
                                     row.names = 1)
# TRRUST: 9396 x 4 matrix, column one is TFs, column two is targets, column 3 is direction of relationship
trrust_gene_target_matrix <- read.table(paste(path, "Input Data Files/trrust_rawdata.human.txt", sep = ""), sep = "\t", header = FALSE)

############################################################
# EXTRACT TFs OF INTEREST FROM SAVED CSV FILES
############################################################
# Get our list of regulatory proteins of interest from the previous files, or keep in
# the local environment -- see "process_mutation_data.R"
# Extract from CSV files:
protein_ids_iprotein <- read.table(paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_iprotein_wConCavity_0.1.csv", sep = ""), 
                          header = TRUE, sep = ",")[,2] # 6214 for BRCA
protein_ids_idomain <- read.table(paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_idomain_wConCavity_0.1.csv", sep = ""), 
                                 header = TRUE, sep = ",")[,2] # 3997 for BRCA
protein_ids_ibindingpos <- read.table(paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_ibindingpos_wConCavity_0.1.csv", sep = ""), 
                                  header = TRUE, sep = ",")[,2] # 810 for BRCA
protein_ids_iprotein_nucacids <- read.table(paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_iprotein_nucacids_wConCavity_0.1.csv", sep = ""), 
                                           header = TRUE, sep = ",")[,2] # 1914 for BRCA
protein_ids_idomain_nucacids <- read.table(paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_idomain_nucacids_wConCavity_0.1.csv", sep = ""), 
                                            header = TRUE, sep = ",")[,2] # 1190 for BRCA
protein_ids_ibindingpos_nucacids <- read.table(paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_ibindingpos_nucacids_wConCavity_0.1.csv", sep = ""), 
                                           header = TRUE, sep = ",")[,2] # 55 for BRCA

#' Quick helper function to convert a vector of SWISSPROT IDs to gene names
#' @param swissprot_ids vector of SWISSPROT IDs to be converted
convert_to_gn <- function(swissprot_ids) {
  return(unique(unlist(lapply(swissprot_ids, function(x) all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'external_gene_name']))))
}
protein_ids_iprotein <- convert_to_gn(protein_ids_iprotein)  # 6063 in BRCA
protein_ids_idomain <- convert_to_gn(protein_ids_idomain)    # 3939 in BRCA
protein_ids_ibindingpos <- convert_to_gn(protein_ids_ibindingpos)  # 843 in BRCA
protein_ids_iprotein_nucacids <- convert_to_gn(protein_ids_iprotein_nucacids)  # 1869 in BRCA 
protein_ids_idomain_nucacids <- convert_to_gn(protein_ids_idomain_nucacids)   # 1174 in BRCA
protein_ids_ibindingpos_nucacids <- convert_to_gn(protein_ids_ibindingpos_nucacids)  # 72 in BRCA

### NOTE: We now also have protein subsets for missense and silent; can also proceed 
# through following protocol with these lists instead (just change filename to append "_missense" or "_silent")


############################################################
# SUBSET MATRIX COLUMNS BY TFs OF INTEREST
############################################################
#' Given three different protein:target matrices (from ENCODE, TRANSFAC, and TRRUST),
#' subsets each by the given protein IDs of interest; returns a list of the subsetted matrices
#' and prints the intersection for each & the total intersection
#' @param encode_mat a protein-target matrix from ENCODE
#' @param transfac_mat a protein-target matrix from TRANSFAC
#' @param trrust_mat a protein-target matrix from TRRUST
#' @param protein_ids a vector of protein names of interest
subset_by_protein_ids <- function(encode_mat, transfac_mat, trrust_mat, protein_ids) {
  # ENCODE
  encode_mat_subset <- encode_mat[,colnames(encode_mat) %fin% protein_ids]
  print(paste(length(unique(colnames(encode_mat_subset))), paste("proteins of interest are found in the ENCODE matrix, or", 
                                             paste(length(unique(colnames(encode_mat_subset))), length(protein_ids), sep = "/"))))
  # TRANSFAC
  transfac_mat_subset <- transfac_mat[,colnames(transfac_mat) %fin% protein_ids]
  print(paste(length(unique(colnames(transfac_mat_subset))), paste("proteins of interest are found in the TRANSFAC matrix, or", 
                                             paste(length(unique(colnames(transfac_mat_subset))), length(protein_ids), sep = "/"))))
  
  # TRRUST
  trrust_mat_subset <- trrust_mat[trrust_mat[,1] %fin% protein_ids,]
  print(paste(length(unique(trrust_mat_subset[,1])), paste("proteins of interest are found in the TRRUST matrix, or", 
                                               paste(length(unique(trrust_mat_subset[,1])), length(protein_ids), sep = "/"))))

  total_union <- unique(c(colnames(encode_mat_subset), colnames(transfac_mat_subset), trrust_mat_subset[,1]))
  print(paste(length(total_union), "total proteins of interest are found in one of the three curated datasets."))
  
  # Return the subsetted dataframes
  return(list("ENCODE" = encode_mat_subset, "TRANSFAC" = transfac_mat_subset, "TRRUST" = trrust_mat_subset))
}
 
iprotein_subset_dfs <- subset_by_protein_ids(encode_gene_attrib_mat, transfac_gene_attrib_mat, 
                                             trrust_gene_target_matrix, protein_ids_iprotein)
  # [1] "86 proteins of interest are found in the ENCODE matrix, or 86/6063"
  # [1] "81 proteins of interest are found in the TRANSFAC matrix, or 81/6063"
  # [1] "332 proteins of interest are found in the TRRUST matrix, or 332/6063"
  # [1] "350 total proteins of interest are found in one of the three curated datasets."
idomain_subset_dfs <- subset_by_protein_ids(encode_gene_attrib_mat, transfac_gene_attrib_mat, 
                                             trrust_gene_target_matrix, protein_ids_idomain)
  # [1] "53 proteins of interest are found in the ENCODE matrix, or 53/3939"
  # [1] "47 proteins of interest are found in the TRANSFAC matrix, or 47/3939"
  # [1] "192 proteins of interest are found in the TRRUST matrix, or 192/3939"
  # [1] "207 total proteins of interest are found in one of the three curated datasets."
ibindingpos_subset_dfs <- subset_by_protein_ids(encode_gene_attrib_mat, transfac_gene_attrib_mat, 
                                             trrust_gene_target_matrix, protein_ids_ibindingpos)
  # [1] "16 proteins of interest are found in the ENCODE matrix, or 16/843"
  # [1] "18 proteins of interest are found in the TRANSFAC matrix, or 18/843"
  # [1] "64 proteins of interest are found in the TRRUST matrix, or 64/843"
  # [1] "66 total proteins of interest are found in one of the three curated datasets."
iprotein_nucacids_subset_dfs <- subset_by_protein_ids(encode_gene_attrib_mat, transfac_gene_attrib_mat, 
                                             trrust_gene_target_matrix, protein_ids_iprotein_nucacids)
  # [1] "64 proteins of interest are found in the ENCODE matrix, or 64/1869"
  # [1] "78 proteins of interest are found in the TRANSFAC matrix, or 78/1869"
  # [1] "236 proteins of interest are found in the TRRUST matrix, or 236/1869"
  # [1] "253 total proteins of interest are found in one of the three curated datasets."
idomain_nucacids_subset_dfs <- subset_by_protein_ids(encode_gene_attrib_mat, transfac_gene_attrib_mat, 
                                             trrust_gene_target_matrix, protein_ids_idomain_nucacids)
  # [1] "42 proteins of interest are found in the ENCODE matrix, or 42/1174"
  # [1] "50 proteins of interest are found in the TRANSFAC matrix, or 50/1174"
  # [1] "139 proteins of interest are found in the TRRUST matrix, or 139/1174"
  # [1] "148 total proteins of interest are found in one of the three curated datasets."
ibindingpos_nucacids_subset_dfs <- subset_by_protein_ids(encode_gene_attrib_mat, transfac_gene_attrib_mat, 
                                             trrust_gene_target_matrix, protein_ids_ibindingpos_nucacids)
  # [1] "5 proteins of interest are found in the ENCODE matrix, or 5/72"
  # [1] "9 proteins of interest are found in the TRANSFAC matrix, or 9/72"
  # [1] "21 proteins of interest are found in the TRRUST matrix, or 21/72"
  # [1] "21 total proteins of interest are found in one of the three curated datasets."

# Now that we have potential downstream targets for some of our proteins, we can make lists for each of 
# these with their potential targets; what is a good data structure for this?


############################################################
# CREATE LISTS OF TARGETS FOR EACH TF OF INTEREST WITH DATA
############################################################
#' Takes a list of protein-target matrices and returns a list of lists 
#' for each column (combined lists of targets for each transcription factor)
#' @param subset_matrices a list of subsetted protein-target matrices, the length
#' is equal to the number of sources (3)
#' @param protein_ids a vector of names of the proteins of interest
get_target_lists <- function(subset_matrices, protein_ids) {
  encode_mat <- subset_matrices[[1]]
  transfac_mat <- subset_matrices[[2]]
  trrust_mat <- subset_matrices[[3]]
  
  target_gene_lists <- lapply(protein_ids, function(id) {
    print(id)
    
    # If we have targets for this regulatory protein ID... extract and compile them
    if(id %fin% colnames(encode_mat) | id %fin% colnames(transfac_mat) | id %fin% trrust_mat[,1]) {
      if(id %fin% colnames(encode_mat)) {
        encode_sub <- encode_mat[,colnames(encode_mat) == id]
        encode_targs <- rownames(encode_mat[which(encode_sub != 0),colnames(encode_mat) == id])
      } else {encode_targs <- c()}
      
      if(id %fin% colnames(transfac_mat)) {
        transfac_sub <- transfac_mat[,colnames(transfac_mat) == id]
        transfac_targs <- rownames(transfac_mat[which(transfac_sub != 0), colnames(transfac_mat) == id])
      } else {transfac_targs <- c()}
      
      if(id %fin% trrust_mat[,1]) {
        trrust_targs <- trrust_mat[trrust_mat[,1] == id, 2]
      } else {trrust_targs <- c()}
      
      all_targs <- unique(c(encode_targs, transfac_targs, trrust_targs))
      if(length(all_targs) > 0) {all_targs}
      else {return(NA)}
    }
    else {return(NA)}
  })
  names(target_gene_lists) <- protein_ids
  target_gene_lists <- target_gene_lists[!is.na(target_gene_lists)]
  return(target_gene_lists)
}


iprotein_targets <- get_target_lists(iprotein_subset_dfs, protein_ids_iprotein)
idomain_targets <- get_target_lists(idomain_subset_dfs, protein_ids_idomain)
ibindingpos_targets <- get_target_lists(ibindingpos_subset_dfs, protein_ids_ibindingpos)
iprotein_targets_nucacids <- get_target_lists(iprotein_nucacids_subset_dfs, protein_ids_iprotein_nucacids)
idomain_targets_nucacids <- get_target_lists(idomain_nucacids_subset_dfs, protein_ids_idomain_nucacids)
ibindingpos_targets_nucacids <- get_target_lists(ibindingpos_nucacids_subset_dfs, protein_ids_ibindingpos_nucacids)


############################################################
# SAVE LISTS TO TEXT FILES
############################################################
# Save these lists to text files that we can parse later
output_path <- paste(path, "Saved Output Data Files/BRCA/Curated TF Data/", sep = "")

#' Takes a regulatory protein-gene targets list and sinks it to a text
#' file using the sink() function
#' @param targets the list of regulatory proteins and their targets
#' @param label the specificity we're looking at ("iprotein", "idomain", etc.)
#' @param path the path we want to write the new text file to
sink_to_text_file <- function(targets, label, path) {
  sink(paste(path, paste(label, "_gene_targets_full.txt", sep = ""), sep = ""))
  print(targets)
  sink()
}

sink_to_text_file(iprotein_targets, "iprotein", output_path)
sink_to_text_file(idomain_targets, "idomain", output_path)
sink_to_text_file(ibindingpos_targets, "ibindingpos", output_path)
sink_to_text_file(iprotein_targets_nucacids, "iprotein_nucacids", output_path)
sink_to_text_file(idomain_targets_nucacids, "idomain_nucacids", output_path)
sink_to_text_file(ibindingpos_targets_nucacids, "ibindingpos_nucacids", output_path)


############################################################
# CONVERT TO UNIFIED DATAFRAME
############################################################
#' Takes a list of the form list(TF name = vector of curated targets),
#' and fills a DF with protein's gene targets (TFs are columns, elements 
#' are targets)
#' @param tf_list the transcription factor-target gene list
compile_targs <- function(tf_list) {
  tot_length <- max(lengths(tf_list))
  targs_df <- data.frame(matrix(nrow = tot_length, ncol = length(tf_list)))
  colnames(targs_df) <- names(tf_list)
  
  for (i in 1:length(tf_list)) {
    # Get the current protein and its list of targets
    tf <- names(tf_list)[i]
    targets <- tf_list[[i]]
    
    # Pad with NAs
    targets <- add_na_vect(targets, tot_length)
    targs_df[,tf] <- targets
  }
  # Return the filled DF
  return(targs_df)
}

iprotein_targs_df <- compile_targs(iprotein_targets)
idomain_targs_df <- compile_targs(idomain_targets)
ibindingpos_targs_df <- compile_targs(ibindingpos_targets)
iprotein_nucacids_targs_df <- compile_targs(iprotein_targets_nucacids)
idomain_nucacids_targs_df <- compile_targs(idomain_targets_nucacids)
ibindingpos_nucacids_targs_df <- compile_targs(ibindingpos_targets_nucacids)

# Write this to a CSV
write.csv(iprotein_targs_df, paste(output_path, "iprotein_curated_targets_df.csv", sep = ""))
write.csv(idomain_targs_df, paste(output_path, "idomain_curated_targets_df.csv", sep = ""))
write.csv(ibindingpos_targs_df, paste(output_path, "ibindingpos_curated_targets_df.csv", sep = ""))
write.csv(iprotein_nucacids_targs_df, paste(output_path, "iprotein_nucacids_curated_targets_df.csv", sep = ""))
write.csv(idomain_nucacids_targs_df, paste(output_path, "idomain_nucacids_curated_targets_df.csv", sep = ""))
write.csv(ibindingpos_nucacids_targs_df, paste(output_path, "ibindingpos_nucacids_curated_targets_df.csv", sep = ""))



#
#
#
#
#
#
############################################################
############################################################
# OPTION TWO: USE RAW ChIP-seq DATA FROM ENCODE2/3
############################################################
############################################################
library(genomation)
library(biomaRt)
library(ChIPpeakAnno)

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/ENCODE ChIP-Seq Files/IDR-Ranked/"

############################################################
# About ENCODE data formats: https://genome.ucsc.edu/FAQ/FAQformat.html#format12

# Genomation Vignette: http://www.bioconductor.org/packages/release/bioc/manuals/genomation/man/genomation.pdf
# NOTE: The IDR ranked files have a different format than normal narrow peak files! Must import & process by hand

# IDR (Irreproducibility Discovery Rate) is a mechanism to assess concordance of peak calls between replicates
# Tutorial for running IDR: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html

# Link to UCSC Table Browser: https://genome.ucsc.edu/cgi-bin/hgTables
# Information about UCSC Genome Browser Table format: https://genome.ucsc.edu/FAQ/FAQformat.html
############################################################


############################################################
### OPTION 1: ChIPpeakAnno PACKAGE ###
############################################################

# Import the BED files
bed_filenames <- list.files(path = path, full.names = TRUE, recursive = FALSE, pattern = ".bed")  # Specify a local path 

# Import the metadata file (links bed file to name of gene in experiment)
metadata <- read.csv(paste(path, "labels.csv", sep = ""), header = TRUE)

# Convert to GRanges & add metadata
gRanges_list <- lapply(1:length(bed_filenames), function(i) toGRanges(bed_filenames[i], format = "BED", header = FALSE))


############################################################
### OPTION 2: MANUAL + BIOMART ###
############################################################
# Import BED files
# bed_filenames <- commandArgs(trailingOnly = TRUE) # Use a list of arguments OR
bed_filenames <- list.files(path = path, full.names = TRUE, recursive = FALSE, pattern = ".bed")  # Specify a local path 

# Import the metadata file (links bed file to name of gene in experiment)
metadata <- read.csv(paste(path, "labels.csv", sep = ""), header = TRUE)

# The headers for normal narrow peaks files, as well as for IDR-ranked files
narrowpeaks_col_headers <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
                             "signalValue", "pValue", "qValue", "summit")
idr_ranked_col_headers <- c(narrowpeaks_col_headers, "globalIDR.-log10", "localID.-log10", 
                            "rep1_chromStart", "rep1_chromEnd", "rep1_signalValue",
                            "rep1_summit", "rep2_chromStart", "rep2_chromEnd", 
                            "rep2_signalValue", "rep2_summit")  
# Note that 'score' in this case is the scaled IDR value
# 'globalIDR' is the value used to calculate the scaled IDR value (above),
  # analogous to a multiple hypothesis correction on a p-value to compute an FDR
# 'localIDR' is akin to the posterior probability of the peak belonging to the
  # irreproducible noise component

# Set the q-value threshold (0.05 tends to be most appropriate for human/ miuse with large genomes)
qval_thres <- 0.05

# Uncomment if using normal narrow peak files
# This function takes a list of bed filenames (normal narrow peak) and reads the
# file to a table
#read_files <- function(bed_filenames) {
#  for (i in 1:length(bed_filenames)) {
#    filename <- bed_filenames[i]
    #bed_track <- readNarrowPeak(file = filename, track.line = FALSE)
#  }
#}

# Create list of the form list(TF1 = ChIP-seq DF, TF2 = ChIP-seq DF, etc.)
protein_chipseq_df_list <- create_prot_chipseq_df_list(bed_filenames, idr_ranked_col_headers, 
                                                       qval_thres, metadata, protein_ids)

# Import file that links gene IDs to position
input_file_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/ID Conversions/"
# geneid_pos_conv <- read.csv(paste(input_file_path, "gene_chrom_num_and_pos_hg17.csv", sep = ""))
#header = TRUE)  # NOTE: THIS IS FOR hg17 (OLD VERSION)
geneid_pos_conv <- read.csv(paste(input_file_path, "gene_chrom_num_and_pos_hg38.csv", sep = ""), header = TRUE)
colnames(geneid_pos_conv)[1] <- "name"

# Import file that links ENST to ENSG 
# enst_ensg_conv <- read.csv(paste(input_file_path, "ensg_enst_conv.csv", sep = ""), header = TRUE)
# Create a mart for ENST to ENSG ID conversions
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Loop through the list of ChIP-Seq dataframes; for each row in each the dataframe,
# call helper function that will identify matching gene ID and add it to the dataframe.
for (i in 1:length(protein_chipseq_df_list)) {
  dataframe <- protein_chipseq_df_list[[i]]
  
  # Update the DF with gene info and add it back into its spot in the list
  protein_chipseq_df_list[[i]] <- update_with_gene_info(dataframe, geneid_pos_conv, mart)
}

# Write this to a file
output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/ChIP-Seq/"
sink(paste(output_path, "chipseq_dataframe_list.txt", sep = ""))
print(protein_chipseq_df_list)
sink()

# Now, for each gene, we have a table updated with their potential targets. Extract these
# in list form and write to a composite table
  # Columns: Genes studied in a ChIP-seq experiment
  # Entries: List of targets
unique_prots <- c()
for (i in 1:length(names(protein_chipseq_df_list))) {
  name_full <- names(protein_chipseq_df_list)[i]
  gene_name <- unlist(strsplit(as.character(name_full), "_", fixed = TRUE))[1]
  if (!(gene_name %in% unique_prots)) {
    unique_prots <- c(unique_prots, gene_name)
  }
}
# Get the max number of targets
#max_num_targets <- 0
#for (i in 1:length(names(protein_chipseq_df_list))) {
#  ensembl_ids <- unlist(protein_chipseq_df_list[[i]]$ensembl_id)
#  ensembl_ids <- lapply()
#}
num_rows <- 100000
chipseq_targs_df <- data.frame(matrix(nrow = num_rows, ncol = length(unique_prots)))   # May have to adjust nrow
colnames(chipseq_targs_df) <- unique_prots
chipseq_targs_df <- compile_chipseq_targs(protein_chipseq_df_list, chipseq_targs_df, num_rows)

# Write this to a CSV
write.csv(chipseq_targs_df, paste(output_path, "chipseq_targets_df.csv", sep = ""))

############################################################
### HELPER FUNCTIONS ###
############################################################

############################################################
# Takes the names of the bed files, the headers for IDR-ranked files, a q-value
# threshold, a metadata file, and protein IDs of interest. 
# Uses a helper function to create and output a list of the 
# form list(TF1 = ChIP-seq DF, TF2 = ChIP-seq DF, etc.) that only
# includes entries with a q-value above the given threshold
# and for proteins of interest
############################################################
create_prot_chipseq_df_list <- function(bed_filenames, idr_ranked_col_headers, qval_thres, metadata, protein_ids) {
  protein_chipseq_df_list <- c()
  
  # Loop through all bedfiles subset ChIP-seq DF by q-value threshold. Extract name
  # of protein and add to list
  for (i in 1:length(bed_filenames)) {
    filename <- bed_filenames[i]
    bed_file_df <- read_bed_file(filename, idr_ranked_col_headers, qval_thres)
    
    # Get the experiment name
    experiment_name <- unlist(strsplit(unlist(strsplit(filename, "/", fixed = TRUE))[12], ".", fixed = TRUE))[1]  # Adjust first index based on length of path
    
    # Use the experiment name to get the protein name
    protein_name <- metadata[metadata$Experiement.Name == experiment_name, 1][1]
    
    # Use the protein name_experiment name as the identifier in the list
    coverage_counter <- 0
    if (protein_name %in% protein_ids) {
      identifier <- paste(protein_name, experiment_name, sep = "_")
      protein_chipseq_df_list[[identifier]] <- bed_file_df
      coverage_counter <- coverage_counter + 1
    }
  }
  print(paste("ChIP-seq IDR-Ranked Coverage:", coverage_counter/length(protein_ids)))
  return(protein_chipseq_df_list)
}

############################################################
# Takes in a bed filename (IDR ranked), along with a list 
# of headers, and reads in the file to a table with appropriate headers. Then,
# filters this table to retain only q-values greater than or equal to the given
# q-value threshold and returns this table.
############################################################
read_bed_file <- function(filename, idr_ranked_col_headers, qval_thres) {
  # Import the bed file as a table and assign headers
  bed_file <- read.table(filename)
  colnames(bed_file) <- idr_ranked_col_headers
  
  # Extract the q-values and convert to normal q-values (not -log10)
  qvalues <- bed_file$qValue
  qvals_untrans <- sapply(qvalues, FUN = function(x) {-(10^x)})
  
  # Apply threshold and subset bed table
  indices_to_keep <- which(qvals_untrans < qval_thres)
  bed_file_subset <- bed_file[indices_to_keep,]
  
  # Return the subsetted table
  return(bed_file_subset)
}

############################################################
# Takes a ChIP-seq IDR-ranked dataframe, along with a conversion document from 
# BioMart that has Ensembl ID, Gene Name, Chromosome #, ChromStart, and ChromEnd.
# Loops through all rows in the IDR-ranked dataframe, and for each peak (with a
# chromosome start and end position) finds a matching gene hit in the conversion
# document. Annotates the dataframe with the Ensembl ID and gene name of the hit,
# and returns the updated dataframe.
############################################################
update_with_gene_info <- function(df, geneid_pos_conv, mart) {
  df_new <- df
  for (i in 1:nrow(df)) {
    # Get the chromosome and start/ end positions
    chrom <- df$chrom[i]
    chromStart <- df$chromStart[i]
    chromEnd <- df$chromEnd[i]
    
    # Subset the conversion doc to only include this chromosome
    geneid_pos_conv_sub <- geneid_pos_conv[geneid_pos_conv$chrom == chrom,]

    # Subset the conversion doc to only include genes that overlap this range
      # OPT 1: Both start and end are in range (gene is encapsulated by this binding range),
      # OPT 2: The start of this gene is in range
      # OPT 3: The end of this gene is in range
      # OPT 4: The gene encapsulates this binding range
    # NOTE: I can define the gene as either the full gene (txStart and txEnd) or the coding
      # region (cdsStart, cdsEnd) (also exon starts and ends)
    geneid_pos_conv_sub <- geneid_pos_conv_sub[((geneid_pos_conv_sub$txStart >= chromStart) & (geneid_pos_conv_sub$txEnd <= chromEnd)) |  # OPT 1
                                                 ((geneid_pos_conv_sub$txStart >= chromStart) & (geneid_pos_conv_sub$txStart <= chromEnd)) |  # OPT 2
                                                 ((geneid_pos_conv_sub$txEnd >= chromStart) & (geneid_pos_conv_sub$txEnd <= chromEnd)) |  # OPT 3
                                                 ((geneid_pos_conv_sub$txStart <= chromStart) & (geneid_pos_conv_sub$txEnd >= chromEnd)),] # OPT 4

    # Add the Ensembl ID(s) and the Swissprot ID(s) to the original dataframe
    if (nrow(geneid_pos_conv_sub) == 0) {
      next
    }
    else {
      # Extract the Ensembl IDs
      ensembl_ids <- as.character(c(geneid_pos_conv_sub$name))
      # Remove the part after the "."
      #ensembl_ids_fixed <- unlist(lapply(1:length(ensembl_ids), function(i) unlist(strsplit(ensembl_ids[i], ".", fixed = TRUE))[1]))
      ensembl_ids_fixed <- getBM(attributes = 'ensembl_gene_id', filters = 'ensembl_transcript_id_version',
                                 values = ensembl_ids, mart = mart)
      #ensembl_gene_ids <- unlist(lapply(1:length(ensembl_ids_fixed), 
                                        #function(i) enst_ensg_conv[enst_ensg_conv$Transcript.stable.ID == ensembl_ids_fixed[i],1]))

      # Sub NA for the Swissprot IDs if there are none listed
      swissprot_ids <- as.character(c(geneid_pos_conv_sub$Gene.name))
      if (length(swissprot_ids) == 0) {swissprot_ids <- NA}
      
      # Add the fixed Ensembl IDs and Swissprot IDs to dataframe
      df_new[i, "ensembl_id"] <- paste(ensembl_ids_fixed, sep = "", collapse = ";")
      df_new[i, "swissprot_id"] <- paste(swissprot_ids, sep = "", collapse = ";")
    }
  }
  # Return the updated dataframe
  return(df_new)
}

############################################################
# Function takes a list of the form list(protein name_experiment name = IDR-Ranked DF),
# along with an empty output dataframe with the columns labeled by the protein name,
# and fills in the DF elements with protein's gene targets
############################################################
compile_chipseq_targs <- function(protein_chipseq_df_list, chipseq_targs_df, num_rows) {
  
  for (i in 1:length(protein_chipseq_df_list)) {
    # Get the current protein and its dataframe
    prot_curr <- unlist(strsplit(names(protein_chipseq_df_list)[i], "_"))[1]
    print(prot_curr)
    df_curr <- protein_chipseq_df_list[[i]]
    
    # Extract the Ensembl IDs 
    ensembl_ids <- unlist(df_curr$ensembl_id)
    # Remove all NA values
    ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]

    # For each Ensembl ID box, get the individual IDs and add to the DF if not 
    # already in there
    unique_ids <- c()
    for (i in 1:length(ensembl_ids)) {
      # Get the new targets from this row
      new_ids <- ensembl_ids[i]
      
      # If there are no targets in this row, skip to next iteration
      if (is.na(new_ids)) {
        next
      }
      
      # Otherwise, split the IDs by semicolons
      new_id_list <- unlist(strsplit(new_ids, split = ";", fixed = TRUE))
      
      # Use unique function to get just gene IDs we haven't seen yet
      unique_ids <- unique(c(unique_ids, new_id_list))
      print(paste("Unique_ids:", unique_ids))
      
    }
    print(length(unique_ids))
    # Pad with NAs
    unique_ids <- add_na_vect(unique_ids, num_rows)
    # Add the list of unique targets for this protein to the dataframe
    chipseq_targs_df[,prot_curr] <- unique_ids
  }
  # Return the filled DF
  return(chipseq_targs_df)
}


############################################################
# Takes a list of IDs, along with the total number of rows in the 
# table, and pads the ID list with NAs to the required length
# then returns
############################################################
add_na_vect <- function(id_list, tot_len) {
  na_vect <- rep_len(NA, tot_len - length(id_list))
  id_list_full <- c(id_list, na_vect)
  return(id_list_full)
}
