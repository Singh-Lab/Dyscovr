############################################################
### Merge Gene Target Data
### Written By: Sara Camilli, August 2020
############################################################
library(TCGAbiolinks)
library(stringr)

# Takes information about regulatory protein:gene target pairings
# from various sources (curated lists, ChIP-seq, ATAC-Seq) and
# merges them into one "downstream target data frame", along with
# proteins of interest that we have no information for (list all
# genes in the human genome as potential targets)
# This data frame will have the following format:
  # Columns: regulatory proteins (SWISSPROT IDs)
  # Entries: target gene names

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"

############################################################
# GET ALL GENES IN THE HUMAN GENOME (POTENTIAL DOWNSTREAM TARGETS)
############################################################
# NOTE: we will only use this set of potential downstream targets for genes that do not have ENCODE data
genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg38") 
genes <- genes[order(genes$start_position),]     # order genes by start position
# Keep only columns of interest
genes <- genes[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position","end_position")]
colnames(genes) <- c("Ensembl.ID", "Gene.Symbol","Chr","Start","End")  # rename columns
# genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
# save(genes_GR,genes,file = "genes_GR.rda", compress = "xz")


############################################################
# GET ALL REGULATORY PROTEINS IN CONSIDERATION
############################################################
prot_path <- paste(main_path, "Mutation/", sep = "")
iprotein_regprots <- read.csv(paste(prot_path, "swissprot_ids_missense_iprotein.csv", sep = ""), header = TRUE, row.names = 1)[,1]
iprotein_nucacids_regprots <- read.csv(paste(prot_path, "swissprot_ids_missense_iprotein_nucacids.csv", sep = ""), header = TRUE, row.names = 1)[,1]
idomain_regprots <- read.csv(paste(prot_path, "swissprot_ids_missense_idomain.csv", sep = ""), header = TRUE, row.names = 1)[,1]
idomain_nucacids_regprots <- read.csv(paste(prot_path, "swissprot_ids_missense_idomain_nucacids.csv", sep = ""), header = TRUE, row.names = 1)[,1]
ibindingpos_regprots <- read.csv(paste(prot_path, "swissprot_ids_missense_ibindingpos.csv", sep = ""), header = TRUE, row.names = 1)[,1]
ibindingpos_nucacids_regprots <- read.csv(paste(prot_path, "swissprot_ids_missense_ibindingpos_nucacids.csv", sep = ""), header = TRUE, row.names = 1)[,1]


############################################################
# RETRIVE THE CURATED TARGETS FOR THESE PROTEINS 
# FROM ENCODE, TRANSFAC, ETC.
############################################################
iprotein_curated_targs <- scan(paste(main_path, "Curated TF Data/iprotein_gene_targets_full.txt", sep = ""), what = "", sep = "\n")
iprotein_nucacids_curated_targs <- scan(paste(main_path, "Curated TF Data/iprotein_nucacids_gene_targets_full.txt", sep = ""), what = "", sep = "\n")
idomain_curated_targs <- scan(paste(main_path, "Curated TF Data/idomain_gene_targets_full.txt", sep = ""), what = "", sep = "\n")
idomain_nucacids_curated_targs <- scan(paste(main_path, "Curated TF Data/idomain_nucacids_gene_targets_full.txt", sep = ""), what = "", sep = "\n")
ibindingpos_curated_targs <- scan(paste(main_path, "Curated TF Data/ibindingpos_gene_targets_full.txt", sep = ""), what = "", sep = "\n")
ibindingpos_nucacids_curated_targs <- scan(paste(main_path, "Curated TF Data/ibindingpos_nucacids_gene_targets_full.txt", sep = ""), what = "", sep = "\n")

#' Process the imported text file into a proper list
#' @param curated_targs a vector containing all the regulatory proteins and their targets,
#' to be converted into list form
convert_to_list <- function(curated_targs) {
  # Extract the names and the entries separately
  names <- curated_targs[which(unlist(lapply(curated_targs, function(x) startsWith(x, "$"))))]
  
  # Compile the entries for each name into a vector
  entries <- unlist(lapply(names, function(name) {
    index <- which(iprotein_curated_targs == name)
    print(index)
    increment <- 1
    next_val <- iprotein_curated_targs[index+increment]
    print(next_val)
    targs <- c()
    while(!startsWith(next_val, '$')) {
      targs <- c(targs, next_val)
      increment <- increment + 1
      next_val <- iprotein_curated_targs[index_name + increment]
    }
    return(targs)
  }))
  
  # Remove the excess '$' characters from names
  names_clean <- unlist(lapply(names, function(x) sub("$", "", x, fixed = TRUE)))
  
  # Remove excess characters and line labels from entries
  entries_clean <- unlist(lapply(entries, function(x) {
    x_new <- str_squish(gsub("[^[:alnum:] ]", "", x))
    x_split <- unlist(strsplit(x_new, " ", fixed = TRUE))
    # Remove any entries that are numbers only
    x_final <- paste(x_split[grepl("\\D", x_split)], collapse = " ")
    return(x_final)
  }))
  
  curated_targs_list <- as.list(entries_clean)
  names(curated_targs_list) <- names_clean
  
  return(curated_targs_list)
}

# Call this function to convert all to lists
iprotein_curated_targs <- convert_to_list(iprotein_curated_targs)
iprotein_nucacids_curated_targs <- convert_to_list(iprotein_ncuacids_curated_targs)
idomain_curated_targs <- convert_to_list(idomain_curated_targs)
idomain_nucacids_curated_targs <- convert_to_list(idomain_nucacids_curated_targs)
ibindingpos_curated_targs <- convert_to_list(ibindingpos_curated_targs)
ibindingpos_nucacids_curated_targs <- convert_to_list(ibindingpos_nucacids_curated_targs)


############################################################
# RETRIVE THE CHIP-EAT TARGETS FOR THESE PROTEINS, IF AVAILABLE
############################################################


############################################################
# RETRIVE THE CHIP-SEQ TARGETS FOR THESE PROTEINS 
############################################################
chipseq_targs_df <- read.csv(paste(main_path, "ChIP-Seq/chipseq_targets_df.csv", sep = ""), header = TRUE)


#TODO: the following two functions will likely need some adjusting as I figure
# out how to format the ChIP-eat and ChIP-seq files for this

############################################################
# MERGE THE DATA INTO ONE UNIFIED DATAFRAME
############################################################
#' Eliminate entries that aren't in the binding proteins list
#' @param targets_df a data frame of regulatory proteins (columns) and their 
#' gene targets (entries), compiled from various sources
#' @param binding_prots a vector of regulatory protein IDs for proteins
#' that bind ligands (of interest)
eliminate_nonoverlapping_prots <- function(targets_df, binding_prots) {
  targets_df_sub <- targets_df[,colnames(targets_df) %fin% binding_prots]
  return(targets_df_sub)
}
curated_gene_targets_df_sub <- eliminate_nonoverlapping_prots(curated_gene_targets_df, binding_prots_w_muts_missense)
#curated_gene_targets_df_sub <- eliminate_nonoverlapping_prots(curated_gene_targets_df, binding_prots_w_muts_in_interact_regs_missense)

chipseq_targs_df_sub <- eliminate_nonoverlapping_prots(chipseq_targs_df, binding_prots_w_muts_missense)
#chipseq_targs_df_sub <- eliminate_nonoverlapping_prots(chipseq_targs_df, binding_prots_w_muts_in_interact_regs_missense)

# Takes the DFs from ChIP-Seq and the curated list and combine the proteins and their targets
gene_targets_df <- merge(curated_gene_targets_df_sub, chipseq_targs_df_sub)  # TODO: specify by.y if merge doesn't perform properly by default


############################################################
# USE FULL GENE LIST AND ENCODE DATA TO GET A SET OF TARGETS
# FOR EACH TRANSCRIPTION FACTOR
############################################################
#' Creates a downstream target data frame, where all proteins 
#' of interest (including those without ChIP-seq experiments/ known targets)
#' have a list of targets to test. Those without any ChIP-seq information
#' will be given the full genome as their list of potential targets.
#' Columns of output DF are proteins and entries are targets
#' @param protein_ids a vector of regulatory proteins of interest
#' @param genes a vector of all genes in the human genome
#' @param genes_targets_df a composite file of specific, known 
#' gene targets (from ENCODE/ ChIP-seq data) for the proteins of interest
create_downstream_target_df <- function(protein_ids, genes, gene_targets_df) {
  # Create the raw DF
  downstream_target_df <- data.frame(matrix(ncol = length(protein_ids), 
                                            nrow = nrow(genes)))
  colnames(downstream_target_df) <- protein_ids
  
  # Loop through all transcription factors, extract their potential downstream targets, and
  # add to data frame
  for (i in 1:length(protein_ids)) {
    # Get their ENCODE data, if they have any
    prot_id <- protein_ids[i]
    target_list <- gene_targets_df[,prot_id]    # TODO: will probably have to do some ID conversion here
    # If there was not an ENCODE/ TRANSFAC entry for this protein, use all potential gene targets
    if (length(target_list) == 0) {
      downstream_target_df[,i] <- genes$Gene.Symbol
    }
    # Otherwise, there was an entry! List only these potential targets in the dataframe, removing any NA
    else {
      downstream_target_df[,i] <- target_list[!is.na(target_list)]  # TODO: add NAs to make same length?
    }
  }
  return(downstream_target_df)
}

# Create the targets dataframe
downstream_targ_df <- create_downstream_target_df(binding_prots_w_muts_missense, genes, gene_targets_df)

# Write to a CSV file
write.csv(downstream_targ_df, paste(main_path, "Linear Model/downstream_targets.csv", sep = ""))