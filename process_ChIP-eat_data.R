############################################################
### Process ChIP-eat Data Files
### Written By: Sara Camilli, March 2021
############################################################

# This file processes ChIP-eat data to create a tissue-type specific list of targets
# for a given transcription factor's files.
# Link to Unibind database with ChIP-eat files: https://unibind.uio.no/

library(dplyr)
library(seqinr)
library(ChIPpeakAnno)
library(EnsDb.Hsapiens.v75)
library(org.Hs.eg.db)

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/ChIP-eat_Data/"
# path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/ChIP-eat_Data/"


############################################################
# IMPORT CHIP-EAT FILES
############################################################
### OPTION 1: BED FILES ### 
bed_files <- list.files(paste(path, "BED Files/TP53/", sep = ""), pattern = ".damo")

# The name of our protein of interest
protein_name <- "TP53"
# protein_name <- "FOXA1"

#' Function takes a vector of bed file names and combines them 
#' into one large table
#' @param path the path to the bed files
#' @param bed_files vector of bed file names
#' @param prot_name the name of the protein we're interested in (bed files
#' will be in a folder under its name)
read_and_combine_bed_files <- function(path, bed_files, prot_name) {
  bed_tables <- lapply(bed_files, function(f) {
    f_name <- paste(path, paste(paste("BED Files/", prot_name, sep = ""), f, sep = "/"), sep = "")
    table <- read.table(f_name, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
    if (ncol(table) == 6) {return(table)}
  })
  comb_bed_table <- distinct(do.call("rbind", bed_tables))
  comb_bed_table <- comb_bed_table[order(comb_bed_table[,1]),]
  return(comb_bed_table)
}

bed_comb_table <- read_and_combine_bed_files(path, bed_files, protein_name)

# Update the column names 
colnames(bed_comb_table) <- c("chrom", "start", "end", "gene_name", "score", "strand")

# Write result to CSV
output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/ChIP-eat/"
# output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/ChIP-eat/"

write.csv(bed_comb_table, paste(output_path, paste(protein_name, "targets_bedfile.csv", sep = "_"), sep = ""))

# Read back
bed_comb_table <- read.csv(paste(output_path, paste(protein_name, "targets_bedfile.csv", sep = "_"), sep = ""), row.names = 1)

#
#
#

### OPTION 2: FASTA FILES ###
fasta_files <- list.files(paste(path, "FASTA Files/", sep = ""), pattern = ".damo")

#' Function takes a vector of FASTA file names and combines
#' all the FASTA files into one 
#' @param path the path to the FASTA files
#' @param fasta_files a vector of the names of the FASTA files
#' @param protein_name the name of the protein we're interested in (FASTA files
#' will be in a folder under its name)
read_and_combine_fasta_files <- function(path, fasta_files, protein_name) {
  fasta_f <- lapply(fasta_files, function(f) {
    f_name <- paste(path, paste(paste("FASTA Files/", protein_name, sep = ""), f, sep = "/"), sep = "")
    fasta <- read.fasta(file = f_name, seqtype = c("AA"), as.string = TRUE)
    return(fasta)
  })
  # Combine into one large FASTA file
  comb_fasta_file <- do.call("c", fasta_f)
  return(comb_fasta_file)
}

comb_fasta_f <- read_and_combine_fasta_files(path, fasta_files, protein_name)

# Write result to FASTA file
write.table(comb_fasta_f, paste(output_path, paste(protein_name, "targets_fasta.fasta")))


############################################################
# ANNOTATE TARGETS USING CHIPSEQANNO PACKAGE
############################################################
#' Function uses the ChIPseqAnno package, with a combined BED file 
#' (all the relevant BED files for a given protein), to create
#' gene target annotations for all peaks. Writes the annotated file and 
#' a list of targets.
#' @param bed_comb_file combined BED file for a given protein
#' @param prot name of protein in question for file, labeling (character)
#' @param output_path path to write outfiles to
get_chipeat_targets <- function(bed_comb_file, prot, output_path) {
  # Get a GRanges object
  peaks <- toGRanges(bed_comb_file, format = "BED", header = TRUE)
  
  # Get annotation data using Homo sapiens reference genome
  annoData <- toGRanges(EnsDb.Hsapiens.v75)
  
  # Annotate the peaks with annotatePeakInBatch
  seqlevelsStyle(peaks) <- seqlevelsStyle(annoData)
  anno <- annotatePeakInBatch(peaks, AnnotationData=annoData)
  
  # Get a pie chart of the overlap features of the peaks
  pie1(table(anno$insideFeature))
  
  # Add gene ID features
  anno <- addGeneIDs(anno, orgAnn="org.Hs.eg.db", 
                     feature_id_type="ensembl_gene_id",
                     IDs2Add=c("symbol"))
  
  # Convert to a DF for easier viewing
  anno_df <- as.data.frame(anno)
  
  # Write this annotation file to a CSV
  write.csv(anno_df, paste(output_path, paste(prot, "_peak_annotation.csv", sep = ""), sep = ""))
  
  # Filter out the rows that have "NA" for the target symbol
  anno_df_filt <- anno_df[!is.na(anno_df$symbol),]
  
  # Get the list of gene targets for this regulatory protein
  targets <- unique(anno_df_filt$symbol)
  length(targets)
  
  # Write these targets to a file
  write(targets, paste(output_path, paste(prot, "_chipeat_targets.txt", sep = ""), sep = ""))
  
  return(targets)
}

targets <- get_chipeat_targets(bed_comb_table, protein_name)

