############################################################
### General Important Functions
### Written By: Sara Camilli, November 2020
############################################################

# This file contains functions that are helpful across files
# and will be used the project

############################################################
#' Use a faster %in% function
library(fastmatch)
`%fin%` <- function(x, table) {
  stopifnot(require(fastmatch))
  fmatch(x, table, nomatch = 0L) > 0L
}
############################################################

############################################################
#' Use a faster as.data.frame function
#' @param a named list l with vectors of equal length
quickdf <- function(l) {
  class(l) <- "data.frame"
  attr(l, "row.names") <- .set_row_names(length(l[[1]]))
  return(l)
}
############################################################


############################################################
# Get all the genes in the human genome and write list to a file
# Get all genes in the human genome and related information
#genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg38") 
#genes <- genes[order(genes$start_position),]     # order genes by start position
# Keep only columns of interest
#genes <- genes[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position","end_position")]
#colnames(genes) <- c("Ensembl.ID", "Gene.Symbol","Chr","Start","End")  # rename columns
#write.table(unique(genes$Ensembl.ID), "C:/Users/sarae/Documents/Mona Lab Work/Summer 2020/Main Project Files/Input Data Files/all_gene_ensembl_ids.txt", 
            #row.names = FALSE, col.names = FALSE, quote = FALSE)   # Use the Ensembl IDs for translation using Uniprot's tool

# To read them back...
#all_genes <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_gene_ensembl_ids.txt", header = FALSE)[,1]
############################################################


############################################################
# Use BioMart to produce a unified table of all the genes and their important attributes
#library(biomaRt)
#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#all_genes_id_conv <- getBM(attributes = c('ensembl_gene_id', 'uniprot_gn_id', 'external_gene_name', 'ensembl_transcript_id', 'ensembl_peptide_id',
                                          #'chromosome_name', 'start_position', 'end_position', 'strand', 'transcript_start', 
                                          #'transcript_end', 'transcription_start_site'),
                           #filters = 'ensembl_gene_id',
                           #values = all_genes,
                           #mart = ensembl)
#write.csv(all_genes_id_conv, "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv")

# To read this table back...
#all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Summer 2020/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)
############################################################


############################################################
#' Helper function that uses the ligand groups document to sort 
#' encountered ligands into buckets. Takes a ligand type, a list 
#' of domain types to fill in, and the ligand groups document.
#' Returns the updated domain types list.
#' @param ligand_type a string representing the type of a given ligand
#' that we'd like to categorize
#' @param domain_types_list a list that we're building that contains
#' the counts of each type of ligand we encounter
#' @param ligand_groups a reference table that can be used to group 
#' obscure ligand labels (subcategories) into their broader category
check_lig_type <- function(ligand_type, domain_types_list, ligand_groups) {
  if (ligand_type == "ALL_") { # ALL
    domain_types_list <- lapply(domain_types_list, function(x) x + 1)
  } else if (ligand_type == "NUCACID_") { # NUCACID
    domain_types_list[["NUCACID"]] <- domain_types_list[["NUCACID"]] + 1
  } else if (str_detect(ligand_type, "DNA")) { # DNA
    domain_types_list[["DNA"]] <- domain_types_list[["DNA"]] + 1
  } else if (ligand_type == "METABOLITE_" | ligand_type == "DRUGLIKE_") { # SM
    domain_types_list[["SM"]] <- domain_types_list[["SM"]] + 1
  } else if (str_detect(ligand_type, "RNA")) { # RNA
    domain_types_list[["RNA"]] <- domain_types_list[["RNA"]] + 1
  } else if (ligand_type == "ION_") { # ION
    domain_types_list[["ION"]] <- domain_types_list[["ION"]] + 1
  } else if (ligand_type == "PEPTIDE_") { # PEPTIDE
    domain_types_list[["PEPTIDE"]] <- domain_types_list[["PEPTIDE"]] + 1
  } else if (ligand_type %fin% ligand_groups$original_ligand_mmCIF_identifier) { # Subcategory
    broad_cat <- ligand_groups[ligand_groups$original_ligand_mmCIF_identifier == ligand_type, 'group_name']
    broad_cat <- unique(broad_cat[!is.na(broad_cat)])
    broad_cat <- unlist(lapply(broad_cat, function(x) str_remove(x, "_")))
    for (i in 1:length(broad_cat)) {
      cat <- broad_cat[i]
      if (cat == "METABOLITE" | cat == "DRUGLIKE") {cat <- "SM"}
      domain_types_list[[cat]] <- domain_types_list[[cat]] + 1
    }
  } else {  
    # print("No identifiable ligand match.")
    return(domain_types_list)
  }
  #print(domain_types_list)
  return(domain_types_list)
}

# Import ligand groups file
#ligand_groups <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/ligand_groups.csv", check.names = FALSE)
# Fix colnames, if needed
#colnames(ligand_groups)[1] <- "group_name"
############################################################
