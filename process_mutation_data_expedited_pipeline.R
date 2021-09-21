############################################################
### Process TCGA Mutation Data File: PIPELINE
### Written By: Sara Camilli, December 2020
############################################################

# This is an expedited pipeline for when we have already processed and saved intermediate files and 
# just want to perform analysis on the results

library(TCGAbiolinks)
library(GenomicRanges)
library(TRONCO)
library(seqinr)
library(stringr)
library(dplyr)
library("RColorBrewer")

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- fread(paste(path, "Input Data Files/all_genes_id_conv.csv", sep = ""), header = TRUE)

############################################################
# IMPORT DATA FROM INTERACDOME, CANBIND, & CONCAVITY
############################################################
# 1. InteracDome

# All Ligands
#interacdome_conf <- fread(paste(path, "Input Data Files/InteracDome/InteracDome_v0.3-confident.csv", sep = ""),
                                  # header = TRUE)
binding_domains_ids_noPF <- fread(paste(path, "Saved Output Data Files/InteracDome/binding_domains_noPF.csv", sep = ""), 
                                  header = TRUE, colClasses = c("x"="character"))[,'x']
binding_domains_ids_noPF <- unlist(binding_domains_ids_noPF)
# Choose which threshold to import
#interacdome_binding_positions_df <- fread(paste(path, "Saved Output Data Files/InteracDome/binding_positions_DF_0.5.csv", sep = ""), header = TRUE)
interacdome_binding_positions_df <- fread(paste(path, "Saved Output Data Files/InteracDome/binding_positions_DF_0_labeled.csv", sep = ""), header = TRUE)
interacdome_domainweights <- fread(paste(path, "Input Data Files/InteracDome/interacdome0.3-pfam31_domainweights-GRCh38.csv", sep = ""))

# Nucleic Acids Only
#interacdome_conf <- fread(paste(path, "Input Data Files/InteracDome/InteracDome_v0.3-confident_nucacids.csv", sep = ""), header = TRUE)
#binding_domains_ids_noPF <- fread(paste(path, "Saved Output Data Files/InteracDome/binding_domains_noPF_nucacids.csv", sep = ""), header = TRUE, colClasses = c("x"="character"))[,'x']
#interacdome_binding_positions_df <- fread(paste(path, "Saved Output Data Files/InteracDome/binding_positions_DF_nucacids_0.csv", sep = ""), header = TRUE)


# 2. CanBind -- choose which threshold to import
# All Ligands
#canbind_df_sub <- fread(paste(path, "Saved Output Data Files/CanBind/canbind_dataframe_labeled_0.5.csv", sep = ""), header = TRUE)
#canbind_swissprot_ids <- unique(canbind_df_sub$Swissprot)
#canbind_binding_site_df <- fread(paste(path, "Saved Output Data Files/CanBind/canbind_binding_ranges_0.5.csv", sep = ""), header = TRUE)

canbind_df_sub <- fread(paste(path, "Saved Output Data Files/CanBind/canbind_dataframe_labeled_0.csv", sep = ""), header = TRUE)
canbind_swissprot_ids <- unique(unlist(lapply(canbind_df_sub$Swissprot, function(x) unlist(strsplit(x, ";")))))
canbind_binding_site_df <- fread(paste(path, "Saved Output Data Files/CanBind/canbind_binding_ranges_0.csv", sep = ""), header = TRUE)

# Nucleic Acids Only
#canbind_df_sub <- fread(paste(path, "Saved Output Data Files/CanBind/canbind_dataframe_labeled_nucacids_0.csv", sep = ""), header = TRUE)
#canbind_swissprot_ids <- unique(canbind_df_sub$Swissprot)
#canbind_binding_site_df <- fread(paste(path, "Saved Output Data Files/CanBind/canbind_binding_ranges_nucacids_0.csv", sep = ""), header = TRUE)

# 3. ConCavity 
concavity_df_sub <- fread(paste(path, "Saved Output Data Files/ConCavity/concavity_dataframe_labeled_0.1.csv", sep = ""), header = TRUE)
concavity_swissprot_ids <- unique(unlist(lapply(concavity_df_sub$Swissprot, function(x) unlist(strsplit(x, ";")))))
concavity_binding_site_df <- fread(paste(path, "Saved Output Data Files/ConCavity/concavity_binding_ranges_0.1.csv", sep = ""), header = TRUE)

# NOTE: ConCavity does not have ligand information; if restriction to nucleic acids only, do not use ConCavity information

############################################################
# IMPORT CLINICAL DATA, MAF FILE, & PROTEOME
############################################################
clin_filename <- paste(path, "Input Data Files/BRCA Data/clinical_data_subset.csv", sep = "")
#clin_filename <- paste(path, "Input Data Files/TCGA Data (ALL)/clinical_data_subset.csv", sep = "")
clinical_df <- fread(clin_filename, header = TRUE)

# Import filtered MAF file
maf_file_df_missense <- fread(paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/maf_file_df_missense.csv", sep = ""), header = TRUE)
#maf_file_df_missense <- fread(paste(path, "Input Data Files/TCGA Data (ALL)/Somatic_Mut_Data/maf_file_df_missense.csv", sep = ""), header = TRUE)

# Imported subsetted proteome file
proteome_subset_missense_filename <- paste(path, "Input Data Files/Proteome/BRCA/proteome_subset_missense.csv", sep = "")
# proteome_subset_missense_filename <- paste(path, "Input Data Files/Proteome/Pan-Cancer/proteome_subset_missense.csv", sep = "")
proteome_subset_missense <- read.fasta(file = proteome_subset_missense_filename, seqtype = c("AA"), as.string = TRUE)

############################################################
# IMPORT DOMAIN FILES
############################################################
# OPT: Import processed domains file from Batch CD-Search if desired
domains_missense <- fread(paste(path, "Input Data Files/BRCA Data/CD-Search/missense_cdsearch_domains_all.csv", sep = ""), header = TRUE,
                             check.names = FALSE)
# domains_missense <- fread(paste(path, "Input Data Files/TCGA Data (ALL)/CD-Search/missense_cdsearch_domains_all.csv", sep = ""), header = TRUE)

# Read in the outputted domains files and protein/Swissprot IDs for each level

# I-Protein
domains_missense_iprotein_sub <- fread(paste(path, "Saved Output Data Files/BRCA/Mutation/iprotein_results_missense.csv", sep = ""), header = TRUE)
# domains_missense_iprotein_sub <- fread(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/iprotein_results_missense.csv", sep = ""), header = TRUE)
protein_ids_missense_iprotein <- fread(paste(path, "Saved Output Data Files/BRCA/Mutation/protein_ids_missense_iprotein.csv", sep = ""))[,2]
swissprot_ids_missense_iprotein <- fread(paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_iprotein.csv", sep = ""))[,2]
#protein_ids_missense_iprotein <- fread(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/protein_ids_missense_iprotein.csv", sep = ""))[,2]
#swissprot_ids_missense_iprotein <- fread(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/swissprot_ids_missense_iprotein.csv", sep = ""))[,2]
maf_subset_df_missense <- fread(paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/maf_file_df_missense_iprotein.csv", sep = ""), header = TRUE)
#maf_subset_df_missense <- fread(paste(path, "Input Data Files/TCGA Data (ALL)/Somatic_Mut_Data/maf_file_df_missense_iprotein.csv", sep = ""), header = TRUE)

# I-Domain
domains_missense_idomain_sub <- fread(paste(path, "Saved Output Data Files/BRCA/Mutation/idomain_results_missense.csv", sep = ""), header = TRUE)
# domains_missense_idomain_sub <- fread(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/idomain_results_missense.csv", sep = ""), header = TRUE)
protein_ids_missense_idomain <- fread(file = paste(path, "Saved Output Data Files/BRCA/Mutation/protein_ids_missense_idomain.csv", sep = ""), header = TRUE)[,2]
swissprot_ids_missense_idomain <- fread(file = paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_idomain.csv", sep = ""), header = TRUE)[,2]
#protein_ids_missense_idomain <- fread(file = paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/protein_ids_missense_idomain.csv", sep = ""), header = TRUE)[,2]
#swissprot_ids_missense_idomain <- fread(file = paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/swissprot_ids_missense_idomain.csv", sep = ""), header = TRUE)[,2]
maf_subset_df_missense_idom <- fread(paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/maf_file_df_missense_idomain.csv", sep = ""), header = TRUE)
#maf_subset_df_missense_idom <- fread(paste(path, "Input Data Files/TCGA Data (ALL)/Somatic_Mut_Data/maf_file_df_missense_idomain.csv", sep = ""), header = TRUE)

# I-Binding Position
domains_missense_ibindingpos_sub <- fread(paste(path, "Saved Output Data Files/BRCA/Mutation/ibindingpos_results_missense.csv", sep = ""), header = TRUE)
#domains_missense_ibindingpos_sub <- fread(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/ibindingpos_results_missense.csv", sep = ""), header = TRUE)
swissprot_ids_missense_ibindingpos <- fread(paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_ibindingpos.csv", sep = ""))[,2]
#swissprot_ids_missense_ibindingpos <- fread(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/swissprot_ids_missense_ibindingpos.csv", sep = ""))
maf_subset_df_missense_ibp <- fread(paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/maf_file_df_missense_ibp.csv", sep = ""))
# maf_subset_df_missense_ibp <- fread(paste(path, "Input Data Files/TCGA Data (ALL)/Somatic_Mut_Data/maf_file_df_missense_ibp.csv", sep = ""))

