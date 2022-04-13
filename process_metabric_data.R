############################################################
### Process METABRIC Data
### Written By: Sara Geraghty, April 2022
############################################################

# METABRIC Data Retrieved from cBioPortal and downloaded as a .tar file, which was unzipped

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/cBioPortal/brca_metabric/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/cBioPortal/brca_metabric/"

output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/METABRIC/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/METABRIC/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


# The goal of this file is to take in METABRIC data files in order to run our model on another
# dataset aside from the TCGA (as a validation of our findings). This file imports each data type of 
# interest and formats them to match the TCGA files, such that they will be able to flow through the 
# various existing pipelines we have for each data type.


############################################################
### PROCESS MUTATION DATA
############################################################
# Import mutation data file
metabric_mutation_df <- read.table(paste0(main_path, "data_mutations.txt"), header = TRUE, sep = "\t")
metabric_MN_mutation_df <- metabric_mutation_df[(metabric_mutation_df$Variant_Classification == "Missense_Mutation") |
                                                  (metabric_mutation_df$Variant_Classification == "Nonsense_Mutation"),]


############################################################
### PROCESS CNA DATA
############################################################
# Import CNA data file (0 is normal, -1 is deletion of 1 copy, -2 is deletion of 2 copies, 2+ is amplification)
# Will have to use bucketing method for this CNA data
metabric_cna_df <- read.table(paste0(main_path, "data_cna.txt"), header = TRUE, sep = "\t")


############################################################
### PROCESS EXPRESSION DATA
############################################################
metabric_expression_df <- read.table(paste0(main_path, "data_mrna_agilent_microarray.txt"), 
                                     header = TRUE, sep = "\t")


############################################################
### PROCESS METHYLATION DATA
############################################################
metabric_methylation_df <- read.table(paste0(main_path, "data_methylation_promoters_rrbs.txt"), 
                                      header = TRUE, sep = "\t")


############################################################
### PROCESS CLINICAL/ SAMPLE DATA
############################################################
metabric_clinical_patient_df <- read.table(paste0(main_path, "data_clinical_patient.txt"), 
                                           header = TRUE, sep = "\t")
metabric_clinical_sample_df <- read.table(paste0(main_path, "data_clinical_sample.txt"), 
                                          header = TRUE, sep = "\t")