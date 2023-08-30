############################################################
### PROCESS PROTEOMIC DATA
### Written By: Sara Geraghty, May 2023
############################################################

library(data.table)
library(purrr)

# This is a file to take all the individual sample proteomic files and merge them into one

# Path to input files
input_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Proteomic_Data/"

# Path to output files
main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


############################################################
### IMPORT PROTEOMIC FILES
############################################################
proteomic_files <- list.files(paste0(input_path, "LGG/"), pattern = "TCGA")   # 435 files
names <- unlist(lapply(proteomic_files, function(f) 
  paste(unlist(strsplit(f, "-", fixed = T))[3:4], collapse = "-")))

new_rows <- lapply(proteomic_files, function(f) 
  fread(paste0(input_path, paste0("LGG/", f)), header = T)[,c('peptide_target', 'protein_expression')])
proteomic_df <- new_rows %>% purrr::reduce(full_join, by = 'peptide_target')
colnames(proteomic_df)[2:ncol(proteomic_df)] <- names

# Write to file
fwrite(proteomic_df, paste0(main_path, "Proteomic/LGG_proteomic.csv"))


############################################################
### VISUALIZE RESULTS USING PROTEOMIC DATA
############################################################
intersecting_patients <- read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids_washu.txt"), header = T)[,1]
proteomic_pats <- unlist(lapply(names, function(n) unlist(strsplit(n, "-", fixed = T))[1]))
proteomic_df_sub <- proteomic_df[,c(1, (which(proteomic_pats %fin% intersecting_patients)+1)), with = F]

mutation_count_df <- fread(paste0(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv"))

# Use get_mutant_patients2 function to get the patients with a mutation in the given driver of interest
idh1_mutant_pats <- get_mutant_patients2(as.data.frame(mutation_count_df), "IDH1")
idh1_mutant_pats_lgg <- intersect(idh1_mutant_pats, colnames(proteomic_df_sub)) # 309 of the 404 patients with an IDH1 mutation are in LGG

# Proteomic target
ensg <- "ENSG00000105281"
symbol <- "SLC1A5"

proteomic_df_sub_targ <- proteomic_df_sub[proteomic_df_sub$peptide_target == "SLC1A5",]