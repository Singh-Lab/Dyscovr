############################################################
### ARCHIVE: PIPELINE FOR ALL MUTATIONS & SILENT MUTATIONS
### Written By: Sara Camilli, November 2020
############################################################

# NOTE: If running through this pipeline again, first change the naming on some of the output files 
# to match the up-to-date nomenclature

maf_file_df_silent <- maf_file_df[maf_file_df$Variant_Classification == "Silent",]
# Narrows to 5673 entries, 4229 unique SWISSPROT IDs, 4362 unique Hugo IDs

maf_swissprot_ids <- maf_file_df$SWISSPROT  #maf_file_df[,69] 
maf_swissprot_ids_silent <- maf_file_df_silent$SWISSPROT

maf_symbols <-  maf_file_df$SYMBOL
maf_symbols_silent <-maf_file_df_silent$SYMBOL

proteome_subset_full <- get_proteome_subset(proteome, maf_swissprot_ids)
proteome_subset_silent <- get_proteome_subset(proteome, maf_swissprot_ids_silent)

proteome_subset_full_symb <- get_proteome_subset(proteome, maf_symbols)
proteome_subset_silent_symb <- get_proteome_subset(proteome, maf_symbols_silent)

length(proteome_subset_full)  # check the number of protein sequences remaining
# The length of this subset: 11956
length(proteome_subset_silent)  # check the number of protein sequences remaining
# The length of this subset: 4222

write.fasta(proteome_subset_full, names = names(proteome_subset_full), as.string = TRUE, 
            file.out = paste(path, "Input Data Files/Proteome/proteome_subset.csv", sep = "")) 
write.fasta(proteome_subset_silent, names = names(proteome_subset_silent), as.string = TRUE, 
            file.out = paste(path, "Input Data Files/Proteome/proteome_subset_silent.csv", sep = ""))

proteome_subset_filename <- paste(path, "Input Data Files/Proteome/proteome_subset.fasta", sep = "")
proteome_subset_silent_filename <- paste(path, "Input Data Files/Proteome/proteome_subset_silent.fasta", sep = "")

proteome_subset_full <- read.fasta(file = proteome_subset_filename, seqtype = c("AA"), as.string = TRUE)
proteome_subset_silent <-read.fasta(file = proteome_subset_silent_filename, seqtype = c("AA"), as.string = TRUE)

prep_cdsearch_files(proteome_subset_full, "full")
prep_cdsearch_files(proteome_subset_silent, "silent")

cdsearch_domains <- get_cdsearch_res("full")
cdsearch_domains_silent <-  get_cdsearch_res("silent")

write.csv(cdsearch_domains, paste(path, "Input Data Files/BRCA Data/CD-Search/full_cdsearch_domains.csv", sep = ""))
write.csv(cdsearch_domains_silent, paste(path, "Input Data Files/BRCA Data/CD-Search/silent_cdsearch_domains.csv", sep = ""))

cdsearch_domains <- read.csv(paste(path, "Input Data Files/BRCA Data/CD-Search/full_cdsearch_domains.csv", sep = ""))
cdsearch_domains_silent <- read.csv(paste(path, "Input Data Files/BRCA Data/CD-Search/silent_cdsearch_domains.csv", sep = ""))

accessions <- cdsearch_domains$Accession
accessions_silent <- cdsearch_domains_silent$Accession

accessions <- trim_accessions(accessions)
accessions_silent <- trim_accessions(accessions_silent)

intersecting_domain_acc <- intersect(accessions, binding_domains_ids_noPF)
length(intersecting_domain_acc)
# 145 total intersecting domains
intersecting_domain_acc_silent <- intersect(accessions_silent, binding_domains_ids_noPF)
length(intersecting_domain_acc_silent)
# 88 total intersecting domains

cdsearch_domains_sub <- subset_by_interacdome_and_canbind(cdsearch_domains,
                                                          intersecting_domain_acc,
                                                          canbind_swissprot_ids)   # 867 unique accessions
cdsearch_domains_silent_sub <- subset_by_interacdome_and_canbind(cdsearch_domains_silent, 
                                                                 intersecting_domain_acc_silent, 
                                                                 canbind_swissprot_ids) # 337 unique accessions

cdsearch_domains_sub <- merge_proteome_and_domains(cdsearch_domains_sub, proteome_subset_full)
cdsearch_domains_silent_sub <- merge_proteome_and_domains(cdsearch_domains_silent_sub, proteome_subset_silent)

length(unique(cdsearch_domains_sub$Accession))  # 867
length(unique(cdsearch_domains_silent_sub$Accession))  # 337

write.csv(cdsearch_domains_sub, file = paste(path, "Saved Output Data Files/BRCA/Mutation/binding_proteins_w_muts.csv", sep = ""))
write.csv(cdsearch_domains_silent_sub, file = paste(path, "Saved Output Data Files/BRCA/Mutation/binding_proteins_w_silent_muts.csv", sep = ""))

cdsearch_domains_sub <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/binding_proteins_w_muts.csv", sep = ""), header = TRUE)
cdsearch_domains_silent_sub <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/binding_proteins_w_silent_muts.csv", sep = ""), header = TRUE)


swissprot_ids <- extract_swissprot_ids(cdsearch_domains_sub)
swissprot_ids_silent <- extract_swissprot_ids(cdsearch_domains_silent_sub)

protein_ids <- extract_protein_names(cdsearch_domains_sub) 
protein_ids_silent <- extract_protein_names(cdsearch_domains_silent_sub)   

length(unique(protein_ids)) # 2097
length(unique(protein_ids_silent)) # 781

write.csv(unique(protein_ids), file = paste(path, "Saved Output Data Files/BRCA/Mutation/unique_binding_prots_w_muts.csv", sep = ""))
write.csv(unique(protein_ids_silent), file = paste(path, "Saved Output Data Files/BRCA/Mutation/unique_binding_prots_w_muts_silent.csv", sep = ""))
write.csv(unique(swissprot_ids), file = paste(path, "Saved Output Data Files/BRCA/Mutation/unique_binding_prots_w_muts_swissprot.csv", sep = ""))
write.csv(unique(swissprot_ids_silent), file = paste(path, "Saved Output Data Files/BRCA/Mutation/unique_binding_prots_w_muts_swissprot_silent.csv", sep = ""))

maf_subset_df <- maf_file_df[maf_file_df$SWISSPROT %in% swissprot_ids,]
maf_subset_df_silent <- maf_file_df_silent[maf_file_df_silent$SWISSPROT %in% swissprot_ids_silent,]

prots_w_mut_in_interact_reg <- get_mutations_in_binding_reg(domain_df = cdsearch_domains_sub,
                                                            maf_df = maf_subset_df,
                                                            canbind_df = canbind_df_sub)
prots_w_mut_in_interact_reg_silent <- get_mutations_in_binding_reg(cdsearch_domains_silent_sub,
                                                                   maf_subset_df_silent,
                                                                   canbind_df_sub)

length(unique(prots_w_mut_in_interact_reg$Accession)) # 415
length(unique(prots_w_mut_in_interact_reg_silent$Accession)) # 131

prots_w_muts_in_interact_reg_UNIQUE <- get_unique_prots(prots_w_mut_in_interact_reg, "nam")
prots_w_muts_in_interact_reg_silent_UNIQUE <- get_unique_prots(prots_w_mut_in_interact_reg_silent, "nam")
print(length(prots_w_muts_in_interact_reg_UNIQUE))
print(length(prots_w_muts_in_interact_reg_silent_UNIQUE))

prots_w_muts_in_interact_reg_swissprot_UNIQUE <- get_unique_prots(prots_w_mut_in_interact_reg, "swiss")
prots_w_muts_in_interact_reg_swissprot_silent_UNIQUE <- get_unique_prots(prots_w_mut_in_interact_reg_silent, "swiss")

write.csv(prots_w_mut_in_interact_reg, file = paste(path, "Saved Output Data Files/BRCA/Mutation/prots_w_mut_in_interact_reg_full.csv", sep = ""))
write.csv(prots_w_mut_in_interact_reg_silent, file = paste(path, "Saved Output Data Files/BRCA/Mutation/prots_w_mut_in_interact_reg_silent.csv", sep = ""))

write.csv(prots_w_muts_in_interact_reg_UNIQUE, file = paste(path, "Saved Output Data Files/BRCA/Mutation/unique_binding_prots_w_muts_in_interact_reg.csv", sep = ""))
write.csv(prots_w_muts_in_interact_reg_silent_UNIQUE, file = paste(path, "Saved Output Data Files/BRCA/Mutation/unique_binding_prots_w_muts_in_interact_reg_silent.csv", sep = ""))
write.csv(prots_w_muts_in_interact_reg_swissprot_UNIQUE, file = paste(path, "Saved Output Data Files/BRCA/Mutation/unique_binding_prots_w_muts_in_interact_reg_swissprot.csv", sep = ""))
write.csv(prots_w_muts_in_interact_reg_swissprot_silent_UNIQUE, file = paste(path, "Saved Output Data Files/BRCA/Mutation/unique_binding_prots_w_muts_in_interact_reg_swissprot_silent.csv", sep = ""))

maf_subset_df <- maf_file_df[maf_file_df$SWISSPROT %in% prots_w_muts_in_interact_reg_swissprot_UNIQUE,]
maf_subset_df_silent <- maf_file_df_silent[maf_file_df_silent$SWISSPROT %in% prots_w_muts_in_interact_reg_swissprot_silent_UNIQUE,]

prots_w_mut_in_interact_reg <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/prots_w_mut_in_interact_reg_full.csv", sep = ""), header = TRUE)
prots_w_mut_in_interact_reg_silent <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/prots_w_mut_in_interact_reg_silent.csv", sep = ""), header = TRUE)



