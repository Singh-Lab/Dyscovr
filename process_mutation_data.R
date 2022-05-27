############################################################
### Process TCGA Mutation Data File
### Written By: Sara Camilli, July 2020
############################################################

library(TCGAbiolinks)
library(GenomicRanges)
library(TRONCO)
library(seqinr)
library(stringr)
library(dplyr)
library("RColorBrewer")

# TCGA Workflow: https://www.bioconductor.org/packages/release/data/experiment/vignettes/TCGAWorkflowData/inst/doc/TCGAWorkflowData.html
# TRONCO Package Vignette: https://bioconductor.org/packages/devel/bioc/vignettes/TRONCO/inst/doc/vignette.pdf
  # NOTE: Can also use TRONCO to import GISTEC CNV data

# OVERALL PROCESSING PIPELINE:

# 1. Use InteracDome, CanBind, and ConCavity files to determine binding 1) domains and 2) positions ("confident" regions); 
# 2. Import a MAF file for cancer type of interest; filter out:
    # 1. all patients that do not have all 4 data types of interest
    # 2. patients who are "hypermutators"
    # 3. all genes that do not have a mutation from that cancer type in their bounds
# 3. For all remaining genes, get their conserved domains using either:
    # a) Batch CD-Search (from NCBI)
    # b) Hmmer (2.32 and 3.0)
    # c) the "DOMAINS" annotation in the MAF file
# 4. I-PROTEIN: Filter out genes that do not have a confident binding domain within its bounds
# 5. I-DOMAIN: Use domain boundaries to filter out genes that do not have a mutation
    # within a confident DNA-binding domain or CanBind/ ConCavity region
# 6. I-BINDING REGION: Use specific binding positions from InteracDome and CanBind to 
    # filter out genes that do not have a mutation within a confident DNA binding position


path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste(path, "Input Data Files/all_genes_id_conv.csv", sep = ""), 
                              header = TRUE)

############################################################
# IMPORT DATA FROM INTERACDOME, CANBIND, & CONCAVITY
############################################################
# 1. InteracDome

# All Ligands
interacdome_conf <- read.csv(paste(path, "Input Data Files/InteracDome/InteracDome_v0.3-confident.csv", sep = ""), 
                             header = TRUE, check.names = FALSE)
binding_domains_ids_noPF <- read.csv(paste(path, "Saved Output Data Files/InteracDome/binding_domains_noPF.csv", sep = ""), 
                                     header = TRUE, check.names = FALSE)[,2]
# Choose which threshold to import
#interacdome_binding_positions_df <- read.csv(paste(path, "Saved Output Data Files/InteracDome/binding_positions_DF_0.5.csv", sep = ""), header = TRUE)
interacdome_binding_positions_df <- read.csv(paste(path, "Saved Output Data Files/InteracDome/binding_positions_DF_0.csv", sep = ""), 
                                             header = TRUE, check.names = FALSE, row.names = 1)

# Nucleic Acids Only
#interacdome_conf <- read.csv(paste(path, "Input Data Files/InteracDome/InteracDome_v0.3-confident_nucacids.csv", sep = ""), header = TRUE, check.names = FALSE)
#binding_domains_ids_noPF <- read.csv(paste(path, "Saved Output Data Files/InteracDome/binding_domains_nucacids_noPF.csv", sep = ""), header = TRUE, check.names = FALSE)[,2]
#interacdome_binding_positions_df <- read.csv(paste(path, "Saved Output Data Files/InteracDome/binding_positions_DF_nucacids_0.csv", sep = ""), header = TRUE, check.names = FALSE)


# 2. CanBind -- choose which threshold to import
# All Ligands
#canbind_df_sub <- read.csv(paste(path, "Saved Output Data Files/CanBind/canbind_dataframe_labeled_0.5.csv", sep = ""), header = TRUE)
#canbind_swissprot_ids <- unique(canbind_df_sub$Swissprot)
#canbind_binding_site_df <- read.csv(paste(path, "Saved Output Data Files/CanBind/canbind_binding_ranges_0.5.csv", sep = ""), header = TRUE)

canbind_df_sub <- read.csv(paste(path, "Saved Output Data Files/CanBind/canbind_dataframe_labeled_0.csv", sep = ""), 
                           header = TRUE)
canbind_swissprot_ids <- unique(unlist(lapply(canbind_df_sub$Swissprot, function(x) 
  unlist(strsplit(x, ";")))))
canbind_binding_site_df <- read.csv(paste(path, "Saved Output Data Files/CanBind/canbind_binding_ranges_0.csv", sep = ""), 
                                    header = TRUE)

# Nucleic Acids Only
#canbind_df_sub <- read.csv(paste(path, "Saved Output Data Files/CanBind/canbind_dataframe_labeled_nucacids_0.csv", sep = ""), header = TRUE)
#canbind_swissprot_ids <- unique(canbind_df_sub$Swissprot)
#canbind_binding_site_df <- read.csv(paste(path, "Saved Output Data Files/CanBind/canbind_binding_ranges_nucacids_0.csv", sep = ""), header = TRUE)

# 3. ConCavity 
concavity_df_sub <- read.csv(paste(path, "Saved Output Data Files/ConCavity/concavity_dataframe_labeled_0.1.csv", sep = ""), 
                             header = TRUE)
concavity_swissprot_ids <- unique(unlist(lapply(concavity_df_sub$Swissprot, function(x) 
  unlist(strsplit(x, ";")))))
concavity_binding_site_df <- read.csv(paste(path, "Saved Output Data Files/ConCavity/concavity_binding_ranges_0.1.csv", sep = ""), 
                                      header = TRUE)

# NOTE: ConCavity does not have ligand information; if restriction to nucleic acids only, do not use ConCavity information


############################################################
# GET DNA MUTATION INFORMATION (MAF FILE) 
############################################################
### GET THE MAF FILENAME
# Upload MAF file of interest (do not merge mutation types to "Mutation")
maf_filename <- paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/TCGA.BRCA.muse.b8ca5856-9819-459c-87c5-94e91aca4032.DR-10.0.somatic.maf", sep = "")
# For METABRIC
# maf_filename <- paste(path, "Input Data Files/BRCA Data/cBioPortal/brca_metabric/data_mutations.txt", sep = "")

# maf_filename <- paste(path, "Input Data Files/Pan-Cancer/Somatic_Mut_Data/TCGA.Aggregate.muse.aggregated.somatic.maf", sep = "")
# maf_file <- import.MAF(maf_filename, sep = "\t", is.TCGA = TRUE, merge.mutation.types = FALSE)
  # This will return the MAF file in TRONCO format; set to.TRONCO to FALSE to import as DF
# Alternatively, upload each individual MAF file for each cancer type to separately filter
  # hypermutators by cancer type
maf_filenames <- list.files(paste(path, "Input Data Files/Pan-Cancer/Somatic_Mut_Data/Individual Pan-Cancer MAF Files/", sep = ""), pattern = ".somatic")

### IMPORT THE MAF FILE USING TRONCO
maf_file_df_unfilt <- import.MAF(maf_filename, sep = "\t", is.TCGA = TRUE, 
                                 merge.mutation.types = FALSE, to.TRONCO = FALSE) 
# For METABRIC
#maf_file_df_unfilt <- import.MAF(maf_filename, sep = "\t", is.TCGA = FALSE, 
                                 #merge.mutation.types = FALSE, to.TRONCO = FALSE) 
  # Can use "filter.fun" param to filter MAF lines using a homemade function
  # BRCA: Number of mutations (rows): 90969, Number of unique Swissprot IDs: 16504; Number of unique symbols: 17990
  # TCGA (All): Number of mutations (rows): X, Number of unique Swissprot IDs: X; Number of unique symbols: X
# If we're using individual MAF files....
maf_file_dfs_unfilt <- lapply(maf_filenames, function(x) import.MAF(paste(path, 
                                                                          paste("Input Data Files/Pan-Cancer/Somatic_Mut_Data/Individual Pan-Cancer MAF Files/", x, sep = ""), sep = ""), 
                                                                    sep = "\t", is.TCGA = TRUE, merge.mutation.types = FALSE, to.TRONCO = FALSE) )

############################################################
# FILTER TO INCLUDE ONLY PATIENTS WITH ALL 4 DATA TYPES
############################################################
# NOTE: This assumes that the clinical file has already been processed; make sure that 
# process_clinical_data.R has been run through properly first!
clin_filename <- paste(path, "Input Data Files/BRCA Data/clinical_data_subset.csv", sep = "")
#clin_filename <- paste(path, "Input Data Files/BRCA Data/cBioPortal/brca_metabric/data_clinical_sample_sub.txt", sep = "")
#clin_filename <- paste(path, "Input Data Files/Pan-Cancer/clinical_data_subset.csv", sep = "")

clinical_df <- read.csv(clin_filename, header = TRUE)
#clinical_df <- read.table(clin_filename, header = TRUE, sep = ",")

# Filter by case ID
maf_file_df <- maf_file_df_unfilt[maf_file_df_unfilt$case_id %fin% clinical_df$case_id,]
#maf_file_df <- maf_file_df_unfilt[maf_file_df_unfilt$Tumor_Sample_Barcode %fin% clinical_df$SAMPLE_ID,]
# This limits to only 173 unique proteins with mutations across these 854 patients

# For pan-cancer with all MAF files, filter each by case ID
maf_file_dfs <- lapply(maf_file_dfs_unfilt, function(df) 
  df[df$case_id %fin% clinical_df$case_id,])


############################################################
# FILTER OUT HYPERMUTATORS
############################################################
# Filter out hypermutators (those who are outliers in the number of mutations they have
# in their tumor sample); Use helper functions (maf_helper_functions.R) 
  # We do this BEFORE subsetting to missense mutations, as we want a sense of how many mutations they have as a whole
hypermut_res <- filter_hypermutators(maf_filename, maf_file_df)
maf_file_df <- hypermut_res[[1]] 
  # BRCA: The filtered maf file still has 30305 entries, 12097 unique Swissprot, 12968 unique names
  # TCGA: The filtered maf file still has X entries, X unique Swissprot, X unique names
  # METABRIC: The filtered maf file still has 6058 entries, 173 unique names
total_mut_count_df <- hypermut_res[[2]]  # The filtered mutation count matrix

write.csv(total_mut_count_df, paste(path, "Saved Output Data Files/BRCA/Mutation/Mutation Count Matrices/total_mut_count_matrix_per_patient.csv", sep = ""))
#write.csv(total_mut_count_df, paste(path, "Saved Output Data Files/BRCA/cBioPortal/METABRIC/Mutation/Mutation Count Matrices/total_mut_count_matrix_per_patient.csv", sep = ""))
#write.csv(total_mut_count_df, paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/Mutation Count Matrices/total_mut_count_matrix_per_patient_aggregate.csv", sep = ""))

# If we're using multiple MAF files, filter hypermutators separately on each and then aggregate
# (This is the case for pan-cancer analysis)
total_mut_count_dfs <- list()
for (i in 1:length(maf_filenames)) {
  maf_filename <- maf_filenames[i]
  cancer_type_label <- unlist(strsplit(maf_filename, split = ".", fixed = TRUE))[2]
  maf_df <- maf_file_dfs[[i]]
  hypermut_res <- filter_hypermutators(paste(path, paste("Input Data Files/Pan-Cancer/Somatic_Mut_Data/Individual Pan-Cancer MAF Files/", 
                                                         maf_filename, sep = ""), sep = ""), maf_df)
  maf_file_dfs[[i]] <- hypermut_res[[1]]
  write.csv(hypermut_res[[2]], paste(path, paste("Saved Output Data Files/Pan-Cancer/Mutation/Mutation Count Matrices/total_mut_count_matrix_per_patient_", 
                                                 paste(cancer_type_label, ".csv", sep = ""), sep = ""), sep = ""))
  total_mut_count_dfs[[cancer_type_label]] <- hypermut_res[[2]]
}
maf_file_df <- data.table::rbindlist(maf_file_dfs, use.names = TRUE)
total_mut_count_df <- do.call(rbind, unname(total_mut_count_dfs))

write.csv(maf_file_df, paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/maf_file_hypermut_filt_indiv.csv", sep = ""))
write.csv(total_mut_count_df, paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/total_mut_count_df_hypermut_filt_indiv.csv", sep = ""))


############################################################
# WRITE TOTAL NUMBER OF MUTATIONS TO CLINICAL FILE
############################################################
clinical_df <- update_clin_with_mut_counts(total_mut_count_df, clinical_df, is_tcga = TRUE)
# clinical_df <- update_clin_with_mut_counts(total_mut_count_df, clinical_df, is_tcga = FALSE)

write.csv(clinical_df, clin_filename)  # Update the saved version with total mutation counts

# If already done, read back the file
clinical_df <- read.csv(clin_filename, header = TRUE)


############################################################
# FILTER OUT NON-MISSENSE/NONSENSE MUTATIONS
############################################################
maf_file_df_missense <- maf_file_df[maf_file_df$Variant_Classification == "Missense_Mutation",]
#maf_file_df_misAndNon <- maf_file_df[(maf_file_df$Variant_Classification == "Missense_Mutation") | 
            #(maf_file_df$Variant_Classification == "Nonsense_Mutation"),]
  # BRCA: Narrows to 16689 entries, 8685 unique SWISSPROT IDs, 8954 unique Hugo IDs
  # TCGA (ALL): Narrows to 300796 entries, 17944 unique SWISSPROT IDs, and 18660 unique Hugo IDs
  # METABRIC: Narrows to 3353 entries, 170 unique Hugo IDs
write.csv(maf_file_df_missense, paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/maf_file_df_missense.csv", sep = ""))
#write.csv(maf_file_df_missense, paste(path, "Input Data Files/BRCA Data/cBioPortal/brca_metabric/data_mutations_missense.csv", sep = ""))
#write.csv(maf_file_df_missense, paste(path, "Input Data Files/Pan-Cancer/Somatic_Mut_Data/maf_file_df_missense.csv", sep = ""))

# If already done, read back the file
maf_file_df_missense <- read.csv(paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/maf_file_df_missense.csv", sep = ""), header = TRUE)
#maf_file_df_missense <- read.csv(paste(path, "Input Data Files/Pan-Cancer/Somatic_Mut_Data/maf_file_df_missense.csv", sep = ""), header = TRUE)


############################################################
# GET PROTEOME INFORMATION
############################################################
# Get all proteins in the human proteome (potentially mutated transcription factors, etc. that
# may impact expression); canonical only
  # read.fasta returns a list of character sequences (when as.string = TRUE)
proteome_filename <- paste(path, "Input Data Files/Proteome/uniprot-reviewed_yes+AND+proteome_up000005640.fasta", sep = "")
proteome <- read.fasta(file = proteome_filename, seqtype = c("AA"), as.string = TRUE)
  # Length of proteome: 20353


############################################################
# FILTER PROTEINS BY PRESENCE OF MUTATIONS
############################################################
# Determine shared ID between the MAF and FASTA files; keep only proteome proteins at the
# intersection of these two files
  # SWISSPROT ID (column 69 in MAF)
maf_swissprot_ids_missense <- maf_file_df_missense$SWISSPROT
#maf_swissprot_ids_missense <- unlist(lapply(maf_file_df_missense$Hugo_Symbol, function(x) 
        #paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'uniprot_gn_id'])), collapse = ";")))

maf_symbols_missense <- maf_file_df_missense$SYMBOL
#maf_symbols_missense <- maf_file_df_missense$Hugo_Symbol


#' Takes a proteome and set of protein Swissprot or Hugo IDs from a given MAF file; 
#' returns a subsetted proteome DF with only overlapping proteins
#' @param proteome a proteome from Uniprot or elsewhere, in FASTA format
#' @param maf_file_ids vector of SWISSPROT or HUGO IDs from a given MAF file
#' @param id_type a string denoting the ID type (either 'swissprot' or 'hugo')
get_proteome_subset <- function(proteome, maf_file_ids, id_type) {
  # Split up the maf file ids, if needed
  maf_file_ids <- unique(unlist(lapply(maf_file_ids, function(x) unlist(strsplit(x, ";", fixed = TRUE)))))

  proteome_subset <- c()
  for (i in 1:length(proteome)) {
    entry <- proteome[i]
    entry_id <- ""
    if (id_type == 'swissprot') {
      entry_id <- unlist(strsplit(attributes(proteome[i])[[1]], "|", fixed = TRUE))[2]
    } else if (id_type == 'hugo') {
      entry_id <- unlist(strsplit(unlist(strsplit(attributes(proteome[i])[[1]], "|", fixed = TRUE))[3], "_", fixed =  TRUE))[1]
    } else {
      print(paste("Unrecognized ID type:", id_type))
    }
    if (entry_id %in% maf_file_ids) {
      proteome_subset <- c(proteome_subset, entry)
    }
  }
  return(proteome_subset)
}

# Swissprot ID subset
proteome_subset_missense <- get_proteome_subset(proteome, maf_swissprot_ids_missense, 'swissprot')
# Symbol subset
proteome_subset_missense_symb <- get_proteome_subset(proteome, maf_symbols_missense, 'hugo')

length(unique(proteome_subset_missense))  # check the number of protein sequences remaining
  # The length of this subset: BRCA - 8675 (11068 for missense + nonsense), TCGA (ALL) - 17899 (18062 for missense + nonsense)
  # 4242 for BRCA for silent, 16725 for Pan-Cancer for silent
  # 160 for METABRIC BRCA missense, 162 for METABRIC BRCA missense + nonsense

# Write proteome subset to a FASTA
write.fasta(proteome_subset_missense, names = names(proteome_subset_missense), as.string = TRUE, 
            file.out = paste(path, "Input Data Files/Proteome/BRCA/proteome_subset_missense.csv", sep = ""))
#write.fasta(proteome_subset_missense, names = names(proteome_subset_missense), as.string = TRUE, 
            #file.out = paste(path, "Input Data Files/Proteome/Pan-Cancer/proteome_subset_missense.csv", sep = ""))

# Alternatively, if we've already done this, just read back proteome subset from FASTA
proteome_subset_missense_filename <- paste(path, "Input Data Files/Proteome/BRCA/proteome_subset_missense.csv", sep = "")
# proteome_subset_missense_filename <- paste(path, "Input Data Files/Proteome/Pan-Cancer/proteome_subset_missense.fasta", sep = "")
proteome_subset_missense <- read.fasta(file = proteome_subset_missense_filename, seqtype = c("AA"), as.string = TRUE)


############################################################
# GET DOMAIN INFORMATION FOR ALL PROTEINS IN THE HUMAN PROTEOME
############################################################
# Keep only proteins that have at least one confident binding domain

# Export the list of remaining proteins back to FASTA format so that it can be uploaded to 
  # Batch CD-Search, a tool that accepts protein sequences in FASTA format (max 4000 entries/ batch) or to HMMER 
  # Batch CD-Search link: https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
  # NOTE: Uniprot's InterPro tool is another option for domain annotation

prep_cdsearch_files(proteome_subset_missense, "missense", "tcga")   # see batch_cdsearch_helper_functions.R
#prep_cdsearch_files(proteome_subset_missense, "missense", "metabric")   # see batch_cdsearch_helper_functions.R


############################################################
# IMPORT RESULTS FROM BATCH CD-SEARCH 
############################################################
domains_missense <- get_cdsearch_res("missense", "tcga")      # see batch_cdsearch_helper_functions.R
#domains_missense <- get_cdsearch_res("missense", "metabric")      # see batch_cdsearch_helper_functions.R

# Write to a file 
write.csv(domains_missense, paste(path, "Input Data Files/BRCA Data/CD-Search/missense_cdsearch_domains_all.csv", sep = ""))
# write.csv(domains_missense, paste(path, "Input Data Files/BRCA Data/CD-Search/missense_cdsearch_domains_all_metabric.csv", sep = ""))
# write.csv(domains_missense, paste(path, "Input Data Files/Pan-Cancer/CD-Search/missense_cdsearch_domains_all.csv", sep = ""))

# If we have already done the above once, simply import the saved domains files
domains_missense <- read.csv(paste(path, "Input Data Files/BRCA Data/CD-Search/missense_cdsearch_domains_all.csv", sep = ""), header = TRUE,
                             check.names = FALSE)
# domains_missense <- read.csv(paste(path, "Input Data Files/Pan-Cancer/CD-Search/missense_cdsearch_domains_all.csv", sep = ""), header = TRUE)


############################################################
# EXTRACT PFAM IDs FROM BATCH CD-SEARCH RESULTS
############################################################
# NOTE: Query #s in the form Q#N - XXXXXXXX, where XXXXXXXX is the sequence ID or 
# first 15 chars of FASTA line
# Get Pfam domains from Batch CD-Search Results

#' Subset domains data frame to only include pfam domains
#' @param domain_df a domain data frame to be subsetted
subset_domains_to_pfam <- function(domain_df) {
  domain_df <- domain_df[grep("pfam", domain_df$Accession, ignore.case = TRUE),]
}
domains_missense <- subset_domains_to_pfam(domains_missense)
accessions_missense <- domains_missense$Accession
  
# Get only the numeric portions of the accessions, and only those from the Pfam database
  # (InteracDome only uses Pfam domains, so we will not get matches for SMART or COG)
accessions_numeric <- regmatches(accessions_missense, regexpr("[[:digit:]]+", accessions_missense)) 


############################################################
# GET DOMAIN OVERLAP BETWEEN BATCH CD-SEARCH AND INTERACDOME
############################################################
# Which of these domains are confidently binding?
# Look for overlap between the domain hits from CD-Search and 
  # the DNA-binding domains from InteracDome
intersecting_domain_acc_missense <- intersect(accessions_numeric, binding_domains_ids_noPF)
length(intersecting_domain_acc_missense)

  # BRCA is 1162 intersecting domains (1250 for missense + nonsense, 895 for silent); TCGA (ALL) is 1453 domains (1437 for missense + nonsense, 1414 for silent)
  # NUC ACIDS ONLY: BRCA is 218 intersecting domains, TCGA (ALL) is 279

  # METABRIC has 131 intersecting domains for missense, 131 for missense + nonsense

############################################################
# I-PROTEIN: SUBSET PROTEOME TO KEEP ONLY PROTEINS WITH EITHER:
  # 1. CONFIDENT BINDING DOMAIN FROM INTERACDOME
  # 2. CONFIDENT BINDING REGION FROM CANBIND
  # 3. CONFIDENT BINDING REGION FROM CONCAVITY
############################################################
#' Subsets the DF from Batch CD-Search to include only proteins that have either 
#' a) a binding domain from InteracDome or b) a binding region from CanBind or 
#' c) a binding region from ConCavity or d) multiple. Returns subset DF.
#' @param domains a data frame of domain information from Batch-CD search
#' @param intersecting_domain_acc a vector of accession numbers for domains that 
#' overlap between Batch CD-Search and InteracDome 
#' @param canbind_swissprot_ids a vector of swissprot IDs with confident interacting 
#' regions from CanBind
#' @param concavity_swissprot_ids OPT: a vector of swissprot IDs with confident 
#' interacting regions from ConCavity
get_iprotein_subset <- function(domains, intersecting_domain_acc, canbind_swissprot_ids, 
                                concavity_swissprot_ids) {
  # Fix the "Query" column name (was i..Query)
  if("ï»¿Query" %in% colnames(domains)) {
    colnames(domains)[which(colnames(domains) == "ï»¿Query")] <- "Query"
  }
  
  # Fix the query column to get rid of query number (we don't really need this)
  domains$Query <- as.character(unlist(lapply(domains$Query, function(x) 
    unlist(strsplit(x, " ", fixed = TRUE))[3])))
  
  # Subset by InteracDome domains
  domains_interacdome_sub <- domains[as.character(regmatches(domains$Accession, 
                                                             gregexpr("[[:digit:]]+",
                                                             domains$Accession)))
                                       %fin% intersecting_domain_acc,]

  # Subset by CanBind proteins with interaction regions
  queries <- domains$Query
  swissprot_ids <- unlist(sapply(queries, FUN = function(x) {unlist(strsplit(x, "|", fixed = TRUE))[2]}))
  matching_indices_canbind <- which(canbind_swissprot_ids %fin% swissprot_ids)
  domains_canbind_sub <- domains[matching_indices_canbind,]
  
  # Subset by Concavity proteins with interaction regions, if provided
  if (!missing(concavity_swissprot_ids)) {
    matching_indices_concavity <- which(concavity_swissprot_ids %fin% swissprot_ids)
    domains_concavity_sub <- domains[matching_indices_concavity,]
    
    # Take the union of these data frames to include all those that remain using either method
    dupl_rows_pt1 <- which(!is.na(match(domains_interacdome_sub$Accession, domains_canbind_sub$Accession)))
    domains_sub_pt1 <- rbind(domains_interacdome_sub, domains_canbind_sub[-dupl_rows_pt1,])
    dupl_rows_pt2 <- !is.na(match(domains_sub_pt1$Accession, domains_concavity_sub$Accession))
    domains_sub_pt2 <- rbind(domains_sub_pt1, domains_concavity_sub[-dupl_rows_pt2,])
    domains_sub <- domains_sub_pt2[!is.na(domains_sub_pt2$Query),]
    #print(nrow(domains_sub)) 
  }
  else {
    dupl_rows <- which(!is.na(match(domains_interacdome_sub$Accession, domains_canbind_sub$Accession)))
    domains_sub <- rbind(domains_interacdome_sub, domains_canbind_sub[-dupl_rows,])
    domains_sub <- domains_sub[!is.na(domains_sub$Query),]
  }

  return(domains_sub)
}


domains_missense_iprotein_sub <- get_iprotein_subset(domains_missense, intersecting_domain_acc_missense, 
                                         canbind_swissprot_ids, concavity_swissprot_ids) 
# No ConCavity version for nucleic acids only
# domains_missense_iprotein_sub <- get_iprotein_subset(domains_missense, intersecting_domain_acc_missense, 
  #canbind_swissprot_ids) 


############################################################
# ADD SEQUENCE INFORMATION TO SUBSETTED PROTEIN DF
############################################################
#' Takes data frame of domain information from Batch-CD search, along with
#' a proteome, and extracts the AA acid sequence for given query, which it adds as
#' an additional column to the data frame. Returns this data frame.
#' @param domains the data frame of domain information from Batch-CD search
#' @param proteome a FASTA-style proteome
merge_proteome_and_domains <- function(domains, proteome) {
  sequence <- c()   # list of sequence strings
  for (i in 1:nrow(domains)) {
    query <- str_remove(domains$Query[i], ">")
    proteome_entry <- proteome[which(attributes(proteome)$names == query)]

    # Get elements of interest from entry
    sequence <- c(sequence, as.character(proteome_entry[1]))
  }
  domains <- cbind(domains, sequence)
  return(domains)
}

domains_missense_iprotein_sub <- merge_proteome_and_domains(domains_missense_iprotein_sub, 
                                                            proteome_subset_missense)


############################################################
# GET UNIQUE PROTEIN IDs AND WRITE TO FILES
############################################################
# Examine how many actual unique accessions we have
length(unique(domains_missense_iprotein_sub$Accession))  
  # All Ligands: BRCA: 1910 (2092 missense + nonsense, 1237 silent); TCGA (ALL): 2767 (2749 missense + nonsense, 2644 silent)
  # Nucleic Acids Only: BRCA: 257; TCGA (ALL): 324

  # METABRIC BRCA: 133 missense (133 missense + nonsense)

# Write this table to a CSV
write.csv(domains_missense_iprotein_sub, file = paste(path, "Saved Output Data Files/BRCA/Mutation/iprotein_results_missense.csv", sep = ""))
# write.csv(domains_missense_iprotein_sub, file = paste(path, "Saved Output Data Files/BRCA/cBioPortal/METABRIC/Mutation/iprotein_results_missense.csv", sep = ""))
# write.csv(domains_missense_iprotein_sub, file = paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/iprotein_results_missense.csv", sep = ""))

# If already done this, read them back from CSV
domains_missense_iprotein_sub <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/iprotein_results_missense.csv", sep = ""), header = TRUE, check.names = FALSE)
# domains_missense_iprotein_sub <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/iprotein_results_missense.csv", sep = ""), header = TRUE, check.names = FALSE)
# domains_missense_iprotein_sub <- read.csv(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/iprotein_results_missense.csv", sep = ""), header = TRUE, check.names = FALSE)

# Extract the Swissprot IDs and general names for all the remaining queries

#' Return all the SWISSPROT IDs from a domain data frame
#' @param domain_df the domain data frame to find IDs in 
extract_swissprot_ids <- function(domain_df) {
  swissprot_ids <- unique(unlist(lapply(domain_df$Query, function(x) 
    unlist(strsplit(x,"|", fixed = TRUE))[2])))
  return(swissprot_ids)
}

#' Return all the protein names from a domain data frame
#' @param domain_df the domain data frame to list the gene names from
extract_protein_names <- function(domain_df) {
  protein_ids <- c()
  for (i in 1:nrow(domain_df)) {
    split_id <- unlist(strsplit(domain_df$Query[i],"|", fixed = TRUE))
    protein_ids <- c(protein_ids, unlist(strsplit(split_id[3], "_", fixed = TRUE))[1])
  }
  return(unique(protein_ids))
}

swissprot_ids_missense_iprotein <- extract_swissprot_ids(domains_missense_iprotein_sub)
protein_ids_missense_iprotein <- extract_protein_names(domains_missense_iprotein_sub)  
length(protein_ids_missense_iprotein)

# All Ligands: BRCA: 6214 (7808 with missense + nonsense, 3040 silent); TCGA (ALL): 11,939 (11,989 with missense + nonsense, 11,348 silent)
# Nuc. Acids Only: BRCA: 1914; TCGA (ALL): 3,644

# METABRIC BRCA: 121 missense (122 for missense + nonsense)

# Write unique protein IDs to a file
write.csv(protein_ids_missense_iprotein, file = paste(path, "Saved Output Data Files/BRCA/Mutation/protein_ids_missense_iprotein.csv", sep = ""))
write.csv(swissprot_ids_missense_iprotein, file = paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_iprotein.csv", sep = ""))
#write.csv(protein_ids_missense_iprotein, file = paste(path, "Saved Output Data Files/BRCA/cBioPortal/METABRIC/Mutation/protein_ids_missense_iprotein.csv", sep = ""))
#write.csv(swissprot_ids_missense_iprotein, file = paste(path, "Saved Output Data Files/BRCA/cBioPortal/METABRIC/Mutation/swissprot_ids_missense_iprotein.csv", sep = ""))
#write.csv(protein_ids_missense_iprotein, file = paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/protein_ids_missense_iprotein.csv", sep = ""))
#write.csv(swissprot_ids_missense_iprotein, file = paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/swissprot_ids_missense_iprotein.csv", sep = ""))

# Read back if necessary
protein_ids_missense_iprotein <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/protein_ids_missense_iprotein.csv", sep = ""))[,2]
swissprot_ids_missense_iprotein <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_iprotein.csv", sep = ""))[,2]
#protein_ids_missense_iprotein <- read.csv(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/protein_ids_missense_iprotein.csv", sep = ""))[,2]
#swissprot_ids_missense_iprotein <- read.csv(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/swissprot_ids_missense_iprotein.csv", sep = ""))[,2]


# Use these IDs to subset the MAF file
maf_subset_df_missense <- maf_file_df_missense[maf_file_df_missense$SWISSPROT %in% swissprot_ids_missense_iprotein,]

# For METABRIC, have to add the SWISSPROT ID and then do the subsetting (does not work properly with Hugo symbol)
maf_file_df_missense$SWISSPROT <- unlist(lapply(maf_file_df_missense$Hugo_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'uniprot_gn_id'])), collapse = ";")))
rows_to_keep <- unlist(lapply(1:nrow(maf_file_df_missense), function(i) {
  prot_ids <- unlist(strsplit(maf_file_df_missense$SWISSPROT[i], ";", fixed = TRUE))
  truth_vect <- unlist(lapply(prot_ids, function(id) ifelse(id %in% swissprot_ids_missense_iprotein, TRUE, FALSE)))
  if(TRUE %in% truth_vect) {return(i)}
  else {return(NA)}
}))
rows_to_keep <- rows_to_keep[!is.na(rows_to_keep)]
maf_subset_df_missense <- maf_file_df_missense[rows_to_keep,]

write.csv(maf_subset_df_missense, paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/maf_file_df_missense_iprotein.csv", sep = ""))
#write.csv(maf_subset_df_missense, paste(path, "Input Data Files/BRCA Data/cBioPortal/brca_metabric/data_mutations_missense_iprotein.csv", sep = ""))
#write.csv(maf_subset_df_missense, paste(path, "Input Data Files/Pan-Cancer/Somatic_Mut_Data/maf_file_df_missense_iprotein.csv", sep = ""))

# If already done, read back the file
maf_subset_df_missense <- read.csv(paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/maf_file_df_missense_iprotein.csv", sep = ""), header = TRUE)
#maf_subset_df_missense <- read.csv(paste(path, "Input Data Files/Pan-Cancer/Somatic_Mut_Data/maf_file_df_missense_iprotein.csv", sep = ""), header = TRUE)


############################################################
# OPTIONAL: ADD PATIENTS AND MUTATION POSITIONS TO DATA FRAME
############################################################
#' Add corresponding patient IDs that have mutations in each protein to the 
#' given I-Protein-level data frame using the information from the MAF file
#' @param iprotein_sub an I-Protein level data frame produced from the above
#' @param maf a MAF file with patient mutational information
#' @param dataset either 'tcga', 'metabric', 'icgc' or 'cptac3'
add_patient_ids <- function(iprotein_sub, maf, dataset) {
  # For each unique protein, get the patients that are mutated in it
  full_prot_set <- unlist(lapply(iprotein_sub$Query, function(x) 
    unlist(strsplit(x, "|", fixed = TRUE))[2]))
  unique_proteins <- unique(full_prot_set)
  
  patients_mut_per_prot <- lapply(unique_proteins, function(p) {
    # Subset the MAF to only the given protein
    maf_sub <- maf[grepl(p, maf$SWISSPROT),]
    
    # Get the patient and sample that have this protein mutated
    unique_patients <- unique(maf_sub$Tumor_Sample_Barcode)
    
    # If TCGA, split this ID apart to get just patient and sample XXXX-XX
    if(dataset == "tcga") {
      unique_patients <- unlist(lapply(unique_patients, function(x) 
        paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))
    }
    # Collapse these patients into a semicolon-separated string
    patient_str <- paste(unique_patients, collapse = ";")
    print(head(patient_str))
    return(patient_str)
  })
  names(patients_mut_per_prot) <- unique_proteins
  
  new_col <- lapply(full_prot_set, function(p) patients_mut_per_prot[[p]])
  print(head(new_col))
  iprotein_sub$Patient <- unlist(new_col)

  return(iprotein_sub)
}

domains_missense_iprotein_sub <- add_patient_ids(domains_missense_iprotein_sub, maf_subset_df_missense, "tcga")
#domains_missense_iprotein_sub <- add_patient_ids(domains_missense_iprotein_sub, maf_subset_df_missense, "metabric")


write.csv(domains_missense_iprotein_sub, file = paste(path, "Saved Output Data Files/BRCA/Mutation/iprotein_results_missense.csv", sep = ""))
#write.csv(domains_missense_iprotein_sub, file = paste(path, "Saved Output Data Files/BRCA/cBioPortal/METABRIC/Mutation/iprotein_results_missense.csv", sep = ""))
#write.csv(domains_missense_iprotein_sub, file = "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Mutation/iprotein_results_missense.csv")


############################################################
# I-DOMAIN: KEEP ONLY PROTEINS WITH A MUTATION ~IN~ A DNA-
# BINDING DOMAIN REGION FROM INTERACDOME OR INTERACTION
# REGION FROM CANBIND/ CONCAVITY
############################################################
# We know the region of the domain in the "To" and "From" columns of the data frame
# Now, we need to check whether, for this particular protein, there are mutations in that
# region in BRCA cancer patients

#' For every protein, see if the "Start" and "End" of any of the mutations falls 
#' within the domain region. Compile all those that do into a new, final DF and return
#' @param domain_df domain data frame from the CD-Batch Search Output, subsetted to 
#' include only InteracDome/ CanBind/ ConCavity domains
#' @param maf_df MAF-style data frame that will have the "Start" and "End" bases 
#' of the mutations of interest
#' @param canbind_binding_site_df a CanBind data frame that has the start and end 
#' positions of each binding site for each protein
#' @param concavity_binding_site_df OPT: a ConCavity data frame that has start 
#' and end positions of each binding site for each protein
get_idomain_subset <- function(domain_df, maf_df, canbind_binding_site_df, concavity_binding_site_df) {
  # Create a DF to hold the final results
  idomain_prot_df <- data.frame(matrix(nrow = 0, ncol = 9))
  
  # Get all the unique proteins in the domain DF
  unique_protein_ids <- unique(unlist(lapply(domain_df$Query, function(x) 
    unlist(strsplit(x, "|", fixed = TRUE))[2])))
  
  for (id in unique_protein_ids) {
    
    print(paste(match(id, unique_protein_ids), paste("/", length(unique_protein_ids))))
    
    # Subset the domain DF, MAF file, CanBind DF, and ConCavity DF (if provided) 
    # to look only at this protein. For each mutation in this protein, check if 
    # it is in a domain or binding region for CanBind and/or ConCavity
    domain_df_sub <- domain_df[grep(id, domain_df$Query),]
    maf_sub <- maf_df[maf_df$SWISSPROT == id,]
    canbind_ranges <- unique(canbind_binding_site_df[grepl(id, canbind_binding_site_df$Swissprot), 
                                                     'Range.of.Binding.Site'])
    if (!missing(concavity_binding_site_df)) {
      concavity_ranges <- unique(concavity_binding_site_df[grepl(id, concavity_binding_site_df$Swissprot), 
                                                           'Range.of.Binding.Site'])
      prot_df <- get_overlap_with_binding_domains(id, domain_df_sub, maf_sub, 
                                                  canbind_ranges, concavity_ranges)
    } else {
      prot_df <- get_overlap_with_binding_domains(id, domain_df_sub, maf_sub, canbind_ranges)
    }
    
    # Bind this protein's data frame to the main DF
    idomain_prot_df <- rbind(idomain_prot_df, prot_df)
  }
  colnames(idomain_prot_df) <- c("Swissprot", "From", "To", "Source", "Accession", 
                                 "Short Name", "Superfamily", "Mut.Pos", "Patient")
  print(idomain_prot_df)
  return(idomain_prot_df)
}

#' Helper function to get the overlap of mutations with 1. conversed binding domains 
#' or 2. CanBind or (optionally) ConCavity binding ranges. Returns a vector of rows 
#' to be added to the new data frame.
#' @param id swissprot ID of current protein
#' @param domain_df_sub the domain DF, subsetted to only this protein
#' @param maf_sub the MAF file, subsetted to only this protein
#' @param canbind_ranges the CanBind binding ranges data frame, subsetted to only this protein
#' @param concavity_ranges the ConCavity binding ranges data frame, subsetted to only this protein
get_overlap_with_binding_domains <- function(id, domain_df_sub, maf_sub, 
                                             canbind_ranges, concavity_ranges) {
  # Iterate over all the mutations in this given protein
  df_full <- lapply(maf_sub$Protein_position, function(x) {
    
    # For each mutation in this protein, get its AA position ("Protein_position", column 55) 
    mut_pos <- unlist(strsplit(x, "/", fixed = TRUE))[1]
    print(mut_pos)
    
    # Make sure this mutation is in the coding region
    if (!mut_pos == "-") {
      # Also get the patient under consideration
      patient_id <- unlist(strsplit(maf_sub$Tumor_Sample_Barcode[which(maf_sub$Protein_position == x)], 
                                    "-", fixed = TRUE))[3]
      
      # First, check if the mutation is in a conserved binding domain 
      domain_rows <- lapply(1:nrow(domain_df_sub), function(i) {
        start_dom <- as.numeric(domain_df_sub[i, 'From'])
        end_dom <- as.numeric(domain_df_sub[i, 'To'])

        if ((as.numeric(mut_pos) >= start_dom) & (as.numeric(mut_pos) <= end_dom)) {
          # Return the filled row
          # Columns to fill: "Swissprot", "From", "To", "Source", "Accession", 
            # "Short Name", "Superfamily", "Mut.Pos", "Patient"
          return(c(id, start_dom, end_dom, "Conserved Domain", domain_df_sub[i, 'Accession'], 
                   domain_df_sub[i, 'Shortname'], domain_df_sub[i, 'Superfamily'], mut_pos, patient_id))
          # Or return(c(id, start_dom, end_dom, "Conserved Domain", domain_df_sub[i, 'Accession'], 
            # domain_df_sub[i, 'Short.name'], domain_df_sub[i, 'Superfamily'], mut_pos, patient_id)) 
            # depending on headers
        }
      })
      #print(domain_rows)
      domain_df <- do.call(rbind, domain_rows)
      
      # Then, check if the mutation is in a binding region for CanBind or ConCavity (if applicable)
      canbind_df <- get_canbind_and_concavity_rows(canbind_ranges, "CanBind", 
                                                   domain_df_sub, mut_pos, patient_id, id)
      #print(paste("CanBind df ", canbind_df))
      
      try({
        concavity_df <- get_canbind_and_concavity_rows(concavity_ranges, "ConCavity", 
                                                       domain_df_sub, mut_pos, patient_id, id)
      }, silent = TRUE)
      #try({print(paste("ConCavity df ", concavity_df))})
      
      # Merge all the rows from these three sources and sort by protein ID
      df <- rbind(domain_df, canbind_df)
      try({df <- rbind(df, concavity_df)}, silent = TRUE)
      
      print(df)
      return(df)
    }
  })
  df_full <- do.call(rbind, df_full)
  return(df_full)
}

#' Helper function that takes the binding ranges for either CanBind or ConCavity, 
#' as well as the information needed to fill in the row, and returns the rows 
#' for mutations that fall in a given binding range
#' @param binding_ranges a set of binding ranges from either CanBind or ConCavity (start:end)
#' @param label a 'source' label (either 'CanBind' or 'ConCavity')
#' @param domain_df_sub a domain DF with information about mutations and domains
#' @param mut_pos the position of the mutation in question
#' @param patient_id the 4-digit patient TCGA ID of the patient in question
#' @param id the swissprot ID of the current protein
get_canbind_and_concavity_rows <- function(binding_ranges, label, domain_df_sub, 
                                           mut_pos, patient_id, id) {
  rows <- lapply(binding_ranges, function(y) {
    start_pos <- as.numeric(unlist(strsplit(y, ":", fixed = TRUE))[1])
    end_pos <- as.numeric(unlist(strsplit(y, ":", fixed = TRUE))[2])
    
    if ((as.numeric(mut_pos) >= start_pos) & (as.numeric(mut_pos) <= end_pos)) {
      # Return the filled row
      # Columns to fill: "Swissprot", "From", "To", "Source", "Accession", "Short Name", "Superfamily", "Mut.Pos", "Patient"
      return(c(id, start_pos, end_pos, label, NA, NA, NA, mut_pos, patient_id))
    }
  })
  df <- do.call(rbind, rows)
  return(df)
}
    

domains_missense_idomain_sub <- get_idomain_subset(domains_missense_iprotein_sub,
                                                   maf_subset_df_missense,
                                                   canbind_binding_site_df, concavity_binding_site_df)
# For nucleic acids, not using ConCavity
#domains_missense_idomain_sub <- get_idomain_subset(domains_missense_iprotein_sub,
#                                                   maf_subset_df_missense,
#                                                   canbind_binding_site_df)

# OPT: Delete any duplicate rows
domains_missense_idomain_sub <- distinct(domains_missense_idomain_sub)


############################################################
# GET UNIQUE PROTEIN IDs AND WRITE TO FILES
############################################################
# Get the number of unique domains in this DF:
length(unique(domains_missense_idomain_sub$Accession)) 

# All Ligands: BRCA: 1130; Pan-Cancer: 2605
# Nuc. Acids Only: BRCA: 151; Pan-Cancer: 310

# Get the Swissprot IDs
swissprot_ids_missense_idomain <- unique(domains_missense_idomain_sub$Swissprot)
print(length(swissprot_ids_missense_idomain)) 

# All Ligands: BRCA: 3997; Pan-Cancer: 11305
# Nuc. Acids Only: BRCA: 1190; Pan-Cancer: 3274

# Write the full table to a CSV
write.csv(domains_missense_idomain_sub, file = paste(path, "Saved Output Data Files/BRCA/Mutation/idomain_results_missense.csv", sep = ""))
#write.csv(domains_missense_idomain_sub, file = paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/idomain_results_missense.csv", sep = ""))

# Write the unique protein list to a CSV
#write.csv(protein_ids_missense_idomain, file = paste(path, "Saved Output Data Files/BRCA/Mutation/protein_ids_missense_idomain.csv", sep = ""))
write.csv(swissprot_ids_missense_idomain, file = paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_idomain.csv", sep = ""))
#write.csv(protein_ids_missense_idomain, file = paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/protein_ids_missense_idomain.csv", sep = ""))
#write.csv(swissprot_ids_missense_idomain, file = paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/swissprot_ids_missense_idomain.csv", sep = ""))

# Use these IDs to further subset the MAF file
maf_subset_df_missense_idom <- maf_subset_df_missense[maf_subset_df_missense$SWISSPROT %in% swissprot_ids_missense_idomain,]
write.csv(maf_subset_df_missense_idom, paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/maf_file_df_missense_idomain.csv", sep = ""))
#write.csv(maf_subset_df_missense_idom, paste(path, "Input Data Files/Pan-Cancer/Somatic_Mut_Data/maf_file_df_missense_idomain.csv", sep = ""))

# Read these files back if needed
# 
# 
# 
maf_subset_df_missense <- read.csv(paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/maf_file_df_missense_idomain.csv", sep = ""), header = TRUE)
#maf_subset_df_missense <- read.csv(paste(path, "Input Data Files/Pan-Cancer/Somatic_Mut_Data/maf_file_df_missense_idomain.csv", sep = ""), header = TRUE)

domains_missense_idomain_sub <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/idomain_results_missense.csv", sep = ""), header = TRUE)
# domains_missense_idomain_sub <- read.csv(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/idomain_results_missense.csv", sep = ""), header = TRUE)

protein_ids_missense_idomain <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/protein_ids_missense_idomain.csv", sep = ""), header = TRUE)
swissprot_ids_missense_idomain <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_idomain.csv", sep = ""), header = TRUE)
#protein_ids_missense_idomain <- read.csv(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/protein_ids_missense_idomain.csv", sep = ""), header = TRUE)
#swissprot_ids_missense_idomain <- read.csv(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/swissprot_ids_missense_idomain.csv", sep = ""), header = TRUE)


############################################################
# I-BINDING POSITION: KEEP ONLY PROTEINS WITH A MUTATION ~IN~ 
# A DNA-BINDING POSITION FROM INTERACDOME OR CANBIND 
############################################################
#' For every protein, see if the any of the mutations fall within a specific 
#' predicted binding position. Compile all those that do into a new, final DF and return
#' @param domain_df the I-Domain data frame containing all the proteins that have 
#' passed through filtering so far
#' @param interacdome_binding_positions_df a binding position DF from InteracDome 
#' (all domains and predicted binding sites)
#' @param canbind_binding_site_df a CanBind data frame that has the positions of 
#' each binding site for each protein
#' @param concavity_df (OPT) a ConCavity DF that has the positions of each binding 
#' site for each protein
get_ibindingpos_subset <- function(domain_df, interacdome_df, canbind_df, concavity_df) {
  swissprot_ids <- unique(domain_df$Swissprot)
  
  # Loop through all proteins 
  prot_dfs <- lapply(1:length(swissprot_ids), function(i) {
      id <- swissprot_ids[i]  # Swissprot ID of current protein
      print(id)
      print(paste(i, paste("/", length(swissprot_ids))))
      
      # Subset domain DF to only look at this protein 
      domain_sub <- domain_df[domain_df$Swissprot == id,]
      
      # Subset the InteracDome, CanBind, and ConCavity DFs to also only look at 
      # this protein and get the AA binding positions. Combine these positions 
      # across sources.
      interacdome_binding_pos <- interacdome_df[grep(id, interacdome_df$Swissprot, 
                                                     ignore.case = FALSE, fixed = TRUE), 'Binding.Pos']
      interacdome_binding_pos <- as.numeric(unlist(lapply(interacdome_binding_pos, 
                                                          function(x) unlist(strsplit(x, split = ",", 
                                                                                      fixed = TRUE)))))
      interacdome_binding_pos <- sort(unique(interacdome_binding_pos), decreasing = FALSE)
      
      canbind_binding_pos <- canbind_df[grep(id, canbind_df$Swissprot, ignore.case = FALSE, 
                                             fixed = TRUE), 'AA_position']
      canbind_binding_pos <- sort(unique(canbind_binding_pos), decreasing = FALSE)
      
      try({
        concavity_binding_pos <- concavity_df[grep(id, concavity_df$Swissprot, 
                                                   ignore.case = FALSE, fixed = TRUE),'AA_position']
        concavity_binding_pos <- sort(unique(concavity_binding_pos), decreasing = FALSE)
      }, silent = TRUE)
      
      # Loop through all patients with mutations in this protein, and check if 
      # those mutations are also in a binding position
      mut_dfs <- lapply(1:nrow(domain_sub), function(j) {
        # Get the mutation position
        mut_pos <- domain_sub[j, 'Mut.Pos']
        print(mut_pos)
        
        # Also get the patient ID
        patient_id <- domain_sub[j, 'Patient']
        
        # Call helper function to check if mutations are in a binding position
        rows_interac <- check_if_in_bp(mut_pos, interacdome_binding_pos, patient_id, "InteracDome", id)
        rows_canbind <- check_if_in_bp(mut_pos, canbind_binding_pos, patient_id, "CanBind", id)
        mut_df <- rbind(rows_interac, rows_canbind)
        try ({
          rows_concavity <- check_if_in_bp(mut_pos, concavity_binding_pos, patient_id, "ConCavity", id)
          mut_df <- rbind(mut_df, rows_concavity)
        }, silent = TRUE) 

        print(mut_df)
        return(mut_df)
      })
      
      if(!is.null(mut_dfs)) {
        prot_df <- do.call(rbind, mut_dfs) 
        return(prot_df)
      }
  })
  
  ibindingpos_df <- do.call(rbind, prot_dfs)
  colnames(ibindingpos_df) <- c("Protein.ID", "Pred.Binding.Pos", "Source", "Mut.Pos", "Patient")
  
  return(as.data.frame(ibindingpos_df))
}

#' Given a mutation and a set of binding positions, checks if the mutation is in 
#' the binding position. If so, creates a new row using the mutation position, the 
#' binding positions, the patient ID, and the name of the tool that listed that binding 
#' position ("label").
#' @param mut_pos the position of the mutation in question
#' @param binding_pos a vector of binding positions within that protein
#' @param patient_id the 4-digit TCGA ID of the patient in question
#' @param label the source of the binding positions ('InteracDome', 'CanBind', or
#' 'ConCavity)
#' @param id the swissprot ID of the protein in question
check_if_in_bp <- function(mut_pos, binding_pos, patient_id, label, id) {
  # As long as this mutation is in the coding region...
  if (!mut_pos == "-") {
    # Is this mutation in the given binding positions?
    if (mut_pos %fin% binding_pos) {
      # Yes! Add this entry, with the mutation position added as well
      row_to_add <- c(id, paste(binding_pos, collapse = ","), label, mut_pos, patient_id)
      print(row_to_add)
      return(row_to_add)
    }
  }
  return(NULL)
}


domains_missense_ibindingpos_sub <- get_ibindingpos_subset(interacdome_binding_positions_df,
                                                           canbind_df_sub, domains_missense_idomain_sub,
                                                           concavity_df_sub)

# For nucleic acids, not using ConCavity
#domains_missense_ibindingpos_sub <- get_ibindingpos_subset(interacdome_binding_positions_df,
#                                                           canbind_df_sub, domains_missense_idomain_sub)

domains_missense_ibindingpos_sub <- distinct(domains_missense_ibindingpos_sub)


############################################################
# GET UNIQUE PROTEIN IDs AND WRITE TO FILES
############################################################
# Get the number of unique proteins in this DF:
swissprot_ids_missense_ibindingpos <- unique(domains_missense_ibindingpos_sub$Protein.ID)
length(swissprot_ids_missense_ibindingpos) 
  # All Ligands: BRCA: 810; Pan- Cancer: 4610
  # Nuc. Acids Only: BRCA: 55; Pan-Cancer: 499
  
# Write the full table to a CSV
write.csv(domains_missense_ibindingpos_sub, file = paste(path, "Saved Output Data Files/BRCA/Mutation/ibindingpos_results_missense.csv", sep = ""))
# write.csv(domains_missense_ibindingpos_sub, file = paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/ibindingpos_results_missense.csv", sep = ""))

# Write the unique protein list to a CSV
write.csv(swissprot_ids_missense_ibindingpos, file = paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_ibindingpos.csv", sep = ""))
#write.csv(swissprot_ids_missense_ibindingpos, file = paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/swissprot_ids_missense_ibindingpos.csv", sep = ""))

# Use these IDs to further subset the MAF file
maf_subset_df_missense_ibp <- maf_subset_df_missense_idom[maf_subset_df_missense_idom$SWISSPROT %in% swissprot_ids_missense_ibindingpos,]
write.csv(maf_subset_df_missense_ibp, paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/maf_file_df_missense_ibp.csv", sep = ""))
#write.csv(maf_subset_df_missense_ibp, paste(path, "Input Data Files/Pan-Cancer/Somatic_Mut_Data/maf_file_df_missense_ibp.csv", sep = ""))

#
#
# Read these files back if needed
domains_missense_ibindingpos_sub <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/ibindingpos_results_missense.csv", sep = ""), header = TRUE)
#domains_missense_ibindingpos_sub <- read.csv(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/ibindingpos_results_missense.csv", sep = ""), header = TRUE)
swissprot_ids_missense_ibindingpos <- read.csv(paste(path, "Saved Output Data Files/BRCA/Mutation/swissprot_ids_missense_ibindingpos.csv", sep = ""))
#swissprot_ids_missense_ibindingpos <- read.csv(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/swissprot_ids_missense_ibindingpos.csv", sep = ""))
maf_subset_df_missense_ibp <- read.csv(paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/maf_file_df_missense_ibp.csv", sep = ""))
# maf_subset_df_missense_ibp <- read.csv(paste(path, "Input Data Files/Pan-Cancer/Somatic_Mut_Data/maf_file_df_missense_ibp.csv", sep = ""))
