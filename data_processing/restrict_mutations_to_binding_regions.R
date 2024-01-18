############################################################
### Limit to Mutations Within Particular Binding Domains or Positions
### ASSOCIATED PUBLICATION INFORMATION
############################################################

library(data.table)
library(GenomicRanges)
library(TRONCO)
library(stringr)
library(dplyr)
library("RColorBrewer")

# Contains additional, unpublished resources to limit mutations considered
# to those that fall in particular ligand-binding regions (domains) or positions

# Local PATH to directory containing MAF data file
PATH <- getwd()
OUTPUT_PATH <- paste0(getwd(), "Output/Mutation/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "Input Data Files/all_genes_id_conv.csv"), 
                              header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
# IMPORT DATA FROM INTERACDOME, CANBIND, & CONCAVITY
############################################################
# 1. InteracDome (https://doi.org/10.1093/nar/gky1224)

# All Ligands
interacdome_conf <- read.csv(paste0(PATH, "Mutation/InteracDome/
                                    InteracDome_v0.3-confident.csv"), 
                             header = T, check.names = F)
binding_domains_ids_noPF <- read.csv(paste0(PATH, "Mutation/InteracDome/
                                            binding_domains_noPF.csv"), 
                                     header = T, check.names = F)[,2]

interacdome_binding_positions_df <- read.csv(paste0(PATH, "Mutation/InteracDome/
                                                    binding_positions_DF_0.csv"), 
                                             header = T, check.names = F, 
                                             row.names = 1)

# Nucleic Acids Only
interacdome_conf <- read.csv(paste0(PATH, "Mutation/InteracDome/
                                    InteracDome_v0.3-confident_nucacids.csv"), 
                             header = T, check.names = F)
binding_domains_ids_noPF <- read.csv(paste0(PATH, "Mutation/InteracDome/
                                            binding_domains_nucacids_noPF.csv"), 
                                     header = T, check.names = F)[,2]
interacdome_binding_positions_df <- read.csv(paste0(PATH, "Mutation/InteracDome/
                                                    binding_positions_DF_nucacids_0.csv"), 
                                             header = T, check.names = F)


# 2. CanBind (10.1093/nar/gkt1305)
# All Ligands

# Threshold 0.5
canbind_df_sub <- read.csv(paste0(PATH, "Mutation/CanBind/
                                  canbind_dataframe_labeled_0.5.csv"), 
                           header = T)
canbind_swissprot_ids <- unique(canbind_df_sub$Swissprot)
canbind_binding_site_df <- read.csv(paste0(PATH, "Mutation/CanBind/
                                           canbind_binding_ranges_0.5.csv"), 
                                    header = T)

# No threshold
canbind_df_sub <- read.csv(paste0(PATH, "Mutation/CanBind/
                                  canbind_dataframe_labeled_0.csv"), 
                           header = T)
canbind_swissprot_ids <- unique(unlist(lapply(canbind_df_sub$Swissprot, function(x) 
  unlist(strsplit(x, ";")))))
canbind_binding_site_df <- read.csv(paste0(PATH, "Mutation/CanBind/
                                           canbind_binding_ranges_0.csv"), 
                                    header = T)

# Nucleic Acids Only
canbind_df_sub <- read.csv(paste0(PATH, "Mutation/CanBind/
                                  canbind_dataframe_labeled_nucacids_0.csv"), 
                           header = T)
canbind_swissprot_ids <- unique(canbind_df_sub$Swissprot)
canbind_binding_site_df <- read.csv(paste0(PATH, "Mutation/CanBind/
                                           canbind_binding_ranges_nucacids_0.csv"), 
                                    header = T)

# 3. ConCavity 
concavity_df_sub <- read.csv(paste0(PATH, "Mutation/ConCavity/
                                    concavity_dataframe_labeled_0.1.csv"), 
                             header = T)
concavity_swissprot_ids <- unique(unlist(lapply(concavity_df_sub$Swissprot, 
                                                function(x) 
                                                  unlist(strsplit(x, ";")))))
concavity_binding_site_df <- read.csv(paste0(PATH, "Mutation/ConCavity/
                                             concavity_binding_ranges_0.1.csv"), 
                                      header = T)

# NOTES: ConCavity does not have ligand information; if restriction to nucleic 
# acids only, do not use ConCavity information


############################################################
# IMPORT PREPROCESSED MAF FILE
############################################################
# Has been filtered for hypermutators in process_mutation_data.R
maf_file_df <- read.csv(paste0(PATH, "Mutation/Individual_PanCancer_MAFs/
                              maf_file_hypermut_filt.csv"))


############################################################
# GET PROTEOME INFORMATION
############################################################
# Get all proteins in the human proteome; canonical only
# read.fasta returns a list of character sequences (when as.string = T)
proteome_filename <- paste0(PATH, "Proteome/uniprot-reviewed_yes+AND+proteome_up000005640.fasta")
proteome <- read.fasta(proteome_filename, seqtype = c("AA"), 
                       as.string = T)
# Length of proteome: 20353


############################################################
# FILTER PROTEINS BY PRESENCE OF MUTATIONS
############################################################
# Determine shared ID between the MAF and FASTA files; keep only proteome 
# proteins at the intersection of these two files
# SWISSPROT ID (column 69 in MAF)
maf_swissprot_ids <- maf_file_df$SWISSPROT

#' Takes a proteome and set of protein Swissprot or Hugo IDs from a given MAF 
#' file; returns a subsetted proteome DF with only overlapping proteins
#' @param proteome a proteome from Uniprot or elsewhere, in FASTA format
#' @param maf_file_ids vector of SWISSPROT or HUGO IDs from a given MAF file
#' @param id_type a string denoting the ID type (either 'swissprot' or 'hugo')
get_proteome_subset <- function(proteome, maf_file_ids, id_type) {
  # Split up the maf file ids, if needed
  maf_file_ids <- unique(unlist(lapply(maf_file_ids, function(x) 
    unlist(strsplit(x, ";", fixed = T)))))
  
  proteome_subset <- c()
  for (i in 1:length(proteome)) {
    entry <- proteome[i]
    entry_id <- ""
    if (id_type == 'swissprot') {
      entry_id <- unlist(strsplit(attributes(proteome[i])[[1]], "|", 
                                  fixed = T))[2]
    } else if (id_type == 'hugo') {
      entry_id <- unlist(strsplit(unlist(strsplit(attributes(proteome[i])[[1]], 
                                                  "|", fixed = T))[3], "_", 
                                  fixed =  T))[1]
    } else {
      print(paste("Unrecognized ID type:", id_type))
    }
    if (entry_id %fin% maf_file_ids) {
      proteome_subset <- c(proteome_subset, entry)
    }
  }
  return(proteome_subset)
}

proteome_subset <- get_proteome_subset(proteome, maf_swissprot_ids, 'swissprot')
length(unique(proteome_subset))  # check the number of protein sequences remaining

# Write proteome subset to a FASTA
write.fasta(proteome_subset, names = names(proteome_subset), as.string = T, 
            file.out = paste0(PATH, "Proteomeproteome_subset.csv"))

# Alternatively, if we've already done this, just read back proteome subset 
# from FASTA
proteome_subset_filename <- paste0(PATH, "Proteome/proteome_subset.csv")
proteome_subset <- read.fasta(proteome_subset_filename, seqtype = c("AA"), 
                              as.string = T)


############################################################
# GET DOMAIN INFORMATION FOR ALL PROTEINS IN THE HUMAN PROTEOME
############################################################
# Keep only proteins that have at least one confident binding domain

# Export the list of remaining proteins back to FASTA format so that it can be 
# uploaded to Batch CD-Search, a tool that accepts protein sequences in FASTA 
# format (max 4000 entries/ batch) or to HMMER 
# Batch CD-Search link: https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
# NOTE: Uniprot's InterPro tool is another option for domain annotation

prep_cdsearch_files(proteome_subset, "missense", "tcga")   # see batch_cdsearch_helper_functions.R
prep_cdsearch_files(proteome_subset, "missense", "metabric")   # see batch_cdsearch_helper_functions.R


############################################################
# IMPORT RESULTS FROM BATCH CD-SEARCH 
############################################################
domains <- get_cdsearch_res("missense", "tcga")      # see batch_cdsearch_helper_functions.R
domains <- get_cdsearch_res("missense", "metabric")      # see batch_cdsearch_helper_functions.R

# Write to a file 
write.csv(domains, paste0(PATH, "CD-Search/missense_cdsearch_domains_all.csv"))
write.csv(domains, paste0(PATH, "CD-Search/missense_cdsearch_domains_all_metabric.csv"))

# If we have already done the above once, simply import the saved domains files
domains <- read.csv(paste0(PATH, "CD-Search/missense_cdsearch_domains_all.csv"), 
                    header = T, check.names = F)
domains <- read.csv(paste0(PATH, "CD-Search/missense_cdsearch_domains_all_metabric.csv"), 
                    header = T, check.names = F)


############################################################
# EXTRACT PFAM IDs FROM BATCH CD-SEARCH RESULTS
############################################################
# NOTE: Query #s in the form Q#N - XXXXXXXX, where XXXXXXXX is the sequence ID 
# or first 15 chars of FASTA line
# Get Pfam domains from Batch CD-Search Results

#' Subset domains data frame to only include pfam domains
#' @param domain_df a domain data frame to be subsetted
subset_domains_to_pfam <- function(domain_df) {
  domain_df <- domain_df[grep("pfam", domain_df$Accession, ignore.case = T),]
}
domains <- subset_domains_to_pfam(domains)
accessions <- domains$Accession

# Get only the numeric portions of the accessions, and only those from the 
# Pfam database (InteracDome only uses Pfam domains, so we will not get matches 
# for SMART or COG)
accessions_numeric <- regmatches(accessions, regexpr("[[:digit:]]+", accessions)) 


############################################################
# GET DOMAIN OVERLAP BETWEEN BATCH CD-SEARCH AND INTERACDOME
############################################################
# Which of these domains are confidently binding? Look for overlap between the 
# domain hits from CD-Search and the DNA-binding domains from InteracDome
intersecting_domain_acc <- intersect(accessions_numeric, 
                                     binding_domains_ids_noPF)
length(intersecting_domain_acc)


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
#' @param canbind_swissprot_ids a vector of swissprot IDs with confident 
#' interacting regions from CanBind
#' @param concavity_swissprot_ids OPT: a vector of swissprot IDs with confident 
#' interacting regions from ConCavity
get_iprotein_subset <- function(domains, intersecting_domain_acc, 
                                canbind_swissprot_ids, 
                                concavity_swissprot_ids) {
  # Fix the "Query" column name (was i..Query)
  if("ï»¿Query" %fin% colnames(domains)) {
    colnames(domains)[which(colnames(domains) == "ï»¿Query")] <- "Query"
  }
  
  # Fix the query column to get rid of query number (we don't really need this)
  domains$Query <- as.character(unlist(lapply(domains$Query, function(x) 
    unlist(strsplit(x, " ", fixed = T))[3])))
  
  # Subset by InteracDome domains
  domains_interacdome_sub <- domains[
    as.character(regmatches(domains$Accession, 
                            gregexpr("[[:digit:]]+", domains$Accession))) %fin% 
      intersecting_domain_acc,]
  
  # Subset by CanBind proteins with interaction regions
  queries <- domains$Query
  swissprot_ids <- unlist(sapply(queries, function(x) 
    {unlist(strsplit(x, "|", fixed = T))[2]}))
  matching_indices_canbind <- which(canbind_swissprot_ids %fin% swissprot_ids)
  domains_canbind_sub <- domains[matching_indices_canbind,]
  
  # Subset by Concavity proteins with interaction regions, if provided
  if (!missing(concavity_swissprot_ids)) {
    matching_indices_concavity <- which(concavity_swissprot_ids %fin% 
                                          swissprot_ids)
    domains_concavity_sub <- domains[matching_indices_concavity,]
    
    # Take the union of these data frames to include all those that remain 
    # using either method
    dupl_rows_pt1 <- which(!is.na(match(domains_interacdome_sub$Accession, 
                                        domains_canbind_sub$Accession)))
    domains_sub_pt1 <- rbind(domains_interacdome_sub, 
                             domains_canbind_sub[-dupl_rows_pt1,])
    dupl_rows_pt2 <- !is.na(match(domains_sub_pt1$Accession, 
                                  domains_concavity_sub$Accession))
    domains_sub_pt2 <- rbind(domains_sub_pt1, 
                             domains_concavity_sub[-dupl_rows_pt2,])
    domains_sub <- domains_sub_pt2[!is.na(domains_sub_pt2$Query),]
  }
  else {
    dupl_rows <- which(!is.na(match(domains_interacdome_sub$Accession, 
                                    domains_canbind_sub$Accession)))
    domains_sub <- rbind(domains_interacdome_sub, 
                         domains_canbind_sub[-dupl_rows,])
    domains_sub <- domains_sub[!is.na(domains_sub$Query),]
  }
  
  return(domains_sub)
}


domains_iprotein_sub <- get_iprotein_subset(domains, intersecting_domain_acc, 
                                            canbind_swissprot_ids, 
                                            concavity_swissprot_ids) 

# No ConCavity version for nucleic acids only
domains_iprotein_sub <- get_iprotein_subset(domains, intersecting_domain_acc, 
                                            canbind_swissprot_ids) 


############################################################
# ADD SEQUENCE INFORMATION TO SUBSETTED PROTEIN DF
############################################################
#' Takes data frame of domain information from Batch-CD search, along with
#' a proteome, and extracts the AA acid sequence for given query, which it adds 
#' as an additional column to the data frame. Returns this data frame.
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

domains_iprotein_sub <- merge_proteome_and_domains(domains_iprotein_sub, 
                                                            proteome_subset)


############################################################
# GET UNIQUE PROTEIN IDs AND WRITE TO FILES
############################################################
# Examine how many actual unique accessions we have
length(unique(domains_iprotein_sub$Accession))  

# Write this table to a CSV
write.csv(domains_iprotein_sub, paste0(PATH, "Mutation/iprotein_results.csv"))

# If already done this, read them back from CSV
domains_iprotein_sub <- read.csv(paste0(PATH, "Mutation/iprotein_results.csv"), 
                                 header = T, check.names = F)


# Extract the Swissprot IDs and general names for all the remaining queries

#' Return all the SWISSPROT IDs from a domain data frame
#' @param domain_df the domain data frame to find IDs in 
extract_swissprot_ids <- function(domain_df) {
  swissprot_ids <- unique(unlist(lapply(domain_df$Query, function(x) 
    unlist(strsplit(x,"|", fixed = T))[2])))
  return(swissprot_ids)
}

#' Return all the protein names from a domain data frame
#' @param domain_df the domain data frame to list the gene names from
extract_protein_names <- function(domain_df) {
  protein_ids <- c()
  for (i in 1:nrow(domain_df)) {
    split_id <- unlist(strsplit(domain_df$Query[i],"|", fixed = T))
    protein_ids <- c(protein_ids, unlist(strsplit(split_id[3], "_", fixed = T))[1])
  }
  return(unique(protein_ids))
}

swissprot_ids_iprotein <- extract_swissprot_ids(domains_iprotein_sub)
protein_ids_iprotein <- extract_protein_names(domains_iprotein_sub)  
length(protein_ids_iprotein)

# Write unique protein IDs to a file
write.csv(protein_ids_iprotein, paste0(PATH, "Mutation/protein_ids_iprotein.csv"))
write.csv(swissprot_ids_iprotein, paste0(PATH, "Mutation/swissprot_ids_iprotein.csv"))

# Read back if necessary
protein_ids_iprotein <- read.csv(paste0(PATH, "Mutation/protein_ids_iprotein.csv"))[,2]
swissprot_ids_iprotein <- read.csv(paste0(PATH, "Mutation/swissprot_ids_iprotein.csv"))[,2]

# Use these IDs to subset the MAF file
maf_subset_df <- maf_file_df[maf_file_df$SWISSPROT %fin% swissprot_ids_iprotein,]

# For METABRIC, have to add the SWISSPROT ID and then do the subsetting (does 
# not work properly with Hugo symbol)
maf_file_df$SWISSPROT <- unlist(lapply(maf_file_df$Hugo_Symbol, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 
                                        'uniprot_gn_id'])), collapse = ";")))
rows_to_keep <- unlist(lapply(1:nrow(maf_file_df), function(i) {
  prot_ids <- unlist(strsplit(maf_file_df$SWISSPROT[i], ";", fixed = T))
  truth_vect <- unlist(lapply(prot_ids, function(id) 
    ifelse(id %fin% swissprot_ids_iprotein, T, F)))
  if(T %fin% truth_vect) {return(i)}
  else {return(NA)}
}))
rows_to_keep <- rows_to_keep[!is.na(rows_to_keep)]
maf_subset_df <- maf_file_df[rows_to_keep,]

write.csv(maf_subset_df, paste0(PATH, "Mutation/maf_file_df_iprotein.csv"))

# If already done, read back the file
maf_subset_df <- read.csv(paste0(PATH, "Mutation/maf_file_df_iprotein.csv"), 
                          header = T, check.names = F)


############################################################
# OPTIONAL: ADD PATIENTS AND MUTATION POSITIONS TO DATA FRAME
############################################################
#' Add corresponding patient IDs that have mutations in each protein to the 
#' given I-Protein-level data frame using the information from the MAF file
#' @param iprotein_sub an I-Protein level data frame produced from the above
#' @param maf a MAF file with patient mutational information
#' @param dataset either 'tcga' or 'metabric'
add_patient_ids <- function(iprotein_sub, maf, dataset) {
  # For each unique protein, get the patients that are mutated in it
  full_prot_set <- unlist(lapply(iprotein_sub$Query, function(x) 
    unlist(strsplit(x, "|", fixed = T))[2]))
  unique_proteins <- unique(full_prot_set)
  
  patients_mut_per_prot <- lapply(unique_proteins, function(p) {
    # Subset the MAF to only the given protein
    maf_sub <- maf[grepl(p, maf$SWISSPROT),]
    
    # Get the patient and sample that have this protein mutated
    unique_patients <- unique(maf_sub$Tumor_Sample_Barcode)
    
    # If TCGA, split this ID apart to get just patient and sample XXXX-XX
    if(dataset == "tcga") {
      unique_patients <- unlist(lapply(unique_patients, function(x) 
        paste(unlist(strsplit(x, "-", fixed = T))[3:4], collapse = "-")))
    }
    # Collapse these patients into a semicolon-separated string
    patient_str <- paste(unique_patients, collapse = ";")
    return(patient_str)
  })
  names(patients_mut_per_prot) <- unique_proteins
  
  new_col <- lapply(full_prot_set, function(p) patients_mut_per_prot[[p]])
  iprotein_sub$Patient <- unlist(new_col)
  
  return(iprotein_sub)
}

domains_iprotein_sub <- add_patient_ids(domains_iprotein_sub, maf_subset_df, 
                                        "tcga")
domains_iprotein_sub <- add_patient_ids(domains_iprotein_sub, maf_subset_df, 
                                        "metabric")


write.csv(domains_iprotein_sub, paste0(PATH, "Mutation/iprotein_results.csv"))
write.csv(domains_iprotein_sub, paste0(PATH, "brca_metabric/Mutation/
                                       iprotein_results.csv"))


############################################################
# I-DOMAIN: KEEP ONLY PROTEINS WITH A MUTATION ~IN~ A DNA-
# BINDING DOMAIN REGION FROM INTERACDOME OR INTERACTION
# REGION FROM CANBIND/ CONCAVITY
############################################################
# We know the region of the domain in the "To" and "From" columns of the DF
# Now, we need to check whether, for this particular protein, there are mutations in that
# region in BRCA cancer patients

#' For every protein, see if the "Start" and "End" of any of the mutations falls 
#' within the domain region. Compile all those that do into a new, final DF and 
#' return
#' @param domain_df domain data frame from the CD-Batch Search Output, subsetted 
#' to include only InteracDome/ CanBind/ ConCavity domains
#' @param maf_df MAF-style data frame that will have the "Start" and "End" bases 
#' of the mutations of interest
#' @param canbind_binding_site_df a CanBind data frame that has the start and 
#' end positions of each binding site for each protein
#' @param concavity_binding_site_df OPT: a ConCavity data frame that has start 
#' and end positions of each binding site for each protein
get_idomain_subset <- function(domain_df, maf_df, canbind_binding_site_df, 
                               concavity_binding_site_df) {
  # Create a DF to hold the final results
  idomain_prot_df <- data.frame(matrix(nrow = 0, ncol = 9))
  
  # Get all the unique proteins in the domain DF
  unique_protein_ids <- unique(unlist(lapply(domain_df$Query, function(x) 
    unlist(strsplit(x, "|", fixed = T))[2])))
  
  for (id in unique_protein_ids) {
    print(paste(match(id, unique_protein_ids), 
                paste("/", length(unique_protein_ids))))
    
    # Subset the domain DF, MAF file, CanBind DF, and ConCavity DF (if provided) 
    # to look only at this protein. For each mutation in this protein, check if 
    # it is in a domain or binding region for CanBind and/or ConCavity
    domain_df_sub <- domain_df[grep(id, domain_df$Query),]
    maf_sub <- maf_df[maf_df$SWISSPROT == id,]
    canbind_ranges <- unique(canbind_binding_site_df[
      grepl(id, canbind_binding_site_df$Swissprot), 'Range.of.Binding.Site'])
    if (!missing(concavity_binding_site_df)) {
      concavity_ranges <- unique(concavity_binding_site_df[
        grepl(id, concavity_binding_site_df$Swissprot), 'Range.of.Binding.Site'])
      prot_df <- get_overlap_with_binding_domains(id, domain_df_sub, maf_sub, 
                                                  canbind_ranges, 
                                                  concavity_ranges)
    } else {
      prot_df <- get_overlap_with_binding_domains(id, domain_df_sub, 
                                                  maf_sub, canbind_ranges)
    }
    
    # Bind this protein's data frame to the main DF
    idomain_prot_df <- rbind(idomain_prot_df, prot_df)
  }
  colnames(idomain_prot_df) <- c("Swissprot", "From", "To", "Source", "Accession", 
                                 "Short Name", "Superfamily", "Mut.Pos", "Patient")
  print(head(idomain_prot_df))
  return(idomain_prot_df)
}

#' Helper function to get the overlap of mutations with 1. conversed binding 
#' domains or 2. CanBind or (optionally) ConCavity binding ranges. Returns a 
#' vector of rows to be added to the new data frame.
#' @param id swissprot ID of current protein
#' @param domain_df_sub the domain DF, subsetted to only this protein
#' @param maf_sub the MAF file, subsetted to only this protein
#' @param canbind_ranges the CanBind binding ranges data frame, subsetted to 
#' only this protein
#' @param concavity_ranges the ConCavity binding ranges data frame, subsetted to 
#' only this protein
get_overlap_with_binding_domains <- function(id, domain_df_sub, maf_sub, 
                                             canbind_ranges, concavity_ranges) {
  # Iterate over all the mutations in this given protein
  df_full <- lapply(maf_sub$Protein_position, function(x) {
    
    # For each mutation in this protein, get its AA position ("Protein_position", 
    # column 55) 
    mut_pos <- unlist(strsplit(x, "/", fixed = T))[1]
    print(mut_pos)
    
    # Make sure this mutation is in the coding region
    if (!mut_pos == "-") {
      # Also get the patient under consideration
      patient_id <- paste(unlist(strsplit(maf_sub$Tumor_Sample_Barcode[
        which(maf_sub$Protein_position == x)], "-", fixed = T))[3:4], 
        collapse = "-")
      
      # First, check if the mutation is in a conserved binding domain 
      domain_rows <- lapply(1:nrow(domain_df_sub), function(i) {
        start_dom <- as.numeric(domain_df_sub[i, 'From'])
        end_dom <- as.numeric(domain_df_sub[i, 'To'])
        
        if ((as.numeric(mut_pos) >= start_dom) & (as.numeric(mut_pos) <= end_dom)) {
          # Return the filled row
          # Columns to fill: "Swissprot", "From", "To", "Source", "Accession", 
          # "Short Name", "Superfamily", "Mut.Pos", "Patient"
          return(c(id, start_dom, end_dom, "Conserved Domain", 
                   domain_df_sub[i, 'Accession'], 
                   domain_df_sub[i, 'Short name'], 
                   domain_df_sub[i, 'Superfamily'], mut_pos, patient_id))
          # Or return(c(id, start_dom, end_dom, "Conserved Domain", 
          # domain_df_sub[i, 'Accession'], domain_df_sub[i, 'Short.name'], 
          # domain_df_sub[i, 'Superfamily'], mut_pos, patient_id)) 
          #, depending on headers
        }
      })
      domain_df <- do.call(rbind, domain_rows)
      
      # Then, check if the mutation is in a binding region for CanBind or 
      # ConCavity (if applicable)
      canbind_df <- get_canbind_and_concavity_rows(canbind_ranges, "CanBind", 
                                                   domain_df_sub, mut_pos, 
                                                   patient_id, id)
      #print(paste("CanBind df ", canbind_df))
      
      try({
        concavity_df <- get_canbind_and_concavity_rows(concavity_ranges, 
                                                       "ConCavity", domain_df_sub, 
                                                       mut_pos, patient_id, id)
      }, silent = T)
      #try({print(paste("ConCavity df ", concavity_df))})
      
      # Merge all the rows from these three sources and sort by protein ID
      df <- rbind(domain_df, canbind_df)
      try({df <- rbind(df, concavity_df)}, silent = T)
      
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
#' @param binding_ranges a set of binding ranges from either CanBind or ConCavity 
#' (start:end)
#' @param label a 'source' label (either 'CanBind' or 'ConCavity')
#' @param domain_df_sub a domain DF with information about mutations and domains
#' @param mut_pos the position of the mutation in question
#' @param patient_id the 4-digit patient TCGA ID of the patient in question
#' @param id the swissprot ID of the current protein
get_canbind_and_concavity_rows <- function(binding_ranges, label, domain_df_sub, 
                                           mut_pos, patient_id, id) {
  rows <- lapply(binding_ranges, function(y) {
    start_pos <- as.numeric(unlist(strsplit(y, ":", fixed = T))[1])
    end_pos <- as.numeric(unlist(strsplit(y, ":", fixed = T))[2])
    
    if ((as.numeric(mut_pos) >= start_pos) & (as.numeric(mut_pos) <= end_pos)) {
      # Return the filled row
      # Columns to fill: "Swissprot", "From", "To", "Source", "Accession", 
      # "Short Name", "Superfamily", "Mut.Pos", "Patient"
      return(c(id, start_pos, end_pos, label, NA, NA, NA, mut_pos, patient_id))
    }
  })
  df <- do.call(rbind, rows)
  return(df)
}


domains_idomain_sub <- get_idomain_subset(domains_iprotein_sub, maf_subset_df,
                                          canbind_binding_site_df, 
                                          concavity_binding_site_df)
# For nucleic acids, not using ConCavity
domains_idomain_sub <- get_idomain_subset(domains_iprotein_sub, maf_subset_df,
                                                   canbind_binding_site_df)

# Delete any duplicate rows
domains_idomain_sub <- distinct(domains_idomain_sub)


############################################################
# GET UNIQUE PROTEIN IDs AND WRITE TO FILES
############################################################
# Get the number of unique domains in this DF:
length(unique(domains_idomain_sub$Accession)) 

# Get the Swissprot IDs
swissprot_ids_idomain <- unique(domains_idomain_sub$Swissprot)
print(length(swissprot_ids_idomain)) 

# Write the full table to a CSV
write.csv(domains_idomain_sub, paste0(PATH, "Mutation/idomain_results.csv"))

# Write the unique protein list to a CSV
write.csv(swissprot_ids_idomain, paste0(PATH, "Mutation/swissprot_ids_idomain.csv"))

# Use these IDs to further subset the MAF file
maf_subset_df_idom <- maf_subset_df[maf_subset_df$SWISSPROT %fin% 
                                      swissprot_ids_idomain,]
write.csv(maf_subset_df_idom, paste0(PATH, "Mutation/maf_file_df_idomain.csv"))


# 
# 
# 
# Read these files back if needed
maf_subset_df <- read.csv(paste0(PATH, "Mutation/maf_file_df_idomain.csv"), 
                          header = T, check.names = F)

domains_idomain_sub <- read.csv(paste0(PATH, "Mutation/idomain_results.csv"), 
                                header = T, check.names = F)

swissprot_ids_idomain <- read.csv(paste0(PATH, "Mutation/swissprot_ids_idomain.csv"), 
                                  header = T, check.names = F)


############################################################
# I-BINDING POSITION: KEEP ONLY PROTEINS WITH A MUTATION ~IN~ 
# A DNA-BINDING POSITION FROM INTERACDOME OR CANBIND 
############################################################
#' For every protein, see if the any of the mutations fall within a specific 
#' predicted binding position. Compile all those that do into a new, final DF and 
#' return
#' @param domain_df the I-Domain data frame containing all the proteins that have 
#' passed through filtering so far
#' @param interacdome_binding_positions_df a binding position DF from InteracDome 
#' (all domains and predicted binding sites)
#' @param canbind_binding_site_df a CanBind data frame that has the positions of 
#' each binding site for each protein
#' @param concavity_df (OPT) a ConCavity DF that has the positions of each 
#' binding site for each protein
get_ibindingpos_subset <- function(domain_df, interacdome_df, canbind_df, 
                                   concavity_df) {
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
                                                   ignore.case = F, fixed = T), 
                                              'Binding.Pos']
    interacdome_binding_pos <- as.numeric(
      unlist(lapply(interacdome_binding_pos, function(x) 
        unlist(strsplit(x, split = ",",  fixed = T)))))
    interacdome_binding_pos <- sort(unique(interacdome_binding_pos), 
                                    decreasing = F)
    
    canbind_binding_pos <- canbind_df[grep(id, canbind_df$Swissprot, 
                                           ignore.case = F, fixed = T), 
                                      'AA_position']
    canbind_binding_pos <- sort(unique(canbind_binding_pos), decreasing = F)
    
    try({
      concavity_binding_pos <- concavity_df[grep(id, concavity_df$Swissprot, 
                                                 ignore.case = F, fixed = T),
                                            'AA_position']
      concavity_binding_pos <- sort(unique(concavity_binding_pos), 
                                    decreasing = F)
    }, silent = T)
    
    # Loop through all patients with mutations in this protein, and check if 
    # those mutations are also in a binding position
    mut_dfs <- lapply(1:nrow(domain_sub), function(j) {
      # Get the mutation position
      mut_pos <- domain_sub[j, 'Mut.Pos']
      print(mut_pos)
      
      # Also get the patient ID
      patient_id <- domain_sub[j, 'Patient']
      
      # Call helper function to check if mutations are in a binding position
      rows_interac <- check_if_in_bp(mut_pos, interacdome_binding_pos, 
                                     patient_id, "InteracDome", id)
      rows_canbind <- check_if_in_bp(mut_pos, canbind_binding_pos, patient_id, 
                                     "CanBind", id)
      mut_df <- rbind(rows_interac, rows_canbind)
      try ({
        rows_concavity <- check_if_in_bp(mut_pos, concavity_binding_pos, 
                                         patient_id, "ConCavity", id)
        mut_df <- rbind(mut_df, rows_concavity)
      }, silent = T) 
      
      print(mut_df)
      return(mut_df)
    })
    
    if(!is.null(mut_dfs)) {
      prot_df <- do.call(rbind, mut_dfs) 
      return(prot_df)
    }
  })
  
  ibindingpos_df <- do.call(rbind, prot_dfs)
  colnames(ibindingpos_df) <- c("Protein.ID", "Pred.Binding.Pos", "Source", 
                                "Mut.Pos", "Patient")
  
  return(as.data.frame(ibindingpos_df))
}

#' Given a mutation and a set of binding positions, checks if the mutation is in 
#' the binding position. If so, creates a new row using the mutation position, 
#' the binding positions, the patient ID, and the name of the tool that listed 
#' that binding position ("label").
#' @param mut_pos the position of the mutation in question
#' @param binding_pos a vector of binding positions within that protein
#' @param patient_id the 4-digit TCGA ID of the patient in question
#' @param label the source of the binding positions ('InteracDome', 'CanBind', 
#' or 'ConCavity')
#' @param id the swissprot ID of the protein in question
check_if_in_bp <- function(mut_pos, binding_pos, patient_id, label, id) {
  # As long as this mutation is in the coding region...
  if (!mut_pos == "-") {
    # Is this mutation in the given binding positions?
    if (mut_pos %fin% binding_pos) {
      # Yes! Add this entry, with the mutation position added as well
      row_to_add <- c(id, paste(binding_pos, collapse = ","), label, 
                      mut_pos, patient_id)
      print(row_to_add)
      return(row_to_add)
    }
  }
  return(NULL)
}


domains_ibindingpos_sub <- get_ibindingpos_subset(domains_idomain_sub,
                                                  interacdome_binding_positions_df,
                                                  canbind_df_sub, concavity_df_sub)

# For nucleic acids, not using ConCavity
domains_ibindingpos_sub <- get_ibindingpos_subset(interacdome_binding_positions_df,
                                                  canbind_df_sub, 
                                                  domains_idomain_sub)

domains_ibindingpos_sub <- distinct(domains_ibindingpos_sub)


############################################################
# GET UNIQUE PROTEIN IDs AND WRITE TO FILES
############################################################
# Get the number of unique proteins in this DF:
swissprot_ids_ibindingpos <- unique(domains_ibindingpos_sub$Protein.ID)
length(swissprot_ids_ibindingpos) 

# Write the full table to a CSV
write.csv(domains_ibindingpos_sub, 
          paste0(PATH, "Mutation/ibindingpos_results.csv"))

# Write the unique protein list to a CSV
write.csv(swissprot_ids_ibindingpos, 
          paste0(PATH, "Mutation/swissprot_ids_ibindingpos.csv"))

# Use these IDs to further subset the MAF file
maf_subset_df_ibp <- maf_subset_df_idom[maf_subset_df_idom$SWISSPROT %fin% 
                                          swissprot_ids_ibindingpos,]
write.csv(maf_subset_df_ibp, paste0(PATH, "Mutation/maf_file_df_ibp.csv"))

#
#
# Read these files back if needed
domains_ibindingpos_sub <- read.csv(paste0(PATH, "Mutation/ibindingpos_results.csv"), 
                                    header = T, check.names = F)
swissprot_ids_ibindingpos <- read.csv(paste0(PATH, "Mutation/swissprot_ids_ibindingpos.csv"))
maf_subset_df_ibp <- read.csv(paste0(PATH, "Mutation/maf_file_df_ibp.csv"))

