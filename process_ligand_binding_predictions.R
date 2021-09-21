############################################################
### Process InteracDome, CanBind, & ConCavity Data
### Written By: Sara Camilli, November 2020
############################################################

library(stringr)
library(dplyr)

# This file takes ligand-binding predictions from three in-house sources:
  # 1. InteracDome "confident" binding domain CSV
  # 2. CanBind track (from PertInInt)
  # 3. ConCavity track (from PertInInt)
# and extracts information about ligand binding a) domains and b) positions

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/"
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)

############################################################
# GET INTERACDOME INTERACTION INFORMATION
############################################################
### OPTION 1: USE INTERACDOME FILE FROM WEBSITE WITH CONFIDENT INTERACTIONS ###
# Upload InteracDome file containing domains with confident DNA interactions
# InteracDome : domain-based interaction pockets
interacdome_conf <- read.csv(paste(path, "Input Data Files/InteracDome/InteracDome_v0.3-confident.csv", sep = ""))   # Remove comment rows before reading in

# OPT: Uncomment to select only DNA-interacting confident domains/ RNA-Interacting Domains; use the nucacids name throughout code below to refer to this DF
interacdome_conf_nucacids <- interacdome_conf[(interacdome_conf$ligand_type == "DNA_" |
                                        interacdome_conf$ligand_type == "DNABASE_" | 
                                        interacdome_conf$ligand_type == "DNABACKBONE" |
                                        interacdome_conf$ligand_type == "RNA_" |
                                        interacdome_conf$ligand_type == "RNABASE_" |
                                        interacdome_conf$ligand_type == "RNABACKBONE_" |
                                        interacdome_conf$ligand_type == "NUCACID_"),]
write.csv(interacdome_conf_nucacids, paste(path, "Input Data Files/InteracDome/InteracDome_v0.3-confident_nucacids.csv", sep = ""))

binding_domains <- unique(interacdome_conf$pfam_id)
binding_domains <- binding_domains[!is.na(binding_domains)]  # Remove any NA values
write.csv(binding_domains, paste(path, "Saved Output Data Files/InteracDome/binding_domains.csv", sep = ""))

# OPT: In order to match these domain names to CD-Search Pfam IDs, remove the "PF" from the ID
binding_domains_ids <- grep('PF', unlist(strsplit(binding_domains, "_", fixed = TRUE)), value = TRUE)
binding_domains_ids_noPF <- unique(unlist(lapply(binding_domains_ids, FUN = function(id) str_remove(id, "PF"))))
write.csv(binding_domains_ids_noPF, paste(path, "Saved Output Data Files/InteracDome/binding_domains_noPF.csv", sep = ""))
# Read back if already done
binding_domains_ids_noPF <- read.csv(paste(path, "Saved Output Data Files/InteracDome/binding_domains_noPF.csv", sep = ""))[,2]

#' Create a data frame that contains the binding positions for each protein domain
#' using input from InteracDome
#' @param interacdome InteracDome file with protein-ligand interaction information
#' @param binding_domains_ids_noPF a vector of ligand-binding domains, cropped to 
#' exclude the "PF" that denotes that the domain is from Pfam
#' @param threshold_interacdome a threshold for the InteracDome binding frequency 
#' to indicate a significant interaction (above threshold)
get_interacdome_binding_pos <- function(interacdome, binding_domains_ids_noPF, threshold_interacdome) {
  binding_positions_df <- data.frame(matrix(ncol = 1, nrow = length(binding_domains_ids_noPF)))
  colnames(binding_positions_df) <- c("Binding.Pos")
  
  for (i in 1:length(binding_domains_ids_noPF)) {
    # Get the binding frequencies for this domain
    domain <- binding_domains_ids_noPF[i]
    binding_freq <- interacdome_conf[unlist(lapply(interacdome$pfam_id, 
                                                   function(x) str_remove(unlist(strsplit(x, "_", fixed = TRUE))[1], "PF"))) == 
                                       domain, 'binding_frequencies']
    # Extract the binding positions that have a frequency greater than the given threshold
    binding_pos <- unlist(lapply(binding_freq, FUN = function(x) which(as.numeric(unlist(strsplit(x, ",", fixed = TRUE)))
                                                                       > threshold_interacdome)))
    # Often, these positions are in the same domain for different ligands. If we have already restricted
    # our domains to just DNA-binding, the InteracDome DF will already be subsetted for this (no additional action required).
    unique_binding_pos <- unique(binding_pos)
    
    # Add these binding positions to the dataframe
    binding_positions_df$Binding.Pos[i] <- paste(unique_binding_pos, collapse = ",")
  }
  rownames(binding_positions_df) <- binding_domains_ids_noPF
  return(binding_positions_df)
}


# Get binding positions for each domain and create a DF for them
binding_positions_df <- get_interacdome_binding_pos(interacdome_conf, binding_domains_ids_noPF, threshold_interacdome)
write.csv(binding_positions_df, paste(path, paste(paste("Saved Output Data Files/InteracDome/binding_positions_DF_", as.character(threshold_interacdome), sep = ""), ".csv", sep = ""), sep = ""))



### OPTION 2: USE PERTININT TRACKFILES ###
interacdome_domainweights <- read.csv(paste(path, "Input Data Files/InteracDome/interacdome0.3-pfam31_domainweights-GRCh38.csv", sep = ""))
interacdome_domsbyprot <- read.csv(paste(path, "Input Data Files/InteracDome/interacdome0.3-pfam31_domsbyprot-GRCh38.csv", sep = ""))

# Subset the domain weights file by the chosen binding frequency threshold
threshold_interacdome <- 0    # binding frequency threshold for individual positions
interacdome_dw_sub <- interacdome_domainweights[interacdome_domainweights$binding_frequency >= threshold_interacdome,]

# OPT: Uncomment to select only DNA-interacting confident domains/ RNA-Interacting Domains; use the nucacids name throughout code below to refer to this DF
interacdome_dw_sub_nucacids <- interacdome_dw_sub[(interacdome_dw_sub$ligand_type == "DNA_" |
                                                   interacdome_dw_sub$ligand_type == "DNABASE_" | 
                                                   interacdome_dw_sub$ligand_type == "DNABACKBONE" |
                                                   interacdome_dw_sub$ligand_type == "RNA_" |
                                                   interacdome_dw_sub$ligand_type == "RNABASE_" |
                                                   interacdome_dw_sub$ligand_type == "RNABACKBONE_" |
                                                   interacdome_dw_sub$ligand_type == "NUCACID_"),]
write.csv(interacdome_dw_sub_nucacids, paste(path, "Input Data Files/InteracDome/interacdome0.3-pfam31_domainweights-GRCh38-nucacids_only.csv", sep = ""))

# Get all the domains 
binding_domains <- unique(interacdome_dw_sub$domain_name)
binding_domains <- binding_domains[!is.na(binding_domains)]  # Remove any NA values; 2152 unique binding domains, 370 unique nucleic acid binding domains
write.csv(binding_domains, paste(path, "Saved Output Data Files/InteracDome/binding_domains.csv", sep = ""))

# Subset the domsbyprot file to only these domains; we'll be using this from here on out
interacdome_domsbyprot_sub <- interacdome_domsbyprot[interacdome_domsbyprot$Pfam_HMM_ID %fin% binding_domains,]

# OPT: In order to match these domain names to CD-Search Pfam IDs, remove the "PF" from the ID
#binding_domains_ids <- grep('PF', unlist(strsplit(binding_domains, "_", fixed = TRUE)), value = TRUE)
binding_domain_ids <- binding_domains[grep("PF", binding_domains)]
binding_domain_ids_noPF <- unlist(lapply(binding_domain_ids, function(x) unlist(strsplit(x, split = "_", fixed = TRUE))[1]))
binding_domain_ids_noPF <- unique(unlist(lapply(binding_domain_ids_noPF, FUN = function(id) str_remove(id, "PF"))))
write.csv(binding_domain_ids_noPF, paste(path, "Saved Output Data Files/InteracDome/binding_domains_noPF.csv", sep = ""))
# Read back if already done
binding_domain_ids_noPF <- read.csv(paste(path, "Saved Output Data Files/InteracDome/binding_domains_noPF.csv", sep = ""))[,2]

#' Create a data frame containing the amino acid binding positions for each protein's domain
#' @param interacdome_dw the domain weights file from InteracDome
#' @param interacdome_dbp the domains-by-protein file from InteracDome
#' @param binding_domain_ids vector of binding domain IDs that bind the particular
#' ligand(s) of interest and exceed a particular binding frequency
get_interacdome_binding_pos <- function(interacdome_dw, interacdome_dbp, binding_domain_ids) {
  # Make a new table that limits the domains to only those that bind the particular 
  # ligand(s) of interest and exceed a particular binding frequency (aka the domains 
  # in the list we're given); keep only the first two columns
  binding_positions_df <- interacdome_dbp[interacdome_dbp$Pfam_HMM_ID %fin% binding_domain_ids, c(1:2)]
  binding_positions_df <- distinct(binding_positions_df)
  
  # For all proteins and all their domains...
  binding_pos_col <- unlist(lapply(1:nrow(binding_positions_df), function(i) {
    
    print(paste(i, paste("/", nrow(binding_positions_df))))
    
    # Get the current protein and domain
    prot_id <- binding_positions_df[i, 'Protein_Sequence_ID']
    domain <- binding_positions_df[i, 'Pfam_HMM_ID']  # we already know that this domain is in our list
    
    # Get the ligand-binding matchstate positions from the domain weights file for 
    # this domain (those that exceed a pre-specified binding frequency threshold)
    matchstate_positions <- unique(interacdome_dw[interacdome_dw$domain_name == domain, 
                                                  'X1.indexed_match_state'])
    matchstate_positions <- matchstate_positions[!is.na(matchstate_positions)]
    print(matchstate_positions)
    
    # Subset the domains by protein file to include only the given protein and domain
    dbp_df_sub <- interacdome_dbp[interacdome_dbp$Protein_Sequence_ID == prot_id,]
    dbp_df_sub <- dbp_df_sub[dbp_df_sub$Pfam_HMM_ID == domain, 'matchstate.AA.0.index.AA.value']
    
    # Get the conversions for all copies of this domain within this protein (concatenate them)
    all_conv <- paste(dbp_df_sub, collapse = ",")
    conversions <- unlist(strsplit(all_conv, split = ",", fixed = TRUE))
    
    # For each matchstate position, convert it to AA position using the conversions
    converted_positions <- unlist(lapply(matchstate_positions, function(matchstate_pos) {
      conv <- conversions[which(unlist(lapply(conversions, function(x) {
        if (unlist(strsplit(x, split = ":", fixed = TRUE))[1] == as.character(matchstate_pos)) {
          return(TRUE)
        } else {return(FALSE)}
      })))]
      aa_pos <- as.numeric(unlist(strsplit(conv, split = ":", fixed = TRUE))[2]) + 1   
      return(as.character(aa_pos))
    }))
    
    # Concatenate the AA positions with commas
    converted_positions_str <- paste(converted_positions, collapse = ",")
    print(converted_positions_str)
    return(converted_positions_str)
  }))
  
  print(binding_pos_col)
  binding_positions_df$Binding.Pos <- binding_pos_col
  
  return(binding_positions_df)
}

# Get AA binding positions for each protein's domain and create a DF for them
binding_positions_df <- get_interacdome_binding_pos(interacdome_dw_sub, 
                                                    interacdome_domsbyprot_sub, 
                                                    binding_domain_ids)
write.csv(binding_positions_df, paste(path, paste(paste("Saved Output Data Files/InteracDome/binding_positions_DF_", 
                                                        as.character(threshold_interacdome), sep = ""), ".csv", sep = ""), sep = ""))

# Add Swissprot IDs
interacdome_prots <- unlist(lapply(unique(binding_positions_df$Protein_Sequence_ID), function(x) unlist(strsplit(x, ".", fixed = TRUE))[1]))
write.csv(interacdome_prots, paste(path, "Saved Output Data Files/InteracDome/interacdome_prots_ensg.csv", sep = ""))  # Write to CSV 

# Order the data frame by the protein name
binding_positions_df <- binding_positions_df[order(binding_positions_df[,'Protein_Sequence_ID']),]

#' Converts a vector of ENSP IDs to SWISSPROT IDs
#' @param prots_ensp a vector of ENSP IDs to convert
#' @param all_genes_id_conv file with multiple IDs for all proteins, from BioMart
get_swissprot_from_ensp <- function(prots_ensp, all_genes_id_conv) {
  swissprot_ids_col <- unlist(lapply(prots_ensp, function(x) {
    ensp <- unlist(strsplit(x, ".", fixed = TRUE))[1]
    paste(all_genes_id_conv[all_genes_id_conv$ensembl_peptide_id == ensp, 'uniprot_gn_id'], collapse = ";")
  }))
  return(swissprot_ids_col)
}

# Add a column of SWISSPROT IDs to this data frame
swissprot_ids_col <- get_swissprot_from_ensp(interacdome_prots, binding_positions_df, all_genes_id_conv)
binding_positions_df$Swissprot <- swissprot_ids_col
write.csv(binding_positions_df, paste(path, "Saved Output Data Files/InteracDome/binding_positions_DF_0_labeled.csv", sep = ""))


############################################################
# GET CANBIND INTERACTION INFORMATION
############################################################
# Get CanBind PertInInt tracks
# CanBind : structurally homologous regions to interaction pockets
canbind_df <- read.table(paste(path, "Input Data Files/CanBind/canbind-biolip-to-ensembl_domainweights-GRCh38.txt", sep = ""), 
                         header = TRUE, sep = "\t")
canbind_df_tracklocs <- read.table(paste(path, "Input Data Files/CanBind/canbind-biolip-to-ensembl_domsbyprot-GRCh38.txt", sep = ""), 
                         header = TRUE, sep = "\t", check.names = FALSE)
colnames(canbind_df_tracklocs)[which(colnames(canbind_df_tracklocs) == "CanBindBindingSite")] <- "EnsemblProtID_BindingSiteNum"

#' Add a column to CanBind DF that has the 1-indexed protein position corresponding to the 
#' given track 0-indexed-matchstate
#' @param df the CanBind data frame
#' @param df_tracklocs the corresponding tracklocs CanBind data frame
add_prot_pos_col <- function(df, df_tracklocs) {
  
  prot_pos_col <- unlist(lapply(1:nrow(df), function(i) {
    matchstate_pos <- df[i,'X1.Indexed.MatchState']
    prot_name <- df[i, 'EnsemblProtID_BindingSiteNum']
    conversions <- unlist(strsplit(df_tracklocs[df_tracklocs$EnsemblProtID_BindingSiteNum == prot_name, 
                                                       'matchstate:AA-0-index:AA-value'], split = ",", fixed = TRUE))
    conv <- conversions[which(unlist(lapply(conversions, function(x) {
        if (unlist(strsplit(x, split = ":", fixed = TRUE))[1] == as.character(matchstate_pos)) {
          return(TRUE)
        } else {return(FALSE)}
     })))]
    prot_pos <- unlist(strsplit(conv, split = ":", fixed = TRUE))[2]    
    return(prot_pos)
  }))
  #print(prot_pos_col)
  df$AA_position <- prot_pos_col 
  df$AA_position <- df$AA_position + 1 # add one because these are 0-indexed
  return(df)
}

canbind_df <- add_prot_pos_col(canbind_df, canbind_df_tracklocs)
write.csv(canbind_df, paste(path, "Input Data Files/CanBind/canbind-full-dataframe-with-AA-pos.csv", sep = ""))

# hist(canbind_df$FracWithin4)    # Make a histogram of the scores

# From this, use a threshold of 0 or 0.5 as a confident interaction score
canbind_thres <- 0
canbind_df_sub <- canbind_df[canbind_df$FracWithin4 >= canbind_thres,]

# Order the dataframe by the protein name
canbind_df_sub <- canbind_df_sub[order(canbind_df_sub[,'EnsemblProtID_BindingSiteNum']),]

# Get a list of proteins with highly scoring binding positions from CanBind
canbind_prots <- unique(canbind_df_sub$EnsemblProtID_BindingSiteNum)
# canbind_prots_ensg <- unlist(lapply(canbind_prots, function(x) unlist(strsplit(x, ".", fixed = TRUE))[1]))
# write.csv(canbind_prots_ensg, paste(path, paste(paste("Saved Output Data Files/CanBind/canbind_prots_ensg_", as.character(canbind_thres), sep = ""), ".csv", sep = ""), sep = ""))  # Write to CSV and upload to website

#' Convert vector of Ensembl IDs to Swissprot IDs 
#' @param prots_ensg vector of ENSG IDs to be converted
#' @param df_sub a CanBind data frame, subsetted by protein name
#' @param all_genes_id_conv file with multiple IDs for all proteins, from BioMart
get_swissprot_from_ensp <- function(prots_ensg, df_sub, all_genes_id_conv) {
  swissprot_ids_col <- unlist(lapply(prots_ensg, function(x) {
    ensp <- unlist(strsplit(x, ".", fixed = TRUE))[1]
    rep(paste(all_genes_id_conv[all_genes_id_conv$ensembl_peptide_id == ensp, 'uniprot_gn_id'], collapse = ";"),
        times = nrow(df_sub[grep(x, df_sub$EnsemblProtID_BindingSiteNum),]))
  }))
  return(swissprot_ids_col)
}
swissprot_ids_col <- get_swissprot_from_ensp(canbind_prots, canbind_df_sub, all_genes_id_conv)
canbind_df_sub$Swissprot <- swissprot_ids_col

# Write this back to an updated CanBind DF
write.csv(canbind_df_sub, paste(path, paste(paste("Saved Output Data Files/CanBind/canbind_dataframe_labeled_", as.character(canbind_thres), sep = ""), ".csv", sep = ""), sep = ""))

# OPT: Subset to only include nucleic acid binding positions
canbind_df_sub_dna <- canbind_df_sub[grep("DNA_", canbind_df_sub$LigandType),]
canbind_df_sub_dnabase <- canbind_df_sub[grep("DNABASE_", canbind_df_sub$LigandType),]
canbind_df_sub_dnabackbone <- canbind_df_sub[grep("DNABACKBONE_", canbind_df_sub$LigandType),]
canbind_df_sub_rna <- canbind_df_sub[grep("RNA_", canbind_df_sub$LigandType),]
canbind_df_sub_rnabase <- canbind_df_sub[grep("RNABASE_", canbind_df_sub$LigandType),]
canbind_df_sub_rnabackbone <- canbind_df_sub[grep("RNABACKBONE_", canbind_df_sub$LigandType),]
canbind_df_sub_nucacids <- canbind_df_sub[grep("NUCACID_", canbind_df_sub$LigandType),]

# Merge these dataframes
canbind_df_sub_tot_nucacids <- rbind(canbind_df_sub_dna, rbind(canbind_df_sub_dnabase, rbind(canbind_df_sub_dnabackbone,
                                                                                             rbind(canbind_df_sub_rna, 
                                                                                                   rbind(canbind_df_sub_rnabase,
                                                                                                         rbind(canbind_df_sub_rnabackbone, 
                                                                                                               canbind_df_sub_nucacids))))))
write.csv(canbind_df_sub_tot_nucacids, paste(path, paste(paste("Saved Output Data Files/CanBind/canbind_dataframe_labeled_nucacids_", as.character(canbind_thres), sep = ""), ".csv", sep = ""), sep = ""))

#' For each protein & binding site, get the range of the binding site 
#' (start to finish, the equivalent of a domain)
#' @param df_sub a CanBind data frame, subsetted by proteins of interest
get_binding_ranges <- function(df_sub) {
  unique_prot_and_bs <- unique(df_sub$EnsemblProtID_BindingSiteNum)
  binding_site_df <- data.frame(matrix(nrow = length(unique_prot_and_bs), ncol = 2))
  colnames(binding_site_df) <- c("Range.of.Binding.Site", "Swissprot")
  rownames(binding_site_df) <- unique_prot_and_bs
  
  for (i in 1:length(unique_prot_and_bs)) {
    # Get the track position range of the binding site
    positions <- as.numeric(df_sub[df_sub$EnsemblProtID_BindingSiteNum == unique_prot_and_bs[i], 'AA_position'])
    range <- c(min(positions), max(positions))
    range_collapsed <- paste(sort(range, decreasing = FALSE), collapse = ":")
    
    # Extract the Swissprot ID
    swissprot <- unique(df_sub[df_sub$EnsemblProtID_BindingSiteNum == unique_prot_and_bs[i], 'Swissprot'])
    swissprot <- swissprot[!is.na(swissprot)]
    
    print(paste(swissprot, range_collapsed, sep = " "))
    
    # Add to DF
    binding_site_df[i,] <- c(range_collapsed, swissprot)
  }
  return(binding_site_df)
}

binding_site_df <- get_binding_ranges(canbind_df_sub)
binding_site_df_nucacids <- get_binding_ranges(canbind_df_sub_tot_nucacids)

# Write this DF to a CSV file
write.csv(binding_site_df, paste(path, paste(paste("Saved Output Data Files/CanBind/canbind_binding_ranges_", as.character(canbind_thres), sep = ""), ".csv", sep = ""), sep = ""))
write.csv(binding_site_df_nucacids, paste(path, paste(paste("Saved Output Data Files/CanBind/canbind_binding_ranges_nucacids_", as.character(canbind_thres), sep = ""), ".csv", sep = ""), sep = ""))





############################################################
# IMPORT CONCAVITY INTERACTION INFORMATION
############################################################
# Get ConCavity PertInInt track and extract domains
# ConCavity : predicted small-molecule binding pockets
concavity_df <- read.table(paste(path, "Input Data Files/ConCavity/concavity_track.track_weights.txt", sep = ""), 
                           header = FALSE, sep = "\t")
concavity_df_tracklocs <- read.table(paste(path, "Input Data Files/ConCavity/concavity_track.track_locations.txt", sep = ""), 
                                     header = FALSE, sep = "\t")
hist(concavity_df[,4], xlim = c(0,1.0))   # Make a histogram of the scores

# OPT: from this, use a threshold of 0.1 as a confident binding site
concavity_threshold <- 0.1
concavity_df_sub <- concavity_df[concavity_df[,'Score'] >= concavity_threshold,]

# Add headers to the file
colnames(concavity_df_sub) <- c("EnsemblProtID_BindingSiteNum", "TrackLabel", "X1.Indexed.MatchState", "Score", "NumInstancesAlways10", "X")
colnames(concavity_df_tracklocs) <- colnames(canbind_df_tracklocs)

# Add AA positions as a column in original ConCavity DF
concavity_df_sub <- add_prot_pos_col(concavity_df_sub, concavity_df_tracklocs)
write.csv(concavity_df_sub, paste(path, "Input Data Files/ConCavity/concavity-dataframe-with-AA-pos-0.1.csv", sep = ""))

# Get a list of proteins with highly scoring binding positions from ConCavity
concavity_prots <- unlist(lapply(unique(concavity_df_sub[,1]), function(x) unlist(strsplit(x, ".", fixed = TRUE))[1]))
write.csv(concavity_prots, paste(path, "Saved Output Data Files/ConCavity/concavity_prots_ensg.csv", sep = ""))  # Write to CSV and upload to website

# Add Swissprot IDs as a column in original ConCavity DF (use helper function from above)
swissprot_ids_col <- get_swissprot_from_ensp(concavity_prots, concavity_df_sub, all_genes_id_conv)
concavity_df_sub$Swissprot <- swissprot_ids_col

# Write this back to an updated CanBind DF
write.csv(concavity_df_sub, paste(path, "Saved Output Data Files/ConCavity/concavity_dataframe_labeled_0.1.csv", sep = ""))

# Additionally, for each protein & binding site, get the range of the binding site 
# (start to finish, the equivalent of a domain); use function from CanBind section
binding_site_df <- get_binding_ranges(concavity_df_sub)

# Write this DF to a CSV file
write.csv(binding_site_df, paste(path, paste(paste("Saved Output Data Files/ConCavity/concavity_binding_ranges_", as.character(concavity_threshold), sep = ""), ".csv", sep = ""), sep = ""))


