############################################################
### CNA Helper Functions
### Written By Sara Geraghty
############################################################

# Helper functions for exploratory analysis and visualization of TCGA CNA files

# Local PATH to directory containing MAF data file
PATH <- getwd()
OUTPUT_PATH <- paste0(getwd(), "Output/CNA/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "Input Data Files/all_genes_id_conv.csv"), 
                              header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
# VISUALIZE CNA EVENT FREQUENCIES AMONG MUTATION CATEGORIES
############################################################
#' Given a CNA data frame, make a barplot and boxplots that
#' visualize the frequencies of each type of CNA event across all
#' genes/ proteins
#' @param cna_df a GISTIC-style CNA data frame, with 0, 1, or -1 values
#' @param label string denoting the level of specificity ("I-Protein", 
#' "I-Domain", "I-Binding Position", or "All Genes")
visualize_cna_freq <- function(cna_df, label) {
  
  # Get the frequency of each type of CNA event (0, 1, -1)
  freq_amp_and_del <- data.frame(matrix(nrow = nrow(cna_df), ncol = 3))
  colnames(freq_amp_and_del) <- c("0", "1", "Neg1")
  rownames(freq_amp_and_del) <- rownames(cna_df)
  
  for (i in 1:nrow(cna_df)) {
    row <- as.numeric(cna_df[i,])
    print(row)
    freq_amp_and_del[i,'0'] <- length(row[row == 0])
    freq_amp_and_del[i,'1'] <- length(row[row == 1])
    freq_amp_and_del[i,'Neg1'] <- length(row[row == -1])
  }
  print(data.matrix(freq_amp_and_del))
  
  # Call helper functions to print stats & construct plots for this 
  # frequency matrix
  print_cna_statistics(cna_df, freq_amp_and_del)
  plot_cna_graphs(freq_amp_and_del, label)
  
  return(freq_amp_and_del)
}


#' Helper function for printing relevant CNA statistics
#' @param cna_df a GISTIC-style CNA data frame, with 0, 1, or -1 values
#' @param freq_amp_and_del a matrix with the frequency of each type of CNA event
print_cna_statistics <- function(cna_df, freq_amp_and_del) {
  
  # Get the average number of non-zero CNA events per patient
  
  # Total number of 1 events across all patients
  totals_1 <- sum(as.numeric(freq_amp_and_del[,2]))  
  # Total number of -1 events across all patients
  totals_neg1 <- sum(as.numeric(freq_amp_and_del[,3])) 
  # Total number of amplification events/ patient
  avg_1 <- totals_1 / ncol(cna_df)  
  # Total number of deletion events per patient
  avg_neg1 <- totals_neg1 / ncol(cna_df) 
  
  print(paste("Average number of amplifications per patient:", avg_1))
  print(paste("Average number of deletions per patient:", avg_neg1))
  
  # Get the average number of non-zero CNA events per gene
  avg_amp <- mean(freq_amp_and_del[,2])
  avg_del <- mean(freq_amp_and_del[,3])
  avg_all <- mean(c(freq_amp_and_del[,2], freq_amp_and_del[,3]))
  print(paste0("Average number of amplifications, deletions, or both per gene:", 
              c(avg_amp, avg_del, avg_all)))
}

#' Helper function for plotting CNA barplots and boxplots
#' Takes a CNA frequency matrix with 3 columns (0, 1, Neg1), where row names are 
#' genes and entries are the average number of neutral events, amplification 
#' events, and deletion events per gene
#' @param cna_freq_matrix a matrix with the frequency of each type of CNA event
#' @param label string denoting the level of specificity ("I-Protein", 
#' "I-Domain", "I-Binding Position", or "All Genes")
plot_cna_graphs <- function(cna_freq_matrix, label) {
  
  # Boxplot of the Total Number of CNA Events 
  boxplot(data.matrix(cna_freq_matrix), main = paste("Number of CNA Events", 
                                                     label, sep = ":"),
          xlab = "Type of cna Event", ylab = paste("Total Number of CNA Events", 
                                                   sep = ":"))
  
  # Individual bar plots for amplifications and deletions
  barplot(sort(cna_freq_matrix[,2], decreasing = T), 
          main = paste(label, " Number of Amplification Events", sep = ":"),
          xlab = "Gene", ylab = "Total Number of Amplification Events", 
          col = "purple")
  barplot(sort(cna_freq_matrix[,3], decreasing = T), 
          main = paste(label, " Number of Deletion Events", sep = ":"),
          xlab = "Gene", ylab = "Total Number of Deletion Events", 
          col = "green")
}


# Run this function
cna_freq_df_iprotein <- visualize_cna_freq(cna_df_iprotein, "I-Protein")
cna_freq_df_idomain <- visualize_cna_freq(cna_df_idomain, "I-Domain")
cna_freq_df_ibindingpos <- visualize_cna_freq(cna_df_ibindingpos, 
                                              "I-Binding Position")
cna_freq_df_allgenes <- visualize_cna_freq(cna_df_allgenes, "All Genes")


############################################################
# VISUALIZE CNA EVENT FREQUENCIES AMONG LIGAND-BINDING
# CATEGORIES
############################################################
#' Given a CNA dataframe and a domains DF, make a barplot and boxplots that
#' visualize the frequencies of each type of CNA event across 
#' each type of ligand-binding domain (DNA, RNA, Peptide, Ion, & SM)
#' @param cna_df a GISTIC-style CNA dataframe, with 0, 1, or -1 values
#' @param domains_df a data frame with information about domains for all 
#' regulatory proteins, produced from the mutation pipelines
#' @param ensg_uniprot_df a conversion file between ENSG and Uniprot IDs, 
#' subsetted for the regulatory proteins of interest
#' @param interacdome_df interacDome output file that predicts the ligands bound
#' by particular domains
#' @param canbind_df_sub canBind output file that predicts the ligands bound by 
#' particular DNA base positions
visualize_cna_freq_by_lig_bind_cat <- function(cna_df, domains_df, 
                                               ensg_uniprot_df, interacdome_df, 
                                               canbind_df_sub) {
  
  # Make empty frequency data frames for each ligand type
  dna_freq_df <- data.frame(matrix(nrow = length(uniprot_ids), ncol = 3))
  rna_freq_df <- data.frame(matrix(nrow = length(uniprot_ids), ncol = 3))
  peptide_freq_df <- data.frame(matrix(nrow = length(uniprot_ids), ncol = 3))
  sm_freq_df <- data.frame(matrix(nrow = length(uniprot_ids), ncol = 3))
  ion_freq_df <- data.frame(matrix(nrow = length(uniprot_ids), ncol = 3))
  
  colnames <- c("0", "1", "Neg1")
  colnames(dna_freq_df) <- colnames 
  colnames(rna_freq_df) <- colnames
  colnames(peptide_freq_df) <- colnames
  colnames(sm_freq_df) <- colnames
  colnames(ion_freq_df) <- colnames
  
  # Set the indices to 1
  dna_index <- 1
  rna_index <- 1
  peptide_index <- 1
  sm_index <- 1
  ion_index <- 1
  
  # Check the ligands that each protein binds, add it to the appropriate tables
  for (i in 1:nrow(ensg_uniprot_df)) {
    id <- ensg_uniprot_df[i, 'Ensg']
    print(id)
    print(paste(i, paste("/", length(uniprot_ids))))
    
    # Get this protein's domains
    uniprot_id <- ensg_uniprot_df[i, 'Uniprot']
    doms <- c()
    if (length(uniprot_id) > 1) {
      uniprot_ids <- unlist(lapply(uniprot_id, function(x) 
        unlist(strsplit(x, split = ";", fixed = T))))
      doms <- unique(unlist(lapply(uniprot_ids, function(x) 
        domains_df[grepl(x, domains_df$Query),'Stripped.Accessions'])))
    } else {
      doms <- unique(domains_df[grepl(uniprot_id, domains_df$Query),
                                'Stripped.Accessions'])
    }
    print(doms)
    
    # Get all the ligands these domains bind
    lig_binding_list <- get_lig_binding_list(doms, uniprot_id, interacdome_df, 
                                             canbind_df_sub)
    
    # Construct the row to add to the data frame
    row_to_add <- get_row_to_add(id, cna_df)
    
    # Add this protein to the appropriate DFs
    if (lig_binding_list[["DNA"]] > 0) {
      dna_freq_df[dna_index,] <- row_to_add
      dna_index <- dna_index + 1
    }
    if(lig_binding_list[["RNA"]] > 0) {
      rna_freq_df[rna_index,] <- row_to_add
      rna_index <- rna_index + 1
    }
    if(lig_binding_list[["SM"]] > 0) {
      sm_freq_df[sm_index,] <- row_to_add
      sm_index <- sm_index + 1
    }
    if(lig_binding_list[["ION"]] > 0) {
      ion_freq_df[ion_index,] <- row_to_add
      ion_index <- ion_index + 1
    }
    if(lig_binding_list[["PEPTIDE"]] > 0) {
      peptide_freq_df[peptide_index,] <- row_to_add
      peptide_index <- peptide_index + 1
    }
  }
  
  # Plot CNA graphs 
  plot_cna_graphs(dna_freq_df, label = "DNA")
  plot_cna_graphs(rna_freq_df, label = "RNA")
  plot_cna_graphs(peptide_freq_df, label = "Peptide")
  plot_cna_graphs(ion_freq_df, label = "Ion")
  plot_cna_graphs(sm_freq_df, label = "Small Molecule")
  
  # Store data frames in a list and return them
  output_list <- list('dna' = dna_freq_df, 'rna' = rna_freq_df, 
                      'ion' = ion_freq_df, 'sm' = sm_freq_df, 
                      'peptide' = peptide_freq_df)
  return(output_list)
}


#' Function takes a set of domains, InteracDome and CanBind DFs, and 
#' creates and returns a list of the count of each ligand they bind
#' @param doms_stripped a vector of domain accessions
#' @param id the ensembl ID of the given protein in question
#' @param interacdome_df interacDome output file that predicts the ligands bound
#' by particular domains
#' @param canbind_df_sub canBind output file that predicts the ligands bound by 
#' particular DNA base positions
get_lig_binding_list <- function(doms_stripped, id, interacdome_df, 
                                 canbind_df_sub) {
  
  lig_binding_list <- list("DNA" = 0, "RNA" = 0, "SM" = 0, "ION" = 0, 
                           "PEPTIDE" = 0, "NUCACID" = 0)
  
  # Get the ligand types this protein binds from CanBind
  if (length(id) > 1) {
    ligand_types_canbind <- unlist(lapply(id, function(x) 
      canbind_df_sub[canbind_df_sub$Swissprot == x, 'LigandType']))
  } else {
    ligand_types_canbind <- canbind_df_sub[canbind_df_sub$Swissprot == id, 
                                           'LigandType']
  }
  ligand_types_canbind_list <- unlist(strsplit(ligand_types_canbind, ",", 
                                               fixed = T))
  
  # For all domains, check the ligand type they bind in InteracDome
  ligand_types_interac <- unlist(lapply(doms_stripped, function(x) 
    interacdome_df[interacdome_df$stripped_pfam_ids == x, 'ligand_type']))
  
  # Get the union of these two
  ligand_types <- union(ligand_types_interac, ligand_types_canbind_list)
  # Taking 'unique' is OK here, since we don't actually care about the number
  # of binding domains for this ligand
  ligand_types <- unique(ligand_types[!is.na(ligand_types)])  # take unique 
  
  # Update the list with counts for these
  if (!length(ligand_types) == 0) {
    for (j in 1:length(ligand_types)) {
      lig_binding_list <- check_lig_type(ligand_types[j], lig_binding_list, 
                                         ligand_groups)
    }
  }
  return(lig_binding_list)
}

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
    broad_cat <- ligand_groups[
      ligand_groups$original_ligand_mmCIF_identifier == ligand_type, 
      'group_name']
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
  return(domain_types_list)
}

# Import ligand groups file
ligand_groups <- read.csv(paste0(PATH, "ligand_groups.csv"), check.names = F)
# Fix colnames, if needed
colnames(ligand_groups)[1] <- "group_name"

############################################################

#' Function that takes an ENSG ID for a regulatory protein, as well as a CNA data 
#' frame, and gets the 0, 1's, and -1's for that protein. Constructs them into a 
#' row to addto the main data frames and returns this row.
#' @param ensg_id the ensembl ID of the regulatory protein in question
#' @param cna_df a GISTIC-style CNA dataframe, with 0, 1, or -1 values
get_row_to_add <- function(ensg_id, cna_df) {
  row <- as.numeric(cna_df[rownames(cna_df) == ensg_id,])
  row_to_add <- c(length(row[row == 0]), length(row[row == 1]), 
                  length(row[row == -1]))
  return(row_to_add)
}

############################################################
# VISUALIZE PROTEINS WITH THE MOST CNA EVENTS OF EACH TYPE
# ACROSS PATIENT CANCER SAMPLE COHORT
############################################################
#' Given a CNA data table, make a plot of the top 100 genes that
#' have CNA events across the most patient cancer samples in the cohort. 
#' We define a CNA event as a deletion or insertion (aka, a copy number 
#' other than two, one per allele)
#' @param cna_df a raw CNA data table, with copy number values per gene
#' (rows are genes, columns are samples)
#' @param OUTPUT_PATH a PATH to where we want the file with the frequencies of 
#' various types of CNA events per gene
#' @param all_genes_id_conv a BioMart file that has gene ID conversions of 
#' various types
get_genes_with_most_cna_events <- function(cna_df, OUTPUT_PATH, 
                                           all_genes_id_conv) {
  
  # Iterate through each ENSG ID and get the number of patients that have 
  # a deletion or amplification in that gene
  new_rows <- lapply(1:nrow(cna_df), function(i) {
    row <- as.numeric(cna_df[i,])
    num_amplif <- length(row[row > 2])
    num_del <- length(row[row < 2])
    tot_num_alt <- num_amplif + num_del
    #print(c(num_amplif, num_del, tot_num_alt))
    
    return(c(num_amplif, num_del, tot_num_alt))
  })
  
  # Bind these new rows into a table 
  tab <- data.frame(do.call(rbind, new_rows))
  colnames(tab) <- c("num.amplifications", "num.deletions", "num.cna.events")
  tab$ensg.id <- cna_df$ensg.id
  
  # Sort the table
  tab <- tab[order(tab[,'num.cna.events'], decreasing = T),]

  # Add the gene name to the table
  tab$gene.name <- unlist(lapply(tab$ensg.id, function(x) 
    paste(unique(all_genes_id_conv[
      all_genes_id_conv$ensembl_gene_id == unlist(strsplit(x, ";", fixed = T))[1], 
      'external_gene_name']), collapse = ";")))
  
  # Write the table to a file
  fwrite(tab, paste0(OUTPUT_PATH, "freq_of_CNA_events_per_gene.csv"))
  
  
  # Plot the most amplified genes
  x <- 15
  barplot(tab$num.cna.events[1:x], 
          main = "Genes with Most Copy Number Alteration Events",
          xlab = "Gene", ylab = "Total Number of CNA Events", col = "lightblue", 
          names.arg = tab$gene.name[1:x], las = 2)
  tab_amp_sorted <- tab[order(tab[,'num.amplifications'], decreasing = T),]
  barplot(tab_amp_sorted$num.amplifications[1:x], 
          main = "Genes with Most Amplification Events",
          xlab = "Gene", ylab = "Total Number of Amplification Events", 
          col = "purple", names.arg = tab_amp_sorted$gene.name[1:x], las = 2)
  tab_del_sorted <- tab[order(tab[,'num.deletions'], decreasing = T),]
  barplot(tab_del_sorted$num.deletions[1:x], 
          main = "Genes with Most Deletion Events",
          xlab = "Gene", ylab = "Total Number of Deletion Events", col = "green", 
          names.arg = tab_del_sorted$gene.name[1:x], las = 2)
}


############################################################
# VISUALIZE FREQUENCY OF SAMPLES WITH X PERCENTAGE OF GENES
# THAT ARE AMPLIFIED
############################################################
#' Creates a histogram pan-cancer or per-cancer that displays the number of 
#' samples that have a copy number of >= 3 across X percentage of genes
#' @param cna_df a copy number alteration DF (samples are columns, genes 
#' are rows)
#' @param panOrPer whether we want to look pan (all cancers aggregated) or 
#' per (each cancer type individually)
#' @param clinical_df if we are looking per cancer, we need a clinical DF to link
#' the sample IDs to their cancer type
visualize_patients_perc_genes_amplif <- function(cna_df, panOrPer, clinical_df) {
  
  cna_df <- apply(cna_df[,2:ncol(cna_df)], MARGIN = c(1,2), FUN = as.integer)
  
  # Create a new T/F data frame that indicates whether the copy # is 3 or more
  cna_df_TF <- apply(cna_df[,2:ncol(cna_df)], MARGIN = c(1,2), function(x) {
    if(!is.na(x)) {
      if(x > 2) {return(T)}
      else {return(F)}
    } else {return(F)}
  })
  
  # Get the percentage of T
  cna_df_percT <- apply(cna_df_TF[,2:ncol(cna_df_TF)], MARGIN = 2, function(y) {
    y <- unlist(y)
    sumT <- length(y[y == T])
    return(sumT / length(y))
  })
  cna_df_percT_df <- as.data.frame(cna_df_percT)
  
  # Plot histogram
  if(panOrPer == "pan") {
    hist(cna_df_percT_df, xlab = "% of Genes with CNA >= 3", 
         ylab = "Number of Samples")
  } else {
    cna_df_percT_df$cancer_type <- unlist(lapply(
      rownames(cna_df_percT_df), function(pat_id) {
        pat_id <- unlist(strsplit(pat_id, "-", fixed = T))[1]
        ct <- unique(unlist(clinical_df[grepl(pat_id, 
                                              clinical_df$case_submitter_id), 
                                        'project_id'])) 
        ct <- unlist(strsplit(ct, "-", fixed = T))[2] 
        if(is.null(ct)) {return(NA)}
        return(ct)
    }))
    cna_df_percT_df <- cna_df_percT_df[!is.na(cna_df_percT_df$cancer_type),]
    ggplot(cna_df_percT_df, aes(cna_df_percT)) + geom_histogram() + 
      labs(x = "% of Genes with CNA >= 3", y = "Num. Patients") + 
      facet_wrap(~cancer_type, scales = "free")
  }
}

clinical_df <- read.csv(paste0(PATH, "clinical_wMutCounts.csv", sep = ""), 
                        header = T, check.names = F)

visualize_patients_perc_genes_amplif(pc_cna_df, "Per", clinical_df)



############################################################
# VISUALIZE PERCENTAGE OF SAMPLES WITH AMPLIFICATIONS IN GIVEN GENE
############################################################
#' Creates a series of pie charts showing the copy number percentages 
#' across a given sample population
#' that have a copy number of >= 3 across X percentage of genes
#' @param cna_df a copy number alteration DF (samples are columns, genes are rows)
visualize_patients_perc_genes_amplif <- function(cna_df) {
  
  vals <- as.numeric(unlist(cna_df[1,]))
  freq_ampl <- lapply(unique(unlist(cna_df[1,])), function(x) 
    return(length(vals[vals == x])))
  names(freq_ampl) <- as.character(unique(unlist(cna_df[1,])))
  freq_ampl_dt <- as.data.frame(freq_ampl, check.names = F)
  print(freq_ampl_dt)
  
  #pie(as.numeric(freq_ampl_dt[1,]), labels = names(freq_ampl_dt))
  print(freq_ampl_dt[,which(colnames(freq_ampl_dt) == 2)])
  
  # Sort into normal/ amplified/ deleted
  bucket_df <- data.frame(
    "Normal" = freq_ampl_dt[,which(colnames(freq_ampl_dt) == 2)],
    "Deleted" = sum(freq_ampl_dt[,which(colnames(freq_ampl_dt) %fin% c(0, 1))]),
    "Amplified" =  sum(freq_ampl_dt[,which(!(colnames(freq_ampl_dt) %fin% 
                                               c(0, 1, 2)))]))
  print(head(bucket_df))
  pie(as.numeric(bucket_df[1,]), labels = names(bucket_df))
}

# In this case, get a CNA DF limited to only patients with NO PIK3CA mutation 
# (our set for which we are evaluating TP53 mutation)
cna_df_brca_noPIK3CAmut <- cna_df_brca[,which(unlist(lapply(
  colnames(cna_df_brca), function(x) 
    unlist(strsplit(x, "-", fixed = T))[1])) %fin% brca_pik3ca_noMut_mn)]
#rownames(cna_df_brca_noPIK3CAmut) <- cna_df[,1]

# Now, we want to limit this to just the CNA status for PIK3CA among this sample 
# population; even if we remove patients with PIK3CA mutation, do they still 
# have a PIK3CA amplification?
cna_df_brca_noPIK3CAmut_pik3ca <- cna_df_brca_noPIK3CAmut[
  rownames(cna_df_brca_noPIK3CAmut) == "ENSG00000121879",]

# Visualize how many of these samples with NO PIK3CA mutation have a PIK3CA 
# amplification -- is it a high %?
visualize_patients_perc_genes_amplif(cna_df_brca_noPIK3CAmut_pik3ca)


### IMPORT NECESSARY FILES ###
# Import domain DF file
domains <- fread(paste0(PATH, "CD-Search/cdsearch_domains_all.csv"), header = T)
# Restrict to only pfam domains and add a column of just the numeric portion
domains <- domains[grepl("pfam", domains$Accession),]  # limit to just pfam IDs
domains$Stripped.Accessions <- as.character(regmatches(domains$Accession, 
                                                       regexpr("[[:digit:]]+",
                                                               domains$Accession)))  # get just the numeric portion

# Import InteracDome
interacdome_df <- fread(paste0(PATH, 
                               "Mutation/InteracDome/binding_domains_noPF.csv"), 
                        header = T)
interac_ids <- interacdome_df[,1]

# Import CanBind
canbind_df_sub <- fread(paste0(PATH, 
                               "Mutation/CanBind/canbind_dataframe_labeled.csv"), 
                        header = T)

# Get the CNA protein IDs and convert to Uniprot IDs
ensg_ids <- rownames(cna_df_allgenes)

ensg_uniprot_df <- data.frame(matrix(nrow = length(ensg_ids), ncol = 2))
colnames(ensg_uniprot_df) <- c('Ensg', "Uniprot")
for (i in 1:length(ensg_ids)) {
  print(i)
  ensg <- ensg_ids[i]
  uniprots <- unique(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == ensg, 
                                       'uniprot_gn_id'])
  uniprots <- paste(uniprots, collapse = ";")
  ensg_uniprot_df[i,] <- c(ensg, uniprots)
}

# Create a file that matches the two for this function (faster than full 
# conversion file)
write.csv(ensg_uniprot_df, paste(OUTPUT_PATH, "cna_uniprot_ensg_conv_df.csv"))

# Or read back, if already done
ensg_uniprot_df <- fread(paste(OUTPUT_PATH, "cna_uniprot_ensg_conv_df.csv", sep = ""), header = T)

# Alternatively, just get the unique Uniprot IDs and save these
uniprot_ids <- unique(uniprot_ids[!is.na(uniprot_ids)])
uniprot_ids <- uniprot_ids[uniprot_ids != ""]
write.csv(uniprot_ids, paste0(OUTPUT_PATH, "cna_uniprot_all_genes.csv"))

# Read CNA files
cna_df_iprotein <- fread(paste0(OUTPUT_PATH, "CNA_DF_iprotein.csv"), header = T)
cna_df_idomain <- fread(paste0(OUTPUT_PATH, "CNA_DF_idomain.csv"), header = T)
cna_df_ibindingpos <- fread(paste0(OUTPUT_PATH, "CNA_DF_ibindingpos.csv"), 
                            header = T)
cna_df_allgenes <- fread(paste0(OUTPUT_PATH, "CNA_DF_AllGenes.csv"), header = T)
cna_df_allgenes <- fread(paste0(OUTPUT_PATH, "CNV_DF_AllGenes_CancerOnly.csv"), 
                         header = T)
cna_df_allgenes <- na.omit(data.frame(cna_df_allgenes))
colnames(cna_df_allgenes)[1] <- 'ensg.id'

### RUN FUNCTION ###
cna_freq_dfs_all_ligands_all_genes <- visualize_cna_freq_by_lig_bind_cat(cna_df_allgenes, 
                                                                         domains, 
                                                                         ensg_uniprot_df, 
                                                                         interacdome_df, 
                                                                         canbind_df_sub)
cna_freq_df_dna <- cna_freq_dfs_all_ligands_all_genes[1]
cna_freq_df_rna <- cna_freq_dfs_all_ligands_all_genes[2]
cna_freq_df_ion <- cna_freq_dfs_all_ligands_all_genes[3]
cna_freq_df_sm <- cna_freq_dfs_all_ligands_all_genes[4]
cna_freq_df_peptide <- cna_freq_dfs_all_ligands_all_genes[5]

# Save results to files
write.csv(cna_freq_df_dna, paste0(OUTPUT_PATH, 
                                  "Freq.Tables/dna_freq_df_all_genes.csv"))
write.csv(cna_freq_df_rna, paste0(OUTPUT_PATH, 
                                  "Freq.Tables/rna_freq_df_all_genes.csv"))
write.csv(cna_freq_df_ion, paste0(OUTPUT_PATH, 
                                  "Freq.Tables/ion_freq_df_all_genes.csv"))
write.csv(cna_freq_df_peptide, paste0(OUTPUT_PATH, 
                                      "Freq.Tables/peptide_freq_df_all_genes.csv"))
write.csv(cna_freq_df_sm, paste0(OUTPUT_PATH, 
                                 "Freq.Tables/sm_freq_df_all_genes.csv"))

# Get genes with most CNA events
get_genes_with_most_cna_events(cna_df_allgenes, OUTPUT_PATH) 

# Limit the CNA DF to just transcription factors
# Import list of human transcription factors, Lambert et al (PMID:29425488)
tf_list <- fread(paste0(PATH, "human_TFs_ensg.txt"), header = F)
# Remove the first two rows
tf_list <- tf_list[3:nrow(tf_list),]
# Limit the CNA DF to just these
cna_df_tfs <- cna_df_allgenes[cna_df_allgenes$ensg.id %fin% 
                                as.character(unlist(tf_list[,1])),] 
fwrite(cna_df_tfs, paste(OUTPUT_PATH, "CNA_DF_TFsOnly.csv", sep = ""))
# Run function
get_genes_with_most_cna_events(cna_df_tfs, OUTPUT_PATH) 
