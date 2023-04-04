############################################################
### Generate Protein Input Tables
### Written By: Sara Geraghty, April 2021
############################################################

# Functions for creating an regulatory protein data frame for each level of specificity (I-Protein, I-Domain,
# and I-Binding Position) in the form of matched Uniprot and ENSG ID inputs. Additional functions for 
# limiting the regulatory protein inputs by frequency/ connectedness in STRING functional network.

# STRING network and confidence threshold
BiocManager::install("STRINGdb")
BiocManager::install("EnsDb.Hsapiens.v86")
#^ This package appears to do the best at gene ID conversion, according to this source: https://shiring.github.io/genome/2016/10/23/AnnotationDbi

library("STRINGdb")
library(igraph)
library("EnsDb.Hsapiens.v86")
# keytypes(EnsDb.Hsapiens.v86)  # list all supported keytypes
library(data.table)
library(dplyr)

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/cBioPortal/METABRIC/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

prot_path <- paste(main_path, "Mutation/", sep = "")


# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


############################################################
# GENERATE PRELIMINARY TABLES
############################################################
# READ IN SWISSPROT IDs FOR EACH LEVEL OF SPECIFICITY
iprotein_prots <- as.character(unlist(fread(paste(prot_path, "swissprot_ids_missense_iprotein.csv", sep = ""), header = TRUE)[,2]))
iprotein_prots_nucacids <- as.character(unlist(fread(paste(prot_path, "swissprot_ids_missense_iprotein_nucacids.csv", sep = ""), header = TRUE)[,2]))

idomain_prots <- as.character(unlist(fread(paste(prot_path, "swissprot_ids_missense_idomain.csv", sep = ""), header = TRUE)[,2]))  # all ligands
idomain_prots_nucacids <- as.character(unlist(fread(paste(prot_path, "swissprot_ids_missense_idomain_nucacids.csv", sep = ""), header = TRUE)[,2]))  # DNA/RNA-binding only

ibindingpos_prots <- as.character(unlist(fread(paste(prot_path, "swissprot_ids_missense_ibindingpos.csv", sep = ""), header = TRUE)[,2]))  # all ligands
ibindingpos_prots_nucacids <- as.character(unlist(fread(paste(prot_path, "swissprot_ids_missense_ibindingpos_nucacids.csv", sep = ""), header = TRUE)[,2]))  # DNA/RNA-binding only

# CONVERT TO ENSG IDs
#' Converts protein uniprot IDs to ensembl IDs
#' @param protein_ids a vector of uniprot IDs to be converted
convert_to_ensembl <- function(protein_ids) {
  return(unlist(lapply(protein_ids, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'ensembl_gene_id']), collapse = ";"))))
}
iprotein_prots_ensg <- convert_to_ensembl(iprotein_prots)
iprotein_prots_nucacids_ensg <- convert_to_ensembl(iprotein_prots_nucacids)
idomain_prots_ensg <- convert_to_ensembl(idomain_prots)
idomain_prots_nucacids_ensg <- convert_to_ensembl(idomain_prots_nucacids)
ibindingpos_prots_ensg <- convert_to_ensembl(ibindingpos_prots)
ibindingpos_prots_nucacids_ensg <- convert_to_ensembl(ibindingpos_prots_nucacids)


#' Combines uniprot and ENSG IDs into a protein IDs data frame, which it writes to a CSV
#' @param protein_ids a vector of uniprot IDs for regulatory proteins
#' @param protein_ids_ensembl a vector of ENSG IDs for regulatory proteins
#' @param label the level of specificity in question ("iprotein", "idomain", etc.)
#' @param prot_path a path for where to write the output data frame
recombine_into_df_and_write <- function(protein_ids, protein_ids_ensembl, label, prot_path) {
  protein_ids_df <- data.frame("swissprot_ids" = protein_ids, "ensg_ids" = protein_ids_ensembl)
  # Remove any that do not have ENSG IDs
  protein_ids_df <- protein_ids_df[!protein_ids_df$ensg_ids == "",]
  # Write this as a CSV file for later use
  write.csv(protein_ids_df, paste(prot_path, paste("Files for Linear Model/", paste(label, "protein_ids_df.csv", sep = "_"), sep = ""), sep = ""))
}

recombine_into_df_and_write(iprotein_prots, iprotein_prots_ensg, "iprotein", prot_path)
recombine_into_df_and_write(iprotein_prots_nucacids, iprotein_prots_nucacids_ensg, "iprotein_nucacids", prot_path)
recombine_into_df_and_write(idomain_prots, idomain_prots_ensg, "idomain", prot_path)
recombine_into_df_and_write(idomain_prots_nucacids, idomain_prots_nucacids_ensg, "idomain_nucacids", prot_path)
recombine_into_df_and_write(ibindingpos_prots, ibindingpos_prots_ensg, "ibindingpos", prot_path)
recombine_into_df_and_write(ibindingpos_prots_nucacids, ibindingpos_prots_nucacids_ensg, "ibindingpos_nucacids", prot_path)


############################################################
# GENERATE FILTERED INPUT DF USING FREQUENCY THRESHOLDS + STRING
############################################################
#' Make regulatory input data frame using a mutational threshold (e.g. >=5%), with an option
#' to limit to known driver genes and to check if these drivers are connected by STRING
#' with a confidence exceeding a given threshold (in this case, keep the more highly mutated of the pair)
#' @param regprot_df a regulatory mutation data frame for a given specificity/ mutation type
#' @param mut_freq_thres a mutation frequency threshold (e.g. 5%) for which to include the 
#' given gene
#' @param abs_mut_count_thres an absolute mutation count threshold (e.g. 5 samples) for 
#' which to include the given gene
#' @param pat_ct_mapping if not NA, then we are either a) applying the thresholds above to each
#' individual cancer type and recombining, or b) writing a separate file per cancer type
#' @param spec_label a specificity label (e.g. "i-protein", "i-domain", or "i-bindingpos")
#' @param patient_ids a set of unique patient/ sample IDs for the given cohort
#' @param driver_df if not NA, use known driver DF to limit to only proteins with known drivers
#' @param all_genes_id_conv a gene ID conversion file from bioMart
#' @param string_db OPT: a STRING network with confidence scores between sets of genes
#' @param string_conf_thres OPT: a confidence threshold above which drivers are considered
#' to be closely related (and only the more highly mutated one is maintained)
#' @param data_type either 'TCGA' or 'METABRIC'
#' @param min_num_cts the minimum number of cancer types in which a given driver must 
#' exceed the mutation frequency threshold in order to be included for pan-cancer runs
create_regulatory_prot_input_df <- function(regprot_df, mut_freq_thres, abs_mut_count_thres, pat_ct_mapping,
                                            spec_label, patient_ids, driver_df, all_genes_id_conv,  
                                            string_db, string_conf_thres, data_type, min_num_cts = 2) {

  # Calculate the minimum number of mutated samples needed to keep a given gene
  mut_count_thres <- mut_freq_thres * length(patient_ids)
  if(mut_count_thres < abs_mut_count_thres) {mut_count_thres <- abs_mut_count_thres}
  print(mut_count_thres)
  
  # Extract the ENSG and Swissprot IDs
  #uniprot_ids <- extract_regprots(spec_label, regprot_df, data_type, patient_ids, mut_count_thres)
  uniprot_ids <- extract_regprots_mutCountMat(regprot_df, data_type, patient_ids, mut_count_thres, all_genes_id_conv)
  
  ensg_ids <- unlist(lapply(uniprot_ids, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id %in% unlist(strsplit(x, ";", fixed = TRUE)), 
                                   'ensembl_gene_id']), collapse = ";")))
  
  output_df <- data.frame("swissprot_ids" = uniprot_ids, "ensg_ids" = ensg_ids)
  print(output_df)
  
  # If we have multiple cancer types, split regprot DF by cancer types, generate thresholds for each, and run
  # separately. In the end, keep only those drivers that are mutated above the given threshold in at least 
  # a given number of cancer types. If we want to create a file per cancer type, do not recombine.
  new_output_df <- output_df
  
  if(length(pat_ct_mapping) > 1) {
    per_ct_thresholds <- lapply(pat_ct_mapping, function(m) {
      count_thres <- mut_freq_thres * length(m) 
      if(count_thres < abs_mut_count_thres) {count_thres <- abs_mut_count_thres}
      return(count_thres)
    })
    names(per_ct_thresholds) <- names(pat_ct_mapping)

    # For each regprot in the above output DF, check whether it exceeds the count thres in enough individual c.t.
    new_output_df_rows <- lapply(1:nrow(output_df), function(i) {
      regprot <- output_df[i, 'swissprot_ids']
      samples <- colnames(regprot_df)[2:ncol(regprot_df)]
      patients <- unlist(lapply(samples, function(s) unlist(strsplit(s, "-", fixed = TRUE))[1]))

      exceeded_per_ct <- unlist(lapply(1:length(per_ct_thresholds), function(j) {
        ct <- names(per_ct_thresholds)[j]
        ct_thres <- as.numeric(unlist(per_ct_thresholds[j]))
        ct_pats <- as.character(unlist(pat_ct_mapping[names(pat_ct_mapping) == ct]))
        num_ct_pats_w_mut <- length(intersect(ct_pats, patients))
        if(num_ct_pats_w_mut > ct_thres) {return(TRUE)}
        else {return(FALSE)}
      }))

      if(length(exceeded_per_ct[exceeded_per_ct == TRUE]) > min_num_cts) {
        return(output_df[i,])
      } else {return(NA)}
    })
    new_output_df_rows <- new_output_df_rows[!is.na(new_output_df_rows)]
    new_output_df <- do.call(rbind, new_output_df_rows)
  }
  print(head(new_output_df))
  
  # If driver_df is not NA, then we want to limit to only driver genes
  if(length(driver_df) > 1) {
    new_output_df <- new_output_df[which(new_output_df$ensg_ids %fin% driver_df$ensembl_gene_id),]
    
    # An optional filter to remove drivers that are closely related to one another
    if(!is.na(string_db)) {
      new_output_df <- filter_string_relatives(new_output_df, string_db, string_conf_thres)
    }
  }
  
  return(new_output_df)
}


#' Helper function to get and return just the regprots from a regprot DF subset
#' @param spec_label a specificity label (e.g. "i-protein", "i-domain", or "i-bindingpos")
#' @param regprot_df a regulatory mutation data frame for a given specificity/ mutation type
#' @param data_type either 'TCGA' or 'METABRIC'
#' @param patient_ids a set of unique patient/ sample IDs for the given cohort
#' @param mut_count_thres mutation count threshold
extract_regprots <- function(spec_label, regprot_df, data_type, patient_ids, mut_count_thres) {
  if(spec_label == "i-protein") {
    unique_regprots <- unique(regprot_df$Swissprot)
    print(length(unique_regprots))
    regprots <- unlist(lapply(unique_regprots, function(regprot) {
      regprot_df_sub <- regprot_df[regprot_df$Swissprot == regprot,]
      patients <- unique(unlist(strsplit(as.character(unlist(regprot_df_sub[, 'Patient'])), ";", fixed = TRUE)))
      # Limit to just patients in the intersecting set
      if(data_type == "TCGA") {
        patients_justID <- unlist(lapply(patients, function(x) 
          unlist(strsplit(x, "-", fixed = TRUE))[1]))
      } else if(data_type == "METABRIC") {
        patients_justID <- patients
      } else {
        print("Unknown data type.")
        patients_justID <- ""
      }
      patients_final <- patients[which(patients_justID %in% patient_ids)]
      if(regprot %in% c("P04637", "P42336", "P01116")) {
        print(regprot)
        print(length(patients_final))
      }
      if (length(patients_final) > mut_count_thres) {
        return(as.character(unique(regprot_df_sub[, 'Swissprot'])))
      } else {return(NA)}
    }))
  } else {
    regprots <- unlist(lapply(1:length(unique(regprot_df$Swissprot)), function(i) {
      regprot <- unique(regprot_df$Swissprot)[i]
      patients <- unlist(regprot_df[regprot_df$Swissprot == regprot, 'Patient'])
      # Limit to just patients in the intersecting set
      if(data_type == "TCGA") {
        patients_justID <- unlist(lapply(patients, function(x) 
          unlist(strsplit(x, "-", fixed = TRUE))[1]))
      } else if(data_type == "METABRIC") {
        patients_justID <- patients
      } else {
        print("Unknown data type.")
        patients_justID <- ""
      }      
      patients_final <- unique(patients[which(patients_justID %in% patient_ids)])
      if (length(patients_final) > mut_count_thres) {return(regprot)}
      else {return(NA)}
    }))
  }
  
  regprots <- unique(regprots[!is.na(regprots)])
  print(head(regprots))
  
  return(regprots)
}

#' Helper function to get and return just the regprots from a mutation count matrix
#' @param regprot_df a regulatory mutation data frame for a given specificity/ mutation type
#' @param data_type either 'TCGA' or 'METABRIC'
#' @param patient_ids a set of unique patient/ sample IDs for the given cohort
#' @param mut_count_thres mutation count threshold
extract_regprots_mutCountMat <- function(mut_count_matrix, data_type, patient_ids, 
                                         mut_count_thres, all_genes_id_conv) {
  
  # Get the number of unique potential regulatory proteins
  unique_regprots <- unique(as.character(mut_count_matrix$Gene_Symbol))
  print(length(unique_regprots))
  
  # Get the patients in the particular cancer type or subtype
  patients <- as.character(unlist(colnames(mut_count_matrix)[!(colnames(mut_count_matrix) %in% c('Gene_Symbol', 'Gene_Symbol.1', "Swissprot"))]))
  # Limit to just patients in the intersecting set
  if(data_type == "TCGA") {
    patients_justID <- unlist(lapply(patients, function(x) 
      unlist(strsplit(x, "-", fixed = TRUE))[1]))
  } else if(data_type == "METABRIC") {
    patients_justID <- patients
  } else {
    print("Unknown data type.")
    patients_justID <- ""
  }
  #patients_final <- patients[which(patients_justID %in% patient_ids)]
  
  # Subset the count matrix to include only these patients
  mut_count_matrix_sub <- mut_count_matrix[, c(1, (which(patients_justID %in% patient_ids)+1), ncol(mut_count_matrix)), with = F]

  # Check how many of these patients have a nonsynonymous mutation in each regprot
  regprots <- unlist(lapply(unique_regprots, function(regprot) {
    mat_regprot_vals <- as.integer(unlist(mut_count_matrix_sub[mut_count_matrix_sub$Gene_Symbol == regprot,]))
    
    num_pats_w_mut <- length(mat_regprot_vals[mat_regprot_vals > 0])

    if (num_pats_w_mut > mut_count_thres) {
      regprot_swissprot <- as.character(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == regprot, 
                                       'uniprot_gn_id']))[1]
        # mut_count_matrix_sub[mut_count_matrix_sub$Gene_Symbol == regprot, 'Swissprot']
      return(as.character(regprot_swissprot))
    } else {return(NA)}
  }))
  
  regprots <- unique(regprots[!is.na(regprots)])
  print(head(regprots))
  
  return(regprots)
}


#' Helper function to see if each pairwise set of driver genes are related to one another
#' in the STRING network with a confidence that exceeds the given threshold. Returns
#' the subsetted data frame
#' @param interactors a list of interaction DFs for each of the drivers in the output DF
#' @param output_df a putative output data frame with columns for swissprot and ensg ids
#' @param string_db a table with STRING relationships and associated confidence scores
#' @param string_conf_thres a confidence threshold above which we eliminate the less 
#' mutated partner in a relationship
#' @param mut_freq_df a mutation frequency DF that has the number of mutations for each
#' gene of interest
#' @param all_genes_id_conv a gene ID conversion document from BioMart
filter_string_relatives <- function(interactors, output_df, string_db, string_conf_thres, 
                                    mut_freq_df, all_genes_id_conv) {
  
  interactors <- lapply(interactors, function(interactor) 
    interactor[interactor$combined_score > string_conf_thres,])
  mut_freq_df <- mut_freq_df[mut_freq_df$Swissprot %in% output_df$swissprot_ids, ]
  
  # Convert swissprot IDs to gene names
  genes <- unlist(lapply(output_df$swissprot_ids, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'external_gene_name']), 
          collapse = ";")))
  empty_ind <- which((genes == "") | (genes == "character(0)"))
  if(length(empty_ind) > 0) {
    genes <- genes[-empty_ind]
    output_df <- output_df[-empty_ind,]
  }
  gene_swissprot_dict <- data.frame(keys = genes, values = output_df$swissprot_ids)
  
  # See if there is any overlap between the genes and the interactors of the other genes
  # If there is, check the mutation frequencies of both
  genes_to_remove <- c()
  
  for (i in 1:length(interactors)) {
    g <- names(interactors)[i]
    g_interactors <- interactors[i]
    
    for (j in (i+1):length(interactors)) {
      g_other <- names(interactors)[j]
      print(paste("Comparing", paste(g, paste("and", g_other))))
      
      if(!((length(g) == 0) | (length(g_other) == 0))) {
        # If they are indeed STRING neighbors with high confidence...
        if((g_other %in% g_interactors$to_gn) | (g_other %in% g_interactors$from_gn)) {
          print("These are high confidence neighbors!")
          
          # Check which of these is more highly mutated (we want to keep this one)
          # Checking for a "dash" is a way of summing all the sample columns and not the labeling columns
          mut_count_g <- sum(mut_freq_df[which(mut_freq_df$Gene_Symbol == g), 
                                         grepl("-", colnames(mut_freq_df))])
          mut_count_g_other <- sum(mut_freq_df[which(mut_freq_df$Gene_Symbol == g_other), 
                                               grepl("-", colnames(mut_freq_df))])
          
          if(mut_count_g > mut_count_g_other) {
            genes_to_remove <- c(genes_to_remove, g_other)
          } else {genes_to_remove <- c(genes_to_remove, g)}
          
        }
      }
    }
  }
  
  genes_to_keep <- setdiff(names(interactors), genes_to_remove)
  genes_to_keep_swissprot <- gene_swissprot_dict[gene_swissprot_dict$keys %in% genes_to_keep, 'values']
  
  output_df_sub <- output_df[output_df$swissprot_ids %in% genes_to_keep_swissprot,]
  
  return(output_df_sub)
}

#' Helper function to get the interactors for each driver gene of interest in the given 
#' table, and writes a file for each driver containing its interactors from STRING
#' @param output_df a putative output data frame with columns for swissprot and ensg ids
#' @param string_db a table with STRING relationships and associated confidence scores
#' @param all_genes_id_conv a gene ID conversion document from BioMart
get_interactors <- function(output_df, string_db, all_genes_id_conv) {
  
  # Convert swissprot IDs to gene names
  genes <- unlist(lapply(output_df$swissprot_ids, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'external_gene_name']), 
          collapse = ";")))
  
  # Get the interactors for each of the genes
  interactors <- lapply(genes, function(g) get_interactions(g, string_db))
  names(interactors) <- genes
  
  
  # Subset each of these interaction lists based on the confidence threshold
  interactors <- lapply(1:length(interactors), function(i) {
    interactor <- interactors[[i]]
    write.csv(interactor, paste0(names(interactors)[i], "_interactors.csv"))
    print(head(interactor))
    return(interactor)
  })
  names(interactors) <- genes
  
  return(interactors)
}

#' Helper function for getting the interactions for a given gene
#' @param gene the gene Hugo symbol/ name
#' @param string_db the STRING network object
get_interactions <- function(gene, string_db) {
  gene_string = string_db$mp( gene )
  
  # Get the neighbors for a given gene
  neighbors <- string_db$get_neighbors(gene_string)
  
  # Convert these from ENSP to ENSG
  neighbors_ensg <- ensembldb::select(EnsDb.Hsapiens.v86, keys=unlist(lapply(neighbors, function(x) 
    unlist(strsplit(x, ".", fixed = TRUE))[2])), 
    keytype = "PROTEINID", columns = c("SYMBOL","PROTEINID","GENEID"))
  
  # Get the interactions for a given gene
  interactions <- distinct(string_db$get_interactions(c(gene_string, neighbors)))
  interactions <- interactions[(interactions$from == gene_string) | (interactions$to == gene_string), ]
  
  interactions$from <- unlist(lapply(interactions$from, function(x)
    unlist(strsplit(x, ".", fixed = TRUE))[2]))
  interactions$to <- unlist(lapply(interactions$to, function(x)
    unlist(strsplit(x, ".", fixed = TRUE))[2]))
  
  # Add the gene name along with the ensg id
  interactions$from_gn <- unlist(lapply(interactions$from, function(x){
    if(grepl(x, gene_string)) {return(gene)}
    symb <- neighbors_ensg[neighbors_ensg$PROTEINID == x, 'SYMBOL']
    if(length(symb) == 0) {symb <- ""}
    return(symb)
  }))
  interactions$to_gn <- unlist(lapply(interactions$to, function(x) {
    if(grepl(x, gene_string)) {return(gene)}
    symb <- neighbors_ensg[neighbors_ensg$PROTEINID == x, 'SYMBOL']
    if(length(symb) == 0) {symb <- ""}
    return(symb)
  }))
  
  print(head(interactions))
  
  return(interactions)
}



############################################################
# GENERATE RELATEDNESS DATA FRAME FROM STRING, FOR PSEUDOGENE GENERATION
############################################################
#' Function to see if each pairwise set of driver genes are related to one another
#' in the STRING network with a confidence that exceeds the given threshold. Returns
#' an output data frame with all the highly confident pairs, for later creation of a 
#' joint 'pseudogene' for each related pair
#' @param interactors a list of interaction DFs for each of the drivers in the output DF
#' @param output_df a putative output data frame with columns for swissprot and ensg ids
#' @param string_conf_thres a confidence threshold above which we eliminate the less 
#' mutated partner in a relationship
#' @param mut_freq_df a mutation frequency DF that has the number of mutations for each
#' gene of interest
#' @param all_genes_id_conv a gene ID conversion document from BioMart
create_related_gene_output_df <- function(interactors, output_df, string_conf_thres, 
                                          mut_freq_df, all_genes_id_conv) {
  interactors <- lapply(interactors, function(interactor) 
    interactor[interactor$combined_score > string_conf_thres,])
  mut_freq_df <- mut_freq_df[mut_freq_df$Swissprot %in% output_df$swissprot_ids, ]
  
  # Convert swissprot IDs to gene names
  genes <- unlist(lapply(output_df$swissprot_ids, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'external_gene_name']), 
          collapse = ";")))
  #print(genes)
  empty_ind <- which((genes == "") | (genes == "character(0)") | (is.na(genes)))
  print(empty_ind)
  if(length(empty_ind) > 0) {
    genes <- genes[-empty_ind]
    output_df <- output_df[-empty_ind,]
  }
  gene_swissprot_dict <- data.frame(keys = genes, values = output_df$swissprot_ids)
  
  # See if there is any overlap between the genes and the interactors of the other genes
  pairing_table <- data.table(matrix(nrow = 0, ncol = 7))
  colnames(pairing_table) <- c("g1_name", "g2_name", "g1_swissprot", "g2_swissprot", 
                               "g1_ensg", "g2_ensg", "string_conf")
  
  for (i in 1:length(genes)) {
    g <- genes[i]
    if(g %in% names(interactors)) {
      g_interactors <- interactors[[which(names(interactors) == g)]]
      print(head(g_interactors))
      
      for (j in (i+1):length(genes)) {
        g_other <- genes[j]
        print(paste("Comparing", paste(g, paste("and", g_other))))
        
        if(!((length(g) == 0) | (length(g_other) == 0))) {
          if(!(is.na(g) | is.na(g_other))) {
            if(!(g == g_other)) {
              # If they are indeed STRING neighbors with high confidence...
              if((g_other %in% g_interactors$to_gn) | (g_other %in% g_interactors$from_gn)) {
                print("These are high confidence neighbors!")
                swissprot_g <- gene_swissprot_dict[gene_swissprot_dict$keys == g, 'values']
                swissprot_g_other <- gene_swissprot_dict[gene_swissprot_dict$keys == g_other, 'values']
                ensg_g <- output_df[output_df$swissprot_ids == swissprot_g, 'ensg_ids']
                ensg_g_other <- output_df[output_df$swissprot_ids == swissprot_g_other, 'ensg_ids']
                string_conf <- g_interactors[(g_interactors$to_gn == g_other) | (g_interactors$from_gn == g_other), 
                                             'combined_score']
                pairing_table <- rbind(pairing_table, data.table("g1_name" = g, "g2_name" = g_other,
                                                                 "g1_swissprot" = swissprot_g, "g2_swissprot" = swissprot_g_other,
                                                                 "g1_ensg" = ensg_g, "g2_ensg" = ensg_g_other, "string_conf" = string_conf))
              }
            }
          }
        }
      }
    }
  }
  
  print(head(pairing_table))
  
  return(pairing_table)
}


# IMPORT NECESSARY FILES
intersecting_patients <- read.table(paste(main_path, "Linear Model/Tumor_Only/intersecting_ids_washu.txt", sep = ""))[,1]

# Limit just to particular subtypes
patient_set_lumA <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.A_patient_ids.txt"), header = TRUE)[,1], intersecting_patients)
patient_set_lumB <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.B_patient_ids.txt"), header = TRUE)[,1], intersecting_patients)
patient_set_basal <- intersect(read.table(paste0(main_path, "Patient Subsets/Basal_patient_ids.txt"), header = TRUE)[,1], intersecting_patients)
patient_set_her2 <- intersect(read.table(paste0(main_path, "Patient Subsets/HER2_patient_ids.txt"), header = TRUE)[,1], intersecting_patients)

# Colon and rectal combined 
patient_set_colorectal <- intersect(read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids_coad_read_washu_new.txt"), header = TRUE)[,1], intersecting_patients)
patient_set_brca_blca_hnsc <- intersect(read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids_brca_blca_hnsc_washu.txt"), header = TRUE)[,1], intersecting_patients)


driver_gene_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/GRCh38_driver_gene_list.tsv", 
                           sep = "\t", header = TRUE, comment.char = "#", skip = 10)


#mut_freq_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_nonsynonymous_IntersectPatientsWashU.csv",
#                        row.names = 1, check.names = FALSE, header = TRUE)
mut_freq_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_nonsynonymous_IntersectPatientsWashU.csv",
                        header = TRUE)
#mutation_regprot_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Linear Model/Tumor_Only/Regprot_Mutation/iprotein_results_nonsynonymous_IntersectPatientsWashU.csv",
                                                        #row.names = 1, check.names = FALSE, header = TRUE)

clinical_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/clinical_data_subset_w_Nonsyn_MutCounts.csv",
                        header = T, check.names = F)

# Generate a sample/patient to cancer type mapping
get_patient_cancer_mapping <- function(clinical_df) {
  specific_types <- unique(clinical_df$project_id)
  patient_cancer_mapping <- lapply(specific_types, function(ct) {
    pats <- clinical_df[grepl(ct, clinical_df$project_id),'case_submitter_id']
    pats_ids <- unlist(lapply(pats, function(pat) 
      unlist(strsplit(pat, "-", fixed = TRUE))[3]))
    return(unique(pats_ids))
  })
  names(patient_cancer_mapping) <- specific_types
  return(patient_cancer_mapping)
}

# Create mapping, if applicable
patient_cancer_mapping <- get_patient_cancer_mapping(clinical_df)

# Call function
regprot_input_df_all_5perc <- create_regulatory_prot_input_df(mut_freq_df, 0.05, 5, patient_cancer_mapping, "i-protein",
                                                              intersecting_patients, NA, all_genes_id_conv, NA, NA, "TCGA")
regprot_input_df_driver_5perc <- create_regulatory_prot_input_df(mut_freq_df, 0.05, 5, patient_cancer_mapping, "i-protein",
                                                                 intersecting_patients, driver_gene_df, all_genes_id_conv, NA, NA, "TCGA")
# Write to file
write.csv(regprot_input_df_driver_5perc, paste0(main_path, "Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.05Freq_drivers_missense.csv"))

# Print out the names of the drivers that were written to a given file
print(paste(unique(unlist(lapply(regprot_input_df_driver_5perc$swissprot_ids, function(id) 
  all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == id, 'external_gene_name']))), collapse = ","))

# Do this per cancer type
regprot_input_df_driver_5perc_list <- lapply(1:length(patient_cancer_mapping), function(i) {
  #ct <- names(patient_cancer_mapping)[i]
  pats <- intersect(intersecting_ids, patient_cancer_mapping[[i]])
  print(length(pats))
  
  regprot_input_df <- create_regulatory_prot_input_df(mut_freq_df, 0.05, 5, NA, 
                                                      "i-protein", pats, driver_gene_df, 
                                                      all_genes_id_conv, NA, NA, "TCGA")
  return(regprot_input_df)
})
names(regprot_input_df_driver_5perc_list) <- names(patient_cancer_mapping)
for(i in 1:length(regprot_input_df_driver_5perc_list)) {
  df <- regprot_input_df_driver_5perc_list[[i]]
  ct <- unlist(strsplit(names(regprot_input_df_driver_5perc_list)[i], "-", fixed = TRUE))[2]
  write.csv(df, paste0(main_path, paste0("Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.05Freq_drivers_nonsynonymous_", paste0(ct, ".csv"))))
}

regprot_input_df_all_5perc_list <- lapply(1:length(patient_cancer_mapping), function(i) {
  #ct <- names(patient_cancer_mapping)[i]
  pats <- patient_cancer_mapping[[i]]
  
  regprot_input_df <- create_regulatory_prot_input_df(mut_freq_df, 0.05, 5, NA, 
                                                      "i-protein", pats, NA, 
                                                      all_genes_id_conv, NA, NA, "TCGA")
  return(regprot_input_df)
})
names(regprot_input_df_all_5perc_list) <- names(patient_cancer_mapping)
for(i in 1:length(regprot_input_df_all_5perc_list)) {
  df <- regprot_input_df_all_5perc_list[[i]]
  ct <- unlist(strsplit(names(regprot_input_df_all_5perc_list)[i], "-", fixed = TRUE))[2]
  write.csv(df, paste0(main_path, paste0("Mutation/Files for Linear Model/iprotein_protein_ids_df_gr0.05Freq_all_nonsynonymous_", paste0(ct, ".csv"))))
}

####################################################################



# Include STRING data
string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=0, 
                          input_directory="", protocol = "http")
string_conf_thres <- 990

regprot_input_df_all_5perc <- filter_string_relatives(regprot_input_df_all_5perc, 
                                                      string_db, string_conf_thres,
                                                      mut_freq_df, all_genes_id_conv, "TCGA")


# Get all the input files and do this for each one
iprotein_input_files_fns <- list.files(paste0(main_path, "Mutation/Files for Linear Model/"), pattern = "iprotein_protein_ids_df_gr")
iprotein_input_files <- lapply(iprotein_input_files_fns, function(x) 
  read.csv(paste0(main_path, paste0("Mutation/Files for Linear Model/", x)), header = TRUE, check.names = FALSE, row.names = 1))

iprotein_input_files_interactors <- lapply(iprotein_input_files, function(f) 
  get_interactors(f, string_db, all_genes_id_conv))

# OPTION 1: filter out closely related relatives, including only the protein in the pair
# with the higher mutation count
new_iprotein_input_files <- lapply(iprotein_input_files, function(f) {
  filtered <- filter_string_relatives(iprotein_input_files_interactors, f, string_db, string_conf_thres, 
                                      mut_freq_df, all_genes_id_conv, "TCGA")
  return(filtered)
})

lapply(1:length(new_iprotein_input_files), function(i) {
  f <- new_iprotein_input_files[[i]]
  n <- unlist(strsplit(iprotein_input_files_fns[i], ".csv", fixed = TRUE))[1]
  write.csv(f, paste0(main_path, paste0("Mutation/Files for Linear Model/", paste0(n, "_string_subset.csv"))))
})

# Read back if already done
iprotein_input_files_interactors_fns <- list.files("C:/Users/sarae/Documents/Mona Lab Work/", pattern = "_interactors.csv")
iprotein_input_files_interactors <- lapply(iprotein_input_files_interactors_fns, function(x)
  read.csv(paste0("C:/Users/sarae/Documents/Mona Lab Work/", x), header = TRUE, check.names = FALSE, row.names = 1))
names <- unlist(lapply(iprotein_input_files_interactors_fns, function(x) unlist(strsplit(x, "_", fixed = TRUE))[1]))
names(iprotein_input_files_interactors) <- names


# OPTION 2: create a file with all of these pairings that can be used to create a joint
# 'pseudogene' that merges the status of both of these genes
related_gene_output_files <- lapply(iprotein_input_files, function(f) {
  related_genes <- create_related_gene_output_df(iprotein_input_files_interactors, f, string_conf_thres, 
                                                 mut_freq_df, all_genes_id_conv)
  return(related_genes)
})

lapply(1:length(related_gene_output_files), function(i) {
  f <- related_gene_output_files[[i]]
  n <- unlist(strsplit(related_gene_output_files[i], ".csv", fixed = TRUE))[1]
  write.csv(f, paste0(main_path, paste0("Mutation/Files for Linear Model/", paste0(n, "_related_string_pairs.csv"))))
})