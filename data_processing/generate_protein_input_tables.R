############################################################
### Generate Protein Input Tables
### PUBLICATION INFORMATION
############################################################

# Functions for creating a protein data frame in the form of matched Uniprot and 
# ENSG ID inputs. Additional functions for limiting the protein inputs by 
# frequency/ connectedness in STRING functional network.

# STRING network and confidence threshold
BiocManager::install("STRINGdb")
BiocManager::install("EnsDb.Hsapiens.v86")
#^ This package appears to do the best at gene ID conversion, 
# according to this source: https://shiring.github.io/genome/2016/10/23/AnnotationDbi

library(dplyr)
library("STRINGdb")
library(data.table)
library("EnsDb.Hsapiens.v86")

# Local PATH to directory containing MAF data file
PATH <- getwd()
OUTPUT_PATH <- paste0(getwd(), "Output/Mutation/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
# ID CONVERSION FUNCTIONS
############################################################
# CONVERT TO ENSG IDs
#' Converts protein uniprot IDs to ensembl IDs
#' @param protein_ids a vector of uniprot IDs to be converted
convert_to_ensembl <- function(protein_ids) {
  return(unlist(lapply(protein_ids, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 
                                   'ensembl_gene_id']), collapse = ";"))))
}

#' Combines uniprot and ENSG IDs into a protein IDs data frame, which it writes 
#' to a CSV file in the given path
#' @param protein_ids a vector of uniprot IDs 
#' @param protein_ids_ensembl a vector of ENSG IDs
#' @param label a label to include in the output filename
#' @param prot_path a path for where to write the output data frame
recombine_into_df_and_write <- function(protein_ids, protein_ids_ensembl, label, 
                                        prot_path) {
  protein_ids_df <- data.frame("swissprot_ids" = protein_ids, 
                               "ensg_ids" = protein_ids_ensembl)
  # Remove any that do not have ENSG IDs
  protein_ids_df <- protein_ids_df[!protein_ids_df$ensg_ids == "",]
  # Write this as a CSV file
  write.csv(protein_ids_df, paste0(prot_path, 
                                   paste0("Dyscovr_Input/", 
                                          paste0(label, "protein_ids_df.csv"))))
}


############################################################
# GENERATE FILTERED INPUT DF USING FREQUENCY THRESHOLDS + STRING
############################################################
#' Make input data frame using a mutational threshold (e.g. >=5%), with an option
#' to limit to known driver genes and to check if these drivers are connected in 
#' STRING with a confidence exceeding a given threshold (in this case, keep the 
#' more highly mutated of the pair)
#' @param protein_df a mutation data frame for a given specificity/ mutation type
#' @param mut_freq_thres a mutation frequency threshold (e.g. 5%) for which to 
#' include the given gene
#' @param abs_mut_count_thres an absolute mutation count threshold (e.g. 5 
#' samples) for which to include the given gene
#' @param pat_ct_mapping if not NA, then we are either a) applying the 
#' thresholds above to each individual cancer type and recombining, or 
#' b) writing a separate file per cancer type
#' @param intersecting_ids, if pat_ct_mapping is not NA, we include intersecting 
#' patient IDs for the given cancer type
#' @param driver_df if not NA, use known driver DF to limit to only proteins 
#' with known cancer driving roles
#' @param all_genes_id_conv a gene ID conversion file from bioMart
#' @param string_db OPT: a STRING network with confidence scores between sets 
#' of genes
#' @param string_conf_thres a confidence threshold above which drivers are 
#' considered to be closely related (and only the more highly mutated one 
#' is maintained)
#' @param data_type either 'TCGA','METABRIC', or 'Chinese_TN'
#' @param top_n option to provide an "n" to limit to only the top n genes above 
#' the given frequency threshold
#' @param min_num_cts the minimum number of cancer types in which a given driver 
#' must exceed the mutation frequency threshold in order to be included for 
#' pan-cancer runs
create_regulatory_prot_input_df <- function(protein_df, mut_freq_thres, 
                                            abs_mut_count_thres, pat_ct_mapping,
                                            intersecting_ids, driver_df, 
                                            all_genes_id_conv, string_db, 
                                            string_conf_thres, data_type, top_n, 
                                            min_num_cts = 2) {
  
  # Calculate the minimum number of mutated samples needed to keep a given gene
  mut_count_thres <- mut_freq_thres * length(intersecting_ids)
  if(mut_count_thres < abs_mut_count_thres) {
    mut_count_thres <- abs_mut_count_thres
  }
  print(paste0("Mut Count Threshold: ", mut_count_thres))
  
  # Extract the ENSG and Swissprot IDs
  uniprot_ids <- extract_protein_ids(protein_df, data_type, intersecting_ids, 
                                              mut_count_thres, all_genes_id_conv,
                                              driver_df, top_n)
  
  if(length(uniprot_ids) > 0) {
    ensg_ids <- unlist(lapply(uniprot_ids, function(x) 
      paste(unique(all_genes_id_conv[
        all_genes_id_conv$uniprot_gn_id %fin% unlist(strsplit(x, ";", fixed = T)), 
                                     'ensembl_gene_id']), collapse = ";")))
    
    output_df <- data.frame("swissprot_ids" = uniprot_ids, "ensg_ids" = ensg_ids)

    # If we have multiple cancer types, split output DF by cancer types, 
    # generate thresholds for each, and run separately. In the end, keep only 
    # those proteins that are mutated above the given threshold in at least 
    # a given number of cancer types. If we want to create a file per cancer 
    # type, do not recombine.
    new_output_df <- output_df
    
    if(length(pat_ct_mapping) > 1) {
      per_ct_thresholds <- lapply(pat_ct_mapping, function(m) {
        count_thres <- mut_freq_thres * length(intersect(m, intersecting_ids)) 
        #if(count_thres < abs_mut_count_thres) {count_thres <- abs_mut_count_thres}
        return(count_thres)
      })
      names(per_ct_thresholds) <- names(pat_ct_mapping)
      
      # Get just the patient IDs (XXXX of XXXX-XX)
      samples <- colnames(protein_df)[!(colnames(protein_df) %fin% 
                                          c("Gene_Symbol", "Swissprot"))]
      patients <- unlist(lapply(samples, function(s) 
        unlist(strsplit(s, "-", fixed = T))[1]))
      
      # For each protein in the above output DF, check whether it exceeds the 
      # count threshold in sufficient individual cancer types
      new_output_df_rows <- lapply(1:nrow(output_df), function(i) {
        protein <- output_df[i, 'swissprot_ids']
        print(protein)
        
        # Check in each cancer type
        exceeded_per_ct <- unlist(lapply(1:length(per_ct_thresholds), function(j) {
          # Get the patients specific to just that cancer type and use to subset
          ct <- names(per_ct_thresholds)[j]
          ct_thres <- as.numeric(unlist(per_ct_thresholds[j]))
          ct_pats <- as.character(unlist(pat_ct_mapping[
            names(pat_ct_mapping) == ct]))
          
          # Subset the protein DF 
          protein_df_sub <- protein_df[, c(1, (which(patients %fin% ct_pats)+1), 
                                           ncol(protein_df)), 
                                       with = F]
          print(dim(protein_df_sub))

          # Get the number of mutations for the given protein in the given 
          # cancer type
          mat_protein_vals <- as.integer(unlist(protein_df_sub[
            grepl(protein, protein_df_sub$Swissprot), ]))
          num_ct_pats_w_mut <- length(mat_protein_vals[mat_protein_vals > 0])
          print(paste(ct, num_ct_pats_w_mut))
          print(num_ct_pats_w_mut / ncol(protein_df_sub))
          
          if(num_ct_pats_w_mut > ct_thres) {return(T)}
          else {return(F)}
        }))
        
        if(length(exceeded_per_ct[exceeded_per_ct == T]) > min_num_cts) {
          return(output_df[i,])
        } else {return(NA)}
      })
      new_output_df_rows <- new_output_df_rows[!is.na(new_output_df_rows)]
      new_output_df <- do.call(rbind, new_output_df_rows)
    }

    # An optional filter to remove drivers that are closely related to one 
    # another using STRING
    if(!is.na(string_db)) {
      new_output_df <- filter_string_relatives(new_output_df, string_db, 
                                               string_conf_thres)
    }
    
    return(new_output_df)
  }
  else {return(NA)}
}


#' Helper function to get and return proteins from a mutation count matrix
#' @param protein_df a mutation data frame for a given specificity/ mutation type
#' @param data_type either 'TCGA', 'METABRIC', or 'Chinese_TN'
#' @param intersecting_ids a set of unique patient/ sample IDs for the given cohort
#' @param mut_count_thres mutation count threshold
#' @param driver_gene_df if not NA, use known driver DF to limit to only proteins 
#' with known cancer driving capacity
#' @param top_n option to provide an "n" to limit to only the top n genes above 
#' the given frequency threshold
extract_protein_ids <- function(mut_count_matrix, data_type, intersecting_ids, 
                                         mut_count_thres, all_genes_id_conv, 
                                         driver_gene_df, top_n) {
  
  # First, make things faster by limiting to proteins with a value of at least 1
  # in at least 5 samples
  mut_count_matrix <- mut_count_matrix[
    which(rowSums(mut_count_matrix[,2:(ncol(mut_count_matrix)-1)] > 0) >= 5),]
  
  # Get the number of unique potential proteins
  unique_proteins <- unique(as.character(mut_count_matrix$Gene_Symbol))
  print(length(unique_proteins))
  
  # Get the patients in the particular cancer type or subtype
  patients <- as.character(unlist(
    colnames(mut_count_matrix)[!(colnames(mut_count_matrix) %fin% 
                                   c('Gene_Symbol', "Swissprot"))]))
  
  # Limit to just patients in the intersecting set
  if(data_type == "TCGA") {
    patients_justID <- unlist(lapply(patients, function(x) 
      unlist(strsplit(x, "-", fixed = T))[1]))
  } else if(data_type == "METABRIC") {
    patients_justID <- patients
  } else {
    print("Unknown data type.")
    patients_justID <- ""
  }

  # Subset the mutation count matrix to include only these patients
  mut_count_matrix_sub <- mut_count_matrix[, c(1, (which(patients_justID %fin% 
                                                     intersecting_ids)+1)), 
                                           with = F]
  print(head(mut_count_matrix_sub))

  # Check how many of these patients have a nonsynonymous mutation in each protein
  proteins_counts <- lapply(unique_proteins, function(protein) {
    mat_protein_vals <- as.integer(unlist(mut_count_matrix_sub[
      mut_count_matrix_sub$Gene_Symbol == protein, ]))
    
    num_pats_w_mut <- length(mat_protein_vals[mat_protein_vals > 0])

    if (num_pats_w_mut > mut_count_thres) {
      print(protein)
      
      if(length(driver_gene_df) > 1) {
        if(protein %fin% driver_gene_df$primary_gene_names) {
          protein_swissprot <- as.character(unlist(
            all_genes_id_conv[all_genes_id_conv$external_gene_name == protein,
                              'uniprot_gn_id']))[1]
          return(data.frame("protein" = as.character(protein_swissprot), 
                            "count" = num_pats_w_mut))
        } else  {return(data.frame("protein" = NA, "count" = NA))}
      } else {
        protein_swissprot <- as.character(unlist(
          all_genes_id_conv[all_genes_id_conv$external_gene_name == protein,
                            'uniprot_gn_id']))[1]
        return(data.frame("protein" = as.character(protein_swissprot), 
                          "count" = num_pats_w_mut))
      }
    } else {return(data.frame("protein" = NA, "count" = NA))}
  })
  
  proteins_counts_df <- do.call(rbind, proteins_counts)
  colnames(proteins_counts_df) <- c("protein", "count")
  
  proteins_counts_df <- proteins_counts_df[
    !(is.na(proteins_counts_df$protein) | 
        (length(proteins_counts_df$protein) == 0)), ]
  proteins_counts_df$count <- as.integer(proteins_counts_df$count)

  if(nrow(proteins_counts_df) > 0) {
    # Restrict to the top N, if needed
    if(!is.na(top_n)) {
      proteins_counts_df <- proteins_counts_df[
        order(proteins_counts_df$count, decreasing = T),]
      proteins_counts_df <- proteins_counts_df[1:top_n,]
    }
    
    proteins <- unique(proteins_counts_df$protein)
    print(head(proteins))
    
    return(proteins)
  } else {return(c())}
}


#' Helper function to see if each pairwise set of driver genes are related to 
#' one another in the STRING network with a confidence that exceeds the given 
#' threshold. Returns the subsetted data frame
#' @param interactors a list of interaction DFs for each of the drivers in the 
#' output DF
#' @param output_df a putative output data frame with columns for swissprot and 
#' ensg ids
#' @param string_db a table with STRING relationships and associated confidence 
#' scores
#' @param string_conf_thres a confidence threshold above which we eliminate the 
#' less mutated partner in a relationship
#' @param mut_freq_df a mutation frequency DF that has the number of mutations 
#' for each gene of interest
#' @param all_genes_id_conv a gene ID conversion document from BioMart
filter_string_relatives <- function(interactors, output_df, string_db, 
                                    string_conf_thres, mut_freq_df, 
                                    all_genes_id_conv) {
  
  interactors <- lapply(interactors, function(interactor) 
    interactor[interactor$combined_score > string_conf_thres,])
  mut_freq_df <- mut_freq_df[mut_freq_df$Swissprot %fin% output_df$swissprot_ids, ]
  
  # Convert swissprot IDs to gene names
  genes <- unlist(lapply(output_df$swissprot_ids, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 
                                   'external_gene_name']), collapse = ";")))
  empty_ind <- which((genes == "") | (genes == "character(0)"))
  if(length(empty_ind) > 0) {
    genes <- genes[-empty_ind]
    output_df <- output_df[-empty_ind,]
  }
  gene_swissprot_dict <- data.frame(keys = genes, 
                                    values = output_df$swissprot_ids)
  
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
        if((g_other %fin% g_interactors$to_gn) | 
           (g_other %fin% g_interactors$from_gn)) {
          print("These are high confidence neighbors!")
          
          # Check which of these is more highly mutated (we want to keep this one)
          # Checking for a "dash" is a way of summing all the sample columns and 
          # not the labeling columns
          mut_count_g <- sum(mut_freq_df[which(
            mut_freq_df$Gene_Symbol == g), grepl("-", colnames(mut_freq_df))])
          mut_count_g_other <- sum(mut_freq_df[which(
            mut_freq_df$Gene_Symbol == g_other), grepl("-", colnames(mut_freq_df))])
          
          if(mut_count_g > mut_count_g_other) {
            genes_to_remove <- c(genes_to_remove, g_other)
          } else {genes_to_remove <- c(genes_to_remove, g)}
          
        }
      }
    }
  }
  
  genes_to_keep <- setdiff(names(interactors), genes_to_remove)
  genes_to_keep_swissprot <- gene_swissprot_dict[gene_swissprot_dict$keys %fin% 
                                                   genes_to_keep, 'values']
  
  output_df_sub <- output_df[output_df$swissprot_ids %fin% 
                               genes_to_keep_swissprot,]
  
  return(output_df_sub)
}

#' Helper function to get the interactors for each driver gene of interest in 
#' the given table, and writes a file for each driver containing its interactors 
#' from STRING
#' @param output_df a putative output data frame with columns for swissprot and 
#' ensg ids
#' @param string_db a table with STRING relationships and associated confidence 
#' scores
#' @param all_genes_id_conv a gene ID conversion document from BioMart
get_interactors <- function(output_df, string_db, all_genes_id_conv) {
  
  # Convert Swissprot IDs to gene names
  genes <- unlist(lapply(output_df$swissprot_ids, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 
                                   'external_gene_name']), collapse = ";")))
  
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
  neighbors_ensg <- ensembldb::select(EnsDb.Hsapiens.v86, 
                                      keys=unlist(lapply(neighbors, function(x) 
    unlist(strsplit(x, ".", fixed = T))[2])), 
    keytype = "PROTEINID", columns = c("SYMBOL","PROTEINID","GENEID"))
  
  # Get the interactions for a given gene
  interactions <- distinct(string_db$get_interactions(c(gene_string, neighbors)))
  interactions <- interactions[(interactions$from == gene_string) | 
                                 (interactions$to == gene_string), ]
  
  interactions$from <- unlist(lapply(interactions$from, function(x)
    unlist(strsplit(x, ".", fixed = T))[2]))
  interactions$to <- unlist(lapply(interactions$to, function(x)
    unlist(strsplit(x, ".", fixed = T))[2]))
  
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
#' Function to see if each pairwise set of driver genes are related to one 
#' another in the STRING network with a confidence that exceeds the given 
#' threshold. Returns an output data frame with all the highly confident pairs, 
#' for later creation of a joint 'pseudogene' for each related pair
#' @param interactors a list of interaction DFs for each of the drivers in 
#' the output DF
#' @param output_df a putative output data frame with columns for swissprot and 
#' ensg ids
#' @param string_conf_thres a confidence threshold above which we eliminate the 
#' less mutated partner in a relationship
#' @param mut_freq_df a mutation frequency DF that has the number of mutations 
#' for each gene of interest
#' @param all_genes_id_conv a gene ID conversion document from BioMart
create_related_gene_output_df <- function(interactors, output_df, 
                                          string_conf_thres, mut_freq_df, 
                                          all_genes_id_conv) {
  interactors <- lapply(interactors, function(interactor) 
    interactor[interactor$combined_score > string_conf_thres,])
  mut_freq_df <- mut_freq_df[mut_freq_df$Swissprot %fin% 
                               output_df$swissprot_ids, ]
  
  # Convert swissprot IDs to gene names
  genes <- unlist(lapply(output_df$swissprot_ids, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 
                                   'external_gene_name']), collapse = ";")))
  empty_ind <- which((genes == "") | (genes == "character(0)") | 
                       (is.na(genes)))
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
    if(g %fin% names(interactors)) {
      g_interactors <- interactors[[which(names(interactors) == g)]]

      for (j in (i+1):length(genes)) {
        g_other <- genes[j]
        print(paste("Comparing", paste(g, paste("and", g_other))))
        
        if(!((length(g) == 0) | (length(g_other) == 0))) {
          if(!(is.na(g) | is.na(g_other))) {
            if(!(g == g_other)) {
              # If they are indeed STRING neighbors with high confidence...
              if((g_other %fin% g_interactors$to_gn) | 
                 (g_other %fin% g_interactors$from_gn)) {
                
                print("These are high confidence neighbors!")
                swissprot_g <- gene_swissprot_dict[gene_swissprot_dict$keys == g, 
                                                   'values']
                swissprot_g_other <- gene_swissprot_dict[
                  gene_swissprot_dict$keys == g_other, 'values']
                ensg_g <- output_df[output_df$swissprot_ids == swissprot_g, 
                                    'ensg_ids']
                ensg_g_other <- output_df[
                  output_df$swissprot_ids == swissprot_g_other, 'ensg_ids']
                string_conf <- g_interactors[(g_interactors$to_gn == g_other) | 
                                               (g_interactors$from_gn == g_other), 
                                             'combined_score']
                pairing_table <- rbind(pairing_table, 
                                       data.table("g1_name" = g, 
                                                  "g2_name" = g_other,
                                                  "g1_swissprot" = swissprot_g, 
                                                  "g2_swissprot" = swissprot_g_other,
                                                  "g1_ensg" = ensg_g, 
                                                  "g2_ensg" = ensg_g_other, 
                                                  "string_conf" = string_conf))
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


############################################################
# IMPORT NECESSARY FILES
############################################################
intersecting_patients <- read.table(
  paste0(PATH, "Dyscovr_Input/intersecting_ids.txt"))[,1]

# Optionally limit just to particular subtypes (in BRCA, for example)
patient_set_lumA <- intersect(read.table(paste0(
  PATH, "Patient_Subsets/Luminal.A_patient_ids.txt"), header = T)[,1], 
  intersecting_patients)
patient_set_lumB <- intersect(read.table(paste0(
  PATH, "Patient_Subsets/Luminal.B_patient_ids.txt"), header = T)[,1], 
  intersecting_patients)
patient_set_basal <- intersect(read.table(paste0(
  PATH, "Patient_Subsets/Basal_patient_ids.txt"), header = T)[,1], 
  intersecting_patients)
patient_set_her2 <- intersect(read.table(paste0(
  PATH, "Patien_Subsets/HER2_patient_ids.txt"), header = T)[,1], 
  intersecting_patients)

driver_gene_df <- read.csv(paste0(PATH, "GRCh38_driver_gene_list.tsv"), sep = "\t", 
                           header = T, comment.char = "#", skip = 10)
driver_gene_df_vogelstein <- driver_gene_df[
  grepl("V", driver_gene_df$cancer_driver_status),]

mut_freq_df <- fread(paste0(PATH, "Dyscovr_Input/Mutation/
                            mut_count_matrix_nonsynonymous_IntersectPatients.csv"),
                     header = T)
clinical_df <- read.csv(paste0(PATH, "clinical_data_subset.csv"),
                        header = T, check.names = F)


############################################################
# CREATE PATIENT-CANCER TYPE MAPPING
############################################################
# Generate a sample/patient to cancer type mapping
get_patient_cancer_mapping <- function(clinical_df) {
  specific_types <- unique(clinical_df$project_id)
  patient_cancer_mapping <- lapply(specific_types, function(ct) {
    pats <- clinical_df[grepl(ct, clinical_df$project_id),'case_submitter_id']
    pats_ids <- unlist(lapply(pats, function(pat) 
      unlist(strsplit(pat, "-", fixed = T))[3]))
    return(unique(pats_ids))
  })
  names(patient_cancer_mapping) <- specific_types
  return(patient_cancer_mapping)
}

# Create mapping, for pan-cancer
patient_cancer_mapping <- get_patient_cancer_mapping(clinical_df)


############################################################
# SAMPLE FUNCTION CALLS
############################################################
# Call function for one cancer type
protein_input_df_all_5perc <- create_regulatory_prot_input_df(
  mut_freq_df, 0.05, 5, patient_cancer_mapping, intersecting_patients, 
  NA, all_genes_id_conv, NA, NA, "TCGA", NA)
protein_input_df_driver_5perc <- create_regulatory_prot_input_df(
  mut_freq_df, 0.05, 5, patient_cancer_mapping, intersecting_patients, 
  driver_gene_df, all_genes_id_conv, NA, NA, "TCGA", NA)

# Write output to file
write.csv(protein_input_df_driver_5perc, paste0(
  PATH, "Dyscovr_Input/Mutation/protein_ids_df_gr0.05Freq_drivers.csv"))

# Print out the names of the drivers that were written to a given file
print(paste(unique(unlist(lapply(
  protein_input_df_driver_5perc$swissprot_ids, function(id) 
    all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == id, 
                      'external_gene_name']))), collapse = ","))

# Call function for each cancer type in a list
protein_input_df_driver_5perc_list <- lapply(
  1:length(patient_cancer_mapping), function(i) {
    pats <- intersect(intersecting_ids, patient_cancer_mapping[[i]])
    print(length(pats))
    protein_input_df <- create_regulatory_prot_input_df(mut_freq_df, 0.05, 5, NA,
                                                        pats, driver_gene_df, 
                                                        all_genes_id_conv, NA, 
                                                        NA, "TCGA", NA)
  return(protein_input_df)
})
names(protein_input_df_driver_5perc_list) <- names(patient_cancer_mapping)

# Write to files
for(i in 1:length(protein_input_df_driver_5perc_list)) {
  df <- protein_input_df_driver_5perc_list[[i]]
  ct <- unlist(strsplit(names(protein_input_df_driver_5perc_list)[i], "-", 
                        fixed = T))[2]
  write.csv(df, paste0(PATH, paste0(
  "Dyscovr_Input/Mutation/iprotein_protein_ids_df_gr0.05Freq_drivers_nonsynonymous_", 
  paste0(ct, ".csv"))))
}

# Print out the names of the drivers that were written to a given file
for(i in 1:length(protein_input_df_driver_5perc_list)) {
  df <- protein_input_df_driver_5perc_list[[i]]
  ct <- unlist(strsplit(names(protein_input_df_driver_5perc_list)[i], "-", 
                        fixed = T))[2]
  if(length(df) > 0) {
    if(nrow(df) > 0) {
      print(ct)
      df <- df[df$swissprot_ids != "",]
      genes <- sort(unique(unlist(lapply(df$swissprot_ids, function(id) 
        all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == id, 
                          'external_gene_name']))))
      if(length(genes) > 30) {
        print(length(genes))
      } else {
        print(paste(genes, collapse = ", "))
      }
    }
  }
}


# Include STRING data
string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=0, 
                          input_directory="", protocol = "http")
string_conf_thres <- 990

protein_input_df_all_5perc <- filter_string_relatives(protein_input_df_all_5perc, 
                                                      string_db, string_conf_thres,
                                                      mut_freq_df, 
                                                      all_genes_id_conv, "TCGA")

# Import all the input files and do this for each one
iprotein_input_files_fns <- list.files(paste0(
  PATH, "Dyscover_Input/Mutation/"), pattern = "protein_ids_df_gr")
iprotein_input_files <- lapply(iprotein_input_files_fns, function(x) 
  read.csv(paste0(PATH, paste0("Dyscover_Input/Mutation/", x)), header = T, 
           check.names = F, row.names = 1))

iprotein_input_files_interactors <- lapply(iprotein_input_files, function(f) 
  get_interactors(f, string_db, all_genes_id_conv))

# OPTION 1: filter out closely related relatives, including only the protein in 
# the pair with the higher mutation count
new_iprotein_input_files <- lapply(iprotein_input_files, function(f) {
  filtered <- filter_string_relatives(iprotein_input_files_interactors, f, 
                                      string_db, string_conf_thres, 
                                      mut_freq_df, all_genes_id_conv, "TCGA")
  return(filtered)
})

lapply(1:length(new_iprotein_input_files), function(i) {
  f <- new_iprotein_input_files[[i]]
  n <- unlist(strsplit(iprotein_input_files_fns[i], ".csv", fixed = T))[1]
  write.csv(f, paste0(PATH, paste0("Dyscover_Input/Mutation/", 
                                   paste0(n, "_string_subset.csv"))))
})

# OPTION 2: create a file with all of these pairings that can be used to create 
# a joint 'pseudogene' that merges the status of both of these genes
related_gene_output_files <- lapply(iprotein_input_files, function(f) {
  related_genes <- create_related_gene_output_df(iprotein_input_files_interactors, 
                                                 f, string_conf_thres, mut_freq_df, 
                                                 all_genes_id_conv)
  return(related_genes)
})

lapply(1:length(related_gene_output_files), function(i) {
  f <- related_gene_output_files[[i]]
  n <- unlist(strsplit(related_gene_output_files[i], ".csv", fixed = T))[1]
  write.csv(f, paste0(PATH, paste0("Dyscover_Input/Mutation/", 
                                   paste0(n, "_related_string_pairs.csv"))))
})