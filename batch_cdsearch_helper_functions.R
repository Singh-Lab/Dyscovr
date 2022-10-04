############################################################
# Helper functions for the running and processing of 
# Batch CD-Search
# Written By: Sara Geraghty, August 2020
############################################################

# Link to Batch CD-Search: https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi

############################################################
### PREP CD-SEARCH FILES
############################################################
#' Function takes a proteome and a label (e.g. "full", "missense", or "silent") 
#' and creates sub-FASTA files that can be uploaded directly to Batch CD-Search's 
#' web interface
#' @param proteome the full human proteome in FASTA format
#' @param label a string for labeling that indicates whether we are looking at 
#' missense mutations, silent mutations, or both ('full', 'missense', 'misAndNon', or 'silent')
#' @param dataset either "tcga", "metabric", "icgc", or "cptac3"
prep_cdsearch_files <- function(proteome, label, dataset) {
  
  #path_cdsearch <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Proteome/BRCA/"
  path_cdsearch <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Proteome/Pan-Cancer/"
  
  # Split the proteome into sub-proteomes for Batch CD-Search if needed
  if (length(proteome) > 3000) {
    proteome_subsets <- split_into_subproteomes(proteome)
    
    # Write a FASTA file for each proteome subset
    for (i in 1:length(proteome_subsets)) {
      prot_subset <- proteome_subsets[i]
      output_fasta_name <- ""
      if(dataset == "tcga") {
        output_fasta_name <- paste0(paste0(paste0(paste0("proteome_subset_", label), "_part"), i), ".fasta")
      } else {
        output_fasta_name <- paste0(paste0(paste0(paste(paste0("proteome_subset_", label), dataset, sep = "_"), "_part"), i), ".fasta")
      }
      
      write.fasta(prot_subset[[1]], names = names(prot_subset[[1]]), as.string = TRUE, file.out = paste0(path_cdsearch, output_fasta_name))
    }
  } else {
    output_filename <- ""
    if(dataset == "tcga") {
      output_filename <- paste(paste("proteome_subset", label, sep = "_"), ".fasta", sep = "")
    } else {
      output_filename <- paste(paste(paste("proteome_subset", label, sep = "_"), dataset, sep = "_"), ".fasta", sep = "")
    }
    write.fasta(proteome, names = names(proteome), as.string = TRUE, file.out = paste(path_cdsearch, output_filename))
  }
}

############################################################
### SPLIT INTO SUB-PROTEOMES FOR BATCH CD-SEARCH
############################################################
#' This function takes a proteome and, based on its length, splits it into 
#' sub-FASTA files for use with Batch CD-Search (accepts max 1000 protein entries,
#' but we use 1000 here to be on the safe side) # UPDATE: this changed from 4K to 1K in 2022
#' @param proteome a proteome of the human genome in FASTA format
split_into_subproteomes <- function(proteome) {
  proteome_subsets <- list()   # A vector of all the resultant proteome subset DFs to return
  
  num_subfiles <- as.integer(length(proteome) / 900) + 1
  
  line_count <-  1
  for (i in 1:num_subfiles) {
    if (!(length(proteome) - line_count < 900)) {
      proteome_subsets[[i]] <- proteome[line_count:(line_count+900)]
    } else {
      proteome_subsets[[i]] <- proteome[line_count:length(proteome)]
    }
    line_count <- line_count + 900
  }
  return(proteome_subsets)
}


############################################################
### GET BATCH CD-SEARCH RESULTS
############################################################
#' Function takes a label for a particular proteome and uploads all the results
#' files from Batch CD-Search for that proteome. It then recombines all the results
#' into a single dataframe, subsets it to include only specific domain hits (rather
#' than superfamilies) and returns the final dataframe
#' NOTE: Batch CD-Search files must have following naming scheme: "<LABEL>_proteome_partX_hitdata.csv" if multipart
#' @param label a string for file retrieval that indicates whether we are looking at 
#' missense mutations, silent mutations, or both ('full', 'missense', or 'silent')
#' @param dataset either "tcga", "metabric", "icgc", or "cptac3"
get_cdsearch_res <- function(label, dataset) {
  #path = "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/CD-Batch Results"
  path = "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/CD-Batch Results"
  
  all_files <- list.files(path)
  subset_files <- all_files[unlist(lapply(all_files, FUN = function(x) {startsWith(x, label)}))]
  if(dataset == "tcga") {subset_files <- subset_files[grepl(paste0(label, "_part"), subset_files)]}
  else {subset_files <- subset_files[grepl(dataset, subset_files)]}
  
  # Loop through files of interest, read into DF, and recombine into a single DF
  full_filename_1 <- paste(path, subset_files[1], sep = "/")
  #print(full_filename_1)
  cd_search_res_df <- read.csv(full_filename_1, header = TRUE, check.names = FALSE)
  if (length(subset_files) > 1) {
    for (i in 2:length(subset_files)) {
      full_filename <- paste(path, subset_files[i], sep = "/")
      res_df <- read.csv(full_filename, header = TRUE, check.names = FALSE)
      print(head(res_df))
      cd_search_res_df <- rbind(cd_search_res_df, res_df)
    }
  }
  cd_search_res_df <- distinct(cd_search_res_df)
  print(cd_search_res_df)
  return(cd_search_res_df)

  # OPT: Keep only specific hits, rather than superfamilies or non-specific 
  # cd_search_res_df_specific_only <- cd_search_res_df[cd_search_res_df$'Hit type' == "specific",]
  # print(nrow(cd_search_res_df_specific_only))
  # print(cd_search_res_df_specific_only)
  # return(cd_search_res_df_specific_only)   
}

#
#
#
#
#
#
#
#
#
#
#
#
#


############################################################
# OPTIONAL MANUAL CODE FOR UNIPROT PROTEOME:
############################################################
# Split into sub-proteome DFs
proteome_subset_part1 <- proteome_subset_full[1:2000]
proteome_subset_part2 <- proteome_subset_full[2001:5000]
proteome_subset_part3 <- proteome_subset_full[5001:8000]
proteome_subset_part4 <- proteome_subset_full[8001:11000]
proteome_subset_part5 <- proteome_subset_full[11001:14000]
proteome_subset_part6 <- proteome_subset_full[14001:length(proteome_subset_full)]

proteome_subset_missense_part1 <- proteome_subset_missense[1:3000]
proteome_subset_missense_part2 <- proteome_subset_missense[3001:6000]
proteome_subset_missense_part3 <- proteome_subset_missense[6001:9000]
proteome_subset_missense_part4 <- proteome_subset_missense[9001:length(proteome_subset_missense)]

proteome_subset_silent_part1 <- proteome_subset_silent[1:3500]
proteome_subset_silent_part2 <- proteome_subset_silent[3501:length(proteome_subset_silent)]

# Write the sub-proteome DFs to FASTA files
fasta_subset_part1 <- write.fasta(proteome_subset_part1, names = names(proteome_subset_part1), 
                                  as.string = TRUE, file.out = "proteome_subset_part1.fasta")
fasta_subset_part2 <- write.fasta(proteome_subset_part2, names = names(proteome_subset_part2), 
                                  as.string = TRUE, file.out = "proteome_subset_part2.fasta")
fasta_subset_part3 <- write.fasta(proteome_subset_part3, names = names(proteome_subset_part3), 
                                  as.string = TRUE, file.out = "proteome_subset_part3.fasta")
fasta_subset_part4 <- write.fasta(proteome_subset_part4, names = names(proteome_subset_part4), 
                                  as.string = TRUE, file.out = "proteome_subset_part4.fasta")
fasta_subset_part5 <- write.fasta(proteome_subset_part5, names = names(proteome_subset_part5), 
                                  as.string = TRUE, file.out = "proteome_subset_part5.fasta")
fasta_subset_part6 <- write.fasta(proteome_subset_part6, names = names(proteome_subset_part6), 
                                  as.string = TRUE, file.out = "proteome_subset_part6.fasta")

fasta_subset_missense_part1 <- write.fasta(proteome_subset_missense_part1, names = names(proteome_subset_missense_part1), 
                                           as.string = TRUE, file.out = "proteome_subset_missense_part1.fasta")
fasta_subset_missense_part2 <- write.fasta(proteome_subset_missense_part2, names = names(proteome_subset_missense_part2), 
                                           as.string = TRUE, file.out = "proteome_subset_missense_part2.fasta")
fasta_subset_missense_part3 <- write.fasta(proteome_subset_missense_part3, names = names(proteome_subset_missense_part3), 
                                           as.string = TRUE, file.out = "proteome_subset_missense_part3.fasta")
fasta_subset_missense_part4 <- write.fasta(proteome_subset_missense_part4, names = names(proteome_subset_missense_part4), 
                                           as.string = TRUE, file.out = "proteome_subset_missense_part4.fasta")

fasta_subset_silent_part1 <- write.fasta(proteome_subset_silent_part1, names = names(proteome_subset_silent_part1), 
                                         as.string = TRUE, file.out = "proteome_subset_silent_part1.fasta")
fasta_subset_silent_part2 <- write.fasta(proteome_subset_silent_part2, names = names(proteome_subset_silent_part2), 
                                         as.string = TRUE, file.out = "proteome_subset_silent_part2.fasta")

# Import the results and recombine into single DFs
# Full Proteome
cd_search_res_pt1 <- read.csv("full_proteome_part1_hitdata.csv", header = TRUE)
cd_search_res_pt2 <- read.csv("full_proteome_part2_hitdata.csv", header = TRUE)
cd_search_res_pt3 <- read.csv("full_proteome_part3_hitdata.csv", header = TRUE)
cd_search_res_pt4 <- read.csv("full_proteome_part4_hitdata.csv", header = TRUE)
cd_search_res_pt5 <- read.csv("full_proteome_part5_hitdata.csv", header = TRUE)
cd_search_res_pt6 <- read.csv("full_proteome_part6_hitdata.csv", header = TRUE)

# Recombine these results into one DF
cd_search_res <- rbind(cd_search_res_pt1, cd_search_res_pt2)
cd_search_res <- rbind(cd_search_res, cd_search_res_pt3)
cd_search_res <- rbind(cd_search_res, cd_search_res_pt4)
cd_search_res <- rbind(cd_search_res, cd_search_res_pt5)
cd_search_res <- rbind(cd_search_res, cd_search_res_pt6)
# Number of rows: 148456

# Missense-Subsetted Proteome
cd_search_res_missense_pt1 <- read.csv("missense_proteome_part1_hitdata.csv", header = TRUE)
cd_search_res_missense_pt2 <- read.csv("missense_proteome_part2_hitdata.csv", header = TRUE)
cd_search_res_missense_pt3 <- read.csv("missense_proteome_part3_hitdata.csv", header = TRUE)
cd_search_res_missense_pt4 <- read.csv("missense_proteome_part4_hitdata.csv", header = TRUE)

cd_search_res_missense <- rbind(cd_search_res_missense_pt1, cd_search_res_missense_pt2)
cd_search_res_missense <- rbind(cd_search_res_missense, cd_search_res_missense_pt3)
cd_search_res_missense <- rbind(cd_search_res_missense, cd_search_res_missense_pt4)
# Number of rows: 128133

# Silent-Subsetted Proteome
cd_search_res_silent_pt1 <- read.csv("silent_proteome_part1_hitdata.csv", header = TRUE)
cd_search_res_silent_pt2 <- read.csv("silent_proteome_part2_hitdata.csv", header = TRUE)
# Number of rows: 86214

cd_search_res_silent <- rbind(cd_search_res_silent_pt1, cd_search_res_silent_pt2)
  
# Subset to get specific hits only
cd_search_res_specific_only <- cd_search_res[cd_search_res$Hit.type == "specific",]
# Number of rows: 61081
cd_search_res_missense_specific_only <- cd_search_res_missense[cd_search_res_missense$Hit.type == "specific",]
# Number of rows: 52896
cd_search_res_silent_specific_only <- cd_search_res_silent[cd_search_res_silent$Hit.type == "specific",]
# Number of rows: 36167
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}