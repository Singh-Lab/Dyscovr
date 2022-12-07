############################################################
# Helper functions for exploratory analysis and visualization
# of MAF files
# Written By: Sara Camilli, July 2020
############################################################
library(maftools)
library(stats)
library(ggplot2)
library(dplyr)
library(tidytext)
library(TCGAbiolinks)
library("RColorBrewer")
library("gplots")
library(data.table)


# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


############################################################

#' SUMMARIZE MAF
#' This function takes a MAF file that has been read by 
#' maftools and analyzes it as a whole using maftools functions;
#' returns gene summary
#' @param maf a MAF file that has been read in by maftools
summarize_maf <- function(maf) {
  # Get summaries of the maf file
  gene_summary <- getGeneSummary(maf)  
  plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  return(gene_summary)

}

############################################################

#' GET PFAM DOMAINS
#' This function takes a MAF file and that has been read by 
#' maftools and an "n" value. Extracts protein and domain 
#' summaries for the top n genes, which it then returns.
#' @param maf a MAF file that has been read in by maftools
#' @param n number of top genes we'd like pfam domains for
get_pfam_domains <- function(maf, n) {
  # Examine the top pfam domains and get summaries using maftools
  maf.pfam = pfamDomains(maf = maf, top = n)   
  maf_protein_summary <- maf.pfam$proteinSummary
  maf_domain_summary <- maf.pfam$domainSummary
  return(maf.pfam)
}

############################################################

#' GET MUTATION COUNT MATRIX
#' This function takes a MAF file that has been read by 
#' maftools and returns a mutation-count matrix
#' THIS ONLY INCLUDES NONSYNONYMOUS VARIANTS unless includeSyn 
#' is specified as FALSE. This only includes primary solid 
#' tumor samples (not normal).
#' @param maf a MAF file that has been read in by maftools
get_mut_count_matrix <- function(maf) {
  #mut_count_matrix <- mutCountMatrix(maf = maf, countOnly = NULL, removeNonMutated = FALSE)
  #mut_count_matrix <- mutCountMatrix(maf = maf, countOnly = "Missense_Mutation")
  #mut_count_matrix <- mutCountMatrix(maf = maf, countOnly = c("Missense_Mutation", 
                                                              #"Nonsense_Mutation"), , removeNonMutated = FALSE)
  mut_count_matrix <- mutCountMatrix(maf = maf, countOnly = c("Missense_Mutation", 
                                                              "Nonsense_Mutation",
                                                              "Nonstop_Mutation",
                                                              "Splice_Site"), removeNonMutated = FALSE)
  #mut_count_matrix <- mutCountMatrix(maf = maf, countOnly = "Silent", 
  #includeSyn = TRUE, removeNonMutated = FALSE)
  #print(mut_count_matrix)
  return(mut_count_matrix)
}

# As an aside--add rows of 0's to this matrix for genes that are not included
allgene_targets <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/allgene_targets.csv", 
                            header = TRUE, row.names = 1, check.names = FALSE)

intersecting_patients_brca <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Linear Model/Tumor_Only/intersecting_ids.txt", 
                                         header = TRUE, row.names = 1, check.names = FALSE)[,1]
intersecting_patients_pc <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Linear Model/Tumor_Only/intersecting_ids.txt", 
                                         header = FALSE, check.names = FALSE)[,1]

#' Function to add genes with no mutations across any sample to the mutation count matrix
#' @param mut_count_matrix a maftools-produced mutation count matrix
#' @param allgene_targets a list of all target genes created from ensembl
add_missing_genes_to_mut_count_mat <- function(mut_count_matrix, allgene_targets) {
  # Get ENSG IDs from gene symbols
  mut_count_mat_ensgs <- unlist(lapply(rownames(mut_count_matrix), function(x) 
    unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'ensembl_gene_id']))[1]))
  
  # Get the genes that are missing
  missing_genes <- setdiff(allgene_targets$ensg, mut_count_mat_ensgs)
  
  # Create a table for them and append to mutation count matrix
  missing_genes_df <- data.table(matrix(nrow = length(missing_genes), ncol = ncol(mut_count_matrix)))
  missing_genes_df[is.na(missing_genes_df)] <- 0
  rownames(missing_genes_df) <- make.names(unlist(lapply(missing_genes, function(x) 
    unique(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == x, 'external_gene_name']))[1])), unique = TRUE)
  colnames(missing_genes_df) <- colnames(mut_count_matrix)
  
  mut_count_matrix_full <- rbind(mut_count_matrix, missing_genes_df)
  
  return(mut_count_matrix_full)
}

#' Function to add patients with no mutation of the the given specificity in any protein
#' @param mut_count_matrix a maftools-produced mutation count matrix
#' @param intersecting_patients a list of patients that have all data types within the given
#' cancer type, before filtering is done
add_missing_patients_to_mut_count_mat <- function(mut_count_matrix, intersecting_patients) {
  colnames(mut_count_matrix) <- unlist(lapply(colnames(mut_count_matrix), function(x)
    paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))
  missing_patients <- setdiff(intersecting_patients, unlist(lapply(colnames(mut_count_matrix), function(x) 
    unlist(strsplit(x, "-", fixed = TRUE))[1])))
  new_pat_df <- data.frame(matrix(nrow = nrow(mut_count_matrix), ncol = length(missing_patients)))
  new_pat_df[is.na(new_pat_df)] <- 0
  colnames(new_pat_df) <- paste0(missing_patients, "-01A")
  
  new_mut_count_mat <- cbind(mut_count_matrix, new_pat_df)
  return(new_mut_count_mat)
}

mut_count_matrix_full <- add_missing_genes_to_mut_count_mat(mut_count_matrix, allgene_targets)
mut_count_matrix_full <- add_missing_patients_to_mut_count_mat(mut_count_matrix_full, intersecting_patients)

write.csv(mut_count_matrix_full, paste0(path, "Saved Output Data Files/BRCA/Mutation/Mutation Count Matrices/mut_count_matrix_nonsyn_nonMutantGenesAdded.csv"))
write.csv(mut_count_matrix_full, paste0(path, "Saved Output Data Files/Pan-Cancer/Mutation/Mutation Count Matrices/mut_count_matrix_nonsynonymous_ALL_inclNonmut.csv"))


############################################################

#' GET AVERAGE # MISSENSE MUTATIONS FOR GIVEN GENE
#' Input a MAF file and corresponding mutation-count matrix.
#' Get the average number of missense mutations for 
#' top genes, E.g. "TP53", "PIK3CA" and produce a Lollipop plot.
#' @param mut_count_matrix a mutation-count matrix produced by maftools
#' @param maf_file the maf file read in by maftools
#' @param top_gene the gene name of the gene we're interested in
#' looking at average number of mutations for
get_avg_missense_muts_per_gene <- function(mut_count_matrix, maf_file, top_gene) {
  index <- which(rownames(mut_count_matrix) == top_gene)
  gene_avg_muts <- mean(as.numeric(mut_count_matrix_missense[index,]))
  lollipopPlot(maf = maf_file, gene = top_gene, showMutationRate = TRUE)
  return(gene_avg_muts)
}

############################################################

#' GET TOTAL NUMBER OF MUTATIONS PER PATIENT
#' Given a mutation count matrix, produces a DF of the total number of 
#' mutations per participant, for each of their samples if applicable.
#' Plots a histogram of the total mutation counts across patients
#' for first sample as well as a boxplot.
#' @param mut_count_matrix a mutation-count matrix produced by maftools
get_num_mutations_per_patient <- function(mut_count_matrix) {
  total_mut_count_df <- data.frame(matrix(ncol = 1, nrow = 0))
  
  for (i in 1:ncol(mut_count_matrix)) {
    id <- colnames(mut_count_matrix)[i]
    new_row <- c(sum(as.numeric(mut_count_matrix[,i])))
    total_mut_count_df <- rbind(total_mut_count_df, new_row)
    rownames(total_mut_count_df)[i] <- id 
  }
  colnames(total_mut_count_df) <- c("Total.Num.Mut")
  
  # Create visualizations
  # hist(total_mut_count_df$Sample.1, xlab = "Total Number of Mutations", ylab = "Number of Patients",
       #main = "Histogram of Total Mutation Count Per Patient",
       #xlim = c(0, max(total_mut_count_df$Sample.1)+20), breaks = 50)
  #boxplot(total_mut_count_df$Sample.1, ylab = "Total Number of Mutations")
  # Appear to be two high-mutational outliers
  #identify(rep(1, length(total_mut_count_df$Sample.1)), total_mut_count_df$Sample.1, 
           #labels = rownames(total_mut_count_df))
  
  return(total_mut_count_df)
}


############################################################

#' GET PATHWAY AND VAF PLOTS
#' Given a maf file and n value, creates following plots:
#' 1. Variant allele frequency plot
#' 2. Oncoplot (top n genes)
#' 3. Oncogenic pathways plot
#' @param maf_file a maf file processed by maftools
#' @param n number of top genes we are interested in for
#' our oncoplot
create_vaf_oncoplots <- function(maf_file, n) {
  plotVaf(maf = maf_file)
  oncoplot(maf = maf_file, top = n)
  OncogenicPathways(maf = maf_file)
}


############################################################

#' VISUALIZE MUTATION DISTRIBUTION
#' Given a maf data frame and corresponding clinical data frame:
#' 1. Reads in the maf using maftools
#' 2. Uses helper function to get a patient ID:cancer type mapping 
#' for all patients in the maf file
#' 3. Gets the mutation count matrix for each cancer type
#' 4. Plots the mutation counts as a histogram per cancer type
#' @param total_mutation_count_matrix a mutation-count matrix produced
#' by maftools
#' @param patient_id_cancer_type_map a clinical data frame produced
#' by the 'make_patient_id_cancer_type_map()' function
visualize_mutation_distrib <- function(total_mutation_count_matrix, patient_id_cancer_type_map) {

  # For every cancer type, get a separate visualization of the mutation
  cancer_types <- unique(patient_id_cancer_type_map[,'Project.ID'])
  print(cancer_types)
  mut_count_matrices <- lapply(cancer_types, function(ct) {
    # Subset the mapping to just this cancer type & get all the patients of this type
    patients <- patient_id_cancer_type_map[patient_id_cancer_type_map$Project.ID == ct, 'Patient.ID'] 
    print(patients)
    # Use these patients to subset the total count matrix file
    #count_matrix_sub <- total_mutation_count_matrix[grepl(paste(patients, collapse = "|"), total_mutation_count_matrix$X),]            #maf[unlist(strsplit(maf$Tumor.Sample.Barcode, "-", fixed = TRUE))[3] %fin% patients,]  # May need to fix this, or use subset
    count_matrix_sub <- total_mutation_count_matrix[grepl(paste(patients, collapse = "|"), rownames(total_mutation_count_matrix)), , drop = FALSE]            #maf[unlist(strsplit(maf$Tumor.Sample.Barcode, "-", fixed = TRUE))[3] %fin% patients,]  # May need to fix this, or use subset
    # Add a column for the cancer type
    count_matrix_sub$Project.ID <- rep(ct, nrow(count_matrix_sub)) 
    return(as.data.frame(count_matrix_sub))
  }) 
  
  # Re-combine these all into one big table for plotting
  total_mut_count_tab <- do.call("rbind", mut_count_matrices)
  print(head(total_mut_count_tab))
  
  return(total_mut_count_tab)
}


#' Helper function that, using a list of patient IDs and a clinical DF from
#' the TCGA, creates a mapping between patient IDs and their cancer type
#' (since cancer type is not listed in aggregated MAF files)
#' @param patient_ids a vector of four-digit patient IDs from the TCGA
#' @param clinical_df a clinical DF produced by the TCGA for the given cohort
make_patient_id_cancer_type_map <- function(patient_ids, clinical_df) {
  
  rows <- lapply(patient_ids, function(id) {
    # Extract the "project ID" for each patient ID and split out the "TCGA-" part
    project_id <- unlist(strsplit(unique(clinical_df[grepl(id, clinical_df$tcga_barcode), 'project_id']), "-", fixed = TRUE))[2]
    if(!length(project_id) == 0) {
      return(c(id, project_id))
    }
  })
  rows <- rows[!is.na(rows)]
  df <- as.data.frame(do.call("rbind", rows))
  colnames(df) <- c("Patient.ID", "Project.ID")
  print(df)
  return(df)
}

# Import needed information
total_mutation_count_matrix <- read.csv(paste(path, "Saved Output Data Files/Pan-Cancer/Mutation/Mutation Count Matrices/total_mut_count_matrix_per_patient.csv", sep = ""))
patient_ids <- unique(unlist(lapply(total_mutation_count_matrix$X, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3])))

# Make the patient ID-cancer type mapping
cancer_patient_id_type_map <- make_patient_id_cancer_type_map(patient_ids, clinical_df)

# Visualize the mutation distribution
total_mut_count_tab <- visualize_mutation_distrib(total_mutation_count_matrix, cancer_patient_id_type_map)

# Plot the total mutation counts as a histogram for each table
ggplot(total_mut_count_tab, aes(x=Total.Num.Mut)) + geom_histogram() +
  labs(x = "Total Mutation Count", y = "Frequency") + facet_wrap(~Project.ID)


### OPTIONAL, IF NEEDED ###
#' Function to get per-cancer sub-MAF files from a pan-cancer aggregated MAF file
#' @param maf a pan-cancer aggregated MAF file to split
#' @param clinical_df an aggregated clinical DF across all cancer types in the TCGA
get_per_cancer_sub_mafs <- function(maf, clinical_df) {
  
  # Get a mutation count matrix for every cancer type; split the maf file into sub-mafs
  # based on patient cancer types
  patient_ids <- unique(unlist(lapply(maf$Tumor_Sample_Barcode, 
                                      function(x) unlist(strsplit(x, "-", fixed = TRUE))[3])))
  patient_id_cancer_type_map <- make_patient_id_cancer_type_map(patient_ids, clinical_df)
  
  # For every cancer type, get a separate visualization of the mutation
  cancer_types <- unique(patient_id_cancer_type_map[,'Project.ID'])
  print(cancer_types)
  lapply(cancer_types, function(ct) {
    # Subset the mapping to just this cancer type & get all the patients of this type
    patients <- patient_id_cancer_type_map[patient_id_cancer_type_map$Project.ID == ct, 'Patient.ID'] 
    # Use these patients to subset the MAF file
    maf_sub <- maf[grepl(patients, maf$Tumor.Sample.Barcode),]            #maf[unlist(strsplit(maf$Tumor.Sample.Barcode, "-", fixed = TRUE))[3] %fin% patients,]  # May need to fix this, or use subset
    
    # Write this MAF to a file
    write.table(maf_sub, paste(path, paste("Saved Output Data Files/Pan-Cancer/Mutation/Sub-MAF Files/maf_file_", paste(ct, ".maf", sep = ""), sep = ""), sep = ""))
  }) 
}

# Call this function
get_per_cancer_sub_mafs(maf_file_df_missense, clinical_df)


############################################################

#' FILTER HYPERMUTATORS
#' Given a maf filename, maf DF, & UUID-Barcode conversion:
#' 1. Reads in the maf using maftools
#' 2. Uses helper function to get mut-count matrix
#' 3. Determines outliers and removes them
#' 4. Returns filtered maf DF and the filtered mutation-count matrix
#' @param maf_filename the name of the maf file, including its path
#' @param maf_df a MAF file processed by maftools
filter_hypermutators <- function(maf_filename, maf_df) {
  
  maf <- read.maf(maf_filename)
  mut_count_matrix <- get_mut_count_matrix(maf)
  total_mut_count_tab <- get_num_mutations_per_patient(mut_count_matrix)
  total_mut_counts <- total_mut_count_tab[,1]
  
  # Get the IQR
  iqr <- IQR(total_mut_counts)
  # Get Q3 of data
  q3 <- as.numeric(quantile(total_mut_counts)[4])
  # Get Q3 + 1.5(IQR)
  thres <- as.numeric(q3 + (1.5 * iqr))

  # Filter the mutation count table by this threshold
  total_mut_count_tab_filt <- subset(total_mut_count_tab, Total.Num.Mut < thres)

  # Get the remaining patients left in the matrix
  remaining_samp <- rownames(total_mut_count_tab_filt)

  # Filter the maf DF to include only these patients & return
  maf_df_filt <- maf_df[maf_df$Tumor_Sample_Barcode %fin% remaining_samp,]

  return(list(maf_df = maf_df_filt, mut_count_df = total_mut_count_tab_filt))
}

############################################################

#' UPDATE CLINICAL FILE WITH TOTAL MUTATION COUNTS
#' Given a total mutation count table (filtered to exclude 
#' hypermutators), and a clinical DF:
#' 1. Extracts all unique case_ids 
#' 2. Removes rows from clinical data frame without matching case ids
#' 3. Adds a column in the clinical DF with the total mutation
#' count for each patient
#' 4. Returns the updated clinical DF
#' @param total_mut_count_tab a mutation-count matrix from maftools that has
#' been filtered to remove hypermutators
#' @param clinical_df a clinical DF from the TCGA or METABRIC
#' @param is_tcga a TRUE/FALSE value indicating if these are TCGA IDs
update_clin_with_mut_counts <- function(total_mut_count_tab, clinical_df, is_tcga) {
  
  # Get the unique patient IDs for our mutation count table
  ids <- rownames(total_mut_count_tab)
  total_mut_counts <- c()
  
  if(is_tcga) {
    ids <- unlist(lapply(ids, function(id) 
      paste(unlist(strsplit(id, split = "-", fixed = TRUE))[1:3], collapse = "-")))
    # Subset the clinical data frame to only these IDs
    clinical_df <- clinical_df[clinical_df$case_submitter_id %fin% ids,]
    total_mut_counts <- unlist(lapply(clinical_df$case_submitter_id, function(case) 
      total_mut_count_tab[grepl(case, rownames(total_mut_count_tab)), 'Total.Num.Mut']))
  } 
  else {
    clinical_df <- clinical_df[clinical_df$SAMPLE_ID %fin% ids,]
    total_mut_counts <- unlist(lapply(clinical_df$SAMPLE_ID, function(case) 
      total_mut_count_tab[grepl(case, rownames(total_mut_count_tab)), 'Total.Num.Mut']))
    
  }

  # Add vector as a new column 
  clinical_df$Total.Num.Muts <- total_mut_counts
  
  return(clinical_df)
}


############################################################

#' FOR A PARTICULAR DRIVER GOI, KEEP ONLY MUTATIONS IN PARTICULAR
#' REGIONS OF INTEREST (E.G. GOI HOTSPOTS)
#' Removes rows of the MAF file with mutations for the given GOI that
#' do not fall into the given vector of hotspot positions
#' @param maf_df a MAF file processed by maftools
#' @param goi the Hugo gene symbol of the gene-of-interest
#' @param vect_of_pos a vector of hotspot positions we'd like to keep
filter_by_region <- function(maf_df, goi, vect_of_pos) {
  maf_df_goi <- maf_df[maf_df$Hugo_Symbol == goi,]
  maf_df_noGoi <- maf_df[maf_df$Hugo_Symbol != goi,]
  
  maf_df_goi_sub <- maf_df_goi[maf_df_goi$Start_Position %fin% vect_of_pos,]
  
  maf_df_comb <- rbind(maf_df_noGoi, maf_df_goi_sub)
  return(maf_df_comb)
}


############################################################

#' VISUALIZE THE NUMBER OF PATIENTS WITH MUTATIONS IN EACH
#' PROTEIN
#' Given a MAF file (filtered as desired), a subsetted domains DF
#' and a label, creates barplots/ boxplots of the frequency of binding among
#' all the patients in the MAF file. Also prints some useful info.
#' @param all_gene_id_conv ID conversion doc from BioMart
#' @param domain_df a subsetted domain data frame that has information
#' about gene domains
#' @param label a string denoting the level of specificity ("I-Protein", 
#' "I-Domain", or "I-Binding Position")
#' @param maf_df a MAF table processed by maftools
visualize_mutation_freqs <- function(all_gene_id_conv, domain_df, label, maf_df) {
  
  # Get the protein IDs from the domain DF
  if (label == "I-Protein") {
    protein_ids <- unique(unlist(lapply(domain_df$Query, function(x) unlist(strsplit(x, "|", fixed = TRUE))[2])))
  } else if (label == "I-Domain") {
    protein_ids <- unique(domain_df$Swissprot)
  } else if (label == "I-Binding Position") {
    protein_ids <- unique(domain_df$Protein.ID)
  }
  print(length(protein_ids))

  # For all protein IDs, see which patient tumors have a mutation in it
  freq_df <- data.frame(matrix(ncol = 2, nrow = length(protein_ids)))
  colnames(freq_df) <- c("ID", "Num.Patients")
  freq_df$ID <- protein_ids
  
  # Construct mutation frequency matrix for visualization
  if (label == "I-Protein") {
    maf_df <- data.frame(maf_df)
    freq_df$Num.Patients <- unlist(lapply(protein_ids, function(x) length(unique(maf_df[grepl(x, maf_df$SWISSPROT), 'Tumor_Sample_Barcode']))))
  } else if (label == "I-Domain") {
    freq_df$Num.Patients <- unlist(lapply(protein_ids, function(x) length(unique(domain_df[grepl(x, domain_df$Swissprot),'Patient']))))
  } else if (label == "I-Binding Position") {
    freq_df$Num.Patients <- unlist(lapply(protein_ids, function(x) length(unique(domain_df[domain_df$Protein.ID == x,'Patient']))))
  } else {
    print(paste("Invalid label:", label))
  }
  
  freq_df <- freq_df[order(freq_df$Num.Patients, decreasing = TRUE),]

  # Create barplots
  # ylab_label <- strsplit(strsplit(label, " ", fixed = TRUE)[1], "-", fixed = TRUE)[2]
  # create_barplot(freq_df, label, ylab_label)

  # Create boxplot
  # create_boxplot(freq_df, label)
  
  # Create histogram
  # create_histogram(freq_df, label)
  
  # Create cumulative histogram
  create_cumulative_histogram(freq_df, label)
    
  # Print some interesting output information
  #print(paste("Percent above 1% frequency:", (length(col_freq[col_freq >1])/length(col_freq))*100))
  #print(paste("Number of proteins with >1% frequency:", length(col_freq[col_freq > 1])))
  print(paste("Number of proteins mutated in patients more than 3 times:", length(freq_df$Num.Patients[freq_df$Num.Patients >3])))
  print(paste("Number of proteins mutated in patients more than 6 times:", length(freq_df$Num.Patients[freq_df$Num.Patients >6])))
  
  freq_df_sub <- subset(freq_df, Num.Patients > 3)
  freq_df_sub$Gene.Name <- unlist(lapply(freq_df_sub$ID, function(x) 
    paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'external_gene_name'])), collapse = ";")))
  print(paste(freq_df_sub$Gene.Name, freq_df_sub$Num.Patients, sep = ":"))
  
  return(freq_df_sub)
}

#' Uses a frequency data frame, along with specificity labels,
#' to plot a barplot of mutation frequency across patients
#' @param freq_df a data frame with each protein and its frequency
#' of mutation across all patients
#' @param label a string label denoting the level of specificity
#' ('I-Protein', 'I-Domain', or 'I-Binding Position')
#' @param ylab_label a string label for the y-axis label ('Protein' 
#' for I-Protein, 'Domain' for I-Domain, or 'Position' for 
#' I-Binding Position)
create_barplot <- function(freq_df, label, ylab_label) {
  barplot(freq_df$Num.Patients, main = paste(label, " Number of Patients with Mutations in Each Protein", sep = ":"), 
          xlab = "Protein", ylab = paste("# of Patient Tumor Samp. with >1 Missense Mut. in Ligand-Binding ", 
                                         ylab_label, sep = ""), col = "blue") #, names.arg = rownames(abs_df))
}

#' Uses a frequency data frame, along with a specificity label,
#' to plot a boxplot of mutation frequency across patients
#' @param freq_df a data frame with each protein and its frequency
#' of mutation across all patients
#' @param label a string label denoting the level of specificity
#' ('I-Protein', 'I-Domain', or 'I-Binding Position')
create_boxplot <- function(freq_df, label) {
  bp = boxplot(freq_df$Freq, main = paste(label, " Frequency of Mutation across patients", sep = ":"), ylab = "Frequency",
               horizontal = TRUE, col = "orange")
  identify(rep(1, length(freq_df$Freq)), freq_df$Freq, labels = rownames(freq_df))
  print(paste("Outliers:", bp$out))
}

#' Uses a frequency data frame, along with a specificity label,
#' to plot a histogram of mutation frequency across patients
#' @param freq_df a data frame with each protein and its frequency
#' of mutation across all patients
#' @param label a string label denoting the level of specificity
#' ('I-Protein', 'I-Domain', or 'I-Binding Position')
create_histogram <- function(freq_df, label) {
  hist(freq_df$Num.Patients, main = paste(label, " Frequency of Mutations across Patient Population", sep = ":"), 
       xlab = "Number of Patients with Mutated Version of Protein", ylab = "Number of Genes", 
       col = brewer.pal(n = 11, name = "RdBu"), xlim = c(2, max(freq_df$Num.Patients)))
}


#' Uses a frequency data frame, along with a specificity label,
#' to plot a cumulative histogram of mutation frequency across patients
#' @param freq_df a data frame with each protein and its frequency
#' of mutation across all patients
#' @param label a string label denoting the level of specificity
#' ('I-Protein', 'I-Domain', or 'I-Binding Position')
create_cumulative_histogram <- function(freq_df, label) {
  cum_vals <- c()
  for (i in 1:max(freq_df$Num.Patients)) {
    cum_vals <- c(cum_vals, length(freq_df$Num.Patients[freq_df$Num.Patients >= i]))
  }
  print(cum_vals)
  
  barplot(height = cum_vals[4:20], main = paste(label, " Frequency of Mutations across Patient Population", sep = ":"), 
          xlab = "Minimum # of Patients with Mutated Version of Protein", ylab = "Number of Proteins", col = 'cornflowerblue',
          border = 'black', space = 0, xpd = FALSE, names.arg = c(4:20), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)
}


# Call function
frequency_df_iprotein <- visualize_mutation_freqs(all_gene_id_conv, domains_missense_iprotein_sub, "I-Protein", maf_subset_df_missense)
frequency_df_idomain <- visualize_mutation_freqs(all_gene_id_conv, domains_missense_idomain_sub, "I-Domain")
frequency_df_ibindingpos <- visualize_mutation_freqs(all_gene_id_conv, domains_missense_ibindingpos_sub, "I-Binding Position")



############################################################
#' GET NUMBER OF REMAINING DOMAINS THAT BIND EACH LIGAND TYPE
#' Function takes subsetted data frame of domain information from 
#' Batch-CD search, along with InteracDome and CanBind DFs and a ligand
#' groups conversion document. Extracts the ligand(s) that each domain is 
#' predicted to bind. Prints the number of binding domains for each ligand 
#' type and creates pie chart & bar chart visualizations.
#' @param domains a data frame with information about protein domains
#' @param interacdome_df an InteracDome data frame with information about
#' the ligands that bind particular protein domains
#' @param canbind_df_sub a CanBind data frame with information about 
#' the ligands that bind particular base positions
#' @param ligand_groups a table matching ligand names to their broader
#' categories
get_num_domains_per_ligand <- function(domains, interacdome_df, canbind_df_sub, ligand_groups) {
  
  # Get the domains from the domain DF
  doms <- domains$Accession
  doms_stripped <- unique(as.character(regmatches(doms, gregexpr("[[:digit:]]+", doms))))  # strip the PF from the name
  print(doms_stripped)
  
  # Create a list to hold counts of the results
  domain_types_list <- list("DNA" = 0, "RNA" = 0, "SM" = 0, "ION" = 0, "PEPTIDE" = 0, "NUCACID" = 0)
  
  # Create columns of stripped IDs for InteracDome for easier comparison
  interac_ids_stripped <- unlist(lapply(interacdome_df$pfam_id, function(x) 
    unlist(strsplit(unlist(strsplit(x, "_", fixed = TRUE))[1], "F", fixed = TRUE))[2]))
  interacdome_df$stripped_pfam_ids <- interac_ids_stripped
  
  # For all domains, check the ligand type they bind
  for (i in 1:length(doms_stripped)) {
    domain_types_list <- get_ligand_types(i, interacdome_df, canbind_df_sub, ligand_groups)
  }
  
  # Make this list into a data frame for visualization, ordered by count
  domain_types_table <- data.frame(matrix(unlist(domain_types_list), nrow = length(domain_types_list),
                                          byrow = TRUE), stringsAsFactors = FALSE)
  domain_types_table$Category <- names(domain_types_list)
  colnames(domain_types_table) <- c('Count', 'Category')
  domain_types_table <- domain_types_table[order(-domain_types_table$Count),]
  print(domain_types_table)
  
  # Make a pie chart
  pie(domain_types_table$Count, labels = domain_types_table$Category)
  # Make a bar chart
  barplot(domain_types_table$Count, names.arg = domain_types_table$Category, 
          xlab = "Ligand Type", ylab = "Number of Binding Domains", col=brewer.pal(n = 6, name = "RdBu"),
          main = "Number of domains that bind each ligand type, in top proteins", ylim = c(0,6000))
  print(paste("Total number of unique domains:", length(doms_stripped)))
}


#' CHECK LIGAND TYPES 
#' A helper function that, given an index i, checks the ligand type bound
#' by the domain in the ith row of the domain DF and updates the given
#' domain types list (which it then returns)
#' @param i the row index of the domain DF 
#' @param interacdome_df an InteracDome data frame with information about
#' the ligands that bind particular protein domains
#' @param canbind_df_sub a CanBind data frame with information about 
#' the ligands that bind particular base positions
#' @param domain_types_list a list containing the count of each ligand type 
#' bound by these proteins' domains
get_ligand_types <- function(i, interacdome_df, canbind_df_sub, domain_types_list) {
  # Check InteracDome
  ligand_types_interac <- interacdome_df[interacdome_df$stripped_pfam_ids == doms_stripped[i], 'ligand_type']
  
  # Check CanBind
  prot_id <- unlist(strsplit(domains$Query[i], "|", fixed = TRUE))[2]
  ligand_types_canbind <- canbind_df_sub[canbind_df_sub$Swissprot == prot_id, 'LigandType']
  ligand_types_canbind_list <- unlist(strsplit(ligand_types_canbind, ",", fixed = TRUE))
  
  # Get the union of these two
  ligand_types <- union(ligand_types_interac, ligand_types_canbind_list)
  ligand_types <- ligand_types[!is.na(ligand_types)]
  
  # Update the list with counts for these
  if (!length(ligand_types) == 0) {
    for (j in 1:length(ligand_types)) {
      domain_types_list <- general_important_functions::check_lig_type(ligand_types[j], 
                                                                       domain_types_list, ligand_groups)
    }
  }
  return(domain_types_list)
}


# Import ligand group category file
ligand_groups <- read.csv(paste(path, "Input Data Files/ligand_groups.csv", sep = ""), check.names = FALSE)
# Fix column names, if needed
colnames(ligand_groups)[1] <- "group_name"

# Call function
get_num_domains_per_ligand(domains_missense_iprotein_sub, interacdome_df, canbind_df_sub, ligand_groups)
get_num_domains_per_ligand(domains_missense_idomain_sub, interacdome_df, canbind_df_sub, ligand_groups)
get_num_domains_per_ligand(domains_missense, interacdome_conf, canbind_df_sub, ligand_groups)  # Also look across all domains, for reference


############################################################

#' GET PROP. OF REMAINING PROTEINS THAT BIND EACH LIGAND TYPE
#' Function takes subsetted data frame of domain information from 
#' Batch-CD search, along with InteracDome and CanBind DFs and a ligand
#' groups conversion document. Do not need to include ConCavity, as it 
#' does not make predictions for specific binding types.
#' For each ligand type, creates a pie chart with the proportion of 
#' remaining proteins predicted to bind that ligand in at least
#' one of its domains, according to InteracDome or CanBind
#' @param domains a data frame with information about protein domains
#' @param interacdome_dw_sub domain weights file from InteracDome
#' @param canbind_df_sub a CanBind data frame with information about 
#' the ligands that bind particular base positions
#' @param ligand_groups a table matching ligand names to their broader
#' categories
#' @param label a string denoting the level of specificity ('I-Protein',
#' 'I-Domain', or 'I-Binding Position')
#' @param ibinding_ids vector of regulatory protein swissprot IDs; needed
#' only for 'I-Domain' and 'I-Binding Position' cases
get_num_proteins_that_bind_each_lig_type <- function(domains, interacdome_dw_sub, canbind_df_sub, 
                                                     ligand_groups, label, ibinding_ids) {
  
  # If I-Protein, we can get the Swissprot IDs of all the proteins we're interested in from
  # the domains DF; otherwise import them as a parameter
  if(label == "I-Protein") {
    proteins <- unique(unlist(lapply(domains$Query, function(x) unlist(strsplit(x, "|", fixed = TRUE))[2])))
  } else {
    proteins <- ibinding_ids
  }
  
  # Create columns of stripped IDs for InteracDome for easier comparison
  interac_ids_stripped <- unlist(lapply(interacdome_dw_sub$domain_name, 
                                        function(x) unlist(strsplit(unlist(strsplit(x, "_", fixed = TRUE))[1],
                                                                    "F", fixed = TRUE))[2]))
  interacdome_dw_sub$stripped_pfam_ids <- interac_ids_stripped
  
  # Save the number of proteins that bind each ligand type for plotting
  num_dna_binders <- 0
  num_rna_binders <- 0
  num_pep_binders <- 0
  num_sm_binders <- 0
  num_ion_binders <- 0
  
  # Loop through all the proteins in the remaining group
  for (i in 1:length(proteins)) {
    prot_id <- proteins[i]
    print(prot_id)
    
    # Create a list to hold counts of the resulting ligands it binds
    domain_types_list <- list("DNA" = 0, "RNA" = 0, "SM" = 0, "ION" = 0, "PEPTIDE" = 0, "NUCACID" = 0)
    
    # Get the domains of this protein and strip the "pfam" from the front 
    if (label == "I-Protein") {domains_sub <- domains[grep(prot_id, domains$Query),]} 
    else {domains_sub <- domains[grep(prot_id, domains$Swissprot),]}
    doms <- domains_sub$Accession
    doms_stripped <- unique(as.character(regmatches(doms, gregexpr("[[:digit:]]+", doms))))
    
    # For all domains, check the ligand type they bind
    for (j in 1:length(doms_stripped)) {
      domain_types_list <- get_ligand_types(j, interacdome_dw_sub, canbind_df_sub, domain_types_list)
    }
    if (domain_types_list$DNA > 0) {num_dna_binders <- num_dna_binders + 1}
    if (domain_types_list$RNA > 0) {num_rna_binders <- num_rna_binders + 1}
    if (domain_types_list$PEPTIDE > 0) {num_pep_binders <- num_pep_binders + 1}
    if (domain_types_list$SM > 0) {num_sm_binders <- num_sm_binders + 1}
    if (domain_types_list$ION > 0) {num_ion_binders <- num_ion_binders + 1}
  }

  print(paste("Number of DNA Binders:", num_dna_binders))
  print(paste("Number of RNA Binders:", num_rna_binders))
  print(paste("Number of Peptide Binders:", num_pep_binders))
  print(paste("Number of SM Binders:", num_sm_binders))
  print(paste("Number of Ion Binders:", num_ion_binders))
  
  counts <- list("DNA" = num_dna_binders, "RNA" = num_rna_binders, "PEPTIDE" = num_pep_binders, 
                 "SM" = num_sm_binders, "ION" = num_ion_binders)
  
  # Make this list into a data frame for visualization, ordered by count
  counts_table <- data.frame(matrix(unlist(counts), nrow = length(counts),
                                         byrow = TRUE), stringsAsFactors = FALSE)
  counts_table$Category <- names(counts)
  colnames(counts_table) <- c('Count', 'Category')
  counts_table <- counts_table[order(-counts_table$Count),]
  print(counts_table)
  
  barplot(counts_table$Count, names.arg = counts_table$Category, 
          xlab = "Ligand Type", ylab = "Number of Proteins", col=brewer.pal(n = 6, name = "RdBu"),
          main = "Number of proteins that bind each ligand type", ylim = c(0,length(proteins)))
  
  total_num_prots <- length(proteins)
  print(paste("Total Number of Proteins:", total_num_prots))
  
  prop_dna_binders <- num_dna_binders / total_num_prots
  prop_rna_binders <- num_rna_binders / total_num_prots
  prop_pep_binders <- num_pep_binders / total_num_prots
  prop_sm_binders <- num_sm_binders / total_num_prots
  prop_ion_binders <- num_ion_binders / total_num_prots

  print(paste("Prop. of DNA Binders:", prop_dna_binders))
  print(paste("Prop. of RNA Binders:", prop_rna_binders))
  print(paste("Prop. of Peptide Binders:", prop_pep_binders))
  print(paste("Prop. of SM Binders:", prop_sm_binders))
  print(paste("Prop. of Ion Binders:", prop_ion_binders))
  
  proportions <- list("DNA" = prop_dna_binders, "RNA" = prop_rna_binders, "PEPTIDE" = prop_pep_binders, 
                      "SM" = prop_sm_binders, "ION" = prop_ion_binders)

  # Make this list into a data frame for visualization, ordered by count
  proportions_table <- data.frame(matrix(unlist(proportions), nrow = length(proportions),
                                          byrow = TRUE), stringsAsFactors = FALSE)
  proportions_table$Category <- names(proportions)
  colnames(proportions_table) <- c('Proportion', 'Category')
  proportions_table <- proportions_table[order(-proportions_table$Proportion),]
  print(proportions_table)
  
  barplot(proportions_table$Proportion, names.arg = proportions_table$Category, 
          xlab = "Ligand Type", ylab = "Proportion of Proteins", col=brewer.pal(n = 6, name = "RdBu"),
          main = "Proportion of proteins that bind each ligand type", ylim = c(0,1))
}

# Import necessary files 
ligand_groups <- read.csv(paste(path, "Input Data Files/ligand_groups.csv", sep = ""), check.names = FALSE)
# Fix colnames, if needed
colnames(ligand_groups)[1] <- "group_name"

# Call function
get_num_proteins_that_bind_each_lig_type(domains_missense_iprotein_sub, interacdome_domainweights, canbind_df_sub, ligand_groups, "I-Protein")
get_num_proteins_that_bind_each_lig_type(domains_missense_idomain_sub, interacdome_domainweights, canbind_df_sub, ligand_groups, "I-Domain", ibinding_ids = swissprot_ids_missense_idomain)
get_num_proteins_that_bind_each_lig_type(domains_missense_idomain_sub, interacdome_domainweights, canbind_df_sub, ligand_groups, "I-Binding Position", ibinding_ids = swissprot_ids_missense_ibindingpos)

get_num_proteins_that_bind_each_lig_type(domains_missense, interacdome_conf, canbind_df_sub, ligand_groups)  # Also look across all domains, for reference


############################################################

#' GET MOST MUTATED DRIVER GENES PER CANCER TYPE
#' Uses a regulatory protein mutation data frame created from process_mutation_data.R
#' as well as a list of known driver genes to output a table of the top N most
#' frequently mutated cancer genes in each cancer type
#' @param cancer_types a vector of cancer types we are interested in getting the
#' top mutated cancer genes for
#' @param clinical_df a clinical data frame to link patients to their cancer type
#' @param regprot_mut_df a regulatory protein mutation data frame of any 
#' given specificity
#' @param driver_gene_df a data frame containing known driver genes from 
#' CGC/Vogelstein
#' @param N the number of top genes to retrieve
#' @param all_genes_id_conv ID conversion data frame from BioMart
get_most_mutated_drivers_per_ct <- function(cancer_types, clinical_df, regprot_mut_df, 
                                            driver_gene_df, N, all_genes_id_conv) {
  
  top_mutated_genes_per_ct <- lapply(cancer_types, function(ct) {
    # Get the patients that are of this cancer type
    ct_patients <- unlist(clinical_df[grepl(ct, clinical_df$project_id), 'case_submitter_id'])
    ct_patients <- unlist(lapply(ct_patients, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))
    
    # Get all the genes that have mutations in at least one of these patients
    regprot_mut_df_sub_patients <- data.frame(matrix(nrow = 0, ncol = 2))
    colnames(regprot_mut_df_sub_patients) <- c("Swissprot", "Patient")
    
    for (i in 1:nrow(regprot_mut_df)) {
      print(paste(i, paste("/", nrow(regprot_mut_df))))
      
      # Get the patients with a mutation here
      pats <- regprot_mut_df[i, 'Patient']
      new_pats <- c()
      overlap <- FALSE
      # Check each patient in this cancer type; if they are listed here, keep them,
      # and discard all other patients
      for(j in ct_patients) {
        if(TRUE %in% grepl(j, pats)) {
          overlap <- TRUE
          new_pats <- unique(c(new_pats, j))
        } 
      }
      # If there are no patients from this cancer type in the given row, eliminate it
      if(overlap == TRUE) {
        regprot_mut_df_sub_patients <- rbind(regprot_mut_df_sub_patients, 
                                             data.frame("Swissprot" = regprot_mut_df$Swissprot[i],
                                                        "Patient" = paste(new_pats, collapse = ";")))
      } 
    }
    
    # Aggregate patients for the same protein and get a count matrix
    count_tab <- data.frame(protein = unique(regprot_mut_df_sub_patients$Swissprot),
                            num_patients_with_mutation = rep(NA, length(unique(regprot_mut_df_sub_patients$Swissprot))))
    for (prot in unique(regprot_mut_df_sub_patients$Swissprot)) {
      pats <- regprot_mut_df_sub_patients[grepl(prot, regprot_mut_df_sub_patients$Swissprot),'Patient']
      unique_pats <- unique(unlist(strsplit(pats, ";", fixed = TRUE)))
      count_tab[count_tab$protein == prot, 'num_patients_with_mutation'] <- length(unique_pats)
    }
    
    print(count_tab)
    
    # Filter to only include genes in the driver gene list
    count_tab$ensg <- unlist(lapply(count_tab$protein, function(x)
      paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x,
                                            'ensembl_gene_id'])), collapse = ";")))
    count_tab_list <- lapply(1:nrow(count_tab), function(i) {
      ensg_ids <- unlist(strsplit(count_tab$ensg[i], ";", fixed = TRUE))
      if(any(ensg_ids %fin% driver_gene_df$ensembl_gene_id)) {return(count_tab[i,])}
    })
    count_tab_drivers <- rbindlist(count_tab_list[!is.null(count_tab_list)])
    
    # Get the names of the top N genes and return
    ordered_count_tab <- count_tab_drivers[order(-count_tab_drivers$num_patients_with_mutation),]
    top_n <- ordered_count_tab[1:10,]
    
    return(top_n)
    
  })
  names(top_mutated_genes_per_ct) <- cancer_types
  
  return(top_mutated_genes_per_ct)
}

cancer_types_of_interest <- c("BRCA", "LGG", "THCA", "PRAD", "HNSC", "LIHC", "LUAD",
                              "BLCA", "UCEC", "STAD", "KIRP", "KIRC", "CESC", "LUSC",
                              "COAD", "SARC", "PCPG", "ESCA", "PAAD", "TGCT")

clinical_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/clinical_wMutCounts.csv",
                        header = TRUE, check.names = FALSE)

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"
regprot_mut_df <- read.csv(paste(path, "Mutation/iprotein_results_missense.csv", sep = ""), 
                           header = TRUE, check.names = FALSE)

driver_gene_df <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/GRCh38_driver_gene_list.tsv",
                           header = TRUE, check.names = FALSE)

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


# Call function
top_drivers_per_ct <- get_most_mutated_drivers_per_ct(cancer_types_of_interest, clinical_df, regprot_mut_df,
                                driver_gene_df, 10, all_genes_id_conv)

# Optionally convert to DF and add the gene name
top_drivers_per_ct_df <- as.data.frame(top_drivers_per_ct)
protein_sub <- top_drivers_per_ct_df[,grepl(".protein", colnames(top_drivers_per_ct_df))]
new_cols <- lapply(1:ncol(protein_sub), function(i) {
  col <- protein_sub[,i]
  names <- unlist(lapply(col, function(x) {paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x,
                                          'external_gene_name'])), collapse = ";")}))
  return(names)
})
new_cols_df <- do.call(cbind, new_cols)
colnames(new_cols_df) <- unlist(lapply(cancer_types_of_interest, function(ct) paste(ct, "geneName", sep = ".")))
top_drivers_per_ct_df <- cbind(top_drivers_per_ct_df, new_cols_df)

# Visualize per cancer type using horizontal bar charts
# Reorganize the data frame first, so that cancer type is a column
top_drivers_per_ct_df2 <- do.call(rbind, top_drivers_per_ct)
top_drivers_per_ct_df2$cancer_type <- rep(cancer_types_of_interest, each = 10)
top_drivers_per_ct_df2$gene_name <- unlist(lapply(top_drivers_per_ct_df2$protein, function(x) {
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'external_gene_name'])), 
        collapse = ";")}))

ggplot(top_drivers_per_ct_df2, aes(y = reorder_within(gene_name, by = num_patients_with_mutation, within = cancer_type, sep = "_"), 
                                   x = num_patients_with_mutation)) + geom_bar(stat = "identity") + 
  labs(x = "Driver Gene", y = "Num. Patients with Mutation") + facet_wrap(~cancer_type, scales = "free")


############################################################
#' DETERMINE CORRELATION BETWEEN MUTATION STATUS OF KEY GENES AND CANCER SUBTYPE IN BRCA
#' Creates a heat map that displays the count of patients that have a given subtype 
#' and a corresponding mutation in a particular gene of interest
#' @param subtype_file contains the subtypes for the given cancer along with the 
#' corresponding patient IDs
#' @param mutation_regprot_df a data frame with the mutation status of patients in 
#' the given cancer type for the genes of interest
#' @param genes_of_interest a vector of genes of interest, swissprot IDs
#' @param totOrFrac either "total" or "fraction" to denote whether we want the total
#' count of mutated patients for each gene, or the fraction of patients that have
#' a mutation in this gene
visualize_mutation_corr_to_subtype <- function(subtype_file, mutation_regprot_df, 
                                               genes_of_interest, totOrFrac) {
  
  # Get all the unique subtypes for BRCA
  unique_subtypes <- unique(unlist(subtype_file$BRCA_Subtype_PAM50))
  unique_subtypes <- unique_subtypes[!(unique_subtypes == "NA")]
  patients <- unlist(subtype_file$patient)
  
  # Create the blank count matrix
  count_mat <- data.frame(matrix(nrow = length(genes_of_interest), ncol = length(unique_subtypes)))
  colnames(count_mat) <- unique_subtypes
  rownames(count_mat) <- genes_of_interest
  count_mat[is.na(count_mat)] <- 0
  
  if(totOrFrac == "fraction") {
    total_counts <- list("Basal" = 0, "Her2" = 0, "LumB" = 0, "LumA" = 0,"Normal" = 0)
  }
  
  #print(count_mat)
  
  for (i in 1:length(genes_of_interest)) {
    gene <- genes_of_interest[i]
    print(gene)
    
    for (j in 1:length(patients)) {
      patient <- patients[j]
      subtype <- as.character(unlist(subtype_file[subtype_file$patient == patient, 'BRCA_Subtype_PAM50']))
      print(subtype)
      
      if (!(subtype == "NA")) {
        pat_id <- unlist(strsplit(patient, "-", fixed = TRUE))[3]
        #print(pat_id)
        
        # If this patient has a mutation in this gene, add to count
        mutation_regprot_df_sub <- mutation_regprot_df[mutation_regprot_df$Swissprot == gene,]
        if(any(grepl(pat_id, mutation_regprot_df_sub$Patient))) {
          count_mat[i, colnames(count_mat) == subtype] <- count_mat[i, colnames(count_mat) == subtype] + 1
        }
        if(totOrFrac == "fraction") {
          total_counts[names(total_counts) == subtype] <- as.integer(total_counts[names(total_counts) == subtype]) + 1
        }
      }
    }
  }
  
  # OPTIONAL: divide them all by the number of patients with that subtype to get 
  # a percentage
  if(totOrFrac == "fraction") {
    print(total_counts)
    for (y in 1:ncol(count_mat)) {
      colnam_y <- colnames(count_mat)[y]
      tot_counts <- total_counts[[colnam_y]]
      count_mat[,y] <- round(count_mat[,y] / tot_counts, digits = 5)
    }
  }
  
  print(count_mat)
  
  # Convert from data frame to matrix
  count_mat <- data.matrix(count_mat)
  labels <- apply(count_mat, c(1,2), as.character)
  
  # Visualize using a heat map with counts displayed
  #col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
  #heatmap(count_mat, scale = "none", col = bluered(100))
  #heatmap(count_mat, scale = "none", col = bluered(100), trace = "none",
            #density.info = "none", cellnote = labels, notecol = "black") # cl
  count_mat_melt <- reshape2::melt(count_mat)
  count_mat_melt$Gene.name <- unlist(lapply(count_mat_melt$Var1, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'external_gene_name']), collapse = ";")))
  colnames(count_mat_melt) <- c("Swissprot", "Subtype", "Mutation.Count", "Gene.name")
  p <- ggplot(count_mat_melt, aes(x = Subtype, y = Gene.name, fill = Mutation.Count)) + 
    geom_tile() + geom_text(aes(label = Mutation.Count), color = "black", size = 6) + 
    theme(axis.title.y=element_blank(), axis.text = element_text(size = 14), 
          axis.title.x = element_text(size = 16)) + scale_fill_gradient(low = "dark orange", high = "yellow", guide = "colorbar")
  print(p)
  
  return(count_mat_melt)
}

# Import subtype file using TCGA Biolinks
brca_subtype <- TCGAquery_subtype(tumor = "BRCA")

# Import regulatory protein mutation data frame
mutation_regprot_df <- fread(paste("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Mutation/iprotein_results_missense.csv", sep = ""),
                             header = TRUE)   # I-Protein, all ligands
mutation_regprot_df$Swissprot <- unlist(lapply(mutation_regprot_df$Query, function(x) 
  unlist(strsplit(x, "|", fixed = TRUE))[2]))

# Get the genes of interest
genes_of_interest_brca <- c("FOXA1", "HUWE1", "PTEN", "AKT1", "FAT1", "TP53", 
                            "MAP2K4", "MAP3K1", "KMT2C", "PIK3CA", "ERBB2")
genes_of_interest_brca_swissprot <- c("P42336", "P04637", "P31749", "Q8NEZ4", "Q13233",
                                      "P60484", "Q14517", "P45985", "P04626", "Q7Z6Z7", "P55317") 

count_matrix <- visualize_mutation_corr_to_subtype(brca_subtype, mutation_regprot_df, 
                                                   genes_of_interest_brca_swissprot,
                                                   "total")
count_matrix_frac <- visualize_mutation_corr_to_subtype(brca_subtype, mutation_regprot_df, 
                                                   genes_of_interest_brca_swissprot,
                                                   "fraction")


############################################################
#' GET A FILE OF PATIENTS THAT POSSESS A GIVEN CANCER SUBTYPE
#' Creates header-less text file that has a simple line-separated list of all the
#' patients (4-letter TCGA ID, ####) that are a member of the provided subtype(s)
#' @param subtype_file contains the subtypes for the given cancer along with the 
#' corresponding patient IDs
#' @param subtypes_of_interest a list of all the subtypes of interest we want to
#' put in the output file
#' @param outpath a path to where we should write the new file
#' @param subtypes_label a character label of what subtype(s) we're looking at in
#' order to appropriately name the new file
get_subtype_patient_list <- function(subtype_file, subtypes_of_interest, outpath,
                                     subtypes_label) {
  output_patient_ids <- c()
  
  for (subtype in subtypes_of_interest) {
    new_ids <- subtype_file[subtype_file$BRCA_Subtype_PAM50 == subtype, "patient"]  #TODO: make this generalizable beyond BRCA
    new_ids_short <- unlist(lapply(new_ids$patient, function(id) unlist(strsplit(id, "-", fixed = TRUE))[3]))
    output_patient_ids <- c(output_patient_ids, new_ids_short)
  }
  
  fwrite(as.data.table(output_patient_ids), paste0(outpath, paste0(subtypes_label, "_patient_ids.txt")))
}


# Import subtype file using TCGA Biolinks
brca_subtype <- TCGAquery_subtype(tumor = "BRCA")

# List the outpath
outpath <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"

# Call function
get_subtype_patient_list(brca_subtype, c("LumA", "LumB"), outpath, "Luminal.A.B")
get_subtype_patient_list(brca_subtype, "LumA", outpath, "Luminal.A")


############################################################
#' DETERMINE THE NUMBER OF FUNCTIONAL COPIES OF A GIVEN GOI FOR EACH PATIENT SAMPLE
#' Given a gene-of-interest, a MAF file, and a CNA file (with absolute CNs), determine
#' for each patient sample the number of functional copies they have of the given GOI
#' @param maf_file a TCGA MAF file, subsetted to only missense and nonsense mutations
#' @param cna_file a TCGA CNA file, from the ASCAT2 pipeline
#' @param patient_ids a set of 4-digit patient ids (XXXX) we are considering
#' @param goi the ensembl (ENSG) ID of a gene of interest
#' @param ts_or_onco a character string denoting whether the goi is a "tumor_suppressor"
#' or an "oncogene"
get_num_functional_copies <- function(maf_file, cna_file, patient_ids, goi, ts_or_onco) {
 
  # Limit the MAF and CNA files to just the gene of interest
  maf_sub <- maf_file[maf_file$Gene == goi,]
  cna_sub <- cna_file[rownames(cna_file) == goi,]
  
  # Limit the MAF and CNA files to just the patients of interest
  maf_sub$Patient_ID <- unlist(lapply(maf_sub$Tumor_Sample_Barcode, function(x)
    unlist(strsplit(x, "-", fixed = TRUE))[3]))
  maf_sub <- maf_sub[maf_sub$Patient_ID %in% patient_ids,]
  
  cna_pats <- unlist(lapply(colnames(cna_sub), function(x)
    unlist(strsplit(x, "-", fixed = TRUE))[1]))
  cna_sub <- cna_sub[,which(cna_pats %in% patient_ids)]
  
  # Create a new DF with all these patients
  output_df <- data.frame("sample.id" = colnames(cna_sub))
  
  # Check what the copy number of the GOI is for each patient
  output_df$Abs.CN <- unlist(lapply(1:nrow(output_df), function(i) {
    samp <- output_df$sample.id[i]
    return(cna_sub[,colnames(cna_sub) == samp])
  }))
  
  #print(output_df)
  
  # Use helper function to get the number of copies with a missense and/or nonsense mutation
  output_df$Num.Missense.Mut.Copies <- get_num_mut_copies(output_df, maf_sub, "Missense_Mutation")
  output_df$Num.Nonsense.Mut.Copies <- get_num_mut_copies(output_df, maf_sub, "Nonsense_Mutation")
  output_df$Num.Mut.Copies <- get_num_mut_copies(output_df, maf_sub, c("Missense_Mutation", "Nonsense_Mutation"))
  
  print(output_df)
  
  output_df$Num.Functional.Copies <- unlist(lapply(1:nrow(output_df), function(i) {
    # If a total deletion (CNA of 0), we have no functional copies -- return 0
    if(output_df[i, 'Abs.CN'] == 0) {return(0)}
    # If we have one copy deleted, then we know we're down at least one functional copy;
    # need to check and see if the other has a mutation (if it's a T.S.)
    else if(output_df[i, 'Abs.CN'] == 1) {
      if (ts_or_onco == "tumor_suppressor") {
        if(output_df[i, 'Num.Mut.Copies'] >= 1) {return(0)}
        else {return(1)}
      # If this is an oncogene, we only need to check for nonsense mutations
      } else {
        if(output_df[i, 'Num.Nonsense.Mut.Copies'] >= 1) {return(0)}
        else {return(1)}
      }
    # If we have at least 2 non-deleted copies, we just have to check for mutation status 
    } else {
      if (ts_or_onco == "tumor_suppressor") {
        # If we have at least 1 mutated allele, subtract this from the total CN
        if(output_df[i, 'Num.Mut.Copies'] >= 1) {
          return(as.integer(output_df[i, 'Abs.CN']) - as.integer(output_df[i, 'Num.Mut.Copies']))
        }
        # If we have no mutated alleles, and no deletions, we have 2 (or more) functional copies
        else {as.integer(output_df[i, 'Abs.CN'])}
      } else {
        # If this is an oncogene, we have to check whether the mutation is nonsense
        # (If missense, likely GOF, and if nonsense, LOF)
        if(output_df[i, 'Num.Nonsense.Mut.Copies'] >= 1) {
          return(as.integer(output_df[i, 'Abs.CN']) - as.integer(output_df[i, 'Num.Nonsense.Mut.Copies']))
        }
        else {as.integer(output_df[i, 'Abs.CN'])}
      }
    }
  }))
  
  return(output_df)
}

#' For each sample, return the allele-specific mutation status for the GOI
#' @param output_df a data frame we are filling in 
#' @param maf_sub a MAF file DF subsetting to a given patient and GOI
#' @param classification either "Missense_Mutation" or "Nonsense_Mutation" to 
#' to indicate the type of mutation we are looked at the mutation status for
get_num_mut_copies <- function(output_df, maf_sub, classification) {
  # Iterate through all the samples and, for each, get the number of mutant copies
  # of a given classification
  num_mut_copies <- unlist(lapply(1:nrow(output_df), function(i) {
    samp <- output_df$sample.id[i]
    
    # Subset the MAF data frame to only this sample and variant type
    if(length(classification) == 1) {
      maf_samp <- maf_sub[(grepl(samp, maf_sub$Tumor_Sample_Barcode)) & 
                            (maf_sub$Variant_Classification == classification),]
    } else {
      maf_samp <- maf_sub[(grepl(samp, maf_sub$Tumor_Sample_Barcode)) & 
                            (maf_sub$Variant_Classification %in% classification),]
    }
    
    num_mutant_copies <- 0
    print(nrow(maf_samp))
    
    # If there's no remaining rows (aka mutations), we have no mutated copies -- return 0
    if (nrow(maf_samp) == 0) {return(num_mutant_copies)}
    
    # Otherwise, look at the individual alleles
    else {
      #print(head(maf_samp))
      mut_pattern <- list("allele1" = 0, "allele2" = 0)
      
      # If there are multiple rows of mutations, iterate through them to see which allele
      # they are specific to
      for (j in 1:nrow(maf_samp)) {
        ref_allele <- maf_samp[j, 'Reference_Allele']
        print(ref_allele)
        # Is the first allele mutated?
        if (maf_samp[j, 'Tumor_Seq_Allele1'] != ref_allele) {
          # Have we already found a mutation of this classification on this allele?
          if(mut_pattern[[1]] != 1) {
            num_mutant_copies <- num_mutant_copies + 1
            mut_pattern[[1]] <- mut_pattern[[1]] + 1
          }
        }
        # Is the second allele mutated?
        if (maf_samp[j, 'Tumor_Seq_Allele2'] != ref_allele) {
          # Have we already found a mutation of this classification on this allele?
          if(mut_pattern[[2]] != 1) {
            num_mutant_copies <- num_mutant_copies + 1
            mut_pattern[[2]] <- mut_pattern[[2]] + 1
          }
        }
      }
      return(num_mutant_copies)
    }
  }))
  return(num_mut_copies)
}


input_path <-  "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/"
maf_df_misAndNon <- read.csv(paste0(input_path, "Somatic_Mut_Data/maf_file_df_missense_nonsense_iprotein.csv"), 
                             header = TRUE, check.names = FALSE)
output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
cna_df <- read.csv(paste0(output_path, "CNV/Gene-level Raw/CNV_DF_AllGenes_CancerOnly.csv"), 
                   header = TRUE, check.names = FALSE, row.names = 1)

patient_ids <- read.table(paste0(output_path, "Linear Model/Tumor_Only/intersecting_ids.txt"), header = TRUE)[,1]

goi <- "ENSG00000141510" # TP53
goi <- "ENSG00000121879" # PIK3CA

ts_or_onco <- "tumor_suppressor"
ts_or_onco <- "oncogene"


# Call function
output_df <- get_num_functional_copies(maf_df_misAndNon, cna_df, patient_ids, goi, ts_or_onco)

# Visualize the distribution of functional copies
hist(output_df$Num.Functional.Copies, breaks = seq(min(output_df$Num.Functional.Copies)-0.5, 
                                                   max(output_df$Num.Functional.Copies)+0.5, by =1),
     xlab = "Num. Functional Copies", main = "TP53 Histogram of Num. Functional Copies (BRCA)")
# Visualize when we log2(F + 1) the values
new_vals <- log2(output_df$Num.Functional.Copies + 1)
hist(new_vals, breaks = seq(min(new_vals), max(new_vals)+0.5, by =1),
     xlab = "Log2(Num. Functional Copies + 1)", main = "TP53 Histogram of Log2(Num. Functional Copies + 1) (BRCA)")


# Create a Venn diagram of patients with both a TP53 deletion and mutation
ggVennDiagram(list("P53 Deletion" = output_df[output_df$Abs.CN %in% c(0,1), 'sample.id'], "P53 Mutation" = output_df[output_df$Num.Mut.Copies >= 1, 'sample.id']), 
              label_alpha = 0, category.names = c("TP53 Deletion", "TP53 Mutation"), set_color = "black",
              set_size = 10, label_size = 8, edge_size = 0) +
  ggplot2::scale_fill_gradient(low="cornsilk1", high = "cadetblue3")


# How many have no functional copies?
barplot(unlist(list("Total P53 KO" = length(output_df[output_df$Num.Functional.Copies == 0, 'sample.id']), 
                   "Partial P53 KO (1 Cp)" = length(output_df[output_df$Num.Functional.Copies == 1, 'sample.id']),
                   ">=2 Functional P53 Cp" = length(output_df[output_df$Num.Functional.Copies >= 2, 'sample.id']))))


# Write this to a file
write.csv(output_df, paste0(output_path, "Linear Model/Tumor_Only/tp53_functional_copies_per_sample.csv"))

#
#
#
#
#
#
#
#
#### SAMPLE MAIN FOR USING THE CODE ABOVE ###

# Use maftools to examine the proportion of patients that possess top mutations
maf_filename <- "TCGA.BRCA.muse.b8ca5856-9819-459c-87c5-94e91aca4032.DR-10.0.somatic.maf"
#maf_filename <- "TCGA.Aggregate.muse.aggregated.somatic.maf"


maf <- read.maf(maf = paste0(main_path, maf_filename))  # Read in maf file using maftools

# Use subsetMAF to also create subsetted MAF objects
subset_maf <- subsetMaf(maf, genes = protein_ids, mafObj = TRUE, isTCGA = TRUE)
subset_maf_missense <- subsetMaf(maf, genes = protein_ids_missense, mafObj = TRUE, 
                                          isTCGA = TRUE, query = "Variant_Classification == 'Missense_Mutation'")
# maftools_subset_maf_silent <- subsetMaf(maftools_maf, genes = protein_ids_silent, mafObj = TRUE, 
# isTCGA = TRUE, query = "Variant_Classification == 'Silent'")

# We want to use this function to summarize the full MAF file and also a subsetted
# MAF file that includes ONLY our transcription factors of interest and patients of interest
full_maf_pfam_genes <- summarize_maf(maf)
full_maf_pfam_doms <- get_pfam_domains(maf, 6)

subset_maf_pfam_genes <- summarize_maf(subset_maf)
subset_maf_pfam_doms <- get_pfam_domains(subset_maf, 6)

subset_maf_pfam_missense_genes <- summarize_maf(subset_maf_missense)
subset_maf_pfam_missense_doms <- get_pfam_domains(subset_maf_missense, 6)

# Get a matrix of all the mutations per gene per sample (for full MAF file and also
# for subsetted MAF files)
mut_count_matrix <- get_mut_count_matrix(maf = maf)
mut_count_matrix_subset <- get_mut_count_matrix(maf = subset_maf)
mut_count_matrix_missense <- get_mut_count_matrix(maf = subset_maf_missense)
# TODO: differentiate between tumor and normal samples

# Save as CSV output
write.csv(mut_count_matrix, file = "mut_count_matrix.csv")
write.csv(mut_count_matrix_subset, file = "mut_count_matrix_subset.csv")
write.csv(mut_count_matrix_missense, file = "mut_count_matrix_subset_missense.csv")
#write.csv(mut_count_matrix_silent, file = "mut_count_matrix_silent.csv")


# Get binding domains for particular proteins
unique(domains_missense_idomain_sub[grep("P60484", domains_missense_idomain_sub$Swissprot),'Accession'])
