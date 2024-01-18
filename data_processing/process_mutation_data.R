############################################################
### Process Mutation Data 
### ASSOCIATED PUBLICATION INFORMATION
############################################################

library(data.table)
library(GenomicRanges)
library(TRONCO)
library(stringr)
library(dplyr)
library("RColorBrewer")
library(maftools)
library(stats)


# TCGA Workflow: https://www.bioconductor.org/packages/release/data/experiment/vignettes/TCGAWorkflowData/inst/doc/TCGAWorkflowData.html
# TRONCO Package Vignette: https://bioconductor.org/packages/devel/bioc/vignettes/TRONCO/inst/doc/vignette.pdf

# Useful information about MAF file formats: 
# GDC specified format with column descriptors: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
# Biostars MAF caveats: https://www.biostars.org/p/69222/

# Import a MAF file for cancer type of interest, or pan-cancer and filter:
# 1. all patients that do not have all data types of interest
# 2. patients who are "hypermutators"

# Local PATH to directory containing MAF data file
PATH <- getwd()
OUTPUT_PATH <- paste0(getwd(), "Output/Mutation/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "Input Data Files/all_genes_id_conv.csv"), 
                              header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
# GET DNA MUTATION INFORMATION (MAF FILE) 
############################################################
# GET THE MAF FILENAME
# Load TCGA MAF file of interest (do not merge mutation types to "Mutation")
maf_filename <- paste0(PATH, "Mutation/Individual_PanCancer_MAFs/
                       TCGA.BRCA.muse.b8ca5856-9819-459c-87c5-94e91aca4032.DR-10.0.somatic.maf")
maf_filename <- paste(PATH, "Mutation/TCGA.Aggregate.muse.aggregated.somatic.maf")
# For METABRIC
maf_filename <- paste0(PATH, "Mutation/brca_metabric/data_mutations.txt")

# For pan-cancer, upload MAF file for each cancer type to separately filter
# hypermutators by cancer type
maf_filenames <- list.files(paste0(PATH, "Mutation/Individual_PanCancer_MAFs/"), 
                            pattern = ".somatic")


# IMPORT THE MAF FILE USING TRONCO
# This will return the MAF file in TRONCO format; set to.TRONCO to F to 
# import as DF
maf_file_df_unfilt <- import.MAF(maf_filename, sep = "\t", is.TCGA = T, 
                                 merge.mutation.types = F, to.TRONCO = F) 
# For METABRIC, set is.TCGA to FALSE
maf_file_df_unfilt <- import.MAF(maf_filename, sep = "\t", is.TCGA = F, 
                                 merge.mutation.types = F, to.TRONCO = F) 
    # Can use "filter.fun" param to filter MAF lines using a homemade function

# If we're using individual MAF files....
maf_file_dfs_unfilt <- lapply(maf_filenames, function(x) 
  import.MAF(paste0(PATH, paste0("Mutation/Individual_PanCancer_MAFs/", x)), 
             sep = "\t", is.TCGA = T, merge.mutation.types = F, to.TRONCO = F))


############################################################
# FILTER TO INCLUDE ONLY PATIENTS WITH ALL 4 DATA TYPES
############################################################
# NOTE: This assumes that the clinical file has already been processed; make sure that 
# process_clinical_data.R has been run through properly first!
clinical_df <- read.csv(paste0(PATH, "Mutation/clinical_data.csv"), header = T)

# For METABRIC
clin_filename <- paste0(PATH, "Mutation/brca_metabric/
                        data_clinical_sample_sub.txt")
clinical_df <- read.table(clin_filename, header = T, sep = ",")

# Filter by case ID
maf_file_df <- maf_file_df_unfilt[maf_file_df_unfilt$case_id %fin% 
                                    clinical_df$case_id,]

# For pan-cancer with all MAF files, filter each by case IDs
maf_file_dfs <- lapply(maf_file_dfs_unfilt, function(df) 
  df[df$case_id %fin% clinical_df$case_id,])


# Write to file, as needed
write.csv(maf_file_df, paste0(PATH, "Mutation/Individual_PanCancer_MAFs/
                              maf_file_hypermut_filt.csv"))


############################################################
# OPTIONAL: KEEP ONLY MUTATIONS IN PARTICULAR POSITIONS FOR 
# A PARTICULAR GENE OF INTEREST
############################################################
#' FOR A PARTICULAR DRIVER, KEEP ONLY MUTATIONS IN PARTICULAR
#' REGIONS OF INTEREST (E.G. GENE-OF-INTEREST (GOI) HOTSPOTS)
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

goi <- "PIK3CA"
vect_of_pos <- c(179234297, 179218303, 179218294, 179203765)
maf_file_df <- filter_by_region(maf_file_df, goi, vect_of_pos)


############################################################
# GENERATE A MUTATION COUNT MATRIX FROM MAF FILES
############################################################
#' GET MUTATION COUNT MATRIX
#' This function takes a MAF file that has been read by 
#' maftools and returns a mutation-count matrix
#' THIS ONLY INCLUDES NONSYNONYMOUS VARIANTS unless includeSyn 
#' is specified as TRUE. This only includes primary solid tumor samples.
#' @param maf a MAF file that has been read in by maftools
get_mut_count_matrix <- function(maf, includeSyn) {
  
  if(includeSyn == T) {
    # An explanation of the difference between a "Splice_Region" variant and a 
    # "Splice_Site" variant: https://github.com/mskcc/vcf2maf/issues/63#issuecomment-236981452
    # Essentially, "splice site" mutations are more likely to alter splicing 
    # (within 2bp of the splice site), while "splice region" mutations are more 
    # uncertain in their effect on splicing (3-8bp from the splice site)
    mut_count_matrix <- mutCountMatrix(maf = maf, 
                                       countOnly = c("Missense_Mutation", 
                                                     "Nonsense_Mutation",
                                                     "Nonstop_Mutation",
                                                     "Splice_Site"), 
                                       removeNonMutated = FALSE)
  } else {
    mut_count_matrix <- mutCountMatrix(maf = maf, countOnly = "Silent", 
                                       includeSyn = TRUE, 
                                       removeNonMutated = FALSE)
  }
  return(mut_count_matrix)
}

# Generate mutation count matrix and cancer patient-ID map
mut_count_matrix <- get_mut_count_matrix(maf_file_df)


############################################################
# MODIFY THE FILTERED MUTATION COUNT MATRIX TO ADD MISSING GENES (NO MUTATIONS 
# IN ANY PATIENT AND PATIENTS (NO MUTATION IN ANY GENE)
############################################################
# Add rows of 0's to this matrix for genes that are not included, because there
# are no patients that have a mutation in this gene

#' Function to add genes with no mutations across any sample to the mutation 
#' count matrix
#' @param mut_count_matrix a maftools-produced mutation count matrix
#' @param allgene_targets a list of all target genes created from Ensembl
add_missing_genes_to_mut_count_mat <- function(mut_count_matrix, 
                                               allgene_targets) {
  # Get ENSG IDs from gene symbols
  mut_count_mat_ensgs <- unlist(lapply(rownames(mut_count_matrix), function(x) 
    unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 
                                    'ensembl_gene_id']))[1]))
  
  # Get the genes that are missing
  missing_genes <- setdiff(allgene_targets$ensg, mut_count_mat_ensgs)
  
  # Create a table for them and append to mutation count matrix
  missing_genes_df <- data.table(matrix(nrow = length(missing_genes), 
                                        ncol = ncol(mut_count_matrix)))
  missing_genes_df[is.na(missing_genes_df)] <- 0
  rownames(missing_genes_df) <- make.names(unlist(lapply(missing_genes, function(x) 
    unique(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == x, 
                                    'external_gene_name']))[1])), unique = TRUE)
  colnames(missing_genes_df) <- colnames(mut_count_matrix)
  
  mut_count_matrix_full <- rbind(mut_count_matrix, missing_genes_df)
  rownames(mut_count_matrix_full) <- c(rownames(mut_count_matrix), 
                                       rownames(missing_genes_df))
  
  return(mut_count_matrix_full)
}

# Import Ensembl file of all human gene IDs
allgene_targets <- read.csv(paste0(PATH, "allgene_targets.csv"), header = T, 
                            row.names = 1, check.names = F)

# Call function
mut_count_matrix <- add_missing_genes_to_mut_count_mat(mut_count_matrix, 
                                                       allgene_targets)


#' Function to add patients with no mutation of the the given specificity in 
#' any protein
#' @param mut_count_matrix a maftools-produced mutation count matrix
#' @param intersecting_patients a list of patients that have all data types 
#' within the given cancer type, before filtering is done
add_missing_patients_to_mut_count_mat <- function(mut_count_matrix, 
                                                  intersecting_patients) {
  
  colnames(mut_count_matrix) <- unlist(lapply(colnames(mut_count_matrix), function(x)
    paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))
  
  missing_patients <- setdiff(intersecting_patients, 
                              unlist(lapply(colnames(mut_count_matrix), function(x) 
                                unlist(strsplit(x, "-", fixed = TRUE))[1])))
  new_pat_df <- data.frame(matrix(nrow = nrow(mut_count_matrix), 
                                  ncol = length(missing_patients)))
  new_pat_df[is.na(new_pat_df)] <- 0
  colnames(new_pat_df) <- paste0(missing_patients, "-01A")  #add a placeholder sample ID
  
  new_mut_count_mat <- cbind(mut_count_matrix, new_pat_df)
  return(new_mut_count_mat)
}

# Call functions
mut_count_matrix_full <- add_missing_genes_to_mut_count_mat(mut_count_matrix, 
                                                            allgene_targets)

# Import the list of unique patient IDs that have all data types 
patient_id_list <- read.table(paste0(PATH, "UNIQUE_PATIENT_IDS.txt"), 
                              header = T)[,1]
mut_count_matrix_full <- add_missing_patients_to_mut_count_mat(mut_count_matrix_full, 
                                                               patient_id_list)

# Write to file
write.csv(mut_count_matrix_full, paste0(PATH, "Mutation/Mutation Count Matrices/
                                        mut_count_matrix_nonsynonymous_ALL_inclNonmut.csv"))



############################################################
# FILTER OUT HYPERMUTATORS
############################################################
#' FILTER HYPERMUTATORS
#' Given a MAF file name and the associated data frame:
#' 1. Reads in the MAF using maftools
#' 2. Uses helper function to get nonsynonymous mutation-count matrix
#' 3. Determines outliers and removes them
#' 4. Returns filtered MAF DF and the filtered mutation-count matrix
#' @param maf_filename the name of the maf file, including its PATH
#' @param maf_df a MAF file processed by maftools
#' @param method either "IQR" or "mean" to denote the method of hypermutator
#' removal
filter_hypermutators <- function(maf_filename, maf_df, method) {
  
  maf <- read.maf(maf_filename)
  
  mut_count_matrix <- get_mut_count_matrix(maf)
  total_mut_count_tab <- get_num_mutations_per_patient(mut_count_matrix)
  total_mut_counts <- total_mut_count_tab[,1]
  
  print(length(total_mut_counts))
  
  # Get the IQR
  if(method == "iqr") {
    iqr <- IQR(total_mut_counts)
    # Get Q3 of data
    q3 <- as.numeric(quantile(total_mut_counts)[4])
    # Get Q3 + 1.5(IQR)
    thres <- as.numeric(q3 + (1.5 * iqr))
    print(thres)
  } else if (method == "mean") {
    mean <- mean(total_mut_counts)
    sd <- sd(total_mut_counts)
    thres <- mean + (sd*3)
    print(thres)
  } else {
    print("Only implemented with 'IQR' and 'mean' methods. Try again with one of
          these inputs.")
    return(NA)
  }

  # Filter the mutation count table by this threshold
  total_mut_count_tab_filt <- subset(total_mut_count_tab, Total.Num.Mut < thres)
  
  # Get the remaining patients left in the matrix
  remaining_samp <- rownames(total_mut_count_tab_filt)
  
  # Filter the MAF DF to include only these patients & return
  maf_df_filt <- maf_df[maf_df$Tumor_Sample_Barcode %fin% remaining_samp,]
  
  # Optional histogram
  #hist(total_mut_count_tab_filt$Total.Num.Mut, xlab = "Number of Variants", 
    #main = "")
  
  return(list(maf_df = maf_df_filt, mut_count_df = total_mut_count_tab_filt))
}

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
  
  return(total_mut_count_df)
}


# Call these functions to filter out hypermutators (those who are outliers in 
# the number of mutations they possess in their tumor sample). 
hypermut_res <- filter_hypermutators(maf_filename, maf_file_df)
maf_file_df <- hypermut_res[[1]] 
total_mut_count_df <- hypermut_res[[2]]  # The filtered mutation count matrix

# If we're using multiple MAF files, filter hypermutators separately on each 
# and then aggregate (this is the case for pan-cancer analysis)
total_mut_count_dfs <- list()
for (i in 1:length(maf_filenames)) {
  maf_filename <- maf_filenames[i]
  cancer_type_label <- unlist(strsplit(maf_filename, split = ".", fixed = T))[2]
  maf_df <- maf_file_dfs[[i]]
  hypermut_res <- filter_hypermutators(paste0(PATH, paste0("Mutation/Individual_PanCancer_MAFs/", 
                                                          maf_filename)), maf_df)
  maf_file_dfs[[i]] <- hypermut_res[[1]]
  write.csv(hypermut_res[[2]], paste0(PATH, paste0("Mutation/Mutation Count Matrices/
                                                  total_mut_count_matrix_per_patient_", 
                                                  paste0(cancer_type_label, ".csv"))))
  total_mut_count_dfs[[cancer_type_label]] <- hypermut_res[[2]]
}
maf_file_df <-rbindlist(maf_file_dfs, use.names = T)
total_mut_count_df <- do.call(rbind, unname(total_mut_count_dfs))


# Write to file
write.csv(total_mut_count_df, paste0(PATH, "Mutation/Mutation Count Matrices/
                                     total_mut_count_matrix_per_patient.csv"))
write.csv(total_mut_count_df, paste0(PATH, "METABRIC/Mutation/
                                     Mutation Count Matrices/
                                     total_mut_count_matrix_per_patient.csv"))

write.csv(maf_file_df, paste0(PATH, "Mutation/maf_file_hypermut_filt_indiv.csv"))
write.csv(total_mut_count_df, paste0(PATH, "Mutation/
                                     total_mut_count_df_hypermut_filt_indiv.csv"))

#' Additional option when we have manually determined thresholds from histograms 
#' of the total number of mutations or some other external source; returns a 
#' filtered mutation count matrix
#' @param mut_count_matrix a pan-cancer mutation count matrix
#' @param pc_hypermutator_thres a table of thresholds for total # of mutations 
#' for each cancer type
#' @param cancer_patient_id_type_map a mapping of patient IDs to cancer type
filter_hypermutators_manual <- function(mut_count_matrix, pc_hypermutator_thres, 
                                        cancer_patient_id_type_map) {
  cts <- unique(cancer_patient_id_type_map$Project.ID)
  mut_count_pats <- unlist(lapply(colnames(mut_count_matrix), function(x) 
    unlist(strsplit(x, "-", fixed = T))[1]))
  
  df_new_list <- lapply(cts, function(ct) {
    ct_pats <- unlist(cancer_patient_id_type_map[
      cancer_patient_id_type_map$Project.ID == ct, 'Patient.ID'])
    mut_ct <- mut_count_matrix[, which(mut_count_pats %in% ct_pats)]
    
    # Get the total number of mutations for each of these patients
    mut_colsums <- colSums(mut_ct)
    
    # Get the threshold for this CT
    thres <- pc_hypermutator_thres[pc_hypermutator_thres$Cancer.Type == ct, 
                                   'Hypermut.Thres']
    
    # Limit to those samples below this threshold
    mut_sub <- mut_ct[, which(mut_colsums < thres)]
    
    return(mut_sub)
  })
  
  mut_count_matrix_sub <- do.call(cbind, df_new_list)
  
  total_mut_count_matrix <- get_num_mutations_per_patient(mut_count_matrix_sub)
  total_mut_count_matrix_full <- visualize_mutation_distrib(
    total_mut_count_matrix, cancer_patient_id_type_map)
  
  # Optional visualization, by cancer type
  #g <- ggplot(total_mut_count_matrix_full, aes(x=Total.Num.Mut)) + 
    #geom_histogram() + labs(x = "Total Mutation Count", y = "Frequency") + 
    #facet_wrap(~Project.ID, scales ="free")
  #print(g)
  
  return(mut_count_matrix_sub)
}


#' Helper function that, using a list of patient IDs and a clinical DF from
#' the TCGA, creates a mapping between patient IDs and their cancer type
#' (since cancer type is not listed in aggregated MAF files)
#' @param patient_ids a vector of four-digit patient IDs from the TCGA
#' @param clinical_df a clinical DF produced by the TCGA for the given cohort
make_patient_id_cancer_type_map <- function(patient_ids, clinical_df) {
  
  rows <- lapply(patient_ids, function(id) {
    # Extract the "project ID" for each patient ID and split out the "TCGA-" part
    project_id <- unlist(strsplit(unique(clinical_df[grepl(
      id, clinical_df$case_submitter_id), 'project_id']), "-", fixed = TRUE))[2]
    if(!length(project_id) == 0) {
      return(c(id, project_id))
    }
  })
  rows <- rows[!is.na(rows)]
  df <- as.data.frame(do.call("rbind", rows))
  colnames(df) <- c("Patient.ID", "Project.ID")
  return(df)
}


#' Given a maf data frame and corresponding clinical data frame, 
#' plots the mutation counts as a histogram per cancer type
#' @param total_mutation_count_matrix a mutation-count matrix produced
#' by maftools
#' @param patient_id_cancer_type_map a clinical data frame produced
#' by the 'make_patient_id_cancer_type_map()' function
visualize_mutation_distrib <- function(total_mutation_count_matrix, 
                                       patient_id_cancer_type_map) {
  
  # For every cancer type, get a separate visualization of the mutation
  cancer_types <- unique(patient_id_cancer_type_map[,'Project.ID'])
  mut_count_matrices <- lapply(cancer_types, function(ct) {
    
    # Subset the mapping to just this cancer type & get all the 
    # patients of this type
    patients <- patient_id_cancer_type_map[
      patient_id_cancer_type_map$Project.ID == ct, 'Patient.ID'] 

    # Use these patients to subset the total count matrix 
    count_matrix_sub <- total_mutation_count_matrix[
      grepl(paste(patients, collapse = "|"), 
            rownames(total_mutation_count_matrix)), , drop = FALSE]
    
    # Add a column for the cancer type
    count_matrix_sub$Project.ID <- rep(ct, nrow(count_matrix_sub)) 
    return(as.data.frame(count_matrix_sub))
  }) 
  
  # Re-combine these all into one big table for plotting
  total_mut_count_tab <- do.call("rbind", mut_count_matrices)
  
  # Plot the total mutation counts as a histogram for each table
  g <- ggplot(total_mut_count_tab, aes(x=Total.Num.Mut)) + geom_histogram() +
    labs(x = "Total Mutation Count", y = "Frequency") + facet_wrap(~Project.ID)
  print(g)

  return(total_mut_count_tab)
}


# Generate patient ID map
patient_ids <- unique(unlist(lapply(colnames(mut_count_matrix), function(x) 
  unlist(strsplit(x, "-", fixed = TRUE))[1])))
cancer_patient_id_type_map <- make_patient_id_cancer_type_map(patient_ids, 
                                                              clinical_df)

# Generate histograms of mutational distributions and use to generate manual
# per-cancer hypermutator thresholds
total_mut_count_tab <- visualize_mutation_distrib(mut_count_matrix, 
                                                  cancer_patient_id_type_map)
write.csv(total_mut_count_tab, paste0(PATH, "Mutation/Mutation Count Matrices/
                                      total_mut_count_matrix_per_patient.csv"))

# Read in created manual or uniform thresholds from the above or other sources
pc_hypermutator_thres_manual <- read.csv(paste0(PATH, "Mutation/
                                                manual_hypermutator_thresholds.csv"))
pc_hypermutator_thres_uniform <- read.csv(paste0(PATH, "Mutation/
                                                 uniform_hypermutator_thresholds.csv"))

# Call function to filter hypermutators using a manual external list
mut_count_matrix_manual_hmFilt <- filter_hypermutators_manual(mut_count_matrix, 
                                                              pc_hypermutator_thres_manual, 
                                                              cancer_patient_id_type_map)
mut_count_matrix_uniform_hmFilt <- filter_hypermutators_manual(mut_count_matrix, 
                                                               pc_hypermutator_thres_uniform, 
                                                               cancer_patient_id_type_map)

############################################################
# WRITE TOTAL NUMBER OF MUTATIONS TO CLINICAL FILE
############################################################
#' Update clinical files with total mutation counts
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

clinical_df <- update_clin_with_mut_counts(total_mut_count_df, clinical_df, 
                                           is_tcga = T)

write.csv(clinical_df, paste0(PATH, "Mutation/clinical_data_w_Nonsyn_MutCounts.csv"))  

