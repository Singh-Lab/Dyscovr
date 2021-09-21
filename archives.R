# Create output dataframe to store the results
# Columns: Patients that have all necessary datatypes 1...j...L
# Rows: Downstream gene targets 1...k...M
# Entries: List of Beta values from output of linear model c(B_0, B_2, ... etc.)
# targets <- downstream_target_df[,i]
# output_dataframe <- data.frame(matrix(ncol = length(patients), nrow = length(targets)))

# Loop through all patients 1..j..L
for (j in 1:length(patients)) {
  # TODO: Extract all the patient-specific information (Sex, Age, and TotalNumMut)
  patient_id <- rownames(patients)[j]
  patient <- patients[j,]
  sex_j <- patient$Sex
  age_j <- patient$Age
  total_mut_j <- patient$Total.Num.Mut
  
  # TODO: Extract all the patient and TF gene specific information (MutStat, MethStat, CNAStat)
  mut_stat_ij <- mutation_df[prot_id, patient_id]
  meth_stat_ij <- methylation_df[prot_id, patient_id]
  cna_stat_ij <- cna_df[prot_id, patient_id]
  
  # Loop through all targets 1...k...M
  for (k in length(targets))
    # TODO: Extract the gene target expression in patient j
    exp_kj <- patient[k]    # This makes no sense yet, fix later
  
  
  
  #
}


# From process_ENCODE_data.R:
protein_id_table <- read.csv("proteins_with_dna_binding_doms_and_muts.csv")







protein_ids <- c()
for (i in 1:nrow(protein_id_table)) {
  split_id <- unlist(strsplit(protein_id_table$Query[i],"|", fixed = TRUE))
  protein_ids <- c(protein_ids, unlist(strsplit(split_id[3], "_", fixed = TRUE))[1])
}
protein_id_table_subset <- read.csv("proteins_with_dna_binding_doms_and_muts_in_doms.csv")
protein_ids_subset <- c()
for (i in 1:nrow(protein_id_table_subset)) {
  split_id <- unlist(strsplit(protein_id_table_subset$Query[i],"|", fixed = TRUE))
  protein_ids_subset <- c(protein_ids_subset, unlist(strsplit(split_id[3], "_", fixed = TRUE))[1])
}

# Keep only unique IDs for this purpose
protein_ids <- unique(protein_ids)
protein_ids_subset <- unique(protein_ids_subset)

cnv_clinical_file_sub <- read.csv("cnv_clinical_data_subset.csv", header = TRUE)
patient_id_list <- cnv_clinical_file_sub$case_id

# Convert patient ids to have periods instead of dashes
patient_id_list_periods <- c()
for (i in 1:length(patient_id_list)) {
  split_id <- unlist(strsplit(patient_id_list[i], split = "-", fixed = TRUE))
  new_id <- paste(split_id, collapse = ".")
  patient_id_list_periods <- c(patient_id_list_periods, new_id)
}

# TODO: Convert these IDs to have periods instead of dashes (what GISTIC uses)


# Rsamtools vignette: http://www.bioconductor.org/packages/2.13/bioc/vignettes/Rsamtools/inst/doc/Rsamtools-Overview.pdf
# Use Rsamtools to process raw bam files
# Rtracklayer vignette: http://127.0.0.1:11563/library/rtracklayer/doc/rtracklayer.pdf
# Use rtracklayer to import bed files and BigWig files
# More about the import function: https://kasperdanielhansen.github.io/genbioconductor/html/rtracklayer_Import.html

# NOTE: BED files are the most commonly used and easiest to process



# Get the names of all the expression files and trim the extension from them
exp_filenames_full <- list.files("C:/Users/sarae/Documents/Mona Lab Work/Summer 2020/Main Project Files/Input Data Files/BRCA Data/Expression_Data")
exp_filenames <- c()
for (i in 1:length(exp_filenames_full)) {
  filename <- exp_filenames_full[i]
  trimmed_fn <- sub(".htseq.counts", "", filename)
  exp_filenames <- c(exp_filenames, trimmed_fn)
}
  
  
curated_gene_targets_list <- convert_to_list(curated_gene_targets)
# curated_gene_targets_list_short <- convert_to_list(curated_gene_targets_short)

# Takes a scanned text file containing a list, each element containing its own vector of values
# Converts this file back to its original list format
convert_to_list <- function(targ_scan) {
  # Separate elements by one or more whitespace
  list <- strsplit(targ_scan, "[[:space:]]+")
  # Extract the first vector element and set it as the list element name
  names(list) <- sapply(list, '[[', 1)
  # Remove the first vector element from each list element
  list <- lapply(list, '[', -1)
  return(list)
}


offset_age <- 8     # Number of columns to add in order to get to the start of the age columns
for (i in 1:length(age_bin_vectors)) {
  input_dataframe[,(i+offset)] <- age_bin_vectors[i]
}


# Takes an in-progress input dataframe with a 'TotalNumMut' column we'd like to fill
# in for all the patients (rows) with their total number of mutations. We have these values
# in a given dataframe (total_mut_df) but the IDs need to be properly matched before
# the column is added. This function matches the IDs, inserts the total number of mutations
# for each patient, and returns the updated DF.
add_total_num_mut <- function(input_df, total_mut_df) {
  # For each row in the input dataframe, we want to find the matching patient in the 
  # total mutation dataframe and insert their value
  for (i in 1:length(input_df)) {
    patient_full <- rownames(input_df)[i]
    patient_sub <- unlist(strsplit(patient_full, "-", fixed = TRUE))[3]   # adjust this if we end up only using 4-digit patient ID in main table
    
    # Find the patient in the mutation DF
    total_mut_df_sub <- total_mut_df[rownames(total_mut_df),]
  }
  
}


# Get the mutation count matrix for the given file
maf <- read.maf(maf_filename)
mut_count_matrix <- get_mut_count_matrix(maf)
colnames(mut_count_matrix) <- unlist(lapply(colnames(mut_count_matrix), 
                                            function(x) paste(unlist(strsplit(x, split = "-", fixed = TRUE))[1:3], 
                                                              collapse = "-")))

# Convert these patient IDs to UUIDs
uuids_mut_count_mat <- unlist(lapply(colnames(mut_count_matrix), 
                                     function(x) uuid_barcode_conv[uuid_barcode_conv$tcga_barcode == x, 2]))
colnames(mut_count_matrix) <- uuids_mut_count_mat  #[!is.na(uuids_mut_count_mat)]
print(mut_count_matrix)






############################################################ 
# Takes in three data tables of information:
# 1. patient_df : Patient data frame with information for patients 1..j..M
# 2. prot_df : Protein data frame with information for proteins 1..i..N
# for patient j
# 3. targ_df : Target data frame with information for target genes 
# 1..k..L of protein i in patient j
# Also takes the numerical values of j, i, and k to know which patient & which protein i
# we are considering at the given time
# Creates a table with all the information needed just for these considerations
# and returns it (columns are inputs to linear model)
############################################################ 
create_df_for_lm(patient_df, prot_df, targ_df, j, i) {
  
  # We're going to build on the target data frame by adding rows of repeating elements
  
  # Extract the row for patient j
  patient_j_row <- patient_df[j,]
  
  # Extract the row for protein i
  prot_i_row <- prot_df[i,]
}




if (length(meth_table_sub > 0)) {
  meth_state_array <- c(meth_state_array, 1)
}
else {
  meth_state_array <- c(meth_state_array, 0)  # otherwise, add 0
}

# Create relational dataframe to hold protein-dependent inputs to linear model for this patient
#prot_dataframe <- data.frame(matrix(ncol = 3, nrow = length(protein_ids)))
#colnames(prot_dataframe) <- c("MutStat_i", "MethStat_i", "CNA_i")
#rownames(prot_dataframe) <- protein_ids

# Create a data frame to hold all the gene k-dependent values for TF i and patient j
# targ_dataframe <- data.frame(matrix(ncol = 4, nrow = nrow(length(targets))))
# colnames(targ_dataframe) <- c("Exp_k", "MutStat_k", "MethStat_k", "CNA_k")
# rownames(targ_dataframe) <- targets

# Add the data and column/row names to the dataframe
# methylation_df <- do.call(cbind, meth_array_list)   # Automatically adds list names as column names
# rownames(methylation_df) <- binding_prots

& tot_num_mut <= 30) {
  tot_mut_df[i,"Tot_Mut_b3"] <- 1      # set this bucket for that row equal to 1
} else if (tot_num_mut > 30 & tot_num_mut <= 40) {
  tot_mut_df[i,"Tot_Mut_b4"] <- 1      # set this bucket for that row equal to 1
} else if (tot_num_mut > 40 & tot_num_mut <= 50) {
  tot_mut_df[i,"Tot_Mut_b5"] <- 1      # set this bucket for that row equal to 1
} else if (tot_num_mut > 50 & tot_num_mut <= 60) {
  tot_mut_df[i,"Tot_Mut_b6"] <- 1      # set this bucket for that row equal to 1
} else if (tot_num_mut > 60 & tot_num_mut <= 70) {
  tot_mut_df[i,"Tot_Mut_b7"] <- 1      # set this bucket for that row equal to 1
} else if (tot_num_mut > 70 & tot_num_mut <= 80) {
  tot_mut_df[i,"Tot_Mut_b8"] <- 1      # set this bucket for that row equal to 1
} else if (tot_num_mut > 80 & tot_num_mut <= 90) {
  tot_mut_df[i,"Tot_Mut_b9"] <- 1      # set this bucket for that row equal to 1
} else if (tot_num_mut > 90) {
  tot_mut_df[i,"Tot_Mut_b10"] <- 1      # set this bucket for that row equal to 1
  
  
  
queries <- cdsearch_domains$Query
swissprot_ids <- unlist(sapply(queries, FUN = function(x) {unlist(strsplit(x, "|", fixed = TRUE))[2]}))
matching_indices <- which(canbind_swissprot_ids %fin% swissprot_ids)
cdsearch_domains_canbind_sub <- cdsearch_domains[matching_indices,]
#print(nrow(cdsearch_domains_canbind_sub))


for (i in 1:length(prot_ids)) {
  ensg <- 
    if (!length(ensg) == 0) {
      ensg_ids <- c(ensg_ids, ensg)
    }
}
unlist(lapply(GISTIC_filenames[2:length(GISTIC_filenames)], function(x) rbind(GISTIC_file, read.table(paste(path, paste("CNV_Data/", x, sep = ""), sep = ""), 
                                                                                                      header = TRUE, sep = "\t", colClasses = "character"))))


update_with_mutation_values <- function(input_dataframe, row_index, target_ensg, patient_id, mutation_df) {
  mutation_df_sub <- mutation_df[grep(target_ensg, mutation_df$Gene),]
  if (nrow(mutation_df_sub) == 0) {  # we didn't find a match
    input_dataframe[row_index, 'MutStat_k'] <- 0
  } else {
    mutation_df_sub <- mutation_df_sub[grep(patient_id, mutation_df$Tumor_Sample_Barcode),]
    if (nrow(mutation_df_sub) == 0) {
      input_dataframe[row_index, 'MutStat_k'] <- 0
    } else {
      input_dataframe[row_index, 'MutStat_k'] <- 1
    }
  }
  return(input_dataframe)
}

# Helper function to update a given input dataframe (at a given row index) with 
# expression information for a particular gene (ENSG ID) and patient 
update_with_expression_values <- function(input_dataframe, row_index, target_ensg, patient_id, expression_df) {
  expression_df_sub <- expression_df[rownames(expression_df) == target_ensg, grep(patient_id, colnames(expression_df))] 
  # OPT: subset to just tumor samples
  expression_value <- expression_df_sub[1, grep("01", colnames(expression_df_sub))] 
  input_dataframe[row_index,'Exp_k'] <- expression_value
  return(input_dataframe)
}

# Helper function to update a given input dataframe (at a given row index) with 
# methylation information for a particular gene (ENSG ID) and patient 
update_with_methylation_values <- function(input_dataframe, row_index, target_ensg, patient_id, methylation_df) {
  methyation_df_sub <- methylation_df[rownames(methylation_df) == target_ensg, grep(patient_id, colnames(methylation_df))]
  if(length(methylation_df_sub) == 0) {
    input_dataframe[row_index,'MethStat_k'] <- NA
  } else {
    methylation_value <- methylation_df_sub[1, grep("01", colnames(methylation_df_sub))] 
    input_dataframe[row_index,'MethStat_k'] <- methylation_value
  }
  return(input_dataframe)
}

# Helper function to update a given input dataframe (at a given row index) with 
# CNA information for a particular gene (ENSG ID) and patient 
update_with_cna_values <- function(input_dataframe, row_index, target_ensg, patient_id, cna_df) {
  cna_value <- cna_df[rownames(cna_df) == target_ensg, grep(patient_id, colnames(cna_df))]
  if(length(cna_value) == 0) {
    input_dataframe[row_index,'CNA_k'] <- NA
  } else {
    #cna_value <- cna_value[1, grep("01", colnames(cna_value))]
    input_dataframe[row_index,'CNA_k'] <- cna_value
  }
  return(input_dataframe)
}


gender <- c()
age_b1 <- c()
age_b2 <- c()
age_b3 <- c()
age_b4 <- c()
age_b5 <- c()
age_b6 <- c()
age_b7 <- c()
age_b8 <- c()
age_b9 <- c()
age_b10 <- c()
race_b1 <- c()
race_b2 <- c()
race_b3 <- c()
race_b4 <- c()
race_b5 <- c()
prior_malg <- c()
treatment_rad <- c()
treatment_pharm <- c()
tot_mut_b1 <- c()
tot_mut_b2 <- c()
tot_mut_b3 <- c()

gender <- c(gender, rep(patient_info[1,1], times = num_rep_times))
age_b1 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
age_b2 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
age_b3 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
age_b4 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
age_b5 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
age_b6 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
age_b7 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
age_b8 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
age_b9 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
age_b10 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
race_b1 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
race_b2 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
race_b3 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
race_b4 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
race_b5 <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
prior_malg <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
treatment_rad <- c(age_b1, rep(patient_info[1,2], times = num_rep_times))
treatment_pharm <- c()
tot_mut_b1 <- c()
tot_mut_b2 <- c()
tot_mut_b3 <- c()

# Create lists for all columns to speed runtime (add full list to dataframe at the end)
mutStat_i <- c()
methStat_i <- c()
cna_i <- c()
exp_k <- c()
mutStat_k <- c()
methStat_k <- c()
cna_k <- c()

patient_characteristics_list <- list("Gender" = c(), "Age_b1" = c(), "Age_b2" = c(), "Age_b3" = c(), "Age_b4" = c(),
                                     "Age_b5" = c(), "Age_b6" = c(), "Age_b7" = c(), "Age_b8" = c(), "Age_b9" = c(),
                                     "Age_b10" = c(), "Race_b1" = c(), "Race_b2" = c(), "Race_b3" = c(), "Race_b4" = c(),
                                     "Race_b5" = c(), "Prior_malig" = c(), "Treatment_rad" = c(), "Treatment_pharm" = c(),
                                     "Tot_mut_b1" = c(), "Tot_mut_b2" = c(), "Tot_mut_b3" = c())

# Loop through all patients j (616)
for (j in 1:2) {    #nrow(patient_dataframe)) {
  
  patient_id <- rownames(patient_dataframe)[j]
  
  # Get this patients clinical information
  patient_info <- patient_dataframe[j,]
  
  # Loop through all transcription factors i of interest (protein_ids)
  #for (i in 1:length(protein_ids)) {
  
  # Get the given protein of interest (protein name)
  #prot_id <- protein_ids[i]
  #print(prot_id)
  
  
  # Get the list of this protein's targets
  #targets <- downstream_target_df[,rownames(downstream_target_df) == prot_id]
  #targets_ensembl <- unique(genes$Ensembl.ID)
  
  # Loop through all targets 1...k...M
  for (k in 1:length(targets_ensembl)) {  # this is 60,727 genes if using full set!
    
    # Get the target ENSG ID
    target_ensg <- targets_ensembl[k]
    
    # Get the target protein ID as well
    #target_protid <- ensg_genename_conv[ensg_genename_conv$Gene.stable.ID == target, 2]
    
    # Extract all the target k-specific information and add it to the input DF
    #given_index <- row_index + (k-1)
    
    # EXPRESSION
    #input_dataframe <- get_expression_values(target_ensg, patient_id, expression_df)
    exp_val <- get_expression_values(target_ensg, patient_id, expression_df)
    if (length(exp_val) > 1) {
      print(paste("Exp val length >1:", exp_val))
    }
    exp_k <- c(exp_k, exp_val)
    # MUTATION
    #input_dataframe <- get_mutation_values(target_ensg, patient_id, mutation_df)
    mut_val <- get_mutation_values(target_ensg, patient_id, mutation_df)
    if (length(mut_val) > 1) {
      print(paste("Mut val length >1:", mut_val))
    }
    mutStat_k <- c(mutStat_k, mut_val)
    # METHYLATION
    #input_dataframe <- get_methylation_values(target_ensg, patient_id, methylation_df)
    meth_val <- get_methylation_values(target_ensg, patient_id, methylation_df)
    if (length(meth_val) > 1) {
      print(paste("Meth val length >1:", meth_val))
    }
    methStat_k <- c(methStat_k, meth_val)
    # CNA
    #input_dataframe <- get_cna_values(target_ensg, patient_id, cna_df)
    cna_val <- get_cna_values(target_ensg, patient_id, cna_df)
    if (length(cna_val) > 1) {
      print(paste("CNA val length >1:", cna_val))
    }
    cna_k <- c(cna_k, cna_val)
  }
  # print(input_dataframe)
  # Now we need to add the protein i information next to every target!
  
  # Get the ENSG version of the protein ID as well
  ensg_id <- uniprot_ensg_conv[uniprot_ensg_conv$From == prot_id, 'To']
  print(ensg_id)
  
  # Extract all the TF gene specific information (MutStat, MethStat, CNAStat)
  #colnames(mutation_df) <- unlist(lapply(colnames(mutation_df), function(x) unlist(strsplit(x, ".", fixed = TRUE))[3]))
  #mut_df_sub <- mutation_df[grep(ensg_id, mutation_df$ensg_ids),]
  #mut_df_sub <- mut_df_sub[grep(patient_id, mutation_df$Tumor_Sample_Barcode),]
  #mut_stat <- ifelse(nrow(mut_df_sub) == 0, 0, 1)
  mut_stat <- get_mutation_values(ensg_id, patient_id, mutation_df)
  #meth_stat <- methylation_df[methylation_df$ensg_ids == ensg_id, grep(patient_id, colnames(mutation_df))]
  meth_stat <- get_methylation_values(ensg_id, patient_id, methylation_df)
  #cna_stat <- cna_df[rownames(cna_df) == ensg_id, grep(patient_id, colnames(mutation_df))]
  cna_stat <- get_cna_values(ensg_id, patient_id, cna_df)
  
  # Add these, along with patient specific information, in next to every target value (need to repeat them)
  mut_stat <- rep(mut_stat, times = length(targets_ensembl))
  mutStat_i <- c(mutStat_i, mut_stat)
  meth_stat <- rep(meth_stat, times = length(targets_ensembl))
  methStat_i <- c(methStat_i, meth_stat)
  cna_stat <- rep(cna_stat, times = length(targets_ensembl))
  cna_i <- c(cna_i, cna_stat)
  #for (rep in 1:length(targets_ensembl)) {
  # Protein inputs 
  #input_dataframe[(rep + row_index),'MutStat_i'] <- mut_stat
  #input_dataframe[(rep + row_index),'MethStat_i'] <- meth_stat
  #input_dataframe[(rep + row_index),'CNA_i'] <- cna_stat
  
  # Patient inputs
  #do.call(function(x) 
  #input_dataframe[rep + row_index, colnames(patient_dataframe)[x]] <- patient_info[x], 
  #args = list(x = 1:length(patient_info)))
  #}
  #print(input_dataframe)
  # Update the row index to the end of the filled section
  #row_index <- row_index + length(targets)
  #}
  
  # This is the number of data points to fill in with repeats of the patient information
  # We want to add to the existing patient characteristics list with the current patient's info, but we
  # need to repeat it enough times to match the length of protein i and gene k (so we get an even table)
  num_rep_times <- length(mut_stat)
  
  #patient_characteristics_list <- lapply(1:length(names(patient_characteristics_list)), 
  #function(x) c(patient_characteristics_list[[x]], rep(patient_info[,x], times = num_rep_times)))
  for (i in 1:length(names(patient_characteristics_list))) {
    patient_characteristics_list[[i]] <- c(patient_characteristics_list[[i]], rep(patient_info[,i], times = num_rep_times))
  }
  print(patient_characteristics_list)
}



patient_characteristics_df <- as.data.frame(do.call(cbind, patient_characteristics_list))
input_dataframe <- cbind(input_dataframe, patient_characteristics_df)    #, by = intersect(names(input_dataframe), names(patients_characteristics_df)))



# TODO: get some statistics on the numbers of NA values (sparsity of this data); potentially eliminate genes that have too many NA values?
# Also check how these models handle NA values
# TODO: run through this in a loop - we want to run an individual linear model for each gene-target pairing across all patients

# Now, we have all the info we need to run our simple linear model function
lm_fit <- NA
if (!is_brca) {
  # Include gender as a covariate
  lm_fit <- lm(Exp_k ~ MutStat_i + MethStat_i + MutStat_k + MethStat_k + CNA_k + CNA_i + Gender +
                 Age_b1 + Age_b2 + Age_b3 + Age_b4 + Age_b5 + Age_b6 + Age_b7 + Age_b8 + 
                 Age_b9 + Age_b10 + Race_b1 + Race_b2 + Race_b3 + Race_b4 + Race_b5 + Prior_malig + 
                 Treatment_rad + Treatment_pharm + Tot_mut_b1 + Tot_mut_b2 + Tot_mut_b3, 
               data = input_dataframe)
} else {
  # If this is breast cancer, we are limiting to only women, so do not include gender as 
  # a covariate
  lm_fit <- lm(formula = (Exp_k ~ MutStat_i + MethStat_i + MutStat_k + MethStat_k + CNA_k + CNA_i + 
                            Age_b1 + Age_b2 + Age_b3 + Age_b4 + Age_b5 + Age_b6 + Age_b7 + Age_b8 + 
                            Age_b9 + Age_b10 + Race_b1 + Race_b2 + Race_b3 + Race_b4 + Race_b5 + Prior_malig + 
                            Treatment_rad + Treatment_pharm + Tot_mut_b1 + Tot_mut_b2 + Tot_mut_b3), 
               data = input_dataframe_no_gen)
  # NOTES: R is dropping some of the variables it predicts to be linearly dependent
  # Some of these categories are unnecessary (can be inferred from the others - read about the dummy variable trap
  # here: https://www.algosome.com/articles/dummy-variable-trap-regression.html)
  
  
}

# Exploratory analysis -- Scatterplots
plot(mutStat_i, exp_k_log, main = "Scatterplot of NF1 mutation status and target gene expression, across BRCA patients")
plot(cna_i, exp_k_log, main = "Scatterplot of NF1 CNA status and target gene expression, across BRCA patients")

# Optional: print a summary of the raw LM results
print(summary(lm_fit))

# Use broom package to put results in a more appealing format
tidy_lmfit <- tidy(lm_fit)
pval <- as.numeric(tidy_lmfit$p.value[3])
std_error <- as.numeric(tidy_lmfit$std.error[2])
tidy_lmfit$qvalues <- qvalue(p = tidy_lmfit$p.value)

# Extract the beta coefficients
coef_lm <- lm_fit$coefficients




if (nrow(meth_df_sub) == 0) {
  print("Methylation: NA")
  return(NA)
} else {
  
  
  
}


### OPTION 2: USE INTERACDOME PERTININT TRACK ###
# NOTE: THIS OPTION APPEARS TO PRODUCE IDENTICAL RESULTS TO THE ABOVE
# interacdome_track <- read.table(paste(path, "Input Data Files/interacdome0.3-pfam31_domainweights-GRCh38.txt", sep = ""), 
# header = TRUE, sep = "\t", check.names = FALSE)
# Remove all values below threshold
# interacdome_track <- interacdome_track[interacdome_track$binding_frequency > threshold_interacdome,]

# Get remaining binding domains
# binding_domains <- unique(interacdome_track$domain_name)
# binding_domains <- binding_domains[!is.na(binding_domains)]  # Remove any NA values
# NOTE: all the same domains as above, just a different number of binding positions
# binding_domains_ids <- grep('PF', unlist(strsplit(binding_domains, "_", fixed = TRUE)), value = TRUE)
# binding_domains_ids_noPF <- unique(unlist(lapply(binding_domains_ids, FUN = function(id) str_remove(id, "PF"))))

# Make a column of stripped domains 
# interacdome_track$stripped_domain_name <- unlist(lapply(interacdome_track$domain_name, function(x) unlist(strsplit(unlist(strsplit(x, "_", fixed = TRUE))[1], "F", fixed = TRUE))[2]))

# binding_positions_df <- data.frame(matrix(ncol = 1, nrow = length(binding_domains_ids_noPF)))
# colnames(binding_positions_df) <- c("Binding.Pos")
# for (i in 1:length(binding_domains_ids_noPF)) {
# Get the binding frequencies for this domain
# domain <- binding_domains_ids_noPF[i]

# Subset interacdome DF to only look at this domain
# interac_sub <- interacdome_track[interacdome_track$stripped_domain_name == domain,]

# Extract the binding positions
# binding_pos <- interac_sub$"1-indexed_match_state"

# Often, these positions are in the same domain for different ligands. If we have already restricted
# our domains to just DNA-binding, the InteracDome DF will already be subsetted for this (no additional action required).
# unique_binding_pos <- unique(binding_pos)

# Add these binding positions to the dataframe
# binding_positions_df$Binding.Pos[i] <- paste(unique_binding_pos, collapse = ",")
# }
# rownames(binding_positions_df) <- binding_domains_ids_noPF
# write.csv(binding_positions_df, paste(path, paste(paste("Saved Output Data Files/InteracDome/binding_positions_DF_track_", as.character(threshold_interacdome), sep = ""), ".csv", sep = ""), sep = ""))


# The following works as well, but is slow with the for-loop; use helper function instead
#swissprot_ids_col <- list()
#for (i in 1:nrow(concavity_df_sub)) {
#  ensg <- unlist(strsplit(concavity_df_sub[i,1], ".", fixed = TRUE))[1]
#  ensg_swissprot_conv_sub <- ensg_swissprot_conv[grepl(ensg, ensg_swissprot_conv$ENSP),]
#  swissprot_ids_col[[i]] <- ensg_swissprot_conv_sub$Entry
#  print(ensg_swissprot_conv_sub$Entry)
#}


get_ibindingpos_subset <- function(swissprot_ids, interacdome_binding_positions_df, maf_df, canbind_df, domain_df, concavity_df) {
  # Create a DF to hold the final results
  ibindingpos_df <- data.frame(matrix(nrow = 0, ncol = 4))
  
  # Fix the query column of the domain DF & add a stripped accessions column
  domain_df$Query <- unlist(lapply(domain_df$Query, function(x) unlist(strsplit(x, "|", fixed = TRUE))[2]))
  domain_df <- domain_df[grepl("pfam", domain_df$Accession),]  # limit to just pfam IDs
  domain_df$Stripped.Accessions <- as.character(regmatches(domain_df$Accession, regexpr("[[:digit:]]+",
                                                                                        domain_df$Accession)))  # get just the numeric portion
  
  # Loop through all proteins 
  for (i in 1:length(swissprot_ids)) {
    id <- swissprot_ids[i]  # Swissprot ID of current protein
    print(id)
    print(paste(i, paste("/", length(swissprot_ids))))
    
    # Subset the MAF file and domain DF to only look at this protein 
    maf_sub <- maf_df[maf_df$SWISSPROT == id,]    
    domain_sub <- domain_df[domain_df$Query == id,]
    
    if (!nrow(domain_sub) == 0) {
      # Use the domain sub to also subset the InteracDome binding positions table (only keep the domains
      # that this particular protein possesses)
      interacdome_binding_positions_df[,1] <- as.character(interacdome_binding_positions_df[,1])
      interacdome_sub <- interacdome_binding_positions_df[unlist(lapply(interacdome_binding_positions_df[,1],
                                                                        function(x) x %fin% domain_sub$Stripped.Accessions)),]
      #print(interacdome_sub)
      
      # Binding positions in InteracDome are currently indexed by the domain, not the full
      # protein length (which makes them uncomparable) -- need to add the "start" position of the domain to 
      # each binding position & consolidate them all into one list
      interacDome_binding_pos <- get_interac_bp(domain_sub, interacdome_sub)
      
      if(!missing(concavity_df)) {
        all_binding_pos <- c(canbind_df[grep(id, canbind_df$Swissprot, ignore.case = FALSE, fixed = TRUE), 'AA_position'], 
                             concavity_df[grep(id, concavity_df$Swissprot, ignore.case = FALSE, fixed = TRUE),'AA_position'], 
                             interacDome_binding_pos)
      } else {
        all_binding_pos <- c(canbind_df[grep(id, canbind_df$Swissprot, ignore.case = FALSE, fixed = TRUE), 'AA_position'], 
                             interacDome_binding_pos)
      }
      
      all_binding_pos <- unique(all_binding_pos)
      
      for (j in 1:nrow(maf_sub)) {
        
        # For each mutation in this protein, get its AA position ("Protein_position", column 55)
        mut_pos <- as.numeric(unlist(strsplit(maf_sub$Protein_position[j], "/", fixed = TRUE))[1])
        print(mut_pos)
        
        # Also get the patient ID
        patient_id <- unlist(strsplit(maf_sub$Tumor_Sample_Barcode[j], "-", fixed = TRUE))[3]
        
        if (!mut_pos == "-") {
          # Is this mutation in a CanBind, ConCavity, or InteracDome position?
          if (mut_pos %fin% all_binding_pos) {
            # Yes! Add this entry, with the mutation position added as well
            row_to_add <- c(id, paste(all_binding_pos, collapse = ","), mut_pos, patient_id)
            ibindingpos_df <- rbind(ibindingpos_df, row_to_add)
          }
        }
      }
    }
  }
  colnames(ibindingpos_df) <- c("Protein.ID", "Pred.Binding.Pos", "Mut.Pos", "Patient")
  return(ibindingpos_df)
}

# Helper function that takes all the binding positions in each domain and, for a given protein,
# adds the start position of the domain to the domain binding position to get the overall binding
# position in the protein. Collapses them all into one vector
get_interac_bp <- function(domain_sub, interacdome_sub) {
  all_binding_pos <- unique(unlist(lapply(interacdome_sub[,1], function(dom) {
    # Get the binding positions for this domain from InteracDome
    print(paste("dom:",dom))
    bps <- interacdome_sub$Binding.Pos[which(interacdome_sub[,1] == dom)]
    print(paste("bps:",bps))
    
    # Subset domain DF to just the domain in question 
    subsetted_dom_df <- domain_sub[domain_sub$Stripped.Accessions == dom,]
    
    # Get the start position(s) for this domain
    start_pos <- subsetted_dom_df$From
    
    if (length(start_pos) > 1) {   # If there is more than one start position of this domain, add positions for all domains
      return(unlist(lapply(start_pos, function(x) unlist(lapply(unlist(strsplit(bps, split = ",", fixed = TRUE)), 
                                                                function(y) as.numeric(y) + as.numeric(x))))))
    } else {    # Otherwise, just add the single start position to all the binding positions in this domain
      return(unlist(lapply(unlist(strsplit(bps, split = ",", fixed = TRUE)), function(x) as.numeric(x) + as.numeric(start_pos))))
    }
    
  })))
  all_binding_pos <- all_binding_pos[!is.na(all_binding_pos)]
  print(all_binding_pos)
  return(all_binding_pos)
}




# Helper function that takes all the binding positions in each domain and, for a given protein,
# adds the start position of the domain to the domain binding position to get the overall binding
# position in the protein. Collapses them all into one vector
get_interac_bp <- function(domain_sub, interacdome_sub) {
  all_binding_pos <- unique(unlist(lapply(interacdome_sub[,1], function(dom) {
    # Get the binding positions for this domain from InteracDome
    print(paste("dom:",dom))
    bps <- interacdome_sub$Binding.Pos[which(interacdome_sub[,1] == dom)]
    print(paste("bps:",bps))
    
    # Subset domain DF to just the domain in question 
    subsetted_dom_df <- domain_sub[domain_sub$Stripped.Accessions == dom,]
    
    # Get the start position(s) for this domain
    start_pos <- subsetted_dom_df$From
    
    if (length(start_pos) > 1) {   # If there is more than one start position of this domain, add positions for all domains
      return(unlist(lapply(start_pos, function(x) unlist(lapply(unlist(strsplit(bps, split = ",", fixed = TRUE)), 
                                                                function(y) as.numeric(y) + as.numeric(x))))))
    } else {    # Otherwise, just add the single start position to all the binding positions in this domain
      return(unlist(lapply(unlist(strsplit(bps, split = ",", fixed = TRUE)), function(x) as.numeric(x) + as.numeric(start_pos))))
    }
    
  })))
  all_binding_pos <- all_binding_pos[!is.na(all_binding_pos)]
  print(all_binding_pos)
  return(all_binding_pos)
}


# From the I-Binding Pos function
# Fix the query column of the domain DF & add a stripped accessions column
domain_df$Query <- unlist(lapply(domain_df$Query, function(x) unlist(strsplit(x, "|", fixed = TRUE))[2]))
domain_df <- domain_df[grepl("pfam", domain_df$Accession),]  # limit to just pfam IDs
domain_df$Stripped.Accessions <- as.character(regmatches(domain_df$Accession, regexpr("[[:digit:]]+",
                                                                                      domain_df$Accession)))  # get just the numeric portion


interacdome_binding_positions_df[,1] <- as.character(interacdome_binding_positions_df[,1])
interacdome_sub <- interacdome_binding_positions_df[unlist(lapply(interacdome_binding_positions_df[,1],
                                                                  function(x) x %fin% domain_sub$Stripped.Accessions)),]

domains <- unique(domain_sub$Stripped.Accessions)
print(paste("domains", domains))
# Also subset the interacdome file to only look at this protein
interacdome_sub <- interacdome_binding_positions_df[grep(id, interacdome_binding_positions_df$Swissprot, 
                                                         ignore.case = TRUE, fixed = TRUE),]
print(interacdome_sub)
interacdome_sub$Stripped.Accessions <- as.character(regmatches(interacdome_sub$Pfam_HMM_ID, regexpr("[[:digit:]]+",
                                                                                                    interacdome_sub$Pfam_HMM_ID)))
interacdome_sub <- interacdome_sub[interacdome_sub$Stripped.Accessions %fin% domains,]
interacDome_binding_pos <- unique(unlist(lapply(interacdome_sub[,'Binding.Pos'], function(x) unlist(strsplit(x, split = ",", fixed = TRUE)))))
print(interacDome_binding_pos)

# Binding positions in InteracDome are currently indexed by the domain, not the full
# protein length (which makes them uncomparable) -- need to add the "start" position of the domain to 
# each binding position & consolidate them all into one list
#interacDome_binding_pos <- get_interac_bp(domain_sub, interacdome_sub)

interacdome_df_sub <- interacdome_df[grep(id, interacdome_df$Swissprot),]
#print(interacdome_df_sub)
canbind_df_sub <- canbind_df[grep(id, canbind_df$Swissprot),]
if(!missing(concavity_df)) {
  concavity_df_sub <- concavity_df[grep(id, concavity_df$Swissprot),]
}

if(!missing(concavity_df)) {
  all_binding_pos <- c(canbind_df[grep(id, canbind_df$Swissprot, ignore.case = FALSE, fixed = TRUE), 'AA_position'], 
                       concavity_df[grep(id, concavity_df$Swissprot, ignore.case = FALSE, fixed = TRUE),'AA_position'], 
                       interacDome_binding_pos)
} else {
  all_binding_pos <- c(canbind_df[grep(id, canbind_df$Swissprot, ignore.case = FALSE, fixed = TRUE), 'AA_position'], 
                       interacDome_binding_pos)
}
all_binding_pos <- sort(unique(all_binding_pos), decreasing = FALSE)


# Extract those unique proteins as a list for study
# Function takes in the newly created dataframe containing the proteins with mutations
# in binding domains
get_unique_prots <- function(domains_missense_idomain_sub, swiss_or_nam) {
  protein_ids_subset <- c()
  for (i in 1:nrow(domains_missense_idomain_sub)) {
    split_id <- unlist(strsplit(domains_missense_idomain_sub$Query[i],"|", fixed = TRUE))
    if (swiss_or_nam == "nam") {
      protein_ids_subset <- c(protein_ids_subset, unlist(strsplit(split_id[3], "_", fixed = TRUE))[1])
    } else {
      protein_ids_subset <- c(protein_ids_subset, unlist(strsplit(split_id[2], "_", fixed = TRUE))[1])
    }
  }
  return(unique(protein_ids_subset))
}

freq_df$Freq <- unlist(lapply(freq_df$Num.Patients, function(x) (x/(length(unique(maf_file_df$Tumor_Sample_Barcode))))*100))
#print(freq_df)

combined_vals_array <- unlist(lapply(1:length(diff_beta_array), function(k) {
  
  # Get the beta value entries
  diff_beta_entry <- diff_beta_array[k]
  diff_beta_bin_entry <- diff_beta_array_binary[k]
  
  # Get the probe number entries
  diff_probe_entry <- diff_probe_array_plus_binary[k]
  diff_probe_abs_entry <- diff_probe_entry[1]
  diff_probe_bin_entry <- diff_probe_entry[2]
  
  # Combine these using semicolons
  
  return(comb_entry)
}))


# Helper function to remove NA rows and rows without gene names/ IDs
# from the unprocessed methylation table
fix_initial_table <- function(meth_table, gene) {
  # Remove NA rows from the table as well as rows that do not have a gene name or ENSG ID; 
  # NOT NECESSARY if we're just subsetting to the gene 
  #meth_table <- meth_table[!(is.na(meth_table$Beta_value) |
  #(meth_table$Gene_Symbol == ".") | 
  #(meth_table$Transcript_ID == ".")),]
  #meth_table <- meth_table %>% filter(!(is.na(Beta_value) | Gene_Symbol == "." | Transcript_ID == "."))
  
  # Subset the table to only look at the given gene
  
  
  return(meth_table)
}


# OPTION 1: USE BIOMART CONVERSION FILE
# biomart_ensg_name_conv <- read.table("ensg_genename_conv.txt", header = TRUE, sep = ",")

# Function takes a bioMart ENSG-Gene Name conversion file, along with a list of
# gene names, and returns those gene name ENSG IDs
#get_ensg_ids <- function(conv_table, gene_names) {
#  ensg_ids <- c()
#  for (i in 1:length(gene_names)) {
#    gene_name <- gene_names[i]
#    if (gene_name %in% conv_table[,'Gene.name']) {
#      index <- which(gene_name %in% conv_table[,'Gene.name'])
#      ensg <- conv_table[index,'Gene.stable.ID']
#      ensg_ids <- c(ensg_ids, ensg)
#    }
#    else {
#      print(paste("Not found in conversion document: ", gene_name))
#    }
#  }
#  return(ensg_ids)
#}


# FROM CNV FILE
# Use helper function to get ENSG IDs
#protein_ids_ensg <- get_ensg_ids(biomart_ensg_name_conv, protein_ids)
#protein_ids_missense_ensg <- get_ensg_ids(biomart_ensg_name_conv, protein_ids_missense)
#protein_ids_subset_ensg <- get_ensg_ids(biomart_ensg_name_conv, protein_ids_subset)
#protein_ids_subset_missense_ensg <- get_ensg_ids(biomart_ensg_name_conv, protein_ids_subset_missense)

# OPTION 2: Use Uniprot's web interface mapping tool
# Link: https://www.uniprot.org/mapping; input swissprot_ids_full
#protein_ids_conv_table <- read.table(paste(path, "ID Conversions/uniprot_to_ensg_conv.txt", sep = ""), header = TRUE, sep = "\t")
# ensg_ids <- protein_ids_conv_table[,2]

#protein_ids_iprotein <- read.csv(paste(output_path, "Mutation/swissprot_ids_missense_iprotein.csv", sep = ""), header = TRUE)[,2]
#protein_ids_idomain <- read.csv(paste(output_path, "Mutation/swissprot_ids_missense_idomain.csv", sep = ""), header = TRUE)[,2]
#protein_ids_ibindingpos <- read.csv(paste(output_path, "Mutation/swissprot_ids_missense_ibindingpos.csv", sep = ""), header = TRUE)[,2]

# Function to convert the Swissprot protein IDs to ENSG IDs
#get_ensg_ids_from_prot_ids <- function(prot_ids, conv_table) {
#  ensg_ids <- unlist(lapply(prot_ids, function(x) ifelse(!length(conv_table[conv_table$From == x, 2]) == 0, 
#                                                         conv_table[conv_table$From == x, 2], NA)))
#  return(ensg_ids[!is.na(ensg_ids)])
#}

#ensg_ids_iprotein <- get_ensg_ids_from_prot_ids(protein_ids_iprotein, protein_ids_conv_table)
#ensg_ids_idomain <- get_ensg_ids_from_prot_ids(protein_ids_idomain, protein_ids_conv_table)
#ensg_ids_ibindingpos <- get_ensg_ids_from_prot_ids(protein_ids_ibindingpos, protein_ids_conv_table)


# FROM METHYLATION FILE

get_average_beta_per_gene_df <- function(filenames, path, matched_patient_IDs, gene_names, label, binding_prots) {
  if (label == "subset") {
    methylation_rows <- get_methylation_rows(filesnames, path, matched_patient_IDs, binding_prots, "beta")
  } else if (label == "full") {
    methylation_rows <- get_methylation_rows(filenames, path, matched_patient_IDs, gene_names, "beta")
  }
  print(methylation_rows)
  # Add the methylation arrays for this patient to the DF
  methylation_df <- as.data.frame(data.table::rbindlist(methylation_rows, use.names = TRUE, fill = TRUE))
  
  # Add patient IDs as column headers
  #colnames(methylation_df) <- matched_patient_IDs
  
  # Add the gene IDs as row names
  if (label == "subset") {
    rownames(methylation_df) <- binding_prots
  } else {
    #rownames(methylation_df) <- gene_names
    rownames(methylation_df) <- gene_names[1:2]
  } 
  print(methylation_df)
  
  
  # Return finished dataframe
  return(methylation_df)
}


# Call this function to get the differential methylation table
avg_beta_per_gene_df <- get_average_beta_per_gene_df(filenames = filenames, 
                                                     path = path, 
                                                     matched_patient_IDs = matched_patient_IDs, 
                                                     gene_names = gene_names,
                                                     label = "full")

# Write to a file
write.csv(avg_beta_per_gene_df, paste(path, "Methylation/Differential Methylation/average_beta_per_gene.csv", sep = ""))

get_absolute_num_markers_per_gene_df <- function(filenames, path, matched_patient_IDs, gene_names, label, binding_prots, conf_thres) {
  if (label == "subset") {
    methylation_rows <- get_methylation_rows(filesnames, path, matched_patient_IDs, binding_prots, "num_probes", conf_thres)
  } else if (label == "full") {
    methylation_rows <- get_methylation_rows(filenames, path, matched_patient_IDs, gene_names, "num_probes", conf_thres)
  }
  print(methylation_rows)
  # Add the methylation arrays for this patient to the DF
  methylation_df <- as.data.frame(data.table::rbindlist(methylation_rows, use.names = TRUE, fill = TRUE))
  
  # Add patient IDs as column headers
  #colnames(methylation_df) <- matched_patient_IDs
  
  # Add the gene IDs as row names
  if (label == "subset") {
    rownames(methylation_df) <- binding_prots
  } else {
    rownames(methylation_df) <- gene_names
    #rownames(methylation_df) <- gene_names[1:2]
  } 
  print(methylation_df)
  
  
  # Return finished dataframe
  return(methylation_df)
}

conf_thres <- 0.8

avg_beta_per_gene_df <- get_average_beta_per_gene_df(filenames = filenames, 
                                                     path = path, 
                                                     matched_patient_IDs = matched_patient_IDs, 
                                                     gene_names = gene_names,
                                                     label = "full")

# Write to a file
write.csv(avg_beta_per_gene_df, paste(path, "Methylation/Differential Methylation/average_beta_per_gene.csv", sep = ""))


get_methylation_rows <- function(filesnames, path, intersecting_ids, gene_names, identifier, thres) {
  
  methylation_rows <- mclapply(1:2, function(i) {    #length(gene_names), function(i) {
    gene <- gene_names[i]
    
    print(paste(i, paste("/", length(gene_names))))
    
    tic("get meth vals all patients")
    
    # Get the methylation values for each patient
    meth_vals_all_patients <- lapply(1:length(intersecting_ids), function(j) {
      
      patient_id <- intersecting_ids[j]
      tumor_marker <- "01A"
      normal_marker <- "11A"
      filenames_patient <- filenames[grepl(patient_id, filenames)]
      filename_tumor <- filenames_patient[grepl("01A", filenames_patient)]
      if (length(filename_tumor) == 0) {
        filename_tumor <- filenames_patient[grepl("01B", filenames_patient)]
        tumor_marker <- "01B"
      }
      filename_normal <- filenames_patient[grepl("11A", filenames_patient)]
      if(length(filename_normal) == 0) {
        filename_normal <- filenames_patient[grepl("11B", filenames_patient)]
        normal_marker <- "11B"
      }
      
      # Extract the methylation tables from these files
      tic("upload tables")
      meth_table_tumor <- data.table::fread(paste(path, paste("Methylation_Data/", filename_tumor, sep = ""),
                                                  sep = ""), sep = "\t", header = TRUE, select = c(Beta_value='numeric', Gene_Symbol='character'))
      meth_table_normal <- data.table::fread(paste(path, paste("Methylation_Data/", filename_normal, sep = ""),
                                                   sep = ""), sep = "\t", header = TRUE, select = c(Beta_value='numeric', Gene_Symbol='character'))
      toc()
      
      # Subset to only the given gene and make sure Beta value isn't NA
      tic("filter tables")
      #meth_table_tumor <- meth_table_tumor[grepl(gene, meth_table_tumor$Gene_Symbol) & !(is.na(meth_table_tumor$Beta_value)),]
      #meth_table_normal <- meth_table_normal[grepl(gene, meth_table_normal$Gene_Symbol) & !(is.na(meth_table_normal$Beta_value)),]
      meth_table_tumor <- meth_table_tumor %>% filter(grepl(gene, meth_table_tumor$Gene_Symbol) & !(is.na(Beta_value)))
      meth_table_normal <- meth_table_normal %>% filter(grepl(gene, meth_table_normal$Gene_Symbol) & !(is.na(Beta_value)))
      # meth_table_tumor <- meth_table_tumor %>% filter(stri_detect_fixed(gene, meth_table_tumor$Gene_Symbol) & !(is.na(Beta_value)))
      # meth_table_normal <- meth_table_normal %>% filter(stri_detect_fixed(gene, meth_table_normal$Gene_Symbol) & !(is.na(Beta_value)))
      toc()
      
      # Extract the differential average beta value for this gene (tumor-normal)
      tic("get values")
      
      if(identifier == "beta") {
        beta_normal <-  mean.default(meth_table_tumor$Beta_value)
        beta_tumor <-  mean.default(meth_table_normal$Beta_value)
        results <- cbind(beta_normal, beta_tumor)
      }
      else if (identifier == "num_probes") {
        abs_normal <- get_num_probes(meth_table_normal, thres)
        abs_tumor <- get_num_probes(meth_table_tumor, thres)
        results <- cbind(abs_normal, abs_tumor)
      }
      else {print("Invalid function.")}
      
      toc()
      
      # Make the results into a mini dataframe
      entry <- as.data.frame(results)
      print(entry)
      colnames(entry) <- c(paste(patient_id, tumor_marker, sep = "_"), paste(patient_id, normal_marker, sep = "_"))
      # Combine all these values into a semicolon separated list 
      
      # Return this value for this gene in this patient
      print(entry)
      return(entry)
    })
    toc()
    # Combine these all into a DF
    tic("combine cols and convert to DF")
    meth_vals_all_patients_df <- as.data.frame(list.cbind(meth_vals_all_patients))  #do.call(cbind, meth_vals_all_patients)
    toc()
    # Return the row for this gene
    print(meth_vals_all_patients_df)
    
    return(meth_vals_all_patients_df)
  })
  return(methylation_rows)
}


# Helper function that, given a tumor and normal methylation table, as well as a 
# gene ID, gets the number of fully methylated probes for each gene and 
# subtracts tumor from normal; returns this difference
get_num_probes <- function(meth_table, beta_thres) {
  
  # Subset table to include only fully methylated probes
  meth_table_sub <- meth_table[meth_table$Beta_value > beta_thres,]
  
  # Get & return the number of probes
  return(nrow(distinct(meth_table_sub)))
}


# TODO: OPTIMIZE THIS FUNCTION TO MAKE IT FASTER (CURRENTLY MUCH TOO SLOW); takes 238.75 seconds for all patients, 2 genes
get_differential_methylation_df <- function(filenames, path, matched_patient_IDs, gene_names, label, binding_prots) {
  upper_beta_threshold <- 0.8   # fully methylated
  
  if (label == "subset") {
    methylation_rows <- get_methylation_rows(filesnames, path, matched_patient_IDs, binding_prots, upper_beta_threshold)
  } else if (label == "full") {
    methylation_rows <- get_methylation_rows(filenames, path, matched_patient_IDs, gene_names, upper_beta_threshold)
  }
  print(methylation_rows)
  # Add the methylation arrays for this patient to the DF
  methylation_df <- as.data.frame(data.table::rbindlist(methylation_rows, use.names = TRUE, fill = TRUE))
  
  # Add patient IDs as column headers
  colnames(methylation_df) <- matched_patient_IDs
  
  # Add the gene IDs as row names
  if (label == "subset") {
    rownames(methylation_df) <- binding_prots
  } else {
    #rownames(methylation_df) <- gene_names
    rownames(methylation_df) <- gene_names[1:2]
  } 
  print(methylation_df)
  
  
  # Return finished dataframe
  return(methylation_df)
}

# Helper function that, given the filenames, the path, the intersecting patient IDs, and 
# a list of gene IDs of interest, iterates through all the genes and creates a list of
# rows consisting of the methylation values of interest for each gene, each entry 
# representing a patient
get_methylation_rows <- function(filesnames, path, intersecting_ids, gene_names, upper_beta_threshold) {
  
  methylation_rows <- mclapply(1:2, function(i) {     #length(gene_names), function(i) {
    gene <- gene_names[i]
    
    print(paste(i, paste("/", length(gene_names))))
    
    tic("get meth vals all patients")
    
    # Get the methylation values for each patient
    meth_vals_all_patients <- lapply(1:length(intersecting_ids), function(j) {
      
      patient_id <- intersecting_ids[j]
      filenames_patient <- filenames[grepl(patient_id, filenames)]
      filename_tumor <- filenames_patient[grepl("01A", filenames_patient)]
      if (length(filename_tumor) == 0) {
        filename_tumor <- filenames_patient[grepl("01B", filenames_patient)]
      }
      filename_normal <- filenames_patient[grepl("11A", filenames_patient)]
      if(length(filename_normal) == 0) {
        filename_normal <- filenames_patient[grepl("11B", filenames_patient)]
      }
      
      print(filename_tumor)
      print(filename_normal)
      
      # Extract the methylation tables from these files
      tic("upload tables")
      meth_table_tumor <- data.table::fread(paste(path, paste("Methylation_Data/", filename_tumor, sep = ""),
                                                  sep = ""), sep = "\t", header = TRUE, select = c(Beta_value='numeric', Gene_Symbol='character'))
      meth_table_normal <- data.table::fread(paste(path, paste("Methylation_Data/", filename_normal, sep = ""),
                                                   sep = ""), sep = "\t", header = TRUE, select = c(Beta_value='numeric', Gene_Symbol='character'))
      #meth_table_tumor <- vroom(paste(path, paste("Methylation_Data/", filename_tumor, sep = ""), sep = ""), 
      #delim = "\t", col_names = TRUE, col_select = c('Beta_value', 'Gene_Symbol'))
      #meth_table_normal <- vroom(paste(path, paste("Methylation_Data/", filename_normal, sep = ""), sep = ""), 
      #delim = "\t", col_names = TRUE, col_select = c('Beta_value', 'Gene_Symbol'))
      toc()
      
      # Subset to only the given gene and make sure Beta value isn't NA
      tic("filter tables")
      #meth_table_tumor <- meth_table_tumor[grepl(gene, meth_table_tumor$Gene_Symbol) & !(is.na(meth_table_tumor$Beta_value)),]
      #meth_table_normal <- meth_table_normal[grepl(gene, meth_table_normal$Gene_Symbol) & !(is.na(meth_table_normal$Beta_value)),]
      meth_table_tumor <- meth_table_tumor %>% filter(grepl(gene, meth_table_tumor$Gene_Symbol) & !(is.na(Beta_value)))
      meth_table_normal <- meth_table_normal %>% filter(grepl(gene, meth_table_normal$Gene_Symbol) & !(is.na(Beta_value)))
      toc()
      
      # Extract the differential average beta value for this gene (tumor-normal)
      tic("get values")
      diff_beta <- extract_avg_beta(meth_table_tumor, meth_table_normal)
      print(diff_beta)
      
      # If the difference exceeds the upper threshold or is less then the lower threshold, 1 (otherwise 0)
      diff_beta_binary <- 0
      if((diff_beta > upper_beta_threshold) || (diff_beta < (-upper_beta_threshold))) {diff_beta_binary <- 1}
      
      # Get the absolute difference in the number of confident methylation probes between tumor and 
      # normal for this gene, as well as whether normal is 0 and tumor is nonzero (or vice versa)
      diff_probe_plus_binary <- get_num_probes(meth_table_tumor, meth_table_normal, upper_beta_threshold)
      diff_probe_abs <- diff_probe_plus_binary[1]
      diff_probe_bin <- diff_probe_plus_binary[2]
      toc()
      
      # Combine all these values into a semicolon separated list 
      comb_entry <- paste(diff_beta, paste(diff_beta_binary, paste(diff_probe_abs, diff_probe_bin, sep = ";"), sep = ";"), sep = ";")
      tic("make df")
      df_to_return <- as.data.frame(comb_entry)
      colnames(df_to_return) <- patient_id
      toc()
      
      # Return this value for this gene in this patient
      print(df_to_return)
      return(df_to_return)
    })
    toc()
    # Combine these all into a DF
    tic("combine cols and convert to DF")
    meth_vals_all_patients_df <- as.data.frame(list.cbind(meth_vals_all_patients))  #do.call(cbind, meth_vals_all_patients)
    toc()
    # Return the row for this gene
    print(meth_vals_all_patients_df)
    
    return(meth_vals_all_patients_df)
  })
  return(methylation_rows)
}


# Helper function that, given a tumor and normal methylation table, as well as a 
# gene ID, gets the average beta for each gene and subtracts tumor from normal; 
# returns this difference in methylation
extract_avg_beta <- function(meth_table_tumor, meth_table_normal) {
  
  # Get the average beta value for each
  avg_beta_tumor <- mean.default(meth_table_tumor$Beta_value)    # mean.default() is much faster than mean() for relatively small vectors
  avg_beta_normal <- mean.default(meth_table_normal$Beta_value)
  
  # Return the difference
  return(avg_beta_tumor - avg_beta_normal)
}

# Helper function that, given a tumor and normal methylation table, as well as a 
# gene ID, gets the number of fully methylated probes for each gene and 
# subtracts tumor from normal; returns this difference
get_num_probes <- function(meth_table_tumor, meth_table_normal, beta_thres) {
  
  # Subset each table to include only fully methylated probes
  meth_table_tumor_sub <- meth_table_tumor[meth_table_tumor$Beta_value > beta_thres,]
  meth_table_normal_sub <- meth_table_tumor[meth_table_normal$Beta_value > beta_thres,]
  
  # Get the number of probes in each
  num_probes_tumor <- nrow(distinct(meth_table_tumor_sub))
  num_probes_normal <- nrow(distinct(meth_table_normal_sub))
  
  diff <- num_probes_tumor - num_probes_normal
  
  # If one of the probe numbers is 0 and the other is nonzero, also return 1 (otherwise 0)
  binary <- ifelse((num_probes_normal == 0 && num_probes_tumor != 0) || (num_probes_normal != 0 && num_probes_tumor == 0), 1, 0) 
  
  # Return the difference and the binary result of whether one probe is 0 & other is nonzero
  return(c(diff, binary))
}

# Call this function to get the differential methylation table
differential_methylation_df <- get_differential_methylation_df(filenames = filenames, 
                                                               path = path, 
                                                               matched_patient_IDs = matched_patient_IDs, 
                                                               gene_names = gene_names,
                                                               label = "full")


# Turn each entry into a character vector
comb_fasta_f_unlisted <- lapply(comb_fasta_f, function(x) paste(unlist(x), collapse = ""))

# Function takes a combined FASTA file for a particular TF and blasts it
# using BLASTN. 
blast_sequences <- function(comb_fasta_f) {
  results <- lapply(comb_fasta_f[1:3], function(entry) {
    entry_char <- paste(">ID-1\n", as.character(entry), sep = "")
    res <- blastSequences(entry_char, database = "blastn", hitListSize = 10, expect = 0.05, as = "data.frame")
    return(res)
  })
  return(results)
}

blasted_seq <- blast_sequences(comb_fasta_f)


# Inner function that applies the BMIQ package to a given sample's Beta values
normalize <- function(samp_beta_vals) {
  # Design is whether the probe is type 1 or type 2
  design <- rep()
  # Use defaults for all the rest of the EM algorithm parameters for this function
  return(BMIQ(samp_beta_vals, design.v = ))
}
methylation_df <- apply(methylation_df, MARGIN = 2, FUN = normalize)


ids_meth_nonmatches <- unlist(lapply(patients_with_matched, function(x) if(!x %fin% colnames(methylation_df)) return(x)))
ids_exp_nonmatches <- unlist(lapply(patients_with_matched, function(x) if(!x %fin% exp_patient_ids) return(x)))
patients_with_matched <- setdiff(patients_with_matched, ids_exp_nonmatches) # 73



list_of_lists <- vector(mode="list", length = ncol(my_matrix))    # will hold all the TF target lists so they can be returned
names(list_of_lists) <- colnames(my_matrix)

for (i in 1:ncol(my_matrix)) {
  
  # Subset the matrix based on the current column and keep all the rows that have a value of "1" or "-1"
  subset_mat <- my_matrix[(my_matrix[,i] == 1) | (my_matrix[,i] == -1), i, drop = FALSE]
  
  # Put the remaining rownames (potential gene targets) into the list and add to the list of lists
  list_of_lists[[i]] <- rownames(subset_mat)
}
return(list_of_lists)

sink(paste(path, "Saved Output Data Files/BRCA/Curated TF Data/encode_gene_targets_short.txt", sep = ""))
print(encode_lists_b)
sink()

sink(paste(path, "Saved Output Data Files/BRCA/Curated TF Data/transfac_gene_targets_full.txt", sep = ""))
print(transfac_lists_a)
sink()

sink(paste(path, "Saved Output Data Files/BRCA/Curated TF Data/transfac_gene_targets_short.txt", sep = ""))
print(transfac_lists_b)
sink()

# Integrate the ENCODE and TRANSFAC targets for full and short
appendList <- function (x, val) 
{
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  for (v in names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
      appendList(x[[v]], val[[v]])
    else c(x[[v]], val[[v]])
  }
  return(x)
}

full_tf_list <- appendList(encode_lists_a, transfac_lists_a)
sink(paste(path, "Saved Output Data Files/BRCA/Curated TF Data/gene_targets_full.txt", sep = ""))
print(full_tf_list)
sink()

short_tf_list <- appendList(encode_lists_b, transfac_lists_b)
sink(paste(path, "Saved Output Data Files/BRCA/Curated TF Data/gene_targets_short.txt", sep = ""))
print(short_tf_list)
sink()



if(!is.data.frame(meth_table_tumor_clean) & !length(meth_table_tumor_clean) == 0) {
  for (i in 1:length(meth_table_tumor_clean)) {
    meth_tab_tum <- meth_table_tumor_clean[[i]]
    print(meth_tab_tum)
    fwrite(meth_tab_tum, paste(output_path, paste(patient_id, paste("-0_", paste(i, "_clean_methylation.csv", sep = ""), sep = ""), sep = ""), sep = ""))
    avg_beta_tab_tum <- get_avg_betas(meth_tab_tum)
    fwrite(avg_beta_tab_tum, paste(output_path, paste(patient_id, paste("-0_", paste(i, "_avg_betas_per_gene.csv", sep = ""), sep = ""), sep = ""), sep = ""))
    num_probes_tab_tumor <- get_abs_num(meth_tab_tum, thres)
    fwrite(num_probes_tab_tumor, paste(output_path, paste(patient_id, paste("-0_", paste(i, "_num_methylated_probes_per_gene.csv", sep = ""), sep = ""), sep = ""), sep = ""))
    if (is.na(avg_beta_tab_tumor)) {avg_beta_tab_tumor <- avg_beta_tab_tum}
    else {
      avg_beta_tab_tumor <- merge(avg_beta_tab_tumor, avg_beta_tab_tum, by.x = "Gene_Symbol", by.y = "Gene_Symbol")
      avg_beta_tab_tumor <- data.frame(Gene_Symbol = avg_beta_tab_tumor$Gene_Symbol, Beta_value = rowMeans(avg_beta_tab_tumor[,-1]))
    }
  }
  #avg_beta_tab_tumor_tmp <- abind(beta_tab_tumor_list, along = length(beta_tab_tumor_list))
} else {
  fwrite(meth_table_tumor_clean, paste(output_path, paste(patient_id, "-0_clean_methylation.csv", sep = ""), sep = ""))
  avg_beta_tab_tumor <- get_avg_betas(meth_table_tumor_clean)
  fwrite(avg_beta_tab_tumor, paste(output_path, paste(patient_id, "-0_avg_betas_per_gene.csv", sep = ""), sep = ""))
  num_probes_tab_tumor <- get_abs_num(meth_table_tumor_clean, thres)
  fwrite(num_probes_tab_tumor, paste(output_path, paste(patient_id, "-0_num_methylated_probes_per_gene.csv", sep = ""), sep = ""))
}
if(!is.data.frame(meth_table_normal_clean) & !length(meth_table_normal_clean) == 0) {
  for (i in 1:length(meth_table_normal_clean)) {
    meth_tab_clean <- meth_table_normal_clean[[i]]
    print(meth_tab_clean)
    fwrite(meth_tab_clean, paste(output_path, paste(patient_id, paste("-11", paste(i, "_clean_methylation.csv", sep = ""), sep = ""), sep = ""), sep = ""))
    avg_beta_tab_norm <- get_avg_betas(meth_tab_clean)
    fwrite(avg_beta_tab_norm, paste(output_path, paste(patient_id, paste("-11_", paste(i, "_avg_betas_per_gene.csv", sep = ""), sep = ""), sep = ""), sep = ""))
    num_probes_tab_normal <- get_abs_num(meth_tab_clean, thres)
    fwrite(num_probes_tab_normal, paste(output_path, paste(patient_id, paste("-11", paste(i, "_num_methylated_probes_per_gene.csv", sep = ""), sep = ""), sep = ""), sep = ""))
    if(is.na(avg_beta_tab_normal)) {avg_beta_tab_normal <- avg_beta_tab_norm}
    else {
      avg_beta_tab_normal <- merge(avg_beta_tab_normal, avg_beta_tab_norm, by.x = "Gene_Symbol", by.y = "Gene_Symbol")
      avg_beta_tab_normal <- data.frame(Gene_Symbol = avg_beta_tab_normal$Gene_Symbol, Beta_value = rowMeans(avg_beta_tab_normal[,-1]))
    }
  }
} else {
  if(!length(meth_table_normal_clean) == 0) {
    fwrite(meth_table_normal_clean, paste(output_path, paste(patient_id, "-11_clean_methylation.csv", sep = ""), sep = ""))
    avg_beta_tab_normal <- get_avg_betas(meth_table_normal_clean)
    fwrite(avg_beta_tab_normal, paste(output_path, paste(patient_id, "-11_avg_betas_per_gene.csv", sep = ""), sep = ""))
    num_probes_tab_normal <- get_abs_num(meth_table_normal_clean, thres)
    fwrite(num_probes_tab_normal, paste(output_path, paste(patient_id, "-11_num_methylated_probes_per_gene.csv", sep = ""), sep = ""))
    
  }
}

#gn_matrix <- data.table::fread(full_fn, header = TRUE, check.names = FALSE, 
#blank.lines.skip = TRUE, fill = TRUE)
# Read in using the vcfR package
#vcf <- read.vcfR(full_fn, verbose = FALSE)


add_library_sizes <- function(sample_df, expression_df) {
  library_sizes <- lapply(rownames(sample_df), function(sample) {
    # Get the column(s) of the patient of interest in the expression data frame
    expression_df_sub <- expression_df[,grepl(patient, colnames(expression_df))]
    
    # If there is more than one sample for this patient...
    if(is.data.frame(expression_df_sub)) {
      # Get the tumor and normal columns separately
      expression_df_tum <- expression_df_sub[,grepl("-0", colnames(expression_df_sub))]
      expression_df_norm <- expression_df_sub[,grepl("-11", colnames(expression_df_sub))]
      
      # If there are more than one tumor or normal samples, average them together
      reads_tum <- c()
      reads_norm <- c()
      if(is.data.frame(expression_df_tum)) {
        reads_tum <- rowMeans(expression_df_tum, na.rm = TRUE)
      } else if (length(expression_df_tum) > 0) {
        reads_tum <- expression_df_tum
      } else {reads_tum <- NA}
      if(is.data.frame(expression_df_norm)) {
        reads_norm <- rowMeans(expression_df_norm, na.rm = TRUE)
      } else if (length(expression_df_norm) > 0) {
        reads_norm <- expression_df_norm
      } else {reads_norm <- NA}
      
      # Get the sum of all the reads for both tumor and normal
      sum_tumor <- sum(reads_tum, na.rm = FALSE)
      sum_normal <- sum(reads_norm, na.rm = FALSE)
      
      return(c(sum_tumor, sum_normal))
    }
    # If there is only one sample for this patient...
    else {
      # Retrieve what kind of sample this is
      patient_id_full <- colnames(expression_df)[grepl(patient, colnames(expression_df))]
      
      sum_tum <- NA
      sum_norm <- NA
      
      # If this is the tumor sample, take the sum
      if (grepl("-0", patient_id_full)) {sum_tum <- sum(expression_df_sub)}
      else {sum_norm <- sum(expression_df_sub)}
      
      return(c(sum_tum, sum_norm))
    }
  })
  
  # Combine all results into new DF and cbind this to the patient DF
  lib_size_df <- do.call(rbind, library_sizes)
  colnames(lib_size_df) <- c("Lib_Size_Tumor", "Lib_Size_Norm")
  patient_df <- cbind(patient_df, lib_size_df)
  
  return(patient_df)
}


# If there is more than one sample for this patient...
if(is.data.frame(samples_norm_df_sub)) {
  # Get the tumor and normal columns separately
  samples_norm_df_tum <- samples_norm_df_sub[,grepl("-0", colnames(samples_norm_df_sub))]
  samples_norm_df_norm <- samples_norm_df_sub[,grepl("-11", colnames(samples_norm_df_sub))]
  
  # If there are more than one tumor or normal samples, average them together
  lib_size_tum <- 0
  lib_size_norm <- 0
  if(is.data.frame(samples_norm_df_tum)) {
    lib_size_tum <- mean.default(samples_norm_df_tum$lib.size)
  } else if (length(samples_norm_df_tum) > 0) {
    lib_size_tum <- samples_norm_df_tum
  } else {lib_size_tum <- NA}
  if(is.data.frame(expression_df_norm)) {
    lib_size_norm <- rowMeans(expression_df_norm, na.rm = TRUE)
  } else if (length(expression_df_norm) > 0) {
    lib_size_norm <- expression_df_norm
  } else {lib_size_norm <- NA}
  
  return(c(lib_size_tum, lib_size_norm))
}
# If there is only one sample for this patient...
else {
  # Retrieve what kind of sample this is
  patient_id_full <- colnames(samples_norm_df_sub)[grepl(patient, colnames(samples_norm_df_sub))]
  
  lib_size_tum <- NA
  lib_size_norm <- NA
  
  # Determine if it is tumor or normal
  if (grepl("-0", patient_id_full)) {lib_size_tum <- samples_norm_df_sub[4]}  #4 is the lib.size column
  else {lib_size_norm <- samples_norm_df_sub[4]}
  
  return(c(lib_size_tum, lib_size_norm))
}




# If there are multiple tumor samples for this patient, average them
if (nrow(patient_sub) > 1) {
  # Get the CPE values, if there are any
  cpe_vals <- patient_sub$CPE
  cpe_vals <- cpe_vals[!is.nan(cpe_vals)]
  if (length(cpe_vals) == 1) {return(cpe_vals)}
  else if (length(cpe_vals) > 1) {return(mean.default(cpe_vals))}
  else {
    vals <- unlist(patient_sub[,3:6])
    vals <- vals[!is.nan(vals)]
    return(median(vals))
  }
}
# Otherwise, just take the value for the one sample
else {
  
}
}
} else {
  return(NA)
}





genotypes <- unlist(lapply(1:nrow(maf), function(i) {
  ref_allele <- maf[i, 'Reference_Allele']
  tumor_allele1 <- maf[i, 'Tumor_Seq_Allele1']
  tumor_allele2 <- maf[i, 'Tumor_Seq_Allele2']
  
  if (ref_allele == "-" | tumor_allele1 == "-" | tumor_allele2 == "-") {return(9)}
  if (ref_allele == tumor_allele1) {
    if (ref_allele == tumor_allele2) {return(2)} 
    else {return(1)}
  } else {
    if (ref_allele == tumor_allele2) {return(1)}
    else {return(0)}
  }
}))

# Get all the unique SNPs
#if ('vcf_region' %in% colnames(maf)) {
#snps <- unique(unlist(lapply(maf$vcf_region, function(x) {
#return(unlist(strsplit(x, ":", fixed = TRUE))[3])
#})))
#} else {
#snps <- unique(maf$dbSNP_RS)
#}

# Get all the patient sample IDs
#samples <- unique(maf$Tumor_Sample_Barcode)

# Loop through all the SNPs and see what each patient is
#genotypes <- lapply(1:length(snps), function(snp) {
# Filter the matrix to only look at this SNP
#if ('vcf_region' %in% colnames(maf)) {maf_sub <- maf_sub[maf_sub$vcf_region == snp,]}
#else {maf_sub <- maf_sub[maf_sub$dbSNP_RS == snp,]}

# Get the IDs of the patients that have this SNP
#pats_with_snp <- which(samples %fin% maf_sub$Tumor_Sample_Barcode)


#})


# METHYLATION ARCHIVES
#filename_tumor <- filenames_patient[grepl("-0", filenames_sample_ids)]
#filename_normal <- filenames_patient[grepl("-11", filenames_sample_ids)]

if(length(filename_tumor) == 1) {
  meth_table_tumor <- fread(paste(path, filename_tumor, sep = ""), 
                            sep = "\t", header = TRUE, 
                            select = c(Beta_value='numeric', Gene_Symbol='character'))
} else {
  meth_table_tumor <- lapply(filename_tumor, function(x) 
    fread(paste(path, x, sep = ""), sep = "\t", header = TRUE, 
          select = c(Beta_value='numeric', Gene_Symbol='character')))
}
if (length(filename_normal) == 1) {
  meth_table_normal <- fread(paste(path, filename_normal, sep = ""), 
                             sep = "\t", header = TRUE, 
                             select = c(Beta_value='numeric', Gene_Symbol='character'))
} else {
  meth_table_normal <- lapply(filename_normal, function(x) 
    fread(paste(path, x, sep = ""), sep = "\t", header = TRUE, 
          select = c(Beta_value='numeric', Gene_Symbol='character')))
}


meth_table_tumor_clean <- clean_methylation_table(meth_table_tumor)
#print("cleaned methylation table tumor:")
#print(head(meth_table_tumor_clean))
meth_table_normal_clean <- clean_methylation_table(meth_table_normal)
#print("cleaned methylation table normal:")
#print(head(meth_table_normal_clean))

fwrite(meth_table_tumor_clean, paste(output_path, paste(patient_id, paste(i, "_clean_methylation.csv", sep = ""), sep = "-"), sep = ""))

fwrite(meth_table_tumor_clean, paste(output_path, paste(patient_id, paste(i, "_clean_methylation.csv", sep = ""), sep = "-"), sep = ""))
fwrite(meth_table_normal_clean, paste(output_path, paste(patient_id, 
                                                         paste("-norm", 
                                                               paste(i, "_clean_methylation.csv", 
                                                                     sep = ""), sep = "_"), sep = "-"), sep = ""))

#' For all rows in the methylation table, get the Beta value for each gene name 
#' in that row and create a new, independent row
#' @param meth_table the methylation table to be cleaned
clean_methylation_table <- function(meth_table) {
  
  # If there are multiple tables, clean each table and return them as a list
  if (!is.data.frame(meth_table)) {   # i.e. if a list of data frames 
    tab <- lapply(meth_table, function(x) {
      if (!length(x) == 0) {
        if (is.data.frame(x) | is.data.table(x)) {return(clean_tab(x))}
        else (return(NA))
      }
    })
    tab <- tab[!is.na(tab)]
    return(tab)
    
    # Otherwise, just clean the one table & return it
  } else {
    if (!length(meth_table) == 0) {
      if (is.data.frame(meth_table) | is.data.table(meth_table)) {return(clean_tab(meth_table))}
    } else {return(NA)}
  }
}


# Take the log of the expression values & library sizes
#linear_model_input_table$ExpStat_k <- unlist(lapply(linear_model_input_table$ExpStat_k), log)
#linear_model_input_table$Lib_Size <- unlist(lapply(linear_model_input_table$Lib_Size), log)
#lm_fit <- .lm.fit(cbind(1, linear_model_input_table[,terms]),
#linear_model_input_table[,"ExpStat_k"],
#offset = linear_model_input_table[,"Lib_Size"])

#linear_model_input_table$MethStat_k <- unlist(lapply(linear_model_input_table$MethStat_k), log)
#linear_model_input_table$Lib_Size <- unlist(lapply(linear_model_input_table$Lib_Size), log)
#lm_fit <- .lm.fit(cbind(1, linear_model_input_table[,terms]),
#linear_model_input_table[,"MethStat_k"],
#offset = linear_model_input_table[,"Lib_Size"])

if (analysis_type == "eQTL") {
  lm_fit <- speedglm::speedlm(formula = (log2(ExpStat_k) ~ MutStat_i + CNAStat_i + MethStat_i +
                                           MutStat_k + CNAStat_k + Age_b3 + Age_b4 + Age_b5 + Age_b6 + Age_b7 + 
                                           Age_b8 + Age_b9 + Prior_malig + Treatment_rad + Treatment_pharm + 
                                           Tot_Mut_b1 + Tot_Mut_b2 + Tumor_purity_b1 + Tumor_purity_b2 + 
                                           Tot_IC_Frac_b1 + Tot_IC_Frac_b2), 
                              offset = log2(Lib_Size),
                              data = linear_model_input_table)
  # TODO: add tumor subtype, PCs, PEER factors, gender, etc. -- get these values from the table?
  
} else if (analysis_type == "meQTL") {
  lm_fit <- lm(formula = (log2(MethStat_k) ~ MutStat_i + CNAStat_i + MethStat_i +
                            MutStat_k + CNAStat_k + Age_b3 + Age_b4 + Age_b5 + Age_b6 + Age_b7 + 
                            Age_b8 + Age_b9 + Prior_malig + Treatment_rad + Treatment_pharm + 
                            Tot_Mut_b1 + Tot_Mut_b2 + Tumor_purity_b1 + Tumor_purity_b2 + 
                            Tot_IC_Frac_b1 + Tot_IC_Frac_b2), 
               offset = log2(Lib_Size),
               data = linear_model_input_table)
  
} else {
  print(paste("Invalid analysis type:", analysis_type))
  return(NA)
}

construct_formula <- function(lm_input_table, analysis_type) {
  colnames_to_incl <- colnames(lm_input_table)
  
  # Note: currently excluding age b1 and age b2 since no patients fall into these categories
  colnames_to_incl <- colnames_to_incl[!((colnames_to_incl == "Age_b1") | (colnames_to_incl == "Age_b2"))]
  
  # Remove the race variable entirely
  colnames_to_incl <- colnames_to_incl[!(grepl("Race", colnames_to_incl))]
  
  # For the bucketed variables, remove the last bucket (we don't need it!)
  bucketed_vars <- c("Tot_Mut_b","Tumor_purity_b", "Tot_IC_Frac_b", "Cancer_type_b")
  vals_to_remove <- unlist(lapply(bucketed_vars, function(var) {
    # Get the last bucket for this variable
    matching_vars <- colnames_to_incl[grepl(var, colnames_to_incl)]
    index <- which(colnames_to_incl == matching_vars[length(matching_vars)])
    return(index)
  }))
  colnames_to_incl <- colnames_to_incl[-vars_to_remove]
  
  # Exclude the library size, as we are using this as an offset instead
  colnames_to_incl <- colnames_to_incl[!(colnames_to_incl == "Lib_Size")]
  
  if (analysis_type == "eQTL") {      
    colnames_to_incl <- colnames_to_incl[!(colnames_to_incl == "ExpStat_k")]
    formula <- paste(colnames_to_incl, collapse = " + ") 
    formula <- paste("log2(ExpStat_k) ~ ", formula, sep = "")
  } else if (analysis_type == "meQTL") {
    colnames_to_incl <- colnames_to_incl[!(colnames_to_incl == "MethStat_k")]
    formula <- paste(colnames_to_incl, collapse = " + ") 
    formula <- paste("log2(MethStat_k) ~ ", formula, sep = "")
  } else {
    print(paste("Invalid analysis type:", analysis_type))
    return(NA)
  }
  return(formula)
}

cna_val <- unique(as.numeric(cna_df[,colnames(cna_df) == sample, with = FALSE]))
#if(length(cna_stat) > 1) {cna_stat <- cna_stat[1] + 1}   # Add 1 as a pseudocount when taking the log
#else {cna_stat <- cna_stat + 1}
# OPT 2: BUCKETED CNA VALUE
# Opt. 2a. Include patients with amplifications
# Opt. 2b. Exclude patients with amplifications
if (length(cna_val) > 1) {
  if (0 %in% cna_val) {cna_stat <- -1}   # deletions
  #else {cna_stat <- 0}   # a. amplifications or no change
  else if (1 %in% cna_val) {cna_stat <- 0}  # b. one copy
  #else {cna_stat <- NA} # b. exclude amplifications
  else {cna_stat <- 1} # at least 2 copies
  #else if (2 %in% cna_val) {cna_stat <- 1}
  #else {cna_stat <- NA}
} else {
  if(length(cna_val) == 0) {cna_stat <- NA}
  else if (is.na(cna_val)) {cna_stat <- NA}
  else if(cna_val == 0) {cna_stat <- -1}
  else if (cna_val == 1) {cna_stat <- 0}
  else {cna_stat <- 1}
  #else if (cna_val == 2) {cna_stat <- 1}
  #else {cna_stat <- NA}
}


# Take the mean of all cancer expression values
#if(is.data.frame(patient_exp_cols)) {
#if (nrow(patient_exp_cols) > 1) {
#patient_exp_cols <- colMeans(patient_exp_cols)
#exp_stat <- log(mean.default(patient_exp_cols))
#} else {
# exp_stat <- log(mean.default(as.numeric(patient_exp_cols[,grepl("-0", colnames(patient_exp_cols))])))
#}
#} else {exp_stat <- log(mean.default(patient_exp_cols))}

get_tnm_methylation <- function(methylation_df, analysis_type, sample) {
  
  ### eQTL ###
  if (analysis_type == "eQTL") {
    
    # Extract the methylation value for this patient in this gene
    meth_val <- methylation_df[,colnames(methylation_df) == patient]
    
    # If there is more than one, take the average
    if (length(meth_val) > 1) {meth_val <- mean.default(meth_val)}
    
    # As long as this methylation value is not NA, either take the raw methylation value
    # or calculate whether the gene is differentially methylated
    if(is.na(meth_val)) {meth_stat <- NA} 
    else {
      if (length(meth_val) == 0) {meth_stat <- NA}
      # OPT 1. (RAW METHYLATION): Does this target gene have a methylation marker in cancer?
      meth_stat <- meth_val
      # OPT 2. (DIFFERENTIAL METHYLATION): Does this target gene have a differential methylation
      # state between tumor and normal?
      #else {
      #meth_stat <- 0
      #meth_beta_thres <- 0.25
      #if (meth_val > meth_beta_thres | meth_val < -meth_beta_thres) {meth_stat <- 1}
      #}
    }
    # Return the methylation stat
    return(meth_stat)
  }
  
  ### meQTL ###
  else if (analysis_type == "meQTL") {
    
    # Get the log fold change in methylation
    meth_cols <- methylation_df[,grepl(patient, colnames(methylation_df))]
    meth_val_tumor <- meth_cols[,grepl("norm", colnames(meth_cols))]
    meth_val_norm <- meth_cols[,grepl("tumor", colnames(meth_cols))]
    log_fc_meth <- log(meth_val_tumor / meth_val_norm)
    
    return(log_fc_meth)
  }
}



# OLD LINEAR_MODEL.R file
############################################################
### Linear Model
### Written By: Sara Camilli, July 2020
############################################################
library(biomaRt)
library(broom)
library(qvalue)
library(parallel)
library(rlang)
library(dplyr)
library(rlist)
library(rockchalk)
library(enrichvs)

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)

############################################################
# GET ALL REGULATORY PROTEINS IN CONSIDERATION
############################################################
# Regulatory proteins will vary depending on the level of specificity in question;
# import only the proteins at the given specificity level of interest under label "protein_ids_df"
prot_path <- paste(main_path, "Mutation/Files for Linear Model", sep = "")

protein_ids_df <- read.csv(paste(prot_path, "iprotein_protein_ids_df", sep = ""), header = TRUE, check.names = FALSE)
# protein_ids_df <- read.csv(paste(prot_path, "iprotein_nucacids_protein_ids_df", sep = ""), header = TRUE, check.names = FALSE)
# protein_ids_df <- read.csv(paste(prot_path, "idomain_protein_ids_df", sep = ""), header = TRUE, check.names = FALSE)
# protein_ids_df <- read.csv(paste(prot_path, "idomain_nucacids_protein_ids_df", sep = ""), header = TRUE, check.names = FALSE)
# protein_ids_df <- read.csv(paste(prot_path, "ibindingpos_protein_ids_df", sep = ""), header = TRUE, check.names = FALSE)
# protein_ids_df <- read.csv(paste(prot_path, "ibindingpos_nucacids_protein_ids_df", sep = ""), header = TRUE, check.names = FALSE)


############################################################
# IMPORT EXPRESSION FILES
############################################################
### TUMOR-NORMAL MATCHED ###
# 1. Unfiltered FPKM
#expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_DF_tumNormMatched.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
#expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_DF_tumNormMatched_meQTL.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)

# 2. Filtered FPKM
expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_filt_DF_tumNormMatched.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_filt_DF_tumNormMatched_meQTL.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)

# 3. TMM
expression_df <- read.csv(paste(main_path, "Expression/expression_tmm_DF_tumNormMatched.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
expression_df <- read.csv(paste(main_path, "Expression/expression_tmm_DF_tumNormMatched_meQTL.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)


### NON-TUMOR-NORMAL MATCHED ###
# FPKM (unfiltered)
#expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_DF.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)  # fpkm-normalized

# FPKM (filtered)
expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_filt_DF.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)  # fpkm-normalized and filtered

# TMM 
#expression_df <- read.csv(paste(main_path, "Expression/tmm_normalized_expression_counts.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
expression_df <- read.csv(paste(main_path, "Expression/expression_tmm_DF_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)

############################################################
# IMPORT METHYLATION FILES
############################################################
### TUMOR-NORMAL MATCHED ###
methylation_df <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_tumNormMatched.csv", sep = ""), header = TRUE, row.names = 1)
methylation_df <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_tumNormMatched_eQTL.csv", sep = ""), header = TRUE, row.names = 1)
methylation_df_dependent <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_tumNormMatched_dependent.csv", sep = ""), header = TRUE, row.names = 1)
methylation_df_dependent <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_tumNormMatched_dependent_eQTL.csv", sep = ""), header = TRUE, row.names = 1)


### NON-TUMOR-NORMAL MATCHED ###
methylation_df <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_cancer_only_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)

############################################################
# IMPORT GENE TARGET MUTATION FILES
############################################################
### TUMOR-NORMAL MATCHED ###
# For eQTL
mutation_targ_df <- read.csv(paste(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense_tumNormMatched.csv", sep = ""), header = TRUE, row.names = 1)
# For meQTL
mutation_targ_df <- read.csv(paste(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense_tumNormMatched_meQTL.csv", sep = ""), header = TRUE, row.names = 1)

### NON-TUMOR-NORMAL MATCHED ###
mutation_targ_df <- read.csv(paste(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
#mutation_targ_df$swissprot <- unlist(lapply(rownames(mutation_targ_df), function(x) paste(unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'uniprot_gn_id']), collapse = ";")))

############################################################
# IMPORT REGULATORY PROTEIN MUTATION FILES
############################################################
### TUMOR-NORMAL MATCHED ###
# For eQTL
mutation_regprot_df <- read.csv(paste(main_path, "Mutation/iprotein_results_missense_tumNormMatched.csv", sep = ""), header = TRUE, row.names = 1)
# For meQTL
mutation_regprot_df <- read.csv(paste(main_path, "Mutation/iprotein_results_missense_tumNormMatched_meQTL.csv", sep = ""), header = TRUE, row.names = 1)

### NON-TUMOR-NORMAL MATCHED ###
mutation_regprot_df <- read.csv(paste(main_path, "Mutation/iprotein_results_missense_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
#mutation_regprot_df$Swissprot <- unlist(lapply(mutation_regprot_df$Query, function(x) unlist(strsplit(x, "|", fixed = TRUE))[2]))

############################################################
# IMPORT CNA FILES
############################################################
### TUMOR-NORMAL MATCHED ###
# For eQTL
cna_df <- read.csv(paste(main_path, "CNV/Gene-level Score/CNV_DF_AllGenes_CancerOnly_tumNormMatched.csv", sep = ""), header = TRUE, row.names = 1)
# For meQTL
cna_df <- read.csv(paste(main_path, "CNV/Gene-level Score/CNV_DF_AllGenes_CancerOnly_tumNormMatched_meQTL.csv", sep = ""), header = TRUE, row.names = 1)

# If rownames are unlabeled
# cna_df_full <- read.csv(paste(main_path, "CNV/Gene-level Score/CNV_DF_AllGenes.csv", sep = ""), header = TRUE, row.names = 1)
# gene_ids <- rownames(cna_df_full)
# rownames(cna_df) <- gene_ids
# rm(cna_df_full)

### NON-TUMOR-NORMAL MATCHED ###
cna_df <- read.csv(paste(main_path, "CNV/Gene-level Raw/CNA_DF_AllGenes_CancerOnly_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
rownames(cna_df) <- unlist(lapply(rownames(cna_df), function(x) unlist(strsplit(x, ".", fixed = TRUE))[1]))

############################################################
# IMPORT PATIENT FILES
############################################################
### TUMOR-NORMAL MATCHED ###
# For eQTL
patient_df <- read.csv(paste(main_path, "Linear Model/patient_dataframe_tumNormMatched.csv", sep = ""), row.names = 1)
# For meQTL
patient_df <- read.csv(paste(main_path, "Linear Model/patient_dataframe_tumNormMatched_meQTL.csv", sep = ""), row.names =  1)


### NON-TUMOR-NORMAL MATCHED ###
patient_df <- read.csv(paste(main_path, "Linear Model/patient_dataframe_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)


############################################################
# SET UP TARGETS FOR PROOF-OF-CONCEPT TESTING
############################################################
sample_protein_uniprot <- "P04637" 
sample_protein_ensg <- "ENSG00000141510"
sample_protein_name <- "P53"

sample_protein_uniprot <- "P55317" 
sample_protein_ensg <- "ENSG00000129514"
sample_protein_name <- "FOXA1"

protein_ids_df_tester <- data.frame(swissprot_ids = sample_protein_uniprot, ensg_ids = sample_protein_ensg)

# Get sample targets data frame
# Option 1: Curated Targets
sample_targets_DF <- read.csv(paste(main_path, "Linear Model/TP53/tp53_curated_targets.csv", sep = ""), row.names = 1, check.names = FALSE)
sample_targets_DF <- read.csv(paste(main_path, "Linear Model/FOXA1/foxa1_curated_targets.csv", sep = ""), row.names = 1, check.names = FALSE)

# Option 2: ChIP-eat Targets
sample_targets_DF <- read.csv(paste(main_path, "Linear Model/TP53/tp53_chipeat_targets.csv", sep = ""), row.names = 1, check.names = FALSE)
sample_targets_DF <- read.csv(paste(main_path, "Linear Model/FOXA1/foxa1_chipeat_targets.csv", sep = ""), row.names = 1, check.names = FALSE)

# Option 3: All Gene Targets
sample_targets_DF <- read.csv(paste(main_path, "Linear Model/allgene_targets.csv", sep = ""), row.names = 1, check.names = FALSE)


# Regardless of method, limit targets to only those that overlap the genes in the expression DF
rows_to_keep <- unlist(lapply(1:length(sample_targets_DF$ensg), function(i) {
  targs <- unlist(strsplit(sample_targets_DF[i,'ensg'], ";", fixed = TRUE))
  if (any(targs %fin% rownames(expression_df))) {return(TRUE)}
  else {return(FALSE)}
}))
sample_targets_DF <- sample_targets_DF[rows_to_keep,]


############################################################
#### OVERVIEW OF LINEAR MODEL
############################################################
### PATIENT CHARACTERISTICS ###
# Gender : Gender of all patients 1..j..M (0 if male, 1 if female) -- ONLY INCLUDE IF NOT BRCA
# Age : Age of all patients 1..j..M in years 
# Buckets: 1 if 0-9, 2 if 10-19, 3 if 20-29, 4 if 30-39, 5 if 40-49, 6 if 50-59, 7 if 60-69, 8 if 70-79, 9 if 80-89, 10 if 90+
# Race : Race/ Ethnicity of all patients 1..j..M 
# Buckets: 1 if White (not Latinx), 2 if White (Latinx), 3 if black or African, 4 if Asian, 5 if other
# Prior_malig : If all patients 1..j..M had a prior malignancy (0 if no, 1 if yes)
# Treatment_rad : If all patients 1..j..M were treated with radiation (0 is not treated, 1 is treated) 
# Treatment_pharm : If all patients 1..j..M were treated with pharmaceutical therapy (0 is not treated, 1 is treated) 
# TotalNumMut : Total number of missense mutations each patient j has (across all patients 1..j..M), integer value

### REGULATORY PROTEIN i CHARACTERISTICS ###
# MutStat_i : The mutation status of protein i (0 if not mutated, 1 if mutated) across all patients 1..j..M
# MethStat_i : The methylation status of protein i (0 if not methylated, 1 if methylated) across all patients 1..j..M  
# Alternative: 0 if not differentially methylated (cancer v. normal), 1 if differentially methylated (cancer v. normal)
# CNA_i : The copy number alteration status of protein i (1 if amplified, -1 if deleted, 0 if no amp. or del. OR copy #) across all patients 1..j..M

### TARGET GENE k CHARACTERISTICS ###
# MutStat_k : The mutation status of gene k (0 if not mutated, 1 if mutated) across all patients 1..j..M
# OPT 1. MethStat_k : The methylation status of protein k (0 if not methylated, 1 if methylated) across all patients 1..j..M  
# Alternative: 0 if not differentially methylated (cancer v. normal), 1 if differentially methylated (cancer v. normal)
# OPT 2. ExpStat_k : The differential expression status of protein i (1 if differentially methylated cancer v. normal, 0 if
# not differentially methylated cancer v. normal)

### DEPENDENT VARIABLE ###
# OPT 1. log(Exp_k,c/Exp_k,n) : the log fold change of the differential expression of target gene k across all patients 1..j..M
# OPT 2. log(Meth_k,c/Meth_k,n) : the log fold change of the differential methylation beta or m-values of target gene k across all patients 1..j..M
# OPT 3. Exp_k,c or Meth_k,c : the cancer expression value (for non-tumor-normal matched runs)

### OFFSETS/ SCALING FACTORS ###
# log(N_j) : log of effective library size (sample-specific scaling factor); the total number of reads for 
# patient j after cross-sample normalization, across all patients 1..j..M
# log(CNA_k,j) : the copy number for target gene k across all patients 1..j..M


############################################################
#### MAIN LINEAR MODEL FUNCTION
############################################################
#' MAIN FUNCTION: This function will run the linear model, taking in the data 
#' frame we've created above.
#' @param protein_ids_df DF of the regulatory proteins of interest, found based 
#' on their DNA-binding regions (columns for swissprot IDs and ENSG IDs)
#' @param downstream_target_df table of regulatory proteins (columns) with gene 
#' targets we're testing for each (entries)
#' @param patient_df table of patient-specific information (sex, age, total number 
#' of mutations, treated or not)
#' @param mutation_df_targ table of mutation counts for each patient, for each 
#' gene in the genome
#' @param mutation_df_regprot a dataframe of which patients have mutations at 
#' a given level of specificity in particular regulatory proteins
#' @param methylation_df table of methylation results (rows are proteins, columns 
#' are patients, entries are methylation values)
#' @param cna_df table of CNA per gene (rows are proteins, columns are patients, 
#' entries are CNA values)
#' @param expression_df table of expression (rows are genes, columns are patients, 
#' entries are expression values)
#' @param analysis_type a string label that reads either "eQTL" or "meQTL" to 
#' determine what kind of model we should run
#' @param tumor_norm_matched a TRUE/FALSE value indicating whether we are looking 
#' at tumor-normal matched data or not
run_linear_model <- function(protein_ids_df, downstream_target_df, patient_df, mutation_df_targ, mutation_df_regprot, 
                             methylation_df, cna_df, expression_df, analysis_type, tumor_norm_matched) {
  
  # We need to get a mini-table for every r_i, t_k combo that we rbind into a master table of results
  results_df_list <- mclapply(1:nrow(protein_ids_df), function(i) {
    
    # Get the given regulatory protein r_i's Swissprot & ENSG IDs
    regprot <- protein_ids_df$swissprot_ids[i]
    regprot_ensg <- unlist(strsplit(protein_ids_df$ensg_ids[i], ";", fixed = TRUE))
    
    # Create a starter table that gets the mutation, CNA, and methylation status 
    # for this regulatory protein and binds it to the patient DF
    starter_df <- fill_regprot_inputs(patient_df, regprot, regprot_ensg, 
                                      mutation_df_regprot, methylation_df, 
                                      cna_df, analysis_type)
    
    ### OPTIONAL ###
    # Randomize the table's protein i mutations across all the patients as a test
    # starter_df$MutStat_i <- sample(starter_df$MutStat_i)
    
    print(starter_df)
    
    # Loop through all the targets for this protein of interest and create a table 
    # for each of them from this starter DF
    regprot_i_results_df_list <- lapply(1:nrow(downstream_target_df), function(k) {
      
      # Get the target t_k's Swissprot & ENSG IDs
      targ <- unlist(strsplit(downstream_target_df$swissprot[k], ";", fixed = TRUE))
      targ_ensg <- unlist(strsplit(downstream_target_df$ensg[k], ";", fixed = TRUE))
      
      # OPT: Filter outlier expression for each group (mutated, unmutated) using a standard 
      # (1.5 x IQR) + Q3 schema
      #starter_df <- filter_expression_df(expression_df, starter_df, targ_ensg)
      
      # Create a full input table to the linear model for this target
      linear_model_input_table <- fill_targ_inputs(starter_df, targ, targ_ensg, 
                                                   mutation_df_targ, methylation_df,
                                                   cna_df, expression_df, analysis_type, 
                                                   tumor_norm_matched)
      
      if(length(linear_model_input_table) == 0) {return(NA)}
      if (is.na(linear_model_input_table)) {return(NA)}
      
      # If not tumor-normal-matched, remove the Lib_Size_Norm column
      if(!tumor_norm_matched) {
        linear_model_input_table <- linear_model_input_table[, -which(colnames(linear_model_input_table) 
                                                                      == 'Lib_Size_Norm')]
      }
      
      # If a row has NA, remove that whole row and predict without it
      #linear_model_input_table_noNA <- na.exclude(linear_model_input_table)
      linear_model_input_table_noNA <- na.omit(linear_model_input_table)
      print(linear_model_input_table_noNA)
      
      # Coerse all columns to be numeric
      linear_model_input_table_final <- as.data.frame(sapply(linear_model_input_table_noNA, as.numeric))
      rownames(linear_model_input_table_final) <- rownames(linear_model_input_table_noNA)
      
      print(linear_model_input_table_final)
      
      # Run the linear model on this input table
      # Note: currently excluding age b1 and age b2 since no patients fall into these categories
      if (analysis_type == "eQTL") {
        #lm_fit <- lm(formula = (LogFCExp_k ~ 0 + MutStat_i + MutStat_k + MethStat_k + CNAStat_k + Age_b3 + 
        #Age_b4 + Age_b5 + Age_b6 + Age_b7 + Age_b8 + Age_b9 + Prior_malig + 
        #Treatment_rad + Treatment_pharm + Tot_Mut_b1 + Tot_Mut_b2), 
        #data = linear_model_input_table_noNA)
        # Using CNAStat_k as an offset rather than as a covariate, adding library size as a covariate
        #print(apply(linear_model_input_table_final, 2, function(x) any(is.na(x))))
        #print(apply(linear_model_input_table_final, 2, function(x) any(is.infinite(x))))
        #print(apply(linear_model_input_table_final, 2, function(x) any(is.nan(x))))
        #print(apply(linear_model_input_table_final, 2, function(x) any(!is.numeric(x))))
        #print(apply(linear_model_input_table_final, 2, function(x) any(is.factor(x))))
        
        print(any(linear_model_input_table_final$CNAStat_k == 0))
        print(any(linear_model_input_table_final$Lib_Size_Tumor == 0))
        
        lm_fit <- lm(formula = (LogFCExp_k ~ MutStat_i + MutStat_k + MethStat_k + Age_b3 + 
                                  Age_b4 + Age_b5 + Age_b6 + Age_b7 + Age_b8 + Age_b9 + Prior_malig + 
                                  Treatment_rad + Treatment_pharm + Tot_Mut_b1 + Tot_Mut_b2), 
                     offset = log(CNAStat_k * Lib_Size_Tumor),
                     data = linear_model_input_table_final)
      } else if (analysis_type == "meQTL") {
        lm_fit <- lm(formula = (LogFCMeth_k ~ 0 + MutStat_i + MutStat_k + LogFCExp_k + CNAStat_k + Age_b3 + 
                                  Age_b4 + Age_b5 + Age_b6 + Age_b7 + Age_b8 + Age_b9 + Prior_malig + 
                                  Treatment_rad + Treatment_pharm + Tot_Mut_b1 + Tot_Mut_b2), 
                     data = linear_model_input_table_final)
      } else {
        print(paste("Invalid analysis type:", analysis_type))
        return(NA)
      }
      # Restrict the output table to just particular columns & just the mutation status of regprot i
      lm_fit <- tidy(lm_fit)
      summary_table <- lm_fit[lm_fit$term == "MutStat_i",]
      
      # Add a column for the regulatory protein and target gene to ID them
      summary_table$R_i <- regprot
      summary_table$T_k <- paste(targ, collapse = ";")
      
      print(summary_table)
      
      # Return this summary table
      return(summary_table)
    })
    
    # Now we have a list of output summary DFs for each of this regulatory protein's targets. Bind them
    # all together.
    #return(do.call("rbind", regprot_i_results_df_list))
    regprot_i_results_df_list <- regprot_i_results_df_list[!is.na(regprot_i_results_df_list)]
    return(as.data.frame(data.table::rbindlist(regprot_i_results_df_list)))
  })
  
  # Now we have a list of the resulR omit ts tables, one table per regulatory protein, with entries for each target. 
  # Bind these all together into one master table.
  #master_df <- do.call("rbind", results_df_list)
  master_df <- as.data.frame(data.table::rbindlist(results_df_list))
  
  # Return this master DF
  return(master_df)
}


# Call this function
master_df <- run_linear_model(protein_ids_df, downstream_target_df, patient_df, mutation_targ_df, mutation_regprot_df, 
                              methylation_df, cna_df, expression_df, "eQTL")
master_df <- run_linear_model(protein_ids_df, downstream_target_df, patient_df, mutation_targ_df, mutation_regprot_df, 
                              methylation_df_dependent, cna_df, expression_df, "meQTL")

# For tester
master_df <- run_linear_model(protein_ids_df_tester, sample_targets_DF, patient_df, mutation_targ_df, mutation_regprot_df,
                              methylation_df, cna_df, expression_df, "eQTL", FALSE)
master_df <- run_linear_model(protein_ids_df_tester, sample_targets_DF, patient_df, mutation_targ_df, mutation_regprot_df,
                              methylation_df_dependent, cna_df, expression_df, "meQTL", FALSE)

master_df <- master_df[order(master_df$p.value, decreasing = FALSE),]

# Write the results to a file
write.csv(master_df, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/iprotein_TP53_output_results_df_TMM_uncorrected.csv", sep = ""))

############################################################
#### FILL IN REGULATORY PROTEIN R_I INPUTS TO TABLE
############################################################
#' Function takes a dataframe of inputs for the given patient p_j, as well as the data files 
#' for regulatory protein r_i (and its IDs), and uses them to construct a partial tabular input
#' to a linear model function for a given regulatory protein r_i of the following form:
#' Columns: Linear model input variables
#' Rows: Patients 
#' 
#' @param patient_df a 'starter DF' that has all the patient characteristics 
#' @param regprot_i_uniprot the uniprot ID of regulatory protein i
#' @param regprot_i_ensg the ensembl ID of regulatory protein i
#' @param mutation_regprot_df the mutation DF for regulatory proteins 
#' @param methylation_df the methylation DF
#' @param cna_df the copy number alteration DF
#' @param analysis_type a string denoting the type of analysis, either 'eQTL' or 'meQTL'
fill_regprot_inputs <- function(patient_df, regprot_i_uniprot, regprot_i_ensg, mutation_regprot_df,
                                methylation_df, cna_df, analysis_type) {
  
  # Filter the dataframes to look at only this regulatory protein
  mutation_regprot_df <- mutation_regprot_df %>% dplyr::filter(Swissprot == regprot_i_uniprot)
  cna_df <- filter_cna_by_ensg(cna_df, regprot_i_ensg)
  methylation_df <- filter_meth_by_ensg(methylation_df, regprot_i_ensg)
  
  # Loop through all patients to get the info on this regulatory protein for each
  regprot_rows <- mclapply(1:nrow(patient_df), function(i) {
    patient <- unlist(strsplit(rownames(patient_df)[i], ".", fixed = TRUE))[1]
    
    # Is this regulatory protein mutated in this patient at the given level of specificity?
    # OPT 1: FOR I-DOMAIN & I-BINDING POSITION:
    #if (nrow(mutation_regprot_df %>% filter(Patient == patient)) > 0) {mut_stat <- 1}
    # OPT 2: FOR I-PROTEIN:
    if(nrow(mutation_regprot_df[grepl(patient, mutation_regprot_df$Patient),]) > 0) {mut_stat <- 1}
    else {mut_stat <- 0}
    
    # OPT 1. RAW CNA VALUE: Does this regulatory protein have a CNA?
    cna_stat <- unique(as.numeric(cna_df[,colnames(cna_df) == patient]))
    if(length(cna_stat) > 1) {cna_stat <- cna_stat[1] + 1}   # Add 1 as a pseudocount when taking the log
    else {cna_stat <- cna_stat + 1}
    # OPT 2: BUCKETED CNA VALUE
    # Opt. 2a. Include patients with amplifications
    # Opt. 2b. Exclude patients with amplifications
    #if (length(cna_val) > 1) {
    #if (0 %in% cna_val) {cna_stat <- -1}   # deletions
    #else {cna_stat <- 0}   # a. amplifications or no change
    #else if (1 %in% cna_val) {cna_stat <- 0}  # b. one copy
    #else {cna_stat <- NA} # b. exclude amplifications
    #else {cna_stat <- 1} # at least 2 copies
    #else if (2 %in% cna_val) {cna_stat <- 1}
    #else {cna_stat <- NA}
    #} else {
    #if(cna_val == 0) {cna_stat <- -1}
    #else if (cna_val == 1) {cna_stat <- 0}
    #else {cna_stat <- 1}
    #else if (cna_val == 2) {cna_stat <- 1}
    #else {cna_stat <- NA}
    #}
    
    # OPT 1. (RAW METHYLATION): Does this regulatory protein have a methylation marker in cancer?
    meth_stat <- mean.default(as.numeric(methylation_df[,grepl(patient, colnames(methylation_df))]))
    # OPT 2. (DIFFERENTIAL METHYLATION): Does this regulatory protein have a differential methylation
    # state between tumor and normal?
    #meth_beta_thres <- 0.25
    #meth_stat <- 0
    #if (analysis_type == "eQTL") {
    #meth_val <- as.numeric(methylation_df[,colnames(methylation_df) == patient])
    #if(length(meth_val) == 0) {meth_stat <- NA}
    #else {
    # If the difference between tumor and normal exceeds a threshold, make it 1
    #if (meth_val > meth_beta_thres | meth_val < -meth_beta_thres) {meth_stat <- 1}
    #}
    #} else {   # meQTL
    #patient_tumor_lab <- paste(patient, ".tumor", sep = "")
    #patient_norm_lab <- paste(patient, ".norm", sep = "")
    #meth_val_tum <- as.numeric(methylation_df[,colnames(methylation_df) == patient_tumor_lab])
    #meth_val_norm <- as.numeric(methylation_df[,colnames(methylation_df) == patient_norm_lab])
    #if (is.na(meth_val_norm) | is.na(meth_val_tum)) {meth_stat <- NA}
    #else {
    # difference <- meth_val_tum - meth_val_norm
    # if (difference > meth_beta_thres | difference < -meth_beta_thres) {meth_stat <- 1}
    #}
    #}
    return(data.frame("MutStat_i" = mut_stat, "CNAStat_i" = cna_stat, "MethStat_i" = meth_stat))
  })
  
  # Bind these rows into the starter DF
  #regprot_i_df <- do.call("rbind", regprot_rows)
  regprot_i_df <- as.data.frame(data.table::rbindlist(regprot_rows))
  colnames(regprot_i_df) <- c("MutStat_i", "CNAStat_i", "MethStat_i")
  starter_df <- cbind(patient_df, regprot_i_df)
  
  # Return this full input DF
  return(starter_df)
}


############################################################
#### FILL IN TARGET T_K INPUTS TO TABLE
############################################################
#' Function takes a dataframe of inputs for the given patient p_j and the regulatory protein r_i,
#' as well as the data files for target t_k (and its ID), and uses them to construct a full tabular input
#' to a linear model function for a given regulatory protein i-target gene k pairing of the following form:
#' Columns: Linear model input variables
#' Rows: Patients 
#' @param starter_df a 'starter DF' that has all the patient and reg. protein r_i characteristics 
#' @param targ_k the swissprot ID of target gene k
#' @param targ_k_ensg the ensembl ID of target gene k
#' @param mutation_targ_df the mutation gene target DF
#' @param methylation_df the methylation DF
#' @param cna_df the copy number alteration DF
#' @param expression_df the gene expression DF
#' @param analysis_type a string denoting the type of analysis, either 'eQTL' or 'meQTL'
#' @param tumor_norm_matched a TRUE/FALSE value indicating whether we are using tumor-normal matched
#' or cancer-only data
fill_targ_inputs <- function(starter_df, targ_k, targ_k_ensg, mutation_targ_df,
                             methylation_df, cna_df, expression_df, analysis_type, tumor_norm_matched) {
  print(targ_k)
  print(targ_k_ensg)
  
  # Begin by limiting the dataframes to just this target gene k
  if (!(targ_k == "" | targ_k_ensg == "")) {
    # Mutation target DF
    mutation_targ_df <- filter_mut_by_uniprot(mutation_targ_df, targ_k) 
    # CNA DF
    cna_df <- filter_cna_by_ensg(cna_df, targ_k_ensg)
    # Methylation DF
    methylation_df <- filter_meth_by_ensg(methylation_df, targ_k_ensg)
    # Expression DF
    exp_rows_to_keep <- unique(unlist(lapply(targ_k_ensg, function(x) grep(x, rownames(expression_df)))))
    expression_df <- expression_df[exp_rows_to_keep,]
    print("Expression DF filt by ENSG ID")
    print(expression_df)
    
    # Continue only if all of these dataframes contain information about this target gene k
    if(!(nrow(expression_df) == 0 | nrow(cna_df) == 0 | nrow(methylation_df) == 0 | nrow(mutation_targ_df) == 0)) {
      # Loop through all patients to get the info for each target t_k column
      target_rows <- mclapply(1:nrow(starter_df), function(i) {
        patient <- unlist(strsplit(rownames(starter_df)[i], ".", fixed = TRUE))[1]
        
        # Is this target mutated? 
        mut_count <- mutation_targ_df[,colnames(mutation_targ_df) == patient]
        # If there is no entry here, then it was not mutated in any patients
        if (length(mut_count) == 0) {
          mut_stat <- 0
        } else if (length(mut_count) == 1) {
          if (mut_count > 0) {mut_stat <- 1} 
          else {mut_stat <- 0}
        } else {
          mut_count <- mean.default(mut_count)
          if(mut_count > 0) {mut_stat <- 1} 
          else {mut_stat <- 0}
        }
        
        # Is this target amplified or deleted?
        cna_stat <- unique(as.numeric(cna_df[,colnames(cna_df) == patient])) 
        # OPT 1. RAW CNA VALUE: Does this regulatory protein have a CNA?
        if(length(cna_stat) > 1) {cna_stat <- cna_stat[1] + 1}
        else {cna_stat <- cna_stat + 1}
        # OPT 2: CNA DELETED OR NOT? 
        # Opt. 2a. Include patients with amplifications
        # Opt. 2b. Exclude patients with amplifications
        #if (length(cna_val) > 1) {
        #if (0 %in% cna_val) {cna_stat <- -1}   # no copies
        #else {cna_stat <- 0}   # amplifications or no change
        #else if (1 %in% cna_val) {cna_stat <- 0}  # no change
        #else {cna_stat <- 1} 
        #else if (2 %in% cna_val) {cna_stat <- 1}
        #else {cna_stat <- NA} # exclude amplications
        #} else {
        #if(length(cna_val) == 0) {cna_stat <- NA}
        #else if (is.na(cna_val)) {cna_stat <- NA}
        #else if (cna_val == 0) {cna_stat <- -1}
        #else if (cna_val == 1) {cna_stat <- 0}
        #else {cna_stat <- 1}
        #else if (cna_val == 2) {cna_stat <- 1}
        #else {cna_stat <- NA}
        #}
        
        # What is the log fold expression change/ log expression of this target in cancer?
        patient_exp_cols <- expression_df[,grepl(patient, colnames(expression_df))]
        # If tumor-normal matched, use helper function
        if(tumor_norm_matched) {exp_stat <- log(handle_exp_cases(patient_exp_cols))}
        # Otherwise just take the mean of all cancer values
        else {
          if(is.data.frame(patient_exp_cols)) {
            if (nrow(patient_exp_cols) > 1) {
              patient_exp_cols <- colMeans(patient_exp_cols)
              exp_stat <- log(mean.default(patient_exp_cols))
            } else {
              exp_stat <- log(mean.default(as.numeric(patient_exp_cols[,grepl("-0", colnames(patient_exp_cols))])))
            }
          } else {exp_stat <- log(mean.default(patient_exp_cols))}
        }
        
        # Convert Inf, -Inf, or NaN values to NA
        if(is.infinite(exp_stat) | is.nan(exp_stat)) {exp_stat <- NA}
        
        # Get the methylation statistic using helper function
        meth_stat <- handle_methylation_cases(methylation_df, analysis_type, patient)
        
        #print(paste("exp_stat", exp_stat))
        #print(paste("mut_stat", mut_stat))
        #print(paste("cna_stat", cna_stat))
        #print(paste("meth_stat", meth_stat))
        
        # Make row for this patient and return
        return(data.frame("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, 
                          "CNAStat_k" = cna_stat, "MethStat_k" = meth_stat))
        
      })
      
      # Bind these rows into the starter DF
      targ_k_df <- as.data.frame(data.table::rbindlist(target_rows))
      print(nrow(targ_k_df))
      
      if (analysis_type == "eQTL") {colnames(targ_k_df) <- c("LogFCExp_k", "MutStat_k", "CNAStat_k", "MethStat_k")}
      else {colnames(targ_k_df) <- c("ExpStat_k", "MutStat_k", "CNAStat_k", "LogFCMeth_k")}
      print(head(targ_k_df))
      full_df <- cbind(starter_df, targ_k_df)
      
      # Return this full input DF
      return(full_df)
    }
  }
  return(NA)
}

############################################################
############################################################
### HELPER FUNCTIONS
############################################################
############################################################

#' Takes an expression DF and a starter DF and returns the starter DF without outlier patients,
#' defined as patients with expression of target gene that exceeds 1.5*IQR + Q3 threshold
#' (separately among patients with a mutation in the regulatory protein and patients without)
#' @param expression_df an expression DF, unfiltered 
#' @param starter_df a starter DF for the regulatory protein of interest i
#' @param ensg the ensembl ID of the target gene k
filter_expression_df <- function(expression_df, starter_df, ensg) {
  patients_mut <- rownames(starter_df[starter_df$MutStat_i == 1,])
  patients_nonmut <- setdiff(rownames(starter_df), patients_mut)
  
  # Filter expression DF to target gene of interest and crop patient IDs
  expression_df_filt <- data.frame()
  if(length(ensg) == 1) {
    expression_df_filt <- expression_df[rownames(expression_df) == ensg,]   # ENSG of target gene
  } else {
    rows_exp_filt <- lapply(ensg, function(e) expression_df[rownames(expression_df) == e,])
    expression_df_filt <- do.call("rbind", rows_exp_filt)
  }
  labels <- unlist(lapply(colnames(expression_df_filt), function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))
  expression_df_upd_labels <- expression_df_filt
  colnames(expression_df_upd_labels) <- labels
  
  # Get the expression values for patients in each mutational group (mutated, non-mutated)
  expression_mut <- list()
  expression_nonmut <- list()
  for (i in 1:ncol(expression_df_upd_labels)) {
    patient <- colnames(expression_df_upd_labels)[i]
    if(patient %fin% patients_mut) {
      expression_mut[[patient]] <- as.numeric(expression_df_upd_labels[,i])
    } else {
      expression_nonmut[[patient]] <- as.numeric(expression_df_upd_labels[,i])
    }
  }
  print(expression_mut)
  print(expression_nonmut)
  
  # Get the IQR, Q3, and expression thresholds for each group
  iqr_mut <- IQR(unlist(expression_mut), na.rm = TRUE)
  iqr_nonmut <- IQR(unlist(expression_nonmut), na.rm = TRUE)
  q3_mut <- as.numeric(quantile(unlist(expression_mut), na.rm = TRUE)[4])
  q3_nonmut <- as.numeric(quantile(unlist(expression_nonmut), na.rm = TRUE)[4])
  threshold_mut <- as.numeric(q3_mut + (1.5 * iqr_mut))
  print(paste("threshold mut:", threshold_mut))
  threshold_nonmut <- as.numeric(q3_nonmut + (1.5 * iqr_nonmut))
  print(paste("threshold nonmut:", threshold_nonmut))
  # TODO: IF THERE ARE MULTIPLE ROWS, SHOULD I TAKE THE IQR SEPARATELY FOR EACH? 
  
  # Get patients that do not exceed threshold 
  if(length(ensg) > 1) {
    filt_express_mult <- function(expression, threshold) {
      filt_expression <- lapply(expression, function(x) {
        if(!any(x > threshold)) {return(x)}})
      filt_expression <- removeNULL(filt_expression)
      return(filt_expression)
    }
    filt_expression_mut <- filt_express_mult(expression_mut, threshold_mut)
    filt_expression_nonmut <- filt_express_mult(expression_nonmut, threshold_nonmut)
  } else {
    filt_expression_mut <- expression_mut[expression_mut <= threshold_mut]
    filt_expression_nonmut <- expression_nonmut[expression_nonmut <= threshold_nonmut]
  }
  patients_to_keep <- c(names(filt_expression_mut), names(filt_expression_nonmut))
  
  # Filter the starter DF based on these patients
  starter_df_filt <- starter_df[rownames(starter_df) %fin% patients_to_keep,]
  
  # Return the filtered starter DF
  print(starter_df_filt)
  print(paste("Nrow original starter", nrow(starter_df)))
  print(paste("Nrow filtered starter", nrow(starter_df_filt)))
  return(starter_df_filt)
}
############################################################

#' Filters a mutation target DF in order to keep only rows that overlap the 
#' gene targets of interest
#' @param mutation_targ_df a mutation target DF to be subsetted
#' @param uniprot_ids a list of uniprot IDs to keep in the DF
filter_mut_by_uniprot <- function(mutation_targ_df, uniprot_ids) {
  mutation_rows_to_keep <- unique(unlist(lapply(uniprot_ids, function(x) grep(x, mutation_targ_df$swissprot))))
  mutation_targ_df <- mutation_targ_df[mutation_rows_to_keep,]
  print("Mutation Targ DF filt by Uniprot ID:")
  print(mutation_targ_df)
  return(mutation_targ_df)
}
############################################################

#' Filters a CNA DF in order to keep only rows that overlap the 
#' gene targets of interest
#' @param cna_df a copy number alteration dataframe to be subsetted
#' @param ensg_ids a list of ensembl IDs to keep in the DF
filter_cna_by_ensg <- function(cna_df, ensg_ids) {
  cna_rows_to_keep <- unique(unlist(lapply(ensg_ids, function(x) grep(x, rownames(cna_df)))))
  cna_df <- cna_df[cna_rows_to_keep,]
  print("CNA filt by ENSG:")
  print(cna_df)
  return(cna_df)
}
############################################################

#' Filters a methylation DF in order to keep only rows that overlap the 
#' gene targets of interest
#' @param methylation_df a methylation dataframe to be subsetted
#' @param ensg_ids a list of ensembl IDs to keep in the DF
filter_meth_by_ensg <- function(methylation_df, ensg_ids) {
  meth_rows_to_keep <- unique(unlist(lapply(ensg_ids, function(x) grep(x, methylation_df$ensg_ids))))
  methylation_df <- methylation_df[meth_rows_to_keep,]
  print("Meth filt by ENSG:") 
  print(methylation_df)
  return(methylation_df)
}
############################################################

#' Given a set of columns from an expression dataframe that align with 
#' a particular patient j and target gene k, calculates and returns the fold 
#' change in expression
#' @param patient_exp_cols the result of subsetting an expression DF to a particular 
#' patient of interest j
handle_exp_cases <- function(patient_exp_cols) {
  fc_exp <- NA
  
  # If we have no patient matches or only 1 unmatched result, set to NA (no tumor-normal matched data)
  if(length(patient_exp_cols) == 0 | length(patient_exp_cols) == 1) {fc_exp <- NA}
  
  # If it's not a dataframe, this was also a sign we had only one column - set to NA (no tumor-normal matched data)
  else if(!is.data.frame(patient_exp_cols)) {fc_exp <- NA}
  
  # If there are 2 columns (one for tumor, one for normal)
  else if(length(patient_exp_cols) == 2) {
    
    # If there's only 1 row (one ENSG match)
    if(nrow(patient_exp_cols) == 1) {
      # Extract the expression value for this gene in this patient for both tumor and normal
      exp_val_tumor <- as.numeric(patient_exp_cols[,grepl("-0", colnames(patient_exp_cols))])
      exp_val_norm <- as.numeric(patient_exp_cols[,grepl("-11", colnames(patient_exp_cols))])
      # Calculate the fold change
      fc_exp <- exp_val_tumor / exp_val_norm
      
      # If there's >1 row (multiple ENSG matches)
    } else if (nrow(patient_exp_cols) > 1) {
      # Get the mean of the values across these different ENSG hits
      exp_val_tumor <- mean.default(as.numeric(unlist(patient_exp_cols[,grepl("-0", colnames(patient_exp_cols))])))
      exp_val_norm <- mean.default(as.numeric(unlist(patient_exp_cols[,grepl("-11", colnames(patient_exp_cols))])))
      # Calculate the fold change
      fc_exp <- exp_val_tumor / exp_val_norm
    }
    else {log_fc_exp <- NA}
    
  }
  # If we have more than one tumor or normal sample; take the average across them
  else if(length(patient_exp_cols) > 2) {
    # Get the mean of the values across these different patient samples
    tumor_exp_vals <- c()
    normal_exp_vals <- c()
    for (i in 1:ncol(patient_exp_cols)) {
      if(grepl("-0", colnames(patient_exp_cols)[i])) {
        tumor_exp_vals <- c(tumor_exp_vals, as.numeric(patient_exp_cols[,i]))
      } else if (grepl("-11", colnames(patient_exp_cols)[i])) {
        normal_exp_vals <- c(normal_exp_vals, as.numeric(patient_exp_cols[,i]))
      }
    }
    exp_val_tumor <- mean.default(as.numeric(tumor_exp_vals))
    exp_val_norm <- mean.default(as.numeric(normal_exp_vals))
    # Calculate the fold change
    fc_exp <- exp_val_tumor / exp_val_norm
  } else {print(paste("Unhandled case:", patient))}
  
  print(paste("Fold change:", fc_exp))
  return(fc_exp)
}
############################################################

#' Given an analysis type ('eQTL' or 'meQTL'), a methylation DF, and a patient j,
#' this function subsets this dataframe to the given patient j and returns either the
#' raw methylation value or the differential methylation value
#' @param methylation_df a methylation DF, subsetted to a particular target gene k 
#' and only cancer samples
#' @param anaylsis_type a string (either 'eQTL' or 'meQTL') that denotes the analysis type
#' @param patient the four-character patient ID for given patient j
handle_methylation_cases <- function(methylation_df, analysis_type, patient) {
  
  ### eQTL ###
  if (analysis_type == "eQTL") {
    
    # Extract the methylation value for this patient in this gene
    meth_val <- methylation_df[,colnames(methylation_df) == patient]
    
    # If there is more than one, take the average
    if (length(meth_val) > 1) {meth_val <- mean.default(meth_val)}
    
    # As long as this methylation value is not NA, either take the raw methylation value
    # or calculate whether the gene is differentially methylated
    if(is.na(meth_val)) {meth_stat <- NA} 
    else {
      if (length(meth_val) == 0) {meth_stat <- NA}
      # OPT 1. (RAW METHYLATION): Does this target gene have a methylation marker in cancer?
      meth_stat <- meth_val
      # OPT 2. (DIFFERENTIAL METHYLATION): Does this target gene have a differential methylation
      # state between tumor and normal?
      #else {
      #meth_stat <- 0
      #meth_beta_thres <- 0.25
      #if (meth_val > meth_beta_thres | meth_val < -meth_beta_thres) {meth_stat <- 1}
      #}
    }
    # Return the methylation stat
    return(meth_stat)
  }
  
  ### meQTL ###
  else if (analysis_type == "meQTL") {
    # If the log fold change in expression is greater than or less than a given 
    # threshold, give 1; otherwise 0
    logfc_exp_binary <- 0
    exp_thres <- 0.5
    if (!is.na(log_fc_exp)) {
      if (log_fc_exp > exp_thres | log_fc_exp < exp_thres) {logfc_exp_binary <- 1}
    } else {logfc_exp_binary <- NA}
    
    # Get the log fold change in methylation
    meth_cols <- methylation_df[,grepl(patient, colnames(methylation_df))]
    meth_val_tumor <- meth_cols[,grepl("norm", colnames(meth_cols))]
    meth_val_norm <- meth_cols[,grepl("tumor", colnames(meth_cols))]
    log_fc_meth <- log(meth_val_tumor / meth_val_norm)
    
    return(log_fc_meth)
  }
}



############################################################
############################################################
### VISUALIZATION FUNCTIONS
############################################################
############################################################

############################################################
#### VISUALIZE BETA VALUE DISTRIBUTION
############################################################
#' Function plots a histogram to visualize the distribution
#' of regulatory protein r_i beta values across all linear models.
#' @param results_table a master DF produced from run_linear_model()
visualize_beta_distrib <- function(results_table) {
  betas <- results_table$estimate
  hist(betas, main = "Histogram of Beta Coefficient Values Across all Reg. Proteins",
       xlab = "Beta Coefficient Value", ylab = "Frequency")
}

# Call this function
visualize_beta_distrib(master_df)


############################################################
#### VISUALIZE P-VALUE DISTRIBUTION
############################################################
#' Function plots a histogram to visualize the distribution
#' of regulatory protein r_i beta values across all linear models.
#' @param results_table a master DF produced from run_linear_model()
visualize_pval_distrib <- function(results_table) {
  pvals <- results_table$p.value[!is.na(results_table$p.value) & !is.infinite(results_table$p.value)]
  hist(pvals, main = "Histogram of p-Values Across all Reg. Proteins",
       xlab = "p-value", ylab = "Frequency")
}

# Call this function
visualize_pval_distrib(master_df)

############################################################
#### VISUALIZE Q-Q PLOT
############################################################
#' Function plots a Q-Q plot to visualize the distribution of p-values
#' and assess whether they come from a uniform distribution
#' @param results_table a master DF produced from run_linear_model()
qqplot_pvals <- function(results_table) {
  qqnorm(results_table$p.value, pch = 1, frame = FALSE)
  qqline(results_table$p.value, col = "steelblue", lwd = 2)
}

# Call this function
qqplot_pvals(master_df)

############################################################
#### VISUALIZE ERROR DISTRIBUTION
############################################################
#' Function plots a histogram of the standard errors produced from
#' all LM runs to assess whether the errors derive from a normal distribution
#' @param results_table a master DF produced from run_linear_model()
visualize_error_distrib <- function(results_table) {
  hist(results_table$std.error, main = "Standard Error Distribution Across All Tests",
       xlab = "Standard Error (SE)", ylab = "Frequency")
}

# Call this function
visualize_error_distrib(master_df)

############################################################
############################################################
#### PEFORM MULTIPLE HYPOTHESIS TESTING CORRECTION
############################################################
############################################################
#' Function takes in an output results table and applies multiple
#' hypothesis testing correction (Storey's q-value correction) to 
#' all p-values in order to add a column of q-values. Returns 
#' the results table with a column for q-values.
#' @param results_table a master DF produced from run_linear_model()
mh_correct <- function(results_table) {
  
  # Get the qvalue object
  qobj <- qvalue(p = results_table$p.value)
  
  # OPT: plot some useful plots & print some useful information
  plot(qobj)
  print(summary(qobj))
  #print(paste("Pi0 (Propr. of true null hypotheses):", qobj$pi0))
  
  qvals <- qobj$qvalues # extract qvalues
  
  results_table$q.value <- qvals # add the qvalues back to the dataframe
  
  # Return the tidied linear model fit with q-values
  return(results_table)
}

# Call this function
master_df_corrected <- mh_correct(master_df)

# Write this to a new file
write.csv(master_df_corrected, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/iprotein_output_results_df_TP53_TMM_corrected.csv", sep = ""))


############################################################
############################################################
#### OBTAIN SIGNIFICANT CORRELATIONS
############################################################
############################################################
#' Function takes in an output results table with q-values and 
#' restricts it to only models that exceed the q-value threshold
#' (are statistically significant correlations). Then ranks the 
#' remaining by q-values and returns the ranked top hits list.
#' @param results_table a master DF produced from run_linear_model() that has q-values added
#' from the mh_correct() function
#' @param qval_thres a threshold for significance for q-values
get_signif_correl <- function(results_table, qval_thres) {
  
  # Limit to only entries that exceed the given qvalue threshold
  results_table_sig <- results_table %>% filter(q.value < qval_thres)
  
  # Sort the table by qvalue
  results_table_sig_ordered <- results_table_sig[order(results_table_sig$q.value, decreasing = FALSE),]
  
  return(results_table_sig_ordered)
}

# Call this function
qval_thres <- 0.05
master_df_sig <- get_signif_correl(master_df_corrected, qval_thres)

# Add a column for target protein name
master_df_sig$T_k.name <- unlist(lapply(master_df_sig$T_k, function(x) paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == unlist(strsplit(x, ";", fixed = TRUE))[1], 'external_gene_name']), collapse = ";")))

# Write these results to a new file
write.csv(master_df_sig, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/iprotein_TP53_rawCNAi_TMM_significant_output_df.csv", sep = ""))                                                                                                                     # 118 (including CNAi and Methi as covar.)

############################################################
############################################################
#### VISUALIZE ENRICHMENT OF TOP HITS IN CGC/ VOGELSTEIN/ etc.
############################################################
############################################################
#' Function to plot the enrichment of cancer genes among the significant target genes
#' @param master_df_sig a master DF produced from run_linear_model() that has q-values
#' and has been thresholded to only those pairings that exceed a significance threshold
#' @param known_cancer_genes_table a data frame containing a compiled list of known cancer genes
#' (CGC, Vogelstein, etc.) with various types of IDs
plot_cgc_enrichment <- function(master_df_sig, known_cancer_genes_table) {
  # Get all the significant target genes 
  significant_target_hits <- master_df_sig$T_k.name
  
  # Fill in a vector of 0 and 1 for each significant target hit to indicate if it is a known cancer gene
  cancer_vect <- unlist(lapply(significant_target_hits, function(x) ifelse(x %fin% known_cancer_genes_table$primary_gene_names, 1, 0)))
  
  # Plot the enrichment (fraction of CGC genes at/ above given rank)
  ranks <- 1:length(significant_target_hits)
  frac_cancer_vect <- c()
  count_of_cancer_genes <- 0
  for (i in 1:length(significant_target_hits)) {
    val_of_curr_tg <- cancer_vect[i]
    count_of_cancer_genes <- count_of_cancer_genes + val_of_curr_tg
    frac <- count_of_cancer_genes / i
    frac_cancer_vect <- c(frac_cancer_vect, frac)
  }
  
  plot(ranks, frac_cancer_vect, main = "Enrichment of Significant Target Genes in Known Cancer Genes",
       xlab = "Rank", ylab = "Fraction of Target Genes that are Known Cancer Genes")
}

# Import a table containing a compiled list of known cancer genes
known_cancer_genes_table <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/GRCh38_driver_gene_list.tsv", sep = "\t",
                                       header = TRUE, check.names = FALSE, comment.char = "#")

# Run function
plot_cgc_enrichment(master_df_sig, known_cancer_genes_table)



create_genotype_file <- function(maf) {
  
  # Get all the unique SNP positions
  snp_positions <- unique(maf$Start_Position)
  
  # Get all the patient sample IDs
  samples <- unique(maf$Tumor_Sample_Barcode)
  
  # Loop through all the SNP positions and get the genotypes for each patient
  # at that position
  genotype_rows <- lapply(snp_positions, function(snp_pos) {
    # Subset the table to SNPs at this position
    maf_sub <- maf[maf$Start_Position == snp_pos,]
    
    # For each patient, if they are not in the DF, give them a 2, if they are,
    # determine if they have 1 or 0 copies of the reference allele
    genotype <- unlist(lapply(samples, function(samp) {
      if (samp %fin% maf_sub$Tumor_Sample_Barcode) {
        samp_row <- maf_sub[maf_sub$Tumor_Sample_Barcode == samp,]
        if(samp_row$Tumor_Seq_Allele1 == samp_row$Reference_Allele) {
          if(samp_row$Tumor_Seq_Allele2 == samp_row$Reference_Allele) {return(2)}
          else {return(1)}
        } else {
          if (samp_row$Tumor_Seq_Allele2 == samp_row$Reference_Allele) {return(1)}
          else{return(0)}
        }
      } else {return(2)}
    }))
    
    # Make the genotype into a string to return
    genotype_str <- paste(genotype, collapse = "")
    
    return(genotype_str)
    #return(as.data.frame(genotype))
  })
  
  # Recombine this list of genotype strings per SNP into a data frame
  #genotype_tab <- do.call(rbind, genotype_rows)
  
  return(genotype_rows)
}

# Import UUID-TCGA Barcode conversion document
uuid_barcode_conv <- read.csv(paste(path, "Input Data Files/BRCA Data/Somatic_Mut_Data/mut_uuid_barcode_conversion.csv", sep = ""))
#uuid_barcode_conv <- read.csv(paste(path, "Input Data Files/TCGA Data (ALL)/Somatic_Mut_Data/mut_uuid_barcode_conversion.csv", sep = ""))
# Alternatively:
# uuid_barcode_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/TCGA Exploratory Analysis/UUID_TCGABarcode_Mapping.csv", 
#header = TRUE, check.names = FALSE, row.names = 1)

uuid_barcode_conv <- uuid_barcode_conv[!duplicated(uuid_barcode_conv$tcga_barcode),]  # delete duplicate rows
uuid_barcode_conv <- uuid_barcode_conv[complete.cases(uuid_barcode_conv),]  # delete any rows with NAs

# Limit to only cancer samples for now
#methylation_df <- methylation_df[, grepl("-0", colnames(methylation_df))]

# Save this file again for future use
#write.csv(methylation_df, paste(main_path, "Methylation/methylation_DF_0.8_cancer_only.csv", sep = ""))

# 2. Differential methylation as a covariate
#differential_methylation_df <- read.csv(paste(main_path, "Methylation/Differential Methylation/differential_avg_beta_per_gene.csv", sep = ""), 
#header = TRUE, row.names = 1, check.names = FALSE)
# differential_methylation_df <- read.csv(paste(main_path, "Methylation/Differential Methylation/differential_num_probes_per_gene.csv", sep = ""), 
# header = TRUE, row.names = 1)
# Make the gene symbols the rownames
#rownames(differential_methylation_df) <- differential_methylation_df$Gene.Symbol

# 3. Differential methylation as a dependent variable
# methylation_df_dependent <- read.csv(paste(main_path, "Methylation/Differential Methylation/average_beta_per_gene.csv", sep = ""), header = TRUE, row.names = 1)
# methylation_df_dependent <- read.csv(paste(main_path, "Methylation/Differential Methylation/differential_num_probes_per_gene.csv", sep = ""), header = TRUE)



# 2. Differential methylation as a dependent variable (meQTL)
patients_with_matched_meQTL <- unique(unlist(lapply(2:length(colnames(methylation_df_dependent)), function(i) 
  unlist(strsplit(colnames(methylation_df_dependent[i]), ".", fixed = TRUE))[1])))
# Get the expression patient names that don't overlap this list
exp_patient_ids <- unlist(lapply(colnames(expression_df), function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))
patients_with_matched_meQTL <- intersect(patients_with_matched_meQTL, exp_patient_ids)
# 87 in BRCA


# Keep only cancer samples
# cna_df <- cna_df[,grepl("-0", colnames(cna_df))] # 374 BRCA samples for GISTIC, 747 for raw
# Remove the "." from the gene ID
#

# Remove sample ID
#colnames(cna_df) <- unlist(lapply(colnames(cna_df), function(x) 
#unlist(strsplit(x, ".", fixed = TRUE))[1]))
#rownames(cna_df) <- gene_ids

# Write this back to a CSV
# write.csv(cna_df, paste(main_path, "CNV/Gene-level Score/CNV_DF_AllGenes_CancerOnly.csv", sep = ""))
#write.csv(cna_df, paste(main_path, "CNV/Gene-level Raw/CNV_DF_AllGenes_CancerOnly.csv", sep = ""))

#
#
#
# Read this back
#cna_df <- read.csv(paste(main_path, "CNV/Gene-level Score/CNV_DF_AllGenes_CancerOnly.csv", sep = ""), header = TRUE, row.names = 1)
cna_df <- read.csv(paste(main_path, "CNV/Gene-level Raw/CNV_DF_AllGenes_CancerOnly.csv", sep = ""), header = TRUE, row.names = 1)


for (i in 1:length(file_uuids)) {
  uuid <- file_uuids[i]
  
  tcga_barcode <- ""
  tryCatch(
    expr = {
      tcga_barcode <- UUIDtoBarcode(uuid, from_type = "case_id")
      uuids_and_barcodes[i,2] <- tcga_barcode[1,2]
    },
    error = function(e) {
      print(paste("Unable to convert UUID ", uuid))
    }
  )
}
# head(uuids_and_barcodes)


source_path <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/"
source(paste(source_path, "general_important_functions.R", sep = ""))


############################################################
# SET UP PARSER ARGUMENTS
############################################################
# Create parser object
parser <- ArgumentParser()

# Specify desired options for the parser

# Name of master DF mutation and CNA files
parser$add_argument("--master_df_mut", default = "output_results_P53_chipeat_iprotein_tmm_rawCNA_methBetaRaw_cibersortTotalFrac_uncorrected_MUT.csv", 
                    type = "character", help = "Name of mutation-specific results file. No path.")
parser$add_argument("--master_df_cna", default = "output_results_P53_chipeat_iprotein_tmm_rawCNA_methBetaRaw_cibersortTotalFrac_uncorrected_CNA.csv", 
                    type = "character", help = "Name of mutation-specific results file. No path.")

# Type of QTL
parser$add_argument("--QTLtype", default = "eQTL", type = "character", 
                    help = "eQTL or meQTL. Default is eQTL.")

# Whether we are using tumor-normal matched or non-matched data
parser$add_argument("--tumNormMatched", default = "FALSE", type = "character",
                    help = "TRUE/FALSE for whether or not we are looking at tumor-normal-matched data. Default is FALSE.")

# Cancer type we are looking like (or pan-cancer)
parser$add_argument("--cancerType", default = "BRCA", type = "character",
                    help = "Type of cancer we're looking at, or pan-cancer. Default is BRCA.")

# Randomization
parser$add_argument("--randomize", default = "FALSE", type = "character", 
                    help = "TRUE/FALSE for whether or not data is randomized. Default is FALSE")

# Whether or not we are running a test on one regulatory protein, and details of 
# the tester protein if so
parser$add_argument("--test", default = "FALSE", type = "character",
                    help = "TRUE/FALSE for whether or not we are running a test with one regulatory protein, vs. many. Default is FALSE.")
parser$add_argument("--tester_name", default = "P53", type = "character",
                    help = "Name of regulatory protein if doing a 1-protein test. Defaults to TP53.")
parser$add_argument("--tester_uniprot_id", default = "P04637", type = "character",
                    help = "Uniprot ID of regulatory protein if doing a 1-protein test. Defaults to TP53 (P04637).")
parser$add_argument("--tester_ensg_id", default = "ENSG00000141510", type = "character",
                    help = "ENSG ID of regulatory protein if doing a 1-protein test. Defaults to TP53 (ENSG00000141510).")

# The decision for if/ how to bucket CNAs and methylation
parser$add_argument("--cna_bucketing", default = "bucket_inclAmp", 
                    type = "character",
                    help = "If/ the type of bucketing we use for CNA values. Default is bucket_inclAmp, but rawCNA and bucket_exclAmp are also options.")
parser$add_argument("--meth_bucketing", default = TRUE, type = "character",
                    help = "If/ the type of bucketing we use for methylation values. Default is FALSE.")

# Label the type of methylation data we are using (Beta, M, or Threshold)
parser$add_argument("--meth_type", default = "Beta", type = "character",
                    help = "The type methylation values we are using. Default is Beta, but other options are M and Threshold_X, where X is the threshold value.")

# What type of test we are running, so we can write the results files to an appropriate directory. 
# Only necessary if --test is F.
parser$add_argument("--run_name", default = "cancer_related_genes", type = "character",
                    help = "Provide a name for the run, if not testing on just one regulatory protein, to write output files to appropriate directory. [default %(default)s]") 

# A name for the group of targets being tested, for labeling the output file.
parser$add_argument("--targets_name", default = "allGenes", type = "character",
                    help = "Provide a name for the group of targets being tested, in order to properly name the output file. [default %(default)s]")

# Whether or not we are including PEER factors and PCs as covariates in the model
parser$add_argument("--include_PEER", default = "TRUE", type = "character",
                    help = "A TRUE/ FALSE value indicating whether or not we are including PEER factors as covariates in our model. Default is TRUE.")
parser$add_argument("--include_pcs", default = "TRUE", type = "character",
                    help = "A TRUE/ FALSE value indicating whether or not we are including principal components as covariates in our model. Default is TRUE.")

# Parse the given arguments
args <- parser$parse_args()


############################################################
### CONVERT STRING LOGICALS TO REAL LOGICALS
############################################################
#' Add a string to bool function which converts a string "boolean" into a boolean type
#' Modified from: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
#' @param v a argparse input value that should be a logical
str2bool <- function(v) {
  if (is.logical(v)) {return(v)}
  if (tolower(v) %in% c('yes', 'true', 't', 'y', '1')) {return(TRUE)}
  else if (tolower(v) %in% c('no', 'false', 'f', 'n', '0')) {return(FALSE)}
  else {
    print("Error. Boolean value expected.")
    return(9)
  }
}

tumNormMatched <- str2bool(args$tumNormMatched)
randomize <- str2bool(args$randomize)
test <- str2bool(args$test)
meth_bucketing <- str2bool(args$meth_bucketing)
include_PEER <- str2bool(args$include_PEER)
include_pcs <- str2bool(args$include_pcs)


# Path to output files
main_path <- paste(source_path, paste("output_files/", args$cancerType, sep = ""), sep = "")
if(tumNormMatched) {main_path <- paste(main_path, "tumor_normal_matched", sep = "")}
else {main_path <- paste(main_path, "tumor_only", sep = "/")}


