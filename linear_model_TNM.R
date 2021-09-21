############################################################
### Linear Model (Tumor-Normal Matched)
### Written By: Sara Camilli, July 2020
############################################################
library(biomaRt)
library(parallel)
library(rlang)
library(dplyr)
library(rlist)
library(rockchalk)
#library(enrichvs)
library(broom)

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

protein_ids_df <- read.csv(paste(prot_path, "iprotein_protein_ids_df", sep = ""), 
                           header = TRUE, check.names = FALSE)
# protein_ids_df <- read.csv(paste(prot_path, "iprotein_nucacids_protein_ids_df", sep = ""), 
                           # header = TRUE, check.names = FALSE)
# protein_ids_df <- read.csv(paste(prot_path, "idomain_protein_ids_df", sep = ""), 
                           # header = TRUE, check.names = FALSE)
# protein_ids_df <- read.csv(paste(prot_path, "idomain_nucacids_protein_ids_df", sep = ""), 
                           # header = TRUE, check.names = FALSE)
# protein_ids_df <- read.csv(paste(prot_path, "ibindingpos_protein_ids_df", sep = ""), 
                            # header = TRUE, check.names = FALSE)
# protein_ids_df <- read.csv(paste(prot_path, "ibindingpos_nucacids_protein_ids_df", sep = ""), 
                            # header = TRUE, check.names = FALSE)


############################################################
# IMPORT EXPRESSION FILES
############################################################
# 1. Unfiltered FPKM
expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_DF_tumNormMatched.csv", sep = ""), 
                          header = TRUE, row.names = 1, check.names = FALSE)
expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_DF_tumNormMatched_meQTL.csv", sep = ""), 
                          header = TRUE, row.names = 1, check.names = FALSE)

# 2. Filtered FPKM
expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_filt_DF_tumNormMatched.csv", sep = ""), 
                          header = TRUE, row.names = 1, check.names = FALSE)
expression_df <- read.csv(paste(main_path, "Expression/expression_fpkm_filt_DF_tumNormMatched_meQTL.csv", sep = ""), 
                          header = TRUE, row.names = 1, check.names = FALSE)

# 3. TMM
expression_df <- read.csv(paste(main_path, "Expression/expression_tmm_DF_tumNormMatched.csv", sep = ""), 
                          header = TRUE, row.names = 1, check.names = FALSE)
expression_df <- read.csv(paste(main_path, "Expression/expression_tmm_DF_tumNormMatched_meQTL.csv", sep = ""), 
                          header = TRUE, row.names = 1, check.names = FALSE)

# 4. Quantile-Normalized
# TODO: Fill

# 5. Rank-Normalized
# TODO: fill

############################################################
# IMPORT METHYLATION FILES
############################################################
methylation_df <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_tumNormMatched.csv", sep = ""), 
                           header = TRUE, row.names = 1, check.names = FALSE)
methylation_df <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_tumNormMatched_eQTL.csv", sep = ""), 
                           header = TRUE, row.names = 1, check.names = FALSE)
methylation_df_dependent <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_tumNormMatched_dependent.csv", sep = ""), 
                                     header = TRUE, row.names = 1, check.names = FALSE)
methylation_df_dependent <- read.csv(paste(main_path, "Methylation/methylation_DF_0.8_tumNormMatched_dependent_eQTL.csv", sep = ""), 
                                     header = TRUE, row.names = 1, check.names = FALSE)


############################################################
# IMPORT GENE TARGET MUTATION FILES
############################################################
# For eQTL
mutation_targ_df <- read.csv(paste(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense_tumNormMatched.csv", sep = ""), 
                             header = TRUE, row.names = 1, check.names = FALSE)
# For meQTL
mutation_targ_df <- read.csv(paste(main_path, "Mutation/Mutation Count Matrices/mut_count_matrix_missense_tumNormMatched_meQTL.csv", sep = ""), 
                             header = TRUE, row.names = 1, check.names = FALSE)


############################################################
# IMPORT REGULATORY PROTEIN MUTATION FILES
############################################################
# For eQTL
mutation_regprot_df <- read.csv(paste(main_path, "Mutation/iprotein_results_missense_tumNormMatched.csv", sep = ""), 
                                header = TRUE, row.names = 1, check.names = FALSE)
# For meQTL
mutation_regprot_df <- read.csv(paste(main_path, "Mutation/iprotein_results_missense_tumNormMatched_meQTL.csv", sep = ""), 
                                header = TRUE, row.names = 1, check.names = FALSE)


############################################################
# IMPORT CNA FILES
############################################################
# For eQTL
cna_df <- read.csv(paste(main_path, "CNV/Gene-level Raw/CNV_DF_AllGenes_CancerOnly_tumNormMatched.csv", sep = ""), 
                   header = TRUE, row.names = 1, check.names = FALSE)
# For meQTL
cna_df <- read.csv(paste(main_path, "CNV/Gene-level Raw/CNV_DF_AllGenes_CancerOnly_tumNormMatched_meQTL.csv", sep = ""), 
                   header = TRUE, row.names = 1, check.names = FALSE)

# If rownames are unlabeled
# cna_df_full <- read.csv(paste(main_path, "CNV/Gene-level Score/CNV_DF_AllGenes.csv", sep = ""), header = TRUE, row.names = 1)
# gene_ids <- rownames(cna_df_full)
# rownames(cna_df) <- gene_ids
# rm(cna_df_full)


############################################################
# IMPORT PATIENT/SAMPLE FILES
############################################################
# For eQTL
patient_df <- read.csv(paste(main_path, "Linear Model/combined_patient_sample_tumNormMatched.csv", sep = ""), 
                       header = TRUE, row.names = 1, check.names = FALSE)
# For meQTL
patient_df <- read.csv(paste(main_path, "Linear Model/combined_patient_sample_tumNormMatched_eQTL.csv", sep = ""), 
                       header = TRUE, row.names =  1, check.names = FALSE)


#########################################################
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
sample_targets_DF <- read.csv(paste(main_path, "Linear Model/TP53/tp53_curated_targets.csv", sep = ""), 
                              row.names = 1, check.names = FALSE)
sample_targets_DF <- read.csv(paste(main_path, "Linear Model/FOXA1/foxa1_curated_targets.csv", sep = ""), 
                              row.names = 1, check.names = FALSE)

# Option 2: ChIP-eat Targets
sample_targets_DF <- read.csv(paste(main_path, "Linear Model/TP53/tp53_chipeat_targets.csv", sep = ""), 
                              row.names = 1, check.names = FALSE)
sample_targets_DF <- read.csv(paste(main_path, "Linear Model/FOXA1/foxa1_chipeat_targets.csv", sep = ""), 
                              row.names = 1, check.names = FALSE)

# Option 3: All Gene Targets
sample_targets_DF <- read.csv(paste(main_path, "Linear Model/allgene_targets.csv", sep = ""), 
                              row.names = 1, check.names = FALSE)


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
  ### PATIENT/ SAMPLE CHARACTERISTICS ###
  # Gender : Gender of all patients 1..j..M (0 if male, 1 if female) -- ONLY INCLUDE IF NOT BRCA
  # Age : Age of all patients 1..j..M in years 
    # Buckets: 1 if 0-9, 2 if 10-19, 3 if 20-29, 4 if 30-39, 5 if 40-49, 6 if 50-59, 7 if 60-69, 8 if 70-79, 9 if 80-89, 10 if 90+
  # Race : Race/ Ethnicity of all patients 1..j..M 
    # Buckets: 1 if White (not Latinx), 2 if White (Latinx), 3 if black or African, 4 if Asian, 5 if other
  # Prior_malig : If all patients 1..j..M had a prior malignancy (0 if no, 1 if yes)
  # Treatment_rad : If all patients 1..j..M were treated with radiation (0 is not treated, 1 is treated) 
  # Treatment_pharm : If all patients 1..j..M were treated with pharmaceutical therapy (0 is not treated, 1 is treated) 
  # TotalNumMut : Total number of missense mutations each patient j has (across all patients 1..j..M), integer value
  # Tumor_purity : An estimate of the 'purity' of the tumor (fraction of tumor cells within tumor sample)
  # Total_IC_Frac : An estimate of the fraction of tumor cells in the sample, from CIBERSORT Abs

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
run_linear_model <- function(protein_ids_df, downstream_target_df, patient_df, mutation_df_targ, mutation_df_regprot, 
                             methylation_df, cna_df, expression_df, analysis_type) {
  
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
                                                   cna_df, expression_df, analysis_type)
      
      if(length(linear_model_input_table) == 0) {return(NA)}
      if (is.na(linear_model_input_table)) {return(NA)}
      
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
        #lm_fit <- lm(formula = (LogFCExp_k ~ MutStat_i + MutStat_k + MethStat_k + CNAStat_k + Age_b3 + 
                                  #Age_b4 + Age_b5 + Age_b6 + Age_b7 + Age_b8 + Age_b9 + Prior_malig + 
                                  #Treatment_rad + Treatment_pharm + Tot_Mut_b1 + Tot_Mut_b2), 
                     #data = linear_model_input_table_noNA)
        # Using CNAStat_k as an offset rather than as a covariate, adding library size as a covariate
        lm_fit <- lm(formula = (LogFCExp_k ~ MutStat_i + CNAStat_i + LogFCMeth_i + MutStat_k + Age_b3 + 
                                  Age_b4 + Age_b5 + Age_b6 + Age_b7 + Age_b8 + Age_b9 + Prior_malig + 
                                  Treatment_rad + Treatment_pharm + Tot_Mut_b1 + Tot_Mut_b2 +
                                  Tumor_purity + Total_Frac_ICD), 
                     offset = log(CNAStat_k * Lib_Size),
                     data = linear_model_input_table_final)
        
      } else if (analysis_type == "meQTL") {
        lm_fit <- lm(formula = (LogFCMeth_k ~ MutStat_i + CNAStat_i + MethStat_i + MutStat_k + LogFCExp_k 
                                + Age_b3 + Age_b4 + Age_b5 + Age_b6 + Age_b7 + Age_b8 + Age_b9 + Prior_malig + 
                                  Treatment_rad + Treatment_pharm + Tot_Mut_b1 + Tot_Mut_b2 + 
                                  Tumor_purity + Total_Frac_ICD), 
                     offset = log(CNAStat_k * Lib_Size),
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
                              methylation_df, cna_df, expression_df, "eQTL")
master_df <- run_linear_model(protein_ids_df_tester, sample_targets_DF, patient_df, mutation_targ_df, mutation_regprot_df,
                              methylation_df_dependent, cna_df, expression_df, "meQTL")


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
  
  # Filter the data frames to look at only this regulatory protein
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
fill_targ_inputs <- function(starter_df, targ_k, targ_k_ensg, mutation_targ_df,
                                  methylation_df, cna_df, expression_df, analysis_type) {
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
    exp_rows_to_keep <- unique(unlist(lapply(targ_k_ensg, function(x) 
      grep(x, rownames(expression_df)))))
    expression_df <- expression_df[exp_rows_to_keep,]
    print("Expression DF filt by ENSG ID")
    print(expression_df)
    
    # Continue only if all of these dataframes contain information about this target gene k
    if(!(nrow(expression_df) == 0 | nrow(cna_df) == 0 | nrow(methylation_df) == 0 | 
         nrow(mutation_targ_df) == 0)) {
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
        exp_stat <- log(handle_exp_cases(patient_exp_cols))
        
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
      
      if (analysis_type == "eQTL") {colnames(targ_k_df) <- c("LogFCExp_k", "MutStat_k", 
                                                             "CNAStat_k", "MethStat_k")}
      else {colnames(targ_k_df) <- c("ExpStat_k", "MutStat_k", "CNAStat_k", "LogFCMeth_k")}
      print(head(targ_k_df))
      full_df <- cbind(starter_df, targ_k_df)
      
      # Return this full input DF
      return(full_df)
    }
  }
  return(NA)
}







