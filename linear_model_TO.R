############################################################
### Linear Model (Tumor-Only): ARCHIVE
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
library(data.table)
library(speedglm)
#library(Rcpp)
#library(RcppEigen)

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)

# Source other files needed
source_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/"
source(paste(source_path, "linear_model_helper_functions.R", sep = ""))
source(paste(source_path, "general_important_functions.R", sep = ""))

############################################################
# GET ALL REGULATORY PROTEINS IN CONSIDERATION
############################################################
# Regulatory proteins will vary depending on the level of specificity in question;
# import only the proteins at the given specificity level of interest under label "protein_ids_df"
prot_path <- paste(main_path, "Mutation/Files for Linear Model/", sep = "")

#protein_ids_df <- fread(paste(prot_path, "iprotein_protein_ids_df.csv", sep = ""), 
                           #header = TRUE)
# protein_ids_df <- fread(paste(prot_path, "iprotein_nucacids_protein_ids_df.csv", sep = ""), 
                            #header = TRUE)
# protein_ids_df <- fread(paste(prot_path, "idomain_protein_ids_df.csv", sep = ""), 
                            #header = TRUE)
# protein_ids_df <- fread(paste(prot_path, "idomain_nucacids_protein_ids_df.csv", sep = ""), 
                            #header = TRUE)
# protein_ids_df <- fread(paste(prot_path, "ibindingpos_protein_ids_df.csv", sep = ""), 
                            #header = TRUE)
# protein_ids_df <- fread(paste(prot_path, "ibindingpos_nucacids_protein_ids_df.csv", sep = ""), 
                            #header = TRUE)

# Import only known cancer-related genes data frame
protein_ids_df <- fread(paste(prot_path, "iprotein_protein_ids_df_cancerRelated.csv", sep = ""),
                        header = TRUE)

# Import only highly deleted TFs
protein_ids_df <- fread(paste(prot_path, "iprotein_protein_ids_df_highlyDelTFs.csv", sep = ""),
                        header = TRUE)


############################################################
# IMPORT EXPRESSION FILES
############################################################
# 1. FPKM (filtered)
expression_df <- fread(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_fpkm_IntersectPatients.csv", sep = ""), 
                          header = TRUE)  # fpkm-normalized and filtered

# 2. TMM 
expression_df <- fread(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_tmm_IntersectPatients.csv", sep = ""), 
                          header = TRUE)

# 3. Quantile-Normalized
expression_df <- fread(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_quantile_norm_IntersectPatients.csv", sep = ""),
                          header = TRUE)

# 4. Rank-Normalized
expression_df <- fread(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_rank_norm_IntersectPatients.csv", sep = ""),
                          header = TRUE)

colnames(expression_df)[1] <- "ensg_id"
# Set the key
#setkey(expression_df, ensg_id)


############################################################
# IMPORT METHYLATION FILES
############################################################
methylation_df <- fread(paste(main_path, "Linear Model/Tumor_Only/Methylation/methylation_0.8_CancerOnly_IntersectPatients.csv", 
                                 sep = ""), header = TRUE)
colnames(methylation_df)[1] <- "gene_name"
# Set the key
#setkey(methylation_df, ensg_id)


############################################################
# IMPORT GENE TARGET MUTATION FILES
############################################################
mutation_targ_df <- fread(paste(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_IntersectPatients.csv", 
                                   sep = ""), header = TRUE)
colnames(mutation_targ_df)[1] <- "gene_name"
mutation_targ_df <- mutation_targ_df[, Swissprot := unlist(lapply(mutation_targ_df$gene_name, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'uniprot_gn_id']), collapse = ";")))]
# Set the key
#setkey(mutation_targ_df, Swissprot)


############################################################
# IMPORT REGULATORY PROTEIN MUTATION FILES
############################################################
mutation_regprot_df <- fread(paste(main_path, "Linear Model/Tumor_Only/Regprot_Mutation/iprotein_results_missense_IntersectPatients.csv", 
                                      sep = ""), header = TRUE)
mutation_regprot_df <- fread(paste(main_path, "Linear Model/Tumor_Only/Regprot_Mutation/iprotein_results_missense_nucacids_IntersectPatients.csv", 
                                      sep = ""), header = TRUE)
mutation_regprot_df <- fread(paste(main_path, "Linear Model/Tumor_Only/Regprot_Mutation/idomain_results_missense_IntersectPatients.csv", 
                                      sep = ""), header = TRUE)
mutation_regprot_df <- fread(paste(main_path, "Linear Model/Tumor_Only/Regprot_Mutation/idomain_results_missense_nucacids_IntersectPatients.csv", 
                                      sep = ""), header = TRUE)
mutation_regprot_df <- fread(paste(main_path, "Linear Model/Tumor_Only/Regprot_Mutation/ibindingpos_results_missense_IntersectPatients.csv", 
                                      sep = ""), header = TRUE)
mutation_regprot_df <- fread(paste(main_path, "Linear Model/Tumor_Only/Regprot_Mutation/ibindingpos_results_missense_nucacids_IntersectPatients.csv", 
                                      sep = ""), header = TRUE)

#mutation_regprot_df$Swissprot <- unlist(lapply(mutation_regprot_df$Query, function(x) 
  #unlist(strsplit(x, "|", fixed = TRUE))[2]))
mutation_regprot_df <- mutation_regprot_df[,5:ncol(mutation_regprot_df)]   # remove first few meaningless columns, if necessary
# Set the key
#setkey(mutation_regprot_df, Swissprot)

############################################################
# IMPORT CNA FILES
############################################################
cna_df <- fread(paste(main_path, "Linear Model/Tumor_Only/CNV/CNA_AllGenes_CancerOnly_IntersectPatients.csv", 
                         sep = ""), header = TRUE)
colnames(cna_df)[1] <- "ensg_id"
#setkey(cna_df, ensg_id)

############################################################
# IMPORT PATIENT/SAMPLE FILES
############################################################
patient_df <- fread(paste(main_path, "Linear Model/Tumor_Only/Patient/combined_patient_sample_IntersectPatients.csv", sep = ""), 
                       header = TRUE)
colnames(patient_df)[1] <- "sample_id"
#setkey(patient_df, sample_id)

############################################################
# SET UP TARGETS FOR PROOF-OF-CONCEPT TESTING
############################################################
sample_protein_uniprot <- "P04637" 
sample_protein_ensg <- "ENSG00000141510"
sample_protein_name <- "P53"

#sample_protein_uniprot <- "P55317" 
#sample_protein_ensg <- "ENSG00000129514"
#sample_protein_name <- "FOXA1"

# A negative control
#sample_protein_uniprot <- "P36956"
#sample_protein_ensg <- "ENSG00000072310"
#sample_protein_name <- "SRBP1"



protein_ids_df_tester <- data.table(swissprot_ids = sample_protein_uniprot, 
                                    ensg_ids = sample_protein_ensg)

# Get sample targets data frame
# Option 1: Curated Targets
#sample_targets_DF <- fread(paste(main_path, "Linear Model/TP53/tp53_curated_targets.csv", sep = ""), 
                           #header = TRUE)
#sample_targets_DF <- fread(paste(main_path, "Linear Model/FOXA1/foxa1_curated_targets.csv", sep = ""), 
                           #header = TRUE)

# Option 2: ChIP-eat Targets
sample_targets_DF <- fread(paste(main_path, "Linear Model/TP53/tp53_chipeat_targets.csv", sep = ""), 
                           header = TRUE)
#sample_targets_DF <- fread(paste(main_path, "Linear Model/FOXA1/foxa1_chipeat_targets.csv", sep = ""), 
                           #header = TRUE)

# Option 3: All Gene Targets
#sample_targets_DF <- fread(paste(main_path, "Linear Model/allgene_targets.csv", sep = ""), 
                           #header = TRUE)

# Option 4: Metabolic Targets
sample_targets_DF <- fread(paste(main_path, "Linear Model/metabolic_targets.csv", sep = ""),
                           header = TRUE)

# Remove rows with any blank elements
sample_targets_DF <- sample_targets_DF[!(sample_targets_DF$swissprot == ""),]

# Remove the first column, if necessary
sample_targets_DF <- sample_targets_DF[,2:3]


# Regardless of method, limit targets to only those that overlap the genes in the expression DF
rows_to_keep <- unlist(lapply(1:length(sample_targets_DF$ensg), function(i) {
  targs <- unlist(strsplit(as.character(sample_targets_DF[i,'ensg']), ";", fixed = TRUE))
  if (any(targs %fin% expression_df$ensg_id)) {return(TRUE)}
  else {return(FALSE)}
}))
sample_targets_DF <- sample_targets_DF[rows_to_keep,]


############################################################
#### OVERVIEW OF LINEAR MODEL
############################################################
  ### PATIENT CHARACTERISTICS ###
  # Gender : Gender of all samples 1..j..M (0 if male, 1 if female) -- ONLY INCLUDE IF NOT BRCA
  # Age : Age of all samples 1..j..M in years 
    # Buckets: 1 if 0-9, 2 if 10-19, 3 if 20-29, 4 if 30-39, 5 if 40-49, 6 if 50-59, 7 if 60-69, 8 if 70-79, 9 if 80-89, 10 if 90+
  # Race : Race/ Ethnicity of all samples 1..j..M 
    # Buckets: 1 if White (not Latinx), 2 if White (Latinx), 3 if black or African, 4 if Asian, 5 if other
  # Prior_malig : If all samples 1..j..M had a prior malignancy (0 if no, 1 if yes)
  # Treatment_rad : If all samples 1..j..M were treated with radiation (0 is not treated, 1 is treated) 
  # Treatment_pharm : If all samples 1..j..M were treated with pharmaceutical therapy (0 is not treated, 1 is treated) 
  # TotalNumMut : Total number of missense mutations each patient j has (across all samples 1..j..M)
  # Tumor_purity : An estimate of the 'purity' of the tumor (fraction of tumor cells within tumor sample) 
    # Buckets: 1 if low purity (0-0.5), 2 if medium purity (0.5-0.75), 3 if high purity (0.75-1.0)
  # Total_IC_Frac : An estimate of the fraction of tumor cells in the sample, from CIBERSORT Abs
    # Buckets: 1 if low ICI (0-0.3), 2 if medium ICI (0.3-0.7), 3 if high ICI (0.7-1.0)
  
  ### REGULATORY PROTEIN i CHARACTERISTICS ###
  # MutStat_i : The mutation status of protein i (0 if not mutated, 1 if mutated) across all tumor samples 1..j..M
  # MethStat_i : The methylation status of protein i (0 if not methylated, 1 if methylated) across all tumor samples 1..j..M  
  # CNA_i : The copy number alteration status of protein i (1 if amplified, -1 if deleted, 0 if no amp. or del. OR copy #) 
      # across all tumor samples 1..j..M
  
  ### TARGET GENE k CHARACTERISTICS ###
  # MutStat_k : The mutation status of gene k (0 if not mutated, 1 if mutated) across all samples 1..j..M
  
  ### DEPENDENT VARIABLE ###
  # log(Exp_k,c) or log(Meth_k,c) : the cancer expression or methylation value 
  
  ### OFFSETS/ SCALING FACTORS ###
  # log(N_j) : log of effective library size (sample-specific scaling factor); the total number of reads for 
  # tumor sample j after cross-sample normalization, across all samples 1..j..M
      # NOTE: This is important only for TMM and FPKM expression, not for rank- or quantile-normalized expression
  # log(CNA_k,j) : the copy number for target gene k across all samples 1..j..M


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
run_linear_model <- function(protein_ids_df, downstream_target_df, patient_df, 
                             mutation_df_targ, mutation_df_regprot, methylation_df, 
                             cna_df, expression_df, analysis_type) {
  
  # We need to get a mini-table for every r_i, t_k combo that we rbind into a master table of results
  results_df_list <- mclapply(1:protein_ids_df[, .N], function(i) {
    
    # Get the given regulatory protein r_i's Swissprot & ENSG IDs
    regprot <- protein_ids_df$swissprot_ids[i]
    regprot_ensg <- unlist(strsplit(protein_ids_df$ensg_ids[i], ";", fixed = TRUE))
    
    print(paste("Regulatory protein", paste(i, paste("/", protein_ids_df[, .N]))))
    
    # Create a starter table that gets the mutation, CNA, and methylation status 
    # for this regulatory protein and binds it to the patient DF
    starter_df <- fill_regprot_inputs(patient_df, regprot, regprot_ensg, 
                                      mutation_df_regprot, methylation_df, 
                                      cna_df)
    
    ### OPTIONAL ###
    # Randomize the table's protein i mutations across all the patients as a test
    # starter_df$MutStat_i <- sample(starter_df$MutStat_i)
    
    #print(starter_df)
    
    # Loop through all the targets for this protein of interest and create a table 
    # for each of them from this starter DF
    regprot_i_results_df_list <- lapply(1:downstream_target_df[, .N], function(k) {
      
      # Get the target t_k's Swissprot & ENSG IDs
      targ <- unlist(strsplit(downstream_target_df$swissprot[k], ";", fixed = TRUE))
      targ_ensg <- unlist(strsplit(downstream_target_df$ensg[k], ";", fixed = TRUE))
      
      print(paste("Target gene", paste(k, paste("/", downstream_target_df[, .N]))))
      
      # OPT: Filter outlier expression for each group (mutated, unmutated) using a standard 
        # (1.5 x IQR) + Q3 schema
      #starter_df <- filter_expression_df(expression_df, starter_df, targ_ensg)
      
      # Create a full input table to the linear model for this target
      linear_model_input_table <- fill_targ_inputs(starter_df, targ, targ_ensg, 
                                                   mutation_df_targ, methylation_df,
                                                   cna_df, expression_df)
      
      if(length(linear_model_input_table) == 0) {return(NA)}
      if (is.na(linear_model_input_table)) {return(NA)}
      
      # If a row has NA, remove that whole row and predict without it
      #linear_model_input_table <- na.exclude(linear_model_input_table)
      linear_model_input_table <- na.omit(linear_model_input_table)
      #print(linear_model_input_table)
      
      if(linear_model_input_table[.N, ] == 0) {return(NA)}
      
      # Convert to matrix (for .lm.fit)
      #linear_model_input_table <- as.matrix(linear_model_input_table)[,,drop=FALSE]

      print(head(linear_model_input_table))
      
      ## OPT: randomize expression across cancer and normal, across all samples
      #linear_model_input_table$ExpStat_k <- sample(linear_model_input_table$ExpStat_k)
      
      # Run the linear model on this input table
      # Note: currently excluding age b1 and age b2 since no patients fall into these categories
      
      #terms <- c("MutStat_i", "CNAStat_i", "MethStat_i", "MutStat_k", "CNAStat_k", 
                 #"Age_b3", "Age_b4", "Age_b5", "Age_b6", "Age_b7","Age_b8", "Age_b9", 
                 #"Prior_malig","Treatment_rad", "Treatment_pharm", "Tot_Mut_b1", 
                 #"Tot_Mut_b2", "Tumor_purity_b1", "Tumor_purity_b2", "Tot_IC_Frac_b1",
                 #"Tot_IC_Frac_b2")
      
      if (analysis_type == "eQTL") {
        lm_fit <- speedglm::speedlm(formula = (log2(ExpStat_k) ~ MutStat_i + CNAStat_i + MethStat_i +
             MutStat_k + CNAStat_k + Age_b3 + Age_b4 + Age_b5 + Age_b6 + Age_b7 + 
               Age_b8 + Age_b9 + Prior_malig + Treatment_rad + Treatment_pharm + 
               Tot_Mut_b1 + Tot_Mut_b2 + Tumor_purity_b1 + Tumor_purity_b2 + 
               Tot_IC_Frac_b1 + Tot_IC_Frac_b2), 
            offset = log2(Lib_Size),
            data = linear_model_input_table)
        # Take the log of the expression values & library sizes
        #linear_model_input_table$ExpStat_k <- unlist(lapply(linear_model_input_table$ExpStat_k), log)
        #linear_model_input_table$Lib_Size <- unlist(lapply(linear_model_input_table$Lib_Size), log)
        #lm_fit <- .lm.fit(cbind(1, linear_model_input_table[,terms]),
                          #linear_model_input_table[,"ExpStat_k"],
                          #offset = linear_model_input_table[,"Lib_Size"])
        
      } else if (analysis_type == "meQTL") {
        lm_fit <- lm(formula = (log2(MethStat_k) ~ MutStat_i + CNAStat_i + MethStat_i +
                                  MutStat_k + CNAStat_k + Age_b3 + Age_b4 + Age_b5 + Age_b6 + Age_b7 + 
                                  Age_b8 + Age_b9 + Prior_malig + Treatment_rad + Treatment_pharm + 
                                  Tot_Mut_b1 + Tot_Mut_b2 + Tumor_purity_b1 + Tumor_purity_b2 + 
                                  Tot_IC_Frac_b1 + Tot_IC_Frac_b2), 
                     offset = log2(Lib_Size),
                     data = linear_model_input_table)
        #linear_model_input_table$MethStat_k <- unlist(lapply(linear_model_input_table$MethStat_k), log)
        #linear_model_input_table$Lib_Size <- unlist(lapply(linear_model_input_table$Lib_Size), log)
        #lm_fit <- .lm.fit(cbind(1, linear_model_input_table[,terms]),
                          #linear_model_input_table[,"MethStat_k"],
                          #offset = linear_model_input_table[,"Lib_Size"])
      } else {
        print(paste("Invalid analysis type:", analysis_type))
        return(NA)
      }
      # Remove the memory from the input tables
      rm(linear_model_input_table)
      gc()
      
      # Tidy the output
      summary_table <- tidy(lm_fit)
      print(lm_fit)
      # Restrict the output table to just particular columns & just the mutation status of regprot i
      #summary_table <- lm_fit[lm_fit$term == "MutStat_i" | lm_fit$term == "CNAStat_i",]
      
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
    return(data.table::rbindlist(regprot_i_results_df_list))
  })
  
  # Now we have a list of the result omit ts tables, one table per regulatory protein, with entries for each target. 
  # Bind these all together into one master table.
  #master_df <- do.call("rbind", results_df_list)
  master_df <- data.table::rbindlist(results_df_list)
  
  # Return this master DF
  return(master_df)
}


############################################################
#### FILL IN REGULATORY PROTEIN R_I INPUTS TO TABLE
############################################################
#' Function takes a data frame of inputs for the given patient p_j, as well as 
#' the data files for regulatory protein r_i (and its IDs), and uses them to 
#' construct a partial tabular input to a linear model function for a given 
#' regulatory protein r_i of the following form:
#' Columns: Linear model input variables
#' Rows: Patients 
#' 
#' @param patient_df a 'starter DF' that has all the patient characteristics 
#' @param regprot_i_uniprot the uniprot ID of regulatory protein i
#' @param regprot_i_ensg the ensembl ID of regulatory protein i
#' @param mutation_regprot_df the mutation DF for regulatory proteins 
#' @param methylation_df the methylation DF
#' @param cna_df the copy number alteration DF
fill_regprot_inputs <- function(patient_df, regprot_i_uniprot, regprot_i_ensg, 
                                mutation_regprot_df, methylation_df, cna_df) {
  
  # Filter the data frames to look at only this regulatory protein
  mutation_regprot_df <- mutation_regprot_df[Swissprot %like% regprot_i_uniprot]
  cna_df <- filter_cna_by_ensg(cna_df, regprot_i_ensg)   #TODO: try and make this faster
  methylation_df <- filter_meth_by_ensg(methylation_df, regprot_i_ensg)  #TODO: try and make this faster
  
  # Loop through all samples to get the info on this regulatory protein for each
  regprot_rows <- mclapply(1:patient_df[, .N], function(i) {
    sample <- patient_df$sample_id[i]
    
    # Is this regulatory protein mutated in this sample at the given level of specificity?
    # OPT 1: FOR I-DOMAIN & I-BINDING POSITION:
    #if (nrow(mutation_regprot_df %>% filter(Patient == patient)) > 0) {mut_stat <- 1}
    # OPT 2: FOR I-PROTEIN:
    mutation_regprot_df_sub <- mutation_regprot_df[Patient %like% sample]
    if(mutation_regprot_df_sub[, .N] > 0) {mut_stat <- 1}
    else {mut_stat <- 0}
    
    # OPT 1. RAW CNA VALUE: Does this regulatory protein have a CNA?
    cna_val<- unique(as.numeric(cna_df[,colnames(cna_df) == sample, with = FALSE]))
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
    
    # Does this regulatory protein have a methylation marker in cancer?
    meth_stat <- mean.default(as.numeric(methylation_df[,grepl(sample, colnames(methylation_df)), 
                                                        with = FALSE]))
    
    #print(paste("Meth stat", meth_stat))
    #print(paste("CNA stat", cna_stat))
    #print(paste("Mut stat", mut_stat))
    
    if(length(cna_stat) < 1) {cna_stat <- NA}
    else if (is.nan(cna_stat)) {cna_stat <- NA}
    else {cna_stat <- cna_stat}
    
    if(length(meth_stat) < 1) {meth_stat <- NA}
    else if (is.nan(meth_stat)) {meth_stat <- NA}
    else {meth_stat <- meth_stat}
    
    return(data.table("MutStat_i" = mut_stat, "CNAStat_i" = cna_stat, "MethStat_i" = meth_stat))
  })
  
  # Bind these rows into the starter DF
  #regprot_i_df <- do.call("rbind", regprot_rows)
  regprot_i_df <- data.table::rbindlist(regprot_rows)
  colnames(regprot_i_df) <- c("MutStat_i", "CNAStat_i", "MethStat_i")
  starter_df <- cbind(patient_df, regprot_i_df)
  
  # Return this full input DF
  return(starter_df)
}


############################################################
#### FILL IN TARGET T_K INPUTS TO TABLE
############################################################
#' Function takes a data frame of inputs for the given patient p_j and the regulatory 
#' protein r_i, as well as the data files for target t_k (and its ID), and uses 
#' them to construct a full tabular input to a linear model function for a given 
#' regulatory protein i-target gene k pairing of the following form:
#' Columns: Linear model input variables
#' Rows: Patients 
#' @param starter_df a 'starter DF' that has all the patient and reg. protein r_i characteristics 
#' @param targ_k the swissprot ID of target gene k
#' @param targ_k_ensg the ensembl ID of target gene k
#' @param mutation_targ_df the mutation gene target DF
#' @param methylation_df the methylation DF
#' @param cna_df the copy number alteration DF
#' @param expression_df the gene expression DF
fill_targ_inputs <- function(starter_df, targ_k, targ_k_ensg, mutation_targ_df,
                             methylation_df, cna_df, expression_df) {
  #print(targ_k)
  #print(targ_k_ensg)
  
  # Begin by limiting the data frames to just this target gene k
  if (!(targ_k == "" | targ_k_ensg == "")) {
    # Mutation target DF
    mutation_targ_df <- filter_mut_by_uniprot(mutation_targ_df, targ_k) 
    # CNA DF
    cna_df <- filter_cna_by_ensg(cna_df, targ_k_ensg)
    # Methylation DF
    methylation_df <- filter_meth_by_ensg(methylation_df, targ_k_ensg)
    # Expression DF
    exp_rows_to_keep <- unique(unlist(lapply(targ_k_ensg, function(x) grep(x, expression_df$ensg_id))))
    expression_df <- expression_df[exp_rows_to_keep,]
    #print("Expression DF filt by ENSG ID")
    #print(expression_df)
    
    # Continue only if all of these data frames contain information about this target gene k
    if(!(expression_df[, .N] == 0 | cna_df[, .N] == 0 | methylation_df[, .N] == 0 | 
         mutation_targ_df[, .N] == 0)) {
      
      # Loop through all samples to get the info for each target t_k column
      target_rows <- mclapply(1:starter_df[, .N], function(i) {
        sample <- starter_df$sample_id[i]
        
        # Is this target mutated? 
        mut_count <- as.numeric(unlist(mutation_targ_df[, colnames(mutation_targ_df) == sample, 
                                                 with = FALSE]))
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
        cna_val <- unique(as.numeric(unlist(cna_df[,colnames(cna_df) == sample, with = FALSE])))
        # OPT 1. RAW CNA VALUE: Does this regulatory protein have a CNA?
        #if(length(cna_stat) > 1) {cna_stat <- cna_stat[1] + 1}
        #else {cna_stat <- cna_stat + 1}
        # OPT 2: CNA DELETED OR NOT? 
        # Opt. 2a. Include patients with amplifications
        # Opt. 2b. Exclude patients with amplifications
        if (length(cna_val) > 1) {
          if (0 %in% cna_val) {cna_stat <- -1}   # no copies
          #else {cna_stat <- 0}   # amplifications or no change
          else if (1 %in% cna_val) {cna_stat <- 0}  # no change
          else {cna_stat <- 1} 
          #else if (2 %in% cna_val) {cna_stat <- 1}
          #else {cna_stat <- NA} # exclude amplications
        } else {
          if(length(cna_val) == 0) {cna_stat <- NA}
          else if (is.na(cna_val)) {cna_stat <- NA}
          else if (cna_val == 0) {cna_stat <- -1}
          else if (cna_val == 1) {cna_stat <- 0}
          else {cna_stat <- 1}
          #else if (cna_val == 2) {cna_stat <- 1}
          #else {cna_stat <- NA}
        }
        
        # What is the log expression of this target in this sample in cancer?
        # Add a pseudocount so we don't take the log of 0
        exp_stat <- as.numeric(unlist(expression_df[, grepl(sample, colnames(expression_df)), 
                                             with = FALSE])) + 1  # use 1, since we are taking the log

        # Take the mean of all cancer expression values
        #if(is.data.frame(patient_exp_cols)) {
          #if (nrow(patient_exp_cols) > 1) {
            #patient_exp_cols <- colMeans(patient_exp_cols)
            #exp_stat <- log(mean.default(patient_exp_cols))
          #} else {
           # exp_stat <- log(mean.default(as.numeric(patient_exp_cols[,grepl("-0", colnames(patient_exp_cols))])))
          #}
        #} else {exp_stat <- log(mean.default(patient_exp_cols))}
    
        
        # Get the methylation statistic 
        meth_stat <- as.numeric(unlist(methylation_df[, grepl(sample, colnames(methylation_df)), 
                                               with = FALSE]))

        #print(paste("exp_stat", exp_stat))
        #print(paste("mut_stat", mut_stat))
        #print(paste("cna_stat", cna_stat))
        #print(paste("meth_stat", meth_stat))

        # Convert Inf, -Inf, or NaN values to NA
        if (length(exp_stat) < 1) {exp_stat <- NA}
        else if (length(exp_stat > 1)) {exp_stat <- mean.default(exp_stat)}  # if there are two ENSG ids, avg them
        else if(is.infinite(exp_stat) | is.nan(exp_stat)) {exp_stat <- NA}
        else {exp_stat <- exp_stat}
        
        if(length(cna_stat) < 1) {cna_stat <- NA}
        else if (is.nan(cna_stat)) {cna_stat <- NA}
        else {cna_stat <- cna_stat}
        
        if(length(meth_stat) < 1) {meth_stat <- NA}
        else if (length(meth_stat) > 1) {meth_stat <- meth_stat[1]}
        else if (is.nan(meth_stat)) {meth_stat <- NA}
        else {meth_stat <- meth_stat}
        
        # Make row for this patient and return
        return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, 
                          "CNAStat_k" = cna_stat, "MethStat_k" = meth_stat))
      })
      
      # Bind these rows into the starter DF
      targ_k_df <- data.table::rbindlist(target_rows)

      colnames(targ_k_df) <- c("ExpStat_k", "MutStat_k", "CNAStat_k", "MethStat_k")
      full_df <- cbind(starter_df, targ_k_df)
      
      # Return this full input DF
      return(full_df)
    }
  }
  return(NA)
}


############################################################
#### CALL FUNCTION
############################################################
# Run gc to free up any loose memory
rm(all_genes_id_conv)
gc()

master_df <- run_linear_model(protein_ids_df, sample_targets_DF, patient_df, 
    mutation_targ_df, mutation_regprot_df, methylation_df, 
    cna_df, expression_df, "eQTL")
#master_df <- run_linear_model(protein_ids_df, downstream_target_df, patient_df, 
    #mutation_targ_df, mutation_regprot_df, methylation_df_dependent, 
    #cna_df, expression_df, "meQTL")

# For tester
master_df <- run_linear_model(protein_ids_df_tester, sample_targets_DF, patient_df, 
                              mutation_targ_df, mutation_regprot_df, methylation_df, 
                              cna_df, expression_df, "eQTL")
#master_df <- run_linear_model(protein_ids_df_tester, sample_targets_DF, patient_df, 
    #mutation_targ_df, mutation_regprot_df, methylation_df_dependent, 
    #cna_df, expression_df, "meQTL")


############################################################
### INITIAL PROCESSING OF RAW FILE 
############################################################
# Order the file by p-value
master_df <- master_df[order(p.value)]

# Write the results to a file
fwrite(master_df, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/iprotein_output_results_TP53_TMM_bucketCNA_iprot_iciTotFrac_uncorrected.csv", sep = ""))
#fwrite(master_df, paste(main_path, "Linear Model/FOXA1/Non-Tumor-Normal Matched/iprotein_output_results_FOXA1_TMM_rawCNA_iprot_iciTotFrac_uncorrected.csv", sep = ""))

# Limit the data frame to just the term of interest (typically either MutStat_i or CNAStat_i)
master_df_mut <- master_df[master_df$term == "MutStat_i",]
master_df_cna <- master_df[master_df$term == "CNAStat_i",]

# Write these to files
fwrite(master_df_mut, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/iprotein_output_results_TP53_TMM_bucketCNA_iprot_iciTotFrac_MUT_uncorrected.csv", sep = ""))
fwrite(master_df_cna, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/iprotein_output_results_TP53_TMM_bucketCNA_iprot_iciTotFrac_CNA_uncorrected.csv", sep = ""))


