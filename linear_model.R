############################################################
### Linear Model 
### Written By: Sara Geraghty, July 2020
############################################################

#!/usr/bin/env Rscript

#library(biomaRt)
library(parallel)
library(rlang)
library(dplyr)
library(broom)
library(data.table)
library(speedglm)
library(argparse)
library(stringr)
library(glmnet)
library(cli)


# Source other files needed
source_path <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/Sara_LinearModel/"
source(paste(source_path, "linear_model_helper_functions.R", sep = ""))
source(paste(source_path, "general_important_functions.R", sep = ""))


############################################################
# SET UP PARSER ARGUMENTS
############################################################
# Create parser object
parser <- ArgumentParser()

# Specify desired options for the parser

# Data set we are running on
parser$add_argument("--dataset", default = "TCGA", type = "character", 
                    help = "TCGA, METABRIC, ICGC, or CPTAC3. Default is TCGA.")

# Type of QTL to run 
parser$add_argument("--QTLtype", default = "eQTL", type = "character", 
                    help = "eQTL or meQTL. Default is eQTL.")


# Whether we are using tumor-normal matched or non-matched data
parser$add_argument("--tumNormMatched", default = "FALSE", type = "character",
                    help = "TRUE/FALSE for whether or not we are looking at tumor-normal-matched data. Default is FALSE.")

# Cancer type we are looking at (or pan-cancer), and patients we are looking at
parser$add_argument("--cancerType", default = "BRCA", type = "character",
                    help = "Type of cancer we're looking at, or PanCancer. Default is BRCA.")
parser$add_argument("--specificTypes", default = "ALL", type = "character",
                    help = "If & only if the cancer type is pan-cancer, denotes which specific cancer types we're including (running per-cancer). If ALL, run pan-cancer.")
parser$add_argument("--clinical_df", default = "clinical.csv", type = "character",
                    help = "If & only if the cancer type is pan-cancer and specific types is not ALL, a clinical data frame name is given to link cancer types to patient IDs.")
parser$add_argument("--patientsOfInterest", default = "", type = "character",
                    help = "Rather given a particular cancer type or types, further restricts the patient pool to particular patients (e.g. from a particular subtype). A filename with a vector of patient IDs (no header) is optionally provided here.")
parser$add_argument("--patientsOfInterestLabel", default = "", type = "character",
                    help = "A character label for the specific patient population we are looking at here (for file naming purposes).")
parser$add_argument("--removeMetastatic", default = "FALSE", type = "character",
                    help = "A TRUE/ FALSE value indicating whether or not we want to exclude metastatic samples from the analysis.")

# Randomization
parser$add_argument("--randomize", default = "FALSE", type = "character", 
                    help = "TRUE/FALSE for whether or not we perform a randomization test. Default is FALSE")

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
parser$add_argument("--incl_nextMutDriver", default = "FALSE", type = "character",
                    help = "Currently implemented only for TP53 and PIK3CA in BRCA. If TRUE, will include mutation status, CNA status, and mutation interaction term of the other (TP53 if PIK3CA, PIK3CA if TP53) in the model as well.")

# The files to use for the other inputs
parser$add_argument("--protein_ids_df", default = "iprotein_protein_ids_df.csv", 
                    type = "character",
                    help = "The name of the protein IDs data frame input (with ENSG and Uniprot IDs). [default %(default)s]")
parser$add_argument("--expression_df", default = "expression_tmm_CancerOnly_IntersectPatients.csv", 
                    type = "character",
                    help = "The name of the expression data frame input. [default %(default)s]")
parser$add_argument("--methylation_df", default = "methylation_bucketed_Beta_CancerOnly_IntersectPatients.csv", 
                    type = "character",
                    help = "The name of the methylation data frame input. [default %(default)s]")
parser$add_argument("--methylation_df_meQTL", default = "methylation_Beta_CancerOnly_IntersectPatients.csv", 
                    type = "character",
                    help = "The name of the methylation data frame input for meQTL, if using a bucketed version for the regulatory protein. [default %(default)s]")
parser$add_argument("--mutation_targ_df", default = "mut_count_matrix_missense_CancerOnly_IntersectPatients.csv", 
                    type = "character",
                    help = "The name of the mutation target data frame input. [default %(default)s]")
parser$add_argument("--mutation_regprot_df", default = "iprotein_results_missense_CancerOnly_IntersectPatients.csv", 
                    type = "character",
                    help = "The name of the mutation regulatory protein data frame input. [default %(default)s]")
parser$add_argument("--cna_df", default = "CNA_AllGenes_CancerOnly_IntersectPatients.csv", 
                    type = "character",
                    help = "The name of the CNA data frame input. [default %(default)s]")
parser$add_argument("--patient_df", default = "combined_patient_sample_cibersort_total_frac_tmm_IntersectPatients.csv", 
                    type = "character",
                    help = "The name of the patient sample data frame input. [default %(default)s]")
parser$add_argument("--target_df", default = "allgene_targets.csv", type = "character",
                    help = "The name of the ChIP-eat or curated targets data frame. [default %(default)s]")
parser$add_argument("--neighboring_cna_df", default = "neighboring_CNA_CancerOnly_IntersectPatients.csv",
                    type = "character",
                    help = "The name of the neighboring CNA data frame input. [default %(default)s]")
parser$add_argument("--functional_copies_df", default = "tp53_functional_copies_per_sample.csv",
                    type = "character",
                    help = "The name of a functional copies DF for a particular tester gene of interest, default TP53.")

# The decision for if/ how to bucket CNAs and methylation
parser$add_argument("--cna_bucketing", default = "bucket_inclAmp", 
                    type = "character",
                    help = "If/ the type of bucketing we use for CNA values. Default is bucket_inclAmp, but rawCNA, bucket_exclAmp, bucket_justAmp, and bucket_justDel are also options.")
parser$add_argument("--meth_bucketing", default = TRUE, type = "character",
                    help = "If/ the type of bucketing we use for methylation values. Default is FALSE.")

# Label the type of methylation data we are using (Beta, M, or Threshold)
parser$add_argument("--meth_type", default = "Beta", type = "character",
                    help = "The type methylation values we are using. Default is Beta, but other options are M and Threshold_X, where X is the threshold value.")

# Label for the explanatory variable we are interested in (which output files to write/ output visualizations to make)

# What type of test we are running, so we can write the results files to an appropriate directory. 
# Only necessary if --test is F.
parser$add_argument("--run_name", default = "cancerRelated", type = "character",
                    help = "Provide a name for the run, if not testing on just one regulatory protein, to write output files to appropriate directory. [default %(default)s]") 

# A name for the group of targets being tested, for labeling the output file.
parser$add_argument("--targets_name", default = "allGenes", type = "character",
                    help = "Provide a name for the group of targets being tested, in order to properly name the output file. [default %(default)s]")

# Whether or not we are including PEER factors and PCs as covariates in the model
parser$add_argument("--num_PEER", default = 0, type = "integer",
                    help = "A value from 0 to 10 indicating the number of PEER factors we are using as covariates in our model. Default is 0.")
parser$add_argument("--num_pcs", default = 2, type = "integer",
                    help = "A value between 0 and 2 indicating the number of principal components we are using as covariates in our model. Default is 2.")
parser$add_argument("--mut_pc_vect", default = "0", type = "character",
                    help = "either '0' for no mutational vectors, or a vector of PCs to include, e.g. 1,3,4. Default is 0.")

# Add a flag for running collinearity diagnostics/ including residuals; can be useful but adds to runtime.
parser$add_argument("--collinearity_diagn", default = "FALSE", type = "character",
                    help = "A TRUE/ FALSE value indicating whether or not we want detailed collinearity diagnostics. Can be useful but adds to runtime; not suggested for large runs. Default is FALSE.")
parser$add_argument("--inclResiduals", default = "FALSE", type = "character",
                    help = "A TRUE/FALSE value indicating whether we want to include residuals in the output DF. Default is FALSE.")


# Add a flag for adding a regularization method
parser$add_argument("--regularization", default = "None", type = "character",
                    help = "The name of the regularization method being used. Currently only implemented for L1 and L2. Default is None.")

# Add a flag for keeping only trans pairings
parser$add_argument("--removeCis", default = "TRUE", type = "character",
                    help = "A TRUE/ FALSE value to indicate whether or not we want to remove cis pairings and look only at trans pairings. Defaults to true, but can be set to false.")

# A flag for using the number of functional copies of regulatory protein per sample, rather than the mutation and CNA status 
parser$add_argument("--useNumFunctCopies", default = "FALSE", type = "character",
                    help = "A TRUE/ FALSE value to indicate whether or not we want to use the allele-specific number of functional copies for the given regulatory tester protein, rather than its separate mutation and CNA status.")
#parser$add_argument("--explan_var_interest", default = "MutStat", type = "character",
                    #help = "The name of the explantory variable we are interested in. If CNAStat, program will automatically add a covariate to account for the precense of a co-amplified or deleted neighboring driver gene.")

# Add a flag for debugging
parser$add_argument("--debug", default = "FALSE", type = "character",
                    help = "A TRUE/ FALSE value indicating whether or not we want detailed printing for debugging purposes. Default is FALSE.")

# Add a flag for including only select covariates in the model, given as a character list with each item separated by a semicolon
parser$add_argument("--select_args", default = "ALL", type = "character",
                    help = "Option to provide a semicolon-separated character string listing all covariates to include, e.g. MutStat_i;CNAStat_i, etc. If this option is chosen, dependent variable must also be listed. [default %(default)s]")
parser$add_argument("--select_args_label", default = "", type = "character",
                    help = "If we are providing a character string for --select_args, this is an option to provide a label for what arguments we are selecting for labeling the output files and visualizations.")


# Parse the given arguments
tryCatch({
  args <- parser$parse_args()
}, error = function(cond) {
  print(cond)
  print(traceback())
})



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
debug <- str2bool(args$debug)
collinearity_diagn <- str2bool(args$collinearity_diagn)
removeCis <- str2bool(args$removeCis)
removeMetastatic <- str2bool(args$removeMetastatic)
useNumFunctCopies <- str2bool(args$useNumFunctCopies)
incl_nextMutDriver <- str2bool(args$incl_nextMutDriver)
inclResiduals <- str2bool(args$inclResiduals)


############################################################
# SET MAIN PATH AND IMPORT GENERALLY NEEDED FILES
############################################################
# Set the paths
input_file_path <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/"


if(args$cancerType == "BRCA") {
  prot_path <- paste(input_file_path, "BRCA/", sep = "")
  targ_path <- paste(input_file_path, "BRCA/target_lists/", sep = "")
  
  if(args$dataset == "METABRIC") {
    main_path <- paste0(input_file_path, "METABRIC/")
  } else {
    if(!tumNormMatched) {
      main_path <- paste(input_file_path, "BRCA/tumor_only/", sep = "")
    } else {
      main_path <- paste(input_file_path, "BRCA/tumor_normal_matched/", sep = "") 
    }
  }
} else {
  prot_path <- paste(input_file_path, "PanCancer/", sep = "")
  targ_path <- paste(input_file_path, "PanCancer/target_lists/", sep = "")
  if(!tumNormMatched) {
    main_path <- paste(input_file_path, "PanCancer/tumor_only/", sep = "")
  } else {
    main_path <- paste(input_file_path, "PanCancer/tumor_normal_matched/", sep = "")
  }
}

# If we are looking at tester targets, we need to go one directory deeper
if(grepl("curated", args$targets_name)) {
  targ_path <- paste(targ_path, "curated_targets/", sep = "")
} else if (grepl("chipeat", args$targets_name)) {
  targ_path <- paste(targ_path, "chipeat_targets/", sep = "")
} else {targ_path <- targ_path}

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- fread("/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/all_genes_id_conv.csv", 
                           header = TRUE)


############################################################
# IMPORT REGULATORY PROTEINS
############################################################
# Regulatory proteins will vary depending on the level of specificity in question;
# import only the proteins at the given specificity level of interest under label "protein_ids_df"
protein_ids_df <- data.frame()
if(!test) {
  protein_ids_df <- fread(paste(prot_path, args$protein_ids_df, sep = ""), 
                          header = TRUE)
} else {
  sample_protein_uniprot <- args$tester_uniprot_id
  sample_protein_ensg <- args$tester_ensg_id
  sample_protein_name <- args$tester_name
  
  print(paste("Tester Protein:", sample_protein_name))
  
  protein_ids_df <- data.table(swissprot_ids = sample_protein_uniprot,
                               ensg_ids = sample_protein_ensg)
}


############################################################
# IMPORT EXPRESSION FILES
############################################################
expression_df <- fread(paste(main_path, paste("expression/", args$expression_df, sep = ""), sep = ""),
                       header = TRUE)

if((colnames(expression_df)[1] == "") & (colnames(expression_df)[2] == "")) {
  expression_df <- expression_df[,2:ncol(expression_df)]
}
if(!("ensg_id" %fin% colnames(expression_df))) {
  colnames(expression_df)[1] <- "ensg_id"
}

if(debug) {
  print("Expression DF")
  print(head(expression_df))
}

# Determine what type of expression this is from the filename
is_rank_or_quant_norm <- FALSE
if(grepl("rank", args$expression_df) | grepl("quantile", args$expression_df)) {is_rank_or_quant_norm <- TRUE}

# Determine whether or not we will need to log-transform the data, or whether this
# has already been done in some way (e.g. by an inverse normal transformation or z-score normalization)
log_expression <- TRUE
if(grepl("transform", args$expression_df) | grepl("zscore", args$expression_df)) {log_expression <- FALSE}


############################################################
# IMPORT METHYLATION FILES
############################################################
methylation_df <- fread(paste(main_path, paste("methylation/", args$methylation_df, sep = ""), sep = ""),
                        header = TRUE)
methylation_df <- methylation_df[,2:ncol(methylation_df)]
colnames(methylation_df)[which(colnames(methylation_df) == "Gene_Symbol")] <- "gene_name"

# If we're using bucketed methylation values AND running an meQTL, we'll need an 
# additional DF with the raw Beta values to use as the outcome variable
if((args$QTLtype == "meQTL") & (grepl("bucketed", args$methylation_df))) {
  methylation_df_meQTL <- fread(paste(main_path, paste("methylation/", args$methylation_df_meQTL, sep = ""), sep = ""),
                                header = TRUE)
  methylation_df_meQTL <- methylation_df_meQTL[,2:ncol(methylation_df_meQTL)]
  colnames(methylation_df_meQTL)[which(colnames(methylation_df_meQTL) == "Gene_Symbol")] <- "gene_name"
  
} else {methylation_df_meQTL <- NA}

if(debug) {
  print("Methylation DF")
  print(head(methylation_df))
  
  if(!is.na(methylation_df_meQTL)) {
    print("Methylation DF meQTL")
    print(head(methylation_df_meQTL))
  }
}


############################################################
# IMPORT GENE TARGET MUTATION FILES
############################################################
mutation_targ_df <- fread(paste(main_path, paste("genetarg_mutation/", args$mutation_targ_df, sep = ""), sep = ""),
                          header = TRUE)
mutation_targ_df <- mutation_targ_df[,2:ncol(mutation_targ_df)]
colnames(mutation_targ_df)[which(colnames(mutation_targ_df) == "Gene_Symbol")] <- "gene_name"
mutation_targ_df <- mutation_targ_df[, Swissprot := unlist(lapply(mutation_targ_df$gene_name, function(x) 
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'uniprot_gn_id'])), collapse = ";")))]

if(debug) {
  print("Mutation Targ DF")
  print(head(mutation_targ_df))
}


############################################################
# IMPORT REGULATORY PROTEIN MUTATION FILES
############################################################
mutation_regprot_df <- fread(paste(main_path, paste("regprot_mutation/", args$mutation_regprot_df, sep = ""), sep = ""),
                             header = TRUE)
mutation_regprot_df <- mutation_regprot_df[,5:ncol(mutation_regprot_df)]   # remove first few meaningless columns, if necessary
mutation_regprot_df <- mutation_regprot_df[!duplicated(mutation_regprot_df),] # remove any duplicate rows

if(debug) {
  print("Mutation Regprot DF")
  print(head(mutation_regprot_df))
}



############################################################
# IMPORT CNA FILES
############################################################
cna_df <- fread(paste(main_path, paste("cna/", args$cna_df, sep = ""), sep = ""),
                header = TRUE)
colnames(cna_df)[1] <- "ensg_id"

if(debug) {
  print("CNA DF")
  print(head(cna_df))
}

#if(args$explan_var_interest == "CNAStat") {
  #neighboring_cna_df <- fread(paste(main_path, paste("cna/", args$neighboring_cna_df, sep = ""), sep = ""),
                              #header = TRUE)
  #neighboring_cna_df <- neighboring_cna_df[neighboring_cna_df$ensg_id %in% protein_ids_df$ensg_id,]  #TODO: fix this to match
#}


############################################################
# IMPORT PATIENT/SAMPLE FILES
############################################################
patient_df <- fread(paste(main_path, paste("patient/", args$patient_df, sep = ""), sep = ""),
                    header = TRUE)
colnames(patient_df)[1] <- "sample_id"

# Restrict the Mut_PC columns to only those given in the argument vector
if(args$dataset == "TCGA") {
  patient_df <- restrict_mut_pc_cols(patient_df, args$mut_pc_vect)
}

if(debug) {
  print("Patient DF")
  print(head(patient_df))
}

# Read if the patient IDs of interest, if provided
patients_of_interest <- ""
if(args$patientsOfInterest != "") {
  patients_of_interest <- as.character(unlist(fread(paste(main_path, paste("patient/groups_of_interest/", 
                                                 args$patientsOfInterest, sep = ""), sep = ""), 
                                header = FALSE)[,1]))
}
if(debug) {
  print("Patients of Interest")
  print(head(patients_of_interest))
}


############################################################
# IMPORT GENE TARGETS DATA FRAME
############################################################
targets_DF <- fread(paste(targ_path, args$target_df, sep = ""), header = TRUE)

# Remove rows with any blank elements
targets_DF <- targets_DF[!(targets_DF$swissprot == ""),]

# Remove the first column, if needed (for all but metabolic targets)
if(ncol(targets_DF) > 2) {
  targets_DF <- targets_DF[,2:3]
}

# Regardless of method, limit targets to only those that overlap the genes in the expression DF
rows_to_keep <- unlist(lapply(1:length(targets_DF$ensg), function(i) {
  targs <- unlist(strsplit(as.character(targets_DF[i,'ensg']), ";", fixed = TRUE))
  if (any(targs %fin% expression_df$ensg_id)) {return(TRUE)}
  else {return(FALSE)}
}))
targets_DF <- targets_DF[rows_to_keep,]

if(debug) {
  print("Targets DF")
  print(head(targets_DF))
}


############################################################
# IMPORT NUM. FUNCTIONAL COPIES DATA FRAME (OPT)
############################################################
# If we are using a data frame with the number of functional copies for the
# given regulatory protein, import that here
if(useNumFunctCopies) {
  num_funct_copies_DF <- fread(paste(main_path, args$functional_copies_df, sep = ""), header = TRUE)
  colnames(num_funct_copies_DF)[1] <- "sample_id"
  
  if(debug) {
    print("Num. Functional Copies DF")
    print(head(num_funct_copies_DF))
  }
}


############################################################
#### OVERVIEW OF LINEAR MODEL
############################################################
### PATIENT CHARACTERISTICS ###
# Gender : Gender of all samples 1..j..M (0 if male, 1 if female) -- ONLY INCLUDE IF NOT BRCA
# Age : Age of all samples 1..j..M in years / 100
# If Buckets: 1 if 0-9, 2 if 10-19, 3 if 20-29, 4 if 30-39, 5 if 40-49, 6 if 50-59, 7 if 60-69, 8 if 70-79, 9 if 80-89, 10 if 90+
# Race : Race/ Ethnicity of all samples 1..j..M -- CURRENTLY EXCLUDING
# Buckets: 1 if White (not Latinx), 2 if White (Latinx), 3 if black or African, 4 if Asian, 5 if other
# Prior_malig : If all samples 1..j..M had a prior malignancy (0 if no, 1 if yes)
# Treatment_rad : If all samples 1..j..M were treated with radiation (0 is not treated, 1 is treated) 
# Treatment_pharm : If all samples 1..j..M were treated with pharmaceutical therapy (0 is not treated, 1 is treated) 
# TotalNumMut : Total number of missense mutations each patient j has (across all samples 1..j..M)
# Tumor_purity : An estimate of the 'purity' of the tumor (fraction of tumor cells within tumor sample) 
# Buckets: 1 if low purity (0-0.5), 2 if medium purity (0.5-0.75), 3 if high purity (0.75-1.0)
# Total_IC_Frac : An estimate of the fraction of tumor cells in the sample, from CIBERSORT Abs
# Buckets: 1 if low ICI (0-0.3), 2 if medium ICI (0.3-0.7), 3 if high ICI (0.7-1.0), for TCGA
# Tumor Type/ Subtype : The type or subtype classification of the tumor sample (buckets equivalent to the number of types/ subtypes)

### REGULATORY PROTEIN i CHARACTERISTICS ###
# MutStat_i : The mutation status of protein i (0 if not mutated, 1 if mutated) across all tumor samples 1..j..M
# MethStat_i : The methylation status of protein i (0 if not methylated, 1 if methylated) across all tumor samples 1..j..M  
# CNA_i : The copy number alteration status of protein i (1 if amplified, -1 if deleted, 0 if no amp. or del. OR copy #) 
# across all tumor samples 1..j..M

# ALTERNATIVE: Rather than MutStat_i and CNA_i, can also use a bulk "Num. Functional Copies" (Num_F_Copies_i) term,
# which is the log2 of the number of functional copies a given sample possesses of the regulatory protein i 
  # Can also bucket (0 = no functional copies, 1 = 1 functional copy, 2 = at least 2 functional copies)

### TARGET GENE k CHARACTERISTICS ###
# MutStat_k : The mutation status of gene k (0 if not mutated, 1 if mutated) across all samples 1..j..M
# CNA_k : the copy number for target gene k across all samples 1..j..M

### DEPENDENT VARIABLE ###
# log2(ExpStat_k,c) or log2(Meth_k,c) : the cancer expression or methylation value (tumor-only)
# log2(ExpStat_k,c/Exp_k,n) or log2(Meth_k,c/ Meth_k,n) : the cancer differential expression or methylation value (tumor-normal matched)

### OFFSETS/ SCALING FACTORS ###
# log2(N_j) : log of effective library size (sample-specific scaling factor); the total number of reads for 
# tumor sample j after cross-sample normalization, across all samples 1..j..M (tumor-only)
# log2(N_j,c / N_j,n) log of differential effective library size (sample-specific scaling factor); the total number of reads for 
# tumor sample j divided by normal sample j after cross-sample normalization, across all matched samples 1..j..M (tumor-normal-matched)
# NOTE: This is important only for TMM and FPKM expression, not for rank- or quantile-normalized expression


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
#' @param methylation_df_meQTL OPTIONAL: table of methylation results (rows are 
#' proteins, columnsare patients, entries are raw methylation values to be used
#' for regulatory protein i in the case of an meQTL and bucketed methylation for t_k)
#' @param cna_df table of CNA per gene (rows are proteins, columns are patients, 
#' entries are CNA values)
#' @param expression_df table of expression (rows are genes, columns are patients, 
#' entries are expression values)
#' @param num_funct_copies_df a data frame with the number of functional copies
#' of a given regulatory protein for each sample of interest 
#' @param neighboring_cna_df a data frame with information about neighboring
#' driver genes on the same chromosome that have matching CNA status to a given gene
#' @param is_rank_or_quant_norm a TRUE/FALSE value indicating whether this is rank-
#' or quantile-normalized data (i.e. whether we need to include library size as an
#' offset)
#' @param log_expression a TRUE/FALSE value indicating whether we need to take the
#' log2(exp + 1), or whether expression has already been transformed to a normal
#' distribution in some way
#' @param analysis_type a string label that reads either "eQTL" or "meQTL" to 
#' determine what kind of model we should run
#' @param tumNormMatched a TRUE/FALSE value indicating whether or not the analysis is 
#' tumor-normal matched
#' @param randomize a TRUE/FALSE value indicating whether or not we are randomizing 
#' expression
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. 
#' Possible values are "bucket_inclAmp", "bucket_exclAmp", "bucket_justAmp", 
#' "bucket_justDel" and "rawCNA"
#' @param meth_bucketing a TRUE/FALSE value indicating whether or not we are bucketing
#' methylation values
#' @param useNumFunctCopies a TRUE/FALSE value indicating whether or not we are using
#' the number of functional copies rather than the mutation/CNA status information
#' @param num_PEER a value from 0-10 indicating the number of PEER factors we 
#' are including as covariates in the model
#' @param num_pcs a value from 0-2 indicating the number of  
#' principal components we are including as covariates in the model 
#' @param debug a TRUE/FALSE value indicating if we are in debug mode and should 
#' use additional prints
#' @param collinearity_diagn a TRUE/FALSE value indicating if we are running 
#' collinearity diagnostics
#' @param outpath a string with the local path for debugging/ collinearity files 
#' we'll be writing from inside the function
#' @param outfn a string with a generic output filename based on the characteristics
#' of the given run; can be modified to write files within this function
#' @param regularization a string with the type of regularization to be used; currently
#' only implemented for "L2" or "None"
#' @param covs_to_incl a vector of covariates that we want to include, or ALL if we
#' want to include all implemented covariates
#' @param removeCis a TRUE/FALSE value indicating whether or not we are eliminating 
#' all cis pairings
#' @param all_genes_id_conv an ID conversion file from BioMart, for use in removing 
#' cis pairings if needed
#' @param incl_nextMutDriver a TRUE/FALSE value indicating, if we are running just on
#' TP53 or PIK3CA in BRCA, whether or not we include mutation/ CNA/ mutation interaction
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3', to denote what columns
#' we are referencing/ specific data types we are using
#' @param inclResiduals a TRUE/FALSE value indicating whether or not we want to include 
#' residuals in our output DF, to evaluate model fit
run_linear_model <- function(protein_ids_df, downstream_target_df, patient_df, 
                             mutation_df_targ, mutation_df_regprot, methylation_df, 
                             methylation_df_meQTL, cna_df, expression_df, num_funct_copies_df,
                             neighboring_cna_df, is_rank_or_quant_norm, log_expression, analysis_type, 
                             tumNormMatched, randomize, cna_bucketing, meth_bucketing, 
                             useNumFunctCopies, num_PEER, num_pcs, 
                             debug, collinearity_diagn, outpath, outfn, regularization, 
                             covs_to_incl, removeCis, all_genes_id_conv, incl_nextMutDriver,
                             dataset, inclResiduals) {
  
  # We need to get a mini-table for every r_i, t_k combo that we rbind into a master table of results
  results_df_list <- lapply(1:protein_ids_df[, .N], function(i) {
    
    # Get the given regulatory protein r_i's Swissprot & ENSG IDs
    regprot <- protein_ids_df$swissprot_ids[i]
    regprot_ensg <- unlist(strsplit(protein_ids_df$ensg_ids[i], ";", fixed = TRUE))
    
    print(paste("Regulatory protein", paste(i, paste("/", protein_ids_df[, .N]))))
    
    # Create a starter table that gets the mutation, CNA, and methylation status 
    # for this regulatory protein and binds it to the patient DF
    starter_df <- fill_regprot_inputs(patient_df = patient_df, 
                                      regprot_i_uniprot = regprot, 
                                      regprot_i_ensg = regprot_ensg, 
                                      mutation_regprot_df = mutation_df_regprot, 
                                      methylation_df = methylation_df, 
                                      cna_df = cna_df,
                                      num_funct_copies_df = num_funct_copies_DF,
                                      useNumFunctCopies = useNumFunctCopies,
                                      cna_bucketing = cna_bucketing, 
                                      meth_bucketing = meth_bucketing, 
                                      tumNormMatched = tumNormMatched, 
                                      incl_nextMutDriver = incl_nextMutDriver,
                                      debug = debug, dataset = dataset)
    
    if(debug) {
      print("Starter DF, with Regprot Inputs")
      print(head(starter_df))
    }
    
    # If remove cis is TRUE, limit the downstream target DF to only trans pairings
    if(removeCis) {
      downstream_target_df <- subset_to_trans_targets(regprot_ensg, downstream_target_df,
                                                      all_genes_id_conv, debug)
    }
    
    # Loop through all the targets for this protein of interest and create a table 
    # for each of them from this starter DF
    regprot_i_results_df_list <- lapply(1:downstream_target_df[, .N], function(k) {
      
      # Get the target t_k's Swissprot & ENSG IDs
      targ <- unlist(strsplit(downstream_target_df$swissprot[k], ";", fixed = TRUE))
      targ_ensg <- unlist(strsplit(downstream_target_df$ensg[k], ";", fixed = TRUE))
      
      print(paste("Target gene", paste(k, paste("/", downstream_target_df[, .N]))))
      
      # OPT: Filter outlier expression for each group (mutated, unmutated) using 
      # a standard (1.5 x IQR) + Q3 schema
      #starter_df <- filter_expression_df(expression_df, starter_df, targ_ensg)
      
      if(length(intersect(targ_ensg, regprot_ensg)) == 0) {
        # Create a full input table to the linear model for this target
        #print(methylation_df_meQTL)
        methylation_df_targ <- methylation_df
        if(!is.null(methylation_df_meQTL)) {
          if(!is.na(methylation_df_meQTL)) {
            methylation_df_targ <- methylation_df_meQTL
          }
        } 
        
        lm_input_table <- fill_targ_inputs(starter_df = starter_df, targ_k = targ, 
                                           targ_k_ensg = targ_ensg, 
                                           mutation_targ_df = mutation_df_targ,
                                           methylation_df = methylation_df_targ, 
                                           cna_df = cna_df, expression_df = expression_df, 
                                           log_expression = log_expression,
                                           cna_bucketing = cna_bucketing,
                                           meth_bucketing = meth_bucketing, 
                                           analysis_type = analysis_type,
                                           tumNormMatched = tumNormMatched, 
                                           debug = debug, dataset = dataset)
        
        
        if(length(lm_input_table) == 0) {return(NA)}
        if (is.na(lm_input_table)) {return(NA)}
        if(lm_input_table[, .N] == 0) {return(NA)}
        
        # If a row has NA, remove that whole row and predict without it
        #lm_input_table <- na.exclude(lm_input_table)
        lm_input_table <- na.omit(lm_input_table)
        
        # Make sure everything is numeric type
        #lm_input_table[,2:ncol(lm_input_table)] <- sapply(lm_input_table[,2:ncol(lm_input_table)], as.numeric) 
        
        if(debug) {
          print("LM Input Table")
          print(head(lm_input_table))
          #new_fn <- str_replace(outfn, "output_results", paste(targ[1], paste(regprot, "lm_input_table", 
                                                                             # sep = "_"), sep = "_"))
          #fwrite(lm_input_table, paste(outpath, paste(new_fn, ".csv", sep = ""), sep = "/"))
        }
        
        # Randomize expression across cancer and normal, across all samples
        if (randomize) {
          if(analysis_type == "eQTL") {lm_input_table$ExpStat_k <- sample(lm_input_table$ExpStat_k)}
          else if (analysis_type == "meQTL") {lm_input_table$MethStat_k <- sample(lm_input_table$MethStat_k)}
          else {print(paste("Analysis type is invalid:", analysis_type))}
        }
        
        formula <- construct_formula(lm_input_table, analysis_type, num_PEER, num_pcs,
                                     cna_bucketing, meth_bucketing, covs_to_incl)
        
        if (regularization == "None") {
          
          if(debug) {
            print(paste("Formula:", formula))
          }
          
          if(inclResiduals) {
            lm_fit <- tryCatch(
              {
                speedglm::speedlm(formula = formula, data = lm_input_table, fitted = TRUE)  # Do not include library size as offset
                #speedglm::speedlm(formula = formula, offset = log2(Lib_Size), 
                #data = lm_input_table)
              }, error = function(cond) {
                message("There was a problem with this linear model run:")
                message(cond)
                print(formula)
                print(head(lm_input_table))
                return(NA)
              })
          } else {
            lm_fit <- tryCatch(
              {
                speedglm::speedlm(formula = formula, data = lm_input_table)  # Do not include library size as offset
                #speedglm::speedlm(formula = formula, offset = log2(Lib_Size), 
                #data = lm_input_table)
              }, error = function(cond) {
                message("There was a problem with this linear model run:")
                message(cond)
                print(formula)
                print(head(lm_input_table))
                return(NA)
              })
          }
          
        } else if (regularization == "L2" | regularization == "ridge" | 
                   regularization == "L1" | regularization == "lasso") {
          
          lm_fit <- run_regularization_model(formula, lm_input_table, type = regularization, debug)
          
        } else {
          lm_fit <- NA
          print("Model only has L1, L2 or None implemented for regularization. Please provide one of these options.")
        }
        
        # Tidy the output
        summary_table <- tidy(lm_fit)
        if(inclResiduals) {
          summary_table$Residual <- lm_input_table$ExpStat_k - predict(lm_fit)
        }
        
        # Restrict the output table to just particular columns, if desired
        #summary_table <- lm_fit[lm_fit$term == "MutStat_i" | lm_fit$term == "CNAStat_i",]
        
        # Add a column for the regulatory protein and target gene to ID them
        summary_table$R_i <- regprot
        summary_table$T_k <- paste(targ, collapse = ";")
        
        if(debug) {
          print("Summary Table")
          print(head(summary_table))
          #new_fn <- str_replace(outfn, "output_results", paste(targ[1], paste(regprot, "summary_table", 
                                                                             # sep = "_"), sep = "_"))
          #new_fn <- paste(new_fn, ".csv", sep = "")
          #fwrite(summary_table, paste(outpath, new_fn, sep = "/"))
        }
        
        # Run collinearity diagnostics 
        if(collinearity_diagn) {
          print("Running Collinearity Diagnostics")
          #fit_lm <- lm(formula = formula, offset = log2(Lib_Size), 
                       #data = lm_input_table)
          fit_lm <- lm(formula = formula, offset = log2(Lib_Size), 
                       data = lm_input_table)
          source(paste(getwd(), "run_collinearity_diagnostics.R", sep = "/"), local = TRUE)
        }
        
        # Remove the memory from the input tables
        rm(lm_input_table)
        gc()
        
        # Return this summary table
        return(summary_table)
      }
      else {return(NA)}
    })
    
    # Now we have a list of output summary DFs for each of this regulatory protein's targets. 
    # Bind them all together.
    #return(do.call("rbind", regprot_i_results_df_list))
    regprot_i_results_df_list <- regprot_i_results_df_list[!is.na(regprot_i_results_df_list)]
    regprot_i_results_df <- as.data.table(rbindlist(regprot_i_results_df_list, fill = TRUE))
    
    if(debug) {
      print("Combined Results DF")
      print(head(regprot_i_results_df))
    }
    
    return(regprot_i_results_df)
  })
        
  # Now we have a list of the results tables, one table per regulatory protein, 
  # with entries for each target. Bind these all together into one master table.
  if(debug) {
    print("Results DF list")
    print(head(results_df_list))
  }
        
  # Check that none of the DF entries are NA or empty data.tables
  results_df_list <- results_df_list[!is.na(results_df_list)]
  results_df_list <- Filter(function(dt) nrow(dt) != 0, results_df_list)
        
  master_df <- rbindlist(results_df_list, use.names = TRUE, fill = TRUE)
  #master_df <- do.call(rbind, results_df_list)
        
  if(debug) {
    print("Master DF")
    print(head(master_df))
  }
        
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
#' @param patient_df a 'starter DF' that has all the patient characteristics 
#' @param regprot_i_uniprot the uniprot ID of regulatory protein i
#' @param regprot_i_ensg the ensembl ID of regulatory protein i
#' @param mutation_regprot_df the mutation DF for regulatory proteins 
#' @param methylation_df the methylation DF
#' @param cna_df the copy number alteration DF
#' @param num_funct_copies_df a data frame with the number of functional copies 
#' of the given regulatory protein (if looking at only 1 regprot)
#' @param useNumFunctCopies a TRUE/FALSE value indicating whether or not we are 
#' using the number of functional copies data
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. Possible
#' values are "bucket_inclAmp", "bucket_exclAmp", and "rawCNA"
#' @param meth_bucketing a TRUE/FALSE value indicating whether or not we are bucketing
#' methylation values
#' @param tumNormMatched a TRUE/FALSE value indicating whether or not the analysis is 
#' tumor-normal matched
#' @param incl_nextMutDriver a TRUE/FALSE value indicating, if we are running just on
#' TP53 or PIK3CA in BRCA, whether or not we include mutation/ CNA/ mutation interaction
#' terms for the other driver gene (e.g. for TP53 if PIK3CA or PIK3CA if TP53). If we
#' are running on TTN (as a control), then we use PIK3CA covariates as a control
#' @param debug a TRUE/FALSE value indicating if we are in debug mode and should 
#' use additional prints
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3', to denote what columns
#' we are referencing/ specific data types we are using
fill_regprot_inputs <- function(patient_df, regprot_i_uniprot, regprot_i_ensg, 
                                mutation_regprot_df, methylation_df, cna_df, num_funct_copies_df,
                                useNumFunctCopies, cna_bucketing, meth_bucketing, tumNormMatched, 
                                incl_nextMutDriver, debug, dataset) {
  
  # Filter the data frames to look at only this regulatory protein
  mutation_regprot_df_sub <- mutation_regprot_df[Swissprot %like% regprot_i_uniprot]
  cna_df_sub <- filter_cna_by_ensg(cna_df, regprot_i_ensg)   #TODO: try and make this faster
  methylation_df_sub <- filter_meth_by_ensg(methylation_df, regprot_i_ensg)  #TODO: try and make this faster
  
  #print(paste("REGPROTi ENSG:", regprot_i_ensg))
  
  # Optionally, add mutation, amplification status covariates for PIK3CA if regprot is TP53
  # and vice versa 
  # Also, optionally add all of these covariates for genes that are not PIK3CA or TP53
  cna_df_pik3ca <- NA
  cna_df_tp53 <- NA
  if(incl_nextMutDriver) {
    if(regprot_i_ensg == "ENSG00000141510") {
      cna_df_pik3ca <- filter_cna_by_ensg(cna_df, "ENSG00000121879")
      mutation_regprot_df_pik3ca <- mutation_regprot_df[Swissprot %like% "P42336"]
    }
    else if(regprot_i_ensg %in% c("ENSG00000121879", "ENSG00000155657")) {
      cna_df_tp53 <- filter_cna_by_ensg(cna_df, "ENSG00000141510")
      mutation_regprot_df_tp53 <- mutation_regprot_df[Swissprot %like% "P04637"]
    }
    else {
      cna_df_pik3ca <- filter_cna_by_ensg(cna_df, "ENSG00000121879")
      mutation_regprot_df_pik3ca <- mutation_regprot_df[Swissprot %like% "P42336"]
      cna_df_tp53 <- filter_cna_by_ensg(cna_df, "ENSG00000141510")
      mutation_regprot_df_tp53 <- mutation_regprot_df[Swissprot %like% "P04637"]
    }
  }

  
  # If tumor-normal matched, limit to just cancer samples (we'll ID the matched normal in each run)
  if(tumNormMatched) {
    patient_df2 <- patient_df[grepl("-0", patient_df$sample_id)]
  } else {
    patient_df2 <- patient_df
  }
  
  # Loop through all tumor samples to get the info on this regulatory protein for each
  regprot_rows <- mclapply(1:patient_df2[, .N], function(i) {
    sample <- patient_df2$sample_id[i]
    if (debug) {
      print(paste("Sample:", sample))
    }
    
    # Is this regulatory protein mutated in this sample at the given level of specificity?
    # OPT 1: FOR I-DOMAIN & I-BINDING POSITION:
    #if (nrow(mutation_regprot_df %>% filter(Patient == sample)) > 0) {mut_stat <- 1}
    # OPT 2: FOR I-PROTEIN:
    mutation_regprot_df_sub <- mutation_regprot_df_sub[Patient %like% sample]
    if(mutation_regprot_df_sub[, .N] > 0) {mut_stat <- 1}
    else {mut_stat <- 0}
    
    # Does this regulatory protein have a CNA in cancer?
    cna_stat <- get_cna_stat(cna_df_sub, sample, cna_bucketing, dataset)
    if(grepl("incl", cna_bucketing)) {cna_stat <- as.integer(cna_stat[[1]])}    
    
    # Does this regulatory protein have a methylation marker in cancer, or differential 
    # methylation between tumor and normal?
    meth_stat <- get_meth_stat(methylation_df_sub, sample, meth_bucketing, tumNormMatched)
    if(meth_bucketing) {
      meth_stat <- as.integer(meth_stat[[1]])
      if(length(meth_stat) > 3) {meth_stat <- meth_stat[1:3]}
    }
    
    # If we are using the functional number of copies info, what is the functional
    # number of copies for this protein in this sample?
    if(useNumFunctCopies) {
      fnc_stat <- num_funct_copies_df[num_funct_copies_df$sample_id == sample, 
                                      'Num.Functional.Copies']
      # Take the log2 + 1 of this value
      fnc_stat <- log2(fnc_stat + 1)
    }
    
    altern_stats <- list(NA,NA,NA,NA,NA,NA)
    
    #' Helper function that will get stats for TP53 or PIK3CA if we are adding them
    #' as covariates as well 
    get_altern_stats <- function(regprot_i_ensg, cna_df_sub, cna_df_altern, bucket_type,
                                 mutation_regprot_df, sample, dataset, mut_stat) {
      
      altern_cna_stat <- get_cna_stat(cna_df_altern, sample, bucket_type, dataset)
      mutation_regprot_df_sub <- mutation_regprot_df[Patient %like% sample]
      if(mutation_regprot_df_sub[, .N] > 0) {
        altern_mut_stat <- 1
      } else {altern_mut_stat <- 0}
      if ((mut_stat == 1) & (altern_mut_stat == 1)) {
        altern_comb_stat <- 1
      } else {altern_comb_stat <- 0}
      
      name <- ""
      if(regprot_i_ensg == "ENSG00000141510") {name <- "PIK3CA."}
      if(regprot_i_ensg == "ENSG00000121879") {name <- "TP53."}
      list_to_ret <- list(altern_cna_stat, altern_mut_stat, altern_comb_stat)
      names(list_to_ret) <- c(paste0(name, "altern_cna_stat"), paste0(name, "altern_mut_stat"),
                              paste0(name, "altern_comb_stat"))
      return(list_to_ret)
    }
    
    if (incl_nextMutDriver) {
      if(regprot_i_ensg == "ENSG00000141510") {
        cna_b <- cna_bucketing
        if(grepl("just", cna_bucketing)) {cna_b <- "bucket_justAmp"}
        altern_stats <- get_altern_stats(regprot_i_ensg, cna_df_sub, cna_df_pik3ca,
                                         cna_b, mutation_regprot_df_pik3ca,
                                         sample, dataset, mut_stat)
        #altern_stats <- get_altern_stats(regprot_i_ensg, cna_df_sub, cna_df_pik3ca,
                                         #"bucket_justAmp", mutation_regprot_df_pik3ca,
                                         #sample, dataset, mut_stat)
      }
      else if(regprot_i_ensg  == "ENSG00000121879") {
        cna_b <- cna_bucketing
        if(grepl("just", cna_bucketing)) {cna_b <- "bucket_justDel"}
        altern_stats <- get_altern_stats(regprot_i_ensg, cna_df_sub, cna_df_tp53,
                                         cna_b, mutation_regprot_df_tp53,
                                         sample, dataset, mut_stat)
        #altern_stats <- get_altern_stats(regprot_i_ensg, cna_df_sub, cna_df_tp53,
                                         #"bucket_justDel", mutation_regprot_df_tp53,
                                         #sample, dataset, mut_stat)
      }
      else {
        cna_b_tp53 <- cna_bucketing
        cna_b_pik3ca <- cna_bucketing
        if(grepl("just", cna_bucketing)) {
          cna_b_pik3ca <- "bucket_justAmp"
          cna_b_tp53 <- "bucket_justDel"
        }
        
        altern_stats_tp53 <- get_altern_stats(regprot_i_ensg, cna_df_sub, cna_df_tp53,
                                         cna_b_tp53, mutation_regprot_df_tp53,
                                         sample, dataset, mut_stat)
        altern_stats_pik3ca <- get_altern_stats(regprot_i_ensg, cna_df_sub, cna_df_pik3ca,
                                         cna_b_pik3ca, mutation_regprot_df_pik3ca,
                                         sample, dataset, mut_stat)
        #altern_stats_tp53 <- get_altern_stats(regprot_i_ensg, cna_df_sub, cna_df_tp53,
                                         #"bucket_justDel", mutation_regprot_df_tp53,
                                         #sample, dataset, mut_stat)
        #altern_stats_pik3ca <- get_altern_stats(regprot_i_ensg, cna_df_sub, cna_df_pik3ca,
                                         #"bucket_justAmp", mutation_regprot_df_pik3ca,
                                         #sample, dataset, mut_stat)
        altern_stats <- c(altern_stats_tp53, altern_stats_pik3ca)
      }
    }
    
    if (debug) {
      print(paste("Meth stat:", meth_stat))
      print(paste("CNA stat:", cna_stat))
      print(paste("Mut stat:", mut_stat))
      if(useNumFunctCopies) {print(paste("FNC stat:", fnc_stat))}
      print("Altern driver stats:")
      print(altern_stats)
    }
    
    outdf <- data.table()
    if(useNumFunctCopies) {
      outdf[, 'FNCStat_i'] <- fnc_stat
      if(!meth_bucketing) {
        outdf[, 'MethStat_i'] <- meth_stat
      } else {
        outdf <- cbind(outdf, data.table("MethStat_i_b1" = meth_stat[1],
                                         "MethStat_i_b2" = meth_stat[2], 
                                         "MethStat_i_b3" = meth_stat[3]))
      }
      if(incl_nextMutDriver) {
        if(!(regprot_i_ensg == "ENSG00000121879")) {
          outdf[, 'TP53_MutStat_i'] <- altern_stats[['TP53.altern_mut_stat']]
          outdf[, 'TP53_DelStat_i'] <- altern_stats[['TP53.altern_cna_stat']]
          outdf[, 'TP53PIK3CA_MutStat_i'] <- altern_stats[['TP53.altern_comb_stat']]
        } 
        if (!(regprot_i_ensg == "ENSG00000141510")) {
          outdf[, 'PIK3CA_MutStat_i'] <- altern_stats[['PIK3CA.altern_mut_stat']]
          outdf[, 'PIK3CA_AmplStat_i'] <- altern_stats[['PIK3CA.altern_cna_stat']]
          outdf[, 'PIK3CAP53_MutStat_i'] <- altern_stats[['PIK3CA.altern_comb_stat']]
        }
      }
      
    } else {
      if(!meth_bucketing) {
        if ((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing))) {
          outdf <- data.table("MutStat_i" = mut_stat, "CNAStat_i" = cna_stat, "MethStat_i" = meth_stat)
        } else {
          outdf <- data.table("MutStat_i" = mut_stat, "CNAStat_i_b1" = cna_stat[1], "CNAStat_i_b2" = cna_stat[2],
                              "CNAStat_i_b3" = cna_stat[3], "MethStat_i" = meth_stat)
        }
        if(incl_nextMutDriver) {
          if(!(regprot_i_ensg == "ENSG00000121879")) {
            outdf[, 'TP53_MutStat_i'] <- altern_stats[['TP53.altern_mut_stat']]
            outdf[, 'TP53_DelStat_i'] <- altern_stats[['TP53.altern_cna_stat']]
            outdf[, 'TP53PIK3CA_MutStat_i'] <- altern_stats[['TP53.altern_comb_stat']]
          } 
          if (!(regprot_i_ensg == "ENSG00000141510")) {
            outdf[, 'PIK3CA_MutStat_i'] <- altern_stats[['PIK3CA.altern_mut_stat']]
            outdf[, 'PIK3CA_AmplStat_i'] <- altern_stats[['PIK3CA.altern_cna_stat']]
            outdf[, 'PIK3CAP53_MutStat_i'] <- altern_stats[['PIK3CA.altern_comb_stat']]
          }
        }

      } else {
        if ((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing))) {
          outdf <- data.table("MutStat_i" = mut_stat, "CNAStat_i" = cna_stat, "MethStat_i_b1" = meth_stat[1],
                              "MethStat_i_b2" = meth_stat[2], "MethStat_i_b3" = meth_stat[3])
        } else {
          outdf <- data.table("MutStat_i" = mut_stat, "CNAStat_i_b1" = cna_stat[1], "CNAStat_i_b2" = cna_stat[2],
                              "CNAStat_i_b3" = cna_stat[3], "MethStat_i_b1" = meth_stat[1], "MethStat_i_b2" = meth_stat[2], 
                              "MethStat_i_b3" = meth_stat[3])
        }
      }
      if(incl_nextMutDriver) {
        if(!(regprot_i_ensg == "ENSG00000121879")) {
          outdf[, 'TP53_MutStat_i'] <- altern_stats[['TP53.altern_mut_stat']]
          outdf[, 'TP53_DelStat_i'] <- altern_stats[['TP53.altern_cna_stat']]
          outdf[, 'TP53PIK3CA_MutStat_i'] <- altern_stats[['TP53.altern_comb_stat']]
        } 
        if (!(regprot_i_ensg == "ENSG00000141510")) {
          outdf[, 'PIK3CA_MutStat_i'] <- altern_stats[['PIK3CA.altern_mut_stat']]
          outdf[, 'PIK3CA_AmplStat_i'] <- altern_stats[['PIK3CA.altern_cna_stat']]
          outdf[, 'PIK3CAP53_MutStat_i'] <- altern_stats[['PIK3CA.altern_comb_stat']]
        }
      }
    }
    return(outdf)
  })
  
  
  # Bind these rows into the starter DF
  if(debug) {print(head(regprot_rows))}
  tryCatch({
    colnam <- colnames(regprot_rows[[1]])
    regprot_i_df <- data.table::rbindlist(regprot_rows)
    colnames(regprot_i_df) <- colnam
    starter_df <- cbind(patient_df, regprot_i_df)
    
  }, error=function(cond){
    print("Unexpected error; regprot rows is not valid: ")
    print(regprot_rows)
    return(starter_df)
  })
  
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
#' @param log_expression a TRUE/ FALSE value indicating whether we take the log2(exp + 1),
#' or simply the expression value given in the DF
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. Possible
#' values are "bucketInclAmp", "bucketExclAmp", and "rawCNA"
#' @param meth_bucketing a TRUE/FALSE value indicating whether or not we are bucketing
#' methylation values
#' @param analysis_type a string label that reads either "eQTL" or "meQTL" to 
#' determine what kind of model we should run
#' @param tumNormMatched a TRUE/FALSE value indicating whether or not the analysis is 
#' tumor-normal matched
#' @param debug a TRUE/FALSE value indicating if we are in debug mode and should 
#' use additional prints
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3', to denote what columns
#' we are referencing/ specific data types we are using
fill_targ_inputs <- function(starter_df, targ_k, targ_k_ensg, mutation_targ_df,
                             methylation_df, cna_df, expression_df, log_expression, cna_bucketing,
                             meth_bucketing, analysis_type, tumNormMatched, debug, dataset) {
  
  if(debug) {
    print(paste("Target:", targ_k))
    print(paste("Target ENSG:", targ_k_ensg))
  }
  
  # Begin by limiting the data frames to just this target gene k
  if (!(targ_k == "" | targ_k_ensg == "")) {
    mutation_targ_df <- filter_mut_by_uniprot(mutation_targ_df, targ_k) 
    cna_df <- filter_cna_by_ensg(cna_df, targ_k_ensg)
    methylation_df <- filter_meth_by_ensg(methylation_df, targ_k_ensg)
    exp_rows_to_keep <- unique(unlist(lapply(targ_k_ensg, function(x) grep(x, expression_df$ensg_id))))
    expression_df <- expression_df[exp_rows_to_keep,]
    
    if(debug) {
      print(paste("dim mutation targ df:", dim(mutation_targ_df)))
      print(paste("dim cna df:", dim(cna_df)))
      print(paste("dim methylation df:", dim(methylation_df)))
      print(paste("dim expression df:", dim(expression_df)))
    }
    
    # Continue only if all of these data frames contain information about this target gene k
    if(!(expression_df[, .N] == 0 | cna_df[, .N] == 0 | methylation_df[, .N] == 0 | 
         mutation_targ_df[, .N] == 0)) {
      
      # Loop through all samples to get the info for each target t_k column
      target_rows <- mclapply(1:starter_df[, .N], function(i) {
        sample <- starter_df$sample_id[i]
        
        # Is this target mutated? 
        mut_stat <- get_mut_stat_targ(mutation_targ_df, sample)
        
        # Is this target amplified or deleted?
        cna_stat <- get_cna_stat(cna_df, sample, cna_bucketing, dataset)
        if(grepl("incl", cna_bucketing)) {cna_stat <- cna_stat[[1]]}    
        
        # Is this target methylated, or differentially methylated?
        if(analysis_type == "meQTL") {meth_bucketing <- FALSE}
        meth_stat <- get_meth_stat(methylation_df, sample, meth_bucketing, tumNormMatched)
        if(meth_bucketing) {
          meth_stat <- as.integer(meth_stat[[1]])
          if(length(meth_stat) > 3) {meth_stat <- meth_stat[1:3]}
        }
        # What is the expression of this target in this sample in cancer?
        exp_stat <- unlist(get_exp_stat(expression_df, sample, tumNormMatched, 
                                        dataset, log_expression))
        
        if(debug) {
          print(paste("exp_stat:", exp_stat))
          print(paste("mut_stat:", mut_stat))
          print(paste("cna_stat:", cna_stat))
          print(paste("meth_stat:", meth_stat))
        }
        
        # Make row for this patient and return
        if(!meth_bucketing) {
          if ((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing))) {
            return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, "CNAStat_k" = cna_stat, 
                              "MethStat_k" = meth_stat))
          } else {
            return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, "CNAStat_k_b1" = cna_stat[1], 
                              "CNAStat_k_b2" = cna_stat[2], "CNAStat_k_b3" = cna_stat[3], "MethStat_k" = meth_stat))
          }
        } else {
          if ((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing))) {
            return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, "CNAStat_k" = cna_stat, 
                              "MethStat_k_b1" = meth_stat[1], "MethStat_k_b2" = meth_stat[2], "MethStat_k_b3" = meth_stat[3]))
          } else {
            return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, "CNAStat_k_b1" = cna_stat[1], 
                              "CNAStat_k_b2" = cna_stat[2], "CNAStat_k_b3" = cna_stat[3], "MethStat_k_b1" = meth_stat[1], 
                              "MethStat_k_b2" = meth_stat[2], "MethStat_k_b3" = meth_stat[3]))
          }
        }
      })
      
      tryCatch({
        # Bind these rows into the starter DF
        colnam <- colnames(target_rows[[1]])
        targ_k_df <- data.table::rbindlist(target_rows)
        colnames(targ_k_df) <- colnam
        full_df <- cbind(starter_df, targ_k_df)
        
      }, error=function(cond){
        print("Unexpected error; target rows is not valid: ")
        print(target_rows)
        full_df <- starter_df
      })
      
      return(full_df)
    }
  } else {
    if(debug) {
      print(paste("Target not found in all DFs:", targ_k))
    }
    return(NA)
  }
}


############################################################
#### CONSTRUCT LM FORMULA
############################################################
#' Given a linear model input table, analysis and cancer type,
#' contructs a formula to pass to the speedlm formula. Removes 
#' unnecessary covariates, etc.
#' @param lm_input_table the LM input table produced from the above function
#' @param analysis_type 'eQTL' or 'meQTL'
#' @param num_PEER a value from 0-10 indicating the number of PEER factors we 
#' are including as covariates in the model
#' @param num_pcs a value from 0-2 indicating the number of  
#' principal components we are including as covariates in the model 
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. Possible
#' values are "bucketInclAmp", "bucketExclAmp", and "rawCNA"
#' @param meth_bucketing a TRUE/FALSE value indicating whether or not we are bucketing
#' methylation values
#' @param covs_to_incl a vector, either c("ALL") to include all covariates, or a vector
#' of all the covariates we want to keep
construct_formula <- function(lm_input_table, analysis_type, num_PEER, num_pcs,
                              cna_bucketing, meth_bucketing, covs_to_incl) {
  colnames_to_incl <- colnames(lm_input_table)[2:ncol(lm_input_table)]

  # Note: currently excluding age b1 and age b2 since no patients fall into these categories
  colnames_to_incl <- colnames_to_incl[!((colnames_to_incl == "Age_b1") | (colnames_to_incl == "Age_b2"))]
  
  # Remove the race variable entirely
  colnames_to_incl <- colnames_to_incl[!(grepl("Race", colnames_to_incl))]
  
  # Remove the PEER factors we are not using
  peer_factors <- colnames_to_incl[grepl("PEER", colnames_to_incl)]
  peer_fact_to_remove <- setdiff(peer_factors, peer_factors[0:num_PEER])
  colnames_to_incl <- colnames_to_incl[!(colnames_to_incl %fin% peer_fact_to_remove)]
  
  # Remove the PCs we are not using
  pcs <- colnames_to_incl[grepl("PC", colnames_to_incl)]
  pcs_to_remove <- setdiff(pcs, pcs[0:num_pcs])
  colnames_to_incl <- colnames_to_incl[!(colnames_to_incl %fin% pcs_to_remove)]  
  
  # For the bucketed variables, remove the last bucket (we don't need it!)
  bucketed_vars <- c("Tot_Mut_b","Tumor_purity_b", "Tot_IC_Frac_b", "Cancer_type_b")
  if(meth_bucketing) {
    bucketed_vars <- c(bucketed_vars, "MethStat_k")
    if(analysis_type == "eQTL") {bucketed_vars <- c(bucketed_vars, "MethStat_i")}
  }
  if(!((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing)))) {
    bucketed_vars <- c(bucketed_vars, c("CNAStat_k", "CNAStat_i"))
  }
  
  vals_to_remove <- unlist(lapply(bucketed_vars, function(var) {
    # Get the last bucket for this variable
    matching_vars <- colnames_to_incl[grepl(var, colnames_to_incl)]
    index <- which(colnames_to_incl == matching_vars[length(matching_vars)])
    return(index)
  }))
  colnames_to_incl <- colnames_to_incl[-vals_to_remove]

  
  # Additionally, if we have any columns that are uniform (the same value in all patients) we'd 
  # like to remove those as well, as long as they are not MutStat_i, CNAStat_i or FNCStat_i
  uniform_val_ind <- which(sapply(lm_input_table, function(x) length(unique(x)) == 1))
  uniform_vals <- colnames(lm_input_table)[uniform_val_ind]
  uniform_vals <- uniform_vals[!((uniform_vals == "MutStat_i") | (uniform_vals == "CNAStat_i") | 
                                   (uniform_vals == "FNCStat_i"))]
  colnames_to_incl <- colnames_to_incl[!(colnames_to_incl %fin% uniform_vals)]
  
  # Exclude the library size, as we are using this as an offset instead
  colnames_to_incl <- colnames_to_incl[!(colnames_to_incl == "Lib_Size")]
  
  # Now, of the variables we have left, include only those given in the covs_to_include parameter
  if(!(covs_to_incl[1] == "ALL")) {
    colnames_to_incl <- unlist(lapply(colnames_to_incl, function(x) {
      if(any(unlist(lapply(covs_to_incl, function(y) grepl(y, x))))) {return(x)}
    }))
  }

  if (analysis_type == "eQTL") {      
    colnames_to_incl <- colnames_to_incl[!(colnames_to_incl == "ExpStat_k")]
    formula <- paste(colnames_to_incl, collapse = " + ") 
    formula <- paste("ExpStat_k ~ ", formula, sep = "")
  } else if (analysis_type == "meQTL") {
    colnames_to_incl <- colnames_to_incl[!(colnames_to_incl == "MethStat_k")]
    formula <- paste(colnames_to_incl, collapse = " + ") 
    formula <- paste("MethStat_k ~ ", formula, sep = "")
  } else {
    print(paste("Invalid analysis type:", analysis_type))
    return(NA)
  }
  return(formula)
}


############################################################
#### CREATE OUTPATH AND OUTFILE NAME 
############################################################
# Create an appropriate file outpath
if(args$cancerType == "PanCancer") {
  if(!(args$specificTypes == "ALL")) {
    cancer_types <- unlist(strsplit(args$specificTypes, ";", fixed = TRUE))
    outpath <- lapply(cancer_types, function(ct) {
      return(create_file_outpath(cancerType = args$cancerType, specificType = ct, test = test, 
                                 tester_name = args$tester_name, run_name = args$run_name,
                                 tumNormMatched = tumNormMatched, QTLtype = args$QTLtype,
                                 dataset = args$dataset))
    })
    names(outpath) <- cancer_types
  } else {
    outpath <- create_file_outpath(cancerType = args$cancerType, specificType = "", test = test, 
                                   tester_name = args$tester_name, run_name = args$run_name,
                                   tumNormMatched = tumNormMatched, QTLtype = args$QTLtype,
                                   dataset = args$dataset)
  }
} else {
  outpath <- create_file_outpath(cancerType = args$cancerType, specificType = "", test = test, 
                                 tester_name = args$tester_name, run_name = args$run_name,
                                 tumNormMatched = tumNormMatched, QTLtype = args$QTLtype,
                                 dataset = args$dataset)
}


if(debug) {
  print(paste("Outpath:", outpath))
}


# Create an appropriate output file name
outfn <- create_output_filename(test = test, tester_name = args$tester_name, run_name = args$run_name,
                                targets_name = args$targets_name, expression_df_name = args$expression_df,
                                cna_bucketing = args$cna_bucketing, mutation_regprot_df_name = args$mutation_regprot_df, 
                                meth_bucketing = meth_bucketing, meth_type = args$meth_type, 
                                patient_df_name = args$patient_df, num_PEER = args$num_PEER, 
                                num_pcs = args$num_pcs, randomize = randomize, covs_to_incl_label = args$select_args_label,
                                patients_to_incl_label = args$patientsOfInterestLabel, removeCis = removeCis,
                                removeMetastatic = removeMetastatic)

if(debug) {
  print(paste("Outfile Name:", outfn))
}


############################################################
# IF RUNNING PER-CANCER WITH SPECIFIC CANCER TYPES, GET 
# THE PATIENTS OF THOSE CANCER TYPES 
############################################################
#' If we are running per-cancer on several individual cancer types, get a mapping
#' between each cancer type and the IDs of the patients classified as that cancer type
#' @param specificTypes a semicolon-separated string of the cancer types we are looking
#' at in this run
#' @param clinical_df a file name for a clinical data frame that can link patient IDs to
#' their respective cancer types
get_patient_cancer_mapping <- function(specificTypes, clinical_df) {
    specific_types <- unlist(strsplit(specificTypes, ";", fixed = TRUE))
    clinical_df <- fread(paste(main_path, paste("clinical/", clinical_df, sep = ""), sep = ""), 
                         header = TRUE)
    patient_cancer_mapping <- lapply(specific_types, function(ct) {
      pats <- clinical_df[grepl(ct, clinical_df$project_id),'case_submitter_id']
      pats_ids <- unlist(lapply(pats$case_submitter_id, function(pat) 
        unlist(strsplit(pat, "-", fixed = TRUE))[3]))
      return(pats_ids)
    })
    names(patient_cancer_mapping) <- specific_types
    return(patient_cancer_mapping)
}

# Call function, if applicable
patient_cancer_mapping <- NA
if(args$cancerType == "PanCancer") {
  if(!(args$specificTypes == "ALL")) {
    patient_cancer_mapping <- get_patient_cancer_mapping(specificTypes = args$specificTypes,
                                                         clinical_df = args$clinical_df)
    if(debug) {print(patient_cancer_mapping)}
  }
}


############################################################
#### CALL FUNCTION
############################################################
# Run gc to free up any loose memory
gc()

# Convert covariates to include from a semicolon-separated character to a vector
tryCatch({
  covs_to_incl <- unlist(strsplit(args$select_args, ";", fixed = TRUE))
}, error = function(cond) {
  print(cond)
  print("Invalid input covariate list. Make sure covariates are semi-colon separated. Running with all covariates.")
  covs_to_incl <- c("ALL")
})

# Run the function itself

# Case 1: If we're running on just one cancer type or PanCancer
if(is.na(patient_cancer_mapping)) {
  
  # If needed, subset these data tables by the given patient IDs of interest using helper functions
  if(patients_of_interest != "") {
    patient_df <- subset_by_intersecting_ids(patients_of_interest, patient_df, FALSE, tumNormMatched, args$dataset)
    mutation_targ_df <- subset_by_intersecting_ids(patients_of_interest, mutation_targ_df, TRUE, tumNormMatched, args$dataset)
    mutation_regprot_df <- subset_regprot_df_by_intersecting_ids(patients_of_interest, mutation_regprot_df, tumNormMatched, args$dataset)
    methylation_df <- subset_by_intersecting_ids(patients_of_interest, methylation_df, TRUE, tumNormMatched, args$dataset)
    if(!is.na(methylation_df_meQTL)) {
      methylation_df_meQTL <- subset_by_intersecting_ids(patients_of_interest, methylation_df_meQTL, TRUE, tumNormMatched, args$dataset)
    } else {methylation_df_meQTL <- NA} 
    cna_df <- subset_by_intersecting_ids(patients_of_interest, cna_df, TRUE, tumNormMatched, args$dataset)
    expression_df <- subset_by_intersecting_ids(patients_of_interest, expression_df, TRUE, tumNormMatched, args$dataset)
    if(useNumFunctCopies) {num_funct_copies_DF <- subset_by_intersecting_ids(patients_of_interest, num_funct_copies_DF, FALSE, tumNormMatched, args$dataset)}
  }
  
  # If needed, exclude metastatic samples
  if(removeMetastatic == TRUE) {
    patient_df <- remove_metastatic_samples(patient_df, FALSE)
    mutation_targ_df <- remove_metastatic_samples(mutation_targ_df, TRUE)
    mutation_regprot_df <- remove_metastatic_samples_regprot(mutation_regprot_df)
    methylation_df <- remove_metastatic_samples(methylation_df, TRUE)
    if(!is.na(methylation_df_meQTL)) {
      methylation_df_meQTL <- remove_metastatic_samples(methylation_df_meQTL, TRUE)
    } else {methylation_df_meQTL <- NA} 
    cna_df <- remove_metastatic_samples(cna_df, TRUE)
    expression_df <- remove_metastatic_samples(expression_df, TRUE)
    if(useNumFunctCopies) {num_funct_copies_DF <- remove_metastatic_samples(num_funct_copies_DF, FALSE)}
  }
    
  master_df <- run_linear_model(protein_ids_df = protein_ids_df, 
                                downstream_target_df = targets_DF, 
                                patient_df = patient_df,
                                mutation_df_targ =  mutation_targ_df,
                                mutation_df_regprot = mutation_regprot_df, 
                                methylation_df = methylation_df, 
                                methylation_df_meQTL = methylation_df_meQTL,
                                cna_df = cna_df,
                                expression_df = expression_df, 
                                is_rank_or_quant_norm = is_rank_or_quant_norm,
                                log_expression = log_expression,
                                analysis_type = args$QTLtype,
                                tumNormMatched = tumNormMatched,
                                randomize = randomize,
                                cna_bucketing = args$cna_bucketing,
                                meth_bucketing = meth_bucketing,
                                num_PEER = args$num_PEER,
                                num_pcs = args$num_pcs,
                                debug = debug,
                                collinearity_diagn = collinearity_diagn,
                                outpath = outpath, outfn = outfn,
                                regularization = args$regularization,
                                covs_to_incl = covs_to_incl,
                                removeCis = removeCis,
                                all_genes_id_conv = all_genes_id_conv,
                                useNumFunctCopies = useNumFunctCopies,
                                num_funct_copies_df = num_funct_copies_DF,
                                incl_nextMutDriver = incl_nextMutDriver,
                                dataset = args$dataset,
                                inclResiduals = inclResiduals)
  
# Case 2: If we're running on multiple cancer types separately
} else if (!(is.na(patient_cancer_mapping))) {
  # Handle the case of running multiple cancer types per-cancer
  # Subset all the input files for each cancer type and run separately
  master_df_list <- lapply(1:length(outpath), function(i) {
    
    # Get the patients for this cancer type
    patient_ids <- unique(unlist(patient_cancer_mapping[[i]]))
    print(patient_ids)
    
    if(patients_of_interest != "") {
      patient_ids <- patient_ids[patient_ids %fin% patients_of_interest]
    }
    
    # Get the outpath for this cancer type
    outpath_i <- outpath[[i]]
    
    # Subset files using patient_ids (using helper function)
    patient_df_sub <- subset_by_intersecting_ids(patient_ids, patient_df, FALSE, tumNormMatched, args$dataset)
    mutation_targ_df_sub <- subset_by_intersecting_ids(patient_ids, mutation_targ_df, TRUE, tumNormMatched, args$dataset)
    mutation_regprot_df_sub <- subset_regprot_df_by_intersecting_ids(patient_ids, mutation_regprot_df, tumNormMatched, args$dataset)
    methylation_df_sub <- subset_by_intersecting_ids(patient_ids, methylation_df, TRUE, tumNormMatched, args$dataset)
    if(!is.na(methylation_df_meQTL)) {
      methylation_df_meQTL_sub <- subset_by_intersecting_ids(patient_ids, methylation_df_meQTL, TRUE, tumNormMatched, args$dataset)
    } else {methylation_df_meQTL_sub <- NA} 
    cna_df_sub <- subset_by_intersecting_ids(patient_ids, cna_df, TRUE, tumNormMatched, args$dataset)
    expression_df_sub <- subset_by_intersecting_ids(patient_ids, expression_df, TRUE, tumNormMatched, args$dataset)
    if(useNumFunctCopies) {num_funct_copies_DF_sub <- subset_by_intersecting_ids(patients_of_interest, num_funct_copies_DF, FALSE, tumNormMatched, args$dataset)}
    
    # If needed, exclude metastatic samples
    if(removeMetastatic == TRUE) {
      patient_df <- remove_metastatic_samples(patient_df, FALSE)
      mutation_targ_df <- remove_metastatic_samples(mutation_targ_df, TRUE)
      mutation_regprot_df <- remove_metastatic_samples_regprot(mutation_regprot_df)
      methylation_df <- remove_metastatic_samples(methylation_df, TRUE)
      if(!is.na(methylation_df_meQTL)) {
        methylation_df_meQTL <- remove_metastatic_samples(methylation_df_meQTL, TRUE)
      } else {methylation_df_meQTL <- NA} 
      cna_df <- remove_metastatic_samples(cna_df, TRUE)
      expression_df <- remove_metastatic_samples(expression_df, TRUE)
      if(useNumFunctCopies) {num_funct_copies_DF <- remove_metastatic_samples(num_funct_copies_DF, FALSE)}
    }
    
    # Run the model with these subsetted files
    master_df <- run_linear_model(protein_ids_df = protein_ids_df, 
                                  downstream_target_df = targets_DF, 
                                  patient_df = patient_df_sub,
                                  mutation_df_targ =  mutation_targ_df_sub,
                                  mutation_df_regprot = mutation_regprot_df_sub, 
                                  methylation_df = methylation_df_sub, 
                                  methylation_df_meQTL = methylation_df_meQTL_sub,
                                  cna_df = cna_df_sub,
                                  expression_df = expression_df_sub, 
                                  is_rank_or_quant_norm = is_rank_or_quant_norm,
                                  log_expression = log_expression,
                                  analysis_type = args$QTLtype,
                                  tumNormMatched = tumNormMatched,
                                  randomize = randomize,
                                  cna_bucketing = args$cna_bucketing,
                                  meth_bucketing = meth_bucketing,
                                  num_PEER = args$num_PEER,
                                  num_pcs = args$num_pcs,
                                  debug = debug,
                                  collinearity_diagn = collinearity_diagn,
                                  outpath = outpath, outfn = outfn,
                                  regularization = args$regularization,
                                  covs_to_incl = covs_to_incl,
                                  removeCis = removeCis,
                                  all_genes_id_conv = all_genes_id_conv,
                                  useNumFunctCopies = useNumFunctCopies,
                                  num_funct_copies_df = num_funct_copies_DF,
                                  incl_nextMutDriver = incl_nextMutDriver,
                                  dataset = args$dataset,
                                  inclResiduals = inclResiduals)
    
    return(master_df)
  })

} else {
  print("Something went wrong with patient-cancer type mapping. Exiting now.")
  master_df <- NA
}



############################################################
#### INITIAL PROCESSING OF RAW FILE & WRITING
############################################################
#' Takes a master data frame produced from the main linear model function and 
#' does some preliminary processing, including ordering by p-value, limiting to
#' just MutStat_i and CNAStat_i terms of interest, and running collinearity 
#' diagnostics if desired
#' @param master_df a master data frame output result from main LM function
#' @param outpath the outpath for the master DF to be written to
#' @param outfn the generic output filename for the master DF
#' @param collinearity_diagn a TRUE/FALSE value indicating whether we are running
#' additional collinearity diagnostics
#' @param debug a TRUE/FALSE value indicating whether or not we are in debug mode
#' @param randomize a TRUE/FALSE value indicating whether or not this was a 
#' randomized run
#' @param useNumFunctCopies a TRUE/FALSE value indicating whether or not we are
#' using the number of functional copies of regprot_i rather than mut/CNA status
process_raw_output_df <- function(master_df, outpath, outfn, collinearity_diagn, 
                                  debug, randomize, useNumFunctCopies) {
  # Order the file by p-value
  print(head(master_df))
  try(master_df <- master_df[order(p.value)])
  #try(master_df <- setorder(master_df, cols = "p.value", order = -1, na.last = TRUE))  
  
  # Write the results to the given file, making the directory if it does not already exist
  outpath <- as.character(unlist(outpath[[1]]))
  print(outpath)
  tryCatch({
    dir.create(outpath, showWarnings = FALSE)
  }, error = function(cond) {
    print(cond)
    print("Invalid outpath. Cannot create specified directory.")
  })
  fwrite(master_df, paste(outpath, paste(outfn, ".csv", sep = ""), sep = "/"))
  
  # Limit the data frame to just the term of interest (typically either MutStat_i or CNAStat_i)
  if(!useNumFunctCopies) {
    master_df_mut <- master_df[master_df$term == "MutStat_i",]
    master_df_cna <- master_df[grepl("CNAStat_i", master_df$term),]
    
    # Write these to files as well 
    fwrite(master_df_mut, paste(outpath, paste(outfn, "_MUT.csv", sep = ""), sep = "/")) 
    fwrite(master_df_cna, paste(outpath, paste(outfn, "_CNA.csv", sep = ""), sep = "/"))  
  } else {
    master_df_fnc <- master_df[master_df$term == "FNCStat_i",]
    fwrite(master_df_fnc, paste(outpath, paste(outfn, "_FNC.csv", sep = ""), sep = "/")) 
  }
  
  
  # Combine the collinearity results and write to a file
  if(collinearity_diagn) {
    collinearity_df <- combine_collinearity_diagnostics(outpath, outfn, debug, randomize)
    # Adjust the output file name for collinearity results
    collin_fn <- str_replace(outfn, "output_results", "collinearity_results")
    if(debug) {print(paste("Collinearity Output FN:", collin_fn))}
    fwrite(collinearity_df, paste(outpath, paste(collin_fn, ".csv", sep = ""), sep = "/"))
  }
  
  return(master_df)
}


if(is.na(patient_cancer_mapping)) {
  outpath <- as.character(unlist(outpath[[1]]))
  outfn <- as.character(unlist(outfn[[1]]))
  master_df <- process_raw_output_df(master_df = master_df, outpath = outpath, outfn = outfn, 
                                     collinearity_diagn = collinearity_diagn, debug = debug,
                                     randomize = randomize, useNumFunctCopies = useNumFunctCopies)
  if(!useNumFunctCopies) {
    master_df_mut <- master_df[master_df$term == "MutStat_i",]
    master_df_cna <- master_df[grepl("CNAStat_i", master_df$term),]
    if((args$tester_name == "TP53") & incl_nextMutDriver) {
      master_df_pik3ca_mut <- master_df[master_df$term == "PIK3CA_MutStat_i",]
      master_df_pik3ca_cna <- master_df[master_df$term == "PIK3CA_AmplStat_i",]
      master_df_p53_pik3ca_mut <- master_df[master_df$term == "PIK3CAP53_MutStat_i",]
    }    
    if((args$tester_name == "PIK3CA") & incl_nextMutDriver) {
      master_df_tp53_mut <- master_df[master_df$term == "TP53_MutStat_i",]
      master_df_tp53_cna <- master_df[master_df$term == "TP53_AmplStat_i",]
      master_df_tp53_pik3ca_mut <- master_df[master_df$term == "TP53PIK3CA_MutStat_i",]
    } 
  } else {
    master_df_fnc <- master_df[master_df$term == "FNCStat_i",]
  }
  
  # Call the file to create output visualizations
  source(paste(source_path, "process_LM_output.R", sep = "")) 
  
} else {
  for (i in 1:length(master_df_list)) {
    master_df <- master_df_list[[i]]
    outpath <- as.character(unlist(outpath[[i]]))
    outfn <- as.character(unlist(outfn[[i]]))
    master_df <- process_raw_output_df(master_df = master_df, outpath = outpath, outfn = outfn, 
                                       collinearity_diagn = collinearity_diagn, debug = debug,
                                       randomize = randomize, useNumFunctCopies = useNumFunctCopies)
    if(!useNumFunctCopies) {
      master_df_mut <- master_df[master_df$term == "MutStat_i",]
      master_df_cna <- master_df[grepl("CNAStat_i", master_df$term),]
      if((args$tester_name == "TP53") & incl_nextMutDriver) {
        master_df_pik3ca_mut <- master_df[master_df$term == "PIK3CA_MutStat_i",]
        master_df_pik3ca_cna <- master_df[master_df$term == "PIK3CA_AmplStat_i",]
        master_df_tp53_pik3ca_mut <- master_df[master_df$term == "PIK3CAP53_MutStat_i",]
      }    
      if((args$tester_name == "PIK3CA") & incl_nextMutDriver) {
        master_df_tp53_mut <- master_df[master_df$term == "TP53_MutStat_i",]
        master_df_tp53_cna <- master_df[master_df$term == "TP53_AmplStat_i",]
        master_df_tp53_pik3ca_mut <- master_df[master_df$term == "TP53PIK3CA_MutStat_i",]
      }    
    } else {
      master_df_fnc <- master_df[master_df$term == "FNCStat_i",]
    }
    
    
    # Call the file to create output visualizations
    source(paste(source_path, "process_LM_output.R", sep = "")) 
  }
}


