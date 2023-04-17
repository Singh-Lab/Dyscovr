############################################################
### Create Linear Model Input Table
### Written By: Sara Geraghty, December 2022
############################################################

# Creates an input table for linear regression in linear_model.R for each gene target
# provided in input tables. Writes all input tables to a directory, where they can be 
# fetched when LM runs.

# This script requires 1 node and 6 CPUs for multithreading.

#!/usr/bin/env Rscript

# Set the library path
library.path <- .libPaths()
#library.path <- "C:/Users/sarae/AppData/Local/R/win-library/4.1"
#library.path <- "C:/Program Files/R/R-4.1.1/library"
#print(library.path)

#library(biomaRt, lib.loc = library.path)
library(parallel, lib.loc = library.path)
library(rlang, lib.loc = library.path)
library(dplyr, lib.loc = library.path)
library(broom, lib.loc = library.path)
library(data.table, lib.loc = library.path)
library(speedglm, lib.loc = library.path)
library(argparse, lib.loc = library.path)
library(stringr, lib.loc = library.path)
library(dqrng, lib.loc = library.path)
library(foreach, lib.loc = library.path)
library(snow, lib.loc = library.path)
library(doParallel, lib.loc = library.path)

# Source other files needed
source_path <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/Sara_LinearModel/"
source(paste(source_path, "general_important_functions2.R", sep = "/"))
source(paste(source_path, "linear_model_helper_functions2.R", sep = "/"))

############################################################
# SET UP PARSER ARGUMENTS
############################################################
# Create parser object
parser <- ArgumentParser()

# Specify desired options for the parser

# Data set we are running on
parser$add_argument("--dataset", default = "TCGA", type = "character", 
                    help = "TCGA, METABRIC, Chinese_TN, ICGC, or CPTAC3. Default is TCGA.")

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
parser$add_argument("--run_query_genes_jointly", default = "TRUE", type = "character",
                    help = "If TRUE, will run all the query regulatory proteins together in the same joint model, rather than running one model per query-target pair.")

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
#parser$add_argument("--mutation_regprot_df", default = "iprotein_results_missense_CancerOnly_IntersectPatients.csv", 
                    #type = "character",
                    #help = "The name of the mutation regulatory protein data frame input. [default %(default)s]")
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

# What type of test we are running, so we can write the results files to an appropriate directory. 
# Only necessary if --test is F.
parser$add_argument("--run_name", default = "top_drivers", type = "character",
                    help = "Provide a name for the run, if not testing on just one regulatory protein, to write output files to appropriate directory. [default %(default)s]") 

# Whether or not we are including PEER factors and PCs as covariates in the model
parser$add_argument("--num_PEER", default = 0, type = "integer",
                    help = "A value from 0 to 10 indicating the number of PEER factors we are using as covariates in our model. Default is 0.")
parser$add_argument("--num_pcs", default = 2, type = "integer",
                    help = "A value between 0 and 2 indicating the number of principal components we are using as covariates in our model. Default is 2.")
parser$add_argument("--mut_pc_vect", default = "0", type = "character",
                    help = "either '0' for no mutational vectors, or a vector of PCs to include, e.g. 1,3,4. Default is 0.")

# Add a flag for keeping only trans pairings
parser$add_argument("--removeCis", default = "TRUE", type = "character",
                    help = "A TRUE/ FALSE value to indicate whether or not we want to remove cis pairings and look only at trans pairings. Defaults to true, but can be set to false.")

# A flag for using the number of functional copies of regulatory protein per sample, rather than the mutation and CNA status 
parser$add_argument("--useNumFunctCopies", default = "FALSE", type = "character",
                    help = "A TRUE/ FALSE value to indicate whether or not we want to use the allele-specific number of functional copies for the given regulatory tester protein, rather than its separate mutation and CNA status.")
#parser$add_argument("--explan_var_interest", default = "MutStat", type = "character",
#help = "The name of the explantory variable we are interested in. If CNAStat, program will automatically add a covariate to account for the presence of a co-amplified or deleted neighboring driver gene.")

# Add a flag for debugging
parser$add_argument("--debug", default = "FALSE", type = "character",
                    help = "A TRUE/ FALSE value indicating whether or not we want detailed printing for debugging purposes. Default is FALSE.")


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
removeCis <- str2bool(args$removeCis)
removeMetastatic <- str2bool(args$removeMetastatic)
useNumFunctCopies <- str2bool(args$useNumFunctCopies)
incl_nextMutDriver <- str2bool(args$incl_nextMutDriver)
runRegprotsJointly <- str2bool(args$run_query_genes_jointly)


############################################################
# SET MAIN PATH AND IMPORT GENERALLY NEEDED FILES
############################################################
# Set the paths
input_file_path <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/"

if(args$cancerType == "BRCA") {
  targ_path <- paste(input_file_path, "BRCA/target_lists/", sep = "")
  
  if(args$dataset == "METABRIC") {
    main_path <- paste0(input_file_path, "METABRIC/")
    prot_path <- paste(input_file_path, "METABRIC/", sep = "")
    
  } else if (args$dataset == "Chinese_TN") {
    main_path <- paste0(input_file_path, "Chinese_TN/")
    prot_path <- paste(input_file_path, "Chinese_TN/", sep = "")
  
  } else {
    if(!tumNormMatched) {
      main_path <- paste(input_file_path, "BRCA/tumor_only/", sep = "")
    } else {
      main_path <- paste(input_file_path, "BRCA/tumor_normal_matched/", sep = "") 
    }
    prot_path <- paste(input_file_path, "BRCA/", sep = "")
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

if(debug) {print(main_path)}

# If we are looking at tester targets, we need to go one directory deeper
#if(grepl("curated", args$targets_name)) {
#  targ_path <- paste(targ_path, "curated_targets/", sep = "")
#} else if (grepl("chipeat", args$targets_name)) {
#  targ_path <- paste(targ_path, "chipeat_targets/", sep = "")
#} else {targ_path <- targ_path}

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
  if(args$cancerType == "PanCancer") {
    if(args$specificTypes == "ALL") {
      protein_ids_df <- fread(paste(prot_path, args$protein_ids_df, sep = ""), 
                              header = TRUE)
    } else {
      types <- unlist(strsplit(args$specificTypes, ";", fixed = T))
      protein_ids_df <- lapply(types, function(ct) {
        new_name <- paste0(prot_path, paste0(args$protein_ids_df, paste0(ct, ".csv")))
        return(fread(new_name, header = TRUE))
      })
      names(protein_ids_df) <- types
    }
  } else {
    protein_ids_df <- fread(paste(prot_path, args$protein_ids_df, sep = ""), 
                            header = TRUE)
  }

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

# Determine what type of expression this is from the file name
is_rank_or_quant_norm <- FALSE
if(grepl("rank", args$expression_df) | grepl("quantile", args$expression_df)) {is_rank_or_quant_norm <- TRUE}

# Determine whether or not we will need to log-transform the data, or whether this
# has already been done in some way (e.g. by an inverse normal transformation or z-score normalization)
log_expression <- TRUE
if(grepl("transform", args$expression_df) | grepl("zscore", args$expression_df) | 
   grepl("sklearn", args$expression_df)) {log_expression <- FALSE}


############################################################
# IMPORT METHYLATION FILES
############################################################
if(args$dataset != "Chinese_TN") {
  methylation_df <- fread(paste(main_path, paste("methylation/", args$methylation_df, sep = ""), sep = ""),
                          header = TRUE)
  #methylation_df <- methylation_df[,2:ncol(methylation_df)]
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
}

############################################################
# IMPORT GENE TARGET MUTATION FILES
############################################################
mutation_targ_df <- fread(paste(main_path, paste("genetarg_mutation/", args$mutation_targ_df, sep = ""), sep = ""),
                          header = TRUE)
#mutation_targ_df <- mutation_targ_df[,2:ncol(mutation_targ_df)]
colnames(mutation_targ_df)[which(colnames(mutation_targ_df) == "Gene_Symbol")] <- "gene_name"
if(!("Swissprot" %in% colnames(mutation_targ_df))) {
  mutation_targ_df <- mutation_targ_df[, Swissprot := unlist(lapply(mutation_targ_df$gene_name, function(x) 
    paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 'uniprot_gn_id'])), collapse = ";")))]
}

if(debug) {
  print("Mutation Targ DF")
  print(head(mutation_targ_df))
}


############################################################
# IMPORT REGULATORY PROTEIN MUTATION FILES
############################################################
#mutation_regprot_df <- fread(paste(main_path, paste("regprot_mutation/", args$mutation_regprot_df, sep = ""), sep = ""),
#                             header = TRUE)
#if(unlist(strsplit(args$mutation_regprot_df, "_", fixed = TRUE))[1] == "iprotein") {
#  mutation_regprot_df <- mutation_regprot_df[,5:ncol(mutation_regprot_df)]   # remove first few meaningless columns, if necessary
#}
#mutation_regprot_df <- mutation_regprot_df[!duplicated(mutation_regprot_df),] # remove any duplicate rows

#if(debug) {
#  print("Mutation Regprot DF")
#  print(head(mutation_regprot_df))
#}



############################################################
# IMPORT CNA FILES
############################################################
cna_df <- fread(paste(main_path, paste("cna/", args$cna_df, sep = ""), sep = ""),
                header = TRUE)
colnames(cna_df)[1] <- "ensg_id"

# Replace any -Inf values with NA
#cna_df[cna_df == -Inf] <- NA

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
targets_DF <- targets_DF[!((targets_DF$swissprot == "") | (targets_DF$ensg == "")),]

# Remove the first column, if needed (for all but metabolic targets)
if(ncol(targets_DF) > 2) {
  targets_DF <- targets_DF[,2:3]
}

# Regardless of method, limit targets to only those that overlap the genes in the expression 
# and methylation DFs
if(args$dataset != "Chinese_TN") {
  methylation_ensg_ids <- unlist(lapply(methylation_df$ensg_ids, function(ids) 
    unlist(strsplit(as.character(ids), ";", fixed = TRUE))))
  rows_to_keep <- unlist(lapply(1:length(targets_DF$ensg), function(i) {
    targs <- unlist(strsplit(as.character(targets_DF[i,'ensg']), ";", fixed = TRUE))
    if (any(targs %fin% expression_df$ensg_id)) {
      if(any(targs %fin% methylation_ensg_ids)) {
        return(TRUE)
      } else {return(FALSE)}
    } else {return(FALSE)}
  }))
  targets_DF <- targets_DF[rows_to_keep,]
}

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
# MAIN FUNCTION FOR CREATING LM INPUT DATA FRAME
############################################################
#' @param protein_ids_df DF of the regulatory proteins of interest, found based 
#' on their DNA-binding regions (columns for swissprot IDs and ENSG IDs)
#' @param downstream_target_df table of regulatory proteins (columns) with gene 
#' targets we're testing for each (entries)
#' @param patient_df table of patient-specific information (sex, age, total number 
#' of mutations, treated or not)
#' @param mutation_df_targ table of mutation counts for each patient, for each 
#' gene in the genome
#' @param methylation_df table of methylation results (rows are proteins, columns 
#' are patients, entries are methylation values)
#' @param methylation_df_meQTL OPTIONAL: table of methylation results (rows are 
#' proteins, columns are patients, entries are raw methylation values to be used
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
#' @param outpath a string with the local path for debugging/ collinearity files 
#' we'll be writing from inside the function
#' @param removeCis a TRUE/FALSE value indicating whether or not we are eliminating 
#' all cis pairings
#' @param removeMetastatic a TRUE/FALSE value indicating whether or not metastatic
#' samples were removed, for file naming
#' @param all_genes_id_conv an ID conversion file from BioMart, for use in removing 
#' cis pairings if needed
#' @param incl_nextMutDriver a TRUE/FALSE value indicating, if we are running just on
#' TP53 or PIK3CA in BRCA, whether or not we include mutation/ CNA/ mutation interaction
#' @param run_query_genes_jointly a TRUE/ FALSE value indicating whether or not we 
#' want to run all query regulatory proteins jointly in the same model, or separately
#' (one model per query-target pair)
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3', to denote what columns
#' we are referencing/ specific data types we are using
#' @param filename_labels a list of labels used for appropriate filename generation
create_lm_input_table <- function(protein_ids_df, downstream_target_df, patient_df, 
                                  mutation_df_targ, methylation_df, methylation_df_meQTL, 
                                  cna_df, expression_df, num_funct_copies_df,
                                  neighboring_cna_df, is_rank_or_quant_norm, log_expression, analysis_type, 
                                  tumNormMatched, randomize, cna_bucketing, meth_bucketing, 
                                  useNumFunctCopies, num_PEER, num_pcs, debug, outpath, 
                                  fixed_lambda_for_rand, removeCis, removeMetastatic, 
                                  all_genes_id_conv, incl_nextMutDriver, run_query_genes_jointly, 
                                  dataset, filename_labels) {
  
  # If we are randomizing expression, do this now per-sample
  if(randomize) {
    if(analysis_type == "eQTL") {
      random_exp_df <- apply(expression_df[,2:ncol(expression_df), with = FALSE], 
                             MARGIN = 2,  function(x) dqrng::dqsample(x))
      expression_df <- cbind(expression_df[,'ensg_id'], data.table(random_exp_df))
    }
    else if (analysis_type == "meQTL") {
      random_methyl_df <- apply(methylation_df[,2:(ncol(methylation_df)-1), with = FALSE], 
                                MARGIN = 2, dqsample)
      methylation_df <- cbind(methylation_df[,'Gene_Symbol'], cbind(data.table(random_methyl_df), 
                                                                    methylation_df[,'ensg_ids']))
    }
    else {print(paste("Analysis type is invalid:", analysis_type))}
  }

  # If we are running individual models per query gene 
  if(!run_query_genes_jointly) {
    protein_input_dfs <- lapply(1:protein_ids_df[, .N], function(i) {
      
      # Get the given regulatory protein r_i's Swissprot & ENSG IDs
      regprot <- protein_ids_df$swissprot_ids[i]
      regprot_ensg <- unlist(strsplit(protein_ids_df$ensg_ids[i], ";", fixed = TRUE))
      
      print(paste("Regulatory protein", paste(i, paste("/", protein_ids_df[, .N]))))
      
      # Create a starter table that gets the mutation, CNA, and methylation status 
      # for this regulatory protein and bind it to the patient DF
      regprot_i_df <- fill_regprot_inputs(patient_df = patient_df,
                                          regprot_i_uniprot = regprot, 
                                          regprot_i_ensg = regprot_ensg, 
                                          mutation_targ_df = mutation_df_targ, 
                                          methylation_df = methylation_df, 
                                          cna_df = cna_df,
                                          num_funct_copies_df = num_funct_copies_DF,
                                          useNumFunctCopies = useNumFunctCopies,
                                          cna_bucketing = cna_bucketing, 
                                          meth_bucketing = meth_bucketing, 
                                          tumNormMatched = tumNormMatched, 
                                          incl_nextMutDriver = incl_nextMutDriver,
                                          run_query_genes_jointly = run_query_genes_jointly,
                                          debug = debug, dataset = dataset)
      starter_df <- cbind(patient_df, regprot_i_df)
      
      if(debug) {
        print("Starter DF, with Regprot Inputs")
        print(head(starter))
      }
      
      # If remove cis is TRUE, limit the downstream target DF to only trans pairings
      if(removeCis) {
        downstream_target_df <- subset_to_trans_targets(regprot_ensg, downstream_target_df,
                                                        all_genes_id_conv, debug)
      }
      
      # Use multithreading to run these jobs in parallel
      num_groups <- 6
      parallel_cluster <- makeCluster(6, type = "SOCK", methods = FALSE, outfile = "")
      
      # Tell R that we want to use these processes to do calculations
      setDefaultCluster(parallel_cluster)
      doParallel::registerDoParallel(parallel_cluster)
      
      list_of_rows <- split(downstream_target_df, f = seq(nrow(downstream_target_df)))
      
      lm_input_tables <- foreach(currentRow = list_of_rows, 
                                 .export = c("methylation_df", "methylation_df_meQTL", "starter_df", "mutation_targ_df",
                                             "cna_df", "expression_df", "log_expression", "cna_bucketing", 
                                             "meth_bucketing", "analysis_type", "tumNormMatched", "debug", 
                                             "dataset", "randomize", "num_PEER", "num_pcs", "removeCis", "removeMetastatic",
                                             "model_type", "run_query_genes_jointly", "regprot", "filename_labels"), 
                                 #.combine = function(x,y) list(x, y),  # don't need a combine function here, since it defaults to returning a list
                                 .inorder = FALSE) %dopar% {
                                   
                                   
         # Source the necessary files for the worker node
         source_path <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/Sara_LinearModel/"
         source(paste(source_path, "general_important_functions2.R", sep = "/"))
         source(paste(source_path, "linear_model_helper_functions2.R", sep = "/"))
         
         # Get the target t_k's Swissprot & ENSG IDs
         targ <- unlist(strsplit(currentRow$swissprot, ";", fixed = TRUE))
         targ_ensg <- unlist(strsplit(currentRow$ensg, ";", fixed = TRUE))
         
         if(debug) {print(targ)}
         
         tryCatch({
           lm_input_table <- fill_targ_inputs(starter_df = starter_df, targ_k = targ, 
                                              targ_k_ensg = targ_ensg, 
                                              mutation_targ_df = mutation_targ_df,
                                              methylation_df = methylation_df, 
                                              cna_df = cna_df, expression_df = expression_df, 
                                              log_expression = log_expression,
                                              cna_bucketing = cna_bucketing,
                                              meth_bucketing = meth_bucketing, 
                                              analysis_type = analysis_type,
                                              tumNormMatched = tumNormMatched, 
                                              debug = debug, dataset = dataset)
         }, error = function(cond) {
           print(cond)
           lm_input_table <- NA
         })
         
         if((length(lm_input_table) == 0) | (lm_input_table[, .N] == 0)) {
           print("Something went wrong; LM input table is empty. Returning NA.")
           lm_input_table <- NA
         }
         
         if(!is.na(lm_input_table)) {
           # First, remove any columns that are entirely NA (e.g., a given driver did not 
           # have CNA or methylation data reported)
           lm_input_table <- lm_input_table[, which(unlist(lapply(lm_input_table, function(x)
             !all(is.na(x))))), with = FALSE]
           
           # Then, if a row has NA, remove that whole row and predict without it
           #lm_input_table <- na.exclude(lm_input_table)
           lm_input_table <- na.omit(lm_input_table)
           
           # Make sure everything is numeric type
           #lm_input_table[,2:ncol(lm_input_table)] <- sapply(lm_input_table[,2:ncol(lm_input_table)], as.numeric) 
           
           if(debug) {
             print("LM Input Table")
             print(head(lm_input_table))
             print(dim(lm_input_table))
           }
           
           # Write this file to the appropriate outpath with the generated file name
           tryCatch({
             lm_input_table_fn <- create_lm_input_table_filename(TRUE, regprot, filename_labels[["run_name"]], 
                                                                 targ, filename_labels[["expression_df_name"]],
                                                                 cna_bucketing, meth_bucketing, 
                                                                 filename_labels[["meth_type"]], 
                                                                 filename_labels[["patient_df_name"]], num_PEER, 
                                                                 num_pcs, randomize, filename_labels[["patients_to_incl_label"]], 
                                                                 removeCis, removeMetastatic)
             
             data.table::fwrite(lm_input_table, paste(outpath, paste0(lm_input_table_fn, ".csv"), sep = "/"))
           }, error = function(cond) {
             print(cond)
             print(dim(lm_input_table))
             print(paste(outpath, paste0(lm_input_table_fn, ".csv"), sep = "/"))
           })
         }
         # Ensure this is what is returned by this iteration of foreach
         lm_input_table <- lm_input_table
      }
      return(lm_input_tables)
    })
    return(protein_input_dfs)
  }
  
  # Otherwise, if we are running all our regulatory proteins jointly in one model
  else {
    
    # Create a starter DF from all the input regulatory proteins
    protein_input_dfs <- lapply(1:protein_ids_df[, .N], function(i) {
      
      # Get the given regulatory protein r_i's Swissprot & ENSG IDs
      regprot <- protein_ids_df$swissprot_ids[i]
      regprot_ensg <- unlist(strsplit(protein_ids_df$ensg_ids[i], ";", fixed = TRUE))
      
      print(paste("Regulatory protein", paste(i, paste("/", protein_ids_df[, .N]))))
      
      # Create a starter table that gets the mutation, CNA, and methylation status 
      # for this regulatory protein and bind it to the patient DF
      regprot_i_df <- fill_regprot_inputs(patient_df = patient_df,
                                          regprot_i_uniprot = regprot, 
                                          regprot_i_ensg = regprot_ensg, 
                                          mutation_targ_df = mutation_df_targ, 
                                          methylation_df = methylation_df, 
                                          cna_df = cna_df,
                                          num_funct_copies_df = num_funct_copies_DF,
                                          useNumFunctCopies = useNumFunctCopies,
                                          cna_bucketing = cna_bucketing, 
                                          meth_bucketing = meth_bucketing, 
                                          tumNormMatched = tumNormMatched, 
                                          incl_nextMutDriver = incl_nextMutDriver,
                                          run_query_genes_jointly = run_query_genes_jointly,
                                          debug = debug, dataset = dataset)
      setnames(regprot_i_df, "MutStat_i", paste0(regprot, "_MutStat_i"))
      setnames(regprot_i_df, "CNAStat_i", paste0(regprot, "_CNAStat_i"))
      setnames(regprot_i_df, "MethStat_i", paste0(regprot, "_MethStat_i"))
      
      # If remove cis is TRUE, limit the downstream target DF to only trans pairings
      #if(removeCis) {
      #  downstream_target_df <- subset_to_trans_targets(regprot_ensg, downstream_target_df,
      #                                                  all_genes_id_conv, debug)
      #}
      
      return(regprot_i_df)
    })
    
    # Bind these together
    regprot_df <- do.call(cbind, protein_input_dfs)
    print(head(regprot_df))
    
    # Add the patient/ sample DF to this 
    starter_df <- cbind(patient_df, regprot_df)
    
    if(debug) {
      print("Starter DF, with Regprot Inputs")
      print(head(starter_df))
    }
    
    # Use multithreading to run these jobs in parallel
    num_groups <- 6
    parallel_cluster <- makeCluster(6, type = "SOCK", methods = FALSE, outfile = "")
    
    # Tell R that we want to use these processes to do calculations
    setDefaultCluster(parallel_cluster)
    doParallel::registerDoParallel(parallel_cluster)
    
    list_of_rows <- split(downstream_target_df, f = seq(nrow(downstream_target_df)))
    if(dataset == "Chinese_TN") {
      methylation_df <- NA
      methylation_df_meQTL <- NA
    }
    
    lm_input_tables <- foreach(currentRow = list_of_rows, 
                              .export = c("methylation_df", "methylation_df_meQTL", "starter_df", "mutation_targ_df",
                                          "cna_df", "expression_df", "log_expression", "cna_bucketing", 
                                          "meth_bucketing", "analysis_type", "tumNormMatched", "debug", 
                                          "dataset", "randomize", "num_PEER", "num_pcs", "removeCis",  
                                          "removeMetastatic", "run_query_genes_jointly", "filename_labels"), 
                              #.combine = function(x,y) list(x, y),  # don't need a combine function here, since it defaults to returning a list
                              .inorder = FALSE) %dopar% {
                                
                                
        # Source the necessary files for the worker node
        source_path <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/Sara_LinearModel/"
        source(paste(source_path, "general_important_functions2.R", sep = "/"))
        source(paste(source_path, "linear_model_helper_functions2.R", sep = "/"))
        
        # Get the target t_k's Swissprot & ENSG IDs
        targ <- unlist(strsplit(currentRow$swissprot, ";", fixed = TRUE))
        targ_ensg <- unlist(strsplit(currentRow$ensg, ";", fixed = TRUE))
        
        if(debug) {print(targ)}
        
        lm_input_table <- fill_targ_inputs(starter_df = starter_df, targ_k = targ, 
                                           targ_k_ensg = targ_ensg, 
                                           mutation_targ_df = mutation_targ_df,
                                           methylation_df = methylation_df, 
                                           cna_df = cna_df, expression_df = expression_df, 
                                           log_expression = log_expression,
                                           cna_bucketing = cna_bucketing,
                                           meth_bucketing = meth_bucketing, 
                                           analysis_type = analysis_type,
                                           tumNormMatched = tumNormMatched, 
                                           debug = debug, dataset = dataset)
        
        if((length(lm_input_table) == 0) | (lm_input_table[, .N] == 0)) {
          print("Something went wrong; LM input table is empty. Returning NA.")
          lm_input_table <- NA
        }
        
        if(!is.na(lm_input_table)) {

          tryCatch({
            # First, remove any columns that are entirely NA or 0 (e.g., a given driver did not 
            # have CNA or methylation data reported)
            lm_input_table <- lm_input_table[, which(unlist(lapply(lm_input_table, function(x)
              !all(is.na(x))))), with = FALSE]
            lm_input_table <- lm_input_table[, colSums(lm_input_table != 0, na.rm = T) > 0, with = FALSE]
            
            # Then, if a row has NA, remove that whole row and predict without it
            #lm_input_table <- na.exclude(lm_input_table)
            lm_input_table <- na.omit(lm_input_table)
            
            # Make sure everything is numeric type
            #lm_input_table[,2:ncol(lm_input_table)] <- sapply(lm_input_table[,2:ncol(lm_input_table)], as.numeric) 
            
            if(debug) {
              print("LM Input Table")
              print(head(lm_input_table))
              print(dim(lm_input_table))
            }
            
            print("LM Input Table")
            print(head(lm_input_table))
            
            # Create an output file name for this LM input table
            lm_input_table_fn <- create_lm_input_table_filename(filename_labels[["test"]], filename_labels[["tester"]], 
                                                                filename_labels[["run_name"]], targ, filename_labels[["expression_df_name"]],
                                                                cna_bucketing, meth_bucketing, filename_labels[["meth_type"]], 
                                                                filename_labels[["patient_df_name"]], num_PEER, 
                                                                num_pcs, randomize, filename_labels[["patients_to_incl_label"]], 
                                                                removeCis, removeMetastatic)
            
            # Write this file to the appropriate outpath with the generated file name
            data.table::fwrite(lm_input_table, file = paste(outpath, paste0(lm_input_table_fn, ".csv"), sep = "/"))
          
          }, error = function(cond) {
            print(cond)
            print(length(lm_input_table_fn))
            print(paste(outpath, paste0(lm_input_table_fn, ".csv"), sep = "/"))
          })
        }
        
        lm_input_table <- lm_input_table
    }
    # tell R that we don't need the processes anymore
    stopCluster(parallel_cluster)
    
    return(lm_input_tables)
  }
}


############################################################
#### FILL IN REGULATORY PROTEIN R_I INPUTS TO TABLE
############################################################
#' Function takes the data files for regulatory protein r_i (and its IDs), 
#' and uses them to  construct a partial tabular input to a linear model 
#' function for a given regulatory protein r_i of the following form:
#' Columns: Linear model input variables
#' Rows: Patients 
#' @param patient_df the starter DF with all of the patient/ sample information
#' @param regprot_i_uniprot the uniprot ID of regulatory protein i
#' @param regprot_i_ensg the ensembl ID of regulatory protein i
#' @param mutation_targ_df the mutation DF for all proteins 
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
#' @param run_query_genes_jointly a TRUE/ FALSE value indicating whether or not we 
#' want to run all query regulatory proteins jointly in the same model, or separately
#' (one model per query-target pair)
#' @param debug a TRUE/FALSE value indicating if we are in debug mode and should 
#' use additional prints
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3', to denote what columns
#' we are referencing/ specific data types we are using
fill_regprot_inputs <- function(patient_df, regprot_i_uniprot, regprot_i_ensg, mutation_targ_df, 
                                methylation_df, cna_df, num_funct_copies_df,
                                useNumFunctCopies, cna_bucketing, meth_bucketing, tumNormMatched, 
                                incl_nextMutDriver, run_query_genes_jointly, debug, dataset) {
  
  # Filter the data frames to look at only this regulatory protein
  mutation_targ_df_sub <- mutation_targ_df[grepl(regprot_i_uniprot, mutation_targ_df$Swissprot),]
  print(dim(mutation_targ_df_sub))
  cna_df_sub <- filter_cna_by_ensg(cna_df, regprot_i_ensg)  
  if(dataset != "Chinese_TN") {
    methylation_df_sub <- filter_meth_by_ensg(methylation_df, regprot_i_ensg) 
  }
   
  #print(paste("REGPROT_i ENSG:", regprot_i_ensg))
  
  # Optionally, add mutation, CNA status covariates for PIK3CA if regprot is TP53 and vice versa 
  if(incl_nextMutDriver & !(run_query_genes_jointly)) {
    tp53_pik3ca_dfs <- get_next_mut_driver_dfs(regprot_i_ensg, cna_df, targ_df)  # TODO: fix this function if we end up needing it again
    cna_df_pik3ca <- tp53_pik3ca_dfs[[1]]
    mutation_regprot_df_pik3ca <- tp53_pik3ca_dfs[[2]]
    cna_df_tp53 <- tp53_pik3ca_dfs[[3]]
    mutation_regprot_df_tp53 <- tp53_pik3ca_dfs[[4]]
  }
  
  # If tumor-normal matched, limit to just cancer samples (we'll ID the matched normal in each run)
  if(tumNormMatched) {patient_df2 <- patient_df[grepl("-0", patient_df$sample_id)]} 
  else {patient_df2 <- patient_df}
  
  # Loop through all tumor samples to get the info on this regulatory protein for each
  regprot_rows <- mclapply(1:patient_df2[, .N], function(i) {
    sample <- patient_df2$sample_id[i]
    if (debug) {print(paste("Sample:", sample))}
    #print(paste("Sample:", sample))
    # Is this regulatory protein mutated in this sample at the given level of specificity?
    # OPT 1: FOR I-DOMAIN & I-BINDING POSITION:
    #if (nrow(mutation_regprot_df %>% filter(Patient == sample)) > 0) {mut_stat <- 1}
    # OPT 2: FOR I-PROTEIN:
    #mutation_regprot_df_sub <- mutation_regprot_df_sub[Patient %like% sample]
    #if(mutation_regprot_df_sub[, .N] > 0) {mut_stat <- 1}
    #else {mut_stat <- 0}
    # OPT 3: FOR MUT COUNT MATRIX
    #print(which(colnames(mutation_targ_df_sub) == sample))
    val <- as.integer(unlist(mutation_targ_df_sub[,which(colnames(mutation_targ_df_sub) == sample), with = F]))
    if(length(val) == 0) {mut_stat <- NA}
    else if(val == 0) {mut_stat <- 0}
    else {mut_stat <- 1}
    
    # Does this regulatory protein have a CNA in cancer?
    cna_stat <- get_cna_stat(cna_df_sub, sample, cna_bucketing, dataset, FALSE, NA)
    #cna_vect_sample <- unlist(cna_df[,colnames(cna_df) == sample, with = F])
    #cna_stat <- get_cna_stat(cna_df_sub, sample, cna_bucketing, dataset, TRUE, cna_vect_sample)
    if(grepl("incl", cna_bucketing)) {cna_stat <- as.integer(cna_stat[[1]])}     
    
    # Does this regulatory protein have a methylation marker in cancer, or differential 
    # methylation between tumor and normal?
    if(dataset != "Chinese_TN") {
      meth_stat <- get_meth_stat(methylation_df_sub, sample, meth_bucketing, tumNormMatched)
      if(meth_bucketing) {
        meth_stat <- as.integer(meth_stat[[1]])
        if(length(meth_stat) > 3) {meth_stat <- meth_stat[1:3]}
      }
    } else {
      meth_stat <- NA
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
    if (incl_nextMutDriver & !(run_query_genes_jointly)) {
      if(regprot_i_ensg == "ENSG00000141510") {
        cna_b <- cna_bucketing
        if(grepl("just", cna_bucketing)) {cna_b <- "bucket_justAmp"}
        altern_stats <- get_altern_stats(regprot_i_ensg, cna_df_sub, cna_df_pik3ca,
                                         cna_b, mutation_regprot_df_pik3ca,
                                         sample, dataset, mut_stat)
      }
      else if(regprot_i_ensg  == "ENSG00000121879") {
        cna_b <- cna_bucketing
        if(grepl("just", cna_bucketing)) {cna_b <- "bucket_justDel"}
        altern_stats <- get_altern_stats(regprot_i_ensg, cna_df_sub, cna_df_tp53,
                                         cna_b, mutation_regprot_df_tp53,
                                         sample, dataset, mut_stat)
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
        altern_stats <- c(altern_stats_tp53, altern_stats_pik3ca)
      }
    }
    
    if (debug) {
      print(paste("Meth stat:", meth_stat))
      print(paste("CNA stat:", cna_stat))
      print(paste("Mut stat:", mut_stat))
      if(useNumFunctCopies) {print(paste("FNC stat:", fnc_stat))}
      #print("Altern driver stats:")
      #print(altern_stats)
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
      if(incl_nextMutDriver & !(run_query_genes_jointly)) {
        outdf <- add_altern_stats_to_outdf(outdf, altern_stats, regprot_i_ensg)
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
          outdf <- add_altern_stats_to_outdf(outdf, altern_stats, regprot_i_ensg)
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
      if(incl_nextMutDriver & !(run_query_genes_jointly)) {
        outdf <- add_altern_stats_to_outdf(outdf, altern_stats, regprot_i_ensg)
      }
    }
    return(outdf)
  })
  
  # Bind these rows together
  if(debug) {print(head(regprot_rows))}
  tryCatch({
    colnam <- colnames(regprot_rows[[1]])
    regprot_i_df <- data.table::rbindlist(regprot_rows)
    colnames(regprot_i_df) <- colnam
    
  }, error=function(cond){
    print("Unexpected error; regprot rows is not valid: ")
    print(regprot_rows)
    return(regprot_i_df)
  })
  
  # Return this bound DF
  return(regprot_i_df)
}



#' Helper function if we are including TP53- and PIK3CA-specific covariates in
#' breast cancer to obtain TP53 and PIK3CA CNA and mutation regulatory protein DFs
#' CURRENTLY DEFUNCT AS OF FEB 2023
#' @param regprot_i_ensg the ensg ID for the current regulatory protein of interest
#' @param cna_df the full CNA input DF, to be subsetted
#' @param mutation_regprot_df the full regulatory protein mutation DF, to be subsetted
get_next_mut_driver_dfs <- function(regprot_i_ensg, cna_df, mutation_regprot_df) {
  cna_df_pik3ca <- NA
  cna_df_tp53 <- NA
  mutation_regprot_df_pik3ca <- NA
  mutation_regprot_df_tp53 <- NA
  
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
  return(list(cna_df_pik3ca, mutation_regprot_df_pik3ca, cna_df_tp53, mutation_regprot_df_tp53))
}


#' Helper function that will get stats for TP53 or PIK3CA if we are adding them
#' as covariates as well 
#' CURRENTLY DEFUNCT AS OF FEB 2023
#' @param regprot_i_ensg the ensg ID for the current regulatory protein of interest
#' @param cna_df_sub the CNA DF, subsetted to the regulatory protein of interest
#' @param cna_df_altern the CNA DF, subsetted to the next mutated driver gene (either 
#' PIK3CA or TP53)
#' @param bucket_type the bucketing type for CNA, for use with the get_cna_stat function
#' @param mutation_regprot_df the mutation regprot DF, unsubsetted
#' @param sample the current sample ID
#' @param dataset the dataset name
#' @param mut_stat the mutation status of the regulatory protein of interest
get_altern_stats <- function(regprot_i_ensg, cna_df_sub, cna_df_altern, bucket_type,
                             mutation_regprot_df, sample, dataset, mut_stat) {
  
  altern_cna_stat <- get_cna_stat(cna_df_altern, sample, bucket_type, dataset, FALSE, NA)
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

#' Helper function to add the alternative stats (e.g. CNA status or mutation status
#' of TP53 and/or PIK3CA) to the output data frame
#' @param outdf the regulatory output df we are currently building
#' @param altern_stats the alternative CNA/ mutation statuses for TP53 and/or PIK3CA
#' @param regprot_i_ensg the ensg ID of the current regulatory protein of interest
add_altern_stats_to_outdf <- function(outdf, altern_stats, regprot_i_ensg) {
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
  return(outdf)
}

############################################################
#### CREATE OUTPATH 
############################################################
# Create an appropriate file outpath
if(args$cancerType == "PanCancer") {
  if(!(args$specificTypes == "ALL")) {
    cancer_types <- unlist(strsplit(args$specificTypes, ";", fixed = TRUE))
    outpath <- lapply(cancer_types, function(ct) {
      return(create_lm_input_file_outpath(cancerType = args$cancerType, specificType = ct, test = test, 
                                 tester_name = args$tester_name, run_name = args$run_name,
                                 tumNormMatched = tumNormMatched, QTLtype = args$QTLtype,
                                 dataset = args$dataset))
    })
    names(outpath) <- cancer_types
  } else {
    outpath <- create_lm_input_file_outpath(cancerType = args$cancerType, specificType = "", test = test, 
                                   tester_name = args$tester_name, run_name = args$run_name,
                                   tumNormMatched = tumNormMatched, QTLtype = args$QTLtype,
                                   dataset = args$dataset)
  }
} else {
  outpath <- create_lm_input_file_outpath(cancerType = args$cancerType, specificType = "", test = test, 
                                 tester_name = args$tester_name, run_name = args$run_name,
                                 tumNormMatched = tumNormMatched, QTLtype = args$QTLtype,
                                 dataset = args$dataset)
}


if(debug) {
  print(paste("Outpath:", outpath))
}

# Note: Output file name will be created dynamically per-gene

## OPTIONAL: If we've already started this run, limit only to genes that are not already
## in the given outpath
# Get the index of where the target gene name will be found in each file name
gene_name_index <- length(unlist(strsplit(args$run_name, "_", fixed = TRUE))) + 3
if(args$patientsOfInterest != "") {
  gene_name_index <- gene_name_index + 1
  outpath <- paste0(outpath, paste0("/subtype_files", args$patientOfInterestLabel))
}
#targets_DF <- limit_to_targets_wo_existing_files(outpath, targets_DF, gene_name_index)

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

# Get the labels needed for dynamic filename creation and add to a list
filename_labels <- list("expression_df_name" = args$expression_df,
                        #"mutation_targ_df_name" = args$mutation_targ_df,
                        "meth_type" = args$meth_type,
                        "patient_df_name" = args$patient_df,
                        "test" = test, "tester_name" = args$tester_name, 
                        "run_name" = args$run_name,
                        "patients_to_incl_label" = args$patientsOfInterestLabel)


# Case 1: If we're running on just one cancer type or PanCancer
if(is.na(patient_cancer_mapping)) {
  
  # If needed, subset these data tables by the given patient IDs of interest using helper functions
  if(patients_of_interest != "") {
    
    patient_df <- subset_by_intersecting_ids(patients_of_interest, patient_df, FALSE, tumNormMatched, args$dataset)
    mutation_targ_df <- subset_by_intersecting_ids(patients_of_interest, mutation_targ_df, TRUE, tumNormMatched, args$dataset)
    #mutation_regprot_df <- subset_regprot_df_by_intersecting_ids(patients_of_interest, mutation_regprot_df, tumNormMatched, args$dataset)
    if(args$dataset != "Chinese_TN") {
      methylation_df <- subset_by_intersecting_ids(patients_of_interest, methylation_df, TRUE, tumNormMatched, args$dataset)
      if(!is.na(methylation_df_meQTL)) {
        methylation_df_meQTL <- subset_by_intersecting_ids(patients_of_interest, methylation_df_meQTL, TRUE, tumNormMatched, args$dataset)
      } else {methylation_df_meQTL <- NA} 
    } else {
      methylation_df <- NA
      methylation_df_meQTL <- NA
    }
    cna_df <- subset_by_intersecting_ids(patients_of_interest, cna_df, TRUE, tumNormMatched, args$dataset)
    expression_df <- subset_by_intersecting_ids(patients_of_interest, expression_df, TRUE, tumNormMatched, args$dataset)
    if(useNumFunctCopies) {num_funct_copies_DF <- subset_by_intersecting_ids(patients_of_interest, num_funct_copies_DF, FALSE, tumNormMatched, args$dataset)}
  }
  
  # If needed, exclude metastatic samples
  if(removeMetastatic == TRUE) {
    patient_df <- remove_metastatic_samples(patient_df, FALSE)
    mutation_targ_df <- remove_metastatic_samples(mutation_targ_df, TRUE)
    #mutation_regprot_df <- remove_metastatic_samples_regprot(mutation_regprot_df)
    if(args$dataset != "Chinese_TN") {
      methylation_df <- remove_metastatic_samples(methylation_df, TRUE)
      if(!is.na(methylation_df_meQTL)) {
        methylation_df_meQTL <- remove_metastatic_samples(methylation_df_meQTL, TRUE)
      } else {methylation_df_meQTL <- NA} 
    } 
    cna_df <- remove_metastatic_samples(cna_df, TRUE)
    expression_df <- remove_metastatic_samples(expression_df, TRUE)
    if(useNumFunctCopies) {num_funct_copies_DF <- remove_metastatic_samples(num_funct_copies_DF, FALSE)}
  }
  
  
  lm_input_tables <- create_lm_input_table(protein_ids_df = protein_ids_df, 
                                          downstream_target_df = targets_DF, 
                                          patient_df = patient_df,
                                          mutation_df_targ =  mutation_targ_df,
                                          #mutation_df_regprot = mutation_regprot_df, 
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
                                          outpath = outpath, 
                                          removeCis = removeCis,
                                          removeMetastatic = removeMetastatic,
                                          all_genes_id_conv = all_genes_id_conv,
                                          useNumFunctCopies = useNumFunctCopies,
                                          num_funct_copies_df = num_funct_copies_DF,
                                          incl_nextMutDriver = incl_nextMutDriver,
                                          run_query_genes_jointly = runRegprotsJointly,
                                          dataset = args$dataset,
                                          filename_labels = filename_labels)
  
# Case 2: If we're running on multiple cancer types separately
} else if (!(is.na(patient_cancer_mapping))) {
  # Handle the case of running multiple cancer types per-cancer
  # Subset all the input files for each cancer type and run separately
  input_df_list <- lapply(1:length(outpath), function(i) {
    
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
    #mutation_regprot_df_sub <- subset_regprot_df_by_intersecting_ids(patient_ids, mutation_regprot_df, tumNormMatched, args$dataset)
    if(args$dataset != "Chinese_TN") {
      methylation_df_sub <- subset_by_intersecting_ids(patient_ids, methylation_df, TRUE, tumNormMatched, args$dataset)
      if(!is.na(methylation_df_meQTL)) {
        methylation_df_meQTL_sub <- subset_by_intersecting_ids(patient_ids, methylation_df_meQTL, TRUE, tumNormMatched, args$dataset)
      } else {methylation_df_meQTL_sub <- NA} 
    } else {
      methylation_df <- NA
      methylation_df_meQTL <- NA
    }
    cna_df_sub <- subset_by_intersecting_ids(patient_ids, cna_df, TRUE, tumNormMatched, args$dataset)
    expression_df_sub <- subset_by_intersecting_ids(patient_ids, expression_df, TRUE, tumNormMatched, args$dataset)
    if(useNumFunctCopies) {num_funct_copies_DF_sub <- subset_by_intersecting_ids(patients_of_interest, num_funct_copies_DF, FALSE, tumNormMatched, args$dataset)}
    
    # If needed, exclude metastatic samples
    if(removeMetastatic == TRUE) {
      patient_df <- remove_metastatic_samples(patient_df, FALSE)
      mutation_targ_df <- remove_metastatic_samples(mutation_targ_df, TRUE)
      #mutation_regprot_df <- remove_metastatic_samples_regprot(mutation_regprot_df)
      if(args$dataset != "Chinese_TN") {
        methylation_df <- remove_metastatic_samples(methylation_df, TRUE)
        if(!is.na(methylation_df_meQTL)) {
          methylation_df_meQTL <- remove_metastatic_samples(methylation_df_meQTL, TRUE)
        } else {methylation_df_meQTL <- NA} 
      } 
      cna_df <- remove_metastatic_samples(cna_df, TRUE)
      expression_df <- remove_metastatic_samples(expression_df, TRUE)
      if(useNumFunctCopies) {num_funct_copies_DF <- remove_metastatic_samples(num_funct_copies_DF, FALSE)}
    }
    
    name <- names(patient_cancer_mapping)[i]
    protein_ids_df_spec <- protein_ids_df[[name]]
    
    # Run the model with these subsetted files
    lm_input_tables <- create_lm_input_table(protein_ids_df = protein_ids_df_spec, 
                                            downstream_target_df = targets_DF, 
                                            patient_df = patient_df_sub,
                                            mutation_df_targ =  mutation_targ_df_sub,
                                            #mutation_df_regprot = mutation_regprot_df_sub, 
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
                                            outpath = outpath_i, 
                                            removeCis = removeCis,
                                            removeMetastatic = removeMetastatic,
                                            all_genes_id_conv = all_genes_id_conv,
                                            useNumFunctCopies = useNumFunctCopies,
                                            num_funct_copies_df = num_funct_copies_DF,
                                            incl_nextMutDriver = incl_nextMutDriver,
                                            run_query_genes_jointly = runRegprotsJointly,
                                            dataset = args$dataset,
                                            filename_labels = filename_labels)
    
    return(lm_input_tables)
  })
  
} else {
  print("Something went wrong with patient-cancer type mapping. Exiting now.")
  lm_input_table <- NA
}


