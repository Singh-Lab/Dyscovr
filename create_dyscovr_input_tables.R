############################################################
### Create Dyscovr Input Tables
### PUBLICATION INFORMATION 
############################################################

# This script requires 1 node and 8 CPUs for multithreading.
#!/usr/bin/env Rscript

# Set the library path
library.path <- .libPaths()

library(dplyr, lib.loc = library.path)
library(snow, lib.loc = library.path)
library(argparse, lib.loc = library.path)
library(data.table, lib.loc = library.path)
library(doParallel, lib.loc = library.path)
library(dqrng, lib.loc = library.path)
library(foreach, lib.loc = library.path)
library(parallel, lib.loc = library.path)
library(fastmatch, lib.loc = library.path)
library(stringr, lib.loc = library.path, quietly = T)

# Creates an input table for Dyscovr linear regression in dyscovr.R for each 
# gene target provided in input tables. Writes all input tables to a directory, 
# where they can be fetched when LM runs.

# Source other files needed
SOURCE_PATH <- "/Genomics/argo/users/scamilli/Dyscovr_Personal"
source(paste(SOURCE_PATH, "general_important_functions.R", sep = "/"))
source(paste(SOURCE_PATH, "dyscovr_helper_functions.R", sep = "/"))


############################################################
# SET UP PARSER ARGUMENTS
############################################################
# Create parser object
parser <- ArgumentParser()

# Specify desired options for the parser

# Data set we are running on
parser$add_argument("--dataset", default = "TCGA", type = "character", 
                    help = "TCGA or METABRIC. Default is TCGA.")

# Cancer type we are looking at (or pan-cancer), and patients we are looking at
parser$add_argument("--cancerType", default = "PanCancer", type = "character",
                    help = "Type of cancer we're looking at, or PanCancer. 
                    Default is PanCancer.")
parser$add_argument("--specificTypes", default = "ALL", type = "character",
                    help = "If and only if the cancer type is PanCancer, 
                    denotes which specific cancer types we're including (running 
                    per each specified cancer type individually). If ALL, 
                    run pan-cancer. Default is ALL.")
parser$add_argument("--clinical_df", default = "clinical.csv", type = "character",
                    help = "If amd only if the cancer type is PanCancer and 
                    specific types is not ALL, a clinical data frame name is 
                    given to link cancer types to patient IDs.")
parser$add_argument("--patientsOfInterest", default = "", type = "character",
                    help = "Further restricts the patient pool to particular 
                    patients (e.g. from a particular subtype or with particular
                    properties). A filename with a vector of patient IDs (no 
                    header) is optionally provided here.")
parser$add_argument("--patientsOfInterestLabel", default = "", type = "character",
                    help = "A character label for the specific patient population 
                    we are looking at (for output file naming purposes).")
parser$add_argument("--removeMetastatic", default = "F", type = "character",
                    help = "A T/ F value indicating whether or not we want 
                    to exclude metastatic samples from the analysis.")

# Randomization validation test
parser$add_argument("--randomize", default = "F", type = "character", 
                    help = "T/F for whether or not we perform a
                    randomization test. Default is F.")

# Whether or not we are running a test on one protein-of-interest ("tester" 
# protein), and details of the tester protein if so
parser$add_argument("--test", default = "F", type = "character",
                    help = "T/F for whether or not we are running a test 
                    with one protein-of-interest. Default is F.")
parser$add_argument("--tester_name", default = "P53", type = "character",
                    help = "Name of regulatory protein if doing a 1-protein test. 
                    Defaults to TP53.")
parser$add_argument("--tester_uniprot_id", default = "P04637", type = "character",
                    help = "Uniprot ID of regulatory protein if doing a 1-protein 
                    test. Defaults to TP53 (P04637).")
parser$add_argument("--tester_ensg_id", default = "ENSG00000141510", 
                    type = "character", help = "ENSG ID of regulatory protein if 
                    doing a 1-protein test. Defaults to TP53 (ENSG00000141510).")
parser$add_argument("--run_query_genes_jointly", default = "T", type = "character",
                    help = "If T, will run all the query regulatory proteins 
                    together in the same joint model, rather than running one 
                    model per query-target pair. Default is T.")

# The files to use for the other data inputs
parser$add_argument("--protein_ids_df", default = "protein_ids_df.csv", 
                    type = "character", help = "The name of the protein IDs data 
                    frame input (with ENSG and Uniprot IDs). [default %(default)s]")
parser$add_argument("--expression_df", 
                    default = "expression_quantile_norm_IntersectPatients.csv", 
                    type = "character", help = "The name of the expression data 
                    frame input. [default %(default)s]")
parser$add_argument("--methylation_df", 
                    default = "methylation_M_CancerOnly_IntersectPatients.csv", 
                    type = "character", help = "The name of the methylation data 
                    frame input. [default %(default)s]")
parser$add_argument("--mutation_df", 
                    default = "mut_count_matrix_nonsynonymous_IntersectPatients.csv", 
                    type = "character", help = "The name of the mutation data 
                    frame input. [default %(default)s]")
parser$add_argument("--cna_df", default = "CNA_AllGenes_IntersectPatients.csv", 
                    type = "character", help = "The name of the CNA data frame 
                    input. [default %(default)s]")
parser$add_argument("--patient_df", 
                    default = "combined_patient_sample_cibersort_total_frac_washu_IntersectPatients.csv", 
                    type = "character", help = "The name of the patient sample 
                    data frame input. [default %(default)s]")
parser$add_argument("--target_df", default = "allgene_targets.csv", 
                    type = "character", help = "The name of the gene targets 
                    data frame. [default %(default)s]")

# The decision for if/ how to bucket CNAs and methylation
parser$add_argument("--cna_bucketing", default = "rawCNA", 
                    type = "character", help = "If/ the type of bucketing we use 
                    for CNA values. Default is rawCNA (non-bucketed representation), 
                    but bucket_inclAmp, bucket_exclAmp, bucket_justAmp, and 
                    bucket_justDel are also options.")
parser$add_argument("--meth_bucketing", default = F, type = "character",
                    help = "If/ the type of bucketing we use for methylation 
                    values. Default is F.")

# Label the type of methylation data we are using (Beta, M, or Threshold)
parser$add_argument("--meth_type", default = "M", type = "character",
                    help = "The type methylation values we are using. Default is 
                    M-value, but Beta is another option.")

# What type of test we are running, so we can write the results files to an 
# appropriate directory. Only necessary if --test is F.
parser$add_argument("--run_name", default = "top_drivers", type = "character",
                    help = "Provide a name for the run, if not testing on just 
                    one protein-of-interest, to write output files to appropriate 
                    directory. [default %(default)s]") 

# An option to specify the number of genotype PCs we include in the model
parser$add_argument("--num_pcs", default = 3, type = "integer",
                    help = "A value between 0 and 3 indicating the number of 
                    principal components we are using as covariates in our model. 
                    Default is 3.")

# Add a flag for keeping only trans pairings
parser$add_argument("--removeCis", default = "T", type = "character",
                    help = "A T/ F value to indicate whether or not we 
                    want to remove cis pairings and look only at trans pairings. 
                    Defaults to T, but can be set to F.")

# Add a flag for debugging
parser$add_argument("--debug", default = "F", type = "character",
                    help = "A T/ F value indicating whether or not we want 
                    detailed printing for debugging purposes. Default is F.")


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
randomize <- str2bool(args$randomize)
test <- str2bool(args$test)
meth_bucketing <- str2bool(args$meth_bucketing)
debug <- str2bool(args$debug)
removeCis <- str2bool(args$removeCis)
removeMetastatic <- str2bool(args$removeMetastatic)
runQueriesJointly <- str2bool(args$run_query_genes_jointly)


############################################################
# SET MAIN PATH AND IMPORT GENERALLY NEEDED FILES
############################################################
# Set the paths
INPUT_FILE_PATH <- "/Genomics/argo/users/scamilli/Dyscovr_Personal/input_files/"

if(args$cancerType == "BRCA") {
  TARG_PATH <- paste0(INPUT_FILE_PATH, "BRCA/target_lists/")
  
  if(args$dataset == "METABRIC") {
    MAIN_PATH <- paste0(INPUT_FILE_PATH, "METABRIC/")

  } else {
    MAIN_PATH <- paste0(INPUT_FILE_PATH, "BRCA/")
  }

} else {
  TARG_PATH <- paste0(INPUT_FILE_PATH, "PanCancer/target_lists/")
  MAIN_PATH <- paste0(INPUT_FILE_PATH, "PanCancer/")
}

if(debug) {print(MAIN_PATH)}

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- fread(paste0(INPUT_FILE_PATH, "all_genes_id_conv.csv"), 
                           header = T)


############################################################
# IMPORT PROTEINS-OF-INTEREST (EG DRIVER GENES)
############################################################
protein_ids_df <- data.frame()
if(!test) {
  if(args$cancerType == "PanCancer") {
    if(args$specificTypes == "ALL") {
      protein_ids_df <- fread(paste0(MAIN_PATH, args$protein_ids_df), header = T)
    } else {
      types <- unlist(strsplit(args$specificTypes, ";", fixed = T))
      protein_ids_df <- lapply(types, function(ct) {
        new_name <- paste0(MAIN_PATH, paste0(args$protein_ids_df, paste0(ct, ".csv")))
        return(fread(new_name, header = T))
      })
      names(protein_ids_df) <- types
    }
  } else {
    protein_ids_df <- fread(paste0(MAIN_PATH, args$protein_ids_df), header = T)
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
expression_df <- fread(paste0(MAIN_PATH, paste0("Expression/", args$expression_df)),
                       header = T)

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
is_rank_or_quant_norm <- F
if(grepl("rank", args$expression_df) | grepl("quantile", args$expression_df)) {
  is_rank_or_quant_norm <- T
}

# Determine whether or not we will need to log-transform the data, or whether 
# this has already been done in some way (e.g. by an inverse normal 
# transformation, z-score normalization, or quantile-normalization)
log_expression <- T
if(grepl("transform", args$expression_df) | grepl("zscore", args$expression_df) | 
   grepl("sklearn", args$expression_df)) {log_expression <- F}


############################################################
# IMPORT METHYLATION FILES
############################################################
methylation_df <- fread(paste0(MAIN_PATH, 
                               paste0("Methylation/", args$methylation_df)),
                        header = T)
colnames(methylation_df)[
  which(colnames(methylation_df) == "Gene_Symbol")] <- "gene_name"

if(debug) {
  print("Methylation DF")
  print(head(methylation_df))
}


############################################################
# IMPORT MUTATION FILES
############################################################
mutation_df <- fread(paste0(MAIN_PATH, paste0("Mutation/", 
                                            args$mutation_df)), header = T)
colnames(mutation_df)[which(colnames(mutation_df) == "Gene_Symbol")] <- "gene_name"
if(!("Swissprot" %fin% colnames(mutation_df))) {
  mutation_df <- mutation_df[, Swissprot := unlist(lapply(
    mutation_df$gene_name, function(x) 
      paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == x, 
                                            'uniprot_gn_id'])), collapse = ";")))]
}

if(debug) {
  print("Mutation DF")
  print(head(mutation_df))
}

############################################################
# IMPORT CNA FILES
############################################################
cna_df <- fread(paste0(MAIN_PATH, paste0("CNA/", args$cna_df)), header = T)
colnames(cna_df)[1] <- "ensg_id"

if(debug) {
  print("CNA DF")
  print(head(cna_df))
}


############################################################
# IMPORT PATIENT/SAMPLE FILES
############################################################
# If we are looking per-cancer, read in each individual patient DF
if((args$cancerType == "PanCancer") & (args$specificTypes != "ALL")) {
  patient_path <- paste0(MAIN_PATH, "Sample/")
  per_cancer_fns <- intersect(list.files(patient_path, pattern = args$patient_df),
			      list.files(patient_path, pattern = "cibersort_total_frac"))  #TODO: change to find a less hacky way to exclude the pan-cancer file
  patient_df <- lapply(per_cancer_fns, function(x) {
    df <- fread(paste0(patient_path, x), header = T)
    colnames(df)[1] <- "sample_id"
    return(df)
  })
  names(patient_df) <- unlist(lapply(per_cancer_fns, function(x) 
    unlist(strsplit(x, "_", fixed = T))[7]))
  
} else {
  patient_df <- fread(paste0(MAIN_PATH, paste0("Sample/", args$patient_df)),
                      header = T)
  colnames(patient_df)[1] <- "sample_id"
}


if(debug) {
  print("Patient-Sample DF")
  print(head(patient_df))
}


# Read if the patient IDs of interest, if provided
patients_of_interest <- ""
if(args$patientsOfInterest != "") {
  patients_of_interest <- as.character(unlist(fread(paste0(
    MAIN_PATH, paste0("Sample/groups_of_interest/", args$patientsOfInterest)), 
    header = F)[,1]))
}

if(debug) {
  print("Patients of Interest")
  print(head(patients_of_interest))
}


############################################################
# IMPORT GENE TARGETS DATA FRAME
############################################################
targets_DF <- fread(paste0(TARG_PATH, args$target_df), header = T)

# Remove rows with any blank elements
targets_DF <- targets_DF[!((targets_DF$swissprot == "") | (targets_DF$ensg == "")),]


# Remove the first column, if needed (row numbers)
if(ncol(targets_DF) > 2) {
  targets_DF <- targets_DF[,2:3]
}

# Regardless of method, limit targets to only those that overlap the genes in the 
# expression and methylation DFs
methylation_ensg_ids <- unlist(lapply(methylation_df$ensg_ids, function(ids) 
  unlist(strsplit(as.character(ids), ";", fixed = T))))
rows_to_keep <- unlist(lapply(1:length(targets_DF$ensg), function(i) {
  targs <- unlist(strsplit(as.character(targets_DF[i,'ensg']), ";", fixed = T))
  if (any(targs %fin% expression_df$ensg_id)) {
    if(any(targs %fin% methylation_ensg_ids)) {
      return(T)
    } else {return(F)}
  } else {return(F)}
}))
targets_DF <- targets_DF[rows_to_keep,]


if(debug) {
  print("Targets DF")
  print(head(targets_DF))
}


############################################################
# CREATE THE PROTEIN-OF-INTEREST/ DRIVER-SPECIFIC DYSCOVR
# INPUT DATA FRAME
############################################################
#' @param protein_ids_df DF of the regulatory proteins of interest, found based 
#' on their DNA-binding regions (columns for Swissprot IDs and ENSG IDs)
#' @param patient_df table of patient- and sample-specific information (sex, age, 
#' total number of mutations, etc.)
#' @param mutation_df table of mutation counts for each patient, for each 
#' gene in the genome
#' @param methylation_df table of methylation results (rows are proteins, columns 
#' are patients, entries are methylation values)
#' @param cna_df table of CNA per gene (rows are proteins, columns are patients, 
#' entries are CNA values)
#' @param filename_labels a list of labels used for appropriate file name 
#' generation
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. 
#' Possible values are "bucket_inclAmp", "bucket_exclAmp", "bucket_justAmp", 
#' "bucket_justDel" and "rawCNA"
#' @param meth_bucketing a T/F value indicating whether or not we are bucketing
#' methylation values
#' @param outpath a string with the local path for writing files
#' @param num_pcs a value from 0-3 indicating the number of genotypic principal 
#' components we are including as covariates 
#' @param removeCis a T/F value indicating whether or not we are eliminating 
#' all cis pairings
#' @param removeMetastatic a T/F value indicating whether or not metastatic
#' samples were removed, for file naming
#' @param dataset either 'TCGA' or 'METABRIC', to denote what columns
#' we are referencing/ specific data types we are using
#' @param run_query_genes_jointly a T/ F value indicating whether or not we 
#' want to run all query regulatory proteins jointly in the same model, or 
#' separately (one model per query-target pair)
#' @param debug a T/F value indicating if we are in debug mode and should 
#' use additional prints
create_regulatory_prot_input_df <- function(protein_ids_df, patient_df, 
                                            mutation_df, methylation_df, 
                                            cna_df, filename_labels,
                                            cna_bucketing, meth_bucketing,
                                            outpath, num_pcs, removeCis, 
                                            removeMetastatic, dataset, 
                                            run_query_genes_jointly, debug) {
  
  # If we are running individual models per query gene 
  if(!run_query_genes_jointly) {
    protein_input_dfs <- lapply(1:protein_ids_df[, .N], function(i) {
      
      # Get the given protein-of-interest's Swissprot & ENSG IDs
      prot <- protein_ids_df$swissprot_ids[i]
      prot_ensg <- unlist(strsplit(protein_ids_df$ensg_ids[i], ";", fixed = T))
      
      print(paste("Driver Protein", paste(i, paste("/", protein_ids_df[, .N]))))
      
      # Create a starter table that contains the mutation, CNA, and methylation 
      # status for this protein-of-interest and binds it to the patient DF
      prot_df <- fill_prot_inputs(patient_df = patient_df,
                                          prot_uniprot = prot, 
                                          prot_ensg = prot_ensg, 
                                          mutation_df = mutation_df, 
                                          methylation_df = methylation_df, 
                                          cna_df = cna_df,
                                          cna_bucketing = cna_bucketing, 
                                          meth_bucketing = meth_bucketing, 
                                          run_query_genes_jointly = 
                                            run_query_genes_jointly,
                                          debug = debug, dataset = dataset)
      starter_df <- cbind(patient_df, prot_df)
      
      starter_df <- check_and_write_starter_df(starter_df, prot, filename_labels, 
                                               cna_bucketing, meth_bucketing, 
                                               num_pcs, randomize, removeCis, 
                                               removeMetastatic,
                                               run_query_genes_jointly, 
                                               outpath, debug)
      
      return(starter_df)
    })
    return(protein_input_dfs)
  } 
  # Otherwise, if we are running all our regulatory proteins jointly in one model
  else {
    # Create a starter DF from all the input regulatory proteins
    protein_input_dfs <- lapply(1:protein_ids_df[, .N], function(i) {
      
      # Get the given protein-of-interest's Swissprot & ENSG IDs
      prot <- protein_ids_df$swissprot_ids[i]
      prot_ensg <- unlist(strsplit(protein_ids_df$ensg_ids[i], ";", fixed = T))
      
      print(paste("Driver Protein", paste(i, paste("/", protein_ids_df[, .N]))))
      
      # Create a starter table that gets the mutation, CNA, and methylation status 
      # for this protein-of-interest and binds it to the patient DF
      prot_df <- fill_prot_inputs(patient_df = patient_df,
                                  prot_uniprot = prot, 
                                  prot_ensg = prot_ensg, 
                                  mutation_df = mutation_df, 
                                  methylation_df = methylation_df, 
                                  cna_df = cna_df,
                                  cna_bucketing = cna_bucketing, 
                                  meth_bucketing = meth_bucketing, 
                                  run_query_genes_jointly = 
                                    run_query_genes_jointly,
                                  debug = debug, dataset = dataset)
      setnames(prot_df, "MutStat_i", paste0(prot, "_MutStat_i"))
      setnames(prot_df, "CNAStat_i", paste0(prot, "_CNAStat_i"))
      setnames(prot_df, "MethStat_i", paste0(prot, "_MethStat_i"))
      
      return(prot_df)
    })
    
    # Bind these together
    prot_df <- do.call(cbind, protein_input_dfs)
    if(debug) {print(head(prot_df))}
    
    # Add the patient/ sample DF to this 
    starter_df <- cbind(patient_df, prot_df)
    starter_df <- check_and_write_starter_df(starter_df, prot, filename_labels, 
                                             cna_bucketing, meth_bucketing, 
                                             num_pcs, removeCis, removeMetastatic, 
                                             run_query_genes_jointly, outpath, 
                                             debug)
    return(starter_df)
  }
}


############################################################
#### FILL PROTEIN-OF-INTEREST/ DRIVER-SPECIFIC VALUES FOR 
#### A GIVEN PROTEIN (HELPER)
############################################################
#' Function takes the IDs and data files for protein-of-interest/driver 
#' and uses them to  construct a partial tabular input to a Dyscovr linear model 
#' of the following form: Columns: Linear model variables, Rows: Patients 
#' @param patient_df the starter DF with all of the patient/ sample information
#' @param prot_uniprot the uniprot ID of given protein 
#' @param prot_ensg the ensembl ID of given protein
#' @param mutation_df the mutation DF for all proteins 
#' @param methylation_df the methylation DF
#' @param cna_df the copy number alteration DF
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. 
#' Possible values are "bucket_inclAmp", "bucket_exclAmp", and "rawCNA"
#' @param meth_bucketing a T/F value indicating whether or not we are bucketing
#' methylation values
#' @param run_query_genes_jointly a T/ F value indicating whether or not we 
#' want to run all query proteins jointly in the same model, or separately
#' (one model per query-target pair)
#' @param debug a T/F value indicating if we are in debug mode and should 
#' use additional prints
#' @param dataset either 'TCGA' or 'METABRIC', to denote what columns we are 
#' referencing/ specific data types we are using
fill_prot_inputs <- function(patient_df, prot_uniprot, prot_ensg, mutation_df, 
                             methylation_df, cna_df, cna_bucketing, 
                             meth_bucketing, run_query_genes_jointly, debug, 
                             dataset) {
  
  # Filter the data frames to look at only this protein-of-interest, using
  # helper functions from dyscovr_helper_functions.R
  mutation_df_sub <- mutation_df[grepl(prot_uniprot, mutation_df$Swissprot),]
  cna_df_sub <- filter_cna_by_ensg(cna_df, prot_ensg)  
  methylation_df_sub <- filter_meth_by_ensg(methylation_df, prot_ensg) 
  
  # Loop through all tumor samples to get the info on this protein for each
  if(debug) {
    print("Num patients in patient DF:")
    print(nrow(patient_df))
  }

  prot_rows <- mclapply(1:patient_df[, .N], function(i) {
    sample <- patient_df$sample_id[i]
    if (debug) {print(paste("Sample:", sample))}
    
    # Is this protein mutated in this sample?
    val <- as.integer(unlist(mutation_df_sub[,which(
      colnames(mutation_df_sub) == sample), with = F]))
    if(length(val) == 0) {mut_stat <- NA}
    else if(val == 0) {mut_stat <- 0}
    else {mut_stat <- 1}
    
    # Does this protein have a CNA in cancer?
    cna_vect_sample <- unlist(cna_df[,colnames(cna_df) == sample, with = F])
    cna_stat <- get_cna_stat(cna_df_sub, sample, cna_bucketing, dataset, T, 
                             cna_vect_sample)
    if(grepl("incl", cna_bucketing)) {cna_stat <- as.integer(cna_stat[[1]])}     
    
    # Does this regulatory protein have a methylation marker in cancer?
    meth_stat <- get_meth_stat(methylation_df_sub, sample, meth_bucketing)
    if(meth_bucketing) {
      meth_stat <- as.integer(meth_stat[[1]])
      if(length(meth_stat) > 3) {meth_stat <- meth_stat[1:3]}
    }
    
    if (debug) {
      print(paste("Meth stat:", meth_stat))
      print(paste("CNA stat:", cna_stat))
      print(paste("Mut stat:", mut_stat))
    }
    
    # Fill an output data frame using these values, depending on how we are
    # representing each variable
    outdf <- data.table()

    if(!meth_bucketing) {
      if ((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing))) {
        outdf <- data.table("MutStat_i" = mut_stat, "CNAStat_i" = cna_stat, 
                            "MethStat_i" = meth_stat)
      } else {
        outdf <- data.table("MutStat_i" = mut_stat, "CNAStat_i_b1" = cna_stat[1], 
                            "CNAStat_i_b2" = cna_stat[2], 
                            "CNAStat_i_b3" = cna_stat[3], 
                            "MethStat_i" = meth_stat)
      }
      
    } else {
      if ((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing))) {
        outdf <- data.table("MutStat_i" = mut_stat, "CNAStat_i" = cna_stat, 
                            "MethStat_i_b1" = meth_stat[1],
                            "MethStat_i_b2" = meth_stat[2], 
                            "MethStat_i_b3" = meth_stat[3])
      } else {
        outdf <- data.table("MutStat_i" = mut_stat, "CNAStat_i_b1" = cna_stat[1], 
                            "CNAStat_i_b2" = cna_stat[2],
                            "CNAStat_i_b3" = cna_stat[3], 
                            "MethStat_i_b1" = meth_stat[1], 
                            "MethStat_i_b2" = meth_stat[2], 
                            "MethStat_i_b3" = meth_stat[3])
      }
    }
    return(outdf)
  })
  
  # Bind these rows together
  tryCatch({
    colnam <- colnames(prot_rows[[1]])
    prot_df <- data.table::rbindlist(prot_rows)
    colnames(prot_df) <- colnam
    
  }, error=function(cond){
    print("Unexpected error; prot rows is not valid: ")
    print(prot_rows)
  })
  
  # Return this bound DF
  return(prot_df)
}


#' Helper function to check and write the starter DF to file
#' @param starter_df a "starter" DF with patient/ sample and protein-of-interest 
#' data
#' @param prot the protein(s) of interest, e.g. driver genes
#' @param filename_labels a list of labels used for appropriate file 
#' name generation
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. 
#' Possible values are "bucket_inclAmp", "bucket_exclAmp", "bucket_justAmp", 
#' "bucket_justDel" and "rawCNA"
#' @param meth_bucketing a T/F value indicating whether or not we are bucketing
#' methylation values
#' @param num_pcs a value from 0-3 indicating the number of  
#' principal components we are including as covariates in the model 
#' @param removeCis a T/F value indicating whether or not we are eliminating 
#' all cis pairings
#' @param removeMetastatic a T/F value indicating whether or not metastatic
#' samples were removed, for file naming
#' @param run_query_genes_jointly a T/ F value indicating whether or not we 
#' want to run all query driver proteins jointly in the same model, or 
#' separately (one model per query-target pair)
#' @param outpath a string with the local path for writing files
#' @param debug a T/F value indicating if we are in debug mode and should use 
#' additional prints
check_and_write_starter_df <- function(starter_df, prot, filename_labels,
                                       cna_bucketing, meth_bucketing, 
                                       num_pcs, removeCis, removeMetastatic,
                                       run_query_genes_jointly, outpath, debug) {
  # Print the data frame, if in debug mode
  if(debug) {
    print("Starter DF, with Driver Protein Inputs")
    print(head(starter_df))
  }
  
  if((length(starter_df) == 0) | (starter_df[, .N] == 0)) {
    print("Something went wrong; LM input table is empty. Returning NA.")
    starter_df <- NA
  }
  
  if(!is.na(starter_df)) {
    # First, remove any columns that are entirely NA (e.g., a given gene did not 
    # have CNA or methylation data reported)
    starter_df <- starter_df[, which(unlist(lapply(starter_df, function(x)
      !all(is.na(x))))), with = F]
    
    # Then, if a row has NA, remove that whole row and predict without it
    starter_df <- na.omit(starter_df)
    
    if(debug) {
      print("Starter DF, Driver Protein Inputs, NA removed")
      print(head(starter_df))
      print(dim(starter_df))
    }
    
    # Write this file to the appropriate outpath with the generated file name,
    # using helper function in dyscovr_helper_functions.R
    tryCatch({
      starter_df_fn <- NA
      
      if(!run_query_genes_jointly) {
        starter_df_fn <- create_driver_input_table_filename(
          T, prot, filename_labels[["run_name"]], cna_bucketing, meth_bucketing, 
          filename_labels[["meth_type"]], filename_labels[["patient_df_name"]],  
          num_pcs, filename_labels[["patients_to_incl_label"]], removeCis, 
          removeMetastatic)
      } else {
        starter_df_fn <- create_driver_input_table_filename(
          F, NA, filename_labels[["run_name"]], cna_bucketing, meth_bucketing,
          filename_labels[["meth_type"]], filename_labels[["patient_df_name"]], 
          num_pcs, filename_labels[["patients_to_incl_label"]], removeCis, 
          removeMetastatic)
      }
      if(debug) {print(starter_df_fn)}
      fwrite(starter_df, paste(outpath, paste0(starter_df_fn, ".csv"), sep = "/"))
      
    }, error = function(cond) {
      print(cond)
      print(dim(starter_df))
      print(paste(outpath, paste0(starter_df_fn, ".csv"), sep = "/"))
    })
  }
  return(starter_df)
}

#
#
#
#
#

############################################################
# MAIN FUNCTION FOR CREATING TARGET-SPECIFIC LM INPUT DATA FRAME
############################################################
#' @param starter_df a "starter" input data frame with patient/ sample and 
#' protein-of-interest (e.g. driver gene) information
#' @param downstream_target_df table of driver genes (columns) with gene 
#' targets we're testing for each (entries)
#' @param mutation_df table of mutation counts for each patient, for each 
#' gene in the genome
#' @param methylation_df table of methylation results (rows are proteins, columns 
#' are patients, entries are methylation values)
#' @param cna_df table of CNA per gene (rows are proteins, columns are patients, 
#' entries are CNA values)
#' @param expression_df table of expression (rows are genes, columns are patients, 
#' entries are expression values)
#' @param is_rank_or_quant_norm a T/F value indicating whether this is rank-
#' or quantile-normalized data (i.e. whether we need to include library size as 
#' an offset)
#' @param log_expression a T/F value indicating whether we need to take the
#' log2(exp + 1), or whether expression has already been transformed to a normal
#' distribution in some way
#' @param randomize a T/F value indicating whether or not we are randomizing 
#' expression
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. 
#' Possible values are "bucket_inclAmp", "bucket_exclAmp", "bucket_justAmp", 
#' "bucket_justDel" and "rawCNA"
#' @param meth_bucketing a T/F value indicating whether or not we are bucketing
#' methylation values
#' @param removeCis a T/F value indicating whether or not we are eliminating 
#' all cis pairings
#' @param removeMetastatic a T/F value indicating whether or not metastatic
#' samples were removed, for file naming
#' @param all_genes_id_conv an ID conversion file from BioMart, for use in removing 
#' cis pairings if needed
#' @param dataset either 'TCGA' or 'METABRIC', to denote what columns
#' we are referencing/ specific data types we are using
#' @param filename_labels a list of labels used for appropriate filename 
#' generation
#' @param outpath a string with the local path for writing files
#' @param debug a T/F value indicating if we are in debug mode and should 
#' use additional prints
create_lm_input_table <- function(starter_df, downstream_target_df, mutation_df, 
                                  methylation_df, cna_df, expression_df, 
                                  is_rank_or_quant_norm, log_expression, 
                                  randomize, cna_bucketing, meth_bucketing, 
                                  removeCis, removeMetastatic, all_genes_id_conv, 
                                  dataset, filename_labels, outpath, debug) {
  
  # If we are randomizing expression, do this now per-sample
  if(randomize) {
    random_exp_df <- apply(expression_df[,2:ncol(expression_df), with = F], 
                           MARGIN = 2,  function(x) dqrng::dqsample(x))
    expression_df <- cbind(expression_df[,'ensg_id'], data.table(random_exp_df))
  }
  
  # If remove cis is T, limit the downstream target DF to only trans pairings
  # using helper function in dyscovr_helper_functions.R
  if(removeCis) {
    downstream_target_df <- subset_to_trans_targets(prot_ensg, 
                                                    downstream_target_df,
                                                    all_genes_id_conv, debug)
  }
  
  # Use multithreading to run these jobs in parallel
  num_groups <- 8
  parallel_cluster <- makeCluster(num_groups, type = "SOCK", methods = F, outfile = "")
  
  # Tell R that we want to use these processes to do calculations
  setDefaultCluster(parallel_cluster)
  doParallel::registerDoParallel(parallel_cluster)
  
  list_of_rows <- split(downstream_target_df, f = seq(nrow(downstream_target_df)))
  
  lm_input_tables <- foreach(
    currentRow = list_of_rows, 
    # Export needed data and functions
    .export = c("methylation_df", "starter_df", "mutation_df", "cna_df", 
                "expression_df", "log_expression", "cna_bucketing", 
                "meth_bucketing", "debug", "dataset", "randomize", "removeCis", 
                "removeMetastatic", "filename_labels",
		"%fin%", "fill_targ_inputs", "create_lm_input_table_filename", 
		"filter_mut_by_uniprot", "filter_cna_by_ensg", "filter_meth_by_ensg",
		"get_mut_stat_targ", "get_cna_stat", "get_meth_stat", "get_exp_stat"), 
    .packages = c("dplyr", "data.table", "dqrng", "parallel", "fastmatch", "stringr"),
    .inorder = F) %dopar% {

      # Source the necessary files for the worker node
      #SOURCE_PATH <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/Dyscovr/"
      #source(paste(SOURCE_PATH, "general_important_functions.R", sep = "/"), local = T)
      #source(paste(SOURCE_PATH, "dyscovr_helper_functions.R", sep = "/"), local = T)
      
      # Get the target t_k's Swissprot & ENSG IDs
      targ <- unlist(strsplit(currentRow$swissprot, ";", fixed = T))[1]
      targ_ensg <- unlist(strsplit(currentRow$ensg, ";", fixed = T))
      
      if(debug) {print(targ)}
      
      tryCatch({
        lm_input_table <- fill_targ_inputs(starter_df = starter_df, 
                                           targ = targ, targ_ensg = targ_ensg, 
                                           mutation_targ_df = mutation_df,
                                           methylation_df = methylation_df,
                                           cna_df = cna_df, 
                                           expression_df = expression_df,
                                           log_expression = log_expression,
                                           cna_bucketing = cna_bucketing,
                                           meth_bucketing = meth_bucketing,
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
        # First, remove any columns that are entirely NA (e.g., a given target 
        # did not have CNA or methylation data reported)
        lm_input_table <- lm_input_table[, which(unlist(lapply(
          lm_input_table, function(x) !all(is.na(x))))), with = F]
        lm_input_table <- lm_input_table[, colSums(lm_input_table != 0, 
                                                   na.rm = T) > 0, with = F]
        
        # Then, if a row has NA, remove that whole row and predict without it
        lm_input_table <- na.omit(lm_input_table)
        
        if(debug) {
          print("LM Input Table") 
          print(head(lm_input_table))
          print(dim(lm_input_table))
        } 
        
        # Write this file to the appropriate outpath with the generated file name
        tryCatch({
          # Create an output file name for this LM input table using helper 
          # function in dyscovr_helper_functions.R
          lm_input_table_fn <- create_lm_input_table_filename(
            targ, filename_labels[["expression_df_name"]], cna_bucketing, 
            meth_bucketing, filename_labels[["meth_type"]], 
            filename_labels[["patient_df_name"]], randomize, 
            filename_labels[["patients_to_incl_label"]], removeCis, 
            removeMetastatic)
          
          # Write this file to appropriate outpath with the generated file name
          fwrite(lm_input_table, paste(outpath, paste0(lm_input_table_fn, ".csv"), 
                                       sep = "/"))
          gc()
        }, error = function(cond) {
          print(cond)
          print(length(lm_input_table_fn))
          print(paste(outpath, paste0(lm_input_table_fn, ".csv"), sep = "/"))
        })
      }
      # Ensure this is what is returned by this iteration of foreach
      lm_input_table <- distinct(lm_input_table)
  }
  # tell R that we don't need the processes anymore
  stopCluster(parallel_cluster)
  
  return(lm_input_tables)
}


############################################################
#### FILL IN TARGET GENE INPUTS TO TABLE
############################################################
#' Function takes a data frame of inputs for the given patient and the protein-
#' of-interest, as well as the data files for target (and its ID), and uses 
#' them to construct a full tabular input to Dyscovr for a given driver i-
#' target gene k pairing of the following form: Columns: Dyscovr input variables;
#' Rows: Patients 
#' @param starter_df a 'starter DF' that has all the patient and 
#' protein-of-interest characteristics 
#' @param targ the swissprot ID of target gene 
#' @param targ_ensg the ensembl ID of target gene
#' @param mutation_targ_df the mutation gene target DF
#' @param methylation_df the methylation DF
#' @param cna_df the copy number alteration DF
#' @param expression_df the gene expression DF
#' @param log_expression a T/ F value indicating whether we take the 
#' log2(exp + 1), or simply the expression value given in the DF
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. 
#' Possible values are "bucketInclAmp", "bucketExclAmp", and "rawCNA"
#' @param meth_bucketing a T/F value indicating whether or not we are bucketing
#' methylation values
#' @param debug a T/F value indicating if we are in debug mode and should 
#' use additional prints
#' @param dataset either 'TCGA' or 'METABRIC' to denote what columns
#' we are referencing/ specific data types we are using
fill_targ_inputs <- function(starter_df, targ, targ_ensg, mutation_targ_df,
                             methylation_df, cna_df, expression_df, 
                             log_expression, cna_bucketing, meth_bucketing, 
                             debug, dataset) {
  
  if(debug) {
    print(paste("Target:", targ))
    print(paste("Target ENSG:", targ_ensg))
  }
  
  # Begin by limiting the data frames to just this target gene using helper 
  # functions in dyscovr_helper_functions.R
  if (!(targ == "" | targ_ensg == "")) {
    
    mutation_targ_df <- filter_mut_by_uniprot(mutation_targ_df, targ) 
    cna_df <- filter_cna_by_ensg(cna_df, targ_ensg)
    methylation_df <- filter_meth_by_ensg(methylation_df, targ_ensg)
    
    exp_rows_to_keep <- unique(unlist(lapply(targ_ensg, function(x) 
      grep(x, expression_df$ensg_id))))
    expression_df <- expression_df[exp_rows_to_keep,]
    
    if(debug) {
      print(paste("Dim mutation targ df:", dim(mutation_targ_df)))
      print(paste("Dim cna df:", dim(cna_df)))
      print(paste("Dim expression df:", dim(expression_df)))
      print(paste("Dim methylation df:", dim(methylation_df)))
      print(!(expression_df[, .N] > 0) & (cna_df[, .N] > 0) & 
                (methylation_df[, .N] > 0) & (mutation_targ_df[, .N] > 0))
    }
    
    # Continue only if all of these data frames contain information about this 
    # target gene
    if(!(expression_df[, .N] > 0) & (cna_df[, .N] > 0) & 
                (methylation_df[, .N] > 0) & (mutation_targ_df[, .N] > 0)) {
      
      print(paste("Target not found in all DFs:", targ))
      return(NA)
    }
      
    # Loop through all samples to get the info for each target t_k column
    target_rows <- lapply(1:starter_df[, .N], function(i) {
      sample <- starter_df$sample_id[i]
      print(sample)
      
      # Is this target mutated? 
      mut_stat <- get_mut_stat_targ(mutation_targ_df, sample)
      
      # Is this target amplified or deleted?
      cna_stat <- get_cna_stat(cna_df, sample, cna_bucketing, dataset, F, NA)
      if(grepl("incl", cna_bucketing)) {cna_stat <- cna_stat[[1]]}    
      
      # Is this target methylated?
      meth_stat <- get_meth_stat(methylation_df, sample, meth_bucketing)
      if(meth_bucketing) {
        meth_stat <- as.integer(meth_stat[[1]])
        if(length(meth_stat) > 3) {meth_stat <- meth_stat[1:3]}
      }
      
      # What is the expression of this target in this sample in cancer?
      exp_stat <- unlist(get_exp_stat(expression_df, sample, dataset, 
                                      log_expression))
      
      if(debug) {
        print(paste("exp_stat:", exp_stat))
        print(paste("mut_stat:", mut_stat))
        print(paste("cna_stat:", cna_stat))
        print(paste("meth_stat:", meth_stat))
      }
      
      # Make row for this patient and return
      if(!meth_bucketing) {
        if ((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing))) {
          return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, 
                            "CNAStat_k" = cna_stat, "MethStat_k" = meth_stat))
        } else {
          return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, 
                            "CNAStat_k_b1" = cna_stat[1], 
                            "CNAStat_k_b2" = cna_stat[2], 
                            "CNAStat_k_b3" = cna_stat[3], 
                            "MethStat_k" = meth_stat))
        }
      } else {
        if ((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing))) {
          return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, 
                            "CNAStat_k" = cna_stat, "MethStat_k_b1" = meth_stat[1], 
                            "MethStat_k_b2" = meth_stat[2], 
                            "MethStat_k_b3" = meth_stat[3]))
        } else {
          return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, 
                            "CNAStat_k_b1" = cna_stat[1], 
                            "CNAStat_k_b2" = cna_stat[2], 
                            "CNAStat_k_b3" = cna_stat[3], 
                            "MethStat_k_b1" = meth_stat[1], 
                            "MethStat_k_b2" = meth_stat[2], 
                            "MethStat_k_b3" = meth_stat[3]))
        }
      }
    })
        
    tryCatch({
      colnam <- colnames(target_rows[[1]])
      targ_df <- rbindlist(target_rows)
      colnames(targ_df) <- colnam
      targ_df <- cbind(starter_df[,'sample_id'], targ_df)
      return(targ_df)
      
    }, error=function(cond){
      print("Unexpected error; target rows is not valid: ")
      print(target_rows)
      return(starter_df)
    })

  } else {
    if(debug) {
      print(paste("Target Uniprot or ENSG ID is missing,", targ))
    }
    return(NA)
  }
}


############################################################
#### CREATE OUTPATH 
############################################################
# Create an appropriate file outpath using helper function from 
# dyscovr_helper_functions.R
if(args$cancerType == "PanCancer") {
  if(!(args$specificTypes == "ALL")) {
    cancer_types <- unlist(strsplit(args$specificTypes, ";", fixed = T))
    outpath <- lapply(cancer_types, function(ct) {
      return(create_lm_input_file_outpath(cancerType = args$cancerType, 
                                          specificType = ct, test = test, 
                                          tester_name = args$tester_name, 
                                          run_name = args$run_name,
                                          dataset = args$dataset,
                                          outpath = SOURCE_PATH))
    })
    names(outpath) <- cancer_types
    
  } else {
    outpath <- create_lm_input_file_outpath(cancerType = args$cancerType, 
                                            specificType = "", test = test, 
                                            tester_name = args$tester_name, 
                                            run_name = args$run_name,
                                            dataset = args$dataset, 
                                            outpath = SOURCE_PATH)
  }
} else {
  outpath <- create_lm_input_file_outpath(cancerType = args$cancerType, 
                                          specificType = "", test = test, 
                                          tester_name = args$tester_name, 
                                          run_name = args$run_name,
                                          dataset = args$dataset, 
                                          outpath = SOURCE_PATH)
}


if(debug) {
  print(paste("Outpath:", outpath))
}

# NOTE: Output file name will be created dynamically per-gene


############################################################
#### OPTIONAL: LIMIT TARGETS FILE TO ONLY THOSE TARGETS WITHOUT A
#### FILE IN THE PROVIDED DIRECTORY
############################################################
#' Limits the target gene DF to only those targets that don't have an existing 
#' file in the outpath provided
#' @param outpath the path to which the LM input files will be written; a list of
#' paths if we are doing this pan-cancer
#' @param targets_DF the data frame with gene targets for LM input file creation
#' @param gene_name_index the numerical index of where the target gene name is 
#' found in a given file name
limit_to_targets_wo_existing_files <- function(outpath, targets_DF, 
                                               gene_name_index) {
  
  # Get all the existing files in this path, if any exist
  files_in_outpath <- list.files(paste0(outpath, "/"), pattern = "lm_input_")
  
  # Check if there are sub-directories in this path that might contain the file
  dirs_in_path <- list.dirs(outpath, recursive = F)
  if(T %fin% unlist(lapply(dirs_in_path, function(d) grepl("part_", d)))) {
    dirs_part <- dirs_in_path[grepl("part_", dirs_in_path)]
    for (d in dirs_part) {
      files_d <- list.files(paste0(d, "/"), pattern = "lm_input_")
      files_in_outpath <- c(files_in_outpath, files_d)
    }
  } 
  genes_in_outpath <- unlist(lapply(files_in_outpath, function(f) 
    unlist(strsplit(f, "_", fixed = T))[gene_name_index]))
  
  rows_to_keep <- unlist(lapply(1:targets_DF[, .N], function(i) {
    ids <- unlist(strsplit(as.character(targets_DF[i, 'swissprot']), ';', 
                           fixed = T))
    if(length(intersect(ids, genes_in_outpath)) > 0) {return(NA)}
    else {return(i)}
  }))
  rows_to_keep <- rows_to_keep[!is.na(rows_to_keep)]
  
  targets_DF <- targets_DF[rows_to_keep,]
  
  print("Targets DF, subsetted")
  print(dim(targets_DF))
  
  return(targets_DF)
}

# Get the index of where the target gene name will be found in each file name
gene_name_index <- length(unlist(strsplit(args$run_name, "_", fixed = T))) + 3
if(args$patientsOfInterest != "") {
  gene_name_index <- gene_name_index + 1
  outpath <- paste0(outpath, paste0("/subtype_files", 
                                    args$patientOfInterestLabel))
}

#targets_DF <- limit_to_targets_wo_existing_files(outpath, targets_DF, 
                                                 #gene_name_index)


############################################################
# IF RUNNING PER-CANCER WITH SPECIFIC CANCER TYPES, GET 
# THE PATIENTS OF THOSE CANCER TYPES 
############################################################
#' If we are running per-cancer on several individual cancer types, get a 
#' mapping between each cancer type and the IDs of the patients classified as 
#' that cancer type
#' @param specificTypes a semicolon-separated string of the cancer types we are 
#' looking at in this run
#' @param clinical_df a file name for a clinical data frame that can link 
#' patient IDs to their respective cancer types
get_patient_cancer_mapping <- function(specificTypes, clinical_df) {
  specific_types <- unlist(strsplit(specificTypes, ";", fixed = T))
  clinical_df <- fread(paste0(MAIN_PATH, paste0("Clinical/", clinical_df)), 
                       header = T)
  patient_cancer_mapping <- lapply(specific_types, function(ct) {
    pats <- clinical_df[grepl(ct, clinical_df$project_id),'case_submitter_id']
    pats_ids <- unlist(lapply(pats$case_submitter_id, function(pat) 
      unlist(strsplit(pat, "-", fixed = T))[3]))
    return(pats_ids)
  })
  names(patient_cancer_mapping) <- specific_types
  return(patient_cancer_mapping)
}

# Call function, if applicable
patient_cancer_mapping <- NA
if(args$cancerType == "PanCancer") {
  if(!(args$specificTypes == "ALL")) {
    patient_cancer_mapping <- get_patient_cancer_mapping(
      specificTypes = args$specificTypes, clinical_df = args$clinical_df)
    if(debug) {print(patient_cancer_mapping)}
  }
}


############################################################
#### CALL FUNCTION TO CREATE LM INPUT TABLES
############################################################
# Run gc to free up any loose memory
gc()

# Get the labels needed for dynamic filename creation and add to a list
filename_labels <- list("expression_df_name" = args$expression_df,
                        "meth_type" = args$meth_type,
                        "patient_df_name" = args$patient_df,
                        "test" = test, "tester_name" = args$tester_name, 
                        "run_name" = args$run_name,
                        "patients_to_incl_label" = args$patientsOfInterestLabel)


# Case 1: If we're running on just one cancer type or PanCancer (all cancers 
# together)
if(is.na(patient_cancer_mapping)) {
  
  # If needed, subset these data tables by the given patient IDs of interest 
  # using helper functions from dyscovr_helper_functions.R
  if(patients_of_interest != "") {
    patient_df <- subset_by_intersecting_ids(patients_of_interest, patient_df, 
                                             F, args$dataset)
    mutation_df <- subset_by_intersecting_ids(patients_of_interest, mutation_df, 
                                              T, args$dataset)
    methylation_df <- subset_by_intersecting_ids(patients_of_interest, 
                                                 methylation_df, T, args$dataset)
    cna_df <- subset_by_intersecting_ids(patients_of_interest, cna_df, T, 
                                         args$dataset)
    expression_df <- subset_by_intersecting_ids(patients_of_interest, 
                                                expression_df, T, args$dataset)
  }
  
  # If desired, exclude metastatic samples (dyscovr_helper_functions.R)
  if(removeMetastatic == T) {
    patient_df <- remove_metastatic_samples(patient_df, F)
    mutation_df <- remove_metastatic_samples(mutation_df, T)
    methylation_df <- remove_metastatic_samples(methylation_df, T)
    cna_df <- remove_metastatic_samples(cna_df, T)
    expression_df <- remove_metastatic_samples(expression_df, T)
  }
  
  # Create protein-of-interest/driver protein starter DF
  starter_df <- create_regulatory_prot_input_df(
    protein_ids_df = protein_ids_df, 
    patient_df = patient_df, 
    mutation_df = mutation_df, 
    methylation_df = methylation_df, 
    cna_df = cna_df, 
    debug = debug, 
    filename_labels = filename_labels,
    cna_bucketing = args$cna_bucketing, 
    meth_bucketing = meth_bucketing, 
    outpath = outpath, 
    num_pcs = args$num_pcs, 
    removeCis = removeCis, 
    removeMetastatic = removeMetastatic, 
    dataset = args$dataset, 
    run_query_genes_jointly = runQueriesJointly)
  
  # Create target gene input tables
  lm_input_tables <- create_lm_input_table(
    starter_df = starter_df, 
    downstream_target_df = targets_DF, 
    mutation_df = mutation_df,
    methylation_df = methylation_df, 
    cna_df = cna_df,
    expression_df = expression_df, 
    is_rank_or_quant_norm = is_rank_or_quant_norm,
    log_expression = log_expression,
    randomize = randomize,
    cna_bucketing = args$cna_bucketing,
    meth_bucketing = meth_bucketing,
    debug = debug,
    outpath = outpath, 
    removeCis = removeCis,
    removeMetastatic = removeMetastatic,
    all_genes_id_conv = all_genes_id_conv,
    dataset = args$dataset,
    filename_labels = filename_labels)
  
  
  # Case 2: If we're running on multiple cancer types individually
} else if (!(is.na(patient_cancer_mapping))) {
  # Subset all the input files for each cancer type and run individually
  input_df_list <- lapply(1:length(outpath), function(i) {
    
    # Get the patients for this cancer type
    patient_ids <- unique(unlist(patient_cancer_mapping[[i]]))
    #print(patient_ids)
    
    # Subset to just patients given in a particular set of interest, if provided
    if(patients_of_interest != "") {
      patient_ids <- patient_ids[patient_ids %fin% patients_of_interest]
    }
    
    # Get the outpath for this cancer type
    outpath_i <- outpath[[i]]
    
    ct <- names(patient_cancer_mapping)[i]
    print(ct)
    print(names(patient_df))
    print(head(patient_df[[1]]))
    patient_df_ct <- patient_df[[which(names(patient_df) == ct)]]

    # Subset files using patient_ids (using helper functions from 
    # dyscovr_helper_functions.R)
    patient_df_sub <- subset_by_intersecting_ids(patient_ids, patient_df_ct, F,
                                                 args$dataset)
    mutation_df_sub <- subset_by_intersecting_ids(patient_ids, mutation_df, T, 
                                                  args$dataset)
    methylation_df_sub <- subset_by_intersecting_ids(patient_ids, methylation_df, 
                                                     T, args$dataset)
    cna_df_sub <- subset_by_intersecting_ids(patient_ids, cna_df, T, args$dataset)
    expression_df_sub <- subset_by_intersecting_ids(patient_ids, expression_df, 
                                                    T, args$dataset)

    # If desired, exclude metastatic samples
    if(removeMetastatic == T) {
      patient_df_sub <- remove_metastatic_samples(patient_df_sub, F)
      mutation_df_sub <- remove_metastatic_samples(mutation_df_sub, T)
      methylation_df_sub <- remove_metastatic_samples(methylation_df_sub, T)
      cna_df_sub <- remove_metastatic_samples(cna_df_sub, T)
      expression_df_sub <- remove_metastatic_samples(expression_df_sub, T)
    }
    
    name <- names(patient_cancer_mapping)[i]
    protein_ids_df_spec <- protein_ids_df[[name]]
    
    # Create protein-of-interest/driver starter DF
    starter_df <- create_regulatory_prot_input_df(
      protein_ids_df = protein_ids_df_spec, 
      patient_df = patient_df_sub, 
      mutation_df = mutation_df_sub, 
      methylation_df = methylation_df_sub, 
      cna_df = cna_df_sub, 
      debug = debug, 
      filename_labels = filename_labels,
      cna_bucketing = args$cna_bucketing, 
      meth_bucketing = meth_bucketing, 
      outpath = outpath_i, 
      num_pcs = args$num_pcs, 
      removeCis = removeCis, 
      removeMetastatic = removeMetastatic, 
      dataset = args$dataset, 
      run_query_genes_jointly = runQueriesJointly)
    
    # Create target gene input tables
    lm_input_tables <- create_lm_input_table(
      starter_df = starter_df, 
      downstream_target_df = targets_DF, 
      mutation_df = mutation_df_sub,
      methylation_df = methylation_df_sub, 
      cna_df = cna_df_sub,
      expression_df = expression_df_sub, 
      is_rank_or_quant_norm = is_rank_or_quant_norm,
      log_expression = log_expression,
      randomize = randomize,
      cna_bucketing = args$cna_bucketing,
      meth_bucketing = meth_bucketing,
      debug = debug,
      outpath = outpath_i, 
      removeCis = removeCis,
      removeMetastatic = removeMetastatic,
      all_genes_id_conv = all_genes_id_conv,
      dataset = args$dataset,
      filename_labels = filename_labels)
    
    return(lm_input_tables)
  })
  
} else {
  print("Something went wrong with patient-cancer type mapping. Exiting now.")
  lm_input_table <- NA
}

