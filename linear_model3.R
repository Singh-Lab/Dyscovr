############################################################
### Linear Model 
### Written By: Sara Geraghty, July 2020
############################################################

#!/usr/bin/env Rscript

# Set the library path
library.path <- .libPaths()
#library.path <- "C:/Users/sarae/AppData/Local/R/win-library/4.1"
#library.path <- "C:/Program Files/R/R-4.1.1/library"
#print(library.path)

#library(biomaRt, lib.loc = library.path)
library(parallel, lib.loc = library.path, quietly = TRUE)
library(rlang, lib.loc = library.path, quietly = TRUE)
library(dplyr, lib.loc = library.path, quietly = TRUE)
library(broom, lib.loc = library.path, quietly = TRUE)
library(data.table, lib.loc = library.path, quietly = TRUE)
library(speedglm, lib.loc = library.path, quietly = TRUE)
library(argparse, lib.loc = library.path, quietly = TRUE)
library(stringr, lib.loc = library.path, quietly = TRUE)
library(dqrng, lib.loc = library.path, quietly = TRUE)
library(foreach, lib.loc = library.path, quietly = TRUE)
library(snow, lib.loc = library.path, quietly = TRUE)
library(doParallel, lib.loc = library.path, quietly = TRUE)
#library(glmnet, lib.loc = library.path)
#library(cli, lib.loc = library.path)
#library(nlme, lib.loc = library.path)
#library(lme4, lib.loc = library.path)
#library(lmerTest, lib.loc = library.path)

# Source other files needed
source_path <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/Sara_LinearModel/"
source(paste(source_path, "general_important_functions2.R", sep = "/"))
source(paste(source_path, "linear_model_helper_functions2.R", sep = "/"))
source(paste(source_path, "run_regularized_models2.R", sep = "/"))
source(paste(source_path, "run_bayesian_lasso.R", sep = "/"))


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

# The file path to the linear model input files from create_lm_input_files.R
parser$add_argument("--input_lm_filepath", default = "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/BRCA/tumor_only/lm_input_tables/top_drivers_0.05/eQTL",
                    type = "character", help = "The path to the input files created from create_lm_input_files.R for the given run")

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
parser$add_argument("--run_query_genes_jointly", default = "TRUE", type = "character",
                    help = "If TRUE, will run all the query regulatory proteins together in the same joint model, rather than running one model per query-target pair.")

# Randomization
parser$add_argument("--randomize", default = "FALSE", type = "character", 
                    help = "TRUE/FALSE for whether or not we perform a randomization test. Default is FALSE")

# Additional input files needed. The directory for the LM input files created from create_lm_input_tables.R can be inferred 
# from the other provided arguments.
parser$add_argument("--expression_df", default = "expression_tmm_CancerOnly_IntersectPatients.csv", 
                    type = "character",
                    help = "The name of the expression data frame input. [default %(default)s]")
parser$add_argument("--cna_df", default = "CNA_DF_CancerOnly_IntersectPatients.csv", 
                    type = "character",
                    help = "The name of the CNA data frame input; only needed if CNA values need adjusting. [default %(default)s]")

#parser$add_argument("--mutation_regprot_df", default = "iprotein_results_missense_CancerOnly_IntersectPatients.csv", 
 #                   type = "character",
  #                  help = "The name of the mutation regulatory protein data frame input. [default %(default)s]")
parser$add_argument("--patient_df", default = "combined_patient_sample_cibersort_total_frac_tmm_IntersectPatients.csv", 
                    type = "character",
                    help = "The name of the patient sample data frame input. [default %(default)s]")


# The decision for if/ how to bucket CNAs and methylation
parser$add_argument("--cna_bucketing", default = "bucket_inclAmp", 
                    type = "character",
                    help = "If/ the type of bucketing we use for CNA values. Default is bucket_inclAmp, but rawCNA, bucket_exclAmp, bucket_justAmp, and bucket_justDel are also options.")
parser$add_argument("--meth_bucketing", default = TRUE, type = "character",
                    help = "If/ the type of bucketing we use for methylation values. Default is FALSE.")

# Label the type of methylation data we are using (Beta, M, or Threshold)
parser$add_argument("--meth_type", default = "Beta", type = "character",
                    help = "The type methylation values we are using. Default is Beta, but other options are M and Threshold_X, where X is the threshold value.")

# Whether or not we are converting absolute CNA values to relative ones
parser$add_argument("--adjust_cna_rel", default = "FALSE", type = "character",
                    help = "Whether or not we are adjusting the CNA values to be log2(CN+1 / N+1) rather than log2(CN+1). Default is FALSE.")

# What type of test we are running, so we can write the results files to an appropriate directory. 
# Only necessary if --test is F.
parser$add_argument("--run_name", default = "cancerRelated", type = "character",
                    help = "Provide a name for the run, if not testing on just one regulatory protein, to write output files to appropriate directory. [default %(default)s]") 

# A name for the group of targets being tested, for labeling the output file.
parser$add_argument("--target_df", default = "metabolic_targets.csv", type = "character",
                    help = "Provide the name of the list of target genes; use to subset the given input DFs.")
parser$add_argument("--targets_name", default = "allGenes", type = "character",
                    help = "Provide a name for the group of targets being tested, in order to properly name the output file. [default %(default)s]")

# Whether or not we are including PEER factors and PCs as covariates in the model; for formula creation
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
parser$add_argument("--correctPerDriver", default = "TRUE", type = "character",
                    help = "A TRUE/FALSE value indicating whether we want to do q-value correction per-driver or across drivers. Default is TRUE (per driver).")

# Add a flag for the type of model to use (e.g. a mixed effects model, rather than a simple linear or regularized linear)
parser$add_argument("--model_type", default = "linear", type = "character",
                    help = "The type of model being used. Use 'mixed' to switch from linear to a mixed effects model.")
# Add a flag for adding a regularization method
parser$add_argument("--regularization", default = "None", type = "character",
                    help = "The name of the regularization method being used. Currently only implemented for L1 and L2. Default is None.")
# If we are using a regularized model (L1 or L2), select the method to evaluate significance.
parser$add_argument("--signif_eval_type", default = "randomization_perSamp", type = "character",
                    help = "If regularization is not None, uses this flag to determine method of evaluating significance. Takes either randomization_perTarg, randomization_perSamp, randomization_predictors, subsampling, or selectiveInference. Default is randomization_predictors.")
# If we are using a regularized model (L1 or L2), and a randomization method for evaluating significance, provide the number of randomizations. Default is 100.
parser$add_argument("--num_randomizations", default = 100, type = "integer",
                    help = "If regularization is not None, and signif_eval_type is some form of randomization, the number of randomizations is provided here.")
parser$add_argument("--fixed_lambda_for_rand", default = "FALSE", type = "character",
                    help = "If regularization is not None, and signif_eval_type is some form of randomization, this TRUE/FALSE value indicates whether we use a cross validation to determine the lambda for each randomization or use a fixed lambda value from the real LASSO.")

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
parser$add_argument("--select_drivers", default = "ALL", type ="character",
                    help = "Option to provide a semicolon-separated character string listing all driver Uniprot IDs to include, if LM input tables contain more driver covariates than is desired.")
parser$add_argument("--path_to_driver_lists", default = "", type ="character",
                    help = "Option to provide a path to files with drivers for individual cancer types. Will only be checked for a path is specificTypes is not ALL.")
parser$add_argument("--driver_file_prefix", default = "", type = "character",
                    help = "If path_to_driver_lists is not empty, this provides a prefix for the driver files we want to access for each cancer.")


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
useNumFunctCopies <- str2bool(args$useNumFunctCopies)
runRegprotsJointly <- str2bool(args$run_query_genes_jointly)
collinearity_diagn <- str2bool(args$collinearity_diagn)
inclResiduals <- str2bool(args$inclResiduals)
fixed_lambda_for_rand <- str2bool(args$fixed_lambda_for_rand)
correct_per_gene <- str2bool(args$correctPerDriver)
adjust_cna_rel <- str2bool(args$adjust_cna_rel)


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

# NOTE: if we opt to run all regulatory proteins jointly in the same model, we repeat these
# three covariates for all given regulatory proteins

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
#' @param list_of_input_dfs a list of input data frames from create_lm_input_tables.R,
#' where the list names are the swissprot IDs of the associated target gene
#' @param expression_df table of expression (rows are genes, columns are patients, 
#' entries are expression values)
#' @param analysis_type a string label that reads either "eQTL" or "meQTL" to 
#' determine what kind of model we should run
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. 
#' Possible values are "bucket_inclAmp", "bucket_exclAmp", "bucket_justAmp", 
#' "bucket_justDel" and "rawCNA"
#' @param meth_bucketing a TRUE/FALSE value indicating whether or not we are bucketing
#' methylation values
#' @param debug a TRUE/FALSE value indicating if we are in debug mode and should 
#' use additional prints
#' @param collinearity_diagn a TRUE/FALSE value indicating if we are running 
#' collinearity diagnostics
#' @param model_type the type of model being used (either linear, the default, or mixed)
#' @param regularization a string with the type of regularization to be used; currently
#' only implemented for "L2" or "None"
#' @param signif_eval_type if regularization is not "None", uses the given type of method
#' to evaluate significance (either some form of "randomization", "subsampling", or "selectiveInference")
#' @param num_randomizations the number of randomizations, if our signif_eval_type 
#' is some form of randomization
#' @param fixed_lambda_for_rand a TRUE/FALSE value indicating whether or not we are using a fixed
#' value of lambda for the randomizations (e.g. from the "real" lasso) or if we are using
#' cross-validation to determine a lambda value
#' @param run_query_genes_jointly a TRUE/ FALSE value indicating whether or not we 
#' want to run all query regulatory proteins jointly in the same model, or separately
#' (one model per query-target pair)
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3', to denote what columns
#' we are referencing/ specific data types we are using
#' @param inclResiduals a TRUE/FALSE value indicating whether or not we want to include 
#' residuals in our output DF, to evaluate model fit
run_linear_model <- function(list_of_input_dfs, expression_df, analysis_type, cna_bucketing, 
                             meth_bucketing, debug, collinearity_diagn, model_type, 
                             regularization, signif_eval_type, num_randomizations,
                             fixed_lambda_for_rand, run_query_genes_jointly, dataset, 
                             inclResiduals) {
  
  # Run through all targets in the downstream target data frame, running a regression model for each
  # Will be bound together into an output master DF
  master_df <- NA
  
  # Create a file to log our optimal lambda values
  #file.create("optimal_lambdas_real.txt")
  
  # Create shuffled input DFs if using an empirical approach to p-value generation
  shuffled_input_dfs <- NA
  if((model_type == "linear") & (regularization %fin% c("L1", "L2", "lasso", "ridge"))) {
    shuffled_input_dfs <- get_shuffled_input_dfs(signif_eval_type, lm_input_dfs[1], expression_df, 
                                                 num_randomizations)
  }
  
  # Use multithreading to run these jobs in parallel
  num_groups <- 6
  # Various types, including 'fork' and 'socket'. Fork is faster and runs on a Linux system (but not Windows)
  parallel_cluster <- makeCluster(6, type = "SOCK", methods = FALSE, outfile = "")
  #parallel_cluster <- makeCluster(6, type = "PSOCK", methods = FALSE, outfile = "")
  #parallel_cluster <- makeCluster(6, type = "FORK", methods = FALSE, outfile = "")
  
  # Tell R that we want to use these processes to do calculations
  setDefaultCluster(parallel_cluster)
  doParallel::registerDoParallel(parallel_cluster)
  
  #list_of_rows <- split(downstream_target_df, f = seq(nrow(downstream_target_df)))
  #print(length(list_of_input_dfs))
  
  master_df <- foreach(i = 1:length(list_of_input_dfs),
                        .export = c("list_of_input_dfs", "cna_bucketing", "meth_bucketing", 
                                    "analysis_type", "debug",  "dataset", "inclResiduals", 
                                    "model_type", "collinearity_diagn", "regularization", 
                                    "signif_eval_type", "num_randomizations", "fixed_lambda_for_rand", 
                                    "shuffled_input_dfs", "run_query_genes_jointly"), 
                        .combine = function(x,y) combine_tab(x, y)) %dopar% {
                        #.final = function(x) setNames(x, names(list_of_input_dfs))) %dopar% {
                          
      # Source the necessary files for the worker node
      source_path <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/Sara_LinearModel/"
      source(paste(source_path, "general_important_functions2.R", sep = "/"))
      source(paste(source_path, "linear_model_helper_functions2.R", sep = "/"))
      if (regularization != "None") {
        source(paste(source_path, "run_regularized_models2.R", sep = "/"))
        source(paste(source_path, "run_bayesian_lasso.R", sep = "/"))
      }
      
      lm_input_table = list_of_input_dfs[[i]]
      
      if(debug) {print(head(lm_input_table))}
      
      # Make sure any -Inf or Inf values are given NA and then omitted
      # If a column is >50% Inf or -Inf, remove that column first
      if(dataset == "METABRIC") {
        lm_input_table <- fix_lm_input_table(lm_input_table)
      }
      
      targ = names(list_of_input_dfs)[i]
      
      # Get the target t_k's Swissprot & ENSG IDs
      #targ <- unlist(strsplit(currentRow$swissprot, ";", fixed = TRUE))
      #targ_ensg <- unlist(strsplit(currentRow$ensg, ";", fixed = TRUE))
      
      if(debug) {print(targ)}
      #if(debug) {print(head(colnames(lm_input_table)))}
      
      if((lm_input_table[, .N] == 0) | (ncol(lm_input_table) == 0)) {
        print(paste("Null input data table for target:", paste0(targ, ". Returning NA.")))
        return(NA)
      }
      swissprot_ids <- unlist(lapply(colnames(lm_input_table)[grepl("MutStat_i", colnames(lm_input_table))], function(x)
        unlist(strsplit(x, "_", fixed = TRUE))[1]))
      input_regprots <- data.table("swissprot" = swissprot_ids)
      formula <- tryCatch(
        {
          construct_formula(lm_input_table, input_regprots, analysis_type, cna_bucketing, 
                            meth_bucketing, run_query_genes_jointly)          
        }, error = function(cond) {
          print(cond)
          return(NA)
        })
      if(debug) {print(paste("Formula:", formula))}
      if(formula == "ExpStat_k ~ ") {
        print("Truncated formula. Something went wrong.")
        print(head(colnames(lm_input_table)))
        print(swissprot_ids)
        return(NA)
      }
                          
      summary_table <- NA
                          
      if ((model_type == "linear") & (regularization == "None")) {
        
        print("Now running linear regression...")
        
        lm_fit <- tryCatch(
          {
            speedglm::speedlm(formula = formula, data = lm_input_table)  # Do not include library size as offset
            #speedglm::speedlm(formula = formula, offset = log2(Lib_Size), data = lm_input_table)
          }, error = function(cond) {
            message(paste("There was a problem with this linear model run:", cond))
            print(formula)
            print("Offending LM Input Table:")
            print(head(lm_input_table))
            print(dim(lm_input_table))
            print(any(is.na(lm_input_table)))
            print(targ)
            return(NA)
          })
        
        if(!is.na(lm_fit)) {
          summary_table <- tidy(lm_fit)
          summary_table <- as.data.table(summary_table[, colSums(is.na(summary_table)) < nrow(summary_table)])
          rm(lm_fit)
        } 
                            
      } else if ((model_type == "linear") & (regularization == "L2" | regularization == "ridge" | 
                                             regularization == "L1" | regularization == "lasso" |
                                             regularization == "bayesian.bgl" | regularization == "bayesian.bglss")) {
        
        # Create a file to hold the lambda values
        summary_table <- run_regularization_model(formula = formula, lm_input_table = lm_input_table, 
                                           type = regularization, debug = debug, meth_bucketing = meth_bucketing, 
                                           cna_bucketing = cna_bucketing, signif_eval_type = signif_eval_type,
                                           shuffled_input_dfs = shuffled_input_dfs, ensg = targ_ensg,
                                           num_randomizations = num_randomizations,
                                           fixed_lambda_for_rand = fixed_lambda_for_rand)

      } else if (model_type == "mixed") {
        # Our random effects variable is our genotype PCs; if we are not including them,
        # run the simple linear model instead
        # TODO: implement for cancer type/ subtype
        if (num_pcs == 0) {
          print("No genotype PCs are being used in regression. Random effects variables are currently only implemented for genotype PCs. Running a linear regression.")
          lm_fit <- tryCatch(
            {
              speedglm::speedlm(formula = formula, data = lm_input_table)  # Do not include library size as offset
              #speedglm::speedlm(formula = formula, offset = log2(Lib_Size), data = lm_input_table)
            }, error = function(cond) {
              message("There was a problem with this linear model run:")
              message(cond)
              return(NA)
            })
        } else {
          # Build a linear mixed model
          lm_fit <- tryCatch(
            {
              lmer(formula = formula, data = lm_input_table, REML = TRUE)      
            }, error = function(cond) {
              message("There was a problem with this linear mixed model run:")
              message(cond)
              return(NA)
            })
        }
        if(!is.na(lm_fit)) {
          summary_table <- as.data.table(tidy(lm_fit))
          rm(lm_fit)
        }
      } else {
        if(is.na(formula)) {
          print("Encountered an error. Formula is NA.")
          return(NA)
        } else {
          print(paste("Model type", paste(model_type, "does not match any of the implemented options. Please provide one of these options.")))
          return(NA)
        }
      }
      
      if(inclResiduals & (!(is.na(summary_table)))) {
        residuals <- as.numeric(lm_input_table$ExpStat_k) - 
          as.numeric(unlist(predict.speedlm(lm_fit, newdata = lm_input_table)))
        summary_table$Residual <- paste(residuals, collapse = ";")
      }
                          
      # Add a column for the regulatory protein and target gene to ID them
      if(!is.na(summary_table)) {
        summary_table$T_k <- paste(targ, collapse = ";")
        if(!run_query_genes_jointly) {
          summary_table$R_i <- regprot
        } else {
          # Remove rows in which the target and the driver are the same
          ri <- unlist(lapply(summary_table$term, function(x)
            unlist(strsplit(x, "_", fixed = TRUE))[1]))
          summary_table$R_i <- ri
          rows_to_keep <- unlist(lapply(1:length(ri), function(i) 
            if(grepl(ri[i], summary_table[i, 'T_k'])) {return(NA)}
            else {return(i)}))
          rows_to_keep <- rows_to_keep[!is.na(rows_to_keep)]
          summary_table <- summary_table[rows_to_keep,]
        }
      }
                            
      # Run collinearity diagnostics 
      if(collinearity_diagn) {
        print("Running Collinearity Diagnostics")
        #fit_lm <- lm(formula = formula, data = lm_input_table)
        fit_lm <- lm(formula = formula, data = lm_input_table)
        source(paste(getwd(), "run_collinearity_diagnostics.R", sep = "/"), local = TRUE)
      }
                
      # Make sure that the summary table for this gene is returned
      #if(length(summary_table) == 0) {summary_table <- NA}
      if(debug) {
        print("Summary Table")
        print(head(summary_table))
      }
      
      # Try to clear up some memory
      gc()
      
      # Return summary table
      summary_table <- summary_table
  }
  print(head(master_df))
  print(dim(master_df))
  
  # tell R that we don't need the processes anymore
  stopCluster(parallel_cluster)
  
  
  # Now we have a list of the results tables, one table per regulatory protein, 
  # with entries for each target. Bind these all together into one master table.
  #if(debug) {
  #  print("Results DF list")
  #  print(head(results_df_list))
  #}
  
  # Check that none of the DF entries are NA or empty data.tables
  #results_df_list <- results_df_list[!is.na(results_df_list)]
  #results_df_list <- Filter(function(dt) nrow(dt) != 0, results_df_list)
  
  #master_df <- rbindlist(results_df_list, use.names = TRUE, fill = TRUE)
  #master_df <- do.call(rbind, results_df_list)  
  #master_df <- na.omit(master_df)
  
  if(debug) {
    print("Master DF")
    print(head(master_df))
  }
        
  # Return this master DF
  return(master_df)
}



############################################################
#### CREATE OUTPATH AND OUTFILE NAME 
############################################################
# Create an appropriate file outpath
if(args$cancerType == "PanCancer") {
  if(!(args$specificTypes == "ALL")) {
    cancer_types <- unlist(strsplit(args$specificTypes, ";", fixed = TRUE))
    outpath <- lapply(cancer_types, function(ct) 
      create_file_outpath(cancerType = args$cancerType, specificType = ct, test = test,
                          tester_name = args$tester_name, run_name = args$run_name,
                          tumNormMatched = tumNormMatched, QTLtype = args$QTLtype,
                          dataset = args$dataset))
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
                                cna_bucketing = args$cna_bucketing, #mutation_regprot_df_name = args$mutation_regprot_df, 
                                meth_bucketing = meth_bucketing, meth_type = args$meth_type, 
                                patient_df_name = args$patient_df[1], num_PEER = args$num_PEER, 
                                num_pcs = args$num_pcs, randomize = randomize, covs_to_incl_label = args$select_args_label,
                                patients_to_incl_label = args$patientsOfInterestLabel, removeCis = FALSE,
                                model_type = args$model_type, removeMetastatic = TRUE, 
                                regularization = args$regularization, signif_eval_type = args$signif_eval_type)[1]

if(debug) {
  print(paste("Outfile Name:", outfn))
}


############################################################
# IMPORT COVARIATES
############################################################
# Convert covariates to include from a semicolon-separated character to a vector
covs_to_incl <- c("ALL")

tryCatch({
  covs_to_incl <- unlist(strsplit(args$select_args, ";", fixed = TRUE))
}, error = function(cond) {
  print(cond)
  print("Invalid input covariate list. Make sure covariates are semi-colon separated. Running with all covariates.")
})


############################################################
# SET PATHS AND IMPORT TARGET GENE DF
############################################################
# If needed, subset this list based on the given targets in the target DF
input_file_path <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/"
main_path <- NA
prot_path <- NA
targ_path <- NA

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
target_df <- fread(paste0(targ_path, args$target_df), header = TRUE)

target_uniprot_ids <- unique(unlist(lapply(target_df$swissprot, function(x)
  unlist(strsplit(x, ";", fixed = TRUE)))))
target_uniprot_ids <- target_uniprot_ids[!is.na(target_uniprot_ids)]


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
#' @param main_path local path to relevant files
#' @param thres a threshold for number of samples to include a cancer type, by default 25
get_patient_cancer_mapping <- function(specificTypes, clinical_df, main_path, thres = 25) {
    specific_types <- unlist(strsplit(specificTypes, ";", fixed = TRUE))
    clinical_df <- fread(paste(main_path, paste("clinical/", clinical_df, sep = ""), sep = ""), 
                         header = TRUE)
    patient_cancer_mapping <- lapply(specific_types, function(ct) {
      pats <- clinical_df[grepl(ct, clinical_df$project_id),'case_submitter_id']
      pats_ids <- unlist(lapply(pats$case_submitter_id, function(pat) 
        unlist(strsplit(pat, "-", fixed = TRUE))[3]))
      if(length(pats_ids) > thres) {return(pats_ids)}
      else {return(NA)}
    })
    names(patient_cancer_mapping) <- specific_types
    patient_cancer_mapping <- patient_cancer_mapping[!is.na(patient_cancer_mapping)]
    return(patient_cancer_mapping)
}

# Call function, if applicable
patient_cancer_mapping <- NA
if(args$cancerType == "PanCancer") {
  if(!(args$specificTypes == "ALL")) {
    patient_cancer_mapping <- get_patient_cancer_mapping(specificTypes = args$specificTypes,
                                                         clinical_df = args$clinical_df,
                                                         main_path = main_path)
    if(debug) {print(patient_cancer_mapping)}
  }
}
print(patient_cancer_mapping)

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
# IMPORT LIST OF PRE-CREATED LM INPUT FILES
############################################################
lm_input_fns <- NA

# Input LM filepath
input_lm_filepath <- args$input_lm_filepath
if(args$patientsOfInterest != "") { 
  input_lm_filepath <- paste(input_lm_filepath, "subtype_files", sep = "/")
}

# Get the index of where the target gene name will be found in each file name
gene_name_index <- length(unlist(strsplit(args$run_name, "_", fixed = TRUE))) + 3

# If we need to adjust CNA values, import the DF
cna_df <- NA
if(adjust_cna_rel) {
  cna_df <- fread(paste(main_path, paste("cna/", args$cna_df, sep = ""), sep = ""),
                  header = TRUE)
  colnames(cna_df)[1] <- "ensg_id"
} 

# If for one cancer type or PanCancer, get the file names
if(length(patient_cancer_mapping) == 1) {
  if(args$patientsOfInterest != "") {
    gene_name_index <- gene_name_index + 1
    lm_input_fns <- list.files(input_lm_filepath, pattern = args$patientsOfInterestLabel, 
                               recursive = F)
  } else {
    lm_input_fns <- list.files(input_lm_filepath, pattern = "lm_input_", recursive = F)
  }
  # Check if there are sub-directories in this path that might contain the file
  dirs_in_path <- list.dirs(input_lm_filepath, recursive = F)
  if(TRUE %in% unlist(lapply(dirs_in_path, function(d) grepl("part_", d)))) {
    dirs_part <- dirs_in_path[grepl("part_", dirs_in_path)]
    gene_name_index <- gene_name_index + 1
    for (d in dirs_part) {
      files_d <- c()
      if(args$patientsOfInterest != "") {
        files_d <- list.files(paste0(d, "/"), pattern = args$patientsOfInterestLabel, 
                              recursive = F)
      } else {
        files_d <- list.files(paste0(d, "/"), pattern = "lm_input_", 
                              recursive = F)
      }
      part <- unlist(strsplit(d, "/", fixed = TRUE))
      the_part <- part[length(part)]
      files_d <- paste(the_part, files_d, sep = "/")
      lm_input_fns <- c(lm_input_fns, files_d)
      
    }
  } 
} else {   # For multiple individual cancer types, get the file names in a list
  spl_fp <- unlist(strsplit(input_lm_filepath, "/", fixed = T))
  index_pc <- which(spl_fp == "PanCancer")
  pt1 <- paste(spl_fp[1:index_pc], collapse = "/")
  pt2 <- paste(spl_fp[(index_pc+1):length(spl_fp)], collapse = "/")
  input_lm_filepath <- c()
  lm_input_fns <- list()
  
  if(args$patientsOfInterest != "") {
    gene_name_index <- gene_name_index + 1
    for (ct in names(patient_cancer_mapping)) {
      new_fp <- paste(pt1, paste(ct, pt2, sep = "/"), sep = "/")
      print(new_fp)
      input_lm_filepath <- c(input_lm_filepath, new_fp)
      new_fns <- list.files(new_fp, pattern = args$patientsOfInterestLabel, 
                            recursive = F)
      lm_input_fns <- append(lm_input_fns, list(new_fns))
    }
  } else {
    for (ct in names(patient_cancer_mapping)) {
      new_fp <- paste(pt1, paste(ct, pt2, sep = "/"), sep = "/")
      print(new_fp)
      input_lm_filepath <- c(input_lm_filepath, new_fp)
      new_fns <- list.files(new_fp, pattern = "lm_input_", 
                            recursive = F)
      lm_input_fns <- append(lm_input_fns, list(new_fns))
    }
  }
  names(lm_input_fns) <- names(patient_cancer_mapping)
  if(debug) {print(lm_input_fns)}
}
print(input_lm_filepath)

# Get the names of all the genes we have LM input files for
gene_names <- NA
if(length(patient_cancer_mapping) == 1) {
  gene_names <- unlist(lapply(lm_input_fns, function(fn) 
    unlist(strsplit(fn, "_", fixed = TRUE))[gene_name_index]))
  if(debug) {print(head(gene_names))}
  
  print("# of input LM DFs")
  print(length(lm_input_fns))
  
} else {
  gene_names <- lapply(lm_input_fns, function(ct_fn) {
    print(paste("# of input LM DFs", length(ct_fn)))
    print(head(ct_fn))
    
    gns <- unlist(lapply(ct_fn, function(fn) 
      unlist(strsplit(fn, "_", fixed = TRUE))[gene_name_index]))
    if(debug) {print(head(gns))}
    return(gns)
  })
  names(gene_names) <- names(patient_cancer_mapping)
}

lm_input_fns_sub <- lm_input_fns
gene_names_sub <- gene_names

if(length(patient_cancer_mapping) == 1) {
  to_keep <- which(gene_names %fin% target_uniprot_ids)
  lm_input_fns_sub <- lm_input_fns[to_keep]
  gene_names_sub <- gene_names[to_keep]
  
  print("# of input LM DFs, subsetted")
  print(length(lm_input_fns_sub))
  
} else {
  for (i in 1:length(lm_input_fns)) {
    input_fns <- lm_input_fns[[i]]
    gns <- gene_names[[i]]
    
    to_keep <- which(gns %fin% target_uniprot_ids)
    lm_input_fns_sub[[i]] <- input_fns[to_keep]
    gene_names_sub[[i]] <- gns[to_keep]
    
    print("# of input LM DFs, subsetted")
    print(length(lm_input_fns_sub[[i]]))
  }
}

# If the number of genes is greater than a particular threshold N (e.g. 8K genes), 
# split into multiple partial runs so that we don't exceed memory constraints
split_target_list_flag <- FALSE
N <- 8000
lm_input_dfs <- NA

if(length(patient_cancer_mapping) == 1) {
  lm_input_dfs <- call_import_lm_file(lm_input_fns_sub, input_lm_filepath, randomize,
                                      gene_names_sub, args$select_drivers, covs_to_incl,
                                      args$num_pcs, args$num_PEER, args$QTLtype, 
                                      cna_df, adjust_cna_rel, N)
  if(length(lm_input_fns_sub) > N) {split_target_list_flag <- TRUE}
  
} else {
  split_target_list_flag <- c()
  lm_input_dfs <- lapply(1:length(lm_input_fns_sub), function(i) {
    ct_fns <- lm_input_fns_sub[[i]]
    ct_gens <- gene_names_sub[[i]]
    ct_filepath <- input_lm_filepath[i]
    print(ct_filepath)
    input_dfs <- call_import_lm_file(ct_fns, ct_filepath, randomize, ct_gens, 
                                    args$select_drivers, covs_to_incl, args$num_pcs, 
                                    args$num_PEER, args$QTLtype, cna_df, adjust_cna_rel, N)
    if(length(ct_fns) > N) {split_target_list_flag <- c(split_target_list_flag, TRUE)}
    else {split_target_list_flag <- c(split_target_list_flag, FALSE)}
    return(input_dfs)
  })
  names(lm_input_dfs) <- names(lm_input_fns_sub)
}



############################################################
#### CALL FUNCTION
############################################################
# Run gc to free up any loose memory
gc()

# Run the function itself

# Case 1: If we're running on just one cancer type or PanCancer
if(length(patient_cancer_mapping) == 1) {
  
  # If needed, subset these data tables by the given patient IDs of interest
  if(patients_of_interest != "") {
    lm_input_dfs <- lapply(lm_input_dfs, function(df) 
      subset_by_intersecting_ids(patients_of_interest, df, FALSE, tumNormMatched, args$dataset))
  }
  
  # Filter samples with outlier expression for each target gene
  #lm_input_dfs <- lapply(lm_input_dfs, function(df) 
    #filter_expression_outlier_samples3(df))
    #collapse_expression_outlier_samples(df))

  master_df <- run_linear_model(list_of_input_dfs = lm_input_dfs, 
                                expression_df = expression_df, 
                                analysis_type = args$QTLtype,
                                cna_bucketing = args$cna_bucketing,
                                meth_bucketing = meth_bucketing,
                                debug = debug,
                                collinearity_diagn = collinearity_diagn,
                                regularization = args$regularization,
                                signif_eval_type = args$signif_eval_type,
                                num_randomizations = args$num_randomizations,
                                fixed_lambda_for_rand = fixed_lambda_for_rand,
                                model_type = args$model_type,
                                run_query_genes_jointly = runRegprotsJointly,
                                dataset = args$dataset,
                                inclResiduals = inclResiduals)
  
  #fwrite(master_df, paste(outpath, paste(outfn, "_part1.csv", sep = ""), sep = "/"))
  
  if(split_target_list_flag) {
    
    # Clear memory from the previous batch
    rm(lm_input_dfs)
    gc()
    
    # Run again on the next batch
    lm_input_dfs <- import_lm_input_fns(lm_input_fns_sub[(N+1):length(lm_input_fns_sub)], input_lm_filepath, 
                                        randomize, gene_names_sub[(N+1):length(lm_input_fns_sub)], 
                                        args$select_drivers, covs_to_incl, args$num_pcs, args$num_PEER, args$QTLtype,
                                        cna_df, adjust_cna_rel)
    
    # If needed, subset these data tables by the given patient IDs of interest
    if(patients_of_interest != "") {
      lm_input_dfs <- lapply(lm_input_dfs, function(df) 
        subset_by_intersecting_ids(patients_of_interest, df, FALSE, tumNormMatched, args$dataset))
    }
    
    # Filter samples with outlier expression for each target gene
    #lm_input_dfs <- lapply(lm_input_dfs, function(df) 
      #filter_expression_outlier_samples3(df))
      #collapse_expression_outlier_samples(df))
    
    master_df_2 <- run_linear_model(list_of_input_dfs = lm_input_dfs, 
                                  expression_df = expression_df, 
                                  analysis_type = args$QTLtype,
                                  cna_bucketing = args$cna_bucketing,
                                  meth_bucketing = meth_bucketing,
                                  debug = debug,
                                  collinearity_diagn = collinearity_diagn,
                                  regularization = args$regularization,
                                  signif_eval_type = args$signif_eval_type,
                                  num_randomizations = args$num_randomizations,
                                  fixed_lambda_for_rand = fixed_lambda_for_rand,
                                  model_type = args$model_type,
                                  run_query_genes_jointly = runRegprotsJointly,
                                  dataset = args$dataset,
                                  inclResiduals = inclResiduals)
    
    #fwrite(master_df_2, paste(outpath, paste(outfn, "_part2.csv", sep = ""), sep = "/"))
    
    master_df <- rbind(master_df, master_df_2)
  }
  
  
# Case 2: If we're running on multiple cancer types separately
} else if (length(patient_cancer_mapping) > 1) {
  # Handle the case of running multiple cancer types per-cancer
  # Subset all the input files for each cancer type and run separately
  master_df_list <- lapply(1:length(outpath), function(i) {
    
    cancer_type <- names(patient_cancer_mapping)[i]
    
    print(paste("Now processing...", cancer_type))
    
    # Get the patients for this cancer type
    if(!is.na(cancer_type)) {
      patient_ids <- unique(unlist(patient_cancer_mapping[[i]]))
      print(head(patient_ids))
      
      if(patients_of_interest != "") {
        patient_ids <- patient_ids[patient_ids %fin% patients_of_interest]
      }
      
      # Subset input DFs to these patients of interest
      lm_input_dfs_ct <- lm_input_dfs[[i]]
      
      lm_input_dfs_sub <- lapply(lm_input_dfs_ct, function(df) {
        df_sub <- subset_by_intersecting_ids(patient_ids, df, FALSE, tumNormMatched, args$dataset)
        if(!(nrow(df_sub) == 0)) {return(df_sub)}
        else {return(NA)}
      })
      lm_input_dfs_sub <- lm_input_dfs_sub[!is.na(lm_input_dfs_sub)]
      if(length(lm_input_dfs_sub) == 0) {return(NA)}
      
      # Filter samples with outlier expression for each target gene
      #lm_input_dfs_sub <- lapply(lm_input_dfs_sub, function(df) 
      #filter_expression_outlier_samples3(df))
      #collapse_expression_outlier_samples(df))
      
      # Subset input DFs to drivers that pass particular threshold
      #if(args$path_to_driver_lists != "") {
      #  driver_filename <- paste0(args$path_to_driver_lists, 
      #                            paste0(args$driver_file_prefix, 
      #                                   paste0(cancer_type, ".csv")))
      #  driver_file <- fread(driver_filename, header = T)
      #  select_driver_vect <- unique(unlist(driver_file$swissprot_ids))
      #  if(debug) {print(select_driver_vect)}
      #  
      #  lm_input_dfs_sub <- lapply(lm_input_dfs_sub, function(df) {
      #    drivers_in_tab <- unique(unlist(lapply(colnames(df)[grepl("MutStat_i", colnames(df))], function(x)
      #      unlist(strsplit(x, "_", fixed = TRUE))[1])))
      #    drivers_to_exclude <- setdiff(drivers_in_tab, select_driver_vect)
      #    drivers_to_exclude_covs <- unlist(lapply(drivers_to_exclude, function(d) 
      #      return(colnames(df)[grepl(d, colnames(df))])))
      #    df_sub <- df[, !(colnames(df) %fin% drivers_to_exclude_covs), with = FALSE]
      #    return(df_sub)
      #  })
      #}
      
      # Get the outpath for this cancer type
      #outpath_i <- outpath[i]
      
      # Run the model with these subsetted files
      master_df <- run_linear_model(list_of_input_dfs = lm_input_dfs_sub, 
                                    expression_df = expression_df, 
                                    analysis_type = args$QTLtype,
                                    cna_bucketing = args$cna_bucketing,
                                    meth_bucketing = meth_bucketing,
                                    debug = debug,
                                    collinearity_diagn = collinearity_diagn,
                                    regularization = args$regularization,
                                    signif_eval_type = args$signif_eval_type,
                                    num_randomizations = args$num_randomizations,
                                    fixed_lambda_for_rand = fixed_lambda_for_rand,
                                    model_type = args$model_type,
                                    run_query_genes_jointly = runRegprotsJointly,
                                    dataset = args$dataset,
                                    inclResiduals = inclResiduals)
      return(master_df)
    } else {return(NA)}
  })
  
  if(TRUE %fin% split_target_list_flag) {
    
    # Clear memory from the previous batch
    rm(lm_input_dfs)
    rm(lm_input_dfs_ct)
    rm(lm_input_dfs_sub)
    gc()
    
    # Run again on the next batch
    lm_input_dfs <- lapply(1:length(lm_input_fns_sub), function(i) {
      ct_fns <- lm_input_fns_sub[[i]]
      ct_gens <- gene_names_sub[[i]]
      input_dfs <- import_lm_input_fns(ct_fns[(N+1):length(ct_fns)], input_lm_filepath, 
                                       randomize, ct_gens[(N+1):length(ct_fns)], 
                                       args$select_drivers, covs_to_incl, args$num_pcs, 
                                       args$num_PEER, args$QTLtype,
                                       cna_df, adjust_cna_rel)
      return(input_dfs)
    })
    names(lm_input_dfs) <- names(lm_input_fns_sub)
    
    master_df_list_2 <- lapply(1:length(outpath), function(i) {
      
      if(split_target_list_flag) {
        cancer_type <- names(patient_cancer_mapping)[i]
        
        # Get the patients for this cancer type
        if(!is.na(cancer_type)) {
          patient_ids <- unique(unlist(patient_cancer_mapping[[i]]))
          print(head(patient_ids))
          
          if(patients_of_interest != "") {
            patient_ids <- patient_ids[patient_ids %fin% patients_of_interest]
          }
          
          # Subset input DFs to these patients of interest
          lm_input_dfs_ct <- lm_input_dfs[[i]]
          lm_input_dfs_sub <- lapply(lm_input_dfs_ct, function(df) {
            df_sub <- subset_by_intersecting_ids(patient_ids, df, FALSE, tumNormMatched, args$dataset)
            if(!(nrow(df_sub) == 0)) {return(df_sub)}
            else {return(NA)}
          })
          lm_input_dfs_sub <- lm_input_dfs_sub[!is.na(lm_input_dfs_sub)]
          if(length(lm_input_dfs_sub) == 0) {return(NA)}
          
          # Filter samples with outlier expression for each target gene
          #lm_input_dfs_sub <- lapply(lm_input_dfs_sub, function(df) 
          #filter_expression_outlier_samples3(df))
          #collapse_expression_outlier_samples(df))
          
          # Subset input DFs to drivers that pass particular threshold
          #if(args$path_to_driver_lists != "") {
          #  driver_filename <- paste0(args$path_to_driver_lists, 
          #                            paste0(args$driver_file_prefix, 
          #                                   paste0(cancer_type, ".csv")))
          #  driver_file <- fread(driver_filename, header = T)
          #  select_driver_vect <- unique(unlist(driver_file$swissprot_ids))
          #  if(debug) {print(select_driver_vect)}
          #  
          #  lm_input_dfs_sub <- lapply(lm_input_dfs_sub, function(df) {
          #    drivers_in_tab <- unique(unlist(lapply(colnames(df)[grepl("MutStat_i", colnames(df))], function(x)
          #      unlist(strsplit(x, "_", fixed = TRUE))[1])))
          #    drivers_to_exclude <- setdiff(drivers_in_tab, select_driver_vect)
          #    drivers_to_exclude_covs <- unlist(lapply(drivers_to_exclude, function(d) 
          #      return(colnames(df)[grepl(d, colnames(df))])))
          #    df_sub <- df[, !(colnames(df) %fin% drivers_to_exclude_covs), with = FALSE]
          #    return(df_sub)
          #  })
          #}
          
          master_df_2 <- run_linear_model(list_of_input_dfs = lm_input_dfs_sub, 
                                          expression_df = expression_df, 
                                          analysis_type = args$QTLtype,
                                          cna_bucketing = args$cna_bucketing,
                                          meth_bucketing = meth_bucketing,
                                          debug = debug,
                                          collinearity_diagn = collinearity_diagn,
                                          regularization = args$regularization,
                                          signif_eval_type = args$signif_eval_type,
                                          num_randomizations = args$num_randomizations,
                                          fixed_lambda_for_rand = fixed_lambda_for_rand,
                                          model_type = args$model_type,
                                          run_query_genes_jointly = runRegprotsJointly,
                                          dataset = args$dataset,
                                          inclResiduals = inclResiduals)
          return(master_df_2)
        } else {return(NA)}
      } else {return(NA)}
    })
    
    master_df_list <- lapply(1:length(master_df_list), function(i) {
      master_df <- master_df_list[[i]]
      master_df_2 <- master_df_list_2[[i]]
      if(length(master_df_2) > 1) {
        return(rbind(master_df, master_df_2))
      } else {return(master_df)}
    })
  }
  #na_ind <- which(is.na(master_df_list))
  #master_df_list <- master_df_list[-na_ind]
  #names(master_df_list) <- names(patient_cancer_mapping)[-na_ind]
  #outpath <- outpath[-na_ind]

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
  if(debug) {print(head(master_df))}
  try(master_df <- master_df[order(p.value)])
  #try(master_df <- setorder(master_df, cols = "p.value", order = -1, na.last = TRUE))  
  
  # Write the results to the given file, making the directory if it does not already exist
  tryCatch({
    dir.create(outpath, showWarnings = FALSE)
  }, error = function(cond) {
    print(cond)
    print("Invalid outpath. Cannot create specified directory.")
  })
  fwrite(master_df, paste(outpath, paste(outfn, ".csv", sep = ""), sep = "/"))
  
  # Limit the data frame to just the term of interest (typically either MutStat_i or CNAStat_i)
  if(!useNumFunctCopies) {
    master_df_mut <- master_df[grepl("MutStat_i", master_df$term),]
    master_df_cna <- master_df[grepl("CNAStat_i", master_df$term),]
    
    # Write these to files as well 
    fwrite(master_df_mut, paste(outpath, paste(outfn, "_MUT.csv", sep = ""), sep = "/")) 
    fwrite(master_df_cna, paste(outpath, paste(outfn, "_CNA.csv", sep = ""), sep = "/"))  
  } else {
    master_df_fnc <- master_df[grepl("FNCStat_i", master_df$term),]
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


if(length(patient_cancer_mapping) == 1) {
  outpath <- as.character(unlist(outpath[[1]]))
  print(outpath)
  outfn <- as.character(unlist(outfn[[1]]))
  master_df <- na.omit(master_df)
  master_df <- process_raw_output_df(master_df = master_df, outpath = outpath, outfn = outfn, 
                                     collinearity_diagn = collinearity_diagn, debug = debug,
                                     randomize = randomize, useNumFunctCopies = useNumFunctCopies)
  if(!useNumFunctCopies) {
    master_df_mut <- master_df[grepl("MutStat_i", master_df$term),]
    master_df_cna <- master_df[grepl("CNAStat_i", master_df$term),]
  } else {
    master_df_fnc <- master_df[master_df$term == "FNCStat_i",]
  }
  
  # Call the file to create output visualizations
  outpath_curr <- outpath
  source(paste(source_path, "process_LM_output.R", sep = "")) 
  
} else {
  for (i in 1:length(master_df_list)) {
    master_df <- na.omit(master_df_list[[i]])
    print(outfn)
    outpath_curr <- as.character(outpath[i])
    print(outpath_curr)
    #outfn <- as.character(outfn[i])
    
    if((length(master_df) > 1) & (!is.na(outpath_curr) & (!is.na(outfn)))) {
      if(('p.value' %fin% colnames(master_df)) & (nrow(master_df) > 0)) {
        master_df <- process_raw_output_df(master_df = master_df, outpath = outpath_curr, outfn = outfn, 
                                           collinearity_diagn = collinearity_diagn, debug = debug,
                                           randomize = randomize, useNumFunctCopies = useNumFunctCopies)
        if(!useNumFunctCopies) {
          master_df_mut <- master_df[grepl("MutStat_i", master_df$term),]
          master_df_cna <- master_df[grepl("CNAStat_i", master_df$term),]
        } else {
          master_df_fnc <- master_df[master_df$term == "FNCStat_i",]
        }
        
        # Call the file to create output visualizations
        ct <- unlist(strsplit(outpath_curr, "/", fixed = T))[13]
        
        source(paste(source_path, "process_LM_output.R", sep = "")) 
      }
    }
    else {
      print(paste("Error: either master DF is empty or outfile name/ path is NA for cancer type",
                  names(master_df_list)[i]))
    }
  }
}


