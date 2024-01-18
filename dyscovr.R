############################################################
### DYSCOVR Linear Regression Model
### PUBLICATION INFORMATION
############################################################

#!/usr/bin/env Rscript

# Set the library path
library.path <- .libPaths()

library(parallel, lib.loc = library.path, quietly = T)
library(rlang, lib.loc = library.path, quietly = T)
library(dplyr, lib.loc = library.path, quietly = T)
library(broom, lib.loc = library.path, quietly = T)
library(data.table, lib.loc = library.path, quietly = T)
library(speedglm, lib.loc = library.path, quietly = T)
library(argparse, lib.loc = library.path, quietly = T)
library(stringr, lib.loc = library.path, quietly = T)
library(dqrng, lib.loc = library.path, quietly = T)
library(foreach, lib.loc = library.path, quietly = T)
library(snow, lib.loc = library.path, quietly = T)
library(doParallel, lib.loc = library.path, quietly = T)
library(fastmatch, lib.loc = library.path, quietly = T)
library(Hmisc, lib.loc = library.path, quietly = T)
library(caret, lib.loc = library.path, quietly = T)
library("gplots", lib.loc = library.path, quietly = TRUE)
library(qvalue, lib.loc = library.path, quietly = TRUE)
library(ggplot2, lib.loc = library.path, quietly = TRUE)
library(olsrr, lib.loc = library.path, quietly = TRUE)
library("RColorBrewer", lib.loc = library.path, quietly = TRUE)
library(glmnet, lib.loc = library.path)
library(gglasso, lib.loc = library.path)

# For Bayesian LASSO
library(statmod, lib.loc = library.path)
library(mvtnorm, lib.loc = library.path)
library(MCMCpack, lib.loc = library.path)
library(SuppDists, lib.loc = library.path)
library(mvnfast, lib.loc = library.path)
library(MBSGS, lib.loc = library.path)
library(MASS, lib.loc = library.path)
library(mgcv, lib.loc = library.path)
library(mnormt, lib.loc = library.path)
library(truncnorm, lib.loc = library.path)

# Source other files needed
SOURCE_PATH <- "/Genomics/argo/users/scamilli/Dyscovr"
source(paste(SOURCE_PATH, "general_important_functions.R", sep = "/"))
source(paste(SOURCE_PATH, "dyscovr_helper_functions.R", sep = "/"))
source(paste(SOURCE_PATH, "perform_collinearity_correction.R", sep = "/"))
source(paste(SOURCE_PATH, "run_regularized_models.R", sep = "/"))
source(paste(SOURCE_PATH, "run_bayesian_lasso.R", sep = "/"))


############################################################
# SET UP PARSER ARGUMENTS
############################################################
# Create parser object
parser <- ArgumentParser()

# Specify desired options for the parser

# Data set we are running on
parser$add_argument("--dataset", default = "TCGA", type = "character", 
                    help = "TCGA or METABRIC. Default is TCGA.")

# The file path to the linear model input files from create_dyscovr_input_tables.R
parser$add_argument("--input_lm_filepath", 
                    default = "/Genomics/grid/users/scamilli/thesis_work/
                    run-model-R/input_files/PanCancer/lm_input_tables/top_drivers_0.05/",
                    type = "character", help = "The path to the input files 
                    created from create_dyscovr_input_tables.R for the given run")

# Cancer type we are looking at (or pan-cancer), and patients we are looking at
parser$add_argument("--cancerType", default = "PanCancer", type = "character",
                    help = "Type of cancer we're looking at, or PanCancer. 
                    Default is PanCancer.")
parser$add_argument("--specificTypes", default = "ALL", type = "character",
                    help = "If and only if the cancer type is PanCancer, denotes 
                    which specific cancer types we're including (if running 
                    per-cancer). If ALL, run pan-cancer across all cancers together.")
parser$add_argument("--clinical_df", default = "clinical.csv", type = "character",
                    help = "If and only if the cancer type is PanCancer and 
                    specific types is not ALL, a clinical data frame name is 
                    given to link cancer types to patient IDs.")
parser$add_argument("--patientsOfInterest", default = "", type = "character",
                    help = "Rather given a particular cancer type or types, 
                    further restricts the patient pool to particular patients 
                    (e.g. from a particular subtype). A filename with a vector 
                    of patient IDs (no header) is optionally provided here.")
parser$add_argument("--patientsOfInterestLabel", default = "", type = "character",
                    help = "A character label for the specific patient population 
                    we are looking at here (for file naming purposes).")

# Whether or not we are running a test on one driver protein, and details of 
# the tester protein if so
parser$add_argument("--test", default = "F", type = "character",
                    help = "T/F for whether or not we are running a test with 
                    one driver protein, vs. many. Default is F.")
parser$add_argument("--tester_name", default = "P53", type = "character",
                    help = "Name of driver protein if doing a 1-protein test. 
                    Defaults to TP53.")
parser$add_argument("--tester_uniprot_id", default = "P04637", type = "character",
                    help = "Uniprot ID of driver protein if doing a 1-protein test. 
                    Defaults to TP53 (P04637).")
parser$add_argument("--tester_ensg_id", default = "ENSG00000141510", 
                    type = "character", help = "ENSG ID of driver protein if 
                    doing a 1-protein test. Defaults to TP53 (ENSG00000141510).")
parser$add_argument("--run_query_genes_jointly", default = "T", type = "character",
                    help = "If T, will run all the query driver proteins together 
                    in the same joint model, rather than running one model per 
                    query-target pair.")

# Randomization
parser$add_argument("--randomize", default = "F", type = "character", 
                    help = "T/F for whether or not we perform a randomization 
                    test. Default is F.")

# Additional input files needed. The directory for the LM input files created 
# from create_dyscovr_input_tables.R can be inferred from the other provided arguments.
parser$add_argument("--expression_df", 
                    default = "expression_quantile_norm_IntersectPatients.csv", 
                    type = "character", help = "The name of the expression data 
                    frame input. [default %(default)s]")
parser$add_argument("--cna_df", default = "CNA_DF_CancerOnly_IntersectPatients.csv", 
                    type = "character", help = "The name of the CNA data frame 
                    input; only needed if CNA values need adjusting. 
                    [default %(default)s]")
parser$add_argument("--patient_df", default = "combined_patient_sample_cibersort_
                    total_frac_washu_IntersectPatients.csv", type = "character",
                    help = "The name of the patient sample data frame input. 
                    [default %(default)s]")

# The decision for if/ how to bucket CNAs and methylation
parser$add_argument("--cna_bucketing", default = "bucket_inclAmp", 
                    type = "character", help = "If/ the type of bucketing we use 
                    for CNA values. Default is bucket_inclAmp, but rawCNA, 
                    bucket_exclAmp, bucket_justAmp, and bucket_justDel are also 
                    options.")
parser$add_argument("--meth_bucketing", default = T, type = "character",
                    help = "If/ the type of bucketing we use for methylation 
                    values. Default is F.")

# Label the type of methylation data we are using (Beta or M)
parser$add_argument("--meth_type", default = "M", type = "character",
                    help = "The type methylation values we are using, either M 
                    or Beta. Default is M.")

# Whether or not we are converting absolute CNA values to relative ones
parser$add_argument("--adjust_cna_rel", default = "F", type = "character",
                    help = "Whether or not we are adjusting the CNA values to be 
                    log2(CN+1 / N+1) rather than log2(CN+1). Default is F.")

# What type of test we are running, so we can write the results files to an 
# appropriate directory. Only necessary if --test is F.
parser$add_argument("--run_name", default = "cancerRelated", type = "character",
                    help = "Provide a name for the run, if not testing on just 
                    one driver protein, to write output files to appropriate 
                    directory. [default %(default)s]") 

# A name for the group of targets being tested, for labeling the output file.
parser$add_argument("--target_df", default = "allgene_targets.csv", 
                    type = "character", help = "Provide the name of the list of 
                    target genes; use to subset the given input DFs. Default is
                    allgene_targets.")
parser$add_argument("--targets_name", default = "allGenes", type = "character",
                    help = "Provide a name for the group of targets being tested, 
                    in order to properly name the output file. 
                    [default %(default)s]")

# Whether or not we are genotypic PCs as covariates in the model; for formula 
# creation
parser$add_argument("--num_pcs", default = 3, type = "integer", help = "A value 
                    between 0 and 3 indicating the number of principal components 
                    we are using as covariates in our model. Default is 3.")

# Add a flag for running detailed collinearity diagnostics/ including residuals; 
# can be useful but adds to run time.
parser$add_argument("--collinearity_diagn", default = "F", type = "character",
                    help = "A T/ F value indicating whether or not we want 
                    detailed collinearity diagnostics. Can be useful but adds to 
                    runtime; not suggested for large runs. Default is F.")
parser$add_argument("--inclResiduals", default = "F", type = "character",
                    help = "A T/F value indicating whether we want to include 
                    residuals in the output DF. Default is F.")
parser$add_argument("--correctPerDriver", default = "T", type = "character",
                    help = "A T/F value indicating whether we want to do q-value 
                    correction per-driver or across drivers. 
                    Default is T (per driver).")

# Similarly, a flag to calculate a Spearman correlation coefficient matrix prior to 
# running regression and either remove variables that are highly correlated, or add 
# an interaction term
parser$add_argument("--collinearity_corr_method", default = "eliminate_vif", 
                    type = "character", help = "Option to specify a method of 
                    dealing with highly correlated variables. Either 'eliminate' 
                    or 'interaction_term' are acceptable inputs. Option to add 
                    '_vif' to each of these to use the VIF score rather than 
                    Spearman correlation (better, but more computationally 
                    intensive). Default is 'eliminate_vif'.")

# Add a flag for the type of model to use (e.g. a mixed effects model, rather 
# than a simple linear or regularized linear)
parser$add_argument("--model_type", default = "linear", type = "character",
                    help = "The type of model being used. Use 'mixed' to switch 
                    from linear to a mixed effects model.")

# Add a flag for adding a regularization method
parser$add_argument("--regularization", default = "None", type = "character",
                    help = "The name of the regularization method being used. 
                    Currently only implemented for L1 and L2. Default is None.")
# If we are using a regularized model (L1 or L2), select the method to evaluate 
# significance.
parser$add_argument("--signif_eval_type", default = "randomization_perSamp", 
                    type = "character", help = "If regularization is not None, 
                    uses this flag to determine method of evaluating significance. 
                    Takes either randomization_perTarg, randomization_perSamp, 
                    randomization_predictors, subsampling, or selectiveInference. 
                    Default is randomization_predictors.")
# If we are using a regularized model (L1 or L2), and a randomization method for 
# evaluating significance, provide the number of randomizations. Default is 100.
parser$add_argument("--num_randomizations", default = 100, type = "integer",
                    help = "If regularization is not None, and signif_eval_type 
                    is some form of randomization, the number of randomizations 
                    is provided here.")
parser$add_argument("--fixed_lambda_for_rand", default = "F", type = "character",
                    help = "If regularization is not None, and signif_eval_type 
                    is some form of randomization, this T/F value indicates 
                    whether we use a cross validation to determine the lambda 
                    for each randomization or use a fixed lambda value from the 
                    real LASSO.")

# Add a flag for debugging
parser$add_argument("--debug", default = "F", type = "character",
                    help = "A T/ F value indicating whether or not we want 
                    detailed printing for debugging purposes. Default is F.")

# Add a flag for including only select covariates in the model, given as a 
# character list with each item separated by a semicolon
parser$add_argument("--select_args", default = "ALL", type = "character",
                    help = "Option to provide a semicolon-separated character 
                    string listing all covariates to include, e.g. MutStat_i, 
                    CNAStat_i, etc. If this option is chosen, dependent variable 
                    must also be listed. [default %(default)s]")
parser$add_argument("--select_args_label", default = "", type = "character",
                    help = "If we are providing a character string for 
                    --select_args, this is an option to provide a label for what 
                    arguments we are selecting for labeling the output files and 
                    visualizations.")
parser$add_argument("--select_drivers", default = "ALL", type ="character",
                    help = "Option to provide a semicolon-separated character 
                    string listing all driver Uniprot IDs to include, if LM input 
                    tables contain more driver covariates than is desired.")
parser$add_argument("--path_to_driver_lists", default = "", type ="character",
                    help = "Option to provide a path to files with drivers for 
                    individual cancer types. Will only be checked for a path is 
                    specificTypes is not ALL.")
parser$add_argument("--driver_file_prefix", default = "", type = "character",
                    help = "If path_to_driver_lists is not empty, this provides 
                    a prefix for the driver files we want to access for 
                    each cancer.")


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
runQueriesJointly <- str2bool(args$run_query_genes_jointly)
collinearity_diagn <- str2bool(args$collinearity_diagn)
inclResiduals <- str2bool(args$inclResiduals)
fixed_lambda_for_rand <- str2bool(args$fixed_lambda_for_rand)
correct_per_gene <- str2bool(args$correctPerDriver)
adjust_cna_rel <- str2bool(args$adjust_cna_rel)


############################################################
#### MAIN DYSCOVR LINEAR MODEL FUNCTION
############################################################
#' MAIN FUNCTION: This function will run the linear model, taking in the data 
#' frames we've created in create_dyscovr_input_tables.R.
#' @param list_of_input_dfs a list of input data frames from 
#' create_dyscovr_input_tables.R, where the list names are the Swissprot IDs of the 
#' associated target gene
#' @param expression_df table of expression (rows are genes, columns are patients, 
#' entries are expression values)
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. 
#' Possible values are "bucket_inclAmp", "bucket_exclAmp", "bucket_justAmp", 
#' "bucket_justDel" and "rawCNA"
#' @param meth_bucketing a T/F value indicating whether or not we are bucketing
#' methylation values
#' @param debug a T/F value indicating if we are in debug mode and should 
#' use additional prints
#' @param collinearity_diagn a T/F value indicating if we are running 
#' collinearity diagnostics
#' @param collinearity_corr_method either 'None', or a method to deal with 
#' correlated variables; either 'eliminate' or 'interaction_term', with an 
#' optional '_vif' on the end to indicate we are using VIF rather than Spearman
#' @param model_type the type of model being used (either linear, the default, 
#' or mixed)
#' @param regularization a string with the type of regularization to be used; 
#' currently only implemented for "L2" or "None"
#' @param signif_eval_type if regularization is not "None", uses the given type 
#' of method to evaluate significance (either some form of "randomization", 
#' "subsampling", or "selectiveInference")
#' @param num_randomizations the number of randomizations, if our signif_eval_type 
#' is some form of randomization
#' @param fixed_lambda_for_rand a T/F value indicating whether or not we are 
#' using a fixed value of lambda for the randomizations (e.g. from the "real" 
#' lasso) or if we are using cross-validation to determine a lambda value
#' @param run_query_genes_jointly a T/ F value indicating whether or not we 
#' want to run all query driver proteins jointly in the same model, or separately
#' (one model per query-target pair)
#' @param dataset either 'TCGA' or 'METABRIC' to denote what columns we are 
#' referencing/ specific data types we are using
#' @param inclResiduals a T/F value indicating whether or not we want to include 
#' residuals in our output DF, to evaluate model fit
#' @param outfn core of an output filename, for use with 
#' run_collinearity_diagnostics.R
#' @param outpath an output filepath, for use with run_collinearity_diagnostics.R
run_linear_model <- function(list_of_input_dfs, expression_df, cna_bucketing, 
                             meth_bucketing, debug, collinearity_diagn, 
                             collinearity_corr_method, model_type, regularization, 
                             signif_eval_type, num_randomizations, 
                             fixed_lambda_for_rand, run_query_genes_jointly, 
                             dataset, inclResiduals, outfn, outpath) {
  
  # Run through all targets in the downstream target data frame, running a 
  # regression model for each. These results will be bound together into a 
  # single output "master" DF
  master_df <- NA
  
  # Create a file to log our optimal lambda values
  #file.create("optimal_lambdas_real.txt")
  
  # If we are using regularization techniques, create shuffled input DFs if 
  # using an empirical approach to p-value generation
  shuffled_input_dfs <- NA
  if((model_type == "linear") & (regularization %fin% 
                                 c("L1", "L2", "lasso", "ridge"))) {
    shuffled_input_dfs <- get_shuffled_input_dfs(signif_eval_type, lm_input_dfs[1], 
                                                 expression_df, num_randomizations)
  }

  # If we are resolving multicollinearity, create a file where we will store
  # variables that were eliminated from each model
  #removed_variables_outfn <- paste(
  #  outpath, paste0(outfn, "_eliminated_variables.csv"), sep = "/")
  #fwrite(data.table("Variables.Removed"), removed_variables_outfn, col.names = F)
  
  # Use multithreading to run these jobs in parallel
  num_groups <- 8
  
  # Various types, including 'fork' and 'socket'. Fork is faster and runs on a 
  # Linux system (but not Windows). Default is SOCK.
  parallel_cluster <- makeCluster(num_groups, type = "SOCK", methods = F, outfile = "")
  #parallel_cluster <- makeCluster(num_groups, type = "FORK", methods = F, outfile = "")
  
  # Tell R that we want to use these processes to do calculations
  setDefaultCluster(parallel_cluster)
  doParallel::registerDoParallel(parallel_cluster)
  
  master_df <- foreach(
    i = 1:length(list_of_input_dfs),
    # Export needed data and functions
    .export = c("list_of_input_dfs", "cna_bucketing", "meth_bucketing", "debug",  
                "dataset", "inclResiduals", "model_type", "collinearity_diagn", 
                "collinearity_corr_method", "regularization", "signif_eval_type", 
                "num_randomizations", "fixed_lambda_for_rand", "shuffled_input_dfs", 
                "run_query_genes_jointly", "outfn", "outpath", 
		"%fin%", "fmatch", "fix_lm_input_table", "correct_collinearity", "correct_collinearity_spearman",
		"get_spearman_combos", "calculate_vif", "combine_collinearity_diagnostics",
  		"construct_formula", "run_regularization_model", "get_groups_for_variables",
		"multiResultClass", "combine_variables_rm_tables", "combine_summary_tables", "combine_res_obj",
		"eval_signif_selectiveInference", "eval_signif_randomization",
		"eval_signif_subsampling", "runBGL", "runBGLSS", "run_regularization_output_in_lm",
		"get_per_tg_randomized_beta_df", "get_per_samp_randomized_beta_df",
		"get_predictor_randomized_beta_df", "get_prop_pval", "get_fitted_norm_pval",
		"selectGroups", "getConfidence", "fitlambdaEM", "calculate_sparsity", "convert_grouptype"),
    .packages = c("dplyr", "data.table", "dqrng", "parallel", "fastmatch", "stringr", "rlang",
		  "broom", "speedglm", "Hmisc", "caret", "gplots", "qvalue", "ggplot2", "olsrr",
		  "RColorBrewer", "glmnet", "gglasso", "statmod", "mvtnorm", "MCMCpack", "SuppDists",
		  "mvnfast", "MBSGS", "MASS", "mgcv", "mnormt", "truncnorm"), 
    .combine = function(x,y) combine_res_obj(x, y)) %dopar% {
      
      # Source the necessary files for the worker node
      #SOURCE_PATH <- "/Genomics/argo/users/scamilli/Dyscovr"
      #source(paste(SOURCE_PATH, "general_important_functions.R", sep = "/"))
      #source(paste(SOURCE_PATH, "dyscovr_helper_functions.R", sep = "/"))
      #source(paste(SOURCE_PATH, "perform_collinearity_correction.R", sep = "/"))
      #if (regularization != "None") {
      #  source(paste(SOURCE_PATH, "run_regularized_models.R", sep = "/"))
      #  source(paste(SOURCE_PATH, "run_bayesian_lasso.R", sep = "/"))
      #}
      
      lm_input_table = list_of_input_dfs[[i]]
      if(debug) {print(head(lm_input_table))}
      
      # Make sure any -Inf or Inf values are given NA and then omitted
      # If a column is >50% Inf or -Inf, remove that column first
      if(dataset %fin% c("METABRIC")) {
        lm_input_table <- fix_lm_input_table(lm_input_table)
      }
      
      targ = names(list_of_input_dfs)[i]
      if(debug) {print(targ)}
      
      if(length(lm_input_table) == 0) {
	print(paste("Null input data table for target:", 
                    paste0(targ, ". Returning NA.")))
      }
      if((lm_input_table[, .N] == 0) | (ncol(lm_input_table) == 0)) {
        print(paste("Null input data table for target:", 
                    paste0(targ, ". Returning NA.")))
        lm_input_table <- NA
      }                  
      # Option to correct for collinear variables using Spearman method
      if(collinearity_corr_method !=  "None") {
        if(!grepl("vif", collinearity_corr_method)) {
          collinearity_res <- correct_collinearity(lm_input_table, 
                                                 method = collinearity_corr_method,
						 target = targ)
                                                 #outfile = removed_variables_outfn)
 	  lm_input_table <- collinearity_res[[1]]
          variables_rm <- collinearity_res[[2]]

          if(debug) {
            print(head(lm_input_table))
          }
        }
      }
      
      swissprot_ids <- unlist(lapply(colnames(lm_input_table)[
        grepl("MutStat_i", colnames(lm_input_table))], function(x)
          unlist(strsplit(x, "_", fixed = T))[1]))
      swissprot_ids <- unique(swissprot_ids[swissprot_ids != ""])
      input_drivers <- data.table("swissprot" = swissprot_ids)
      print(head(input_drivers))      

      # Generate formula dynamically
      formula <- tryCatch({
        construct_formula(lm_input_table, input_drivers, cna_bucketing, 
				meth_bucketing, run_query_genes_jointly) 
        }, error = function(cond) {
          print(cond) 
          return(NA)
      })
      print(formula)
      
      # If we are using a VIF method to correct for collinearity, do this now 
      # after we have created the formula for the model fit
      if(grepl("vif", collinearity_corr_method)) {
        tryCatch({
          collinearity_res <- correct_collinearity(lm_input_table, formula,
                                                 method = collinearity_corr_method,
						 target = targ)
                                                 #outfile = removed_variables_outfn)
          lm_input_table <- collinearity_res[[1]]
          variables_rm <- collinearity_res[[2]]

          }, error = function(cond) {
            print(cond)
            if(grepl("str2lang", cond)) {
              print(paste("FORMULA:", formula))
              print("Head of LM input table:")
              print(head(lm_input_table))
            }
            print("Unable to correct collinearity. Returning NA.")
            lm_input_table <- NA
          })
        if(debug) {
          print(head(lm_input_table))
        }
        
        # Adjust the formula to reflect the new input table 
        if(!is.na(lm_input_table)) {
          formula <- tryCatch({
            construct_formula(lm_input_table, input_drivers, cna_bucketing, 
				meth_bucketing, run_query_genes_jointly)
            }, error = function(cond) {
              print(cond)
              return(NA)
          })
        }
        
        if(debug) {print(paste("Formula:", formula))}
        
        # Ensure the formula looks correct
        if(formula == "ExpStat_k ~ ") {
          print("Truncated formula. Something went wrong.")
          print(head(colnames(lm_input_table)))
          print(swissprot_ids)
          return(NA)
        }
        
        summary_table <- NA
        lm_fit <- NA
        
        # Run regression based on model type and desired regularization
        if ((model_type == "linear") & (regularization == "None")) {
          print("Now running linear regression...")
          
          if(!is.na(lm_input_table)) {
            lm_fit <- tryCatch({
              speedglm::speedlm(formula = formula, data = lm_input_table)  
            }, error = function(cond) {
              message(paste("There was a problem with this linear model run:", 
                            cond))
              print(formula)
              print("Offending LM Input Table:")
              print(head(lm_input_table))
              print(dim(lm_input_table))
              print(any(is.na(lm_input_table)))
              print(targ)
              return(NA)
            })
          }
          if(!is.na(lm_fit)) {
            summary_table <- tidy(lm_fit)
            summary_table <- as.data.table(
              summary_table[, colSums(is.na(summary_table)) < nrow(summary_table)])
            rm(lm_fit)
          } 
        } else if ((model_type == "linear") & (regularization == "L2" | 
                                               regularization == "ridge" | 
                                               regularization == "L1" | 
                                               regularization == "lasso" |
                                               regularization == "bayesian.bgl" | 
                                               regularization == "bayesian.bglss")) {
          if(!is.na(lm_input_table)) {
            summary_table <- run_regularization_model(
              formula = formula, lm_input_table = lm_input_table, 
              type = regularization, debug = debug, meth_bucketing = meth_bucketing, 
              cna_bucketing = cna_bucketing, signif_eval_type = signif_eval_type, 
              shuffled_input_dfs = shuffled_input_dfs, ensg = targ_ensg,
              num_randomizations = num_randomizations,
              fixed_lambda_for_rand = fixed_lambda_for_rand)
           }
        } else if (model_type == "mixed") {
          
          # Our random effects variable is our genotype PCs; if we are not 
          # including them, run the simple linear model instead
          if(!is.na(lm_input_table)) {
            if (num_pcs == 0) {
              print("No genotype PCs are being used in regression. 
                    Random effects variables are currently only implemented for 
                    genotype PCs. Running a linear regression.")
              tryCatch({
                lm_fit <- speedglm::speedlm(formula = formula, 
                                            data = lm_input_table) 
                }, error = function(cond) {
                  message("There was a problem with this linear model run:")
                  message(cond)
                })
              } else {
                # Build a linear mixed model
                tryCatch({
                  lm_fit <-  lmer(formula = formula, data = lm_input_table, 
                                  REML = T)
                  }, error = function(cond) {
                    message("There was a problem with this linear mixed model run:")
                    message(cond) 
                  })
              }
            }
          if(!is.na(lm_fit)) {
            summary_table <- as.data.table(tidy(lm_fit))
            rm(lm_fit)
          }
        } else {
          if(is.na(formula)) {
            print("Encountered an error. Formula is NA.")
            summary_table <- NA
          } else {
            print(paste("Model type", paste(model_type, "does not match any of 
                                            the implemented options. Please 
                                            provide one of these options.")))
            summary_table <- NA
          }
        }
        if(inclResiduals & (!(is.na(summary_table)))) {
          residuals <- as.numeric(lm_input_table$ExpStat_k) - 
            as.numeric(unlist(predict.speedlm(lm_fit, newdata = lm_input_table)))
          summary_table$Residual <- paste(residuals, collapse = ";")
        }
        
        # Add a column for the driver protein and target gene to ID them
        if(!is.na(summary_table)) {
          summary_table$T_k <- paste(targ, collapse = ";")
          # Remove rows in which the target and the driver are the same, prior
          # to MHT correction
          ri <- unlist(lapply(summary_table$term, function(x)
            unlist(strsplit(x, "_", fixed = T))[1]))
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
      #  print("Running Collinearity Diagnostics")
      #  fit_lm <- lm(formula = formula, data = lm_input_table)
      #  source(paste(getwd(), "run_collinearity_diagnostics.R", sep = "/"), 
      #         local = T)
	 print("Collinearity diagnostics currently in repair.")
      }
      
      # Make sure that the summary table for this gene is returned
      if(debug) {
        print("Summary Table")
        print(head(summary_table))
      }
      
      # Try to clear up some memory
      gc()
                         
      # Return summary table
      final_results <- multiResultClass(summary_table, variables_rm)
      return(final_results)
  }
  
  # tell R that we don't need the processes anymore
  stopCluster(parallel_cluster)
  
  # Add target names to the DF containing removed variables
  #variables_rm_df <- fread(removed_variables_outfn, header = T)
  #variables_rm_df <- cbind(data.table("Target" = names(list_of_input_dfs), 
   #                                   variables_rm_df))
  #fwrite(variables_rm_df, removed_variables_outfn)
  
  if(debug) {
    print("Dyscovr Master DF")
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
    cancer_types <- unlist(strsplit(args$specificTypes, ";", fixed = T))
    outpath <- lapply(cancer_types, function(ct) 
      create_file_outpath(cancerType = args$cancerType, specificType = ct, 
                          test = test, tester_name = args$tester_name, 
                          run_name = args$run_name, dataset = args$dataset,
			  outpath = SOURCE_PATH))
    names(outpath) <- cancer_types
  } else {
    outpath <- create_file_outpath(cancerType = args$cancerType, specificType = "", 
                                   test = test, tester_name = args$tester_name, 
                                   run_name = args$run_name, dataset = args$dataset,
				   outpath = SOURCE_PATH)
  }
} else {
  outpath <- create_file_outpath(cancerType = args$cancerType, specificType = "", 
                                 test = test, tester_name = args$tester_name, 
                                 run_name = args$run_name, dataset = args$dataset,
				 outpath = SOURCE_PATH)
}

if(debug) {print(paste("Outpath:", outpath))}

# Create an appropriate output file name
outfn <- create_output_filename(test = test, tester_name = args$tester_name, 
                                run_name = args$run_name, 
                                targets_name = args$targets_name, 
                                expression_df_name = args$expression_df,
                                cna_bucketing = args$cna_bucketing, 
                                meth_bucketing = meth_bucketing, 
                                meth_type = args$meth_type, 
                                patient_df_name = args$patient_df[1], 
                                num_pcs = args$num_pcs, 
                                randomize = randomize, 
                                covs_to_incl_label = args$select_args_label,
                                patients_to_incl_label = args$patientsOfInterestLabel, 
                                removeCis = F, model_type = args$model_type, 
                                removeMetastatic = T, 
                                regularization = args$regularization, 
                                signif_eval_type = args$signif_eval_type)[1]

if(debug) {print(paste("Outfile Name:", outfn))}


############################################################
# IMPORT COVARIATES
############################################################
# Convert covariates to include from a semicolon-separated character to a vector
covs_to_incl <- c("ALL")

tryCatch({
  covs_to_incl <- unlist(strsplit(args$select_args, ";", fixed = T))
}, error = function(cond) {
  print(cond)
  print("Invalid input covariate list. Make sure covariates are semi-colon 
        separated. Running with all covariates.")
})


############################################################
# SET PATHS AND IMPORT TARGET GENE DF
############################################################
# If needed, subset this list based on the given targets in the target DF
INPUT_PATH <- "/Genomics/argo/users/scamilli/Dyscovr/input_files/"
MAIN_PATH <- NA
TARG_PATH <- NA

if(args$cancerType == "BRCA") {
  TARG_PATH <- paste0(INPUT_PATH, "BRCA/target_lists/")
  
  if(args$dataset == "METABRIC") {
    MAIN_PATH <- paste0(INPUT_PATH, "METABRIC/")
  } else {
    MAIN_PATH <- paste0(INPUT_PATH, "BRCA/") 
  }
} else {
  TARG_PATH <- paste0(INPUT_PATH, "PanCancer/target_lists/")
  MAIN_PATH <- paste0(INPUT_PATH, "PanCancer/")
}

target_df <- fread(paste0(TARG_PATH, args$target_df), header = T)

target_uniprot_ids <- unique(unlist(lapply(target_df$swissprot, function(x)
  unlist(strsplit(x, ";", fixed = T)))))
target_uniprot_ids <- target_uniprot_ids[!is.na(target_uniprot_ids)]


############################################################
# IF RUNNING PER-CANCER WITH SPECIFIC CANCER TYPES, GET THE PATIENTS 
# OF THOSE CANCER TYPES 
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
  clinical_df <- fread(paste0(MAIN_PATH, paste0("clinical/", clinical_df)), 
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
      specificTypes = args$specificTypes, clinical_df = args$clinical_df,
      main_path = main_path)
    if(debug) {print(patient_cancer_mapping)}
  }
}

# Read if the patient IDs of interest, if provided
patients_of_interest <- ""
if(args$patientsOfInterest != "") {
  patients_of_interest <- as.character(unlist(
    fread(paste0(MAIN_PATH, paste0("patient/groups_of_interest/", 
                                   args$patientsOfInterest)), header = F)[,1]))
}

if(debug) {
  print("Patients of Interest")
  print(head(patients_of_interest))
}


############################################################
# IMPORT LIST OF PRE-CREATED LM INPUT TARGET FILES
############################################################
lm_input_fns <- NA
driver_input_fn <- NA

# Input LM filepath
input_lm_filepath <- args$input_lm_filepath

# Get the index of where the target gene name will be found in each file name
#gene_name_index <- length(unlist(strsplit(args$run_name, "_", fixed = T))) + 3
gene_name_index <- 3

# If we need to adjust CNA values, import the DF
cna_df <- NA
if(adjust_cna_rel) {
  cna_df <- fread(paste0(MAIN_PATH, paste0("cna/", args$cna_df)), header = T)
  colnames(cna_df)[1] <- "ensg_id"
} 

# For multiple individual cancer types, get the file names in a list
if((args$cancerType == "PanCancer") & (args$specificTypes != "ALL")) {
  spl_fp <- unlist(strsplit(input_lm_filepath, "/", fixed = T))
  index_pc <- which(spl_fp == "PanCancer")
  pt1 <- paste(spl_fp[1:index_pc], collapse = "/")
  pt2 <- paste(spl_fp[(index_pc+1):length(spl_fp)], collapse = "/")
  
  input_lm_filepath <- c()
  lm_input_fns <- list()
  driver_input_fn <- list()
  
  for (ct in names(patient_cancer_mapping)) {
    new_fp <- paste(pt1, paste(ct, pt2, sep = "/"), sep = "/")
    input_lm_filepath <- c(input_lm_filepath, new_fp)
    new_fns <- list.files(new_fp, pattern = "lm_input_", recursive = F)
    if(!(args$patientsOfInterestLabel %fin% c("", "male", "female"))) {
      gene_name_index <- gene_name_index + 1
      new_fns <- intersect(new_fns, 
                           list.files(new_fp, 
                                      pattern = args$patientsOfInterestLabel, 
                                      recursive = F))
    }
    lm_input_fns <- append(lm_input_fns, list(new_fns))
    new_driver_fn <- list.files(new_fp, pattern = "driver_input_", 
                                recursive = F)
    driver_input_fn <- append(driver_input_fn, new_driver_fn)
  }
  names(lm_input_fns) <- names(patient_cancer_mapping)
  if(debug) {print(lm_input_fns)}
  
# If for one cancer type or PanCancer, get the file names
} else {   
  lm_input_fns <- unique(list.files(input_lm_filepath, pattern = "lm_input_", 
                                    recursive = F))
  driver_input_fn <- list.files(input_lm_filepath, pattern = "driver_input_", 
                                recursive = F)
  
  if(!(args$patientsOfInterestLabel %fin% c("", "male", "female"))) {
    gene_name_index <- gene_name_index + 1
    lm_input_fns <- intersect(lm_input_fns, 
                              list.files(input_lm_filepath, 
                                         pattern = args$patientsOfInterestLabel, 
                                         recursive = F))
  } 
  # Check if there are sub-directories in this path that might contain the file
  dirs_in_path <- list.dirs(input_lm_filepath, recursive = F)
  if(T %fin% unlist(lapply(dirs_in_path, function(d) grepl("part_", d)))) {
    dirs_part <- dirs_in_path[grepl("part_", dirs_in_path)]
    gene_name_index <- gene_name_index + 1
    for (d in dirs_part) {
      files_d <- list.files(paste0(d, "/"), pattern = "lm_input_", recursive = F)
      if(!(args$patientsOfInterestLabel %fin% c("", "male", "female"))) {
        files_d <- intersect(files_d, 
                             list.files(paste0(d, "/"), 
                                        pattern = args$patientsOfInterestLabel,
                                        recursive = F))
      } 
      part <- unlist(strsplit(d, "/", fixed = T))
      the_part <- part[length(part)]
      files_d <- paste(the_part, files_d, sep = "/")
      lm_input_fns <- unique(c(lm_input_fns, files_d))
    }
  } 
}

# Get the names of all the genes we have LM input files for
gene_names <- NA
if((args$cancerType == "PanCancer") & (args$specificTypes != "ALL")) {
  gene_names <- lapply(lm_input_fns, function(ct_fn) {
    print(paste("# of Input LM DFs", length(ct_fn)))

    gns <- unlist(lapply(ct_fn, function(fn) 
      unlist(strsplit(fn, "_", fixed = T))[gene_name_index]))
    if(debug) {print(head(gns))}
    return(gns)
  })
  names(gene_names) <- names(patient_cancer_mapping)
  
} else {
  gene_names <- unlist(lapply(lm_input_fns, function(fn) 
    unlist(strsplit(fn, "_", fixed = T))[gene_name_index]))
  if(debug) {print(head(gene_names))}
  
  print("# of Input LM DFs")
  print(length(lm_input_fns))
}

lm_input_fns_sub <- lm_input_fns
gene_names_sub <- gene_names

if((args$cancerType == "PanCancer") & (args$specificTypes != "ALL")) {
  for (i in 1:length(lm_input_fns)) {
    input_fns <- lm_input_fns[[i]]
    gns <- gene_names[[i]]
    
    to_keep <- which(gns %fin% target_uniprot_ids)
    lm_input_fns_sub[[i]] <- input_fns[to_keep]
    gene_names_sub[[i]] <- gns[to_keep]
    
    print("# of Input LM DFs, subsetted")
    print(length(lm_input_fns_sub[[i]]))
  }
} else {
  to_keep <- which(gene_names %fin% target_uniprot_ids)
  lm_input_fns_sub <- lm_input_fns[to_keep]
  gene_names_sub <- gene_names[to_keep]
  
  print("# of Input LM DFs, subsetted")
  print(length(lm_input_fns_sub))
}


# If the number of genes is greater than a particular threshold N (e.g. 8K genes), 
# split into multiple partial runs so that we don't exceed memory constraints
split_target_list_flag <- F
N <- 8000
lm_input_dfs <- NA

if((args$cancerType == "PanCancer") & (args$specificTypes != "ALL")) {
  split_target_list_flag <- c()
  lm_input_dfs <- list()
  for (i in 1:length(lm_input_fns_sub)) {
    ct_fns <- lm_input_fns_sub[[i]]
    ct_gens <- gene_names_sub[[i]]
    ct_filepath <- input_lm_filepath[i]
    ct_driver_input_fn <- driver_input_fn[[i]]
    input_dfs <- call_import_lm_file(
      ct_fns, ct_filepath, ct_driver_input_fn, randomize, ct_gens, 
      args$select_drivers, covs_to_incl, args$num_pcs, cna_df, 
      adjust_cna_rel, N, debug)
    if(length(ct_fns) > N) {split_target_list_flag <- c(split_target_list_flag, T)}
    else {split_target_list_flag <- c(split_target_list_flag, F)}
    lm_input_dfs <- append(lm_input_dfs, list(input_dfs))
  }
  names(lm_input_dfs) <- names(lm_input_fns_sub)
  
} else {
  lm_input_dfs <- call_import_lm_file(lm_input_fns_sub, input_lm_filepath, 
                                      driver_input_fn, randomize, gene_names_sub, 
                                      args$select_drivers, covs_to_incl,
                                      args$num_pcs, cna_df, adjust_cna_rel,
                                      N, debug)
  if(length(lm_input_fns_sub) > N) {split_target_list_flag <- T}
}


############################################################
#### CALL FUNCTION
############################################################
# Run gc to free up any loose memory
gc()

# Run the function itself

# Case 1: If we're running on just one cancer type or PanCancer
if(!((args$cancerType == "PanCancer") & (args$specificTypes != "ALL"))) {
  
  # If needed, subset these data tables by the given patient IDs of interest
  if(patients_of_interest != "") {
    lm_input_dfs <- lapply(lm_input_dfs, function(df) 
      subset_by_intersecting_ids(patients_of_interest, df, F, args$dataset))
  }
  
  master_res <- run_linear_model(list_of_input_dfs = lm_input_dfs, 
                                expression_df = expression_df, 
                                cna_bucketing = args$cna_bucketing,
                                meth_bucketing = meth_bucketing,
                                debug = debug,
                                collinearity_diagn = collinearity_diagn,
                                collinearity_corr_method = args$collinearity_corr_method,
                                regularization = args$regularization,
                                signif_eval_type = args$signif_eval_type,
                                num_randomizations = args$num_randomizations,
                                fixed_lambda_for_rand = fixed_lambda_for_rand,
                                model_type = args$model_type,
                                run_query_genes_jointly = runQueriesJointly,
                                dataset = args$dataset,
                                inclResiduals = inclResiduals,
                                outfn = outfn, outpath = outpath)
  master_df <- master_res$summary_table
  variable_rm_df <- master_res$variables_rm
  
  if(split_target_list_flag) {
    
    # Clear memory from the previous batch
    rm(lm_input_dfs)
    gc()
    
    # Run again on the next batch
    lm_input_dfs <- import_lm_input_fns(
      lm_input_fns_sub[(N+1):length(lm_input_fns_sub)], input_lm_filepath, 
      driver_input_fn, randomize, gene_names_sub[(N+1):length(lm_input_fns_sub)], 
      args$select_drivers, covs_to_incl, args$num_pcs, cna_df, adjust_cna_rel)
    
    # If needed, subset these data tables by the given patient IDs of interest
    if(patients_of_interest != "") {
      lm_input_dfs <- lapply(lm_input_dfs, function(df) 
        subset_by_intersecting_ids(patients_of_interest, df, F, args$dataset))
    }
    
    master_res_2 <- run_linear_model(list_of_input_dfs = lm_input_dfs, 
                                    expression_df = expression_df, 
                                    cna_bucketing = args$cna_bucketing,
                                    meth_bucketing = meth_bucketing,
                                    debug = debug,
                                    collinearity_diagn = collinearity_diagn,
                                    collinearity_corr_method = args$collinearity_corr_method,
                                    regularization = args$regularization,
                                    signif_eval_type = args$signif_eval_type,
                                    num_randomizations = args$num_randomizations,
                                    fixed_lambda_for_rand = fixed_lambda_for_rand,
                                    model_type = args$model_type,
                                    run_query_genes_jointly = runQueriesJointly,
                                    dataset = args$dataset,
                                    inclResiduals = inclResiduals,
                                    outfn = outfn, outpath = outpath)
    master_df_2 <- master_res_2$summary_table
    variable_rm_df_2 <- master_res_2$variables_rm

    master_df <- rbind(master_df, master_df_2)
    variable_rm_df <- rbind(variable_rm_df, variable_rm_df_2)
    master_res <- multiResultClass(master_df, variable_rm_df)
  }
  
  
# Case 2: If we're running on multiple cancer types separately
} else if ((args$cancerType == "PanCancer") & (args$specificTypes != "ALL")) {
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
        df_sub <- subset_by_intersecting_ids(patient_ids, df, F, args$dataset)
        if(!(nrow(df_sub) == 0)) {return(df_sub)}
        else {return(NA)}
      })
      lm_input_dfs_sub <- lm_input_dfs_sub[!is.na(lm_input_dfs_sub)]
      if(length(lm_input_dfs_sub) == 0) {return(NA)}
      
      # Subset input DFs to drivers that pass particular threshold
      if(args$path_to_driver_lists != "") {
        driver_filename <- paste0(args$path_to_driver_lists, 
                                  paste0(args$driver_file_prefix, 
                                         paste0(cancer_type, ".csv")))
        driver_file <- fread(driver_filename, header = T)
        
        select_driver_vect <- unique(unlist(driver_file$swissprot_ids))
        
        if(debug) {
          print("Selected Drivers:")
          print(select_driver_vect)
        }
        
        lm_input_dfs_sub <- lapply(lm_input_dfs_sub, function(df) {
          drivers_in_tab <- unique(unlist(lapply(colnames(df)[
            grepl("MutStat_i", colnames(df))], function(x)
              unlist(strsplit(x, "_", fixed = T))[1])))
          drivers_to_exclude <- setdiff(drivers_in_tab, select_driver_vect)
          drivers_to_exclude_covs <- unlist(lapply(drivers_to_exclude, function(d) {
            if(d != "") {
              return(colnames(df)[grepl(d, colnames(df), fixed = T)])
            } else {return(NA)}
          }))
          drivers_to_exclude_covs <- drivers_to_exclude_covs[!is.na(
            drivers_to_exclude_covs)]
          df_sub <- df[, !(colnames(df) %fin% drivers_to_exclude_covs), with = F]
          return(df_sub)
        })
      }
      lm_input_dfs_sub <- lm_input_dfs_sub[!is.na(lm_input_dfs_sub)]
      
      # Run the model with these subsetted files
      master_res <- run_linear_model(list_of_input_dfs = lm_input_dfs_sub, 
                                    expression_df = expression_df, 
                                    cna_bucketing = args$cna_bucketing,
                                    meth_bucketing = meth_bucketing,
                                    debug = debug,
                                    collinearity_diagn = collinearity_diagn,
                                    collinearity_corr_method = 
                                      args$collinearity_corr_method,
                                    regularization = args$regularization,
                                    signif_eval_type = args$signif_eval_type,
                                    num_randomizations = args$num_randomizations,
                                    fixed_lambda_for_rand = fixed_lambda_for_rand,
                                    model_type = args$model_type,
                                    run_query_genes_jointly = runQueriesJointly,
                                    dataset = args$dataset,
                                    inclResiduals = inclResiduals,
                                    outfn = outfn, outpath = outpath[i])
      return(master_res)
    } else {return(NA)}
  })
  
  if(T %fin% split_target_list_flag) {
    
    # Clear memory from the previous batch
    rm(lm_input_dfs)
    rm(lm_input_dfs_ct)
    rm(lm_input_dfs_sub)
    gc()
    
    # Run again on the next batch
    lm_input_dfs <- lapply(1:length(lm_input_fns_sub), function(i) {
      ct_fns <- lm_input_fns_sub[[i]]
      ct_gens <- gene_names_sub[[i]]
      ct_filepath <- input_lm_filepath[i]
      ct_driver_input_fn <- driver_input_fn[[i]]
      input_dfs <- import_lm_input_fns(ct_fns[(N+1):length(ct_fns)], 
                                       ct_filepath, ct_driver_input_fn, 
                                       randomize, ct_gens[(N+1):length(ct_fns)], 
                                       args$select_drivers, covs_to_incl, 
                                       args$num_pcs, cna_df, adjust_cna_rel)
      return(input_dfs)
    })
    names(lm_input_dfs) <- names(lm_input_fns_sub)
    
    master_df_list_2 <- lapply(1:length(outpath), function(i) {
      
      if(split_target_list_flag[i]) {
        cancer_type <- names(patient_cancer_mapping)[i]
        
        # Get the patients for this cancer type
        if(!is.na(cancer_type)) {
          patient_ids <- unique(unlist(patient_cancer_mapping[[i]]))

          if(patients_of_interest != "") {
            patient_ids <- patient_ids[patient_ids %fin% patients_of_interest]
          }
          
          # Subset input DFs to these patients of interest
          lm_input_dfs_ct <- lm_input_dfs[[i]]
          lm_input_dfs_sub <- lapply(lm_input_dfs_ct, function(df) {
            df_sub <- subset_by_intersecting_ids(patient_ids, df, F, args$dataset)
            if(!(nrow(df_sub) == 0)) {return(df_sub)}
            else {return(NA)}
          })
          lm_input_dfs_sub <- lm_input_dfs_sub[!is.na(lm_input_dfs_sub)]
          if(length(lm_input_dfs_sub) == 0) {return(NA)}
          
          # Subset input DFs to drivers that pass particular threshold
          if(args$path_to_driver_lists != "") {
            driver_filename <- paste0(args$path_to_driver_lists, 
                                      paste0(args$driver_file_prefix, 
                                             paste0(cancer_type, ".csv")))
            driver_file <- fread(driver_filename, header = T)
            select_driver_vect <- unique(unlist(driver_file$swissprot_ids))
            if(debug) {print(select_driver_vect)}
            
            lm_input_dfs_sub <- lapply(lm_input_dfs_sub, function(df) {
              drivers_in_tab <- unique(unlist(lapply(colnames(df)[
                grepl("MutStat_i", colnames(df))], function(x)
                  unlist(strsplit(x, "_", fixed = T))[1])))
              drivers_to_exclude <- setdiff(drivers_in_tab, select_driver_vect)
              drivers_to_exclude_covs <- unlist(lapply(
                drivers_to_exclude, function(d) {
                  if(d != "") {
                    return(colnames(df)[grepl(d, colnames(df), fixed = T)])
                  } else {return(NA)}
              }))
              drivers_to_exclude_covs <- drivers_to_exclude_covs[
                !is.na(drivers_to_exclude_covs)]
              df_sub <- df[, !(colnames(df) %fin% drivers_to_exclude_covs), 
                           with = F]
              return(df_sub)
            })
          }
          
          master_res_2 <- run_linear_model(list_of_input_dfs = lm_input_dfs_sub, 
                                          expression_df = expression_df, 
                                          cna_bucketing = args$cna_bucketing,
                                          meth_bucketing = meth_bucketing,
                                          debug = debug,
                                          collinearity_diagn = collinearity_diagn,
                                          collinearity_corr_method = 
                                            args$collinearity_corr_method,
                                          regularization = args$regularization,
                                          signif_eval_type = args$signif_eval_type,
                                          num_randomizations = args$num_randomizations,
                                          fixed_lambda_for_rand = fixed_lambda_for_rand,
                                          model_type = args$model_type,
                                          run_query_genes_jointly = runQueriesJointly,
                                          dataset = args$dataset,
                                          inclResiduals = inclResiduals,
                                          outfn = outfn, outpath = outpath[i])
          return(master_res_2)
        } else {return(NA)}
      } else {return(NA)}
    })
    
    master_df_list <- lapply(1:length(master_df_list), function(i) {
      master_res <- master_df_list[[i]]
      master_res_2 <- master_df_list_2[[i]]
      if(length(master_res_2) > 1) {
        master_df <- master_res$summary_table
        master_df_2 <- master_res_2$summary_table
        variables_rm_df <- master_res$variables_rm
        variables_rm_df_2 <- master_res_2$variables_rm
        return(multResultClass(rbind(master_df, master_df_2),
                               rbind(variables_rm_df, variables_rm_df_2)))
      } else {return(master_res)}
    })
  }
  
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
#' @param master_df a master data frame output result from Dyscovr
#' @param outpath the outpath for the master DF to be written to
#' @param outfn the generic output filename for the master DF
#' @param collinearity_diagn a T/F value indicating whether we are running
#' additional collinearity diagnostics
#' @param debug a T/F value indicating whether or not we are in debug mode
#' @param randomize a T/F value indicating whether or not this was a 
#' randomized run
process_raw_output_df <- function(master_df, outpath, outfn, collinearity_diagn, 
                                  debug, randomize) {
  # Order the file by p-value
  if(debug) {print(head(master_df))}
  try(master_df <- master_df[order(p.value)])

  # Write the results to the given file, making the directory if it does not 
  # already exist
  tryCatch({
    dir.create(outpath, showWarnings = F)
  }, error = function(cond) {
    print(cond)
    print("Invalid outpath. Cannot create specified directory.")
  })
  fwrite(master_df, paste(outpath, paste0(outfn, ".csv"), sep = "/"))
  
  # Limit the data frame to just the term of interest (typically either 
  # MutStat_i or CNAStat_i)
  master_df_mut <- master_df[grepl("MutStat_i", master_df$term),]
  #master_df_cna <- master_df[grepl("CNAStat_i", master_df$term),]
  
  # Write these to files as well 
  fwrite(master_df_mut, paste(outpath, paste0(outfn, "_MUT.csv"), sep = "/")) 
  #fwrite(master_df_cna, paste(outpath, paste0(outfn, "_CNA.csv"), sep = "/"))  
  
  # Combine the collinearity results and write to a file
  if(collinearity_diagn) {
    collinearity_df <- combine_collinearity_diagnostics(outpath, outfn, 
                                                        debug, randomize)
    # Adjust the output file name for collinearity results
    collin_fn <- str_replace(outfn, "output_results", "collinearity_results")
    if(debug) {print(paste("Collinearity Output FN:", collin_fn))}
    fwrite(collinearity_df, paste(outpath, paste0(collin_fn, ".csv"), sep = "/"))
  }
  return(master_df)
}


if((args$cancerType == "PanCancer") & (args$specificTypes != "ALL")) {
  for (i in 1:length(master_df_list)) {
    master_df <- na.omit(master_df_list[[i]]$summary_table)
    outpath_curr <- as.character(outpath[i])

    if((length(master_df) > 1) & (!is.na(outpath_curr) & (!is.na(outfn)))) {
      if(('p.value' %fin% colnames(master_df)) & (nrow(master_df) > 0)) {
        master_df <- process_raw_output_df(master_df = master_df, 
                                           outpath = outpath_curr, outfn = outfn, 
                                           collinearity_diagn = collinearity_diagn, 
                                           debug = debug, randomize = randomize)
        master_df_mut <- master_df[grepl("MutStat_i", master_df$term),]
        #master_df_cna <- master_df[grepl("CNAStat_i", master_df$term),]

        # Call the file to create output visualizations
        ct <- unlist(strsplit(outpath_curr, "/", fixed = T))[13]

        variables_rm_df <- na.omit(master_df_list[[i]]$variables_rm)
        removed_variables_outfn <- paste(outpath_curr, paste0(outfn, "_eliminated_variables.csv"), 
					 sep = "/")
	fwrite(variables_rm_df, paste0(removed_variables_outfn, ".csv"))
        
        source(paste(SOURCE_PATH, "process_dyscovr_output.R", sep = "/")) 
      }
    }
    else {
      print(paste("Error: either master DF is empty or outfile name/ path is NA 
                  for cancer type", names(master_df_list)[i]))
    }
  }
} else {
  outpath <- as.character(unlist(outpath[[1]]))
  outfn <- as.character(unlist(outfn[[1]]))
  master_df <- na.omit(master_res$summary_table)
  master_df <- process_raw_output_df(master_df = master_df, outpath = outpath, 
                                     outfn = outfn, 
                                     collinearity_diagn = collinearity_diagn, 
                                     debug = debug,
                                     randomize = randomize)
  master_df_mut <- master_df[grepl("MutStat_i", master_df$term),]
  #master_df_cna <- master_df[grepl("CNAStat_i", master_df$term),]

  variables_rm_df <- na.omit(master_res$variables_rm)
  removed_variables_outfn <- paste(outpath, paste0(outfn, "_eliminated_variables.csv"), 
					 sep = "/")
  fwrite(variables_rm_df, paste0(removed_variables_outfn, ".csv"))
  
  # Call the file to create output visualizations
  outpath_curr <- outpath
  source(paste(SOURCE_PATH, "process_dyscovr_output.R", sep = "/")) 
}

