############################################################
### CREATE PEER INPUTS 
### Written By: Sara Geraghty, August 2021
############################################################

# This file creates the necessary input files to run the PEER framework on 
# command line. 

#!/usr/bin/env Rscript

library(peer)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
suppressPackageStartupMessages(library(argparse))

# Link to tutorial: https://github.com/PMBio/peer/wiki/Tutorial
# Sample code from the GTEx consortium: https://vatlab.github.io/sos-docs/doc/examples/RNASeqGTEx.html

# Main path
main_path <- "/Genomics/grid/users/scamilli/thesis_work/run-peer/"

#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

############################################################
### SET UP PARSER FOR INPUTS
############################################################
# Create parser object
parser <- ArgumentParser()

# Specify desired options for the parser

# Expression DF
parser$add_argument("-e", "--expression_df", default = "tmm_normalized_expression_counts.csv", 
                    type = "character", help = "The name of the expression matrix to input.")

# Covariate DF
parser$add_argument("--covariate_df", default = "combined_patient_sample_DF_cibersort_abs.csv", type = "character",
                    help = "The name of the covariate matrix to input.")

# Logical options to specify whether or not we are including the covariate matrix, 
# adding the mean as another option, adjusting the defaults, etc.
parser$add_argument("--incl_CovMat", default = "TRUE", type = "character",
                    help = "A true/false value indicating whether or not we are including the covariate matrix. [default %(default)s]")
parser$add_argument("--incl_Mean", default = "FALSE", type = "character",
                    help = "A true/false value indicating whether or not we are including the expression mean as an additional covariate. [default %(default)s]")
parser$add_argument("--top10Kgenes", default = "FALSE", type = "character",
                    help = "A true/false value indicating whether or not we are including only the top 10K expressed genes in the analysis (like GTEx). [default %(default)s]")
parser$add_argument("--cancerType", default = "BRCA", type = "character",
                    help = "The cancer type we are looking at (or, pan-cancer for all cancer types. [default %(default)s]")

# Arguments for results writing
parser$add_argument("--expType", default = "TMM", type = "character",
                    help = "The type of expression (for results writing). TMM, TPM, R-N, Q-N, Counts, or FPKM. [default %(default)s]")
parser$add_argument("--covType", default = "CIBERSORT_Abs", type = "character",
                    help = "The type of covariate matrix (for results writing). TIMER, CIBERSORT, CIBERSORT_Abs, CIBERSORT_Abs_Tot_Frac, or CIBERSORT_Abs_Tot_Frac_Bucket [default %(default)s]")


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

incl_CovMat <- str2bool(args$incl_CovMat)
incl_Mean <- str2bool(args$incl_Mean)
top10Kgenes <- str2bool(args$top10Kgenes)

print(incl_CovMat)
print(incl_Mean)
print(top10Kgenes)


############################################################
### IMPORT AND INVERT EXPRESSION MATRIX
############################################################
expression_df <- read.csv(paste(main_path, paste("expression_matrices/", 
                                                 paste(args$cancerType, args$expression_df, 
                                                       sep = "/"), sep = ""), sep = ""),
                          header = TRUE, check.names = FALSE, row.names = 1)

# Invert matrix so samples are rows and genes are columns
expression_df_t <- t(expression_df)
print(dim(expression_df_t))

# Potentially keep only the top X genes 
#' Function to keep the top genes, from GTEx consortium
#' @param exp_mat expression matrix, where genes are the columns and samples are rows
#' @param num the number of top genes to keep
getTopGenes <- function(exp_mat, num = 10000) {
  if (ncol(exp_mat) <= num) {
    return(exp_mat)
  } else {
    top.index <- whichpart(total.expr <- colSums(exp_mat, na.rm = TRUE), num)
    return(exp_mat[,top.index])
  }
}

whichpart <- function(x, n) {
  nx <- length(x)   # the number of genes
  p <- nx-n         # the number of genes - 10K
  xp <- sort(x, partial=p)[p]  # we get what the value is of the "# of genes - n"st element is
  which(x > xp)     # get all those colSums greater than this (e.g. top 10K)
}

if(top10Kgenes) {
  print("getting top 10K genes")
  expression_df_t <- getTopGenes(expression_df_t)
}

############################################################
### READ IN COVARIATE MATRIX
############################################################
# Covariate matrix has samples as rows and covariates as columns
# Include sample and patient covariates

#covar_matrix <- read.csv(paste(main_path, "Linear Model/combined_patient_sample_DF_cibersort_total_frac.csv", sep = ""))
if(incl_CovMat) {
  covar_matrix <- read.csv(paste(main_path, paste("covariate_matrices/", paste(args$cancerType, args$covariate_df, sep = "/"), sep = ""), sep = ""),
                           header = TRUE, check.names = FALSE, row.names = 1)
  print(dim(covar_matrix))
}


############################################################
### RUN PEER IN R 
############################################################

# Create the model
print("creating model")
model = PEER()

# Set the observed data
print("adding expression matrix")
PEER_setPhenoMean(model,as.matrix(expression_df_t))
print(dim(PEER_getPhenoMean(model)))


# Automatically add an additional factor to account for mean expression
# This seems to only be useful when we are looking at non-normalized expression
if(incl_Mean) {
  print("adding mean")
  PEER_setAdd_mean(model, TRUE)
}

if(incl_CovMat) {
  print("setting covariate matrix")
  PEER_setCovariates(model, as.matrix(covar_matrix))
}

# num_confounders <- 10
# num_confounders <- ncol(expression_df_t) * .25

#' Function to get the number of confounders from GTEX
#' @param ss sample size
getNumPeer <- function(ss) {
  if (ss<150) return (min(15, ceiling(ss / 5)))
  else if (ss >=150 && ss < 250) return(30)
  else return(35)
}
print("calculating number of confounders")
num_confounders <- getNumPeer(nrow(expression_df_t))

PEER_setNk(model, num_confounders)
PEER_getNk(model)   # should equal num_confounders

print("beginning to update model")
PEER_update(model) # Will converge after a certain number of iterations

############################################################
### WRITE TO FILES AND CREATE PLOTS
############################################################
print("model is complete!")
factors = PEER_getX(model)
factors <- factors[,colSums(is.na(factors))<nrow(factors)]
print(dim(factors))

weights = PEER_getW(model)
weights <- weights[,colSums(is.na(weights))<nrow(weights)]
print(dim(weights))

precision = PEER_getAlpha(model)
precision <- precision[,colSums(is.na(precision))<nrow(precision)]
print(dim(precision))

residuals = PEER_getResiduals(model)
residuals <- residuals[,colSums(is.na(residuals))<nrow(residuals)]
print(dim(residuals))

# Write results to files
print("writing results to files")
if(incl_CovMat) {cov_mat <- args$covType} else {cov_mat <- "noCov"}
if(incl_Mean) {mean <- "inclMean"} else {mean <- "noMean"}
if(top10Kgenes) {top10k <- "top10k"} else {top10k <- "allGenes"}
footer <- paste(args$expType, paste(cov_mat, paste(mean, top10k, sep = ", "), sep = ", "), sep = ", ")

c <- paste0("InferredCov", 1:ncol(factors))

colnames(factors) <- c
rownames(factors) <- rownames(expression_df_t)

rownames(precision) <- c
colnames(precision) <- "Alpha"
precision <- as.data.frame(precision)

precision$Relevance <- 1.0 / precision$Alpha

rownames(residuals) <- rownames(expression_df_t)
colnames(residuals) <- colnames(expression_df_t)

fn_factors <-  paste(main_path, paste("results/", paste(args$cancerType, paste("peer_factors (", 
                               paste(footer, ").csv", sep = ""), sep = ""), sep = "/"), sep = ""), sep = "")
write.csv(t(factors), fn_factors)

fn_prec <-  paste(main_path, paste("results/", paste(args$cancerType, paste("peer_precision (", 
                               paste(footer, ").csv", sep = ""), sep = ""), sep = "/"), sep = ""), sep = "")
write.csv(precision, fn_prec)

fn_residuals <-  paste(main_path, paste("results/", paste(args$cancerType, paste("peer_residuals (", 
                               paste(footer, ").csv", sep = ""), sep = ""), sep = "/"), sep = ""), sep = "")
write.csv(residuals, fn_residuals)


# Plot the results
print("observing output")

fn_prec <- paste(main_path, paste("results/", paste(args$cancerType, paste("Precision (", 
                               paste(footer, ").png", sep = ""), sep = ""), sep = "/"), sep = ""), sep = "")
png(fn_prec, width = 400, height = 400)
plot(precision)
dev.off()

fn_weights <- paste(main_path, paste("results/", paste(args$cancerType, paste("Weights (", 
                               paste(footer, ").png", sep = ""), sep = ""), sep = "/"), sep = ""), sep = "")
png(fn_weights, width = 400, height = 400)
plot(weights)
dev.off()

fn_factors <- paste(main_path, paste("results/", paste(args$cancerType, paste("Factors (", 
                               paste(footer, ").png", sep = ""), sep = ""), sep = "/"), sep = ""), sep = "")
png(fn_factors, width = 400, height = 400)
plot(factors)
dev.off()

fn_residuals <- paste(main_path, paste("results/", paste(args$cancerType, paste("Residuals (", 
                              paste(footer, ").png", sep = ""), sep = ""), sep = "/"), sep = ""), sep = "")
png(fn_residuals, width = 400, height = 400)
plot(residuals)
dev.off()


############################################################
### OPTIONAL INFERENCE PARAMS
############################################################
# Set the max number of iterations (default is 1000)
#PEER_setNmax_iterations(model, 100)

# Set the tolerance (default:  (bound=0.001, variance=0.00001))
#PEER_setTolerance(model, 1)
#PEER_setVarTolerance(model, 0.1)

# Specify the alpha and beta of the gamma prior
# PEER uses uninformative priors on weight precision and noise precision by default
# (Alpha a = 0.001, Alpha b = 0.1, Eps a = 0.1, Eps b = 10)
#PEER_setPriorAlpha(model,0.001,0.1)
#PEER_setPriorEps(model,0.1,10.)


# Optional Add-Ons (To Overlay)

# Import clinical DFs 
brca_clinical_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/clinical.csv",
                          header = TRUE)
#pc_clinical_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/clinical.csv",
                        #header = TRUE)


############################################################
### ADD PATIENT INFO TO OVERLAY
############################################################
#' Extract information about age, gender, and race and add as rows for each patient
#' to the factors DF
#' @param factors the factors DF output from PEER
#' @param clin_df the clinical DF for the cohort of interest
add_patient_info <- function(factors, clin_df) {
  
  # Add three additional rows to the factor DF
  df <- data.frame(matrix(nrow = 3, ncol = ncol(factors)))
  colnames(df) <- colnames(factors)
  rownames(df) <- c("age", "race", "gender")
  factors <- rbind(factors, df)
  
  for (i in 1:ncol(factors)) {
    samp <- colnames(factors)[i]
    patient <- paste(unlist(strsplit(samp, "-", fixed = TRUE))[1:3], collapse = "-")
    print(patient)
    
    patient_row <- clin_df[clin_df$case_submitter_id == patient, ]

    age <-  unique(patient_row[,'age_at_index'])
    race <- unique(patient_row[, 'race'])
    gender <- unique(patient_row[ , 'gender']) 
    
    factors[rownames(factors) == 'age', i] <- age
    factors[rownames(factors) == 'race', i] <- race
    factors[rownames(factors) == 'gender', i] <- gender
  }
  return(factors)
}

# Read back if necessary
factors <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/PEER/peer_factors (FPKM_TO, noCov, noMean, allGenes).csv", 
                    header = TRUE, check.names = FALSE, row.names = 1)
# If necessary, remove the TARGET, other project samples
factors <- factors[,grepl("TCGA", colnames(factors))]

# Call function
factors_new <- t(add_patient_info(factors, brca_clinical_df))
#factors_new <- t(add_patient_info(factors, pc_clinical_df))

# Write to file
fwrite(factors_new, "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/PEER/peer_factors (R-N, noCov, noMean, top10K, withPatChar).csv")
#fwrite(factors_new, "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/PEER/peer_factors (R-N, noCov, noMean, top10K, withPatChar).csv")

############################################################
### CREATE PLOTS WITH DEMOGRAPHIC INFO
############################################################
factors_new <- as.data.frame(factors_new)
factors_new$age <- as.integer(factors_new$age)
factors_new$InferredCov1 <- as.numeric(factors_new$InferredCov1)
factors_new$InferredCov2 <- as.numeric(factors_new$InferredCov2)
factors_new$names <- rownames(factors_new)

factors_new$race <- unlist(lapply(factors_new$race, function(r) if(r == "") 
{return("not reported")} else {return(r)}))

# Tumor or normal
factors_new$tumOrNorm <- unlist(lapply(rownames(factors_new), function(id) {
  if(grepl("-0", id, fixed = TRUE)) {return("Tumor")}
  else {return("Normal")}
}))

# Make everything numeric
factors_new <- cbind(apply(factors_new[,1:35], MARGIN = 2, function(x) as.numeric(x)),
                     factors_new[,36:ncol(factors_new)])

ggplot(factors_new, aes(x = InferredCov1, y = InferredCov2, color = age)) + 
  geom_point() #+ theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  #geom_text(aes(label= ifelse(InferredCov1 > quantile(InferredCov1, 0.95),
                              #as.character(names),'')),hjust=0,vjust=0)

ggplot(factors_new, aes(x = InferredCov1, y = InferredCov2, color = race)) + 
  geom_point() #+ theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  #geom_text(aes(label= ifelse(InferredCov1 > quantile(InferredCov1, 0.95),
                              #as.character(names),'')),hjust=0,vjust=0)

ggplot(factors_new, aes(x = InferredCov1, y = InferredCov2, color = tumOrNorm)) + 
  geom_point() #+ theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
#geom_text(aes(label= ifelse(InferredCov1 > quantile(InferredCov1, 0.95),
#as.character(names),'')),hjust=0,vjust=0)

# Directly correlate these factors with age and race

# Check for normality of variables
shapiro.test(factors_new$age)
shapiro.test(factors_new$InferredCov1)
shapiro.test(factors_new$InferredCov2)

# Calculate correlations
cor_age_InferredCov1 <- cor.test(factors_new$InferredCov1, factors_new$age, method = "spearman")
cor_age_InferredCov2 <- cor.test(factors_new$InferredCov2, factors_new$age, method = "spearman")

# Plot correlations
ggscatter(factors_new, x = "age", y = "InferredCov1", add = "reg.line", conf.int = TRUE,
          cor.coeff = TRUE, cor.method = "spearman", xlab = "Age", ylab = "InferredCov1") + 
  stat_cor(method = "spearman")
ggscatter(factors_new, x = "age", y = "InferredCov2", add = "reg.line", conf.int = TRUE,
          cor.coeff = TRUE, cor.method = "spearman", xlab = "Age", ylab = "InferredCov2") + 
  stat_cor(method = "spearman")

# For race, do an ANOVA and also pairwise comparisons
#ggscatter(factors_new, x = "race", y = "InferredCov1", add = "reg.line", conf.int = TRUE,
#cor.coeff = TRUE, cor.method = "spearman", xlab = "Race", ylab = "InferredCov1")
ggboxplot(factors_new, x = "race", y = "InferredCov1", palette = "jco", ylab = "InferredCov1") + 
  stat_compare_means(method = "anova")
compare_means(InferredCov1 ~ race, data = factors_new)
my_comparisons <- list(c("white", "black or african american"), c("white", "not reported"),
                       c("asian", "black or african american"), c("asian", "not reported"),
                       c("black or african american", "not reported"))
ggboxplot(factors_new, x = "race", y = "InferredCov1", palette = "jco", ylab = "InferredCov1") + 
  stat_compare_means(comparisons = my_comparisons) #+ stat_compare_means(label.y = 50)


#ggscatter(factors_new, x = "race", y = "InferredCov2", add = "reg.line", conf.int = TRUE,
#cor.coeff = TRUE, cor.method = "spearman", xlab = "Race", ylab = "InferredCov2")
ggboxplot(factors_new, x = "race", y = "InferredCov2", palette = "jco", ylab = "InferredCov2") + 
  stat_compare_means(method = "anova")
compare_means(InferredCov2 ~ race, data = factors_new)
my_comparisons <- list(c("white", "black or african american"), c("asian", "black or african american"), 
                       c("black or african american", "not reported"))
ggboxplot(factors_new, x = "race", y = "InferredCov2", palette = "jco", ylab = "InferredCov2") + 
  stat_compare_means(comparisons = my_comparisons) #+ stat_compare_means(label.y = 50)


# Make a wrapped plot of all of the Inferred Covariates
covs <- unlist(lapply(1:35, function(i) paste("InferredCov", i, sep = "")))
ggscatter(factors_new, x = "age", y = covs, add = "reg.line", conf.int = TRUE,
          cor.coeff = TRUE, cor.method = "spearman", xlab = "Age") +
  stat_cor(method = "spearman") + facet_wrap()
  

############################################################
### OPTIONAL: CREATE TUMOR/ NORMAL COVARIATE MATRICES FOR 
# PEER
############################################################
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

expression_df_fpkm <- read.csv(paste(main_path, "Expression/expression_fpkm_gn_DF.csv", sep = ""),
                       header = TRUE, check.names = FALSE, row.names = 1)
expression_df_qn <- read.csv(paste(main_path, "Expression/expression_quantile_norm_DF.csv", sep = ""),
                       header = TRUE, check.names = FALSE, row.names = 1)
expression_df_rn <- read.csv(paste(main_path, "Expression/expression_rank_norm_DF.csv", sep = ""),
                       header = TRUE, check.names = FALSE, row.names = 1)
expression_df_tmm <- read.csv(paste(main_path, "Expression/tmm_normalized_expression_counts.csv", sep = ""),
                       header = TRUE, check.names = FALSE, row.names = 1)

#' For the given expression matrix, determines if each sample is tumor or normal 
#' and creates/returns a two-column data frame with each sample ID and the tumor/
#' normal status
#' @param expression_df a normalized expression DF (columns are samples, rows are
#' ENSG IDs)
create_tumNorm_covMat <- function(expression_df) {
  new_rows <- lapply(colnames(expression_df), function(nam) {
    status <- ifelse(grepl("-0", nam), 0, 1)
    new_df <- data.frame(status)
    colnames(new_df) <- "Status"
    rownames(new_df) <- nam
    return(new_df)
  })
  # Bind all the rows together
  cov_mat <- do.call(rbind, new_rows)
  print(head(cov_mat))
  return(cov_mat)
}

# Call function
covMat_fpkm <- create_tumNorm_covMat(expression_df_fpkm)
covMat_qn <- create_tumNorm_covMat(expression_df_qn)
covMat_rn <- create_tumNorm_covMat(expression_df_rn)
covMat_tmm <- create_tumNorm_covMat(expression_df_tmm)

# Write to files
write.csv(covMat_fpkm, paste(main_path, "PEER/Covariate_Matrices/fpkm_tumNorm_covs.csv", sep = ""))
write.csv(covMat_qn, paste(main_path, "PEER/Covariate_Matrices/qn_tumNorm_covs.csv", sep = ""))
write.csv(covMat_rn, paste(main_path, "PEER/Covariate_Matrices/rn_tumNorm_covs.csv", sep = ""))
write.csv(covMat_tmm, paste(main_path, "PEER/Covariate_Matrices/tmm_tumNorm_covs.csv", sep = ""))


############################################################
### OPTIONAL: GET CANCER SAMPLE ONLY VERSIONS OF EXP. MATRICES
############################################################
# Get cancer sample-only expression matrices for the non-tumor-normal-matched cases
expression_df_fpkm_to <- expression_df_fpkm[,grepl("-0", colnames(expression_df_fpkm))]
expression_df_qn_to <- expression_df_qn[,grepl("-0", colnames(expression_df_qn))]
expression_df_rn_to <- expression_df_rn[,grepl("-0", colnames(expression_df_rn))]
expression_df_tmm_to <- expression_df_tmm[,grepl("-0", colnames(expression_df_tmm))]

# Write to files
write.csv(expression_df_fpkm_to, paste(main_path, "Expression/expression_fpkm_gn_DF_TO.csv", sep = ""))
write.csv(expression_df_qn_to, paste(main_path, "Expression/expression_quantile_norm_DF_TO.csv", sep = ""))
write.csv(expression_df_rn_to, paste(main_path, "Expression/expression_rank_norm_DF_TO.csv", sep = ""))
write.csv(expression_df_tmm_to, paste(main_path, "Expression/expression_tmm_DF_TO.csv", sep = ""))

