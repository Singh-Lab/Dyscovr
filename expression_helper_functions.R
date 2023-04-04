############################################################
### Expression Helper Functions
### Written By: Sara Camilli, April 2021
############################################################

# This file contains helper functions for the visualization of expression profiles of 
# particular genes under different conditions

library(dplyr)
library(edgeR)
library(data.table)

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Expression/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Expression/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)

############################################################
### DEFINE A REGULATORY PROTEIN OF INTEREST
############################################################
tp53 <- "P04637" # TP53
tp53_ensg <- "ENSG00000141510"

foxa1 <- "P55317" #FOXA1
foxa1_ensg <- "ENSG00000129514"

pik3ca <- "P42336"
pik3ca_ensg <- "ENSG00000121879"

ttn <- "Q8WZ42"
ttn_ensg <- "ENSG00000155657"

kras <- "P01116"
kras_ensg <- "ENSG00000133703"

############################################################
### IMPORT FILES NECESSARY TO CREATE A STARTER DF FOR THIS
### PROTEIN OF INTEREST
############################################################
# Non-Tumor-Normal Matched Files
methylation_df <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Methylation/methylation_M_CancerOnly_IntersectPatients.csv", sep = ""), 
                           header = TRUE, row.names = 1, check.names = FALSE)
mutation_targ_df <- read.csv(paste(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_nonsense_CancerOnly_IntersectPatients.csv", sep = ""), 
                             header = TRUE, row.names = 1, check.names = FALSE)
mutation_regprot_df <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Regprot_Mutation/iprotein_results_nonsynonymous_IntersectPatients.csv", sep = ""), 
                                header = TRUE, row.names = 1, check.names = FALSE)
cna_df <- read.csv(paste(main_path, "Linear Model/Tumor_Only/CNV/CNA_AllGenes_CancerOnly_IntersectPatients.csv", sep = ""), 
                   header = TRUE, row.names = 1, check.names = FALSE)
rownames(cna_df) <- unlist(lapply(rownames(cna_df), function(x) unlist(strsplit(x, ".", fixed = TRUE))[1]))
patient_df <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Patient/combined_patient_sample_cibersort_total_frac_tmm_IntersectPatients.csv", sep = ""), 
                       header = TRUE, row.names = 1, check.names = FALSE)

methylation_df <- as.data.table(methylation_df)
mutation_targ_df <- as.data.table(mutation_targ_df)
mutation_regprot_df <- as.data.table(mutation_regprot_df)
cna_df$ensg_id <- rownames(cna_df)
cna_df <- as.data.table(cna_df)
patient_df$sample_id <- rownames(patient_df)
patient_df <- as.data.table(patient_df)

# See "fill_regprot_inputs" function in linear_model.R file
starter_df_foxa1 <- fill_regprot_inputs(patient_df, foxa1, foxa1_ensg, mutation_regprot_df, methylation_df, 
                                  cna_df, "eQTL")
starter_df_tp53 <- fill_regprot_inputs(patient_df, tp53, tp53_ensg, mutation_regprot_df, methylation_df, 
                                        cna_df, "rawCNA", FALSE, FALSE, FALSE)

############################################################
### DEFINE A TARGET GENE OF INTEREST
############################################################
ensg <- "ENSG00000073584"


############################################################
### IMPORT EXPRESSION DATAFRAMES
############################################################
expression_df_fpkm <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_fpkm_CancerOnly_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
expression_df_tmm <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_tmmSDGr1filtByExpr_CancerOnly_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
expression_df_quantile_norm <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_quantile_norm_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
expression_df_rank_norm <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_rank_norm_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)

############################################################
### CREATE OUTLIER FILTERED FPKM AND TMM DATAFRAMES
############################################################
# See 'filter_expression_df' function in linear_model_helper_functions.R
  # Will split the patient groups based on whether or not they have a mutation in the given
  # target gene, and then filter outliers in each group using a 1.5*IQR + Q3 schema
starter_df_fpkm_filt_foxa1 <- filter_expression_df(expression_df_fpkm, starter_df_foxa1, ensg)
starter_df_tmm_filt_foxa1 <- filter_expression_df(expression_df_tmm, starter_df_foxa1, ensg)

starter_df_fpkm_filt_tp53 <- filter_expression_df(expression_df_fpkm, starter_df_tp53, ensg)
starter_df_tmm_filt_tp53 <- filter_expression_df(expression_df_tmm, starter_df_tp53, ensg)

# Get the filtered expression DF from the filtered starter DF
patients_fpkm_filt_foxa1 <- rownames(starter_df_fpkm_filt_foxa1)
patients_tmm_filt_foxa1 <- rownames(starter_df_tmm_filt_foxa1)

patients_fpkm_filt_tp53 <- rownames(starter_df_fpkm_filt_tp53)
patients_tmm_filt_tp53 <- rownames(starter_df_tmm_filt_tp53)

patient_ids_fpkm <- unlist(lapply(colnames(expression_df_fpkm), function(x) 
  unlist(strsplit(x, "-", fixed = TRUE))[3]))
patient_ids_tmm <- unlist(lapply(colnames(expression_df_tmm), function(x) 
  unlist(strsplit(x, "-", fixed = TRUE))[3]))

cols_to_keep_fpkm_foxa1 <- which(patient_ids_fpkm %fin% patients_fpkm_filt_foxa1)
cols_to_keep_tmm_foxa1 <- which(patient_ids_tmm %fin% patients_tmm_filt_foxa1)
cols_to_keep_fpkm_tp53 <- which(patient_ids_fpkm %fin% patients_fpkm_filt_tp53)
cols_to_keep_tmm_tp53 <- which(patient_ids_tmm %fin% patients_tmm_filt_tp53)

expression_df_fpkm_filt_foxa1 <- expression_df_fpkm[,cols_to_keep_fpkm_foxa1]
expression_df_tmm_filt_foxa1 <- expression_df_tmm[,cols_to_keep_tmm_foxa1]
expression_df_fpkm_filt_tp53 <- expression_df_fpkm[,cols_to_keep_fpkm_tp53]
expression_df_tmm_filt_tp53 <- expression_df_tmm[,cols_to_keep_tmm_tp53]


############################################################
### VISUALIZE EXPRESSION OF TARGET GENE
############################################################
#' Visualize the expression of a particular target gene using a scatterplot
#' across all patients
#' @param expression_df a data frame with expression values (columns are patient IDs-sample IDs, 
#' rows are ENSG IDs)
#' @param ensg the ENSG ID of the target gene of interest
visualize_tg_expression <- function(expression_df, ensg) {
  plot(as.numeric(expression_df[expression_df$ensg_id == ensg,]), xlab = "Patient", ylab = "Expression")
  abline(h = mean(as.numeric(expression_df[expression_df$ensg_id == ensg,])), col = "red")
}


############################################################
### VISUALIZE EXPRESSION OF TARGET GENE, BY CANCER V. NORMAL
############################################################
#' Get expression profiles for patient cancer samples vs. normal samples
#' and visualize using a box plot
#' @param expression_df a data frame with expression values (columns are patient IDs-sample IDs, 
#' rows are ENSG IDs)
#' @param ensg the ENSG ID of the target gene of interest
visualize_tg_expression_cvn <- function(expression_df, ensg) {
  expression_cancer <- as.numeric(expression_df[rownames(expression_df) == ensg, grepl("-0", colnames(expression_df))])
  expression_normal <- as.numeric(expression_df[rownames(expression_df) == ensg, grepl("-11", colnames(expression_df))])
  plot(expression_cancer, xlab = "Patient", ylab = "Expression", col = "red")
  points(expression_normal, col = "blue")
  legend(500, 3000, legend = c("Cancer", "Normal"), col = c("red", "blue"), lty = 1)  # Adjust placement of legend 
  
  # Make box plot versions
  dataList <- list("Cancer" = expression_cancer, "Normal" = expression_normal)
  boxplot(dataList, ylab = "Expression")
  
  #simple t-test
  ttest_res <- t.test(expression_normal, expression_cancer)
  print(ttest_res$p.value)
}


############################################################
### VISUALIZE EXPRESSION OF TARGET GENE, BY MUTATED/ UNMUTATED
### IN REGULATORY PROTEIN OF INTEREST
############################################################
#' Get a vector of patients that have a mutation in this given regulatory protein 
#' @param mutation_regprot_df a data frame with information about mutations in given protein domains
#' @param regprot_name the SWISSPROT ID of the regulatory protein of interest
get_mutant_patients <- function(mutation_regprot_df, regprot_name) {
  regprot_rows <- which(unlist(lapply(mutation_regprot_df$Query, function(x) 
    unlist(strsplit(unlist(strsplit(x, "|", fixed = TRUE))[3], "_", fixed = TRUE))[1] == regprot_name)))
  mutant_patients <- unique(unlist(strsplit(mutation_regprot_df[regprot_rows, "Patient"], ";", fixed = TRUE)))
  return(mutant_patients)
}

mutation_regprot_df <- as.data.frame(mutation_regprot_df)
foxa1_mutant_patients <- get_mutant_patients(mutation_regprot_df, "FOXA1")
tp53_mutant_patients <- get_mutant_patients(mutation_regprot_df, "P53")
pik3ca_mutant_patients <- get_mutant_patients(mutation_regprot_df, "PK3CA")
kras_mutant_patients <- get_mutant_patients(mutation_regprot_df, "KRAS")

#' Get a vector of patients that have a deletion in this given regulatory protein 
#' @param cna_df a data frame with information about CNAs in given genes
#' @param regprot_ensg the ENSG ID of the regulatory protein of interest
get_deleted_patients <- function(cna_df, regprot_ensg) {
  cnas <- unlist(cna_df[cna_df$ensg_id == regprot_ensg, ])
  to_keep <- which(cnas < 2)
  samps <- colnames(cna_df)[to_keep]
  return(samps)
}

tp53_deleted_patients <- get_deleted_patients(cna_df, tp53_ensg)

tp53_disrupted_patients <- unique(c(tp53_deleted_patients, tp53_mutant_patients))

get_amplified_patients <- function(cna_df, regprot_ensg) {
  cnas <- unlist(cna_df[cna_df$ensg_id == regprot_ensg, ])
  to_keep <- which(cnas > 2)
  samps <- colnames(cna_df)[to_keep]
  return(samps)
}

pik3ca_amplified_patients <- get_amplified_patients(cna_df, pik3ca_ensg) 

pik3ca_disrupted_patients <- unique(c(pik3ca_amplified_patients, pik3ca_mutant_patients))


# Limit the expression data frame to just cancer samples
express_df_cancer_counts <- expression_df_counts[,!(grepl("-11", colnames(expression_df_counts)))]
express_df_cancer_fpkm <- expression_df_fpkm[,!(grepl("-11", colnames(expression_df_fpkm)))]
express_df_cancer_tmm <- expression_df_tmm[,!(grepl("-11", colnames(expression_df_tmm)))]
express_df_cancer_quant_norm <- expression_df_quantile_norm[,!(grepl("-11", colnames(expression_df_quantile_norm)))]
express_df_cancer_rank_norm <- expression_df_rank_norm[,!(grepl("-11", colnames(expression_df_rank_norm)))]

express_df_cancer_fpkm_filt_tp53 <- expression_df_fpkm_filt_tp53[,!(grepl("-11", colnames(expression_df_fpkm_filt_tp53)))]
express_df_cancer_tmm_filt_tp53 <- expression_df_tmm_filt_tp53[,!(grepl("-11", colnames(expression_df_tmm_filt_tp53)))]
express_df_cancer_fpkm_filt_foxa1 <- expression_df_fpkm_filt_foxa1[,!(grepl("-11", colnames(expression_df_fpkm_filt_foxa1)))]
express_df_cancer_tmm_filt_foxa1 <- expression_df_tmm_filt_foxa1[,!(grepl("-11", colnames(expression_df_tmm_filt_foxa1)))]

#' Get the expression in each of the groups (mutated/ unmutated)
#' @param express_df_cancer a data frame with expression values (columns are patient IDs-sample IDs, 
#' rows are ENSG IDs), subsetted to only cancer samples
#' @param mutant_patients a vector of all the patients that have a mutation in the given
#' regulatory protein
#' @param ensg the ENSG ID of the target gene of interest
get_expression_by_mut_group <- function(express_df_cancer, mutant_patients, ensg) {
  #colnames(express_df_cancer) <- unlist(lapply(colnames(express_df_cancer), function(x) unlist(strsplit(x, "-", fixed = TRUE))[1]))
  expression_mutants <- as.numeric(express_df_cancer[rownames(express_df_cancer) == ensg, colnames(express_df_cancer) %in% mutant_patients])
  expression_normal <- as.numeric(express_df_cancer[rownames(express_df_cancer) == ensg, !(colnames(express_df_cancer) %in% mutant_patients)])
  return(list("mutants" = expression_mutants, "normal" = expression_normal))
}

get_expression_by_mut_group_incl_norm <- function(expression_df, mutant_patients, ensg) {
  expression_normal <- as.numeric(expression_df[rownames(expression_df) == ensg, grepl("-11", colnames(expression_df))])
  express_df_cancer <- expression_df[,!(grepl("-11", colnames(expression_df)))]
  #colnames(express_df_cancer) <- unlist(lapply(colnames(express_df_cancer), function(x) paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))
  print(head(express_df_cancer))
  expression_mutants <- as.numeric(express_df_cancer[rownames(express_df_cancer) == ensg, colnames(express_df_cancer) %in% mutant_patients])
  expression_wt <- as.numeric(express_df_cancer[rownames(express_df_cancer) == ensg, !(colnames(express_df_cancer) %in% mutant_patients)])
  return(list("normal" = expression_normal, "cancer_wt" = expression_wt, "cancer_mutants" = expression_mutants))
}


#' Get the expression in each of the functional copy groups 
#' @param express_df_cancer a data frame with expression values (columns are patient IDs-sample IDs, 
#' rows are ENSG IDs), subsetted to only cancer samples
#' @param funct_copy_dict a two-column data frame with sample ID and # of functional copies
#' @param ensg the ENSG ID of the target gene of interest
#' @param ts_or_onco a character value indicating whether or not this is a "tumor_suppressor" 
#' or an "oncogene"
get_expression_by_fnc <- function(express_df_cancer, funct_copy_dict, ensg, ts_or_onco) {
  funct_copy_list <- list()
  num_funct_copy_options <- length(unique(funct_copy_dict$Num.Functional.Copies))
  for (i in 1:num_funct_copy_options) {
    samp_ids <- funct_copy_dict[funct_copy_dict$Num.Functional.Copies == (i-1), 'sample.id']
    funct_copy_list[[i]] <- as.numeric(express_df_cancer[rownames(express_df_cancer) == ensg, colnames(express_df_cancer) %in% samp_ids])
  }
  if(ts_or_onco == "tumor_suppressor") {
    grouped_list <- list("Full KO" = funct_copy_list[[1]], "Partial KO" = funct_copy_list[[2]],
                         "Normal/Amp" = unlist(funct_copy_list[3:length(funct_copy_list)]))
    #grouped_list <- list("Full/Partial KO" =  c(funct_copy_list[[1]], funct_copy_list[[2]]),
                          #"Normal/Amp" = unlist(funct_copy_list[3:length(funct_copy_list)]))
  } else if (ts_or_onco == "oncogene") {
    grouped_list <- list("Full or Partial KO" = c(funct_copy_list[[1]], funct_copy_list[[2]]), "Normal" = funct_copy_list[[3]],
                         "Amplification" = unlist(funct_copy_list[4:length(funct_copy_list)]))
    
  } else {print("Must supply either 'tumor_suppressor' or 'oncogene'.")}
  return(grouped_list)
}

expression_tmm <- get_expression_by_fnc(expression_df_tmm, output_df[,c("sample.id", "Num.Functional.Copies")], ensg, "tumor_suppressor")
boxplot(expression_tmm, ylab = "Expression (FPKM)", main = paste("CTSA", paste("Expression By", paste("TP53", "FNC Status"))))


### MAKE VISUALIZATIONS ###
#' Makes a boxplot to show the expression in each group (mutated/unmutated) for each 
#' type of expression normalization
#' @param mutant_patients vector of patients that have a mutation in the regulatory
#' protein of interest 
#' @param ensg the ENSG ID of the target gene of interest
#' @param regprot_name the gene name of the regulatory protein of interest
#' @param express_df_cancer_counts the expression data frame, limited to cancer patients,
#' measuring expression as counts
#' @param express_df_cancer_fpkm the expression data frame, limited to cancer patients,
#' measuring expression as FPKM
#' @param express_df_cancer_fpkm_filt the expression data frame, limited to cancer patients,
#' measuring expression as FPKM, filtered to exclude expression outliers
#' @param express_df_cancer_tmm the expression data frame, limited to cancer patients,
#' measuring expression as TMM
#' @param express_df_cancer_tmm_filt the expression data frame, limited to cancer patients,
#' measuring expression as TMM, filtered to exclude expression outliers
#' @param express_df_cancer_quant_norm the expression data frame, limited to cancer patients,
#' measuring expression as quantile-normalized expression
#' @param express_df_cancer_rank_norm the expression data frame, limited to cancer patients,
#' measuring expression as rank-normalized expression
make_expression_by_mut_group_vis <- function(mutant_patients, ensg, tg_name, regprot_name, 
                                             express_df_cancer_counts, express_df_cancer_fpkm,
                                             express_df_cancer_fpkm_filt, express_df_cancer_tmm, 
                                             express_df_cancer_tmm_filt, express_df_cancer_quant_norm, 
                                             express_df_cancer_rank_norm) {
  
  # COUNTS
  expression_counts <- get_expression_by_mut_group(express_df_cancer_counts, mutant_patients, ensg)
  expression_counts_mutants <- expression_counts[[1]]
  expression_counts_normal <- expression_counts[[2]]
  dataList_counts <- list("No Mutation" = expression_counts_normal, "Mutation" = expression_counts_mutants)
  boxplot(dataList_counts, ylab = "Expression (Counts)", main = paste(tg_name, paste("Expression By", paste(regprot_name, "Mutation Status"))))
  
  # FPKM
  expression_fpkm <- get_expression_by_mut_group(express_df_cancer_fpkm, mutant_patients, ensg)
  expression_fpkm_mutants <- expression_fpkm[[1]]
  expression_fpkm_normal <- expression_fpkm[[2]]
  dataList_fpkm <- list("No Mutation" = expression_fpkm_normal, "Mutation" = expression_fpkm_mutants)
  boxplot(dataList_fpkm, ylab = "Expression (FPKM)", main = paste(tg_name, paste("Expression By", paste(regprot_name, "Mutation Status"))))
  
  # FPKM, Outliers Filtered
  expression_fpkm_filt <- get_expression_by_mut_group(express_df_cancer_fpkm_filt, mutant_patients, ensg)
  expression_fpkm_filt_mutants <- expression_fpkm_filt[[1]]
  expression_fpkm_filt_normal <- expression_fpkm_filt[[2]]
  dataList_fpkm_filt <- list("No Mutation" = expression_fpkm_filt_normal, "Mutation" = expression_fpkm_filt_mutants)
  boxplot(dataList_fpkm_filt, ylab = "Expression (FPKM, Outliers Filtered)", 
          main = paste(tg_name, paste("Expression By", paste(regprot_name, "Mutation Status"))))
  
  # TMM
  expression_tmm <- get_expression_by_mut_group(express_df_cancer_tmm, mutant_patients, ensg)
  expression_tmm_mutants <- expression_tmm[[1]]
  expression_tmm_normal <- expression_tmm[[2]]
  dataList_tmm <- list("No Mutation" = expression_tmm_normal, "Mutation" = expression_tmm_mutants)
  boxplot(dataList_tmm, ylab = "Expression (TMM)", main = paste(tg_name, paste("Expression By", paste(regprot_name, "Mutation Status"))))
  
  # TMM, Outliers Filtered
  expression_tmm_filt <- get_expression_by_mut_group(express_df_cancer_tmm_filt, mutant_patients, ensg)
  expression_tmm_filt_mutants <- expression_tmm_filt[[1]]
  expression_tmm_filt_normal <- expression_tmm_filt[[2]]
  dataList_tmm_filt <- list("No Mutation" = expression_tmm_filt_normal, "Mutation" = expression_tmm_filt_mutants)
  boxplot(dataList_tmm_filt, ylab = "Expression (TMM, Outliers Filtered)", main = paste(tg_name, paste("Expression By", paste(regprot_name, "Mutation Status"))))
  
  # Quantile-Normalized
  expression_quantile_norm <- get_expression_by_mut_group(express_df_cancer_quant_norm, mutant_patients, ensg)
  expression_quantile_norm_mutants <- expression_quantile_norm[[1]]
  expression_quantile_norm_normal <- expression_quantile_norm[[2]]
  dataList_quantile_norm <- list("No Mutation" = expression_quantile_norm_normal, "Mutation" = expression_quantile_norm_mutants)
  boxplot(dataList_quantile_norm, ylab = "Expression (Quantile-Normalized)", main = paste(tg_name, paste("Expression By", paste(regprot_name, "Mutation Status"))))
  
  # Rank-Normalized
  expression_rank_norm <- get_expression_by_mut_group(express_df_cancer_rank_norm, mutant_patients, ensg)
  expression_rank_norm_mutants <- expression_rank_norm[[1]]
  expression_rank_norm_normal <- expression_rank_norm[[2]]
  dataList_rank_norm <- list("No Mutation" = expression_rank_norm_normal, "Mutation" = expression_rank_norm_mutants)
  boxplot(dataList_rank_norm, ylab = "Expression (Rank-Normalized)", main = paste(tg_name, paste("Expression By", paste(regprot_name, "Mutation Status"))))
}

# Sample call
make_expression_by_mut_group_vis(tp53_mutant_patients, ensg, "ZNF763", "TP53", 
                                             express_df_cancer_counts, express_df_cancer_fpkm,
                                             express_df_cancer_fpkm_filt_tp53, express_df_cancer_tmm, 
                                             express_df_cancer_tmm_filt_tp53, express_df_cancer_quant_norm, 
                                             express_df_cancer_rank_norm)


#' Visualization of a paired analysis of expression in cancer / expression in normal
#' for each driver mutational category. Also calculates and displays a paired 
#' Wilcoxon test between these mutational categories.
#' @param mutant_patients vector of patients that have a mutation in the regulatory
#' protein of interest 
#' @param ensg the ENSG ID of the target gene of interest
#' @param targ_name the name of the target gene of interest
#' @param regprot_name the gene name of the regulatory protein of interest
#' @param expression_df the expression data frame, limited to only patients in the 
#' given subtype of interest
make_paired_expr_by_mut_group_vis <- function(mutant_patients, ensg, targ_name, 
                                              regprot_name, expression_df) {
  expression_df <- expression_df[rownames(expression_df) == ensg, ]
  
  samples_normal <- colnames(expression_df)[grepl("-11", colnames(expression_df))]
  
  patients_normal <- unlist(lapply(samples_normal, function(x) unlist(strsplit(x, "-", fixed = TRUE))[1]))
  
  samples_cancer <- colnames(expression_df)[grepl("-0", colnames(expression_df))]
  patients_cancer <- unlist(lapply(samples_cancer, function(x) unlist(strsplit(x, "-", fixed = TRUE))[1]))
  
  intersecting_pats <- intersect(patients_normal, patients_cancer)
  print(length(intersecting_pats))
  
  samples_normal_sub <- unlist(lapply(samples_normal, function(x) 
    ifelse(unlist(strsplit(x, "-", fixed = TRUE))[1] %in% intersecting_pats, x, NA)))
  samples_normal_sub <- samples_normal_sub[!is.na(samples_normal_sub)]
  
  samples_cancer_sub <- unlist(lapply(samples_cancer, function(x) 
    ifelse(unlist(strsplit(x, "-", fixed = TRUE))[1] %in% intersecting_pats, x, NA)))
  samples_cancer_sub <- samples_cancer_sub[!is.na(samples_cancer_sub)]
  
  
  cancer_norm_ratios_mut <- c()
  cancer_norm_ratios_wt <- c()
  paired_vals_df <- data.frame(matrix(nrow = length(intersecting_pats), ncol = 3))
  colnames(paired_vals_df) <- c("cancer", "normal", "status")
  
  for(i in 1:length(samples_normal_sub)) {
    samp_norm <- samples_normal_sub[i]
    pat <- unlist(strsplit(samp_norm, "-", fixed = TRUE))[1]
    samp_cancer <- unique(samples_cancer_sub[grepl(pat, samples_cancer_sub)])
    if(length(samp_cancer) > 1) {samp_cancer <- samp_cancer[grepl("-01", samp_cancer)][1]}

    expression_normal <- as.numeric(expression_df[, colnames(expression_df) == samp_norm])
    expression_cancer <- as.numeric(expression_df[, colnames(expression_df) == samp_cancer])
    
    if(samp_cancer %in% mutant_patients) {
      cancer_norm_ratios_mut <- c(cancer_norm_ratios_mut, (expression_cancer / expression_normal))
      paired_vals_df[i,] <- c(expression_cancer, expression_normal, "Mut")
    } else {
      cancer_norm_ratios_wt <- c(cancer_norm_ratios_wt, (expression_cancer / expression_normal))
      paired_vals_df[i,] <- c(expression_cancer, expression_normal, "WT")
    }
  }
  
  dataList <- list("TP53 WT (Cancer / Normal)" = cancer_norm_ratios_wt,
                   "TP53 Mut. (Cancer / Normal)" = cancer_norm_ratios_mut)
  
  boxplot(dataList, ylab = "Cancer to Normal TMM Expression Ratio", xlab = paste(regprot_name, "Mutation Status"),
          main = paste(targ_name, paste("Tumor-to-Normal Expression Ratio By", paste(regprot_name, "Mutation Status"))))
  abline(h = 1, col = "red")
  
  wilcox_res <- wilcox.test(x = cancer_norm_ratios_wt, y = cancer_norm_ratios_mut, 
                            paired = FALSE, alternative = "less")
  print(wilcox_res)
  text(1, max(c(cancer_norm_ratios_wt, cancer_norm_ratios_mut)), 
       paste("Wilcoxon p-value: ", round(wilcox_res$p.value, digits = 4)))
  
  print(head(paired_vals_df))
  paired_wilcox_res <- wilcox.test(as.numeric(paired_vals_df$cancer), 
                                   as.numeric(paired_vals_df$normal), 
                                  paired = TRUE, alternative = "greater")
  print(paired_wilcox_res)
  
  paired_wilcox_res_wt <- wilcox.test(as.numeric(paired_vals_df[paired_vals_df$status == "WT", 'cancer']), 
                                      as.numeric(paired_vals_df[paired_vals_df$status == "WT", 'normal']), 
                                      paired = TRUE, alternative = "greater")
  print(paired_wilcox_res_wt)
  
  paired_wilcox_res_mut <- wilcox.test(as.numeric(paired_vals_df[paired_vals_df$status == "Mut", 'cancer']), 
                                      as.numeric(paired_vals_df[paired_vals_df$status == "Mut", 'normal']), 
                                      paired = TRUE, alternative = "greater")
  print(paired_wilcox_res_mut)
}

colnames(expression_df_tmm) <- unlist(lapply(colnames(expression_tmm), function(x) 
  paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))
make_paired_expr_by_mut_group_vis(tp53_mutant_patients, ensg, "SHMT2", "TP53", expression_df_tmm)


# For performing analyses on particular subtypes
patient_set <- read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt"), header = TRUE)[,1]
patient_set_lt20pamp <- intersect(read.table(paste0(main_path, "Patient Subsets/LessThan20PercAmp_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumA <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.A_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumB <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.B_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumAB <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.A.B_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_basal <- intersect(read.table(paste0(main_path, "Patient Subsets/Basal_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_her2 <- intersect(read.table(paste0(main_path, "Patient Subsets/HER2_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_normLike <- intersect(read.table(paste0(main_path, "Patient Subsets/Normal-like_patient_ids.txt"), header = TRUE)[,1], patient_set)

expression_df_pats <- unlist(lapply(colnames(expression_df), function(x)
  unlist(strsplit(x, "-", fixed = TRUE))[1]))
expression_df_lumA <- expression_df[, which(expression_df_pats %in% patient_set_lumA)]
expression_df_lumB <- expression_df[, which(expression_df_pats %in% patient_set_lumB)]
expression_df_basal <- expression_df[, which(expression_df_pats %in% patient_set_basal)]
expression_df_her2 <- expression_df[, which(expression_df_pats %in% patient_set_her2)]