############################################################
### PERFORM META-ANALYSIS
### Written By: Sara Geraghty, July 2020
############################################################

# Functions to perform a meta-analysis on the results from
# linear_model.R and combine Betas from runs in different cancer
# types or subtypes

# install.packages("meta", "metasens")
library(meta)
library(metasens)
library(metafor)
library(ggplot2)


# Path to output files
main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

# Path to where output figures should be saved
output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/BRCA/Linear Model/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/Pan-Cancer/Linear Model/"


# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


############################################################
############################################################
#### PERFORM A META-ANALYSIS TO AGGREGATE SIGNAL
############################################################
############################################################
# Use the package "meta" for perform a meta-analysis, treating the various subtype
# results as results from separate studies
# Link to documentation: https://cran.r-project.org/web/packages/meta/meta.pdf

#' Create an input data frame of the type expected by meta
#' @param list_of_master_dfs a list of the master DFs with results we are interested
#' in combining in our meta-analysis. These are already subsetted to only have results
#' for one regulatory protein of interest. List names are the names of the subgroups.
#' @param vect_of_n a vector the same length as the above list, with the number of 
#' samples for each subgroup (n)
#' @param targ the external gene name of the target gene of interest
create_meta_input_tab <- function(list_of_master_dfs, vect_of_n, targ) {
  meta_input_tab <- data.frame("study" = names(list_of_master_dfs), 
                               "n_patients" = vect_of_n)
  # Get the Beta estimate and t-statistic value for each of the subgroups, for
  # the given target
  meta_input_tab$beta <- lapply(list_of_master_dfs, function(x) 
    x[x$T_k.name == targ, 'estimate']) 
  meta_input_tab$tstat <- lapply(list_of_master_dfs, function(x) 
    x[x$T_k.name == targ, 'statistic']) 
  meta_input_tab$std.error <- lapply(list_of_master_dfs, function(x) 
    x[x$T_k.name == targ, 'std.error']) 
  
  return(meta_input_tab)
}

master_df_list_p53 <- list("LumA"= lumA_allgenes_mut_p53, "LumB" = lumB_allgenes_mut_p53,
                           "Basal" = basal_allgenes_mut_p53, "HER2" = her2_allgenes_mut_p53)

patient_set <- read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt"), header = TRUE)[,1]
patient_set_lt20pamp <- intersect(read.table(paste0(main_path, "Patient Subsets/LessThan20PercAmp_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumA <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.A_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_lumB <- intersect(read.table(paste0(main_path, "Patient Subsets/Luminal.B_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_basal <- intersect(read.table(paste0(main_path, "Patient Subsets/Basal_patient_ids.txt"), header = TRUE)[,1], patient_set)
patient_set_her2 <- intersect(read.table(paste0(main_path, "Patient Subsets/HER2_patient_ids.txt"), header = TRUE)[,1], patient_set)

patient_set_pc <- read.table(paste0(main_path, "Linear Model/Tumor_Only/intersecting_ids.txt"), header = TRUE)[,1]

blca_subtype <- TCGAquery_subtype(tumor = "BLCA")
hnsc_subtype <- TCGAquery_subtype(tumor = "HNSC")
blca_subtype$justPat <- unlist(lapply(blca_subtype$patient, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))
hnsc_subtype$patient <- as.character(hnsc_subtype$patient)
hnsc_subtype$justPat <- unlist(lapply(hnsc_subtype$patient, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))

patient_set_blca <- intersect(patient_set_pc, blca_subtype$justPat)
patient_set_hnsc <- intersect(patient_set_pc, hnsc_subtype$justPat)


n_of_subtypes <- c(length(patient_set_lumA), length(patient_set_lumB), length(patient_set_basal), length(patient_set_her2))
n_of_subtypes <- c(length(patient_set), length(patient_set_blca), length(patient_set_hnsc))

meta_input_tab_p53 <- create_meta_input_tab(master_df_list_p53, n_of_subtypes, "DDB2")
meta_input_tab_p53 <- create_meta_input_tab(master_df_list_p53, n_of_subtypes, "SRD5A1")


#' Perform the meta-analysis using a random-effects model
#' @param meta_input_tab an input table for meta-analysis as created by helper function
#' @param targ the external gene name of the target gene of interest
perform_meta_analysis <- function(meta_input_tab, targ) {
  rma_model <- rma(yi = as.numeric(unlist(meta_input_tab$tstat)), 
                   sei = as.numeric(unlist(meta_input_tab$std.error)))
  print(summary(rma_model))
  
  # Make a forest plot 
  forest(rma_model, slab = meta_input_tab$study)
  
  # Add details
  text(-10, -1.5, pos=4, cex=0.75, bquote(paste("RE Model (Q = ", 
                                                .(formatC(rma_model$QE, digits=2, format="f")), ", df = ", .(rma_model$k - rma_model$p),
                                                ", p = ", .(formatC(ma_model$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                                .(formatC(rma_model$I2, digits=1, format="f")), "%)")))
  # Add title
  grid.text(paste(targ, "Meta-Analysis"), .5, .8, gp=gpar(cex=1.5))
  
  return(rma_model)
}

perform_meta_analysis(meta_input_tab_p53, "DDB2")
perform_meta_analysis(meta_input_tab_p53, "SRD5A1")


#' Perform meta-analysis across a large gene set, getting a p-value for each gene of interest
#' based on the given cancer types or subtypes.
#' @param list_of_master_dfs a list of the master DFs with results we are interested
#' in combining in our meta-analysis. These are already subsetted to only have results
#' for one regulatory protein of interest. List names are the names of the subgroups.
#' @param vect_of_n a vector the same length as the above list, with the number of 
#' samples for each subgroup (n)
#' @param vect_of_targets a vector of gene targets we are interested in doing a 
#' meta-analysis on
perform_meta_analysis_across_genes <- function(list_of_master_dfs, vect_of_n, 
                                               vect_of_targets) {
  output_df <- data.frame("gene.name" = vect_of_targets)
  
  out <- lapply(vect_of_targets, function(targ) {
    input_tab <- create_meta_input_tab(list_of_master_dfs, vect_of_n, targ)
    res <- perform_meta_analysis(input_tab, targ)
    pval <- res$pval
    estimate <- res$beta
    return(data.frame("p.value" = pval, "estimate" = estimate))
  })
  out_df <- do.call(rbind, out)
  output_df <- cbind(output_df, out_df)
  
  return(output_df)
}

target_gene_set <- intersect(lumA_metabol_p53$T_k.name, intersect(lumB_metabol_p53$T_k.name,
                                                                  intersect(basal_metabol_p53$T_k.name, her2_metabol_p53$T_k.name)))
target_gene_set <- intersect(brca_metabol_p53$T_k.name, intersect(blca_metabol_p53$T_k.name,
                                                                  hnsc_metabol_p53$T_k.name))

brca_subtypes_meta_analysis_res <- perform_meta_analysis_across_genes(master_df_list_p53, n_of_subtypes,
                                                                      target_gene_set)
cancerTypes_meta_analysis_res <- perform_meta_analysis_across_genes(master_df_list_p53, n_of_subtypes,
                                                                    target_gene_set)

# Order by pval
brca_subtypes_meta_analysis_res <- brca_subtypes_meta_analysis_res[order(brca_subtypes_meta_analysis_res$p.value),]
cancerTypes_meta_analysis_res <- cancerTypes_meta_analysis_res[order(cancerTypes_meta_analysis_res$p.value),]

# Add q-values
brca_subtypes_meta_analysis_res$q.value <- qvalue(brca_subtypes_meta_analysis_res$p.value)$qvalue
cancerTypes_meta_analysis_res$q.value <- qvalue(cancerTypes_meta_analysis_res$p.value)$qvalue

# Visualize pval distribution
hist(brca_subtypes_meta_analysis_res$p.value, main = "", xlab = "P-Value")
hist(cancerTypes_meta_analysis_res$p.value, main = "", xlab = "P-Value")

# Write to a file
write.csv(brca_subtypes_meta_analysis_res, paste0(main_path, "Linear Model/TP53/Tumor_Only/eQTL/Meta-Analysis Results/brca_subtypes_metabolism_meta_analysis.csv"))
write.csv(cancerTypes_meta_analysis_res, paste0(main_path, "Linear Model/TP53/Tumor_Only/eQTL/Meta-Analysis Results/brca_blca_hnsc_metabolism_meta_analysis.csv"))

