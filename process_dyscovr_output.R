############################################################
### Process Dyscovr Output
### PUBLICATION INFORMATION 
############################################################

# Given an output 'master file' from Dyscovr, this file performs the following 
# functions:
  # Basic Output Visualization (Beta Distribution, P-Value Distribution,
  # Q-Q Plot, and SE Distribution)
  # Multiple hypothesis testing correction
  # Addition of interpretable gene names

#!/usr/bin/env Rscript

# install.packages("gplots")

#library("gplots", lib.loc = library.path, quietly = TRUE)
#library(broom, lib.loc = library.path, quietly = TRUE)
#library(qvalue, lib.loc = library.path, quietly = TRUE)
#library(ggplot2, lib.loc = library.path, quietly = TRUE)
#library(stringr, lib.loc = library.path, quietly = TRUE)
#library(dplyr, lib.loc = library.path, quietly = TRUE)
#library(olsrr, lib.loc = library.path, quietly = TRUE)
#library("RColorBrewer", lib.loc = library.path, quietly = TRUE)


############################################################
# SET PATH
############################################################
# Path to where output figures should be saved
OUTPUT_VIS_PATH <- "/Genomics/argo/users/scamilli/Dyscovr/output_visualizations/"
if(args$dataset != "TCGA") {OUTPUT_VIS_PATH <- paste0(OUTPUT_VIS_PATH, 
                                                      paste0(args$dataset, "/"))}
OUTPUT_VIS_PATH <- paste0(OUTPUT_VIS_PATH, args$cancerType)

if(test) {
  OUTPUT_VIS_PATH <- paste(OUTPUT_VIS_PATH, args$tester_name, sep = "/")
} else {
  OUTPUT_VIS_PATH <- paste(OUTPUT_VIS_PATH, args$run_name, sep = "/")
}
if ((args$cancerType == "PanCancer") & (!(args$specificTypes == "ALL"))) {
  OUTPUT_VIS_PATH <- paste(OUTPUT_VIS_PATH, ct, sep = "/")
}
print(paste("Output Visualization Path:", OUTPUT_VIS_PATH)) 

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(SOURCE_PATH, "/input_files/all_genes_id_conv.csv"), 
			header = T)

# Adjust the outfn for visualizations
outfn_vis <- paste(unlist(strsplit(outfn, "_", fixed = TRUE))[2:length(
  unlist(strsplit(outfn, "_", fixed = TRUE)))], collapse = "_")


############################################################
#### PEFORM MULTIPLE HYPOTHESIS TESTING CORRECTION
############################################################
#' Function takes in an output results table and applies multiple hypothesis 
#' testing correction (Storey's q-value correction) to all p-values in order to 
#' add a column of q-values. Returns the results table with a column for q-values.
#' @param results_table a combined output DF produced from Dyscovr
#' @param per_gene a T/F value indicating whether or not we are doing the 
#' correction per-driver (T) or together across all drivers (F)
#' @param fn_qvalvis a file name with path for a q-value visualization
#' produced from the qvalue package
#' @param fn_qvalsum a file name with path for a q-value summary 
#' with useful information produced from the qvalue package
mh_correct <- function(results_table, per_gene, fn_qvalvis, fn_qvalsum) {
  
  # Get the qvalue object
  qobj <- NA
  if(length(results_table$p.value) < 100) {
    # With a small number of pvalues we may not be able to accurately estimate 
    # pi0, so we set to 1 (the equivalent of B-H correction)
    qobj <- qvalue(p = results_table$p.value, pi0 = 1)
    if(per_gene) {
      print("Too few p-values to do the correction per-driver.")
    }
    qvals <- qobj$qvalues # extract qvalues
    
    results_table$q.value <- qvals # add the qvalues back to the data frame
    
  } else {
    if(per_gene) {
      unique_drivers <- unique(results_table$R_i)
      list_of_corrected_tabs <- lapply(unique_drivers, function(d) {
        res_tab_sub <- results_table[results_table$R_i == d, ]
        tryCatch({
          qobj <- qvalue(p = res_tab_sub$p.value)
        }, error = function(cond) {
          if(grepl("The estimated pi0 <= 0", as.character(cond))) {
            qobj <- qvalue(p = res_tab_sub$p.value, lambda = 0.05)
          } else {qobj <- NA}
        })
        qvals <- qobj$qvalues # extract qvalues
        res_tab_sub$q.value <- qvals # add the qvalues back to the data frame
        return(res_tab_sub)
      })
      results_table <- do.call(rbind, list_of_corrected_tabs)
      results_table <- results_table[order(results_table$q.value),]
      
    } else {
      tryCatch({
        qobj <- qvalue(p = results_table$p.value)
      }, error = function(cond) {
        if(grepl("The estimated pi0 <= 0", as.character(cond))) {
          qobj <- qvalue(p = results_table$p.value, lambda = 0.05)
        } else {qobj <- NA}
      })
      
      # Plot some useful plots & print some useful information
      png(fn_qvalvis, width = 450, height = 450, type="cairo")
      plot(qobj)
      dev.off()
      
      png(fn_qvalsum, width = 450, height = 450, type="cairo")
      print(summary(qobj))
      dev.off()
      #print(paste("Pi0 (Propr. of true null hypotheses):", qobj$pi0))
      
      qvals <- qobj$qvalues # extract qvalues
      
      results_table$q.value <- qvals # add the qvalues back to the data frame
    }
  }
  # Return the tidied linear model fit with q-values
  return(results_table)
}

# Call this function
if("p.value" %fin% colnames(master_df_mut)) {
  fn_mut_qvalvis <- paste(OUTPUT_VIS_PATH, 
                          paste0("Q-ValueVisualiz_", 
                                 paste0(paste0(outfn_vis, "_MUT"), ").png")), 
                          sep = "/")
  fn_mut_qvalsum <- paste(OUTPUT_VIS_PATH, 
                          paste0("Q-ValueSummary_", 
                                 paste0(paste0(outfn_vis, "_MUT"), ").png")), 
                          sep = "/")

  tryCatch({
    master_df_mut_corrected <- mh_correct(master_df_mut, correct_per_gene, 
                                          fn_mut_qvalvis, fn_mut_qvalsum)
  }, error = function(cond) {print(cond)})
  
} else {
  master_df_mut_corrected <- master_df_mut
}

############################################################
#### ADD GENE NAMES 
############################################################
#' Given a master data frame result from Dyscovr, add a column for target gene 
#' name and driver gene name. Return the updated data frame.
#' @param master_df a data.table object produced from the Dyscovr model
#' @param all_genes_id_conv a bioMart file with conversions between different 
#' gene ID types
#' @param runQueriesJointly  a TRUE/ FALSE value to indicate if we are running 
#' all driver genes jointly in the same model
add_targ_regprot_gns <- function(master_df, all_genes_id_conv, runQueriesJointly ) {
  master_df$T_k.name <- unlist(lapply(master_df$T_k, function(x) 
    paste(unique(unlist(all_genes_id_conv[
      all_genes_id_conv$uniprot_gn_id == unlist(strsplit(x, ";", fixed = TRUE))[1], 
      'external_gene_name'])), collapse = ";")))
  
  # Add a column for the driver gene name
  unique_ri <- data.frame("R_i" = unique(master_df$R_i))
  unique_ri$R_i.name <- unlist(lapply(unique_ri$R_i, function(x) 
    paste(unique(all_genes_id_conv[
      all_genes_id_conv$uniprot_gn_id == unlist(strsplit(x, ";", fixed = TRUE))[1], 
      'external_gene_name']), collapse = ";")))
  master_df$R_i.name <- unlist(lapply(master_df$R_i, function(x) 
    unique_ri[unique_ri$R_i == x, 'R_i.name']))
  
  return(master_df)
}

tryCatch({
  # Call this function
  print(paste("Run driver jointly:", runQueriesJointly))
  master_df_mut_corrected <- add_targ_regprot_gns(master_df_mut_corrected, 
                                                  all_genes_id_conv, 
                                                  runQueriesJointly)

  # Write this to a new file
  outfn <- str_replace(outfn, "uncorrected", "corrected") 
  print(paste("NEW OUTPUT FILENAME:", outfn))
  #print(paste(outpath_curr, paste(outfn, paste("_MUT", ".csv", sep = ""), 
  #sep = ""), sep = "/"))
  fwrite(master_df_mut_corrected, paste(outpath_curr, 
                                        paste0(outfn, 
                                               paste0("_MUT", ".csv")), 
                                        sep = "/"))
}, error = function(cond) {print(cond)})


############################################################
############################################################
#### BASIC VISUALIZTION OF OUTPUT
############################################################
############################################################

############################################################
#### VISUALIZE BETA VALUE DISTRIBUTION
############################################################
#' Function plots a histogram to visualize the distribution
#' of driver protein R_i beta values across all models.
#' @param results_table a master DF produced from Dyscovr
visualize_beta_distrib <- function(results_table) {
  betas <- results_table$estimate
  hist(betas, main = "Histogram of Mutation Coefficient Values Across all Drivers",
       xlab = "Coefficient Value", ylab = "Frequency", col = "darkseagreen2")
}

# Create specific filenames
fn_mut <- paste(OUTPUT_VIS_PATH, 
                paste0("BetaDistrib_", 
                       paste0(paste0(outfn_vis, "_MUT"), ").png")), sep = "/")
# print(paste("FN MUT:", fn_mut))

# Call this function & save output
tryCatch({
  png(fn_mut, width = 450, height = 350, type="cairo")
  visualize_beta_distrib(master_df_mut)
  dev.off()
}, error = function(cond) {print(cond)})


############################################################
#### VISUALIZE P-VALUE DISTRIBUTION
############################################################
#' Function plots a histogram to visualize the distribution
#' of combined driver p-values across all models.
#' @param results_table a master DF produced from Dyscovr
visualize_pval_distrib <- function(results_table) {
  pvals <- results_table$p.value[!is.na(results_table$p.value) & 
                                   !is.infinite(results_table$p.value)]
  hist(pvals, main = "Histogram of p-Values",
       xlab = "p-value", ylab = "Frequency", col = "blueviolet")
}

# Create specific filenames
if("p.value" %fin% colnames(master_df_mut)) {
  fn_mut <- paste(OUTPUT_VIS_PATH, 
                  paste0("P-ValDistrib_", 
                         paste0(paste0(outfn_vis, "_MUT"), ").png")), sep = "/")

  # Call this function & save output
  tryCatch({
    png(fn_mut, width = 450, height = 350, type="cairo")
    visualize_pval_distrib(master_df_mut)
    dev.off()
  }, error = function(cond) {print(cond)})
}


############################################################
#### VISUALIZE Q-Q PLOT
############################################################
#' Function plots a Q-Q plot to visualize the distribution of p-values
#' and assess whether they come from a uniform distribution
#' @param results_table a master DF produced from Dyscovr
qqplot_pvals <- function(results_table) {
  qqnorm(results_table$p.value, pch = 1, frame = FALSE)
  qqline(results_table$p.value, col = "steelblue", lwd = 2)
}
# NOTE: currently not called 


############################################################
#### VISUALIZE ERROR DISTRIBUTION
############################################################
#' Function plots a histogram of the standard errors produced from
#' all LM runs to assess whether the errors derive from a normal distribution
#' @param results_table a master DF produced from Dyscovr
visualize_error_distrib <- function(results_table) {
  hist(results_table$std.error, 
       main = "Standard Error Distribution Across All Tests",
       xlab = "Standard Error (SE)", ylab = "Frequency")
}
# NOTE: currently not called 