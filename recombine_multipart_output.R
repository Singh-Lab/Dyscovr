#!/usr/bin/env Rscript

# Set the library path
library.path <- .libPaths()

library(qvalue, lib.loc = library.path, quietly = TRUE)
library(data.table, lib.loc = library.path, quietly = TRUE)
library(argparse, lib.loc = library.path, quietly = TRUE)

# This function takes multipart Dyscovr output and recombines it into one,
# redoing the per-gene multiple hypothesis testing correction and ordering by 
# q-value.

# Create parser object
parser <- ArgumentParser()

# Specify desired options for the parser

# Path to Dyscovr output files
parser$add_argument("--path", default = "", type = "character", 
                    help = "Set the path to the output files to be merged.")

# Input file name for the first part of the recombined Dyscovr output
parser$add_argument("--fn", default = "", type = "character", 
                    help = "The name of the first part of the output (e.g. the filename of pt1")

# Parse the given arguments
tryCatch({
  args <- parser$parse_args()
}, error = function(cond) {
  print(cond)
  print(traceback())
})


############################################################
#### IMPORT MULTIPART DYSCOVR OUTPUT AND RECOMBINE
############################################################
# List all files that match the pattern given in the first file
fn_split <- unlist(strsplit(args$fn, "pt1", fixed = T))
fns_all <- intersect(list.files(args$path, pattern = fn_split[1]),
                   list.files(args$path, pattern = fn_split[2]))
files <- lapply(fns_all, function(f) fread(paste0(args$path, f)))

# Bind together the output files
recombined_output <- do.call(rbind, files)

############################################################
#### PEFORM MULTIPLE HYPOTHESIS TESTING CORRECTION
############################################################
#' Function takes in an output results table and applies multiple hypothesis 
#' testing correction (Storey's q-value correction) to all p-values in order to 
#' add a column of q-values. Returns the results table with a column for q-values.
#' @param results_table a combined output DF produced from Dyscovr
#' @param per_gene a T/F value indicating whether or not we are doing the 
#' correction per-driver (T) or together across all drivers (F)
mh_correct <- function(results_table, per_gene) {
  
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
      
      qvals <- qobj$qvalues # extract qvalues
      
      results_table$q.value <- qvals # add the qvalues back to the data frame
    }
  }
  # Return the tidied linear model fit with q-values
  return(results_table)
}

# Call function
recombined_output_corrected <- mh_correct(recombined_output, correct_per_gene = T)

# Reorder by new q-value
recombined_output_corrected <- recombined_output_corrected[order(recombined_output_corrected$q.value),]

# Write to file
new_outfn <- paste(fn_split[1], fn_split[2])
fwrite(recombined_output_corrected, paste0(args$path, new_outfn))