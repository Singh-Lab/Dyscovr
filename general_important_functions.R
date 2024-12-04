############################################################
### General Important Functions
### Written by Sara Geraghty, Princeton University
### https://www.biorxiv.org/content/10.1101/2024.11.20.624509v1
############################################################

# This file contains functions that are helpful across files
# and will be used the project

############################################################
#' Use a faster %in% function
#library(fastmatch, lib.loc = library.path)
`%fin%` <- function(x, table) {
  #stopifnot(require(fastmatch))
  fmatch(x, table, nomatch = 0L) > 0L
}
############################################################

############################################################
#' Add a string to bool function which converts a string "boolean" into a boolean 
#' type. Modified from: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
#' @param v a argparse input value that should be a logical
str2bool <- function(v) {
  if (is.logical(v)) {return(v)}
  if (tolower(v) %fin% c('yes', 'T', 't', 'y', '1', 'true')) {return(T)}
  else if (tolower(v) %fin% c('no', 'F', 'f', 'n', '0', 'false')) {return(F)}
  else {
    print("Error. Boolean value expected.")
    return(9)
  }
}
############################################################

