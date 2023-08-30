############################################################
### Collinearity Diagnostics, Model Fit, and Variable Contribution
### Written By: Sara Geraghty, September 2021
############################################################

# Uses the olsrr package and others to interrogate collinearity diagnostics,
# model fit, and variable contribution of linear model. The linear regression
# is passed in as an lm object along with the path for the visualizations and 
# files to be written to.

library(olsrr)
library(stringr)

# Link to olsrr vignette: https://cran.r-project.org/web/packages/olsrr/vignettes/regression_diagnostics.html


############################################################
### GET VIF, TOLERANCE, & CONDITION INDEX
############################################################
#' A function that, given a linear model fit and information about the fit, 
#' creates a VIF table and an eigenvalue condition index table and writes
#' them both to the given path.
#' @param fit_lm the linear model fit, in the form of an lm object produced from lm()
#' @param path the path for the output files
#' @param fn the generic filename to modify for the output files
#' @param name_targ the name of the target gene from the fit
get_diag_tables <- function(fit_lm, path, fn, name_targ) {
  
  ols_diagn_tabs <- ols_coll_diag(fit_lm)
  vif_tab <- ols_diagn_tabs[[1]]
  eig_cindex_tab <- ols_diagn_tabs[[2]]
  
  fn <- str_replace(fn, "output_results", name_targ[1])
  fn <- paste(path, fn, sep = "/")
  fwrite(vif_tab, paste(fn, "vif_tab.csv", sep = "_"))
  fwrite(eig_cindex_tab, paste(fn, "eig_cindex_tab.csv", sep = "_"))
}

get_diag_tables(fit_lm, outpath, outfn, targ)


############################################################
### MODEL FIT ASSESSMENT
############################################################
#' A function that, given a linear model fit and information about the fit, 
#' runs various model fit assessments, including creating a residual spread plot,
#' correlations, and an observed vs. predicted plot
#' @param fit_lm the linear model fit, in the form of an lm object produced from lm()
#' @param path the path for the output files
#' @param fn the generic filename to modify for the output visualizations
#' @param name_targ the name of the target gene from the fit
assess_model_fit <- function(fit_lm, path, fn, name_targ) {
  outpath_vis <- str_replace(outpath, "output_files", "output_visualizations")
  fn <- str_replace(fn, "output_results", name_targ)
  fn <- paste(path, fn, sep = "/")
  
  # Residual fit spread plot
  pdf(file = paste(fn, "resid_fit_spread.pdf", sep = "_"))
  ols_plot_resid_fit(fit_lm)
  dev.off()
  
  # Correlations 
  #corr <- ols_correlations(fit_lm)
  #fwrite(corr, paste(fn, "correlations.csv", sep = "_"))
  
  # Observed vs. fitted values (to assess fit of model)
  pdf(file = paste(fn, "observ_fit.pdf", sep = "_"))
  ols_plot_obs_fit(fit_lm)
  dev.off()
}

assess_model_fit(fit_lm, outpath, outfn, targ)


############################################################
### VARIABLE CONTRIBUTIONS
############################################################
#' A function that, given a linear model fit and information about the fit, 
#' runs various variable contribution assessments, including creating an added 
#' variable plot and residual plus component plot
#' @param fit_lm the linear model fit, in the form of an lm object produced from lm()
#' @param path the path for the output files
#' @param name_targ the name of the target gene from the fit
assess_variable_contributions <- function(fit_lm, path, fn, name_targ) {
  fn <- str_replace(fn, "output_results", name_targ)
  path <- str_replace(path, "output_files", "output_visualizations")
  fn <- paste(path, fn, sep = "/")
  
  # Variable contribution plots
  pdf(file = paste(fn, "added_variable.pdf", sep = "_"))
  ols_plot_added_variable(fit_lm)
  dev.off()
  
  # Residual component plots
  pdf(file = paste(fn, "residual_components.pdf", sep = "_"))
  ols_plot_comp_plus_resid(fit_lm)
  dev.off()
}

#assess_variable_contributions(fit_lm, tmp_path, regprot, targ)

