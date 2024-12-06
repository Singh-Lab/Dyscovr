############################################################
### Run Regularized Models
### Written by Sara Geraghty, Princeton University
### https://www.biorxiv.org/content/10.1101/2024.11.20.624509v1
############################################################

# This file contains a set of functions that will allow the Dyscovr model to 
# perform an  L1 (LASSO) or L2 (Ridge) regularization, evaluate model fit, and 
# calculate associated significance metrics.

library.path <- .libPaths()

#library(glmnet, lib.loc = library.path)
#library(gglasso, lib.loc = library.path)
#library(stringr, lib.loc = library.path)
#library(dqrng, lib.loc = library.path)


############################################################
#' If we are doing a regularized version of our model, run the details here
#' using functions from glmnet/gglasso packages
#' @param formula the string version of the formula we would like to use for our
#' regression model
#' @param lm_input_table the input table with values corresponding to all 
#' variables in the formula
#' @param type a string indicating the type of regularization (ridge, lasso, or 
#' bayesian)
#' @param debug a T/ F value indicating whether or not we are in debug 
#' mode
#' @param meth_bucketing a T/F value indicating whether or not methylation 
#' was bucketed; only used in the case of group lasso
#' @param cna_bucketing a T/F value indicating whether or not CNA was
#' bucketed; only used in the case of group lasso
#' @param signif_eval_type a character string giving the type of method we will 
#' use to evaluate significance of our Betas. Options include, "randomization",
#' "subsampling", or "selectiveInference"
#' @param shuffled_input_dfs a list of shuffled input DFs for doing the 
#' randomization_predictor method of p-value generation, or shuffled expression 
#' DFs for doing randomization_perSamp or randomization_perTarg
#' @param ensg the ENSG ID of the current target gene, such that we can subset 
#' the DF appropriately
#' @param num_randomizations if we are doing randomizations per target or per 
#' sample, this provides the number of randomizations to perform
#' @param fixed_lambda_for_rand a T/F value indicating whether or not we 
#' are using a fixed value of lambda for the randomizations (e.g. from the "real" 
#' lasso) or if we are using cross-validation to determine a lambda value
run_regularization_model <- function(formula, lm_input_table, type, debug, 
                                     meth_bucketing, cna_bucketing, 
                                     signif_eval_type, shuffled_input_dfs, ensg, 
                                     num_randomizations, fixed_lambda_for_rand) {
  
  lm_input_table <- as.data.frame(na.omit(lm_input_table))
  
  # Get the data of interest using the formula
  spl_formula <- unlist(strsplit(formula, " ~ ", fixed = T))
  y_var <- trimws(spl_formula[1], which = "both")
  y_data <- as.matrix(lm_input_table[,colnames(lm_input_table) == y_var])
  x_vars <- unlist(strsplit(trimws(spl_formula[2], which = "both"), " + ", 
                            fixed = T))
  x_data <- as.matrix(lm_input_table[,colnames(lm_input_table) %fin% x_vars])
  
  # Use a cross validation glmnet to get the best lambda value; alpha = 0 
  # indicates that we are doing L2-regularization
  #lambdas <- 10^seq(2, -3, by = -.1)
  lambdas <- 10^seq(2, -3, by = -.2)
  best_model <- NA
  optimal_lambda <- NA
  
  # Ridge regression
  if((type == "ridge") | (type == "L2")) {
    tryCatch({
      ridge.cv <- cv.glmnet(x_data, y_data, alpha = 0, lambda = lambdas)
      optimal_lambda <- ridge.cv$lambda.min
      #if(debug) {print(paste("Optimal Lambda:", optimal_lambda))}
      # Run the model with the best lambda
      best_model <- glmnet(x_data, y_data, alpha = 0, family = "gaussian", 
                           lambda = optimal_lambda)
    }, error = function(cond) {best_model <- NA})
    
    # Lasso
  } else if((type == "lasso") | (type == "L1")) {
    tryCatch({
      
      # Scale and center the x and y data
      x_data <- scale.default(x_data, center = T, scale = T)
      y_data <- y_data - mean(y_data)
      
      # Use a group lasso using the gglasso package, to account
      # for the dummy variables (bucketed variables)
      v.group <- get_groups_for_variables(x_vars, meth_bucketing, cna_bucketing, 
                                          "label_per_var")
      
      # Get the optimal lambda using cross-validation
      #lasso.cv <- cv.glmnet(x_data, y_data, alpha = 1, lambda = lambdas, 
      #standardize = T, nfolds = 5)
      # Default number of folds is 5
      lasso.cv <- cv.gglasso(x_data, y_data, group = v.group, lambda = lambdas, 
                             loss = "ls")
      optimal_lambda <- lasso.cv$lambda.min
      print(paste("Optimal Lambda:", optimal_lambda))
      
      # Record this value
      #write.table(optimal_lambda, file = "optimal_lambdas_real.txt", 
      #           row.names = F, append = T)
      
      # Run the model with the best lambda
      #best_model <- glmnet(x_data, y_data, alpha = 1, family = "gaussian",
      #lambda = optimal_lambda, standardize = F, intercept = F)
      best_model <- gglasso(x_data, y_data, lambda = optimal_lambda, 
                            group = v.group, loss = "ls", intercept = F)
      
      # Get a p-value or some form of significance measure for these results,
      # as specified in the argument
      if(signif_eval_type == "selectiveInference") {
        best_model <- eval_signif_selectiveInference(best_model, x_data, y_data,
                                                     optimal_lambda)
        
      } else if(grepl("randomization", signif_eval_type)) {
        if(fixed_lambda_for_rand) {
          best_model <- eval_signif_randomization(
            best_model = best_model, sample_ids = lm_input_table$sample_id, 
            x_data = x_data, x_vars = x_vars, y_data = y_data, v.group = v.group, 
            lambdas = optimal_lambda, debug = debug, 
            shuffled_input_dfs = shuffled_input_dfs, ensg = ensg, 
            num_randomizations = num_randomizations, 
            signif_eval_type = signif_eval_type)
          
        } else {
          # Use only ~half the number of possible lambdas as for the real LASSO, 
          # for improved runtime
          lambdas <- sort(lambdas)
          lambdas_sub <- c(lambdas[c(T, F)], 100)
          
          best_model <- eval_signif_randomization(
            best_model = best_model, sample_ids = lm_input_table$sample_id, 
            x_data = x_data, x_vars = x_vars, y_data = y_data, v.group = v.group, 
            lambdas = lambdas_sub, debug = debug, 
            shuffled_input_dfs = shuffled_input_dfs, ensg = ensg, 
            num_randomizations = num_randomizations, 
            signif_eval_type = signif_eval_type)
        }
        
        
      } else if(signif_eval_type == "subsampling") {
        # Use only ~half the number of possible lambdas as for the real LASSO, 
        # for improved runtime
        lambdas <- sort(lambdas)
        lambdas_sub <- c(lambdas[c(T, F)], 100)
        
        best_model <- eval_signif_subsampling(best_model, x_data, y_data, 
                                              v.group, lambdas_sub, debug)
        
      } else {
        print("Error: only selectiveInference, randomization, and subsampling are 
        valid inputs for type of significance evaluation. Please try again.")
        best_model <- NA
      }
      
    }, error = function(cond) {
      print(paste("error:", cond))
      best_model <- NA
    })
    
    
    # Bayesian Lasso
  } else if (type %fin% c("bayesian.bgl", "bayesian.bglss")) {
    
    x_data <- scale.default(x_data, center = T, scale = T)
    
    # Y data needs to be a vector
    y_data <- as.numeric(unlist(y_data))
    y_data <- y_data - mean(y_data)
    
    # Use a group lasso using the gglasso package, to account
    # for the dummy variables (bucketed variables)
    v.group <- NA
    if(type == "bayesian.bglss") {
      v.group <- get_groups_for_variables(x_vars, meth_bucketing, cna_bucketing, 
                                          "group_sizes")
    } else {
      v.group <- get_groups_for_variables(x_vars, meth_bucketing, cna_bucketing, 
                                          "label_per_var")
    }

    # Call external Bayesian functions from run_bayesian_lasso.R
    tryCatch({
      if(type == "bayesian.bgl") {
        bgl_res <- runBGL(x_data, y_data, groups = v.group, scale = F)
        best_model <- data.table("term" = colnames(x_data), 
                                 "estimate" = bgl_res$Sparse_Beta, 
                                 "p.value" = bgl_res$Confidence)
      } else {
        bglss_res <- runBGLSS(x_data, y_data, group_sizes = v.group, 
                              scale = F)
        best_model <- data.table("term" = colnames(x_data), 
                                 "estimate" = bglss_res$pos_median, 
                                 "p.value" = bglss_res$Confidence)
      }
    }, error = function(cond) {
      print(paste("error:", cond))
      best_model <- NA
    }
  )
    
  } else {
    print(paste("Only implemented for ridge, lasso, bayesian.bgl and bayesian.bglss;", 
                paste(type, "is not implemented. Proceeding with ridge.")))
    tryCatch({
      ridge.cv <- cv.glmnet(x_data, y_data, alpha = 0, lambda = lambdas)
      optimal_lambda <- ridge.cv$lambda.min
      if(debug) {print(paste("Optimal Lambda:", optimal_lambda))}
    }, error = function(cond) {best_model <- NA})
    
    # Run the model with the best lambda
    best_model <- glmnet(x_data, y_data, alpha = 0, lambda = optimal_lambda)
  }
  
  # Option to return an input data frame that will then be 
  # put back into a multiple linear regression model with LASSO-selected x variables
  #input_df <- run_regularization_output_in_lm(best_model, x_data, y_data)
  #return(input_df)
  
  if(debug) {print(head(best_model))}
  
  return(best_model)
}

# NOTE: getting an object not found error for formula_new; adjusted according to 
# this post: https://stackoverflow.com/questions/8218196/object-not-found-error-when-passing-model-formula-to-another-function


#' Evaluation metrics function from plural sight guide (R-squared and RMSE)
#' https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r
#' @param model a regression model
#' @param df the input DF for the regression
#' @param predications the results from stat's "predict" function for the model 
#' & input DF
#' @param target the target or y variable
eval_metrics <- function(model, df, predictions, target){
  resids = df[,target] - predictions
  resids2 = resids**2
  N = length(predictions)
  r2 = as.character(round(summary(model)$r.squared, 2))
  adj_r2 = as.character(round(summary(model)$adj.r.squared, 2))
  print(paste("R-squared:", adj_r2)) #Adjusted R-squared
  print(paste("RMSE:", as.character(round(sqrt(sum(resids2)/N), 2)))) #RMSE
}


#' Helper function to get the categorical variables that we have bucketed (e.g. 
#' one-hot encoded) for group lasso, and assign them to groups. Returns a vector 
#' of group assignments for all variables in x_data.
#' @param x_vars a vector of the names of the x variables (covariates)
#' @param meth_bucketing a T/F value indicating whether or not methylation 
#' was bucketed
#' @param cna_bucketing a T/F value indicating whether or not CNA was
#' bucketed
#' @param type_of_grouping either "label_per_var" or "group_sizes"; label_per_var 
#' indicates that we return a vector of the length of X_vars, with a number for 
#' each position indicating that variables group membership; group_sizes 
#' indicates that we give a vector of the size of each group, assuming that 
#' group members are contiguous in the X data input
get_groups_for_variables <- function(x_vars, meth_bucketing, cna_bucketing, 
                                     type_of_grouping) {
  
  # Obtain the bucketed, categorical variables
  bucketed_vars <- c("Tot_Mut_b","Tumor_purity_b", "Tot_IC_Frac_b", "Cancer_type_b")
  #if(meth_bucketing) {
  #  bucketed_vars <- c(bucketed_vars, "MethStat_k")
  #  if(analysis_type == "eQTL") {bucketed_vars <- c(bucketed_vars, "MethStat_i")}
  #}
  #if(!((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing)))) {
  #  bucketed_vars <- c(bucketed_vars, "CNAStat_k", "CNAStat_i")
  #}
  
  # If we'd like for all driver covariates to be selected together in group 
  # lasso (e.g. for driver X, X_MutStat_i, X_CNAStat_i, and X_MethStat_i would 
  # be a group), implement that here
  unique_driver_ids <- unique(unlist(lapply(x_vars, function(x) {
    if(grepl("MutStat_i|CNAStat_i|MethStat_i", x)) {
      id <- unlist(strsplit(x, "_", fixed = T))[1]
      return(id)
    } else {return(NA)}
  })))
  unique_driver_ids <- unique_driver_ids[!is.na(unique_driver_ids)]
  bucketed_vars <- c(bucketed_vars, unique_driver_ids)
  
  # Assign these variables to numerical groups (non-bucketed variables are given
  # their own group number)
  grp_index <- 1
  curr_grp <- ""
  group_assignments <- c()
  
  for(i in 1:length(x_vars)) {
    x <- unlist(strsplit(x_vars[i], "_", fixed = T))[1]
    # This is a bucketed variable
    if(T %fin% unlist(lapply(bucketed_vars, function(var) 
      ifelse(grepl(x, var), T, F)))) {
        # If it's a bucketed variable, but not in the current group, advance the
        # group index
        if ((curr_grp != x) & (curr_grp != "")) {
          grp_index <- grp_index + 1
      } 
      # Change the group assignment and assign the group index
      curr_grp <- x
      group_assignments <- c(group_assignments, grp_index)
    } else {
      # Advance and assign the group index
      grp_index <- grp_index + 1
      group_assignments <- c(group_assignments, grp_index)
      curr_grp <- x
    }
  }
  
  if(type_of_grouping == "group_sizes") {
    unique_groups <- unique(group_assignments)
    new_group_assignments <- sapply(unique_groups, function(g) 
      length(group_assignments[group_assignments == g]))
    return(new_group_assignments)
  }
  #print(x_vars)
  #print(group_assignments)
  
  return(group_assignments)
}


#' If we want to re-run the output from LASSO using a simple, multiple linear 
#' regression model, create a new input DF that has only the covariates selected
#' by LASSO
#' Concept adapted from: https://stats.stackexchange.com/questions/410173/lasso-regression-p-values-and-coefficients
#' @param best_model the model output from glmnet or gglasso
#' @param x_data the columns of the input data frame for the full set of 
#' x-variables
#' @param y_data the column of the input data frame corresponding to the
#' outcome variable
run_regularization_output_in_lm <- function(best_model, x_data, y_data) {
  # Adjust the input DF to pipe the selected covariates back into a LM to get 
  # p-values
  
  x_data_sub <- NA
  if(!is.na(best_model)) {
    
    # Get the non-zero covariates from the the model (selected by group lasso)
    covs <- as.matrix(coef(best_model))
    keep_x <- rownames(covs)[covs != 0]
    keep_x <- keep_x[!keep_x == "(Intercept)"]
    
    # Add these to the standard model coefficients, using unique to ensure no 
    # duplicates
    std_coeff <- unlist(lapply(rownames(covs), function(x) 
      ifelse(grepl("MutStat_i|CNAStat_i|MethStat_i", x), NA, x)))
    std_coeff <- std_coeff[!is.na(std_coeff) & !std_coeff == "(Intercept)"]
    if(debug) {
      print("Std. Coefficients")
      print(head(std_coeff))
    }
    keep_x <- unique(c(keep_x, std_coeff))
    # Subset the input data frame to reflect the standard covariates and the 
    # selected drivers covs
    x_data_sub <- as.data.frame(x_data[, colnames(x_data) %fin% keep_x])
    if(ncol(x_data_sub) == 1) {colnames(x_data_sub) <- as.character(keep_x)}
    
    if(debug) {
      print("X Data, subsetted")
      print(head(x_data_sub))
    }
  }
  
  input_df <- NA
  
  if(!(is.null(colnames(x_data_sub)) | is.na(best_model))) {
    input_df <- as.data.frame(cbind(y_data, x_data_sub))
    names(input_df)[1] <- "ExpStat_k"
    if(debug) {
      print("New input data frame, post-lasso regression")
      print(head(input_df))
    }
  } 
  return(input_df)
}

#' Evaluate the significance of the learned Betas using the
#' selectiveInference package's fixedLassoInf (a subsampling
#' method with a fixed lambda value)
#' @param best_model the LASSO model output from gglasso or glmnet
#' @param x_data the x-data matrix
#' @param y_data the y-data matrix
#' @param optimal_lambda the learned optimal value of lambda from
#' cross-validation
eval_signif_selectiveInference <- function(best_model, x_data, y_data, 
                                           optimal_lambda) {
  # Get p-values for LASSO using the selectiveInference package in R
  # The -1 when getting the Beta is to remove the intercept Beta
  # Need to divide lambda by n because glmnet uses a different objective function 
  # than selectiveInference
  betas <- coef(best_model, x = x_data, y = y_data, 
                s = (optimal_lambda / nrow(x_data)), exact = T)[-1]
  if(length(betas[betas != 0]) > 0) {
    pvals_and_intervals <- fixedLassoInf(x_data, y_data, beta = betas, 
                                         lambda = optimal_lambda)
    pvals <- pvals_and_intervals$pv
    print("P-values")
    print(pvals)
    intervals <- pvals_and_intervals$ci
    sds <- pvals_and_intervals$sd
    
    # Get the p-values for everything, assigning 1 to those variables that 
    # were not selected
    pvals_full <- NA
    intervals_full <- NA
    sds_full <- NA
    index <- 1
    for(i in 1:length(betas)) {
      b <- betas[i]
      if(b == 0) {
        pvals_full[i] <- NA
        intervals_full[i] <- NA
        sds_full[i] <- NA
      }
      else {
        pvals_full[i] <- pvals[index]
        intervals_full[i] <- paste0(intervals[index, 1], intervals[index, 2], 
                                    collapse = "-")
        sds_full <- sds[index]
        index <- index + 1
      }
    }
    print("P-values, post addition of NA's for non-selected variables")
    if(debug) {print(head(pvals_full))}
    
    # Construct an output DF using all the measures we are interested in
    best_model_res <- as.data.frame(best_model$beta)
    colnames(best_model_res) <- c("estimate")
    best_model_res$term <- rownames(best_model_res)
    best_model_res$p.value <- pvals_full
    best_model_res$sds <- sds_full
    best_model_res$conf.int <- intervals_full
    best_model <- best_model_res
  } else {best_model <- NA}
  
  return(best_model)
}

#' A method to generate p-values via a randomization scheme: y_data is 
#' randomized a fixed number of times (e.g. 100) and Betas are learned using the 
#' same group LASSO scheme as for the real data. These Betas become a 
#' distribution for which to perform a one-sample t-test and generate a p-value. 
#' This p-value represents the probability that this Beta is not significantly 
#' different from a null Beta, e.g. that which is generated when there should 
#' absolutely be no correlation between x and y.
#' @param best_model the LASSO model output from gglasso or glmnet
#' @param sample_ids the sample IDs, in the order of the LM input table
#' @param x_data the x-data matrix
#' @param x_vars x-variables from the formula
#' @param y_data the y-data matrix
#' @param v.group a vector of groups within the group LASSO framework 
#' @param lambdas a vector of lambda possibilities for cross-validation; if 
#' using a fixed lambda, this is a vector of length one with a fixed lambda value
#' @param debug a T/ F value indicating whether or not we are in debug mode
#' @param shuffled_input_dfs a list of shuffled input DFs for doing the 
#' randomization_predictor method of p-value generation or for the per-sample 
#' expression method
#' @param ensg the ENSG ID of the target gene of interest, for subsetting the 
#' expression DF appropriately
#' @param num_randomizations the number of times we will randomize our y-data 
#' in order to generate our null distributions for each x
#' @param signif_eval_type the type of randomization we'll be doing 
#' (per_predictor, per_target, or per_sample)
eval_signif_randomization <- function(best_model, sample_ids, x_data, x_vars, 
                                      y_data, v.group, lambdas, debug, 
                                      shuffled_input_dfs, ensg, 
                                      num_randomizations, signif_eval_type) {
  
  # Do the randomization per-target a given number of times and get a 
  # num_randomizations x ncol(x_data) table of random Betas 
  if(signif_eval_type == "randomization_predictors") {
    random_beta_df <- get_predictor_randomized_beta_df(shuffled_input_dfs, 
                                                       sample_ids, x_data, x_vars, 
                                                       y_data, v.group, lambdas)
  } else if (signif_eval_type == "randomization_perTarg") {
    random_beta_df <- get_per_tg_randomized_beta_df(x_data, y_data, v.group, 
                                                    lambdas, num_randomizations)
  } else if (signif_eval_type == "randomization_perSamp") {
    random_beta_df <- get_per_samp_randomized_beta_df(x_data, sample_ids, v.group, 
                                                      lambdas, shuffled_input_dfs, 
                                                      ensg, num_randomizations, 
                                                      debug)
  } else {print(paste("Invalid significance evaluation type:", signif_eval_type))}
  
  if(debug) {print(head(random_beta_df))}
  
  # Get the betas from the model with the real data
  betas_best_model <- as.data.table(t(data.table(as.numeric(best_model$beta))))
  colnames(betas_best_model) <- colnames(x_data)
  if(debug) {print(head(betas_best_model))}
  
  # Run this normal model 9 more times with 90% of the data, in order to get variance
  #betas_more_best_models <- lapply(1:9, function(i) {
  #  
  # Randomly sample 90% of the data (to avoid overfitting)
  #  rows_sampled <- sample(1:nrow(x_data), size = (0.9 * nrow(x_data)))
  #  x_data_sub <- x_data[rows_sampled,]
  #  y_data_sub <- y_data[rows_sampled,]
  
  #  lasso.cv <- cv.gglasso(x_data_sub, y_data_sub, group = v.group, lambda = lambdas, 
  #                         loss = "ls")
  #  optimal_lambda <- lasso.cv$lambda.min
  
  # Run the model with the best lambda
  #  n_best_model <- gglasso(x_data, y_data, lambda = optimal_lambda, 
  #                          group = v.group, loss = "ls")
  #  betas_n <- t(data.table(as.numeric(n_best_model$beta)))
  #  colnames(betas_n) <- colnames(x_data)
  
  #  return(betas_n)
  #})
  #real_betas_df <- rbind(as.data.table(betas_best_model), 
  #                       rbindlist(lapply(betas_more_best_models, as.data.table)))
  #print(head(real_betas_df))
  
  # For each independent variable in the model, compare its Beta to the
  # randomized Betas using a two-sided, one-sample t-test, and return the 
  # t-statistic and p-value
  tstat_pval <- lapply(1:ncol(betas_best_model), function(i) {
    betas_random <- as.numeric(unlist(random_beta_df[,i, with=F]))
    if(debug) {print(paste("Mean of Random Betas:", mean(betas_random)))}

    term <- colnames(betas_best_model)[i]
    best_beta <- as.numeric(unlist(betas_best_model[
      ,colnames(betas_best_model) == term, with = F]))
    
    if(debug) {print(term)}
    
    # Other methods of getting p-values
    #p.prop <- get_prop_pval(as.data.table(betas_best_model), term, i, 
    # random_beta_df, nrow(x_data))
    #p.oneSamp <- t.test(as.numeric(unlist(random_beta_df[,i, with=F])), 
    #mu = best_beta, alternative = "two.sided")
    #p.mapToNorm <- get_fitted_norm_pval(betas_best_model, term, i, random_beta_df)
    
    # Empirical p-value is (r+1)/(n+1), where n is the # of simulations and r is the # of replicates with
    # a Beta at least as large as the one observed
    # Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC379178/#:~:text=However%2C%20Davison%20and%20Hinkley%20(1997,of%20the%20same%20random%20variable.
    betas_random_abs <- abs(betas_random)
    r <- length(betas_random_abs[betas_random_abs >= abs(best_beta)])
    p.empirical <- (r + 1) / (num_randomizations + 1)
    
    return(data.frame("p.value" = p.empirical, "estimate" = best_beta, #"statistic" = tstat, 
                      "term" = term)) #"p.oneSamp" = p.oneSamp$p.value, 
                                      #"p.fittedNorm" = p.mapToNorm, "p.twoSamp" = pval))
  })
  output_df <- rbindlist(tstat_pval)
  
  if(debug) {print(head(output_df))}
  return(output_df)
}


#' Helper function to get per-target gene randomized Beta data frame
#' @param x_data the x-data matrix
#' @param y_data the y-data matrix
#' @param v.group a vector of groups within the group LASSO framework 
#' @param lambdas a vector of lambda possibilities for cross-validation, or one 
#' optimal lambda val
#' @param num_randomizations the number of times we will randomize our y-data 
#' in order to generate our null distributions for each x
get_per_tg_randomized_beta_df <- function(x_data, y_data, v.group, lambdas, 
                                          num_randomizations) {
  
  output_rows <- lapply(1:num_randomizations, function(i) {
    print(paste(i, paste("/", num_randomizations)))
    y_data_random <- dqrng::dqsample(y_data)
    optimal_lambda <- NA
    if(length(lambdas) == 1) {
      optimal_lambda <- lambdas
    } else {
      lasso.cv <- cv.gglasso(x_data, y_data_random, group = v.group, 
                             lambda = lambdas, loss = "ls")
      optimal_lambda <- lasso.cv$lambda.min
    }
    # Run the model with the best lambda
    rand_model <- gglasso(x_data, y_data_random, lambda = optimal_lambda, 
                          group = v.group, loss = "ls")
    betas_rand <- t(data.table(as.numeric(rand_model$beta)))
    colnames(betas_rand) <- colnames(x_data)
    return(betas_rand)
  })
  
  random_beta_df <- rbindlist(lapply(output_rows, as.data.table), 
                              use.names = T)
  
  return(random_beta_df)
}

#' Helper function to get per-sample randomized Beta data frame
#' @param x_data the x-data matrix
#' @param sample_ids a vector of sample IDs in the x-data matrix
#' @param v.group a vector of groups within the group LASSO framework 
#' @param lambdas a vector of lambda possibilities for cross-validation, or one 
#' optimal lambda val
#' @param shuffled_input_dfs a list of shuffled input DFs for doing the 
#' randomization_predictor method of p-value generation or for the per-sample 
#' expression method
#' @param ensg the ensg id of the target gene, for subsetting the expression DF
#' @param num_randomizations the number of times we will randomize our y-data 
#' in order to generate our null distributions for each x
#' @param debug a T/F value indicating if we are in debug mode
get_per_samp_randomized_beta_df <- function(x_data, sample_ids, v.group, lambdas, 
                                            shuffled_input_dfs, ensg, 
                                            num_randomizations, debug) {
  
  expression_df <- shuffled_input_dfs[[1]]
  index <- which(grepl(ensg, expression_df$ensg_id))
  
  output_rows <- lapply(shuffled_input_dfs, function(df) {
    y_data_rand <- t(as.matrix(as.numeric(
      df[index, which(colnames(df) %fin% sample_ids), with = F])))
    if(debug) {print(head(y_data_rand))}
    
    optimal_lambda <- NA
    if(length(lambdas) == 1) {optimal_lambda <- lambdas} 
    else {
      lasso.cv <- cv.gglasso(x_data, y_data_rand, group = v.group, 
                             lambda = lambdas, loss = "ls")
      optimal_lambda <- lasso.cv$lambda.min
    }
    
    # Run the model with the best lambda
    rand_model <- gglasso(x_data, y_data_rand, lambda = optimal_lambda, 
                          group = v.group, loss = "ls")
    betas_rand <- t(data.table(as.numeric(rand_model$beta)))
    colnames(betas_rand) <- colnames(x_data)
    return(betas_rand)
  })
  
  random_beta_df <- rbindlist(lapply(output_rows, as.data.table), 
                              use.names = T)
  
  return(random_beta_df)
}


#' Helper function to get per-driver/ predictor randomized Beta data frame
#' @param shuffled_input_dfs the shuffled x-data matrices
#' @param sample_ids the sample IDs, to subset the shuffled dfs
#' @param x_data the x-data matrix (has the real target-level features)
#' @param x_vars the x-variables from the formula
#' @param y_data the y-data matrix
#' @param v.group a vector of groups within the group LASSO framework 
#' @param lambdas a vector of lambda possibilities for cross-validation, or one 
#' optimal lambda val
get_predictor_randomized_beta_df <- function(shuffled_input_dfs, sample_ids, 
                                             x_data, x_vars, y_data, v.group, 
                                             lambdas) {
  # Keep only the non-driver columns of the original input DF
  cols_to_keep <- setdiff(colnames(x_data), colnames(shuffled_input_dfs[[1]]))
  x_data_sub <- x_data[,which(colnames(x_data) %fin% cols_to_keep)]
  
  output_rows <- lapply(shuffled_input_dfs, function(input_df) {
    #print(paste(i, paste("/", num_randomizations)))
    input_df <- input_df[input_df$sample_id %fin% sample_ids, 
                         colnames(input_df) != "sample_id", with = F]
    input_df <- input_df[,colnames(input_df) %fin% x_vars, with = F]
    x_data_full <- cbind(as.matrix(input_df), x_data_sub)
    #print(head(x_data_full))
    rows_na <- which(rowSums(is.na(x_data_full)) > 0)
    rows_keep <- setdiff(1:nrow(x_data_full), rows_na)
    x_data_full <- as.matrix(x_data_full[rows_keep,])
    y_data <- as.matrix(y_data[rows_keep,])
    
    optimal_lambda <- NA
    if(length(lambdas) == 1) {
      optimal_lambda <- lambdas
    } else {
      lasso.cv <- cv.gglasso(x_data_full, y_data, group = v.group, 
                             lambda = lambdas, loss = "ls")
      optimal_lambda <- lasso.cv$lambda.min
    }
    
    # Run the model with the best lambda
    rand_model <- gglasso(x_data_full, y_data, lambda = optimal_lambda, 
                          group = v.group, loss = "ls")
    betas_rand <- t(data.table(as.numeric(rand_model$beta)))
    colnames(betas_rand) <- colnames(x_data)
    return(betas_rand)
  })
  
  #fn <- paste0("opt_lambdas_", paste0(runif(1,1,50), ".png"))
  #png(fn, width = 450, height = 450)
  #hist(optimal_lambdas, xlab = "optimal lambda for randomizations", ylab = "frequency")
  #dev.off()
  
  random_beta_df <- rbindlist(lapply(output_rows, as.data.table), use.names = T)
  
  return(random_beta_df)
}

#' Get proportions-based p-value, simply using the randomized distribution and 
#' pnorm built-in function. Returns the p-value.
#' Useful videos: https://www.youtube.com/watch?v=Xz_yMO0iHF8; https://www.khanacademy.org/math/ap-statistics/xfb5d8e68:inference-categorical-proportions/xfb5d8e68:carrying-out-test-proportion/v/calculating-a-z-statistic-in-a-significance-test
#' @param betas_best_model the Beta values for the best model
#' @param term the current term we want to calculate a p-value for
#' @param random_beta_df the randomized Beta values
#' @param i the index for our random beta DF
#' @param n the number of samples
get_prop_pval <- function(betas_best_model, term, i, random_beta_df, n) {
  # We are treating these like proportions, so we want to take the absolute value, 
  # then do a one-tailed test (is p_hat significantly greater than p?)
  p_hat <- abs(as.numeric(unlist(betas_best_model[
    ,colnames(betas_best_model) == term, with = F])))
  if(!(p_hat == 0)) {
    p <- abs(mean(as.numeric(unlist(random_beta_df[,i, with=F]))))
    sd <- p_hat * (1-p_hat)
    se <- sqrt(sd / n)
    z <- (p_hat - p) / se
    #p.prop <- pnorm(z, mean = p, sd = sd, lower.tail = F)    # multiply by 2 for two-tailed test
    p.prop <- pnorm(z, mean = p, sd = sd, lower.tail = F) 
  } else {p.prop <- 1}
  
  return(p.prop)
}

#' Get fitted normal distribution p-value, using the fitdistrplus package in R 
#' Reference: https://www.biostars.org/p/108166/
#' @param betas_best_model the Beta values for the best model
#' @param term the current term we want to calculate a p-value for
#' @param random_beta_df the randomized Beta values
#' @param i the index for our random beta DF
get_fitted_norm_pval <- function(betas_best_model, term, i, random_beta_df) {
  p_hat <- as.numeric(unlist(betas_best_model[,colnames(betas_best_model) == term, 
                                              with = F]))
  rand_betas <- as.numeric(unlist(random_beta_df[,i, with=F]))
  
  tryCatch({
    fit <- fitdist(rand_betas, "norm")  # by maximum likelihood
    
    # Map p_hat to the fitted distribution
    p_value <- pnorm(p_hat, mean = fit$estimate[['mean']], 
                     sd = fit$estimate[['sd']])
    print(paste("P-value:", p_value))
    
    return(p_value)
    
  }, error = function(cond) {
    print(cond)
    return(NA)
  })
  
}

#' This function uses a method called stability selection to get a "stability 
#' score" using a subsampling approach. We repeat our group LASSO B times on a 
#' random subset of half the data, and in every run, check which features are 
#' chosen in the top L. After the B subsamples, we assign a score to each feature 
#' based on how often it was selected as part of L. B tends to be between 100-500, 
#' L typically about 5, and a stability score >= 0.6 tends to be a good 
#' "significance" threshold. NOTE: this function extends this by using the Betas 
#' from the subsamples to perform a one-sample t-test and generate a p-value. 
#' Both the stability score and the p-value are reported and returned.
#' @param best_model the LASSO model output from gglasso or glmnet
#' @param x_data the x-data matrix
#' @param y_data the y-data matrix
#' @param v.group a vector of groups within the group LASSO framework 
#' @param lambdas a vector of lambda possibilities for cross-validation
#' @param debug a T/ F value indicating whether or not we are in debug mode
#' @param B the number of subsampling events
#' @param L the size of the top feature set 
eval_signif_subsampling <- function(best_model, x_data, y_data, v.group, lambdas, 
                                    debug, B = 100, L = 5) {
  # Do the subsampling a given number of times (B) and get a B x num_features
  # table of sampled Betas
  N <- nrow(x_data) / 2
  
  output_rows <- mclapply(1:B, function(i) {
    # Get random row numbers for a subsample
    rows_in_subsample <- sample(1:nrow(x_data), size = N)
    x_data_sub <- x_data[rows_in_subsample,]
    y_data_sub <- y_data[rows_in_subsample]
    
    lasso.cv <- cv.gglasso(x_data_sub, y_data_sub, group = v.group, 
                           lambda = lambdas, loss = "ls")
    optimal_lambda <- lasso.cv$lambda.min
    # Run the model with the best lambda
    sub_model <- gglasso(x_data_sub, y_data_sub, lambda = optimal_lambda, 
                         group = v.group, loss = "ls")
    betas_subsamp <- t(data.table(as.numeric(sub_model$beta)))
    colnames(betas_subsamp) <- colnames(x_data)
    return(betas_subsamp)
  })
  subsampled_beta_df <- rbindlist(lapply(output_rows, as.data.table), 
                                  use.names = T)
  colnames(subsampled_beta_df) <- colnames(x_data)
  
  #if(debug) {print(head(subsampled_beta_df))}
  
  # Get the top L terms per subsample, excluding nuisance variables
  topL_persubsample <- lapply(1:nrow(subsampled_beta_df), function(i) {
    # Get the top L terms for each subsample
    sorted_df <- sort(abs(subsampled_beta_df[i,]), decreasing=T)
    # Eliminate nuisance variables, if desired
    sorted_df <- sorted_df[, which(grepl("Stat_i", colnames(sorted_df))), 
                           with = F]
    # Remove those with Beta = 0
    sorted_df <- sorted_df[, colSums(sorted_df != 0) > 0, with = F]
    # If we still have columns (at least one non-nuisance variable was non-zero), 
    # return it; o.w. NA
    if(ncol(sorted_df) > 0) {
      # Pad with NA if number of remaining columns is less than L
      if(ncol(sorted_df) < L) {
        na_df <- data.table(matrix(nrow = 1, ncol = (L - ncol(sorted_df))))
        sorted_df <- cbind(sorted_df, na_df)
      }
      return(t(as.data.frame(colnames(sorted_df)[1:L])))
    } else {return(NA)}
  })
  topL_persubsample <- topL_persubsample[!is.na(topL_persubsample)]
  topL_persubsample <- rbindlist(lapply(topL_persubsample, as.data.frame))
  
  #best_model_tab <- data.table(as.numeric(unlist(best_model$beta)))
  
  
  # For each independent variable in the model, compute and return the stability 
  # score and the fraction of times that it was non-zero
  ss_and_fracNon0 <- lapply(1:ncol(subsampled_beta_df), function(i) {
    
    term <- colnames(subsampled_beta_df)[i]
    
    # Compute the stability score (fraction of times this term is present in the 
    # top L)
    rows_w_term <- unlist(lapply(1:nrow(topL_persubsample), function(i) {
      row <- topL_persubsample[i,]
      if(term %fin% row) {return(i)}
      else {return(NA)}
    }))
    rows_w_term <- rows_w_term[!is.na(rows_w_term)]
    ss <- nrow(topL_persubsample[rows_w_term,]) / nrow(subsampled_beta_df)
    
    # Compute the number/ fraction of times it was nonzero
    betas_col <- as.numeric(unlist(subsampled_beta_df[,i, with = F]))
    frac_non0 <- length(betas_col[betas_col > 0])
    #frac_non0 <- length(betas_col[betas_col > 0]) / length(betas_col)
    
    # Get the actual estimate
    est <- as.numeric(unlist(best_model$beta[rownames(best_model$beta) == term,]))
    
    return(data.frame("estimate" = est, "term" = term, "stability.score" = ss, 
                      "frac.nonzero" = frac_non0))
  })
  output_df <- rbindlist(ss_and_fracNon0)
  
  if(debug) {print(head(output_df))}
  return(output_df)
}
