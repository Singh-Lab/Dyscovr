############################################################
### PERFORM COLLINEARITY COLLECTION
### Written by Sara Geraghty, Princeton University
### https://www.biorxiv.org/content/10.1101/2024.11.20.624509v1
############################################################

# Helper functions to remove covariates that are correlated/ collinear

############################################################

#' Calculates a Spearman correlation matrix among all variables in an input 
#' table and, according to the provided method, strategically excludes them if 
#' they are significantly correlated 
#' @param lm_input_table a single LM input table to linear model
#' @param formula a LM formula, for using a VIF metric of correction
#' @param method the way to handle collinear variables; either "eliminate" or 
#' "interaction_term", with an optional '_vif' or on the end to indicate we are 
#' using VIF rather than Spearman/ Jaccard
#' @param outfile the combined path and filename of the file where we are 
#' dynamically storing the variables that we remove for each target gene
#' @param vif_thres a vif threshold, default is 5
correct_collinearity <- function(lm_input_table, formula, method, target, #outfile, 
                                 vif_thres = 5) {
  
  lm_input_matrix <- as.matrix(lm_input_table[,!(colnames(lm_input_table) %fin% 
                                                   c("ExpStat_k", "sample_id")), 
                                              with = F])
  variables_removed <- c()
  
  # Remove variables that are entirely 0 
  lm_input_matrix <- lm_input_matrix[, colSums(lm_input_matrix != 0) > 0]
  
  # 1. SPEARMAN CORRELATION METHOD
  
  # Remove cancer subtypes that are correlated to a driver mutation status using 
  # Spearman. Use the Hmisc package to calculate a correlation matrix, get 
  # the p-values
  corr_res <- correct_collinearity_spearman(lm_input_matrix, method)
  lm_input_matrix <- corr_res[[1]]
  variables_removed <- c(variables_removed, as.character(unlist(corr_res[[2]])))
  
  # 2. VIF/ TOLERANCE COLLINEARITY METHOD
  
  # Make new formula
  form_spl <- unlist(strsplit(formula, " ~ ", fixed = T))
  form_spl_x <- unlist(strsplit(form_spl[2], " + ", fixed = T))
  form_spl_x <- form_spl_x[form_spl_x %fin% colnames(lm_input_matrix)]
  formula_new <- paste(as.character(form_spl[1]), 
                       paste(form_spl_x, collapse = " + "), sep = " ~ ")
  
  lm_fit <- NA
  lm_data <- cbind(lm_input_table[,c("sample_id", "ExpStat_k"), with = F], 
                   lm_input_matrix)
  tryCatch({
    lm_fit <- lm(data = lm_data, formula = formula_new)
  }, error = function(cond) {
    print(formula_new)
    print(cond)
    print(lm_data)
  })
  
  # Calculate VIF
  vif_res <- calculate_vif(lm_fit, lm_data, formula_new)
  if(!(is.na(vif_res[[2]]))) {
    variables_removed <- c(variables_removed, as.character(vif_res[[2]]))
  }
  vif_res_new <- vif_res[[1]]
  
  # Get the names of the variables with a VIF score that exceeds some threshold
  #ind_sig_vals <- which(vif_res_new > vif_thres)
  #names_all <- unlist(strsplit(unlist(strsplit(
  #  formula_new, " ~ ", fixed = T))[2], " + ", fixed = T))
  #names <- names_all[ind_sig_vals]
  sig_vals <- vif_res_new[vif_res_new > vif_thres]
  names <- names(sig_vals)

  
  # Keep only the variables that are *not* bucketed (these will automatically 
  # have high VIF scores because of their linearity with the other buckets)
  names <- names[!grepl("_b", names)]
  
  # Exclude other variables that exceed this VIF threshold
  if(length(names) > 0) {
    
    # Check if there are driver mutation variables here - we want to keep these 
    # and make sure their VIF goes down
    if(T %fin% grepl("MutStat_i", names)) {
      mutstat_vars <- names[grepl("MutStat_i", names)]
      
      while(length(mutstat_vars) > 0) {
        names <- names[!(names %fin% mutstat_vars)]
        if(length(names) == 0) {break}
        lm_input_matrix <- lm_input_matrix[,!(colnames(lm_input_matrix) %fin% 
                                                names)]
        variables_removed <- c(variables_removed, names)
        
        # Make new formula
        form_spl <- unlist(strsplit(formula_new, " ~ ", fixed = T))
        form_spl_x <- unlist(strsplit(form_spl[2], " + ", fixed = T))
        form_spl_x <- form_spl_x[form_spl_x %fin% colnames(lm_input_matrix)]
        formula_new <- paste(as.character(form_spl[1]), 
                             paste(form_spl_x, collapse = " + "), sep = " ~ ")
        
        # Fit the new model
        lm_fit <- NA
        lm_data <- cbind(lm_input_table[,c("sample_id", "ExpStat_k"), with = F], 
                         lm_input_matrix)
        tryCatch({
          lm_fit <- lm(data = lm_data, formula = formula_new)
        }, error = function(cond) {
          print(formula_new)
          print(cond)
          print(lm_data)
        })
        
        # Recalculate the VIF scores, and make sure that the new MutStat_i 
        # variable's VIF is lower than the threshold (if not, repeat until it is)
        vif_res <- calculate_vif(lm_fit, lm_data, formula_new)
        if(!(is.na(vif_res[[2]]))) {
          variables_removed <- c(variables_removed, as.character(vif_res[[2]]))
        }
        vif_res_new <- vif_res[[1]]

        #ind_sig_vals <- which(vif_res_new > vif_thres)
        #names_all <- unlist(strsplit(unlist(strsplit(
          #formula_new, " ~ ", fixed = T))[2], " + ", fixed = T))
        #names <- names_all[ind_sig_vals]
	names <- names(vif_res_new[vif_res_new > vif_thres])
        names <- names[!grepl("_b", names)]
        mutstat_vars <- names[grepl("MutStat_i", names)]
      }
      
    } else {
      lm_input_matrix <- lm_input_matrix[,!(colnames(lm_input_matrix) %fin% names)]
      variables_removed <- c(variables_removed, names)
    }
  }
  
  print("Variables to be removed due to collinearity:")
  print(variables_removed)
  print(list(paste(as.character(unlist(variables_removed)), collapse = ",")))
  
  # Write removed variables to file
  variables_rm_df <- data.table('Target' = target,
                                'Variables.Rm' = paste(as.character(unlist(variables_removed)),
                                                          collapse = ","))
  
  lm_input_matrix <- cbind(lm_input_table[,c("sample_id", "ExpStat_k"), 
                                          with = F], lm_input_matrix)
  return(list(as.data.frame(lm_input_matrix), variables_rm_df))
}


############################################################

#' Helper function to correct collinearity using a Spearman correlation/ Jaccard 
#' similarity method Spearman is for computing correlation between continuous 
#' variables; Jaccard is for measuring similarity between binary variables or 
#' continuous/binary combinations
correct_collinearity_spearman <- function(lm_input_matrix, method) {
  
  variables_removed <- c()
  
  # Use the Hmisc package to calculate a correlation matrix, get the p-values
  corr_mat <- as.data.frame(rcorr(lm_input_matrix, type = "spearman")$r)
  corr_mat_p <- as.data.frame(rcorr(lm_input_matrix, type = "spearman")$P)
  
  combos <- get_spearman_combos(corr_mat, corr_mat_p)
  
  # Select which variables to keep, if we are eliminating collinear variables, or
  # create interaction terms for them. Always eliminate the subtype variable 
  # first. Keep the target's covariates over the others.
  for(c in combos) {
    c_spl <- unlist(strsplit(c, ":", fixed = T)) 
    var1 <- as.character(c_spl[1])
    var2 <- as.character(c_spl[2])
    
    # Ignore the bucketed variables or variables for the same driver (we expect 
    # these to be collinear, by definition)
    if((!is.null(var1)) & (!is.null(var2))) {
      var1_spl <- unlist(strsplit(var1, "_", fixed = T))
      var2_spl <- unlist(strsplit(var2, "_", fixed = T))
      if(!(var1_spl[1] == var2_spl[1])) {
        if(!(grepl("PC", var1) & grepl("PC", var2))) {
          
          # Regardless, remove cancer subtypes that are based on mutation status 
          # of a driver
          if((var1_spl[length(var1_spl)] == "i") & (var2_spl[1] == "Cancer")) {
            lm_input_matrix <- lm_input_matrix[,colnames(lm_input_matrix) != var2]
            variables_removed <- c(variables_removed, var2)
          }
          if((var2_spl[length(var2_spl)] == "i") & (var1_spl[1] == "Cancer")) {
            lm_input_matrix <- lm_input_matrix[,colnames(lm_input_matrix) != var1]
            variables_removed <- c(variables_removed, var1)
          }
          
          if(!grepl("vif", method)) {
            
            # Elimination Approach
            if(method == "eliminate") {
              # We can't predict a factor well if it is correlated to the 
              # methylation, CNA, or mutation status of that target gene itself, 
              # even though this is very much in the realm of possibility
              if(var1_spl[length(var1_spl)] == "k") {
                if(var2_spl[length(var2_spl)] != "k") {
                  if(var2 %fin% colnames(lm_input_matrix)) {
                    lm_input_matrix <- lm_input_matrix[,colnames(lm_input_matrix) 
                                                       != var2]
                    variables_removed <- c(variables_removed, var2)
                  }
                }
              }
              if(var2_spl[length(var2_spl)] == "k") {
                if(var1_spl[length(var1_spl)] != "k") {
                  if(var1 %fin% colnames(lm_input_matrix)) {
                    lm_input_matrix <- lm_input_matrix[,colnames(lm_input_matrix) 
                                                       != var1]
                    variables_removed <- c(variables_removed, var1)
                  }
                }
              }
              # NOTE: If both of these terms are for nuisance variables like age 
              # and gender, keep them
            
            # Interaction Term Approach
            } else  {
              
              if(method != "interaction_term") {
                print("Only 'eliminate' and 'interaction_term' are implemented 
                methods for dealing with variable collinearity. Using interaction 
                      term method.")
              }
              if(((var1_spl[length(var1_spl)] == "k") & 
                  (var2_spl[length(var2_spl)] != "k")) | 
                 ((var1_spl[length(var1_spl)] != "k") & 
                  (var2_spl[length(var2_spl)] == "k"))) {
                new_var <- paste(var1, paste(var2, "interac", sep = "_"), 
                                 sep = "_")

                vals <- unlist(lapply(1:nrow(lm_input_matrix), function(i) {
                  e1 <- unique(as.numeric(
                    lm_input_matrix[i,colnames(lm_input_matrix) == var1])) 
                  e2 <- unique(as.numeric(
                    lm_input_matrix[i,colnames(lm_input_matrix) == var2])) 
                  return(ifelse(((e1 > 0) & (e2 > 0)) | 
                                  ((e1 < 0) & (e2 < 0)), 1, 0))
                }))
                
                lm_input_matrix <- cbind(lm_input_matrix, vals)
                colnames(lm_input_matrix)[ncol(lm_input_matrix)] <- new_var
              }
            } 
          }
        }
      } 
    } 
  }
  
  return(list(lm_input_matrix, unique(variables_removed)))
}

############################################################

#' Helper function to get combinations of significantly correlated variables 
#' from a Spearman correlation matrix
#' @param corr_mat spearman correlation matrix
#' @param corr_mat_p spearman correlation pvalue matrix
#' @param sig_thres the p-value significance threshold for considering
#' a pair of variables to be correlated (default is 1E-05)
#' @param corr_thres the absolute value of the minimum spearman correlation to 
#' qualify two variables as being correlated (default is 0.7)
get_spearman_combos <- function(corr_mat, corr_mat_p, sig_thres = 1e-05, 
                                corr_thres = 0.7) {
  
  names <- rownames(corr_mat)
  num_rows <- nrow(corr_mat)
  
  combos <- unique(unlist(lapply(1:ncol(corr_mat), function(i) {
    curr_var <- colnames(corr_mat)[i]
    curr_vals <- corr_mat[i:num_rows, i]
    curr_pvals <- corr_mat_p[i:num_rows, i]
    curr_names <- names[i:num_rows]
    
    # Remove NA/NAN values (assume these are not significant)
    na_nan_ind <- which(is.na(curr_pvals) | is.nan(curr_pvals))
    curr_pvals <- curr_pvals[-na_nan_ind]
    curr_names <- curr_names[-na_nan_ind]
    curr_vals <- curr_vals[-na_nan_ind]

    if(any(curr_pvals < sig_thres)) {
      corr_ind <- intersect(which(curr_pvals < sig_thres), 
                            which(abs(curr_vals) > corr_thres))
      corr_vars <- unlist(curr_names[corr_ind])
      if(length(corr_vars) > 0) {
        tuples <- lapply(corr_vars, function(v) 
          return(paste(curr_var, v, sep = ":")))
        return(tuples)
      } else {return(NA)}
    }
  })))
  combos <- combos[!is.na(combos)]

  return(combos)
}

############################################################

#' Calculate VIF scores using either the caret package
#' @param lm_fit the fitted linear regression
#' @param lm_input_table the LM input table, in case we have aliased coefficients
#' @param formula the LM formula, in case we have aliased coefficients
calculate_vif <- function(lm_fit, lm_input_table, formula) {
  
  variables_removed <- NA
  
  tryCatch({
    vif_res <- car::vif(lm_fit, type = "predictor")
    
  }, error=function(cond) {
    # Handle the case of aliased coefficients
    if(grepl("aliased coefficients", cond)) {
      boolF <- F
      variables_removed <- c()
      i <- 1
      while(boolF == F) {
        lm_input_table_sub <- lm_input_table[, !(colnames(lm_input_table) %fin% 
                                                   c("ExpStat_k", "sample_id")), 
                                             with = F]
        
        # Create a correlation DF for input variables and ID those that are 
        # perfectly correlated (aliased)
        cor_df <- cor(lm_input_table_sub)
        cor_df[row(cor_df) == col(cor_df)] <- NA
        ind <- rbind(which(cor_df == 1, arr.ind = T), 
                     which(cor_df == -1, arr.ind = T))
        
        if(is.null(ind)) {return(NA)}
        
        # Remove aliased variables
        variables_removed <- c(variables_removed, unique(rownames(ind)))
        
        # Recreate formula and input table without these variables
        formula_spl <- unlist(strsplit(formula, " ~ ", fixed = T))
        formula_spl_x <- unlist(strsplit(formula_spl[2], " + ", fixed = T))
        formula_spl_x <- formula_spl_x[!(formula_spl_x %fin% variables_removed)]
        formula_new <- paste(as.character(formula_spl[1]), 
                             paste(formula_spl_x, collapse = " + "), sep = " ~ ")
        
        lm_input_table <- lm_input_table[,!(colnames(lm_input_table) %fin% 
                                              variables_removed), with = F]
        lm_fit <- lm(data = lm_input_table, formula = formula_new)
        
        try({
          # Re-run VIF 
          vif_res <- car::vif(lm_fit, type = "predictor")
          boolF <- T
        })
        i <- i + 1
        
        # Quit if we have reached an endless loop
        if(i > 5) {
          return(NA)
        }
      }
      
    } else {
      print(cond)
      return(NA)
    }
  })
  
  return(list(vif_res, variables_removed))
}


############################################################

#' A function to combine calculated collinearity diagnostic statistics, 
#' specifically the tolerance, VIF, eigenvalue, and condition index values for 
#' each covariate, and returns these aggregate statistics.
#' @param path a path to the folder with all of the collinearity output files
#' @param outfn the file name for the output files, which has all the details of 
#' the current run
#' @param debug if we are in debug mode and should add additional prints
#' @param randomize whether we have randomized this run or not
combine_collinearity_diagnostics <- function(path, outfn, debug, randomize) {
  vif_tabs <- list.files(path, pattern = "vif_tab")
  eig_cindex_tabs <- list.files(path, pattern = "eig_cindex_tab")
  
  # Get only the tables that match the rest of the run as well
  outfn_spl <- unlist(strsplit(outfn, "_", fixed = T))
  outfn_sub <- paste(outfn_spl[3:length(outfn_spl)], collapse = "_")
  
  vif_tabs <- vif_tabs[grepl(outfn_sub, vif_tabs)]
  eig_cindex_tabs <- eig_cindex_tabs[grepl(outfn_sub, eig_cindex_tabs)]
  
  if(randomize) {
    vif_tabs <- vif_tabs[grepl("RANDOMIZED", vif_tabs)]
    eig_cindex_tabs <- eig_cindex_tabs[grepl("RANDOMIZED", eig_cindex_tabs)]
  } else {
    vif_tabs <- vif_tabs[!grepl("RANDOMIZED", vif_tabs)]
    eig_cindex_tabs <- eig_cindex_tabs[!grepl("RANDOMIZED", eig_cindex_tabs)]
  }
  
  
  if(debug) {
    print("VIF Tabs")
    print(head(vif_tabs))
    print("Eig Cindex Tabs")
    print(head(eig_cindex_tabs))
  }
  
  tmp <- fread(paste(path, vif_tabs[1], sep = "/"), header = T)
  results_table <- data.frame(matrix(ncol= 5, nrow = nrow(tmp)+1))
  colnames(results_table) <- c("Covariate", "Tolerance", "VIF", "Eigenvalue", 
                               "Condition_Index")
  results_table$Covariate <- c(tmp$Variables, "--")
  
  if(debug) {print(head(results_table))}
  
  for (i in 1:length(vif_tabs)) {
    combo <- paste(unlist(strsplit(vif_tabs[i], "_", fixed = T))[1:2], 
                   collapse = "_")
    if(debug) {print(paste("TargGene and Regprot Combo:", combo))}
    vif_tab <- fread(paste(path, vif_tabs[i], sep = "/"), header = T)
    tryCatch({
      eig_cindex_tab <- fread(paste(path, eig_cindex_tabs[
        grepl(combo, eig_cindex_tabs)], sep = "/"), header = T)
    }, error=function(cond){
      print(cond)
      tab_name <- eig_cindex_tabs[grepl(combo, eig_cindex_tabs)]
      print("Eig cindex tab name:", tab_name)
      eig_cindex_tab <- fread(paste(path, tab_name[[1]], sep = "/"), 
                              header = T)
    })
    
    for(i in 1:nrow(results_table)) {
      var <- results_table$Covariate[i]
      
      if (!var == "--") {
        tolerance <- as.numeric(vif_tab[vif_tab$Variable == var, 'Tolerance'])
        if(debug) {print(paste("Tolerance", tolerance))}
        if(length(tolerance) > 0) {
          results_table[i, 'Tolerance'] <- mean(c(as.numeric(
            results_table[i, 'Tolerance']), tolerance), na.rm = T)
        } 
        vif <- as.numeric(vif_tab[vif_tab$Variable == var, 'VIF'])
        if(length(vif) > 0) {
          results_table[i, 'VIF'] <- mean(c(as.numeric(
            results_table[i, 'VIF']), vif), na.rm = T)
        }
        
      } 
      
      results_table[i, 'Eigenvalue'] <- mean(c(as.numeric(
        results_table[i, 'Eigenvalue']), as.numeric(
          eig_cindex_tab$Eigenvalue[i])), na.rm = T)
      results_table[i, 'Condition_Index'] <- mean(c(as.numeric(
        results_table[i, 'Condition_Index']), as.numeric(
          eig_cindex_tab$`Condition Index`[i])), na.rm = T)
    }
  }
  return(results_table)
}
