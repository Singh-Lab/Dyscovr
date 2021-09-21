############################################################
### Linear Model Helper Functions
### Written By: Sara Geraghty, July 2020
############################################################

# This file contains a set of helper functions that will be called by the linear
# model (from linear_model.R or linear_model.R) while it is running.


############################################################

#' Filters a mutation target DF in order to keep only rows that overlap the 
#' gene targets of interest
#' @param mutation_targ_df a mutation target DF to be subsetted
#' @param uniprot_ids a list of uniprot IDs to keep in the DF
filter_mut_by_uniprot <- function(mutation_targ_df, uniprot_ids) {
  mutation_rows_to_keep <- unique(unlist(lapply(uniprot_ids, function(x) {
    #print(x)
    return(grep(x, mutation_targ_df$Swissprot))
  })))
  #print(mutation_rows_to_keep)
  mutation_targ_df <- mutation_targ_df[mutation_rows_to_keep,]
  #print("Mutation Targ DF filt by Uniprot ID:")
  #print(mutation_targ_df)
  return(mutation_targ_df)
}
############################################################

#' Filters a CNA DF in order to keep only rows that overlap the 
#' gene targets of interest
#' @param cna_df a copy number alteration data frame to be subsetted
#' @param ensg_ids a list of ensembl IDs to keep in the DF
filter_cna_by_ensg <- function(cna_df, ensg_ids) {
  cna_rows_to_keep <- unique(unlist(lapply(ensg_ids, function(x) 
    grep(x, cna_df$ensg_id))))
  cna_df <- cna_df[cna_rows_to_keep,]
  #print("CNA filt by ENSG:")
  #print(cna_df)
  return(cna_df)
}
############################################################

#' Filters a methylation DF in order to keep only rows that overlap the 
#' gene targets of interest
#' @param methylation_df a methylation data frame to be subsetted
#' @param ensg_ids a list of ensembl IDs to keep in the DF
filter_meth_by_ensg <- function(methylation_df, ensg_ids) {
  meth_rows_to_keep <- unique(unlist(lapply(ensg_ids, function(x) 
    grep(x, methylation_df$ensg_ids))))
  methylation_df <- methylation_df[meth_rows_to_keep,]
  #print("Meth filt by ENSG:") 
  #print(methylation_df)
  return(methylation_df)
}

############################################################

#' Given a parameter for if/how we are bucketing CNA, extracts and 
#' processes the CNA value from the CNA data frame
#' @param cna_df a CNA data frame of interest
#' @param sample the sample ID of the current sample
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. 
#' Possible values are "bucketInclAmp", "bucketExclAmp", and "rawCNA"
get_cna_stat <- function(cna_df, sample, cna_bucketing) {
  
  if(cna_bucketing == "rawCNA") {
    cna_stat <- unique(as.integer(unlist(cna_df[,colnames(cna_df) == sample, with = FALSE])))
    if (length(cna_stat) > 1) {cna_stat <- cna_stat[1]}
    cna_stat <- log2(cna_stat + 1)
    if(length(cna_stat) < 1) {cna_stat <- NA}
    if (is.nan(cna_stat)) {cna_stat <- NA}
  } else {
    str_buckets <- as.character(unlist(cna_df[,colnames(cna_df) == sample, with = FALSE]))
    if(!(length(str_buckets) == 0)) {
      buckets_list <- unlist(strsplit(unlist(strsplit(str_buckets, "(", fixed = TRUE))[2], ")", fixed = TRUE))[1]
      cna_stat <- list(as.integer(unlist(strsplit(buckets_list, ",", fixed = TRUE))))
    } else {cna_stat <- list(NA, NA, NA)}
    print(cna_stat)
  }  
  return(cna_stat)
}

############################################################

#' Given a mutation target DF and a sample of interest, extracts
#' the 0 or 1 mutation statistic to indicate non-mutated/mutated
#' @param mutation_targ_df a CNA data frame of interest
#' @param sample the sample ID of the current sample
get_mut_stat_targ <- function(mutation_targ_df, sample) {
  mut_count <- as.numeric(unlist(mutation_targ_df[, grepl(sample, colnames(mutation_targ_df)), 
                                                  with = FALSE]))
  # If there is no entry here, then it was not mutated in any patients
  if (length(mut_count) == 0) {mut_stat <- 0} 
  else if (length(mut_count) == 1) {
    if (mut_count > 0) {mut_stat <- 1} 
    else {mut_stat <- 0}
  } else {
    mut_count <- mean.default(mut_count)
    if(mut_count > 0) {mut_stat <- 1} 
    else {mut_stat <- 0}
  }
  return(mut_stat)
}


############################################################

#' Given a parameter for if we are bucketing methylation, extracts and 
#' processes the methylation value from the methylaton data frame
#' @param methylation_df a methylation data frame of interest
#' @param sample the sample ID of the current sample
#' @param meth_bucketing a T/F value indicating if we are bucketing 
#' methylation values. 
#' @param tumNormMatched a TRUE/FALSE value indicating whether or not the analysis is 
#' tumor-normal matched
get_meth_stat <- function(methylation_df, sample, meth_bucketing, tumNormMatched) {
  
  if(!tumNormMatched) {
    if(meth_bucketing) { 
      str_buckets <- as.character(unlist(methylation_df[,grepl(sample, colnames(methylation_df)), with = FALSE]))
      #print(str_buckets)
      buckets_list <- unlist(strsplit(unlist(strsplit(str_buckets, "(", fixed = TRUE))[2], ")", fixed = TRUE))[1]
      meth_stat <- list(as.integer(unlist(strsplit(buckets_list, ",", fixed = TRUE))))
      #print(meth_stat)
    } else {  
      meth_stat <- mean.default(as.numeric(methylation_df[,grepl(sample, colnames(methylation_df)), 
                                                          with = FALSE]))
      if(length(meth_stat) < 1) {meth_stat <- NA}
      if (length(meth_stat) > 1) {meth_stat <- meth_stat[1]}
      if (is.nan(meth_stat)) {meth_stat <- NA}
    }
  } else {
    meth_stat <- get_tnm_methylation(methylation_df, sample, meth_bucketing)
    print(meth_stat)
    if(length(meth_stat) < 1) {meth_stat <- NA}
    if (length(meth_stat) > 1) {meth_stat <- meth_stat[1]}
    if (is.nan(meth_stat)) {meth_stat <- NA}
  }
  return(meth_stat)
}

############################################################

#' Given an expression and variable indicating if we are looking tumor-normal-
#' matched, extracts and processes the expression or fold change expression 
#' @param expression_df an expression data frame of interest
#' @param sample the sample ID of the current sample
#' @param tumNormMatched a TRUE/FALSE value indicating whether or not the analysis is 
#' tumor-normal matched
get_exp_stat <- function(expression_df, sample, tumNormMatched) {
  if(!tumNormMatched) {
    # Add a pseudocount so we don't take the log of 0 (1 is typically used for logged values)
    exp_stat <- log2(as.numeric(unlist(expression_df[, grepl(sample, colnames(expression_df)), 
                                                     with = FALSE])) + 1)
  } else {
    exp_stat <- get_tnm_fc(expression_df, sample)
  }
  
  if (length(exp_stat) < 1) {exp_stat <- NA}
  if (length(exp_stat > 1)) {exp_stat <- mean.default(exp_stat)}  # if there are two ENSG ids, avg them
  if(is.infinite(exp_stat) | is.nan(exp_stat)) {exp_stat <- NA}
  
  return(exp_stat)
}


############################################################

### NOTE: FOR FOR TUMOR-NORMAL MATCHED ONLY ###
#' Given a methylation DF, a variable to denote if we are bucketing
#' methylation, and a sample j, this function subsets this data frame 
#' to the given sample j and returns the differential methylation value
#' (bucketed or not). Buckets:
#' Bucket 1: Lower methylation in tumor (meth_logfc < -0.5)
#' Bucket 2: Similar methylation in tumor and normal (-0.5 <= meth_logfc < 0.5)
#' Bucket 3: Higher methylation in tumor (meth_logfc >= 0.5)
#' @param methylation_df a methylation DF, subsetted to a particular target gene k 
#' @param sample the four-character patient ID for given sample j
#' @param meth_bucketing a T/F value indicating if we are bucketing 
#' methylation values. 
get_tnm_methylation <- function(methylation_df, sample, meth_bucketing) {
  
  meth_logfc <- get_tnm_fc(methylation_df, sample)
  
  # If we are bucketing differential methylation, do so here before we return
  if(meth_bucketing) {
    lower_bound <- -0.5
    higher_bound <- 0.5
    
    # Determine the bucket that this falls into
    if(meth_logfc < lower_bound) {return(c(1,0,0))}
    else if (meth_logfc > upper_bound) {return(c(0,0,1))}
    else {return(c(0,1,0))}
    
  } else {return(meth_logfc)}
}

############################################################
#' A generic function to calculate the tumor-normal-matched fold change in 
#' either expression or methylation, given the appropriate data frame and 
#' the current sample ID. If there are multiple normal samples for a given
#' patient, it averages across them.
#' @param df the relevant expression or methylation DF
#' @param sample the sample ID of the sample in question
get_tnm_fc <- function(df, sample) {
  
  # Use the sample ID to subset the data frame
  val <- df[, grepl(sample, colnames(df))]
  
  # Get any corresponding normal samples as well; start by getting the patient ID from the sample
  patient <- unlist(strsplit(sample, "-", fixed = TRUE))[1]
  patient_n <- paste(patient, "-1", sep = "")
  val <- cbind(val, df[, grepl(patient_n, colnames(df))])
  if (colnames(val)[1] == "val") {colnames(val)[1] <- sample}
  
  # Create our empty fold change variable
  fc <- NA
  
  # If we have no patient matches or only 1 unmatched result, set to NA (no tumor-normal 
  # matched data)
  if(length(val) == 0 | length(val) == 1) {return(NA)}
  
  # If it's not a data frame, this was also a sign we had only one column - set 
  # to NA (no tumor-normal matched data)
  else if(!is.data.frame(val)) {return(NA)}
  
  # If there are 2 columns (one for tumor, one for normal)
  else if(length(val) == 2) {
    
    # If there's only 1 row (one ENSG match)
    if(nrow(val) == 1) {
      # Extract the value for this gene in this patient for both tumor and normal
      val_tumor <- as.numeric(val[,grepl("-0", colnames(val))])
      val_norm <- as.numeric(val[,grepl("-1", colnames(val))])
      
      # Calculate the fold change
      fc <- val_tumor / val_norm
      
      # If there's >1 row (multiple ENSG matches)
    } else if (nrow(val) > 1) {
      # Get the mean of the values across these different ENSG hits
      val_tumor <- mean.default(as.numeric(unlist(val[,grepl("-0", colnames(val))])))
      val_norm <- mean.default(as.numeric(unlist(val[,grepl("-11", colnames(val))])))
      
      # Calculate the fold change
      fc <- val_tumor / val_norm
    } else {fc <- NA}
  }
  # If we have more than one normal sample; take the average across them
  else if(length(val) > 2) {
    # Get the mean of the values across these different normal samples
    normal_vals <- c()
    val_norm <- val[,grepl("-1", colnames(val))]
    for (i in 1:ncol(val_norm)) {
      normal_vals <- c(normal_vals, as.numeric(val_norm[,i]))
    }
    val_norm <- mean.default(as.numeric(normal_vals))
    
    # Calculate the fold change
    fc <- val_tumor / val_norm
    
  } else {print(paste("Unhandled case:", sample))}
  
  print(paste("Fold change:", fc))
  
  # Take the log of the fold change, with a pseudocount
  logfc <- log2(fc + 1)
  
  return(logfc)
}


############################################################

#' Takes an expression DF and a starter DF and returns the starter DF without outlier patients,
#' defined as patients with expression of target gene that exceeds 1.5*IQR + Q3 threshold
#' (separately among patients with a mutation in the regulatory protein and patients without)
#' @param expression_df an expression DF, unfiltered 
#' @param starter_df a starter DF for the regulatory protein of interest i
#' @param ensg the ensembl ID of the target gene k
filter_expression_df <- function(expression_df, starter_df, ensg) {
  patients_mut <- rownames(starter_df[starter_df$MutStat_i == 1,])
  patients_nonmut <- setdiff(rownames(starter_df), patients_mut)
  
  # Filter expression DF to target gene of interest and crop patient IDs
  expression_df_filt <- data.frame()
  if(length(ensg) == 1) {
    expression_df_filt <- expression_df[rownames(expression_df) == ensg,]   # ENSG of target gene
  } else {
    rows_exp_filt <- lapply(ensg, function(e) expression_df[rownames(expression_df) == e,])
    expression_df_filt <- do.call("rbind", rows_exp_filt)
  }
  labels <- unlist(lapply(colnames(expression_df_filt), function(x) 
    unlist(strsplit(x, "-", fixed = TRUE))[3]))
  expression_df_upd_labels <- expression_df_filt
  colnames(expression_df_upd_labels) <- labels
  
  # Get the expression values for patients in each mutational group (mutated, non-mutated)
  expression_mut <- list()
  expression_nonmut <- list()
  for (i in 1:ncol(expression_df_upd_labels)) {
    patient <- colnames(expression_df_upd_labels)[i]
    if(patient %fin% patients_mut) {
      expression_mut[[patient]] <- as.numeric(expression_df_upd_labels[,i])
    } else {
      expression_nonmut[[patient]] <- as.numeric(expression_df_upd_labels[,i])
    }
  }
  print(expression_mut)
  print(expression_nonmut)
  
  # Get the IQR, Q3, and expression thresholds for each group
  iqr_mut <- IQR(unlist(expression_mut), na.rm = TRUE)
  iqr_nonmut <- IQR(unlist(expression_nonmut), na.rm = TRUE)
  q3_mut <- as.numeric(quantile(unlist(expression_mut), na.rm = TRUE)[4])
  q3_nonmut <- as.numeric(quantile(unlist(expression_nonmut), na.rm = TRUE)[4])
  threshold_mut <- as.numeric(q3_mut + (1.5 * iqr_mut))
  print(paste("threshold mut:", threshold_mut))
  threshold_nonmut <- as.numeric(q3_nonmut + (1.5 * iqr_nonmut))
  print(paste("threshold nonmut:", threshold_nonmut))
  # TODO: IF THERE ARE MULTIPLE ROWS, SHOULD I TAKE THE IQR SEPARATELY FOR EACH? 
  
  # Get patients that do not exceed threshold 
  if(length(ensg) > 1) {
    filt_express_mult <- function(expression, threshold) {
      filt_expression <- lapply(expression, function(x) {
        if(!any(x > threshold)) {return(x)}})
      filt_expression <- removeNULL(filt_expression)
      return(filt_expression)
    }
    filt_expression_mut <- filt_express_mult(expression_mut, threshold_mut)
    filt_expression_nonmut <- filt_express_mult(expression_nonmut, threshold_nonmut)
  } else {
    filt_expression_mut <- expression_mut[expression_mut <= threshold_mut]
    filt_expression_nonmut <- expression_nonmut[expression_nonmut <= threshold_nonmut]
  }
  patients_to_keep <- c(names(filt_expression_mut), names(filt_expression_nonmut))
  
  # Filter the starter DF based on these patients
  starter_df_filt <- starter_df[rownames(starter_df) %fin% patients_to_keep,]
  
  # Return the filtered starter DF
  print(starter_df_filt)
  print(paste("Nrow original starter", nrow(starter_df)))
  print(paste("Nrow filtered starter", nrow(starter_df_filt)))
  return(starter_df_filt)
}
