############################################################
### Dyscovr Helper Functions
### PUBLICATION INFORMATION
############################################################

library.path <- .libPaths()

#library(data.table, lib.loc = library.path, quietly = T)
#library(dqrng, lib.loc = library.path, quietly = T)
#library(stringr, lib.loc = library.path, quietly = T)

# This file contains a set of helper functions that will be called by the 
# create_dyscovr_input_tables.R or dyscovr.R while they are running.

############################################################

#' Filters a mutation DF in order to keep only rows that overlap the given 
#' genes of interest
#' @param mutation_targ_df a mutation target DF to be subsetted
#' @param uniprot_ids a list of uniprot IDs to keep in the DF
filter_mut_by_uniprot <- function(mutation_targ_df, uniprot_ids) {
  mutation_rows_to_keep <- unique(unlist(lapply(uniprot_ids, function(x) {
    return(grep(x, mutation_targ_df$Swissprot))
  })))
  mutation_targ_df <- mutation_targ_df[mutation_rows_to_keep,]
  return(mutation_targ_df)
}

############################################################

#' Filters a CNA DF in order to keep only rows that overlap the given
#' genes of interest
#' @param cna_df a copy number alteration data frame to be subsetted
#' @param ensg_ids a list of ensembl IDs to keep in the DF
filter_cna_by_ensg <- function(cna_df, ensg_ids) {
  cna_rows_to_keep <- unique(unlist(lapply(ensg_ids, function(x) 
    grep(x, cna_df$ensg_id))))
  cna_df <- cna_df[cna_rows_to_keep,]
  return(cna_df)
}

############################################################

#' Filters a methylation DF in order to keep only rows that overlap the given
#' genes of interest
#' @param methylation_df a methylation data frame to be subsetted
#' @param ensg_ids a list of ensembl IDs to keep in the DF
filter_meth_by_ensg <- function(methylation_df, ensg_ids) {
  meth_rows_to_keep <- unique(unlist(lapply(ensg_ids, function(x) 
    grep(x, methylation_df$ensg_ids))))
  methylation_df <- methylation_df[meth_rows_to_keep,]
  return(methylation_df)
}


############################################################

#' Given a parameter for if/how we are bucketing CNA, extracts and processes the 
#' CNA value from the CNA data frame
#' @param cna_df a CNA data frame of interest
#' @param sample the sample ID of the current sample
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. 
#' Possible values are "bucket_inclAmp", "bucket_exclAmp", "bucket_justAmp",
#' "bucket_justDel", and "rawCNA"
#' @param dataset either 'TCGA' or 'METABRIC'
#' @param use_relative a T/F value indicating whether or not we are
#' taking into account the mean CN value for a given patient 
#' (e.g. log2(CN+1 / N+1) rather than log2(CN+1))
#' @param cna_vals_sample if use_relative is T, this is the full set of copy 
#' numbers across all genes for the given sample, in order to calculate the mean 
#' CN for that sample; otherwise NA
get_cna_stat <- function(cna_df, sample, cna_bucketing, dataset, use_relative, 
                         cna_vals_sample) {
  
  cna_stat <- unique(as.integer(unlist(cna_df[,colnames(cna_df) == sample, 
                                              with = F])))
  
  # "Raw" CNA: Log2(CNA Stat + 1) or Log2(CNA Stat + 1 / Mean CNA + 1)
  if(cna_bucketing == "rawCNA") {

    if (length(cna_stat) > 1) {cna_stat <- cna_stat[1]}
    if(length(cna_stat) < 1) {cna_stat <- NA}
    if (is.nan(cna_stat)) {cna_stat <- NA}
    
    if(use_relative) {
      # Get rid of the top 10% and bottom 10% of CNA values for this sample
      ten_perc <- 0.1 * length(cna_vals_sample)
      cna_vals_sample <- sort(cna_vals_sample)
      cna_vals_sample_sub <- cna_vals_sample[
        ten_perc:(length(cna_vals_sample)-ten_perc)]
      
      # Find the mean of the remaining values
      mean_cn_samp <- mean(cna_vals_sample_sub)
      
      # Calculate the new cna_stat using this value
      cna_stat <- log2((cna_stat + 1) / (mean_cn_samp + 1))
      
    } else {
      cna_stat <- log2(cna_stat + 1)
    }
    
  # 1 if amplified ("justAmp") or deleted ("justDel"), 0 if 2+ copies
  } else if ((cna_bucketing == "bucket_justAmp") | 
             (cna_bucketing == "bucket_justDel")) {
    
    if (length(cna_stat) > 1) {cna_stat <- cna_stat[1]}
    
    if(length(cna_stat) < 1) {cna_stat <- NA}
    else if (is.nan(cna_stat) | is.na(cna_stat)) {cna_stat <- NA}
    else {
      if (cna_bucketing == "bucket_justAmp") {
        if(dataset == 'TCGA') {
          if(cna_stat > 2) {cna_stat <- 1}
          else {cna_stat <- 0}
        } else if (dataset == "METABRIC") {
          if(cna_stat > 0) {cna_stat <- 1}
          else {cna_stat <- 0}
        } else {
          print("Currently only implemented for TCGA and METABRIC.")
          return(NA)
        }
      } else {
        if(dataset == "TCGA") {
          # Include partial deletions
          if(cna_stat < 2) {cna_stat <- 1}
          # Include only full deletions
          else {cna_stat <- 0}
          
        } else if (dataset == "METABRIC") {
          # Include partial deletions
          if(cna_stat < 0) {cna_stat <- 1}
          # Include only full deletions
          else {cna_stat <- 0}
        } else {
          print("Currently only implemented for TCGA and METABRIC.")
          return(NA)
        }
      }
    }
    
  # 3 buckets for deleted, normal, and amplified
  } else {
    
    if(dataset == "TCGA") {
      
      str_buckets <- as.character(unlist(cna_df[,colnames(cna_df) == sample, 
                                                with = F]))
      if(!(length(str_buckets) == 0)) {
        buckets_list <- unlist(strsplit(unlist(strsplit(
          str_buckets, "(", fixed = T))[2], ")", fixed = T))[1]
        cna_stat <- list(as.integer(unlist(strsplit(buckets_list, ",", 
                                                    fixed = T))))
      } else {cna_stat <- list(NA, NA, NA)}
      
    } else if (dataset == "METABRIC") {
      
      if (length(cna_stat) > 1) {cna_stat <- cna_stat[1]}
      if(length(cna_stat) < 1) {cna_stat <- list(NA, NA, NA)}
      else if (is.nan(cna_stat) | is.na(cna_stat)) {cna_stat <- list(NA, NA, NA)}
      else {
        if(cna_stat < 0) {cna_stat <- list(1,0,0)}
        else if (cna_stat > 0) {cna_stat <- list(0,0,1)}
        else if (cna_stat == 0) {cna_stat <- list (0,1,0)}
        else {print("Invalid CNA value. Returning NA."); return(NA)}
      }
    } else {
      print("Currently only implemented for TCGA and METABRIC.")
      return(NA)
    }
  }  
  return(cna_stat)
}


############################################################

#' Given a mutation target DF and a sample of interest, extracts the 0 or 1 
#' mutation statistic to indicate non-mutated/mutated
#' @param mutation_targ_df a mutation data frame of interest
#' @param sample the sample ID of the current sample
get_mut_stat_targ <- function(mutation_targ_df, sample) {
  mut_count <- as.numeric(unlist(mutation_targ_df[, grepl(
    sample, colnames(mutation_targ_df)), with = F]))
  
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

#' Given a parameter for if we are bucketing methylation, extracts and processes 
#' the methylation value from the methylaton data frame
#' @param methylation_df a methylation data frame of interest
#' @param sample the sample ID of the current sample
#' @param meth_bucketing a T/F value indicating if we are bucketing 
#' methylation values. 
get_meth_stat <- function(methylation_df, sample, meth_bucketing) {
  
  if(meth_bucketing) { 
    meth_buckets <- methylation_df[,grepl(sample, colnames(methylation_df)), 
                                   with = F]
    str_buckets <- as.character(unlist(meth_buckets))
    if(length(str_buckets) > 1) {
      if(length(unique(str_buckets)) == 1) {str_buckets <- unique(str_buckets)}
      else {str_buckets <- str_buckets[1]}
    }

   if(!(length(str_buckets) == 0)) {
      buckets_list <- unlist(strsplit(unlist(strsplit(
        str_buckets, "(", fixed = T))[2], ")", fixed = T))[1]
      meth_stat <- list(as.integer(unlist(strsplit(buckets_list, ",", 
                                                   fixed = T))))
    } else {meth_stat <- list(NA, NA, NA)}

  } else {  
    meth_stat <- unlist(methylation_df[,grepl(sample, colnames(methylation_df)), 
                                       with = F])
    if (length(meth_stat) > 1) {meth_stat <- mean.default(as.numeric(meth_stat))}
    if(length(meth_stat) < 1) {meth_stat <- NA}
    if (is.nan(meth_stat)) {meth_stat <- NA}
  }
  
  return(as.numeric(meth_stat))
}

############################################################

#' Given an expression DF and a given cancer sample, extracts and processes the 
#' expression or fold change expression for that sample and protein-of-interest
#' @param expression_df an expression data frame of interest
#' @param sample the sample ID of the current sample
#' @param dataset either 'TCGA' or 'METABRIC'
#' @param log_expression a T/F value indicating whether or not we need to 
#' take the log2(exp + 1), or whether the expression values have already been 
#' transformed to a normal distribution
get_exp_stat <- function(expression_df, sample, dataset, log_expression) {
  exp_stat <- as.numeric(unlist(expression_df[, grepl(
    sample, colnames(expression_df)), with = F]))
  
  if((dataset == "TCGA") & log_expression) {
    # Add a pseudocount so we don't take the log of 0 (1 is typically used for 
    # logged values)
    exp_stat <- log2(exp_stat + 1)
  }

  if (length(exp_stat) < 1) {exp_stat <- NA}
  # if there are two ENSG IDs, average them
  if (length(exp_stat > 1)) {exp_stat <- mean.default(exp_stat)}  
  if(is.infinite(exp_stat) | is.nan(exp_stat)) {exp_stat <- NA}
  
  return(exp_stat)
}


############################################################

#' Create the appropriate output file path (add the cancer type, the 
#' protein-of-interest group or tester protein name, etc)
#' @param cancerType the name of the type of cancer, or "PanCancer" for all 
#' cancers
#' @param specificType if we are looking pan-cancer, but for specific cancer 
#' types individually, the specific cancer type we are looking at is given
#' @param test a T/ F value indicating whether this is a test on a single
#' regulatory protein
#' @param tester_name if this is a test, the name of the tester regulatory protein
#' @param run_name string value indicating the type of run if we are not running 
#' a test on a single regulatory protein 
#' @param dataset either 'TCGA' or 'METABRIC'
#' @param outpath a starter outpath from which to build
create_lm_input_file_outpath <- function(cancerType, specificType, test, 
                                         tester_name, run_name, dataset, 
                                         outpath) {
  
  outpath <- paste(outpath, "input_files", sep = "/")
  
  if(dataset != "TCGA") {outpath <- paste(outpath, dataset, sep = "/")} 
  else {
    outpath <- paste(outpath, cancerType, sep = "/")
    if(!((specificType == "") | (specificType == "ALL"))) {
      outpath <- paste(outpath, specificType, sep = "/")
    }
  }
  outpath <- paste(outpath, "lm_input_tables", sep = "/")
  
  if(test) {
    outpath <- paste(outpath, tester_name, sep = "/")
  } else {outpath <- paste(outpath, run_name, sep = "/")}
  
  print(paste("File Outpath:", outpath))
  
  return(outpath)
}

############################################################

#' Create the appropriate output file name for the Dyscovr DF result based on the
#' conditions of the given run
#' @param test a T/ F value indicating whether this is a test on a single
#' protein-of-interest/ driver gene
#' @param tester_name if this is a test, the name of the tester protein
#' @param run_name string value indicating the type of run, if we are not running 
#' a test on a single protein
#' @param cna_bucketing a string description of how we are bucketing/ handling 
#' CNA values
#' @param meth_bucketing a T/ F value indicating if we are in some way 
#' bucketing our methylation value
#' @param meth_type 'Beta' or 'M' depending on the type of methylation we are
#' using
#' @param patient_df_name the string name used to import the patient DF
#' @param num_pcs an integer value from 0 to 2 indicating the number of principal
#' components we are using
#' @param randomize a T/ F value indicating whether or not we are 
#' randomizing our dependent variable (expression)
#' @param patients_to_incl_label a label indicating the subpopulation we are 
#' restricting to, if any
#' @param removeCis a T/ F value indicating whether or not we have 
#' removed cis gene pairings 
#' @param removeMetastatic a T/ F value indicating whether or not we have 
#' removed metastatic samples from the analysis
create_driver_input_table_filename <- function(test, tester_name, run_name, 
                                                cna_bucketing, meth_bucketing, 
                                                meth_type, patient_df_name, 
                                                num_pcs, patients_to_incl_label, 
                                                removeCis, removeMetastatic) {
  outfn <- "driver_input"
  
  # Add the name of the patient population of interest, if there is one
  if(!(patients_to_incl_label == "")) {
    outfn <- paste(outfn, patients_to_incl_label, sep = "_")
  }
  
  # Add the name of the driver protein/ group of driver proteins
  if(test) {
    outfn <- paste(outfn, tester_name, sep = "_")
  } else {outfn <- paste(outfn, run_name, sep = "_")} 
  
  # Add the decisions for CNA bucketing, mutation restriction, methyl. bucketing
  cna_bucketing_lab <- "rawCNA"
  if (cna_bucketing != "rawCNA") {cna_bucketing_lab <- paste0("CNA", cna_bucketing)}
  
  if(!(is.na(meth_type) | (meth_type == "NA"))) {
    meth_b_lab <- "Raw"
    if(meth_bucketing) {meth_b_lab <- "Bucketed"}
    meth_label <- paste0("meth", paste0(meth_type, meth_b_lab))
    outfn <- paste0(outfn, paste0(cna_bucketing_lab, meth_label))
    
  } else {
    outfn <- paste(outfn, cna_bucketing_lab, sep = "_")
  }
  
  # Add PCs if needed
  if(!(num_pcs == 0))  {outfn <- paste(outfn, paste0(num_pcs, "PCs"), sep = "_")}
  
  # If we are removing cis comparisons, add a label to denote this
  if(removeCis) {outfn <- paste(outfn, "RmCis", sep = "_")}
  
  # If we are NOT removing metastatic samples, add a label to denote this
  if(!removeMetastatic) {outfn <- paste(outfn, "notRmMetast", sep = "_")}
  
  print(paste("Driver Outfile Name", outfn))
  
  return(outfn)
}


############################################################

#' Create the appropriate output file name for the master DF result based on the
#' conditions of the given run
#' @param target_uniprot a Uniprot ID of the target gene for the given DF
#' @param expression_df_name the used to import the expression DF
#' @param cna_bucketing a string description of how we are bucketing/ handling 
#' CNA values
#' @param meth_bucketing a T/ F value indicating if we are in some way 
#' bucketing our methylation value
#' @param meth_type 'Beta' or 'M' depending on the type of methylation we are
#' using
#' @param patient_df_name the string name used to import the patient DF
#' @param randomize a T/ F value indicating whether or not we are 
#' randomizing our dependent variable (expression or methylation)
#' @param patients_to_incl_label a label indicating the subpopulation we are 
#' restricting to, if any
#' @param removeCis a T/ F value indicating whether or not we have 
#' removed cis gene pairings 
#' @param removeMetastatic a T/ F value indicating whether or not we have 
#' removed metastatic samples from the analysis
create_lm_input_table_filename <- function(target_uniprot, expression_df_name, 
                                           cna_bucketing, meth_bucketing, 
                                           meth_type, patient_df_name, randomize, 
                                           patients_to_incl_label, 
                                           removeCis, removeMetastatic) {
  outfn <- "lm_input"
  
  # Add the name of the patient population of interest, if there is one
  if(!(patients_to_incl_label == "")) {
    outfn <- paste(outfn, patients_to_incl_label, sep = "_")
  }
  
  # Add the name of the target 
  outfn <- paste(outfn, target_uniprot[1], sep = "_")
  
  # Add the decisions for expression normalization, CNA bucketing, mutation 
  # restriction, methylation bucketing
  expr_norm <- unlist(strsplit(expression_df_name, "_", fixed = T))[2]
  
  cna_bucketing_lab <- "rawCNA"
  if (cna_bucketing != "rawCNA") {cna_bucketing_lab <- paste0("CNA", cna_bucketing)}
  
  if(!(is.na(meth_type) | (meth_type == "NA"))) {
    meth_b_lab <- "Raw"
    if(meth_bucketing) {meth_b_lab <- "Bucketed"}
    meth_label <- paste0("meth", paste0(meth_type, meth_b_lab))
    
    outfn <- paste(outfn, paste(expr_norm, paste(cna_bucketing_lab, meth_label, 
                                                 sep = "_"), sep = "_"), sep = "_")
  } else {
    outfn <- paste(outfn, paste(expr_norm, cna_bucketing_lab, sep = "_"), sep = "_")
  }
  
  # If we are removing cis comparisons, add a label to denote this
  if(removeCis) {outfn <- paste(outfn, "RmCis", sep = "_")}
  
  # If we randomized expression, add "_RANDOMIZED" to the end of the file name
  if (randomize) {outfn <- paste(outfn, "_RANDOMIZED", sep = "")}
  
  print(paste("Outfile Name", outfn))
  
  return(outfn)
}


############################################################

#' Create the appropriate output file path (add the cancer type, the driver 
#' protein group, etc)
#' @param cancerType the name of the type of cancer, or "PanCancer" for all 
#' cancers
#' @param specificType if we are looking pan-cancer, but for specific cancer 
#' types individually, the specific cancer type we are looking at is given
#' @param test a T/ F value indicating whether this is a test on a single
#' regulatory protein
#' @param tester_name if this is a test, the name of the tester regulatory 
#' protein
#' @param run_name string value indicating the type of run if we are not running 
#' a test on a single regulatory protein 
#' @param dataset either 'TCGA' or 'METABRIC'
#' @param outpath a starter outpath from which to build
create_file_outpath <- function(cancerType, specificType, test, tester_name, 
                                run_name, dataset, outpath) {
  outpath <- paste(outpath, "output_files", sep = "/")
   
  if(dataset != "TCGA") {outpath <- paste(outpath, dataset, sep = "/")}
  outpath <- paste(outpath, cancerType, sep = "/")
  
  if(test) {
    outpath <- paste(outpath, tester_name, sep = "/")
  } else {outpath <- paste(outpath, run_name, sep = "/")}
  
  if(!(specificType == "")) {outpath <- paste(outpath, specificType, sep = "/")}
  
  print(paste("File outpath:", outpath))
  
  return(outpath)
}


############################################################

#' Create the appropriate output file name for the Dyscovr DF result based on the
#' conditions of the given run
#' @param test a T/ F value indicating whether this is a test on a single
#' protein-of-interest/ driver
#' @param tester_name if this is a test, the name of the tester protein
#' @param run_name string value indicating the type of run if we are not running 
#' a test on a single tester protein
#' @param targets_name a string name of the category of targets we are testing
#' @param expression_df_name the string name used to import the expression DF
#' @param cna_bucketing a string description of how we are bucketing/ handling 
#' CNA values
#' @param meth_bucketing a T/ F value indicating if we are in some way 
#' bucketing our methylation value
#' @param meth_type 'Beta' or 'M' depending on the type of methylation we are
#' using
#' @param patient_df_name the string name used to import the patient DF
#' @param num_pcs an integer value from 0 to 3 indicating the number of 
#' principal components we are using
#' @param randomize a T/ F value indicating whether or not we are 
#' randomizing our dependent variable (expression or methylation)
#' @param covs_to_incl_label a label indicating what covariates we've included 
#' in this run
#' @param patients_to_incl_label a label indicating the subpopulation we are 
#' restricting to, if any
#' @param model_type either "linear" or "mixed" (to indicate a linear mixed 
#' effects model)
#' @param removeCis a T/ F value indicating whether or not we have 
#' removed cis gene pairings 
#' @param removeMetastatic a T/ F value indicating whether or not we have 
#' removed metastatic samples from the analysis
#' @param regularization either "None", "L1" or "Lasso", or "L2" or "Ridge" to 
#' indicate if regularization was performed
#' @param signif_eval_type if regularization is not "None", provides additional 
#' info on how we generated significance metrics
create_output_filename <- function(test, tester_name, run_name, targets_name, 
                                   expression_df_name, cna_bucketing, 
                                   meth_bucketing, meth_type, patient_df_name, 
                                   num_pcs, randomize, covs_to_incl_label, 
                                   patients_to_incl_label, model_type, removeCis, 
                                   removeMetastatic, regularization, 
                                   signif_eval_type) {
  outfn <- "res"
  
  # Add the name of the patient population of interest, if there is one
  if(!(patients_to_incl_label == "")) {
    outfn <- paste(outfn, patients_to_incl_label, sep = "_")
  }
  
  # Add the name of the driver gene/ group of driver genes
  if(test) {outfn <- paste(outfn, tester_name, sep = "_")}
  else {outfn <- paste(outfn, run_name, sep = "_")} 
  
  # Add the name of the target group
  outfn <- paste(outfn, targets_name, sep = "_")
  
  # Add the decisions for expression normalization, CNA bucketing, mutation 
  # restriction, methylation bucketing, etc 
  expr_norm <- unlist(strsplit(expression_df_name, "_", fixed = T))[2]
  
  cna_bucketing_lab <- "rawCNA"
  if (cna_bucketing != "rawCNA") {cna_bucketing_lab <- paste0("CNA", cna_bucketing)}
  
  if(!(is.na(meth_type) | (meth_type == "NA"))) {
    meth_b_lab <- "Raw"
    if(meth_bucketing) {
      meth_b_lab <- "Bucketed"
    }
    meth_label <- paste0("meth", paste0(meth_type, meth_b_lab))
    
    outfn <- paste(outfn, paste(expr_norm, paste(cna_bucketing_lab, meth_label, 
                                                 sep = "_"), sep = "_"), sep = "_")
  } else {
    outfn <- paste(outfn, paste(expr_norm, cna_bucketing_lab, sep = "_"), sep = "_")
  }
  
  # Add PCs if needed
  if(!(num_pcs == 0))  {outfn <- paste(outfn, paste0(num_pcs, "PCs"), sep = "_")}
  
  # If we are removing cis comparisons, add a label to denote this
  if(removeCis) {outfn <- paste(outfn, "RmCis", sep = "_")}
  
  # If we are NOT removing metastatic samples, add a label to denote this
  if(!removeMetastatic) {outfn <- paste(outfn, "notRmMetast", sep = "_")}
  
  # Linear or mixed model
  if(!model_type == "linear") {outfn <- paste(outfn, "mixed", sep = "_")}
  
  # If we have regularized in some way, denote this
  if(!(regularization == "None")) {
    if (regularization %fin% c("L1", "lasso")) {
      if(grepl("randomization", signif_eval_type)) {
        if(signif_eval_type == "randomization_predictors") {
          outfn <- paste(outfn, "L1rand_pred", sep = "_")
        } else {outfn <- paste(outfn, "L1rand", sep = "_")}
      }
      else if(signif_eval_type == "subsampling") {
        outfn <- paste(outfn, "L1subsamp", sep = "_")
      }
      else if (signif_eval_type == "selectiveInference") {
        outfn <- paste(outfn, "L1si", sep = "_")
      }
      else {
        print("Error; unrecognized argument for method of evaluating LASSO 
              significance.")
        outfn <- paste(outfn, "L1", sep = "_")
      }
    }
    if (regularization %fin% c("L2", "ridge")) {
      outfn <- paste(outfn, "L2", sep = "_")
    }
    if (regularization == "bayesian.bgl") {
      outfn <- paste(outfn, "bayesian.bgl", sep = "_")
    }
    if (regularization == "bayesian.bglss") {
      outfn <- paste(outfn, "bayesian.bglss", sep = "_")
    }
  }
  
  # If we've restricted the number of covariates, add that label
  if(!(covs_to_incl_label == "")) {
    outfn <- paste(outfn, covs_to_incl_label, sep = "_")
  }
  
  # Add "uncorrected" to the file name so we know MHT correction has not yet 
  # been performed on these results
  outfn <- paste0(outfn, "_uncorrected")
  
  # If we randomized expression, add "_RANDOMIZED" to the end of the file name
  if (randomize) {outfn <- paste0(outfn, "_RANDOMIZED")}
  
  print(paste("Outfile Name", outfn))
  
  return(outfn)
}


############################################################

#' Short function for combining the output tibbles produced for each row in 
#' list_of_rows; used to recombine output of foreach, which are given as 
#' multiResultClass objects
#' @param x_obj the first output object of results
#' @param y_obj the second output object of results
#' @param expected_col_num the expected number of columns, in case we need to 
#' turn 2 NA inputs into a table accepted in the re-combining process
combine_res_obj <- function(x_obj, y_obj, expected_col_num = 7) {
  
  # First, recombine the summary tables
  x_tab <- x_obj$summary_table
  y_tab <- y_obj$summary_table

  comb_dt <- combine_summary_tables(x_tab, y_tab)

  # Next, recombine the removed variables object
  x_vars_rm <- x_obj$variables_rm
  y_vars_rm <- y_obj$variables_rm

  comb_vars_rm_dt <- combine_variables_rm_tables(x_vars_rm, y_vars_rm)

  return(multiResultClass(comb_dt, comb_vars_rm_dt))
}


#' Helper function for recombining summary tables
#' @param x_tab summary table from first target
#' @param y_tab summary_table from second target
combine_summary_tables <- function(x_tab, y_tab) {
  # Check if length of either table is 0; if so, handle like an NA table
  if((length(x_tab) == 0) | (length(y_tab) == 0)) {
    if (length(x_tab) == 0) {
      x_tab <- as.data.table(t(rep(NA, times = expected_col_num)))
    }
    if (length(y_tab) == 0) {
      y_tab <- as.data.table(t(rep(NA, times = expected_col_num)))
    }
  }
  # Check if the length is 1 (if so, it is NA)
  else if(length(x_tab) == 1 | length(y_tab) == 1) {
    # Check if either table is NA
    if (is.na(x_tab)) {
      if(is.na(y_tab)) {
        x_tab <- as.data.table(t(rep(NA, times = expected_col_num)))
        y_tab <- as.data.table(t(rep(NA, times = expected_col_num)))
      } else {
        x_tab <- as.data.table(t(rep(NA, times = ncol(y_tab))))
        colnames(x_tab) <- colnames(y_tab)
      }
    } else {
      if(is.na(y_tab)) {
        y_tab <- as.data.table(t(rep(NA, times = ncol(x_tab))))
        colnames(y_tab) <- colnames(x_tab)
      }
    }
  }
  # Otherwise, the length is greater than 1
  else {
    # Remove any NA or NAN values
    if(T %fin% c(is.na(x_tab), is.na(y_tab))) {
      x_tab <- na.omit(x_tab)
      y_tab <- na.omit(y_tab)
    }
    
    # Again, check if removing NAs/NaNs has set the length of either table is 0; 
    # if so, handle like an NA table
    if((length(x_tab) == 0) | (length(y_tab) == 0)) {
      if (length(x_tab) == 0) {
        x_tab <- as.data.table(t(rep(NA, times = expected_col_num)))
      }
      if (length(y_tab) == 0) {
        y_tab <- as.data.table(t(rep(NA, times = expected_col_num)))
      }
      
    } else {
      # Check if the number of rows is 0
      #print(dim(x_tab))
      #print(dim(y_tab))
      if(x_tab[, .N] == 0) {
        x_tab <- rbind(x_tab, as.data.table(t(rep(NA, times = ncol(x_tab)))), 
                       use.names = F)
      }
      if(y_tab[, .N] == 0) {
        y_tab <- rbind(y_tab, as.data.table(t(rep(NA, times = ncol(y_tab)))), 
                       use.names = F)
      }
    }
  }

  comb_dt <- tryCatch({
    rbindlist(list(x_tab, y_tab), use.names = T, fill = T)
  },error=function(cond){
    print(cond)
    return(x_tab)
  })
  return(comb_dt)
}


#' Helper function for recombining tables of variables that have been removed
#' @param x_vars_rm variables removed table from first target
#' @param y_vars_rm variables removed table from second target
#' @param expected_col_num the number of expected columns (2)
combine_variables_rm_tables <- function(x_vars_rm, y_vars_rm, expected_col_num = 2) {
  # Check if length of either table is 0; if so, handle like an NA table
  if((length(x_vars_rm) == 0) | (length(y_vars_rm) == 0)) {
    if (length(x_vars_rm) == 0) {
      x_vars_rm <- as.data.table(t(rep(NA, times = expected_col_num)))
    }
    if (length(y_vars_rm) == 0) {
      y_vars_rm <- as.data.table(t(rep(NA, times = expected_col_num)))
    }
  }
  # Check if the length is 1 (if so, it is NA)
  else if(length(x_vars_rm) == 1 | length(y_vars_rm) == 1) {
    # Check if either table is NA
    if (is.na(x_vars_rm)) {
      if(is.na(y_vars_rm)) {
        x_vars_rm <- as.data.table(t(rep(NA, times = expected_col_num)))
        y_vars_rm <- as.data.table(t(rep(NA, times = expected_col_num)))
      } else {
        x_tab <- as.data.table(t(rep(NA, times = ncol(y_tab))))
        colnames(x_vars_rm) <- colnames(y_vars_rm)
      }
    } else {
      if(is.na(y_vars_rm)) {
        y_vars_rm <- as.data.table(t(rep(NA, times = ncol(x_vars_rm))))
        colnames(y_vars_rm) <- colnames(x_vars_rm)
      }
    }
  }
  # Otherwise, the length is greater than 1
  else {
    # Remove any NA or NAN values
    if(T %fin% c(is.na(x_vars_rm), is.na(y_vars_rm))) {
      x_vars_rm <- na.omit(x_vars_rm)
      y_vars_rm <- na.omit(y_vars_rm)
    }
    
    # Again, check if removing NAs/NaNs has set the length of either table is 0; 
    # if so, handle like an NA table
    if((length(x_vars_rm) == 0) | (length(y_vars_rm) == 0)) {
      if (length(x_vars_rm) == 0) {
        x_vars_rm <- as.data.table(t(rep(NA, times = expected_col_num)))
      }
      if (length(y_vars_rm) == 0) {
        y_vars_rm <- as.data.table(t(rep(NA, times = expected_col_num)))
      }
      
    } else {
      # Check if the number of rows is 0
      #print(dim(x_vars_rm))
      #print(dim(y_vars_rm))
      if(x_vars_rm[, .N] == 0) {
        x_vars_rm <- rbind(x_vars_rm, as.data.table(t(rep(NA, times = ncol(x_vars_rm)))), 
                       use.names = F)
      }
      if(y_vars_rm[, .N] == 0) {
        y_vars_rm <- rbind(y_vars_rm, as.data.table(t(rep(NA, times = ncol(y_vars_rm)))), 
                       use.names = F)
      }
    }
  }

  comb_dt <- tryCatch({
    rbindlist(list(x_vars_rm, y_vars_rm), use.names = T, fill = T)
  },error=function(cond){
    print(cond)
    return(x_vars_rm)
  })
  return(comb_dt)
}


############################################################

#' A class which holds multiple results for each foreach loop iteration.
#' Each loop iteration populates two properties: $summary_table and $variables_rm.
#' @param summary_table the output summary table results from linear regression
#' for the given target gene
#' @param variables_rm a table with the variables removed for that target gene
#' during multicollinearity correction
multiResultClass <- function(summary_table = NULL, variables_rm = NULL) {
  me <- list(summary_table = summary_table, variables_rm = variables_rm)

  ## Set the name for the class
  class(me) <- append(class(me), "multiResultClass")
  return(me)
}


############################################################

#' A function to get a list of randomized input DFs (either per-sample or 
#' per-target expression DFs, or shuffled starter DFs) that will be input to our 
#' regularized regression model to generate an empirical p-value
#' @param signif_eval_type either "randomization_predictors", 
#' "randomization_perTarg", or "randomization_perSamp" to denote what type of 
#' significance evaluation we are doing, and therefore what data frame we will 
#' need to randomize and along what axis
#' @param lm_input_df for the predictor randomization; a data frame with 
#' covariates for the drivers 
#' @param expression_df for the per-sample and per-target expression 
#' randomization; an expression data frame (first column is ENSG ID)
#' @param num_randomizations the number of randomizations we are performing 
#' (will equal the length of the output list)
get_shuffled_input_dfs <- function(signif_eval_type, lm_input_df, expression_df, 
                                   num_randomizations) {
  
  shuffled_input_dfs <- NA
  
  # If we are doing a regularized regression model with a randomization approach 
  # to generating a p-value, we want to pre-shuffle the input data frame 
  # per-driver and save these DFs as a list
  if (signif_eval_type == "randomization_predictors") {
    swissprot_ids <- unique(unlist(lapply(colnames(lm_input_df)[
      grepl("MutStat_i", colnames(lm_input_df))], function(x)
        unlist(strsplit(x, "-", fixed = T))[1])))
    
    driver_indices <- lapply(swissprot_ids, function(s) 
      which(grepl(s, colnames(starter_df))))
    remaining_indices <- setdiff(1:ncol(starter_df), unlist(driver_indices))
    
    non_driver_starter_df <- starter_df[, remaining_indices, with = F]
    
    shuffled_input_dfs <- lapply(1:num_randomizations, function(i) {
      indiv_driver_shuffled_df <- lapply(1:length(driver_indices), function(j) {
        driver_i <- driver_indices[[j]]
        shuffled_cols <- starter_df[dqsample(1:nrow(starter_df)), driver_i, 
                                    with = F]
        return(shuffled_cols)
      })
      shuffled_df <- do.call("cbind", c(indiv_driver_shuffled_df, 
                                        non_driver_starter_df))
      return(shuffled_df)
    })
    
  } else if (signif_eval_type == "randomization_perSamp") {
    expression_df_sub <- expression_df[,colnames(expression_df) %fin% 
                                         lm_input_df$sample_id, with = F]
    
    shuffled_input_dfs <- lapply(1:num_randomizations, function(i) {
      expression_df_rand <- apply(expression_df_sub, MARGIN = 2,  function(x) 
        dqrng::dqsample(x))
      expression_df_rand[is.na(expression_df_rand)] <- 0
      
      expression_df_rand <- cbind(expression_df$ensg_id, expression_df_rand)
      colnames(expression_df_rand)[1] <- "ensg_id"
      
      return(as.data.table(expression_df_rand))
    })
    
  } else if (signif_eval_type == "randomization_perTarg") {
    # NOTE: would be slower and more memory inefficient to do it this way.
    # Not currently implemented.
    print("randomization_perTarg is not currently implemented. Please try again
    with a different evaluation type.")
    
  } else {
    print(paste("error: invalid significance evaluation type", signif_eval_type))
  }
  
  return(shuffled_input_dfs)
}


############################################################

#' Make sure any -Inf or Inf values in the LM input table are given NA and then 
#' omitted. If a column is >25% Inf or -Inf, remove that column first
#' @param lm_input_table the LM input table for a given target gene
fix_lm_input_table <- function(lm_input_table) {
  
  nrow <- lm_input_table[, .N]
  
  # Check for columns that are mostly -Inf or Inf
  cols_to_keep <- unlist(lapply(1:ncol(lm_input_table), function(i) {
    col <- unlist(lm_input_table[,i, with = F])
    if(length(col[(col == Inf) | (col == -Inf) | 
                  (col == "Inf") | (col == "-Inf")]) > (0.25 * nrow)) {
      return(NA)
    } else {return(i)}
  }))
  cols_to_keep <- cols_to_keep[!is.na(cols_to_keep)]
  
  lm_input_table_sub <- lm_input_table[, cols_to_keep, with = F]
  
  # Then, remove all remaining rows with Inf or -Inf values
  lm_input_table_sub[(lm_input_table_sub == -Inf) | 
                       (lm_input_table_sub == "-Inf") |
                       (lm_input_table_sub == Inf) | 
                       (lm_input_table_sub == "Inf")] <- NA
  lm_input_table_sub <- na.omit(lm_input_table_sub)
  
  return(lm_input_table_sub)
}

############################################################

#' Import all the Dyscovr input tables from a list. 
#' Given a vector of linear model input file names (without path attached)
#' this function reads them into a list of data tables that are returned
#' @param lm_input_fns_sub vector of LM input files
#' @param input_lm_filepath the path to the linear model input files
#' @param driver_input_fn the filename for the regulatory protein(s) portion of 
#' the LM input DF
#' @param randomize T/F value to indicate whether or not we are randomizing 
#' the target variable
#' @param gene_names_sub vector of the genes associated with each of the LM 
#' input files
#' @param select_drivers a given semicolon-separated string with the drivers we 
#' wish to include in each DF
#' @param select_args a vector of covariates we plan to include in the 
#' consequent model
#' @param num_pcs the number of genotypic PCs to include
#' @param cna_df the copy number data frame, in case we want to adjust the CNA 
#' columns
#' @param adjust_cna_rel a T/F value indicating whether or not we are adjusting
#' the CNA columns to make them relative
import_lm_input_fns <- function(lm_input_fns_sub, input_lm_filepath,
                                driver_input_fn, randomize, gene_names_sub, 
                                select_drivers, select_args, num_pcs, 
                                cna_df, adjust_cna_rel) {
  
  # Get a threshold for number of samples
  # Set threshold to 25, as suggested in this paper: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0229345
  thres <- 25
  print(paste0("Minimum number of samples with data for a given gene:", thres))
  
  num_wo_yvar <- 0
  
  cna_samp_dict <- NA
  if(adjust_cna_rel) {cna_samp_dict <- get_cna_samp_dict(cna_df)}
  
  # Import the protein-of-interest input DF
  driver_df <- fread(paste(input_lm_filepath, driver_input_fn, sep = "/"), 
                      header = T)
  
  lm_input_dfs <- lapply(lm_input_fns_sub, function(fn) {
    
    # Import file
    if(!is.na(fn)) {
      file <- fread(paste(input_lm_filepath, fn, sep = "/"), header = T)
      
      # Filter and randomize, if needed
        if("ExpStat_k" %fin% colnames(file)) {
          if(nrow(file) >= thres) {
            if(randomize) {file$ExpStat_k <- dqsample(file$ExpStat_k)}
            
            # Add the driver/ clinical portions
            file <- merge(driver_df, file, by = 'sample_id', all.x = F, all.y = T)
            
            if(adjust_cna_rel) {file <- adjust_cna_relative(file, cna_samp_dict)}
            
            # Do some preliminary filtering of columns we do not need
            file_filt <- filter_input_file(file, select_drivers, select_args, 
                                           num_pcs)
            return(file_filt)
            
          } else {
            return(NA)
          }

      } else {
        num_wo_yvar <- num_wo_yvar + 1
        return(NA)
      }
    }
  })
  if(num_wo_yvar > 0) {
    print(paste("Number of genes without a y-variable column, also discarded: ",  
                num_wo_yvar))
  }
  
  names(lm_input_dfs) <- gene_names_sub
  return(lm_input_dfs)
}


############################################################

#' Helper function to call the import_lm_input_fns() function, given the size of 
#' the target set
#' @param lm_input_fns_sub the names of the input filenames to import (a vector)
#' @param input_lm_filepath the path where the input files are found
#' @param driver_input_fn the name of the driver protein(s) part of the LM 
#' input tables
#' @param randomize T/F value, whether we are randomizing expression
#' @param gene_names_sub the names of all the genes we are looking to import
#' @param select_drivers the names of particular drivers we'd like to look at
#' @param select_args particular covariates to include in the model
#' @param num_pcs the numerical number of genotypic PCs to include
#' @param cna_df the copy number data frame, in case we want to adjust the CNA 
#' columns
#' @param adjust_cna_rel a T/F value indicating whether or not we are adjusting
#' the CNA columns to make them relative
#' @param N threshold for splitting into two segments
#' @param debug a T/F value to indicate whether we're in debug mode
call_import_lm_file <- function(lm_input_fns_sub, input_lm_filepath, 
                                driver_input_fn, randomize, gene_names_sub, 
                                select_drivers, select_args, num_pcs, cna_df, 
                                adj_cna_rel, N, debug) {
  
  if(length(lm_input_fns_sub) > N) {
    lm_input_dfs <- import_lm_input_fns(lm_input_fns_sub[1:N], input_lm_filepath, 
                                        driver_input_fn, randomize, 
                                        gene_names_sub[1:N], select_drivers, 
                                        select_args, num_pcs, cna_df, 
                                        adj_cna_rel)
  } else {
    lm_input_dfs <- import_lm_input_fns(lm_input_fns_sub, input_lm_filepath, 
                                        driver_input_fn, randomize, 
                                        gene_names_sub, select_drivers, 
                                        select_args, num_pcs, cna_df, 
                                        adj_cna_rel)
  }
  
  if(debug) {
    print(paste0("Number of genes discarded:", 
                 length(lm_input_dfs[is.na(lm_input_dfs)])))
  }

  lm_input_dfs <- lm_input_dfs[!is.na(lm_input_dfs)]
  
  if(debug) {
    print("Length imported input DFs")
    print(length(lm_input_dfs))
  }
  
  return(lm_input_dfs)
}


############################################################

#' Adjust the CNA columns to be relative, rather than absolute
#' For all CNAStat columns (non-bucketed), adjust the value to 
#' be the relative CNA rather than the absolute CNA; in other words,
#' rather than log2(CN + 1), it would be log2(CN+1) / log2((CN+1) / (N+1)),
#' where N is the mean CN value for the middle 80% of the patient's CN values
#' @param file the file to adjust
#' @param cna_samp_dict a CN-sample dictionary that contains the N for each 
#' sample
adjust_cna_relative <- function(file, cna_samp_dict) {
  
  cna_cols <- which(grepl("CNA", colnames(file)))
  samples <- file$sample_id
  
  for (i in cna_cols) {
    cna_vals <- as.integer(unlist(file[,i, with = F]))
    new_cna_vals <- unlist(lapply(1:length(cna_vals), function(j) {
      cn <- cna_vals[j]
      # Undo the log2
      cn <- 2^cn
      # Divide by N+1 and take the log2 of that
      n <- cna_samp_dict[[samples[j]]]
      return(log2(cn/(n+1)))
    }))
    file[,i] <- new_cna_vals
  }
  return(file)
}


############################################################

#' Get a copy number sample dictionary (named list) that, for each sample, has 
#' the mean copy number for that sample, across the middle 80% of their genes
#' @param cna_df full copy number DF for all samples
get_cna_samp_dict <- function(cna_df) {
  
  samples <- colnames(cna_df)
  samples <- samples[!(samples %fin% c("Gene_Symbol", "ensg_id", "ensg_ids", 
                                       "Swissprot"))]
  ten_perc <- 0.1 * nrow(cna_df)
  
  dictionary <- lapply(samples, function(samp) {
    cna_vals_samp <- as.integer(unlist(cna_df[,colnames(cna_df) == samp, 
                                              with = F]))
    
    # Get rid of the top 10% and bottom 10% of CNA values for this sample
    cna_vals_samp <- sort(cna_vals_samp)
    cna_vals_samp_sub <- cna_vals_samp[ten_perc:(length(cna_vals_samp)-ten_perc)]
    
    # Find the mean of the remaining values
    mean_cn_samp <- mean(cna_vals_samp_sub)
    return(mean_cn_samp)
  })
  names(dictionary) <- samples
  
  return(dictionary)
}


############################################################

#' Filter out columns we don't need from the Dyscovr LM input table given the 
#' selected drivers, arguments, number of PCs, etc. 
#' @param lm_input_table the LM input table for a given target gene
#' @param select_drivers vector of driver Uniprot IDs to keep in model
#' @param covs_to_incl a vector, either c("ALL") to include all covariates, or a 
#' vector of all the covariates we want to keep
#' @param num_pcs a value from 0-3 indicating the number of  
#' principal components we are including as covariates in the model 
filter_input_file <- function(lm_input_table, select_drivers, covs_to_incl, 
                              num_pcs) {
  
  colnames_to_incl <- colnames(lm_input_table)
  
  # Excluding age_b1 and age_b2 since no patients fall into these categories
  colnames_to_incl <- colnames_to_incl[!((colnames_to_incl == "Age_b1") | 
                                           (colnames_to_incl == "Age_b2"))]
  
  # Remove the PCs we are not using
  pcs <- colnames_to_incl[grepl("PC", colnames_to_incl)]
  pcs_to_remove <- setdiff(pcs, pcs[0:num_pcs])
  colnames_to_incl <- colnames_to_incl[!(colnames_to_incl %fin% pcs_to_remove)]  
  
  # Remove the driver genes we are not using
  if(select_drivers != "ALL") {
    select_driver_vect <- unlist(strsplit(select_drivers, ";", fixed = T))
    drivers_in_tab <- unlist(lapply(
      colnames_to_incl[grepl("MutStat_i", colnames_to_incl)], function(x)
        unlist(strsplit(x, "_", fixed = T))[1]))
    drivers_to_exclude <- setdiff(drivers_in_tab, select_driver_vect)
    drivers_to_exclude <- drivers_to_exclude[drivers_to_exclude != ""]
    drivers_to_exclude_covs <- unlist(lapply(drivers_to_exclude, function(d) 
      return(colnames_to_incl[grepl(d, colnames_to_incl)])))
    colnames_to_incl <- colnames_to_incl[!(colnames_to_incl %fin% 
                                             drivers_to_exclude_covs)]
  }
  
  # Additionally, if we have any columns that are uniform (the same value in all 
  # patients) we'd like to remove those as well, as long as they are not 
  # MutStat_i or CNAStat_i 
  uniform_val_tf <- sapply(lm_input_table, function(x) length(unique(x)) == 1)
  if(T %fin% uniform_val_tf) {
    uniform_val_ind <- which(uniform_val_tf)
    uniform_vals <- colnames(lm_input_table)[uniform_val_ind]
    uniform_vals <- uniform_vals[!(grepl("MutStat_i", uniform_vals) | 
                                     grepl("CNAStat_i", uniform_vals))]
    colnames_to_incl <- colnames_to_incl[!(colnames_to_incl %fin% uniform_vals)]
  }
  
  # Now, of the variables we have left, include only those given in the 
  # covs_to_include parameter
  if(!(covs_to_incl[1] == "ALL")) {
    colnames_to_incl <- unlist(lapply(colnames_to_incl, function(x) {
      if(any(unlist(lapply(covs_to_incl, function(y) grepl(y, x))))) {return(x)}
    }))
  }
  lm_input_table_sub <- lm_input_table[, c(1, which(colnames(lm_input_table) %fin% 
                                                      colnames_to_incl)), with = F]
  
  return(lm_input_table_sub)
}


############################################################

#' Helper function to remove samples with outlier expression values in the given 
#' gene using standard IQR-based thresholds (Q1 and Q3)
#' @param lm_input_table the linear model input table
filter_expression_outlier_samples <- function(lm_input_table) {
  
  if(!is.na(lm_input_table)) {
    if(!lm_input_table[, .N] == 0) {
      # Get the expression values
      expr_vals <- unlist(as.numeric(lm_input_table$ExpStat_k))
      
      # Get the IQR
      iqr <- IQR(expr_vals)
      
      # Get Q1 and Q3 of data
      q1 <- as.numeric(quantile(expr_vals)[2])
      q3 <- as.numeric(quantile(expr_vals)[4])
      
      # Get Q1 - 1.5(IQR), lower limit, and Q3 + 1.5(IQR), upper limit
      lower_thres <- as.numeric(q1 - (1.5 * iqr))
      upper_thres <- as.numeric(q3 + (1.5 * iqr))
      
      lm_input_table_sub <- lm_input_table[
        (lm_input_table$ExpStat_k < upper_thres) & 
          (lm_input_table$ExpStat_k > lower_thres), ]
      #print(paste("Number of samples removed due to extreme expression outliers: ", 
      #(nrow(lm_input_table) - nrow(lm_input_table_sub))))
      
      return(lm_input_table_sub)
    }
  }
  return(NA)
}

#' Helper function to remove samples with outlier expression values in the given 
#' gene using a standard deviation threshold (> 3 SDs from the mean value)
#' @param lm_input_table the linear model input table
filter_expression_outlier_samples2 <- function(lm_input_table) {
  
  if(!is.na(lm_input_table)) {
    if(!lm_input_table[, .N] == 0) {
      # Get the expression values
      expr_vals <- unlist(as.numeric(lm_input_table$ExpStat_k))
      
      # Get the mean and standard deviation
      mean <- mean(expr_vals)
      std_dev <- sd(expr_vals)
      
      # Get the upper and lower bounds
      lower_thres <- as.numeric(mean - (std_dev * 3))
      upper_thres <- as.numeric(mean + (std_dev * 3))
      
      lm_input_table_sub <- lm_input_table[
        (lm_input_table$ExpStat_k < upper_thres) & 
          (lm_input_table$ExpStat_k > lower_thres), ]
      #print(paste("Number of samples removed due to extreme expression outliers: ", 
      #(nrow(lm_input_table) - nrow(lm_input_table_sub))))
      
      return(lm_input_table_sub)
    }
  }
  return(NA)
}


#' Helper function to remove samples with outlier expression values in the given 
#' gene using a hard cutoff of values 4 and -4, since with quantile-
#' normalization and per-sample normal transformation all samples have the same 
#' distribution and range of expression values. Values that were originally zero 
#' will be below -4 after the normalization procedure and should be removed.
#' @param lm_input_table the linear model input table
filter_expression_outlier_samples3 <- function(lm_input_table) {
  
  if(!is.na(lm_input_table)) {
    if(!lm_input_table[, .N] == 0) {
      
      lm_input_table_sub <- lm_input_table[abs(lm_input_table$ExpStat_k) > 4, ]
      #print(paste("Number of samples removed due to extreme expression outliers: ", 
      #(nrow(lm_input_table) - nrow(lm_input_table_sub))))
      
      return(lm_input_table_sub)
    }
  }
  return(NA)
}

#' Helper function to alter samples with outlier expression values in the given 
#' gene using a hard cutoff of values 4 and -4, since with quantile-normalization 
#' and per-sample normal transformation all samples have the same distribution 
#' and range of expression values. Values that were originally zero will be 
#' below -4 after the normalization procedure. Takes any values above or below 4 
#' and truncate them to be between -4 and 4 (if < -4, be -4, if >4, be 4).
#' @param lm_input_table the linear model input table
collapse_expression_outlier_samples <- function(lm_input_table) {
  
  if(!is.na(lm_input_table)) {
    if(!lm_input_table[, .N] == 0) {
      
      lm_input_table$ExpStat_k[which(lm_input_table$ExpStat_k < -4)] <- -4
      lm_input_table$ExpStat_k[which(lm_input_table$ExpStat_k > 4)] <- 4
      #print(paste("Number of samples removed due to extreme expression outliers: ", 
      #(nrow(lm_input_table) - nrow(lm_input_table_sub))))
      
      return(lm_input_table)
    }
  }
  return(NA)
}


############################################################

#' Construct LM input formula for Dyscovr given a linear model input table and
#' cancer type, to pass to the speedlm. Removes unnecessary covariates, etc.
#' @param lm_input_table the LM input table produced from the above function
#' @param protein_ids_df the input DF of all drivers/ regulatory proteins
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. 
#' Possible values are "bucketInclAmp", "bucketExclAmp", and "rawCNA"
#' @param meth_bucketing a T/F value indicating whether or not we are bucketing
#' methylation values
#' @param runQueriesJointly a T/ F value indicating whether or not we are 
#' running all our drivers/ regulatory proteins jointly in one model
construct_formula <- function(lm_input_table, protein_ids_df, cna_bucketing, 
                              meth_bucketing, runQueriesJointly) {
  
  print("Now constructing formula...")
  
  # Exclude sample ID and mutation status of driver genes (covariates of interest)
  colnames_to_incl <- colnames(lm_input_table)[!(colnames(lm_input_table) == 
                                                   "sample_id")]
  colnames_to_incl <- colnames_to_incl[!(colnames_to_incl == "_MutStat_i")]
  #if(length(colnames_to_incl) < 3) {print(lm_input_table)}
  
  # For the bucketed variables, remove the last bucket (we don't need it!)
  bucketed_vars <- c("Tot_Mut_b","Tumor_purity_b", "Tot_IC_Frac_b", 
                     "Cancer_type_b")
  
  if(meth_bucketing) {
    bucketed_vars <- c(bucketed_vars, "MethStat_k")
    if(!runQueriesJointly) {
      bucketed_vars <- c(bucketed_vars, "MethStat_i")
    } else {
      meth_stat_i <- unlist(lapply(protein_ids_df$swissprot, function(x) 
        paste0(x, "_MethStat_i")))
      bucketed_vars <- c(bucketed_vars, meth_stat_i)
    }
  }
  if(!((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing)))) {
    bucketed_vars <- c(bucketed_vars, "CNAStat_k")
    if(!runQueriesJointly) {
      bucketed_vars <- c(bucketed_vars, "CNAStat_i")
    } else {
      cna_stat_i <- unlist(lapply(protein_ids_df$swissprot, function(x) 
        paste0(x, "_CNAStat_i")))
      bucketed_vars <- c(bucketed_vars, cna_stat_i)
    }
  }
  
  vals_to_remove <- unlist(lapply(bucketed_vars, function(var) {
    # Get the last bucket for this variable
    matching_vars <- colnames_to_incl[grepl(var, colnames_to_incl)]
    index <- which(colnames_to_incl == matching_vars[length(matching_vars)])
    return(index)
  }))
  colnames_to_incl <- colnames_to_incl[-vals_to_remove]
  
  # Remove these variables and construct formula from the remaining variables
  colnames_to_incl <- colnames_to_incl[!(colnames_to_incl == "ExpStat_k")]
  formula <- paste(colnames_to_incl, collapse = " + ") 
  formula <- paste("ExpStat_k ~", formula)
  
  return(formula)
}


############################################################

#' Given a list of  patient IDs (XXXX, e.g. A0WY), subsets a given data frame 
#' with column/row names that contain sample IDs (XXXX-XXX, e.g. A0WY-01A) to 
#' include only those given patients
#' @param patients a vector of patient IDs
#' @param df the data frame to subset using the patient ids
#' @param colNames a T/F value indicating whether the patient/sample IDs are
#' found on the rows or the columns of the DF (T is columns, F is rows)
#' @param dataset either 'TCGA' or 'METABRIC'
subset_by_intersecting_ids <- function(patients, df, colNames, dataset) {
  
  df <- as.data.frame(df)
  df_adj <- NA
  
  column_names_to_keep <- c("ensg_id", "Swissprot", "swissprot", 
                            "gene_name", "ensg_ids")
  
  if(colNames == T) {
    # Keep only intersecting patients
    label_col_indices <- which(colnames(df) %fin% column_names_to_keep)
    just_patients <- colnames(df)
    if(dataset == "TCGA") {
      just_patients <- unlist(lapply(colnames(df), function(x) 
        unlist(strsplit(x, "-", fixed = T))[1]))
    }
    df_adj <- df[, c(label_col_indices, which(just_patients %fin% patients))]

    # Keep only cancer samples
    if (dataset == "TCGA") {
      df_adj_label_col_indices <- which(colnames(df_adj) %fin% 
                                          column_names_to_keep)
      df_adj <- df_adj[, c(df_adj_label_col_indices, 
                           which(grepl("-0", colnames(df_adj), fixed = T)))]
      
      # Remove any duplicate samples
      df_adj <- df_adj[, !(grepl("-1", colnames(df_adj), fixed = T))]
    }
    
  } else {
    if(length(df) > 0) {
      if(nrow(df) > 0) {
        if(length(unlist(strsplit(as.character(df$sample_id[1]), "-", 
                                  fixed = T))) == 4) {
          df$sample_id <- unlist(lapply(df$sample_id, function(x) 
            paste(unlist(strsplit(x, "-", fixed = T))[3:4], collapse = "-")))
        }
        
        # Keep only intersecting patients
        just_patients <- df$sample_id
        if(dataset == "TCGA") {
          just_patients <- unlist(lapply(df$sample_id, function(x) 
            unlist(strsplit(x, "-", fixed = T))[1]))
        }
        df_adj <- df[which(just_patients %fin% patients),]
      }
    }
  }
  return(as.data.table(df_adj))
}


############################################################

#' Removes metastatic samples from a given data frame (sample ID will be
#' -06 or -07 to indicate metastasis) - only for TCGA
#' @param df the data frame to remove metastatic samples from
#' @param colNames a T/F value indicating whether the patient/sample IDs are
#' found on the rows or the columns of the DF (T is columns, F is rows)
remove_metastatic_samples <- function(df, colNames) {
  df <- as.data.frame(df)
  
  column_names_to_keep <- c("ensg_id", "Swissprot", "swissprot", "gene_name", 
                            "ensg_ids", "Patient", "patient")
  
  if(colNames == T) {
    # Keep only primary patients
    label_col_indices <- which(colnames(df) %fin% column_names_to_keep)
    just_samples <- unlist(lapply(colnames(df), function(x) 
      unlist(strsplit(x, "-", fixed = T))[2]))
    df_adj <- df[, unique(c(label_col_indices, 
                            intersect(which(just_samples != "06"), 
                                      which(just_samples != "07"))))]

  } else {
    if(length(unlist(strsplit(as.character(df$sample_id[1]), "-", 
                              fixed = T))) == 4) {
      df$sample_id <- unlist(lapply(df$sample_id, function(x) 
        paste(unlist(strsplit(x, "-", fixed = T))[3:4], collapse = "-")))
    }
    
    # Keep only primary patients
    just_samples <- unlist(lapply(df$sample_id, function(x) 
      unlist(strsplit(x, "-", fixed = T))[2]))
    df_adj <- df[intersect(which(just_samples != "06"),
                           which(just_samples != "07")),]
  }
  return(as.data.table(df_adj))
}


############################################################

#' Subset the given targets to only genes that are trans to the current driver
#' protein (found more than X (10E6 default) BP away on the same chromosome or 
#' on another chromosome entirely.
#' @param driver the regulatory protein in question, ENSG ID
#' @param target_df the full data table of targets we want to subset
#' @param all_genes_id_conv a conversion file from BioMart with information about
#' chromosome number and position
#' @param X the number of bases within which we will consider a gene 'cis'
#' @param debug a T/F value indicating whether or not we are in debug mode
subset_to_trans_targets <- function(driver, target_df, all_genes_id_conv, 
                                    debug, X = 10E6) {
  
  # Get the chromosome, starting position, and ending position (in BP) for the 
  # given driver
  chr_driver <- NA
  start_position_driver <- NA
  end_position_driver <- NA
  
  for (r in driver) {
    chr_driver <- c(chr_driver, as.numeric(unlist(all_genes_id_conv[
      all_genes_id_conv$ensembl_gene_id == r, 'chromosome_name'])))
    start_position_driver <- c(start_position_driver, as.numeric(
      unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == r, 
                               'start_position'])))
    end_position_driver <- c(end_position_driver, as.numeric(
      unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == r, 
                               'end_position'])))
  }
  
  chr_driver <- unique(chr_driver[!is.na(chr_driver)])
  start_position_driver <- unique(
    start_position_driver[!is.na(start_position_driver)])
  end_position_driver <- unique(
    end_position_driver[!is.na(end_position_driver)])
  
  if(debug) {
    print(paste("Chr_driver", chr_driver))
    print(paste("Start_pos_driver", start_position_driver))
    print(paste("End_pos_driver", end_position_driver))
  }
  
  # For each target, check if it is trans or not
  trans_target_df_rows <- lapply(1:target_df[, .N], function(i) {
    
    # Get the target
    t <- unlist(strsplit(target_df$ensg[i], ";", fixed = T))
    
    # Get it's chromosome
    if (length(t) > 1) {t <- t[1]}
    chr_targ <- unique(as.numeric(unlist(all_genes_id_conv[
      all_genes_id_conv$ensembl_gene_id == t, 'chromosome_name'])))
    
    # If not the same chromosome, return it
    if((length(chr_driver) == 0) | (length(chr_targ) == 0)) {
      return(target_df[i,])
    }
    # NAs are induced by coercion for X & Y chromosomes
    if(is.na(chr_targ)) {return(target_df[i,])}  
    if(chr_driver == "") {return(target_df[i,])}
    if (!(chr_targ == chr_driver)) {return(target_df[i,])}
    
    # If on the same chromosome, check the distance apart
    else {
      start_position_targ <- unique(as.numeric(unlist(all_genes_id_conv[
        all_genes_id_conv$ensembl_gene_id == t, 'start_position'])))
      end_position_targ <- unique(as.numeric(unlist(all_genes_id_conv[
        all_genes_id_conv$ensembl_gene_id == t, 'end_position'])))
      
      if(debug) {
        print(paste("Chr_targ", chr_targ))
        print(paste("Start_pos_targ", start_position_targ))
        print(paste("End_pos_targ", end_position_targ))
      }
      
      distance_min <- min(abs(start_position_targ - end_position_driver),
                          abs(start_position_driver - end_position_targ),
                          abs(start_position_targ - start_position_driver),
                          abs(end_position_driver - end_position_targ))
      if (distance_min > X) {return(target_df[i,])}
    }
  })
  
  trans_target_df <- do.call(rbind, trans_target_df_rows)
  
  return(trans_target_df)
}


