############################################################
### Linear Model Helper Functions
### Written By: Sara Geraghty, July 2020
############################################################

# This file contains a set of helper functions that will be called by the linear
# model (from linear_model.R or linear_model.R) while it is running.

library.path <- .libPaths()

library(stringr, lib.loc = library.path, quietly = TRUE)
library(data.table, lib.loc = library.path, quietly = TRUE)
library(broom, lib.loc = library.path, quietly = TRUE)

############################################################

#' Given a patient/ sample DF (as a data.table object), as well
#' as a vector of PCs to include, restricts the patient/sample DF
#' to only these PCs
#' @param patient_df a data.table object of imported patient/sample DF
#' @param mut_pc_vect a vector of integers denoting mutation matrix PCs
#' to include, given as a character (e.g. "1,2,4")
restrict_mut_pc_cols <- function(patient_df, mut_pc_vect) {
  # Convert the PC vector to be numeric
  mut_pc_vect <- as.numeric(unlist(strsplit(mut_pc_vect, ",", fixed = TRUE)))
  
  # Identify the PC columns within the patient DF
  patient_df_pc_cols <- which(grepl("Mut_PC", colnames(patient_df)))
  
  if(length(patient_df_pc_cols) > 0) {
    # Ensure the largest given PC is within the DF (discard any that are larger than what we have)
    if(TRUE %in% unlist(lapply(mut_pc_vect, function(x) ifelse(x > length(patient_df_pc_cols), TRUE, FALSE)))) {
      print("Some mutation matrix PCs given exceed the number of PCs provided in patient DF. Restricting PCs
          to those within given range.")
      mut_pc_vect <- mut_pc_vect[mut_pc_vect <= length(patient_df_pc_cols)]
    }
    
    # Remove any PC columns we don't want to keep
    cols_to_remove <- unlist(lapply(1:length(patient_df_pc_cols), function(i) {
      if(!(i %in% mut_pc_vect)) {return(patient_df_pc_cols[i])}
    }))
    patient_df <- patient_df[, (cols_to_remove) := NULL]
  }
  
  return(patient_df)
}

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
#' Possible values are "bucket_inclAmp", "bucket_exclAmp", "bucket_justAmp",
#' "bucket_justDel", and "rawCNA"
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3'
get_cna_stat <- function(cna_df, sample, cna_bucketing, dataset) {

  # Log2(CNA Stat + 1)
  if(cna_bucketing == "rawCNA") {
    cna_stat <- unique(as.integer(unlist(cna_df[,colnames(cna_df) == sample, with = FALSE])))
    if (length(cna_stat) > 1) {cna_stat <- cna_stat[1]}
    cna_stat <- log2(cna_stat + 1)
    if(length(cna_stat) < 1) {cna_stat <- NA}
    if (is.nan(cna_stat)) {cna_stat <- NA}
  
  # 1 if amplified ("justAmp") or deleted ("justDel"), 0 if 2+ copies
  } else if ((cna_bucketing == "bucket_justAmp") | (cna_bucketing == "bucket_justDel")) {
    cna_stat <- unique(as.integer(unlist(cna_df[,colnames(cna_df) == sample, with = FALSE])))
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
        } else {print("Currently only implemented for TCGA and METABRIC."); return(NA)}
      } else {
        if(dataset == "TCGA") {
          # Include partial deletions
          if(cna_stat < 2) {cna_stat <- 1}
          # Include only full deletions
          #if(cna_stat == 0) {cna_stat <- 1}
          else {cna_stat <- 0}
        } else if (dataset == "METABRIC") {
          # Include partial deletions
          if(cna_stat < 0) {cna_stat <- 1}
          # Include only full deletions
          #if(cna_stat == -2) {cna_stat <- 1}
          else {cna_stat <- 0}
        } else {print("Currently only implemented for TCGA and METABRIC."); return(NA)}
      }
    }
  
  # 3 buckets for deleted, normal, and amplified
  } else {
    if(dataset == "TCGA") {
      str_buckets <- as.character(unlist(cna_df[,colnames(cna_df) == sample, with = FALSE]))
      if(!(length(str_buckets) == 0)) {
        buckets_list <- unlist(strsplit(unlist(strsplit(str_buckets, "(", fixed = TRUE))[2], ")", fixed = TRUE))[1]
        cna_stat <- list(as.integer(unlist(strsplit(buckets_list, ",", fixed = TRUE))))
      } else {cna_stat <- list(NA, NA, NA)}
    } else if (dataset == "METABRIC") {
      cna_stat <- unique(as.integer(unlist(cna_df[,colnames(cna_df) == sample, with = FALSE])))
      if (length(cna_stat) > 1) {cna_stat <- cna_stat[1]}
      if(length(cna_stat) < 1) {cna_stat <- list(NA, NA, NA)}
      else if (is.nan(cna_stat) | is.na(cna_stat)) {cna_stat <- list(NA, NA, NA)}
      else {
        if(cna_stat < 0) {cna_stat <- list(1,0,0)}
        else if (cna_stat > 0) {cna_stat <- list(0,0,1)}
        else if (cna_stat == 0) {cna_stat <- list (0,1,0)}
        else {print("Invalid CNA value. Returning NA."); return(NA)}
      }
    } else {print("Currently only implemented for TCGA and METABRIC."); return(NA)}
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
      meth_buckets <- methylation_df[,grepl(sample, colnames(methylation_df)), with = FALSE]
      str_buckets <- as.character(unlist(meth_buckets))
      if(length(str_buckets) > 1) {
        if(length(unique(str_buckets)) == 1) {str_buckets <- unique(str_buckets)}
        else {str_buckets <- str_buckets[1]}
      }
      print(str_buckets)
      if(!(length(str_buckets) == 0)) {
        buckets_list <- unlist(strsplit(unlist(strsplit(str_buckets, "(", fixed = TRUE))[2], ")", fixed = TRUE))[1]
        meth_stat <- list(as.integer(unlist(strsplit(buckets_list, ",", fixed = TRUE))))
      } else {meth_stat <- list(NA, NA, NA)}
      #print(meth_stat)
    } else {  
      meth_stat <- unlist(methylation_df[,grepl(sample, colnames(methylation_df)), with = FALSE])
      if (length(meth_stat) > 1) {meth_stat <- mean.default(as.numeric(meth_stat))}
      if(length(meth_stat) < 1) {meth_stat <- NA}
      if (is.nan(meth_stat)) {meth_stat <- NA}
    }
  } else {
    meth_stat <- get_tnm_methylation(methylation_df, sample, meth_bucketing)
    print(meth_stat)
    if(length(meth_stat) < 1) {meth_stat <- NA}
    if (length(meth_stat) > 1) {meth_stat <- meth_stat[1]}
    if (is.nan(meth_stat)) {meth_stat <- NA}
  }
  return(as.numeric(meth_stat))
}

############################################################

#' Given an expression and variable indicating if we are looking tumor-normal-
#' matched, extracts and processes the expression or fold change expression 
#' @param expression_df an expression data frame of interest
#' @param sample the sample ID of the current sample
#' @param tumNormMatched a TRUE/FALSE value indicating whether or not the analysis is 
#' tumor-normal matched
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3'
#' @param log_expression a TRUE/FALSE value indicating whether or not we need to take 
#' the log2(exp + 1), or whether the expression values have already been transformed
#' to a normal distribution
get_exp_stat <- function(expression_df, sample, tumNormMatched, dataset, log_expression) {
  if(!tumNormMatched) {
    exp_stat <- as.numeric(unlist(expression_df[, grepl(sample, colnames(expression_df)), 
                                                with = FALSE]))
    if((dataset == "TCGA") & log_expression) {
      # Add a pseudocount so we don't take the log of 0 (1 is typically used for logged values)
      exp_stat <- log2(exp_stat + 1)
    }
  
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
#' Create the appropriate output file path (add the cancer type, the regulatory protein group 
#' or tester name, tumor-normal matched or tumor only, and the analysis type (eQTL or meQTL)
#' @param cancerType the name of the type of cancer, or "PanCancer" for all cancers
#' @param specificType if we are looking pan-cancer, but for specific cancer types 
#' individually, the specific cancer type we are looking at is given
#' @param test a TRUE/ FALSE value indicating whether this is a test on a single
#' regulatory protein
#' @param tester_name if this is a test, the name of the tester regulatory protein
#' @param run_name string value indicating the type of run if we are not running 
#' a test on a single regulatory protein 
#' @param tumNormMatched a TRUE/ FALSE value indicating whether or not this is a 
#' tumor-normal matched run
#' @param QTL_type either "eQTL" or "meQTL" depending on the type of run we are running
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3'
create_lm_input_file_outpath <- function(cancerType, specificType, test, tester_name, 
                                run_name, tumNormMatched, QTLtype, dataset) {
  outpath <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files"
  
  if(dataset != "TCGA") {
    outpath <- paste(outpath, dataset, sep = "/")
    
  } else {
    outpath <- paste(outpath, cancerType, sep = "/")
    if(!(specificType == "")) {outpath <- paste(outpath, specificType, sep = "/")}
    
    if(tumNormMatched) {
      outpath <- paste(outpath, "tumor_normal_matched", sep = "/")
    } else {outpath <- paste(outpath, "tumor_only", sep = "/")}
  }

  outpath <- paste(outpath, "lm_input_tables", sep = "/")
  
  if(test) {
    outpath <- paste(outpath, tester_name, sep = "/")
  } else {outpath <- paste(outpath, run_name, sep = "/")}
  
  outpath <- paste(outpath, QTLtype, sep = "/")
  print(paste("File outpath:", outpath))
  
  return(outpath)
}


############################################################
#' Create the appropriate output file name for the master DF result based on the
#' conditions of the given run
#' @param test a TRUE/ FALSE value indicating whether this is a test on a single
#' regulatory protein
#' @param tester_name if this is a test, the name of the tester regulatory protein
#' @param run_name string value indicating the type of run if we are not running 
#' a test on a single regulatory protein
#' @param target_uniprot a Uniprot ID of the target gene for the given DF
#' @param expression_df_name the string name used to import the expression DF
#' @param cna_bucketing a string description of how we are bucketing/ handling 
#' CNA values
#' @param meth_bucketing a TRUE/ FALSE value indicating if we are in some way 
#' bucketing our methylation value
#' @param meth_type 'Beta' or 'M' depending on the type of methylation we are
#' using
#' @param patient_df_name the string name used to import the patient DF
#' @param num_PEER an integer value from 0 to 10 indicating the number of PEER 
#' factors we are using
#' @param num_pcs an integer value from 0 to 2 indicating the number of principal
#' components we are using
#' @param randomize a TRUE/ FALSE value indicating whether or not we are randomizing
#' our dependent variable (expression or methylation)
#' @param patients_to_incl_label a label indicating the subpopulation we are restricting 
#' to, if any
#' @param removeCis a TRUE/ FALSE value indicating whether or not we have removed
#' cis gene pairings 
#' @param removeMetastatic a TRUE/ FALSE value indicating whether or not we have removed 
#' metastatic samples from the analysis
create_lm_input_table_filename <- function(test, tester_name, run_name, target_uniprot, expression_df_name,
                                   cna_bucketing, meth_bucketing, meth_type, patient_df_name, num_PEER, 
                                   num_pcs, randomize, patients_to_incl_label, removeCis, removeMetastatic) {
  outfn <- "lm_input"
  
  # Add the name of the patient population of interest, if there is one
  if(!(patients_to_incl_label == "")) {
    outfn <- paste(outfn, patients_to_incl_label, sep = "_")
  }
  
  # Add the name of the regulatory protein/ group of regulatory proteins
  if(test) {
    outfn <- paste(outfn, tester_name, sep = "_")
  } else {outfn <- paste(outfn, run_name, sep = "_")} 
  
  # Add the name of the target 
  outfn <- paste(outfn, target_uniprot[1], sep = "_")
  
  # Add the decisions for expression normalization, CNA bucketing, mutation restriction, 
  # methylation bucketing, and immune cell deconvolution bucketing 
  expr_norm <- unlist(strsplit(expression_df_name, "_", fixed = TRUE))[2]
  cna_bucketing_lab <- "rawCNA"
  if (cna_bucketing != "rawCNA") {cna_bucketing_lab <- paste("CNA", cna_bucketing, sep = "")}
  #mutation_restr <- unlist(strsplit(mutation_regprot_df_name, "_", fixed = TRUE))[1]
  meth_b_lab <- "Raw"
  if(meth_bucketing) {meth_b_lab <- "Bucketed"}
  meth_label <- paste("meth", paste(meth_type, meth_b_lab, sep = ""), sep = "")
  ici_label_spl <- unlist(strsplit(patient_df_name, "_", fixed = TRUE))
  ici_label <- ici_label_spl[4]
  if(ici_label_spl[5] == "total") {ici_label <- paste(ici_label, "TotalFrac", sep = "")}
  if(ici_label_spl[5] == "abs") {ici_label <- paste(ici_label, "Abs", sep = "")}
  
  outfn <- paste(outfn, paste(expr_norm, paste(cna_bucketing_lab, paste(meth_label, ici_label, sep = "_"), 
                                                                     sep = "_"), sep = "_"), sep = "_")
  # paste(mutation_restr, , sep = "_")
  
  # Add PEER factors/ PCs if needed
  if(!(num_PEER == 0)) {outfn <- paste(outfn, paste(num_PEER, "PEER", sep = ""), sep = "_")}
  if(!(num_pcs == 0))  {outfn <- paste(outfn, paste(num_pcs, "PCs", sep = ""), sep = "_")}
  
  # If we are NOT removing cis comparisons, add a label to denote this
  if(removeCis) {outfn <- paste(outfn, "RmCis", sep = "_")}
  
  # If we are NOT removing metastatic samples, add a label to denote this
  if(!removeMetastatic) {outfn <- paste(outfn, "notRmMetast", sep = "_")}
  
  # If we randomized expression, add "_RANDOMIZED" to the end of the file name
  if (randomize) {outfn <- paste(outfn, "_RANDOMIZED", sep = "")}
  
  print(paste("Outfile Name", outfn))
  
  return(outfn)
}


############################################################
#' Create the appropriate output file path (add the cancer type, the regulatory protein group 
#' or tester name, tumor-normal matched or tumor only, and the analysis type (eQTL or meQTL)
#' @param cancerType the name of the type of cancer, or "PanCancer" for all cancers
#' @param specificType if we are looking pan-cancer, but for specific cancer types 
#' individually, the specific cancer type we are looking at is given
#' @param test a TRUE/ FALSE value indicating whether this is a test on a single
#' regulatory protein
#' @param tester_name if this is a test, the name of the tester regulatory protein
#' @param run_name string value indicating the type of run if we are not running 
#' a test on a single regulatory protein 
#' @param tumNormMatched a TRUE/ FALSE value indicating whether or not this is a 
#' tumor-normal matched run
#' @param QTL_type either "eQTL" or "meQTL" depending on the type of run we are running
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3'
create_file_outpath <- function(cancerType, specificType, test, tester_name, 
                                run_name, tumNormMatched, QTLtype, dataset) {
  outpath <- "/Genomics/grid/users/scamilli/thesis_work/run-model-R/output_files"
  
  if(dataset != "TCGA") {outpath <- paste(outpath, dataset, sep = "/")}
  
  outpath <- paste(outpath, cancerType, sep = "/")
  
  if(test) {
    outpath <- paste(outpath, tester_name, sep = "/")
  } else {outpath <- paste(outpath, run_name, sep = "/")}
  
  if(dataset == "TCGA") {
    if(tumNormMatched) {
      outpath <- paste(outpath, "tumor_normal_matched", sep = "/")
    } else {outpath <- paste(outpath, "tumor_only", sep = "/")}
  }
  
  outpath <- paste(outpath, QTLtype, sep = "/")
  
  if(!(specificType == "")) {outpath <- paste(outpath, specificType, sep = "/")}
  
  print(paste("File outpath:", outpath))
  
  return(outpath)
}


############################################################
#' Create the appropriate output file name for the master DF result based on the
#' conditions of the given run
#' @param test a TRUE/ FALSE value indicating whether this is a test on a single
#' regulatory protein
#' @param tester_name if this is a test, the name of the tester regulatory protein
#' @param run_name string value indicating the type of run if we are not running 
#' a test on a single regulatory protein
#' @param targets_name a string name of the category of targets we are testing
#' @param expression_df_name the string name used to import the expression DF
#' @param cna_bucketing a string description of how we are bucketing/ handling 
#' CNA values
#' @param meth_bucketing a TRUE/ FALSE value indicating if we are in some way 
#' bucketing our methylation value
#' @param meth_type 'Beta' or 'M' depending on the type of methylation we are
#' using
#' @param patient_df_name the string name used to import the patient DF
#' @param num_PEER an integer value from 0 to 10 indicating the number of PEER 
#' factors we are using
#' @param num_pcs an integer value from 0 to 2 indicating the number of principal
#' components we are using
#' @param randomize a TRUE/ FALSE value indicating whether or not we are randomizing
#' our dependent variable (expression or methylation)
#' @param covs_to_incl_label a label indicating what covariates we've included in this run
#' @param patients_to_incl_label a label indicating the subpopulation we are restricting 
#' to, if any
#' @param model_type either "linear" or "mixed" (to indicate a linear mixed effects model)
#' @param removeCis a TRUE/ FALSE value indicating whether or not we have removed
#' cis gene pairings 
#' @param removeMetastatic a TRUE/ FALSE value indicating whether or not we have removed 
#' metastatic samples from the analysis
#' @param regularization either "None", "L1" or "Lasso", or "L2" or "Ridge" to indicate
#' if regularization was performed
#' @param signif_eval_type if regularization is not "None", provides additional info on how
#' we generated significance metrics
create_output_filename <- function(test, tester_name, run_name, targets_name, expression_df_name,
                                   cna_bucketing, meth_bucketing, meth_type, patient_df_name, 
                                   num_PEER, num_pcs, randomize, covs_to_incl_label, patients_to_incl_label, 
                                   model_type, removeCis, removeMetastatic, regularization, signif_eval_type) {
  outfn <- "res"
  
  # Add the name of the patient population of interest, if there is one
  if(!(patients_to_incl_label == "")) {
    outfn <- paste(outfn, patients_to_incl_label, sep = "_")
  }
  
  # Add the name of the regulatory protein/ group of regulatory proteins
  if(test) {
    outfn <- paste(outfn, tester_name, sep = "_")
  } else {outfn <- paste(outfn, run_name, sep = "_")} 
  
  # Add the name of the target group
  outfn <- paste(outfn, targets_name, sep = "_")
  
  # Add the decisions for expression normalization, CNA bucketing, mutation restriction, 
  # methylation bucketing, and immune cell deconvolution bucketing 
  expr_norm <- unlist(strsplit(expression_df_name, "_", fixed = TRUE))[2]
  cna_bucketing_lab <- "rawCNA"
  if (cna_bucketing != "rawCNA") {cna_bucketing_lab <- paste("CNA", cna_bucketing, sep = "")}
  #mutation_restr <- unlist(strsplit(mutation_regprot_df_name, "_", fixed = TRUE))[1]
  meth_b_lab <- "Raw"
  if(meth_bucketing) {
    meth_b_lab <- "Bucketed"
  }
  meth_label <- paste("meth", paste(meth_type, meth_b_lab, sep = ""), sep = "")
  ici_label_spl <- unlist(strsplit(patient_df_name, "_", fixed = TRUE))
  ici_label <- ici_label_spl[4]
  if(ici_label_spl[5] == "total") {ici_label <- paste(ici_label, "TotalFrac", sep = "")}
  if(ici_label_spl[5] == "abs") {ici_label <- paste(ici_label, "Abs", sep = "")}
  
  outfn <- paste(outfn, paste(expr_norm, paste(cna_bucketing_lab, paste(meth_label, ici_label, sep = "_"), 
                                                                     sep = "_"), sep = "_"), sep = "_")
  # paste(mutation_restr, , sep = "_")
  
  # Add PEER factors/ PCs if needed
  if(!(num_PEER == 0)) {outfn <- paste(outfn, paste(num_PEER, "PEER", sep = ""), sep = "_")}
  if(!(num_pcs == 0))  {outfn <- paste(outfn, paste(num_pcs, "PCs", sep = ""), sep = "_")}
  
  # If we are NOT removing cis comparisons, add a label to denote this
  if(removeCis) {outfn <- paste(outfn, "RmCis", sep = "_")}
  
  # If we are NOT removing metastatic samples, add a label to denote this
  if(!removeMetastatic) {outfn <- paste(outfn, "notRmMetast", sep = "_")}
  
  if(!model_type == "linear") {outfn <- paste(outfn, "mixed", sep = "_")}
  
  # If we have regularized in some way, denote this
  if(!(regularization == "None")) {
    if (regularization %in% c("L1", "lasso")) {
      if(grepl("randomization", signif_eval_type)) {
        if(signif_eval_type == "randomization_predictors") {
          outfn <- paste(outfn, "L1rand_pred", sep = "_")
        } else {outfn <- paste(outfn, "L1rand", sep = "_")}
      }
      else if(signif_eval_type == "subsampling") {outfn <- paste(outfn, "L1subsamp", sep = "_")}
      else if (signif_eval_type == "selectiveInference") {{outfn <- paste(outfn, "L1si", sep = "_")}}
      else {
        print("Error; unrecognized argument for method of evaluating LASSO significance.")
        outfn <- paste(outfn, "L1", sep = "_")
      }
    }
    if (regularization %in% c("L2", "ridge")) {outfn <- paste(outfn, "L2", sep = "_")}
    if (regularization == "bayesian.bgl") {outfn <- paste(outfn, "bayesian.bgl", sep = "_")}
    if (regularization == "bayesian.bglss") {outfn <- paste(outfn, "bayesian.bglss", sep = "_")}
  }
  
  # If we've restricted the number of covariates, add that label
  if(!(covs_to_incl_label == "")) {outfn <- paste(outfn, covs_to_incl_label, sep = "_")}
  
  # Add "uncorrected" to the file name so we know MHT correction has not yet been done on these results
  outfn <- paste(outfn, "_uncorrected", sep = "")
  
  # If we randomized expression, add "_RANDOMIZED" to the end of the file name
  if (randomize) {outfn <- paste(outfn, "_RANDOMIZED", sep = "")}
  
  print(paste("Outfile Name", outfn))
  
  return(outfn)
}

############################################################
#' Short function for combining the output tibbles produced for each row in list_of_rows
#' Used to recombine output of foreach
#' @param x_tab the first output table of results
#' @param y_tab the second output table of results
#' @param expected_col_num the expected number of columns, in case we need to turn
#' 2 NA inputs into a table accepted in the re-combining process
combine_tab <- function(x_tab, y_tab, expected_col_num = 7) {
  
  # Check if length of either table is 0; if so, handle like an NA table
  if((length(x_tab) == 0) | (length(y_tab) == 0)) {
    if (length(x_tab) == 0) {x_tab <- as.data.table(t(rep(NA, times = expected_col_num)))}
    if (length(y_tab) == 0) {y_tab <- as.data.table(t(rep(NA, times = expected_col_num)))}
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
    if(TRUE %in% c(is.na(x_tab), is.na(y_tab))) {
      x_tab <- na.omit(x_tab)
      y_tab <- na.omit(y_tab)
    }
    
    # Again, check if removing NAs/NaNs has set the length of either table is 0; 
    # if so, handle like an NA table
    if((length(x_tab) == 0) | (length(y_tab) == 0)) {
      if (length(x_tab) == 0) {x_tab <- as.data.table(t(rep(NA, times = expected_col_num)))}
      if (length(y_tab) == 0) {y_tab <- as.data.table(t(rep(NA, times = expected_col_num)))}
      
    } else {
      # Check if the number of rows is 0
      if(x_tab[, .N] == 0) {
        x_tab <- rbind(x_tab, as.data.table(t(rep(NA, times = ncol(x_tab)))), use.names = F)
      }
      if(y_tab[, .N] == 0) {
        y_tab <- rbind(y_tab, as.data.table(t(rep(NA, times = ncol(y_tab)))), use.names = F)
      }
    }
  }
  #print(dim(y_tab))

  comb_dt <- tryCatch({
    rbindlist(list(x_tab, y_tab), use.names = TRUE, fill = TRUE)
    #rbind(x_tab, y_tab)
    },error=function(cond){
      print(cond)
      return(x_tab)
    })

  return(comb_dt)
}


############################################################
#' A function to get a list of randomized input DFs (either per-sample or per-target expression
#' DFs, or shuffled starter DFs) that will be input to our regularized regression model 
#' to generate an empirical p-value
#' @param signif_eval_type either "randomization_predictors", "randomization_perTarg", or
#' "randomization_perSamp" to denote what type of significance evaluation we are doing,
#' and therefore what data frame we will need to randomize and along what axis
#' @param lm_input_df for the predictor randomization; a data frame with covariates for the drivers 
#' @param expression_df for the per-sample and per-target expression randomization; an expression
#' data frame (first column is ENSG ID)
#' @param num_randomizations the number of randomizations we are performing (will equal the 
#' length of the output list)
get_shuffled_input_dfs <- function(signif_eval_type, lm_input_df, expression_df, num_randomizations) {
  
  shuffled_input_dfs <- NA
  
  # If we are doing a regularized regression model with a randomization approach to 
  # generating a p-value, we want to pre-shuffle the input data frame per-driver and save 
  # these DFs as a list
  if (signif_eval_type == "randomization_predictors") {
    swissprot_ids <- unique(unlist(lapply(colnames(lm_input_df)[grepl("MutStat_i", colnames(lm_input_df))], function(x)
      unlist(strsplit(x, "-", fixed = TRUE))[1])))
    driver_indices <- lapply(swissprot_ids, function(s) which(grepl(s, colnames(starter_df))))
    remaining_indices <- setdiff(1:ncol(starter_df), unlist(driver_indices))
    non_driver_starter_df <- starter_df[, remaining_indices, with = FALSE]
    shuffled_input_dfs <- lapply(1:num_randomizations, function(i) {
      indiv_driver_shuffled_df <- lapply(1:length(driver_indices), function(j) {
        driver_i <- driver_indices[[j]]
        shuffled_cols <- starter_df[dqsample(1:nrow(starter_df)), driver_i, with = FALSE]
        return(shuffled_cols)
      })
      shuffled_df <- do.call("cbind", c(indiv_driver_shuffled_df, non_driver_starter_df))
      return(shuffled_df)
    })
  } else if (signif_eval_type == "randomization_perSamp") {
    expression_df_sub <- expression_df[,colnames(expression_df) %fin% lm_input_df$sample_id, with = FALSE]
    
    shuffled_input_dfs <- lapply(1:num_randomizations, function(i) {
      expression_df_rand <- apply(expression_df_sub, MARGIN = 2,  function(x) dqrng::dqsample(x))
      expression_df_rand[is.na(expression_df_rand)] <- 0
      
      expression_df_rand <- cbind(expression_df$ensg_id, expression_df_rand)
      colnames(expression_df_rand)[1] <- "ensg_id"

      return(as.data.table(expression_df_rand))
    })
    
  } else if (signif_eval_type == "randomization_perTarg") {
    #expression_df_sub <- cbind(expression_df$ensg_id, expression_df[,colnames(expression_df) %fin% 
                                                                      #starter_df$sample_id, with = FALSE])
    #colnames(expression_df_sub)[1] <- "ensg_id"
    
    #shuffled_input_dfs <- lapply(1:num_randomizations, function(i) {
    #  expression_df_rand <- apply(expression_df[,2:ncol(expression_df), with = FALSE], 
                                  #MARGIN = 1,  function(x) dqrng::dqsample(x))
    #  expression_df_rand[is.na(expression_df_rand)] <- 0
    #  return(expression_df_rand)
    #})

    # NOTE: would be slower and more memory inefficient to do it this way. Do the randomizations per target gene.
  } else {
    print(paste("error: invalid significance evaluation type", signif_eval_type))
  }
  
  return(shuffled_input_dfs)
}


############################################################
#### FILL IN TARGET T_K INPUTS TO TABLE
############################################################
#' Function takes a data frame of inputs for the given patient p_j and the regulatory 
#' protein r_i, as well as the data files for target t_k (and its ID), and uses 
#' them to construct a full tabular input to a linear model function for a given 
#' regulatory protein i-target gene k pairing of the following form:
#' Columns: Linear model input variables
#' Rows: Patients 
#' @param starter_df a 'starter DF' that has all the patient and reg. protein r_i characteristics 
#' @param targ_k the swissprot ID of target gene k
#' @param targ_k_ensg the ensembl ID of target gene k
#' @param mutation_targ_df the mutation gene target DF
#' @param methylation_df the methylation DF
#' @param cna_df the copy number alteration DF
#' @param expression_df the gene expression DF
#' @param log_expression a TRUE/ FALSE value indicating whether we take the log2(exp + 1),
#' or simply the expression value given in the DF
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. Possible
#' values are "bucketInclAmp", "bucketExclAmp", and "rawCNA"
#' @param meth_bucketing a TRUE/FALSE value indicating whether or not we are bucketing
#' methylation values
#' @param analysis_type a string label that reads either "eQTL" or "meQTL" to 
#' determine what kind of model we should run
#' @param tumNormMatched a TRUE/FALSE value indicating whether or not the analysis is 
#' tumor-normal matched
#' @param debug a TRUE/FALSE value indicating if we are in debug mode and should 
#' use additional prints
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3', to denote what columns
#' we are referencing/ specific data types we are using
fill_targ_inputs <- function(starter_df, targ_k, targ_k_ensg, mutation_targ_df,
                             methylation_df, cna_df, expression_df, log_expression, cna_bucketing,
                             meth_bucketing, analysis_type, tumNormMatched, debug, dataset) {
  
  #print(paste("Target:", targ_k))
  
  if(debug) {
    print(paste("Target:", targ_k))
    print(paste("Target ENSG:", targ_k_ensg))
  }
  
  # Begin by limiting the data frames to just this target gene k
  if (!(targ_k == "" | targ_k_ensg == "")) {
    mutation_targ_df <- filter_mut_by_uniprot(mutation_targ_df, targ_k) 
    cna_df <- filter_cna_by_ensg(cna_df, targ_k_ensg)
    methylation_df <- filter_meth_by_ensg(methylation_df, targ_k_ensg)
    exp_rows_to_keep <- unique(unlist(lapply(targ_k_ensg, function(x) grep(x, expression_df$ensg_id))))
    expression_df <- expression_df[exp_rows_to_keep,]
    
    if(debug) {
      print(paste("dim mutation targ df:", dim(mutation_targ_df)))
      print(paste("dim cna df:", dim(cna_df)))
      print(paste("dim methylation df:", dim(methylation_df)))
      print(paste("dim expression df:", dim(expression_df)))
    }
    
    # Continue only if all of these data frames contain information about this target gene k
    if(!(expression_df[, .N] == 0 | cna_df[, .N] == 0 | methylation_df[, .N] == 0 | 
         mutation_targ_df[, .N] == 0)) {
      
      # Loop through all samples to get the info for each target t_k column
      target_rows <- lapply(1:starter_df[, .N], function(i) {
        sample <- starter_df$sample_id[i]
        
        # Is this target mutated? 
        mut_stat <- get_mut_stat_targ(mutation_targ_df, sample)
        
        # Is this target amplified or deleted?
        cna_stat <- get_cna_stat(cna_df, sample, cna_bucketing, dataset)
        if(grepl("incl", cna_bucketing)) {cna_stat <- cna_stat[[1]]}    
        
        # Is this target methylated, or differentially methylated?
        if(analysis_type == "meQTL") {meth_bucketing <- FALSE}
        meth_stat <- get_meth_stat(methylation_df, sample, meth_bucketing, tumNormMatched)
        if(meth_bucketing) {
          meth_stat <- as.integer(meth_stat[[1]])
          if(length(meth_stat) > 3) {meth_stat <- meth_stat[1:3]}
        }
        # What is the expression of this target in this sample in cancer?
        exp_stat <- unlist(get_exp_stat(expression_df, sample, tumNormMatched, 
                                        dataset, log_expression))
        
        if(debug) {
          print(paste("exp_stat:", exp_stat))
          print(paste("mut_stat:", mut_stat))
          print(paste("cna_stat:", cna_stat))
          print(paste("meth_stat:", meth_stat))
        }
        
        # Make row for this patient and return
        if(!meth_bucketing) {
          if ((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing))) {
            return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, "CNAStat_k" = cna_stat, 
                              "MethStat_k" = meth_stat))
          } else {
            return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, "CNAStat_k_b1" = cna_stat[1], 
                              "CNAStat_k_b2" = cna_stat[2], "CNAStat_k_b3" = cna_stat[3], "MethStat_k" = meth_stat))
          }
        } else {
          if ((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing))) {
            return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, "CNAStat_k" = cna_stat, 
                              "MethStat_k_b1" = meth_stat[1], "MethStat_k_b2" = meth_stat[2], "MethStat_k_b3" = meth_stat[3]))
          } else {
            return(data.table("ExpStat_k" = exp_stat, "MutStat_k" = mut_stat, "CNAStat_k_b1" = cna_stat[1], 
                              "CNAStat_k_b2" = cna_stat[2], "CNAStat_k_b3" = cna_stat[3], "MethStat_k_b1" = meth_stat[1], 
                              "MethStat_k_b2" = meth_stat[2], "MethStat_k_b3" = meth_stat[3]))
          }
        }
      })
      
      tryCatch({
        # Bind these rows into the starter DF
        colnam <- colnames(target_rows[[1]])
        targ_k_df <- data.table::rbindlist(target_rows)
        colnames(targ_k_df) <- colnam
        full_df <- cbind(starter_df, targ_k_df)
        
      }, error=function(cond){
        print("Unexpected error; target rows is not valid: ")
        print(target_rows)
        full_df <- starter_df
      })
      
      return(full_df)
    } else {
      if(debug) {
        print(paste("Target not found in all DFs:", targ_k))
      }
      return(NA)
    }
  } else {
    if(debug) {
      print(paste("Target Uniprot or ENSG ID is missing,", targ_k))
    }
    return(NA)
  }
}


############################################################
#### LIMIT TARGETS FILE TO ONLY THOSE TARGETS WITHOUT A
#### FILE IN THE PROVIDED DIRECTORY
############################################################
#' Limits the target gene DF to only those targets that don't have an 
#' existing file in the outpath provided
#' @param outpath the path to which the LM input files will be written
#' @param targets_DF the data frame with gene targets for LM input file creation
#' @param gene_name_index the numerical index of where the target gene name is found
#' in a given file name
limit_to_targets_wo_existing_files <- function(outpath, targets_DF, gene_name_index) {
  # Get all the existing files in this path, if any exist
  files_in_outpath <- list.files(paste0(outpath, "/"), pattern = "lm_input_")
  
  # Check if there are sub-directories in this path that might contain the file
  dirs_in_path <- list.dirs(outpath, recursive = F)
  if(TRUE %in% unlist(lapply(dirs_in_path, function(d) grepl("part_", d)))) {
    dirs_part <- dirs_in_path[grepl("part_", dirs_in_path)]
    for (d in dirs_part) {
      files_d <- list.files(paste0(d, "/"), pattern = "lm_input_")
      files_in_outpath <- c(files_in_outpath, files_d)
    }
  } 
  genes_in_outpath <- unlist(lapply(files_in_outpath, function(f) 
    unlist(strsplit(f, "_", fixed = TRUE))[gene_name_index]))
  
  rows_to_keep <- unlist(lapply(1:targets_DF[, .N], function(i) {
    ids <- unlist(strsplit(as.character(targets_DF[i, 'swissprot']), ';', fixed = TRUE))
    if(length(intersect(ids, genes_in_outpath)) > 0) {return(NA)}
    else {return(i)}
  }))
  rows_to_keep <- rows_to_keep[!is.na(rows_to_keep)]

  targets_DF <- targets_DF[rows_to_keep,]
  print("Targets DF, subsetted")
  print(dim(targets_DF))
  
  return(targets_DF)
}

############################################################
#### ADJUST THE LINEAR MODEL INPUT TABLE TO EXCLUDE -INF OR INF VALS
############################################################
#' Make sure any -Inf or Inf values are given NA and then omitted
#' If a column is >25% Inf or -Inf, remove that column first
#' @param lm_input_table the LM input table for a given target gene
fix_lm_input_table <- function(lm_input_table) {
  
  nrow <- lm_input_table[, .N]
  
  # Check for columns that are mostly -Inf or Inf
  cols_to_keep <- unlist(lapply(1:ncol(lm_input_table), function(i) {
    col <- unlist(lm_input_table[,i, with = FALSE])
    if(length(col[(col == Inf) | (col == -Inf) | 
                  (col == "Inf") | (col == "-Inf")]) > (0.25 * nrow)) {
      return(NA)
    } else {return(i)}
  }))
  cols_to_keep <- cols_to_keep[!is.na(cols_to_keep)]
  
  lm_input_table_sub <- lm_input_table[, cols_to_keep, with = FALSE]
  
  # Then, remove all remaining rows with Inf or -Inf values
  lm_input_table_sub[(lm_input_table_sub == -Inf) | (lm_input_table_sub == "-Inf") |
                   (lm_input_table_sub == Inf) | (lm_input_table_sub == "Inf")] <- NA
  lm_input_table_sub <- na.omit(lm_input_table_sub)
  
  return(lm_input_table_sub)
}

############################################################
#### IMPORT ALL THE LINEAR MODEL INPUT TABLES FROM A LIST
############################################################
#' Given a vector of linear model input file names (without path attached)
#' this function reads them into a list of data tables that are returned
#' @param lm_input_fns_sub vector of LM input files
#' @param input_lm_filepath the path to the linear model input files
#' @param randomize TRUE/FALSE value to indicate whether or not we are randomizing 
#' the target variable
#' @param gene_names_sub vector of the genes associated with each of the LM input files
#' @param select_drivers a given semicolon-separated string with the drivers we 
#' wish to include in each DF
#' @param select_args a vector of covariates we plan to include in the consequent model
#' @param num_pcs the number of genotypic PCs to include
#' @param num_PEER the number of PEER factors to include
#' @param qtl_type the type of model (eQTL or meQTL)
import_lm_input_fns <- function(lm_input_fns_sub, input_lm_filepath, randomize, gene_names_sub,
                                select_drivers, select_args, num_pcs, num_PEER, qtl_type) {
  
  # Get a threshold for number of samples
  # Set threshold to 25, as suggested in this paper: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0229345
  thres <- 25
  print(paste0("Minimum number of samples with data for a given gene:", thres))
  
  num_wo_yvar <- 0
  
  lm_input_dfs <- lapply(lm_input_fns_sub, function(fn) {
    
    # Import file
    if(!is.na(fn)) {
      file <- fread(paste(input_lm_filepath, fn, sep = "/"), header = TRUE)
      
      # Filter and randomize, if needed
      if(qtl_type == "eQTL") {
        if("ExpStat_k" %fin% colnames(file)) {
          if(nrow(file) >= thres) {
            if(randomize) {file$ExpStat_k <- dqsample(file$ExpStat_k)}
            # Do some preliminary filtering of columns we do not need
            file_filt <- filter_input_file(file, select_drivers, select_args, 
                                           num_pcs, num_PEER)
            #print(head(file_filt))
            return(file_filt)
          } else {
            num_wo_yvar <- num_wo_yvar + 1
            return(NA)
          }
        } else {return(NA)}
        
      } else if(qtl_type == "meQTL") {
        if("MethStat_k" %fin% colnames(file)) {
          if(nrow(file) >= thres) {
            if(randomize) {file$MethStat_k <- dqsample(file$MethStat_k)}
            # Do some preliminary filtering of columns we do not need
            file_filt <- filter_input_file(file, select_drivers, select_args, 
                                           num_pcs, num_PEER)
            #print(head(file_filt))
            return(file_filt)
          } else {return(NA)}
        } else {
          num_wo_yvar <- num_wo_yvar + 1
          return(NA)
        }
      } else {
        print("QTL type is invalid. Only eQTL and meQTL are implemented. Returning NA.")
        return(NA)
      }
    }
  })
  
  print(paste("Number of genes without a y-variable column, also discarded: ",  num_wo_yvar))
  
  names(lm_input_dfs) <- gene_names_sub
  return(lm_input_dfs)
}


############################################################
#### FILTER THE LM INPUT TABLE TO EXCLUDE COLUMNS WE DON'T NEED
############################################################
#' Filter out columns we don't need from the LM input table given the selected
#' drivers, arguments, number of PCs/ PEER facors, and race. 
#' @param lm_input_table the LM input table for a given target gene
#' @param select_drivers vector of driver Uniprot IDs to keep in model
#' @param covs_to_incl a vector, either c("ALL") to include all covariates, or a vector
#' of all the covariates we want to keep
#' @param num_PEER a value from 0-10 indicating the number of PEER factors we 
#' are including as covariates in the model
#' @param num_pcs a value from 0-2 indicating the number of  
#' principal components we are including as covariates in the model 
filter_input_file <- function(lm_input_table, select_drivers, covs_to_incl, num_pcs, num_PEER) {

  colnames_to_incl <- colnames(lm_input_table)
  
  # Note: currently excluding age b1 and age b2 since no patients fall into these categories
  colnames_to_incl <- colnames_to_incl[!((colnames_to_incl == "Age_b1") | (colnames_to_incl == "Age_b2"))]
  
  # Remove the race variable entirely
  colnames_to_incl <- colnames_to_incl[!(grepl("Race", colnames_to_incl))]

  # Remove the PEER factors we are not using
  peer_factors <- colnames_to_incl[grepl("PEER", colnames_to_incl)]
  peer_fact_to_remove <- setdiff(peer_factors, peer_factors[0:num_PEER])
  colnames_to_incl <- colnames_to_incl[!(colnames_to_incl %fin% peer_fact_to_remove)]
  
  # Remove the PCs we are not using
  pcs <- colnames_to_incl[grepl("PC", colnames_to_incl)]
  pcs_to_remove <- setdiff(pcs, pcs[0:num_pcs])
  colnames_to_incl <- colnames_to_incl[!(colnames_to_incl %fin% pcs_to_remove)]  
  
  # Exclude the library size
  colnames_to_incl <- colnames_to_incl[!(colnames_to_incl == "Lib_Size")]
  
  if(select_drivers != "ALL") {
    select_driver_vect <- unlist(strsplit(select_drivers, ";", fixed = TRUE))
    drivers_in_tab <- unlist(lapply(colnames_to_incl[grepl("MutStat_i", colnames_to_incl)], function(x)
      unlist(strsplit(x, "_", fixed = TRUE))[1]))
    drivers_to_exclude <- setdiff(drivers_in_tab, select_driver_vect)
    drivers_to_exclude_covs <- unlist(lapply(drivers_to_exclude, function(d) 
      return(colnames_to_incl[grepl(d, colnames_to_incl)])))
    colnames_to_incl <- colnames_to_incl[!(colnames_to_incl %fin% drivers_to_exclude_covs)]
  }
  
  # Additionally, if we have any columns that are uniform (the same value in all patients) we'd 
  # like to remove those as well, as long as they are not MutStat_i, CNAStat_i or FNCStat_i
  uniform_val_tf <- sapply(lm_input_table, function(x) length(unique(x)) == 1)
  if(TRUE %fin% uniform_val_tf) {
    uniform_val_ind <- which(uniform_val_tf)
    uniform_vals <- colnames(lm_input_table)[uniform_val_ind]
    uniform_vals <- uniform_vals[!(grepl("MutStat_i", uniform_vals) | grepl("CNAStat_i", uniform_vals) | 
                                     grepl("FNCStat_i", uniform_vals))]
    colnames_to_incl <- colnames_to_incl[!(colnames_to_incl %fin% uniform_vals)]
  }
  
  # Now, of the variables we have left, include only those given in the covs_to_include parameter
  if(!(covs_to_incl[1] == "ALL")) {
    colnames_to_incl <- unlist(lapply(colnames_to_incl, function(x) {
      if(any(unlist(lapply(covs_to_incl, function(y) grepl(y, x))))) {return(x)}
    }))
  }
  lm_input_table_sub <- lm_input_table[, c(1, which(colnames(lm_input_table) %fin% colnames_to_incl)), with = FALSE]
  
  return(lm_input_table_sub)
}


#' Helper function to remove samples with outlier expression values in the given gene
#' using standard IQR-based thresholds (Q1 and Q3)
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
      
      lm_input_table_sub <- lm_input_table[(lm_input_table$ExpStat_k < upper_thres) & 
                                             (lm_input_table$ExpStat_k > lower_thres), ]
      #print(paste("Number of samples removed due to extreme expression outliers: ", 
                  #(nrow(lm_input_table) - nrow(lm_input_table_sub))))
      
      return(lm_input_table_sub)
    }
  }
  return(NA)
}

#' Helper function to remove samples with outlier expression values in the given gene
#' using a standard deviation threshold (> 3 SDs from the mean value)
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
      
      lm_input_table_sub <- lm_input_table[(lm_input_table$ExpStat_k < upper_thres) & 
                                             (lm_input_table$ExpStat_k > lower_thres), ]
      #print(paste("Number of samples removed due to extreme expression outliers: ", 
                  #(nrow(lm_input_table) - nrow(lm_input_table_sub))))
      
      return(lm_input_table_sub)
    }
  }
  return(NA)
}


#' Helper function to remove samples with outlier expression values in the given gene
#' using a hard cutoff of values 4 and -4, since with quantile-normalization and per-sample
#' normal transformation all samples have the same distribution and range of expression values.
#' Values that were originally zero will be below -4 after the normalization procedure and 
#' should be removed.
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

#' Helper function to alter samples with outlier expression values in the given gene
#' using a hard cutoff of values 4 and -4, since with quantile-normalization and per-sample
#' normal transformation all samples have the same distribution and range of expression values.
#' Values that were originally zero will be below -4 after the normalization procedure. Take any
#' values above or below 4 and truncate them to be between -4 and 4 (if < -4, be -4, if >4, be 4).
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
#### CONSTRUCT LM FORMULA
############################################################
#' Given a linear model input table, analysis and cancer type,
#' contructs a formula to pass to the speedlm formula. Removes 
#' unnecessary covariates, etc.
#' @param lm_input_table the LM input table produced from the above function
#' @param protein_ids_df the input DF of all drivers/ regulatory proteins
#' @param analysis_type 'eQTL' or 'meQTL'
#' @param cna_bucketing a string indicating if/how we are bucketing CNA values. Possible
#' values are "bucketInclAmp", "bucketExclAmp", and "rawCNA"
#' @param meth_bucketing a TRUE/FALSE value indicating whether or not we are bucketing
#' methylation values
#' @param runRegprotsJointly a TRUE/ FALSE value indicating whether or not we are 
#' running all our drivers/ regulatory proteins jointly in one model
construct_formula <- function(lm_input_table, protein_ids_df, analysis_type, 
                              cna_bucketing, meth_bucketing, runRegprotsJointly) {
  
  print("Now constructing formula...")

  # Exclude sample ID
  colnames_to_incl <- colnames(lm_input_table)[2:ncol(lm_input_table)]
  #print(length(colnames_to_incl))
  if(length(colnames_to_incl) < 3) {print(lm_input_table)}
  
  # For the bucketed variables, remove the last bucket (we don't need it!)
  bucketed_vars <- c("Tot_Mut_b","Tumor_purity_b", "Tot_IC_Frac_b", "Cancer_type_b")
  if(meth_bucketing) {
    bucketed_vars <- c(bucketed_vars, "MethStat_k")
    if(!runRegprotsJointly) {
      if(analysis_type == "eQTL") {bucketed_vars <- c(bucketed_vars, "MethStat_i")}
    } else {
      meth_stat_i <- unlist(lapply(protein_ids_df$swissprot, function(x) paste0(x, "_MethStat_i")))
      bucketed_vars <- c(bucketed_vars, meth_stat_i)
    }
  }
  if(!((cna_bucketing == "rawCNA") | (grepl("just", cna_bucketing)))) {
    bucketed_vars <- c(bucketed_vars, "CNAStat_k")
    if(!runRegprotsJointly) {
      bucketed_vars <- c(bucketed_vars, "CNAStat_i")
    } else {
      cna_stat_i <- unlist(lapply(protein_ids_df$swissprot, function(x) paste0(x, "_CNAStat_i")))
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
  
  if (analysis_type == "eQTL") {      
    colnames_to_incl <- colnames_to_incl[!(colnames_to_incl == "ExpStat_k")]
    formula <- paste(colnames_to_incl, collapse = " + ") 
    formula <- paste("ExpStat_k ~ ", formula, sep = "")
  } else if (analysis_type == "meQTL") {
    colnames_to_incl <- colnames_to_incl[!(colnames_to_incl == "MethStat_k")]
    formula <- paste(colnames_to_incl, collapse = " + ") 
    formula <- paste("MethStat_k ~ ", formula, sep = "")
  } else {
    print(paste("Invalid analysis type:", analysis_type))
    return(NA)
  }
  
  return(formula)
}


############################################################
#' A function to combine the collinearity diagnostic statistics,
#' specifically the tolerance, VIF, eigenvalue, and condition index
#' values for each covariate, and returns these aggregate statistics.
#' @param path a path to the folder with all of the collinearity
#' output files
#' @param outfn the file name for the output files, which has all the details of 
#' the current run
#' @param debug if we are in debug mode and should add additional prints
#' NOTE: any condition index over 30 indicates strong collinearity
#' @param randomize whether we have randomized this run or not
combine_collinearity_diagnostics <- function(path, outfn, debug, randomize) {
  vif_tabs <- list.files(path, pattern = "vif_tab")
  eig_cindex_tabs <- list.files(path, pattern = "eig_cindex_tab")
  
  
  # Get only the tables that match the rest of the run as well
  outfn_spl <- unlist(strsplit(outfn, "_", fixed = TRUE))
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
  
  tmp <- fread(paste(path, vif_tabs[1], sep = "/"), header = TRUE)
  results_table <- data.frame(matrix(ncol= 5, nrow = nrow(tmp)+1))
  colnames(results_table) <- c("Covariate", "Tolerance", "VIF", "Eigenvalue", "Condition_Index")
  results_table$Covariate <- c(tmp$Variables, "--")
  
  if(debug) {print(head(results_table))}
    
  for (i in 1:length(vif_tabs)) {
    combo <- paste(unlist(strsplit(vif_tabs[i], "_", fixed = TRUE))[1:2], collapse = "_")
    if(debug) {print(paste("TargGene and Regprot Combo:", combo))}
    vif_tab <- fread(paste(path, vif_tabs[i], sep = "/"), header = TRUE)
    tryCatch({
      eig_cindex_tab <- fread(paste(path, eig_cindex_tabs[grepl(combo, eig_cindex_tabs)], sep = "/"), 
                              header = TRUE)
    }, error=function(cond){
      print(cond)
      tab_name <- eig_cindex_tabs[grepl(combo, eig_cindex_tabs)]
      print("Eig cindex tab name:", tab_name)
      eig_cindex_tab <- fread(paste(path, tab_name[[1]], sep = "/"), 
                              header = TRUE)
    })
      
      
    for(i in 1:nrow(results_table)) {
      var <- results_table$Covariate[i]
      
      if (!var == "--") {
        tolerance <- as.numeric(vif_tab[vif_tab$Variable == var, 'Tolerance'])
        if(debug) {print(paste("Tolerance", tolerance))}
        if(length(tolerance) > 0) {
          results_table[i, 'Tolerance'] <- mean(c(as.numeric(results_table[i, 'Tolerance']), tolerance),
                                                na.rm = TRUE)
        } 
        vif <- as.numeric(vif_tab[vif_tab$Variable == var, 'VIF'])
        if(length(vif) > 0) {
          results_table[i, 'VIF'] <- mean(c(as.numeric(results_table[i, 'VIF']), vif), na.rm = TRUE)
        }
        
      } 
      
      results_table[i, 'Eigenvalue'] <- mean(c(as.numeric(results_table[i, 'Eigenvalue']), 
                                               as.numeric(eig_cindex_tab$Eigenvalue[i])),
                                             na.rm = TRUE)
      results_table[i, 'Condition_Index'] <- mean(c(as.numeric(results_table[i, 'Condition_Index']), 
                                                    as.numeric(eig_cindex_tab$`Condition Index`[i])), 
                                                  na.rm = TRUE)
    }
  }
  return(results_table)

}


############################################################
#' Given a list of  patient IDs (XXXX, e.g. A0WY), it subsets a 
#' given data frame with column/row names that contain sample IDs (XXXX-XXX, e.g. A0WY-01A)
#' to include only those given patients
#' @param patients a vector of patient IDs
#' @param df the data frame to subset using the patient ids
#' @param colNames a TRUE/FALSE value indicating whether the patient/sample IDs are
#' found on the rows or the columns of the DF (TRUE is columns, FALSE is rows)
#' @param tumNormMatched a TRUE/FALSE value indicating whether we are looking at 
#' tumor-normal matched samples (if FALSE, restrict to only cancer samples)
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3'
subset_by_intersecting_ids <- function(patients, df, colNames, tumNormMatched, dataset) {
  df <- as.data.frame(df)
  df_adj <- NA

  column_names_to_keep <- c("ensg_id", "Swissprot", "swissprot", "gene_name", "ensg_ids")
  
  if(colNames == TRUE) {
    # Keep only intersecting patients
    label_col_indices <- which(colnames(df) %fin% column_names_to_keep)
    just_patients <- colnames(df)
    if(dataset == "TCGA") {
      just_patients <- unlist(lapply(colnames(df), function(x) 
        unlist(strsplit(x, "-", fixed = TRUE))[1]))
    }
    df_adj <- df[, c(label_col_indices, which(just_patients %fin% patients))]
    #print(head(df_adj))
    
    # Keep only cancer samples, if TCGA
    if (dataset == "TCGA") {
      if(!tumNormMatched) {
        df_adj_label_col_indices <- which(colnames(df_adj) %fin% column_names_to_keep)
        df_adj <- df_adj[, c(df_adj_label_col_indices, which(grepl("-0", colnames(df_adj), 
                                                                   fixed = TRUE)))]
      }
      
      # Remove any duplicate samples
      df_adj <- df_adj[, !(grepl("-1", colnames(df_adj), fixed = TRUE))]
    }
    
  } else {
    if(length(df) > 0) {
      if(nrow(df) > 0) {
        if(length(unlist(strsplit(as.character(df$sample_id[1]), "-", fixed = TRUE))) == 4) {
          df$sample_id <- unlist(lapply(df$sample_id, function(x) 
            paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))
        }
        # Keep only intersecting patients
        just_patients <- df$sample_id
        if(dataset == "TCGA") {
          just_patients <- unlist(lapply(df$sample_id, function(x) 
            unlist(strsplit(x, "-", fixed = TRUE))[1]))
        }
        df_adj <- df[which(just_patients %fin% patients),]
        
        # Keep only cancer samples
        #if(dataset == "TCGA") {
        #  if(!tumNormMatched) {
        #    df_adj <- df_adj[which(grepl("-0", df_adj$sample_id, fixed = TRUE)),]
        #  }
        #  
        #  # Remove any duplicate samples
        #  df_adj <- df_adj[!(grepl('-1', df_adj$sample_id, fixed = TRUE)),]
        #}
      }
    }
  }
  return(as.data.table(df_adj))
}



############################################################
#' Given a list of intersecting patient IDs (XXXX, e.g. A0WY), it subsets a 
#' given regulatory protein data frame with row entries that contain semicolon-separated
#' sample IDs (XXXX-XXX, e.g. A0WY-01A) to include only those given patients
#' @param intersecting_patients a vector of intersecting patient IDs
#' @param mutation_regprot_df the data frame to subset using the patient ids
#' @param tumNormMatched a TRUE/FALSE value indicating whether we are looking at 
#' tumor-normal matched samples (if FALSE, restrict to only cancer samples)
#' @param dataset either 'TCGA', 'METABRIC', 'ICGC', or 'CPTAC3'
subset_regprot_df_by_intersecting_ids <- function(patients, mutation_regprot_df, 
                                                  tumNormMatched, dataset) {
  
  # Adjust each row to remove patients not in the intersecting patients vector
  regprot_df_new_patient_labels <- lapply(mutation_regprot_df$Patient, function(x) {
    # Split apart the semicolon separated sample IDs
    spl_patients <- unlist(strsplit(x, ";", fixed = TRUE))
    spl_patients_adj <- spl_patients
    # Extract just the 4-digit patient ID, if TCGA
    if(dataset == "TCGA") {
      spl_patients_adj <- unlist(lapply(spl_patients, function(x) 
        unlist(strsplit(x, "-", fixed = TRUE))[1]))
    }
    # Check if this patient is in the intersecting patients list and return TRUE if so, F o.w.
    matching_patient_ind <- unlist(lapply(spl_patients_adj, function(y) {
      if(y %fin% patients) {
        return(TRUE)
      } else {return(FALSE)}
    }))
    # If there are overlapping patients in this row, then we want to keep this row
    if(TRUE %in% matching_patient_ind) {
      samps_to_keep <- spl_patients[matching_patient_ind]
      
      # Keep only cancer samples, if TCGA
      if(dataset == "TCGA") {
        if(!tumNormMatched) {
          samp_ids <- unlist(lapply(samps_to_keep, function(x) unlist(strsplit(x, "-", fixed = TRUE))[2]))
          samps_to_keep <- samps_to_keep[unlist(lapply(samp_ids, function(x) startsWith(x, "0")))]
        }
        
        # Keep only non-duplicate samples
        samps_to_keep <- samps_to_keep[!grepl(".1", samps_to_keep, fixed = TRUE)]
      }
      
      # If we still have samples, return them in this row 
      if(length(samps_to_keep) > 0) {
        return(paste(samps_to_keep, collapse = ";"))
      } else {return(NA)}
    } else {return(NA)}
  })
  
  # Remove NA
  mutation_regprot_df_sub <- mutation_regprot_df[unlist(lapply(regprot_df_new_patient_labels, function(x) 
    ifelse(is.na(x), FALSE, TRUE))),]
  # Update the DF with the new intersecting patient labels
  mutation_regprot_df_sub$Patient <- regprot_df_new_patient_labels[!is.na(regprot_df_new_patient_labels)]
  mutation_regprot_df_sub$Patient <- unlist(lapply(mutation_regprot_df_sub$Patient, function(x) return(unlist(x))))
  
  return(mutation_regprot_df_sub)
}


############################################################
#' Removes metastatic samples from a given data frame (sample ID will be
#' -06 or -07 to indicate metastasis) - only for TCGA
#' @param df the data frame to remove metastatic samples from
#' @param colNames a TRUE/FALSE value indicating whether the patient/sample IDs are
#' found on the rows or the columns of the DF (TRUE is columns, FALSE is rows)
remove_metastatic_samples <- function(df, colNames) {
  df <- as.data.frame(df)
  
  column_names_to_keep <- c("ensg_id", "Swissprot", "swissprot", "gene_name", "ensg_ids", "Patient", "patient")
  
  if(colNames == TRUE) {
    # Keep only primary patients
    label_col_indices <- which(colnames(df) %fin% column_names_to_keep)
    just_samples <- unlist(lapply(colnames(df), function(x) 
      unlist(strsplit(x, "-", fixed = TRUE))[2]))
    df_adj <- df[, unique(c(label_col_indices, intersect(which(just_samples != "06"), 
                                                         which(just_samples != "07"))))]
    #print(head(df_adj))
    
  } else {
    if(length(unlist(strsplit(as.character(df$sample_id[1]), "-", fixed = TRUE))) == 4) {
      df$sample_id <- unlist(lapply(df$sample_id, function(x) 
        paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))
    }
    # Keep only primary patients
    just_samples <- unlist(lapply(df$sample_id, function(x) 
      unlist(strsplit(x, "-", fixed = TRUE))[2]))
    df_adj <- df[intersect(which(just_samples != "06"),
                           which(just_samples != "07")),]
    #print(head(df_adj))
  }
  return(as.data.table(df_adj))
}


############################################################
#' Removes metastatic samples from a given regprot data frame (sample ID will be
#' -06 or -07 to indicate metastasis) - only for TCGA
#' @param mutation_regprot_df the regprot data frame to remove metastatic samples from
remove_metastatic_samples_regprot <- function(mutation_regprot_df) {

  # Adjust each row to remove metastatic samples
  regprot_df_new_patient_labels <- lapply(mutation_regprot_df$Patient, function(x) {
    # Split apart the semicolon separated sample IDs
    spl_patients <- unlist(strsplit(x, ";", fixed = TRUE))
    # Extract just the 2-digit sample ID
    samp_ids <- unlist(lapply(spl_patients, function(x) 
      unlist(strsplit(x, "-", fixed = TRUE))[2]))
    # Check if this sample is primary and return TRUE if so, F o.w.
    primary_samp_ind <- unlist(lapply(samp_ids, function(y) {
      if((y != "06") & (y != "07")) {
        return(TRUE)
      } else {return(FALSE)}
    }))
    # If there are primary samples in this row, then we want to keep this row
    if(TRUE %in% primary_samp_ind) {
      samps_to_keep <- spl_patients[primary_samp_ind]
      
      # If we still have samples, return them in this row 
      if(length(samps_to_keep) > 0) {
        return(paste(samps_to_keep, collapse = ";"))
      } else {return(NA)}
    } else {return(NA)}
  })
  
  # Remove NA
  mutation_regprot_df_sub <- mutation_regprot_df[unlist(lapply(regprot_df_new_patient_labels, function(x) 
    ifelse(is.na(x), FALSE, TRUE))),]
  # Update the DF with the new intersecting patient labels
  mutation_regprot_df_sub$Patient <- regprot_df_new_patient_labels[!is.na(regprot_df_new_patient_labels)]
  mutation_regprot_df_sub$Patient <- unlist(lapply(mutation_regprot_df_sub$Patient, function(x) return(unlist(x))))
  
  return(mutation_regprot_df_sub)
}

############################################################
#' Subset the given targets to only genes that are trans to the current regulatory 
#' protein (found more than X (10E6 default) BP away on the same chromosome or on another 
#' chromosome entirely.
#' @param regprot the regulatory protein in question, ENSG ID
#' @param target_df the full data table of targets we want to subset
#' @param all_genes_id_conv a conversion file from BioMart with information about
#' chromosome number and position
#' @param X the number of bases within which we will consider a gene 'cis'
#' @param debug a TRUE/FALSE value indicating whether or not we are in debug mode
subset_to_trans_targets <- function(regprot, target_df, all_genes_id_conv, debug, X = 10E6) {
  
  # Get the chromosome, starting position, and ending position (in BP) for the 
  # regulatory protein
  chr_regprot <- NA
  start_position_regprot <- NA
  end_position_regprot <- NA
  
  for (r in regprot) {
    chr_regprot <- c(chr_regprot, as.numeric(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == r,
                                                               'chromosome_name'])))
    start_position_regprot <- c(start_position_regprot, 
                                as.numeric(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == r, 
                                                             'start_position'])))
    end_position_regprot <- c(end_position_regprot, 
                              as.numeric(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == r, 
                                                           'end_position'])))
  }
  
  chr_regprot <- unique(chr_regprot[!is.na(chr_regprot)])
  start_position_regprot <- unique(start_position_regprot[!is.na(start_position_regprot)])
  end_position_regprot <- unique(end_position_regprot[!is.na(end_position_regprot)])
  
  if(debug) {
    print(paste("Chr_regprot", chr_regprot))
    print(paste("Start_pos_regprot", start_position_regprot))
    print(paste("End_pos_regprot", end_position_regprot))
  }

  # For each target, check if it is trans or not
  trans_target_df_rows <- lapply(1:target_df[, .N], function(i) {
    
    # Get the target
    t <- unlist(strsplit(target_df$ensg[i], ";", fixed = TRUE))
    
    # Get it's chromosome
    if (length(t) > 1) {t <- t[1]}
    chr_targ <- unique(as.numeric(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == t, 
                                                    'chromosome_name'])))
    
    # If not the same chromosome, return it
    if((length(chr_regprot) == 0) | (length(chr_targ) == 0)) {return(target_df[i,])}
    if(is.na(chr_targ)) {return(target_df[i,])}              # NAs are induced by coersion for X & Y chrom.
    if(chr_regprot == "") {return(target_df[i,])}
    if (!(chr_targ == chr_regprot)) {return(target_df[i,])}
    
    # If on the same chromosome, check the distance apart
    else {
      start_position_targ <- unique(as.numeric(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == t, 
                                                                 'start_position'])))
      end_position_targ <- unique(as.numeric(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == t, 
                                                               'end_position'])))
      
      if(debug) {
        print(paste("Chr_targ", chr_targ))
        print(paste("Start_pos_targ", start_position_targ))
        print(paste("End_pos_targ", end_position_targ))
      }
      
      distance_min <- min(abs(start_position_targ - end_position_regprot),
                          abs(start_position_regprot - end_position_targ),
                          abs(start_position_targ - start_position_regprot),
                          abs(end_position_regprot - end_position_targ))
      if (distance_min > X) {return(target_df[i,])}
    }
  })
 
  trans_target_df <- do.call(rbind, trans_target_df_rows)
  
  return(trans_target_df)
}


############################################################

#' Takes an expression DF and a starter DF and returns the starter DF without outlier patients,
#' defined as patients with expression of target gene that exceeds 1.5*IQR + Q3 threshold
#' (separately among patients with a mutation in the regulatory protein and patients without)
#' @param expression_df an expression DF, unfiltered 
#' @param starter_df a starter DF for the regulatory protein of interest i
#' @param ensg the ensembl ID of the target gene k
filter_expression_df_outliers <- function(expression_df, starter_df, ensg) {
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

