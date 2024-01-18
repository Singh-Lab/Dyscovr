############################################################
### Process Gene Expression (RNA-Seq) Data: FPKM and Counts (Archival)
### ASSOCIATED PUBLICATION INFORMATION
############################################################

library(stringr)
library(edgeR)
library(dplyr)
library(tibble)
library(parallel)
library(gdata)
library(matrixStats)
library(data.table)
library("gplots")

# Archival functions for processing of FPKM-normalized expression data and HT-Seq
# count expression data; alternative normalization methods (e.g. rank-
# normalization, z-score normalization, etc.)

# Local path to directory containing gene expression data files
PATH <- getwd()
OUTPUT_PATH <- paste0(getwd(), "Output/Expression/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)

# Format of output data frame: 
# Rows : Genes (ENSG ID)
# Columns : Sample ID (TCGA ID, ex TCGA-AR-A10Z-01A, where "AR" is cohort ID, 
# "A10Z" is patient ID, "01A" is sample ID); METABRIC ID, e.g. MB-XXXX or 
# MTS-XXXX
# Entries : expression value (may be either normalized or raw counts depending 
# on the type of data we are using)


############################################################
### CREATE COMPILED EXPRESSION DATAFRAME FROM HT-SEQ COUNT 
### DATA OR FPKM-NORMALIZED DATA
############################################################
#' Set up the compiled output data frame using a sample file
#' @param path the local path to the subdirectories of the cohort of interest
#' @param sample_sheet_subset a GDC sample sheet subsetted to only patients 
#' of interest
#' @param counts_or_fpkm a string that denotes that we're looking at 
#' counts ("counts") or FPKM normalization ("fpkm")
create_output_df <- function(path, sample_sheet_subset, counts_or_fpkm) {
  
  patient_filenames <- sample_sheet_subset$File.Name
  if(counts_or_fpkm == "counts") {
    samp_filename <- paste0("Expression/Counts/", patient_filenames[1])
  }
  else if(counts_or_fpkm == "fpkm") {
    samp_filename <- paste0("Expression/FPKM/", patient_filenames[1])}
  else {
    print("Unknown RNA-seq data format, please try again.")
    return(0)
  }
  
  samp_file <- read.table(paste0(path, samp_filename), header = F, sep = "\t")
  gene_ids <- strip_gene_ids(samp_file$V1)
  sample_ids <- sample_sheet_subset$Sample.ID
  
  # Loop through the file names and extract the expression information 
  # for each patient
  expression_list <- mclapply(1:length(patient_filenames), function(i) {
    # Read in the expression file for a given patient
    exp_filename <- patient_filenames[i]
    exp_file <- 0
    
    if(counts_or_fpkm == "fpkm") {
      exp_filename <- sub(".htseq.counts", ".FPKM.txt", exp_filename)
      if(exp_filename %fin% list.files(paste0(path, "Expression/FPKM/"))) {
        exp_file <- fread(paste0(path, paste0("Expression/FPKM/", exp_filename)), 
                                      header = F, sep = "\t")
      }
    } else if (counts_or_fpkm == "counts") {
      if(exp_filename %fin% list.files(paste0(path, "Expression/Counts/"))) {
        exp_file <- fread(paste0(path, paste0("Expression/Counts/", exp_filename)), 
                                      header = F, sep = "\t")
      }
    }
    else {
      print("Unknown RNA-seq data format, please try again.")
      return(0)
    }
    # Fix the gene IDs
    if(!exp_file == 0) {counts <- as.numeric(as.data.frame(exp_file)[,2])}
    else {counts <- rep(NA, times = length(gene_ids))}

    print(paste(i, paste("/", length(patient_filenames))))
    
    return(counts_df)
  })
  expression_dataframe <- as.data.frame(do.call(cbind, expression_list))
  rownames(expression_dataframe) <- gene_ids
  colnames(expression_dataframe) <- sample_ids
  
  return(expression_dataframe)
}

#' Helper function that takes a vector of ENSG IDs and strips the version from 
#' them (the period and the bit after)
#' @param gene_ids a vector of ENSG IDs to be stripped
strip_gene_ids <- function(gene_ids) {
  stripped_ids <- unlist(lapply(gene_ids, FUN = function(id) 
    unlist(strsplit(id, split = ".", fixed = T))[1]))
  return(stripped_ids)
}

expression_df_counts <- create_output_df(PATH, sample_sheet_subset, "counts")
expression_df_fpkm <- create_output_df(PATH, sample_sheet_subset, "fpkm")


# Remove any columns that are entirely NA
expression_df_counts <- expression_df_counts[,which(
  unlist(lapply(expression_df_counts, function(x) !all(is.na(x)))))]
expression_df_fpkm <- expression_df_fpkm[,which(
  unlist(lapply(expression_df_fpkm, function(x) !all(is.na(x)))))]

# Remove the rows that do not start with ENSG
expression_df_counts <- expression_df_counts[grepl(
  "ENSG", rownames(expression_df_counts)),]
expression_df_fpkm <- expression_df_fpkm[grepl(
  "ENSG", rownames(expression_df_fpkm)),]


# Save the expression data frame to a CSV
output_filename_counts <- paste0(OUTPUT_PATH, 
                                 "Expression/expression_counts_DF.csv")
output_filename_fpkm <- paste0(OUTPUT_PATH, "Expression/expression_fpkm_DF.csv")

write.csv(expression_df_counts, file = output_filename_counts, row.names = T)
write.csv(expression_df_fpkm, file = output_filename_fpkm, row.names = T)

#
#
#
# Read back if necessary
expression_df_counts <- read.csv(output_filename_counts, header = T, 
                                 row.names = 1, check.names = F)
expression_df_fpkm <- read.csv(output_filename_fpkm, header = T, 
                               row.names = 1, check.names = F)


############################################################
### FILTER BY DIFFERENTIAL EXPRESSION ACROSS PATIENTS ###
############################################################
#' Filter out genes that have less than a standard deviation of 
#' 1 across all patients in the cohort (e.g. do not significantly vary)
#' @param expression_df compiled expression data frame from create_output_df()
filter_by_sd <- function(expression_df) {
  
  # Get the standard deviation across each row of the DF
  expression_matrix <- as.matrix(expression_df)
  sd_vector <- rowSds(expression_matrix, na.rm = T)
  
  # Limit the matrix to only those rows with a SD less than 1
  expression_matrix_filt <- expression_matrix[which(sd_vector > 1),]
  
  return(as.data.frame(expression_matrix_filt))
}

# Call function
expression_df_filt <- filter_by_sd(expression_df)

# Write to file
write.csv(expression_df_filt, paste0(OUTPUT_PATH, "expression_fpkm_DF_SDGr1.csv"))


############################################################
### FILTER BY MEAN/MEDIAN EXPRESSION ###
############################################################
# Visualize read count distributions
hist(rowMedians(as.matrix(expression_df_counts)), 
     main = "Histogram of Median Read Counts Per Gene", 
     xlab = "Median # of Reads per Gene, Across All Patients")
hist(rowMeans(expression_df_counts), 
     main = "Histogram of Average Read Counts Per Gene", 
     xlab = "Average # of Reads per Gene, Across All Patients")
hist(rowSums(expression_df_counts), 
     main = "Histogram of Total Read Counts Per Gene", 
     xlab = "Sum of Reads per Gene, Across All Patients")


# Filter overall by a given mean or median read count
X <- 10
expression_df_counts_minMedX <- expression_df_counts[rowMedians(as.matrix(
  expression_df_counts)) > X,]
expression_df_counts_minAvgX <- expression_df_counts[rowMeans(
  expression_df_counts) > X,]


# Import the clinical DF from the TCGA
clinical_df <- read.csv(paste0(PATH, "clinical_data.csv"), 
                        header = T, check.names = F)


#' Given an expression data frame and a clinical data frame, keeps only genes 
#' that, per cancer type, exceed a mean or median of X
#' @param patient_cancer_mapping a list of all the 4-digit TCGA patient IDs per 
#' cancer type
#' @param expression_df a merged, pan-cancer expression count DF
#' @param medOrAvg "Med" or "Avg" to indicate if we are using the median 
#' or average
#' @param X a median or average below which we will drop the gene
filter_per_ct <- function(patient_cancer_mapping, expression_df, medOrAvg, X) {
  
  expression_df_filt_cols <- lapply(1:length(patient_cancer_mapping), function(i) {
    patients <- patient_cancer_mapping[[i]]
    cols_to_keep <- unlist(lapply(1:ncol(expression_df), function(j) {
      colnam <- colnames(expression_df)[j]
      pat <- unlist(strsplit(colnam, "-", fixed = T))[3]
      return(ifelse(pat %fin% patients, T, F))
    }))
    expression_df_sub_pats <- expression_df[,cols_to_keep]
    
    if (medOrAvg == "Med") {
      new_df <- expression_df_sub_pats[rowMedians(as.matrix(
        expression_df_sub_pats)) > X, ]
      
    } else if (medOrAvg == "Avg") {
      new_df <- expression_df_sub_pats[rowMeans(as.matrix(
        expression_df_sub_pats)) > X,]
      
    } else {
      print(paste("Invalid analysis type:", medOrAvg))
      return(0)
    }
    new_df$ensg <- rownames(new_df)
    new_df <- as.data.table(new_df)
    return(new_df)
  })
  
  expression_df_counts_filt <- expression_df_filt_cols[[1]]
  for(i in 2:length(expression_df_filt_cols)) {
    expression_df_counts_filt <- merge(expression_df_counts_filt, 
                                       expression_df_filt_cols[[i]], 
                                       all = T, by = 'ensg')
  }
  
  return(expression_df_counts_filt)
}

expression_df_counts_tcga <- expression_df_counts[,grepl(
  "TCGA", colnames(expression_df_counts))]


# Create patient-cancer type matching
specific_types <- unique(clinical_df$project_id)
patient_cancer_mapping <- lapply(specific_types, function(ct) {
  pats <- clinical_df[grepl(ct, clinical_df$project_id), 'case_submitter_id']
  pats_ids <- unique(unlist(lapply(pats, function(pat) 
    unlist(strsplit(pat, "-", fixed = T))[3])))
  return(pats_ids)
})
names(patient_cancer_mapping) <- specific_types

expression_df_counts_minAvg10_perCt <- as.data.frame(filter_per_ct(
  patient_cancer_mapping, expression_df_counts_tcga, "Avg", 10))
expression_df_counts_minMed10_perCt <- as.data.frame(filter_per_ct(
  patient_cancer_mapping, expression_df_counts_tcga, "Med", 10))


write.csv(expression_df_counts_minAvg10_perCt, 
          paste0(OUTPUT_PATH, "expression_counts_DF_mean_Gr10_perCt.csv"))
write.csv(expression_df_counts_minMed10_perCt, 
          paste0(OUTPUT_PATH, "expression_counts_DF_mean_MedGr10_perCt.csv"))


############################################################
### QUNATILE- OR RANK-NORMALIZE EXPRESSION ###
############################################################
#' Quantile normalize an expression data frame and return a normalized version
#' (ranked on a per-gene basis)
#' Function taken from https://bioinformatics.stackexchange.com/questions/6863/how-to-quantile-normalization-on-rna-seq-counts
#' A useful paper on types of quantile normalization for gene expression data: 
#' https://www.nature.com/articles/s41598-020-72664-6#:~:text=The%20quantile%20normalization%20(QN)%20procedure,rank%20with%20this%20average%20value.
#' @param df a raw count expression data frame
quantile_normalize <- function(df){
  df_rank <- apply(df, MARGIN = 2, function(y) rank(y, ties.method="average", 
                                                    na.last = NA))
  df_sorted <- as.data.frame(apply(df, MARGIN = 2, sort))
  row_means <- rowMeans(df_sorted)
  
  index_to_mean <- function(my_index, my_means){
    return(my_means[my_index])
  }
  
  df_final <- as.data.frame(apply(df_rank, MARGIN = 2, function(y) {
    index_to_mean(y, my_means = row_means)
  }))
  rownames(df_final) <- rownames(df)
  return(df_final)
}

#' Apply an inverse normal transformation per gene such that, for each gene, its
#' expression across samples are matched to the quantiles of a standard normal 
#' distribution. Apply to quantile-normalized DF
#' @param df a quantile-normalized expression DF
inverse_normal_transform <- function(df) {
  new_df <- apply(df, MARGIN = 2, function(x) 
    qnorm((rank(x,na.last=NA)-0.5)/sum(!is.na(x))))
  return(new_df)
}


#' Rank-normalizes a given expression count DF
#' @param df the expression count DF to be normalized
rank_normalize <- function(df) {
  df_rank <- as.data.frame(apply(df, MARGIN = 2, function(y) {
    return(rank(y, ties.method = "average") / length(y))
  }))
  return(df_rank)
}


# For any given sample, visualize the distribution of rank-normalized values (should be uniform)
hist(expression_df_rank_norm$`TCGA-AC-A6IV-01A`, main = "", 
     xlab = "Expression of Patient A6IV-01A")

# Count the number of duplicates in each column
duplicate_info_list <- apply(expression_df_counts, MARGIN = 2, function(y) {
  y <- as.integer(y)
  print(paste("Num duplicated", sum(duplicated(y))))
  print(paste("Num 0's", length(y[y==0])))
  print(paste("Num 1's", length(y[y==1])))
  print(paste("Most common values", names(sort(summary(as.factor(y)), 
                                               decreasing = T)[1:10])))
  return(data.frame("Num. Duplicated" = sum(duplicated(y)), 
                    "Num. 0s" = length(y[y==0]),
                    "Num. 1s" = length(y[y==1]), "Most Common Values" = 
                      paste(names(sort(summary(as.factor(y)), 
                                       decreasing = T)[1:10]), collapse = ",")))
})
duplicate_info <- rbindlist(duplicate_info_list)


############################################################
### Z-SCORE NORMALIZATION: RPKM/ FPKM STARTING VALUES ###
############################################################
#' Perform a z-score normalization of the given expression DF
#' Steps: if FPKM values, subtract the overall average gene abundance from the 
#' RPKM-normalized expression for each gene, and divide that result by the 
#' standard deviation (SD) of all of the measured counts across all samples.
#' @param expression_df a gene by sample expression DF with RPKM count values
zscore_normalize <- function(expression_df) {
  exp_rowMeans <- rowMeans(as.matrix(expression_df))
  exp_rowSds <- rowSds(as.matrix(expression_df))
  
  new_cols <- lapply(1:ncol(expression_df), function(i) {
    zscores <- unlist(lapply(1:nrow(expression_df), function(j) {
      gene_rpkm <- expression_df[j, i]
      gene_mean <- exp_rowMeans[j]
      gene_sd <- exp_rowSds[j]
      
      zscore <- (gene_rpkm - gene_mean) / gene_sd
      return(zscore)
    }))
    return(zscores)
  })
  new_df <- do.call(cbind, new_cols)
  colnames(new_df) <- colnames(expression_df)
  rownames(new_df) <- rownames(expression_df)
  
  return(new_df)
}


############################################################
### FILTER COMPILED EXPRESSION DATAFRAME: FPKM ONLY ###
############################################################
#' Process table to remove genes with negligible expression (FPKM)
#' @param expression_df the FPKM expression DF to filter
#' @param thres the threshold for FPKM values to keep (typically 1)
filter_expression_df <- function(expression_df, thres) {
  # Filter out genes with < threshold avg. FPKM across all patients
  expression_df_filt <- expression_df[rowMeans(expression_df) > thres,]
  return(expression_df_filt)
}

fpkm_thres <- 1
expression_df_filt <- filter_expression_df(expression_df_fpkm, fpkm_thres)


# Create a version with gene names rather than ENSG IDs
new_rownames <- unlist(lapply(rownames(expression_df_filt), function(ensg) {
  gn <- paste(unique(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == ensg, 
                                       'external_gene_name']), collapse = ";")
  if (!(gn == "") | (!(grepl("Y_RNA", gn)))) {return(gn)}
  else {return(NA)}
}))
indic_to_delete <- which(is.na(new_rownames))
expression_df_filt <- expression_df_filt[-indic_to_delete,]
rownames(expression_df_filt) <- new_rownames[-indic_to_delete]

# Write to file
write.csv(expression_df_filt, 
          paste0(OUTPUT_PATH, "expression_fpkm_filt_gn_DF.csv"))
