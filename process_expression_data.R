############################################################
### Process Expression Data (RNA-Seq)
### Written By: Sara Geraghty, August 2020
############################################################

library(stringr)
library(edgeR)
library(dplyr)
library(tibble)
library(parallel)
library(TCGAbiolinks)

# This file processes the various expression data files for each patient to:
# Combine the files across various patients into a single data frame that we can 
# easily access in our linear model

# Format of output dataframe: 
# Rows : Genes (ENSG ID)
# Columns : TCGA Sample ID (ex TCGA-AR-A10Z-01A, where "AR" is cohort ID, "A10Z" 
# is patient ID, "01A" is sample ID)
# Entries : expression value (may be either normalized or raw counts depending on the
# type of data we are using)

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/"
#path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/"

output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Expression/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Expression/"

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


############################################################
############################################################
### IMPORT PATIENT TCGA IDS & SAMPLE SHEETS ###
############################################################
############################################################
# Unique patient TCGA IDs
patient_tcga_ids <- read.table(paste(path, "unique_brca_patient_ids.txt", sep = ""), header =  TRUE)[,1]
# patient_tcga_ids <- read.table(paste(path, "unique_patient_ids.txt", sep = ""), header =  TRUE)[,2]

# IF NEEDED FOR PAN-CANCER: Merge multiple files
#' This function merges the sample sheets for pan-cancer expression data, if needed
#' @param path the path to the sample sheets
merge_ss_files <- function(path) {
  sample_sheet1 <- read.csv(paste(path, paste("Expression_Data/Counts/", "gdc_sample_sheet.pt1.csv", sep = ""), sep = ""))
  sample_sheet2 <- read.csv(paste(path, paste("Expression_Data/Counts/", "gdc_sample_sheet.pt2.csv", sep = ""), sep = ""))
  merged_ss <- rbind(sample_sheet1, sample_sheet2)
  return(merged_ss)
}
# Run function
sample_sheet <- merge_ss_files(path)

# Write this back to file 
write.csv(sample_sheet, paste(path, "Expression_Data/Counts/gdc_sample_sheet_total.csv", sep = ""))

# Import the GDC sample sheet for file conversions
sample_sheet <- read.csv(paste(path, "Expression_Data/Counts/gdc_sample_sheet.csv", sep = ""), header = TRUE)
#sample_sheet <- read.csv(paste(path, "Expression_Data/FPKM/gdc_sample_sheet.csv", sep = ""), header = TRUE)
sample_sheet$File.Name <- unlist(lapply(sample_sheet$File.Name, function(x) 
  paste(unlist(strsplit(x, ".", fixed = TRUE))[1:3], collapse = ".")))

#' Subset sample sheet to only patients of interest
#' @param patient_tcga_ids vector of the 4-digit patient TCGA IDs of interest
#' @param sample_sheet the GDC sample sheet we'd like to subset by patients
#' @param version a string (either "original" or "edgeR") which indicates the 
#' sample sheet format
subset_samp_sheet <- function(patient_tcga_ids, sample_sheet, version) {
  indices_of_int <- unlist(lapply(patient_tcga_ids, function(x) {
    if (version == "original") {return(which(grepl(x, sample_sheet$Sample.ID)))}
    else if (version == "edgeR") {return(which(grepl(x, sample_sheet$description)))}
    else {
      print(paste("Unknown version,", version))
      return(NA)
    }
  }))
  sample_sheet_subset <- sample_sheet[indices_of_int,]
  return(sample_sheet_subset)
}

sample_sheet_subset <- subset_samp_sheet(patient_tcga_ids, sample_sheet, "original")

############################################################
############################################################
### FILTERING & NORMALIZATION: COUNTS ONLY ###
############################################################
############################################################

############################################################
### PREPARE INPUT FILES FOR USE WITH EDGER ###
############################################################
# Import a separate version in the format required by the edgeR package
targets <- read.csv(paste(path, "Expression_Data/Counts/gdc_sample_sheet_for_edgeR.csv", sep = ""), 
                    header = TRUE)
targets$files <- unlist(lapply(targets$files, function(x) 
  paste(unlist(strsplit(x, ".", fixed = TRUE))[1:3], collapse = ".")))

targets$files <- unlist(lapply(targets$files, function(x) 
  paste(path, paste("Expression_Data/Counts/", x, sep = ""), sep = "")))

# Subset edgeR-formatted version using the above helper function
targets_sub <- subset_samp_sheet(patient_tcga_ids, targets, "edgeR")

# If needed, add ".txt" at the end or remove ".counts"
# targets_sub$files <- unlist(lapply(targets_sub$files, function(x) 
  #paste(x, ".txt", sep = "")))
# targets_sub$files <- unlist(lapply(targets_sub$files, function(x) 
  #paste(unlist(strsplit(x, ".", fixed = TRUE))[1:2], collapse = "."))

targets_sub$files <- unlist(lapply(targets_sub$files, function(x) {
   if (endsWith(x, ".htseq.counts.txt")) {
     begin <- unlist(strsplit(x, ".", fixed = TRUE))[1]
     return(paste(begin, "htseq.counts", sep = "."))
   } else {return(x)}
 }))

# Limit data frame to only those entries that overlap the files we have
split_fn <- unlist(lapply(targets_sub$files, function(x) unlist(strsplit(x, "/", fixed = TRUE))[11]))
targets_sub2 <- targets_sub[which(unlist(lapply(split_fn, function(x) 
  x %fin% list.files(paste(path, "Expression_Data/Counts/", sep = ""), 
                     pattern = ".htseq")))),]

# If pan-cancer, split apart the file into individual cancer types
#' Takes a targets object (GDC sample sheet in edgeR-accepted format) and breaks
#' it into a list of targets objects, each containing only the samples from the same
#' cancer type
#' @param targets the targets object created from importing a GDC sample sheet
#' in edgeR-accepted format
#' @param dictionary a two-column file that relates each TCGA sample ID to its 
#' corresponding cancer type
split_by_ct <- function(targets, dictionary) {
  # Add a temporary column with the cancer types that we can use to split apart
  # the data frame
  targets$cancer.types <- unlist(lapply(targets$description, function(samp_id) {
    proj_id <- dictionary[dictionary$Sample.ID == samp_id, 'Project.ID']
    ct <- unlist(strsplit(proj_id, "-", fixed = TRUE))[2]
    return(ct)
  }))
  
  # Split the data frame based on the cancer type
  sub_targets <- lapply(unique(targets$cancer.types), function(ct) {
    sub_t <- targets[targets$cancer.types == ct, c('files', 'group', 'description')]
    return(sub_t)
  })
  
  # Add the cancer type name to the list
  names(sub_targets) <- unique(targets$cancer.types)
  
  return(sub_targets)
}

dictionary <- sample_sheet[,c('Sample.ID', 'Project.ID')]
list_of_ct_targets <- split_by_ct(targets_sub2, dictionary)

# Remove any duplicate rows
list_of_ct_targets <- lapply(list_of_ct_targets, distinct)


############################################################
### USE EDGER FOR COUNT NORMALIZATION (TMM) ###
############################################################
#' Takes a subsetted gene targets expression data frame (counts) and, using the
#' edgeR package, TMM normalizes the counts. Returns the TMM-normalized
#' expression data frame.
#' @param targets a count expression data frame, subsetted to only include
#' patients of interest 
normalize_counts_edgeR <- function(targets) {

  d <- readDGE(targets)
  
  # Remove MetaTags
  meta_tags <- grep("^__", rownames(d))
  d <- d[-meta_tags,]
  print(d$samples)
  print(dim(d))
  
  # If doing in two pieces, recombine into one
  #d_full <- c(d, d2)
  
  # Filter out tags that are not expressed in at least four libraries
  keep <- rowSums(cpm(d) > 1) >= 4
  d <- d[keep,]
  # BRCA: went from 60482x840 to 24109x840

  # Reset the library sizes
  d$samples$lib.size <- colSums(d$counts)
  
  # Look at the different sample groups
  print(levels(d$samples$group))
  
  # Reset the colnames to the sample names
  colnames(d$counts) <- d$samples$description
  
  # Perform TMM normalization
  d <- calcNormFactors(d, method = "TMM")
  # Alternatively, use the upper-quartile normalization method ("upperquartile") or the relative log
  # expression normalization method ("RLE")
  
  return(d)
}

# Call this function
d <- normalize_counts_edgeR(targets_sub)
# Or, if we are looking pan-cancer, apply to each list
d_list <- lapply(list_of_ct_targets, normalize_counts_edgeR)
  
#' DGEList Function to get the TMM-normalized counts using TMM factors
cpm.DGEList <- function(y, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25, ...) {
  lib.size <- y$samples$lib.size
  if(normalized.lib.sizes) lib.size <- lib.size * y$samples$norm.factors
  cpm.default(y$counts, lib.size = lib.size, log = log, prior.count = prior.count)
}

tmm <- cpm(d)
# Or, if we are looking pan-cancer, apply to each item in list
tmm_list <- lapply(d_list, cpm)

# Optionally remove the dot from the ENSG IDs
rownames(tmm) <- unlist(lapply(rownames(tmm), function(x) 
  unlist(strsplit(x, ".", fixed = TRUE))[1]))
# Or, if we are looking pan-cancer, apply to each item in list
tmm_list <- lapply(tmm_list, function(tmm) {
  rownames(tmm) <- unlist(lapply(rownames(tmm), function(x) 
    unlist(strsplit(x, ".", fixed = TRUE))[1]))
  return(tmm)
})

############################################################
### OPT: ESTIMATE & VISUALIZE DISPERSION ###
############################################################
#' Estimate and visualize dispersion
#' @param d edgeR output data structure to visualize
visualize_dispersion <- function(d) {
  d <- estimateCommonDisp(d, verbose = TRUE)
  # BRCA: Disp = 0.92427 , BCV = 0.9614
  # Pan-Cancer: 
  d <- estimateTagwiseDisp(d, trend="none")
  plotBCV(d, cex=0.4)
}
visualize_dispersion(d)
# Or, if we are looking pan-cancer, apply to each item in list
lapply(d_list, visualize_dispersion)


############################################################
### WRITE NORMALIZED RESULTS TO FILES ###
############################################################
# Write the results to a file to be saved 
write.csv(d$samples, paste(output_path, "tmm_normalized_expression_samples.csv", sep = ""))
write.csv(d$counts, paste(output_path, "tmm_unnormalized_expression_counts.csv", sep = ""))
write.csv(tmm, paste(output_path, "tmm_normalized_expression_counts.csv", sep = ""))

# Or, if looking pan-cancer, write the results of each cancer type to a separate file
for (i in 1:length(d_list)) {
  d <- d_list[[i]]
  tmm <- tmm_list[[i]]
  write.csv(d$samples, paste(output_path, paste(names(d_list)[i], "tmm_normalized_expression_samples.csv", sep = "_"), sep = ""))
  write.csv(d$counts, paste(output_path, paste(names(d_list)[i], "tmm_unnormalized_expression_counts.csv", sep = "_"), sep = ""))
  write.csv(tmm, paste(output_path, paste(names(d_list)[i], "tmm_normalized_expression_counts.csv", sep = "_"), sep = ""))
}

# Also, merge TMM-normalized expression and samples into "master files" and write these as well
tmm_list <- lapply(tmm_list, as.data.frame)
tmm_list <- lapply(tmm_list, tibble::rownames_to_column)
tmm_normalized_full <- tmm_list %>% reduce(full_join, by = 'rowname')
rownames(tmm_normalized_full) <- tmm_normalized_full$rowname
tmm_normalized_full <- tmm_normalized_full[,2:ncol(tmm_normalized_full)]
write.csv(tmm_normalized_full, paste(output_path, "ALL_CT_tmm_normalized_expression_counts.csv", sep = ""))


samples_list <- lapply(d_list, function(d) return(as.data.frame(d$samples)))
samples_full <- do.call(rbind, samples_list)
write.csv(samples_full, paste(output_path, "ALL_CT_tmm_normalized_expression_samples.csv", sep = ""))


############################################################
############################################################
### CREATE COMPILED EXPRESSION DATAFRAME ###
############################################################
############################################################
#' Set up the compiled output data frame using a sample file
#' @param path the local path to the subdirectories of the cohort of interest
#' @param sample_sheet_subset a GDC sample sheet subsetted to only patients of interest
#' @param counts_or_fpkm a string that denotes that we're looking at counts ("counts")
#' or FPKM normalization ("fpkm")
create_output_df <- function(path, sample_sheet_subset, counts_or_fpkm) {
  
  patient_filenames <- sample_sheet_subset$File.Name
  if(counts_or_fpkm == "counts") {samp_filename <- paste("Expression_Data/Counts/", 
                                                         patient_filenames[1], sep = "")}
  else if(counts_or_fpkm == "fpkm") {samp_filename <- paste("Expression_Data/FPKM/", 
                                                            patient_filenames[1], sep = "")}
  else {
    print("Unknown RNA-seq data format, please try again.")
    return(0)
  }
  samp_file <- read.table(paste(path, samp_filename, sep = ""), header = FALSE, sep = "\t")
  gene_ids <- strip_gene_ids(samp_file$V1)
  sample_ids <- sample_sheet_subset$Sample.ID
  
  # Loop through the file names and extract the expression information for each patient
  expression_list <- mclapply(1:length(patient_filenames), function(i) {
    # Read in the expression file for a given patient
    exp_filename <- patient_filenames[i]
    exp_file <- 0
    
    if(counts_or_fpkm == "fpkm") {
      exp_filename <- sub(".htseq.counts", ".FPKM.txt", exp_filename)
      if(exp_filename %fin% list.files(paste(path, "Expression_Data/FPKM/", sep = ""))) {
        exp_file <- data.table::fread(paste(path, paste("Expression_Data/FPKM/", 
                                                        exp_filename, sep = ""), sep = ""), 
                                      header = FALSE, sep = "\t")
      }
    } else if (counts_or_fpkm == "counts") {
      if(exp_filename %fin% list.files(paste(path, "Expression_Data/Counts/", sep = ""))) {
        exp_file <- data.table::fread(paste(path, paste("Expression_Data/Counts/", 
                                                        exp_filename, sep = ""), sep = ""), 
                                      header = FALSE, sep = "\t")
      }
    }
    else {
      print("Unknown RNA-seq data format, please try again.")
      return(0)
    }
    # Fix the gene IDs
    if(!exp_file == 0) {counts <- as.numeric(as.data.frame(exp_file)[,2])}
    else {counts <- rep(NA, times = length(gene_ids))}
    print(counts)

    print(paste(i, paste("/", length(patient_filenames))))
    
    return(counts)
  })
  expression_dataframe <- as.data.frame(do.call(cbind, expression_list))
  rownames(expression_dataframe) <- gene_ids
  colnames(expression_dataframe) <- sample_ids
  
  return(expression_dataframe)
}

#' HELPER: Takes a vector of ENSG IDs and strips the version from them
#' (the period and the bit after)
#' @param gene_ids a vector of ENSG IDs to be stripped
strip_gene_ids <- function(gene_ids) {
  stripped_ids <- unlist(lapply(gene_ids, FUN = function(id) 
    unlist(strsplit(id, split = ".", fixed = TRUE))[1]))
  return(stripped_ids)
}

expression_dataframe_counts <- create_output_df(path, sample_sheet_subset, "counts")
expression_dataframe_fpkm <- create_output_df(path, sample_sheet_subset, "fpkm")

# Remove any columns that are entirely NA
expression_dataframe_counts <- expression_dataframe_counts[,which(unlist(lapply(expression_dataframe_counts, 
                                                                                function(x) !all(is.na(x)))))]
expression_dataframe_fpkm <- expression_dataframe_fpkm[,which(unlist(lapply(expression_dataframe_fpkm, 
                                                                            function(x) !all(is.na(x)))))]
# Leaves 6825 columns

# Save the expression data frame to a CSV
output_filename_counts <- paste(output_path, "expression_counts_DF.csv", sep = "")
output_filename_fpkm <- paste(output_path, "expression_fpkm_DF.csv", sep = "")

write.csv(expression_dataframe_counts, file = output_filename_counts, row.names = TRUE)
write.csv(expression_dataframe_fpkm, file = output_filename_fpkm, row.names = TRUE)
#
#
#
# Read back if necessary
expression_dataframe_counts <- read.csv(output_filename_counts, header = TRUE, row.names = 1, 
                                      check.names = FALSE, )
expression_dataframe_fpkm <- read.csv(output_filename_fpkm, header = TRUE, row.names = 1, 
                                 check.names = FALSE, )


############################################################
############################################################
### QUNATILE- OR RANK-NORMALIZE EXPRESSION: COUNTS ONLY ###
############################################################
############################################################
# OPT: COUNTS DATA ONLY, Rank or quantile normalize the expression data frame
#' Rank or quantile normalize the expression data frame and return a rank-normalized version
#' (ranked on a per-gene basis)
#' Function taken from https://bioinformatics.stackexchange.com/questions/6863/how-to-quantile-normalization-on-rna-seq-counts
#' @param df a raw count expression data frame
quantile_normalize <- function(df){
  df_rank <- apply(df, MARGIN = 2, function(y) rank(y, ties.method="random"))
  df_sorted <- data.frame(apply(df, MARGIN = 2, sort))
  df_mean <- apply(df_sorted, MARGIN = 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- as.data.frame(apply(df_rank, MARGIN = 2, function(y) {
    index_to_mean(y, my_mean = df_mean)
  }))
  rownames(df_final) <- rownames(df)
  return(df_final)
}

expression_df_quantile_norm <- quantile_normalize(expression_dataframe_counts)

# Write to DF
write.csv(expression_df_quantile_norm, paste(output_path, "expression_quantile_norm_DF.csv", sep = ""))

#' Rank-normalizes a given expression count DF
#' @param df the expression count DF to be normalized
rank_normalize <- function(df) {
  df_rank <- as.data.frame(apply(df, MARGIN = 2, function(y) {
    return(rank(y, ties.method = "random") / length(y))
  }))
  return(as.data.frame(df_rank))
}

expression_df_rank_norm <- rank_normalize(expression_dataframe_counts)

# Write to DF
write.csv(expression_df_rank_norm, paste(output_path, "expression_rank_norm_DF.csv", sep = ""))


############################################################
############################################################
### FILTER COMPILED EXPRESSION DATAFRAME: FPKM ONLY ###
############################################################
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
expression_dataframe_filt <- filter_expression_df(expression_dataframe, fpkm_thres)
  # BRCA: Filters to 16,236 genes

# Save the expression data frame to a CSV
output_filename_fpkm_filt <- paste(output_path, "expression_fpkm_filt_DF.csv", sep = "")
write.csv(expression_dataframe_filt, file = output_filename_fpkm_filt, row.names = TRUE)


# Create a version with gene names rather than ENSG IDs
new_rownames <- unlist(lapply(rownames(expression_dataframe_filt), function(ensg) {
  gn <- paste(unique(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == ensg, 
                                       'external_gene_name']), collapse = ";")
  if (!(gn == "") | (grepl("Y_RNA", gn))) {return(gn)}
  else {return(NA)}
}))
indic_to_delete <- which(is.na(new_rownames))
expression_df_filt <- expression_df_filt[-indic_to_delete,]
rownames(expression_dataframe_filt) <- new_rownames[!is.na(new_rownames)]
write.csv(expression_dataframe_filt, paste(output_path, "expression_fpkm_filt_gn_DF.csv", 
                                           sep = ""))


############################################################
############################################################
### GET A TUMOR-ONLY VERSION OF AN EXPRESSION DF ###
############################################################
############################################################
#' In case a tumor-only version is needed (e.g., for PEER),
#' this function takes an expression DF and limits it to just cancer
#' samples
#' @param expression_df
limit_to_tumor_only <- function(expression_df) {
  cols_to_keep <- which(grepl("-0", colnames(expression_df)))
  expression_df <- expression_df[, cols_to_keep]
  return(expression_df)
}

# Call function
expression_df_tumor_only <- limit_to_tumor_only(expression_df_filt)

write.csv(expression_df_tumor_only, paste(output_path, "expression_fpkm_filt_gn_DF_TO.csv", sep = ""))


#
#
#
#
#
#
#
#
#
#
#
#
#
#
# ALTERNATIVE: USE TCGABiolinks functions for normalization and preprocessing

# Download
exp_query <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      legacy = FALSE,
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts",
                      data.format = "txt",
                      experimental.strategy = "RNA-Seq",
                      access = "open")
GDCdownload(exp_query)
brca_exp_df <- GDCprepare(exp_query, summarizedExperiment = TRUE)

# Import the clinical data
brca_clin <- GDCquery_clinic("TCGA-BRCA", "Clinical")

# Preprocessing
dataPrep_brca <- TCGAanalyze_Preprocessing(object = brca_exp_df,
                                           cor.cut = 0.6,
                                           datatype = "HTSeq - Counts",
                                           filename = "BRCA_IlluminaHiSeq_RNASeq.png")

# Normalization
dataNorm_brca <- TCGAanalyze_Normalization(tabDF = dataPrep_brca, geneInfo = geneInfoHT,
                                           #method = "gcContent")
                                           method = "geneLength")

# Filtering
dataFilt_brca <- TCGAanalyze_Filtering(tabDF = dataNorm_brca, method = "quantile", 
                                       qnt.cut = 0.25)

save(dataFilt_brca, file = paste(output_path, "BRCA_Norm_IlluminaHiSeq.rda"))


