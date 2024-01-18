############################################################
### Process Gene Expression (RNA-Seq) Data
### ASSOCIATED PUBLICATION INFORMATION
############################################################

library(data.table)
library(edgeR)
library(purrr)
library(dplyr)
library(tibble)
library(parallel)
library(gdata)
library(matrixStats)
library(data.table)
library("gplots")

# Process the STAR count expression data files for each TCGA sample; combine the 
# files across various patients into a single data frame

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

# edgeR User's Manual: https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf


############################################################
### IMPORT PATIENT TCGA IDS & SAMPLE SHEETS ###
############################################################
# Import the GDC sample sheet for file conversions; if pan-cancer, use a merged
# sample sheet
sample_sheet <- read.csv(paste0(PATH, "Expression/gdc_sample_sheet.tsv"), 
                         header = T, check.names = F, sep = "\t")
# Use a sample sheet specifically formatted for use with edgeR
# (re-formatted to have 3 columns: files (filenames), description (TCGA IDs),
# and group (Primary or Metastasis))
sample_sheet_edgeR <- data.frame("files" = sample_sheet$`File Name`, 
                                 "description" = sample_sheet$`Sample ID`,
                                 "group" = sample_sheet$`Sample Type`)
write.csv(sample_sheet_edgeR, paste0(PATH, "Expression/
                                     gdc_sample_sheet_for_edgeR.csv"))
# Read back, as needed
sample_sheet_edgeR <- read.csv(paste0(PATH, "Expression/
                                      gdc_sample_sheet_for_edgeR.csv"), 
                               header = T, check.names = F)

# Import the list of unique patient IDs that have all data types 
patient_id_list <- read.table(paste0(PATH, "UNIQUE_PATIENT_IDS.txt"), 
                              header = T)[,1]

#' Subset sample sheet to only patients of interest
#' @param patient_tcga_ids vector of the 4-digit patient TCGA IDs of interest
#' @param sample_sheet the GDC sample sheet we'd like to subset by patients
#' @param version a string (either "original" or "edgeR") which indicates the 
#' sample sheet format
subset_samp_sheet <- function(patient_tcga_ids, sample_sheet, version) {
  indices_of_int <- unlist(lapply(patient_tcga_ids, function(x) {
    if (version == "original") {
      return(which(grepl(x, sample_sheet$`Sample ID`)))
    }
    else if (version == "edgeR") {
      return(which(grepl(x, sample_sheet$description)))
    }
    else {
      print(paste("Unknown version,", version))
      return(NA)
    }
  }))
  sample_sheet_subset <- sample_sheet[indices_of_int,]
  return(sample_sheet_subset)
}

sample_sheet_subset <- subset_samp_sheet(patient_tcga_ids, sample_sheet, 
                                         "original")
sample_sheet_subset_edgeR <- subset_samp_sheet(patient_tcga_ids, 
                                               sample_sheet_edgeR, "edgeR")



############################################################
### LIMIT TO CANCER SAMPLES ONLY ###
############################################################
# *Proceed with edgeR-formatted files for count data*

# For pan-cancer, combine "Primary Tumor" & "Additional - New Primary" together
sample_sheet_subset_edgeR[which(sample_sheet_subset_edgeR$group == 
                                  "Additional - New Primary"), 
                          'group'] <- "Primary Tumor"

# Remove the solid tissue normal samples
sample_sheet_subset_edgeR_to <- sample_sheet_subset_edgeR[!(sample_sheet_subset_edgeR$group %fin% 
                                                              c('Solid Tissue Normal')), ] 


############################################################
### SPLIT PAN-CANCER SAMPLE SHEET BY CANCER TYPES ###
############################################################
#' If pan-cancer, split apart the file into individual cancer types
#' Takes a targets object (GDC sample sheet in edgeR-accepted format) and breaks
#' it into a list of targets objects, each containing only the samples from 
#' the same cancer type
#' @param targets the targets object created from importing a GDC sample sheet
#' in edgeR-accepted format
#' @param dictionary a two-column file that relates each TCGA sample ID to its 
#' corresponding cancer type
split_by_ct <- function(targets, dictionary) {
  # Add a temporary column with the cancer types that we can use to split apart
  # the data frame
  targets$cancer.types <- unlist(lapply(targets$description, function(samp_id) {
    proj_id <- dictionary[dictionary$`Sample ID` == samp_id, 'Project ID']
    ct <- unlist(strsplit(proj_id, "-", fixed = T))[2]
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

dictionary <- sample_sheet[,c('Sample ID', 'Project ID')]
list_of_ct_targets <- split_by_ct(sample_sheet_subset_edgeR_to, dictionary)

# Remove any duplicate rows
list_of_ct_targets <- lapply(list_of_ct_targets, distinct)


############################################################
### PRE-PROCESS STAR FILES TO CREATE COUNT MATRICES ###
############################################################
LOCAL_PATH <- paste0(PATH, "Expression/")

expression_files <- list.files(LOCAL_PATH, pattern = "star_gene_counts", 
                               recursive = F)

# Subset files to include only necessary columns and remove header lines
for (name in expression_files) {
  print(name)
  fn <- paste0(LOCAL_PATH, name)
  file <- fread(fn, header = F, sep = "\t", check.names = F, skip = 6, 
                select = c(1,4))
  fwrite(file, paste0(LOCAL_PATH, name), col.names = F, sep = "\t")
}


############################################################
### FILTERING & NORMALIZATION USING EDGER ###
############################################################
#' Takes a subsetted gene targets expression data frame (counts) and, using the
#' edgeR package, filters and TMM normalizes the counts. Returns the 
#' TMM-normalized expression data frame.
#' @param targets a count expression data frame, subsetted to only include
#' patients of interest 
#' @param local_path local path to the directory with the expression files
normalize_counts_edgeR <- function(targets, local_path) {
  
  d <- readDGE(targets, path = local_path, columns = c(1,2))
  print(dim(d$counts))

  # Use edgeR's filterByExpr() function to filter reads
  # Function documentation: https://rdrr.io/bioc/edgeR/man/filterByExpr.html
  keep <- filterByExpr(cpm(d), group = "group") 
  d <- d[keep,]
  
  # Reset the library sizes
  d$samples$lib.size <- colSums(d$counts)
  
  # Reset the colnames to the sample names
  colnames(d$counts) <- d$samples$description
  
  # Perform TMM normalization
  d <- calcNormFactors(d, method = "TMM")
  
  return(d)
}


#' DGEList function to get the TMM-normalized counts using TMM factors
cpm.DGEList <- function(y, normalized.lib.sizes=T, log=F, prior.count=0.25, ...) {
  lib.size <- y$samples$lib.size
  if(normalized.lib.sizes) lib.size <- lib.size * y$samples$norm.factors
  cpm.default(y$counts, lib.size = lib.size, log = log, prior.count = prior.count)
}


# For a single cancer type, call this function
d <- normalize_counts_edgeR(sample_sheet_subset_edgeR_to, LOCAL_PATH)

# Or, if we are looking pan-cancer, apply to each in a list
d_list <- lapply(list_of_ct_targets, function(x) {
  tryCatch({
    d <- normalize_counts_edgeR(x, paste0(PATH, "Expression_Data/Counts/"))
    return(d)
  }, error = function(cond) {
    print(cond)
  })
})


tmm <- cpm(d)
# Or, if we are looking pan-cancer, apply to each item in list
tmm_list <- lapply(d_list, cpm)


# Strip the version numbers from the ENSG IDs
rownames(tmm) <- unlist(lapply(rownames(tmm), function(x) 
  unlist(strsplit(x, ".", fixed = T))[1]))
# Or, if we are looking pan-cancer, apply to each item in list
tmm_list <- lapply(tmm_list, function(tmm) {
  rownames(tmm) <- unlist(lapply(rownames(tmm), function(x) 
    unlist(strsplit(x, ".", fixed = T))[1]))
  return(tmm)
})


############################################################
### OPTIONAL: ESTIMATE & VISUALIZE DISPERSION ###
############################################################
#' Estimate and visualize dispersion
#' @param d edgeR output data structure to visualize
visualize_dispersion <- function(d) {
  
  d <- estimateDisp(d)
  plotBCV(d, cex=0.4)
  
  return(d)
}

visualize_dispersion(d)
# Or, if we are looking pan-cancer, apply to each item in list
lapply(d_list, visualize_dispersion)


############################################################
### WRITE NORMALIZED RESULTS TO FILES ###
############################################################
# For a single cancer type
write.csv(d$samples, paste0(OUTPUT_PATH, "samples.csv"))
write.csv(d$counts, paste0(OUTPUT_PATH, "expression_counts_DF_edgeRfilt_TO.csv"))
write.csv(tmm, paste0(OUTPUT_PATH, "tmm_normalized_expression_counts.csv"))

# Pan-cancer
for (i in 1:length(d_list)) {
  d <- d_list[[i]]
  tmm <- tmm_list[[i]]
  write.csv(d$samples, paste0(OUTPUT_PATH, paste0(names(d_list)[i], 
                                                  "samples.csv")))
  write.csv(d$counts, paste0(OUTPUT_PATH, 
                             paste0(names(d_list)[i], 
                                    "expression_counts_DF_edgeRfilt_TO.csv")))
  write.csv(tmm, paste0(OUTPUT_PATH, 
                        paste0(names(d_list)[i], 
                               "tmm_normalized_expression_counts.csv")))
}


# Merge per-cancer count files into one pan-cancer count file
d_list_counts <- lapply(d_list, function(d) as.data.frame(d$counts))
d_list_counts <- lapply(d_list_counts, tibble::rownames_to_column)
d_counts_full <- d_list_counts %>% purrr::reduce(full_join, by = 'rowname')
rownames(d_counts_full) <- d_counts_full$rowname
d_counts_full <- d_counts_full[,2:ncol(d_counts_full)]
write.csv(d_counts_full, paste0(OUTPUT_PATH, 
                               "ALL_CT_expression_counts_DF_edgeRfilt_TO.csv"))


# Merge per-cancer TMM-normalized files into one pan-cancer count file
tmm_list <- lapply(tmm_list, as.data.frame)
tmm_list <- lapply(tmm_list, tibble::rownames_to_column)
tmm_normalized_full <- tmm_list %>% purrr::reduce(full_join, by = 'rowname')
rownames(tmm_normalized_full) <- tmm_normalized_full$rowname
tmm_normalized_full <- tmm_normalized_full[,2:ncol(tmm_normalized_full)]
write.csv(tmm_normalized_full, paste(OUTPUT_PATH, 
                                     "ALL_CT_tmm_normalized_expression.csv"))

# Merge per-cancer sample files into one pan-cancer count file
samples_list <- lapply(d_list, function(d) return(as.data.frame(d$samples)))
samples_full <- do.call(rbind, samples_list)
write.csv(samples_full, paste0(OUTPUT_PATH, 
                               "ALL_CT_tmm_normalized_expression_samples.csv"))

# Also write tumor-only sample file
samples_full_to <- samples_full[grepl("-0", samples_full$description),]
write.csv(samples_full_to, paste(OUTPUT_PATH, 
                                 "ALL_CT_tmm_normalized_expression_samples_TO.csv"))


### RUN QUANTILE-NORMALIZATION PIPELINE ON PREPROCESSED COUNT FILES USING 
### QUANTILE_NORMALIZE_EXPRESSION.PY 



############################################################
### POST-PROCESSING OF QUANTILE-NORMALIZED FILES ###
############################################################
# Recombine expression files that have been individually quantile-normalized 
# per cancer type into one pan-cancer file
qn_expression_fns <- list.files(OUTPUT_PATH, 
                                pattern = "quantile_normalized.csv")
qn_expression_dfs <- lapply(qn_expression_fns, function(x) {
  f <- fread(paste0(OUTPUT_PATH, x), header = T)
  colnames(f)[1] <- "ensg_id"
  return(f)
})

qn_expression_df_full <- qn_expression_dfs %>% purrr::reduce(full_join, 
                                                             by = 'ensg_id')
rownames(qn_expression_df_full) <- qn_expression_df_full$ensg_id
qn_expression_df_full <- qn_expression_df_full[, which(colnames(
  qn_expression_df_full) != 'ensg_id'), with = F]
write.csv(qn_expression_df_full, 
          paste0(OUTPUT_PATH, "ALL_CT_expression_quantile_normalized_sklearn.csv"))
