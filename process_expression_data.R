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
library(gdata)
library(matrixStats)
library(data.table)
library("gplots")

# This file processes the various expression data files for each patient to:
# Combine the files across various patients into a single data frame that we can 
# easily access in our linear model

# Format of output dataframe: 
# Rows : Genes (ENSG ID)
# Columns : TCGA Sample ID (ex TCGA-AR-A10Z-01A, where "AR" is cohort ID, "A10Z" 
# is patient ID, "01A" is sample ID)
# Entries : expression value (may be either normalized or raw counts depending on the
# type of data we are using)

# edgeR User's Manual: https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/"
#path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/"

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
patient_tcga_ids <- read.table(paste(path, "unique_brca_patient_ids_2.txt", sep = ""), header =  TRUE)[,1]
# patient_tcga_ids <- read.table(paste(path, "unique_patient_ids_2.txt", sep = ""), header =  TRUE)[,'x']

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
sample_sheet <- read.csv(paste(path, "Expression_Data/Counts/gdc_sample_sheet.csv", sep = ""), header = TRUE, 
                         check.names = FALSE, row.names = 1)
#sample_sheet <- read.csv(paste(path, "Expression_Data/FPKM/gdc_sample_sheet.csv", sep = ""), header = TRUE, 
#check.names = FALSE, row.names = 1)
sample_sheet$`File Name` <- unlist(lapply(sample_sheet$`File Name`, function(x) 
  paste(unlist(strsplit(x, ".", fixed = TRUE))[1:3], collapse = ".")))

#' Subset sample sheet to only patients of interest
#' @param patient_tcga_ids vector of the 4-digit patient TCGA IDs of interest
#' @param sample_sheet the GDC sample sheet we'd like to subset by patients
#' @param version a string (either "original" or "edgeR") which indicates the 
#' sample sheet format
subset_samp_sheet <- function(patient_tcga_ids, sample_sheet, version) {
  indices_of_int <- unlist(lapply(patient_tcga_ids, function(x) {
    if (version == "original") {return(which(grepl(x, sample_sheet$`Sample ID`)))}
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
sample_sheet_subset <- subset_samp_sheet(patient_tcga_ids, sample_sheet, "edgeR")


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
                    header = TRUE, check.names = FALSE)
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

# OPTIONAL: also limit DF to only cancer samples
targets_sub3 <- targets_sub2[targets_sub2$group != 'Solid Tissue Normal',]
#targets_sub3 <- targets_sub3[!grepl("Normal", targets_sub2$group),]

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
list_of_ct_targets <- split_by_ct(targets_sub3, dictionary)

# Remove any duplicate rows
list_of_ct_targets <- lapply(list_of_ct_targets, distinct)


############################################################
### OPT: ADJUST GROUPS TO GROUP BY MUTATION STATUS ###
### OF A GIVEN PROTEIN OF INTEREST ###
############################################################
#' Takes a subsetted edgeR targets DF and changes the grouping
#' to be by the mutation status of a given protein of interest
#' (e.g., TP53) so that we can look at DE genes between the 
#' mutated and nonmutated sample groups.
#' @param targets_sub a subsetted edgeR targets DF
#' @param prot a protein whose mutation status we are 
#' interested in (Hugo Symbol)
#' @param mutation_df a mutation count DF that correlates patient
#' ID to the mutation status of the given protein of interest
group_by_mutation_status <- function(targets_sub, prot, mutation_df) {
  
  # Get the patients with/ without a mutation in the given protein
  mutation_df_sub <- mutation_df[rownames(mutation_df) == prot,]
  patients_mut <- colnames(mutation_df_sub[,(mutation_df_sub[1,] > 0)])
  #patients_nonmut <- setdiff(colnames(mutation_df_sub), patients_mut)
  
  # Add these new groups to the edgeR targets DF under the 'group' label 
  # (patients with and without mutation)
  targets_sub_patients <- unlist(lapply(targets_sub$description, function(x)
    paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))
  targets_sub$group <- unlist(lapply(targets_sub_patients, function(x) 
    ifelse((x %fin% patients_mut), 1, 0)))
  
  return(targets_sub)
}

# Define the protein of interest
prot_of_interest <- "TP53"  # TP53
# prot_of_interest <- "PIK3CA"  # PIK3CA

# Import the mutation count DF
mut_count_matrix <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Mutation/Mutation Count Matrices/mut_count_matrix_missense.csv",
                         header = TRUE, check.names = FALSE, row.names = 1)
# Adjust the column names to keep just the patient ID and not the sample ID
colnames(mut_count_matrix) <- unlist(lapply(colnames(mut_count_matrix), function(x)
  paste(unlist(strsplit(x, "-", fixed = TRUE))[3:4], collapse = "-")))

targets_sub3_groupedByTP53Mut <- group_by_mutation_status(targets_sub3, prot_of_interest, mut_count_matrix)

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
  

  # Filter out tags that are not expressed in at least four libraries
  #keep <- rowSums(cpm(d) > 1) >= 4
  
  # Alternatively, use edgeR's filterByExpr() function
  # Function documentation: https://rdrr.io/bioc/edgeR/man/filterByExpr.html
  keep <- filterByExpr(cpm(d), group = "group") 
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
  #d <- estimateCommonDisp(d, verbose = TRUE)
  # BRCA: Disp = 0.92427 , BCV = 0.9614
  # Pan-Cancer: 
  #d <- estimateTagwiseDisp(d, trend="none")
  #plotBCV(d, cex=0.4)
  
  d <- estimateDisp(d)
  plotBCV(d, cex=0.4)
  
  return(d)
}
visualize_dispersion(d)
# Or, if we are looking pan-cancer, apply to each item in list
lapply(d_list, visualize_dispersion)


############################################################
### OPT: ESTIMATE & VISUALIZE NAIVE DIFF. GENE EXPRESSION ###
############################################################
#' Estimate and visualize top differentially expressed genes,
#' without using a linear modeling approach. Examine the 
#' differences between the DE genes by mutation status of a 
#' given protein (e.g. TP53).  
#' @param d edgeR output data structure to visualize
visualize_de_by_mut_status <- function(d) {
  
  d <- visualize_dispersion(d)
  
  # Get the DE genes across groups (1 is mutated, 0 is not mutated)
  et <- exactTest(d)
  topTags <- topTags(et, n = 15000, adjust.method = "BH", p.value = 0.05)
  print(topTags)
  de1 <- decideTestsDGE(et, adjust.method = "BH", p.value = 0.05)
  print(summary(de1))
  
  # Adjust the topTags table
  rownames(topTags$table) <- unlist(lapply(rownames(topTags$table), function(x) 
    unlist(strsplit(as.character(x), ".", fixed = TRUE))[1]))
  print(head(topTags$table))
  topTags$table$gene.name <- unlist(lapply(rownames(topTags$table), function(ensg)
    paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == ensg, 'external_gene_name'])), collapse = ";")))
  
  return(topTags)
}


topHits_tp53 <- visualize_de_by_mut_status(d_tp53)
# topHits_pik3ca <- visualize_de_by_mut_status(d_pik3ca)

# Or, if we are looking pan-cancer, apply to each item in list
lapply(d_list, visualize_de_by_mut_status)


#' Visualize the top DE genes by mutation status of multiple proteins of interest using a heat map
#' @param topHits_prot1 a formal class TopTags object with the top hits from the first protein of interest
#' @param topHits_prot2 a formal class TopTags object with the top hits from the second protein of interest
#' @param prot1_name the external gene name of the first protein of interest
#' @param prot2_name the external gene name of the second protein of interest
#' @param subset_genes an optional subset of genes that we are interested in looking at expression
#' effects of (e.g. metabolic genes)
create_de_by_mut_status_heatmap <- function(topHits_prot1, topHits_prot2, prot1_name, prot2_name, subset_genes) {
  
  # Get the tables from the TopTags object
  topHits_prot1_tab <- topHits_prot1$table
  topHits_prot2_tab <- topHits_prot2$table
  
  # Limit them to just the gene.name and the logFC
  topHits_prot1_tab <- topHits_prot1_tab[,c("logFC", "gene.name")]
  topHits_prot2_tab <- topHits_prot2_tab[,c("logFC", "gene.name")]
  
  # Adjust the logFC column name
  colnames(topHits_prot1_tab)[which(colnames(topHits_prot1_tab) == "logFC")] <- paste(prot1_name, "logFC", sep = ".")
  colnames(topHits_prot2_tab)[which(colnames(topHits_prot2_tab) == "logFC")] <- paste(prot2_name, "logFC", sep = ".")
  
  # Remove any gene.name entries that are blank
  topHits_prot1_tab <- topHits_prot1_tab[-which(topHits_prot1_tab$gene.name == ""),]
  topHits_prot2_tab <- topHits_prot2_tab[-which(topHits_prot2_tab$gene.name == ""),]
  
  # Merge them together by gene name
  topHits_tab <- merge(topHits_prot1_tab, topHits_prot2_tab, by = "gene.name", all = TRUE)
  
  # Deal with duplicates by adding .X to them (where X is a number)
  duplicates <- unique(topHits_tab$gene.name[which(duplicated(topHits_tab$gene.name))])
  for (d in duplicates) {
    indices_of_d <- which(topHits_tab$gene.name == d)
    for (i in 1:length(indices_of_d)) {
      index <- indices_of_d[i]
      topHits_tab$gene.name[index] <- paste(d, i, sep = ".")
    }
  }
  
  # Make the gene.name a rowname rather than a column for plotting
  rownames(topHits_tab) <- topHits_tab$gene.name
  topHits_tab <- topHits_tab[, which(colnames(topHits_tab) != "gene.name")]
  
  
  # Turn into a matrix
  topHits_matrix <- data.matrix(topHits_tab)
  
  # Plot a default heatmap
  #col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
  #heatmap(matrix, scale = "none", col = col)
  
  # For the purposes of the heatmap, convert all NAs to 0
  topHits_matrix[is.na(topHits_matrix)] <- 0
  
  print(head(topHits_matrix))
  
  # If "subset_genes" is given, limit the matrix to just these genes 
  try(topHits_matrix <- topHits_matrix[which(rownames(topHits_matrix) %in% subset_genes),],
      silent = TRUE)
  
  # Plot an enhanced heatmap
  heatmap.2(topHits_matrix, scale = "none", col = bluered(100), trace = "none",
            density.info = "none", labCol = "", dendrogram = c("row"), 
            add.expr = text(x = seq_along(colnames(topHits_matrix)), 
                            y = -2, srt = 0, labels = colnames(topHits_matrix), 
                            xpd = NA, cex = 2, pos = 1))
  
  return(topHits_matrix)
}

tp53_gn <- "TP53"
pik3ca_gn <- "PIK3CA"

metabolic_genes <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/metabolism_gene_list_ensg.txt",
                              header = TRUE)[,1]
metabolic_genes_gn <- unique(unlist(lapply(metabolic_genes, function(ensg)
  paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == ensg, 'external_gene_name'])), collapse = ";"))))


logFC_matrix <- create_de_by_mut_status_heatmap(topHits_tp53, topHits_pik3ca, tp53_gn, pik3ca_gn)
logFC_matrix_metabolism <- create_de_by_mut_status_heatmap(topHits_tp53, topHits_pik3ca, tp53_gn, pik3ca_gn, metabolic_genes_gn)


############################################################
#### OPT: ESTIMATE SPEARMAN CORRELATION OF LOGFC ###
############################################################
#' Compute and print the Spearman correlation of the Betas and of the T-statistic,
#' given two groups of interest
#' @param results_table the output matrix from "create_de_by_mut_status_heatmap" function
#' @param ri_1 the external gene name of the first protein of interest
#' @param ri_2 the external gene name of the second protein of interest
compute_and_print_spearman <- function(results_table, ri_1, ri_2) {
  
  target_genes <- unique(rownames(results_table))
  
  grp1_logFC <- results_table[,1]
  grp2_logFC <- results_table[,2]
  
  # Get logFC spearman
  logFC_spearman <- cor.test(grp1_logFC, grp2_logFC, method = "spearman")
  logFC_spearman_stat <- as.numeric(logFC_spearman$estimate)
  logFC_spearman_pval <- logFC_spearman$p.value

  # Print the results
  print(paste("Spearman results for", paste(ri_1, paste("and", ri_2))))
  print(paste("logFC correlation of", paste(logFC_spearman_stat, paste(", p-value of", logFC_spearman_pval))))

  # Create a plot to visualize the correlations
  plot(grp1_logFC, grp2_logFC, pch = 19, col = "lightblue", #main = "logFC Spearman Correlation",
       xlab = paste(ri_1, "logFC"), ylab = paste(ri_2, "logFC"))
  abline(lm(grp2_logFC ~ grp1_logFC), col = "red", lwd = 3)
  text(labels = paste("Correlation:", paste(round(logFC_spearman_stat, 6), 
                                            paste(", p-value:", round(logFC_spearman_pval, 6)))), 
       x = max(grp1_logFC, na.rm = TRUE)-sd(grp1_logFC, na.rm = TRUE)*5, 
       y = max(grp2_logFC, na.rm = TRUE)-sd(grp2_logFC, na.rm = TRUE), col = "black")
}

ri_1 <- "TP53"
ri_2 <- "PIK3CA"

# Call function
compute_and_print_spearman(logFC_matrix, ri_1, ri_2)
compute_and_print_spearman(logFC_matrix_metabolism, ri_1, ri_2)


# ASIDE: Compute & visualize overlap between naive DE genes and those predicted from the model
tp53_pik3ca_master_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/output_results_P53_PIK3CA_metabolicTargs_iprotein_tmm_rawCNA_methMRaw_cibersortTotalFrac_allButCancerType_uncorrected_MUT.csv", 
                                  header = TRUE, check.names = FALSE)

# Use processing functions from process_LM_output_manual.R to get top hits from model
tp53_pik3ca_master_df <- mh_correct(tp53_pik3ca_master_df)
tp53_pik3ca_master_df <- add_targ_regprot_gns(tp53_pik3ca_master_df, all_genes_id_conv)
tp53_pik3ca_master_df_sig <- get_signif_correl(tp53_pik3ca_master_df, qval_thres = 0.1)
tophits_tp53 <- tp53_pik3ca_master_df_sig[tp53_pik3ca_master_df_sig$R_i.name == "TP53", 'T_k.name']
tophits_pik3ca <- tp53_pik3ca_master_df_sig[tp53_pik3ca_master_df_sig$R_i.name == "PIK3CA", 'T_k.name']

# Make Venn diagrams
vd_tp53 <- venn.diagram(list(tophits_tp53, logFC_matrix_metabolism$TP53.logFC[logFC_matrix_metabolism$TP53.logFC > 0]),
             category.names = c("Model Hits", "Naive DE Hits"), filename = NULL, output = TRUE,
             lwd = 2, lty = 'blank', fill = c("red", "blue"), cex = 2, fontface = "bold",
             fontfamily = "sans", cat.cex = 2, cat.fontface = "bold",cat.fontfamily = "sans")
grid::grid.draw(vd_tp53)

vd_pik3ca <- venn.diagram(list(tophits_pik3ca, logFC_matrix_metabolism$PIK3CA.logFC[logFC_matrix_metabolism$PIK3CA.logFC > 0]),
                        category.names = c("Model Hits", "Naive DE Hits"), filename = NULL, output = TRUE,
                        lwd = 2, lty = 'blank', fill = c("red", "blue"), cex = 2, fontface = "bold",
                        fontfamily = "sans", cat.cex = 2, cat.fontface = "bold", cat.fontfamily = "sans")
grid::grid.draw(vd_pik3ca)

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
tmm_normalized_full <- tmm_list %>% purrr::reduce(full_join, by = 'rowname')
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
    
    return(counts_df)
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

# Remove the rows that do not start with ENSG
expression_dataframe_counts <- expression_dataframe_counts[grepl("ENSG", rownames(expression_dataframe_counts)),]
expression_dataframe_fpkm <- expression_dataframe_fpkm[grepl("ENSG", rownames(expression_dataframe_fpkm)),]

# OPT: use edgeR to filter out any genes that are lowly expressed
counts_matrix <- as.matrix(expression_dataframe_counts)
keep <- filterByExpr(counts_matrix, group = )
counts_matrix <- counts_matrix[keep,]

expression_dataframe_counts <- as.data.frame(counts_matrix)

# Remove genes that are not expressed in any cell
#expression_dataframe_counts <- expression_dataframe_counts[rowSums(expression_dataframe_counts) > 0,]
#expression_dataframe_fpkm <- expression_dataframe_fpkm[rowSums(expression_dataframe_fpkm) > 0,]

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
                                      check.names = FALSE)
expression_dataframe_fpkm <- read.csv(output_filename_fpkm, header = TRUE, row.names = 1, 
                                 check.names = FALSE)

############################################################
############################################################
### OPT: FILTER BY DIFFERENTIAL EXPRESSION ACROSS PATIENTS ###
############################################################
############################################################
#' Filter out genes that have less than a standard deviation of 
#' 1 across all patients in the cohort (e.g. do not significantly 
#' vary)
#' @param expression_df
filter_by_sd <- function(expression_df) {
  
  # Get the standard deviation across each row of the DF
  expression_matrix <- as.matrix(expression_df)
  sd_vector <- rowSds(expression_matrix, na.rm = TRUE)
  
  print(head(sd_vector))
  
  # Limit the matrix to only those rows with a SD less than 1
  expression_matrix_filt <- expression_matrix[which(sd_vector > 1),]
  
  return(as.data.frame(expression_matrix_filt))
}

expression_dataframe_filt <- filter_by_sd(expression_dataframe)


############################################################
############################################################
### FILTER BY MEAN/MEDIAN EXPRESSION: COUNTS ONLY ###
############################################################
############################################################

# Visualize Read Count Distributions
hist(rowMedians(as.matrix(expression_dataframe_counts)), main = "Histogram of Median Read Counts Per Gene", 
     xlab = "Median # of Reads per Gene, Across All Patients")
hist(rowMeans(expression_dataframe_counts), main = "Histogram of Average Read Counts Per Gene", 
     xlab = "Average # of Reads per Gene, Across All Patients")
hist(rowSums(expression_dataframe_counts), main = "Histogram of Total Read Counts Per Gene", 
     xlab = "Sum of Reads per Gene, Across All Patients")

# Filter overall by a given mean or median read count
X <- 10
expression_dataframe_counts_minMedX <- expression_dataframe_counts[rowMedians(as.matrix(expression_dataframe_counts)) > X,]
expression_dataframe_counts_minAvgX <- expression_dataframe_counts[rowMeans(expression_dataframe_counts) > X,]

# Filter per cancer type
clinical_df <- read.csv(paste(path, "clinical_data_subset_2.csv", sep = ""), 
                        header = TRUE, check.names = FALSE)

#' Given an expression data frame and a clinical data frame, keeps only genes that,
#' per cancer type, exceed a mean or median of X
#' @param patient_cancer_mapping a list of all the 4-digit TCGA patient IDs per cancer type
#' @param expression_df a merged, pan-cancer expression count DF
#' @param medOrAvg "Med" or "Avg" to indicate if we are using the median or average
#' @param X a median or average below which we will drop the gene
filter_per_ct <- function(patient_cancer_mapping, expression_df, medOrAvg, X) {
  
  expression_df_filt_cols <- lapply(1:length(patient_cancer_mapping), function(i) {
    patients <- patient_cancer_mapping[[i]]
    cols_to_keep <- unlist(lapply(1:ncol(expression_df), function(j) {
      colnam <- colnames(expression_df)[j]
      pat <- unlist(strsplit(colnam, "-", fixed = TRUE))[3]
      return(ifelse(pat %fin% patients, TRUE, FALSE))
    }))
    expression_df_sub_pats <- expression_df[,cols_to_keep]
    
    if (medOrAvg == "Med") {
      new_df <- expression_df_sub_pats[rowMedians(as.matrix(expression_df_sub_pats)) > X, ]
      
    } else if (medOrAvg == "Avg") {
      new_df <- expression_df_sub_pats[rowMeans(as.matrix(expression_df_sub_pats)) > X,]
      
    } else {
      print(paste("Invalid analysis type:", medOrAvg))
      return(0)
    }
    new_df$ensg <- rownames(new_df)
    new_df <- as.data.table(new_df)
    return(new_df)
  })
  
  print(head(expression_df_filt_cols))
  #expression_dataframe_counts_filt <- do.call(cbindX, expression_df_filt_cols)
  #expression_dataframe_counts_filt <- expression_df_filt_cols %>% reduce(full_join, by = 'rowname')
  #expression_dataframe_counts_filt <- Reduce(function(x, y) merge(x, y, all=TRUE, by='row.names'), 
                                             #expression_df_filt_cols)
  expression_dataframe_counts_filt <- expression_df_filt_cols[[1]]
  for(i in 2:length(expression_df_filt_cols)) {
    expression_dataframe_counts_filt <- merge(expression_dataframe_counts_filt, 
                                               expression_df_filt_cols[[i]], 
                                              all = TRUE, by = 'ensg')
  }
  
  return(expression_dataframe_counts_filt)
}

expression_dataframe_counts_tcga <- expression_dataframe_counts[,grepl("TCGA", colnames(expression_dataframe_counts))]

# Create patient-cancer type matching
specific_types <- unique(clinical_df$project_id)
patient_cancer_mapping <- lapply(specific_types, function(ct) {
  pats <- clinical_df[grepl(ct, clinical_df$project_id),'case_submitter_id']
  pats_ids <- unique(unlist(lapply(pats, function(pat) unlist(strsplit(pat, "-", fixed = TRUE))[3])))
  return(pats_ids)
})
names(patient_cancer_mapping) <- specific_types

expression_dataframe_counts_minAvg10_perCt <- as.data.frame(filter_per_ct(patient_cancer_mapping, 
                                                            expression_dataframe_counts_tcga,
                                                            "Avg", 10))
expression_dataframe_counts_minMed10_perCt <- as.data.frame(filter_per_ct(patient_cancer_mapping,
                                                            expression_dataframe_counts_tcga,
                                                            "Med", 10))

# Remove any entirely NA rows


write.csv(expression_dataframe_counts_minAvg10_perCt, paste0(output_path, "expression_counts_DF_mean_Gr10_perCt.csv"))
write.csv(expression_dataframe_counts_minMed10_perCt, paste0(output_path, "expression_counts_DF_mean_MedGr10_perCt.csv"))


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
  df_sorted <- as.data.frame(apply(df, MARGIN = 2, sort))
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

#expression_df_quantile_norm <- quantile_normalize(expression_dataframe_counts)
expression_df_quantile_norm <- quantile_normalize(expression_dataframe_counts_minMedX)
expression_df_quantile_norm <- quantile_normalize(expression_dataframe_counts_minAvgX)

# Write to DF
write.csv(expression_df_quantile_norm, paste(output_path, "expression_quantile_norm_DF.csv", sep = ""))


#' Rank-normalizes a given expression count DF
#' @param df the expression count DF to be normalized
rank_normalize <- function(df) {
  df_rank <- as.data.frame(apply(df, MARGIN = 2, function(y) {
    return(rank(y, ties.method = "average") / length(y))
  }))
  return(df_rank)
}

expression_df_rank_norm <- rank_normalize(expression_dataframe_counts)
expression_df_rank_norm <- rank_normalize(expression_dataframe_counts_minMedX)
expression_df_rank_norm <- rank_normalize(expression_dataframe_counts_minAvgX)

# Write to DF
write.csv(expression_df_rank_norm, paste(output_path, "expression_rank_norm_DF.csv", sep = ""))

########################
# For any given sample, visualize the distribution of rank-normalized values (should be uniform)
hist(expression_df_rank_norm$`TCGA-AC-A6IV-01A`, main = "", xlab = "Expression of Patient A6IV-01A")

# Aside: count the number of duplicates in each column
duplicate_info_list <- apply(expression_dataframe_counts, MARGIN = 2, function(y) {
  print(head(y))
  y <- as.integer(y)
  print(paste("Num duplicated", sum(duplicated(y))))
  print(paste("Num 0's", length(y[y==0])))
  print(paste("Num 1's", length(y[y==1])))
  print(paste("Most common values", names(sort(summary(as.factor(y)), decreasing = TRUE)[1:10])))
  return(data.frame("Num. Duplicated" = sum(duplicated(y)), "Num. 0s" = length(y[y==0]),
                    "Num. 1s" = length(y[y==1]), "Most Common Values" = 
                      paste(names(sort(summary(as.factor(y)), decreasing = TRUE)[1:10]), collapse = ",")))
})
duplicate_info <- rbindlist(duplicate_info_list)



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
  if (!(gn == "") | (!(grepl("Y_RNA", gn)))) {return(gn)}
  else {return(NA)}
}))
indic_to_delete <- which(is.na(new_rownames))
# indic_to_delete <- which(new_rownames == "" | new_rownames == "Y_RNA")
expression_df_filt <- expression_df_filt[-indic_to_delete,]
rownames(expression_df_filt) <- new_rownames[-indic_to_delete]
write.csv(expression_df_filt, paste(output_path, "expression_fpkm_filt_gn_DF.csv", 
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


