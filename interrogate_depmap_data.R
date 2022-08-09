###################################################################################
### INTERROGATE DEPMAP DATA
### By: Sara Geraghty, Dec. 2021
###################################################################################

library(ggplot2)
library(reshape2)
library(ggpubr)
library(rstatix)
library(broom)
library(ggrepel)
library(rstatix)

# Functions to investigate drug sensitivity, knockout (dependency) information, and other phenotypes 
# in the DepMap cell lines based on cell line mutation/ copy number status

# Also includes a similar section on interrogating drug sensitivity data from PharmacoDB

###################################################################################
# DRUG SENSITIVITY ANALYSIS (DEPMAP)
###################################################################################
# Import the drug sensitivity and CCLE mutation files
metformin_sensitivity <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/depmap_metformin.csv", 
                                  header = TRUE, check.names = FALSE)
metformin_sensitivity_brca <- metformin_sensitivity[grepl("Breast", metformin_sensitivity$Lineage),]

methotrexate_sensitivity <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/methotrexate_sensitivity_primaryPRISM.csv", 
                                     header = TRUE, check.names = FALSE)
methotrexate_sensitivity <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/methotrexate_sensitivity_secondaryPRISM.csv", 
                                     header = TRUE, check.names = FALSE)
methotrexate_sensitivity <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/methotrexate_sensitivity_auc.csv", 
                                     header = TRUE, check.names = FALSE)
methotrexate_sensitivity_brca <- methotrexate_sensitivity[grepl("Breast", methotrexate_sensitivity$Lineage),]

fu_sensitivity_brca <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/5-fluorouracil_sensitivity_brca_primaryPRISM.csv", 
                             header = TRUE, check.names = FALSE)

ccle_drug_sensitivity <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/DepMap CCLE/drug_sensitivity_screen.csv", 
                                  header = TRUE, check.names = FALSE)

ccle_mutations <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/DepMap CCLE/CCLE_mutations.csv", 
                           header = TRUE, check.names = FALSE)
ccle_missense_mutations <- ccle_mutations[ccle_mutations$Variant_Classification == "Missense_Mutation",]


#' Analyze differential sensitivity between mutated/ nonmutated cell lines (of a particular given gene)
#' using a boxplot and corresponding Wilcoxon test
#' @param sensitivity_file a DepMap sensitivity file for a drug of interest, subsetted to the cell lines of interest
#' @param mutation_file a DepMap mutation file
#' @param gene the gene whose mutation status is of interest (Hugo Symbol)
analyze_differential_sensitivity <- function(sensitivity_file, mutation_file, gene) {
  # Subset the mutation file to the gene of interest
  mutation_file_sub <- mutation_file[mutation_file$Hugo_Symbol == gene,]
  
  # Get the cell lines with a mutation in this gene
  cell_lines_mut <- unique(mutation_file_sub$DepMap_ID)
  
  # Add this information to the drug sensitivity file
  sensitivity_file$Mut <- unlist(lapply(sensitivity_file$`Depmap ID`, function(id) {
    if(id %in% cell_lines_mut) {return(1)}
    else {return(0)}
  }))
  
  # Plot as a boxplot
  sensitivity_mutants <- unlist(sensitivity_file[sensitivity_file$Mut == 1, 2])
  sensitivity_normal <- unlist(sensitivity_file[sensitivity_file$Mut == 0, 2])
  boxplot(list("Mutation" = sensitivity_mutants, "No Mutation" = sensitivity_normal), 
          ylab = "Drug sensitivity", main = paste("Sensitivity by", paste(gene, "Missense Mutation Status")))
  
  # Wilcoxon
  print(wilcox.test(sensitivity_mutants, y = sensitivity_normal))
}

# Call function
analyze_differential_sensitivity(metformin_sensitivity_brca, ccle_missense_mutations, "TP53")
analyze_differential_sensitivity(metformin_sensitivity_brca, ccle_missense_mutations, "PIK3CA")

analyze_differential_sensitivity(methotrexate_sensitivity_brca, ccle_missense_mutations, "TP53")
analyze_differential_sensitivity(methotrexate_sensitivity_brca, ccle_missense_mutations, "PIK3CA")

analyze_differential_sensitivity(fu_sensitivity_brca, ccle_missense_mutations, "TP53")
analyze_differential_sensitivity(fu_sensitivity_brca, ccle_missense_mutations, "PIK3CA")


#' Analyze differential sensitivity across two genes, considering cases where cell lines have mutations
#' exclusively in one or the other, both, or neither. Plots a boxplot and performs a corresponding 
#' ANOVA test
#' @param sensitivity_file a DepMap sensitivity file for a drug of interest, subsetted to the cell lines of interest
#' @param mutation_file a DepMap mutation file
#' @param gene1 the first gene whose mutation status is of interest (Hugo Symbol)
#' @param gene2 the second gene whose mutation status is of interest (Hugo Symbol)
analyze_differential_sensitivity_twoGenes <- function(sensitivity_file, mutation_file, gene1, gene2) {
  # Subset the mutation file to the genes of interest
  mutation_file_gene1 <- mutation_file[mutation_file$Hugo_Symbol == gene1,]
  mutation_file_gene2 <- mutation_file[mutation_file$Hugo_Symbol == gene2,]
  
  # Get the cell lines with a mutation in each gene
  cell_lines_mut_gene1 <- unique(mutation_file_gene1$DepMap_ID)
  cell_lines_mut_gene2 <- unique(mutation_file_gene2$DepMap_ID)
  
  # Get the cell line overlap
  cell_lines_mut_gene1_gene2 <- intersect(cell_lines_mut_gene1, cell_lines_mut_gene2)
  
  # Add this information to the drug sensitivity file
  sensitivity_file$MutG1 <- unlist(lapply(sensitivity_file$`Depmap ID`, function(id) {
    if(id %in% cell_lines_mut_gene1) {return(1)}
    else {return(0)}
  }))
  sensitivity_file$MutG2 <- unlist(lapply(sensitivity_file$`Depmap ID`, function(id) {
    if(id %in% cell_lines_mut_gene2) {return(1)}
    else {return(0)}
  }))
  
  # ANOVA
  table_for_anova <- data.frame("sensitivity" = sensitivity_file[,2], "status" = unlist(lapply(1:nrow(sensitivity_file), function(i) {
    if(sensitivity_file$MutG1[i] & sensitivity_file$MutG2[i]) {return("Both")}
    else if(!sensitivity_file$MutG1[i] & !sensitivity_file$MutG2[i]) {return("Neither")}
    else if (sensitivity_file$MutG1[i] & !sensitivity_file$MutG2[i]) {return(paste(as.character(gene1), "Mut Only"))}
    else if (!sensitivity_file$MutG1[i] & sensitivity_file$MutG2[i]) {return(paste(as.character(gene2), "Mut Only"))}
    else {return(NA)}
  })))
  aov_res <- aov(sensitivity ~ status, data = table_for_anova)
  print(summary(aov_res))
  
  boxplot(sensitivity ~ status, data = table_for_anova, ylab = "Drug sensitivity", 
          main = "Sensitivity to 5-FU in BRCA Cell Lines, by Missense Mutation Status")
  
  return(table_for_anova)
 
}

# Call function
tab <- analyze_differential_sensitivity_twoGenes(metformin_sensitivity_brca, ccle_missense_mutations, "TP53", "PIK3CA")
tab <- analyze_differential_sensitivity_twoGenes(methotrexate_sensitivity_brca, ccle_missense_mutations, "TP53", "PIK3CA")
tab <- analyze_differential_sensitivity_twoGenes(fu_sensitivity_brca, ccle_missense_mutations, "TP53", "PIK3CA")


# Aside: Test for a significant difference between two genes of interest
ttest_res <- t.test(x = tab[tab$status == "TP53 Mut Only", 'sensitivity'], y = tab[tab$status == "PIK3CA Mut Only", 'sensitivity'])


###################################################################################
# DRUG SENSITIVITY ANALYSIS (PHARMACODB)
###################################################################################
# PharmacoDB Web Portal: https://pharmacodb.pmgenomics.ca/
# AAC: Area Above Curve (referring to a dose-response curve); DSS: Drug Sensitivity Score (normalized against dose)

# Import PharmacoDB data downloaded from Web Portal

# 5-FU Data
fu_sensitivity <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/PharmacoDB/breast_5-FU_dose_response.csv", 
                              header = TRUE, check.names = FALSE)
fu_stats <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/PharmacoDB/breast_5-FU_stat_table.csv", 
                     header = TRUE, check.names = FALSE)

# Methotrxate Data
methotrexate_sensitivity <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/PharmacoDB/breast_Methotrexate_dose_response.csv", 
                           header = TRUE, check.names = FALSE)
methotrexate_stats <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/PharmacoDB/breast_Methotrexate_stat_table.csv", 
                     header = TRUE, check.names = FALSE)

# Import CCLE's mutation data for BRCA cell lines in TP53 and PIK3CA
ccle_mutations <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Mutation_21Q4_Public_P53_PIK3CA.csv", 
                           header = TRUE, check.names = FALSE)

# Get overlapping cell lines
print(length(intersect(fu_stats$`Cell Line`, ccle_mutations$cell_line_display_name)))  #21
fu_stats_sub <- fu_stats[fu_stats$`Cell Line` %in% ccle_mutations$cell_line_display_name,]

print(length(intersect(methotrexate_stats$`Cell Line`, ccle_mutations$cell_line_display_name)))  #20
methotrexate_stats_sub <- methotrexate_stats[methotrexate_stats$`Cell Line` %in% ccle_mutations$cell_line_display_name,]



#' Analyze differential sensitivity across two genes, considering cases where cell lines have mutations
#' exclusively in one or the other, both, or neither. Plots a boxplot and performs a corresponding 
#' ANOVA test. Modified for PharmacoDB files.
#' @param sensitivity_file a PharmacoDB sensitivity file for a drug of interest, subsetted to the cell lines of interest
#' @param mutation_file a CCLE mutation file
#' @param gene1 the first gene whose mutation status is of interest (Hugo Symbol)
#' @param gene2 the second gene whose mutation status is of interest (Hugo Symbol)
analyze_differential_sensitivity_twoGenes_pharmacoDB <- function(sensitivity_file, mutation_file, gene1, gene2) {
  
  # Get the cell lines with a mutation in each gene
  cell_lines_mut_gene1 <- unique(mutation_file[mutation_file[ncol(mutation_file)-1] == 1,
                                                             'cell_line_display_name'])
  cell_lines_mut_gene2 <- unique(mutation_file[mutation_file[ncol(mutation_file)] == 1,
                                               'cell_line_display_name'])
  
  # Get the cell line overlap
  cell_lines_mut_gene1_gene2 <- intersect(cell_lines_mut_gene1, cell_lines_mut_gene2)
  
  # Add this information to the drug sensitivity file
  sensitivity_file$MutG1 <- unlist(lapply(sensitivity_file$`Cell Line`, function(id) {
    if(id %in% cell_lines_mut_gene1) {return(1)}
    else {return(0)}
  }))
  sensitivity_file$MutG2 <- unlist(lapply(sensitivity_file$`Cell Line`, function(id) {
    if(id %in% cell_lines_mut_gene2) {return(1)}
    else {return(0)}
  }))
  print(head(sensitivity_file))
  
  # ANOVA
  table_for_anova <- data.frame("sensitivity" = sensitivity_file[,'DSS1 (arb.)'], 
                                "status" = unlist(lapply(1:nrow(sensitivity_file), function(i) {
    if(sensitivity_file$MutG1[i] & sensitivity_file$MutG2[i]) {return("Both")}
    else if(!sensitivity_file$MutG1[i] & !sensitivity_file$MutG2[i]) {return("Neither")}
    else if (sensitivity_file$MutG1[i] & !sensitivity_file$MutG2[i]) {return(paste(as.character(gene1), "Mut Only"))}
    else if (!sensitivity_file$MutG1[i] & sensitivity_file$MutG2[i]) {return(paste(as.character(gene2), "Mut Only"))}
    else {return(NA)}
  })))
  print(table_for_anova)
  aov_res <- aov(sensitivity ~ status, data = table_for_anova)
  print(summary(aov_res))
  
  boxplot(sensitivity ~ status, data = table_for_anova, ylab = "Drug sensitivity", 
          main = "Sensitivity to Methotrexate in BRCA Cell Lines, by Missense Mutation Status")
  
  return(table_for_anova)
  
}

# Call function
tab <- analyze_differential_sensitivity_twoGenes_pharmacoDB(fu_stats_sub, ccle_mutations, "TP53", "PIK3CA")
tab <- analyze_differential_sensitivity_twoGenes_pharmacoDB(methotrexate_stats_sub, ccle_mutations, "TP53", "PIK3CA")


###################################################################################
# CRISPR KO DEPENDENCY DATA, INDIV. GENES
###################################################################################
mthfd1L_crispr <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/MTHFD1L CRISPR (DepMap 21Q4 Public+Score Chronos).csv", header= TRUE)
mthfd1L_crispr_brca <- mthfd1L_crispr[mthfd1L_crispr$Lineage == "Breast",]

expression_cell_lines_tp53_mut <- mthfd1L_crispr_brca[mthfd1L_crispr_brca$Depmap.ID %in% intersect(cell_lines_tp53_mut, mthfd1L_crispr_brca$Depmap.ID), 'Expression.21Q4.Public']
expression_cell_lines_pik3ca_mut <- mthfd1L_crispr_brca[mthfd1L_crispr_brca$Depmap.ID %in% intersect(cell_lines_pik3ca_mut, mthfd1L_crispr_brca$Depmap.ID), 'Expression.21Q4.Public']

boxplot(list("TP53 Mut" = expression_cell_lines_tp53_mut, "PIK3CA Mut" = expression_cell_lines_pik3ca_mut), 
         ylab = "Expression of MTHFD1L", main = "MTHFD1 Expression, by Missense Mutation Status")

cell_lines_tp53_mut_ko <- intersect(cell_lines_tp53_mut, mthfd1L_crispr_brca$Depmap.ID)
cell_lines_pik3ca_mut_ko <- intersect(cell_lines_pik3ca_mut, mthfd1L_crispr_brca$Depmap.ID)
overlap <- intersect(cell_lines_tp53_mut_ko, cell_lines_pik3ca_mut_ko)

cell_lines_tp53_mut_ko_excl <- setdiff(cell_lines_tp53_mut_ko, overlap)
cell_lines_pik3ca_mut_ko_excl <- setdiff(cell_lines_pik3ca_mut_ko, overlap)

expression_cell_lines_tp53_mut_excl <- mthfd1L_crispr_brca[mthfd1L_crispr_brca$Depmap.ID %in% cell_lines_tp53_mut_ko_excl, 'Expression.21Q4.Public']
expression_cell_lines_pik3ca_mut_excl <- mthfd1L_crispr_brca[mthfd1L_crispr_brca$Depmap.ID %in% cell_lines_pik3ca_mut_ko_excl, 'Expression.21Q4.Public']

boxplot(list("TP53 Mut Excl" = expression_cell_lines_tp53_mut_excl, "PIK3CA Mut Excl" = expression_cell_lines_pik3ca_mut_excl), 
         ylab = "Expression of MTHFD1L", main = "MTHFD1 Expression, by Missense Mutation Status")

ko_score_tp53_mut <- mthfd1L_crispr_brca[mthfd1L_crispr_brca$Depmap.ID %in% cell_lines_tp53_mut_ko, 'CRISPR..DepMap.21Q4.Public.Score..Chronos.']
ko_score_pik3ca_mut <- mthfd1L_crispr_brca[mthfd1L_crispr_brca$Depmap.ID %in% cell_lines_pik3ca_mut_ko, 'CRISPR..DepMap.21Q4.Public.Score..Chronos.']

boxplot(list("TP53 Mut" = ko_score_tp53_mut, "PIK3CA Mut" = ko_score_pik3ca_mut), 
         ylab = "Expression of MTHFD1L", main = "MTHFD1 KO Score, by Missense Mutation Status")

ko_score_tp53_mut_excl <- mthfd1L_crispr_brca[mthfd1L_crispr_brca$Depmap.ID %in% cell_lines_tp53_mut_ko_excl, 'CRISPR..DepMap.21Q4.Public.Score..Chronos.']
ko_score_pik3ca_mut_excl <- mthfd1L_crispr_brca[mthfd1L_crispr_brca$Depmap.ID %in% cell_lines_pik3ca_mut_ko_excl, 'CRISPR..DepMap.21Q4.Public.Score..Chronos.']

boxplot(list("TP53 Mut Excl" = ko_score_tp53_mut_excl, "PIK3CA Mut Excl" = ko_score_pik3ca_mut_excl), 
         ylab = "Expression of MTHFD1L", main = "MTHFD1 KO Score, by Missense Mutation Status")


###################################################################################
# CRISPR KO DEPENDENCY DATA, ALL SIGNIF. GENES
###################################################################################

### IMPORT DEPMAP PORTAL FEATURES ###
# Link to DepMap downloads: https://depmap.org/portal/download/

# Use DepMap features to create an aggregate file for all BRCA cell lines and for a particular gene set of interest
signif_hits_crispr_data <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/CRISPR_(DepMap_21Q4_Public+Score,_Chronos)_subsetted.csv",
                                header = TRUE, check.names = FALSE)
signif_hits_rnai_data <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/RNAi_(Achilles+DRIVE+Marcotte,_DEMETER2)_subsetted.csv",
                                  header = TRUE, check.names = FALSE)
#signif_hits_expression_data <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Expression_21Q4_Public_subsetted.csv",
                                        #header = TRUE, check.names = FALSE)
# For all genes
signif_hits_expression_data <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Expression_22Q1_Public.csv", header = TRUE, check.names = FALSE)

# Get the mutations for TP53 & PIK3CA among these cell lines
tp53_pik3ca_mutations <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Mutation_21Q4_Public_P53_PIK3CA.csv",
                                  header = TRUE, check.names = FALSE)

# Get the relative copy number changes for TP53 & PIK3CA among these cell lines
tp53_pik3ca_cnas <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Copy_Number_21Q4_Public_subsetted.csv",
                             header = TRUE, check.names = FALSE)
# Adjust column names
colnames(tp53_pik3ca_cnas)[ncol(tp53_pik3ca_cnas)-1] <- "TP53.CNA"
colnames(tp53_pik3ca_cnas)[ncol(tp53_pik3ca_cnas)] <- "PIK3CA.CNA"

# Get the absolute copy number for TP53 & PIK3CA among these cell lines
#tp53_pik3ca_cnas <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Copy_Number_(Absolute)_subsetted.csv",
#header = TRUE, check.names = FALSE)

### IMPORT CELL MODEL PASSPORT FEATURES ###
# Link to Cell Model Passport downloads: https://cellmodelpassports.sanger.ac.uk/downloads

depmap_priority_scores <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Cell Model Passports/depmap-priority-scores.csv",
                                   header = TRUE, check.names = FALSE)
depmap_priority_scores_brca <- depmap_priority_scores[grepl("Breast", depmap_priority_scores$`analysis name`),]

cell_passport_mut <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Cell Model Passports/mutations_all_20211208.csv",
                              header = TRUE, check.names = FALSE)

cell_passport_cna_gistic <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Cell Model Passports/cnv_gistic_20191101.csv",
                                     header = TRUE, check.names = FALSE)
cell_passport_cna_picnic <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Cell Model Passports/cnv_abs_copy_number_picnic_20191101.csv",
                                     header = TRUE, check.names = FALSE)

cell_passport_fpkm <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Cell Model Passports/rnaseq_fpkm_20211124.csv",
                               header = TRUE, check.names = FALSE)
cell_passport_tpm <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Cell Model Passports/rnaseq_tpm_20211124.csv",
                              header = TRUE, check.names = FALSE)
cell_passport_counts <- read.csv("/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/Cell Model Passports/rnaseq_read_count_20211124.csv",
                                 header = TRUE, check.names = FALSE)


#' Adjust the knockout and expression data frames to "melt" them so that we can make boxplots from them,
#' as well as adding the mutation and CNA data
#' @param signif_hits_df a DepMap knockout or expression data frame, subsetted to only include genes that are 
#' significant hits from a particular run of the LM
#' @param mutation_df a DepMap mutation data frame
#' @param cna_df a DepMap CNA data frame
#' @param genes_of_interest a vector of the Hugo IDs of gene whose mutation/ CNA status is of interest
#' @param del_or_amp_vect a vector of "deletion" and "amplification" strings corresponding to which
#' we are interested for each given gene (should be same length as genes_of_interest)
adjust_depmap_df <- function(signif_hits_df, mutation_df, cna_df, genes_of_interest, del_or_amp_vect) {
  # Melt data frame to get it ready for boxplots
  signif_hits_df_melt <- melt(signif_hits_df)
  
  # Adjust column names
  colnames(signif_hits_df_melt)[which(colnames(signif_hits_df_melt) == "variable")] <- "gene"
  
  # Make the "value" column numeric, the "gene" column a character
  signif_hits_df_melt$value <- as.numeric(unlist(signif_hits_df_melt$value))
  signif_hits_df_melt$gene <- as.character(unlist(signif_hits_df_melt$gene))
  
  # Add mutation and CNA status of genes of interest to the data frame
  signif_hits_df_melt <- merge(signif_hits_df_melt, mutation_df[,c("depmap_id", genes_of_interest)], 
                               by = "depmap_id")
  genes_of_interest_cna <- unlist(lapply(genes_of_interest, function(x) paste0(x, ".CNA")))
  signif_hits_df_melt <- merge(signif_hits_df_melt, cna_df[,c("depmap_id", genes_of_interest_cna)], 
                               by = "depmap_id")
  
  # Make the mutation columns factors
  mutation_cols <- which(colnames(signif_hits_df_melt) %in% genes_of_interest)
  signif_hits_df_melt[, mutation_cols] <- lapply(signif_hits_df_melt[, mutation_cols], as.factor)
  
  
  # Bucket the CNA columns and also make them factors
  cna_cols <- which(colnames(signif_hits_df_melt) %in% genes_of_interest_cna)
  signif_hits_df_melt[, cna_cols] <- lapply(1:length(cna_cols), function(i) {
    col <- signif_hits_df_melt[, cna_cols[i]]
    bucketed_col <- bucket_cna(col, del_or_amp_vect[i])
    bucketed_col <- as.factor(bucketed_col)
    return(bucketed_col)
  })
  
  return(signif_hits_df_melt)
}

#' Bucket CNA values, given whether we want to prioritize deletions or amplifications
#' @param cna_vals vector of CNA values 
#' @param deletionOrAmp either "deletion" or "amplification" to indicate which we are prioritizing
bucket_cna <- function(cna_vals, deletionOrAmp) {
  bucketed_cna_vals <- unlist(lapply(cna_vals, function(x) {
    if(deletionOrAmp == "deletion") {
      # we use 3 here because there is a pseudocount of 1
      if (x < log2(3)) {return(1)}
      else {return(0)}
    } else {
      if (x > log2(3)) {return(1)}
      else {return(0)}
    }
  }))
  return(bucketed_cna_vals)
}


signif_hits_crispr_data <- adjust_depmap_df(signif_hits_crispr_data, tp53_pik3ca_mutations, tp53_pik3ca_cnas,
                                            c("TP53", "PIK3CA"), c("deletion", "amplification"))
signif_hits_rnai_data <- adjust_depmap_df(signif_hits_rnai_data, tp53_pik3ca_mutations, tp53_pik3ca_cnas,
                                            c("TP53", "PIK3CA"), c("deletion", "amplification"))
signif_hits_expression_data <- adjust_depmap_df(signif_hits_expression_data, tp53_pik3ca_mutations, tp53_pik3ca_cnas,
                                            c("TP53", "PIK3CA"), c("deletion", "amplification"))



# Subset to just the top hits for TP53 & PIK3CA specifically
tp53_top_hits <- c("MTHFD1L", "SHMT2", "GSTP1", "TYMS", "PSPH", "CDO1", "MTHFD1", "CSAD", "AHCY", 
                   "MTHFR", "MTR", "FTCD")
# additional hits from just TP53/ PIK3CA run: "SARDH", "DLD", "MAT2A", "GLDC"
tp53_top_hits_specific <- c(tp53_top_hits, c("SARDH", "DLD", "MAT2A", "GLDC"))

pik3ca_top_hits <- c("TYMS", "MTHFD1L", "PSPH", "GLDC")
# additional hits from just TP53/ PIK3CA run: "SHMT2", "FTCD", "MAT1A", "MTHFR", "AHCY", "CHKB", "MAT2A", "BHMT2", "GCLM"
pik3ca_top_hits_specific <- c(pik3ca_top_hits, c("SHMT2", "FTCD", "MAT1A", "MTHFR", "AHCY", "CHKB", "MAT2A", "BHMT2", "GCLM"))


#' Given a list of top hits, subset the given DF to this list of hits
#' @param signif_hits_df a DepMap knockout or expression data frame, subsetted to only include genes that are 
#' significant hits from a particular run of the LM
#' @param sig_hits_goi a vector of significant hits for only a particular gene of interest
subset_to_gene_tophits <- function(signif_hits_df, sig_hits_goi) {
  
  signif_hits_df <- signif_hits_df[signif_hits_df$gene %fin% sig_hits_goi,]
  
  return(signif_hits_df)
}


signif_hits_crispr_data_tp53 <- subset_to_gene_tophits(signif_hits_crispr_data, tp53_top_hits_specific)
signif_hits_rnai_data_tp53 <- subset_to_gene_tophits(signif_hits_rnai_data, tp53_top_hits_specific)
signif_hits_expression_data_tp53 <- subset_to_gene_tophits(signif_hits_expression_data, tp53_top_hits_specific)

signif_hits_crispr_data_pik3ca <- subset_to_gene_tophits(signif_hits_crispr_data, pik3ca_top_hits_specific)
signif_hits_rnai_data_pik3ca <- subset_to_gene_tophits(signif_hits_rnai_data, pik3ca_top_hits_specific)
signif_hits_expression_data_pik3ca <- subset_to_gene_tophits(signif_hits_expression_data, pik3ca_top_hits_specific)


#' Run a significance test for each combo, given an input DF and either
#' "mut" or "cna
#' @param df input DepMap data frame, melted & subsetted
#' @param goi gene of interest, Hugo ID
#' @param mut_or_cna either "mut", to signify we want to test dependency by mutation status, 
#' or "cna" to indicate we want to test dependency by CNA status
run_ttest <- function(df, goi, mut_or_cna) {
  stat_test <- NA
  
  if(mut_or_cna == "mut") {
    stat.test <- df %>%
      group_by(gene) %>%
      #wilcox.test(value ~ TP53) %>%
      t_test(reformulate(goi, response = "value", intercept = FALSE)) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")
  } else {
    goi_cna <- paste0(goi, ".CNA")
    stat.test <- df %>%
      group_by(gene) %>%
      #wilcox.test(value ~ TP53) %>%
      t_test(reformulate(goi_cna, response = "value", intercept = FALSE)) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")
  }
  
  return(stat.test)
}

ttest_crispr_tp53_mut <- run_ttest(signif_hits_crispr_data_tp53, "TP53", "mut")
ttest_rnai_tp53_mut <- run_ttest(signif_hits_rnai_data_tp53, "TP53", "mut")
ttest_expression_tp53_mut <- run_ttest(signif_hits_expression_data_tp53, "TP53", "mut")

ttest_crispr_tp53_cna <- run_ttest(signif_hits_crispr_data_tp53, "TP53", "cna")
ttest_rnai_tp53_cna <- run_ttest(signif_hits_rnai_data_tp53, "TP53", "cna")
ttest_expression_tp53_cna <- run_ttest(signif_hits_expression_data_tp53, "TP53", "cna")

ttest_crispr_pik3ca_mut <- run_ttest(signif_hits_crispr_data_pik3ca, "PIK3CA", "mut")
ttest_rnai_pik3ca_mut <- run_ttest(signif_hits_rnai_data_pik3ca, "PIK3CA", "mut")
ttest_expression_pik3ca_mut <- run_ttest(signif_hits_expression_data_pik3ca, "PIK3CA", "mut")

ttest_crispr_pik3ca_cna <- run_ttest(signif_hits_crispr_data_pik3ca, "PIK3CA", "cna")
ttest_rnai_pik3ca_cna <- run_ttest(signif_hits_rnai_data_pik3ca, "PIK3CA", "cna")
ttest_expression_pik3ca_cna <- run_ttest(signif_hits_expression_data_pik3ca, "PIK3CA", "cna")


#' Plot boxplot of each gene's dependency across BRCA cell lines, by gene of interest mutation status
#' @param signif_hits_df a DepMap knockout or expression data frame, subsetted to only include genes that are 
#' significant hits from a particular run of the LM
#' @param gene_of_interest a Hugo ID for a given gene of interest
#' @param mut_or_cna either "mut", to signify we want to plot dependency by mutation status, or "cna" to
#' indicate we want to plot dependency by CNA status
#' @param ttest_res optional: can include the results of a t-test with the significance values on plot
make_dependency_boxplot <- function(signif_hits_df, gene_of_interest, mut_or_cna, ttest_res) {
  if(mut_or_cna == "mut") {
    bxp <- ggplot(signif_hits_df, aes_string(x = "gene", y = "value", fill = gene_of_interest)) + geom_boxplot() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
  } else if (mut_or_cna == "cna") {
    gene_of_interest_cna <- paste0(gene_of_interest, ".CNA")
    bxp <- ggplot(signif_hits_df, aes_string(x = "gene", y = "value", fill = gene_of_interest_cna)) + geom_boxplot() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
  } else {
    print(paste(mut_or_cna, "is not defined. Please try again with either 'mut' or 'cna'."))
  }
  print(bxp)

  # If ttest_res is defined, add the p-values onto the plot
  try({
    ttest_res <- ttest_res %>% add_xy_position(x = "gene", dodge = 0.8)
    print(ttest_res)
    bxp + stat_pvalue_manual(ttest_res,  label = "p.adj.signif", tip.length = 0) + 
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
    
  }, silent = TRUE)
}

make_dependency_boxplot(signif_hits_crispr_data_tp53, "TP53", "mut", ttest_crispr_tp53_mut)
make_dependency_boxplot(signif_hits_rnai_data_tp53, "TP53", "mut", ttest_rnai_tp53_mut)
make_dependency_boxplot(signif_hits_expression_data_tp53, "TP53", "mut", ttest_expression_tp53_mut)

make_dependency_boxplot(signif_hits_crispr_data_tp53, "TP53", "cna", ttest_crispr_tp53_cna)
make_dependency_boxplot(signif_hits_rnai_data_tp53, "TP53", "cna", ttest_rnai_tp53_cna)
make_dependency_boxplot(signif_hits_expression_data_tp53, "TP53", "cna", ttest_expression_tp53_cna)

make_dependency_boxplot(signif_hits_crispr_data_pik3ca, "PIK3CA", "mut", ttest_crispr_pik3ca_mut)
make_dependency_boxplot(signif_hits_rnai_data_pik3ca, "PIK3CA", "mut", ttest_rnai_pik3ca_mut)
make_dependency_boxplot(signif_hits_expression_data_pik3ca, "PIK3CA", "mut", ttest_expression_pik3ca_mut)

make_dependency_boxplot(signif_hits_crispr_data_pik3ca, "PIK3CA", "cna", ttest_crispr_pik3ca_cna)
make_dependency_boxplot(signif_hits_rnai_data_pik3ca, "PIK3CA", "cna", ttest_rnai_pik3ca_cna)
make_dependency_boxplot(signif_hits_expression_data_pik3ca, "PIK3CA", "cna", ttest_expression_pik3ca_cna)


# Split this by predicted upregulation/ predicted downregulation, so we can see how closely DepMap aligns
pred_upregulat_tp53 <- c("AHCY", "DLD", "FTCD", "GLDC", "GSTP1", "MTHFD1", "MTHFD1L", "PSPH", "SHMT2", "TYMS")
pred_downregulat_tp53 <- c("CDO1", "CSAD", "MAT2A", "MTHFR", "MTR", "SARDH")

pred_upregulat_pik3ca <- c("BHMT2", "MTHFR")
pred_downregulat_pik3ca <- c("AHCY", "CHKB", "GLDC", "FTCD", "MAT1A", "MTHFDL1", "PSPH", "SHMT2", "TYMS")

make_dependency_boxplot(signif_hits_expression_data_tp53[signif_hits_expression_data_tp53$gene %in% pred_upregulat_tp53,], "TP53", "mut")
make_dependency_boxplot(signif_hits_expression_data_tp53[signif_hits_expression_data_tp53$gene %in% pred_downregulat_tp53,], "TP53", "mut")
make_dependency_boxplot(signif_hits_expression_data_pik3ca[signif_hits_expression_data_pik3ca$gene %in% pred_upregulat_pik3ca,], "PIK3CA", "mut")
make_dependency_boxplot(signif_hits_expression_data_pik3ca[signif_hits_expression_data_pik3ca$gene %in% pred_downregulat_pik3ca,], "PIK3CA", "mut")

make_dependency_boxplot(signif_hits_expression_data_tp53[signif_hits_expression_data_tp53$gene %in% pred_upregulat_tp53,], "TP53", "cna")
make_dependency_boxplot(signif_hits_expression_data_tp53[signif_hits_expression_data_tp53$gene %in% pred_downregulat_tp53,], "TP53", "cna")
make_dependency_boxplot(signif_hits_expression_data_pik3ca[signif_hits_expression_data_pik3ca$gene %in% pred_upregulat_pik3ca,], "PIK3CA", "cna")
make_dependency_boxplot(signif_hits_expression_data_pik3ca[signif_hits_expression_data_pik3ca$gene %in% pred_downregulat_pik3ca,], "PIK3CA", "cna")


make_dependency_boxplot(signif_hits_crispr_data_tp53[signif_hits_crispr_data_tp53$gene %in% pred_upregulat_tp53,], "TP53", "mut")
make_dependency_boxplot(signif_hits_crispr_data_tp53[signif_hits_crispr_data_tp53$gene %in% pred_downregulat_tp53,], "TP53", "mut")
make_dependency_boxplot(signif_hits_crispr_data_pik3ca[signif_hits_crispr_data_pik3ca$gene %in% pred_upregulat_pik3ca,], "PIK3CA", "mut")
make_dependency_boxplot(signif_hits_crispr_data_pik3ca[signif_hits_crispr_data_pik3ca$gene %in% pred_downregulat_pik3ca,], "PIK3CA", "mut")

make_dependency_boxplot(signif_hits_crispr_data_tp53[signif_hits_crispr_data_tp53$gene %in% pred_upregulat_tp53,], "TP53", "cna")
make_dependency_boxplot(signif_hits_crispr_data_tp53[signif_hits_crispr_data_tp53$gene %in% pred_downregulat_tp53,], "TP53", "cna")
make_dependency_boxplot(signif_hits_crispr_data_pik3ca[signif_hits_crispr_data_pik3ca$gene %in% pred_upregulat_pik3ca,], "PIK3CA", "cna")
make_dependency_boxplot(signif_hits_crispr_data_pik3ca[signif_hits_crispr_data_pik3ca$gene %in% pred_downregulat_pik3ca,], "PIK3CA", "cna")


make_dependency_boxplot(signif_hits_rnai_data_tp53[signif_hits_rnai_data_tp53$gene %in% pred_upregulat_tp53,], "TP53", "mut")
make_dependency_boxplot(signif_hits_rnai_data_tp53[signif_hits_rnai_data_tp53$gene %in% pred_downregulat_tp53,], "TP53", "mut")
make_dependency_boxplot(signif_hits_rnai_data_pik3ca[signif_hits_rnai_data_pik3ca$gene %in% pred_upregulat_pik3ca,], "PIK3CA", "mut")
make_dependency_boxplot(signif_hits_rnai_data_pik3ca[signif_hits_rnai_data_pik3ca$gene %in% pred_downregulat_pik3ca,], "PIK3CA", "mut")

make_dependency_boxplot(signif_hits_rnai_data_tp53[signif_hits_rnai_data_tp53$gene %in% pred_upregulat_tp53,], "TP53", "cna")
make_dependency_boxplot(signif_hits_rnai_data_tp53[signif_hits_rnai_data_tp53$gene %in% pred_downregulat_tp53,], "TP53", "cna")
make_dependency_boxplot(signif_hits_rnai_data_pik3ca[signif_hits_rnai_data_pik3ca$gene %in% pred_upregulat_pik3ca,], "PIK3CA", "cna")
make_dependency_boxplot(signif_hits_rnai_data_pik3ca[signif_hits_rnai_data_pik3ca$gene %in% pred_downregulat_pik3ca,], "PIK3CA", "cna")


# Split cell lines by whether they are derived from primary or metastatic, to see if there is a difference
cell_line_sample_info <- read.csv(paste0("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/cell_line_sample_info.csv"),
                                  header = TRUE, check.names = FALSE)
brca_cl_info <- cell_line_sample_info[grepl("Breast", cell_line_sample_info$primary_disease),]

# Get the cell lines we have mutation/ expression data for
cls_with_dat <- unique(signif_hits_crispr_data$depmap_id)

# Subset the sample info DF to only the cell lines of interest
brca_cl_info_sub <- brca_cl_info[brca_cl_info$DepMap_ID %fin% cls_with_dat,]

# Add a "primary" or "metastatic" label to the signif_hits data frames

#' Adds a label for each cell line denoting whether the cell line is derived from a primary or metastatic site
#' @param sample_info a data frame correlating cell line ID to site
#' @param signif_hits_df a data frame to add the site information to
add_site <- function(sample_info, signif_hits_df) {
  sites_vect <- unlist(lapply(signif_hits_df$depmap_id, function(id) {
    lab <- sample_info[sample_info$DepMap_ID == id, 'primary_or_metastasis']
    if(length(lab) == 0) {lab <- NA}
    return(lab)
  }))
  signif_hits_df$site <- sites_vect
  return(signif_hits_df)
}

signif_hits_expression_data_tp53 <- add_site(brca_cl_info_sub, signif_hits_expression_data_tp53)
signif_hits_expression_data_pik3ca <- add_site(brca_cl_info_sub, signif_hits_expression_data_pik3ca)
signif_hits_crispr_data_tp53 <- add_site(brca_cl_info_sub, signif_hits_crispr_data_tp53)
signif_hits_crispr_data_pik3ca <- add_site(brca_cl_info_sub, signif_hits_crispr_data_pik3ca)
signif_hits_rnai_data_tp53 <- add_site(brca_cl_info_sub, signif_hits_rnai_data_tp53)
signif_hits_rnai_data_pik3ca <- add_site(brca_cl_info_sub, signif_hits_rnai_data_pik3ca)


# Segregate them into separate primary and metastatic files, leaving out the NA CLs
signif_hits_expression_data_primary_tp53 <- signif_hits_expression_data_tp53[signif_hits_expression_data_tp53$site == "Primary",]
signif_hits_expression_data_metastatic_tp53 <- signif_hits_expression_data_tp53[signif_hits_expression_data_tp53$site == "Metastasis",]

signif_hits_expression_data_primary_pik3ca <- signif_hits_expression_data_pik3ca[signif_hits_expression_data_pik3ca$site == "Primary",]
signif_hits_expression_data_metastatic_pik3ca <- signif_hits_expression_data_pik3ca[signif_hits_expression_data_pik3ca$site == "Metastasis",]


make_dependency_boxplot(signif_hits_expression_data_primary_tp53[signif_hits_expression_data_primary_tp53$gene %in% pred_upregulat_tp53,], "TP53", "mut")
make_dependency_boxplot(signif_hits_expression_data_metastatic_tp53[signif_hits_expression_data_metastatic_tp53$gene %in% pred_downregulat_tp53,], "TP53", "mut")
make_dependency_boxplot(signif_hits_expression_data_pik3ca[signif_hits_expression_data_pik3ca$gene %in% pred_upregulat_pik3ca,], "PIK3CA", "mut")
make_dependency_boxplot(signif_hits_expression_data_pik3ca[signif_hits_expression_data_pik3ca$gene %in% pred_downregulat_pik3ca,], "PIK3CA", "mut")



# Code to test expression of a few cell lines of interest
cl_00 <- "ACH-000349"
cl_01 <- "ACH-001388"
cl_10 <- "ACH-000288"
ids_of_interest <- c(cl_00, cl_01, cl_10)

pred_up_tp53_down_pik3ca <- intersect(pred_upregulat_tp53, pred_downregulat_pik3ca)
pred_down_tp53_up_pik3ca <- intersect(pred_downregulat_tp53, pred_upregulat_pik3ca)
pred_up_tp53_noEff_pik3ca <- unique(c(setdiff(pred_upregulat_tp53, c(pred_upregulat_pik3ca, pred_downregulat_pik3ca))))

#' Make comparative barplot of the expression of 3 cell lines in various mutational categories
#' @param signif_hits_df a DepMap expression data frame, subsetted to only include genes that are 
#' significant hits from a particular run of the LM
make_expression_barplot <- function(signif_hits_df) {
  ggplot(signif_hits_df, aes_string(x = "gene", y = "value", fill = "depmap_id")) + geom_bar(position = "dodge", stat= "identity") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}

make_expression_barplot(signif_hits_expression_data[(signif_hits_expression_data$gene %in% pred_up_tp53_down_pik3ca) &
                                                                   (signif_hits_expression_data$depmap_id %in% ids_of_interest),])
make_expression_barplot(signif_hits_expression_data[(signif_hits_expression_data$gene %in% pred_down_tp53_up_pik3ca) &
                                                      (signif_hits_expression_data$depmap_id %in% ids_of_interest),])
make_expression_barplot(signif_hits_expression_data[(signif_hits_expression_data$gene %in% pred_up_tp53_noEff_pik3ca) &
                                                      (signif_hits_expression_data$depmap_id %in% ids_of_interest),])

#' For the given target genes, does a Wilcoxon test to see if there
#' is significant DE between the GOI mutated and non-mutated DepMap groups.
#' Prints the target gene and the p-value of the Wilcoxon test, along with 
#' the directionality of the effect
#' @param expression_df an expression data frame with mutational status appended
#' @param goi the Hugo Symbol of a gene of interest
wilcox_of_targs <- function(expression_df, goi) {
  
  targs <- unique(expression_df$gene)
  print(targs)
  pvals <- c()
  
  for (t in targs) {
    print(paste("Target:", t))
    
    exp_mut <- as.numeric(expression_df[(expression_df$gene == t) & 
                                          (expression_df[, goi] == 1), 'value'])
    exp_noMut <- as.numeric(expression_df[(expression_df$gene == t) & 
                                            (expression_df[, goi] == 0), 'value'])
    
    wilcox_res <- wilcox.test(exp_mut, exp_noMut, exact = FALSE)   # use default two-sided alternative
    
    pval <- wilcox_res$p.value
    print(paste("Wilcoxon p-value:", pval))
    pvals <- c(pvals, pval)
    
    if(mean(exp_mut) > mean(exp_noMut)) {print("Upregulation")}
    else {print("Downregulation")}
  }
  return(list("Targets" = targs, "P.value" = pvals))
}

pvalues_depmap <- wilcox_of_targs(signif_hits_expression_data_primary_tp53, "TP53")

pvalues$Signif <- unlist(lapply(pvalues$P.value, function(x) ifelse(x < 0.05, 1, 0)))
pvalues$ChipEat <- unlist(lapply(pvalues$Targets, function(x) ifelse(x %in% tp53_chipeat_targs, 1, 0)))

jac_res <- clujaccard(pvalues$ChipEat, pvalues$Signif, zerobyzero = NA)

###################################################################################
# INVESTIGATE SLOPES USING A SIMPLE LM OF DEPMAP VAL ~ MUT/ AMP STATUS
###################################################################################
#' Run a simple linear model relating DepMap score to the mutation or CNA
#' status of TP53/ PIK3CA. Creates a volcano plot with Beta of fitted LM on the x-axis
#' and -log(pval) on the y-axis
#' @param depmap_table a depmap table with CRISPR or RNAi KO scores for given gene set 
#' @param mut_or_cna a string that reads either "mut" or "cna" to indicate which we are 
#' looking to relate to depmap 
#' @param gene_name the name of the gene to investigate for mutations/ CNAs (currently
#' implemented for either 'TP53' or 'PIK3CA')
depmap_lm <- function(depmap_table, mut_or_cna, gene_name) {
    input_df <- data.frame("Depmap" = depmap_table$value, "Gene" = depmap_table$gene)
    if(gene_name == 'TP53') {
        if (mut_or_cna == "mut") {
            input_df$status <- depmap_table$TP53
        } else {
            input_df$status <- depmap_table$TP53.CNA
        }
    } else if (gene_name == "PIK3CA") {
        if (mut_or_cna == "mut") {
            input_df$status <- depmap_table$PIK3CA
        } else {
            input_df$status <- depmap_table$PIK3CA.CNA
        }
    } else {
        print("error. implemented only for TP53 and PIK3CA.")
    }
    print(input_df)
    master_lm_res <- data.frame()
    
    for (i in 1:length(unique(input_df$Gene))) {
        gene <- unique(input_df$Gene)[i]
        input_df_gene <- input_df[input_df$Gene == gene,]
        lm_res <- lm(Depmap ~ status, data = input_df_gene)
        lm_res_summary <- summary(lm_res)
        lm_res_tidy <- as.data.frame(tidy(lm_res_summary))
        lm_res_tidy$Gene <- gene
        master_lm_res <- rbind(master_lm_res, lm_res_tidy)
    }
    print(master_lm_res)
    
    betas <- master_lm_res$estimate
    p.vals <- master_lm_res$p.value
    genes <- master_lm_res$Gene
    plt_input <- data.frame("Betas" = betas, "p.val" = p.vals, "Gene" = genes,
                            "Neg.Log_p.val" = unlist(lapply(p.vals, function(x) 
        return(-log(x, base = 10)))))
    # Add differential dependency
    plt_input$diffDep <- "NONE"
    plt_input$diffDep[plt_input$Betas > 0 & plt_input$p.val < 0.05] <- "MORE DEP."
    plt_input$diffDep[plt_input$Betas < 0 & plt_input$p.val < 0.05] <- "LESS DEP."

    plt_input$diffDepLabel <- NA
    plt_input$diffDepLabel[plt_input$diffDep != "NONE"] <- plt_input$Gene[plt_input$diffDep != "NONE"]
    
    print(plt_input)
    
    plt <- ggplot(plt_input, aes(x = Betas, y = Neg.Log_p.val, col = diffDep, label = diffDepLabel)) + geom_point() +
        geom_hline(yintercept = -log10(0.05), col="red") + theme_minimal() +
        geom_text_repel() + scale_color_manual(values=c("blue", "black", "red")) +
        ggtitle(paste(gene_name, mut_or_cna))
    print(plt)
}

depmap_lm(signif_hits_crispr_data_tp53, "mut", "TP53")
depmap_lm(signif_hits_crispr_data_tp53, "cna", "TP53")
depmap_lm(signif_hits_crispr_data_pik3ca, "mut", "PIK3CA")
depmap_lm(signif_hits_crispr_data_pik3ca, "cna", "PIK3CA")

depmap_lm(signif_hits_crispr_data, "mut", "TP53")
depmap_lm(signif_hits_crispr_data, "cna", "TP53")
depmap_lm(signif_hits_crispr_data, "mut", "PIK3CA")
depmap_lm(signif_hits_crispr_data, "cna", "PIK3CA")

