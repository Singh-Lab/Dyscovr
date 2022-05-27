# Interrogate CPTAC Data from cBioPortal (Script)

library(ggplot2)
library(reshape2)
library(ggpubr)
library(rstatix)
library(broom)
library(ggrepel)
library(rstatix)
library(fpc)

cptac_brca_mutation_data <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/cBioPortal/brca_cptac_2020/data_mutations.txt", 
                                     sep = "\t", header = TRUE, check.names = FALSE)
cptac_brca_missense_nonsense_mutation_data <- cptac_brca_mutation_data[(cptac_brca_mutation_data$Variant_Classification == "Missense_Mutation") |
                                                                         (cptac_brca_mutation_data$Variant_Classification == "Nonsense_Mutation"),]

cptac_brca_fpkm_data <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/cBioPortal/brca_cptac_2020/data_mrna_seq_fpkm.txt", 
                                 sep = "\t", header = TRUE, check.names = FALSE)
cptac_brca_fpkm_data_complete <- cptac_brca_fpkm_data[rowSums(is.na(cptac_brca_fpkm_data[,2:ncol(cptac_brca_fpkm_data)])) == 0,]

cptac_brca_proteomic_data <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/cBioPortal/brca_cptac_2020/data_protein_quantification.txt", 
                                      sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
rownames(cptac_brca_proteomic_data) <- unlist(lapply(rownames(cptac_brca_proteomic_data), function(x)
  unlist(strsplit(x, "|", fixed = TRUE))[1]))

cptac_brca_tp53_mutation <- cptac_brca_missense_nonsense_mutation_data[cptac_brca_missense_nonsense_mutation_data$Hugo_Symbol == "TP53",]
tp53_mut_patients <- unique(cptac_brca_tp53_mutation$Tumor_Sample_Barcode)
length(tp53_mut_patients)
tp53_not_mut_patients <- setdiff(unique(cptac_brca_missense_nonsense_mutation_data$Tumor_Sample_Barcode), tp53_mut_patients)

cptac_brca_pik3ca_mutation <- cptac_brca_missense_nonsense_mutation_data[cptac_brca_missense_nonsense_mutation_data$Hugo_Symbol == "PIK3CA",]
pik3ca_mut_patients <- unique(cptac_brca_pik3ca_mutation$Tumor_Sample_Barcode)
length(pik3ca_mut_patients)
pik3ca_not_mut_patients <- setdiff(unique(cptac_brca_missense_nonsense_mutation_data$Tumor_Sample_Barcode), pik3ca_mut_patients)

cptac_brca_clinical_sample <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/cBioPortal/brca_cptac_2020/data_clinical_sample.txt", 
                                       sep = "\t", header = TRUE, check.names = FALSE, row.names = 1, comment.char = "#")
cptac_lumA_samples <- cptac_brca_clinical_sample[cptac_brca_clinical_sample$PAM50 == "LumA", 'SAMPLE_ID']
cptac_lumB_samples <- cptac_brca_clinical_sample[cptac_brca_clinical_sample$PAM50 == "LumB", 'SAMPLE_ID']
cptac_basal_samples <- cptac_brca_clinical_sample[cptac_brca_clinical_sample$PAM50 == "Basal", 'SAMPLE_ID']
cptac_her2_samples <- cptac_brca_clinical_sample[cptac_brca_clinical_sample$PAM50 == "Her2", 'SAMPLE_ID']
cptac_normLike_samples <- cptac_brca_clinical_sample[cptac_brca_clinical_sample$PAM50 == "Normal-like", 'SAMPLE_ID']


# MTHFD1L + TP53
ensg <- "ENSG00000120254"
symbol <- "MTHFD1L"

expression_mutants <- as.numeric(cptac_brca_fpkm_data_complete[cptac_brca_fpkm_data_complete$Hugo_Symbol == symbol, 
                                                               colnames(cptac_brca_fpkm_data_complete) %in% tp53_mut_patients])
expression_normal <- as.numeric(cptac_brca_fpkm_data_complete[cptac_brca_fpkm_data_complete$Hugo_Symbol == symbol, 
                                                              !colnames(cptac_brca_fpkm_data_complete) %in% tp53_mut_patients])
tp53_expression_by_mut_group <- list("Mutation" = expression_mutants, "No Mutation" = expression_normal)
boxplot(tp53_expression_by_mut_group, ylab = "Expression (FPKM)", 
        main = "MTHFD1L Expression By TP53 Mutation Status")

wilcox.test(expression_mutants, y = expression_normal)

# TYMS + PIK3CA
ensg <- "ENSG00000176890"
symbol <- "TYMS"

expression_mutants <- as.numeric(unlist(cptac_brca_fpkm_data[cptac_brca_fpkm_data$Hugo_Symbol == symbol, 
                                                             colnames(cptac_brca_fpkm_data) %in% pik3ca_mut_patients]))
expression_normal <- as.numeric(unlist(cptac_brca_fpkm_data[cptac_brca_fpkm_data$Hugo_Symbol == symbol, 
                                                            !colnames(cptac_brca_fpkm_data) %in% pik3ca_mut_patients]))
pik3ca_expression_by_mut_group <- list("Mutation" = expression_mutants, "No Mutation" = expression_normal)
boxplot(pik3ca_expression_by_mut_group, ylab = "Expression (FPKM)", 
        main = "TYMS Expression By PIK3CA Mutation Status")

wilcox.test(expression_mutants, y = expression_normal)

# TP53 + others
symbol <- "SHMT2"
symbol <- "GSTP1"

# Do this on a larger scale across all genes of interest

# Subset to just the top hits for TP53 & PIK3CA specifically
tp53_top_hits <- c("MTHFD1L", "SHMT2", "GSTP1", "TYMS", "PSPH", "CDO1", "MTHFD1", "CSAD", "AHCY", 
                   "MTHFR", "MTR", "FTCD")
# additional hits from just TP53/ PIK3CA run: "SARDH", "DLD", "MAT2A", "GLDC"
tp53_top_hits_specific <- c(tp53_top_hits, c("SARDH", "DLD", "MAT2A", "GLDC"))

pik3ca_top_hits <- c("TYMS", "MTHFD1L", "PSPH", "GLDC")
# additional hits from just TP53/ PIK3CA run: "SHMT2", "FTCD", "MAT1A", "MTHFR", "AHCY", "CHKB", "MAT2A", "BHMT2", "GCLM"
pik3ca_top_hits_specific <- c(pik3ca_top_hits, c("SHMT2", "FTCD", "MAT1A", "MTHFR", "AHCY", "CHKB", "MAT2A", "BHMT2", "GCLM"))


#' Given a list of top hits, subset the given DF to this list of hits
#' @param signif_hits_df a CPTAC expression data frame
#' @param sig_hits_goi a vector of significant hits for only a particular gene of interest
subset_to_gene_tophits <- function(signif_hits_df, sig_hits_goi) {
  
  signif_hits_df <- signif_hits_df[signif_hits_df$Hugo_Symbol %fin% sig_hits_goi,]
  
  return(signif_hits_df)
}

cptac_brca_fpkm_data_complete_tp53 <- subset_to_gene_tophits(cptac_brca_fpkm_data_complete, tp53_top_hits_specific)
cptac_brca_fpkm_data_complete_pik3ca <- subset_to_gene_tophits(cptac_brca_fpkm_data_complete, pik3ca_top_hits_specific)


#' Get the patients that have a mutated version of each gene of interest and add to the expression DF
#' @param expression_df a CPTAC expression DF 
#' @param mutation_df a CPTAC mutation DF
#' @param goi the Hugo Symbol of a gene of interest
fix_expression_df <- function(expression_df, mutation_df, goi) {
  
  exp_df_m <- melt(expression_df)
  colnames(exp_df_m) <- c("Hugo_Symbol", "Patient_ID", "FPKM")
  
  print(exp_df_m)
  
  # Add mutation status of the GOI
  status <- unlist(lapply(exp_df_m$Patient_ID, function(id) {
    if(id %in% mutation_df$Tumor_Sample_Barcode) {return(1)}
    else {return(0)}
  }))
  exp_df_m[,4] <- status
  lab <- paste0(goi, "_Mut")
  colnames(exp_df_m)[4] <- lab
  
  return(exp_df_m)
}

cptac_brca_fpkm_data_complete_tp53 <- fix_expression_df(cptac_brca_fpkm_data_complete_tp53, 
                                                        cptac_brca_tp53_mutation, "TP53")
cptac_brca_fpkm_data_complete_pik3ca <- fix_expression_df(cptac_brca_fpkm_data_complete_pik3ca, 
                                                          cptac_brca_pik3ca_mutation, "PIK3CA")

cptac_brca_fpkm_data_complete_tp53$TP53_Mut <- as.factor(cptac_brca_fpkm_data_complete_tp53$TP53_Mut)
cptac_brca_fpkm_data_complete_pik3ca$PIK3CA_Mut <- as.factor(cptac_brca_fpkm_data_complete_pik3ca$PIK3CA_Mut)


# Plot boxplot of each gene's expression across CPTAC samples, by gene of interest mutation status
ggplot(cptac_brca_fpkm_data_complete_tp53, aes_string(x = "Hugo_Symbol", y = "FPKM", fill = "TP53_Mut")) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(cptac_brca_fpkm_data_complete_pik3ca, aes_string(x = "Hugo_Symbol", y = "FPKM", fill = "PIK3CA_Mut")) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Split this by predicted upregulation/ predicted downregulation, so we can see how closely DepMap aligns
pred_upregulat_tp53 <- c("AHCY", "DLD", "FTCD", "GLDC", "GSTP1", "MTHFD1", "MTHFD1L", "PSPH", "SHMT2", "TYMS")
pred_downregulat_tp53 <- c("CDO1", "CSAD", "MAT2A", "MTHFR", "MTR", "SARDH")

pred_upregulat_pik3ca <- c("BHMT2", "MTHFR")
pred_downregulat_pik3ca <- c("AHCY", "CHKB", "GLDC", "FTCD", "MAT1A", "MTHFDL1", "PSPH", "SHMT2", "TYMS")



ggplot(cptac_brca_fpkm_data_complete_tp53[cptac_brca_fpkm_data_complete_tp53$Hugo_Symbol %in% pred_upregulat_tp53,], 
       aes_string(x = "Hugo_Symbol", y = "FPKM", fill = "TP53_Mut")) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(cptac_brca_fpkm_data_complete_tp53[cptac_brca_fpkm_data_complete_tp53$Hugo_Symbol %in% pred_downregulat_tp53,], 
       aes_string(x = "Hugo_Symbol", y = "FPKM", fill = "TP53_Mut")) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(cptac_brca_fpkm_data_complete_pik3ca[cptac_brca_fpkm_data_complete_pik3ca$Hugo_Symbol %in% pred_upregulat_pik3ca,], 
       aes_string(x = "Hugo_Symbol", y = "FPKM", fill = "PIK3CA_Mut")) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(cptac_brca_fpkm_data_complete_pik3ca[cptac_brca_fpkm_data_complete_pik3ca$Hugo_Symbol %in% pred_downregulat_pik3ca,], 
       aes_string(x = "Hugo_Symbol", y = "FPKM", fill = "PIK3CA_Mut")) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#' Get the patients that have a mutated version of each gene of interest and add to the proteomic DF
#' @param proteomic_df a CPTAC proteomic DF 
#' @param mutation_df a CPTAC mutation DF
#' @param goi the Hugo Symbol of a gene of interest
fix_proteomic_df <- function(proteomic_df, mutation_df, goi) {
  
  proteomic_df$gene.name <- rownames(proteomic_df)
  proteomic_df_m <- melt(proteomic_df)
  colnames(proteomic_df_m) <- c("Hugo_Symbol", "Patient_ID", "Prot.Quant")
  
  print(head(proteomic_df_m))
  
  # Add mutation status of the GOI
  status <- unlist(lapply(proteomic_df_m$Patient_ID, function(id) {
    if(id %in% mutation_df$Tumor_Sample_Barcode) {return(1)}
    else {return(0)}
  }))
  proteomic_df_m[,4] <- status
  lab <- paste0(goi, "_Mut")
  colnames(proteomic_df_m)[4] <- lab
  
  proteomic_df_m <- na.omit(proteomic_df_m)
  
  return(proteomic_df_m)
}

cptac_brca_proteomic_data_m <- fix_proteomic_df(cptac_brca_proteomic_data,  cptac_brca_tp53_mutation, "TP53")
cptac_brca_proteomic_data_m <- fix_proteomic_df(cptac_brca_proteomic_data,  cptac_brca_pik3ca_mutation, "PIK3CA")


#' For the given target genes, does a Wilcoxon test to see if there
#' is significant DE between the GOI mutated and non-mutated CPTAC-3 groups.
#' Prints the target gene and the p-value of the Wilcoxon test, along with 
#' the directionality of the effect
#' @param cptac_df a "complete" FPKM data frame with mutational status appended
#' @param goi the Hugo Symbol of a gene of interest
#' @param label either "FPKM" or "Prot.Quant" to indicate what kind of expression
#' we are looking at
wilcox_of_targs <- function(cptac_df, goi, label) {
  
  targs <- unique(cptac_df$Hugo_Symbol)
  pvals <- c()
  
  for (t in targs) {
    print(paste("Target:", t))
    
    colname <- paste(goi, "Mut", sep = "_")
    exp_mut <- as.numeric(cptac_df[(cptac_df$Hugo_Symbol == t) & (cptac_df[,colname] == 1), label])
    exp_noMut <- as.numeric(cptac_df[(cptac_df$Hugo_Symbol == t) & (cptac_df[,colname] == 0), label])
    
    wilcox_res <- wilcox.test(exp_mut, exp_noMut)   # use default two-sided alternative
    
    pval <- wilcox_res$p.value
    print(paste("Wilcoxon p-value:", pval))
    pvals <- c(pvals, pval)
    
    if(mean(exp_mut) > mean(exp_noMut)) {print("Upregulation")}
    else {print("Downregulation")}
  }
  return(list("Targets" = targs, "P.value" = pvals))
}

pvalues <- wilcox_of_targs(cptac_brca_fpkm_data_complete_tp53, "TP53", "FPKM")
pvalues <- wilcox_of_targs(cptac_brca_proteomic_data_m, "TP53", "Prot.Quant")

pvalues <- wilcox_of_targs(cptac_brca_fpkm_data_complete_pik3ca, "PIK3CA", "FPKM")
pvalues <- wilcox_of_targs(cptac_brca_proteomic_data_m, "PIK3CA", "Prot.Quant")


pvalues$Signif <- unlist(lapply(pvalues$P.value, function(x) ifelse(x < 0.05, 1, 0)))
pvalues$ChipEat <- unlist(lapply(pvalues$Targets, function(x) ifelse(x %in% tp53_chipeat_targs, 1, 0)))

pvalues <- as.data.frame(pvalues)
pvalues <- pvalues[order(pvalues$P.value),]

jac_res <- clujaccard(pvalues$ChipEat, pvalues$Signif, zerobyzero = NA)


# Do this for individual subtypes 
cptac_brca_fpkm_data_complete_tp53_lumA <- cptac_brca_fpkm_data_complete_tp53[cptac_brca_fpkm_data_complete_tp53$Patient_ID %in% 
                                                                                cptac_lumA_samples,]
cptac_brca_fpkm_data_complete_tp53_lumB <- cptac_brca_fpkm_data_complete_tp53[cptac_brca_fpkm_data_complete_tp53$Patient_ID %in% 
                                                                                cptac_lumB_samples,]
cptac_brca_fpkm_data_complete_tp53_basal <- cptac_brca_fpkm_data_complete_tp53[cptac_brca_fpkm_data_complete_tp53$Patient_ID %in% 
                                                                                cptac_basal_samples,]
cptac_brca_fpkm_data_complete_tp53_her2 <- cptac_brca_fpkm_data_complete_tp53[cptac_brca_fpkm_data_complete_tp53$Patient_ID %in% 
                                                                                cptac_her2_samples,]

cptac_brca_proteomic_data_m_lumA <- cptac_brca_proteomic_data_m[cptac_brca_proteomic_data_m$Patient_ID %in%
                                                                  cptac_lumA_samples,]
cptac_brca_proteomic_data_m_lumB <- cptac_brca_proteomic_data_m[cptac_brca_proteomic_data_m$Patient_ID %in%
                                                                  cptac_lumB_samples,]
cptac_brca_proteomic_data_m_basal <- cptac_brca_proteomic_data_m[cptac_brca_proteomic_data_m$Patient_ID %in%
                                                                  cptac_basal_samples,]
cptac_brca_proteomic_data_m_her2 <- cptac_brca_proteomic_data_m[cptac_brca_proteomic_data_m$Patient_ID %in%
                                                                  cptac_her2_samples,]

pvalues <- wilcox_of_targs(cptac_brca_fpkm_data_complete_tp53_lumA, "TP53", "FPKM")
pvalues <- wilcox_of_targs(cptac_brca_proteomic_data_m_lumA, "TP53", "Prot.Quant")
