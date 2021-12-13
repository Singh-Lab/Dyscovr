Interrogate CPTAC Data from cBioPortal (Script)

cptac_brca_mutation_data <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/cBioPortal/brca_cptac_2020/data_mutations.txt", sep = "\t", header = TRUE, check.names = FALSE)
cptac_brca_missense_mutation_data <- cptac_brca_mutation_data[cptac_brca_mutation_data$Variant_Classification == "Missense_Mutation",]

cptac_brca_fpkm_data <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/cBioPortal/brca_cptac_2020/data_mrna_seq_fpkm.txt", sep = "\t", header = TRUE, check.names = FALSE)
cptac_brca_fpkm_data_complete <- cptac_brca_fpkm_data[rowSums(is.na(cptac_brca_fpkm_data[,2:ncol(cptac_brca_fpkm_data)])) == 0,]

cptac_brca_tp53_mutation <- cptac_brca_missense_mutation_data[cptac_brca_missense_mutation_data$Hugo_Symbol == "TP53",]
tp53_mut_patients <- unique(cptac_brca_tp53_mutation$Tumor_Sample_Barcode)
length(tp53_mut_patients)
[1] 25
tp53_not_mut_patients <- setdiff(unique(cptac_brca_missense_mutation_data$Tumor_Sample_Barcode), tp53_mut_patients)

cptac_brca_pik3ca_mutation <- cptac_brca_missense_mutation_data[cptac_brca_missense_mutation_data$Hugo_Symbol == "PIK3CA",]
pik3ca_mut_patients <- unique(cptac_brca_pik3ca_mutation$Tumor_Sample_Barcode)
length(pik3ca_mut_patients)
[1] 35
pik3ca_not_mut_patients <- setdiff(unique(cptac_brca_missense_mutation_data$Tumor_Sample_Barcode), pik3ca_mut_patients)

# MTHFD1L + TP53
ensg <- "ENSG00000120254"
symbol <- "MTHFD1L"

expression_mutants <- as.numeric(cptac_brca_fpkm_data_complete[cptac_brca_fpkm_data_complete$Hugo_Symbol == symbol, colnames(cptac_brca_fpkm_data_complete) %in% tp53_mut_patients])
expression_normal <- as.numeric(cptac_brca_fpkm_data_complete[cptac_brca_fpkm_data_complete$Hugo_Symbol == symbol, !colnames(cptac_brca_fpkm_data_complete) %in% tp53_mut_patients])
tp53_expression_by_mut_group <- list("Mutation" = expression_mutants, "No Mutation" = expression_normal)
boxplot(tp53_expression_by_mut_group, ylab = "Expression (FPKM)", 
        main = "MTHFD1L Expression By TP53 Mutation Status")

wilcox.test(expression_mutants, y = expression_normal)

# TYMS + PIK3CA
ensg <- "ENSG00000176890"
symbol <- "TYMS"

expression_mutants <- as.numeric(unlist(cptac_brca_fpkm_data[cptac_brca_fpkm_data$Hugo_Symbol == symbol, colnames(cptac_brca_fpkm_data) %in% pik3ca_mut_patients]))
expression_normal <- as.numeric(unlist(cptac_brca_fpkm_data[cptac_brca_fpkm_data$Hugo_Symbol == symbol, !colnames(cptac_brca_fpkm_data) %in% pik3ca_mut_patients]))
pik3ca_expression_by_mut_group <- list("Mutation" = expression_mutants, "No Mutation" = expression_normal)
boxplot(pik3ca_expression_by_mut_group, ylab = "Expression (FPKM)", 
        main = "TYMS Expression By PIK3CA Mutation Status")

wilcox.test(expression_mutants, y = expression_normal)

# TP53 + others
symbol <- "SHMT2"
symbol <- "GSTP1"
