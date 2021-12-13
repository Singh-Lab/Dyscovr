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

###################################################################################
# DRUG SENSITIVITY ANALYSIS
###################################################################################
expression_metformin <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/depmap_metformin.csv", header = TRUE, check.names = FALSE)
expression_brca_metformin <- expression_metformin[grepl("Breast", expression_metformin$Lineage),]

ccle_mutations <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/DepMap CCLE/CCLE_mutations.csv", header = TRUE, check.names = FALSE)
ccle_missense_mutations <- ccle_mutations[ccle_mutations$Variant_Classification == "Missense_Mutation",]
ccle_tp53_missense_mutations <- ccle_missense_mutations[ccle_missense_mutations$Hugo_Symbol == "TP53",]
cell_lines_tp53_mut <- unique(ccle_tp53_missense_mutations$DepMap_ID)

# Match mutation status to cell line
expression_brca_metformin$TP53_Mut <- unlist(lapply(expression_brca_metformin$`Depmap ID`, function(id) {
     if(id %in% cell_lines_tp53_mut) {return(TRUE)}
     else {return(FALSE)}
}))


# Plot as a boxplot
sensitivity_mutants <- unlist(expression_brca_metformin[expression_brca_metformin$TP53_Mut == TRUE, 2])
sensitivity_normal <- unlist(expression_brca_metformin[expression_brca_metformin$TP53_Mut == FALSE, 2])
boxplot(list("Mutation" = sensitivity_mutants, "No Mutation" = sensitivity_normal), 
         ylab = "Drug sensitivity", main = "Sensitivity to Metformin in BRCA Cell Lines, by TP53 Missense Mutation Status")

# Wilcoxon
wilcox.test(sensitivity_mutants, y = sensitivity_normal)

ccle_drug_sensitivity <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/DepMap CCLE/drug_sensitivity_screen.csv", header = TRUE, check.names = FALSE)
fu_sensitiv_brca <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/DepMap/5fluorouracil_sensitivity_brca.csv", header = TRUE, check.names = FALSE)
fu_sensitiv_brca$TP53_Mut <- unlist(lapply(fu_sensitiv_brca$`DepMap ID`, function(id) {
     if(id %in% cell_lines_tp53_mut) {return(TRUE)}
     else {return(FALSE)}
}))

tp53_sensitivity_mutants <- unlist(fu_sensitiv_brca[fu_sensitiv_brca$TP53_Mut == TRUE, 2])
tp53_sensitivity_normal <- unlist(fu_sensitiv_brca[fu_sensitiv_brca$TP53_Mut == FALSE, 2])
boxplot(list("Mutation" = tp53_sensitivity_mutants, "No Mutation" = tp53_sensitivity_normal), 
         ylab = "Drug sensitivity", main = "Sensitivity to 5-FU in BRCA Cell Lines, by TP53 Missense Mutation Status")

# Wilcoxon
wilcox.test(tp53_sensitivity_mutants, y = tp53_sensitivity_normal)

fu_sensitiv_brca$PIK3CA_Mut <- unlist(lapply(fu_sensitiv_brca$`DepMap ID`, function(id) {
    if(id %in% cell_lines_pik3ca_mut) {return(TRUE)}
    else {return(FALSE)}    
}))
pik3ca_sensitivity_mutants <- unlist(fu_sensitiv_brca[fu_sensitiv_brca$PIK3CA_Mut == TRUE, 2])
pik3ca_sensitivity_normal <- unlist(fu_sensitiv_brca[fu_sensitiv_brca$PIK3CA_Mut == FALSE, 2])
boxplot(list("Mutation" = pik3ca_sensitivity_mutants, "No Mutation" = pik3ca_sensitivity_normal), 
         ylab = "Drug sensitivity", main = "Sensitivity to 5-FU in BRCA Cell Lines, by PIK3CA Missense Mutation Status")

# Wilcoxon
wilcox.test(pik3ca_sensitivity_mutants, y = pik3ca_sensitivity_normal)

tp53_pik3ca_sensitivity_mutants <- unlist(fu_sensitiv_brca[(fu_sensitiv_brca$TP53_Mut == TRUE) & (fu_sensitiv_brca$PIK3CA_Mut == TRUE), 2])
tp53_pik3ca_sensitivity_normal <- unlist(fu_sensitiv_brca[(fu_sensitiv_brca$TP53_Mut == FALSE) & (fu_sensitiv_brca$PIK3CA_Mut == FALSE), 2])

tp53_excl_sensitivity_mutants <- setdiff(tp53_sensitivity_mutants, tp53_pik3ca_sensitivity_mutants)
pik3ca_excl_sensitivity_mutants <- setdiff(pik3ca_sensitivity_mutants, tp53_pik3ca_sensitivity_mutants)

boxplot(list("TP53 Mut Only" = tp53_excl_sensitivity_mutants, "PIK3CA Mut Only" = pik3ca_excl_sensitivity_mutants,
              "Both" = tp53_pik3ca_sensitivity_mutants, "Neither" = tp53_pik3ca_sensitivity_normal), 
         ylab = "Drug sensitivity", main = "Sensitivity to 5-FU in BRCA Cell Lines, by Missense Mutation Status")

# ANOVA
table_for_anova <- data.frame("sensitivity" = fu_sensitiv_brca$`5-fluorouracil sensitivity`, "status" = unlist(lapply(1:nrow(fu_sensitiv_brca), function(i) {
    if(fu_sensitiv_brca$TP53_Mut[i] & fu_sensitiv_brca$PIK3CA_Mut[i]) {return("Both")}
    else if(!fu_sensitiv_brca$TP53_Mut[i] & !fu_sensitiv_brca$PIK3CA_Mut[i]) {return("Neither")}
    else if (fu_sensitiv_brca$TP53_Mut[i] & !fu_sensitiv_brca$PIK3CA_Mut[i]) {return("TP53.Mut")}
    else if (!fu_sensitiv_brca$TP53_Mut[i] & fu_sensitiv_brca$PIK3CA_Mut[i]) {return("PIK3CA.Mut")}
    else {return(NA)}
})))
aov_res <- aov(sensitivity ~ status, data = table_for_anova)
summary(aov_res)


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

# Use DepMap features to create an aggregate file for all BRCA cell lines and for a particular gene set of interest
curr_dir <- getwd()
signif_hits_crispr_data <- read.csv(paste0(curr_dir, "/2020-2021/Research Data/DepMap/CRISPR_(DepMap_21Q4_Public+Score,_Chronos)_subsetted.csv"),
                                header = TRUE, check.names = FALSE)
signif_hits_rnai_data <- read.csv(paste0(curr_dir, "/2020-2021/Research Data/DepMap/RNAi_(Achilles+DRIVE+Marcotte,_DEMETER2)_subsetted.csv"),
                                  header = TRUE, check.names = FALSE)
signif_hits_expression_data <- read.csv(paste0(curr_dir, "/2020-2021/Research Data/DepMap/Expression_21Q4_Public_subsetted.csv"),
                                        header = TRUE, check.names = FALSE)

# Melt data frames to get them ready for boxplots
signif_hits_crispr_data_melt <- melt(signif_hits_crispr_data)
signif_hits_rnai_data_melt <- melt(signif_hits_rnai_data)
signif_hits_expression_data_melt <- melt(signif_hits_expression_data)

# Adjust column names
colnames(signif_hits_crispr_data_melt)[which(colnames(signif_hits_crispr_data_melt) == "variable")] <- "gene"
colnames(signif_hits_rnai_data_melt)[which(colnames(signif_hits_rnai_data_melt) == "variable")] <- "gene"
colnames(signif_hits_expression_data_melt)[which(colnames(signif_hits_expression_data_melt) == "variable")] <- "gene"

# Make the "value" column numeric, the "gene" column a character
signif_hits_crispr_data_melt$value <- as.numeric(unlist(signif_hits_crispr_data_melt$value))
signif_hits_rnai_data_melt$value <- as.numeric(unlist(signif_hits_rnai_data_melt$value))
signif_hits_expression_data_melt$value <- as.numeric(unlist(signif_hits_expression_data_melt$value))

signif_hits_crispr_data_melt$gene <- as.character(unlist(signif_hits_crispr_data_melt$gene))
signif_hits_rnai_data_melt$gene <- as.character(unlist(signif_hits_rnai_data_melt$gene))
signif_hits_expression_data_melt$gene <- as.character(unlist(signif_hits_expression_data_melt$gene))

# Get the mutations for TP53 & PIK3CA among these cell lines
tp53_pik3ca_mutations <- read.csv(paste0(curr_dir, "/2020-2021/Research Data/DepMap/Mutation_21Q4_Public_subsetted.csv"),
                                  header = TRUE, check.names = FALSE)

# Get the copy number changes for TP53 & PIK3CA among these cell lines
tp53_pik3ca_cnas <- read.csv(paste0(curr_dir, "/2020-2021/Research Data/DepMap/Copy_Number_21Q4_Public_subsetted.csv"),
                             header = TRUE, check.names = FALSE)
# Adjust column names
colnames(tp53_pik3ca_cnas)[ncol(tp53_pik3ca_cnas)-1] <- "TP53.CNA"
colnames(tp53_pik3ca_cnas)[ncol(tp53_pik3ca_cnas)] <- "PIK3CA.CNA"

# Add TP53 & PIK3CA mutation status & CNA status to the above data sets
signif_hits_crispr_data_melt <- merge(signif_hits_crispr_data_melt, tp53_pik3ca_mutations[,c("depmap_id", "TP53", "PIK3CA")], by = "depmap_id")
signif_hits_rnai_data_melt <- merge(signif_hits_rnai_data_melt, tp53_pik3ca_mutations[,c("depmap_id", "TP53", "PIK3CA")], by = "depmap_id")
signif_hits_expression_data_melt <- merge(signif_hits_expression_data_melt, tp53_pik3ca_mutations[,c("depmap_id", "TP53", "PIK3CA")], by = "depmap_id")

signif_hits_crispr_data_melt <- merge(signif_hits_crispr_data_melt, tp53_pik3ca_cnas[,c("depmap_id", "TP53.CNA", "PIK3CA.CNA")], by = "depmap_id")
signif_hits_rnai_data_melt <- merge(signif_hits_rnai_data_melt, tp53_pik3ca_cnas[,c("depmap_id", "TP53.CNA", "PIK3CA.CNA")], by = "depmap_id")
signif_hits_expression_data_melt <- merge(signif_hits_expression_data_melt, tp53_pik3ca_cnas[,c("depmap_id", "TP53.CNA", "PIK3CA.CNA")], by = "depmap_id")

# Make the TP53/ PIK3CA mutation columns factors
signif_hits_crispr_data_melt$TP53 <- as.factor(unlist(signif_hits_crispr_data_melt$TP53))
signif_hits_rnai_data_melt$TP53 <- as.factor(unlist(signif_hits_rnai_data_melt$TP53))
signif_hits_expression_data_melt$TP53 <- as.factor(unlist(signif_hits_expression_data_melt$TP53))

signif_hits_crispr_data_melt$PIK3CA <- as.factor(unlist(signif_hits_crispr_data_melt$PIK3CA))
signif_hits_rnai_data_melt$PIK3CA <- as.factor(unlist(signif_hits_rnai_data_melt$PIK3CA))
signif_hits_expression_data_melt$PIK3CA <- as.factor(unlist(signif_hits_expression_data_melt$PIK3CA))

# Make the TP53/ PIK3CA CNA columns numeric
#signif_hits_crispr_data_melt$TP53.CNA <- as.numeric(unlist(signif_hits_crispr_data_melt$TP53.CNA))
#signif_hits_rnai_data_melt$TP53.CNA <- as.numeric(unlist(signif_hits_rnai_data_melt$TP53.CNA))
#signif_hits_expression_data_melt$TP53.CNA <- as.numeric(unlist(signif_hits_expression_data_melt$TP53.CNA))

#signif_hits_crispr_data_melt$PIK3CA.CNA <- as.numeric(unlist(signif_hits_crispr_data_melt$PIK3CA.CNA))
#signif_hits_rnai_data_melt$PIK3CA.CNA <- as.numeric(unlist(signif_hits_rnai_data_melt$PIK3CA.CNA))
#signif_hits_expression_data_melt$PIK3CA.CNA <- as.numeric(unlist(signif_hits_expression_data_melt$PIK3CA.CNA))

# Alternatively, convert CNA to be bucketed (<1.098 is a deletion, 1.098-1.1 is around normal, > 1.1 is amplification)
bucket_cna <- function(cna_vals, deletionOrAmp) {
    bucketed_cna_vals <- unlist(lapply(cna_vals, function(x) {
        if(deletionOrAmp == "deletion") {
            if (x < 1.098) {return(0)}
            else {return(1)}
        } else {
            if (x <= 1.099) {return(0)}
            else {return(1)}
        }
    }))
    return(bucketed_cna_vals)
}

signif_hits_crispr_data_melt$TP53.CNA <- as.factor(bucket_cna(signif_hits_crispr_data_melt$TP53.CNA, "deletion"))
signif_hits_rnai_data_melt$TP53.CNA <- as.factor(bucket_cna(signif_hits_rnai_data_melt$TP53.CNA, "deletion"))
signif_hits_expression_data_melt$TP53.CNA <- as.factor(bucket_cna(signif_hits_expression_data_melt$TP53.CNA, "deletion"))

signif_hits_crispr_data_melt$PIK3CA.CNA <- as.factor(bucket_cna(signif_hits_crispr_data_melt$PIK3CA.CNA, "amplification"))
signif_hits_rnai_data_melt$PIK3CA.CNA <- as.factor(bucket_cna(signif_hits_rnai_data_melt$PIK3CA.CNA, "amplification"))
signif_hits_expression_data_melt$PIK3CA.CNA <- as.factor(bucket_cna(signif_hits_expression_data_melt$PIK3CA.CNA, "amplification"))


# Subset to just the top hits for TP53 & PIK3CA specifically
tp53_top_hits <- c("MTHFD1L", "SHMT2", "GSTP1", "TYMS", "PSPH", "CDO1", "MTHFD1", "CSAD", "AHCY", 
                   "MTHFR", "MTR", "FTCD")
# additional hits from just TP53/ PIK3CA run: "SARDH", "DLD", "MAT2A", "GLDC"
tp53_top_hits_specific <- c(tp53_top_hits, c("SARDH", "DLD", "MAT2A", "GLDC"))

pik3ca_top_hits <- c("TYMS", "MTHFD1L", "PSPH", "GLDC")
# additional hits from just TP53/ PIK3CA run: "SHMT2", "FTCD", "MAT1A", "MTHFR", "AHCY", "CHKB", "MAT2A", "BHMT2", "GCLM"
pik3ca_top_hits_specific <- c(pik3ca_top_hits, c("SHMT2", "FTCD", "MAT1A", "MTHFR", "AHCY", "CHKB", "MAT2A", "BHMT2", "GCLM"))
    
signif_hits_crispr_data_melt_tp53 <- signif_hits_crispr_data_melt[signif_hits_crispr_data_melt$gene %in% tp53_top_hits,]
signif_hits_rnai_data_melt_tp53 <- signif_hits_rnai_data_melt[signif_hits_rnai_data_melt$gene %in% tp53_top_hits,]
signif_hits_expression_data_melt_tp53 <- signif_hits_expression_data_melt[signif_hits_expression_data_melt$gene %in% tp53_top_hits,]

signif_hits_crispr_data_melt_pik3ca <- signif_hits_crispr_data_melt[signif_hits_crispr_data_melt$gene %in% pik3ca_top_hits,]
signif_hits_rnai_data_melt_pik3ca <- signif_hits_rnai_data_melt[signif_hits_rnai_data_melt$gene %in% pik3ca_top_hits,]
signif_hits_expression_data_melt_pik3ca <- signif_hits_expression_data_melt[signif_hits_expression_data_melt$gene %in% pik3ca_top_hits,]

# Run a significance test for each combo
run_ttest <- function(df) {
    stat.test <- df %>%
        group_by(gene) %>%
        #wilcox.test(value ~ TP53) %>%
        t_test(value ~ TP53) %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")
    return(stat.test)
}

ttest_crispr_tp53 <- run_ttest(signif_hits_crispr_data_melt_tp53)
ttest_rnai_tp53 <- run_ttest(signif_hits_rnai_data_melt_tp53)
ttest_expression_tp53 <- run_ttest(signif_hits_expression_data_melt_tp53)

# Plot boxplots of each gene's dependency across BRCA cell lines, by TP53 or PIK3CA mutation status
crispr_bxp_tp53 <- ggplot(signif_hits_crispr_data_melt_tp53, aes(x = gene, y = value, fill = TP53)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
rnai_bxp_tp53 <- ggplot(signif_hits_rnai_data_melt_tp53, aes(x = gene, y = value, fill = TP53)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
express_bxp_tp53 <- ggplot(signif_hits_expression_data_melt_tp53, aes(x = gene, y = value, fill = TP53)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

crispr_bxp_pik3ca <- ggplot(signif_hits_crispr_data_melt_pik3ca, aes(x = gene, y = value, fill = PIK3CA)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
rnai_bxp_pik3ca <- ggplot(signif_hits_rnai_data_melt_pik3ca, aes(x = gene, y = value, fill = PIK3CA)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
express_bxp_pik3ca <- ggplot(signif_hits_expression_data_melt_pik3ca, aes(x = gene, y = value, fill = PIK3CA)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Add p-values onto the box plots
ttest_crispr_tp53 <- ttest_crispr_tp53 %>% add_xy_position(x = "gene", dodge = 0.8)
ttest_rnai_tp53 <- ttest_rnai_tp53 %>% add_xy_position(x = "gene", dodge = 0.8)
ttest_expression_tp53 <- ttest_expression_tp53 %>% add_xy_position(x = "gene", dodge = 0.8)

crispr_bxp_tp53 + stat_pvalue_manual(ttest_crispr_tp53,  label = "p", tip.length = 0) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
rnai_bxp_tp53 + stat_pvalue_manual(ttest_rnai_tp53,  label = "p", tip.length = 0) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
express_bxp_tp53 + stat_pvalue_manual(ttest_expression_tp53,  label = "p", tip.length = 0) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Plot boxplots of each gene's dependency across BRCA cell lines, by TP53 or PIK3CA CNA status
crispr_bxp_tp53_cna <- ggplot(signif_hits_crispr_data_melt_tp53, aes(x = gene, y = value, fill = TP53.CNA)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_discrete(name = "TP53 CNA", labels = c("Deletion", "Normal/ Amplif."))
rnai_bxp_tp53_cna <- ggplot(signif_hits_rnai_data_melt_tp53, aes(x = gene, y = value, fill = TP53.CNA)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_discrete(name = "TP53 CNA", labels = c("Deletion", "Normal/ Amplif."))
express_bxp_tp53_cna <- ggplot(signif_hits_expression_data_melt_tp53, aes(x = gene, y = value, fill = TP53.CNA)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_discrete(name = "TP53 CNA", labels = c("Deletion", "Normal/ Amplif."))

crispr_bxp_pik3ca_cna <- ggplot(signif_hits_crispr_data_melt_pik3ca, aes(x = gene, y = value, fill = PIK3CA.CNA)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_discrete(name = "PIK3CA CNA", labels = c("Normal/Deletion", "Amplification"))
rnai_bxp_pik3ca_cna <- ggplot(signif_hits_rnai_data_melt_pik3ca, aes(x = gene, y = value, fill = PIK3CA.CNA)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_discrete(name = "PIK3CA CNA", labels = c("Normal/Deletion", "Amplification"))
express_bxp_pik3ca_cna <- ggplot(signif_hits_expression_data_melt_pik3ca, aes(x = gene, y = value, fill = PIK3CA.CNA)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_discrete(name = "PIK3CA CNA", labels = c("Normal/Deletion", "Amplification"))


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
    plt_input$diffDep[plt_input$Betas > 0 & plt_input$p.val < 0.05] <- "UP"
    plt_input$diffDep[plt_input$Betas < 0 & plt_input$p.val < 0.05] <- "DOWN"

    plt_input$diffDepLabel <- NA
    plt_input$diffDepLabel[plt_input$diffDep != "NONE"] <- plt_input$Gene[plt_input$diffDep != "NONE"]
    
    print(plt_input)
    
    plt <- ggplot(plt_input, aes(x = Betas, y = Neg.Log_p.val, col = diffDep, label = diffDepLabel)) + geom_point() +
        geom_hline(yintercept = -log10(0.05), col="red") + theme_minimal() +
        geom_text_repel() + scale_color_manual(values=c("blue", "black", "red")) +
        ggtitle(paste(gene_name, mut_or_cna))
    print(plt)
}

depmap_lm(signif_hits_crispr_data_melt_tp53, "mut", "TP53")
depmap_lm(signif_hits_crispr_data_melt_tp53, "cna", "TP53")
depmap_lm(signif_hits_crispr_data_melt_pik3ca, "mut", "PIK3CA")
depmap_lm(signif_hits_crispr_data_melt_pik3ca, "cna", "PIK3CA")

depmap_lm(signif_hits_crispr_data_melt, "mut", "TP53")


