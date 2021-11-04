##############################################################################
### MAFTOOLS EXPLORATORY ANALYSIS ###
##############################################################################

# This file uses maftools package to examine and plot characteristics of a
# given MAF file

# Link to Tutorial/ Vignette: https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html

BiocManager::install("maftools")
library(maftools)
library(TCGAbiolinks)
library("gplots")
library(ggplot2)

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/Somatic_Mut_Data/"

##############################################################################
### BRCA DATA ###
##############################################################################

# Read in maf file
maf_filename <- "TCGA.BRCA.muse.b8ca5856-9819-459c-87c5-94e91aca4032.DR-10.0.somatic.maf"
maf <- read.maf(maf = paste0(path, maf_filename)) # clinicalData = mut_clin_file)

# Print summaries of the maf
print(maf)  # general summary
print(getSampleSummary(maf))  # sample summary
print(getGeneSummary(maf))    # gene summary

# Visualize maf file
plotmafSummary(maf = maf, rmOutlier = FALSE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# Remove some of the top mutated proteins that are not driver genes
maf_dt <- fread(paste0(path, maf_filename), header= TRUE)
maf_df <- as.data.frame(maf_dt)
maf_df_sub <- maf_df[maf_df$Hugo_Symbol %in% c("TP53", "PIK3CA", "KMT2C", "MAP2K4", "FOXA1", "FAT1", 
                                               "ERBB2", "PTEN", "HUWE1", "MAPK1", "AKT1"),]
maf_df_sub <- maf_df[!(maf_df$Hugo_Symbol %in% c("TTN", "MUC16", "HMCN1", "CDH1", "USH2A", "MUC4", 
                                                 "RYR2", "FLG", "SYNE1", "MUC17", "SPTA1", "NCOR1", 
                                                 "NEB", "ZFHX4", "RYR3", "DMD", "SYNE2", "CSMD3", "FAT3",
                                                 "PKHD1L1", "DNAH11", "LRP2", "MDN1", "XIRP2", "ABCA13",
                                                 "ADGRV1", "APOB", "CSMD1", "CCDC168", "OBSCN", 
                                                 "PCLO", "DNAH17", "CACNA1E", "PRUNE2", "KIAA1109",
                                                 "LRP1B", "AHNAK", "MYCBP2", "PRKDC", "DNAH9", "RELN", 
                                                 "NF1", "CSMD2", "UTRN", "MXRA5", "HERC2", "FLG2")),]
maf_sub <- read.maf(maf_df_sub)
#write.table(maf_df_sub, paste0(path, "TCGA.BRCA.muse.b8ca5856-9819-459c-87c5-94e91aca4032.DR-10.0.somatic.topDrivers.maf"), 
            #sep = "\t", quote = FALSE)
plotmafSummary(maf = maf_sub, rmOutlier = FALSE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)


# Draw an oncoplot/ waterfall plot for the maf file
oncoplot(maf = maf, top = 10)

# Lollipop plot for PIK3CA and TP53, top mutated genes
lollipopPlot(maf = maf, gene = 'PIK3CA', showMutationRate = TRUE)
lollipopPlot(maf = maf, gene = 'TP53', showMutationRate = TRUE)
lollipopPlot(maf = maf, gene = "CDH1", showMutationRate = TRUE)

# Compare mutation load against other TCGA cohorts
maf.mutload = tcgaCompare(maf = maf, cohortName = 'BRCA')

# Rainfall plot for detecting hyper-mutated loci
rainfallPlot(maf = maf, detectChangePoints = TRUE, pointSize = 0.6)

# Variant Allele Frequencies Plot (clonal status of top mutated genes)
plotVaf(maf = maf)

# Somatic Interactions Visualization
  # Uses pair-wise Fisher's Exact Test to detect significant pairs of genes
somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1))

# ID Cancer Driver Genes Using Frequency-Based Algorithm
brca.sig = oncodrive(maf = maf, minMut = 5, pvalMethod = 'zscore')
head(brca.sig)
  # TOP 6: AKT1, KRAS, PIK3CA, DGKB, PAXIP1, GRXCR1
plotOncodrive(res = brca.sig, fdrCutOff = 0.05, useFraction = TRUE, bubbleSize = 1, labelSize = 0.7)

# Adding and summarizing Pfam domains
brca.pfam = pfamDomains(maf = maf, top = 6, labelSize = 0.7)
brca.pfam$proteinSummary[,1:7, with = FALSE]  # Protein summary
brca.pfam$domainSummary[,1:3, with = FALSE]   # Domain summary

# Store the protein and domain summaries as objects - these will be useful to us!
brca_protein_summary <- brca.pfam$proteinSummary
brca_domain_summary <- brca.pfam$domainSummary

# Oncogenic pathways
OncogenicPathways(maf = maf)

# NOTE: easy to subset the maf using subsetMaf function (?subsetMaf)
maf_subset_missense <- subsetMaf(maf = maf, query = "Variant_Classification == 'Missense_Mutation'")
#maf_subset_silent <- subsetMaf(maf = maf, query = "Variant_Classification == 'Silent'")
#maf_subset_mis_and_sil <- subsetMaf(maf = maf, query = "Variant_Classification == 'Missense_Mutation' |
                                    #Variant_Classification == 'Silent'")


##############################################################################
### GBM DATA ###
##############################################################################
# Read in maf file
gbm_maf_filename <- "TCGA.GBM.muse.59a84472-27d4-497c-8f37-8bc447ff9374.DR-10.0.somatic.maf"
gbm_maf <- read.maf(maf = gbm_maf_filename) 

# Print summaries of the maf
print(gbm_maf)  # general summary
print(getSampleSummary(gbm_maf))  # sample summary
print(getGeneSummary(gbm_maf))    # gene summary

# Visualize maf file
plotmafSummary(maf = gbm_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# Draw an oncoplot/ waterfall plot for the maf file
oncoplot(maf = gbm_maf, top = 10)

# Lollipop plot for PIK3CA and TP53, top mutated genes
lollipopPlot(maf = gbm_maf, gene = 'PTEN', showMutationRate = TRUE)
lollipopPlot(maf = gbm_maf, gene = 'TP53', showMutationRate = TRUE)
lollipopPlot(maf = gbm_maf, gene = "TTN", showMutationRate = TRUE)

# Compare mutation load against other TCGA cohorts
maf.mutload = tcgaCompare(maf = gbm_maf, cohortName = 'GBM')

# Rainfall plot for detecting hyper-mutated loci
rainfallPlot(maf = gbm_maf, detectChangePoints = TRUE, pointSize = 0.6)

# Variant Allele Frequencies Plot (clonal status of top mutated genes)
plotVaf(maf = gbm_maf)

# Somatic Interactions Visualization
# Uses pair-wise Fisher's Exact Test to detect significant pairs of genes
somaticInteractions(maf = gbm_maf, top = 25, pvalue = c(0.05, 0.1))

# ID Cancer Driver Genes Using Frequency-Based Algorithm
gbm.sig = oncodrive(maf = gbm_maf, minMut = 5, pvalMethod = 'zscore')
head(gbm.sig)
# TOP 6: AKT1, KRAS, PIK3CA, DGKB, PAXIP1, GRXCR1
plotOncodrive(res = gbm.sig, fdrCutOff = 0.05, useFraction = TRUE, bubbleSize = 1, labelSize = 0.7)

# Adding and summarizing Pfam domains
gbm.pfam = pfamDomains(maf = gbm_maf, top = 7, labelSize = 0.7)
gbm.pfam$proteinSummary[,1:7, with = FALSE]  # Protein summary
gbm.pfam$domainSummary[,1:3, with = FALSE]   # Domain summary

# Store the protein and domain summaries as objects - these will be useful to us!
gbm_protein_summary <- gbm.pfam$proteinSummary
gbm_domain_summary <- gbm.pfam$domainSummary

# Oncogenic pathways
OncogenicPathways(maf = gbm_maf)

# NOTE: easy to subset the maf using subsetMaf function (?subsetMaf)
maf_subset_missense <- subsetMaf(maf = maf, query = "Variant_Classification == 'Missense_Mutation'")
#maf_subset_silent <- subsetMaf(maf = maf, query = "Variant_Classification == 'Silent'")
#maf_subset_mis_and_sil <- subsetMaf(maf = maf, query = "Variant_Classification == 'Missense_Mutation' |
#Variant_Classification == 'Silent'")

