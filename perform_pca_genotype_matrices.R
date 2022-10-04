############################################################
### Take PCA of Genotype Matrices
### Written By: Sara Geraghty, June 2021
############################################################

# This file imports the unmasked MAF files for all patients, across cancer types.
# This is protected data in the GDC ('Aggregated Somatic Mutation', WXS, processed
# using the MUSE pipeline).

# More information about GDC pre-processing and mutation calling: 
# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#masked-somatic-aggregation-workflow

# To perform PCA, this file uses the two built-in R methods to perform PCA
# princomp() : spectral decomposition (examines the covariances/ correlations between variables)
# prcomp() : singular value decomposition (examines the covariances/ correlations between individuals)
# NOTE: prcomp() is the preferred method due to slightly better numerical accuracy

# FOR PCA STUFF
# SeqArray package tutorial: https://bioc.ism.ac.jp/packages/3.5/bioc/vignettes/SeqArray/inst/doc/R_Integration.html

#library(SeqArray)
#library(factoextra)   # for visualization of PCA results
library(TRONCO)
library(Rcpp)
#library(vcfR)
#library(pcadapt)
#library(gdsfmt)
#library(SNPRelate) 
library(data.table)
library(dplyr)
library(ggpubr)
library(reshape2)

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Genotype_Mut_Data/"


############################################################
### MAKE USE OF PCs THAT HAVE ALREADY BEEN CALLED ON GENOTYPE 
### FILES
############################################################
# Link to file downloads: https://gdc.cancer.gov/about-data/publications/CCG-AIM-2020
# Link to associated Carrot-Zhang Paper: https://www.cell.com/cancer-cell/pdf/S1535-6108(20)30211-7.pdf

# 1. Washington University
washu_pcs <- read.table(paste0(main_path, "WashU_PCA_ethnicity_assigned.tsv"), 
                        header = TRUE, check.names = FALSE, sep = "\t")

# 2. UCSF 
ucsf_pcs <- read.csv(paste0(main_path, "UCSF_Ancestry_Calls.csv"), 
                     header = TRUE, check.names = FALSE)

# 3. The Broad Institute
broad_pcs <- read.csv(paste0(main_path, "Broad_ancestry_PCA.txt"), 
                     header = TRUE, check.names = FALSE)


############################################################
### MANUALLY CREATE FILES FOR EIGENSTRAT's
### smartpca FUNCTION
############################################################

# Information about the input file formats: https://reich.hms.harvard.edu/software/InputFileFormats

# Three necessary file types: 
  # 1. Genotype file (1 line per SNP, each with one character per individual, no spaces)
    # 0 means no copies of reference allele, 1 means one copy of ref. allele, 
    # 2 means two copies, 9 is missing data
    # SNPS MUST MATCH THOSE FROM THE SNP FILE, IN SAME ORDER
  # 2. SNP file (1 SNP per line, with four columns)
    # SNP_ID (ex. rs0000)
    # Chromosome_Num (use 23 for X chromosome)
    # Genetic_Position (Morgans or centiMorgans; can use 0 for everything for smartpca)
    # Physical_Position (bases)
    # Previous letter
    # New letter
  # 3. Indiv. file (1 line per individual, with three columns)
    # Sample_ID (ex. SAMPLE0, any type of ID will do)
    # Gender (M (male), F (female), or U (unknown))
    # Status (Case, Control, or Ignore; can also use your own labels or not put anything)

# Currently using MAF file, which captures mutations -- do I want to use VCF files 
# for normal patients instead to capture normal variation?

# Define outpath
outpath <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Eigensoft Files/"
#outpath <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Eigensoft Files/"


############################################################
#' Given a MAF file and corresponding clinical file from the TCGA, creates a 
#' ".indiv" file with the above format to be used in Eigensoft's smartpca function.
#' @param maf an input maf file, either pan-cancer or for a particular cancer type
#' @param clinical_file an input clinical file from the TCGA that corresponds to
#' the same patient set as the maf file
#' @param ct_name the name of the cancer type in question, is "allCancers" if 
#' looking pan-cancer
create_indiv_file <- function(maf, clinical_file, ct_name) {
  
  # Get the duplicated patient IDs
  unique_tcga_tumor_barcodes <- unique(maf$Tumor_Sample_Barcode)
  patient_ids_from_tsbs <- unlist(lapply(unique_tcga_tumor_barcodes, function(x) 
    unlist(strsplit(x, "-", fixed = TRUE))[3]))
  indices_of_duplicates <- which(duplicated(patient_ids_from_tsbs))
  ids_of_duplicates <- unique(patient_ids_from_tsbs[indices_of_duplicates])
  
  # Get the gender information and ID for each sample
  new_row_list <- lapply(1:nrow(maf), function(i) {
    # Get sample ID
    full_id <- maf[i, 'Tumor_Sample_Barcode']
    samp_id <- paste(unlist(strsplit(full_id, "-", fixed = TRUE))[1:3], 
                     collapse = "-")
    patient_id <- unlist(strsplit(samp_id, "-", fixed = TRUE))[3]
    
    # Get the gender corresponding to this ID
    gender <- unique(unlist(clinical_file[clinical_file$case_submitter_id == samp_id, 'gender']))
    
    # Handle duplicate patient cases
    if (patient_id %fin% ids_of_duplicates) {
      # When it comes to this patient, is this the first, second, third, etc. time we're 
      # encountering one of their samples?
      all_patient_samp_ids <- unique(maf$Tumor_Sample_Barcode[grepl(patient_id, maf$Tumor_Sample_Barcode)])
      which_samp <- which(full_id == all_patient_samp_ids)
      if(!which_samp == 1) {samp_id <- paste(samp_id, which_samp, sep = "-")}
    }
    
    print(samp_id)
    print(gender)
    print(ct_name)
    
    # Create a new row and return
    new_row <- data.frame(samp_id, gender, ct_name)
    
    return(new_row)
  })
  print(head(new_row_list))
  
  # Bind into a table for easier writing
  tab <- rbindlist(new_row_list)
  
  # Remove duplicated entries
  tab <- distinct(tab)
  
  return(tab)
}

############################################################

#' Given a MAF file, creates a ".gen" file of the format given above to be used
#' in Eigensoft's smartpca function
#' @param maf an input maf file, either pan-cancer or for a particular cancer type
create_genotype_file <- function(maf, snp_tab) {
  
  # Get all unique samples 
  samples <- unique(maf$Tumor_Sample_Barcode)
  
  # For each of the SNPs in the SNP file, get the genotypes of that SNP for 
  # every sample
  genotype_rows <- lapply(snp_tab$bp_position, function(bp) {
    # Subset the maf file to SNPs at this position
    maf_sub <- maf[maf$Start_Position == bp,]
    
    # For each patient, if they are not in the DF, give them a 2, if they are,
    # determine if they have 1 or 0 copies of the reference allele
    genotype <- unlist(lapply(samples, function(samp) {
      if (samp %fin% maf_sub$Tumor_Sample_Barcode) {
        samp_row <- maf_sub[maf_sub$Tumor_Sample_Barcode == samp,]
        if(samp_row$Tumor_Seq_Allele1 == samp_row$Reference_Allele) {
          if(samp_row$Tumor_Seq_Allele2 == samp_row$Reference_Allele) {return(2)}
          else {return(1)}
        } else {
          if (samp_row$Tumor_Seq_Allele2 == samp_row$Reference_Allele) {return(1)}
          else{return(0)}
        }
      } else {return(2)}
    }))
    
    # Make the genotype into a string to return
    genotype_str <- paste(genotype, collapse = "")
    
    return(genotype_str)
  })
  
  # Recombine this list of genotype strings per SNP into a data frame
  #genotype_tab <- do.call(rbind, genotype_rows)
  
  return(genotype_rows)
}


############################################################

#' Given a MAF file, creates a ".snp" file of the format given above to be used
#' in Eigensoft's smartpca function
#' @param maf an input maf file, either pan-cancer or for a particular cancer type
create_snp_file <- function(maf) {

  new_row_list <- lapply(1:nrow(maf), function(i) {
    
    # Get the information from the "vcf_region" column, if it exists
    # In the format: chrx:XXXXXXXX:rsYYYYYYYYY:C:T"
    if ('vcf_region' %in% colnames(maf)) {
      vcf_info <- maf[i, 'vcf_region']
      vcf_info_spl <- unlist(strsplit(vcf_info, ":", fixed = TRUE))
      
      # Extract individual necessary components
      chrom_num <- unlist(strsplit(vcf_info_spl[1], "r", fixed = TRUE))[2]
      bp_position <- vcf_info_spl[2]
      rsID <- vcf_info_spl[3]
      prior_base <- vcf_info_spl[4]
      new_base <- vcf_info_spl[5]
      
    } else {
      # Remove mitochondrial SNPs
      if (!((maf[i, 'Chromosome'] == "") | (maf[i, 'Chromosome'] == "M"))) {
        chrom_num <- unlist(strsplit(maf[i, 'Chromosome'], "r", fixed = TRUE))[2]
        bp_position <- maf[i, 'Start_Position']
        rsID <- maf[i, 'dbSNP_RS']
        # Replace unknown SNPs with a generic name so that we can retain them
        prior_base <- maf[i, 'Reference_Allele']
        new_base <- maf[i, 'Allele']
      }
    }
    
    if(rsID == "novel" | rsID == "" | rsID == ".") {
      rsID <- paste("rs", bp_position, sep = "")
    }   # is this necessary? Or should we eliminate them entirely?
    
    # Return as a list
    new_row <- data.frame(rsID, chrom_num, 0, bp_position, prior_base, new_base)
    
    return(new_row)
  })
  print(head(new_row_list))
  
  # Bind into a table for easier writing
  tab <- rbindlist(new_row_list)
  
  # Remove any duplicate SNPs found in multiple patients
  tab <- tab[!duplicated(tab)]
  
  
  return(tab)
}



#' Runs the above three functions to generate and write the necessary files for 
#' eigensoft's smartpca function
#' @param maf_f a mutation (.maf) file for the given cohort
#' @param ct_name a string with the name of the cancer type, or "allCancers" if 
#' we are looking at pan-cancer
#' @param clinical_file an associated clinical file from the TCGA
#' @param outpath a path to the location to write the eigensoft files to
get_eigensoft_files <- function(maf_f, ct_name, clinical_file, outpath) {
  
  # Create all three file types for this cancer type, & write each to a file at the outpath given
  snp_tab <- create_snp_file(maf_f)
  fwrite(snp_tab, paste(outpath, paste(ct_name, "_snp_file.snp", sep = ""), sep = ""), sep = " ")
  
  indiv_tab <- create_indiv_file(maf_f, clinical_file, ct_name)
  fwrite(indiv_tab, paste(outpath, paste(ct_name, "_indiv_file.indiv", sep = ""), sep = ""), sep = " ")
  
  genotype_tab <- create_genotype_file(maf_f, snp_tab)
  fwrite(genotype_tab, paste(outpath, paste(ct_name, "_genotype_file.eigenstratgeno", sep = ""), sep = ""), sep = "\n")
  
}


# BRCA
brca_matrix_fn <- paste(main_path, "MAF_files/TCGA.BRCA.muse.cd0d1636-ae95-4141-94b3-8218eb3b8b25.DR-10.0.protected.maf", sep = "")
brca_matrix <- import.MAF(brca_matrix_fn, sep = "\t", is.TCGA = TRUE, 
                       merge.mutation.types = FALSE, to.TRONCO = FALSE) 

# P-C
genotype_matrices_fns <- list.files(paste(main_path, "MAF_files/", sep = ""), pattern = ".muse")
genotype_matrix_list <- lapply(genotype_matrices_fns, function(fn) {
  full_fn <- paste(main_path, paste("MAF_files/", fn, sep = ""), sep = "")
  gn_matrix <- import.MAF(full_fn, sep = "\t", is.TCGA = TRUE, 
                          merge.mutation.types = FALSE, to.TRONCO = FALSE) 
  return(gn_matrix)
})

# For P-C, merge all the MAF files into one
genotype_matrix_pc <- do.call(rbind, genotype_matrix_list)


# Get the clinical files of interest
# TUMOR-NORMAL MATCHED
brca_clinical_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/clinical_data_subset.csv",
                          header = TRUE)
# NON-TUMOR-NORMAL MATCHED
brca_clinical_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/clinical.csv",
                          header = TRUE)
  
  
pc_clinical_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/clinical_data_subset.csv",
                        header = TRUE)


# Limit the MAF files to only include patients in the clinical DFs
brca_matrix$tcga_id <- unlist(lapply(brca_matrix$Tumor_Sample_Barcode, function(x) 
  paste(unlist(strsplit(x, "-", fixed = TRUE))[1:3], collapse = "-")))
#brca_matrix_sub <- brca_matrix[brca_matrix$tcga_id %fin% unique(brca_clinical_df$tcga_barcode),]
brca_matrix_sub <- brca_matrix[brca_matrix$tcga_id %fin% unique(brca_clinical_df$case_submitter_id),]

genotype_matrix_pc$tcga_id <- unlist(lapply(genotype_matrix_pc$Tumor_Sample_Barcode, function(x) 
  paste(unlist(strsplit(x, "-", fixed = TRUE))[1:3], collapse = "-")))
genotype_matrix_pc_sub <- genotype_matrix_pc[genotype_matrix_pc$tcga_id %fin% unique(pc_clinical_df$tcga_barcode),]


# Call this function

# BRCA
get_eigensoft_files(brca_matrix_sub, "BRCA", brca_clinical_df, outpath)

# Pan-Cancer
get_eigensoft_files(genotype_matrix_pc_sub, "allCancers", pc_clinical_df, outpath, TRUE)



############################################################
# IF NEEDED: add extra individuals to .indiv file to correspond to multiple samples
# per patient in the MAF file
############################################################
# Convert genotype to data frame
genotype_tab_mat <- as.matrix(genotype_tab)
genotype_tab_df_rows <- lapply(1:nrow(genotype_tab_mat), function(i) {
  row <- unlist(genotype_tab_mat[i,])
  return(as.numeric(unlist(strsplit(row, "", fixed = TRUE))))
})
genotype_tab_df <- do.call(rbind, genotype_tab_df_rows)

# Get the duplicated IDs
unique_tcga_tumor_barcodes <- unique(brca_matrix_sub$Tumor_Sample_Barcode)
patient_ids_from_tsbs <- unlist(lapply(unique_tcga_tumor_barcodes, function(x) unlist(strsplit(x, "-", fixed = TRUE))[3]))
which(duplicated(patient_ids_from_tsbs))
indices_of_duplicates <- which(duplicated(patient_ids_from_tsbs))
ids_of_duplicates <- unique(patient_ids_from_tsbs[indices_of_duplicates])

# Add fake columns with these ids (another with a -1 at the end)
add_additional_samples <- function(indiv_tab, ids_of_duplicates, brca_matrix_sub) {

  # Loop through all the IDs that have more than one tumor sample (and thus potential >1 genotype)
  for (i in 1:length(ids_of_duplicates)) {
    id <- ids_of_duplicates[i]
    print(id)
    
    # Get unique ids
    unique_ids <- unique(brca_matrix_sub[grepl(id, brca_matrix_sub$Tumor_Sample_Barcode), 
                           'Tumor_Sample_Barcode'])
    print(unique_ids)

    # Find the ID that came before these in the .indiv file so we can insert it 
    # in the proper place
    min_indic_of_other_samps <- unlist(lapply(unique_ids[2:length(unique_ids)], function(x) {
      indices <- which(x == brca_matrix_sub$Tumor_Sample_Barcode)
      return(min(indices))
    }))
    print(min_indic_of_other_samps)
    
    # For each of these IDs, get the ID before and insert an additional row after that one
    for (i in 1:length(min_indic_of_other_samps)) {
      min_index <- min_indic_of_other_samps[i]
      id_before <- paste(unlist(strsplit(brca_matrix_sub[min_index-1, 'Tumor_Sample_Barcode'], "-", fixed = TRUE))[1:3], 
                         collapse = "-")
      print(id_before)
      
      # Get the index of THIS ID in the indiv file
      prior_id_index_indiv <- which(grepl(id_before, indiv_tab$samp_id))
      print(prior_id_index_indiv)
      
      # Insert new row for additional sample HERE
      row_to_copy <- indiv_tab[grepl(id, indiv_tab$samp_id), ]
      if(nrow(row_to_copy) > 1) {row_to_copy <- row_to_copy[1,]}
      full_id <- as.character(unlist(row_to_copy[, 'samp_id']))
      gender <- as.character(unlist(row_to_copy[, 'gender']))
      ct <- as.character(unlist(row_to_copy[, 'ct_name']))
      
      new_row <- data.frame(paste(full_id, i, sep = "-"), gender, ct)
      colnames(new_row) <- c("samp_id", "gender", "ct_name")
      print(new_row)
      
      indiv_tab <- insertRow2(indiv_tab, new_row, (prior_id_index_indiv+1))
    }
  }
  return(indiv_tab)
}

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  return(existingDF)
}

insertRow2 <- function(existingDF, newrow, r) {
  existingDF <- rbind(existingDF,newrow)
  existingDF <- existingDF[order(c(1:(nrow(existingDF)-1),r-0.5)),]
  row.names(existingDF) <- 1:nrow(existingDF)
  return(existingDF)  
}

# Call function
indiv_tab_extraRows <- add_additional_samples(indiv_tab, ids_of_duplicates, brca_matrix_sub)

# Write to file 
write.table(indiv_tab_extraRows, paste(outpath, "BRCA_indiv_file.indiv", sep = ""), sep = " ", 
            col.names = FALSE, quote = FALSE, row.names = FALSE)



############################################################
### VISUALIZE PCA RESULTS
############################################################
pcs <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/smartpca/BRCA_output_eigenvect.pca",
                  header = TRUE)
# If needed
rownames(pcs) <- pcs[,1]

# Make a scatter plot of the first 2 PCs
ggplot(pcs, aes(x = V2, y = V3)) + geom_point() + xlab("PC1") + ylab("PC2") 

# Identify the names of outliers
ggplot(pcs, aes(x = V2, y = V3)) + geom_point() + 
  geom_text(aes(label= ifelse(V3 > quantile(V3, 0.95), as.character(V1),'')), hjust=1, vjust=0)

# Overlay clinical information
# Import clinical DFs 
brca_clinical_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/clinical.csv",
                          header = TRUE)
brca_clinical_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/clinical_data_subset.csv",
                          header = TRUE)
#pc_clinical_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/clinical_data_subset.csv",
#header = TRUE)

#' Extract information about age, gender, and race and add as rows for each patient
#' to the PCs DF
#' @param pcs the pca output from eigensoft's smartpca
#' @param clin_df the clinical DF for the cohort of interest
add_patient_info <- function(pcs, clin_df) {
  
  # Add three additional columns to the PCA DF
  new_rows <- lapply(rownames(pcs), function(id) {
    
    # Get the age, race, and gender info for this patient # TODO: add subtype as well?
    patient_row <- clin_df[clin_df$case_submitter_id == id, ]
    
    age <-  unique(patient_row[,'age_at_index'])
    race <- unique(patient_row[, 'race'])
    gender <- unique(patient_row[ , 'gender']) 
    
    return(c(age, race, gender))
  })

  new_df <- do.call(rbind, new_rows)
  pcs <- cbind(pcs, new_df)
  
  return(pcs)
}

# Call function
pcs_new <- add_patient_info(pcs, brca_clinical_df)
pcs_new[,'race'] <- as.character(pcs_new[,'race'])
pcs_new[,'age_at_index'] <- as.integer(pcs_new[,'age_at_index'])

# Write to file
fwrite(pcs_new, "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/smartpca/BRCA_output_withClinInfo.csv")

# Make a scatter plot of the first 2 PCs colored by clinical characteristic
ggplot(pcs_new, aes(x = V2, y = V3, color = age_at_index)) + geom_point() + 
  xlab("PC1") + ylab("PC2") 

pcs_new$race <- unlist(lapply(pcs_new$race, function(r) {
    if(r == "" | r == "character(0)") {return("not reported")}
    else {return(r)}
  }))
ggplot(pcs_new, aes(x = V2, y = V3, color = race)) + geom_point() + 
  xlab("PC1") + ylab("PC2") 


# Calculate the Pearson correlation between variables and plot
# NOTE: can also use Spearman or Kendall, just input "spearman" or "kendall" into 'method' param)

# Add a numerical representation of race for association analysis
pcs_new$race_numerical <- unlist(lapply(pcs_new$race, function(r) {
  if(r == "white") {return(0)}
  #else if (r == "black or african american") {return(1)}
  #else if (r == "asian") {return(2)}
  #else {return(3)}
  else {return(1)}
}))

#pcs_new$race <- unlist(lapply(pcs_new$race, function(r) if(r == "") 
  #{return("not reported")} else {return(r)}))

# Check for normality of variables
shapiro.test(pcs_new$age_at_index)
shapiro.test(pcs_new$race_numerical)
shapiro.test(pcs_new$V2)
shapiro.test(pcs_new$V3)

# Calculate correlations
cor_age_pc1 <- cor.test(pcs_new$V2, pcs_new$age_at_index, method = "spearman")
cor_race_pc1 <- cor.test(pcs_new$V2, pcs_new$race_numerical, method = "spearman")
cor_age_pc2 <- cor.test(pcs_new$V3, pcs_new$age_at_index, method = "spearman")
cor_race_pc2 <- cor.test(pcs_new$V3, pcs_new$race_numerical, method = "spearman")

# Plot correlations
ggscatter(pcs_new, x = "age_at_index", y = "V2", add = "reg.line", conf.int = TRUE,
          cor.coeff = TRUE, cor.method = "spearman", xlab = "Age", ylab = "PC1")
ggscatter(pcs_new, x = "age_at_index", y = "V3", add = "reg.line", conf.int = TRUE,
          cor.coeff = TRUE, cor.method = "spearman", xlab = "Age", ylab = "PC2")

# For race, do an ANOVA and also pairwise comparisons
#ggscatter(pcs_new, x = "race_numerical", y = "V2", add = "reg.line", conf.int = TRUE,
          #cor.coeff = TRUE, cor.method = "spearman", xlab = "Race", ylab = "PC1")
ggboxplot(pcs_new, x = "race", y = "V2", palette = "jco", ylab = "PC1") + stat_compare_means(method = "anova")
compare_means(V2 ~ race, data = pcs_new)
my_comparisons <- list(c("white", "black or african american"), c("white", "not reported"),
                       c("asian", "black or african american"), c("asian", "not reported"),
                       c("black or african american", "not reported"))
ggboxplot(pcs_new, x = "race", y = "V2", palette = "jco", ylab = "PC1") + 
  stat_compare_means(comparisons = my_comparisons) #+ stat_compare_means(label.y = 50)


#ggscatter(pcs_new, x = "race_numerical", y = "V3", add = "reg.line", conf.int = TRUE,
          #cor.coeff = TRUE, cor.method = "spearman", xlab = "Race", ylab = "PC2")
ggboxplot(pcs_new, x = "race", y = "V3", palette = "jco", ylab = "PC2") + stat_compare_means(method = "anova")
compare_means(V3 ~ race, data = pcs_new)
my_comparisons <- list(c("white", "black or african american"), c("asian", "black or african american"), 
                        c("black or african american", "not reported"))
ggboxplot(pcs_new, x = "race", y = "V3", palette = "jco", ylab = "PC2") + 
  stat_compare_means(comparisons = my_comparisons) #+ stat_compare_means(label.y = 50)



############################################################
### PERFORM PCA USING AFFYMETRIX SNP 6.0 GENOTYPE FILES
#############################################################

# Modeled off of code from: https://bioinfo.uth.edu/LungCancerSubtypes/code_for_pca_with_TCGA_data.txt
# Gives detailed R code for performing PCA with TCGA genotype files
# To run R in this cluster, first run: module add R

# Additional useful PLINK vignettes: 
# Processing 1000 Genomes reference data for ancestry estimation: https://meyer-lab-cshl.github.io/plinkQC/articles/Genomes1000.html
# Ancestry estimation based on reference samples of known ethnicities: https://meyer-lab-cshl.github.io/plinkQC/articles/AncestryCheck.html
# Tutorial:Produce PCA bi-plot for 1000 Genomes Phase III - Version 2: https://www.biostars.org/p/335605/

#path <- "/home/scamilli/Breast_Cancer/SNP 6.0 Genotype Arrays/"
path <- "/home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/"

library(reshape2)

## 1. FORMAT GENOTYPE FILES FROM THE TCGA LEGACY ARCHIVE, EXTRACT NORMAL BLOOD SAMPLES ##

# Load all the sample SNP 6.0 files
all.brca <- list.files(path, recursive = TRUE, pattern = "*.birdseed.data.txt$")
# Or, get all TCGA and subset by the names in the SDRF file

# Load the SDRF file that contains the naming information for all the samples
brca.names <- read.delim(paste0(path, "broad.mit.edu_BRCA.Genome_Wide_SNP_6.sdrf.txt"), header = TRUE)
#all.brca <- all.tcga[all.tcga %in% brca.names$Derived.Array.Data.Matrix.File.1]

# Split the file names so that they match the names in the SDRF file
#all.brca.split <- colsplit(all.brca, "/", c("folder_name", "sample_name"))
brca.final.matched <- brca.names[brca.names$Derived.Array.Data.Matrix.File.1 %in% all.brca, c(1,2,31)]
brca.split.barcode <- colsplit(brca.final.matched$Comment..TCGA.Barcode., "-", 
                               c("project", "TSS", "participant", "sample", "portion", "plate", "center"))
brca.final.matched$patient <- unlist(lapply(brca.final.matched$Comment..TCGA.Barcode., function(x) 
  unlist(strsplit(x, "-", fixed = TRUE))[3]))

# Extract the germline blood only files (is this necessary?)
brca.germline.blood <- brca.final.matched[grepl("10", brca.split.barcode$sample),]
# Only 244 BRCA blood samples; 272 total BRCA normal samples (including matched tissue normal)

# To get duplicates (those with normal blood and a tumor):
brca.germline.blood <- unique(brca.split.barcode[grepl("10", brca.split.barcode$sample), 'participant'])
brca.germline.all.norm <- unique(brca.split.barcode[grepl("10", brca.split.barcode$sample) | grepl("11", brca.split.barcode$sample), 'participant'])
brca.tumor <- unique(brca.split.barcode[grepl("01", brca.split.barcode$sample, fixed = TRUE), 'participant'])

brca.intersecting.ids <- read.table("/home/scamilli/Breast_Cancer/intersecting_ids.txt", header = TRUE, row.names = 1)[,1]
#length(intersect(brca.duplicates, brca.intersecting.ids))

brca.germline.blood.intersect <- intersect(brca.germline.blood, brca.intersecting.ids)
brca.germline.all.intersect <- intersect(brca.germline.all.norm, brca.intersecting.ids)
brca.tumor.intersect <- intersect(brca.tumor, brca.intersecting.ids)

length(intersect(brca.germline.all.intersect, brca.tumor.intersect))
# only 33 patients with one of each

brca.final.matched.germline.blood <- brca.final.matched[brca.final.matched$patient %in% brca.germline.blood.intersect,]
brca.final.matched.germline.all <- brca.final.matched[brca.final.matched$patient %in% brca.germline.all.intersect,]
brca.final.matched.tumor <- brca.final.matched[brca.final.matched$patient %in% brca.tumor.intersect,]
#print(brca.final.matched.duplicates[brca.final.matched.duplicates$patient == "A3TN", c('Derived.Array.Data.Matrix.File.1', 'Comment..TCGA.Barcode.')])

#write.table(intersect(unique(brca.duplicates), brca.intersecting.ids), paste0(path, "brca.duplicates.txt"))
#write.table(intersect(unique(brca.duplicates.incl.tissue.norm), brca.intersecting.ids), paste0(path, "brca.duplicates.incl.tissue.norm.txt"))
#write.table(intersect(unique(brca.tumor), brca.intersecting.ids), paste0(path, "brca.tumor.txt"))

# Get the pairs of file names for each ID with both tumor and normal
pairs.all <- lapply(brca.germline.all.intersect, function(x) 
  unlist(brca.final.matched[brca.final.matched$patient == x, 'Derived.Array.Data.Matrix.File.1']))
names(pairs.all) <- brca.germline.all.intersect
# Get the # of entries with two files
to_keep <- unlist(lapply(pairs.all, function(x) ifelse(length(x) > 1, TRUE, FALSE)))
pairs.all.two.files <- pairs.all[to_keep]
pairs.all.two.files <- lapply(pairs.all.two.files, function(x) as.data.frame(x))
names(pairs.all.two.files) <- names(pairs.all)[to_keep]

# Remove any with 3
pairs.all.two.files.sub <- lapply(pairs.all.two.files, function(x) {
  if(nrow(x) > 2) {
    # Get the starting char for all three
    start_chars <- unlist(lapply(x[,1], function(k) 
      unlist(strsplit(k, "_", fixed = TRUE))[1]))
    print(start_chars)
    # Get the most common start char and keep these files
    u_start_char <- unique(start_chars)
    dominant_start_char <- u_start_char[which.max(tabulate(match(start_chars, u_start_char)))]
    print(dominant_start_char)
    x_sub <- x[x[,1] == dominant_start_char,]
    return(x_sub)
  }
  else {return(x)}
})
pairs.all.two.files.sub <- pairs.all.two.files.sub[unlist(lapply(pairs.all.two.files.sub, function(x) 
  ifelse(length(x) > 0, TRUE, FALSE)))]
pairs.all.two.files.df <- do.call(cbind, pairs.all.two.files.sub)
colnames(pairs.all.two.files.df) <- names(pairs.all.two.files.sub)
pairs.all.two.files.df.t <- t(pairs.all.two.files.df)
colnames(pairs.all.two.files.df.t) <- c("File1", "File2")
write.csv(pairs.all.two.files.df.t, paste0(path, "brca.tum.norm.matched.files.csv"))

#' Given two filenames, checks if the lines (without confidence for SNP calling) are
#' identical. Returns the number of mismatching lines and the percentage of the total.
#' @param f1 the first file (either tumor or normal)
#' @param f2 the second file (either tumor or normal)
check_percentage_identical <- function(f1, f2) {

  file1 <- read.table(f1, header = FALSE, skip = 2)
  file2 <- read.table(f2, header = FALSE, skip = 2)
  
  file1 <- file1[,1:2]
  file2 <- file2[,1:2]
  
  colnames(file1) <- c("SNP", "Call")
  colnames(file2) <- c("SNP", "Call")
  
  file3 <- merge(file1, file2, by = "SNP")
  print(head(file3))
  print(dim(file3))
  
  n_differing_rows <- nrow(file3[file3$Call.x != file3$Call.y,])
  print(paste("Num. Differing Rows:", n_differing_rows))
  
  tot_r <- nrow(file1)
  perc_diff <- (n_differing_rows / tot_r) * 100
  print(paste("Percentage Difference", perc_diff))
  
  return(perc_diff)
}

percentage_differences <- unlist(lapply(1:nrow(brca.tum.norm.matched), function(x) 
  check_percentage_identical(brca.tum.norm.matched[x,1], brca.tum.norm.matched[x,2])))

# Plot a histogram of the percentage differences
pdf("percentage.diff.brca.pdf", width = 450, height = 450)
plot(hist(percentage_differences, xlab = "Percentage Discrepancy in SNP Calls (T vs. N)", ylab = "Frequency"))
dev.off()

print(range(percentage_differences))

percentage_differences_df <- data.frame("Patient" = rownames(brca.tum.norm.matched),
                                        "Perc.Discrepancy" = percentage_differences)

# Add the tumor mutational burden for each of these patients
# Import the clinical DF with total mutation count
brca.clinical.df <- read.csv(paste0(path, "brca_clinical_data_subset_w_Nonsyn_MutCounts.csv"),
                             header = TRUE, check.names = FALSE, row.names = 1)
brca.clinical.df$patient <- unlist(lapply(brca.clinical.df$tcga_barcode, function(x)
  unlist(strsplit(x, "-", fixed = TRUE))[3]))
percentage_differences_df$Tot.Mut.Count <- unlist(lapply(percentage_differences_df$Patient, function(p) {
  return(unique(brca.clinical.df[brca.clinical.df$patient == p, 'Total.Num.Muts']))
}))

# Make a scatter plot of the percentage discrepancy vs. the total nonsynonymous somatic
# mutational burden to see if there is a correlation
lm_res = lm(Tot.Mut.Count ~ Perc.Discrepancy, data = percentage_differences_df)
r2 = round(summary(lm_res)$adj.r.squared, digits = 4)

pdf("percentage.diff.vs.tmb.brca.pdf", width = 4, height = 4)
plot(x = as.numeric(percentage_differences_df$Perc.Discrepancy), 
     y = as.numeric(percentage_differences_df$Tot.Mut.Count), 
     xlab = "Percentage Discrepancy in SNP Calls (T vs. N)", 
     ylab = "Total # Nonsynonymous Somatic Mutations", pch = 19)
#abline(a = 0, b = 10, col = "red")
abline(lm_res, col = "red")
text(x = 2, y = 60, labels = paste("R2:", r2))
dev.off()


# (In Linux) diff -q KEBAB_p_TCGASNP_226_227_N_GenomeWideSNP_6_H07_1151588.birdseed.data.txt KEBAB_p_TCGASNP_226_227_N_GenomeWideSNP_6_E02_1151484.birdseed.data.txt


# Load the clinical file from the TCGA
brca.clin <- read.delim('nationwidechildrens.org_clinical_patient_brca.txt', skip = 1)
brca.clin.to.use <- brca.clin[-1,] # remove the first row for file matching

# Write this to a new file
write.table(brca.germline.blood, file = paste0(path, "blood_normal_BRCA.tsv"), col.names = TRUE,
            row.names = FALSE, quote = FALSE, sep = "\t")


## 2. EXTRACT FILE NAMES AND FORMAT LOW-CONFIDENCE CALLS USING THE COMMAND LINE ##

# Extract just the third column from the files we just created (sample) using command line
# cut -f3 blood_normal_BRCA.tsv > blood_normal_copy_BRCA.tsv 
# mkdir BRCA_TCGA_normal_blood

# Put the following in a script called: BRCA_blood_genotype_copy.sh
# This will move all the BRCA birdseed files to this subdirectory
#!/bin/bash

#while read line
#do
#   find BRCA -iname "$line" -exec cp '{}' BRCA_TCGA_normal_blood \;
#done < BRCA/blood_normal_copy_BRCA.tsv
#chmod u+x BRCA_blood_genotype_copy.sh 

# Then, run the script using SLURM or bash
#sbatch BRCA_blood_genotype_copy.sh
#bash BRCA_blood_genotype_copy.sh

## 3. REPLACE BIRDSEED CALLS WITH CONFIDENCE < 0.1 THRESHOLD WITH -9 ##
#for f in BRCA_TCGA_normal_blood/*.-9_conversion.txt
#   do sed 's/ //g' $f | awk '{if ($3>0.1){$2=-9}}{print $0}' | sed 's/ /\t/g' > $f.-9_conversion.txt
#done


## 4. COMBINE ALL GENOTYPE CALLS TOGETHER PER CANCER TYPE ##

# Put the following in a script called: combine_reformatted_tcga_blood_genotype_BRCA.sh
#!/bin/bash

# Create temporary directories to hold results
#tmp=$(mktemp)
#tmp2=$(mktemp)
#tmp3=$(mktemp)

#echo $tmp;

#for file in BRCA_TCGA_normal_blood/*.birdseed.data.txt
#do
#   if [ -s "$tmp" ]
#   then
#       cut -f 1,2 "$file" > "$tmp3"
#       join "$tmp" "$tmp3" > "$tmp2"
#   else
#       cut -f 1,2 "$file" > "$tmp3"
#       cp "$tmp3" "$tmp2"
#   fi
#   cp "$tmp2" "$tmp"
#done

#chmod u+x combine_reformatted_tcga_blood_genotype_BRCA.sh

# Then, run the script on SLURM
# sbatch combine_reformatted_tcga_blood_genotype_BRCA.sh
# bash combine_reformatted_tcga_blood_genotype_BRCA.sh

# Copy the tmp files to a new location
# Note the location of the tmp files (e.g. /tmp/tmp.60ypaUcH1Q) and use that in the following command
#cp -i /tmp/tmp.60ypaUcH1Q ./reformatted_files/BRCA_combined_reformatted_genotype_files.txt

# Ours: /tmp/tmp.JZSySVbhch


## 5. USE R TO CONVERT COMBINED GENOTYPE FILES TO PLINK ALLELE CODED FORMAT OF A,C,G,T and 0 0 for missing calls ##
# module add R

brca.combined <- read.table(paste0(path, "reformatted_files/BRCA_combined_reformatted_genotype_files.txt"))
brca.combined <- brca.combined[-c(1,2),] # remove the top two header lines

# or brca.combined <- read.table(paste0(path, "BRCA_files_for_plink/reformatted_files/BRCA_combined_reformatted_genotype_files2.txt"), skip = 2)

# Download the Affymetrix mapping file from http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping
# NOTE - LINK IS DOWN; downloaded from here: https://www.thermofisher.com/us/en/home/life-science/microarray-analysis/microarray-data-analysis/genechip-array-annotation-files.html
# Import the mapping file and skip the header info
affymap <- read.csv("/home/scamilli/GenomeWideSNP_6.na35.annot.csv", header = TRUE, skip = 18)

# Merge together some of the attributes from the TCGA data and the affymap in order to create a map file for PLINK
brca.combined.merge <- merge(affymap[, c(1:5)], brca.combined, by.x = c("Probe.Set.ID"), by.y = c("V1")) 
brca.map.file.for.plink <- brca.combined.merge[ , c(3,2,4)] 

# Extract sample names to create the ped file for PLINK
brca.sample.names <- as.data.frame(paste("BRCA_sample", seq(1:ncol(brca.combined[ ,-1])), sep = "")) 
# Another way:
#brca.sample.names <- read.table(file = paste0(path, "reformatted_files/BRCA_combined_reformatted_genotype_files2.txt"), header = F, nrows = 1)
#brca.sample.names <- as.character(unlist(brca.sample.names))
#brca.sample.names <- brca.sample.names[grepl("GenomeWideSNP_6", brca.sample.names)]
#brca.sample.names <- as.data.frame(brca.sample.names)

# Remove probe IDs without SNP positions; make the new map and genotype files for PLINK
pos.to.remove <- which(brca.map.file.for.plink$Physical.Position == "---")
brca.map.file.new <- brca.map.file.for.plink[-pos.to.remove, ]

brca.combined.merge.2 <- merge(affymap[, c(1:5, 9, 10)], brca.combined, 
                               by.x = c("Probe.Set.ID"), by.y = c("V1"))

# Remove bad probes without SNP positions
brca.combined.merge.3 <- brca.combined.merge.2[-pos.to.remove, ]
rm(brca.combined.merge.2)

# Recode the 0,1,2 alleles from Birdseed to the actual bases for PLINK
brca.ped <- ifelse(brca.combined.merge.3 == 0, paste(brca.combined.merge.3$Allele.A, brca.combined.merge.3$Allele.A, sep = " "),
                   ifelse(brca.combined.merge.3 == 1, paste(brca.combined.merge.3$Allele.A, brca.combined.merge.3$Allele.B, sep = " "),
                          ifelse(brca.combined.merge.3 == 2, paste(brca.combined.merge.3$Allele.B, brca.combined.merge.3$Allele.B, sep = " "),
                                 paste("0", "0", sep = " "))))

# Remove the non-sample stuff from the new ped file
brca.ped.for.plink <- brca.ped[, -c(1:7)]

# ANOTHER WAY
#brca.ped <- apply(brca.combined.merge.3[,8:ncol(brca.combined.merge.3)], MARGIN = c(1,2), function(x) {
#  if(x == 0) {paste(brca.combined.merge.3$Allele.A, brca.combined.merge.3$Allele.A, sep = " ")}
#  else if (x == 1) {paste(brca.combined.merge.3$Allele.A, brca.combined.merge.3$Allele.B, sep = " ")}
#  else if (x == 2) {paste(brca.combined.merge.3$Allele.B, brca.combined.merge.3$Allele.B, sep = " ")}
#  else {paste("0", "0", sep = " ")}
#})

# Transpose the file for use with PLINK
brca.ped.for.plink.t <- t(brca.ped.for.plink)

# Write the new files for use with PLINK
write.table(brca.ped.for.plink.t, file = "BRCA_genotype_file_for_plink.ped", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(brca.map.file.new, file = "BRCA_genotype_file_for_plink.map", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(brca.sample.names, file = "BRCA_sample_names.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


## 6. COMBINE THE NEWLY CODED GENOTYPE FILES TOGETHER USING PLINK, IF RUNNING MORE THAN ONE CANCER TYPE ##
# Command line in Linux, with PLINK software installed

# Paste the sample names to the ped files; change the name to match
#paste -d" " BRCA_sample_names.txt BRCA_genotype_file_for_plink.ped > BRCA_genotype_file_for_plink_with_sample_names.ped
#cp BRCA_genotype_file_for_plink.map BRCA_genotype_file_for_plink_with_sample_names.map 

# Make the bfiles from the map and PED files using PLINK
# The --no-fid --no-parents, etc. flags mark that these .fam and .ped files lack family ID, parental ID, sex, and phenotype columns
# Added --map3 option to denote that there are only 3 columns
#./plink --file BRCA_genotype_file_for_plink_with_sample_names --no-fid --no-parents --no-sex --no-pheno --make-bed --out BRCA_genotype_file_for_plink_with_sample_names
# FOR ME SPECIFICALLY: ./plink-1.07-x86_64/plink --noweb --file BRCA_genotype_file_for_plink_with_sample_names --no-fid --no-parents --no-sex --no-pheno --make-bed --map3 --out BRCA_genotype_file_for_plink_with_sample_names
# Can also add the --maf 0.05 to filter those SNPs below a given minor allele frequency threshold
#./plink-1.07-x86_64/plink --noweb --file BRCA_genotype_file_for_plink_with_sample_names --no-fid --no-parents --no-sex --no-pheno --make-bed --map3 --out BRCA_genotype_file_for_plink_with_sample_names


# OPTIONAL: If we are doing this separately for different cancer types in the TCGA, use the mergelist function to 
# combine them all
# Create a file called tcga_genotype_files_for_merging_in_plink.txt containing the following:
#BRCA_genotype_file_for_plink_with_sample_names
#LUAD_genotype_file_for_plink_with_sample_names
#LUSC_genotype_file_for_plink_with_sample_names
#HNSC_genotype_file_for_plink_with_sample_names
#BLCA_genotype_file_for_plink_with_sample_names
#GBM_genotype_file_for_plink_with_sample_names
#LGG__genotype_file_for_plink_with_sample_names
# ... (include all necessary files)

# Run mergelist
#./plink --merge-list tcga_genotype_files_for_merging_in_plink.txt --make-bed --out all_tcga_genotypes_merged

# Remove potential triallelic SNPs using missnp
#./plink --file BRCA_genotype_file_for_plink_with_sample_names --no-fid --no-parents --no-sex --no-pheno --exclude all_tcga_genotypes_merged-merge.missnp --make-bed --out BRCA_genotype_file_for_plink_with_sample_names_tri_allelic_removed
# <Repeat this command for all cancer types of interest, then run mergelist again >

# Contents of tcga_genotype_files_for_merging_in_plink_tri_allelic_removed.txt file:
#BRCA_genotype_file_for_plink_with_sample_names_tri_allelic_removed
#LUAD_genotype_file_for_plink_with_sample_names_tri_allelic_removed
#LUSC_genotype_file_for_plink_with_sample_names_tri_allelic_removed
#HNSC_genotype_file_for_plink_with_sample_names_tri_allelic_remove
#BLCA_genotype_file_for_plink_with_sample_names_tri_allelic_removed
#GBM_genotype_file_for_plink_with_sample_names_tri_allelic_removed
#LGG_genotype_file_for_plink_with_sample_names_tri_allelic_removed
# ... (include all necessary files)

# Run mergelist
#./plink --merge-list tcga_genotype_files_for_merging_in_plink_tri_allelic_removed.txt --make-bed --out all_tcga_genotypes_merged_tri_allelic_removed


## 7. TRIM THE SNPs IN LD PRIOR TO PCA ##
# Use an R2 value of 0.5

# Make a set of SNPs for the LD trim from the merged file, using PLINK's write-snplist command
#./plink-1.07-x86_64/plink --noweb --bfile all_tcga_genotypes_merged_tri_allelic_removed --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.1  --write-snplist --out all_tcga_genotypes_merged_tri_allelic_removed_for_ld_trim
# Actual: ./plink-1.07-x86_64/plink --noweb --bfile /home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/BRCA_files_for_plink/BRCA_genotype_file_for_plink_with_sample_names --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.1 --write-snplist --out /home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/BRCA_files_for_plink/BRCA_genotype_file_for_plink_with_sample_names_for_ld_trim

# Run the LD trim on these SNPs using PLINK formatted 1000 Genomes European file from https://vegas2.qimrberghofer.edu.au/g1000p3_EUR.tar.gz
#./plink-1.07-x86_64/plink --noweb --bfile g1000p3_EUR --indep-pairwise 1000 1 0.5 --extract all_tcga_genotypes_merged_tri_allelic_removed_for_ld_trim.snplist --out all_tcga_genotypes_merged_tri_allelic_removed_ld_snps_pruned
# Actual: ./plink-1.07-x86_64/plink --noweb --bfile /home/scamilli/g1000p3_EUR/g1000p3_EUR --indep-pairwise 1000 1 0.5 --extract /home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/BRCA_files_for_plink/BRCA_genotype_file_for_plink_with_sample_names_for_ld_trim.snplist --out /home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/BRCA_files_for_plink/BRCA_genotype_file_for_plink_with_sample_names_ld_snps_pruned
# --indep-pairwise: window size (in kb), step size, and r^2 threshold

# NOTE: This takes a really long time! Run in the background with nohup, since there is no slurm job management system on gen-singhtmp cluster
# i.e. https://www.serverwatch.com/guides/detach-processes-with-disown-and-nohup/
#nohup ./plink-1.07-x86_64/plink --noweb --bfile /home/scamilli/g1000p3_EUR/g1000p3_EUR --indep-pairwise 1000 1 0.5 --extract /home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/BRCA_files_for_plink/BRCA_genotype_file_for_plink_with_sample_names_for_ld_trim.snplist --out /home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/BRCA_files_for_plink/BRCA_genotype_file_for_plink_with_sample_names_ld_snps_pruned


# UPDATE (9-8-22): We might not actually need to use a reference file; also, LD-pruning using Plink 1.9 is 
# much faster. New command:
#nohup ./plink-1.9-x86_64/plink --noweb --bfile /home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/BRCA_files_for_plink/BRCA_genotype_file_for_plink_with_sample_names --keep /home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/BRCA_files_for_plink/BRCA_genotype_file_for_plink_with_sample_names.fam --extract /home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/BRCA_files_for_plink/BRCA_genotype_file_for_plink_with_sample_names_for_ld_trim.snplist --indep-pairwise 50 5 0.5 --out /home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/BRCA_files_for_plink/BRCA_genotype_file_for_plink_with_sample_names_ld_snps_pruned_no_ref


## 10. RUN PCA USING THE PRUNE.IN SNPS, REMOVING OUTLIERS ##
# NOTE: Plink version 1.07 does not have a built-in PCA function; need to use version 1.9 for this
#./plink --bfile all_tcga_genotypes_merged_tri_allelic_removed --extract all_tcga_genotypes_merged_tri_allelic_removed_ld_snps_pruned.prune.in --pca --out pca_all_tcga_genotypes_merged_no_tri_using_ld_pruned_snps
# Actual: ./plink-1.9-x86_64/plink --noweb --bfile /home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/BRCA_files_for_plink/BRCA_genotype_file_for_plink_with_sample_names --extract /home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/BRCA_files_for_plink/BRCA_genotype_file_for_plink_with_sample_names_ld_snps_pruned_no_ref.prune.in --pca --out pca_brca_genotypes_using_ld_pruned_snps

# Look at the plot and create a file with a list of the samples that are outliers and need removal, e.g. brca_samples_to_exclude_pca.txt
# Then, re-run the PCA without the outliers
#./plink-1.9-x86_64/plink --noweb -bfile all_tcga_genotypes_mergedtri_allelic_removed --remove brca_samples_to_exclude_pca.txt --extract all_tcga_genotypes_merged_white_only_tri_allelic_removed_ld_snps_pruned.prune.in --pca --out pca_all_tcga_genotypes_merged_no_tri_using_ld_pruned_snps_brca_outliers_excluded

# VIEW THE PCA RESULTS!



# NOTE: If we decide to use a 1000 Genomes reference file other than the European one, we can create it like this:
# Adapted from 'Processing 1000 Genomes reference data for ancestry estimation': https://meyer-lab-cshl.github.io/plinkQC/articles/Genomes1000.html
# https://dougspeed.com/1000-genomes-project/

#mkdir 1000Grefdir
#cd 1000Grefdir/

#pgen=https://www.dropbox.com/s/e5n8yr4n7y91fyp/all_hg38.pgen.zst?dl=1
#pvar=https://www.dropbox.com/s/cy46f1c8yutd1h4/all_hg38.pvar.zst?dl=1
#sample=https://www.dropbox.com/s/3j9zg103fi8cjfs/hg38_corrected.psam?dl=1

#wget $pgen
#mv 'all_hg38.pgen.zst?dl=1' all_hg38.pgen.zst
#./plink-1.07-x86_64/plink --noweb --zst-decompress all_hg38.pgen.zst > all_hg38.pgen
  # this failed for me: used this alternative: unzstd /home/scamilli/1000Grefdir/all_hg38.pgen.zst

#wget $pvar
#mv 'all_hg38.pvar.zst?dl=1' all_hg38.pvar.zst
# unzstd /home/scamilli/1000Grefdir/all_hg38.pvar.zst

#wget $sample
#mv 'hg38_corrected.psam?dl=1' all_hg38.psam

#./plink-1.07-x86_64/plink --noweb --pfile $1000Grefdir/all_hg38 --max-alleles 2 --make-bed --out $1000Grefdir/all_hg38

# Different option to restrict to autosomal SNPs with MAF>0.01 (and excluding duplicates and SNPs with name ".")
#echo "." > exclude.snps
#./plink-1.07-x86_64/plink --noweb --make-bed --out /home/scamilli/1000Grefdir/all_hg38 --pgen /home/scamilli/1000Grefdir/all_hg38.pgen --pvar /home/scamilli/1000Grefdir/all_hg38.pvar --psam /home/scamilli/1000Grefdir/all_hg38.psam --maf 0.01 --autosome --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all --exclude /home/scamilli/1000Grefdir/exclude.snps


# OR: like this: 
# Adapted from https://www.biostars.org/p/335605/
# Link to 1000 Genomes data file: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/


############################################################
### VISUALIZE PCA RESULTS
############################################################
local_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/smartpca or plink/Non-Tumor-Normal Matched/"
pcs <- read.table(paste0(local_path, "pca_brca_genotypes_using_ld_pruned_snps_no_ref.eigenvec"),
                  header = FALSE, sep = " ")

# Get the TCGA barcodes for the patient and replace the rownames with these
barcode_mapping <- read.csv(paste0(local_path, "blood_normal_BRCA.tsv"), header = TRUE, check.names = FALSE,
                            sep = "\t")
pcs$tcga_barcode <- unlist(lapply(pcs[,1], function(x) barcode_mapping[grepl(x, barcode_mapping$Derived.Array.Data.Matrix.File.1), 
                                                                       'Comment..TCGA.Barcode.']))

# If needed
rownames(pcs) <- unlist(lapply(pcs$tcga_barcode, function(x) paste(unlist(strsplit(x, "-", fixed = TRUE))[1:3], collapse = "-")))
pcs <- pcs[,3:ncol(pcs)]
pcs <- pcs[,-ncol(pcs)]
colnames(pcs) <- paste0("PC", 1:ncol(pcs))

# Make a scatter plot of the first 2 PCs
ggplot(pcs, aes(x = PC1, y = PC2)) + geom_point() + xlab("PC1") + ylab("PC2") 

# Identify the names of outliers
ggplot(pcs, aes(x = PC1, y = PC2)) + geom_point() + 
  geom_text(aes(label= ifelse(PC2 > quantile(PC2, 0.95), as.character(PC1),'')), hjust=1, vjust=0)


# Overlay clinical information
# Import clinical DFs 
#brca_clinical_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/clinical.csv",
                         # header = TRUE)
brca_clinical_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/clinical_data_subset.csv",
                          header = TRUE)
#pc_clinical_df <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/clinical_data_subset.csv",
#header = TRUE)

# Subset the patients to those in the clinical DF
pcs_sub <- pcs[rownames(pcs) %in% brca_clinical_df$tcga_barcode,]

#' Extract information about age, gender, and race and add as rows for each patient
#' to the PCs DF
#' @param pcs the pca output from eigensoft's smartpca
#' @param clin_df the clinical DF for the cohort of interest
add_patient_info <- function(pcs, clin_df) {
  
  # Add three additional columns to the PCA DF
  new_rows <- lapply(rownames(pcs), function(id) {
    
    # Get the age, race, and gender info for this patient # TODO: add subtype as well?
    patient_row <- clin_df[clin_df$tcga_barcode == id, ]
    
    age <-  unique(patient_row[,'age_at_index'])
    race <- unique(patient_row[, 'race'])
    gender <- unique(patient_row[ , 'gender']) 
    
    return(c(age, race, gender))
  })
  
  new_df <- do.call(rbind, new_rows)
  pcs <- cbind(pcs, new_df)
  
  return(pcs)
}

# Call function
pcs_new <- add_patient_info(pcs_sub, brca_clinical_df)
pcs_new[,'race'] <- as.character(pcs_new[,'race'])
pcs_new[,'age_at_index'] <- as.integer(pcs_new[,'age_at_index'])

# Write to file
fwrite(pcs_new, paste0(local_path, "pca_brca_genotypes_using_ld_pruned_snps_no_ref_with_clinical_info.eigenvect.csv"))

# Make a scatter plot of the first 2 PCs colored by clinical characteristic
ggplot(pcs_new, aes(x = PC1, y = PC2, color = age_at_index)) + geom_point() + 
  xlab("PC1") + ylab("PC2") 

pcs_new$race <- unlist(lapply(pcs_new$race, function(r) {
  if(r == "" | r == "character(0)") {return("not reported")}
  else {return(r)}
}))
ggplot(pcs_new, aes(x = PC1, y = PC2, color = race)) + geom_point() + 
  xlab("PC1") + ylab("PC2") 


# Calculate the Pearson correlation between variables and plot
# NOTE: can also use Spearman or Kendall, just input "spearman" or "kendall" into 'method' param)

# Add a numerical representation of race for association analysis
pcs_new$race_numerical <- unlist(lapply(pcs_new$race, function(r) {
  if(r == "white") {return(0)}
  #else if (r == "black or african american") {return(1)}
  #else if (r == "asian") {return(2)}
  #else {return(3)}
  else {return(1)}
}))

#pcs_new$race <- unlist(lapply(pcs_new$race, function(r) if(r == "") 
#{return("not reported")} else {return(r)}))

# Check for normality of variables
shapiro.test(pcs_new$age_at_index)
shapiro.test(pcs_new$race_numerical)
shapiro.test(pcs_new$PC1)
shapiro.test(pcs_new$PC2)

# Calculate correlations
cor_age_pc1 <- cor.test(pcs_new$PC1, pcs_new$age_at_index, method = "spearman")
cor_race_pc1 <- cor.test(pcs_new$PC1, pcs_new$race_numerical, method = "spearman")
cor_age_pc2 <- cor.test(pcs_new$PC2, pcs_new$age_at_index, method = "spearman")
cor_race_pc2 <- cor.test(pcs_new$PC2, pcs_new$race_numerical, method = "spearman")

# Plot correlations
ggscatter(pcs_new, x = "age_at_index", y = "PC1", add = "reg.line", conf.int = TRUE,
          cor.coeff = TRUE, cor.method = "spearman", xlab = "Age", ylab = "PC1")
ggscatter(pcs_new, x = "age_at_index", y = "PC2", add = "reg.line", conf.int = TRUE,
          cor.coeff = TRUE, cor.method = "spearman", xlab = "Age", ylab = "PC2")

# For race, do an ANOVA and also pairwise comparisons
#ggscatter(pcs_new, x = "race_numerical", y = "V2", add = "reg.line", conf.int = TRUE,
#cor.coeff = TRUE, cor.method = "spearman", xlab = "Race", ylab = "PC1")
ggboxplot(pcs_new, x = "race", y = "PC1", palette = "jco", ylab = "PC1") + stat_compare_means(method = "anova")
compare_means(PC1 ~ race, data = pcs_new)
my_comparisons <- list(c("white", "black or african american"), c("white", "not reported"),
                       c("asian", "black or african american"), c("asian", "not reported"),
                       c("black or african american", "not reported"))
ggboxplot(pcs_new, x = "race", y = "PC1", palette = "jco", ylab = "PC1") + 
  stat_compare_means(comparisons = my_comparisons) #+ stat_compare_means(label.y = 50)


#ggscatter(pcs_new, x = "race_numerical", y = "V3", add = "reg.line", conf.int = TRUE,
#cor.coeff = TRUE, cor.method = "spearman", xlab = "Race", ylab = "PC2")
ggboxplot(pcs_new, x = "race", y = "PC2", palette = "jco", ylab = "PC2") + stat_compare_means(method = "anova")
compare_means(PC2 ~ race, data = pcs_new)
my_comparisons <- list(c("white", "black or african american"), c("asian", "black or african american"), 
                       c("black or african american", "not reported"))
ggboxplot(pcs_new, x = "race", y = "PC2", palette = "jco", ylab = "PC2") + 
  stat_compare_means(comparisons = my_comparisons) #+ stat_compare_means(label.y = 50)



# Create a scree plot to visualize the variance explained by each PC
eigenvals <- read.table(paste0(local_path, "pca_brca_genotypes_using_ld_pruned_snps_no_ref.eigenval"),
                        header = FALSE, sep = " ")[,1]
plot(x = seq(1:length(eigenvals)), y = as.numeric(eigenvals), type = "o", xlab = "Principal Component", ylab = "Variance")

