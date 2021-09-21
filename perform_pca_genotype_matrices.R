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

main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/Genotype_Mut_Data/"


############################################################
### IMPORT GENOTYPE MATRIX FILES
############################################################

# OPTION 1: MAF FILES
genotype_matrices_fns <- list.files(paste(main_path, "MAF_files/", sep = ""), pattern = ".muse")
genotype_matrix_list <- lapply(genotype_matrices_fns, function(fn) {
  full_fn <- paste(main_path, paste("MAF_files/", fn, sep = ""), sep = "")
  gn_matrix <- import.MAF(full_fn, sep = "\t", is.TCGA = TRUE, 
                          merge.mutation.types = FALSE, to.TRONCO = FALSE) 
  return(gn_matrix)
})

# OPTION 2: VCF FILES
genotype_matrices_fns <- list.files(paste(main_path, "VCF_files/", sep = ""), pattern = ".vcf")
genotype_matrix_list <- lapply(genotype_matrices_fns, function(fn) {
  full_fn <- paste(main_path, paste("VCF_files/", fn, sep = ""), sep = "")
  
  uuid <- unlist(strsplit(fn, split = ".", fixed = TRUE))[1]
  
  # Use the SeqArray package to convert VCF to GDS file type (recommended)
  gn_matrix_gds <- SeqArray::seqVCF2GDS(full_fn, 
                                        out = paste(main_path, paste("GDS_files/", paste(uuid,
                                                                   ".gds", sep = ""), sep = ""), sep = ""))  
  # Alternative: Use SNPRelate package to convert VCF to GDS
  #SNPRelate::snpgdsVCF2GDS_R(full_fn, paste(main_path, paste("GDS_files/", 
                                                             #paste(uuid, ".gds", sep = ""), sep = ""), sep = ""))
  #print(snpgdsSummary(paste(main_path, paste("GDS_files/", 
                                             #paste(uuid, ".gds", sep = ""), sep = ""), sep = "")))
  
  gn_matrix_gds <- seqOpen(paste(main_path, paste("GDS_files/", 
                                                  paste(uuid, ".gds", sep = ""), sep = ""), sep = ""))

  return(gn_matrix_gds)
})

#TODO: do I have to merge the VCF files before or after I convert to GDS?


############################################################
### CALCULATE ALLELE FREQUENCIES
############################################################
#' Given an imported GDS file, calculates the allele frequencies
calculate_af <- function(file) {
  af <- seqApply(file, "genotype", as.is = "double", margin = "by.variant",
                 FUN = CalcAlleleFreq)
  return(af)
}

#' This function in C++ is taken directly from the SeqArray() tutorial 
#' described above
cppFunction("
    double CalcAlleleFreq(IntegerVector x)
    {
        int len=x.size(), n=0, n0=0;
        for (int i=0; i < len; i++)
        {
            int g = x[i];
            if (g != NA_INTEGER)
            {
                n++;
                if (g == 0) n0++;
            }
        }
        return double(n0) / n;
    }")

allele_frequencies <- lapply(genotype_matrix_list, calculate_af)


############################################################
### PERFORM PCA ON GENOTYPE MATRICES
############################################################
# Covariance variable with an initial value
s <- 0

seqApply(f, "$dosage", function(x) {
  p <- 0.5 * mean(x, na.rm=TRUE)      # allele frequency
  g <- (x - 2*p) / sqrt(p*(1-p))      # normalized by allele frequency
  g[is.na(g)] <- 0                    # correct missing values
  s <<- s + (g %o% g)                 # update the cov matrix in the parent environment
}, margin="by.variant", .progress=TRUE)

# Scaled by the number of samples over the trace
s <- s * (nrow(s) / sum(diag(s)))

# Eigen-decomposition
eig <- eigen(s)


############################################################
### VISUALIZE PCA RESULTS
############################################################
plot(eig$vectors[,1], eig$vectors[,2], xlab="PC 1", ylab="PC 2")


############################################################
### SAVE PCA RESULTS TO FILES
############################################################






############################################################
### ALTERNATIVE: MANUALLY CREATE FILES FOR EIGENSTRAT's
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
