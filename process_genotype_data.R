############################################################
### PROCESS GENOTYPE DATA
### Written By: Sara Geraghty, June 2023
############################################################
# More information about GDC pre-processing and mutation calling: 
# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#masked-somatic-aggregation-workflow

# Processes Affymetrix genotype SNP6.0 files in order to extract germline variants that may module the 
# downstream expression changes of somatic alterations

library(reshape2)
library(data.table)

# Adjust this path for the cluster
all_genes_id_conv <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", 
                           header = TRUE)

# Limiting to variants within 1MB of the transcription start site of a protein coding gene
MIN_DIST <- 1000000  # 1MB

############################################################
## 1.PREPROCESSED GENOTYPE FILES FROM SAYAMAN ET AL. 
## Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8414660/
## Github: https://github.com/rwsayaman/TCGA_PanCancer_Genotyping_Imputation/blob/main/QualityControlAnalysis/README.md
## Link to Downloads: https://gdc.cancer.gov/about-data/publications/CCG-AIM-2020

# Have QC-corrected, stranded data for both HRC (haplotype reference consortium)
# and 1000G project. Palindromic SNPs removed prior to stranding. Both in VCF format.

# NOTE: Stranding is essentially comparing genotype SNPs to SNPs from the HRC or 
# the 1000K Genomes project, and removing SNPs with differing alleles, SNPs with > 0.2 
# allele frequency difference, and SNPs not in the reference panel

# $module add R/4.2.1
# $R

path <- "/home/scamilli/Preprocessed_Genotype_Sayaman/"

# Format of the QC-corrected data from the 1000G project: VCF format, one file per chromosome
# Column names are samples with Birdseed IDs, rows are SNPs, entries are 0/1, 1/1, or ./.

############################################################
### PART ONE: RESTRUCTURE EACH CHROMOSOME FILE TO USE: THE TCGA
### ID OF THE GIVEN PATIENT AND TO BE EITHER 0 (no variant), 1 (variant),
### or NA (no info); alternative would be to use the vcfR package
# NOTE: Stored in gen-singhtmp.princeton.edu cluster as restructure_genotype_files.R
############################################################
birdseed_tcga_id_map <- fread(paste0(path, "information_file_composition/Map_TCGAPatientID_BirdseedFileID.txt"),
                              header = F, check.names = F)[,c(1,3), with=F]
colnames(birdseed_tcga_id_map) <- c("Birdseed", "TCGA")

# For each chromosome file, create a copy of the file that uses TCGA IDs 
# rather than Birdseed IDs and changes the variant information to be integer type
# either 0 or 1 or NA
chrom_files <- list.files(paste0(path, "Controlled_Data/1000G_Stranded_vcf/"), pattern = ".vcf")

# First 9 columns are informational, and first 5 are the most important (CHROM, POS, ID, REF, ALT)
chrom_files_new <- lapply(chrom_files, function(f) {
  print(paste("Now processing...", f))
  
  file <- fread(paste0(path, paste0("Controlled_Data/1000G_Stranded_vcf/", f)), header = T)
  colnames(file)[1] <- "CHROM"
  new_colnames <- unlist(lapply(colnames(file)[10:ncol(file)], function(c) {
    spl_c <- unlist(strsplit(c, "_", fixed=T))
    ind <- which(spl_c == spl_c[1])[2]
    spl_c <- paste(spl_c[1:(ind-1)], collapse = "_")
    tcga_c <- as.character(birdseed_tcga_id_map[birdseed_tcga_id_map$Birdseed == spl_c, 'TCGA'])
    return(tcga_c)
  }))
  colnames(file) <- c(colnames(file)[1:9], new_colnames)
  
  for (j in 10:ncol(file)) set(file, j = j, value = convert_to_integer(file[[j]]))
  
  new_fn <- unlist(strsplit(unlist(strsplit(f, ".", fixed = T))[4], "_", fixed = T))
  new_fn <- paste0(new_fn[length(new_fn)], "_preprocessed.csv")
  fwrite(file, paste0(path, paste0("Controlled_Data/1000G_Stranded_vcf/", new_fn)))
  
  return(file)
})

convert_to_integer <- function(val) {
  #print(head(val))
  new_vals <- unlist(lapply(val, function(v) {
    if(v == "0/0") {return(as.integer(0))}   # homozygous reference
    else if (v == "0/1") {return(as.integer(1))} # heterozygous
    else if (v == "1/1") {return(as.integer(2))}  # homozygous alternate
    else {return(NA)}
  }))
  return(new_vals)
}

# $Rscript restructure_genotype_files.R

############################################################
### PART TWO: RESTRICT FILES TO GERMLINE VARIANTS WITHIN 1 MB
### OF A TRANSCRIPTION START SITE (TSS) AND ADD GENE ANNOTATION
############################################################
# Limiting to variants within 1MB of the transcription start site of a protein coding gene
MIN_DIST <- 1000000  # 1MB

all_genes_id_conv <- fread(paste0(path, "Controlled_Data/1000G_Stranded_vcf/all_genes_id_conv.csv"), header = TRUE)
all_genes_id_conv_sub <- all_genes_id_conv[all_genes_id_conv$uniprot_gn_id != "",]
chrom_files <- list.files(paste0(path, "Controlled_Data/1000G_Stranded_vcf/"), pattern = "preprocessed.csv")

# Get the TSS of all human genes
transcription_start_sites <- unique(all_genes_id_conv_sub$transcription_start_site)
max_pos <- max(all_genes_id_conv_sub$end_position)

# Get the ranges of positions within 1MB of a TSS
#ranges_accept_pos <- unique(unlist(lapply(transcription_start_sites, function(tss) {
#  rang_up <- max(1, as.integer(tss)-MIN_DIST)
#  rang_down <- min(as.integer(tss)+MIN_DIST, max_pos)
#  return(c(rang_up:as.integer(tss), (as.integer(tss)+1):rang_down))
#})))
start_pos <- max(1, as.integer(transcription_start_sites[1])-MIN_DIST)
end_pos <- min(as.integer(transcription_start_sites[length(transcription_start_sites)])+MIN_DIST, max_pos)
curr_pos <- start_pos

while(curr_pos < end_pos) {
  curr_end_pos <- curr_pos + MIN_DIST
  if(!any((transcription_start_sites >= curr_pos) & (transcription_start_sites < curr_end_pos))) {
    print("REGION WITH NO TSS")
  }
  curr_pos <- curr_pos + MIN_DIST
}

# There are no 1MB regions without any TSS -- do we want to use something smaller than this?

# Associate each SNP with its closest gene
lapply(chrom_files, function(f) {
  print(paste("Now processing...", f))

  file <- fread(paste0(path, paste0("Controlled_Data/1000G_Stranded_vcf/", f)), header = T)
  colnames(file)[1] <- "CHROM"
  
  gene_info <- lapply(file$POS, function(pos) {
    return(all_genes_id_conv_sub[which.min(abs(all_genes_id_conv_sub$transcription_start_site - pos)),
                          c('ensembl_gene_id', 'uniprot_gn_id', 'external_gene_name')])
  })
  gene_info_dt <- rbindlist(gene_info)
  colnames(gene_info_dt) <- c("ENSG", "UNIPROT", "GENE_NAME")
  
  file_new <- cbind(file[,1:5, with=F], cbind(gene_info_dt, file[,(6:ncol(file)),with=F]))
  fwrite(file_new, paste0(path, paste0("Controlled_Data/1000G_Stranded_vcf/", f)))
})


############################################################
### PART THREE: FOR EACH GENE, GET ALL THE SNPs PROXIMAL TO IT; 
### CREATE ONE DF FOR ALL CHROMOSOMES WITH THE SNPs ASSOC. TO EACH GENE
############################################################
# TODO: currently assoc. each SNP with one gene; could a SNP be assoc. with multiple genes?

# Format is a gene/SNP by patient DF


















#
#
#
#
#
#
#


############################################################
### ARCHIVAL CODE
############################################################
# Import the mapping between TCGA IDs and Birdseed IDs
id_mapping <- fread(paste0(path, "OpenAccess_data/Map_TCGAPatientID_BirdseedFileID.txt"), 
                    header = F)[,c(1,3), with = F]
colnames(id_mapping) <- c("birdseed", "tcga")

# Output table format: sample(column) x snp status(row)
# Column names are sample IDs, row names are SNPs
output_snp_df <- data.frame(matrix(ncol = nrow(id_mapping), nrow = 0))
colnames(output_snp_df) <- id_mapping$tcga

# Also maintain a mapping between target genes and SNP IDs (list);
# names of list will be target genes
snp_list <- list(rep(NA, times = length(unique(all_genes_id_conv$ensembl_gene_id))))
names(snp_list) <- unique(all_genes_id_conv$ensembl_gene_id)

# For each chromosome, get the target genes on that chromosome and their start and end positions.
# Then, for each, look for a SNP position(s) that is found nearby that target gene (<1MB).
# If there is, record whether each individual contains that SNP. 
for (chrom in 1:22) {
  
  # Get target genes for this chromosome
  target_genes <- all_genes_id_conv[all_genes_id_conv$chromosome_name == chrom, 
                                    c('ensembl_gene_id', 'start_position', 'end_position')]
  
  # Import that chromosome's genotype data file
  fn <- "TCGAGenotyping_ZivLabWhitelist10128.v6_clean.rmPalindromic.1000G_chr"
  geno_file <- fread(paste0(path, paste0(fn, paste0(chrom, ".vcf"))), header = T,
                     skip = 6)  # Skip the comment characters
  # The first 9 columns are: CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
  # followed by the sample IDs; we are interested in the position and whether each
  # sample has that SNP
  
  # Keep only SNPs that are within 1MB of a target gene from our set
  geno_file_sub <- geno_file[]

  gc()
}


############################################################
## 2. FORMAT GENOTYPE FILES FROM THE TCGA LEGACY ARCHIVE ##

#path <- "/home/scamilli/Breast_Cancer/SNP 6.0 Genotype Arrays/"
path <- "/home/scamilli/Pan-Cancer/Genotype_Mut_Data/SNP6_Genotype_files/"

# Modeled off of code from: https://bioinfo.uth.edu/LungCancerSubtypes/code_for_pca_with_TCGA_data.txt
# Gives detailed R code for performing PCA with TCGA genotype files
# To run R in this cluster, first run: module add R

# Additional useful PLINK vignettes: 
# Processing 1000 Genomes reference data for ancestry estimation: https://meyer-lab-cshl.github.io/plinkQC/articles/Genomes1000.html
# Ancestry estimation based on reference samples of known ethnicities: https://meyer-lab-cshl.github.io/plinkQC/articles/AncestryCheck.html
# Tutorial:Produce PCA bi-plot for 1000 Genomes Phase III - Version 2: https://www.biostars.org/p/335605/

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
# brca.germline.blood <- brca.final.matched[grepl("10", brca.split.barcode$sample),]
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
