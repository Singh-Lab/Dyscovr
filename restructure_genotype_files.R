library(reshape2)
library(data.table)

path <- "/home/scamilli/Preprocessed_Genotype_Sayaman/"

# Format of the QC-corrected data from the 1000G project: VCF format, one file per chromosome
# Column names are samples with Birdseed IDs, rows are SNPs, entries are 0/1, 1/1, or ./.

############################################################
### PART ONE: RESTRUCTURE EACH CHROMOSOME FILE TO USE: THE TCGA
### ID OF THE GIVEN PATIENT AND TO BE EITHER 0 (no variant), 1 (variant),
### or NA (no info); alternative would be to use the vcfR package
############################################################
birdseed_tcga_id_map <- fread(paste0(path, "information_file_composition/Map_TCGAPatientID_BirdseedFileID.txt"),
                              header = F, check.names = F)[,c(1,3), with=F]
colnames(birdseed_tcga_id_map) <- c("Birdseed", "TCGA")

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

# For each chromosome file, create a copy of the file that uses TCGA IDs 
# rather than Birdseed IDs and changes the variant information to be integer type
# either 0 or 1 or NA
chrom_files <- list.files(paste0(path, "Controlled_Data/1000G_Stranded_vcf/"), pattern = ".vcf")

# First 9 columns are informational, and first 5 are the most important (CHROM, POS, ID, REF, ALT)
lapply(chrom_files, function(f) {
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
})

