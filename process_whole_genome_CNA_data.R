############################################################
### Process Whole-Genome CNA Data 
### Written By: Sara Camilli, February 2021
############################################################

# This file processes the whole-genome CNA (asgatNGS Files) data file(s) in order 
# to merge across all patients of interest

# Format of output dataframe: 
# Rows : Patient TCGA ID (4-digit) OPT: (-Sample ID (01A is tumor, 11A is normal))
# Columns : Chromosome, Start, End, Copy Number, Major
# Entries : CNA value (positive represents an amplification event, negative represents deletion event,
# 0 represents no amplification or deletion)

#library(copynumber)
library(dplyr)
library(rlang)
library(rlist)
library(ggplot2)
library(DescTools)

path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/CC Data/"
#path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/TCGA Data (ALL)/"

output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/CC/CNV/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/CNV/"


############################################################
### INPUT THE WHOLE-GENOME CNA DATA
############################################################
filenames <- list.files(paste(path, "WGS_CNV/", sep = ""), pattern = ".copy_number_variation.seg")
  # 154 files for CC

############################################################
### VISUALIZE THE WHOLE-GENOME CNA DATA
############################################################
# Start with chromosome 21 across all patients
# Get the length of this chromosome manually
chrom <- "chr21"
start_val <- 1 
end_val <- 20000000

start_val <- 17000000
end_val <- 22000000

start_val <- 10000000
end_val <- 15000000
length_chrom <- end_val - start_val

#' Visualization Function 1: per-nucleotide, calculates the frequency that this 
#' nucleotide is part of a CNA event (frequency in # of patients)
#' @param filenames a vector of the filenames for all patients in the given group
#' @param chrom a chromosome of interest
#' @param length_chrom the length of the chromosome or segment, in base pairs
#' @param end_val the end of the segment we are looking at, in bp
#' @param start_val the start of the segment we are looking at
#' @param label "amplification", for CN >=6, "deletion", for CN == 0, or 
#' "copy number loss", for CN == 1
get_cna_events_per_bp <- function(filenames, chrom, length_chrom, end_val, start_val, label) {
  
  # We need to create a per-nucleotide frequency table that we can plot 
  # Rows: nucleotide positions
  # Column: frequency (# of patients)
  freq_table <- data.frame(matrix(nrow = length_chrom + 1, ncol = 2))
  colnames(freq_table) <- c("BP", "Freq")
  rownames(freq_table) <- 1:(length_chrom + 1)
  
  # Loop through all patients
  freq_vectors <- lapply(1:length(filenames), function(i) {
    print(paste(i, paste("/", length(filenames))))
    
    # Import the CNA table for this patient
    cna_table <- data.table::fread(paste(path, paste("WGS_CNV/", filenames[i], sep = ""), sep = ""), sep = "\t", 
                                   header = TRUE, select = c(Chromosome='character', 
                                                             Start='integer', 
                                                             End='integer', 
                                                             Copy_Number='integer'), 
                                   strip.white = TRUE)
    
    # Limit to just the chromosome of interest, and only affecting the segment of interest
    # Either the whole CNA is in range, the start val is in range, or the end val is in range
    cna_table_sub <- cna_table %>% filter((Chromosome == chrom) & (Start < end_val) 
                                          & (End > start_val)) 
    
    # Filter based on the kind of CNA event we are looking at
    if (label == "amplification") {
      cna_table_sub <- cna_table_sub %>% filter(Copy_Number >= 6)
    } else if (label == "deletion") {
      cna_table_sub <- cna_table_sub %>% filter(Copy_Number == 0)
    } else if (label == "copy number loss") {
      cna_table_sub <- cna_table_sub %>% filter(Copy_Number == 1)
    } else {
      print(paste("Inappropriate label", label))
      return(NA)
    }
    
    # Make a binary 0/1 vector for each position in this chromosome, to show whether or
    # not it was amplified in the given patient
    freq_vect <- rep(0, times = length_chrom + 1)
    
    if (nrow(cna_table_sub) > 0) {
      # Get all the unique positions with a CNA event in this chromosome
      positions_to_add <- unique(unlist(lapply(1:nrow(cna_table_sub), function(j) {
        
        #print(cna_table_sub[j,])
        
        # Get the start and end positions of this CNA range
        start_pos <- as.integer(cna_table_sub[j, 'Start'])
        end_pos <- as.integer(cna_table_sub[j, 'End'])
        print(start_pos)
        print(end_pos)
        
        # Make sure this range overlaps our segment
        #if(length(intersect(start_pos:end_pos, start_val:end_val)) == 0) {return(NA)}
        if(!(TRUE %in% (start_pos:end_pos %overlaps% start_val:end_val))) {return(NA)}

        # If they extend outside our segment of interest, restrict them to the segment
        if(start_pos < start_val) {start_pos <- start_val}
        if(end_pos > end_val) {end_pos <- end_val}
        
        return(start_pos:end_pos)
      })))
      #print(head(positions_to_add))
      
      # Update these positions so that they are 1-indexed to our segment and we 
      # can put them in the existing vector
      positions_to_add <- positions_to_add[!is.na(positions_to_add)] - (start_val - 1)   
                                                            # the -1 is to prevent 0-indexing
      print(positions_to_add)
      
      # Update these positions as '1' in the vector
      freq_vect[positions_to_add] <- 1
      
    }
    print(length(freq_vect))
    return(freq_vect)
  })
  
  # Get the sum of all the '1' values across all patients (columns) in this data 
  # frame for each BP (row) and add to DF
  freq_table$BP <- start_val:end_val
  freq_table$Freq <- rowSums(list.cbind(freq_vectors), na.rm = TRUE) 
  
  print(freq_table)
  return(freq_table)
}


# Run the function
chrom_freq_table <- get_cna_events_per_bp(filenames, chrom, length_chrom, end_val, 
                                          start_val, "amplification")
chrom_freq_table <- get_cna_events_per_bp(filenames, chrom, length_chrom, end_val, 
                                          start_val, "deletion")
chrom_freq_table <- get_cna_events_per_bp(filenames, chrom, length_chrom, end_val, 
                                          start_val, "copy number loss")


# Write the frequency table result to a file
write.csv(chrom_freq_table, paste(output_path, paste("frequency_cnas_per_bp_", 
                                                     paste(chrom, ".csv", sep = ""), sep = ""), sep = ""))

# Read back if needed
chrom_freq_table <- read.csv(paste(output_path, paste("frequency_cnas_per_bp_", 
                                                      paste(chrom, ".csv", sep = ""), sep = ""), sep = ""), header = TRUE)

# Visualize
barplot(chrom_freq_table$Freq, xlab = "Base Pair Position", ylab = "Frequency (Num. of Patients)", 
     main = paste("Number of Patients with a CNA Event Affecting Each Base, ", chrom, sep = ""), space = 0)
# Add a smooth curve
lines(density(as.numeric(chrom_freq_table$Freq)), col = "blue", lwd = 2)
lines(density(as.numeric(chrom_freq_table$Freq), adjust = 2), lty = "dotted", col = "darkgreen", lwd = 2)

# Alternative visualization using ggplot2
ggplot(data = chrom_freq_table, aes(x=BP, y=Freq)) + geom_line() + 
  labs(title = paste("Number of Patients with a Deletion CNA Event Affecting Each Base, ", chrom, sep = ""),
       x = "Base Pair Position", y = "Frequency (Num. of Patients)")

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


chrom <- "chr21"
#length_chrom <- 47000000
length_chrom_half1 <- length_chrom / 2

#start_val <- 17000000
#end_val <- 22000000

start_val <- 10000000
end_val <- 15000000
length_chrom <- end_val - start_val

#' Visualization Function 2: per-nucleotide, finds the avg. number of copies across patients
#' @param filenames a vector of the filenames for all patients in the given group
#' @param chrom a chromosome of interest
#' @param length_chrom the length of the chromosome or segment, in base pairs
#' @param end_val the end of the segment we are looking at, in bp
#' @param start_val the start of the segment we are looking at
get_avg_cna_events_per_bp <- function(filenames, chrom, length_chrom, end_val, start_val) {
  # Loop through all patients
  cna_count_tables <- lapply(1:length(filenames), function(i) {
    print(paste(i, paste("/", length(filenames))))
    
    # Import the CNA table for this patient
    cna_table <- data.table::fread(paste(path, paste("WGS_CNV/", filenames[i], sep = ""), sep = ""), 
                                   sep = "\t", header = TRUE, select = c(Chromosome='character', 
                                                                         Start='numeric', 
                                                                         End='numeric',
                                                                         Copy_Number='numeric',
                                                                         Major_Copy_Number='numeric', 
                                                                         Minor_Copy_Number='numeric'), 
                                   strip.white = TRUE)
    
    # Limit to just the chromosome of interest, and only affecting the segment of interest
    # Either the whole CNA is in range, the start val is in range, or the end val is in range
    cna_table_sub <- cna_table %>% filter((Chromosome == chrom) & (Start < end_val) & (End > start_val))                      
    print(cna_table_sub)
    
    if (nrow(cna_table_sub) > 0) {
      # Make an integer matrix for each position in this chromosome, to show the number 
      # of times it was amplified in the given patient
      num_cnas_table <- matrix(nrow = length_chrom + 1, ncol = 3)
      #num_cnas_table[,1] <- start_pos:end_pos
      num_cnas_table[,1:3] <- 0
      
      # Update all the positions for this patient in each CNA range that overlaps 
      # our range of interest
      for (j in 1:nrow(cna_table_sub)) {
        # Get the start and end positions of this CNA range
        start_pos <- as.numeric(cna_table_sub[j, 'Start'])
        end_pos <- as.numeric(cna_table_sub[j, 'End'])
        
        # Make sure this range overlaps our segment
        #if(length(intersect(start_pos:end_pos, start_val:end_val)) == 0) {return(NA)}
        if(!(TRUE %in% (start_pos:end_pos %overlaps% start_val:end_val))) {return(NA)}
        
        # If they extend outside our segment of interest, restrict them to the segment
        if(start_pos < start_val) {start_pos <- start_val}
        if(end_pos > end_val) {end_pos <- end_val}
        
        positions <- start_pos:end_pos
        
        # Update them so they are 1-indexed to our segment
        positions <- positions - (start_val - 1)
        
        # Extract the corresponding number of amplifications for this position
        copy_num <- as.numeric(cna_table_sub[j, 'Copy_Number'])
        major_copy_num <- as.numeric(cna_table_sub[j, 'Major_Copy_Number'])
        minor_copy_num <- as.numeric(cna_table_sub[j, 'Minor_Copy_Number'])
        print(c(copy_num, major_copy_num, minor_copy_num))
        
        # Use these values to update the data frame
        num_cnas_table[positions,1] <- num_cnas_table[positions,1] + copy_num
        num_cnas_table[positions,2] <- num_cnas_table[positions,2] + major_copy_num
        num_cnas_table[positions,3] <- num_cnas_table[positions,3] + minor_copy_num
          #c(copy_num, major_copy_num, minor_copy_num)
      }
      return(num_cnas_table)
    }
  })
  
  # We need to create a per-nucleotide CNA amplification average table that we can plot 
  # Rows: nucleotide positions
  # Column: average number of amplification events at this position, across patients
  avg_cna_table <- data.frame(matrix(nrow = length_chrom+1, ncol = 4))
  colnames(avg_cna_table) <- c("BP", "Avg.Copy.Num", "Avg.Major.Copy.Num", "Avg.Minor.Copy.Num")
  rownames(avg_cna_table) <- 1:(length_chrom+1)
  
  # Get the average of all the amplification values across all patients in this 
  # data frame for each BP (row) and add to DF
  #print(cna_count_tables)
  
  # Get copy numbers into a new table
  avg_cna_table$BP <- start_val:end_val
  avg_cna_table$Avg.Copy.Num <- rowMeans(do.call("cbind", lapply(cna_count_tables, function(x) {x[,1]})))   
  avg_cna_table$Avg.Major.Copy.Num <- rowMeans(do.call("cbind", lapply(cna_count_tables, function(x) {x[,2]})))
  avg_cna_table$Avg.Minor.Copy.Num <- rowMeans(do.call("cbind", lapply(cna_count_tables, function(x) {x[,3]})))
  
  
  print(avg_cna_table)
  return(avg_cna_table)
}


# Run the function
chrom_avg_cna_table <- get_avg_cna_events_per_bp(filenames, chrom, length_chrom, 
                                                 end_val, start_val)

# Write the frequency table result to a file
write.csv(chrom_avg_cna_table, paste(output_path, paste("average_num_cnas_per_bp_", 
                                                        paste(chrom, "_BP10-15.csv", sep = ""), 
                                                        sep = ""), sep = ""))

# Read back if needed
chrom_avg_cna_table <- read.csv(paste(output_path, 
                                      paste("average_num_cnas_per_bp_", 
                                            paste(chrom, "_BP10-15.csv", sep = ""), sep = ""), sep = ""), 
                                header = TRUE)

# Visualize
barplot(chrom_avg_cna_table$Avg.Copy.Num, xlab = "Base Pair Position", 
        ylab = "Avg. Number of CNA Events", main = paste("Avg. Number of CNA Events Affecting Each Base, Across Patients, ", chrom, sep = ""), 
        space = 0)
# Add a smooth curve
lines(density(as.numeric(chrom_avg_cna_table$Avg.Copy.Num)), col = "blue", lwd = 2)
lines(density(as.numeric(chrom_avg_cna_table$Avg.Copy.Num), adjust = 2), lty = "dotted", 
      col = "darkgreen", lwd = 2)

# Alternative visualization using ggplot2
ggplot(data = chrom_avg_cna_table, aes(x=BP, y=Avg.Copy.Num)) + geom_line(size = 2) + 
  labs(title = paste("Avg. Number of CNA Events Affecting Each Base, Across Patients, ", chrom, sep = ""),
       x = "Base Pair Position", y = "Avg. Number of CNA Events")

ggplot(data = chrom_avg_cna_table, aes(x=BP, y=Avg.Major.Copy.Num)) + geom_line(size = 2) + 
  labs(title = paste("Avg. Number of Major CNA Events Affecting Each Base, Across Patients, ", chrom, sep = ""),
       x = "Base Pair Position", y = "Avg. Number of Major CNA Events")

ggplot(data = chrom_avg_cna_table, aes(x=BP, y=Avg.Minor.Copy.Num)) + geom_line(size = 2) + 
  labs(title = paste("Avg. Number of Minor CNA Events Affecting Each Base, Across Patients, ", chrom, sep = ""),
       x = "Base Pair Position", y = "Avg. Number of Minor CNA Events")

# Plot all 3 on the same graph
ggplot(data = chrom_avg_cna_table, aes(x=BP))  +
  geom_line(aes(y = Avg.Minor.Copy.Num), color = "darkgreen") + 
  geom_line(aes(y = Avg.Major.Copy.Num), color = "steelblue") + 
  geom_line(aes(y = Avg.Copy.Num), color = "darkred") +
  labs(title = paste("Avg. Number of CNA Events Affecting Each Base, Across Patients, ", chrom, sep = ""),
       x = "Base Pair Position", y = "Avg. Number of Minor CNA Events")

