############################################################
# Code to create Suppl. Figure 1 Visualizations
# Written by Sara Geraghty
# PUBLICATION INFORMATION
############################################################

library(data.table)
library(ggplot2)
library(ggsci)
library(cowplot)

# Local PATH to directory containing Dyscovr output files
PATH <-  paste0(getwd(), "Output/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)

############################################################
### IMPORT PAN-CANCER OUTPUT FILE(S)
############################################################
outfn <- "res_top_0.05_allGenes_quantile_rawCNA_methMRaw_3PCs_Nonsyn.Drivers.Vogel.elim.vif5.sp0.7_corrected_MUT.csv"
pc_allGenes <- read.csv(paste0(PATH, outfn))

# Associated eliminated variables from multicollinearity
elim_vars_fn <- "res_top_0.05_orig_allGenes_quantile_rawCNA_methMRaw_3PCs_elim.vif5.sp0.7_uncorrected_eliminated_variables.csv"
pc_elimVars <- read.csv(paste0(PATH, elim_vars_fn))

# Randomized output
outfn_rand <- "res_top_0.05_allGenes_quantile_rawCNA_methMRaw_3PCs_Nonsyn.Drivers.Vogel.elim.vif5.sp0.7_corrected_MUT_RANDOMIZED.csv"
pc_allGenes_rand <- read.csv(paste0(PATH, outfn_rand))

# Full output, with all covariates
outfn_full <- "res_top_0.05_allGenes_quantile_rawCNA_methMRaw_3PCs_Nonsyn.Drivers.Vogel.elim.vif5.sp0.7_uncorrected.csv"
pc_allGenes_full <- read.csv(paste0(PATH, outfn_full))

### PER-CANCER ###
perCancer_fns <- intersect(list.files(path = PATH, recursive = T,
                                         pattern = "_corrected_MUT"), 
                              intersect(list.files(path = PATH, recursive = T,
                                                   pattern = "allGenes_"),
                                        list.files(path = PATH, recursive = T,
                                                   pattern = "Nonsyn.Drivers.Vogel.elim.vif.5")))
perCancer <- lapply(perCancer_fns, function(f) 
  fread(paste0(PATH, f), header = T))
names(perCancer) <- unlist(lapply(perCancer, function(x)
  unlist(strsplit(x, "/", fixed = T))[1]))


############################################################
### PART A: REAL V. RANDOMIZED HISTOGRAM
############################################################
#' Plot a multi-layer histogram showing the p-value distributions for 
#' a particular driver (or all drivers, if NA), both "real" and "randomized" 
#' versions (H0 vs. HA) with significant hits (defined as being below a given 
#' q-value threshold) highlighted
#' @param real_master_df a master DF for the "real" run (HA)
#' @param random_master_df a master DF for the "randomized" run (H0)
#' @param goi a driver gene-of-interest
#' @param qval_thres a q-value threshold for significance
plot_pvalue_histograms_real_and_random <- function(real_master_df, 
                                                   random_master_df, 
                                                   goi, qval_thres) {
  
  pvals_real <- real_master_df[, c('p.value', 'T_k.name')]
  pvals_random <- random_master_df[, c('p.value', 'T_k.name')]
  
  # Subset to the given driver gene, if one is specified
  if(!is.na(goi)) {
    pvals_real <- real_master_df[real_master_df$R_i.name == goi, 
                                 c('p.value', 'T_k.name')]
    pvals_random <- random_master_df[random_master_df$R_i.name == goi, 
                                     c('p.value', 'T_k.name')]
  } 
  
  pvals <- merge(pvals_real, pvals_random, by = "T_k.name")
  colnames(pvals)[colnames(pvals) == "p.value.x"] <- "Real"
  colnames(pvals)[colnames(pvals) == "p.value.y"] <- "Random"
  pvals_m <- melt(pvals)
  
  # Get the p-value corresponding to q-value threshold (if we want to add a 
  # dotted line corresponding to this threshold)
  #max_qval_under_thres <- max(real_master_df[real_master_df$q.value < qval_thres, 
  #                                           'q.value'])
  #pval_thres <- real_master_df[real_master_df$q.value == max_qval_under_thres, 
  #                            'p.value']
  
  p <- ggplot(pvals_m, aes(x = value, fill = variable)) + 
    geom_histogram(alpha = 0.4, position = "identity", bins = 50) + 
    theme_minimal() + xlab("P-Value") + ylab("Frequency") + 
    labs(fill = "") + theme(axis.title = element_text(face = "bold", size = 14),
                            axis.text = element_text(face = "bold", size = 12),
                            legend.text = element_text(size=12), 
                            legend.position = c(0.8, 0.8),
                            legend.title = element_text(face = "bold", size = 14)) +
    scale_fill_nejm() 
  # To add a dotted line for the p-value threshold
  #geom_vline(xintercept = pval_thres, linetype='dotted', linewidth = 1) + 
  
  p
}

# Call function
plot_pvalue_histograms_real_and_random(pc_allGenes, pc_allGenes_rand, NA, 0.01)


############################################################
### PART B: FREQUENCY OF TOP COVARIATES
############################################################
# Remove the (Intercept) terms
pc_allGenes_full <- pc_allGenes_full[pc_allGenes_full$term != "(Intercept)",]

# Of the top X remaining tests, get what the terms are and plot the proportions
x <- 500
pc_allGenes_full_topx <- pc_allGenes_full[1:x,]
terms <- unique(pc_allGenes_full$term)
terms_counts <- unlist(lapply(terms, function(x) 
  nrow(pc_allGenes_full_topx[pc_allGenes_full_topx$term == x,])))
terms_counts_df <- data.frame('Term' = terms, 'Freq' = terms_counts)

# Visualize this as a Barplot 
ggplot(data = terms_counts_df[1:5,], aes(x = reorder(Term, -Freq), 
                                         y = Freq, fill = Term)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  theme_minimal() + scale_fill_nejm() +
  xlab("Covariate") + ylab(paste0("\n", paste0("Frequency Among Top ", x))) +
  theme(axis.text = element_text(size = 16, face = "bold"), 
        axis.text.x = element_text(angle = 45, vjust =1, hjust=1), 
        axis.title=element_text(size=18,face="bold"), 
        legend.title=element_text(size=14, face="bold"),
        legend.text = element_text(size=12))


############################################################
### PART C: FREQUENCY VARIABLES REMOVED FROM COLLINEARITY ANALYSIS (VIF)
############################################################
#' Uses the files output from the correct_collinearity functions to calculate 
#' the number of times each variable was removed due to collinearity
#' @param elim_vars_df a data frame with the variables removed from the 
#' regression for each target gene
calculate_collinearity_removal_freq <- function(elim_vars_df) {
  elim_vars <- unlist(lapply(elim_vars_df$Variables.Removed, function(x) 
    unlist(strsplit(x, ",", fixed = T))))
  elim_vars_table <- table(elim_vars)
  
  # Adjust table
  elim_vars_table_sub <- as.data.frame(elim_vars_table[
    !grepl("Cancer_type", names(elim_vars_table))])
  elim_vars_table_sub$Freq <- as.numeric(elim_vars_table_sub$Freq)
  elim_vars_table_sub$elim_vars <- as.factor(elim_vars_table_sub$elim_vars)
  
  # Make more readable
  elim_vars_table_sub$elim_vars <- unlist(lapply(elim_vars_table_sub$elim_vars, function(x) {
    x <- as.character(x)
    if(grepl("_k", x)) {return(paste0("Target.", unlist(strsplit(x, "_", fixed = T))[1]))}
    else if(grepl("_i", x)) {
      spl <- unlist(strsplit(x, "_", fixed = T))
      d <- unique(pc_allGenes[pc_allGenes$R_i == spl[1], 'R_i.name'])
      return(paste(d, spl[2], sep = "."))
    } else {return(gsub("_", ".", x))}
  }))
  
  return(elim_vars_table_sub)
}

# Call function
elim_vars_table <- calculate_collinearity_removal_freq(pc_elimVars)

# Create barplot
ggplot(elim_vars_table, aes(x = reorder(elim_vars, -Freq), y = Freq)) + 
  geom_bar(position = "dodge", width = 0.95, stat = "identity", 
           show.legend = FALSE, color = "black", fill = "#20854EFF") + 
  xlab("Eliminated Variable") + ylab("Frequency") +
  theme(axis.text = element_text(face="bold", size = 16, angle = 45, 
                                 vjust = 0.5, hjust = 1), 
        axis.title=element_text(size=18, face="bold"), 
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = 'white'))


############################################################
### PART D: COWPLOT OF NONSYNONYMOUS MUTATION FREQUENCY
############################################################
# Import necessary files, including clinical DF (for patient-cancer type mapping)
# and nonsynonymous mutation matrix
clinical_df <- read.csv(paste0(PATH, "clinical_data.csv"),
                        header = T, check.names = F)
mut_count_matrix <- read.csv(paste0(PATH, "Dyscovr_Input/mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv"),
                             header = T, check.names = F)

# Generate a sample/patient to cancer type mapping
get_patient_cancer_mapping <- function(clinical_df) {
  specific_types <- unique(clinical_df$project_id)
  patient_cancer_mapping <- lapply(specific_types, function(ct) {
    pats <- clinical_df[grepl(ct, clinical_df$project_id),'case_submitter_id']
    pats_ids <- unlist(lapply(pats, function(pat) 
      unlist(strsplit(pat, "-", fixed = TRUE))[3]))
    return(unique(pats_ids))
  })
  names(patient_cancer_mapping) <- specific_types
  return(patient_cancer_mapping)
}

# Create mapping
patient_cancer_mapping <- get_patient_cancer_mapping(clinical_df)
names(patient_cancer_mapping) <- unlist(lapply(names(patient_cancer_mapping), function(x)
  unlist(strsplit(x, "-", fixed = T))[2]))


# Identify unique drivers across the various cancer types
unique_drivers <- unique(unlist(lapply(perCancer, function(x) 
  unique(x$R_i.name)))) # 47 drivers
unique_drivers_uniprot <- unique(unlist(lapply(perCancer, function(x) 
  unique(x$R_i))))

#' Get the per-cancer frequencies of each of these drivers
#' @param unique_drivers a vector with the names of each unique driver gene
#' @param unique_drivers_uniprot a vector with the Uniprot IDs of each unique
#' driver gene
#' @param mut_count_matrix a mutation count matrix across patients
#' @param patient_cancer_mapping a list of cancer types, with each list entry
#' containing a vector of patient IDs for all patients belonging to that cancer
#' type
get_driver_freq_dfs <- function(unique_drivers, unique_drivers_uniprot, 
                                mut_count_matrix, patient_cancer_mapping) {
  
  driver_dfs <- lapply(1:length(unique_drivers), function(d_i) {
    d <- unique_drivers[d_i]
    d_uniprot <- unique_drivers_uniprot[d_i]

    colnames(mut_count_matrix) <- unlist(lapply(colnames(mut_count_matrix), function(x)
      unlist(strsplit(x, "-", fixed = T))[1]))
    mut_count_matrix_d <- mut_count_matrix[mut_count_matrix$Gene_Symbol == d,]
    
    ct_dfs <- lapply(1:length(patient_cancer_mapping), function(i) {
      ct <- names(patient_cancer_mapping)[i]
      patients_ct <- unlist(colnames(mut_count_matrix_d))[
        2:(ncol(mut_count_matrix_d)-1)]
      
      if(ct != "PC") {
        patients_ct <- unlist(patient_cancer_mapping[
          names(patient_cancer_mapping) == ct])
      }
      mutation_df_ct <- mut_count_matrix_d[, which(colnames(mut_count_matrix_d) %fin% 
                                                    patients_ct)]
      freq <- as.data.frame(length(which(as.integer(mutation_df_ct) > 0)) / 
                              (length(which(as.integer(mutation_df_ct) == 0)) + 
                                 length(which(as.integer(mutation_df_ct) > 0))))
      
      colnames(freq) <- "Freq"
      freq$Driver <- d
      freq$Cancer.Type <- ct
      
      return(freq)
    })
    
    ct_dfs <- ct_dfs[!is.na(ct_dfs)]
    ct_df_full <- do.call(rbind, ct_dfs)
    return(ct_df_full)
  })
  
  return(driver_dfs)
}

driver_dfs <- get_driver_freq_dfs(unique_drivers, unique_drivers_uniprot,
                                  mut_count_matrix, patient_cancer_mapping)

# Combine these into one data frame for plotting
driver_freq_df <- do.call(rbind, driver_dfs)
colnames(driver_freq_df)[which(colnames(driver_freq_df) == "Freq")] <- "Frequency"
# Remove those rows with frequency zero
driver_freq_df <- driver_freq_df[driver_freq_df$Frequency > 0,]

# Identify the number of cancer types in which each driver gene is present, for
# visual sorting purposes
num_cts <- unlist(lapply(unique(driver_freq_df$Driver), function(d)
  return(nrow(driver_freq_df[driver_freq_df$Driver == d,]))))
names(num_cts) <- unique(driver_freq_df$Driver)
driver_freq_df$Num.CTs <- unlist(lapply(driver_freq_df$Driver, function(d)
  num_cts[names(num_cts) == d]))
# Also identify the absolute total frequency of mutation of each driver pan-cancer,
# as another option for sorting
driver_freq_df$Total.Freq <- unlist(lapply(driver_freq_df$Driver, function(d)
  sum(driver_freq_df[driver_freq_df$Driver == d, 'Frequency'])))

# Optionally limit just to those with frequency >= 5%
driver_freq_df_5perc <- driver_freq_df
driver_freq_df_5perc <- driver_freq_df_5perc[driver_freq_df_5perc$Frequency > 0.05,] 

# Create cowplot visualization
ggplot(driver_freq_df_5perc, aes(x=Cancer.Type, y = reorder(Driver, Num.CTs), 
                           color = Frequency, size = Frequency)) + 
  geom_point() + cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.title = element_text(face = "bold", size = 14), 
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12), axis.ticks = element_blank(),
        axis.line  = element_blank()) +
  ylab('Driver Gene') + xlab("\nCancer Type") + 
  scale_color_gradient(low = "#20854EFF", high = "#FFDC91FF", 
                       name = 'Nonsynonymous\nMutation Frequency') +
  background_grid() #,oob = scales::squish)

# Separate into cancer types with > 75 samples and those with less
cancer_types_lt75 <- c("CHOL", "DLBC", "GBM", "KICH", "LAML", "OV", "READ", 
                       "SARC", "STAD", "TCGT", "THYM", "UCS", "UVM")
cancer_types_gr75 <- c("ACC", "BLCA", "BRCA", "CESC", "COAD", "ESCA", "HNSC",
                       "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC",
                       "MESO", "PAAD", "PRAD", "PCPG", "THCA", "UCEC")

ggplot(driver_freq_df_5perc, aes(x=Cancer.Type, y = reorder(Driver, Total.Freq), 
                           color = Frequency, size = Frequency)) + 
  geom_point() + cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        axis.title = element_text(face = "bold", size = 14), 
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12), axis.ticks = element_blank(),
        axis.line  = element_blank()) +
  ylab('Driver Gene') + xlab("\nCancer Type") + 
  scale_color_gradient(low = "#20854EFF", high = "#FFDC91FF", 
                       name = 'Nonsynonymous\nMutation Frequency') + #, 
  scale_x_discrete(limits = c(cancer_types_gr75, cancer_types_lt75)) +
  background_grid() #,oob = scales::squish)


