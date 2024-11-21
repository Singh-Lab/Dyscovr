############################################################
# Code to create Suppl. Figure 3 Visualizations
# Written by Sara Geraghty
# PUBLICATION INFORMATION
############################################################

library(data.table)
library(ggplot2)
library(cowplot)
library(wCorr)

############################################################
### IMPORT PAN-CANCER OUTPUT FILE(S)
############################################################
outfn <- "res_top_0.05_allGenes_quantile_rawCNA_methMRaw_3PCs_Nonsyn.Drivers.Vogel.elim.vif5.sp0.7_corrected_MUT.csv"
pc_allGenes <- read.csv(paste0(PATH, outfn))

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
### PART A: COWPLOT OF NONSYNONYMOUS MUTATION FREQUENCY
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


############################################################
### PART B: Pairwise (Weighted) Spearman Correlation Barplot
############################################################
#' Take the Spearman correlation, weighted Spearman, dot product, or cosine
#' similarity of the coefficients fit in each pairwise set of cancer types at
#' or above some q-value threshold, or a HG test of the hits that are shared
#' at that threshold. Plot barplot containing all driver genes.
#' @param perCancer list of output data frames from Dyscovr per-cancer, named by
#' the cancer type
#' @param qval_thres a q-value threshold for significance
#' @param method either 'Spearman', 'WeightedSpearman', 'DP', 'CS', or 'HG' to 
#' indicate what kind of test we are performing
#' @param drivers a vector of the driver gene names that are highly mutated in
#' at least 3 cancer types/ that we'd like to consider
create_perdriver_percancer_barplot <- function(perCancer, qval_thres, method,
                                               drivers) {
  
  if((method == "Spearman") | (method == "WeightedSpearman") | (method == "DP") | 
     (method == "CS")) {
    
    # Use helper function to get the correlation values for each driver (pairwise
    # between each set of cancer types in which that driver is frequently mutated)
    vals_per_driver <- calculate_correlations_spearman_etc(perCancer, qval_thres, 
                                                           method, drivers)
    input_df <- melt(vals_per_driver)
    
    ## SPEARMAN ## 
    if(method == "Spearman") {
      colnames(input_df) <- c("Spearman", "Driver")
      
      g <- ggplot(input_df, aes(x = Driver, y = Spearman, fill = Driver)) + 
        geom_boxplot() + theme_minimal() +
        scale_fill_manual(values = c("TP53" = "#0072B5FF", "PIK3CA" = "#BC3C29FF", 
                                     "KRAS" = "#20854EFF", "CTNNB1" = "#E18727FF", 
                                     "FBXW7" = "palevioletred3")) +
        geom_hline(yintercept = 0, linetype = "dashed") + 
        ylab("Pairwise Spearman Correlation") + 
        scale_x_discrete(limits = rev) + 
        theme(axis.text = element_text(face = "bold", size = 12),
              axis.title = element_text(face = "bold", size = 14),
              axis.ticks = element_blank(), axis.line = element_blank(),
              panel.border = element_blank(), legend.position = "none",
              panel.grid.major = element_line(color='#eeeeee'))
      print(g)
    }
    ## WEIGHTED SPEARMAN ## 
    if(method == "WeightedSpearman") {
      colnames(input_df) <- c("WeightedSpearman", "Driver")
      print(head(input_df))
      g <- ggplot(input_df, aes(x = Driver, y = WeightedSpearman, fill = Driver)) + 
        geom_boxplot() + theme_minimal() +
        scale_fill_manual(values = c("TP53" = "#0072B5FF", "PIK3CA" = "#BC3C29FF", 
                                     "KRAS" = "#20854EFF", "CTNNB1" = "#E18727FF", 
                                     "FBXW7" = "palevioletred3")) +
        geom_hline(yintercept = 0, linetype = "dashed") + 
        ylab("Pairwise Weighted Spearman Correlation") + 
        scale_x_discrete(limits = rev) + 
        theme(axis.text = element_text(face = "bold", size = 12),
              axis.title = element_text(face = "bold", size = 14),
              axis.ticks = element_blank(), axis.line = element_blank(),
              panel.border = element_blank(), legend.position = "none",
              panel.grid.major = element_line(color='#eeeeee'))
      print(g)
    }
    ## DOT PRODUCT ## 
    if(method == "DP") {
      colnames(input_df) <- c("Dot.Product", "Driver")
      print(head(input_df))
      g <- ggplot(input_df, aes(x = Driver, y = Dot.Product, fill = Driver)) + 
        geom_boxplot() + theme_minimal() +
        scale_fill_manual(values = c("TP53" = "#0072B5FF", "PIK3CA" = "#BC3C29FF", 
                                     "KRAS" = "#20854EFF", "CTNNB1" = "#E18727FF", 
                                     "FBXW7" = "palevioletred3")) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        ylab("Pairwise Dot Product") + 
        scale_x_discrete(limits = rev) + 
        theme(axis.text = element_text(face = "bold", size = 12),
              axis.title = element_text(face = "bold", size = 14),
              axis.ticks = element_blank(), axis.line = element_blank(),
              panel.border = element_blank(), legend.position = "none",
              panel.grid.major = element_line(color='#eeeeee'))
      print(g)
    }
    ## COSINE SIMILARITY ## 
    if (method == "CS") {
      colnames(input_df) <- c("Cosine.Similarity", "Driver")
      print(head(input_df))
      g <- ggplot(input_df, aes(x = Driver, y = Cosine.Similarity, fill = Driver)) + 
        geom_boxplot() + theme_minimal() +
        scale_fill_manual(values = c("TP53" = "#0072B5FF", "PIK3CA" = "#BC3C29FF", 
                                     "KRAS" = "#20854EFF", "CTNNB1" = "#E18727FF", 
                                     "FBXW7" = "palevioletred3")) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        ylab("Pairwise Cosine Similarity") + 
        scale_x_discrete(limits = rev) + 
        theme(axis.text = element_text(face = "bold", size = 12),
              axis.title = element_text(face = "bold", size = 14),
              axis.ticks = element_blank(), axis.line = element_blank(),
              panel.border = element_blank(), legend.position = "none",
              panel.grid.major = element_line(color='#eeeeee'))
      print(g)
    }
    
    ## HYPERGEOMETRIC TEST ## 
  } else if (method == "HG") {
    
    hg_per_driver <- calculate_hg_per_driver(perCancer, qval_thres, method, drivers)
    
    input_df <- melt(hg_per_driver)
    colnames(input_df) <- c("HG.pvalue", "Driver")
    print(head(input_df))
    ggplot(input_df, aes(x = Driver, y = HG.pvalue)) + geom_boxplot()
    
  } else {
    print("Only implemented for methods 'Spearman', 'DP', 'CS', 'WeightedSpearman',
    and 'HG'. Please try again with one of these inputs.")
  }
}

#' Helper function to subset each output data frame in the list by driver d and
#' by q-value; keep only the necessary columns (estimate and target name)
#' @param perCancer list of output data frames from Dyscovr per-cancer
#' @param d the name of the given driver
#' @param qval_thres a q-value threshold for significance
#' @param min_hits the minimum number of hits for a given driver in a given 
#' cancer type
subset_output_list <- function(perCancer, d, qval_thres, min_hits = 15) {
  perCancer_sub <- lapply(1:length(perCancer), function(i) {
    x <- perCancer[[i]]
    ct <- names(perCancer)[i]
    
    if(d %fin% x$R_i.name) {
      # Keep only if there are at least 15 hits for this driver
      if(nrow(x[(x$R_i.name == d) & (x$q.value < 0.2),]) < min_hits) {return(NA)}
      
      est_newname <- paste0('estimate.', ct)
      colnames(x)[which(colnames(x) == 'estimate')] <- est_newname
      
      # Uncomment to use the T-statistic instead
      #est_newname <- paste0('statistic.', ct)
      #colnames(x)[which(colnames(x) == 'statistic')] <- est_newname
      
      est_newname_qval <- paste0('q.value.', ct)
      colnames(x)[which(colnames(x) == 'q.value')] <- est_newname_qval
      
      cols_to_keep <- c(est_newname, 'T_k.name', est_newname_qval)
      
      if(!is.na(qval_thres)) {
        if(nrow(x[(x$q.value < qval_thres) & (x$R_i.name == d),]) > 0) {
          return(x[(x$q.value < qval_thres) & (x$R_i.name == d), 
                   which(colnames(x) %fin% cols_to_keep), with=F])
        } else {return(NA)}
      } else {return(x[x$R_i.name == d, which(colnames(x) %fin% cols_to_keep), 
                       with=F])}
    } else {return(NA)}
  })
  perCancer_sub <- perCancer_sub[!is.na(perCancer_sub)]
  names(perCancer_sub) <- unlist(lapply(perCancer_sub, function(x) 
    return(unlist(strsplit(colnames(x)[1], ".", fixed = T))[2])))
  
  return(perCancer_sub)
}

#' Calculates the correlation (according to the given method) of the pairwise 
#' mutation coefficients for each driver gene 
#' @param perCancer list of output data frames from Dyscovr per-cancer, named by
#' the cancer type
#' @param qval_thres a q-value threshold for significance
#' @param method either 'Spearman', 'WeightedSpearman', 'DP', 'CS', or 'HG' to 
#' indicate what kind of test we are performing
#' @param drivers a vector of the driver gene names that are highly mutated in
#' at least 3 cancer types/ that we'd like to consider
calculate_correlations_spearman_etc <- function(perCancer, qval_thres, method,
                                                drivers) {
  
  vals_per_driver <- lapply(drivers, function(d) {
    perCancer_sub <- subset_output_list(perCancer, d, qval_thres)
    
    # Merge estimates by gene
    if(length(perCancer_sub) > 1) {
      perCancer_merged <- Reduce(function(...) 
        merge(..., by = 'T_k.name', all = T), perCancer_sub)
      
      # Keep only rows that have at least two non-NA values
      #perCancer_merged <- perCancer_merged[rowSums(!is.na(
      #  perCancer_merged)) > 1,]
      
      if(nrow(perCancer_merged) > 9) {
        
        ## SPEARMAN ##
        if(method == "Spearman") {
          perCancer_merged <- perCancer_merged[,!which(grepl("q.value", 
                                                             colnames(perCancer_merged))),
                                               with=F]
          # Generate correlation matrix using the Hmisc package
          corr_mat <- as.data.frame(rcorr(x = as.matrix(
            perCancer_merged[,2:ncol(perCancer_merged)]), 
            type = "spearman")$r)
          vals <- unlist(lapply(1:(ncol(corr_mat)-1), function(i) {
            return(corr_mat[(i+1):nrow(corr_mat),i])
          }))
          return(vals)
          
        } else {
          perCancer_qvals <- perCancer_merged[,c(1, which(grepl("q.value",
                                                                colnames(perCancer_merged)))),
                                              with=F]
          perCancer_merged <- perCancer_merged[,!which(grepl("q.value",
                                                             colnames(perCancer_merged))),
                                               with=F]
          coeff_vects <- lapply(2:ncol(perCancer_merged), function(i) {
            vals <- unlist(perCancer_merged[,i,with=F])
            return(as.numeric(vals))
          })
          qval_vects <- lapply(2:ncol(perCancer_qvals), function(i) {
            vals <- unlist(perCancer_qvals[,i,with=F])
            return(as.numeric(vals))
          })
          
          # Take the pairwise CS/ DP/ Weighted Spearman of each vector
          vals <- c()
          for (i in 1:(length(coeff_vects)-1)) {
            for (j in (i+1):length(coeff_vects)) {
              v1 <- coeff_vects[[i]]
              v2 <- coeff_vects[[j]]
              q1 <- qval_vects[[i]]
              q2 <- qval_vects[[j]]
              ## WEIGHTED SPEARMAN ##
              if (method == "WeightedSpearman") {
                vals_to_rm <- unique(c(which(is.na(v1) | is.nan(v1)), 
                                       which(is.na(v2) | is.nan(v2))))
                if(length(vals_to_rm) > 0) {
                  v1 <- v1[-vals_to_rm]
                  v2 <- v2[-vals_to_rm]
                  q1 <- q1[-vals_to_rm]
                  q2 <- q2[-vals_to_rm]
                }
                w <- q1 * q2
                #w <- pmin(q1, q2, na.rm=T)
                w <- unlist(lapply(w, function(l) 1-l))
                
                # Use the WCorr package
                spearman <- weightedCorr(v1, v2, method = "Spearman", 
                                         weights = w)
                vals <- c(vals, spearman)
              }
              ## DOT PRODUCT ##
              if (method == "DP") {
                dp <- as.numeric(sum(v1*v2, na.rm = T))
                vals <- c(vals, dp)
                
                ## COSINE SIMILARITY ##
              } else {
                vals_to_rm <- unique(c(which(is.na(v1)), which(is.na(v2))))
                if(length(vals_to_rm) > 0) {
                  v1 <- v1[-vals_to_rm]
                  v2 <- v2[-vals_to_rm]
                }
                cos <- cosine(v1, v2)
                vals <- c(vals, cos)
              }
            }
          }
          return(vals)
        }
      } else {return(NA)}
    } else{return(NA)}
  })
  names(vals_per_driver) <- drivers
  vals_per_driver <- vals_per_driver[!is.na(vals_per_driver)]
  
  return(vals_per_driver)
}

#' Calculates the Hypergeometric enrichment of the pairwise sets of mutation
#' coefficients for each driver gene
#' @param perCancer list of output data frames from Dyscovr per-cancer, named by
#' the cancer type
#' @param qval_thres a q-value threshold for significance
#' @param method either 'Spearman', 'WeightedSpearman', 'DP', 'CS', or 'HG' to 
#' indicate what kind of test we are performing
#' @param drivers a vector of the driver gene names that are highly mutated in
#' at least 3 cancer types/ that we'd like to consider
calculate_hg_per_driver <- function(perCancer, qval_thres, method, drivers) {
  
  hg_per_driver <- lapply(drivers, function(d) {
    targets_by_ct <- lapply(1:length(perCancer), function(i) {
      x <- perCancer[[i]]
      if(!is.na(qval_thres)) {
        if(nrow(x[(x$q.value < qval_thres) & (x$R_i.name == d),]) > 0) {
          return(as.character(unlist(x[(x$q.value < qval_thres) & (x$R_i.name == d), 
                                       'T_k.name', with=F])))
        } else {return(NA)}
      } else {return(as.character(unlist(x[x$R_i.name == d, 'T_k.name', with=F])))}
    })
    names(targets_by_ct) <- names(perCancer)
    targets_by_ct <- targets_by_ct[!is.na(targets_by_ct)]
    targets_by_ct <- targets_by_ct[lengths(targets_by_ct) != 0]
    
    if(length(targets_by_ct) > 1) {
      # Take the pairwise HG of each set of cancer types
      hg <- c()
      for (i in 1:(length(targets_by_ct)-1)) {
        len_hits1 <- length(targets_by_ct[[i]])
        ct1 <- names(targets_by_ct)[i]
        dfi <- perCancer[[ct1]]
        
        for (j in (i+1):length(targets_by_ct)) {
          # Calculate whether overlap is significant using an approx. of an 
          # exact hypergeometric test
          len_hits2 <- length(targets_by_ct[[j]])
          target_overlap <- intersect(targets_by_ct[[i]], 
                                      targets_by_ct[[j]])
          len_ol <- length(target_overlap)
          
          ct2 <- names(targets_by_ct)[j]
          dfj <- perCancer[[ct2]]
          n <- length(intersect(as.character(unlist(dfi[dfi$R_i.name == d, 
                                                        'T_k.name'])), 
                                as.character(unlist(dfj[dfj$R_i.name == d, 
                                                        'T_k.name']))))
          print(n)
          
          # n is size of urn; hits1 is the number of genes from set 1. Say we 
          # were to draw the number of balls the size of hits2 (# of genes from 
          # set 2), w/o replacement. M of these balls are also found in hits2. 
          # What are the odds that M >= len_ol?
          phyper_res <- phyper(len_ol-1, len_hits1, n-len_hits1, len_hits2, 
                               lower.tail = F, log.p = F)
          print(paste("p-value:", phyper_res))
          hg <- c(hg, phyper_res)
        }
      }
      return(hg)
    } else {return(NA)}
  })
  names(hg_per_driver) <- drivers
  hg_per_driver <- hg_per_driver[!is.na(hg_per_driver)]
  
  return(hg_per_driver)
}

# Call function
top5_drivers <- c("TP53", "PIK3CA", "KRAS", "CTNNB1", "FBXW7")
create_perdriver_percancer_barplot(perCancer, NA, 'Spearman', top5_drivers)


# PART C: SEE CODE FOR FIGURE 3, PART B