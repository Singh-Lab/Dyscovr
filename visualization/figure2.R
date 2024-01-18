############################################################
# Code to Create Figure 2 Visualizations
# Written by Sara Geraghty
# PUBLICATION INFORMATION
############################################################

library(data.table)
library(ggplot2)
library(UpSetR)
library(Hmisc)
library(reshape2)
library(wCorr)
library(broom)

# Local PATH to directory containing Dyscovr output files
PATH <-  paste0(getwd(), "Output/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
### IMPORT PER-CANCER OUTPUT FILE(S)
############################################################
perCancer_fns <- intersect(list.files(path = PATH, recursive = T,
                                         pattern = "_corrected_MUT"), 
                              intersect(list.files(path = PATH, recursive = T,
                                                   pattern = "allGenes_"),
                                        list.files(path = PATH, recursive = T,
                                                   pattern = "Nonsyn.Drivers.Vogel.elim.vif.5")))
perCancer <- lapply(perCancer_fns, function(f) 
  fread(paste0(PATH, f), header = T))
names(perCancer) <- unlist(lapply(pc_perCancer, function(x)
  unlist(strsplit(x, "/", fixed = T))[1]))


############################################################
### PART A: PERCENTAGE BARPLOT OF DRIVER HITS PER-CANCER
############################################################
#' Create driver bar with the percentage of hits for each driver gene in each
#' cancer type, with number of hits at the particular significance threshold 
#' labeled
#' @param list_of_master_dfs list of output master DFs with target gene names 
#' and q-values, names are cancer types or subtypes
#' @param tophit_thres a percentage or q-value threshold (e.g. 0.05) within 
#' which to consider a hit "significant"
#' @param perc_or_qval_or_ss whether the threshold is a percentage, a q-value or 
#' a stability score
#' @param all_genes_id_conv a conversion table to convert R_i Uniprot IDs to 
#' intelligible gene names
create_driver_bar <- function(list_of_master_dfs, tophit_thres, 
                              perc_or_qval_or_ss, all_genes_id_conv) {
  
  freq.tables <- lapply(1:length(list_of_master_dfs), function(i) {
    master_df <- list_of_master_dfs[[i]]
    cancer_type <- names(list_of_master_dfs)[[i]]
    
    # Create the frequency table
    freq.table <- generate_frequency_input_table(master_df, tophit_thres, 
                                                 perc_or_qval_or_ss, 
                                                 all_genes_id_conv)
    if(length(freq.table) < 2) {return(NA)}
    
    # Add a column of percentages
    sum <- sum(freq.table$Frequency)
    if(sum < 10) {return(NA)}
    freq.table$Percentage <- unlist(lapply(freq.table$Frequency, function(x) 
      (x/sum)*100))
    
    # Add the cancer type
    freq.table$Cancer_Type <- cancer_type
  
    return(freq.table)
  })
  freq.tables <- freq.tables[!is.na(freq.tables)]
  freq.table <- do.call(rbind, freq.tables)
  freq.table[is.na(freq.table)] <- 0

  # Make the drivers that make up less than 5% of the results across all cancer 
  # types "other"
  unique_drivers <- unique(as.character(freq.table$Driver))
  drivers_to_rm <- unlist(lapply(unique_drivers, function(d) {
    # Get the total % across all cancer types
    total_perc <- sum(freq.table[freq.table$Driver == d, 'Percentage'])
    total_frac <- sum(freq.table[freq.table$Driver == d, 'Frequency'])
    if((total_perc < 10) & (total_frac < 100)) {return(d)}
    else {return(NA)}
  }))
  drivers_to_rm <- drivers_to_rm[!is.na(drivers_to_rm)]
  freq.table$Driver <- unlist(lapply(freq.table$Driver, as.character))
  freq.table[which(freq.table$Driver %fin% drivers_to_rm), 'Driver'] <- "Other"
  freq.table$Driver <- as.factor(freq.table$Driver)

  # Adjust the order of the drivers based on total frequency + # of cancer types
  freq.per.driver <- lapply(unique(freq.table$Driver), function(d) 
    sum(freq.table[freq.table$Driver == d, 'Frequency']))
  names(freq.per.driver) <- unique(freq.table$Driver)
  freq.per.driver.df <- melt(as.data.frame(freq.per.driver))
  colnames(freq.per.driver.df) <- c("Driver", "Total.Frequency")
  
  # Add the number of cancer types as well
  freq.per.driver.df$Num.CTs <- unlist(lapply(freq.per.driver.df$Driver, function(d)
    length(unique(freq.table[freq.table$Driver == d, 'Cancer_Type']))))
  
  freq.per.driver.df <- freq.per.driver.df[with(freq.per.driver.df, 
                                                order(Num.CTs, Total.Frequency, 
                                                      decreasing = T)),]
  freq.per.driver.df$Rank <- 1:nrow(freq.per.driver.df)
  freq.per.driver.df[freq.per.driver.df$Driver == 'Other', 'Rank'] <- max(
    freq.per.driver.df$Rank) + 1
  freq.per.driver.df <-  freq.per.driver.df[order(freq.per.driver.df$Rank),]

  freq.table <- merge(freq.table, freq.per.driver.df, by = "Driver", all = T)
  freq.table <- freq.table[order(freq.table$Rank),]
  
  freq.table$Frequency.Label <- unlist(lapply(1:nrow(freq.table), function(i) {
    perc <- freq.table[i, "Percentage"]
    freq <- freq.table[i, "Frequency"]
    return(ifelse(perc > 2.5, freq, ""))
  }))
  freq.table$Driver <- as.factor(freq.table$Driver)
  freq.table$Cancer_Type <- as.factor(freq.table$Cancer_Type)
  
  # Add the rank of each cancer type based on the total frequency of the top drivers
  drivers <- unique(as.character(freq.table$Driver))
  remaining_cancer_types <- unique(as.character(freq.table$Cancer_Type))
  ranks <- list()
  for(d in drivers) {
    print(remaining_cancer_types)
    if(length(remaining_cancer_types) == 0) {break}
    sub <- freq.table[(freq.table$Driver == d) & 
                        (freq.table$Cancer_Type %fin% remaining_cancer_types),]
    new_ranks <- c()
    if(length(ranks) > 0) {
      new_ranks <- rank(desc(sub$Percentage)) + max(as.numeric(ranks))
    } else {new_ranks <- rank(desc(sub$Percentage))}
    names(new_ranks) <- sub$Cancer_Type
    ranks <- c(ranks, new_ranks)
    remaining_cancer_types <- setdiff(remaining_cancer_types, sub$Cancer_Type)
  }
  freq.table$rank_from_driver <- unlist(lapply(freq.table$Cancer_Type, function(ct)
    ranks[names(ranks) == ct]))
  
  # Create plot
  g <- ggplot(freq.table, aes(x = fct_reorder(Cancer_Type, rank_from_driver), 
                         y = Percentage, 
                         fill = fct_reorder(Driver, Rank, .desc = T), 
                         label = Frequency.Label)) + 
    geom_bar(position="stack", stat="identity", color = "black") +
    theme_minimal() + xlab("Cancer Type") + 
    ylab(paste0("Percentage of Hits q<", tophit_thres)) +
    scale_fill_manual(values = distinct_colors,
                      breaks = as.character(freq.per.driver.df$Driver),
                      name = "Driver") +
    geom_text(size = 3, position = position_stack(vjust = 0.5)) +
    theme(axis.title = element_text(face = "bold", size = 14), 
          axis.text = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12))
  g
}

#' Helper function to generate the input frequency table for drivers, given 
#' appropriate significance thresholds
#' @param master_df output master DF with target gene names and q-values
#' @param tophit_thres a percentage or q-value threshold (e.g. 0.05) within 
#' which to  consider a hit "significant" 
#' @param perc_or_qval_or_ss whether the threshold is a percentage, a q-value,
#' or a stability score
#' @param all_genes_id_conv a conversion table to convert R_i Uniprot IDs to 
#' intelligible gene names
generate_frequency_input_table <- function(master_df, tophit_thres, 
                                           perc_or_qval_or_ss, 
                                           all_genes_id_conv) {
  if(perc_or_qval_or_ss == "perc") {
    # Get the top n% of hits
    master_df_topn <- master_df[1:(tophit_thres * nrow(master_df)),]
    if('q.value' %in% colnames(master_df)) {
      print(paste("Q-Value at", 
                  paste(tophit_thres, 
                        paste("threshold:", master_df[tophit_thres*nrow(master_df), 
                                                      'q.value']))))
    }
    if('stability.score' %in% colnames(master_df)) {
      print(paste("Stability score at", 
                  paste(tophit_thres, 
                        paste("threshold:", master_df[tophit_thres*nrow(master_df), 
                                                      'stability.score']))))
    }
  } else if (perc_or_qval_or_ss == "qval") {
    # Get the n hits below the qvalue threshold
    master_df_topn <- master_df[master_df$q.value < tophit_thres,]
    print(paste("Number of hits below the q-value threshold:", 
                nrow(master_df_topn)))
  } else if (perc_or_qval_or_ss == "ss") {
    # Get the n hits greater than or equal to the stability score threshold
    master_df_topn <- master_df[master_df$stability.score >= tophit_thres,]
    print(paste("Number of hits above or equal to threshold:", 
                nrow(master_df_topn)))
  } else {
    print("The only possible thresholding options are perc or qval. 
          Please try again.")
    return(NA)
  }
  
  # Get the unique drivers, and create a frequency table of the number of 
  # times that each driver appears
  if(nrow(master_df_topn) > 0) {
    unique_drivers <- unique(master_df_topn$term)
    unique_drivers <- unique_drivers[!is.na(unique_drivers)]
    
    freq.table <- as.data.frame(lapply(unique_drivers, function(x)
      nrow(master_df_topn[master_df_topn$term == x,])))
    unique_driver_names <- unlist(lapply(unique_drivers, function(d) {
      driver_uniprot <- unlist(strsplit(d, "_", fixed = TRUE))[1]
      driver_gn <- unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == 
                                              driver_uniprot, 
                                            'external_gene_name'])
      return(driver_gn)
    }))
    names(freq.table) <- unique_driver_names
    
    freq.table.m <- melt(freq.table)
    colnames(freq.table.m) <- c("Driver", "Frequency")
    freq.table.m <- freq.table.m[order(freq.table.m$Frequency, decreasing = TRUE),]
    
    return(freq.table.m)
  }
  else {return(NA)}
}

# Call function
create_driver_bar(perCancer, 0.2, "qval", all_genes_id_conv)


############################################################
### PART B: UpSet Plots per Driver Gene
############################################################
#' Create an UpSet plot for a given gene-of-interest (goi) to visualize the
#' shared hits between various cancer types
#' @param list_of_master_dfs a named list of output DFs from Dyscovr, where 
#' names correspond to cancer types 
#' @param goi the Uniprot ID of the driver gene-of-interest
#' @param qval_thres q-value threshold for significance
#' @param top_n an optional field to limit to only the top n hits 
plot_overlap_between_subgroups_upset <- function(list_of_master_dfs, goi, 
                                                 qval_thres, top_n) {
  
  sig_hits <- lapply(1:length(list_of_master_dfs), function(i) {
    m <- list_of_master_dfs[[i]]
    if(!(goi %fin% m$R_i)) {return(NA)}
    
    sig_hits_tmp <- as.character(unlist(m[(m$R_i == goi) & 
                                            (m$q.value < qval_thres), 
                                          'T_k.name']))
    if(!is.na(top_n)) {
      if(length(sig_hits_tmp) > top_n) {sig_hits_tmp <- sig_hits_tmp[1:top_n]}
    }
    return(sig_hits_tmp)
  })
  names(sig_hits) <- names(list_of_master_dfs)
  sig_hits <- sig_hits[!is.na(sig_hits)]
  sig_hits <- sig_hits[lengths(sig_hits) != 0]
  
  ylab <- paste(goi, paste("Hit Intersections, q <", qval_thres))
  
  # Software is limited to only 5 groups, to restrict to just the cancer types
  # with the greatest absolute number of hits
  sig_hits_top5 <- sig_hits[order(lengths(sig_hits), decreasing = T)][1:5]
  
  # Set color scheme
  main_bar_col <- c("#0072B5FF")
  sets_bar_col <- c("#FFDC91FF")
  matrix_col <- c("#E18727FF")
  shade_col <- c("wheat4")
  
  if(!is.na(top_n)) {ylab <- paste0(ylab, 
                                    paste(", top", 
                                          paste(top_n, "hits per cancer type")))}
  # Create UpSet plot
  plt <- upset(fromList(sig_hits_top5),
               order.by = 'freq', point.size = 3.5, line.size = 2, 
               mainbar.y.label = paste0("Intersection Size, q<", qval_thres), 
               sets.x.label = "Number of Hits", 
               text.scale = c(2, 2, 2, 2, 2, 2),
               main.bar.color = main_bar_col,
               sets.bar.color = sets_bar_col,
               matrix.color = matrix_col,
               shade.color = shade_col)
  print(plt)
}

# Call function
plot_overlap_between_subgroups_upset(perCancer, "P04637", 0.2, NA) # TP53
plot_overlap_between_subgroups_upset(perCancer, "P42336", 0.2, NA) # PIK3CA
plot_overlap_between_subgroups_upset(perCancer, "P01116", 0.2, NA) # KRAS
plot_overlap_between_subgroups_upset(perCancer, "O75874", 0.2, NA) # IDH1

# NOTE: Annotations of shared target genes were added manually in a graphics editor


############################################################
### PART C: Pairwise Weighted Spearman Correlation Barplot
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
create_perdriver_percancer_barplot(perCancer, 0.2, 'WeightedSpearman', top5_drivers)


############################################################
### PART D: Spearman Correlation between TCGA-BRCA and METABRIC
############################################################   
#' Create a Spearman barplot that, when given two corresponding lists of master 
#' DFs (named according to the R_i), computes the Spearman correlation of the 
#' mutation coefficient (Beta) values. Plots a barplot of the results.
#' @param source1_master_list a list of master DFs from a certain source 1
#' @param source2_master_list a list of master DFs (for the same genes) from a 
#' certain source 2. Two source lists should be the same length.
#' @param qval_thres optional q-value threshold for significance
#' @param label a character label of the data sets for plotting
create_spearman_barplot <- function(source1_master_list, source2_master_list, 
                                    qval_thres, label) {
  
  spearman_correlations <- lapply(1:length(source1_master_list), function(i) {
    df1 <- source1_master_list[[i]]
    df2 <- source2_master_list[[i]]
    
    if(!(is.na(qval_thres))) {
      df1 <- df1[df1$q.value < qval_thres,]
      df2 <- df2[df2$q.value < qval_thres,]
    }
    
    df_merged <- merge(df1, df2, by = "T_k.name", all = FALSE)
    
    betas_df1 <- df_merged$estimate.x
    betas_df2 <- df_merged$estimate.y
    
    neglogpvals_df1 <- -log2(df_merged$p.value.x)
    neglogpvals_df2 <- -log2(df_merged$p.value.y)
    
    neglogpvals_wdir_df1 <- unlist(lapply(1:length(neglogpvals_df1), function(i) {
      beta <- betas_df1[i]
      pval <- neglogpvals_df1[i]
      return(ifelse(beta < 0, (-1)*pval, pval))
    }))
    neglogpvals_wdir_df2 <- unlist(lapply(1:length(neglogpvals_df2), function(i) {
      beta <- betas_df2[i]
      pval <- neglogpvals_df2[i]
      return(ifelse(beta < 0, (-1)*pval, pval))
    }))
    
    #spearman <- tidy(cor.test(neglogpvals_wdir_df1, neglogpvals_wdir_df2, 
    #                          method = "spearman", use = "pairwise"))
    # Try a weighted Spearman
    w <- df_merged$p.value.x * df_merged$p.value.y
    #w <- pmin(df_merged$p.value.x, df_merged$p.value.y, na.rm=T)
    w <- unlist(lapply(w, function(l) 1-l))
    
    # Use the WCorr package
    spearman <- weightedCorr(betas_df1, betas_df2, method = "Spearman", 
                             weights = w)
    
    return(spearman)
  })
  spearman_cor_vals <- as.numeric(unlist(spearman_correlations))
  #spearman_cor_vals <- as.numeric(unlist(lapply(spearman_correlations, function(x) 
    #x$estimate)))
  spearman_correlations_df <- data.frame("Gene" = names(source1_master_list), 
                                         "Spearman" = spearman_cor_vals)
  #print("P-Values")
  #print(unlist(lapply(spearman_correlations, function(x) 
    #as.numeric(x$p.value))))
  #print("Correlations")
  #print(as.numeric(unlist(lapply(spearman_correlations, function(x) 
    #x$estimate))))
  
  ylab <- paste("Weighted SCC,", label)
  if(!is.na(qval_thres)) {
    ylab <- paste(ylab, paste0("(q <", paste0(qval_thres, ")")))
  }
  
  p <- ggplot(spearman_correlations_df, aes(x = reorder(Gene, -Spearman, mean), 
                                            y = Spearman, fill = Gene)) + 
    geom_bar(stat = "identity", show.legend = FALSE) + scale_fill_nejm() + 
    theme_minimal() + xlab("Gene") + ylab(ylab) +
    theme(axis.text = element_text(face="bold", size = 12), 
          axis.title=element_text(size=14,face="bold"))
  
  print(p)
}

tcga_brca <- perCancer[["BRCA"]]
source1_master_list <- list("TP53" = tcga_brca[tcga_brca$R_i.name == "TP53",], 
                            "PIK3CA" = tcga_brca[tcga_brca$R_i.name == "PIK3CA",])

# Import METABRIC results
metabric_res <- read.csv(paste0(PATH, "METABRIC/res_top_0.05_allGenes_df_rawCNA_methMRaw_Nonsyn.Driver.Vogel.sklearn.elim.vif.5_corrected_MUT.csv"))

source2_master_list <- list("TP53" = metabric_res[metabric_res$R_i.name == "TP53",], 
                            "PIK3CA" = metabric_res[metabric_res$R_i.name == "PIK3CA",])

# Call function
create_spearman_barplot(source1_master_list, source2_master_list, 0.2, 
                        "TCGA-BRCA & METABRIC") 
