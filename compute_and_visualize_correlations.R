############################################################
### COMPUTE AND VISUALIZE CORRELATIONS
### Written By: Sara Geraghty, July 2022
############################################################

# Includes use and visualization of Spearman and Pearson correlations
# to examine the output from linear_model.R


# NEJM color palatte: https://nanx.me/ggsci/reference/pal_nejm.html

library(ggplot2)
library(stringr)
library(dplyr)
library(VennDiagram)
library("RColorBrewer")

# Path to output files
main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

# Path to where output figures should be saved
output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/BRCA/Linear Model/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/Pan-Cancer/Linear Model/"


# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


############################################################
############################################################
#### SPEARMAN CORRELATION OF BETAS AND T-STATISTICS
############################################################
############################################################
#' Compute and print the Spearman correlation of the Betas and of the T-statistic,
#' given two groups of interest
#' @param results_table the output master DF from the linear model
#' @param ri_1 the external gene name of the first protein of interest
#' @param ri_2 the external gene name of the second protein of interest
compute_and_print_spearman <- function(results_table, ri_1, ri_2) {
  if((ri_1 %fin% results_table$R_i.name) & (ri_2 %fin% results_table$R_i.name)) {
    
    target_genes <- setdiff(unique(results_table$T_k.name), c(ri_1, ri_2))
    
    # Mini functions to get Betas/t-statistics for each target gene
    #' @param results_table the master DF 
    #' @param target_genes the list of target genes
    #' @param ri the given regulatory protein
    #' @param type "Betas" or "t-statistics" to indicate what value we are returning
    get_values <- function(results_table, target_genes, ri, type) {
      vals <- unlist(lapply(target_genes, function(tg) {
        if(type == "Betas") {
          est <- results_table[(results_table$R_i.name == ri) & (results_table$T_k.name == tg), "estimate"]
        } else {
          est <- results_table[(results_table$R_i.name == ri) & (results_table$T_k.name == tg), "statistic"]
        }
        if (length(est) == 0) {est <- 0}  # To ensure the lengths of the Beta vectors are the same
        if(length(est) > 1) {est <- est[1]}
        return(est)
      }))
    }
    
    grp1_Betas <- get_values(results_table, target_genes, ri_1, "Betas")
    grp2_Betas <- get_values(results_table, target_genes, ri_2, "Betas")
    
    print(length(grp1_Betas))
    print(length(grp2_Betas))
    
    # Get Betas spearman
    betas_spearman <- cor.test(grp1_Betas, grp2_Betas, method = "spearman")
    betas_spearman_stat <- as.numeric(betas_spearman$estimate)
    betas_spearman_pval <- betas_spearman$p.value
    
    # Print the results
    print(paste("Spearman results for", paste(ri_1, paste("and", ri_2))))
    print(paste("Beta correlation of", paste(betas_spearman_stat, paste(", p-value of", betas_spearman_pval))))
    
    # Create a plot to visualize the correlations
    plot(grp1_Betas, grp2_Betas, pch = 19, col = alpha("lightblue", 0.4), #main = "Betas Spearman Correlation",
         xlab = paste(ri_1, "Betas"), ylab = paste(ri_2, "Betas"))
    abline(lm(grp2_Betas ~ grp1_Betas), col = "red", lwd = 3)
    text(labels = paste("Correlation:", paste(round(betas_spearman_stat, 6), 
                                              paste(", p-value:", round(betas_spearman_pval, 6)))), 
         x = max(grp1_Betas, na.rm = TRUE)-sd(grp1_Betas, na.rm = TRUE)*3, 
         y = max(grp2_Betas, na.rm = TRUE)-sd(grp2_Betas, na.rm = TRUE), col = "black")
    
    
  } else {print("Error. Provided regprots are not in the given master DF.")}
}

ri_1 <- "TP53"
ri_2 <- "PIK3CA"

# Call function
compute_and_print_spearman(master_df, ri_1, ri_2)



#' Compute and print Spearman, given two separate results tables
#' @param results_table1 first output master DF from the linear model
#' @param results_table2 second output master DF from the linear model
#' @param ri_1 the external gene name of the first protein of interest
#' @param ri_2 the external gene name of the second protein of interest
compute_and_print_spearman_multDF <- function(results_table1, results_table2, ri_1, ri_2) {
  if((ri_1 %fin% results_table1$R_i.name) & (ri_2 %fin% results_table2$R_i.name)) {
    
    target_genes <- intersect(unique(results_table1$T_k.name), unique(results_table2$T_k.name))
    
    # Mini functions to get Betas/t-statistics for each target gene
    #' @param results_table the master DF 
    #' @param target_genes the list of target genes
    #' @param ri the given regulatory protein
    #' @param type "Betas" or "t-statistics" to indicate what value we are returning
    get_values <- function(results_table, target_genes, ri, type) {
      vals <- unlist(lapply(target_genes, function(tg) {
        if(type == "Betas") {
          est <- results_table[(results_table$R_i.name == ri) & (results_table$T_k.name == tg), "estimate"]
        } else {
          est <- results_table[(results_table$R_i.name == ri) & (results_table$T_k.name == tg), "statistic"]
        }
        if (length(est) == 0) {est <- 0}  # To ensure the lengths of the Beta vectors are the same
        return(est)
      }))
    }
    
    grp1_Betas <- get_values(results_table1, target_genes, ri_1, "Betas")
    grp2_Betas <- get_values(results_table2, target_genes, ri_2, "Betas")
    
    # Get Betas spearman
    betas_spearman <- cor.test(grp1_Betas, grp2_Betas, method = "spearman")
    betas_spearman_stat <- as.numeric(betas_spearman$estimate)
    betas_spearman_pval <- betas_spearman$p.value
    
    # Print the results
    print(paste("Spearman results for", paste(ri_1, paste("and", ri_2))))
    print(paste("Beta correlation of", paste(betas_spearman_stat, paste(", p-value of", betas_spearman_pval))))
    
    # Create a plot to visualize the correlations
    # "#BC3C29FF", "#0072B5FF"
    plot(grp1_Betas, grp2_Betas, pch = 19, col = alpha("lightblue", 0.4), #main = "Betas Spearman Correlation",
         xlab = "", ylab = "", bty="n")
    abline(lm(grp2_Betas ~ grp1_Betas), col = "red", lwd = 3)
    text(labels = paste("Correlation:", paste(round(betas_spearman_stat, 6), 
                                              paste(", p-value:", round(betas_spearman_pval, 6)))), 
         x = max(grp1_Betas, na.rm = TRUE)-sd(grp1_Betas, na.rm = TRUE)*3, 
         y = max(grp2_Betas, na.rm = TRUE)-sd(grp2_Betas, na.rm = TRUE), 
         col = "black", font=2)
    mtext(side=1, line=3, paste(ri_1, "Betas"), font=2, cex=1.2)
    mtext(side=2, line=3, paste(ri_2, "Betas"), font=2, cex=1.2)
    
    
    
  } else {print("Error. Provided regprots are not in the given master DF.")}
}

ri_1 <- "TP53"
ri_2 <- "PIK3CA"

# Call function
compute_and_print_spearman_multDF(master_df1, master_df2, ri_1, ri_2)


#' Compute and print the Spearman correlation of the T-statistics for two separate
#' subpopulations
#' @param results_table_subpop1 the output master DF from the linear model for 
#' the first subgroup
#' #' @param results_table_subpop2 the output master DF from the linear model for 
#' the second subgroup
#' @param ri the external gene name of the protein of interest
#' @param subpop1 the name of the first subgroup
#' @param subpop2 the name of the second subgroup
#' @param sub_to_top if given, subsets to the genes with smallest X and largest X
#' t-statistics from the two subpopulations
compute_and_print_spearman_subpops <- function(results_table_subpop1, results_table_subpop2, 
                                               ri, subpop1, subpop2, sub_to_top) {
  if((ri %fin% results_table_subpop1$R_i.name) & (ri %fin% results_table_subpop2$R_i.name)) {
    
    if(is.infinite(sub_to_top)) {
      target_genes <- intersect(unique(results_table_subpop1$T_k.name), 
                                unique(results_table_subpop2$T_k.name))
    } else {
      # Use a helper function to subset the target genes
      target_genes <- subset_targ_genes(results_table_subpop1, results_table_subpop2, 
                                        ri, sub_to_top)
    }
    
    # Mini functions to get t-statistics for each target gene
    #' @param results_table the master DF 
    #' @param target_genes the list of target genes
    #' @param ri the given regulatory protein
    get_values <- function(results_table, target_genes, ri) {
      vals <- unlist(lapply(target_genes, function(tg) {
        est <- results_table[(results_table$R_i.name == ri) & (results_table$T_k.name == tg), "statistic"]
        if (length(est) == 0) {est <- 0}  # To ensure the lengths of the Beta vectors are the same
        return(est)
      }))
    }
    
    subpop1_tstats <- get_values(results_table_subpop1, target_genes, ri)
    subpop2_tstats <- get_values(results_table_subpop2, target_genes, ri)
    
    # Get Spearman
    spearman <- cor.test(subpop1_tstats, subpop2_tstats, method = "spearman")
    spearman_stat <- as.numeric(spearman$estimate)
    spearman_pval <- spearman$p.value
    
    # Print the results
    print(paste("Spearman results for", paste(subpop1, paste("and", paste(subpop2, paste("in", ri))))))
    print(paste("T-statistic correlation of", paste(spearman_stat, paste(", p-value of", spearman_pval))))
    
    # Create a plot to visualize the correlations
    plot(subpop1_tstats, subpop2_tstats, pch = 19, col = alpha("lightblue", 0.4), #main = "Betas Spearman Correlation",
         xlab = paste(subpop1, paste(ri, "T-statistics")), ylab = paste(subpop2, paste(ri, "T-statistics")))
    abline(lm(subpop1_tstats ~ subpop2_tstats), col = "red", lwd = 3)
    text(labels = paste("Correlation:", paste(round(spearman_stat, 6), 
                                              paste(", p-value:", round(spearman_pval, 6)))), 
         x = max(subpop1_tstats, na.rm = TRUE)-sd(subpop1_tstats, na.rm = TRUE)*3, 
         y = max(subpop2_tstats, na.rm = TRUE)-sd(subpop2_tstats, na.rm = TRUE), col = "black")
    
    
  } else {print("Error. Provided regprot not in the given master DF.")}
}

#' Helper function for subsetting the target genes to the top X/ bottom X t-statistics 
#' in the two subpopulations
#' @param results_table_subpop1 the output master DF from the linear model for 
#' the first subgroup
#' @param results_table_subpop2 the output master DF from the linear model for 
#' the second subgroup
#' @param ri the external gene name of the protein of interest
#' @param sub_to_top if given, subsets to the genes with smallest X and largest X
#' t-statistics from the two subpopulations
subset_targ_genes <- function(results_table_subpop1, results_table_subpop2, ri, sub_to_pop) {
  
  # First, subset to the regprot of interest
  results_table_subpop1_sub <- results_table_subpop1[results_table_subpop1$R_i.name == ri,]
  results_table_subpop2_sub <- results_table_subpop2[results_table_subpop2$R_i.name == ri,]
  
  # Next, get the top X and bottom X t-statistics 
  results_table_subpop1_sub <- arrange(results_table_subpop1_sub, desc(estimate))
  results_table_subpop2_sub <- arrange(results_table_subpop2_sub, desc(estimate))
  
  top100_subpop1 <- results_table_subpop1_sub[1:sub_to_pop, 'T_k.name']
  bottom100_subpop1 <- results_table_subpop1_sub[(nrow(results_table_subpop1_sub)-sub_to_pop):nrow(results_table_subpop1_sub), 
                                                 'T_k.name']
  top100_subpop2 <- results_table_subpop2_sub[1:sub_to_pop, 'T_k.name']
  bottom100_subpop2 <- results_table_subpop2_sub[(nrow(results_table_subpop2_sub)-sub_to_pop):nrow(results_table_subpop2_sub), 
                                                 'T_k.name']
  
  top_genes <- unique(c(top100_subpop1, bottom100_subpop1, top100_subpop2, bottom100_subpop2))
  
  return(top_genes)
}

ri <- "TP53"
ri <- "PIK3CA"

# Call function
compute_and_print_spearman_subpops(lumA_allgenes_mut, lumB_allgenes_mut, ri, 
                                   "LumA", "LumB", 100)


############################################################
############################################################
#### CREATE A BARPLOT TO COMPARE MULTIPLE SPEARMAN CORRELATIONS
############################################################
############################################################
#' Create a Spearman barplot that, when given two corresponding lists of master DFs
#' (named according to the R_i), computes the Spearman correlation of the Beta values.
#' Plots a barplot of the Beta values.
#' @param source1_master_list a list of master DFs from a certain source 1
#' @param source2_master_list a list of master DFs (for the same genes) from a certain
#' source 2. Two source lists should be the same length.
create_spearman_barplot <- function(source1_master_list, source2_master_list) {
  
  spearman_correlations <- lapply(1:length(source1_master_list), function(i) {
    df1 <- source1_master_list[[i]]
    df2 <- source2_master_list[[i]]
    
    df_merged <- merge(df1, df2, by = "T_k.name", all = FALSE)
    
    betas_df1 <- df_merged$estimate.x
    betas_df2 <- df_merged$estimate.y
    
    spearman <- tidy(cor.test(betas_df1, betas_df2, method = "spearman", use = "pairwise"))
    
    return(spearman)
  })
  spearman_cor_vals <- as.numeric(unlist(lapply(spearman_correlations, function(x) x$estimate)))
  spearman_correlations_df <- data.frame("Gene" = names(source1_master_list), 
                                         "Spearman" = spearman_cor_vals)
  print("P-Values")
  print(unlist(lapply(spearman_correlations, function(x) as.numeric(x$p.value))))
  print("Correlations")
  print(as.numeric(unlist(lapply(spearman_correlations, function(x) x$estimate))))
  
  p <- ggplot(spearman_correlations_df, aes(x = reorder(Gene, -Spearman, mean), 
                                            y = Spearman, fill = Gene)) + 
    geom_bar(stat = "identity", show.legend = FALSE) + scale_color_nejm() + 
    theme_minimal() + xlab("\nGene") + ylab("Spearman Correlation\n") +
    theme(axis.text = element_text(face="bold", size = 12), 
          axis.title=element_text(size=14,face="bold"))
  
  return(p)
}

source1_master_list <- list("TP53" = allgenes_p53_incl2GenoPCs, "PIK3CA" = allgenes_pik3ca_incl2GenoPCs,
                            "KMT2C" = allgenes_kmt2c, "GATA3" = allgenes_gata3)
source2_master_list <- list("TP53" = allgenes_p53_metabric, "PIK3CA" = allgenes_pik3ca_metabric,
                            "KMT2C" = allgenes_kmt2c_metabric, "GATA3" = allgenes_gata3_metabric)

create_spearman_barplot(source1_master_list, source2_master_list)


############################################################
############################################################
#### CALCULATE PEARSON'S CORRELATION COEFFICIENT (PCC) 
#### BETWEEN THE SIGNIFICANT HITS FOR DIFFERENT RUNS
############################################################
############################################################
#' Given two sets of top hit genes from two separate master DFs,
#' computes, prints, and returns the PCC value
#' @param master_df1 the first master DF
#' @param master_df2 the second master DF
#' @param qval_thres a threshold for the q-value, above which values are 
#' considered significant
compute_pcc <- function(master_df1, master_df2, qval_thres) {
  
  # Subset to only significant hits
  master_df1_sig <- master_df1[master_df1$q.value < qval_thres, 'T_k.name']
  master_df2_sig <- master_df2[master_df2$q.value < qval_thres, 'T_k.name']
  
  rank_df_master1 <- data.frame("Gene" = master_df1_sig, "Rank.Master1" = 1:length(master_df1_sig))
  rank_df_master2 <- data.frame("Gene" = master_df2_sig, "Rank.Master2" = 1:length(master_df2_sig))
  
  print(head(rank_df_master1))
  
  rank_df <- merge(rank_df_master1, rank_df_master2, by = "Gene")
  
  print(head(rank_df))
  
  pcc <- cor.test(rank_df$Rank.Master1, rank_df$Rank.Master2, method = "pearson", use = "pairwise")
  print(pcc)
  return(pcc$p.value)
  
  # Get a rank-based representation of each value
  #unique_targs <- unique(c(master_df1$T_k.name, master_df2$T_k.name))
  #list1_ranking <- lapply(unique_targs, function(t) which(master_df1$T_k.name == t))
  #list2_ranking <- lapply(unique_targs, function(t) which(master_df2$T_k.name == t))
  #names(list1_ranking) <- unique_targs
  #names(list2_ranking) <- unique_targs
  
  #print(head(list1_ranking))
  #print(head(list2_ranking))
  
  #print(head(as.numeric(unlist(list1_ranking, use.names = FALSE))))
  
  #pcc <- cor(as.numeric(unlist(list1_ranking, use.names = FALSE)), 
  #as.numeric(unlist(list2_ranking, use.names = FALSE)), 
  #method = "pearson")
  #print(pcc)
  #return(pcc$p.value)
}

