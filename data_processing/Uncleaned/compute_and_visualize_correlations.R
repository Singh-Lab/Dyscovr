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
library(SuppDists)
library(broom)
library(ggsci)

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
    plot(grp1_Betas, grp2_Betas, pch = 19, col = alpha("#0072B599", 0.4), #main = "Betas Spearman Correlation",
         xlab = paste(ri_1, "Betas"), ylab = paste(ri_2, "Betas"))
    abline(lm(grp2_Betas ~ grp1_Betas), col = "#BC3C29FF", lwd = 3)
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
#' @param beta_thres a Beta value that the target gene must exceed to be included, optional
#' @param qval_thres a q-value that the target gene must be below to be included, optional
compute_and_print_spearman_multDF <- function(results_table1, results_table2, ri_1, ri_2,
                                              beta_thres, qval_thres) {
  
  # Optionally limit to genes with a qval that is below some threshold
  if(!is.na(qval_thres)) {
    results_table1 <- results_table1[results_table1$q.value < qval_thres,]
    results_table2 <- results_table2[results_table2$q.value < qval_thres,]
  }
  
  if((ri_1 %fin% results_table1$R_i.name) & (ri_2 %fin% results_table2$R_i.name)) {
    
    target_genes <- intersect(unique(results_table1[results_table1$R_i.name == ri_1, 'T_k.name']), 
                              unique(results_table2[results_table2$R_i.name == ri_2, 'T_k.name']))
    
    grp1_Betas <- get_values(results_table1, target_genes, ri_1, "negLogPval", beta_thres)
    grp2_Betas <- get_values(results_table2, target_genes, ri_2, "negLogPval", beta_thres)
    print(head(grp1_Betas))
    print(head(grp2_Betas))
    
    # Remove extreme outliers 
    #iqr <- IQR(c(grp1_Betas, grp2_Betas))
    # Get Q3 of data
    #q3 <- as.numeric(quantile(c(grp1_Betas, grp2_Betas))[4])
    # Get Q3 + 1.5(IQR)
    #thres <- as.numeric(q3 + (1.5 * iqr))
    #print(thres)
    #ind <- c(which((grp1_Betas > thres) | (grp2_Betas > thres)))
    
    #grp1_Betas <- grp1_Betas[-ind]
    #grp2_Betas <- grp2_Betas[-ind]
    
    #na_ind <- unique(c(which(is.na(grp1_Betas)), which(is.na(grp1_Betas))))
    #if(length(na_ind) > 0) {
    #  grp1_Betas <- grp1_Betas[-na_ind]
    #  grp2_Betas <- grp2_Betas[-na_ind]
    #}

    #print(length(grp2_Betas))
    
    # Get Betas spearman
    betas_spearman <- cor.test(grp1_Betas, grp2_Betas, method = "spearman", exact = F)
    betas_spearman_stat <- as.numeric(betas_spearman$estimate)
    betas_spearman_pval <- betas_spearman$p.value
    #betas_spearman_precise_pval <- pSpearman()
    print(betas_spearman)
    
    # Print the results
    print(paste("Spearman results for", paste(ri_1, paste("and", ri_2))))
    print(paste("Correlation of", paste(betas_spearman_stat, paste(", p-value of", betas_spearman_pval))))
    
    # Create a plot to visualize the correlations
    # "#BC3C29FF", "#0072B5FF"
    plot(grp1_Betas, grp2_Betas, pch = 19, col = alpha("#0072B599", 0.4), #main = "Betas Spearman Correlation",
         xlab = "", ylab = "", bty="n")
    abline(lm(grp2_Betas ~ grp1_Betas), col = "#BC3C29FF", lwd = 3)
    
    # Get R2 and P-value as well
    lm_res <- summary(lm(grp2_Betas ~ grp1_Betas))
    r2 <- lm_res$r.squared
    pval <- as.numeric(lm_res$coefficients[,4])[2]
    print(paste("R2:", paste(r2, paste(", p-value:", pval))))
    
    #text(labels = paste("Correlation:", paste(round(betas_spearman_stat, 4), 
    #                                          paste(", p-value:", format(betas_spearman_pval, digits = 6, scientific = TRUE)))), 
    #     x = max(grp1_Betas, na.rm = TRUE)-sd(grp1_Betas, na.rm = TRUE)*3, 
    #     y = max(grp2_Betas, na.rm = TRUE)-sd(grp2_Betas, na.rm = TRUE), 
    #     col = "black", font=2)
    
    r2_round <- round(r2, 4)
    pval_formatted <- format(pval, digits = 2, scientific = TRUE)
    label <- bquote(R^2 ~ "=" ~ .(r2_round) ~ ", p =" ~ .(pval_formatted))
    text(labels = label, 
         x = max(grp1_Betas, na.rm = TRUE)-sd(grp1_Betas, na.rm = TRUE)*3, 
         y = max(grp2_Betas, na.rm = TRUE), #-sd(grp2_Betas, na.rm = TRUE), 
         col = "black", font=2)
    mtext(side=1, line=3, paste(ri_1, "T-Statistic, TCGA-BRCA"), font=2, cex=1.2)
    mtext(side=2, line=3, paste(ri_2, "T-Statistic, METABRIC"), font=2, cex=1.2)
    
    
    
  } else {print("Error. Provided regprots are not in the given master DF.")}
}

# Mini function to get Betas/t-statistics for each target gene
#' @param results_table the master DF 
#' @param target_genes the list of target genes
#' @param ri the given regulatory protein
#' @param type "Betas" or "t-statistics" to indicate what value we are returning
get_values <- function(results_table, target_genes, ri, type, beta_thres) {
  results_table_ri <- results_table[results_table$R_i.name == ri,]
  vals <- c()
  if(type == "Betas") {
    vals <- as.numeric(unlist(lapply(target_genes, function(tg) 
      results_table_ri[results_table_ri$T_k.name == tg, 'estimate'])))
  }
  else if(type == "spearman") {
    vals <- as.numeric(unlist(lapply(target_genes, function(tg) 
      results_table_ri[results_table_ri$T_k.name == tg, 'statistic'])))
  } else if (type == "negLogPval") {
    vals <- as.numeric(unlist(lapply(target_genes, function(tg) {
      pval <- -log2(results_table_ri[results_table_ri$T_k.name == tg, 'p.value'])
      beta <- results_table_ri[results_table_ri$T_k.name == tg, 'estimate']
      return(ifelse(beta < 0, (pval*-1), pval))
    })))
  } else {
    print("Unimplemented type. Returning NA.")
    vals <- rep(NA, times = length(target_genes))
  }
  
  names(vals) <- target_genes
  if(!is.na(beta_thres)) {
    vals <- unlist(lapply(vals, function(est) ifelse(abs(est < beta_thres), NA, est)))
    vals <- vals[!is.na(vals)]
  }
  return(vals)
}

ri_1 <- "TP53"
ri_2 <- "PIK3CA"

# Call function
compute_and_print_spearman_multDF(master_df1, master_df2, ri_1, ri_2, NA, NA)

#' Helper to apply this function systematically across multiple data types
#' @param list_of_data_tables1 the first list of master DFs
#' @param list_of_data_tables2 the first list of master DFs
#' @param ri_1 the external gene name of the first protein of interest
#' @param ri_2 the external gene name of the second protein of interest
#' @param beta_thres a Beta value that the target gene must exceed to be included, optional
#' @param qval_thres a q-value that the target gene must be below to be included, optional
compute_and_print_spearman_multDF_multCTs <- function(list_of_data_tables1, list_of_data_tables2, 
                                                      ri_1, ri_2, beta_thres, qval_thres) {
  intersecting_names <- intersect(names(list_of_data_tables1), names(list_of_data_tables2))
  print(intersecting_names)
  
  spearman_res_list <- lapply(1:length(intersecting_names), function(i) {
    name <- intersecting_names[i]
    print(name)
    results_table1 <- list_of_data_tables1[[which(names(list_of_data_tables1) == name)]]
    results_table2 <- list_of_data_tables2[[which(names(list_of_data_tables2) == name)]]
    
    if(length(results_table1) > 0) {results_table1 <- as.data.frame(results_table1)}
    else {return(NA)}
    if(length(results_table2) > 0) {results_table2 <- as.data.frame(results_table2)}
    else {return(NA)}

    # Optionally limit to genes with a qval that is below some threshold
    if(!is.na(qval_thres)) {
      results_table1_sub <- results_table1[results_table1$q.value < qval_thres,]
      results_table2_sub <- results_table2[results_table2$q.value < qval_thres,]
    } else {
      results_table1_sub <- results_table1
      results_table2_sub <- results_table2
    }
    
    if((ri_1 %fin% results_table1_sub$R_i.name) & (ri_2 %fin% results_table2_sub$R_i.name)) {
      
      target_genes <- intersect(unique(results_table1_sub[results_table1_sub$R_i.name == ri_1, 'T_k.name']), 
                                unique(results_table2_sub[results_table2_sub$R_i.name == ri_2, 'T_k.name']))

      grp1_Betas <- get_values(results_table1_sub, target_genes, ri_1, "T-statistic", beta_thres)
      grp2_Betas <- get_values(results_table2_sub, target_genes, ri_2, "T-statistic", beta_thres)
      
      na_ind <- unique(c(which(is.na(grp1_Betas)), which(is.na(grp1_Betas))))
      if(length(na_ind) > 0) {
        grp1_Betas <- grp1_Betas[-na_ind]
        grp2_Betas <- grp2_Betas[-na_ind]
      }
      
      #print(length(grp2_Betas))
      
      # Get Betas spearman
      betas_spearman <- NA
      tryCatch({
        betas_spearman <- cor.test(grp1_Betas, grp2_Betas, method = "spearman", exact = F)
      }, error=function(cond){
        print(cond)
      })
      if(!(length(betas_spearman) == 1)) {
        betas_spearman_stat <- as.numeric(betas_spearman$estimate)
        betas_spearman_pval <- betas_spearman$p.value
        #betas_spearman_precise_pval <- pSpearman()
        print(name)
        print(paste("Spearman results for", paste(ri_1, paste("and", ri_2))))
        print(paste("T-Statistic correlation of", paste(betas_spearman_stat, paste(", p-value of", betas_spearman_pval))))
        
        return(data.frame("Cancer.Type" = rep(name, times = length(grp1_Betas)), "Grp1.Betas" = grp1_Betas,
                          "Grp2.Betas" = grp2_Betas, "Pval" = rep(betas_spearman_pval, times = length(grp1_Betas))))
      } else {return(NA)}
    }
  })
  spearman_res_list <- spearman_res_list[!is.na(spearman_res_list)]
  spearman_res_df_full <- do.call(rbind, spearman_res_list)
  
  g <- ggplot(spearman_res_df_full, aes(x=Grp1.Betas, y=Grp2.Betas)) + geom_point() +
    geom_smooth(method=lm) + labs(x = paste(ri_1, "T-statistic, Male"), y = paste(ri_1, "T-statistic, Female")) + 
    facet_wrap(~Cancer.Type, scales ="free")
  print(g)
  
  return(spearman_res_df_full)
}


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
#' @param qval_thres optional q-value threshold for significance
create_spearman_barplot <- function(source1_master_list, source2_master_list, qval_thres) {
  
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
    
    #spearman <- tidy(cor.test(betas_df1, betas_df2, method = "spearman", use = "pairwise"))
    spearman <- tidy(cor.test(neglogpvals_wdir_df1, neglogpvals_wdir_df2, 
                              method = "spearman", use = "pairwise"))
    
    return(spearman)
  })
  spearman_cor_vals <- as.numeric(unlist(lapply(spearman_correlations, function(x) x$estimate)))
  spearman_correlations_df <- data.frame("Gene" = names(source1_master_list), 
                                         "Spearman" = spearman_cor_vals)
  print("P-Values")
  print(unlist(lapply(spearman_correlations, function(x) as.numeric(x$p.value))))
  print("Correlations")
  print(as.numeric(unlist(lapply(spearman_correlations, function(x) x$estimate))))
  
  ylab <- "Spearman Corr, METABRIC & TCGA-BRCA"
  if(!is.na(qval_thres)) {ylab <- paste(ylab, paste0("(q <", paste0(qval_thres, ")")))}
  p <- ggplot(spearman_correlations_df, aes(x = reorder(Gene, -Spearman, mean), 
                                            y = Spearman, fill = Gene)) + 
    geom_bar(stat = "identity", show.legend = FALSE) + scale_fill_nejm() + 
    theme_minimal() + xlab("Gene") + ylab(ylab) +
    theme(axis.text = element_text(face="bold", size = 12), 
          axis.title=element_text(size=14,face="bold"))
  
  return(p)
}

source1_master_list <- list("TP53" = allgenes_p53_incl2GenoPCs, "PIK3CA" = allgenes_pik3ca_incl2GenoPCs,
                            "KMT2C" = allgenes_kmt2c, "GATA3" = allgenes_gata3)
source2_master_list <- list("TP53" = allgenes_p53_metabric, "PIK3CA" = allgenes_pik3ca_metabric,
                            "KMT2C" = allgenes_kmt2c_metabric, "GATA3" = allgenes_gata3_metabric)

create_spearman_barplot(source1_master_list, source2_master_list, NA) 

# Spearman correlation table
#' @param source1_master_list a list of master DFs from a certain source 1
#' @param source2_master_list a list of master DFs (for the same genes) from a certain
#' source 2. Two source lists should be the same length.
#' @param qval_thres optional q-value threshold for significance
generate_spearman_correlation_table <- function(source1_master_list, source2_master_list, qval_thres) {
  spearman_correlation_tables <- lapply(1:length(source1_master_list), function(i) {
    df1 <- source1_master_list[[i]]
    df2 <- source2_master_list[[i]]
    
    #df1_sig_up <- df1[(df1$q.value < qval_thres) & (df1$estimate > 0), c('estimate', 'p.value', 'T_k.name')]
    #df1_sig_down <- df1[(df1$q.value < qval_thres) & (df1$estimate < 0), c('estimate', 'p.value', 'T_k.name')]
    df1_sig <- df1[(df1$q.value < qval_thres), c('estimate', 'p.value', 'T_k.name')]
    df1_nonsig <- df1[(df1$q.value >= qval_thres), c('estimate', 'p.value', 'T_k.name')]
    #df1_list <- list(df1_sig_up, df1_sig_down, df1_nonsig)
    df1_list <- list(df1_sig, df1_nonsig)
    
    #df2_sig_up <- df2[(df2$q.value < qval_thres) & (df2$estimate > 0), c('estimate', 'p.value', 'T_k.name')]
    #df2_sig_down <- df2[(df2$q.value < qval_thres) & (df2$estimate < 0), c('estimate', 'p.value', 'T_k.name')]
    df2_sig <- df2[(df2$q.value < qval_thres), c('estimate', 'p.value', 'T_k.name')]
    df2_nonsig <- df2[(df2$q.value >= qval_thres), c('estimate', 'p.value', 'T_k.name')]
    #df2_list <- list(df2_sig_up, df2_sig_down, df2_nonsig)
    df2_list <- list(df2_sig, df2_nonsig)
    
    #results_table <- data.frame("a" = c(1,4,7), "b" = c(2,5,8), "c" = c(3,6,9))
    results_table <- data.frame("a" = c(1,3), "b" = c(2,4))

    get_spearman <- function(x) {
      print(x)
      df1 <- NA
      df2 <- NA
      
      if(x < 3) {
        df1 <- df1_list[[1]] 
        df2 <- df2_list[[x]] 
      } else {
        df1 <- df1_list[[2]] 
        df2 <- df2_list[[x-2]] 
      }
      
      #if(x < 4) {
      #  df1 <- df1_list[[1]] 
      #  df2 <- df2_list[[x]] 
      #} else if((x > 3) & (x < 7)) {
      #  df1 <- df1_list[[2]] 
      #  df2 <- df2_list[[x-3]] 
      #} else {
      #  df1 <- df1_list[[3]] 
      #  df2 <- df2_list[[x-6]] 
      #}
      
      df_merged <- merge(df1, df2, by = "T_k.name", all = FALSE)
      print(head(df_merged))

      if(nrow(df_merged) == 0) {return(0)}
      
      neglogpvals_wdir_df1 <- unlist(lapply(1:nrow(df_merged), function(i) {
        beta <- df_merged$estimate.x[i]
        pval <- -log2(df_merged$p.value.x[i])
        return(ifelse(beta < 0, (-1)*pval, pval))
      }))
      neglogpvals_wdir_df2 <- unlist(lapply(1:nrow(df_merged), function(i) {
        beta <- df_merged$estimate.y[i]
        pval <- -log2(df_merged$p.value.y[i])
        return(ifelse(beta < 0, (-1)*pval, pval))
      }))
      #print(head(neglogpvals_wdir_df1))
      #print(head(neglogpvals_wdir_df2))
      spearman <- tidy(cor.test(neglogpvals_wdir_df1, neglogpvals_wdir_df2, 
                                method = "spearman", use = "pairwise"))
      print(spearman)
      print(spearman$p.value)
      return(spearman$p.value)
    }
    
    results_table <- apply(results_table, MARGIN = c(2,1), FUN = get_spearman)
    #rownames(results_table) <- c("Sig.Up", "Sig.Down", "Neither")
    #colnames(results_table) <- c("Sig.Up", "Sig.Down", "Neither")
    rownames(results_table) <- c("Sig", "Not Sig")
    colnames(results_table) <- c("Sig", "Not Sig")
    
    print(names(source1_master_list)[i])
    print(results_table)
    return(results_table)
  })
  return(spearman_correlation_tables)
}

pval_spearman_tables <- generate_spearman_correlation_table(list("TP53" = pc_allGenes[pc_allGenes$R_i.name == "TP53",], "PIK3CA" = pc_allGenes[pc_allGenes$R_i.name == "PIK3CA",]), 
                                    list("TP53" = metabric_allGenes[metabric_allGenes$R_i.name == "TP53",], "PIK3CA" = metabric_allGenes[metabric_allGenes$R_i.name == "PIK3CA",]), 
                                   0.01) 
names(pval_spearman_tables) <- c("TP53", "PIK3CA")


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

compute_pcc(top_drivers_0.05[["BRCA"]][top_drivers_0.05[["BRCA"]]$R_i.name == "TP53",], metabric_allGene[metabric_allGene$R_i.name == "TP53",], 0.1)
