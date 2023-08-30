############################################################
### Process Linear Model Output
### Written By: Sara Geraghty, May 2021
############################################################

# Given an output 'master file' from a linear model run, this
# file performs the following functions:
# 1. Basic Output Visualization (Beta Distribution, P-Value Distribution,
# Q-Q Plot, and SE Distribution)
# 2. Perform multiple hypothesis testing correction
# 4. Obtain significant correlations (q-value thresholding)


# BiocManager::install("pathview")
# BiocManager::install("RCy3")
library(TCGAbiolinks)
library(broom)
library(qvalue)
library(ggplot2)
library(stringr)
library(dplyr)
library(VennDiagram)
library("RColorBrewer")
library("pathview")
library("RCy3")
library(igraph)
library(scales)
library(reshape2)
library(ggsci)
library(parallel)
library(forcats)

# NEJM color palatte: https://nanx.me/ggsci/reference/pal_nejm.html


# Path to output files
main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
#main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"

# Path to where output figures should be saved
output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/BRCA/Linear Model/"
#output_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Output Visualizations/Pan-Cancer/Linear Model/"


# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)


## Read in the master DF output files from multiple cancer types
top_drivers_0.05_fns <- intersect(list.files(path = "C:/Users/sarae/top_0.05/", recursive = T,
                                   pattern = "_corrected_MUT"), intersect(
                                  list.files(path = "C:/Users/sarae/top_0.05/", recursive = T,
                                             pattern = "metabolicTargs"),
                                  list.files(path = "C:/Users/sarae/top_0.05/", recursive = T,
                                             pattern = "Nonsyn.Drivers.Vogel.Top5")))
top_drivers_0.05_fns <- top_drivers_0.05_fns[!grepl("Archives", top_drivers_0.05_fns)]
top_drivers_0.1_fns <- intersect(list.files(path = "C:/Users/sarae/top_0.1/", recursive = T,
                                   pattern = "_corrected_MUT"),intersect(
                                     list.files(path = "C:/Users/sarae/top_0.1/", recursive = T,
                                                pattern = "metabolicTargs"),
                                     list.files(path = "C:/Users/sarae/top_0.1/", recursive = T,
                                                pattern = "Nonsyn.Drivers.Vogel")))
top_drivers_0.15_fns <- intersect(list.files(path = "C:/Users/sarae/top_0.15/", recursive = T,
                                            pattern = "_corrected_MUT"),intersect(
                                              list.files(path = "C:/Users/sarae/top_0.15/", recursive = T,
                                                         pattern = "metabolicTargs"),
                                              list.files(path = "C:/Users/sarae/top_0.15/", recursive = T,
                                                         pattern = "Nonsyn.Drivers.Vogel")))
top_drivers_0.05 <- lapply(top_drivers_0.05_fns, function(f) 
  fread(paste0("C:/Users/sarae/top_0.05/", f), header = T))
names(top_drivers_0.05) <- unlist(lapply(top_drivers_0.05_fns, function(x)
  unlist(strsplit(x, "/", fixed = T))[1]))
top_drivers_0.1 <- lapply(top_drivers_0.1_fns, function(f) 
  fread(paste0("C:/Users/sarae/top_0.1/", f), header = T))
names(top_drivers_0.1) <- unlist(lapply(top_drivers_0.1_fns, function(x)
  unlist(strsplit(x, "/", fixed = T))[1]))
top_drivers_0.15 <- lapply(top_drivers_0.15_fns, function(f) 
  fread(paste0("C:/Users/sarae/top_0.15/", f), header = T))
names(top_drivers_0.15) <- unlist(lapply(top_drivers_0.15_fns, function(x)
  unlist(strsplit(x, "/", fixed = T))[1]))

# If each cancer type has multiple files, combine them
top_drivers_0.05 <- lapply(unique(names(top_drivers_0.05)), function(ct)
  do.call(rbind, top_drivers_0.05[names(top_drivers_0.05) == ct]))
names(top_drivers_0.05) <- unique(unlist(lapply(top_drivers_0.05_fns, function(x)
  unlist(strsplit(x, "/", fixed = T))[1])))
# Redo the MHT correction
top_drivers_0.05 <- lapply(top_drivers_0.05, function(x) mh_correct(x, T))
top_drivers_0.05 <- lapply(top_drivers_0.05, function(x) x[order(x$q.value),])

# Write these back
lapply(1:length(top_drivers_0.05), function(i) {
  df <- top_drivers_0.05[[i]]
  ct <- names(top_drivers_0.05)[i]
  write.csv(df, paste0("C:/Users/sarae/top_0.05/", 
                       paste(ct, "res_top_0.05_allGenes_quantile_rawCNA_methMRaw_3PCs_Nonsyn.Drivers.Vogel.elim.vif.5_corrected_MUT.csv", sep = "/")))
})

# If we've already done this 
top_drivers_0.05_fns <- top_drivers_0.05_fns[!grepl("pt", top_drivers_0.05_fns)]
top_drivers_0.05 <- lapply(top_drivers_0.05_fns, function(f) 
  fread(paste0("C:/Users/sarae/top_0.05/", f), header = T))
names(top_drivers_0.05) <- unlist(lapply(top_drivers_0.05_fns, function(x)
  unlist(strsplit(x, "/", fixed = T))[1]))


top_drivers_0.05_fns_metabol_top5only <- top_drivers_0.05_fns_metabol_top5only[(grepl("BLCA", top_drivers_0.05_fns_metabol_top5only)) | (grepl("HNSC", top_drivers_0.05_fns_metabol_top5only)) | grepl("COAD", top_drivers_0.05_fns_metabol_top5only) | (grepl("ESCA", top_drivers_0.05_fns_metabol_top5only)) | grepl("LGG", top_drivers_0.05_fns_metabol_top5only) | grepl("LUSC", top_drivers_0.05_fns_metabol_top5only) | grepl("UCEC", top_drivers_0.05_fns_metabol_top5only)]
top_drivers_0.05_fns_metabol_top5only <- top_drivers_0.05_fns_metabol_top5only[!grepl("Full_Vogel", top_drivers_0.05_fns_metabol_top5only)]
top_drivers_0.05_fns_metabol_top5only <- top_drivers_0.05_fns_metabol_top5only[!grepl("male", top_drivers_0.05_fns_metabol_top5only)]
top_drivers_0.05_metabol_top5only <- lapply(top_drivers_0.05_fns_metabol_top5only, function(f) 
  fread(paste0("C:/Users/sarae/top_0.05/", f), header = T))
names(top_drivers_0.05_metabol_top5only) <- unlist(lapply(top_drivers_0.05_fns_metabol_top5only, function(x)
  unlist(strsplit(x, "/", fixed = T))[1]))

# Make simple version for user viewing
top_drivers_0.05_simple <- lapply(top_drivers_0.05, function(x) 
  x[, c("estimate", "std.error", "statistic", "p.value", "q.value", "T_k.name", "R_i.name")])
top_drivers_0.05_simple <- lapply(top_drivers_0.05_simple, function(x) 
  do.call(cbind, list(x[,"R_i.name"], x[,"T_k.name"], x[, c("estimate", "std.error", "statistic", "p.value", "q.value")])))
top_drivers_0.05_simple <- lapply(top_drivers_0.05_simple, function(x) {
  colnames(x) <- c("Driver_Name", "Target_Name", "Mutation_Coefficient", "Standard_Error", 
                   "TStatistic", "pvalue", "qvalue")
  return(x)
})
lapply(1:length(top_drivers_0.05_simple), function(i) {
  fn <- paste0("TCGA_", paste0(names(top_drivers_0.05_simple)[i], "_Output.csv"))
  write.csv(top_drivers_0.05_simple[[i]], paste0("C:/Users/sarae/top_0.05/", fn))
})

############################################################
############################################################
#### BASIC VISUALIZTION OF OUTPUT
############################################################
############################################################

############################################################
#### VISUALIZE BETA VALUE DISTRIBUTION
############################################################
#' Function plots a histogram to visualize the distribution
#' of regulatory protein r_i beta values across all linear models.
#' @param results_table a master DF produced from run_linear_model()
visualize_beta_distrib <- function(results_table) {
  betas <- results_table$estimate
  hist(betas, main = "Histogram of Beta Coefficient Values",
       xlab = "Beta Coefficient Value", ylab = "Frequency")
}

# Call this function & save output
fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/Mutation, eQTL/Beta Distribution (I-Protein, TMM, bucketCNA, iciTotFrac).png", sep = "")
png(fn, width = 450, height = 350)
visualize_beta_distrib(master_df_mut)
dev.off()

fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/CNA, eQTL/Beta Distribution (I-Protein, TMM, bucketCNA, iciTotFrac).png", sep = "")
png(fn, width = 450, height = 350)
visualize_beta_distrib(master_df_cna)
dev.off()


############################################################
#### VISUALIZE P-VALUE DISTRIBUTION
############################################################
#' Function plots a histogram to visualize the distribution
#' of regulatory protein r_i beta values across all linear models.
#' @param results_table a master DF produced from run_linear_model()
visualize_pval_distrib <- function(results_table) {
  #pvals <- results_table$p.value[!is.na(results_table$p.value) & 
                                   #!is.infinite(results_table$p.value)]
  #hist(pvals, main = "",
       #xlab = "p-value", ylab = "Frequency", col = "blueviolet")
  
  ggplot(results_table, aes(x = p.value)) + 
    geom_histogram(alpha = 1, position = "identity", bins = 50, fill = "#0072B5FF") + 
    theme_minimal() + xlab("P-Value") + ylab("Frequency") + 
    labs(fill = "") + theme(axis.title = element_text(face = "bold", size = 14),
                            axis.text = element_text(face = "bold", size = 12),
                            legend.text = element_text(size=12),
                            legend.title = element_text(face = "bold", size = 14)) +
    geom_vline(xintercept = 0.05, linetype='dashed', size = 1)+
    annotate("text", y = 100, x = 0.01, label = "p < 0.05", hjust = -0.5, size = 5, fontface = "bold")
  
}

# Call this function & save output
fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/Mutation, eQTL/P-Value Distribution (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 450, height = 350)
visualize_pval_distrib(master_df_mut)
dev.off()

fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/CNA, eQTL/P-Value Distribution (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 450, height = 350)
visualize_pval_distrib(master_df_cna)
dev.off()

############################################################
#### VISUALIZE Q-Q PLOT
############################################################
#' Function plots a Q-Q plot to visualize the distribution of p-values
#' and assess whether they come from a uniform distribution
#' @param results_table a master DF produced from run_linear_model()
qqplot_pvals <- function(results_table) {
  qqnorm(results_table$p.value, pch = 1, frame = FALSE)
  qqline(results_table$p.value, col = "#0072B5FF", lwd = 2)
}

#' Function plots a Q-Q plot to visualize the distribution of standard error valeus
#' and assess whether they come from a uniform distribution
#' @param results_table a master DF produced from run_linear_model()
qqplot_stderror <- function(results_table) {
  qqnorm(results_table$std.error, pch = 1, frame = FALSE)
  qqline(results_table$std.error, col = "#0072B5FF", lwd = 2)
}


# Call this function
#qqplot_pvals(master_df)


############################################################
#### VISUALIZE ERROR DISTRIBUTION
############################################################
#' Function plots a histogram of the standard errors produced from
#' all LM runs to assess whether the errors derive from a normal distribution
#' @param results_table a master DF produced from run_linear_model()
visualize_error_distrib <- function(results_table) {
  hist(results_table$std.error, main = "Standard Error Distribution Across All Tests",
       xlab = "Standard Error (SE)", ylab = "Frequency")
}

# Call this function
#visualize_error_distrib(master_df)


############################################################
############################################################
#### PEFORM MULTIPLE HYPOTHESIS TESTING CORRECTION
############################################################
############################################################
#' Function takes in an output results table and applies multiple
#' hypothesis testing correction (Storey's q-value correction) to 
#' all p-values in order to add a column of q-values. Returns 
#' the results table with a column for q-values.
#' @param results_table a master DF produced from run_linear_model()
#' @param per_gene a T/F value indicating whether or not we are doing the correction
#' per-R_i (T) or together across all R_i (F)
mh_correct <- function(results_table, per_gene) {
  
  # Get the qvalue object
  qobj <- NA
  if(length(results_table$p.value) < 100) {
    # With a small number of pvalues we may not be able to accurately estimate pi0,
    # so we set to 1 (the equivalent of B-H correction)
    qobj <- qvalue(p = results_table$p.value, pi0 = 1)
    
    # OPT: plot some useful plots & print some useful information
    plot(qobj)
    print(summary(qobj))
    #print(paste("Pi0 (Propr. of true null hypotheses):", qobj$pi0))
    
    qvals <- qobj$qvalues # extract qvalues
    
    results_table$q.value <- qvals # add the qvalues back to the data frame
    
  } else {
    if(per_gene) {
      unique_drivers <- unique(results_table$R_i)
      list_of_corrected_tabs <- lapply(unique_drivers, function(d) {
        res_tab_sub <- results_table[results_table$R_i == d, ]
        qobj <- qvalue(p = res_tab_sub$p.value)
        qvals <- qobj$qvalues # extract qvalues
        res_tab_sub$q.value <- qvals # add the qvalues back to the data frame
        return(res_tab_sub)
      })
      results_table <- do.call(rbind, list_of_corrected_tabs)
      results_table <- results_table[order(results_table$q.value),]
    }
    else {
      qobj <- qvalue(p = results_table$p.value)
      # OPT: plot some useful plots & print some useful information
      plot(qobj)
      print(summary(qobj))
      #print(paste("Pi0 (Propr. of true null hypotheses):", qobj$pi0))
      
      qvals <- qobj$qvalues # extract qvalues
      
      results_table$q.value <- qvals # add the qvalues back to the data frame
    }
  }
  # Return the tidied linear model fit with q-values
  return(results_table)
}

# Call this function
master_df_mut_corrected <- mh_correct(master_df_mut, T)
master_df_cna_corrected <- mh_correct(master_df_cna, T)

#Q-Value Visualization (I-Protein, FPKM, rawCNA, methBeta, iciTotFrac)

top_drivers_0.05 <- lapply(top_drivers_0.05, function(x) if(nrow(x) > 0) mh_correct(x, T))
top_drivers_0.1 <- lapply(top_drivers_0.1, function(x) if(nrow(x) > 0) mh_correct(x, T))

############################################################
############################################################
#### ADD GENE NAMES TO FILE
############################################################
############################################################
#' Given a master data frame result from model, add a column for target gene name
#' and regulatory protein gene name. Return the updated data frame.
#' @param master_df_sig a data.table object produced from the linear_model.R function
#' @param all_genes_id_conv a bioMart file with conversions between different gene ID types
add_targ_regprot_gns <- function(master_df_sig, all_genes_id_conv) {
  master_df_sig$T_k.name <- unlist(mclapply(master_df_sig$T_k, function(x) 
    paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == unlist(strsplit(x, ";", fixed = TRUE))[1], 
                                   'external_gene_name'])), collapse = ";")))
  
  # Add a column for the regulatory protein name
  unique_ri <- data.frame("R_i" = unique(master_df_sig$R_i))
  unique_ri$R_i.name <- unlist(lapply(unique_ri$R_i, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == unlist(strsplit(x, ";", fixed = TRUE))[1], 
                                   'external_gene_name']), collapse = ";")))
  master_df_sig$R_i.name <- unlist(lapply(master_df_sig$R_i, function(x) 
    unique_ri[unique_ri$R_i == x, 'R_i.name']))
  
    return(master_df_sig)
}

# Call this function
master_df_mut_corrected <- add_targ_regprot_gns(master_df_mut_corrected, all_genes_id_conv)
master_df_cna_corrected <- add_targ_regprot_gns(master_df_cna_corrected, all_genes_id_conv)

# Write this to a new file
fwrite(master_df_mut_corrected, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/eQTL/output_results_P53_chipeat_iprotein_TMM_rawCNA_cibersortTotFrac_corrected_MUT.csv", sep = ""))
fwrite(master_df_cna_corrected, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/eQTL/output_results_P53_chipeat_iprotein_TMM_rawCNA_cibersortTotFrac_corrected_CNA.csv", sep = ""))

top_drivers_0.05 <- lapply(top_drivers_0.05, function(x) if(nrow(x) > 0) add_targ_regprot_gns(x, all_genes_id_conv))
top_drivers_0.1 <- lapply(top_drivers_0.1, function(x) if(nrow(x) > 0) add_targ_regprot_gns(x, all_genes_id_conv))

############################################################
############################################################
#### OBTAIN SIGNIFICANT CORRELATIONS
############################################################
############################################################
#' Function takes in an output results table with q-values and 
#' restricts it to only models that exceed the q-value threshold
#' (are statistically significant correlations). Then ranks the 
#' remaining by q-values and returns the ranked top hits list.
#' @param results_table a master DF produced from run_linear_model() that has q-values added
#' from the mh_correct() function
#' @param qval_thres a threshold for significance for q-values
get_signif_correl <- function(results_table, qval_thres) {
  
  # Limit to only entries that exceed the given qvalue threshold
  results_table_sig <- results_table %>% filter(q.value < qval_thres)
  
  # Sort the table by qvalue
  results_table_sig_ordered <- results_table_sig[order(results_table_sig$q.value, 
                                                       decreasing = FALSE),]
  
  return(results_table_sig_ordered)
}

# Call this function
qval_thres <- 0.1
master_df_mut_sig <- get_signif_correl(master_df_mut_corrected, qval_thres)
master_df_cna_sig <- get_signif_correl(master_df_cna_corrected, qval_thres)

# Write these results to a new file
fwrite(master_df_mut_sig, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/eQTL/significant_output_P53_chipeat_iprotein_tmm_rawCNA_methBetaRaw_cibersortTotFrac_MUT.csv", sep = ""))                                                                                                                     # 118 (including CNAi and Methi as covar.)
fwrite(master_df_cna_sig, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/eQTL/significant_output_P53_chipeat_iprotein_tmm_rawCNA_methBetaRaw_cibersortTotFrac_CNA.csv", sep = ""))                                                                                                                     # 118 (including CNAi and Methi as covar.)


############################################################
############################################################
#### BAR CHARTS TO SHOW NUMBERS OF SIGNIFICANT HITS
#### FOR EACH OF QUERY GENES, ACROSS MULTIPLE DFs
############################################################
############################################################
#' Create barplot with the number of missense/nonsense and silent significant hits 
#' for drivers and non-drivers
#' @param list_of_master_dfs a list of master DFs, with the names of the list corresponding
#' to the names of the regulatory proteins of interest; alternatively, could be one combined master DF
#' @param q_thres a q-value threshold for significance
#' @param list_of_silent_dfs OPT: a corresponding list of silent DFs to include in the chart
create_num_sig_hits_barplot <- function(list_of_master_dfs, q_thres, list_of_silent_dfs) {
  # Get the number of significant hits for each master DF
  df <- list_of_master_dfs[[1]]
  unique_drivers <- unique(df$R_i.name)
  sig_hits_df <- data.frame("Gene" = unique_drivers, 
                            "Num.Signif.Hits.q" = rep(0, times = length(unique_drivers)))
  colnames(sig_hits_df)[2] <- paste0(colnames(sig_hits_df)[2], q_thres)
  sig_hits_df[,2] <- unlist(lapply(unique_drivers, function(d) {
      return(nrow(df[(df$q.value < q_thres) & (df$R_i.name == d),]))
  }))
  
  print(head(sig_hits_df))
  
  if(!is.na(list_of_silent_dfs)) {
    sig_hits_df_silent <- data.frame("Gene" = names(list_of_silent_dfs), 
                                     "Num.Signif.Hits.q" = rep(0, times = length(list_of_silent_dfs)))
    colnames(sig_hits_df_silent)[2] <- paste0(colnames(sig_hits_df_silent)[2], q_thres)
    sig_hits_df_silent[,2] <- unlist(lapply(list_of_silent_dfs, function(df) {
      return(nrow(df[df$q.value < q_thres,]))
    }))
    sig_hits_df$Mutation.Type <- "Nonsynonymous"
    sig_hits_df_silent$Mutation.Type <- "Synonymous"
    
    sig_hits_df <- rbind(sig_hits_df, sig_hits_df_silent)
    
    p <- ggplot(sig_hits_df, aes(fill = Mutation.Type, x = reorder(Gene, -Num.Signif.Hits.q0.1, mean), 
                                 y = Num.Signif.Hits.q0.1)) + 
      scale_fill_manual(values = c("#FFDC91FF", "#20854EFF")) +
      geom_bar(position = "dodge", stat = "identity", color = "black") + theme_minimal() + #scale_color_nejm() +
      xlab("\nGene") + ylab(paste0("\n", paste0("Num. Signif. Hits, q < ", q_thres))) +
      theme(axis.text = element_text(face="bold", size = 16), 
            axis.title=element_text(size=18,face="bold"))
    
  } else {
    p <- ggplot(sig_hits_df, aes(y = Num.Signif.Hits.q0.1, fill = Gene, 
                                 x = reorder(Gene, -Num.Signif.Hits.q0.1, mean))) + 
      geom_bar(position = "dodge", width = 0.95, stat = "identity", show.legend = FALSE, color = "black") + 
      scale_fill_manual(values = c("#BC3C29FF", "#0072B5FF", "gray", "#20854EFF", "#FFDC91FF", 
                                              "#20854EFF", "#FFDC91FF", "#BC3C29FF", "#0072B5FF", 
                                              "#20854EFF", "#FFDC91FF", "#BC3C29FF", "#0072B5FF", "gray", "gray")) + #scale_color_nejm() +
      xlab("Gene") + ylab(paste0("\n", paste0("Number of hits (q < ", paste0(q_thres, ")\n")))) +
      #theme_minimal() +
      theme(axis.text = element_text(face="bold", size = 16), 
            axis.title=element_text(size=18, face="bold"), panel.grid.major = element_blank(),
            panel.background = element_rect(fill = 'white'))
  }
  
  return(p)
}


list_of_master_dfs <- list("TP53" = allgenes_p53, "PIK3CA" = allgenes_pik3ca,
                           "TTN" = allgenes_ttn, "FOXA1" = allgenes_foxa1, 
                           "SF3B1" = allgenes_sf3b1, "USH2A" = allgenes_ush2a)
# "GATA3" = allgenes_gata3, "KMT2C" = allgenes_kmt2c

create_num_sig_hits_barplot(list_of_master_dfs, 0.2, NA)
create_num_sig_hits_barplot(list_of_master_dfs, 0.1, NA)


#' Simpler version for one DF
qval_thres <- 0.01
sig_hits_df <- master_df[master_df$q.value < qval_thres,]
sig_hits_df_freq <- melt(table(sig_hits_df$R_i.name))
colnames(sig_hits_df_freq) <- c("Driver", "Num.Hits")
sig_hits_df_freq <- sig_hits_df_freq[order(sig_hits_df_freq$Num.Hits),]
ggplot(sig_hits_df_freq, aes(y = Num.Hits, fill = Driver, 
                             x = reorder(Driver, -Num.Hits, mean))) + 
  geom_bar(position = "dodge", width = 0.95, stat = "identity", show.legend = FALSE, color = "black") + 
  scale_fill_manual(values = c("#FFDC91FF", "#20854EFF", "#BC3C29FF", "#0072B5FF", "gray")) + #scale_color_nejm() +
  xlab("Driver") + ylab(paste0("\n", paste0("Number of hits (q < ", paste0(qval_thres, ")")))) +
  theme(axis.text = element_text(face="bold", size = 16), 
        axis.title=element_text(size=18, face="bold"), panel.grid.major = element_blank(),
        panel.background = element_rect(fill = 'white'))


#' Creates a stacked bar plot to show the number of significant hits
#' (at given q-value threshold) found for each query gene, for 
#' each particular subtype of interest
#' @param qval_dict_list a list of qvalue dictionaries (as output above),
#' one for each subtype of interest. Each entry should have the same 
#' query genes.
#' @param thres a qvalue threshold for significance
create_stacked_sig_hits_bar_plot <- function(qval_dict_list, thres) {
  sub_dfs <- lapply(1:length(qval_dict_list), function(i) {
    subtype_entry <- qval_dict_list[[i]]
    prot_dfs <- lapply(1:length(subtype_entry), function(j) {
      prot_entry <- subtype_entry[[j]]
      # Turn this list into a DF
      prot_df <- data.frame("Query.Prot" = names(subtype_entry)[j],
                            "Num.Sig.Hits" = length(prot_entry$Names[which(prot_entry$Qvals < thres)]))
      #prot_df$Query.Prot <- rep(names(prot_entry)[j], times = ncol(prot_df))
      return(prot_df)
    })
    combined_prot_df <- do.call(rbind, prot_dfs)
    print(head(combined_prot_df))
    # Add the name of the subtype
    print(names(qval_dict_list)[i])
    combined_prot_df$Subtype <- rep(names(qval_dict_list)[i], 
                                    times = nrow(combined_prot_df))
    return(combined_prot_df)
  }) 
  # Combined these into one 
  full_df <- do.call(rbind, sub_dfs)
  #print(head(full_df))
  
  # Create a stacked bar plot from this
  ggplot(full_df, aes(fill=Query.Prot, y=Num.Sig.Hits, x=Subtype)) + 
    geom_bar(position="stack", stat="identity")
  
  return(full_df)
}

#' Perform q-value correction on each R_i individually 
#' @param master_df
get_qval_dict <- function(master_df) {
  
  master_df <- master_df[!(is.na(master_df$p.value)),]
  
  qval_dict <- lapply(1:length(unique(master_df$R_i.name)), function(i) {
    r_i <- unique(master_df$R_i.name)[i]
    #print(r_i)
    pvals <- master_df[master_df$R_i.name == r_i, 'p.value']
    #print(pvals)
    qobj <- NA
    if(length(pvals) < 100) {
      # With a small number of pvalues we may not be able to accurately estimate pi0,
      # so we set to 1 (the equivalent of B-H correction)
      qobj <- qvalue(p = pvals, pi0 = 1)
    } else {
      qobj <- qvalue(p = pvals)
    }
    
    # OPT: plot some useful plots & print some useful information
    #plot(qobj)
    #print(summary(qobj))
    #print(paste("Pi0 (Propr. of true null hypotheses):", qobj$pi0))
    
    qvals <- qobj$qvalues
    return(list("Names" = master_df[master_df$R_i.name == r_i, 'T_k.name'],
                "Qvals" = qvals))
  })
  names(qval_dict) <- unique(master_df$R_i.name)
  return(qval_dict)
}


#' Use this qval to print the number of significant hits for each R_i
#' @param qval_dict a q-value dictionary as output above (Names and Qvals
#' for each query protein)
print_num_sig_hits <- function(qval_dict, thres) {
  for(i in 1:length(qval_dict)) {
    print(paste("Protein name:", names(qval_dict)[i]))
    print(paste("Number of significant hits, q <", thres))
    curr_entry <- qval_dict[[i]]
    print(length(curr_entry$Names[which(curr_entry$Qvals < thres)]))
  }
}


# Create a list of qval dictionaries for each subtype
qval_dict_lumA <- get_qval_dict(master_df_lumA)
qval_dict_lumB <- get_qval_dict(master_df_lumB)
qval_dict_basal <- get_qval_dict(master_df_basal)
qval_dict_her2 <- get_qval_dict(master_df_her2)

qval_dict_list <- list("LumA" = qval_dict_lumA, "LumB" = qval_dict_lumB,
                       "Basal" = qval_dict_basal, "HER2" = qval_dict_her2)

combined_df <- create_stacked_sig_hits_bar_plot(qval_dict_list, 0.1)

driver_genes <- c("TP53", "PIK3CA", "KMT2C")
combined_df$Status <- unlist(lapply(combined_df$Query.Prot, function(g) 
  ifelse(g %in% driver_genes, "Driver", "Non-Driver")))

nb = length(unique(combined_df$Query.Prot))
nm = length(unique(combined_df$Status))
colors = apply(expand.grid(seq(70,40,length=nm), 100, seq(15,375,length=nb+1)[1:nb]), 1, 
               function(x) hcl(x[3],x[2],x[1]))
colors[3] <- "#E8909C"

ggplot(combined_df, aes(fill=interaction(Query.Prot, Status), y=Num.Sig.Hits, x=Subtype)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=colors) + theme_classic()


# Examine significant hit overlap
p53_top_allgenes <- list("LumA" = qval_dict_lumA_allgenes$TP53$Names[which(qval_dict_lumA_allgenes$TP53$Qvals < 0.2)],
                         "LumB" = qval_dict_lumB_allgenes$TP53$Names[which(qval_dict_lumB_allgenes$TP53$Qvals < 0.2)],
                         "Basal" = qval_dict_basal_allgenes$TP53$Names[which(qval_dict_basal_allgenes$TP53$Qvals < 0.2)],
                         "HER2" = qval_dict_her2_allgenes$TP53$Names[which(qval_dict_her2_allgenes$TP53$Qvals < 0.2)])
plt <- venn.diagram(p53_top_allgenes, category.names = c("LumA", "LumB", "Basal", "HER2"), 
                    filename = NULL, output = TRUE, lwd = 2, lty = 'blank', 
                    fill = c("red", "blue", "green", "yellow"), cex = 2, fontface = "bold",
                    fontfamily = "sans", cat.cex = 2, cat.fontface = "bold",cat.fontfamily = "sans")
                    #cat.default.pos = "outer", cat.pos = c(180, 180, 180, 180)) #, cat.fontfamily = "sans", rotation = 1)
grid::grid.draw(plt)



############################################################
############################################################
#### PLOT OVERLAP BETWEEN MUTATION & CNA TOP HITS
############################################################
############################################################
#' Plots the overlap in significant top target gene hits for mutation 
#' and CNA results. 
#' @param master_df_mut_sig mutation master DF, with gene names added and 
#' thresholded for significance
#'@param master_df_cna_sig CNA master DF, with gene names added and 
#' thresholded for significance
plot_tophit_overlap <- function(master_df_mut_sig, master_df_cna_sig) {
  #overlap_genes <- intersect(master_df_mut_sig$T_k.name, master_df_cna_sig$T_k.name)
  #print(overlap_genes)
  
  # Alternatively, look at only combinations of overlap
  pairs_mut_sig <- unlist(lapply(1:nrow(master_df_mut_sig), function(i)
    paste(master_df_mut_sig[i, 'R_i.name'], master_df_mut_sig[i, 'T_k.name'], sep = ":")))
  pairs_cna_sig <- unlist(lapply(1:nrow(master_df_cna_sig), function(i)
    paste(master_df_cna_sig[i, 'R_i.name'], master_df_cna_sig[i, 'T_k.name'], sep = ":")))
  print(paste("Length intersect:", length(intersect(pairs_mut_sig, pairs_cna_sig))))
  
  pairs_list <- list("Mutation" = pairs_mut_sig, "CNA" = pairs_cna_sig)
  
  # Plot a Venn Diagram
  #myCol <- brewer.pal(2, "Pastel2")
  plt <- venn.diagram(pairs_list, category.names = c("Mutation", "CNA"), filename = NULL, output = TRUE,
                      lwd = 2, lty = 'blank', fill = c("red", "blue"), cex = 2, fontface = "bold",
                      fontfamily = "sans", cat.cex = 2, cat.fontface = "bold",cat.fontfamily = "sans",
                      cat.default.pos = "outer", cat.pos = c(180, 180)) #, cat.fontfamily = "sans", rotation = 1)
  grid::grid.draw(plt)
}

# Call this function and write to PNG file
fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/Top Gene Hit Overlap (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 450, height = 350)
plot_tophit_overlap(master_df_mut_sig, master_df_cna_sig)
dev.off()


#' Looks at only top hits for each given master DF (only one regprot), below 
#' a given q-value threshold, and makes a Venn diagram to show overlap, along with 
#' reporting whether there is more overlap than might be expected given
#' the size of the total target pool
#' @param master_df1 the first master DF
#' @param master_df2 the second master DF
#' @param qval_thres the threshold below which we will consider the pairing significant
#' @param size_targgene_pool the number of target genes we are looking at
#' @param master_df1_label a label for what the first master DF is
#' @param master_df2_label a label for what the second master DF is
get_statistical_overlap_of_top_hits <- function(master_df1, master_df2, qval_thres, size_targgene_pool,
                                    master_df1_label, master_df2_label) {
  # Limit the master DFs to only significant hits
  master_df1_sig <- master_df1[master_df1$q.value < qval_thres,]
  master_df2_sig <- master_df2[master_df2$q.value < qval_thres,]
  
  # Plot a Venn diagram of the intersection
  myCol <- brewer.pal(3, "Pastel1")
  grid.newpage()
  v <- venn.diagram(x = list(master_df1_sig$T_k.name, master_df2_sig$T_k.name), 
                    filename = NULL, #paste(master_df1_label, paste(master_df2_label, "BRCA.LumAB.png", sep = "_"), sep = "_"),
                    #lwd = 2, lty = 'blank', col = c("#440154ff", '#21908dff'),  
                    cex = 3, fill = myCol[1:2],
                    fontface = "bold", fontfamily = "sans",
                    cat.cex = 1.5, cat.fontface = "bold", #cat.default.pos = "outer", 
                    cat.fontfamily = "sans", #rotation = 1,
                    category.names = c(master_df1_label, master_df2_label),
                    hyper.test = TRUE, lower.tail = FALSE)
  grid.draw(v)
  
  # Get the intersection between the top genes
  intersecting_genes <- intersect(master_df1_sig$T_k.name, master_df2_sig$T_k.name)
  num_intersecting_genes <- length(intersecting_genes)
  print(paste("Number of Significant Mutation Hits:", nrow(master_df1_sig)))
  print(paste("Number of Significant CNA Hits:", nrow(master_df2_sig)))
  print(paste("Length Overlap:", num_intersecting_genes))
  
  # Use the hypergeometric distribution to determine if the overlap is more
  # than we would expect by chance
  print(phyper(q = num_intersecting_genes - 1, m = nrow(master_df1_sig),
               n = size_targgene_pool - nrow(master_df1_sig), k = nrow(master_df2_sig),
               lower.tail = FALSE))
}

master_df_tp53_mut <- read.csv(paste0(main_path, "Linear Model/TP53/Tumor_Only/eQTL/output_results_LumAB_P53_metabolicTargs_iprotein_tmm_CNAbucket_justDel_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_uncorrected_MUT.csv"), header = TRUE, check.names = FALSE)
master_df_tp53_cna <- read.csv(paste0(main_path, "Linear Model/TP53/Tumor_Only/eQTL/output_results_LumAB_P53_metabolicTargs_iprotein_tmm_CNAbucket_justDel_methMRaw_cibersortTotalFrac_rmCis_allButCancerType_uncorrected_CNA.csv"), header = TRUE, check.names = FALSE)

# Call function
# Num. metabolic targets in Recon3D: 1837
get_statistical_overlap_of_top_hits(master_df_tp53_mut, master_df_tp53_cna, 0.1, 1837, "TP53 Mutation", "TP53 Deletion")


############################################################
############################################################
#### VISUALIZE THE RELATIVE RANKS OF TOP HITS FOR TWO REGPROTS
############################################################
############################################################
#' Creates a line graph with two lines - one for each regulatory protein
#' of interest. Plots the rank (of the top hit from the master DF) vs.
#' the cumulative sum to that rank of hits that were for the given regprot.
#' @param master_df a DF produced as output from the linear model, MHT corrected
#' @param regprot1 the gene name of the first regulatory protein of interest
#' @param regprot2 the gene name of the second regulatory protein of interest
visualize_tophit_relative_ranks <- function(master_df, regprot1, regprot2) {
  
  regprot1_hits <- unlist(lapply(1:nrow(master_df), function(i) 
    ifelse(master_df$R_i.name[i] == regprot1, 1, 0)))
  regprot2_hits <- unlist(lapply(1:nrow(master_df), function(i) 
    ifelse(master_df$R_i.name[i] == regprot2, 1, 0)))
  
  regprot1_cumul <- unlist(lapply(1:length(regprot1_hits), function(i) {
    return(sum(as.numeric(regprot1_hits[1:i])))
  }))
  regprot2_cumul <- unlist(lapply(1:length(regprot2_hits), function(i) {
    return(sum(as.numeric(regprot2_hits[1:i])))
  }))
  
  df_tmp <- data.frame(regprot1 = regprot1_cumul, regprot2 = regprot2_cumul)
  df_tmp$rank <- 1:nrow(df_tmp)
  
  plot(df_tmp$rank, df_tmp[,1], col = "orange", lty = 1, xlab = "rank", 
       ylab = "cumulative # of top hits")
  points(df_tmp$rank, df_tmp[,2], col = "darkgreen", lty = 1)
  
  legend(x = "bottomright", y=NULL, legend = c(regprot1, regprot2), fill = c("orange", "darkgreen"))
}

visualize_tophit_relative_ranks(master_df_mut_corrected, "TP53", "PIK3CA")


#' Creates a line graph with two lines - one for each regulatory protein
#' of interest. Plots the rank (of the top hit from the master DF) vs.
#' the Beta to that rank of hits that were for the given regprot.
#' @param master_df a DF produced as output from the linear model, MHT corrected
#' @param regprot1 the gene name of the first regulatory protein of interest
#' @param regprot2 the gene name of the second regulatory protein of interest
visualize_tophit_relative_Betas <- function(master_df, regprot1, regprot2) {
  
  regprot1_hits <- unlist(lapply(1:nrow(master_df), function(i) 
    ifelse(master_df$R_i.name[i] == regprot1, 1, 0)))
  regprot2_hits <- unlist(lapply(1:nrow(master_df), function(i) 
    ifelse(master_df$R_i.name[i] == regprot2, 1, 0)))
  
  get_beta_by_rank <- function(regprot_hits, master_df) {
    prev_Beta <- 0
    regprot_cumul <- c()
    for (i in 1:length(regprot_hits)) {
      if(regprot_hits[i] == 0) {
        regprot_cumul <- c(regprot_cumul, prev_Beta)
      } else {
        beta <- master_df$estimate[i]
        regprot_cumul <- c(regprot_cumul, beta)
        prev_Beta <- beta
      }
    }
    return(regprot_cumul)
  }
  
  regprot1_cumul <- get_beta_by_rank(regprot1_hits, master_df)
  regprot2_cumul <- get_beta_by_rank(regprot2_hits, master_df)
  
  df_tmp <- data.frame(regprot1 = regprot1_cumul, regprot2 = regprot2_cumul)
  df_tmp$rank <- 1:nrow(df_tmp)
  
  plot(df_tmp$rank, df_tmp[,1], col = "orange", lty = 1, xlab = "rank", 
       ylab = "Relative Betas")
  points(df_tmp$rank, df_tmp[,2], col = "darkgreen", lty = 1)
  abline(a = 0, b = 0)
  
  legend(x = "bottomright", y=NULL, legend = c(regprot1, regprot2), fill = c("orange", "darkgreen"))
}

visualize_tophit_relative_Betas(master_df_mut_corrected, "TP53", "PIK3CA")


############################################################
############################################################
#### VISUALIZE ON RELEVANT PATHWAYS
############################################################
############################################################
#' Visualize our results on pathways using the Pathview package
#' Link to vignette: https://pathview.r-forge.r-project.org/pathview.pdf
#' @param master_df_corrected output master DF with target gene names and q-values
#' @param pathway_id the Pathview pathway ID ("0XXXX")
#' @param label a label for visualization
visualize_pathway <- function(master_df_corrected, pathway_id, label) {
  gene_data <- master_df_corrected$estimate
  names(gene_data) <- master_df_corrected$T_k.name
  
  # Run using Pathview
  pv.out <- pathview(gene.data = gene_data, pathway.id = pathway_id,
                     species = "hsa", out.suffix = label, kegg.native = T,
                     gene.idtype = "SYMBOL")
}

pathway_tp53 <- "04115"
pathway_breastcancer <- "05224"
pathway_cancerCarbonMetabolism <- "05230"
pathway_cancerCholineMetabolism <- "05231"
pathway_pi3k_akt1 <- "04151"

visualize_pathway(master_df_corrected, pathway_tp53, "cancerRelated_metabolicTargs")
visualize_pathway(master_df_corrected, pathway_breastcancer, "cancerRelated_metabolicTargs")
visualize_pathway(master_df_corrected, pathway_cancerCarbonMetabolism, "cancerRelated_metabolicTargs")
visualize_pathway(master_df_corrected, pathway_cancerCholineMetabolism, "cancerRelated_metabolicTargs")
visualize_pathway(master_df_corrected, pathway_pi3k_akt1, "topDrivers_allGeneTargs")


############################################################
############################################################
####  DETERMINE WHAT COVARIATES ARE IMPORTANT
############################################################
############################################################
# Remove the (Intercept) terms
master_df <- master_df[master_df$term != "(Intercept)",]

# Of the top X remaining tests, get what the terms are and plot the proportions
x <- 1000
master_df_topx <- master_df[1:x,]
terms <- unique(master_df_topx$term)
terms_counts <- unlist(lapply(terms, function(x) nrow(master_df_topx[master_df_topx$term == x,])))
terms_counts_df <- data.frame('Term' = terms, 'freq' = terms_counts)
pie(terms_counts_df$freq, labels = terms_counts_df$Term, main = paste("Most Significant Covariates (Top", paste(x, "from All Tests)")))

# Visualize this as a barplot instead
ggplot(data = terms_counts_df[1:7,], aes(x = reorder(Term, -freq), y = freq, fill = Term)) + 
  geom_bar(position = "dodge", stat = "identity") + theme_minimal() + scale_fill_nejm() +
  xlab("Covariate") + ylab(paste0("\n", paste0("Frequency Among Top ", x))) +
  theme(axis.text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust =1, hjust=1), 
        axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14, face="bold"),
        legend.text = element_text(size=12))

# Merge bucketed variables and fix names
unique_terms <- unique(unlist(lapply(terms_counts_df$Term, function(x) {
  if(grepl("_b", x)) {
    spl_x <- unlist(strsplit(x, "_", fixed = T))
    return(paste(spl_x[1:(length(spl_x)-1)], collapse = "_"))
  } else {return(x)}
})))
unique_terms_counts <- unlist(lapply(unique_terms, function(t) 
  sum(terms_counts_df[grepl(t, terms_counts_df$Term), 'freq'])))
names(unique_terms_counts) <- unique_terms
unique_terms_df <- as.data.frame(unique_terms_counts)

# Manually adjust names
unique_terms_df$Term <- c("Methylation Status (Target)", "Frac of Immune Cell Infiltration", 
                               "Cancer Type + Subtype", "CNA Status (Target)", "TP53 CNA Status") #, 
                               #"PIK3CA CNA Status", "Gender")
colnames(unique_terms_df)[1] <- "Frequency"

# Visualize
ggplot(data = unique_terms_df, aes(x = reorder(Term, -Frequency), y = Frequency, fill = Term)) + 
  geom_bar(position = "dodge", stat = "identity") + theme_minimal() + scale_fill_nejm() +
  xlab("Covariate") + ylab(paste0("\n", paste0("Frequency Among Top ", x))) +
  theme(axis.text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust =1, hjust=1), 
        axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14, face="bold"),
        legend.text = element_text(size=12), legend.position = c(0.75,0.75))

# Get the proportions of all the tests with q-value <0.05
master_df_sig <- master_df[master_df$p.value < 0.05,]
terms <- unique(master_df_sig$term)
terms_counts <- unlist(lapply(terms, function(x) nrow(master_df_sig[master_df_sig$term == x,])))
terms_counts_df <- data.frame('term' = terms, 'freq' = terms_counts)
pie(terms_counts_df$freq, labels = terms_counts_df$term, main = "Categories of Significant Covariates (All Tests)")


############################################################
############################################################
#### FOR ANALYSES WITH MULTIPLE REGULATORY PROTEINS OF 
#### INTEREST, CREATE PIE GRAPHS SHOWING THE FREQUENCY OF
#### HITS FOR EACH DRIVER IN THE TOP N% OF HITS
############################################################
############################################################
#' Creates pie charts of number of hits in the top n % for 
#' each driver gene/ regulatory protein (R_i) of interest
#' @param master_df output master DF with target gene names and q-values
#' @param tophit_thres a percentage or q-value threshold (e.g. 0.05) within which to 
#' consider a hit "significant" or important enough for inclusion in pie chart
#' @param perc_or_qval_or_ss whether the threshold is a percentage or a q-value or a stability score
#' @param all_genes_id_conv a conversion table to convert R_i uniprot IDs to 
#' intelligible gene names
create_driver_pie <- function(master_df, tophit_thres, perc_or_qval_or_ss, all_genes_id_conv) {
  
  freq.table.m <- generate_frequency_input_table(master_df, tophit_thres, perc_or_qval_or_ss, all_genes_id_conv)
  
  ggplot(freq.table.m, aes(x = "", y = Frequency, fill = Driver)) + geom_col(color = "black") +
    coord_polar(theta = "y") + theme_void() + 
    geom_text(aes(label = Frequency), position = position_stack(vjust = 0.5),
              fontface = "bold", size = 6) +
    scale_fill_manual(values = c("TP53" = "#0072B5FF", "EP300" = "#7876B1FF", "PIK3CA" = "#BC3C29FF", "KRAS" = "#20854EFF",
                                 "CSMD3" = "#FFDC91FF", "IDH1" = "#EE4C97FF", "CIC" = "beige", "CTNNB1" = "cornflowerblue", 
                                 "PPP2R1A" = "lightpink", "FLG" = "gray", "KMT2C" = "seagreen2", "GATA3" = "lightpink", "PDE4DIP" = "#E18727FF", "MAP3K1" ="gray",
                                 "FBXW7" = "#E18727FF", "NSD1" = "mediumaquamarine", "SETD2" = "palevioletred4", "MET" = "lightpink", "STK11" = "wheat3",
                                 "SMAD4" = "lightyellow", "HRAS" = "mediumorchid3", "SPOP" = "#FFDC91FF", "BRAF" = "seagreen2", "NRAS" = "lightgray", "ARID1A" = "palevioletred2",
                                 "FGFR3" = "seagreen3", "KDM6A" = "lightblue", "CASP8" = "wheat2", "NOTCH1" = "lightpink3", "PTEN" = "red3", "SMARCA4" = "cyan", "NFE2L2" = "orange",
                                 "TSC1" = "cornflowerblue", "RB1" = "lightblue", "STAG2" = "cyan", "CREBBP" = "red2", "ATM" = "lightyellow", "PBRM1" = "pink2")) +
    theme(legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12)) #, legend.position = "bottom")
  
}


#' Helper function to generate the input frequency table for drivers, given appropriate
#' significance thresholds
#' @param master_df output master DF with target gene names and q-values
#' @param tophit_thres a percentage or q-value threshold (e.g. 0.05) within which to 
#' consider a hit "significant" or important enough for inclusion in pie chart
#' @param perc_or_qval_or_ss whether the threshold is a percentage or a q-value or a stability score
#' @param all_genes_id_conv a conversion table to convert R_i uniprot IDs to 
#' intelligible gene names
generate_frequency_input_table <- function(master_df, tophit_thres, 
                                           perc_or_qval_or_ss, all_genes_id_conv) {
  if(perc_or_qval_or_ss == "perc") {
    # Get the top n% of hits
    master_df_topn <- master_df[1:(tophit_thres*nrow(master_df)),]
    if('q.value' %in% colnames(master_df)) {
      print(paste("Q-Value at", paste(tophit_thres, paste("threshold:", master_df[tophit_thres*nrow(master_df), 'q.value']))))
    }
    if('stability.score' %in% colnames(master_df)) {
      print(paste("Stability score at", paste(tophit_thres, paste("threshold:", master_df[tophit_thres*nrow(master_df), 'stability.score']))))
    }
  } else if (perc_or_qval_or_ss == "qval") {
    # Get the n hits below the qvalue threshold
    master_df_topn <- master_df[master_df$q.value < tophit_thres,]
    print(paste("Number of hits below the q-value threshold:", nrow(master_df_topn)))
  } else if (perc_or_qval_or_ss == "ss") {
    # Get the n hits greater than or equal to the stability score threshold
    master_df_topn <- master_df[master_df$stability.score >= tophit_thres,]
    print(paste("Number of hits above or equal to threshold:", nrow(master_df_topn)))
  } else {
    print("The only possible thresholding options are perc or qval. Please try again.")
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
      driver_gn <- unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == driver_uniprot, 
                                            'external_gene_name'])
      return(driver_gn)
    }))
    names(freq.table) <- unique_driver_names
    print(freq.table)
    
    #pie(as.integer(freq.table[1,]), labels = colnames(freq.table))
    
    freq.table.m <- melt(freq.table)
    colnames(freq.table.m) <- c("Driver", "Frequency")
    freq.table.m <- freq.table.m[order(freq.table.m$Frequency, decreasing = TRUE),]
    print(freq.table.m)
    
    return(freq.table.m)
  }
  else {return(NA)}
}
#"#BC3C29FF", "#E18727FF", "#FFDC91FF", "khaki1", "seagreen2", "#20854EFF", "mediumaquamarine",
#"cyan2", "skyblue1", "#0072B5FF",  "darkblue", "#7876B1FF", "mediumorchid", "plum1",
#"lightpink", "#EE4C97FF", "palevioletred4", "beige", "wheat3", "saddlebrown", "gray", "black"

# Call function
create_driver_pie(master_df, 0.2, "qval", all_genes_id_conv)


#' Create driver bar -- a more aesthetic alternative to the pie charts (above)
#' @param list_of_master_dfs list of output master DFs with target gene names and q-values,
#' names are cancer types or subtypes
#' @param tophit_thres a percentage or q-value threshold (e.g. 0.05) within which to 
#' consider a hit "significant" or important enough for inclusion in pie chart
#' @param perc_or_qval_or_ss whether the threshold is a percentage or a q-value or a stability score
#' @param all_genes_id_conv a conversion table to convert R_i uniprot IDs to 
#' intelligible gene names
create_driver_bar <- function(list_of_master_dfs, tophit_thres, perc_or_qval_or_ss, all_genes_id_conv) {
  
  freq.tables <- lapply(1:length(list_of_master_dfs), function(i) {
    master_df <- list_of_master_dfs[[i]]
    cancer_type <- names(list_of_master_dfs)[[i]]
    
    # Create the frequency table
    freq.table <- generate_frequency_input_table(master_df, tophit_thres, perc_or_qval_or_ss, all_genes_id_conv)
    if(length(freq.table) < 2) {return(NA)}
    
    # Add a column of percentages
    sum <- sum(freq.table$Frequency)
    if(sum < 10) {return(NA)}
    freq.table$Percentage <- unlist(lapply(freq.table$Frequency, function(x) (x/sum)*100))
    
    # Add the cancer type
    freq.table$Cancer_Type <- cancer_type
    
    # Adjust the column names to include the cancer type
    #colnames(freq.table)[2:3] <- c(paste0("Freq.", cancer_type), paste0("Perc.", cancer_type))
    return(freq.table)
  })
  #freq.table <- Reduce(function(x, y) merge(x, y, by = "Driver", all = T), freq.tables)
  freq.tables <- freq.tables[!is.na(freq.tables)]
  freq.table <- do.call(rbind, freq.tables)
  freq.table[is.na(freq.table)] <- 0
  #print(freq.table)
  
  # Adjust the order of the drivers based on total frequency + number of cancer types
  freq.per.driver <- lapply(unique(freq.table$Driver), function(d) 
    sum(freq.table[freq.table$Driver == d, 'Frequency']))
  names(freq.per.driver) <- unique(freq.table$Driver)
  freq.per.driver.df <- melt(as.data.frame(freq.per.driver))
  colnames(freq.per.driver.df) <- c("Driver", "Total.Frequency")
  
  # Add the number of cancer types as well
  freq.per.driver.df$Num.CTs <- unlist(lapply(freq.per.driver.df$Driver, function(d)
    length(unique(freq.table[freq.table$Driver == d, 'Cancer_Type']))))
  
  freq.per.driver.df <- freq.per.driver.df[with(freq.per.driver.df, order(Num.CTs, Total.Frequency, decreasing = T)),]
  freq.per.driver.df$Rank <- 1:nrow(freq.per.driver.df)
  print(freq.per.driver.df)
  
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
    print(d)
    print(remaining_cancer_types)
    if(length(remaining_cancer_types) == 0) {break}
    sub <- freq.table[(freq.table$Driver == d) & 
                        (freq.table$Cancer_Type %fin% remaining_cancer_types),]
    new_ranks <- c()
    if(length(ranks) > 0) {
      new_ranks <- rank(desc(sub$Percentage)) + max(as.numeric(ranks))
    } else {new_ranks <- rank(desc(sub$Percentage))}
    names(new_ranks) <- sub$Cancer_Type
    print(new_ranks)
    ranks <- c(ranks, new_ranks)
    print(ranks)
    remaining_cancer_types <- setdiff(remaining_cancer_types, sub$Cancer_Type)
  }
  #print(ranks)
  freq.table$rank_from_driver <- unlist(lapply(freq.table$Cancer_Type, function(ct)
    ranks[names(ranks) == ct]))
  
  #freq.table.list <- c()
  #for(i in 1:length(unique(freq.table$Driver))) {
  #  d <- unique(freq.table$Driver)[i]
  #  sub <- freq.table[freq.table$Driver == d,]
  #  sub$rank_within_driver <- rank(sub$Percentage) 
  #  freq.table.list[[i]] <- sub
  #}
  #freq.table.upd <- do.call(rbind, freq.table.list)
  #freq.table.upd <- freq.table.upd[order(freq.table.upd$Rank),]
  #lvls <- names(sort(tapply(freq.table$Driver == "TP53", freq.table$Cancer_Type, mean), decreasing = T))
  #print(lvls)
  #print(freq.table)
  
  # Option to not show the legend for small, unreadable boxes
  
  
  ggplot(freq.table, aes(x = fct_reorder(Cancer_Type, rank_from_driver), y = Percentage, 
                         fill = fct_reorder(Driver, Rank, .desc = T), label = Frequency.Label)) + 
    geom_bar(position="stack", stat="identity", color = "black") +
    theme_minimal() + xlab("Cancer Type") + ylab(paste0("Percentage of Hits q<", tophit_thres)) +
    scale_fill_manual(values = c("TP53" = "#0072B5FF", "PIK3CA" = "#BC3C29FF", "KRAS" = "#20854EFF", "IDH1" = "#FFDC91FF", 
                                 "CTNNB1" = "#E18727FF", "BRAF" = "#7876B1FF", "NSD1" = "#EE4C9799", "SPOP" = "cyan3",
                                 "NRAS" = "mediumaquamarine", "NFE2L2" = "beige", "NOTCH1" = "slategray1", 
                                 "FBXW7" = "palevioletred3", "FGFR3" = "darkolivegreen3", "STK11" = "orange4", 
                                 "CASP8" = "#6F99AD99", "TSC1" = "cornflowerblue", "RB1" = "yellow3", "SMARCA4" = "thistle2",
                                 "HRAS" = "honeydew2", "SETD2" = "tan", "CIC" = "coral3", "PTEN" = "plum3",
                                 "FGFR2" = "palegreen2", "ZNF521" = "lightsalmon", "MET" = "seagreen3",
                                 "SMAD4" = "slateblue3", "NF1" = "orchid4", "CREBBP" = "paleturquoise2", "PBRM1" = "lightcoral",
                                 "EGFR" = "navy", "ATM" = "white", "STAG2" = "greenyellow", "GNAS" = "gainsboro", 
                                 "EP300" = "darkseagreen", "PPP2R1A" = "deepskyblue"),
                      breaks = as.character(freq.per.driver.df$Driver),
                      name = "Driver") +
    geom_text(size = 3, position = position_stack(vjust = 0.5)) +
    theme(axis.title = element_text(face = "bold", size = 14), 
          axis.text = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12))
  
}

create_driver_bar(top_drivers_0.05, 0.2, "qval", all_genes_id_conv)



# OLD COLOR SCHEME:
scale_fill_manual(values = c("TP53" = "#0072B5FF", "EP300" = "mediumorchid1", "PIK3CA" = "#BC3C29FF", "KRAS" = "#20854EFF",
                             "CSMD3" = "#FFDC91FF", "IDH1" = "#EE4C97FF", "CIC" = "beige", "CTNNB1" = "cornflowerblue", 
                             "PPP2R1A" = "lightpink", "FLG" = "gray", "GATA3" = "lightpink", "PDE4DIP" = "orange", "MAP3K1" ="gray4",
                             "FBXW7" = "mediumaquamarine", "NSD1" = "#E18727FF", "SETD2" = "palevioletred3", "MET" = "#7876B1FF", "STK11" = "palevioletred2",
                             "SMAD4" = "lightsalmon", "HRAS" = "mediumorchid3", "SPOP" = "seagreen3", "BRAF" = "#FFDC91FF", "NRAS" = "lightgray", "ARID1A" = "black",
                             "FGFR3" = "seagreen1", "KDM6A" = "lightblue", "CASP8" = "wheat2", "NOTCH1" = "lightpink3", "PTEN" = "red4", "SMARCA4" = "cyan2", "NFE2L2" = "yellow2",
                             "TSC1" = "palevioletred4", "RB1" = "lightblue", "STAG2" = "cyan3", "CREBBP" = "red2", "ATM" = "orange3", "GNAS" = "yellow3",
                             "ZNF521" = "tan", "NF1" = "gray", "FGFR2" = "coral3", "EGFR" = "darkblue", "PBRM1" = "darkolivegreen2")) #"KMT2C" = "seagreen2"
                  


############################################################
############################################################
#### RECOMBINE P-VALUES USING FISHER'S OR STOUFFER's METHODS
#### FOR TESTS ACROSS MULTIPLE SUBTYPES OR CANCER TYPES
############################################################
############################################################
#' Using the metap package, recombine the p-values from across 
#' cancer types of subtypes to explore global effects
#' @param list_of_master_dfs a list of data frames, with the names of the list
#' corresponding to the cancer type or subtype; should already be subsetted to
#' only include one R_i
recombine_pvals <- function(list_of_master_dfs) {
  # Get all the unique target genes 
  unique_targets <- unique(unlist(lapply(list_of_master_dfs, function(x) 
    unlist(x$T_k.name))))
  
  # For each of these, get the p-values from across the various cancer types or subtypes
  new_pvals <- unlist(lapply(unique_targets, function(targ) {
    print(targ)
    pvals <- unlist(lapply(list_of_master_dfs, function(df) {
      if(targ %in% df$T_k.name) {
        return(df[df$T_k.name == targ, 'p.value'])
      }
    }))
    # Do the correction using these p-values
    pval_new <- metap::sumlog(pvals, log.p = FALSE)$p
    print(pval_new)
    return(pval_new)
  }))
  
  # Create a new data frame with all targets and their p-values
  new_df <- data.frame("T_k.name" = unique_targets, "p.value" = new_pvals)
  
  # Multiple hypothesis testing correction
  new_df$q.value <- qvalue(new_df$p.value)$qvalues
  
  return(new_df)
}

# Call function 
new_master_0.05 <- recombine_pvals(list_of_master_dfs_0.05)

new_master_0.05 <- new_master_0.05[order(new_master_0.05$q.value),]

new_master_0.05_tp53$R_i.name <- "TP53"

############################################################
############################################################
#### DETERMINE WHICH VARIABLES WERE REMOVED MOST OFTEN WHEN
#### CORRECTING FOR COLLINEARITY
############################################################
############################################################
#' Uses the files output from the correct_collinearity functions to calculate 
#' the number of times each variable was removed due to collinearity
#' @param elim_vars_df a data frame with the variables removed from the regression 
#' for each target gene
calculate_collinearity_removal_freq <- function(elim_vars_df) {
  elim_vars <- unlist(lapply(elim_vars_df$Variables.Removed, function(x) 
    unlist(strsplit(x, ",", fixed = T))))
  print(unique(elim_vars))
  #elim_vars <- unlist(lapply(elim_vars, function(x) 
    #ifelse(grepl("Cancer_type", x), "Cancer_type", x)))
  
  elim_vars_table <- table(elim_vars)
  print(elim_vars_table)

}


############################################################
############################################################
#### FOR A GIVEN PAN-CANCER HIT, DETERMINE WHICH INDIVIDUAL
#### CANCER TYPES IT IS ALSO SIGNIFICANT IN (AND DIRECTIONALITY)
############################################################
############################################################
#' For a given pan-cancer hit, see which individual cancer types it is 
#' significant in (for potential experiments). Returns a data frame with
#' this information and prints it to the console.
#' @param top_drivers_0.05 a list of master DFs from individual cancer types,
#' with name corresponding to the cancer type
#' @param gois list of gene names that are pan-cancer targets of a given driver
#' @param driver the gene name of the driver for which this gene is a hit
#' @param qval_thres a q-value threshold for determining hit significance
get_per_cancer_analysis_of_pc_hit <- function(top_drivers_0.05, gois, driver, qval_thres) {
  pvals_and_betas <- lapply(1:length(top_drivers_0.05), function(i) {
    df <- top_drivers_0.05[[i]]
    ct <- names(top_drivers_0.05)[i]
    
    if(driver %fin% df$R_i.name) {
      df_sub <- df[(df$T_k.name %fin% gois) & (df$R_i.name == driver), ]
      df_sub_sig <- df_sub[df_sub$q.value < qval_thres]
      if(nrow(df_sub_sig) > 0) {
        to_return <- df_sub_sig[, c('T_k.name', 'R_i.name', 'estimate', 'q.value')]
        to_return$cancer_type <- ct
        return(to_return)
      } else {return(NA)}
    } else {return(NA)}
  })
  pvals_and_betas <- pvals_and_betas[!is.na(pvals_and_betas)]
  
  pvals_and_betas_df <- do.call(rbind, pvals_and_betas)
  
  print(paste("For driver and target combinations for", driver))
  print(pvals_and_betas_df)
  
  return(pvals_and_betas_df)
}

kras_metabolic_pc_hits <- c("NT5E", "TBXAS1", "SLC27A1", "FPGS", "CANT1", "GALNT10", "ALDH3B1", "FAXDC2",
                            "SAT1", "DPYSL2", "AGPAT2", "GMDS", "SRD5A3", "SLC18A1", "SLC26A9", 
                            "B4GALNT3", "ACOX2", "SLC14A1", "ATP6V0A4", "MGLL", "PIP5K1C", "PLCB3",
                            "SLC12A2", "CYP3A5", "ST3GAL5", "PFKP", "CPS1", "MTMR14", "CHPF2",
                            "ACO1", "ABCC3", "AQP5", "GNE", "PDHB", "CA8", "SLCO4A1", "SLC3A1",
                            "SLC37A1", "SLC35C1", "DOLPP1")
pik3ca_metabolic_pc_hits <- c("PIK3R1", "PIK3R3", "CYP4X1", "PLCE1", "SGPL1", "PPIP5K2",
                              "HADH", "CRAT", "MYO5B", "IDH2", "HMGCR", "PIGN",
                              "GPX8", "ASRGL1", "SLC37A1", "UGT2B11", "SPTLC2")
idh1_metabolic_pc_hits <- c("INPPL1", "CYP2U1", "CHST14", "B4GALT1", "SLC14A1", "MTHFD2",
                              "AHCYL1", "LPL", "UCP2", "IMPDH2", "DERA", "CBR4", "SLC12A4",
                            "SLC11A1", "ABCA1", "SLC35B2")
tp53_metabolic_pc_hits <-c("ENO1", "TK1", "SLC35A2", "G6PD", "RRM2", "PRPS1", 
                           "GMPS", "PGK1", "PDHA1", "PRDX1", "CAD", "SMS", "SRM", "IL4I1", "MTHFD2", "ME1",
                           "HPRT1", "RPIA", "ABCD1", "SUV39H1", "GLO1", "PPAT", "PPT1",
                           "MTHFD1L", "CHKA", "PGD", "SLC25A19", "DEGS1", "GAPDH", "PSAT1",
                           "SPHK1", "CHST11", "RPE", "FTL", "SLC29A4", "ACYP1", "NOL9",
                           "PYCR1", "AK4", "PDE6D", "PIGU", "B4GALNT1", "CTSA", 
                           "GNPDA1", "GFPT2", "DTYMK", "PFKFB4", "PIGW", "CYP27B1",
                           "SMYD3", "TXNRD1", "OSER1", "SLC7A5", "LBR", "SLC38A1",
                           "B4GALT5", "PTPDC1", "B3GNT4", "PHEX", "ST6GALNAC5", "PTDSS1",
                           "HSD17B10", "CTPS1", "IMPDH1", "PHGDH", "BCAT1", "UMPS",
                           "SMPD4", "FUCA2", "PKM", "MTHFD1L", "ENO2", "ALDOA", "GPI", 
                           "HK1", "GART", "PFAS", "NME3", "CTPS1", "CTPS2", "H6PD", "MTR", "MTHFR", "FGFR3")

get_per_cancer_analysis_of_pc_hit(top_drivers_0.05, kras_metabolic_pc_hits, "KRAS", 0.2)
get_per_cancer_analysis_of_pc_hit(top_drivers_0.05, pik3ca_metabolic_pc_hits, "PIK3CA", 0.2)
get_per_cancer_analysis_of_pc_hit(top_drivers_0.05, idh1_metabolic_pc_hits, "IDH1", 0.2)
get_per_cancer_analysis_of_pc_hit(top_drivers_0.05, tp53_metabolic_pc_hits, "TP53", 0.2)

# Get the genes from a KEGG GSEA file:
print.table(unlist(strsplit(results_gsea_pik3ca_tbl_sig[results_gsea_pik3ca_tbl_sig$Description == "Valine, leucine and isoleucine degradation", 'core_enrichment'], "/", fixed = T)), 
            quote=F, row.names = F)
# Convert to gene names using the following: https://www.genome.jp/kegg/tool/conv_id.html
tmp <- master_df[(master_df$T_k.name %in% valine_leucine_genes) & (master_df$q.value < 0.01),]
print(tmp[order(tmp$T_k.name),])


# ASIDE: JOSH GRANT FIGURE #
josh_df <- read.csv("C:/Users/sarae/Documents/plotTable_interestingZF_bindingEnrichment_MP.csv", header = TRUE)
josh_df <- read.csv("C:/Users/sarae/Documents/plotTable_interestingZF_bindingEnrichment_CCCP.csv", header = TRUE)

josh_df$zf <- factor(josh_df$zf, levels = josh_df$zf[order(-log(josh_df$fisher.p))])

g <- ggplot(josh_df, aes(x = zf, y = -log(fisher.p), fill = zf)) + 
  geom_bar(stat = 'identity', width = 0.95, show.legend = FALSE, color = "black") + 
  theme_minimal() + scale_fill_manual(values = c("#BC3C29FF", "#0072B5FF", "#20854EFF")) +
  labs(y = "-log(p)", x = "Gene") + #title = "Enrichment for binding (M Phase)") + 
  theme(axis.text = element_text(size = 16, face = "bold"), title = element_text(size = 18),
        axis.title = element_text(size = 18, face = "bold"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'))

g



