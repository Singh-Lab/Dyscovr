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

# NEJM color palatte: https://nanx.me/ggsci/reference/pal_nejm.html


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
  pvals <- results_table$p.value[!is.na(results_table$p.value) & 
                                   !is.infinite(results_table$p.value)]
  hist(pvals, main = "",
       xlab = "p-value", ylab = "Frequency", col = "blueviolet")
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
  qqline(results_table$p.value, col = "steelblue", lwd = 2)
}

#' Function plots a Q-Q plot to visualize the distribution of standard error valeus
#' and assess whether they come from a uniform distribution
#' @param results_table a master DF produced from run_linear_model()
qqplot_stderror <- function(results_table) {
  qqnorm(results_table$std.error, pch = 1, frame = FALSE)
  qqline(results_table$std.error, col = "steelblue", lwd = 2)
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
mh_correct <- function(results_table) {
  
  # Get the qvalue object
  qobj <- NA
  if(length(results_table$p.value) < 100) {
    # With a small number of pvalues we may not be able to accurately estimate pi0,
    # so we set to 1 (the equivalent of B-H correction)
    qobj <- qvalue(p = results_table$p.value, pi0 = 1)
  } else {
    qobj <- qvalue(p = results_table$p.value)
  }
  
  # OPT: plot some useful plots & print some useful information
  plot(qobj)
  print(summary(qobj))
  #print(paste("Pi0 (Propr. of true null hypotheses):", qobj$pi0))
  
  qvals <- qobj$qvalues # extract qvalues
  
  results_table$q.value <- qvals # add the qvalues back to the data frame
  
  # Return the tidied linear model fit with q-values
  return(results_table)
}

# Call this function
master_df_mut_corrected <- mh_correct(master_df_mut)
master_df_cna_corrected <- mh_correct(master_df_cna)

#Q-Value Visualization (I-Protein, FPKM, rawCNA, methBeta, iciTotFrac)


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
  master_df_sig$T_k.name <- unlist(lapply(master_df_sig$T_k, function(x) 
    paste(unique(unlist(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == unlist(strsplit(x, ";", fixed = TRUE))[1], 
                                   'external_gene_name'])), collapse = ";")))
  
  # Add a column for the regulatory protein name
  master_df_sig$R_i.name <- unlist(lapply(master_df_sig$R_i, function(x) 
    paste(unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == unlist(strsplit(x, ";", fixed = TRUE))[1], 
                                   'external_gene_name']), collapse = ";")))
  return(master_df_sig)
}

# Call this function
master_df_mut_corrected <- add_targ_regprot_gns(master_df_mut_corrected, all_genes_id_conv)
master_df_cna_corrected <- add_targ_regprot_gns(master_df_cna_corrected, all_genes_id_conv)

# Write this to a new file
fwrite(master_df_mut_corrected, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/eQTL/output_results_P53_chipeat_iprotein_TMM_rawCNA_cibersortTotFrac_corrected_MUT.csv", sep = ""))
fwrite(master_df_cna_corrected, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/eQTL/output_results_P53_chipeat_iprotein_TMM_rawCNA_cibersortTotFrac_corrected_CNA.csv", sep = ""))


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
#' to the names of the regulatory proteins of interest
#' @param q_thres a q-value threshold for significance
#' @param list_of_silent_dfs OPT: a corresponding list of silent DFs to include in the chart
create_num_sig_hits_barplot <- function(list_of_master_dfs, q_thres, list_of_silent_dfs) {
  # Get the number of significant hits for each master DF
  sig_hits_df <- data.frame("Gene" = names(list_of_master_dfs), 
                            "Num.Signif.Hits.q" = rep(0, times = length(list_of_master_dfs)))
  colnames(sig_hits_df)[2] <- paste0(colnames(sig_hits_df)[2], q_thres)
  sig_hits_df[,2] <- unlist(lapply(list_of_master_dfs, function(df) {
    return(nrow(df[df$q.value < q_thres,]))
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
      scale_fill_manual(values = c("#BC3C29FF", "#0072B5FF", "gray")) + #scale_color_nejm() +
      xlab("Gene") + ylab(paste0("\n", paste0("Number of hits (q < ", paste0(q_thres, ")")))) +
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
x <- 500
master_df_topx <- master_df[1:x,]
terms <- unique(master_df_topx$term)
terms_counts <- unlist(lapply(terms, function(x) nrow(master_df_topx[master_df_topx$term == x,])))
terms_counts_df <- data.frame('term' = terms, 'freq' = terms_counts)
pie(terms_counts_df$freq, labels = terms_counts_df$term, main = paste("Most Significant Covariates (Top", paste(x, "from All Tests)")))


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
  
  pie(as.integer(freq.table[1,]), labels = colnames(freq.table))
  
}


create_driver_pie(master_df, 0.05, "qval", all_genes_id_conv)

















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



