############################################################
# Code to Create Figure 2 Visualizations
# Written by Sara Geraghty
# PUBLICATION INFORMATION
############################################################

library(data.table)
library(ggplot2)
library(KSgeneral)
library(dorothea)
library(STRINGdb)
library(igraph)
library(GOSemSim)
library(ggrepel)
library(ggsci)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(org.Hs.eg.db)

keytypes(org.Hs.eg.db)

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

outfn_wttn <- "res_top_0.05_allGenes_quantile_rawCNA_methMRaw_3PCs_Nonsyn.All.Vogel.elim.vif5.sp0.7_corrected_MUT.csv"
pc_allGenes_wttn <- read.csv(paste0(PATH, outfn_wttn))


############################################################
### PART A: BARPLOT OF HITS AT Q < 0.01
############################################################
qval_thres <- 0.01
pc_allGenes_sig <- pc_allGenes_wttn[pc_allGenes_wttn$q.value < qval_thres,]
pc_allGenes_sig_freq <- melt(table(pc_allGenes_sig$R_i.name))
colnames(pc_allGenes_sig_freq) <- c("Driver", "Num.Hits")
pc_allGenes_sig_freq <- pc_allGenes_sig_freq[order(pc_allGenes_sig_freq$Num.Hits),]
ggplot(pc_allGenes_sig_freq, aes(y = Num.Hits, fill = Driver, 
                             x = reorder(Driver, -Num.Hits, mean))) + 
  geom_bar(position = "dodge", width = 0.95, stat = "identity", 
           show.legend = F, color = "black") + 
  scale_fill_manual(values = c("#FFDC91FF", "#20854EFF", "#BC3C29FF", 
                               "#0072B5FF", "gray")) + 
  xlab("Driver") + 
  ylab(paste0("\n", paste0("Number of hits (q < ", paste0(qval_thres, ")")))) +
  theme(axis.text = element_text(face="bold", size = 16), 
        axis.title=element_text(size=18, face="bold"), 
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = 'white'))

  
############################################################
### PART B-C: BARPLOT OF ENRICHMENT IN CGC AND TF TARGETS
############################################################  
#' Use a one-sided Fisher's Exact Test/ Hypergeometric (HG) test or a 
#' Kolmogorov-Smirnov (K-S) test to test for statistical enrichment
#' @param master_df a master DF produced from Dyscovr that has q-values
#' @param known_targs a vector of known targets for a gene of interest (master
#' DF should already be subsetted to just this gene); ex. ChIP-eat targets, STRING
#' @param type the type of test, either "fisher's exact", "k-s", or "both"
#' @param qval_thres q-value threshold for HG enrichment
#' @param top_n option to specify a top n hits, rather than a q-value threshold
compute_statistical_enrichment <- function(master_df, known_targs, type, 
                                           qval_thres, top_n) {
  final_res <- NA
  
  if ((type == "fisher's exact") | (type == "both")) {
    sig_targets <- c()
    nonsig_targets <- c()
    
    if(!is.na(qval_thres)) {
      sig_targets <- unique(master_df[master_df$q.value < qval_thres, 
                                      'T_k.name'])
      nonsig_targets <- unique(master_df[master_df$q.value > qval_thres, 
                                         'T_k.name'])
    } else if(!is.na(top_n)) {
      sig_targets <- unique(master_df[1:top_n, 'T_k.name'])
      nonsig_targets <- unique(master_df[(top_n+1):nrow(master_df), 'T_k.name'])
    } else {
      print("Must specify a q-value threshold or top N in order to perform a 
            Fischer's exact test. Please try again.")
      return(NA)
    }
    targets <- unique(master_df$T_k.name)
    
    num_sig_targets <- length(sig_targets)
    num_nonsig_targets <- length(nonsig_targets)
    
    known_targs_i <- intersect(known_targs, targets)
    notknown_targs_i <- setdiff(targets, known_targs)
    
    num_known_targets <- length(known_targs_i)
    num_notknown_targets <- length(notknown_targs_i)
    
    num_known_sig_targets <- length(intersect(sig_targets, known_targs_i))
    num_known_nonsig_targets <- length(intersect(nonsig_targets, known_targs_i))
    num_notknown_sig_targets <- length(intersect(sig_targets, notknown_targs_i))
    num_notknown_nonsig_targets <- length(intersect(nonsig_targets, notknown_targs_i))
    
    contigency_table <- matrix(c(num_known_sig_targets, 
                                 num_known_nonsig_targets,
                                 num_notknown_sig_targets, 
                                 num_notknown_nonsig_targets),
                               nrow = 2, ncol = 2)
    print(contigency_table)
    
    fisher_res <- fisher.test(contigency_table, alternative = "greater")
    print(fisher_res$p.value)
    
    final_res <- fisher_res
  }
  if ((type == "k-s") | (type == "both")) {
    # Create a uniform distribution the same length as the number of genes
    uniform <- 1:length(master_df$T_k.name)
    
    # Get the ranks of the known targets within our ranked list
    known_targs <- setdiff(known_targs, setdiff(known_targs, master_df$T_k.name))
    ranks_of_targs <- unlist(lapply(known_targs, function(x) {
      if(x %fin% master_df$T_k.name) {return(min(which(master_df$T_k.name == x,)))}
      else {return(NA)}
    }))
    ranks_of_targs <- ranks_of_targs[!(is.infinite(ranks_of_targs) | 
                                         (is.na(ranks_of_targs)))]

    # Perform the K-S test, with the uniform distribution
    ks_res <- ks.test(ranks_of_targs, uniform, alternative = "greater")
    print(ks_res$p.value)
    
    # Alternative to print the K-S statistic exactly, using the KSgeneral package
    # (Note that this is now a two-sided test)
    if(ks_res$p.value == 0) {
      print("Printing KSgeneral two-sided K-S p-value:")
      ks_res_exact <- disc_ks_test(ranks_of_targs, ecdf(uniform), exact = T,
                                   alternative = "greater")
      print(ks_res_exact$p.value)
      ks_res <- ks_res_exact
    }
    
    if(type == "both") {final_res <- list(final_res, ks_res)}
    else {final_res <- ks_res}
  }
  if (!(type %in% c("fisher's exact", "k-s", "both"))) {
    print("Only implemented for 'Fisher's exact' or 'k-s' tests.")
    return(0)
  }
  return(final_res)
}

### CANCER GENE SETS ###
# Import a table containing a compiled list of known cancer genes
# CGC genes are from Cancer Gene Census genes (from COSMIC) downloaded January 21, 2019
# COSMIC: https://www.sanger.ac.uk/tool/cosmic/
# Vogelstein genes are from  Vogelstein et al. (Science 2013), 
# doi: 10.1126/science.1235122, Tables S2A, S2B, S3A, S3B, S3C and S4
known_cancer_genes_table <- read.table(paste0(PATH, "Validation_Files/GRCh38_driver_gene_list.tsv"), 
                                       sep = "\t", header = T, check.names = F, 
                                       comment.char = "#", skip = 10)
known_cancer_genes <- known_cancer_genes_table$primary_gene_names
vogelstein_genes <- known_cancer_genes_table[
  grepl("V", known_cancer_genes_table$cancer_driver_status), 'primary_gene_names']
cgc_genes <- known_cancer_genes_table[
  grepl("C", known_cancer_genes_table$cancer_driver_status), 'primary_gene_names']


### TF EFFECTORS FROM DOROTHEA ###
# https://bioconductor.org/packages/release/data/experiment/vignettes/dorothea/inst/doc/dorothea.html
dorothea_net <- dorothea::dorothea_hs

# PIK3CA
foxo_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("FOXO1|FOXO3|FOXO4|FOXO6", dorothea_net$tf) , 'target']))  #907
srebp_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("SREBF", dorothea_net$tf) , 'target'])) #31
myc_dorothea_targets <- unique(unlist(dorothea_net[
  dorothea_net$tf == "MYC", 'target'])) #386
hif1a_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("HIF1A", dorothea_net$tf) , 'target'])) #148
atf4_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("ATF4", dorothea_net$tf) , 'target'])) #16
nrf2_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("NFE2L2", dorothea_net$tf) , 'target'])) #11

pik3ca_tf_targets <- unique(c(foxo_dorothea_targets, srebp_dorothea_targets, 
                              myc_dorothea_targets, hif1a_dorothea_targets, 
                              atf4_dorothea_targets, nrf2_dorothea_targets)) # 1411
compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == "PIK3CA",], 
                               pik3ca_tf_targets, "both", 0.01, NA)

# KRAS
ets1_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("ETS1", dorothea_net$tf) , 'target']))  #145
ets2_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("ETS2", dorothea_net$tf) , 'target']))  #30
etv1_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("ETV1", dorothea_net$tf) , 'target']))  #29
elk1_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("ELK1", dorothea_net$tf) , 'target']))  #31
creb135_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("CREB1|CREB3|CREB5", dorothea_net$tf) , 'target']))  #1779
foxo1_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("FOXO1", dorothea_net$tf) , 'target']))  #43
myc_dorothea_targets <- unique(unlist(dorothea_net[
  dorothea_net$tf == "MYC", 'target'])) #386

kras_tf_targets <- unique(c(ets1_dorothea_targets, ets2_dorothea_targets, 
                            etv1_dorothea_targets, elk1_dorothea_targets, 
                            foxo1_dorothea_targets, myc_dorothea_targets, 
                            creb135_dorothea_targets)) # 2591
compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == "KRAS",], 
                               kras_tf_targets, "both", 0.01, NA)

# IDH1
atf3_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("ATF3", dorothea_net$tf) , 'target']))  #1304
jun_dorothea_targets <- unique(unlist(dorothea_net[
  grepl("JUN", dorothea_net$tf), 'target']))  #176
sox8_dorothea_targets <- unique(unlist(dorothea_net[
  dorothea_net$tf == "SOX8", 'target']))  #802
mycn_dorothea_targets <- unique(unlist(dorothea_net[
  dorothea_net$tf == "MYCN", 'target']))  #13
hey1_dorothea_targets <- unique(unlist(dorothea_net[
  dorothea_net$tf == "HEY1", 'target']))  #422
nr2f2_dorothea_targets <- unique(unlist(dorothea_net[
  dorothea_net$tf == "NR2F2", 'target'])) #18

idh1_tf_targets <- unique(c(atf3_dorothea_targets, jun_dorothea_targets, 
                            sox8_dorothea_targets, mycn_dorothea_targets,
                            hey1_dorothea_targets, nr2f2_dorothea_targets))
compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == "IDH1",], 
                               idh1_tf_targets, "both", 0.01, NA)


#' Create enrichment barplot of -log10(p-values) from enrichment test for each
#' driver gene
#' @param pc_allGenes results DF from Dyscovr, pan-cancer
#' @param drivers vector of the names of the driver genes of interest
#' @param target_set a list of 1, with the entry corresponding to a vector of 
#' names of genes in the target set. Alternatively, if we have a different
#' target set per driver, a list of target sets with names corresponding
#' to driver genes
#' @param label_enrichment_type a label for the type of enrichment test (e.g. "K-S")
#' @param label_target_set a label for the name of the target gene set (e.g. CGC)
#' @param qval_thres a q-value threshold for significance, only relevant for 
create_enrichment_barplot <- function(pc_allGenes, drivers, label_enrichment_type,
                                      target_set, label_target_set, qval_thres) {
  enrichments <- lapply(drivers, function(d) {
    print(d)
    
    e <- NA
    # One target set for all drivers (e.g. CGC genes)
    if(length(target_set) == 1) {
      e <- compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == d,], 
                                          unlist(target_set), 
                                          tolower(label_enrichment_type), 
                                          qval_thres)
    # One target set per driver (e.g. TF targets)
    } else {
      target_set_d <- target_set[[which(names(target_set) == d)]]
      e <- compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == d,], 
                                          target_set_d, tolower(label_enrichment_type), 
                                          qval_thres)
    }

    p <- -log10(e$p.value)
    return(p)
  })
  names(enrichments) <- drivers
  
  input <- melt(as.data.frame(enrichments))
  colnames(input) <- c("Driver", "negLog10pval")
  
  p <- ggplot(input, aes(y = negLog10pval, fill = Driver, 
                    x = reorder(Driver, negLog10pval, mean))) + 
    geom_bar(position = "dodge", width = 0.95, stat = "identity", 
             show.legend = F, color = "black") + 
    scale_fill_manual(values = c("#0072B5FF", "#BC3C29FF", "#20854EFF", "#FFDC91FF")) + 
    #scale_fill_nejm() + 
    xlab("Driver") +
    ylab(paste("-log10(pval) of", 
                paste(label_enrichment_type, 
                      paste("Enrichment,\n", paste(label_target_set, "Genes"))))) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black") + 
    coord_flip() + #theme_minimal() +
    scale_y_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text = element_text(face="bold", size = 14), 
          axis.title=element_text(size=16, face="bold"), 
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = 'white'))
  print(p)
}

drivers <- c("TP53", "PIK3CA", "KRAS", "IDH1")

# Call function (CGC)
create_enrichment_barplot(pc_allGenes, drivers, "K-S", list(cgc_genes), "CGC", 0.01)

# Call function (DoRothEA targets)
tf_targets <- list("TP53" = tp53_dorothea_targets, 
                         "PIK3CA" = pik3ca_tf_targets,
                         "KRAS" = kras_tf_targets,
                         "IDH1" = idh1_tf_targets)
create_enrichment_barplot(pc_allGenes, drivers, "K-S", tf_targets, 
                          "DoRothEA TF Targets (Pooled)", 0.01)


############################################################
### PART D: TP53 CUMULATIVE FRACTION OF KNOWN TARGETS
############################################################  
#' Function to plot the enrichment of various target groups among the significant 
#' gene hits, together on one plot, where the enrichment in the given source is 
#' denoted by color and line form
#' @param master_df a master DF produced from Dyscovr that has q-values
#' @param target_sets_list a list of each of the target sets we are interested 
#' in doing enrichment on, with each set's name/ label as the name for the 
#' list item
#' @param goi a string denoting the gene-of-interest, for labeling graph
#' @param thres a threshold for the number of hit genes we show in the full plot; 
#' if NA then we plot enrichment across all gene targets
#' @param cutout_thres a threshold for the inset/cutout enrichment plot 
#' ('zoomed' in on)
#' @param silent_df OPTIONAL: an additional, corresponding DF from Dyscovr with 
#' silent mutation results
plot_combined_enrichment <- function(master_df, target_sets_list, goi, thres, 
                                     cutout_thres, silent_df) {
  
  master_df <- master_df[order(master_df$q.value),]
  if(length(silent_df) > 1) {
    silent_df <- silent_df[order(silent_df$q.value),]
  }
  
  # For each source, create a 0 and 1 vector for each  target hit to indicate if  
  # it is a a hit in that source
  target_set_vector_list <- lapply(target_sets_list, function(source) {
    source <- setdiff(source, setdiff(source, master_df$T_k.name))
    vect <- unlist(lapply(1:nrow(master_df), function(i) {
      x <- master_df[i, 'T_k.name']
      return(ifelse(x %fin% source, 1, 0))
    }))
    return(vect)
  })
  
  target_set_vector_list_silent <- NA
  if(length(silent_df) > 1) {
    names(target_set_vector_list) <- paste0("MN.", names(target_sets_list))
    silent_df$T_k.name <- as.factor(silent_df$T_k.name)
    target_set_vector_list_silent <- lapply(target_sets_list, function(source) {
      vect <- unlist(lapply(1:nrow(silent_df), function(i) {
        x <- silent_df[i, 'T_k.name']
        return(ifelse(x %fin% source, 1, 0))
      }))
      return(vect)
    })
    names(target_set_vector_list_silent) <- paste0("S.", names(target_sets_list))
  } else {
    names(target_set_vector_list) <- names(target_sets_list)
  }
  
  # Use helper function to get fraction of genes at given rank
  target_set_fraction_df <- get_frac_of_genes_at_rank(target_set_vector_list, 
                                                      master_df)
  target_set_fraction_df_silent <- NA
  if(length(silent_df) > 1) {
    target_set_fraction_df_silent <- get_frac_of_genes_at_rank(
      target_set_vector_list_silent, silent_df)
    target_set_fraction_df <- merge(target_set_fraction_df, 
                                    target_set_fraction_df_silent, by = "Rank")
  }
  
  target_sets_fraction_df_m <- melt(target_set_fraction_df, "Rank")
  colnames(target_sets_fraction_df_m) <- c("Rank", "Source", "Frac")
  target_sets_fraction_df_m$Mutation.Type <- unlist(lapply(
    target_sets_fraction_df_m$Source, function(x) {
      val <- "Nonsynonymous"
      if(grepl("S.", x, fixed = TRUE)) {val <- "Synonymous"}
      return(val)
  }))
  
  # Create inset plot
  p2 <- ggplot(target_sets_fraction_df_m[target_sets_fraction_df_m$Rank %fin% 
                                           1:cutout_thres,], 
               mapping = aes(x = Rank, y = Frac, color = Source)) + 
    geom_line(aes(linetype=Mutation.Type, alpha = Mutation.Type), size = 1.25) +
    scale_x_continuous(limits = c(1,cutout_thres)) + 
    scale_color_nejm() + scale_alpha_manual(values=c(0.9,0.9,0.4,0.4)) +
    scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#BC3C29FF", "#0072B5FF")) + 
    theme(legend.position = "none", axis.title.x = element_blank(), 
          axis.title.y = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x = element_text(size = 14, face = "bold"), 
          axis.text.y =  element_text(size = 14, face = "bold"))
  
  p <- ggplot(target_sets_fraction_df_m[target_sets_fraction_df_m$Rank %fin% 
                                          1:thres,], 
              mapping = aes(x = Rank, y = Frac, color = Source)) + 
    geom_line(aes(linetype=Mutation.Type, alpha = Mutation.Type), size = 1.25) +
    theme(axis.text = element_text(size = 16, face = "bold"), 
          axis.title = element_text(size = 18, face = "bold"), 
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = 'white'), 
          legend.text=element_text(size=14), 
          legend.title = element_text(size = 16)) +
    xlab("Gene Rank") + ylab("Fraction in Set of Known Targets") + 
    scale_x_continuous(limits = c(1,thres)) +
    scale_color_nejm() +  scale_alpha_manual(values=c(0.9,0.9,0.4,0.4)) +
    scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#BC3C29FF", "#0072B5FF"),
                       labels = c(names(target_sets_list), 
                                  rep("", times = length(names(target_sets_list))))) + 
    annotation_custom(ggplotGrob(p2), xmin =(thres-round(thres/1.5)), 
                      xmax = thres, ymin = 0.30, ymax = 1.0) +
    geom_rect(data=target_sets_fraction_df_m[1:cutout_thres,], 
              aes(xmin = 1, xmax = cutout_thres, ymin = 0, ymax = 1, 
                  fill = "gray"), color = NA, fill = alpha("gray", .01))
  return(p)
}

#' Get the fraction of genes at or above the given rank, for each set
#' @param target_set_vector_list a named list of each source with hit/not 
#' binary vector
#' @param master_df master df with ordered top hits
get_frac_of_genes_at_rank <- function(target_set_vector_list, master_df) {
  target_set_fraction_list <- lapply(target_set_vector_list, function(vect) {
    frac_nw_vect <- c()
    count_of_nw_genes <- 0
    for (i in 1:nrow(master_df)) {
      val_of_curr_tg <- vect[i]
      count_of_nw_genes <- count_of_nw_genes + val_of_curr_tg
      frac <- count_of_nw_genes / i
      frac_nw_vect <- c(frac_nw_vect, frac)
    }
    return(frac_nw_vect)
  })
  
  names(target_set_fraction_list) <- names(target_set_vector_list)
  target_sets_fraction_df <- as.data.frame(target_set_fraction_list)
  
  target_sets_fraction_df$Rank <- 1:nrow(master_df)
  
  return(target_sets_fraction_df)
}


### CURATED TARGETS ###
# TP53 target genes from Fischer et al., 2017: https://pubmed.ncbi.nlm.nih.gov/28288132/
tp53_curated_targets <- read.csv(paste0(PATH, "Validation_Files/TP53_Targets_Fischer_2017_Oncogene.csv"),
                                 header = T, check.names = F)[,1]

### DoRothEA  ###
# https://bioconductor.org/packages/release/data/experiment/vignettes/dorothea/inst/doc/dorothea.html
dorothea_net <- dorothea::dorothea_hs

tp53_dorothea_targets <- unlist(dorothea_net[dorothea_net$tf == "TP53", 'target'])

# Combine into list
tp53_select_target_set_list <- list("Curated_Fischer_2017" = tp53_curated_targets,
                                    "DoRothEA" = tp53_dorothea_targets)

# Call function
plot_combined_enrichment(pc_allGenes[pc_allGenes$R_i.name == "TP53",],
                         tp53_select_target_set_list, "TP53", 500, 50, NA)


############################################################
### PART E: PER-DRIVER NETWORK OVERLAY
############################################################  
# Useful tutorial: https://mr.schochastics.net/material/netVizR/

#' Generate results graph display for a driver of interest 
#' This graphical display takes confidence interaction relationships from 
#' a network like STRING or HumanBase, and displays these relationships
#' as a graph. The graph will only include hits with a q-value below the
#' given threshold. The edges and genes are colored by whether activation or 
#' repression is predicted by Dyscovr when the driver is mutated. 
#' @param master_df an output DF from Dyscovr
#' @param goi the name of the driver gene of interest
#' @param qval_thres a q-value threshold, below which a target-driver
#' pairing is considered significant
#' @param network a table with the network information from either 
#' HumanBase or STRING
#' @param network_label either "STRING" or "HumanBase" to denote 
#' what type of network this is
#' @param conf_thres a confidence threshold for a connection between two 
#' genes in the given network
#' @param top_n if not NA, gives a value of the top N most significant targets 
#' to use
create_graphical_representation <- function(master_df, goi, qval_thres, 
                                            network, network_label, conf_thres, 
                                            top_n) {
  
  # Subset the master DF to the given goi and q-value threshold
  master_df_sub <- master_df[(master_df$R_i.name == goi) & 
                               (master_df$q.value < qval_thres),]
  if(!is.na(top_n)) {master_df_sub <- master_df_sub[1:top_n, ]}
  nodes <-  c(goi, unique(master_df_sub$T_k.name))
  
  network_model_targs <- NA
  
  # Subset the network to the significant targets, plus the GOI
  if(network_label == "STRING") {
    network_model_targs <- network[(network$node1_name %fin% nodes) & 
                                     (network$node2_name %fin% nodes),]
  } else if (network_label == "HumanBase") {
    network_model_targs <- network[(network$node1_name %fin% nodes) & 
                                     (network$node2_name %fin% nodes),]
    
  } else {
    print("Error in given network label. Only implemented for STRING and 
          HumanBase networks.")
  }
  
  # Remove cases where nodes 1 and 2 are the same, or rows where the values in 
  # the row are duplicated in a different order. Divide all confidence scores by 
  # 1000 for STRING
  network_model_targs <- subset(network_model_targs, node1_name != node2_name & 
                                  !duplicated(cbind(pmin(node1_name, node2_name), 
                                                    pmax(node1_name, node2_name))))
  if(network_label == "STRING") {
    network_model_targs$combined_score <- network_model_targs$combined_score / 1000
  }
  
  # For each of the targets, get its a) directionality (up- or down-regulated), 
  # b) confidence score in the given network, in relation to the driver 
  edge_rows <- lapply(1:nrow(network_model_targs), function(i) {
    conf <- network_model_targs[i, 'combined_score']   
    if(conf > conf_thres) {
      g1 <- network_model_targs[i, 'node1_name']
      g2 <- network_model_targs[i, 'node2_name']
      return(data.frame("src" = g1, "target" = g2, 
                        "conf" = as.numeric(conf)))
    } else {return (NA)}
  })
  edge_rows <- edge_rows[!is.na(edge_rows)]
  edge_table <- do.call(rbind, edge_rows)
  colnames(edge_table) <- c("src", "target", "conf")

  unique_nodes <- unique(c(edge_table$src, edge_table$target))
  node_table <- as.data.frame(unique_nodes)
  colnames(node_table) <- "nodes"
  node_table$dir <- as.factor(unlist(lapply(unique_nodes, function(n) {
    dir <- NA
    if(n == goi) {dir <- 0}
    else {
      dir <- ifelse(as.numeric(master_df_sub[master_df_sub$T_k.name == n, 
                                             'estimate'])[1] > 0, 1, -1)
    }
    return(dir)
  })))
  node_table$node_size <- unlist(lapply(node_table$dir, function(d) {
    if(d == 0) {return(2)}
    else {return(1)}
  }))
  print(node_table)

  # Convert to ggraph format
  network_table <- tbl_graph(nodes = node_table, edges = edge_table, directed = F)
  
  # Create ggraph
  driver_targ_network <- ggraph(network_table, layout = "fr") +   # other options: stress
    geom_edge_link0(edge_colour = "darkgray", edge_width = 1) +  
    geom_node_point(aes(fill = dir, size = node_size), shape = 21) +    
    scale_size(range = c(5,10), guide = 'none') + 
    geom_node_label(aes(label = nodes), repel = TRUE, show.legend = F, 
                    label.size = NA, label.padding = 0.05) +
    theme_graph() +
    scale_fill_manual(values = c("#0072B5FF", "lightgreen", "#E18727FF"), 
                      labels = c("Pred. Downregulation", "Driver", 
                                 "Pred. Upregulation"),
                      name = "Direction of Regulation") +
    theme(legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12), legend.position="bottom") 

  # Render the network
  show(driver_targ_network)
  
  return(node_table)
}

#' Helper function to adjust the STRING network file for use with visual rendering
#' @param string_nw an unadjusted STRING network file downloaded from STRING
#' https://string-db.org/
#' @param string_nw_info a helper file from STRING with additional gene annotations
#' @param all_genes_id_conv ID conversion file from BioMart
adjust_string_nw_files <- function(string_nw, string_nw_info, all_genes_id_conv) {

  colnames(string_nw) <- c("node1", "node2", "combined_score")
  string_nw$index <- 1:nrow(string_nw)
  
  # Add info for any missing IDs
  missing_ids <- setdiff(string_nw_info$string_protein_id, 
                         unique(c(string_nw$node1, string_nw$node2)))
  missing_info <- lapply(missing_ids, function(id) {
    new_id <- unlist(strsplit(id, ".", fixed = T))[2]
    print(new_id)
    if(new_id %fin% all_genes_id_conv$ensembl_peptide_id) {
      id_row <- all_genes_id_conv[all_genes_id_conv$ensembl_peptide_id == new_id,]
      print(id_row)
      name <- unique(unlist(id_row[,'external_gene_name']))
      return(data.frame("string_protein_id" = id, "preferred_name" = name, 
                        "protein_size" = NA, "annotation" = NA))
    } else {return(NA)}
  })
  missing_info <- missing_info[!is.na(missing_info)]
  missing_info_df <- do.call(rbind, missing_info)
  string_nw_info <- rbind(string_nw_info, missing_info_df)
  
  # Convert any ENSG IDs to external gene names
  string_nw_info$preferred_name <- unlist(lapply(string_nw_info$preferred_name, function(n) {
    new_n <- n
    if(grepl("ENSG", n, fixed = T)) {
      new_n <- paste(unique(all_genes_id_conv[all_genes_id_conv$ensembl_gene_id == n, 
                                              'external_gene_name']), collapse = ";")
    }
    return(new_n)
  }))
  
  missing_ids_from_run <- setdiff(pc_allGenes$T_k.name, 
                                  string_nw_info$preferred_name)
  missing_info_from_run <- lapply(missing_ids_from_run, function(name) {
    ensp <- unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == name, 
                                     'ensembl_peptide_id'])
    ensp <- paste0("9606.", ensp)
    return(data.frame("string_protein_id" = ensp, "preferred_name" = name, 
                      "protein_size" = NA, "annotation" = NA))
  })
  missing_info_from_run_df <- do.call(rbind, missing_info_from_run)
  missing_info_from_run_df <- missing_info_from_run_df[
    !(missing_info_from_run_df$string_protein_id == "9606."),]
  string_nw_info <- rbind(string_nw_info, missing_info_from_run_df)
  
  string_nw_pt1 <- merge(string_nw[, c('node1', 'combined_score', 'index')], 
                         string_nw_info[, c('string_protein_id', 'preferred_name')], 
                         by.x = 'node1', by.y = 'string_protein_id')
  colnames(string_nw_pt1) <- c("node1", "combined_score", "index", "node1_name")
  string_nw_pt1 <- string_nw_pt1[order(string_nw_pt1$index),]
  
  string_nw_pt2 <- merge(string_nw[, c('node2', 'combined_score', 'index')], 
                         string_nw_info[, c('string_protein_id', 'preferred_name')], 
                         by.x = 'node2', by.y = 'string_protein_id')
  colnames(string_nw_pt2) <- c("node2", "combined_score", "index", "node2_name")
  string_nw_pt2 <- string_nw_pt2[order(string_nw_pt2$index),]
  
  string_nw_full <- merge(string_nw_pt1, string_nw_pt2, 
                          by = c('combined_score', 'index'))
  
  string_nw_full$node1 <- unlist(lapply(string_nw_full$node1, function(x) 
    unlist(strsplit(x, ".", fixed = T))[2]))
  string_nw_full$node2 <- unlist(lapply(string_nw_full$node2, function(x) 
    unlist(strsplit(x, ".", fixed = T))[2]))
  
  return(string_nw_full)
}

# Import the String network files and adjust (only need to do this once)
string_nw <- fread(paste0(PATH, "Validation_Files/string.9606.protein.links.v11.5.txt"),
                   header = T)
string_nw_info <- read.table(paste0(PATH, "Validation_Files/string.9606.protein.info.v11.5.txt"), 
                             header = T, sep = "\t")
string_nw_full <- adjust_string_nw_files(string_nw, string_nw_info, 
                                         all_genes_id_conv)
fwrite(string_nw_full, paste0(
  PATH, "Validation_Files/string.9606.protein.links.v11.5.namesAdded.txt"))

# Read back, once completed
string_nw_full <- fread(paste0(PATH, "Validation_Files/string.9606.protein.links.v11.5.namesAdded.txt"))

# Call function
create_graphical_representation(pc_allGenes[pc_allGenes$R_i.name == "TP53",], 
                                "TP53", 0.01, string_nw_full, "STRING", 0.4, NA)


############################################################
### PART F: PER-DRIVER VOLCANO PLOTS
############################################################ 
#' Create a volcano plot for a given driver in a particular cancer type or 
#' subtype, with significant hits highlighted and labeled
#' @param master_df a data frame produced by Dyscovr
#' @param goi a driver gene-of-interest uniprot ID
#' @param qval_thres a q-value threshold for significance
create_volcano_plot <- function(master_df, goi, qval_thres) {
  
  # Subset to the GOI 
  master_df <- master_df[master_df$R_i == goi,]
  master_df <- na.omit(master_df)
  
  # Get the log2(Beta), which we are considering here to be like log2(fold change)
  log2_beta <- unlist(lapply(master_df$estimate, function(e) { 
    if((is.nan(e)) | (length(e) == 0)) {return(0)}
    else {return(e)}
  }))

  # Get the -log10(qvalue)
  neg_log10_qval <- unlist(lapply(master_df$q.value, function(q) { 
    if((is.nan(q)) | (length(q) == 0)) {return(-log10(1))}
    else {return(-log10(q))}
  }))

  up_or_down <- unlist(lapply(1:nrow(master_df), function(i) {
    if(master_df[i, 'q.value'] < qval_thres) {
      beta_sign <- ifelse(master_df[i, 'estimate'] > 0, 1, 0) 
      if(beta_sign == 1) {return("up")}
      else {return("down")}
    } else {return("ns")}
  }))
  
  gene_names <- unlist(lapply(1:nrow(master_df), function(i) 
    ifelse(-log10(master_df$q.value[i]) > 5, master_df$T_k.name[i], NA)))

  # Put this into a DF for plotting
  plot_df <- data.frame("log2_beta" = log2_beta, "neg_log10_qval" = neg_log10_qval,
                        "up_or_down" = up_or_down, "gene" = gene_names)
  plot_df[which(plot_df$neg_log10_qval == Inf), 'neg_log10_qval'] <- max(
    plot_df$neg_log10_qval[plot_df$neg_log10_qval != Inf]) + 5
  
  xlim_min <- -1
  xlim_max <- 1
  
  # Optionally squish each point outside the given x-axis limits
  #plot_df[which(plot_df$log2_beta > xlim_max), 'log2_beta'] <- xlim_max
  #plot_df[which(plot_df$log2_beta < xlim_min), 'log2_beta'] <- xlim_min
  
  # Set color and alpha params
  cols <- c("up" = "#E18727FF", "down" = "#0072B5FF", "ns" = "grey") 
  alphas <- c("up" = 0.95, "down" = 0.95, "ns" = 0.5)
  
  p <- ggplot(plot_df, aes(x = log2_beta, y = neg_log10_qval, fill = up_or_down, 
                           alpha = up_or_down, label = gene)) + 
    geom_point(shape = 21, size = 2.5) + xlab("Mutation Coefficient") + 
    ylab("-log10(q-value)") + theme_minimal() + 
    ylim(0, max(plot_df$neg_log10_qval[plot_df$neg_log10_qval != Inf]) + 1) +
    geom_hline(yintercept = -log10(qval_thres), linetype = "dashed") +
    scale_fill_manual(values = cols) + scale_alpha_manual(values = alphas) + 
    theme(legend.position = "none", 
          axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 12)) + 
    geom_text_repel(min.segment.length = unit(0.05, "lines"), size = 4, 
                    force = 3.5, box.padding = 0.6, colour = "black") + 
    xlim(xlim_min, xlim_max)
  p
}

create_volcano_plot(pc_allGenes[pc_allGenes$R_i.name == "TP53",], "P04637", 0.01)
create_volcano_plot(pc_allGenes[pc_allGenes$R_i.name == "PIK3CA",], "P42336", 0.01)
create_volcano_plot(pc_allGenes[pc_allGenes$R_i.name == "KRAS",], "P01116", 0.01)
create_volcano_plot(pc_allGenes[pc_allGenes$R_i.name == "IDH1",], "O75874", 0.01)


############################################################
### SUPPLEMENTAL TABLE 2: TP53 ENRICHMENTS IN OTHER EXTERNAL DATASETS
############################################################ 
### TRRUST ###
# https://www.grnpedia.org/trrust/
trrust_df <- read.table(paste0(PATH, "Validation_Files/trrust_rawdata.human.txt"),
                        header = F, check.names = F)
trrust_df_tp53 <- trrust_df[(trrust_df[,1] == "TP53") | (trrust_df[,2] == "TP53"),]
tp53_trrust_targets <- unique(c(trrust_df_tp53$V1, trrust_df_tp53$V2))
tp53_trrust_targets <- tp53_trrust_targets[tp53_trrust_targets != "TP53"]

compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == "TP53"], 
                               tp53_trrust_targets, "k-s", 0.01)

### hTFtarget ###
# http://bioinfo.life.hust.edu.cn/hTFtarget
tf_target_df <- read.table(paste0(PATH, "Validation_Files/TF-Target-information.txt"),
                           header = T, check.names = F, sep = "\t")
tf_target_df_tp53 <- tf_target_df[grepl("P53", tf_target_df$TF),]
tf_target_tp53 <- unique(tf_target_df_tp53$target)

compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == "TP53"], 
                               tf_target_tp53, "k-s", 0.01)

### Reactome ###
# https://reactome.org/PathwayBrowser/
# Transcriptional Regulation by TP53 (R-HSA-3700989)
tp53_transcriptional_reg_pw <- read.csv(paste0(
  PATH, "Validation_Files/Transcriptional Regiulation by TP53_R-HSA-3700989.csv"), 
                                        header = T, check.names = F)
tp53_transcriptional_reg_pw_gns <- unlist(lapply(tp53_transcriptional_reg_pw$MoleculeName, function(x) 
  unlist(strsplit(x, " ", fixed = TRUE))[2]))
tp53_transcriptional_reg_pw_gns <- unique(tp53_transcriptional_reg_pw_gns[
  tp53_transcriptional_reg_pw_gns != "TP53"])

compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == "TP53"], 
                               tp53_transcriptional_reg_pw_gns, "k-s", 0.01)

# TP53 Regulates Metabolic Genes (R-CFA-5628897)
tp53_regulates_metabolism_pw <- read.csv(paste0(
  nw_path, "Validation_Files/TP53 Regulates Metabolic Genes_R-HSA-5628897.csv"), 
                                         header = TRUE, check.names = FALSE)
tp53_regulates_metabolism_pw_gns <- unlist(lapply(tp53_regulates_metabolism_pw$MoleculeName, function(x) 
  unlist(strsplit(x, " ", fixed = TRUE))[2]))
tp53_regulates_metabolism_pw_gns <- unique(tp53_regulates_metabolism_pw_gns[tp53_regulates_metabolism_pw_gns != "TP53"])

compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == "TP53"], 
                               tp53_regulates_metabolism_pw_gns, "k-s", 0.01)

### KEGG ###
# https://www.gsea-msigdb.org/gsea/msigdb/cards/KEGG_P53_SIGNALING_PATHWAY
# KEGG P53 Signaling Pathway (hsa04115)
tp53_kegg_pathway <- read.table(paste0(
  PATH, "Validation_Files/KEGG_P53_SIGNALING_PATHWAY.v7.5.1.txt"), sep = "\t", header = T)
tp53_kegg_pathway_genes <- unlist(strsplit(tp53_kegg_pathway$KEGG_P53_SIGNALING_PATHWAY[19], 
                                           ",", fixed = T))

compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == "TP53"], 
                               tp53_kegg_pathway_genes, "k-s", 0.01)

### HUMAN BASE ###
# Link to download: https://hb.flatironinstitute.org/download
global_nw <- fread(paste0(input_path, "HumanBase/global_top/global_top"), header = F)
colnames(global_nw) <- c("entrez1", "entrez2", "posterior_prob")

# Limit to posterior probability of at least 0.5 (default on website)
global_nw_0.5 <- global_nw[global_nw$posterior_prob > 0.5, ]

# Change Entrez IDs to gene names
mapping <- as.data.frame(bitr(unique(c(global_nw_0.5$entrez1, global_nw_0.5$entrez2)), 
                              fromType = "ENTREZID", toType = "SYMBOL", 
                              OrgDb=org.Hs.eg.db, drop = T))
global_nw_0.5$gene_name1 <- unlist(lapply(global_nw_0.5$entrez1, function(x) {
  symb <- mapping[mapping$ENTREZID == x, 'SYMBOL']
  if(length(symb) == 0) {return(NA)}
  else{return(symb)}
}))
global_nw_0.5$gene_name2 <- unlist(lapply(global_nw_0.5$entrez2, function(x) {
  symb <- mapping[mapping$ENTREZID == x, 'SYMBOL']
  if(length(symb) == 0) {return(NA)}
  else{return(symb)}
}))

fwrite(global_nw_0.5, paste0(input_path, "HumanBase/global_top/global_top_0.5_idconv.csv"))

# Get driver interactors from HumanBase 
tp53_hb_interactors <- unique(unlist(global_nw_0.5[(global_nw_0.5$gene_name1 == "TP53") |
                                                  (global_nw_0.5$gene_name2 == "TP53"), 
                                                c("gene_name1", "gene_name2")]))  #262
compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == "TP53"], 
                               tp53_hb_interactors, "k-s", 0.01)

### INTACT ###
# Website: https://www.ebi.ac.uk/intact/home
# Search TP53 and download set of interactors
tp53_intact_interactors <- read.table(paste0(PATH, "Validation_Files/IntAct/tp53_interactors.txt"),
                                      header = T, sep = "\t")
tp53_intact_interactors[,1] <- unlist(lapply(tp53_intact_interactors[,1], function(x) 
  unlist(strsplit(x, ":", fixed = T))[2]))
tp53_intact_interactors[,2] <- unlist(lapply(tp53_intact_interactors[,2], function(x) 
  unlist(strsplit(x, ":", fixed = T))[2]))
tp53_intact_interactors <- unique(unlist(tp53_intact_interactors[(tp53_intact_interactors[,1] == "P04637") |
                                                     (tp53_intact_interactors[,2] == "P04637"), c(1,2)]))
tp53_intact_interactors <- tp53_intact_interactors[tp53_intact_interactors != "P04637"]

# Convert to gene names, or modify enrichment function
tp53_intact_interactors <- unlist(lapply(tp53_intact_interactors, function(ia) 
  unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == ia, 'external_gene_name'])))

compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == "TP53"], 
                               tp53_intact_interactors, "k-s", 0.01)  

### BIOGRID ###
# Website: https://thebiogrid.org/
# Search TP53 and download set of interactors

tp53_biogrid_interactors <- read.table(paste0(PATH, "Validation_Files/BioGrid/TP53/BIOGRID-GENE-113010-4.4.233.tab3.txt"),
                                      header = T, sep = "\t")
tp53_biogrid_interactors <- unique(unlist(c(tp53_biogrid_interactors$`Systematic Name Interactor A`,
                                            tp53_biogrid_interactors$`Systematic Name Interactor B`)))
tp53_biogrid_interactors <- tp53_biogrid_interactors[tp53_biogrid_interactors != "TP53"]
compute_statistical_enrichment(pc_allGenes[pc_allGenes$R_i.name == "TP53"], 
                               tp53_biogrid_interactors, "k-s", 0.01)  

