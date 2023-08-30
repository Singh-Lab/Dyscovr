############################################################
### CREATE ENRICHMENT VISUALIZATIONS
### Written By: Sara Geraghty, July 2022
############################################################

library(grid)
#library(ggpubr)
library("STRINGdb")
library(igraph)
library(EnsDb.Hsapiens.v86)
library("dorothea")
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
#### COMPUTE STATISTICAL ENRICHMENT USING A HYPERGEOMETRIC
#### (ONE-SIDED FISHER'S EXACT TEST) AND/OR A K-S TEST
############################################################
############################################################
#' Use a one-sided Fisher's exact test/ Hypergeometric test or a Kolmogorov-Smirnov
#' test to test for statistical enrichment
#' @param master_df a master DF produced from run_linear_model() that has q-values
#' @param known_targs a vector of known targets for a gene of interest (master
#' DF should already be subsetted to just this gene); ex. ChIP-eat targets, HumanBase
#' @param type the type of test, either "fisher's exact" or "k-s"
#' @param qval_thres q-value threshold for HG enrichment
compute_statistical_enrichment <- function(master_df, known_targs, type, qval_thres) {
  final_res <- NA
  
  if ((type == "fisher's exact") | (type == "both")) {
    sig_targets <- unique(master_df[master_df$q.value < qval_thres, 'T_k.name'])
    nonsig_targets <- unique(master_df[master_df$q.value > qval_thres, 'T_k.name'])
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
    
    contigency_table <- matrix(c(num_known_sig_targets, num_known_nonsig_targets,
                                 num_notknown_sig_targets, num_notknown_nonsig_targets),
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
    ranks_of_targs <- ranks_of_targs[!(is.infinite(ranks_of_targs) | (is.na(ranks_of_targs)))]
    #print(TRUE %in% duplicated(ranks_of_targs))
    
    # Perform the K-S test, with the uniform distribution
    ks_res <- ks.test(ranks_of_targs, uniform, alternative = "greater")
    
    #print(ks_res)
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

# Call function
compute_statistical_enrichment(allgenes_p53, tp53_nw_targs, "both", 0.1)


input <- melt(as.data.frame(list("TP53" = -log10(1.571079e-06), "PIK3CA" = -log10(0.00031161), 
                                 "KRAS" = -log10(0.009104977), "IDH1" = -log10(0.7308022))))
colnames(input) <- c("Driver", "negLog10pval")

ggplot(input, aes(y = negLog10pval, fill = Driver, x = reorder(Driver, negLog10pval, mean))) + 
  geom_bar(position = "dodge", width = 0.95, stat = "identity", show.legend = FALSE, color = "black") + 
  scale_fill_manual(values = c("#0072B5FF", "#BC3C29FF", "#20854EFF", "#FFDC91FF")) + #scale_color_nejm() +
  xlab("Driver") + ylab(paste0("-log10(pval) of HG Enrichment,\n Top 100 STRING Neighbors (q<", paste0(0.001, ")\n"))) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black") + coord_flip() + #theme_minimal() +
  scale_y_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text(face="bold", size = 14), 
        axis.title=element_text(size=16, face="bold"), panel.grid.major = element_blank(),
        panel.background = element_rect(fill = 'white'))

#' Given a set of STRING/ HumanBase confidence values and p-values values from the model,
#' calculate a spearman correlation between the confidence values and 1-(p.value)
#' @param master_df the master DF
#' @param nw_table a string/humanbase table for the gene with targets and confidence values 0-1
#' @param goi the name of the regulatory gene-of-interest
compute_network_spearman <- function(master_df, nw_table, goi) {
  
  # Convert the p-values to 1-p.value
  conv_pvals <- data.frame("gene.name" = master_df$T_k.name, "conv.pval" = 1 - (master_df$p.value))
  nw_table <- nw_table[nw_table$node1 == goi,]
  string_confvals <- data.frame("gene.name" = nw_table$node2, 
                                "conf" = nw_table$combined_score)
  comb_df <- merge(conv_pvals, string_confvals, by = "gene.name", all = FALSE)
  print(dim(comb_df))
  
  spearman <- cor.test(comb_df$conv.pval, comb_df$conf, method = "spearman", use = "pairwise")
  print(spearman)
  return(spearman$p.value)
  
}

compute_network_spearman(allgenes_p53_tcga, tp53_neighbors_conf500, "TP53")


############################################################
############################################################
#### PLOT COMBINED ENRICHMENT ACROSS MULTIPLE GENE SET SOURCES
############################################################
############################################################
#' Function to plot the enrichment of various target groups among the significant gene hits,
#' together on one plot, where the enrichment in the given source is denoted by color and line form
#' @param master_df a master DF produced from run_linear_model() that has q-values
#' @param target_sets_list a list of each of the target sets we are interested in doing enrichment on,
#' with each set's name/ label as the name for the list item
#' @param goi a string denoting the gene-of-interest, for labeling graph
#' @param thres a threshold for the number of hit genes we show in the full plot; in NA then we plot
#' enrichment across all gene targets
#' @param cutout_thres a threshold for the inset/cutout enrichment plot ('zoomed' in on)
#' @param silent_df OPT: an additional, corresponding DF with silent mutation results
plot_combined_enrichment <- function(master_df, target_sets_list, goi, thres, cutout_thres, silent_df) {
  
  # For each source, create a 0 and 1 vector for each  target hit to indicate if it 
  # is a a hit in that source
  target_set_vector_list <- lapply(target_sets_list, function(source) {
    source <- setdiff(source, setdiff(source, master_df$T_k.name))
    vect <- unlist(lapply(1:nrow(master_df), function(i) {
      x <- master_df[i, 'T_k.name']
      return(ifelse(x %fin% source, 1, 0))
    }))
    return(vect)
  })
  
  target_set_vector_list_silent <- NA
  if(!is.na(silent_df)) {
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
  
  #print(target_set_vector_list)
  #print(target_set_vector_list_silent)
  
  # Use helper function to get fraction of genes at given rank
  target_set_fraction_df <- get_frac_of_genes_at_rank(target_set_vector_list, master_df)
  target_set_fraction_df_silent <- NA
  if(!is.na(silent_df)) {
    target_set_fraction_df_silent <- get_frac_of_genes_at_rank(target_set_vector_list_silent,
                                                               silent_df)
    target_set_fraction_df <- merge(target_set_fraction_df, target_set_fraction_df_silent, by = "Rank")
  }
  
  
  target_sets_fraction_df_m <- melt(target_set_fraction_df, "Rank")
  colnames(target_sets_fraction_df_m) <- c("Rank", "Source", "Frac")
  target_sets_fraction_df_m$Mutation.Type <- unlist(lapply(target_sets_fraction_df_m$Source, function(x) {
    val <- "Nonsynonymous"
    if(grepl("S.", x, fixed = TRUE)) {val <- "Synonymous"}
    return(val)
  }))
  
  print(head(target_sets_fraction_df_m))
  print(unique(target_sets_fraction_df_m$Source))
  
  # Create inset plot
  p2 <- ggplot(target_sets_fraction_df_m[target_sets_fraction_df_m$Rank %in% 1:cutout_thres,], 
               mapping = aes(x = Rank, y = Frac, color = Source)) + 
    geom_line(aes(linetype=Mutation.Type), size = 1.25) +
    #geom_point(size = 1) + 
    scale_x_continuous(limits = c(1,cutout_thres)) + 
    #scale_color_nejm() + 
    #scale_alpha_manual(values=c(0.9,0.9,0.4,0.4)) +
    #scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#BC3C29FF", "#0072B5FF")) + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.minor = element_blank(), axis.text.x = element_text(size = 14, face = "bold"), 
          axis.text.y =  element_text(size = 14, face = "bold"))
  
  p <- ggplot(target_sets_fraction_df_m[target_sets_fraction_df_m$Rank %in% 1:thres,], 
              mapping = aes(x = Rank, y = Frac, color = Source)) + 
    geom_line(aes(linetype=Mutation.Type), size = 1.25) +
    #theme_minimal() +
    theme(axis.text = element_text(size = 16, face = "bold"), 
          axis.title = element_text(size = 18, face = "bold"), #panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), #panel.border = element_blank(),
          panel.background = element_rect(fill = 'white'), legend.text=element_text(size=14), 
          legend.title = element_text(size = 16)) +
    #geom_point(size = 1) + 
    #ggtitle(paste("Fraction of Model-Prioritized Target Genes in", paste(goi, "Sourced Targets"))) +
    xlab("Gene Rank") + ylab("Fraction in Set of Known Targets") + 
    scale_x_continuous(limits = c(1,thres)) +
    #scale_color_nejm() + 
    #scale_alpha_manual(values=c(0.9,0.9,0.4,0.4)) +
    #scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#BC3C29FF", "#0072B5FF"),
                       #labels = c(names(target_sets_list), rep("", times = length(names(target_sets_list))))) + 
    annotation_custom(ggplotGrob(p2), xmin =(thres-round(thres/1.5)), xmax = thres, ymin = 0.30, ymax = 1.0) +
    geom_rect(data=target_sets_fraction_df_m[1:cutout_thres,], 
              aes(xmin = 1, xmax = cutout_thres, ymin = 0, ymax = 1, fill = "gray"),
              color = NA, fill = alpha("gray", .01))
  
  #vp <- viewport(width = 0.4, height = 0.4, x = 0.8, y = 0.2)
  #print(p)
  #print(p, vp = vp)
  
  return(p)
}


#' Get the fraction of genes at or above the given rank, for each set
#' @param target_set_vector_list a named list of each source with hit/not binary vector
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

# Create a list of all the gene sets we want to look for enrichment in
tp53_target_set_list <- list("ChIP-eat" = tp53_chipeat_targets, "Compiled.Cancer" = known_cancer_genes,
                             "STRING" = tp53_string_nw_targs, "HumanBase" = tp53_nw_targs, 
                             "Curated.Fischer.2017" = tp53_curated_targets)
tp53_target_set_list <- list("Compiled.Cancer" = known_cancer_genes,
                             "STRING" = tp53_string_nw_targs, "HumanBase" = tp53_nw_targs, 
                             "Curated.Fischer.2017" = tp53_curated_targets, 
                             "Curated.TRRUST" = tp53_trrust_targets, "Curated.TF-Target" = tf_target_tp53_full)

pik3ca_target_set_list <- list("Compiled.Cancer" = known_cancer_genes,
                               "STRING" = pik3ca_string_nw_targs, "HumanBase" = pik3ca_nw_targs,
                               "Curated.Cizkova.2017" = pik3ca_curated_targs)

# For grant
tp53_select_target_set_list <- list("Curated.Fischer.2017" = tp53_curated_targets, "KEGG.hsa04115" = tp53_kegg_pathway_genes,
                                    "TRRUST" = tp53_trrust_targets, "DoRothEA" = tp53_dorothea_targets)
tp53_select_target_set_list <- list("Curated.Fischer.2017" = tp53_curated_targets,
                                    "DoRothEA" = tp53_dorothea_targets)
tp53_select_target_set_list_cna <- list("STRING" = tp53_string_nw_targs,
                                        "TRRUST" = tp53_trrust_targets_upstr)
pik3ca_select_target_set_list <- list("Curated.Cizkova.2017" = pik3ca_curated_targs, "KEGG.hsa04151" = pik3ca_kegg_pathway_genes)


# Call function
plot_combined_enrichment(allgenes_p53, tp53_target_set_list, "TP53", 500, 50, NA)
plot_combined_enrichment(allgenes_pik3ca, pik3ca_target_set_list, "PIK3CA", 500, 50, NA)



#' Function to plot the enrichment of a target gene set among the significant gene hits from multiple R_is,
#' together on one plot, where the enrichment of targets from the given GOI is denoted by color and line form
#' @param list_of_master_dfs a named list of master DF produced from run_linear_model(), each of which
#' has gene names and q-values. Names are the name of the GOI for that master DF
#' @param target_gene_set a vector of targets of interest
#' @param thres a threshold for the number of hit genes we show in the full plot; in NA then we plot
#' enrichment across all gene targets
#' @param cutout_thres a threshold for the inset/cutout enrichment plot ('zoomed' in on)
#' @param name_of_target_set a name for the target set of genes for labeling graph
plot_combined_enrichment_mult_genes <- function(list_of_master_dfs, target_gene_set, 
                                                thres, cutout_thres, name_of_target_set) {
  
  # For each GOI, create a 0 and 1 vector for each  target hit to indicate if it 
  # is a a hit in the given source
  goi_vector_list <- lapply(list_of_master_dfs, function(goi_master) {
    vect <- unlist(lapply(1:nrow(goi_master), function(i) {
      x <- goi_master[i, 'T_k.name']
      return(ifelse(x %fin% target_gene_set, 1, 0))
    }))
    return(vect)
  })
  
  #' Get the fraction of genes at or above the given rank, for each GOI
  goi_fraction_list <- lapply(1:length(goi_vector_list), function(list_index) {
    goi_vector <- goi_vector_list[[list_index]]
    frac_vect <- c()
    count_of_match_genes <- 0
    for (i in 1:length(goi_vector)) {
      val_of_curr_tg <- goi_vector[i]
      count_of_match_genes <- count_of_match_genes + val_of_curr_tg
      frac <- count_of_match_genes / i
      frac_vect <- c(frac_vect, frac)
    }
    frac_vect_df <- as.data.frame(frac_vect)
    colnames(frac_vect_df) <- names(list_of_master_dfs)[list_index]
    frac_vect_df$Rank <- 1:length(goi_vector)
    return(frac_vect_df)
  })
  min_df_nrow <- min(unlist(lapply(goi_fraction_list, function(x) return(as.numeric(nrow(x))))))
  goi_fraction_list <- lapply(goi_fraction_list, function(df) return(df[1:min_df_nrow,]))
  print(head(goi_fraction_list))
  goi_fraction_df <- Reduce(function(df1, df2) merge(df1, df2, by = "Rank", all = FALSE), goi_fraction_list)
  #goi_fraction_df <- merge(goi_fraction_list, by = "Rank")
  print(head(goi_fraction_df))
  
  goi_fraction_df_m <- melt(goi_fraction_df, "Rank")
  colnames(goi_fraction_df_m) <- c("Rank", "GOI", "Frac")
  print(head(goi_fraction_df_m))
  print(unique(goi_fraction_df_m$GOI))
  
  # Create inset plot
  p2 <- ggplot(goi_fraction_df_m[goi_fraction_df_m$Rank %in% 1:cutout_thres,], 
               mapping = aes(x = Rank, y = Frac, color = GOI)) + geom_line(aes(linetype=GOI), size = 1.25) +
    #geom_point(size = 1) + 
    scale_color_nejm() + 
    scale_x_continuous(limits = c(1,cutout_thres)) +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.minor = element_blank())
  
  p <- ggplot(goi_fraction_df_m[goi_fraction_df_m$Rank %in% 1:thres,], 
              mapping = aes(x = Rank, y = Frac, color = GOI)) + geom_line(aes(linetype=GOI), size = 1.25) +
    #geom_point(size = 1) + 
    scale_color_nejm() + 
    ggtitle(paste("Fraction of Model-Prioritized Target Genes in", name_of_target_set)) +
    xlab("Gene Rank") + ylab("Fraction of Top Hit Genes in Source Gene Set") + 
    scale_x_continuous(limits = c(1,thres)) +
    theme(axis.text.x =  element_text(size = 12), axis.text.y =  element_text(size = 12),
          axis.title = element_text(size = 13, face = "bold"), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.border = element_blank()) + theme_minimal() +
    annotation_custom(ggplotGrob(p2), xmin =(thres-round(thres/1.5)), xmax = thres, ymin = 0.30, ymax = 1.0) +
    geom_rect(data=goi_fraction_df_m[1:cutout_thres,], 
              aes(xmin = 1, xmax = cutout_thres, ymin = 0, ymax = 1, fill = "gray"),
              color = NA, fill = alpha("gray", .01))
  
  return(p)
  
}

# Call function
list_of_master_dfs <- list("TP53" = allgenes_p53_tcga, "PIK3CA" = allgenes_pik3ca_tcga, 
                           "TTN" = allgenes_ttn_tcga)
plot_combined_enrichment_mult_genes(list_of_master_dfs, breast_biomarkers, 500, 50, "Breast Biomarkers")


############################################################
############################################################
#### PLOT BAR PLOTS OF INDIVIDUAL ENRICHMENT IN NETWORKS,  
#### CANCER GENES, AND ChIP-EAT TARGETS
############################################################
############################################################

#' Function to plot the enrichment of network targets among the significant target genes
#' @param master_df a master DF produced from run_linear_model() that has q-values
#' @param network_targs a vector of gene targets from HumanBase/ STRING for the given protein
#' of interest in the given tissue
#' @param goi a string denoting the gene-of-interest, for labeling graph
#' @param bin_size if binning, the bin size (0 indicates no binning)
plot_network_enrichment <- function(master_df, network_targs, goi, bin_size) {
  # Combine all the target genes above the given threshold into one vector
  target_hits <- unlist(master_df$T_k.name)
  print(length(target_hits))
  
  # Fill in a vector of 0 and 1 for each significant target hit to indicate if it 
  # is a HumanBase
  nw_vect <- unlist(lapply(target_hits, function(x) {
    return(ifelse(x %fin% network_targs, 1, 0))
  }))
  
  # Plot the enrichment (fraction of network genes at/ above given rank)
  ranks <- 1:length(target_hits)
  frac_nw_vect <- c()
  count_of_nw_genes <- 0
  for (i in 1:length(target_hits)) {
    val_of_curr_tg <- nw_vect[i]
    count_of_nw_genes <- count_of_nw_genes + val_of_curr_tg
    frac <- count_of_nw_genes / i
    frac_nw_vect <- c(frac_nw_vect, frac)
  }
  
  if(bin_size == 0) {
    df <- data.frame("Rank" = as.factor(ranks[1:100]), "Frac.of.Targets" = frac_nw_vect[1:100])
    p <- ggplot(df, mapping = aes(x = Rank, y = Frac.of.Targets)) + geom_point(fill="#69b3a2", size = 2) + 
      #ggtitle(paste("Enrichment of Significant Target Genes in", paste(goi, "Network Targets"))) +
      xlab("Gene Rank") + ylab("Enrichment of Network Genes") + 
      scale_x_discrete(breaks = c(1, 25, 50, 75, 100), labels = c("1", "25", "50", "75", "100")) +
      theme(axis.text.x =  element_text(size = 12), axis.text.y =  element_text(size = 12),
            axis.title = element_text(size = 13, face = "bold"))
    
  } else {
    spl_vals <- split(frac_nw_vect, ceiling(seq_along(frac_nw_vect) / bin_size))
    head(spl_vals)
    
    spl_vals_df <- data.frame("vals" = frac_nw_vect)
    grp_labs <- rep(1:length(spl_vals), each = bin_size)
    num_grps <- length(unique(grp_labs))
    while(length(grp_labs) > nrow(spl_vals_df)) {grp_labs <- grp_labs[1:length(grp_labs)-1]}
    spl_vals_df$grp <- as.integer(grp_labs)
    print(head(spl_vals_df))
    
    
    spl_vals_means <- data.frame("grp" = as.integer(1:num_grps), 
                                 "mean" = as.numeric(lapply(1:num_grps, function(i) 
                                   mean(as.numeric(spl_vals_df[spl_vals_df$grp == i, 'vals'])))))
    print(head(spl_vals_means))
    
    p <- ggplot(spl_vals_df[spl_vals_df$grp %in% 1:15,], mapping = aes(x = as.factor(grp), y = vals)) + geom_boxplot(fill="#69b3a2") + 
      geom_point(data = spl_vals_means[spl_vals_means$grp %in% 1:15,], mapping = aes(x = grp, y = mean), color="red") +
      geom_line(data = spl_vals_means[spl_vals_means$grp %in% 1:15,], mapping = aes(x = grp, y = mean, group=1)) +
      xlab(paste("Group #, Bin Size:", bin_size)) + ylab(paste("Fraction of", paste(goi, "Network Targets"))) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"))
    
  }
  return(p)
}

# Call function
plot_network_enrichment(allgene_sixMostMut_p53, tp53_string_nw_targs, "TP53", 10)

# Create a stacked plot with both enrichments on the same plot
tp53_hb_plt <- plot_network_enrichment(allgene_sixMostMut_p53, tp53_nw_targs, "TP53", 10)
tp53_string_plt <- plot_network_enrichment(allgene_sixMostMut_p53, tp53_string_nw_targs, "TP53", 10)
ggarrange(tp53_hb_plt, tp53_string_plt, labels = c("A", "B"), ncol = 1, nrow = 2)


#' Function to plot the enrichment of ChIP-eat targets among the significant target genes
#' @param master_df a master DF produced from run_linear_model() that has q-values
#' @param chipeat_targs a vector of known ChIP-eat targets for a gene of interest (master
#' DF should already be subsetted to just this gene)
#' @param goi a string denoting the gene-of-interest, for labeling graph
#' @param bin_size if binning, the bin size (0 indicates no binning)
plot_chipseq_enrichment <- function(master_df, chipeat_targs, goi, bin_size) {
  # Get all the target genes 
  target_hits <- unlist(master_df$T_k.name)
  print(length(target_hits))
  print(length(chipeat_targs))
  
  # Fill in a vector of 0 and 1 for each significant target hit to indicate if it 
  # is a ChIP-eat target
  chipeat_vect <- unlist(lapply(target_hits, function(x) {
    return(ifelse(x %fin% chipeat_targs, 1, 0))
  }))
  
  # Plot the enrichment (fraction of ChIP-eat genes at/ above given rank)
  ranks <- 1:length(target_hits)
  frac_chipeat_vect <- c()
  count_of_chipeat_genes <- 0
  for (i in 1:length(target_hits)) {
    val_of_curr_tg <- chipeat_vect[i]
    count_of_chipeat_genes <- count_of_chipeat_genes + val_of_curr_tg
    frac <- count_of_chipeat_genes / i
    frac_chipeat_vect <- c(frac_chipeat_vect, frac)
  }
  
  if(bin_size == 0) {
    plot(ranks, frac_chipeat_vect, pch = 16, 
         main = paste("Enrichment of Significant Target Genes in", paste(goi, "ChIP-eat Targets")),
         xlab = "Rank", ylab = "Fraction of Target Genes that are ChIP-eat Targets")
  } else {
    #num_bins <- ceiling(length(ranks) / bin_size)
    spl_vals <- split(frac_chipeat_vect, ceiling(seq_along(frac_chipeat_vect) / bin_size))
    head(spl_vals)
    
    spl_vals_df <- data.frame("vals" = frac_chipeat_vect)
    grp_labs <- rep(1:length(spl_vals), each = bin_size)
    num_grps <- length(unique(grp_labs))
    while(length(grp_labs) > nrow(spl_vals_df)) {grp_labs <- grp_labs[1:length(grp_labs)-1]}
    spl_vals_df$grp <- as.integer(grp_labs)
    print(head(spl_vals_df))
    
    
    spl_vals_means <- data.frame("grp" = as.integer(1:num_grps), 
                                 "mean" = as.numeric(lapply(1:num_grps, function(i) 
                                   mean(as.numeric(spl_vals_df[spl_vals_df$grp == i, 'vals'])))))
    print(head(spl_vals_means))
    
    ggplot(spl_vals_df[spl_vals_df$grp %in% 1:75,], mapping = aes(x = as.factor(grp), y = vals)) + geom_boxplot(fill="#69b3a2") + 
      geom_point(data = spl_vals_means[spl_vals_means$grp %in% 1:75,], mapping = aes(x = grp, y = mean), color="red") +
      geom_line(data = spl_vals_means[spl_vals_means$grp %in% 1:75,], mapping = aes(x = grp, y = mean, group=1)) +
      xlab(paste("Group #, Bin Size:", bin_size)) + ylab("Fraction of Targets that are ChIP-eat Targets") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
}

#' Function to plot the enrichment of cancer genes among the significant target genes
#' @param master_df_sig a master DF produced from run_linear_model() that has q-values
#' and has been thresholded to only those pairings that exceed a significance threshold
#' @param known_cancer_genes_table a data frame containing a compiled list of known 
#' cancer genes (CGC, Vogelstein, etc.) with various types of IDs
#' @param tfcancer_df a data frame from TFcancer, which has lists of transcription
#' factors shown to be associated with cancer from curated literature
plot_cancer_enrichment <- function(master_df_sig, known_cancer_genes_table, tfcancer_df) {
  # Get all the significant target genes 
  significant_target_hits <- unlist(master_df_sig$T_k.name)
  print(significant_target_hits)
  
  # Fill in a vector of 0 and 1 for each significant target hit to indicate if it 
  # is a known cancer gene
  cancer_vect <- unlist(lapply(significant_target_hits, function(x) {
    val1 <- ifelse(x %fin% known_cancer_genes_table$primary_gene_names, 1, 0)
    val2 <- ifelse(x %fin% tfcancer_df$gene, 1, 0)
  }))
  
  # Plot the enrichment (fraction of CGC genes at/ above given rank)
  ranks <- 1:length(significant_target_hits)
  frac_cancer_vect <- c()
  count_of_cancer_genes <- 0
  for (i in 1:length(significant_target_hits)) {
    val_of_curr_tg <- cancer_vect[i]
    count_of_cancer_genes <- count_of_cancer_genes + val_of_curr_tg
    frac <- count_of_cancer_genes / i
    frac_cancer_vect <- c(frac_cancer_vect, frac)
  }
  
  plot(ranks, frac_cancer_vect, main = "Enrichment of Significant Target Genes in Known Cancer Genes",
       xlab = "Rank", ylab = "Fraction of Target Genes that are Known Cancer Genes")
}

# Call function
fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/Mutation, eQTL/Enrichment in Known Cancer Genes (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 650, height = 350)
plot_cancer_enrichment(master_df_mut_sig, known_cancer_genes_table, tfcancer_df)
dev.off()

fn <- paste(output_path, "TP53 (Test)/Non-Tumor-Normal-Matched/CNA, eQTL/Enrichment in Known Cancer Genes (I-Protein, TMM, bucketCNA, iciTotFrac).png", 
            sep = "")
png(fn, width = 650, height = 350)
plot_cancer_enrichment(master_df_cna_sig, known_cancer_genes_table, tfcancer_df)
dev.off()


############################################################
############################################################
#### PERFORM GSEA USING STRINGDB
############################################################
############################################################
#' Function to get the enrichment using the STRINGdB package, which has functions 
#' to calculate enrichment in alternative sets (such as PubMed, etc.)
#' @param master_df a master DF produced from run_linear_model() that has q-values
#' @param qval_thres a q-value threshold for what we consider a significant model hit
#' @param top_n alternative to look at the top N genes at the given significance threshold
#' @param background if we ran the model on a subset of genes to begin with (e.g. 
#' metabolic genes), provide this gene set as a background for enrichment (gene symbols)
get_stringdb_enrichment <- function(master_df, qval_thres, top_n, background) {
  # Instantiate the string DB object 
  string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=0, 
                            input_directory="")
  
  # Set the background, if needed
  if(!is.na(background)) {
    background <- as.data.frame(background)
    colnames(background) <- "gene"
    background <- string_db$map(background, "gene", removeUnmappedRows = TRUE)
    
    string_db <- string_db$set_background(background$STRING_id)
  }
  
  # Calculate the enrichment for all sources  
  hits <- master_df[master_df$q.value < qval_thres, 'T_k.name']
  stringdb_enrichment <- string_db$get_enrichment(hits[1:top_n])
  print(head(stringdb_enrichment, n=20))
  
  return(stringdb_enrichment)
}

# Call function
stringdb_enrichment_tp53 <- get_stringdb_enrichment(allgenes_p53_panCts_inclBRCAsubtypes, 0.2, NA)
stringdb_enrichment_pik3ca <- get_stringdb_enrichment(allgenes_pik3ca_panCts_inclBRCAsubtypes, 0.2, NA)


############################################################
############################################################
#### CALCULATE ENRICHMENT USING NETWORK CONFIDENCE SCORES
#### EMPIRICAL P-VALUES
############################################################
############################################################
#' Use STRING network and associated confidence scores to generate a 
#' distribution of scores for the given GOI. See where each significant 
#' target falls on this distribution and use that to perform a hypothesis 
#' test and generate a p-value, which is returned.
#' @param string_nw a STRING network downloaded from the STRING website, or processed
#' using the STRINGdb package
#' @param goi the Uniprot/Swissprot ID for a given driver of interest
#' @param goi_name the Hugo symbol/ gene name for the given driver of interest
#' @param master_df an output data frame with results for the given GOI
#' @param qval_thres a q-value threshold for signficance (e.g. 0.1)
#' @param num_trials the number of times to draw from the confidences to generate
#' a distribution (e.g. 100)
generate_empirical_p_string <- function(string_nw, goi, goi_name, master_df, 
                                        qval_thres, num_trials) {
  # Subset the data frame to only the given GOI
  master_df_goi <- NA
  if("R_i" %in% colnames(master_df)) {
    master_df_goi <- master_df[grepl(goi, master_df$R_i),]
  } else {
    master_df_goi <- master_df[grepl(goi, master_df$term),]
  }
  
  # Subset the data frame using the q-value threshold
  master_df_goi_sub <- master_df_goi[master_df_goi$q.value < qval_thres,]
  
  # Subset the STRING network to only the given GOI
  #string_nw_goi <- string_nw[((string_nw$node1 == goi_name) | (string_nw$node2 == goi_name)),]  # will already be subsetted by our max # of interactions

  # Get the number of significant hits for this GOI
  num_sig_hits <- nrow(master_df_goi_sub)
  print(paste("Num. Significant Hits:", num_sig_hits))

  # We will draw this number of hits a given number of times to generate a distribution
  if(num_sig_hits > 0) {
    distr_vals <- unlist(lapply(1:num_trials, function(i) {
      v <- mean(sample(string_nw$combined_score, size = num_sig_hits, replace = TRUE))
      return(v)
    }))
    hist(distr_vals, main = "", xlab = paste("\nRandomized STRING Conf. Scores for", 
                                              paste0(goi_name, paste0("\nSamp. Size = # Hits at q<", 
                                                                 paste0(qval_thres, paste0(", # Trials = ", num_trials))))),
         ylab = "Frequency Across Trials", col = "#0072B599") #xlim = c(min(distr_vals), max(distr_vals)))
     
    # Now, get the confidences for each of our actual hits
    real_conf_scores <- unlist(lapply(master_df_goi_sub$T_k.name, function(x) {
      print(x)
      if(x %fin% unique(c(string_nw$node1_name, string_nw$node2_name))) {
        return(as.numeric(unlist(string_nw[(string_nw$node1_name == x) | (string_nw$node2_name == x), 
                             'combined_score'])))
      } else {return(NA)}
    }))
      #as.numeric(string_nw_goi[(string_nw_goi$node2 == x) | (string_nw_goi$node1 == x), 
      #'combined_score'])}))
    
    real_conf_scores_mean <- mean(real_conf_scores[!is.na(real_conf_scores)])
    print(real_conf_scores_mean)
    
    # Plot where we fall on this distribution
    abline(v = real_conf_scores_mean, col = "#BC3C29FF")
    
    # Perform a simple one-sample t-test (assume that we have a random sample
    # of values from a theoretical pool of genes that are called significant by 
    # our model for having "real" relationships to our GOI). We compare this against 
    # the distribution of confidences for EVERY gene in relationship to our GOI. 
    # Is the mean confidences for our "real" targets greater than the general mean
    # confidence score for our GOI?
    ttest_res <- t.test(real_conf_scores, mu = mean(distr_vals), alternative = "greater")
    print(ttest_res)
    
    # Add p-value to histogram
    x_val <- real_conf_scores_mean
    if(real_conf_scores_mean > max(distr_vals)) {x_val <- max(distr_vals) - 50}
    text(x = x_val, y  = round(num_trials/4.5), cex=0.8,
         labels = paste0("One-sided t-test p = ", round(ttest_res$p.value, digits = 4)))
    text(x = x_val, y  = round(num_trials/4.5)-10, cex=0.8,
         labels = paste0("Mean of Real Conf. Vals ", round(real_conf_scores_mean, digits = 4)))
    
    return(ttest_res$p.value)
  }
  else {
    print("No significant hits for the given GOI at the given significance threshold.")
    return(NA)
  }
}


#' Process STRINGdb network to get interactions for the given GOI with the combined score
#' @param goi the gene name of the GOI
#' @param all_genes_id_conv a gene ID conversion file from BioMart
process_stringdb_network_for_goi <- function(goi, all_genes_id_conv) {
  # Download the STRINGdb network
  options(download.file.method="libcurl")
  string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 0, 
                            input_directory = "")
  gene = string_db$mp(tolower(goi))
  
  # Get the neighbors for a given gene
  neighbors <- string_db$get_neighbors(tolower(goi))
  interactions <- string_db$get_interactions(c(gene, neighbors))
  
  # Subset interactions that only include gene of interest
  interactions$from <- unlist(lapply(interactions$from, function(x)
    unlist(strsplit(x, ".", fixed = TRUE))[2]))
  interactions$to <- unlist(lapply(interactions$to, function(x)
    unlist(strsplit(x, ".", fixed = TRUE))[2]))
  ensp <- unique(all_genes_id_conv[all_genes_id_conv$external_gene_name == goi, 
                                   'ensembl_peptide_id'])
  interactions <- interactions[(interactions$from %in% ensp) | 
                                 (interactions$to %in% ensp),]
  # Add additional gene symbols
  interactions$from_name <- unlist(lapply(interactions$from, function(x) {
    new_id <- unique(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_peptide_id == x, 
                                              'external_gene_name']))
    if(length(new_id) > 1) {new_id <- new_id[1]}
    if (length(new_id) == 0) {new_id <- NA}
    return(new_id)
  }))
  interactions$to_name <- unlist(lapply(interactions$to, function(x) {
    new_id <- unique(unlist(all_genes_id_conv[all_genes_id_conv$ensembl_peptide_id == x, 
                                              'external_gene_name']))
    if(length(new_id) > 1) {new_id <- new_id[1]}
    if (length(new_id) == 0) {new_id <- NA}
    return(new_id)
  }))
  
  return(na.omit(distinct(interactions)))
}

# Import the String network file
string_nw_full <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/STRING/string.9606.protein.links.v11.5.namesAdded.txt")

# Limit to given GOI
goi <- "P04637"
goi_name <- "TP53"

# goi_name <- "PIK3CA"

# Alternatively, get and process the STRING db network
interactions <- process_stringdb_network_for_goi(goi, all_genes_id_conv)

# Call function
generate_empirical_p_string(string_nw, goi, goi_name, allgenes_topDrivers_0.025_ucsf_noLasso,
                           0.2, 1000)
#generate_empirical_p_string(interactions, goi, goi_name, allgenes_topDrivers_0.025_ucsf_noLasso,
                           # 0.2, 1000)


############################################################
############################################################
#### IMPORT THE NECESSARY GENE SETS
############################################################
############################################################

### CANCER GENE SETS ###
# Import a table containing a compiled list of known cancer genes
known_cancer_genes_table <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/GRCh38_driver_gene_list.tsv", sep = "\t",
                                       header = TRUE, check.names = FALSE, comment.char = "#", skip = 10)
known_cancer_genes <- known_cancer_genes_table$primary_gene_names
vogelstein_genes <- known_cancer_genes_table[grepl("V", known_cancer_genes_table$cancer_driver_status), 'primary_gene_names']
curated_cancer_genes <- known_cancer_genes_table[grepl("V|K|L|B", known_cancer_genes_table$cancer_driver_status), 'primary_gene_names']
cgc_genes <- known_cancer_genes_table[grepl("C", known_cancer_genes_table$cancer_driver_status), 'primary_gene_names']

### TFCANCER TABLE ###
# http://lcbb.swjtu.edu.cn/tfcancer/
tfcancer_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/"
# FOR BRCA
tfcancer_df <- read.csv(paste(tfcancer_path, "BRCA Data/Validation_Files/BRCA_TFcancer.csv", sep = ""),
                        header = TRUE, check.names = FALSE)

### CHIP-EAT/ eCLIP TARGETS ###
# Import TP53-specific ChIP-eat targets from Unibind
tp53_chipeat_targets <- read.table(paste0(main_path, "ChIP-eat/tp53_chipeat_targets.txt"), header = FALSE)[,1]
# SF3B1 eCLIP targets from POSTAR3
sf3b1_eclip_targets <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/POSTAR3_CLIPdb_SF3B1_targets.csv", 
                                header = TRUE, check.names = FALSE)
sf3b1_eclip_targets_numBindSites_gr100 <- sf3b1_eclip_targets[sf3b1_eclip_targets$`Binding site records` > 100, 'Target gene symbol']
sf3b1_eclip_targets_numBindSites_gr125 <- sf3b1_eclip_targets[sf3b1_eclip_targets$`Binding site records` > 125, 'Target gene symbol']


### CURATED TARGETS FROM PAPERS ###
# TP53 target genes from Fischer et al., 2017: https://pubmed.ncbi.nlm.nih.gov/28288132/
tp53_curated_targets <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Validation_Files/TP53_Targets_Fischer_2017_Oncogene.csv",
                                 header = FALSE, check.names = FALSE)[,1]
# PIK3CA target genes from Cizkova et al., 2010: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3012715/
pik3ca_curated_targs <- c("WNT5A", "TCF7L2", "MSX2", "TNFRSF11B", "SEC14L2", "MSX2", 
                          "TFAP2B", "NRIP3", "CYP4Z1", "CYPZ2P", "SLC40A1", "LTF", "LIMCH1")


### TF-TARGET INFORMATION FROM TOOLS/ DATABASES ###
# TRRUST
# https://www.grnpedia.org/trrust/
trrust_df <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Validation_Files/trrust_rawdata.human.txt",
                        header = FALSE, check.names = FALSE)
trrust_df_tp53 <- trrust_df[(trrust_df[,1] == "TP53") | (trrust_df[,2] == "TP53"),]
tp53_trrust_targets <- unique(c(trrust_df_tp53$V1, trrust_df_tp53$V2))
tp53_trrust_targets <- tp53_trrust_targets[tp53_trrust_targets != "TP53"]  # 195 targets

# Or, look just at cases where TP53 is upstream (activation or repression)
trrust_df_tp53_upstr <- trrust_df[trrust_df[,1] == "TP53",]
tp53_trrust_targets_upstr <- unique(trrust_df_tp53$V2)
tp53_trrust_targets_upstr <- tp53_trrust_targets_upstr[tp53_trrust_targets_upstr != "TP53"]  # 163 targets


# hTFtarget
# http://bioinfo.life.hust.edu.cn/hTFtarget
tf_target_df <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Validation_Files/TF-Target-information.txt",
                           header = TRUE, check.names = FALSE, sep = "\t")
tf_target_df_breast <- tf_target_df[grepl("breast", tf_target_df$tissue),]
tf_target_df_tp53 <- tf_target_df_breast[grepl("P53", tf_target_df_breast$TF),]
tf_target_tp53 <- tf_target_df_tp53$target
tf_target_tp53  <- tf_target_tp53[tf_target_tp53 != "TP53"]


# DoRothEA 
# https://bioconductor.org/packages/release/data/experiment/vignettes/dorothea/inst/doc/dorothea.html

dorothea_net <- dorothea::dorothea_hs
tp53_dorothea_targets <- unlist(dorothea_net[dorothea_net$tf == "TP53", 'target'])

nfkb_dorothea_targets <- unlist(dorothea_net[grepl("NFKB", dorothea_net$tf), 'target'])
foxo_dorothea_targets <- unlist(dorothea_net[grepl("FOXO", dorothea_net$tf), 'target'])
creb_dorothea_targets <- unlist(dorothea_net[grepl("CREB", dorothea_net$tf), 'target'])
gata3_dorothea_targets <- unlist(dorothea_net[grepl("GATA3", dorothea_net$tf), 'target'])
foxa1_dorothea_targets <- unlist(dorothea_net[grepl("FOXA1", dorothea_net$tf), 'target'])


### GENETIC PATHWAYS ###
# Import the network interaction files for the genes of interest
nw_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/Network_Data/"
nw_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/"


### REACTOME PATHWAYS ###
pi3k_akt_signaling_reactome_pw <- read.csv(paste0(nw_path, "Reactome/PI3K-AKT Signaling in Cancer_R-HSA-2219528.csv"), 
                                           header = TRUE, check.names = FALSE)
pi3k_akt_signaling_reactome_pw_gns <- unlist(lapply(pi3k_akt_signaling_reactome_pw$MoleculeName, function(x) 
  unlist(strsplit(x, " ", fixed = TRUE))[2]))
pi3k_akt_signaling_reactome_pw_gns <- unique(pi3k_akt_signaling_reactome_pw_gns[pi3k_akt_signaling_reactome_pw_gns != "PIK3CA"])

pi3k_constitutative_signal_reactome_pw <- read.csv(paste0(nw_path, "Reactome/Constituative Signaling by Aberrant PI3K in Cancer_R-HSA-2219530.csv"), 
                                                   header = TRUE, check.names = FALSE)
pi3k_constitutative_signal_reactome_pw_gns <- unlist(lapply(pi3k_constitutative_signal_reactome_pw$MoleculeName, function(x) 
  unlist(strsplit(x, " ", fixed = TRUE))[2]))
pi3k_constitutative_signal_reactome_pw_gns <- unique(pi3k_constitutative_signal_reactome_pw_gns[pi3k_constitutative_signal_reactome_pw_gns != "PIK3CA"])

wnt_signaling_reactome_pw <- read.csv(paste0(nw_path, "Reactome/Wnt Signaling_R-HSA-195721.csv"), 
                                      header = TRUE, check.names = FALSE)
wnt_signaling_reactome_pw_gns <- unlist(lapply(wnt_signaling_reactome_pw$MoleculeName, function(x) 
  unlist(strsplit(x, " ", fixed = TRUE))[2]))
wnt_signaling_reactome_pw_gns <- unique(wnt_signaling_reactome_pw_gns[wnt_signaling_reactome_pw_gns != "PIK3CA"])


mtorc_signaling_reactome_pw <- read.csv(paste0(nw_path, "Reactome/mTORC Signaling_R-HSA-165159.csv"), 
                                        header = TRUE, check.names = FALSE)
mtorc_signaling_reactome_pw_gns <- unlist(lapply(mtorc_signaling_reactome_pw$MoleculeName, function(x) 
  unlist(strsplit(x, " ", fixed = TRUE))[2]))
mtorc_signaling_reactome_pw_gns <- unique(mtorc_signaling_reactome_pw_gns[mtorc_signaling_reactome_pw_gns != "PIK3CA"])

tp53_regulates_metabolism_pw <- read.csv(paste0(nw_path, "Reactome/TP53 Regulates Metabolic Genes_R-HSA-5628897.csv"), 
                                         header = TRUE, check.names = FALSE)
tp53_regulates_metabolism_pw_gns <- unlist(lapply(tp53_regulates_metabolism_pw$MoleculeName, function(x) 
  unlist(strsplit(x, " ", fixed = TRUE))[2]))
tp53_regulates_metabolism_pw_gns <- unique(tp53_regulates_metabolism_pw_gns[tp53_regulates_metabolism_pw_gns != "TP53"])
#tp53_regulates_metabolism_pw_gns <- setdiff(tp53_regulates_metabolism_pw_gns, metabol_p53$T_k.name)

tp53_transcriptional_reg_pw <- read.csv(paste0(nw_path, "Reactome/Transcriptional Regiulation by TP53_R-HSA-3700989.csv"), 
                                         header = TRUE, check.names = FALSE)
tp53_transcriptional_reg_pw_gns <- unlist(lapply(tp53_transcriptional_reg_pw$MoleculeName, function(x) 
  unlist(strsplit(x, " ", fixed = TRUE))[2]))
tp53_transcriptional_reg_pw_gns <- unique(tp53_transcriptional_reg_pw_gns[tp53_transcriptional_reg_pw_gns != "TP53"])


# PIK3CA's TF targets are primarily NFKB, CREB, and FOXO
nfkb_trrust_targets <- unique(trrust_df[trrust_df[,1] == "NFKB1", 'V2'])
nfkb_trrust_targets <- nfkb_trrust_targets[nfkb_trrust_targets != "NFKB1"]  # 301 targets
foxo_trrust_targets <- unique(trrust_df[trrust_df[,1] == "FOXO1", 'V2'])
foxo_trrust_targets <- foxo_trrust_targets[foxo_trrust_targets != "FOXO1"]  # 24 targets
pik3ca_downstream_targets <- unique(c(nfkb_trrust_targets, foxo_trrust_targets))


### KEGG PATHWAYS ###
tp53_kegg_pathway <- read.table(paste0(nw_path, "KEGG/KEGG_P53_SIGNALING_PATHWAY.v7.5.1.txt"), sep = "\t", header = TRUE)
tp53_kegg_pathway_genes <- unlist(strsplit(tp53_kegg_pathway$KEGG_P53_SIGNALING_PATHWAY[19], ",", fixed = TRUE))

pik3ca_kegg_pathway <- read.table(paste0(nw_path, "KEGG/KEGG_PI3K-AKT_SIGNALING_PATHWAY.v7.5.1.txt"), sep = "\t", header = FALSE)
pik3ca_kegg_pathway_genes <- unlist(lapply(1:nrow(pik3ca_kegg_pathway), function(r) {
  row <- pik3ca_kegg_pathway$V1[r]
  gene <- unlist(strsplit(unlist(strsplit(row, "(RefSeq) ", fixed = TRUE))[2], ", ", fixed = TRUE))[1]
  return(gene)
}))
pik3ca_upstr_genes <- c("EGF", "GRB2", "SOS1", "SOS2", "HRAS", "TLR2", "TLR4", "IRS1", "RAF1", "MEK", "ERK", "RAC1", 
                        "SYK", "CD19", "JAK1", "JAK2", "JAK3", "ITGA1", "ITGA2", "ITGA2B", "ITGA3", "ITGA4", "ITGA5",
                        "ITGA6", "ITGA7", "ITGA8", "ITGA9", "ITGA10", "ITGA11", "ITGAV", "ITGAB1", "ITGB3", "ITGB4", 
                        "ITGB5", "ITGB6", "ITGB7", "ITGB8", "FAK", "MAGI1", "MAGI2", "PTEN", "KRAS", "THEM4", "PDPK1", 
                        "NRAS", "HSP90AA1", "HSP90AB1", "HSP90B1", "IGH", "PIK3AP1", "MTCP1", "TCL1A", "TCL1B", "PHLPP2",
                        "PHLPP1", "CRTC2", "CDC37", "GNB5", "GNB1", "GNB2", "GNB3", "GNG3", "GNG4", "GNG5", "GNG7", "GNG10",
                        "GNG11", "GNGT1", "GNGT2", "GNG13", "GNG2", "GNG12", "GNGB4", "GNG8", "LPAR6", "CHRM1", "CHRM2", 
                        "LPAR1", "F2R", "LPAR3", "LPAR4", "LPAR5", "LPAR2", "PTK2", "LAMC3", "CHAD", "COL1A1", "COL1A2",
                        "COL2A1", "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6", "COL6A1", "COL6A2", "COL6A3",
                        "COL9A1", "COL9A2", "COL9A3", "COMP", "COL6A6", "LAMB4", "FN1", "COL6A5", "LAMA1", "TNC", "IBSP",
                        "LAMA2", "LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMB3", "LAMC1", "LAMC2", "RELN", "TNN",
                        "SPP1", "THBS1", "THBS1", "THBS3", "TNR", "TNXB", "VTN", "VWF", "CSF3", "CSH1", "CSH2", "EPO",
                        "GH1", "GH2", "IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNA10", "IFNA13",
                        "IFNA14", "IFNA16", "IFNA17", "IFNA21", "IFNB1", "IL2", "IL3", "IL4", "IL6", "IL7", "OSM", "PRL",
                        "CSF1", "EFNA1", "EFNA2", "EFNA3", "EFNA4", "EFNA5", "EREG", "FGF1", "FGF2", "FGF3", "FGF4", "FGF5",
                        "FGF6", "FGF7", "FGF8", "FGF9", "FGF10", "VEGFD", "FLT3LG", "FGF20", "FGF21", "FGF22", "ANGPT1", 
                        "ANGPT2", "HGF", "IGF1", "IGF2", "INS", "AREG", "KITLG", "NGF", "NTF3", "NTF4", "ANGPT4", "PDGFA",
                        "PDGFB", "PGF", "PDGFC", "BDNF", "TGFA", "VEGFA", "VEGFB", "VEGFC", "PDGFD", "FGF23", "FGF18",
                        "FGF17", "FGF16", "FGF19", "CSF1R", "EGFR", "EPHA2", "ERBB2", "ERBB3", "ERBB4", "FGFR1", "FGFR3",
                        "FGFR2", "FGFR4", "FLT1", "FLT3", "FLT4", "IGF1R", "INSR", "KDR", "KIT", "MET", "NGFR", "NTRK1",
                        "NTRK2", "PDGFRA", "PDGFRB", "TEK")

pik3ca_kegg_pathway_genes_upstr <- pik3ca_kegg_pathway_genes[!(pik3ca_kegg_pathway_genes %fin% pik3ca_upstr_genes)] # limits to a set of 120 genes

pik3ca_kegg_pathway_genes_downstrDNA <- c("PCK1", "PCK2", "G6PC1", "G6PC2", "G6PC3", "CCND1", "CDK2", "CDK4", "CDK6", 
                                          "CDKN1B", "RBL2", "FASLG", "BCL2L11", "BCL2L1", "BCL2", "MCL1", "MYB", "TP53", 
                                          "AKT", "PDK1")  # 20 genes

kras_kegg_pathway <- read.csv(paste0(nw_path, "KEGG/KEGG_MAPK_SIGNALING_PATHWAY.v2022.1.csv"), 
                              header = FALSE, check.names = FALSE)
kras_kegg_pathway_genes <- as.character(unlist(strsplit(kras_kegg_pathway$V2[20], ",", fixed = TRUE)))



### HUMANBASE TARGETS ###
# https://hb.flatironinstitute.org/

goi <- "TP53"
goi <- "PIK3CA"

tp53_mamm_epi_nw <- read.table(paste0(nw_path, "hb_mammary_epithelium_tp53_network.txt"), 
                               header = TRUE, sep = ",")
tp53_mamm_gl_nw <- read.table(paste0(nw_path, "hb_mammary_gland_tp53_network.txt"), 
                              header = TRUE, sep = ",")
pik3ca_mamm_epi_nw <- read.table(paste0(nw_path, "hb_mammary_epithelium_pik3ca_network.txt"), 
                                 header = TRUE, sep = ",")
pik3ca_mamm_gl_nw <- read.table(paste0(nw_path, "hb_mammary_gland_pik3ca_network.txt"), 
                                header = TRUE, sep = ",")

tp53_colon_nw <- read.table(paste0(nw_path, "HumanBase/hb_colon_tp53_network.txt"), 
                            header = TRUE, sep = ",")

tp53_global_nw <- read.table(paste0(nw_path, "tp53_humanbase_global_top50.txt"), header = TRUE, sep = ",")
pik3ca_global_nw <- read.table(paste0(nw_path, "pik3ca_humanbase_global_top50.txt"), header = TRUE, sep = ",")

# Get rid of any white space
tp53_mamm_epi_nw <- as.data.frame(lapply(tp53_mamm_epi_nw, trimws))
tp53_mamm_gl_nw <- as.data.frame(lapply(tp53_mamm_gl_nw, trimws))
pik3ca_mamm_epi_nw <- as.data.frame(lapply(pik3ca_mamm_epi_nw, trimws))
pik3ca_mamm_gl_nw <- as.data.frame(lapply(pik3ca_mamm_gl_nw, trimws))

tp53_colon_nw <- as.data.frame(lapply(tp53_colon_nw, trimws))


tp53_global_nw <- as.data.frame(lapply(tp53_global_nw, trimws))
pik3ca_global_nw <- as.data.frame(lapply(pik3ca_global_nw, trimws))

# Limit the interactions to just gois
tp53_mamm_epi_nw <- tp53_mamm_epi_nw[(tp53_mamm_epi_nw$GENE1 == "TP53") | (tp53_mamm_epi_nw$GENE2 == "TP53"),]
tp53_mamm_gl_nw <- tp53_mamm_gl_nw[(tp53_mamm_gl_nw$GENE1 == "TP53") | (tp53_mamm_gl_nw$GENE2 == "TP53"),]
pik3ca_mamm_epi_nw <- pik3ca_mamm_epi_nw[(pik3ca_mamm_epi_nw$GENE1 == "TP53") | (pik3ca_mamm_epi_nw$GENE2 == "PIK3CA"),]
pik3ca_mamm_gl_nw <- pik3ca_mamm_gl_nw[(pik3ca_mamm_gl_nw$GENE1 == "TP53") | (pik3ca_mamm_gl_nw$GENE2 == "PIK3CA"),]

#tp53_colon_nw <- tp53_colon_nw[(tp53_colon_nw$GENE1 == "TP53") | (tp53_colon_nw$GENE2 == "TP53"),]


# Limit the interactions to just those above a specific confidence threshold
conf_thres <- 0.5

tp53_mamm_epi_nw_sig <- tp53_mamm_epi_nw[tp53_mamm_epi_nw$WEIGHT > conf_thres,]
tp53_mamm_gl_nw_sig <- tp53_mamm_gl_nw[tp53_mamm_gl_nw$WEIGHT > conf_thres,]
pik3ca_mamm_epi_nw_sig <- pik3ca_mamm_epi_nw[pik3ca_mamm_epi_nw$WEIGHT > conf_thres,]
pik3ca_mamm_gl_nw_sig <- pik3ca_mamm_gl_nw[pik3ca_mamm_gl_nw$WEIGHT > conf_thres,]

# Get a union of the significant targets for each goi
tp53_nw_targs <- unique(c(tp53_mamm_epi_nw_sig$GENE1, tp53_mamm_epi_nw_sig$GENE2, 
                          tp53_mamm_gl_nw_sig$GENE1, tp53_mamm_gl_nw_sig$GENE2))
tp53_nw_targs <- tp53_nw_targs[tp53_nw_targs != "TP53"]
pik3ca_nw_targs <- unique(c(pik3ca_mamm_epi_nw_sig$GENE1, pik3ca_mamm_epi_nw_sig$GENE2, 
                            pik3ca_mamm_gl_nw_sig$GENE1, pik3ca_mamm_gl_nw_sig$GENE2))
pik3ca_nw_targs <- pik3ca_nw_targs[pik3ca_nw_targs != "PIK3CA"]

tp53_colon_nw_targs <- unique(tp53_colon_nw$SYMBOL)
tp53_colon_nw_targs <- tp53_colon_nw_targs[tp53_colon_nw_targs != "TP53"]

### STRING ###
# https://string-db.org/

# Import STRING network targets (top 100)
tp53_string_nw_targs <- read.table(paste0(nw_path, "STRING/TP53_interactors_string_top100.txt"),
                                   header = FALSE)[,1]
pik3ca_string_nw_targs <- read.table(paste0(nw_path, "STRING/PIK3CA_interactors_string_top100.txt"),
                                     header = FALSE)[,1]
kras_string_nw_targs <- unique(unlist(read.csv(paste0(nw_path, "STRING/kras_string_interactions_top100.csv"),
                                   header = T, check.names = F)[,1:2]))
kras_string_nw_targs <- kras_string_nw_targs[kras_string_nw_targs != "KRAS"]
idh1_string_nw_targs <- unique(unlist(read.csv(paste0(nw_path, "STRING/idh1_string_interactions_top100.csv"),
                                               header = T, check.names = F)[,1:2]))
idh1_string_nw_targs <- idh1_string_nw_targs[idh1_string_nw_targs != "IDH1"]

# Import the top 500 STRING network targets
tp53_string_nw_targs_top500 <- read.table(paste0(nw_path, "STRING/tp53_string_interactions_top500.csv"),
                                          header = TRUE, sep = ",")
tp53_string_nw_targs_top500 <- unique(c(tp53_string_nw_targs_top500$node1, tp53_string_nw_targs_top500$node2))
tp53_string_nw_targs_top500 <- tp53_string_nw_targs_top500[tp53_string_nw_targs_top500 != "TP53"]
pik3ca_string_nw_targs_top500 <- read.table(paste0(nw_path, "STRING/pik3ca_string_interactions_top500.csv"),
                                            header = FALSE, sep = ",")
pik3ca_string_nw_targs_top500 <- unique(c(pik3ca_string_nw_targs_top500$V1, pik3ca_string_nw_targs_top500$V2))
pik3ca_string_nw_targs_top500 <- pik3ca_string_nw_targs_top500[pik3ca_string_nw_targs_top500 != "PIK3CA"]


# Other STRING target lists (from "process_string_data.R")
tp53_neighbors_dist2_weight5 <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project FIles/Saved Output Data Files/Pan-Cancer/Network_Data/string_tp53_neighbors_dist2_weight5.txt",
                                           header = TRUE)[,'SYMBOL']
pik3ca_neighbors_dist2_weight5 <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Network_Data/string_pik3ca_neighbors_dist2_weight5.txt",
                                             header = TRUE)[,'SYMBOL']

# Note: 100 = 0.1, 500 = 0.5, etc.
tp53_neighbors_conf100 <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Network_Data/string_tp53_neighbors_confThres100.txt",
                                     header = TRUE)[,'SYMBOL']
pik3ca_neighbors_conf100 <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Network_Data/string_pik3ca_neighbors_confThres100.txt",
                                       header = TRUE)[,'SYMBOL']

tp53_neighbors_conf500 <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Network_Data/string_tp53_neighbors_confThres500.csv",
                                   header = TRUE, check.names = FALSE)
tp53_neighbors_conf500 <- unique(c(tp53_neighbors_conf500$`#node1`, tp53_neighbors_conf500$node2))
tp53_neighbors_conf500 <- tp53_neighbors_conf500[tp53_neighbors_conf500 != "TP53"]

pik3ca_neighbors_conf500 <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Network_Data/string_pik3ca_neighbors_confThres500.csv",
                                     header = TRUE, check.names = FALSE)
pik3ca_neighbors_conf500 <- unique(c(pik3ca_neighbors_conf500$`#node1`, pik3ca_neighbors_conf500$node2))
pik3ca_neighbors_conf500 <- pik3ca_neighbors_conf500[pik3ca_neighbors_conf500 != "PIK3CA"]


tp53_neighbors_top500 <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/STRING/tp53_string_interactions_top500.csv",
                                  header = TRUE, check.names = FALSE)


### Human Protein Atlas Co-Expression Cluster ###
# https://www.proteinatlas.org/
pik3ca_humanProtAtlas_cluster <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Validation_Files/pik3ca_humanProtAtlas_expr_cluster.csv", 
                                          header = TRUE, check.names = FALSE)[, 'Gene']


### BREAST DEGS AND BIOMARKERS ###
# Breast DEGs (tumor vs. normal) from Pranavathiyani et al., 2019, Genes and Diseases: https://www.sciencedirect.com/science/article/pii/S2352304218300801?via%3Dihub#ec-research-data
breast_degs <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/Validation_Files/breast_degs_pranavathiyani_2009_genesAndDiseases.csv",
                        header = TRUE, check.names = FALSE)[, 'DEG']
breast_degs_p53 <- breast_degs[breast_degs != "TP53"]
breast_degs_pik3ca <- breast_degs[breast_degs != "PIK3CA"]

# Breast biomarker genes from meta-analysis, Abba et al., 2010, Biomarker Insights: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2978930/#SD2
breast_biomarkers <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/BRCA Data/Validation_Files/breast_top_biomarkers_abba_2010_biomarkerInsights.csv",
                              header = TRUE, check.names = FALSE)[, 'Gene']
breast_biomarkers_p53 <- breast_biomarkers[breast_biomarkers != "TP53"]
breast_biomarkers_pik3ca <- breast_biomarkers[breast_biomarkers != "PIK3CA"]





























### ADDITIONAL FUNCTIONS FOR TFCANCER TABLE PROCESSING ###


# FOR PAN-CANCER (must combine files)
#' Creates a TFcancer file that is the combination of all the individual cancer
#' type TFcancer files (33 cancer types). Adds a column for the cancer type.
#' @param tfcancer_path local path to the TF cancer files
merge_tfcancer_dfs <- function(tfcancer_path) {
  # List all the files in the path
  tfcancer_files <- list.files(paste(tfcancer_path, "TCGA Data (ALL)/Validation_Files/TFcancer/"))
  
  # Read in the first DF 
  tfcancer_df <- read.csv(paste(tfcancer_path, tfcancer_files[1], sep = ""), header = TRUE,
                          check.names = FALSE)
  tfcancer_df$cancer_type <- rep(unlist(strsplit(tfcancer_files[1], ".", fixed = TRUE))[1], 
                                 nrow(tfcancer_df))
  # Merge the rest into this DF with a column for cancer type
  for (i in 2:length(tfcancer_files)) {
    new_tfcancer_df <- read.csv(paste(tfcancer_path, tfcancer_files[i], sep = ""), header = TRUE,
                                check.names = FALSE)
    new_tfcancer_df$cancer_type <- rep(unlist(strsplit(tfcancer_files[i], ".", fixed = TRUE))[1], 
                                       nrow(new_tfcancer_df))
    tfcancer_df <- rbind(tfcancer_df, new_tfcancer_df)
  }
  return(tfcancer_df)
}

# Call function for pan-cancer
#tfcancer_df <- merge_tfcancer_dfs(tfcancer_path)

# Write this back to a CSV
#fwrite(tfcancer_df, paste(path, "Validation_Files/TFcancer/ALL_CANCERS.csv", sep = ""))


#' Function adjusts the TFcancer data frame to 1. make all the
#' gene names uppercase; 2. separate out the semicolon-separated
#' gene target names into individual rows
#' @param tfcancer_df the imported data frame with the information
#' from the TFcancer database
adjust_tfcancer_df <- function(tfcancer_df) {
  # First, break up all the comma-separated targets into new rows
  new_dfs <- lapply(1:nrow(tfcancer_df), function(i) {
    # Get the target genes into a vector with whitespace stripped
    targ_genes <- unlist(strsplit(stringr::str_replace_all(tfcancer_df[i, 'gene'], fixed(" "), ""), 
                                  ";", fixed = TRUE))
    targ_genes <- targ_genes[!is.na(targ_genes)]
    print(targ_genes)
    
    # For each, create a data frame row with everything duplicated from the original
    # DF except for the gene column. Return a new data frame with all these rows.
    samp_df <- tfcancer_df[i,]
    samp_df$gene <- targ_genes[1]
    for (j in 2:length(targ_genes)) {
      samp_df <- rbind(samp_df, tfcancer_df[i,])
      samp_df[j, 'gene'] <- targ_genes[j]
    }
    return(samp_df)
  })
  
  # Re-combine all these partial DFs back into one
  new_tfcancer_df <- do.call(rbind, new_dfs)
  
  return(new_tfcancer_df)
}

# Call this function
#new_tfcancer_df <- adjust_tfcancer_df(tfcancer_df)

# Write to a new file

# For BRCA
#fwrite(new_tfcancer_df, paste(tfcancer_path, "BRCA Data/Validation_Files/BRCA_TFcancer_adj.csv", sep = ""))
# For pan-cancer
#fwrite(new_tfcancer_df, paste(path, "Validation_Files/TFcancer/ALL_CANCERS_adj.csv", sep = ""))


#' Add columns to the master results data frame to signify 
#' whether the prediction aligns with something found in the
#' TFcancer database. If the TF in question is found in the
#' database, and the gene in question is found to be 'regulated
#' by' that gene, then the regulation type is given, along with
#' all the corresponding information. Similarly, if a gene in 
#' question is found to be in the database, and the TF in 
#' question is found to be 'targeted by' that gene, the 
#' information will also be given. Adjusted master DF is returned.
#' @param master_df_sig the master results DF to be annotated
#' @param tfcancer_df the TFcancer data frame with the information
#' to be used for annotation
annotate_master_w_tfcancer <- function(master_df_sig, tfcancer_df) {
  # Look at each top hit individually, creating a new DF to cbind
  # to the master DF
  tfcancerdf_list <- lapply(1:nrow(master_df_sig), function(i) {
    #tf <- master_df_sig[i, 'R_i.name']
    tf <- "p53"
    gene <- master_df_sig[i, 'T_k.name']
    
    # Subset the TFcancer database to see if the TF or gene are in there
    tf_sub <- tfcancer_df[tfcancer_df$tf == tf,]
    gene_sub <- tfcancer_df[tfcancer_df$gene == gene,]
    
    # Row entries to fill in
    df_to_return <- data.frame(tf_present = FALSE, gene_present= FALSE, any_relationship = FALSE,
                               characteristics = "", regulation = "", hallmark = "", 
                               original_text = "")
    
    if (!(nrow(tf_sub) == 0)) {df_to_return$tf_present <- TRUE}
    if (!(nrow(gene_sub) == 0)) {df_to_return$gene_present <- TRUE}
    
    both_sub <- tf_sub[grepl(gene, tf_sub$gene),]
    if (!nrow(both_sub) == 0) {
      df_to_return$any_relationship <- TRUE
      df_to_return$characteristics <- both_sub$characteristics
      df_to_return$regulation <- both_sub$regulation
      df_to_return$hallmark <- both_sub$hallmark
      df_to_return$original_text <- both_sub$original_text
    }
    
    return(df_to_return)
  })
  
  # Rbind all these DFs together
  tfcancer_res_df <- do.call(rbind, tfcancerdf_list)
  
  # Add this DF to the original master results DF and return
  master_df_sig <- cbind(master_df_sig, tfcancer_res_df)
  
  return(master_df_sig)
}

# Run function
#master_df_sig <- annotate_master_w_tfcancer(master_df_sig, tfcancer_df)

# Write to a file
#fwrite(master_df_sig, paste(main_path, "Linear Model/TP53/Non-Tumor-Normal Matched/iprotein_significant_output_TP53_TMM_rawCNA_iprot_iciTotFrac.csv", sep = ""))





