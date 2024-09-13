############################################################
# Code to Create Figure 4 Visualizations
# Written by Sara Geraghty
# PUBLICATION INFORMATION
############################################################

library(data.table)
library(ggplot2)

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


########################################################################
### DEPMAP CO-DEPENDENCY
########################################################################
# Import co-dependent list
kbtbd2_codep <- read.csv(paste0(PATH, "DepMap/KBTBD2's Top 100 Codependencies for CRISPR (DepMap Public 23Q4+Score Chronos).csv"),
                         header = T, check.names = F)

# Import insulin signaling genes from KEGG/MKEGG
# Downloaded from: https://www.gsea-msigdb.org/gsea/msigdb/cards/KEGG_INSULIN_SIGNALING_PATHWAY
insulin_signaling_genes <- read.csv(paste0(PATH, "KEGG/KEGG_INSULIN_SIGNALING_PATHWAY.v2023.2.Hs.csv"), 
                                    header = T, check.names = F)[17,2]
insulin_signaling_genes <- unlist(strsplit(insulin_signaling_genes, ",", fixed = T)) 
#insulin_signaling_genes <- insulin_signaling_genes[insulin_signaling_genes != "PIK3CA"] # 136 genes

# Import insulin secretion genes (hsa04911), which is significantly downregulated
# Downloaded from: https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+path:hsa04911
insulin_secretion_table <- read.table(paste0(PATH, "KEGG/KEGG_INSULIN_SECRETION_hsa04911.txt"),
                                      sep = "\t", header = F)
insulin_secretion_ids_kegg <- unlist(lapply(insulin_secretion_table[,1], function(entry) {
  entry_spl <- unlist(strsplit(entry, " ", fixed = T))
  entry_spl <- entry_spl[!(entry_spl == "")]
  entry <- entry_spl[grepl("hsa:", entry_spl)]
  entry_id <- unlist(strsplit(entry, ":", fixed = T))[2]
  return(entry_id)
}))
mapping_kegg <- as.data.frame(bitr_kegg(insulin_secretion_ids_kegg, 
                                        fromType = "kegg", toType = "uniprot", 
                                        drop = T, organism = "hsa"))
mapping <- as.data.frame(bitr(mapping_kegg$uniprot, fromType = "UNIPROT", 
                              toType = "SYMBOL", drop = T, OrgDb = "org.Hs.eg.db"))
insulin_secretion_genes <- unique(mapping$SYMBOL)

# Import KEGG PI3K signaling genes
pi3k_signaling_pathway <- read.table(paste0(PATH, "KEGG/KEGG_PI3K-AKT_SIGNALING_PATHWAY.v7.5.1.txt"),
                                 header = F, sep = '\t')
pi3k_kegg_pathway_genes <- unlist(lapply(1:nrow(pi3k_signaling_pathway), function(r) {
  row <- pi3k_signaling_pathway$V1[r]
  gene <- unlist(strsplit(unlist(strsplit(row, "(RefSeq) ", fixed = TRUE))[2], ", ", fixed = TRUE))[1]
  return(gene)
}))

# Option 1: Bar chart of top N co-dependent genes from DepMap (CRISPR or RNAi),
# colored by membership in relevant pathways

#' Create a descending bar plot of top N co-dependent genes from DepMap
#' colored by membership in relevant pathways
#' @param codep_df a codependency data frame from DepMap
#' @param pw_of_interest a vector of genes in the pathway of interest
#' @param pw_of_interest_name the name of the pathway of interest for labeling
#' @param N the number of top co-dependent genes to display
#' @param goi the name of the gene of interest, for labeling
create_codependency_barplot <- function(codep_df, pw_of_interest, pw_of_interest_name,
                                        N, goi) {
  # Limit to just the top N co-dependent genes (already sorted by abs(codep_score))
  codep_df_sub <- codep_df[1:N,]
  
  # Add a column to denote whether each gene is in the pathway of interest
  codep_df_sub$in_pw <- unlist(lapply(codep_df_sub$Gene, function(g) 
    ifelse(g %fin% pw_of_interest, "Yes", "No")))
  
  # Create barplot
  g <- ggplot(codep_df_sub, aes(x = reorder(Gene, -abs(Correlation)), 
                                y = Correlation, fill = as.factor(in_pw))) +
    geom_bar(stat = "identity", position = position_dodge()) +
    xlab(paste(goi, paste("Co-Dependent Genes, Top", N))) + ylab("Correlation") + 
    labs(fill = paste("In", pw_of_interest_name)) +
    theme_minimal() + #scale_fill_nejm() +
    scale_fill_manual(values = c("No" = "gray", "Yes" = "#20854EFF")) +
    theme(axis.text.x=element_text(size=14, angle = 90), 
          axis.title.x=element_text(size=14, face="bold"), 
          axis.title.y=element_text(size=14, face="bold"), axis.text.y=element_text(size=14),
          legend.title = element_text(size=14, face="bold"), 
          legend.text = element_text(size=14), legend.position = "bottom") 
  print(g)
}

# Call function
create_codependency_barplot(kbtbd2_codep, insulin_signaling_genes, 
                            "Insulin Signaling Pathway", 50, "KBTBD2")

# KBTBD2 STRING Interactors conf > 0.4
kbtbd2_string_interactors <- c("BCAS1", "LSM5", "CUL3", "AP5B1", "A1CF", "FAM53C",
                               "SYPL2", "UNCX", "COP1", "FMR1", "SPATA33")
intersect(kbtbd2_codep$Gene, kbtbd2_string_interactors) # no intersection
# KBTBD2 HumanBase Interactors  (top 50)
kbtbd2_hb_interactors <- read.table(paste0(PATH, "HumanBase/kbtbd2_network_interactors.txt"),
                                    header = T, sep = ",")[,'SYMBOL']
intersect(kbtbd2_codep$Gene, kbtbd2_hb_interactors) # only HNRNPM


# Option 2: STRING sub-network for KBTBD2
#' Generate results graph display for a target of interest 
#' This graphical display takes confidence interaction relationships from 
#' STRING, and displays these relationships as a graph. 
#' @param goi the name of the gene of interest
#' @param network a table with the network information from STRING
#' @param conf_thres a confidence threshold for a connection between two 
#' genes in the given network
#' @param top_n if not NA, gives a value of the top N co-dependencies to use
#' @param codep_df a codependency data frame from DepMap for the given goi
create_graphical_representation_target <- function(goi, network, conf_thres, top_n,
                                                   codep_df) {
  if(!is.na(top_n)) {codep_df <- codep_df[1:top_n,]}
  genes <- c(goi, codep_df$Gene)
  
  # Subset the network to the significant target and its interactors
  network_model_targs <- network[(network$node1_name %fin% genes) & 
                                   (network$node2_name %fin% genes),]
  
  # Remove cases where nodes 1 and 2 are the same, or rows where the values in 
  # the row are duplicated in a different order. Divide all confidence scores by 1000 
  network_model_targs <- subset(network_model_targs, node1_name != node2_name & 
                                  !duplicated(cbind(pmin(node1_name, node2_name), 
                                                    pmax(node1_name, node2_name))))
  network_model_targs$combined_score <- network_model_targs$combined_score / 1000

  # For each of the targets, get its a) directionality (correlated or anti-correlated), 
  # b) confidence score in the given network, in relation to the goi
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
  print(head(edge_table))
  
  unique_nodes <- unique(c(edge_table$src, edge_table$target))
  node_table <- as.data.frame(unique_nodes)
  colnames(node_table) <- "nodes"
  node_table$dir <- as.factor(unlist(lapply(unique_nodes, function(n) {
    dir <- NA
    if(n == goi) {dir <- "2"}
    else {
      if(n %fin% codep_df$Gene) {
        dir <- ifelse(as.numeric(codep_df[codep_df$Gene == n, 
                                          'Correlation'])[1] > 0, "1", "-1")
      } else {
        dir <- "0"
      }
    }
    return(dir)
  })))
  print(head(node_table))
  node_table$node_size <- unlist(lapply(node_table$dir, function(d) {
    if(length(d) == 0) {return(NA)}
    if(d == "2") {return(2)}
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
    geom_node_label(aes(label = nodes), repel = T, show.legend = F, 
                    label.size = NA, label.padding = 0.05) +
    theme_graph() +
    scale_fill_manual(values = c("-1" = "#0072B5FF", "2" = "#20854EFF",  
                                 "1" = "#E18727FF"),  #"0" = "lightgray", 
                      labels = c("Neg. Correlation", 
                                 "Pos. Correlation", "Gene-of-Interest"), #"Not in DepMap Co-Dependecy Top 100"
                      name = "Direction of\nCo-Dependency\nCorrelation") +
    theme(legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12), legend.position="bottom") 
  
  # Render the network
  show(driver_targ_network)
  
  return(node_table)
}

string_nw_full <- fread(paste0(PATH, "Validation_Files/string.9606.protein.links.v11.5.namesAdded.txt"))

# Call function
create_graphical_representation_target("KBTBD2", string_nw_full, 0.4, NA, 
                                       kbtbd2_codep)


########################################################################
### SURVIVAL CURVES 
########################################################################
#' Create a survival curve for a given driver gene (and, optionally, expression
#' of a predicted target gene) using TCGA clinical data
#' @param clinical_df a clinical DF from the TCGA for a certain patient set
#' @param driver_name the gene name of the driver of interest
#' @param mutation_count_matrix a mutation count matrix from maftools for TCGA data
#' @param target_name the ENSG ID of a target gene of interest
#' @param expression_df an expression DF (quantile-normalized)
#' @param target_label the gene name of the target gene, for labeling
create_survival_curves <- function(clinical_df, driver_name, mutation_count_matrix,
                                   target_gene, expression_df, target_label) {
  clinical_df_survival <- clinical_df[,colnames(clinical_df) %fin% 
                                        c("case_submitter_id", "vital_status",
                                          "days_to_death", "days_to_last_follow_up")]
  clinical_df_survival <- distinct(clinical_df_survival[which(!is.na(
    clinical_df_survival$vital_status)),])
  #table(clinical_df_survival$vital_status)   #11532 alive, 4182 dead, 14 not reported 
  
  clinical_df_survival$status <- ifelse(clinical_df_survival$vital_status == 'Alive', 
                                        0, 1)
  clinical_df_survival$patient <- unlist(lapply(clinical_df_survival$case_submitter_id, function(x)
    unlist(strsplit(x, "-", fixed = T))[3]))
  
  clinical_df_survival$days <- ifelse(clinical_df_survival$status == 0,
                                      clinical_df_survival$days_to_last_follow_up,
                                      clinical_df_survival$days_to_death)
  
  #print(clinical_df_survival)
  
  # Merge the driver mutation status information with the clinical information
  if(!("patient" %fin% colnames(mutation_count_matrix))) {
    mutation_count_matrix_driver <- melt(mutation_count_matrix[
      mutation_count_matrix$Gene_Symbol == driver_name,])[,c("variable", "value")]
    colnames(mutation_count_matrix_driver) <- c('patient', 'mutation_status')
    #print(head(mutation_count_matrix_driver))
    mutation_count_matrix_driver$patient <- unlist(lapply(mutation_count_matrix_driver$patient, function(x)
      unlist(strsplit(as.character(x), "-", fixed = T))[1]))
    mutation_count_matrix_driver$mutation_status <- unlist(lapply(
      mutation_count_matrix_driver$mutation_status, function(x) ifelse(x > 0, 1, 0)))
  } else {
    mutation_count_matrix_driver <- mutation_count_matrix
  }
  
  survival_df <- merge(mutation_count_matrix_driver, 
                       clinical_df_survival[,c("patient","days","status")],
                       by = "patient")
  
  # Optionally add in target gene expression information as well
  expression_df_target <- NA
  if(!is.na(target_gene)) {
    # Get patients that are have high expression (e.g. top 1/3) or low expression
    # (e.g. bottom 1/3) for the given target gene
    expression_df_target <- as.data.frame(t(expression_df[expression_df$ensg_id == 
                                                            target_gene,]))
    expression_df_target$patient <- rownames(expression_df_target)
    expression_df_target <- expression_df_target[2:nrow(expression_df_target),]
    colnames(expression_df_target)[1] <- "expression"
    expression_df_target <- expression_df_target[order(expression_df_target$expression),]
    onethird <- ceiling(nrow(expression_df_target) / 3)
    bottom_third <- expression_df_target[1:onethird,]
    top_third <- expression_df_target[(onethird*2):nrow(expression_df_target),]
    
    # Add an "expression status" term for patients that fall in these buckets
    expression_df_target$exp_status <- unlist(lapply(expression_df_target$patient, function(id) {
      if(id %fin% top_third$patient) {return(1)}
      else if (id %fin% bottom_third$patient) {return(0)}
      else {return(NA)}
    }))
    expression_df_target <- na.omit(expression_df_target)
    expression_df_target$patient <- unlist(lapply(expression_df_target$patient, function(x)
      unlist(strsplit(as.character(x), "-", fixed = T))[1]))
    print(head(expression_df_target))
  }
  if(!is.na(target_gene)) {
    survival_df <- merge(survival_df, expression_df_target[,c('patient', 'exp_status')],
                         by = "patient", all = F)
  }
  survival_df <- survival_df[survival_df$days != "'--",]
  survival_df$days <- as.numeric(survival_df$days)
  survival_df$treated <- unlist(lapply(survival_df$response, function(x) 
    ifelse(is.na(x), 0, 1))) 
  print(head(survival_df))
  
  # Use the survminer package to create the curves
  fit <- NA
  if(!is.na(target_gene)) {
    fit <- survfit(Surv(time = days, event = status) ~ mutation_status + exp_status, #+ treated, 
                   data = survival_df)
    print(head(fit))
    p <- ggsurvplot(fit, data = survival_df, 
                    palette = c("#0072B5FF", "#FFDC91FF", "#BC3C29FF", "#20854EFF"),
                    pval = T, ggtheme = theme_minimal(), 
                    legend.title = paste(driver_name, 
                                         paste("Mutation Status +\n", 
                                               paste(target_label, "Expression Status"))), 
                    legend = c(0.75,0.75), # "bottom",
                    legend.labs = c("No Mutation, Low Expression", 
                                    "No Mutation, High Expression", 
                                    "Mutation, Low Expression", 
                                    "Mutation, High Expression"), 
                    font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                    font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                    font.legend = c(10, "bold", "black")) + 
      xlab("Time (Days)") + ylab("Survival Probability")
    print(p)
    
  } else {
    fit <- survfit(Surv(time = days, event = status) ~ mutation_status, #+ treated, 
                   data = survival_df)
    p <- ggsurvplot(fit, data = survival_df, 
                    palette = c("#0072B5FF", "#BC3C29FF"),
                    pval = T, ggtheme = theme_minimal(), 
                    legend.title = paste(driver_name, "Mutation Status"), 
                    legend = "bottom", legend.labs = c("No Mutation", "Mutation"), 
                    font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                    font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                    font.legend = c(12, "bold", "black")) + 
      xlab("Time (Days)") + ylab("Survival Probability")
    print(p)
  }
  
  return(survival_df)
}

# Import pan-cancer clinical DF
clinical_df <- read.csv(paste0(PATH, "clinical.csv"), 
                        header = T, check.names = F)
# Import the pan-cancer mutation count matrix
mutation_count_matrix <- read.csv(paste0(PATH, "Mutation/mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv"), 
                                  header = T, check.names = F)
# Import the pan-cancer expression matrix
expression_df_qn <- read.csv(paste0(PATH, "Expression/expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv"), 
                             header = T, check.names = F)

# Call function (e.g. KRAS and NT5E, PIK3CA and KBTBD2)
kras_survival_res <- create_survival_curves(clinical_df, "KRAS", mutation_count_matrix,
                                            "ENSG00000135318", expression_df_qn, "NT5E")
pik3ca_survival_res <- create_survival_curves(clinical_df, "PIK3CA", mutation_count_matrix,
                                              "ENSG00000170852", expression_df_qn, "KBTBD2")


########################################################################
### GROWTH CURVES
########################################################################
mcf7_zero_um_control <- c(84.59252, 92.76691, 96.03666)
mcf7_0.1_um_control <- c(98.93485, 79.68789, 93.73297)
mcf7_1_um_control <- c(61.55561, 64.37949, 64.52811)
mcf7_10_um_control <- c(23.95343, 28.41219, 44.53802)

mcf7_zero_um_sikbtbd2 <- c(33.31682, 32.49938, 31.75625)
mcf7_0.1_um_sikbtbd2 <- c(37.3297, 36.28932, 34.80307)
mcf7_1_um_sikbtbd2 <- c(23.73049, 35.32326, 27.81769)
mcf7_10_um_sikbtbd2 <- c(15.25886, 22.46718, 17.71117)

doses <- c(rep(0, times = 3), rep(0.1, times = 3), 
           rep(1, times = 3), rep(10, times = 3))

input_growthcurve_mcf7 <- data.frame(
  "Label" = c(rep("Control", times = 12), rep("siKBTBD2", times = 12)),
  "Dosage" = rep(doses, times = 2),
  "Growth" = c(mcf7_zero_um_control, mcf7_0.1_um_control, 
               mcf7_1_um_control, mcf7_10_um_control,
               mcf7_zero_um_sikbtbd2, mcf7_0.1_um_sikbtbd2,
               mcf7_1_um_sikbtbd2, mcf7_10_um_sikbtbd2)
)
input_growthcurve_mcf7$Dosage <- as.factor(input_growthcurve_mcf7$Dosage)

mean_df <- input_growthcurve_mcf7 %>% group_by(Label, Dosage) %>% 
  summarize(mean = mean(Growth), se = sd(Growth)/sqrt(3))

ggplot(mean_df, aes(fill = Label, y = mean, x = Dosage)) +
  geom_bar(position="dodge", stat = "identity") + scale_fill_nejm() + 
  theme_minimal() + xlab("Dosage (um) of RLY-2608") + ylab("Relative Cell Growth") +
  labs(col = "Group") + 
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), 
                position = position_dodge(width=0.9), width = 0.2)
