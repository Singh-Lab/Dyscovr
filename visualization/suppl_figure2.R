############################################################
# Code to create Suppl. Figure 2 Visualizations
# Written by Sara Geraghty
# https://www.biorxiv.org/content/10.1101/2024.11.20.624509v1
############################################################

library(data.table)
library(ggplot2)
library(STRINGdb)
library(igraph)
library(GOSemSim)
library(ggrepel)
library(ggsci)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(org.Hs.eg.db)

# Local PATH to directory containing Dyscovr output files
PATH <-  paste0(getwd(), "Output/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


############################################################
### PART A-C: PER-DRIVER NETWORK OVERLAY
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
string_nw <- fread(paste0(PATH, "Validation_Files/9606.protein.links.v12.0.txt"),
                   header = T)
string_nw_info <- read.table(paste0(PATH, "Validation_Files/9606.protein.info.v12.0.txt"), 
                             header = T, sep = "\t")
string_nw_full <- adjust_string_nw_files(string_nw, string_nw_info, 
                                         all_genes_id_conv)
fwrite(string_nw_full, paste0(
  PATH, "Validation_Files/9606.protein.links.v12.0.namesAdded.txt"))

# Read back, once completed
string_nw_full <- fread(paste0(PATH, "Validation_Files/9606.protein.links.v12.0.namesAdded.txt"))

# Call function
create_graphical_representation(pc_allGenes[pc_allGenes$R_i.name == "PIK3CA",], 
                                "PIK3CA", 0.01, string_nw_full, "STRING", 0.4, 100)
create_graphical_representation(pc_allGenes[pc_allGenes$R_i.name == "KRAS",], 
                                "KRAS", 0.01, string_nw_full, "STRING", 0.4, 100)
create_graphical_representation(pc_allGenes[pc_allGenes$R_i.name == "IDH1",], 
                                "IDH1", 0.01, string_nw_full, "STRING", 0.4, 100)
