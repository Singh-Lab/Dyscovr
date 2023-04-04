###########################################################################################################
# PROCESS STRING DATA
# Written by: Sara Geraghty, April 2022
###########################################################################################################

# Script to process String database interactions for a given gene of interest and extract
# network interactions for gene set enrichment analysis

###########################################################################################################
# OPTION 1: USE STRING NETWORK FROM WEB INTERFACE
# Download String networks from website: https://string-db.org/cgi/network?taskId=bSXrcGoHJTto&sessionId=beGJRBX9wxZ4
###########################################################################################################
# Import the String network file
string_nw <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/STRING/string_interactions.tsv",
                      sep = "\t", header = TRUE, check.names = FALSE)
colnames(string_nw)[1] <- "node1"

# Limit to given GOI
goi <- "TP53"
# goi <- "PIK3CA"
string_nw_goi <- string_nw[(string_nw$node1 == goi),]  # will already be subsetted by our max # of interactions

# Visualize the scores for this GOI
string_nw_goi <- string_nw_goi[order(string_nw_goi$combined_score),]
hist(string_nw_goi$combined_score)

# Get the list of genes that made the cutoff for top interactions with our GOI
interactors <- string_nw_goi$node2

# Write these interactors to a file
write.table(interactors, paste0(getwd(), paste0("/Desktop/", paste0(goi, "_interactors_string_top100.txt"))),
            quote = FALSE, row.names = FALSE)


###########################################################################################################
# OPTION 2: USE STRINGDB R PACKAGE
###########################################################################################################
BiocManager::install("STRINGdb")
BiocManager::install("EnsDb.Hsapiens.v86")
#^ This package appears to do the best at gene ID conversion, according to this source: https://shiring.github.io/genome/2016/10/23/AnnotationDbi

library("STRINGdb")
library(igraph)
library("EnsDb.Hsapiens.v86")
# keytypes(EnsDb.Hsapiens.v86)  # list all supported keytypes

string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=0, 
                          input_directory="", protocol = "http")

# NOTE: STRINGdb also offers a way to compute enrichment in a given gene set in GO, KEGG, REACTOME, 
# PubMed, UniProt Keywords, and PFAM/INTROPRO/SMART domains
# enrichment <- string_db$get_enrichment( hits ); head(enrichment, n = 20)

tp53 = string_db$mp( "tp53" )
pik3ca = string_db$mp( "pik3ca" )

# Get the neighbors for a given gene
tp53_neighbors <- string_db$get_neighbors(tp53)
pik3ca_neighbors <- string_db$get_neighbors(pik3ca)

# Convert these from ENSP to ENSG
tp53_neighbors_ensg <- ensembldb::select(EnsDb.Hsapiens.v86, 
                                         keys=unlist(lapply(tp53_neighbors, function(x) 
                                           unlist(strsplit(x, ".", fixed = TRUE))[2])), 
                                         keytype = "PROTEINID", columns = c("SYMBOL","PROTEINID","GENEID"))
pik3ca_neighbors_ensg <- ensembldb::select(EnsDb.Hsapiens.v86, 
                                         keys=unlist(lapply(pik3ca_neighbors, function(x) 
                                           unlist(strsplit(x, ".", fixed = TRUE))[2])), 
                                         keytype = "PROTEINID", columns = c("SYMBOL","PROTEINID","GENEID"))


# Get the interactions for a given gene
tp53_interactions <- string_db$get_interactions(c(tp53, tp53_neighbors))
pik3ca_interactions <- string_db$get_interactions(c(pik3ca, pik3ca_neighbors))


#' Given an 'interactions' DF from STRINGDb, a gene-of-interest, and a threshold
#' for the combined score from STRING (can be 0), subsets the interactions DF
#' to only those distinct interactions exceeding the score and adjusts the ENSP
#' IDs to ENSG and gene symbol IDs
#' @param interactions an interactions DF from STRINGDb
#' @param goi a gene-of-interest symbol
#' @param combined_score_thres a threshold for the combined score from STRING
adjust_interaction_df <- function(interactions, goi, combined_score_thres) {
  
  # Get only distinct interactions (eliminate duplicate rows)
  interactions <- distinct(interactions)
  
  # Fix the ENSP labels
  interactions$from <- unlist(lapply(interactions$from, function(x) 
    unlist(strsplit(x, ".", fixed = TRUE))[2]))
  interactions$to <- unlist(lapply(interactions$to, function(x) 
    unlist(strsplit(x, ".", fixed = TRUE))[2]))
  
  # Get all the ENSP variants for our GOI
  conv <- ensembldb::select(EnsDb.Hsapiens.v86, keys=c(goi), keytype = "SYMBOL", 
                            columns = c("GENEID","PROTEINID"))
  protein_ids <- unique(conv$PROTEINID)
  protein_ids <- protein_ids[!is.na(protein_ids)]
  
  # Select only interactions that our GOI is a part of 
  interactions_goi <- interactions[(interactions$from %in% protein_ids) | 
                                     (interactions$to %in% protein_ids), ]
  
  # Select only interactions above a given combined score threshold
  interactions_goi_sub <- interactions_goi[interactions_goi$combined_score > combined_score_thres, ]
  
  # Extract this info and switch ENSP IDs to ENSG IDs/ gene symbols
  interactions_geneNames_to <- ensembldb::select(EnsDb.Hsapiens.v86, 
                                              keys=interactions_goi_sub$to, 
                                              keytype = "PROTEINID", columns = c("SYMBOL","GENEID"))
  interactions_geneNames_from <- ensembldb::select(EnsDb.Hsapiens.v86, 
                                                 keys=interactions_goi_sub$from, 
                                                 keytype = "PROTEINID", columns = c("SYMBOL","GENEID"))
  interactions_geneNames_to <- interactions_geneNames_to[, c("SYMBOL", "GENEID")]
  interactions_geneNames_from <- interactions_geneNames_from[, c("SYMBOL", "GENEID")]
  interactions_geneNames <- distinct(rbind(interactions_geneNames_to, interactions_geneNames_from))
  #colnames(interactions_geneNames_to) <- paste("_TO", colnames(interactions_geneNames_to))
  #colnames(interactions_geneNames_from) <- paste("_FROM", colnames(interactions_geneNames_from))
  #interactions_geneNames <- cbind(cbind(interactions_geneNames_to, interactions_geneNames_from), 
                                  #interactions_goi_sub$combined_score)
  
  # Return this DF of interactions
  return(interactions_geneNames)
}

tp53_interactions_noThres <- adjust_interaction_df(tp53_interactions, "TP53", 0)
tp53_interactions_200 <- adjust_interaction_df(tp53_interactions, "TP53", 200)

pik3ca_interactions_200 <- adjust_interaction_df(pik3ca_interactions, "PIK3CA", 200)


# Alternatively, can play directly with the STRING igraph object
igraph <- string_db$get_graph()

# Determine basic qualities of the graph
is_directed(igraph)   # FALSE
is_weighted(igraph)   # FALSE

# Determine attributes of vertices and edges
unique(names(vertex_attr(igraph))) # 'name', in the form of '9606.ENSP00000XXXXXX', where 9696 is the human species ID
unique(names(edge_attr(igraph))) # 'combined_score', in the form of an integer value

# Convert the combined scores to be weights. Notes that the combined score is the confidence score (e.g. 0.7) multiplied by 1000
# to make it an integer (0.7 = 700). Thus, the larger the score, the higher the confidence. In order to find shortest paths, we must
# invert these values (1000-X). This would change a low score of 0.3 to 0.7, and a high score of 0.8 to 0.2, making it easier to use built-in
# functions that consider weight the cost of the path
E(igraph)$weight <- 1000 - edge_attr(igraph)$combined_score
is_weighted(igraph)   # TRUE
head(E(igraph)$weight)

# Get all genes within a certain fixed distance (e.g. combined score of 350) from our GOI
#tp53_paths <- all_simple_paths(igraph, from = "9606.ENSP00000269305", to = V(igraph), cutoff = 2)   # weights will automatically be used

#' A simple algorithm to collect and record neighbors for the given origin, with set distance 
#' or weight maximum.
#' @param igraph the igraph object
#' @param origin the name of the origin node
#' @param dist_cutoff (optional) a maximum distance away from origin node
#' @param weight_cutoff (optional) a maximum total weigh from the origin node
get_neighborhood <- function(igraph, origin, dist_cutoff, weight_cutoff) {
  goi_subgraph <- igraph
  
  # If there is a distance cutoff, start by getting all neighbors within the cutoff
  if(!is.na(dist_cutoff)) {
    goi_subgraph <- make_ego_graph(goi_subgraph, order = dist_cutoff, nodes = origin)[[1]]
    closest_neighbors <- unlist(lapply(V(goi_subgraph)$name, function(v) 
      unlist(strsplit(v, ".", fixed = TRUE))[2]))
  }
  
  # Then, get all neighbors within the weight cutoff, if it exists
  if(!is.na(weight_cutoff)) {
    #goi_subgraph <- graph.neighborhood(subgraph.edges(igraph, eids = which(E(igraph)$weight < weight_cutoff)), 
                                           #order = 1)
    #res_matrix <- matrix(nrow = 1, ncol = length(unlist(V(goi_subgraph))))
    nonorigin_verts <- setdiff(V(goi_subgraph)$name, origin)
    distances <- distances(goi_subgraph, v = origin, to = nonorigin_verts, 
                           weights = E(goi_subgraph)$weight, algorithm = c("dijkstra"))
    
    # Limit to vertices whose shortest path distance is less than the given weight cutoff
    distances <- as.data.frame(distances)
    print(range(distances[1,]))
    distances_sub <- distances[, which(as.numeric(distances[1,]) < weight_cutoff)]
    
    # Identify these neighbors
    closest_neighbors <- unlist(lapply(colnames(distances_sub), function(x) 
      unlist(strsplit(x, ".", fixed = TRUE))[2]))
  }
  
  # Convert closest neighbors to ENSG/ Symbols
  res <- ensembldb::select(EnsDb.Hsapiens.v86, keys = closest_neighbors, keytype = "PROTEINID", 
                            columns = c("GENEID", "SYMBOL"))

  #head(goi_subgraph)
  #plot(goi_subgraph[[1]], edge.label = E(goi_subgraph[[1]]$weight))
  
  #return(goi_subgraph)
  return(res)
}


tp53_neighborhood <- get_neighborhood(igraph, origin = "9606.ENSP00000269305", 
                                      dist_cutoff = 2, weight_cutoff = 5)
pik3ca_neighborhood <- get_neighborhood(igraph, origin = "9606.ENSP00000269305", 
                                      dist_cutoff = 2, weight_cutoff = 5)

###########################################################################################################
###########################################################################################################
### FLOW THROUGH METABOLIC SUBNETWORKS
###########################################################################################################
###########################################################################################################
# Useful tutorial: https://rpubs.com/HWH/913747

all_genes_id_conv <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/all_genes_id_conv.csv", header = TRUE)

# Restrict the network to just the genes in Recon3D and our driver genes of interest
metabolic_genes <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/metabolic_targets.csv", 
                            header = TRUE, check.names = FALSE)
metabolic_genes_uniprot <- unlist(lapply(metabolic_genes$swissprot, function(s) 
  unlist(strsplit(s, ";", fixed = TRUE))))
driver_genes <- unique(master_df$R_i)

genes_of_interest <- unique(c(driver_genes, metabolic_genes_uniprot))

metabolic_genes_names <- unlist(lapply(metabolic_genes_uniprot, function(x) 
  unique(all_genes_id_conv[all_genes_id_conv$uniprot_gn_id == x, 'external_gene_name'])))
driver_genes <- unique(master_df$R_i.name)
genes_of_interest_names <- unique(c(metabolic_genes_names, driver_genes))

###########################################################################################################
# SET UP STRING NETWORK
###########################################################################################################
string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=0, 
                          input_directory="", protocol = "http")

# Import the Uniprot 2 STRING conversion file to get STRING IDs
# Conversion files from STRING can be found here: https://string-db.org/mapping_files/
id_map <- read.table("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/STRING/human.uniprot_2_string.2018.tsv",
                     sep = "\t", header = F)
colnames(id_map) <- c("species_id", "uniprot_gn_id", "ensp_id", "unknown", "string_id")
genes_of_interest_stringIDs <- unlist(lapply(genes_of_interest, function(g) 
  id_map[grepl(g, id_map$uniprot_gn_id), 'string_id']))

# Generate a subnetwork using these IDs
string_metabolic_subnw <- string_db$get_subnetwork(genes_of_interest_stringIDs)

# Visualize this network
string_db$plot_network(genes_of_interest_stringIDs, required_score = 400)
subnw_summary <- string_db$get_summary(genes_of_interest_stringIDs, required_score = 400)
# SUMMARY: 1096 proteins, 6256 interactions, 4019 expected interactions (p-value = 0)

# ALTERNATIVE: USE STRING TABLE 
string_nw_full <- fread("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Input Data Files/Pan-Cancer/Network_Data/STRING/string.9606.protein.links.v11.5.namesAdded.txt")

string_nw_metabolism <- string_nw_full[(string_nw_full$node1_name %in% genes_of_interest_names) |
                                         (string_nw_full$node2_name %in% genes_of_interest_names),]
string_nw_metabolism_conf <- string_nw_metabolism[string_nw_metabolism$combined_score > 400,]

# Remove entries that are identical in both directions
string_nw_metabolism_conf <- subset(string_nw_metabolism_conf, node1_name != node2_name & 
                                !duplicated(cbind(pmin(node1_name, node2_name), 
                                                  pmax(node1_name, node2_name))))
string_nw_metabolism_conf_names <- string_nw_metabolism_conf[, c("node1_name", "node2_name")]

# Infer directionality in some way? Or assume bidirectionality, after signal leaves source node (undirected network)

# Reform into a table that can be used as input to igraph
string_nw_metabolism_conf_forigraph <- data.frame("from" = string_nw_metabolism_conf$node1_name,
                                                  "to" = string_nw_metabolism_conf$node2_name,
                                                  "conf" = string_nw_metabolism_conf$combined_score)
graph <- igraph::graph_from_data_frame(string_nw_metabolism_conf_forigraph, directed = FALSE,
                                       vertices = string_nw_metabolism_conf_names)
print(graph, e=TRUE, v=TRUE)


###########################################################################################################
# SET UP HUMANBASE NETWORK
###########################################################################################################
# Mammary Epithelium
# Link to network download: https://hb.flatironinstitute.org/download ('Top edges')

# Already processed to include gene names and be above confidence of 0.4
hb_mammary_epithelium_conf <- read.csv(paste0(nw_path, "HumanBase/mammary_epithelium_0.4conf.csv"), 
                                       header = T, check.names = F)

# Restrict to metabolic genes and drivers
hb_mammary_epi_metabolism <- hb_mammary_epithelium_conf[(hb_mammary_epithelium_conf$node1_name %in% genes_of_interest_names) |
                                         (hb_mammary_epithelium_conf$node2_name %in% genes_of_interest_names),]

# Remove entries that are identical in both directions
hb_mammary_epi_metabolism <- subset(hb_mammary_epi_metabolism, node1_name != node2_name & 
                                      !duplicated(cbind(pmin(node1_name, node2_name), 
                                                        pmax(node1_name, node2_name))))
hb_mammary_epi_metabolism_names <- hb_mammary_epi_metabolism[, c("node1_name", "node2_name")]

# Infer directionality in some way? Or assume bidirectionality, after signal leaves source node (undirected network)

# Reform into a table that can be used as input to igraph
hb_mammary_epi_metabolism_forigraph <- data.frame("from" = hb_mammary_epi_metabolism$node1_name,
                                                  "to" = hb_mammary_epi_metabolism$node2_name,
                                                  "conf" = hb_mammary_epi_metabolism$combined_score)
graph <- igraph::graph_from_data_frame(hb_mammary_epi_metabolism, directed = FALSE,
                                       vertices = hb_mammary_epi_metabolism_names)
print(graph, e=TRUE, v=TRUE)


###########################################################################################################
# 
###########################################################################################################

