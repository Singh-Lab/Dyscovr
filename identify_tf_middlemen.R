##############################################################################
# IDENTIFY THE "MIDDLEMAN" TRANSCRIPTION FACTOR
# Created By: Sara Geraghty, Jan. 2023
##############################################################################

# Use DoRothEA to check for enrichment in a given driver's target gene set of 
# the targets of a known TF. Provides an ordered table of candidate TFs that 
# may be a "middleman" linking a mutated driver gene to dysregulated downstream 
# target genes.

# DoRothEA 
# https://bioconductor.org/packages/release/data/experiment/vignettes/dorothea/inst/doc/dorothea.html
# Contains 1333 TFs with targets compiled from curation and various data sources

library(dorothea)
library("STRINGdb")

# Download the DoRothEA network

dorothea_net <- dorothea::dorothea_hs

# Local path
main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/"
main_path <- "C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/"


#' Using the DoRothEA network, creates a table of TF candidates
#' prioritized by their likelihood of being a "middleman" in the interaction
#' between a mutated driver gene and the given ordered list of dysregulated 
#' target genes. Uses both HG and KS tests to determine enrichment, borrowing
#' from function in 'create_enrichment_visualizations.R'.
#' @param dorothea_net the dorothea network of TF/target pairs
#' @param master_df an output DF from my model
#' @param goi a driver gene of interest for enrichment (gene name)
#' @param qval_thres a q-value threshold for HG test enrichment
create_candidate_middlemen_table <- function(dorothea_net, master_df, goi, qval_thres) {
  
  master_df <- master_df[master_df$R_i.name == goi,]
  
  unique_tfs <- unique(dorothea_net$tf)
  
  list_of_partial_dfs <- lapply(unique_tfs, function(tf) {
    tf_targs <- unlist(dorothea_net[dorothea_net$tf == tf, 'target'])
    enrichment_res <- compute_statistical_enrichment(master_df, tf_targs, "both", qval_thres)
    hg_pval <- enrichment_res[[1]]$p.value
    ks_pval <- enrichment_res[[2]]$p.value
    
    return(data.frame("TF" = tf, "Num.TF.Targets" = length(tf_targs),
                      "Hypergeom.Pval" = hg_pval, "K-S.Pval" = ks_pval))
  })
  
  resulting_df <- do.call(rbind, list_of_partial_dfs)
  resulting_df <- resulting_df[order(resulting_df$Hypergeom.Pval),]
}

candidate_middleman_table_pik3ca <- create_candidate_middlemen_table(dorothea_net, 
                                                                     top_drivers_tcga_0.05_pik3ca, 
                                                                     "PIK3CA", 0.01)
candidate_middleman_table_kras <- create_candidate_middlemen_table(dorothea_net, 
                                                                     top_drivers_tcga_0.05_kras, 
                                                                     "KRAS", 0.01)
candidate_middleman_table_idh1 <- create_candidate_middlemen_table(dorothea_net, 
                                                                   top_drivers_tcga_0.05_idh1, 
                                                                   "IDH1", 0.01)

#' Check if these TFs are expressed above a baseline level in the samples with
#' a mutation in the given driver, in order to narrow this list. Add the 
#' absolute difference in expression between mutant/non-mutant samples, as well
#' as the percentage difference
#' @param expression_df an expression DF across all samples of interest
#' @param mutant_patients a vector of patient IDs with a mutation in the given GOI
#' @param candidate_mm_df a data frame of candidate TF middlemen, from
#' 'create_candidate_middlemen_table' function
#' @param all_genes_id_conv ID conversion table from BioMart
add_differential_expression_data <- function(expression_df, mutant_patients,
                                             candidate_mm_df, all_genes_id_conv) {
  
  new_dfs <- lapply(1:nrow(candidate_mm_df), function(i) {
    tf <- candidate_mm_df[i, 'TF']
    print(paste("Now processing", tf))
    tf_ensg <- unique(unlist(all_genes_id_conv[all_genes_id_conv$external_gene_name == tf, 'ensembl_gene_id']))[1]
    
    expr_by_mut_group <- get_expression_by_mut_group(expression_df, mutant_patients, tf_ensg)
    avg_expr_muts <- mean(expr_by_mut_group$mutants, na.rm = T)
    avg_expr_nonmuts <- mean(expr_by_mut_group$normal, na.rm = T)
    
    differential_mut_vs_nonmut <- avg_expr_muts - avg_expr_nonmuts

    perc_diff <- (differential_mut_vs_nonmut / ((avg_expr_muts + avg_expr_nonmuts) / 2)) * 100

    expr_df <- data.frame("Diff.Expr.by.Mut.Status" = differential_mut_vs_nonmut,
                          "Perc.Expr.Change.by.Mut.Status" = perc_diff,
                          "Overall.Avg.Expr.In.Mut.Group" = avg_expr_muts)
    new_df <- cbind(candidate_mm_df[i, ], expr_df)
    
    return(new_df)
  })
  
  full_df <- do.call(rbind, new_dfs)
  
  return(full_df)
}


### HELPER FUNCTIONS FROM 'EXPRESSION_HELPER_FUNCTIONS.R'
get_mutant_patients <- function(mutation_regprot_df, regprot_name) {
  regprot_rows <- which(unlist(lapply(mutation_regprot_df$Query, function(x) 
    unlist(strsplit(unlist(strsplit(x, "|", fixed = TRUE))[3], "_", fixed = TRUE))[1] == regprot_name)))
  mutant_patients <- unique(unlist(strsplit(mutation_regprot_df[regprot_rows, "Patient"], ";", fixed = TRUE)))
  return(mutant_patients)
}

get_mutant_patients2 <- function(mutation_count_matrix, regprot_name) {
  regprot_row <- mutation_count_matrix[mutation_count_matrix$Gene_Symbol == regprot_name, ]
  ind <- which(as.integer(regprot_row) > 0)
  return(colnames(regprot_row)[ind])
}


get_expression_by_mut_group <- function(express_df_cancer, mutant_patients, ensg) {
  #colnames(express_df_cancer) <- unlist(lapply(colnames(express_df_cancer), function(x) 
    #unlist(strsplit(x, "-", fixed = TRUE))[1]))
  #express_df_cancer <- express_df_cancer[,2:ncol(express_df_cancer)]
  expression_mutants <- as.numeric(express_df_cancer[express_df_cancer$ensg_id == ensg, colnames(express_df_cancer) %fin% mutant_patients])
  expression_normal <- as.numeric(express_df_cancer[express_df_cancer$ensg_id == ensg, !(colnames(express_df_cancer) %fin% mutant_patients)])
  return(list("mutants" = expression_mutants, "normal" = expression_normal))
}


# Import expression DF (TMM normalization is all positive values, which is easier to interpret)
expression_df_quantile_norm <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
expression_df_tmm_norm <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Expression/expression_tmm_filtByExpr_SDGr1_CancerOnly_IntersectPatients.csv", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
expression_df_tmm_norm <- expression_df_tmm_norm[,!(colnames(expression_df_tmm_norm) == "ensg_id.1")]

# Import mutation regprot DF
mutation_regprot_df <- read.csv(paste(main_path, "Linear Model/Tumor_Only/Regprot_Mutation/iprotein_results_nonsynonymous_IntersectPatients.csv", sep = ""), 
                                header = TRUE, row.names = 1, check.names = FALSE)
mutation_count_df <- read.csv(paste(main_path, "Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv", sep = ""), 
                              header = TRUE, check.names = FALSE)


#pik3ca_mutant_patients <- get_mutant_patients(mutation_regprot_df, "PK3CA")
#cdh1_mutant_patients <- get_mutant_patients(mutation_regprot_df, "CADH1")

pik3ca_mutant_patients <- get_mutant_patients2(mutation_count_df, "PIK3CA")
#cdh1_mutant_patients <- get_mutant_patients2(mutation_count_df, "CDH1")
kras_mutant_patients <- get_mutant_patients2(mutation_count_df, "KRAS")
idh1_mutant_patients <- get_mutant_patients2(mutation_count_df, "IDH1")

# Call function
candidate_middleman_table_pik3ca <- add_differential_expression_data(expression_df_quantile_norm, pik3ca_mutant_patients,
                                                                     candidate_middleman_table_pik3ca, all_genes_id_conv)
candidate_middleman_table_cdh1 <- add_differential_expression_data(expression_df_quantile_norm, cdh1_mutant_patients,
                                                                     candidate_middleman_table_cdh1, all_genes_id_conv)
candidate_middleman_table_kras <- add_differential_expression_data(expression_df_quantile_norm, kras_mutant_patients,
                                                                   candidate_middleman_table_kras, all_genes_id_conv)

#' Phosphorylation can either activate or deactivate a TF. If the TF is activated, we
#' expect to see its targets (which overlap the drivers's hits) be upregulated. If the TF
#' is deactivated, we expect to see its targets (which overlap the driver's hits) be
#' downregulated. This function gets number of the TF's targets that both a) overlap
#' the drivers significant hits, and b) are either up- or down-regulated according to 
#' the model. 
#' @param master_df an output DF from my model
#' @param candidate_mm_df a data frame of candidate TF middlemen, from
#' 'create_candidate_middlemen_table' function
#' @param dorothea_net the dorothea network of TF/target pairs
#' @param goi a driver gene of interest for enrichment (gene name)
#' @param qval_thres q-value threshold for determining the GOI's significant hits
get_up_and_down_regulated_targets <- function(master_df, candidate_mm_df, 
                                              dorothea_net, goi, qval_thres) {
  
  master_df <- master_df[master_df$R_i.name == goi,]
  sig_hits <- master_df[master_df$q.value < qval_thres, 'T_k.name']
  print(head(sig_hits))
  
  # Instantiate the string DB object 
  string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=0, 
                            input_directory="")
  
  new_dfs <- lapply(1:nrow(candidate_mm_df), function(i) {
    tf <- candidate_mm_df[i, 'TF']
    print(paste("Now processing", tf))
    
    tf_targs <- unlist(dorothea_net[dorothea_net$tf == tf, 'target'])
    intersection <- intersect(as.character(unlist(sig_hits)), tf_targs)
    print(length(intersection))
    
    if(length(intersection) > 0) {
      intersection_betas <- master_df[master_df$T_k.name %fin% intersection, 'estimate']
      num_upregulated <- length(intersection_betas[intersection_betas > 0])
      num_downregulated <- length(intersection_betas[intersection_betas < 0])
      
      intersection_names <- master_df[master_df$T_k.name %fin% intersection, 'T_k.name']
      stringdb_enrichment <- string_db$get_enrichment(intersection_names)
      enriched_pws <- paste(stringdb_enrichment[stringdb_enrichment$category == "RCTM", 'description'], 
                            collpase = ";")
      
      up_down_df <- data.frame("Num.Upreg.Overlapping.Targs" = num_upregulated,
                               "Num.Downreg.Overlapping.Targs" = num_downregulated,
                               "Enriched.RCTM.PWs" = enriched_pws)
      new_df <- cbind(distinct(candidate_mm_df[i, ]), up_down_df)
      
    } else {
      new_df <- cbind(candidate_mm_df[i, ], data.frame("Num.Upreg.Overlapping.Targs" = NA,
                                 "Num.Downreg.Overlapping.Targs" = NA, "Enriched.RCTM.PWs" = NA))
    }
    return(new_df)
  })
  
  full_df <- do.call(rbind, new_dfs)
  full_df$K.S.Qval <- qvalue(full_df$K.S.Pval)$qvalues
  full_df <- full_df[order(full_df$K.S.Qval),]
  
  return(full_df)
}

candidate_middleman_table_pik3ca <- get_up_and_down_regulated_targets(top_drivers_tcga_0.05_pik3ca, candidate_middleman_table_pik3ca,
                                                                     dorothea_net, "PIK3CA", 0.2)
candidate_middleman_table_cdh1 <- get_up_and_down_regulated_targets(top_drivers_tcga_0.05_cdh1, candidate_middleman_table_cdh1, 
                                                                    dorothea_net, "CDH1", 0.2)

write.csv(candidate_middleman_table_pik3ca, paste0(main_path, "Linear Model/Middleman_TF/pik3ca_middlemen_q0.1_qn.csv"))



# Get the actual intersecting genes for particular cases of interest
pik3ca_top_hits <- top_drivers_tcga_0.05_pik3ca[top_drivers_tcga_0.05_pik3ca$q.value < 0.2, 'T_k.name']

myt1l_targets <- unlist(dorothea_net[dorothea_net$tf == 'MYT1L', 'target'])
myt_targets <- unlist(dorothea_net[dorothea_net$tf == 'MYT', 'target'])
st18_targets <- unlist(dorothea_net[dorothea_net$tf == 'ST18', 'target'])
srebf2_targets <- unlist(dorothea_net[dorothea_net$tf == 'SREBF2', 'target'])

myt1l_intersect <- intersect(pik3ca_top_hits, myt1l_targets)
myt_intersect <- intersect(pik3ca_top_hits, myt_targets)
st18_intersect <- intersect(pik3ca_top_hits, st18_targets)
srebf2_intersect <- intersect(pik3ca_top_hits, srebf2_targets)

# CDH1
cdh1_top_hits <- top_drivers_tcga_0.05_cdh1[top_drivers_tcga_0.05_cdh1$q.value < 0.2, 'T_k.name']

print(paste(CREBZF_intersect, collapse = ", "))
