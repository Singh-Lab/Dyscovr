#############################################################################
### SYNTHETIC LETHALITY ANALYSIS PAN-CANCER, DEPMAP
### Written by Sara Geraghty, Princeton University
### https://www.biorxiv.org/content/10.1101/2024.11.20.624509v1
#############################################################################

# Run pan-cancer regression synthetic lethality analysis across multiple drivers 
# on the cluster using SLURM job submission 

#!/usr/bin/env Rscript

# Set the library path
library.path <- .libPaths()

library(data.table, lib.loc = library.path, quietly = T)
library(broom, lib.loc = library.path, quietly = T)
library(dplyr, lib.loc = library.path, quietly = T)
library(qvalue, lib.loc = library.path, quietly = T)
library(metap, lib.loc = library.path, quietly = T)
library(fastmatch, lib.loc = library.path, quietly = T)
library(car, lib.loc = library.path, quietly = T)
library(multcomp, lib.loc = library.path, quietly = T)
library(EmpiricalBrownsMethod, lib.loc = library.path, quietly = T)

# Set input paths (PATH_IN is location to Dyscovr output files, PATH_DEPMAP is location to DepMap input files)
PATH_IN <- "/Genomics/argo/users/scamilli/Dyscovr_Personal/output_files/PanCancer/top_0.05_orig/"
PATH_DEPMAP <- "/Genomics/argo/users/scamilli/Dyscovr_Personal/input_files/PanCancer/DepMap/"

# Set output path where synthetic lethality output will be written
PATH_OUT <- paste0(PATH_IN, "synthetic_lethal_output/")

`%fin%` <- function(x, table) {
  #stopifnot(require(fastmatch))
  fmatch(x, table, nomatch = 0L) > 0L
}

############################################################
### IMPORT PAN- and PER-CANCER OUTPUT FILE(S)
############################################################
# Pan-Cancer
outfn <- "res_top_0.05_allGenes_quantile_rawCNA_methMRaw_3PCs_Nonsyn.Drivers.Vogel.elim.vif5.sp0.7_corrected_MUT.csv"
pc_allGenes <- read.csv(paste0(PATH_IN, outfn), header = T, check.names = F, 
                        row.names = 1)

# Per-Cancer
perCancer_fns <- intersect(list.files(path = PATH_IN, recursive = T,
                                      pattern = "_corrected_MUT"), 
                           intersect(list.files(path = PATH_IN, recursive = T,
                                                pattern = "allGenes_"),
                                     list.files(path = PATH_IN, recursive = T,
                                                pattern = "Nonsyn.Drivers.Vogel.elim.vif")))
perCancer <- lapply(perCancer_fns, function(f) 
  fread(paste0(PATH_IN, f), header = T))
names(perCancer) <- unlist(lapply(perCancer_fns, function(x)
  unlist(strsplit(x, "/", fixed = T))[1]))

# Ensure we only keep the 19 cancer types we are actively considering
cts <- c("ACC", "BLCA", "BRCA", "CESC", "COAD", "ESCA", "HNSC", "KIRC", "KIRP",
         "LGG", "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "PCPG", "PRAD", "THCA", "UCEC")
perCancer <- perCancer[names(perCancer) %fin% cts]


############################################################
### IMPORT DEPMAP FILES
############################################################
# DepMap files can be downloaded directly from their website using custom 
# (https://depmap.org/portal/download/custom/) or full (https://depmap.org/portal/download/all/)
# download webpages

# Import CRISPRi, expression, and mutation data
crispr <- read.csv(paste0(PATH_DEPMAP, "CRISPR_(DepMap_Public_23Q2+Score,_Chronos).csv"), 
                   header = T, check.names = F, row.names = 1)

# Quantile-normalized expression (call quantile_normalize_expression.py on the downloaded expression matrix)
#expression <- read.csv(paste0(PATH_DEPMAP, "Expression_Public_23Q2.csv"), 
#                                        header = T, check.names = F)
expression_qn <- read.csv(paste0(PATH_DEPMAP, "expression_quantile_normalized_sklearn.csv"), 
                          header = T, check.names = F, row.names = 1)

mutations <- read.csv(paste0(PATH_DEPMAP, "OmicsSomaticMutations.csv"),
                      header = T, check.names = F)

# Get mutations specifically for driver genes of interest in the analysis
types <- c("MISSENSE", "NONSENSE", "SPLICE_SITE", "NONSTOP")
mutations_idh1 <- mutations[(mutations$HugoSymbol == "IDH1") & 
                              (mutations$VariantInfo %fin% types),]
mutations_kras <- mutations[(mutations$HugoSymbol == "KRAS") & 
                              (mutations$VariantInfo %fin% types),]
mutations_tp53 <- mutations[(mutations$HugoSymbol == "TP53") & 
                              (mutations$VariantInfo %fin% types),]
mutations_pik3ca <- mutations[(mutations$HugoSymbol == "PIK3CA") & 
                                (mutations$VariantInfo %fin% types),]

# Cell line metadata information
cell_line_sample_info <- read.csv(paste0(PATH_DEPMAP, "cell_line_sample_info.csv"),
                                  header = T, check.names = F)

# Create a mapping between TCGA cancer type and cell line primary disease
cell_line_cancer_type_mapping <- list(
  "ACC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Adrenal Cancer", 'DepMap_ID'],
  "BLCA" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Bladder Cancer", 'DepMap_ID'],
  "BRCA" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Breast Cancer", 'DepMap_ID'],
  "CESC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Cervical Cancer", 'DepMap_ID'],
  "COAD" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Colon/Colorectal Cancer", 'DepMap_ID'],
  "ESCA" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Esophageal Cancer", 'DepMap_ID'],
  "HNSC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Head and Neck Cancer", 'DepMap_ID'],
  "KICH" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Kidney Cancer", 'DepMap_ID'],
  "KIRC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Kidney Cancer", 'DepMap_ID'],
  "KIRP" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Kidney Cancer", 'DepMap_ID'],
  "LGG" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Brain Cancer", 'DepMap_ID'],
  "LIHC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Liver Cancer", 'DepMap_ID'],
  "LUAD" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Lung Cancer", 'DepMap_ID'],
  "LUSC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Lung Cancer", 'DepMap_ID'],
  "MESO" = NA,
  "PAAD" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Pancreatic Cancer", 'DepMap_ID'],
  "PRAD" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Prostate Cancer", 'DepMap_ID'],
  "SARC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Sarcoma", 'DepMap_ID'],
  "STAD" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Gastric Cancer", 'DepMap_ID'],
  "THCA" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Thyroid Cancer", 'DepMap_ID'],
  "UCEC" = cell_line_sample_info[cell_line_sample_info$primary_disease == "Endometrial/Uterine Cancer", 'DepMap_ID']
)

############################################################
### PREPROCESS DEPMAP FILES
############################################################
#' Adjust the knockout and expression data frames to merge them with the mutation
#' data and "melt" them so that we can more easily make plots
#' @param depmap_df a DepMap knockout or expression data frame
#' @param mutation_df a DepMap mutation data frame
#' @param genes_of_interest a vector of the Hugo IDs of genes whose mutation
#' status is of interest
adjust_depmap_df <- function(depmap_df, mutation_df, genes_of_interest) {
  
  # Melt data frame to get it ready for plots
  depmap_df_melt <- melt(as.matrix(depmap_df))
  #print(head(depmap_df_melt))
  
  # Adjust column names
  depmap_id_col <- ifelse(grepl("ACH", depmap_df_melt[1,1]), 1, 2)
  colnames(depmap_df_melt)[depmap_id_col] <- "depmap_id"
  gene_id_col <- ifelse(depmap_id_col == 1, 2, 1)
  colnames(depmap_df_melt)[gene_id_col] <- "gene"
  #print(head(depmap_df_melt))
  
  # Make the "value" column numeric, the "gene" column a character
  depmap_df_melt$value <- as.numeric(unlist(depmap_df_melt$value))
  depmap_df_melt$gene <- as.character(unlist(depmap_df_melt$gene))
  
  # Add mutation status of genes of interest 
  unique_cls <- unique(intersect(mutation_df$ModelID, depmap_df_melt$depmap_id))
  #print(length(unique_cls))
  
  # Get the cell lines that have a mutation
  unique_drivers <- unique(mutation_df$HugoSymbol)
  mutation_df_driver_cols <- as.data.frame(lapply(unique_drivers, function(driver) {
    cls <- mutation_df[mutation_df$HugoSymbol == driver, 'ModelID']
    return(as.factor(unlist(lapply(unique_cls, function(cl) 
      ifelse(cl %fin% cls, 1, 0)))))
  }))
  colnames(mutation_df_driver_cols) <- unique_drivers
  mutation_df_new <- cbind(data.frame("depmap_id" = unique_cls), 
                           mutation_df_driver_cols)
  #print(head(mutation_df_new))
  
  depmap_df_melt <- merge(depmap_df_melt, mutation_df_new, by = "depmap_id")
  #depmap_df_melt[is.na(depmap_df_melt)] <- 0
    
  return(depmap_df_melt)
}

mutations_drivers <- do.call(rbind, list(mutations_idh1, mutations_kras, 
                                         mutations_pik3ca, mutations_tp53))
crispr <- adjust_depmap_df(crispr, mutations_drivers,  
                           c("TP53", "PIK3CA", "KRAS", "IDH1"))
expression_qn <- adjust_depmap_df(expression_qn, mutations_drivers, 
                                  c("TP53", "PIK3CA", "KRAS", "IDH1"))

print("Num unique cell lines with all data:")
print(length(intersect(crispr$depmap_id, expression_qn$depmap_id)))


############################################################
### PRIMARY SYNTHETIC LETHALITY FUNCTIONS
############################################################
#' Explore the relationship between the expression of a driver gene (mutated or
#' unmutated) and cell viability when each of a given set of targets is knocked
#' out via CRISPRi using a simple linear regression model across all cancer types
#' with a one-hot encoded term for cancer type
#' @param crispr a DepMap data frame with CRISPRi knockout data
#' @param expression a DepMap data frame with gene expression data
#' @param driver_mut the name of the driver gene of interest
#' @param perCancer the named list of Dyscovr output files, limited to
#' cancer types containing the provided driver gene
#' @param qval_thres a q-value threshold
#' @param use_pc_hits a T/F value indicating whether we are limiting to targets
#' that are statistically significant pan-cancer
#' @param pc_hits if restricting to PC hits, a vector of PC hits
#' @param use_drugbank a T/F value indicating whether we are limiting to targets
#' with a drug in DrugBank'
#' @param drugbank_gns if restricting by DrugBank, a vector of gene names with 
#' known associated drugs from the database
#' @param use_dependency_check a T/F value indicating whether we are checking for
#' differential dependency by mutation status prior to looking for a correlation
#' between unmutated expression of driver and target dependency
#' @param cell_line_cancer_type_mapping a mapping between depmap IDs and cancer type
#' @param terms_of_interest a vector of terms we'd like to include as variables
#' (options include 'mutation' and 'expression')
#' @param recombine either "None", "linearHypothesis", or "multcomp" to indicate 
#' whether or not we are recombining the driver mutation and expression status using 
#' the car package's linearHypothesis function or the multcomp package's multcomp
#' function
#' @param thres_cls an integer denoting the minimum number of cell lines in a 
#' given cancer type that do not have a mutation in the given driver
perform_lm_analysis_crosscancers <- function(crispr, expression, driver_mut, 
                                             perCancer, qval_thres, use_pc_hits, 
                                             pc_hits, use_drugbank, drugbank_gns, 
                                             use_dependency_check, 
                                             cell_line_cancer_type_mapping,
                                             terms_of_interest, recombine,
                                             thres_cls = 15) {
  
  tophits <- unique(unlist(lapply(perCancer, function(master_df) {
    hits <- as.character(unlist(master_df[(master_df$q.value < qval_thres) & 
                                            (master_df$R_i.name == driver_mut),
                                          'T_k.name']))
    return(hits)
  })))
  
  # Limit CRISPR and expression DepMap data to just cancer types in which 
  # driver is mutated
  cts <- names(perCancer)
  depmap_ids <- unlist(lapply(1:length(cell_line_cancer_type_mapping), function(i) {
    if(names(cell_line_cancer_type_mapping)[i] %fin% cts) {
      return(unlist(cell_line_cancer_type_mapping[[i]]))
    } else {return(NA)}
  }))
  depmap_ids <- depmap_ids[!is.na(depmap_ids)]
  
  crispr_sub <- crispr[crispr$depmap_id %fin% depmap_ids,]
  expression_sub <- expression[expression$depmap_id %fin% depmap_ids,]
  
  # Add cancer type
  new_cols <- lapply(1:length(cell_line_cancer_type_mapping), function(i) {
    ct <- names(cell_line_cancer_type_mapping)[i]
    print(ct)
    depmaps_ct <- unlist(cell_line_cancer_type_mapping[[i]])
    intersect_depmaps <- intersect(expression_sub$depmap_id, 
                                   intersect(crispr_sub$depmap_id, depmaps_ct))
    return(data.table("depmap_id" = intersect_depmaps, 
                      "cancer_type" = rep(ct, times = length(intersect_depmaps))))
  })
  ct_df <- distinct(do.call(rbind, new_cols))
  print(head(ct_df))
  
  crispr_sub <- left_join(crispr_sub, ct_df, by = "depmap_id", 
                          relationship = "many-to-many")
  expression_sub <- left_join(expression_sub, ct_df, by = "depmap_id", 
                              relationship = "many-to-many")
  #print(head(expression_sub))
  #print(head(crispr_sub))
  #print(nrow(expression_sub[expression_sub[,driver_mut] == 1,]))
    
  # Ensure there are at least a certain number (e.g. 15) lines that do not have
  # a mutation in the given driver AND that do have a mutation in the given driver
  if(nrow(expression_sub[which(expression_sub[,driver_mut] == 0),]) < thres_cls) {
    print(paste("There are fewer than", 
                paste(thres_cls, 
                      paste("cell lines that lack a driver mutation in", driver_mut))))
    return(NA)
  }
  if(nrow(expression_sub[which(expression_sub[,driver_mut] == 1),]) < thres_cls) {
    print(paste("There are fewer than", 
                paste(thres_cls, 
                      paste("cell lines that have a driver mutation in", driver_mut))))
    return(NA)
  }
  
  # If desired, limit to DrugBank genes
  if(use_drugbank) {
    tophits <- tophits[tophits %fin% drugbank_gns]
  }
  
  # If desired, limit to pan-cancer hits
  if(use_pc_hits) {
    tophits <- tophits[tophits %fin% pc_hits]
  }
  print(head(tophits))
  
  # If desired, limit be differential dependency; particularly, significantly  
  # lower viability in the unmutated case
  if(use_dependency_check) {
    wilcoxon_res <- data.table("gene" = tophits)
    wilcoxon_res$pval <- unlist(lapply(wilcoxon_res$gene, function(g) {
      dep_mut <- crispr_sub[(crispr_sub$gene == g) & 
                              (crispr_sub[,driver_mut] == 1), 'value']
      dep_noMut <- crispr_sub[(crispr_sub$gene == g) &
                                (crispr_sub[,driver_mut] == 0), 'value']
      if((length(dep_noMut) < 1) | (length(dep_mut) < 1)) {return(NA)}
      if(is.na(dep_noMut) | is.na(dep_mut)) {return(NA)}
      res <- wilcox.test(dep_mut, dep_noMut, alternative = "greater", exact = F)
      return(res$p.value)
    }))
    wilcoxon_res$qval <- qvalue(wilcoxon_res$pval)$qvalues
    wilcoxon_res <- wilcoxon_res[order(wilcoxon_res$pval),]
    signif_hits <- wilcoxon_res[wilcoxon_res$qval < qval_thres, 'gene']  
    tophits <- tophits[tophits %fin% signif_hits]
  }
  
  tophits_depmap <- NA
  
  # Perform regression
  if(length(tophits) > 0) {
    expression_sub <- expression_sub[expression_sub$gene %fin% c(driver_mut),]
    expression_driver <- expression_sub[expression_sub$gene == driver_mut, ]
    
    driver_col <- as.integer(which(colnames(expression_driver) == driver_mut))
    expression_driver <- cbind(expression_driver[, c(driver_col)],
                               expression_driver[, c("depmap_id", "gene", 
                                                     "value", "cancer_type")])
    colnames(expression_driver)[1] <- "driver_mut_status"
    colnames(expression_driver)[
      which(colnames(expression_driver) == "value")] <- "expression_val"
    
    lm_dfs <- lapply(tophits, function(gene) {
      print(gene)
      if(!(gene %fin% crispr_sub$gene)) {return(NA)}
      crispr_sub <- crispr_sub[crispr_sub$gene %fin% c(gene),]
      crispr_gene <- crispr_sub[crispr_sub$gene == gene, c("depmap_id", "gene", 
                                                           "value", "cancer_type")]
      colnames(crispr_gene)[
        which(colnames(crispr_gene) == "value")] <- "crispr_val"
      df_list <- list(crispr_gene, expression_driver)
      #print(df_list)
      
      data_input <- merge(crispr_gene, expression_driver, 
                          by = c("depmap_id", "cancer_type"), all = T)
      data_input$driver_mut_status <- as.factor(data_input$driver_mut_status)
      #print(head(data_input))
      
      if(nrow(data_input) < 5) {return(NA)}
      
      lm_res <- NA
      tryCatch({
      	formula <- create_formula(terms_of_interest, data_input)
      	lm_res <- speedglm::speedlm(formula = formula, data = data_input)
        summary_table <- tidy(lm_res)
        print(summary_table)

        summary_table <- as.data.table(
          summary_table[, colSums(is.na(summary_table)) < nrow(summary_table)])
        	
      	if(recombine == "linearHypothesis") {
      		model <- lm(formula, data = data_input)
      		lm_res <- linearHypothesis(model, 
      			"driver_mut_status1 + expression_val = 0")
      		summary_table <- tidy(lm_res)
      		print(summary_table)
      		#pval <- as.numeric(attr(res, "value"))   #extract just the value of the linear function
      		#print(pval)
      	} 
      	else if (recombine == "multcomp") {
      		model <- lm(formula, data = data_input) 
      		lm_res <- glht(model, linfct = diag(length(coef(model))))  #linfct is a matrix defining the linear functions (hypotheses) to be tested
      		summary_table1 <- tidy(lm_res)
      		print(summary_table1)
      	} 
	else if (recombine == "browns") {
                data_input_t <- apply(t(data_input[,c("expression_val", "driver_mut_status")]), 2, as.numeric)
		colnames(data_input_t) <- unlist(data_input$depmap_id)
                rownames(data_input_t) <- c("expression_val", "driver_mut_status")
		data_input_t <- data_input_t[rowSums(is.na(data_input_t)) != ncol(data_input_t), ]
		data_input_t <- data_input_t[, colSums(is.na(data_input_t)) != nrow(data_input_t)]
		#data_input_t <- na.omit(data_input_t)
		#print(data_input_t)
		pvalues <- as.numeric(unlist(summary_table[summary_table$term %fin% c("expression_val", "driver_mut_status1"), 
					'p.value']))
		#print(pvalues)
		browns_res <- empiricalBrownsMethod(data_matrix = data_input_t, p_values = pvalues)
		print(browns_res)
		#pvals_new <- as.numeric(browns_res[["P_test"]])
		summary_table[2:3, 'p.value'] <- rep(browns_res, times = 2)
	}
	else if (recombine == "kosts") {
                data_input_t <- apply(t(data_input[,c("expression_val", "driver_mut_status")]), 2, as.numeric)
		colnames(data_input_t) <- unlist(data_input$depmap_id)
                rownames(data_input_t) <- c("expression_val", "driver_mut_status")
		data_input_t <- data_input_t[rowSums(is.na(data_input_t)) != ncol(data_input_t), ]
		data_input_t <- data_input_t[, colSums(is.na(data_input_t)) != nrow(data_input_t)]
		#data_input_t <- na.omit(data_input_t)
		pvalues <- as.numeric(unlist(summary_table[summary_table$term %fin% c("expression_val", "driver_mut_status1"), 
					'p.value']))
		#print(pvalues)
		kosts_res <- kostsMethod(data_matrix = data_input_t, p_values = pvalues)
		print(kosts_res)
		#pvals_new <- as.numeric(kosts_res[["P_test"]])
		summary_table[2:3, 'p.value'] <- rep(kosts_res, times = 2)

	}
	else if (recombine == "None") {
		summary_table <- summary_table 
	}
      	else {
      		print(paste(recombine, "is not currently implemented."))
      		return(NA)
      	}

      }, error = function(cond) {
        print(cond)
      })
      
      #print(summary_table)
      
      if(length(summary_table) > 1) {
        summary_table <- summary_table[summary_table$term != "(Intercept)",]
        lm_res_terms <- summary_table$term
        lm_res_stats <- as.numeric(summary_table$estimate)
        lm_res_pvals <- summary_table$p.value
        
        if(recombine == "multcomp") {
        	summary_table1 <- summary_table1[2:nrow(summary_table1),]
        	lm_res_pvals <- as.numeric(summary_table1$adj.p.value)
        	lm_res_stats <- as.numeric(summary_table1$estimate)
        }

        print(lm_res_pvals)
        return(data.frame("stat" = lm_res_stats, "pval" = lm_res_pvals, 
                          "term" = lm_res_terms,
                          "gene" = rep(gene, times = length(lm_res_stats))))
      } else {return(NA)}
    })
    
    lm_dfs <- lm_dfs[!is.na(lm_dfs)]
    lm_res <- do.call(rbind, lm_dfs)
    
    lm_res <- na.omit(lm_res)
    if(length(lm_res) == 0) {return(NA)}
    if(nrow(lm_res) == 0) {return(NA)}
    
    unique_terms <- unique(lm_res$term)
    adj_dfs <- lapply(unique_terms, function(term) {
      lm_res_sub <- lm_res[lm_res$term == term,]
      #lm_res_sub$padj <- p.adjust(lm_res_sub$pval, "BH")
      lm_res_sub$qvalue <- qvalue(lm_res_sub$pval)$qvalues
      return(lm_res_sub)
    })
    lm_res_sig <- do.call(rbind, adj_dfs)
    #lm_res_sig <- lm_res_sig[order(lm_res_sig$padj),]
    lm_res_sig <- lm_res_sig[order(lm_res_sig$qvalue),]
    
    tophits_depmap <- lm_res_sig[(lm_res_sig$pval < 1), ] #& 
    					# (lm_res_sig$stat > 0),]
    if(length(tophits_depmap) == 0) {return(NA)}
    if(nrow(tophits_depmap) == 0) {return(NA)}

    # Uncomment to create associated plot
    #dependency_by_expression_analysis_plot(crispr_sub, expression_sub, 
    #                                       tophits_depmap[!(grepl("cancer_type", 
    #                                                              tophits_depmap$term)), 
    #                                                        'gene'], 
    #                                       driver_mut, 0)
    
    return(tophits_depmap)
  } else {return(NA)}
}


#' Construct linear regression formula using a set of variables indicating which
#' terms we are including in a given model
#' @param terms_of_interest a vector of driver terms of interest to include in the
#' model, possible options are 'mutation' and 'expression'
#' @param data_input a data frame with input to the regression, each column
#' containing values for a particular term
create_formula <- function(terms_of_interest, data_input) {
  formula <- "crispr_val ~ "
  
  # Driver mutation status
  if(('mutation' %fin% terms_of_interest) & 
     (length(unique(data_input$driver_mut_status)) > 1)) {
    formula <- paste0(formula, 'driver_mut_status ')
  }
    
  # Driver expression
  if('expression' %fin% terms_of_interest) {
    formula <- paste(formula, 'expression_val ', sep = "+ ")
  }
    
  # Cancer type
  # if(length(unique(data_input$cancer_type)) > 1) {
  if('cancer_type' %fin% colnames(data_input)) {
    formula <- paste(formula, 'cancer_type', sep = "+ ")
  }
  
  print(formula)
  if(formula == "crispr_val ~ ") { return(NA) }
  else {return(formula)}
}


# Recombine p-values for targets across modalities
#' Using the metap package, recombine the p-values for target gene results for
#' both driver expression and driver mutation status
#' @param results a data frame of results from get_synthetic_lethals
#' without subsetting to significant hits- change p-value threshold to 1)
recombine_pvals_byterm_fisher <- function(results) {
  unique_terms <- unique(results$term)
  unique_targets <- unique(results$gene)
  
  # For each of these, get the p-values from across the various cancer types 
  new_pvals <- lapply(unique_targets, function(targ) {
    print(targ) 
    results_targ <- results[results$gene == targ, ]
    if(nrow(results_targ) == 1) {
      return(data.frame("pval" = results_targ$pval, "gene" = targ))}
    pvals <- results_targ$pval
    
    # Do the correction using these p-values
    pval_new <- metap::sumlog(pvals, log.p = F)$p
    return(data.frame("pval" = pval_new, "gene" = targ))
  })
  
  new_df <- do.call(rbind, new_pvals)
  
  # Multiple hypothesis testing correction
  new_df$qval <- qvalue(new_df$pval)$qvalues
  #new_df$p.adj <- p.adjust(new_df$pval, "BH")
  
  return(new_df)
}


############################################################
### CALL FUNCTIONS
############################################################
# Run this across multiple drivers in one call
drivers_of_interest <- c("TP53", "PIK3CA", "KRAS")

lapply(drivers_of_interest, function(driver_name) {
  perCancer_driver <- lapply(perCancer, function(x) {
    if(driver_name %fin% x$R_i.name) return(x)
    else {return(NA)}
  })
  perCancer_driver <- perCancer_driver[!is.na(perCancer_driver)]
  
  # Get pan-cancer hits at a given threshold
  pc_qval_thres <- 0.2
  pc_hits_driver <- pc_allGenes[(pc_allGenes$R_i.name == driver_name) &
                                  (pc_allGenes$q.value < pc_qval_thres), 
                                'T_k.name']
    
  ## No p-value recombination ##
  lm_results_panq0.2_perq0.2 <- perform_lm_analysis_crosscancers(
    crispr, expression_qn, driver_name, perCancer_driver, 0.2, T, 
    pc_hits_driver, F, NA, F, cell_line_cancer_type_mapping,
    c('mutation', 'expression'), "None")
  lm_results_panq0.2_perq0.2 <- lm_results_panq0.2_perq0.2[!grepl("cancer_type", 
                                                                  lm_results_panq0.2_perq0.2$term),]
  write.csv(lm_results_panq0.2_perq0.2, 
            paste0(PATH_OUT, paste0(driver_name, 
                                "/synthetic_lethals_pancancer_regression_qn_panq0.2_perq0.2.csv")))

  ## Fisher method p-value recombination ##
  lm_results_panq0.2_perq0.2 <- recombine_pvals_byterm_fisher(lm_results_panq0.2_perq0.2)
  lm_results_panq0.2_perq0.2 <- lm_results_panq0.2_perq0.2[order(lm_results_panq0.2_perq0.2$qval),]
  write.csv(lm_results_panq0.2_perq0.2, 
            paste0(PATH_OUT, paste0(driver_name, 
                                "/synthetic_lethals_pancancer_regression_qn_panq0.2_perq0.2_fishersp.csv")))
  
  ## Linear hypothesis testing p-value recombination ##
  lm_results_panq0.2_perq0.2_linHyp <- perform_lm_analysis_crosscancers(
    crispr, expression_qn, driver_name, perCancer_driver, 0.2, T, 
    pc_hits_driver, F, NA, F, cell_line_cancer_type_mapping,
    c('mutation', 'expression'), "linearHypothesis")
  write.csv(lm_results_panq0.2_perq0.2_linHyp, 
           paste0(PATH_OUT, paste0(driver_name, 
                               "/synthetic_lethals_pancancer_regression_qn_panq0.2_perq0.2_linHyp.csv")))
  
  ## Simultaneous inference p-value recombination ##                        
  lm_results_panq0.2_perq0.2_multcomp <- perform_lm_analysis_crosscancers(
    crispr, expression_qn, driver_name, perCancer_driver, 0.2, T, 
    pc_hits_driver, F, NA, F, cell_line_cancer_type_mapping,
    c('mutation', 'expression'), "multcomp")
  lm_results_panq0.2_perq0.2_multcomp <- lm_results_panq0.2_perq0.2_multcomp[!grepl("cancer_type", 
                                                                  lm_results_panq0.2_perq0.2_multcomp$term),]
  write.csv(lm_results_panq0.2_perq0.2_multcomp, 
           paste0(PATH_OUT, paste0(driver_name, 
                               "/synthetic_lethals_pancancer_regression_qn_panq0.2_perq0.2_multcomp.csv")))

  ## Brown's method p-value recombination ##  
  lm_results_panq0.2_perq0.2_browns <- perform_lm_analysis_crosscancers(
    crispr, expression_qn, driver_name, perCancer_driver, 0.2, T, 
    pc_hits_driver, F, NA, F, cell_line_cancer_type_mapping,
    c('mutation', 'expression'), "browns")
  lm_results_panq0.2_perq0.2_browns <- lm_results_panq0.2_perq0.2_browns[!grepl("cancer_type", 
                                                                  lm_results_panq0.2_perq0.2_browns$term),]
  write.csv(lm_results_panq0.2_perq0.2_browns, 
           paste0(PATH_OUT, paste0(driver_name, 
                               "/synthetic_lethals_pancancer_regression_qn_panq0.2_perq0.2_browns.csv")))

  ## Kost's method p-value recombination ##  
  lm_results_panq0.2_perq0.2_kosts <- perform_lm_analysis_crosscancers(
    crispr, expression_qn, driver_name, perCancer_driver, 0.2, T, 
    pc_hits_driver, F, NA, F, cell_line_cancer_type_mapping,
    c('mutation', 'expression'), "kosts")
  lm_results_panq0.2_perq0.2_kosts <- lm_results_panq0.2_perq0.2_kosts[!grepl("cancer_type", 
                                                                  lm_results_panq0.2_perq0.2_kosts$term),]
  write.csv(lm_results_panq0.2_perq0.2_kosts, 
           paste0(PATH_OUT, paste0(driver_name, 
                               "/synthetic_lethals_pancancer_regression_qn_panq0.2_perq0.2_kosts.csv")))

})
