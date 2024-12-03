############################################################
# Code to Create Figure 5 Visualizations
# Written by Sara Geraghty
# https://www.biorxiv.org/content/10.1101/2024.11.20.624509v1
############################################################

library(data.table)
library(ggplot2)
library(survival)
library(survminer)

# Local PATH to directory containing Dyscovr output files
PATH <-  paste0(getwd(), "Output/")

# Generalized ID conversion table from BiomaRt
all_genes_id_conv <- read.csv(paste0(PATH, "all_genes_id_conv.csv"), header = T)

# Source general, widely-used functions for speed enhancements
source(general_important_functions.R)


########################################################################
### PART B: DEPMAP DRUG SENSITIVITY SPEARMAN WITH KBTBD2 KO 
########################################################################
# Correlate drug sensitivity to target gene dependency 
# AUC(CTD^2) -- better quality and more cell lines
drug_sensitivity <- read.csv(paste0(PATH, "DepMap/Drug_sensitivity_AUC_(CTD^2).csv"),
                             header = T, check.names = F)
colnames(drug_sensitivity)[1] <- "depmap_id"

# Read in CRISPR data and reorganize
crispr <- read.csv(paste0(PATH, "DepMap/CRISPR_(DepMap_Public_23Q2+Score,_Chronos).csv"),
                             header = T, check.names = F)
colnames(crispr)[1] <- "depmap_id"
crispr_melt <- melt(crispr)
#print(head(depmap_df_melt))

# Adjust column names
depmap_id_col <- ifelse(grepl("ACH", crispr_melt[1,1]), 1, 2)
colnames(crispr_melt)[depmap_id_col] <- "depmap_id"
gene_id_col <- ifelse(depmap_id_col == 1, 2, 1)
colnames(crispr_melt)[gene_id_col] <- "gene"
#print(head(depmap_df_melt))

# Make the "value" column numeric, the "gene" column a character
crispr_melt$value <- as.numeric(unlist(crispr_melt$value))
crispr_melt$gene <- as.character(unlist(crispr_melt$gene))

# Annotate file with GOI dependency for each cell line
goi <- "KBTBD2"
crispr_sub <- crispr_melt[crispr_melt$gene == goi, c("depmap_id", "value")]
colnames(crispr_sub) <- c("depmap_id", paste0(goi, "_dependency"))
drug_sensitivity <- merge(crispr_sub, drug_sensitivity, by = "depmap_id")
drug_sensitivity <- drug_sensitivity[,c(1,2,10:ncol(drug_sensitivity))]

#' Calculate correlation between GOI dependency and sensitivity to each drug
#' @param drug_sensitivity a merged DF with drug sensitivity and GOI dependency data
#' @param corr_type either 'spearman' or 'pearson' to indicate the type of correlation
#' @param qval_thres a q-value threshold for significance
#' @param target_label the name of the target gene of interest
calc_corr_dep_sens <- function(drug_sensitivity, corr_type, qval_thres, target_label) {
  dep_vals <- drug_sensitivity[,2]
  print(head(dep_vals))
  
  corr_dfs <- lapply(3:ncol(drug_sensitivity), function(i) {
    drug <- colnames(drug_sensitivity)[i]
    print(drug)
    
    sens_vals <- drug_sensitivity[,i]
    corr_res <- NA
    
    if(tolower(corr_type) == "spearman") {
      corr_res <- cor.test(dep_vals, sens_vals, method = "spearman")
    } else if(tolower(corr_type) == "pearson") {
      corr_res <- cor.test(dep_vals, sens_vals, method = "pearson")
    } else {
      print("Only implemented for Spearman and Pearson. Please choose one.")
    }
    outdf <- data.frame("drug" = drug, "corr" = as.numeric(unlist(corr_res$estimate)), 
                        "pval" = as.numeric(unlist(corr_res$p.value)))
    return(outdf)
  })
  corr_df <- do.call(rbind, corr_dfs)
  
  # MHT correction
  corr_df$qval <- qvalue(corr_df$pval)$qvalues
  corr_df$neglog10qval <- unlist(lapply(corr_df$qval, function(q) -log10(q)))
  #corr_df$neglog10pval <- unlist(lapply(corr_df$pval, function(p) -log10(p)))
  
  corr_df$up_or_down <- unlist(lapply(1:nrow(corr_df), function(i) {
    if(corr_df[i, 'qval'] < qval_thres) {
      beta_sign <- ifelse(corr_df[i, 'corr'] > 0, 1, 0) 
      if(beta_sign == 1) {return("pos")}
      else {return("neg")}
    } else {return("ns")}
  }))
  
  corr_df$drug_label <- unlist(lapply(1:nrow(corr_df), function(i) 
    ifelse(corr_df$qval[i] < (qval_thres*0.3), 
           unlist(strsplit(corr_df$drug[i], " (", fixed = T))[1], NA)))
  print(head(corr_df))
  
  # Create visualization
  cols <- c("pos" = "#BC3C29FF", "neg" = "#0072B5FF", "ns" = "grey") 
  alphas <- c("pos" = 0.95, "neg" = 0.95, "ns" = 0.5)
  
  p <- ggplot(corr_df, aes(x = corr, y = neglog10qval, fill = up_or_down, 
                           alpha = up_or_down, label = drug_label)) + 
    geom_point(shape = 21, size = 3) + 
    xlab(paste(corr_type, "Correlation,", 
               paste(target_label, "CRISPRi Dependency +\n Drug Sensitivity (AUC)"))) + 
    ylab("-log10(q-value)") + theme_minimal() + 
    geom_hline(yintercept = -log10(qval_thres), linetype = "dashed") +
    scale_fill_manual(values = cols) + scale_alpha_manual(values = alphas) + 
    theme(legend.position = "none", 
          axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 12)) + 
    geom_text_repel(min.segment.length = unit(0.1, "lines"), size = 5, 
                    force = 3.5, box.padding = 0.4) 
  print(p)
  
  return(corr_df[order(corr_df$qval),])
}

# Call function
corr_df <- calc_corr_dep_sens(drug_sensitivity, "Spearman", 0.2, goi)


########################################################################
### PART C-D: SURVIVAL CURVES 
########################################################################
#' Create a survival curve for a given driver gene (and, optionally, expression
#' of a predicted target gene) using TCGA clinical data
#' @param clinical_df a clinical DF from the TCGA for a certain patient set
#' @param driver_name the gene name of the driver of interest
#' @param mutation_count_matrix a mutation count matrix from maftools for TCGA data
#' @param target_name the ENSG ID of a target gene of interest
#' @param expression_df an expression DF (quantile-normalized)
#' @param target_label the gene name of the target gene, for labeling
#' @param patient_sample_df a patient/sample combined DF with regression variables
#' for TCGA clinical data
#' @param surv_curve_type either log-rank (the default) or a Cox proportional 
#' hazards model (adjusted survival)
#' @param method a method for selecting "high" and "low" patients for target gene
#' expression, either 'znorm', 'thirds', or quantile'
create_survival_curves <- function(clinical_df, driver_name, mutation_count_matrix,
                                   target_gene, expression_df, target_label,
                                   patient_sample_df, surv_curve_type = "log-rank",
                                   method = "znorm") {
  
  clinical_df_survival <- clinical_df[,colnames(clinical_df) %fin% 
                                        c("case_submitter_id", "vital_status",
                                          "days_to_death", "days_to_last_follow_up",
                                          "project_id")]
  clinical_df_survival <- distinct(clinical_df_survival[which(!is.na(
    clinical_df_survival$vital_status)),])
  #table(clinical_df_survival$vital_status)   #11532 alive, 4182 dead, 14 not reported 
  
  clinical_df_survival$status <- ifelse(clinical_df_survival$vital_status == 'Alive', 0, 1)
  clinical_df_survival$patient <- unlist(lapply(clinical_df_survival$case_submitter_id, function(x)
    unlist(strsplit(x, "-", fixed = T))[3]))
  clinical_df_survival$cancer_type <- unlist(lapply(clinical_df_survival$project_id, function(x)
    unlist(strsplit(x, "-", fixed = T))[2]))
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
    
    # Uncomment to exclusively look at mutated patients
    #mutation_count_matrix_driver <- mutation_count_matrix_driver[
    #  mutation_count_matrix_driver$mutation_status > 0,]
    
    mutation_count_matrix_driver$patient <- unlist(lapply(mutation_count_matrix_driver$patient, function(x)
      unlist(strsplit(as.character(x), "-", fixed = T))[1]))
    
    # Binarize mutation status
    mutation_count_matrix_driver$mutation_status <- unlist(lapply(
      mutation_count_matrix_driver$mutation_status, function(x) ifelse(x > 0, 1, 0)))
    
  } else {
    mutation_count_matrix_driver <- mutation_count_matrix
  }
  
  survival_df <- merge(mutation_count_matrix_driver, 
                       clinical_df_survival[,c("patient","days","status","cancer_type")],
                       by = "patient", all = F)
  
  # Optionally add in target gene expression information as well
  expression_df_target <- NA
  if(!is.na(target_gene)) {
    expression_df_target <- get_target_expr_categorization(expression_df, 
                                                           target_gene, method)
    survival_df <- merge(survival_df, expression_df_target[,c('patient', 'exp_status')],
                         by = "patient", all = F)
  }
  
  # Add patient and sample characteristics, if using a Cox model
  if(surv_curve_type == "cox") {
    confounders <- colnames(patient_sample_df)[grepl("Tumor_purity", 
                                                     colnames(patient_sample_df))]
    confounders <- c(confounders, colnames(patient_sample_df)[grepl("Tot_IC_Frac", 
                                                                    colnames(patient_sample_df))])
    confounders <- c(confounders, colnames(patient_sample_df)[grepl("PC", 
                                                                    colnames(patient_sample_df))])
    confounders <- c(confounders, colnames(patient_sample_df)[grepl("Tot_Mut", 
                                                                    colnames(patient_sample_df))])
    confounders <- c(confounders, colnames(patient_sample_df)[grepl("Treatment", 
                                                                    colnames(patient_sample_df))])
    confounders <- c(confounders, colnames(patient_sample_df)[grepl("Cancer_type", 
                                                                    colnames(patient_sample_df))])
    confounders <- c(confounders, "Age", "Gender", "Prior_malig", "patient")
    patient_sample_df_sub <- patient_sample_df[,colnames(patient_sample_df) %fin% confounders]
    
    # Remove columns that are entirely zero
    if(any(grepl("Cancer_type", colnames(patient_sample_df_sub)))) {
      patient_sample_df_ct <- patient_sample_df_sub[,grepl("Cancer_type", colnames(patient_sample_df_sub))]
      print(colSums(patient_sample_df_ct))
      patient_sample_df_ct <- patient_sample_df_ct[, which(colSums(patient_sample_df_ct) != 0)]
      patient_sample_df_sub <- cbind(patient_sample_df_sub[,!grepl("Cancer_type", colnames(patient_sample_df_sub))],
                                     patient_sample_df_ct)
    }
    
    patient_sample_df_sub$patient <- unlist(lapply(patient_sample_df_sub$patient, function(x)
      unlist(strsplit(x, "-", fixed = T))[1]))
    
    survival_df <- merge(survival_df, patient_sample_df_sub, by = "patient", all = F)
  }
  
  survival_df <- survival_df[survival_df$days != "'--",]
  survival_df$days <- as.numeric(survival_df$days)
  print(head(survival_df))
  
  # Use the survminer package to create the curves
  fit <- NA
  if(!is.na(target_gene)) {
    if(surv_curve_type == "log-rank") {
      fit <- survfit(Surv(time = days, event = status) ~ mutation_status + exp_status, #+ treated, 
                     data = survival_df)
      #fit <- survfit(Surv(time = days, event = status) ~ exp_status, #+ treated, 
      #               data = survival_df)
      print(head(fit))
      
      p <- ggsurvplot(fit, data = survival_df, 
                      palette = c("#0072B5FF", "#FFDC91FF", "#BC3C29FF", "#20854EFF"),
                      #palette = c("#0072B5FF", "#BC3C29FF"),
                      pval = T, ggtheme = theme_minimal(), 
                      legend.title = paste(driver_name, 
                                           paste("Mutation Status +\n", 
                                                 paste(target_label, "Expression Status"))), 
                      #legend.title = paste(target_label, "Expression Status"), 
                      legend = c(0.8,0.8), # "bottom",
                      legend.labs = c("No Mutation, Low Expression", 
                                      "No Mutation, High Expression", 
                                      "Mutation, Low Expression", 
                                      "Mutation, High Expression"), 
                      font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                      font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                      font.legend = c(10, "bold", "black")) + 
        xlab("Time (Days)") + ylab("Survival Probability")
      print(p)
      
      # Cox proportional hazards
    } else {
      formula <- "Surv(time = days, event = status) ~ mutation_status + exp_status + Age + Gender + Prior_malig + Treatment_rad + Treatment_pharm + PC1 + PC2 + PC3 + Tumor_purity_b1 + Tumor_purity_b2 + Tumor_purity_b3 + Tot_IC_Frac_b1 + Tot_IC_Frac_b2 + Tot_IC_Frac_b3"
      if(any(grepl("Cancer_type", colnames(survival_df)))) {
        formula <- paste(formula, paste(
          colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
          collapse = " + "), sep = " + ")
      }
      res.cox <- coxph(as.formula(formula), data = survival_df)  #+ cancer_type
      print(res.cox)
      p <- ggadjustedcurves(res.cox, data = survival_df, variable = "exp_status",
                            legend.title = paste(target_label, "Expression Status"), 
                            palette = c("#0072B5FF", "#BC3C29FF"),
                            pval = T, ggtheme = theme_minimal(), 
                            legend = "bottom", legend.labs = c("Low Expression", "High Expression"), 
                            font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                            font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                            font.legend = c(12, "bold", "black"))
      print(p)
    }
  } else {
    if(surv_curve_type == "log-rank") {
      fit <- survfit(Surv(time = days, event = status) ~ mutation_status, #+ treated, 
                     data = survival_df)
      print(head(fit))
      
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
      
      # Cox proportional hazards
    } else {
      formula <- "Surv(time = days, event = status) ~ mutation_status + Age + Gender + Prior_malig + Treatment_rad + Treatment_pharm + PC1 + PC2 + PC3 + Tumor_purity_b1 + Tumor_purity_b2 + Tumor_purity_b3 + Tot_IC_Frac_b1 + Tot_IC_Frac_b2 + Tot_IC_Frac_b3"
      if(any(grepl("Cancer_type", colnames(survival_df)))) {
        formula <- paste(formula, paste(
          colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
          collapse = " + "), sep = " + ")
      }
      res.cox <- coxph(as.formula(formula), data = survival_df)  #+ cancer_type
      print(res.cox)
      p <- ggadjustedcurves(res.cox, data = survival_df, variable = "mutation_status",
                            legend.title = paste(driver_name, "Mutation Status"), 
                            palette = c("#0072B5FF", "#BC3C29FF"),
                            pval = T, ggtheme = theme_minimal(), 
                            legend = "bottom", legend.labs = c("No Mutation", "Mutation"), 
                            font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                            font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                            font.legend = c(12, "bold", "black"))
      print(p)
    }
    
  }
  
  return(survival_df)
}

#' Helper function to add expression categorization by patient for the given 
#' target gene; type of categorization is specified
#' @param expression_df a quantile-normalized expression data frame across genes
#' and patient samples
#' @param target_gene the name of the target gene of interest
#' @param method the method of categorization, either 'thirds' for taking the top
#' and bottom thirds, 'quantile' to take those below Q1 or above Q3, or 'znorm' 
#' to take >1.5 sigma and <-1.5 sigma for a z-score normalized expression matrix
get_target_expr_categorization <- function(expression_df, target_gene, method) {
  
  expression_df_target <- as.data.frame(t(expression_df[expression_df$ensg_id == 
                                                          target_gene,]))
  expression_df_target$patient <- rownames(expression_df_target)
  expression_df_target <- expression_df_target[2:nrow(expression_df_target),]
  colnames(expression_df_target)[1] <- "expression"
  
  # Get patients that are have high expression (e.g. top 1/3) or low expression
  # (e.g. bottom 1/3) for the given target gene
  if(method == "thirds") {
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
    
    # Get patients that are above the third quantile or below the first quantile
  } else if(method == "quantile") {
    expr_vals <- as.numeric(unlist(expression_df_target$expression))
    
    # Get Q1 and Q3 of data
    q1 <- as.numeric(quantile(expr_vals)[2])
    q3 <- as.numeric(quantile(expr_vals)[4])
    lt_q1 <- expression_df_target[expression_df_target$expression < q1, ]
    gr_q3 <- expression_df_target[expression_df_target$expression > q3, ]
    
    # Add an "expression status" term for patients that fall in these buckets
    expression_df_target$exp_status <- unlist(lapply(expression_df_target$patient, function(id) {
      if(id %fin% gr_q3$patient) {return(1)}
      else if (id %fin% lt_q1$patient) {return(0)}
      else {return(NA)}
    }))
    
  } else if (method == "znorm") {
    expr_vals <- as.numeric(unlist(expression_df_target$expression))
    
    sigma1.5 <- sd(expr_vals) * 1.5
    sigma <- sd(expr_vals)
    mean <- mean(expr_vals)
    
    low <- expression_df_target[expression_df_target$expression < (mean - sigma), ]
    high <- expression_df_target[expression_df_target$expression > (mean + sigma), ]
    
    expression_df_target$exp_status <- unlist(lapply(expression_df_target$patient, function(id) {
      if(id %fin% high$patient) {return(1)}
      else if (id %fin% low$patient) {return(0)}
      else {return(NA)}
    }))
    
  } else {
    print("Currently only implemented for methods 'thirds', 'quantile', and 'znorm'. 
          Please try again with one of these methods specified.")
  }
  
  expression_df_target <- na.omit(expression_df_target)
  expression_df_target$patient <- unlist(lapply(expression_df_target$patient, function(x)
    unlist(strsplit(as.character(x), "-", fixed = T))[1]))
  print(head(expression_df_target))
  
  return(expression_df_target)
}


# Import pan-cancer clinical DF
clinical_df <- read.csv(paste0(PATH, "clinical.csv"), header = T, check.names = F)

# Import the pan-cancer mutation count matrix
mutation_count_matrix <- read.csv(paste0(PATH, "Mutation/mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv"), 
                                  header = T, check.names = F)
# Import the pan-cancer expression matrix
expression_df_qn <- read.csv(paste0(PATH, "Expression/expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv"), 
                             header = T, check.names = F)
expression_df_znorm <- cbind(expression_df_qn[,1], 
                             as.data.frame(scale(expression_df_qn[,2:ncol(expression_df_qn)])))
colnames(expression_df_znorm)[1] <- "ensg_id"

# Import patient sample DF
patient_sample_df <- read.csv(paste0(PATH, "Patient/combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv"),
                              header = T, check.names = F)
colnames(patient_sample_df)[1] <- "patient"


# Call function (e.g. PIK3CA and KBTBD2)
pik3ca_survival_res <- create_survival_curves(clinical_df, "PIK3CA", mutation_count_matrix,
                                              "ENSG00000170852", expression_df_znorm, 
                                              "KBTBD2", patient_sample_df, 
                                              method = "znorm", surv_curve_type = "cox")


# Call function for breast cancer, specifically
clinical_df_brca <- clinical_df[clinical_df$project_id == "TCGA-BRCA",]
brca_patients <- unlist(lapply(clinical_df_brca$case_submitter_id, function(x)
  unlist(strsplit(x, "-", fixed = T))[3]))
mutation_pats <- unlist(lapply(colnames(mutation_count_matrix)[2:ncol(mutation_count_matrix)], function(x)
  unlist(strsplit(x, "-", fixed = T))[1]))
expression_pats <- unlist(lapply(colnames(expression_df_znorm)[2:ncol(expression_df_znorm)], function(x)
  unlist(strsplit(x, "-", fixed = T))[1]))
patient_sample_df_sub <- patient_sample_df[,!grepl("Cancer_type", colnames(patient_sample_df))]
pik3ca_survival_res <- create_survival_curves(clinical_df_brca, "PIK3CA", 
                                              mutation_count_matrix[,c(1, (which(mutation_pats %fin% brca_patients)+1))],
                                              "ENSG00000170852", 
                                              expression_df_znorm[,c(1, (which(expression_pats %fin% brca_patients)+1))], 
                                              "KBTBD2", patient_sample_df_sub,
                                              method = "znorm", surv_curve_type = "cox")
