########################################################################
### SURVIVAL CURVES AND DRUG RESPONSE ANALYSIS
### Written by Sara Geraghty, Princeton University
### https://www.biorxiv.org/content/10.1101/2024.11.20.624509v1
########################################################################

# Plot survival curves and perform drug sensitivity analysis using TCGA clinical data

library(ggplot2)
library(TCGAbiolinks)
library(data.table)
library(survival)
library(survminer)
library(dplyr)

PATH <- getwd()

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
    patient_sample_df_ct <- patient_sample_df_sub[,grepl("Cancer_type", colnames(patient_sample_df_sub))]
    print(colSums(patient_sample_df_ct))
    patient_sample_df_ct <- patient_sample_df_ct[, which(colSums(patient_sample_df_ct) != 0)]
    patient_sample_df_sub <- cbind(patient_sample_df_sub[,!grepl("Cancer_type", colnames(patient_sample_df_sub))],
                                   patient_sample_df_ct)
    
    patient_sample_df_sub$patient <- unlist(lapply(patient_sample_df_sub$patient, function(x)
      unlist(strsplit(x, "-", fixed = T))[1]))
    
    survival_df <- merge(survival_df, patient_sample_df_sub, by = "patient", all = F)
  }

  survival_df <- survival_df[survival_df$days != "'--",]
  survival_df$days <- as.numeric(survival_df$days)
  #survival_df$treated <- unlist(lapply(survival_df$response, function(x) 
  #  ifelse(is.na(x), 0, 1))) 
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
      formula <- paste(formula, paste(
        colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
        collapse = " + "), sep = " + ")
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
      formula <- paste(formula, paste(
        colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
        collapse = " + "), sep = " + ")
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
#' @param meth   
                      legend = "bottom", legend.labs = c("No Mutation", "Mutation"), 
                      font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                      font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                      font.legend = c(12, "bold", "black")) + 
        xlab("Time (Days)") + ylab("Survival Probability")
      print(p)
      
      # Cox proportional hazards
    } else {
      formula <- "Surv(time = days, event = status) ~ mutation_status + Age + Gender + Prior_malig + Treatment_rad + Treatment_pharm + PC1 + PC2 + PC3 + Tumor_purity_b1 + Tumor_purity_b2 + Tumor_purity_b3 + Tot_IC_Frac_b1 + Tot_IC_Frac_b2 + Tot_IC_Frac_b3"
      formula <- paste(formula, paste(
        colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
        collapse = " + "), sep = " + ")
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
#' @param expression_df a quantile-normalized expression data fram  across genes
#' and patient samples
#' @param target_gene the name of the target gene of interest
#' @param meth   
                      legend = "bottom", legend.labs = c("No Mutation", "Mutation"), 
                      font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                      font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                      font.legend = c(12, "bold", "black")) + 
        xlab("Time (Days)") + ylab("Survival Probability")
      print(p)
      
      # Cox proportional hazards
    } else {
      formula <- "Surv(time = days, event = status) ~ mutation_status + Age + Gender + Prior_malig + Treatment_rad + Treatment_pharm + PC1 + PC2 + PC3 + Tumor_purity_b1 + Tumor_purity_b2 + Tumor_purity_b3 + Tot_IC_Frac_b1 + Tot_IC_Frac_b2 + Tot_IC_Frac_b3"
      formula <- paste(formula, paste(
        colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
        collapse = " + "), sep = " + ")
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
#' @param expression_df a quantile-normalized expression data fram  across genes
#' and patient samples
#' @param target_gene the name of the target gene of interest
#' @param meth   
                      legend = "bottom", legend.labs = c("No Mutation", "Mutation"), 
                      font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                      font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                      font.legend = c(12, "bold", "black")) + 
        xlab("Time (Days)") + ylab("Survival Probability")
      print(p)
      
      # Cox proportional hazards
    } else {
      formula <- "Surv(time = days, event = status) ~ mutation_status + Age + Gender + Prior_malig + Treatment_rad + Treatment_pharm + PC1 + PC2 + PC3 + Tumor_purity_b1 + Tumor_purity_b2 + Tumor_purity_b3 + Tot_IC_Frac_b1 + Tot_IC_Frac_b2 + Tot_IC_Frac_b3"
      formula <- paste(formula, paste(
        colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
        collapse = " + "), sep = " + ")
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
#' @param expression_df a quantile-normalized expression data fram  across genes
#' and patient samples
#' @param target_gene the name of the target gene of interest
#' @param meth   
                      legend = "bottom", legend.labs = c("No Mutation", "Mutation"), 
                      font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                      font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                      font.legend = c(12, "bold", "black")) + 
        xlab("Time (Days)") + ylab("Survival Probability")
      print(p)
      
      # Cox proportional hazards
    } else {
      formula <- "Surv(time = days, event = status) ~ mutation_status + Age + Gender + Prior_malig + Treatment_rad + Treatment_pharm + PC1 + PC2 + PC3 + Tumor_purity_b1 + Tumor_purity_b2 + Tumor_purity_b3 + Tot_IC_Frac_b1 + Tot_IC_Frac_b2 + Tot_IC_Frac_b3"
      formula <- paste(formula, paste(
        colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
        collapse = " + "), sep = " + ")
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
#' @param expression_df a quantile-normalized expression data fram  across genes
#' and patient samples
#' @param target_gene the name of the target gene of interest
#' @param meth   
                      legend = "bottom", legend.labs = c("No Mutation", "Mutation"), 
                      font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                      font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                      font.legend = c(12, "bold", "black")) + 
        xlab("Time (Days)") + ylab("Survival Probability")
      print(p)
      
      # Cox proportional hazards
    } else {
      formula <- "Surv(time = days, event = status) ~ mutation_status + Age + Gender + Prior_malig + Treatment_rad + Treatment_pharm + PC1 + PC2 + PC3 + Tumor_purity_b1 + Tumor_purity_b2 + Tumor_purity_b3 + Tot_IC_Frac_b1 + Tot_IC_Frac_b2 + Tot_IC_Frac_b3"
      formula <- paste(formula, paste(
        colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
        collapse = " + "), sep = " + ")
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
#' @param expression_df a quantile-normalized expression data fram  across genes
#' and patient samples
#' @param target_gene the name of the target gene of interest
#' @param meth   
                      legend = "bottom", legend.labs = c("No Mutation", "Mutation"), 
                      font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                      font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                      font.legend = c(12, "bold", "black")) + 
        xlab("Time (Days)") + ylab("Survival Probability")
      print(p)
      
      # Cox proportional hazards
    } else {
      formula <- "Surv(time = days, event = status) ~ mutation_status + Age + Gender + Prior_malig + Treatment_rad + Treatment_pharm + PC1 + PC2 + PC3 + Tumor_purity_b1 + Tumor_purity_b2 + Tumor_purity_b3 + Tot_IC_Frac_b1 + Tot_IC_Frac_b2 + Tot_IC_Frac_b3"
      formula <- paste(formula, paste(
        colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
        collapse = " + "), sep = " + ")
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
#' @param expression_df a quantile-normalized expression data fram  across genes
#' and patient samples
#' @param target_gene the name of the target gene of interest
#' @param meth   
                      legend = "bottom", legend.labs = c("No Mutation", "Mutation"), 
                      font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                      font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                      font.legend = c(12, "bold", "black")) + 
        xlab("Time (Days)") + ylab("Survival Probability")
      print(p)
      
      # Cox proportional hazards
    } else {
      formula <- "Surv(time = days, event = status) ~ mutation_status + Age + Gender + Prior_malig + Treatment_rad + Treatment_pharm + PC1 + PC2 + PC3 + Tumor_purity_b1 + Tumor_purity_b2 + Tumor_purity_b3 + Tot_IC_Frac_b1 + Tot_IC_Frac_b2 + Tot_IC_Frac_b3"
      formula <- paste(formula, paste(
        colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
        collapse = " + "), sep = " + ")
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
#' @param expression_df a quantile-normalized expression data fram  across genes
#' and patient samples
#' @param target_gene the name of the target gene of interest
#' @param meth   
                      legend = "bottom", legend.labs = c("No Mutation", "Mutation"), 
                      font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                      font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                      font.legend = c(12, "bold", "black")) + 
        xlab("Time (Days)") + ylab("Survival Probability")
      print(p)
      
      # Cox proportional hazards
    } else {
      formula <- "Surv(time = days, event = status) ~ mutation_status + Age + Gender + Prior_malig + Treatment_rad + Treatment_pharm + PC1 + PC2 + PC3 + Tumor_purity_b1 + Tumor_purity_b2 + Tumor_purity_b3 + Tot_IC_Frac_b1 + Tot_IC_Frac_b2 + Tot_IC_Frac_b3"
      formula <- paste(formula, paste(
        colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
        collapse = " + "), sep = " + ")
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
#' @param expression_df a quantile-normalized expression data fram  across genes
#' and patient samples
#' @param target_gene the name of the target gene of interest
#' @param meth   
                      legend = "bottom", legend.labs = c("No Mutation", "Mutation"), 
                      font.main = c(16, "bold", "black"), font.x = c(14, "bold", "black"), 
                      font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), 
                      font.legend = c(12, "bold", "black")) + 
        xlab("Time (Days)") + ylab("Survival Probability")
      print(p)
      
      # Cox proportional hazards
    } else {
      formula <- "Surv(time = days, event = status) ~ mutation_status + Age + Gender + Prior_malig + Treatment_rad + Treatment_pharm + PC1 + PC2 + PC3 + Tumor_purity_b1 + Tumor_purity_b2 + Tumor_purity_b3 + Tot_IC_Frac_b1 + Tot_IC_Frac_b2 + Tot_IC_Frac_b3"
      formula <- paste(formula, paste(
        colnames(survival_df)[grepl("Cancer_type", colnames(survival_df))], 
        collapse = " + "), sep = " + ")
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
#' @param expression_df a quantile-normalized expression data framace="bold"), 
          legend.text = element_text(size=14), legend.position = "bottom") 
  print(g)
  
  
  # Chi-squared googness-of-fit test
  
  # Difference in overall response vs. disease?
  print(chisq.test(c(sum(input_df[(input_df$Response.Aggr == "Response") & (input_df$Mutation.Status == 0), 'Frequency']),
                     sum(input_df[(input_df$Response.Aggr == "Response") & (input_df$Mutation.Status == 1), 'Frequency'])), 
                   p = c(0.5, 0.5))$p.value)
  print(chisq.test(c(sum(input_df[(input_df$Response.Aggr == "Disease") & (input_df$Mutation.Status == 0), 'Frequency']),
                     sum(input_df[(input_df$Response.Aggr == "Disease") & (input_df$Mutation.Status == 1), 'Frequency'])), 
                   p = c(0.5, 0.5))$p.value)
  print(chisq.test(c(sum(input_df[(input_df$Response.Aggr == "Response") & (input_df$Mutation.Status == 1), 'Frequency']),
                     sum(input_df[(input_df$Response.Aggr == "Disease") & (input_df$Mutation.Status == 1), 'Frequency'])), 
                   p = c(0.5, 0.5))$p.value)
  
  # Fisher's exact test
  conting_table_response <- data.frame("No.Mut" = sum(input_df[(input_df$Response.Aggr == "Response") & (input_df$Mutation.Status == 0), 'Frequency']),
                                       "Mut" = sum(input_df[(input_df$Response.Aggr == "Response") & (input_df$Mutation.Status == 1), 'Frequency']))
  conting_table_disease <- data.frame("No.Mut" = sum(input_df[(input_df$Response.Aggr == "Disease") & (input_df$Mutation.Status == 0), 'Frequency']),
                                      "Mut" = sum(input_df[(input_df$Response.Aggr == "Disease") & (input_df$Mutation.Status == 1), 'Frequency']))
  conting_table <- rbind(conting_table_response, conting_table_disease)
  rownames(conting_table) <- c("Response", "Disease")
  
  test <- fisher.test(conting_table)
  print("Fisher p-value:")
  print(test$p.value)   #0.0536
  
  return(input_df)
}

vin_input_df <- create_drug_response_barplot(mutation_count_matrix_tp53, "vin.response", "\nTUBB-Targeting", "TP53")
vin_input_df <- create_drug_response_barplot(mutation_count_matrix_tp53, "vin.response", "\nTUBB-Targeting", "TP53")
etop_input_df <- create_drug_response_barplot(mutation_count_matrix_tp53, "etop.response", "\nEtoposide", "TP53")
irino_input_df <- create_drug_response_barplot(mutation_count_matrix_idh1, "irino.response", "\nIrinotexan", "IDH1")