########################################################################
### SURVIVAL CURVES AND DRUG RESPONSE ANALYSIS
########################################################################

# Plot survival curves and perform drug sensitivity analysis using TCGA clinical data

library(ggplot2)
library(TCGAbiolinks)
library(data.table)
library(survival)
library(survminer)

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
    mean <- mean(expr_vals)
    
    low <- expression_df_target[expression_df_target$expression < (mean - sigma1.5), ]
    high <- expression_df_target[expression_df_target$expression > (mean + sigma1.5), ]
    
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
clinical_df <- read.csv(paste0(PATH, "Input Data Files/Pan-Cancer/clinical_data_subset.csv"), 
                        header = T, check.names = F)

# Import the pan-cancer mutation count matrix
mutation_count_matrix <- read.csv(paste0(PATH, "GeneTarg_Mutation/mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv"), 
                                  header = T, check.names = F)

# Import the pan-cancer expression matrix
expression_df_qn <- read.csv(paste0(PATH, "Expression/expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU_orig.csv"), 
                             header = T, check.names = F)

# Import the combined patient/sample clinical data matrix
patient_sample_df <- read.csv(paste0(PATH, "Patient/combined_patient_sample_cibersort_total_frac_inclSubtypes_IntersectPatientsWashU.csv"),
                              header = T, check.names = F)
colnames(patient_sample_df)[1] <- "patient"

# Call function (e.g. KRAS and NT5E, or PIK3CA and KBTBD2)
kras_survival_res <- create_survival_curves(clinical_df, "KRAS", mutation_count_matrix,
                                            "ENSG00000135318", expression_df_qn, "NT5E",
                                            patient_sample_df)
pik3ca_survival_res <- create_survival_curves(clinical_df, "PIK3CA", mutation_count_matrix,
                                              "ENSG00000170852", expression_df_qn, "KBTBD2",
                                              patient_sample_df)

# Limit just to specific cancer types
pik3ca_relevant_cancer_types <- c("TCGA-BLCA", "TCGA-BRCA", "TCGA-COAD", "TCGA-ESCA", 
                                  "TCGA-GBM", "TCGA-HNSC", "TCGA-LGG", "TCGA-LUSC", 
                                  "TCGA-UCEC")
clinical_df_pik3ca_relevant <- clinical_df[clinical_df$project_id %fin% 
                                             pik3ca_relevant_cancer_types,]
pik3ca_relevant_patients <- unique(unlist(lapply(clinical_df_pik3ca_relevant$case_submitter_id, function(id)
  unlist(strsplit(id, "-", fixed = T))[3])))
mutation_patients <- unlist(lapply(colnames(mutation_count_matrix)[2:ncol(mutation_count_matrix)], function(x)
  unlist(strsplit(x, "-", fixed = T))[1]))
expression_patients <- unlist(lapply(colnames(expression_df_qn)[2:ncol(expression_df_qn)], function(x)
  unlist(strsplit(x, "-", fixed = T))[1]))

pik3ca_survival_res <- create_survival_curves(clinical_df_pik3ca_relevant, "PIK3CA", 
                                              mutation_count_matrix[c(1, (which(mutation_patients %fin% pik3ca_relevant_patients)+1))],
                                              "ENSG00000170852", 
                                              expression_df_qn[c(1, (which(expression_patients %fin% pik3ca_relevant_patients)+1))], 
                                              "KBTBD2", patient_sample_df)

# Z-score normalize the expression matrix (transform each variable value 
# by subtracting its mean and dividing by its SD, resulting in a new variable with
# mean of 0 and SD of 1, with each value representing # of SDs from the mean)
expression_df_znorm <- cbind(expression_df_qn[,1], 
                             as.data.frame(scale(expression_df_qn[,2:ncol(expression_df_qn)])))
colnames(expression_df_znorm)[1] <- "ensg_id"
hist(as.numeric(expression_df_qn[expression_df_qn$ensg_id == "ENSG00000170852", 
                                 2:ncol(expression_df_qn)]), 
     main = "", xlab = "Expression", col = "cornflowerblue")
hist(as.numeric(expression_df_znorm[expression_df_znorm$ensg_id == "ENSG00000170852", 
                                    2:ncol(expression_df_znorm)]), 
     main = "", xlab = "Expression", col = "cornflowerblue")

# Run survival analyses using z-score normalized expression, keeping patients that
# are above 1.5 sigma from the mean or below -1.5 sigma
pik3ca_survival_res <- create_survival_curves(clinical_df, "PIK3CA", mutation_count_matrix,
                                              "ENSG00000170852", expression_df_znorm, "KBTBD2",
                                              patient_sample_df, method = "znorm")
expression_patients_znorm <- unlist(lapply(colnames(expression_df_znorm)[
  2:ncol(expression_df_znorm)], function(x)
    unlist(strsplit(x, "-", fixed = T))[1]))

# Call function
pik3ca_survival_res <- create_survival_curves(clinical_df_pik3ca_relevant, "PIK3CA", 
                                              mutation_count_matrix[c(1, (which(mutation_patients %fin% pik3ca_relevant_patients)+1))],
                                              "ENSG00000170852", 
                                              expression_df_znorm[c(1, (which(expression_patients_znorm %fin% pik3ca_relevant_patients)+1))], 
                                              "KBTBD2", patient_sample_df,
                                              method = "znorm")

# Use Cox proportional hazards model
pik3ca_survival_res <- create_survival_curves(clinical_df_pik3ca_relevant, "PIK3CA", 
                                              mutation_count_matrix[c(1, (which(mutation_patients %fin% pik3ca_relevant_patients)+1))],
                                              "ENSG00000170852", 
                                              expression_df_znorm[c(1, (which(expression_patients_znorm %fin% pik3ca_relevant_patients)+1))], 
                                              "KBTBD2", patient_sample_df, 
                                              method = "znorm", surv_curve_type = "cox")
pik3ca_survival_res <- create_survival_curves(clinical_df, "PIK3CA", 
                                              mutation_count_matrix, "ENSG00000170852", 
                                              expression_df_znorm, "KBTBD2", 
                                              patient_sample_df, method = "znorm", 
                                              surv_curve_type = "cox")

# Alternatively, create a survival curve individually per cancer type
pik3ca_patients_perct <- lapply(pik3ca_relevant_cancer_types, function(ct) {
  clinical_df_ct <- clinical_df[clinical_df$project_id == ct,]
  patients <- unique(unlist(lapply(clinical_df_ct$case_submitter_id, function(id)
    unlist(strsplit(id, "-", fixed = T))[3])))
  return(patients)
})
names(pik3ca_patients_perct) <- pik3ca_relevant_cancer_types

surv_res_perct <- lapply(1:length(pik3ca_patients_perct), function(i) {
  ct <- names(pik3ca_patients_perct)[i]
  pats <- pik3ca_patients_perct[[i]]
  print(ct)
  
  try({
    surv_res <- create_survival_curves(clinical_df_pik3ca_relevant, "PIK3CA", 
                                       mutation_count_matrix[c(1, (which(mutation_patients %fin% pats)+1))],
                                       "ENSG00000170852", 
                                       expression_df_znorm[c(1, (which(expression_patients_znorm %fin% pats)+1))], 
                                       "KBTBD2", patient_sample_df, method = "znorm")
    return(surv_res)
  })
  return(NA)
})
names(surv_res_perct) <- pik3ca_relevant_cancer_types

#
#
#
#
#


### DRUG RESPONSE ###
# Data source: http://lifeome.net/supp/drug_response/
drug_response <- read.table(paste0(PATH, "Pan-Cancer/Drug_Response/drug_response.txt"), 
                            header = T, sep = "\t")

# How many patients in each cancer type were treated with drugs?
patients <- distinct(drug_response[,c(1,2)])
cancer_freq_table <- melt(table(patients$cancers))
colnames(cancer_freq_table) <- c("Cancer.Type", "Num.Unique.Patients")
cancer_freq_table <- cancer_freq_table[order(cancer_freq_table$Num.Unique.Patients, decreasing = T),]
cancer_freq_table$Num.Unique.Patients <- as.numeric(cancer_freq_table$Num.Unique.Patients)
ggplot(cancer_freq_table, aes(x = reorder(Cancer.Type, -Num.Unique.Patients), 
                              y = Num.Unique.Patients)) + 
  geom_bar(position="stack", stat="identity", fill = "#0072B5FF") + theme_minimal() +
  xlab("TCGA Cancer Type") + ylab("Number of Unique Patients") +
  theme(axis.text = element_text(face="bold", size = 16), axis.text.x = element_text(angle = 45, vjust =1, hjust=1), 
        axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14, face="bold"),
        legend.text = element_text(size=12))

# In how many cancer types was each drug used?
drugs <- distinct(drug_response[,c(1,3)])
drug_freq_table <- melt(table(drugs$drug.name))
colnames(drug_freq_table) <- c("Drug.Name", "Num.Cancer.Types")
drug_freq_table <- drug_freq_table[order(drug_freq_table$Num.Cancer.Types, decreasing = T),]
drug_freq_table$Num.Cancer.Types <- as.numeric(drug_freq_table$Num.Cancer.Types)
drug_freq_table$Drug.Name <- as.character(gsub("\\\\", ".", gsub(">", ".", gsub("<", ".", drug_freq_table$Drug.Name))))
ggplot(drug_freq_table, aes(x = reorder(Drug.Name, -Num.Cancer.Types), 
                            y = Num.Cancer.Types)) + 
  geom_bar(position="stack", stat="identity", fill = "#BC3C29FF") + theme_minimal() +
  xlab("Drug Name") + ylab("Number of Cancer Types Treated") +
  theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 45, vjust =1, hjust=1), 
        axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14, face="bold"),
        legend.text = element_text(size=12))

# Histogram of the frequency values
ggplot(drug_freq_table, aes(x = Num.Cancer.Types)) + geom_histogram(position = "identity", bins = 5, fill = "#BC3C29FF") + 
  theme_minimal() + xlab("Number of Cancer Types Treated") + ylab("Frequency") +
  theme(axis.text = element_text(size = 16, face="bold"), axis.text.x = element_text(angle = 45, vjust =1, hjust=1), 
        axis.title=element_text(size=18,face="bold"))

# Get a table of the top drugs with their cancer types
drugs_w_cancer_types <- data.table("drug.name" = unique(drugs$drug.name),
                                   "cancer.types" = unlist(lapply(unique(drugs$drug.name), function(drug)
                                     paste(unique(drugs[drugs$drug.name == drug, 'cancers']), collapse = ", "))))

# Call function
tp53_survival_res <- create_survival_curves(clinical_df, "TP53", 
                                            mutation_count_matrix, drug_response)

# Plot drug response by mutation status 
mutation_count_matrix_driver <- melt(mutation_count_matrix[mutation_count_matrix$Gene_Symbol == driver_name,])[,c("variable", "value")]
colnames(mutation_count_matrix_driver) <- c('patient', 'mutation_status')
#print(head(mutation_count_matrix_driver))
mutation_count_matrix_driver$patient <- unlist(lapply(mutation_count_matrix_driver$patient, function(x)
  unlist(strsplit(as.character(x), "-", fixed = T))[1]))
mutation_count_matrix_driver$mutation_status <- unlist(lapply(mutation_count_matrix_driver$mutation_status, function(x)
  ifelse(x > 0, 1, 0)))

# Merge with drug treatment and response information
drug_response_cancer <- drug_response[drug_response$cancers %in% c("LGG", "LUAD"),]
drug_response_drug <- drug_response[drug_response$drug.name == "Pemetrexed",]

drug_response_drug <- drug_response[drug_response$drug.name %fin% c("Vinblastine", "Vincristine", "Vinorelbine"),]

mutation_count_matrix_driver$vin.response <- unlist(lapply(mutation_count_matrix_driver$patient, function(x) {
  print(x)
  val <- unique(drug_response_drug[grepl(x, drug_response_drug$patient.arr), 'response'])
  if(length(val) == 0) {return(NA)}
  if(length(val) > 1) {return(val[1])}
  else{return(val)}
}))

#' Create a bar plot for patient drug response by mutation status and conduct a 
#' Chi-squared goodness of fit test and Fisher's exact test
#' @param mutation_count_matrix_driver driver mutation status for each sample
#' @param col_name the name of the mutation count matrix that has the drug response information
#' @param lab a label with the name of the drug
#' @param goi the name of the driver gene-of-interest
create_drug_response_barplot <- function(mutation_count_matrix_driver, col_name, lab, goi) {
  input_df_nomut <- melt(table(mutation_count_matrix_driver[mutation_count_matrix_driver$mutation_status == 0,
                                                            col_name]))
  colnames(input_df_nomut) <- c("Response", "Frequency")
  input_df_nomut$Mutation.Status <- rep(0, times = nrow(input_df_nomut))
  input_df_mut <- melt(table(mutation_count_matrix_driver[mutation_count_matrix_driver$mutation_status == 1,
                                                          col_name]))
  colnames(input_df_mut) <- c("Response", "Frequency")
  input_df_mut$Mutation.Status <- rep(1, times = nrow(input_df_mut))
  input_df <- rbind(input_df_nomut, input_df_mut)
  input_df$Response <- factor(input_df$Response, levels = c("Complete Response", "Partial Response",
                                                            "Stable Disease", "Clinical Progressive Disease"))
  input_df$Response.Aggr <- unlist(lapply(input_df$Response, function(x) ifelse(grepl("Response", x), "Response", "Disease")))
  input_df$Response.Aggr <- factor(input_df$Response.Aggr, levels = c("Response", "Disease"))
  input_df$Response.Aggr <- as.character(input_df$Response.Aggr)
  
  g <- ggplot(input_df, aes(x = Response.Aggr, y = Frequency, fill = as.factor(Mutation.Status))) +
    geom_bar(stat = "identity", position = position_dodge()) +
    xlab(paste(lab, "Drug Response")) + ylab("Frequency") + labs(fill = paste(goi, "Mutation Status")) +
    theme_minimal() + scale_fill_nejm() +
    theme(axis.text.x=element_text(size=14), axis.title.x=element_text(size=14, face="bold"), 
          axis.title.y=element_text(size=14, face="bold"), axis.text.y=element_text(size=14),
          legend.title = element_text(size=14, face="bold"), 
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