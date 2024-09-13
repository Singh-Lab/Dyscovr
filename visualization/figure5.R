############################################################
# Code to Create Figure 5 Visualizations
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


########################################################################
### PART B: DEPMAP DRUG SENSITIVITY SPEARMAN WITH KBTBD2 KO 
########################################################################
# Correlate drug sensitivity to target gene dependency 
drug_sensitivity <- read.csv(paste0(PATH, "DepMap/Drug_sensitivity_(PRISM_Re
                                    purposing_Primary_Screen).csv"),
                             header = T, check.names = F)
# AUC(CTD^2) -- better quality and more cell lines
drug_sensitivity <- read.csv(paste0(PATH, "DepMap/Drug_sensitivity_AUC_(CTD^2).csv"),
                             header = T, check.names = F)
colnames(drug_sensitivity)[1] <- "depmap_id"

# Annotate file with GOI dependency for each cell line
goi <- "KBTBD2"
crispr_sub <- crispr[crispr$gene == goi, c("depmap_id", "value")]
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
### PART F: GROWTH CURVES (USING DATA FROM TXT FILE)
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
