############################################################
# Code to create Suppl. Figure 5 Visualizations
# Written by Sara Geraghty
# PUBLICATION INFORMATION
############################################################

library(data.table)
library(ggplot2)
library(ggsci)

############################################################
### IMPORT PAN-CANCER OUTPUT FILE(S)
############################################################
outfn <- "res_top_0.05_allGenes_quantile_rawCNA_methMRaw_3PCs_Nonsyn.Drivers.Vogel.elim.vif5.sp0.7_corrected_MUT.csv"
pc_allGenes <- read.csv(paste0(PATH, outfn))

# Associated eliminated variables from multicollinearity
elim_vars_fn <- "res_top_0.05_orig_allGenes_quantile_rawCNA_methMRaw_3PCs_elim.vif5.sp0.7_uncorrected_eliminated_variables.csv"
pc_elimVars <- read.csv(paste0(PATH, elim_vars_fn))

# Full output, with all covariates
outfn_full <- "res_top_0.05_allGenes_quantile_rawCNA_methMRaw_3PCs_Nonsyn.Drivers.Vogel.elim.vif5.sp0.7_uncorrected.csv"
pc_allGenes_full <- read.csv(paste0(PATH, outfn_full))


############################################################
### PART A: FREQUENCY OF TOP COVARIATES
############################################################
# Remove the (Intercept) terms
pc_allGenes_full <- pc_allGenes_full[pc_allGenes_full$term != "(Intercept)",]

# Of the top X remaining tests, get what the terms are and plot the proportions
x <- 500
pc_allGenes_full_topx <- pc_allGenes_full[1:x,]
terms <- unique(pc_allGenes_full$term)
terms_counts <- unlist(lapply(terms, function(x) 
  nrow(pc_allGenes_full_topx[pc_allGenes_full_topx$term == x,])))
terms_counts_df <- data.frame('Term' = terms, 'Freq' = terms_counts)

# Visualize this as a Barplot 
ggplot(data = terms_counts_df[1:5,], aes(x = reorder(Term, -Freq), 
                                         y = Freq, fill = Term)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  theme_minimal() + scale_fill_nejm() +
  xlab("Covariate") + ylab(paste0("\n", paste0("Frequency Among Top ", x))) +
  theme(axis.text = element_text(size = 16, face = "bold"), 
        axis.text.x = element_text(angle = 45, vjust =1, hjust=1), 
        axis.title=element_text(size=18,face="bold"), 
        legend.title=element_text(size=14, face="bold"),
        legend.text = element_text(size=12))


############################################################
### PART B: FREQUENCY VARIABLES REMOVED FROM COLLINEARITY ANALYSIS (VIF)
############################################################
#' Uses the files output from the correct_collinearity functions to calculate 
#' the number of times each variable was removed due to collinearity
#' @param elim_vars_df a data frame with the variables removed from the 
#' regression for each target gene
calculate_collinearity_removal_freq <- function(elim_vars_df) {
  elim_vars <- unlist(lapply(elim_vars_df$Variables.Removed, function(x) 
    unlist(strsplit(x, ",", fixed = T))))
  elim_vars_table <- table(elim_vars)
  
  # Adjust table
  elim_vars_table_sub <- as.data.frame(elim_vars_table[
    !grepl("Cancer_type", names(elim_vars_table))])
  elim_vars_table_sub$Freq <- as.numeric(elim_vars_table_sub$Freq)
  elim_vars_table_sub$elim_vars <- as.factor(elim_vars_table_sub$elim_vars)
  
  # Make more readable
  elim_vars_table_sub$elim_vars <- unlist(lapply(elim_vars_table_sub$elim_vars, function(x) {
    x <- as.character(x)
    if(grepl("_k", x)) {return(paste0("Target.", unlist(strsplit(x, "_", fixed = T))[1]))}
    else if(grepl("_i", x)) {
      spl <- unlist(strsplit(x, "_", fixed = T))
      d <- unique(pc_allGenes[pc_allGenes$R_i == spl[1], 'R_i.name'])
      return(paste(d, spl[2], sep = "."))
    } else {return(gsub("_", ".", x))}
  }))
  
  return(elim_vars_table_sub)
}

# Call function
elim_vars_table <- calculate_collinearity_removal_freq(pc_elimVars)

# Create barplot
ggplot(elim_vars_table, aes(x = reorder(elim_vars, -Freq), y = Freq)) + 
  geom_bar(position = "dodge", width = 0.95, stat = "identity", 
           show.legend = FALSE, color = "black", fill = "#20854EFF") + 
  xlab("Eliminated Variable") + ylab("Frequency") +
  theme(axis.text = element_text(face="bold", size = 16, angle = 45, 
                                 vjust = 0.5, hjust = 1), 
        axis.title=element_text(size=18, face="bold"), 
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = 'white'))
