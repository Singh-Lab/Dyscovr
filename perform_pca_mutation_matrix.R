#################################################################
### TAKE PCA OF MUTATION MATRIX
### Written By: Sara Geraghty, March 2022
#################################################################

library(tidyverse)
library(ggplot2)

mut_count_mat <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_nonsense_CancerOnly_IntersectPatients.csv", 
                          header= TRUE, check.names = FALSE, row.names = 1)
#mut_count_mat <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_nonsense_CancerOnly_IntersectPatients.csv", 
                          #header= TRUE, check.names = FALSE, row.names = 1)
#mut_count_mat <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Linear Model/Tumor_Only/GeneTarg_Mutation/mut_count_matrix_missense_nonsense_BRCA_BLCA_HNSC_CancerOnly_IntersectPatients.csv", 
                          #header= TRUE, check.names = FALSE, row.names = 1)



rownames(mut_count_mat) <- mut_count_mat$Gene_Symbol
#mut_count_mat <- mut_count_mat[mut_count_mat$Gene_Symbol != "TP53",]  # to exclude TP53, or other geme

mutation_regprot_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/BRCA/Linear Model/Tumor_Only/Regprot_Mutation/iprotein_results_missense_nonsense_CancerOnly_IntersectPatients.csv", 
                                header = TRUE, check.names = FALSE)
#mutation_regprot_df <- read.csv("C:/Users/sarae/Documents/Mona Lab Work/Main Project Files/Saved Output Data Files/Pan-Cancer/Linear Model/Tumor_Only/Regprot_Mutation/iprotein_results_missense_nonsense_CancerOnly_IntersectPatients.csv", 
                                #header = TRUE, check.names = FALSE)

unique_patients <- unique(unlist(lapply(mutation_regprot_df$Patient, function(x) unlist(strsplit(x, ";", fixed = TRUE)))))

mut_count_mat <- mut_count_mat[,c(1, (which(colnames(mut_count_mat)[2:ncol(mut_count_mat)] %in%
                                             unique_patients)+1))]

mut_count_mat_t <- t(mut_count_mat[,2:ncol(mut_count_mat)])
mut_count_mat_t <- apply(mut_count_mat_t, 2, as.numeric)
rownames(mut_count_mat_t) <- colnames(mut_count_mat[2:ncol(mut_count_mat)])

# Perform TCA
pca_t <- prcomp(mut_count_mat_t)
pc1.5 <- pca_t$x[,1:5]
pc1.5 <- as.data.frame(pc1.5)
pc1.5$PIK3CA.Mut <- unlist(lapply(rownames(pc1.5), function(x) 
  mut_count_mat_t[ rownames(mut_count_mat_t) == x, 
                          colnames(mut_count_mat_t) == "PIK3CA"]))
pc1.5$TP53.Mut <- unlist(lapply(rownames(pc1.5), function(x) 
  mut_count_mat_t[ rownames(mut_count_mat_t) == x, 
                   colnames(mut_count_mat_t) == "TP53"]))
pc1.5$Name <- rownames(pc1.5)


mut_count_mat_t_binary <- apply(mut_count_mat_t, c(1,2), function(x) {
  if(x >=1) {return(1)}
  else{return(0)}
})
mut_count_mat_t_binary <- as.data.frame(mut_count_mat_t_binary)
pca_t_bin <- prcomp(mut_count_mat_t_binary)
pc1.5_bin <- pca_t_bin$x[,1:5]
pc1.5_bin <- as.data.frame(pc1.5_bin)

pc1.5_bin$PIK3CA.Mut <- unlist(lapply(rownames(pc1.5_bin), function(x) 
  mut_count_mat_t_binary[ rownames(mut_count_mat_t_binary) == x, 
                          colnames(mut_count_mat_t_binary) == "PIK3CA"]))
pc1.5_bin$TP53.Mut <- unlist(lapply(rownames(pc1.5_bin), function(x) 
  mut_count_mat_t_binary[ rownames(mut_count_mat_t_binary) == x, 
                          colnames(mut_count_mat_t_binary) == "TP53"]))
pc1.5_bin$PIK3CA.Mut <- as.logical(pc1.5_bin$PIK3CA.Mut)
pc1.5_bin$TP53.Mut <- as.logical(pc1.5_bin$TP53.Mut)

pc1.5_bin$Name <- rownames(pc1.5_bin)

var_explained <- pca_t$sdev^2/sum(pca_t$sdev^2)

pc1.5 <- tidyr::unite(pc1.5, "TP53_PIK3CA", TP53.Mut, PIK3CA.Mut, remove = F)
pc1.5_bin <- tidyr::unite(pc1.5_bin, "TP53_PIK3CA", TP53.Mut, PIK3CA.Mut, remove = F)


pc1.5 %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2, color = TP53_PIK3CA)) + geom_point(size=2) +
  theme_bw(base_size=16) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
             y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top") +
  geom_text(aes(label= ifelse(PC2 > quantile(PC2, 0.99),
                              as.character(Name),'')),hjust=1,vjust=0) +
  geom_text(aes(label= ifelse(PC1 > quantile(PC1, 0.99),
                              as.character(Name),'')),hjust=0,vjust=0)



