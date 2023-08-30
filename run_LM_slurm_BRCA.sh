#!/bin/bash
#SBATCH --qos=1day
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=54GB
#SBATCH --job-name=lm_topDrivers
#SBATCH --mail-user=scamilli@princeton.edu
#SBATCH --mail-type=fail,time_limit

module add R/4.0.3

# SAMPLE RUN: TP53 and all targets/ ChIP-eat targets (NON-TNM)
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "TRUE" --incl_nextMutDriver "TRUE" --expression_df "expression_quantile_norm_edgeRfilt_SDGr1_inverseNtransform_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatients.csv" --mutation_targ_df "iprotein_mut_count_matrix_nonsynonymous_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_nonsynonymous_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_tmm_normAge_IntersectPatients.csv" --target_df "glycolysis_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --targets_name "glycolysisTargs" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 0 --mut_pc_vect "0" --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;PIK3CA_MutStat_i;PIK3CA_CNAStat_i;PIK3CAP53_MutStat_i;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "PIK3CA_covs.inclCT.invNTr.Nonsyn" --removeCis "TRUE" --dataset "TCGA" --inclResiduals "FALSE"

#--patientsOfInterest "LT20PAmp.NotTP53andPIK3CAMNmut_patient_ids.txt" --patientsOfInterestLabel "LT20PAmp.noMutOverlap" \
#--patientsOfInterest "HER2_patient_ids.txt" --patientsOfInterestLabel "HER2" \
#--patientsOfInterest "LessThan20PercAmp_patient_ids.txt" --patientsOfInterestLabel "LT20PAmp" \


# SAMPLE RUN: PIK3CA and all targets/ ChIP-eat targets (NON-TNM)
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "TRUE" --tester_name "PIK3CA" --tester_uniprot_id "P42336" --tester_ensg_id "ENSG00000121879" --incl_nextMutDriver "TRUE" \
#--expression_df "expression_tmmSDGr1filtByExpr_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatients.csv" --mutation_targ_df "iprotein_mut_count_matrix_nonsynonymous_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_nonsynonymous_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_tmm_normAge_IntersectPatients.csv" --target_df "metabolic_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --targets_name "metabolicTargs" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 0 --mut_pc_vect "0" --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;TP53_MutStat_i;TP53_CNAStat_i;TP53PIK3CA_MutStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "TP53_Covs.inclCT.Nonsyn" --removeCis "TRUE" --dataset "TCGA"

#--functional_copies_df "pik3ca_functional_copies_per_sample.csv"
#--patientsOfInterest "LessThan20PercAmp_patient_ids.txt" --patientsOfInterestLabel "LT20PAmp" \

# SAMPLE RUN: IDH1 and metabolic targets (NON-TNM)
#Rscript linear_model.R --QTLtype "meQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "TRUE" --tester_name "IDH1" --tester_uniprot_id "O75874" --tester_ensg_id "ENSG00000138413" \
#--expression_df "expression_tmmSDGr1filtByExpr_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatients.csv" --mutation_targ_df "mut_count_matrix_missense_nonsense_CancerOnly_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_missense_nonsense_CancerOnly_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_tmm_normAge_IntersectPatients.csv" --target_df "allgene_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --targets_name "allGenes" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 0 --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "inclCT" --removeCis "TRUE" --dataset "TCGA"

# SAMPLE RUN: Other Genes and all targets (NON-TNM) 
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "TRUE" --tester_name "TTN" --tester_uniprot_id "Q8WZ42" --tester_ensg_id "ENSG00000155657" --incl_nextMutDriver "TRUE" \
#--expression_df "expression_quantile_norm_edgeRfilt_SDGr1_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatients.csv" --mutation_targ_df "mut_count_matrix_missense_nonsense_CancerOnly_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_missense_nonsense_CancerOnly_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_tmm_normAge_IntersectPatients.csv" --target_df "allgene_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --targets_name "allGenes" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 2 --mut_pc_vect "0" --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;TP53_MutStat_i;TP53_CNAStat_i;PIK3CA_MutStat_i;PIK3CA_CNAStat_i;TP53PIK3CA_MutStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "TP53_Covs.inclCT.MN.edgeRfilt.SDGr1" --removeCis "TRUE" --dataset "TCGA"

# SAMPLE RUN: Other Genes and all targets (NON-TNM) 
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" --removeMetastatic "TRUE" --incl_nextMutDriver "TRUE" \
#--patientsOfInterest "Basal_patient_ids.txt" --patientsOfInterestLabel "Basal" \
#--test "TRUE" --tester_name "MUC16" --tester_uniprot_id "Q8WXI7" --tester_ensg_id "ENSG00000181143" \
#--expression_df "expression_quantile_norm_edgeRfilt_SDGr1_sklearn_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatients.csv" --mutation_targ_df "iprotein_mut_count_matrix_nonsynonymous_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_nonsynonymous_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_tmm_normAge_IntersectPatients.csv" --target_df "metabolic_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --targets_name "metabolicTargs" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 0 --mut_pc_vect "0" --debug "TRUE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;TP53_MutStat_i;TP53_CNAStat_i;TP53PIK3CA_MutStat_i;PIK3CA_MutStat_i;PIK3CA_CNAStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "TP53_PIK3CA_Covs.inclCT.sklearn.Nonsyn" --removeCis "TRUE" --dataset "TCGA" --run_query_genes_jointly "FALSE"

# RBP SF3B1 (mutated in 11); USH2A is highly mutated but not cancer relevant (negative control)
#--test "TRUE" --tester_name "SF3B1" --tester_uniprot_id "O75533" --tester_ensg_id "ENSG00000115524" \
#--test "TRUE" --tester_name "USH2A" --tester_uniprot_id "O75445" --tester_ensg_id "ENSG00000042781" \
#--test "TRUE" --tester_name "MYC" --tester_uniprot_id "P01106" --tester_ensg_id "ENSG00000136997" \
#--test "TRUE" --tester_name "KMT2C" --tester_uniprot_id "Q8NEZ4" --tester_ensg_id "ENSG00000055609" \
#--test "TRUE" --tester_name "GATA3" --tester_uniprot_id "P23771" --tester_ensg_id "ENSG00000107485"  \
#--test "TRUE" --tester_name "FOXA1" --tester_uniprot_id "P55317" --tester_ensg_id "ENSG00000129514" \
#--test "TRUE" --tester_name "MUC16" --tester_uniprot_id "Q8WXI7" --tester_ensg_id "ENSG00000181143" \

#--patientsOfInterest "LessThan20PercAmp_patient_ids.txt" --patientsOfInterestLabel "LT20PAmp" \

# SAMPLE RUN: Cancer Related Genes/ Highly Deleted TFs and Metabolic Targets (NON-TNM); with regularization
#Rscript linear_model2.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "FALSE" --protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_drivers_nonsynonymous.csv" --incl_nextMutDriver "FALSE" \
#--expression_df "expression_quantile_norm_edgeRfilt_SDGr1_sklearn_CancerOnly_IntersectPatientsWashU.csv" \
#--methylation_df "methylation_M_CancerOnly_IntersectPatientsWashU.csv" --methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatientsWashU.csv" \
#--mutation_targ_df "iprotein_mut_count_matrix_nonsynonymous_IntersectPatientsWashU.csv" --mutation_regprot_df "iprotein_results_nonsynonymous_IntersectPatientsWashU.csv" \
#--cna_df "CNA_AllGenes_CancerOnly_IntersectPatientsWashU.csv" --patient_df "combined_patient_sample_cibersort_total_frac_tmm_normAge_IntersectPatientsWashU.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "top_drivers" \
#--targets_name "allGenes" --num_PEER 0 --num_pcs 3 --debug "TRUE" --collinearity_diagn "FALSE" --regularization "bayesian.bglss" --useNumFunctCopies "FALSE" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "inclCT.sklearn.Nonsyn.0.05.WashU" --removeCis "TRUE" --dataset "TCGA" --run_query_genes_jointly "TRUE" --model_type "linear" 


# SAMPLE RUN: Cancer Related Genes/ Highly Deleted TFs and Metabolic Targets (NON-TNM); without regularization
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "FALSE" --protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_drivers_nonsynonymous.csv" --incl_nextMutDriver "FALSE" \
#--expression_df "expression_quantile_norm_edgeRfilt_SDGr1_sklearn_CancerOnly_IntersectPatientsWashU.csv" \
#--methylation_df "methylation_M_CancerOnly_IntersectPatientsWashU.csv" --methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatientsWashU.csv" \
#--mutation_targ_df "iprotein_mut_count_matrix_nonsynonymous_IntersectPatientsWashU.csv" --mutation_regprot_df "iprotein_results_nonsynonymous_IntersectPatientsWashU.csv" \
#--cna_df "CNA_AllGenes_CancerOnly_IntersectPatientsWashU.csv" --patient_df "combined_patient_sample_cibersort_total_frac_tmm_normAge_IntersectPatientsWashU.csv" \
#--target_df "metabolic_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "top_drivers" \
#--targets_name "metabolicTargs" --num_PEER 0 --num_pcs 3 --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" --useNumFunctCopies "FALSE" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "inclCT.sklearn.Nonsyn.0.05.WashU.original" --removeCis "TRUE" --dataset "TCGA" --run_query_genes_jointly "TRUE" --model_type "linear" 

#--patientsOfInterest "Luminal.A_patient_ids.txt" --patientsOfInterestLabel "LumA" \
#--patientsOfInterest "LessThan20PercAmp_patient_ids.txt" --patientsOfInterestLabel "LT20PAmp" \
#--patientsOfInterest "NotTP53andPIK3CAmut_patient_ids.txt" --patientsOfInterestLabel "NoMutOverlap" \


# SAMPLE RUN: Cancer Related Genes/ Highly Deleted TFs and Metabolic Targets (NON-TNM)
Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --run_query_genes_jointly "TRUE" \
--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_edgeRfilt_sklearn_CancerOnly_IntersectPatientsWashU.csv" \
--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/BRCA/tumor_only/lm_input_tables/top_drivers_0.05_relCNA/eQTL" \
--patient_df "combined_patient_sample_cibersort_total_frac_IntersectPatientsWashU.csv" --run_name "top_drivers_0.05_relCNA" \
--target_df "allgene_targets.csv" --targets_name "allGenes" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
--num_PEER 0 --num_pcs 3 --debug "TRUE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
--select_args_label "inclCT.sklearn.Nonsyn.0.05.WashU" 

#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/BRCA/tumor_only/lm_input_tables/top_drivers_0.025/eQTL" \
#--patientsOfInterest "HER2_patient_ids.txt" --patientsOfInterestLabel "HER2" \
#--mutation_regprot_df "iprotein_results_nonsynonymous_IntersectPatientsWashU.csv" \


# FILES THAT WE HAVE TO CHOOSE FROM:
###########################################################################################
# Protein IDs DF
###########################################################################################
# I-Protein: "iprotein_protein_ids_df.csv"
# I-Protein (Nucleic Acids): "iprotein_nucacids_protein_ids_df.csv"
# I-Protein (Cancer-Related Genes): "iprotein_protein_ids_df_cancerRelated.csv"
# I-Protein (Highly Deleted TFs): "iprotein_protein_ids_df_highlyDelTFs.csv"

# I-Domain: "idomain_protein_ids_df.csv"
# I-Domain (Nucleic Acids): "idomain_nucacids_protein_ids_df.csv"

# I-Binding Position: "ibindingpos_protein_ids_df.csv"
# I-Binding Position (Nucleic Acids)" "ibindingpos_nucacids_protein_ids_df.csv"


###########################################################################################
# Expression DF
###########################################################################################
# FPKM: "expression_fpkm_CancerOnly_IntersectPatients.csv"  # fpkm-normalized and filtered
# TMM: "expression_tmm_CancerOnly_IntersectPatients.csv"
# Quantile-Normalized: "expression_quantile_norm_CancerOnly_IntersectPatients.csv"
# Rank-Normalized: "expression_rank_norm_CancerOnly_IntersectPatients.csv"

###########################################################################################
# Methylation DF
###########################################################################################
# Beta Values, Cancer Only: "methylation_Beta_CancerOnly_IntersectPatients.csv"
# M-Values, Cancer Only: "methylation_M_CancerOnly_IntersectPatients.csv"

# Beta (Bucketed, 3 Buckets), Cancer Only: "methylation_bucketed_Beta_CancerOnly_IntersectPatients.csv"
# M (Bucketed, 3 Buckets), Cancer Only: "methylation_bucketed_M_CancerOnly_IntersectPatients.csv"

# Threshold (Beta > 0.8), Cancer Only: "methylation_Beta_0.8_CancerOnly_IntersectPatients.csv"
# Threshold (Beta > 0.5), Cancer Only: "methylation_Beta_0.5_CancerOnly_IntersectPatients.csv"

###########################################################################################
# Mutation Target DF
###########################################################################################
# Missense: "mut_count_matrix_missense_IntersectPatients.csv"

###########################################################################################
# Regulatory Protein Mutation DF
###########################################################################################
# I-Protein: "iprotein_results_missense_IntersectPatients.csv"
# I-Protein (Nucleic Acids): "iprotein_results_missense_nucacids_IntersectPatients.csv"
# I-Domain: "idomain_results_missense_IntersectPatients.csv"
# I-Domain (Nucleic Acids): "idomain_results_missense_nucacids_IntersectPatients.csv"
# I-Binding Position: "ibindingpos_results_missense_IntersectPatients.csv"
# I-Binding Position (Nucleic Acids): "ibindingpos_results_missense_nucacids_IntersectPatients.csv"

###########################################################################################
# CNA DF
###########################################################################################
# All Genes: "CNA_AllGenes_IntersectPatients.csv"
# All Genes (Cancer Samples Only): "CNA_AllGenes_CancerOnly_IntersectPatients.csv"
# All Genes (Cancer Samples Only), Bucketed, Including Amplifications: "CNA_AllGenes_Bucketed_InclAmp_CancerOnly_IntersectPatients.csv"
# TFs Only: "CNA_TFsOnly_IntersectPatients.csv"

###########################################################################################
# Patient DF
###########################################################################################
# CIBERSORT Abs I.C.D.: 
# 1. "combined_patient_sample_cibersort_abs_tmm_IntersectPatients.csv"
# 2. "combined_patient_sample_cibersort_abs_fpkm_IntersectPatients.csv"
# 3. "combined_patient_sample_cibersort_abs_fpkm_top10k_IntersectPatients.csv"
# 4. "combined_patient_sample_cibersort_abs_qn_top10k_IntersectPatients.csv"
# 5. "combined_patient_sample_cibersort_abs_rn_IntersectPatients.csv"
# 6. "combined_patient_sample_cibersort_abs_rn_top10k_IntersectPatients.csv"

# TIMER I.C.D.: 
# 1. "combined_patient_sample_timer_tmm_IntersectPatients.csv"
# 2. "combined_patient_sample_timer_fpkm_IntersectPatients.csv"
# 3. "combined_patient_sample_timer_fpkm_top10k_IntersectPatients.csv"
# 4. "combined_patient_sample_timer_qn_top10k_IntersectPatients.csv"
# 5. "combined_patient_sample_timer_rn_IntersectPatients.csv"
# 6. "combined_patient_sample_timer_rn_top10k_IntersectPatients.csv"

# CIBERSORT Total Frac. I.C.D.: 
# 1. "combined_patient_sample_cibersort_total_frac_tmm_IntersectPatients.csv"
# 2. "combined_patient_sample_cibersort_total_frac_fpkm_IntersectPatients.csv"
# 3. "combined_patient_sample_cibersort_total_frac_fpkm_top10k_IntersectPatients.csv"
# 4. "combined_patient_sample_cibersort_total_frac_qn_top10k_IntersectPatients.csv"
# 5. "combined_patient_sample_cibersort_total_frac_rn_IntersectPatients.csv"
# 6. "combined_patient_sample_cibersort_total_frac_rn_top10k_IntersectPatients.csv"

###########################################################################################
# Sample Protein Info
###########################################################################################
#sample_protein_uniprot <- "P04637" 
#sample_protein_ensg <- "ENSG00000141510"
#sample_protein_name <- "P53"

#sample_protein_uniprot <- "P55317" 
#sample_protein_ensg <- "ENSG00000129514"
#sample_protein_name <- "FOXA1"

# Negative Controls
#sample_protein_uniprot <- "O75445"
#sample_protein_ensg <- "ENSG00000042781"
#sample_protein_name <- "USH2A"

#sample_protein_uniprot <- "Q92736"
#sample_protein_ensg <- "ENSG00000198626"
#sample_protein_name <- "RYR2"

###########################################################################################
# Targets DF
###########################################################################################
# Option 1: Curated Targets
# TP53 Curated Targets: "tp53_curated_targets.csv"
# FOXA1 Curated Targets: "foxa1_curated_targets.csv"

# Option 2: Target Lists (Generic)
# Oncogenes and Tumor Suppressors Retrieved from: https://cancerres.aacrjournals.org/content/canres/suppl/2012/01/23/0008-5472.CAN-11-2266.DC1/T3_74K.pdf
# Oncogenes: "oncogene_targets.csv"
# Tumor-suppressors: "tumor_suppressor_targets.csv"
# Metabolic targets: "metabolic_targets.csv"
# All gene targets: "allgene_targets.csv"

# Option 2: ChIP-eat Targets
# TP53 ChIP-eat Targets: "tp53_chipeat_targets.csv"
# FOXA1 ChIP-eat Targets: "foxa1_chipeat_targets.csv"

# Option 3: All Gene Targets: "allgene_targets.csv"

# Option 4: Metabolic Targets: "metabolic_targets.csv"
