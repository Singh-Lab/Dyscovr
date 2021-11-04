#!/bin/bash
#SBATCH --mem=8192
#SBATCH --qos=1wk
#SBATCH --nodes=1
#SBATCH --job-name=lm_run_cancerRelated_noPEER_metabolicPWtargs_tmm_rawCNA_M_eQTL_allButCancerType
#SBATCH --mail-user=scamilli@princeton.edu
#SBATCH --mail-type=fail,time_limit


#SAMPLE RUN: TP53 and all targets (TNM)
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "TRUE" --cancerType "BRCA" --randomize "FALSE" \
#--test "TRUE" --expression_df "expression_tmm_DF_tnm.csv" --methylation_df "methylation_DF_Beta_tnm.csv" \
#--mutation_targ_df "mut_count_matrix_missense_tnm.csv" --mutation_regprot_df "iprotein_results_missense_tnm.csv" \
#--cna_df "CNV_DF_AllGenes_tnm.csv" --patient_df "combined_patient_sample_DF_cibersort_total_frac_tmm_tnm.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "Beta" \
#--targets_name "allGenes" --num_PEER 10 --num_pcs 0 --debug "FALSE" --collinearity_diagn "TRUE" --regularization "None"

# SAMPLE RUN: TP53 and all targets/ ChIP-eat targets (NON-TNM)
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" \
#--test "TRUE" --expression_df "expression_tmm_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatients.csv" --mutation_targ_df "mut_count_matrix_missense_CancerOnly_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_missense_CancerOnly_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_tmm_normAge_IntersectPatients.csv" --target_df "allgene_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --targets_name "allGenes" --num_PEER 0 \
#--num_pcs 0 --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut" \
#--select_args_label "allButCancerType"

# SAMPLE RUN: IDH4 and metabolic targets (NON-TNM)
#Rscript linear_model.R --QTLtype "meQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" \
#--test "TRUE" --tester_name "IDH1" --tester_uniprot_id "O75874" --tester_ensg_id "ENSG00000138413" \
#--expression_df "expression_tmm_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatients.csv" --mutation_targ_df "mut_count_matrix_missense_CancerOnly_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_missense_CancerOnly_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_tmm_normAge_IntersectPatients.csv" --target_df "metabolic_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --targets_name "metabolicTargs" --num_PEER 0 \
#--num_pcs 0 --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut" \
#--select_args_label "allButCancerType"


# SAMPLE RUN: USH2A and all targets (NON-TNM) -- NEGATIVE CONTROL
#Rscript linear_model.R --QTLtype "meQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" \
#--test "TRUE" --tester_name "USH2A" --tester_uniprot_id "O75445" --tester_ensg_id "ENSG00000042781" \
#--expression_df "expression_tmm_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_bucketed_Beta_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_Beta_CancerOnly_IntersectPatients.csv" --mutation_targ_df "mut_count_matrix_missense_CancerOnly_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_missense_CancerOnly_IntersectPatients.csv" --cna_df "CNA_AllGenes_Bucketed_InclAmp_CancerOnly_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_tmm_IntersectPatients.csv" --target_df "allgene_targets.csv" \
#--cna_bucketing "bucket_inclAmp" --meth_bucketing "TRUE" --meth_type "Beta" --targets_name "allGenes" --num_PEER 0 --num_pcs 0 \
#--debug "FALSE" --collinearity_diagn "TRUE" --regularization "None" --select_args "MutStat_i;ExpStat_k"

# SAMPLE RUN: Cancer Related Genes/ Highly Deleted TFs and Metabolic Targets (NON-TNM)
Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "TRUE" \
--test "FALSE" --protein_ids_df "iprotein_protein_ids_df_cancerRelated.csv" --expression_df "expression_tmm_CancerOnly_IntersectPatients.csv" \
--methylation_df "methylation_M_CancerOnly_IntersectPatients.csv" --methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatients.csv" \
--mutation_targ_df "mut_count_matrix_missense_CancerOnly_IntersectPatients.csv" --mutation_regprot_df "iprotein_results_missense_CancerOnly_IntersectPatients.csv" \
--cna_df "CNA_AllGenes_CancerOnly_IntersectPatients.csv" --patient_df "combined_patient_sample_cibersort_total_frac_tmm_normAge_IntersectPatients.csv" \
--target_df "metabolic_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "cancerRelated" \
--targets_name "metabolicTargs" --num_PEER 0 --num_pcs 0 --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut" \
--select_args_label "allButCancerType"



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
