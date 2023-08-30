#!/bin/bash
#SBATCH --qos=1day
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=24GB
#SBATCH --job-name=lm_topDrivers_metabric
#SBATCH --mail-user=scamilli@princeton.edu
#SBATCH --mail-type=fail,time_limit

module add R/4.0.3


Rscript linear_model3.R --dataset "METABRIC" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --run_query_genes_jointly "TRUE" \
--test "FALSE" --randomize "FALSE" --expression_df "expression_df_log2intensity_filtBySD0.2_IntersectPatients.csv" \
--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/METABRIC/lm_input_tables/top_0.05/eQTL" \
--patient_df "patient_sample_df_cibersort_total_frac_IntersectPatients.csv" --run_name "top_0.05" \
--target_df "allgene_targets.csv" --targets_name "allGenes" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
--num_PEER 0 --num_pcs 0 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
--select_args_label "Nonsyn.Driver.Vogel.elim.vif.5" --select_drivers "P04637;P42336;Q13233;P23771;Q5VU43;P46531" 

#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/BRCA/tumor_only/lm_input_tables/top_drivers_0.025/eQTL" \
#--mutation_regprot_df "iprotein_mutation_df_nonsynonymous_IntersectPatients.csv" \


## OLD ##

# SAMPLE RUN: TP53 and all targets/ ChIP-eat targets (NON-TNM)
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" --removeMetastatic "FALSE" \
#--test "TRUE" --incl_nextMutDriver "TRUE" --expression_df "expression_log2intensity_filtBySD0.2_IntersectPatients.csv" --methylation_df "methylation_df_RRBS_logit_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_df_RRBS_logit_IntersectPatients.csv" --mutation_targ_df "mutation_count_df_missense_nonsense_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_mutation_df_missense_nonsense_IntersectPatients.csv" --cna_df "cna_df_allgenes_IntersectPatients.csv" \
#--patient_df "patient_sample_df_cibersort_total_frac_wTMB_IntersectPatients.csv" --target_df "allgene_targets.csv" \
#--cna_bucketing "bucket_justDel" --meth_bucketing "FALSE" --meth_type "logit" --targets_name "allGenes" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 0 --mut_pc_vect "0" --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;PIK3CA_MutStat_i;PIK3CA_CNAStat_i;PIK3CAP53_MutStat_i;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "PIK3CA_covs.inclCT.MN" --removeCis "TRUE" --dataset "METABRIC"

#--patientsOfInterest "LT20PAmp.NotTP53andPIK3CAMNmut_patient_ids.txt" --patientsOfInterestLabel "LT20PAmp.noMutOverlap" \
#--patientsOfInterest "Luminal.B_patient_ids.txt" --patientsOfInterestLabel "LumB" \
#--patientsOfInterest "LessThan20PercAmp_patient_ids.txt" --patientsOfInterestLabel "LT20PAmp" \

# SAMPLE RUN: PIK3CA and all targets/ ChIP-eat targets (NON-TNM)
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" --removeMetastatic "FALSE" \
#--test "TRUE" --tester_name "PIK3CA" --tester_uniprot_id "P42336" --tester_ensg_id "ENSG00000121879" --incl_nextMutDriver "TRUE" \
#--expression_df "expression_log2intensity_filtBySD0.2_IntersectPatients.csv" --methylation_df "methylation_df_RRBS_logit_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_df_RRBS_logit_IntersectPatients.csv" --mutation_targ_df "mutation_count_df_missense_nonsense_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_mutation_df_missense_nonsense_IntersectPatients.csv" --cna_df "cna_df_allgenes_IntersectPatients.csv" \
#--patient_df "patient_sample_df_cibersort_total_frac_wTMB_IntersectPatients.csv" --target_df "allgene_targets.csv" \
#--cna_bucketing "bucket_justAmp" --meth_bucketing "FALSE" --meth_type "logit" --targets_name "allGenes" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 0 --mut_pc_vect "0" --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;TP53_MutStat_i;TP53_AmplStat_i;TP53PIK3CA_MutStat_i;CNAStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "TP53_Covs.inclCT.MN" --removeCis "TRUE" --dataset "METABRIC"

#--functional_copies_df "pik3ca_functional_copies_per_sample.csv"

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
#--select_args_label "inclCT" --removeCis "TRUE"

# SAMPLE RUN: Other Genes and all targets (NON-TNM) 
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" --removeMetastatic "FALSE" --incl_nextMutDriver "TRUE" \
#--test "TRUE" --tester_name "FOXA1" --tester_uniprot_id "P55317" --tester_ensg_id "ENSG00000129514" \
#--expression_df "expression_df_log2intensity_filtBySD0.2_IntersectPatients.csv" --methylation_df "methylation_df_RRBS_logit_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_df_RRBS_logit_IntersectPatients.csv" --mutation_targ_df "mutation_count_df_missense_nonsense_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_mutation_df_missense_nonsense_IntersectPatients.csv" --cna_df "cna_df_allgenes_IntersectPatients.csv" \
#--patient_df "patient_sample_df_cibersort_total_frac_wTMB_IntersectPatients.csv" --target_df "allgene_targets.csv" \
#--cna_bucketing "bucket_justDel" --meth_bucketing "FALSE" --meth_type "logit" --targets_name "allGenes" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 0 --mut_pc_vect "0" --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;TP53_MutStat_i;TP53_AmplStat_i;TP53PIK3CA_MutStat_i;CNAStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "TP53_Covs.inclCT.MN" --removeCis "TRUE" --dataset "METABRIC"

#--patientsOfInterest "LessThan20PercAmp_patient_ids.txt" --patientsOfInterestLabel "LT20PAmp" \

#--test "TRUE" --tester_name "SF3B1" --tester_uniprot_id "O75533" --tester_ensg_id "ENSG00000115524" \  # O
#--test "TRUE" --tester_name "USH2A" --tester_uniprot_id "O75445" --tester_ensg_id "ENSG00000042781" \  # C
#--test "TRUE" --tester_name "KMT2C" --tester_uniprot_id "Q8NEZ4" --tester_ensg_id "ENSG00000055609" \  # TS
#--test "TRUE" --tester_name "GATA3" --tester_uniprot_id "P23771" --tester_ensg_id "ENSG00000107485" \ # B
#--test "TRUE" --tester_name "TTN" --tester_uniprot_id "Q8WZ42" --tester_ensg_id "ENSG00000155657" \    # C
#--test "TRUE" --tester_name "FOXA1" --tester_uniprot_id "P55317" --tester_ensg_id "ENSG00000129514" \  # TS

# SAMPLE RUN: Cancer Related Genes/ Highly Deleted TFs and Metabolic Targets (NON-TNM)
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --randomize "FALSE" --removeMetastatic "TRUE" \
#--patientsOfInterest "Luminal.B.NotTP53andPIK3CAMNmut_patient_ids.txt" --patientsOfInterestLabel "LumB.NoMutOverlap.MN" \
#--test "FALSE" --protein_ids_df "iprotein_protein_ids_df_sixMostMut.csv" --expression_df "expression_tmmSDGr1filtByExpr_CancerOnly_IntersectPatients.csv" \
#--methylation_df "methylation_M_CancerOnly_IntersectPatients.csv" --methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatients.csv" \
#--mutation_targ_df "mut_count_matrix_missense_nonsense_CancerOnly_IntersectPatients.csv" --mutation_regprot_df "iprotein_results_missense_nonsense_CancerOnly_IntersectPatients.csv" \
#--cna_df "CNA_AllGenes_CancerOnly_IntersectPatients.csv" --patient_df "combined_patient_sample_cibersort_total_frac_tmm_normAge_IntersectPatients.csv" \
#--target_df "matt_metabolic_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "sixMostMut" \
#--targets_name "1CTargs" --num_PEER 0 --num_pcs 0 --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" --useNumFunctCopies "FALSE" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut" \
#--select_args_label "allButCancerType" --removeCis "TRUE"

#--patientsOfInterest "Luminal.A_patient_ids.txt" --patientsOfInterestLabel "LumA" \
#--patientsOfInterest "LessThan20PercAmp_patient_ids.txt" --patientsOfInterestLabel "LT20PAmp" \
#--patientsOfInterest "NotTP53andPIK3CAmut_patient_ids.txt" --patientsOfInterestLabel "NoMutOverlap" \



