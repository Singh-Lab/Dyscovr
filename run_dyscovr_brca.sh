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

## RUN DYSCOVR SPECIFICALLY ON TCGA-BRCA (TEST CASE) ##
Rscript dyscovr.R --dataset "TCGA" --cancerType "BRCA" --run_query_genes_jointly "TRUE" \
--test "FALSE" --randomize "FALSE" --removeMetastatic "TRUE" \
--expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU.csv" \
--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/BRCA/dyscovr_input_tables/top_drivers_0.05" \
--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
--run_name "top_drivers_0.05" --target_df "allgene_targets.csv" --targets_name "allGenes" \
--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --num_pcs 3 \
--debug "FALSE" --model_type "linear" --collinearity_corr_method "eliminate_vif" \
--collinearity_diagn "FALSE" --regularization "None"  \
--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
--select_args_label "WashUPCs.elim.vif5.sp0.7" --removeCis "TRUE" --dataset "TCGA" --inclResiduals "FALSE"


## RUN DYSCOVR SPECIFICALLY ON PIK3CA (TEST CASE) ##
Rscript dyscovr.R --cancerType "BRCA" --run_query_genes_jointly "TRUE" \
--test "TRUE" --tester_name "PIK3CA" --tester_uniprot_id "P42336" --tester_ensg_id "ENSG00000121879" \
--randomize "FALSE" --removeMetastatic "TRUE" \
--expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU.csv" \
--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/BRCA/dyscovr_input_tables/top_drivers_0.05" \
--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
--target_df "glycolysis_targets.csv" --targets_name "glycolysisTargs"\
--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --num_pcs 3 \
--debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" \
--collinearity_corr_method "eliminate_vif" --regularization "None" \
--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
--select_args_label "WashUPCs.elim.vif5.sp0.7" --removeCis "TRUE" --dataset "TCGA" --inclResiduals "FALSE"

# To test on a particular subtype, for example HER2
--patientsOfInterest "HER2_patient_ids.txt" --patientsOfInterestLabel "HER2" \


## RUN DYSCOVR WITH REGULARIZATION ##
Rscript dyscovr.R --cancerType "BRCA" --run_query_genes_jointly "TRUE"\
--test "FALSE" --randomize "FALSE" --removeMetastatic "TRUE" \
--protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_drivers_nonsynonymous.csv" \
--expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU.csv" \
--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/BRCA/dyscovr_input_tables/top_drivers_0.05" \
--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
--target_df "allgene_targets.csv" --targets_name "allGenes" \
--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --num_pcs 3 \
--run_name "top_drivers_0.05" --debug "FALSE" --model_type "linear"  \
--collinearity_diagn "FALSE" --collinearity_corr_method "eliminate_vif" \
--regularization "bayesian.bglss"\
--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
--select_args_label "WashUPCs;BGLSS.elim.vif5.sp0.7" --dataset "TCGA" --inclResiduals "FALSE"

