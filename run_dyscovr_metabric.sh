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

## RUN DYSCOVR ON THE METABRIC COHORT ##
Rscript linear_model3.R --dataset "METABRIC" --cancerType "BRCA" --run_query_genes_jointly "TRUE" \
--test "FALSE" --randomize "FALSE" --expression_df "expression_df_quantile_norm_sklearn_IntersectPatients.csv" \
--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/METABRIC/lm_input_tables/top_0.05" \
--patient_df "patient_sample_df_cibersort_total_frac_IntersectPatients.csv" \
--run_name "top_0.05" --target_df "allgene_targets.csv" --targets_name "allGenes" \
--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --num_pcs 0 \
--debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" \
--collinearity_corr_method "eliminate_vif" --regularization "None"  \
--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
--select_args_label "WashUPCs.elim.vif5.sp0.7" \
--select_drivers "P04637;P42336;Q13233;P23771;Q5VU43;P46531" 