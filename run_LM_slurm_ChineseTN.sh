#!/bin/bash
#SBATCH --qos=1day
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=24GB
#SBATCH --job-name=lm_topDrivers_chineseTN
#SBATCH --mail-user=scamilli@princeton.edu
#SBATCH --mail-type=fail,time_limit

module add R/4.0.3

Rscript linear_model3.R --dataset "Chinese_TN" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --run_query_genes_jointly "TRUE" \
--test "FALSE" --randomize "FALSE" --expression_df "expression_df_FPKM_filtBySD1_uniformHypermutRm_IntersectPatients.csv" \
--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/Chinese_TN/lm_input_tables/top_0.05/eQTL" \
--patient_df "patient_sample_df_uniformHypermutRm_IntersectPatients.csv" --run_name "top_0.05" \
--target_df "allgene_targets.csv" --targets_name "allGenes" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
--num_PEER 0 --num_pcs 0 --debug "TRUE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
--select_args "ExpStat_k;MutStat_i;CNAStat_i;CNAStat_k;MutStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac" \
--select_args_label "Nonsyn.Driver.Vogel" --select_drivers "P04637;P42336" 
