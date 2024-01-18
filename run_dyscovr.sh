#!/bin/bash
#SBATCH --qos=1day
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80GB
#SBATCH --job-name=lm_topDrivers
#SBATCH --mail-user=scamilli@princeton.edu
#SBATCH --mail-type=fail,time_limit

module add R/4.0.5

# If needed:
# install.packages(c("speedglm", "Hmisc", "caret", "gglasso", "statmod", "mvtnorm", "MCMCpack", "SuppDists", "mvnfast", "MBSGS", "MASS", "mgcv", "mnormt", "truncnorm", "glmnet", "stringr", "rlang", "broom", "qvalue", "gplots", "ggplot2", "olsrr", "RColorBrewer"), dependencies = TRUE, INSTALL_opts = '--no-lock')
# BiocManager::install('qvalue')

## RUN DYSCOVR PAN-CANCER ##
#Rscript dyscovr.R --dataset "TCGA" --cancerType "PanCancer" --specificTypes "ALL" \
#--clinical_df "clinical.csv" --run_query_genes_jointly "T" --test "F" \
#--randomize "F" --expression_df "expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/argo/users/scamilli/Dyscovr/input_files/PanCancer/lm_input_tables/top_0.05" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--run_name "top_0.05" --target_df "metabolic_targets.csv" --targets_name "metabolicTargs" \
#--cna_bucketing "rawCNA" --meth_bucketing "F" --meth_type "M" --num_pcs 3 \
#--debug "F" --model_type "linear" --collinearity_diagn "F" \
#--collinearity_corr_method "eliminate_vif" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "elim.vif5.sp0.7" --select_drivers "P04637;P42336;P01116;O75874" 

# To add TTN, add 'Q8WZ42' to the list of select_drivers

## RUN DYSCOVR PAN-CANCER ON ALL GENE TARGETS ##
## For memory purposes, allgene_targets.csv has been broken into 5 parts; run
## each separately and bind results together afterwards. Will need to repeat 
## q-value correction and sorting.
Rscript dyscovr.R --dataset "TCGA" --cancerType "PanCancer" --specificTypes "ALL" \
--clinical_df "clinical.csv" --run_query_genes_jointly "T" --test "F" \
--randomize "F" --expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU_orig.csv" \
--input_lm_filepath "/Genomics/argo/users/scamilli/Dyscovr/input_files/PanCancer/lm_input_tables/top_0.05_orig" \
--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
--run_name "top_0.05_orig" --target_df "allgene_targets_pt5.csv" --targets_name "allGenespt5" \
--cna_bucketing "rawCNA" --meth_bucketing "F" --meth_type "M" --num_pcs 3 \
--debug "F" --model_type "linear" --collinearity_diagn "F" \
--collinearity_corr_method "eliminate_vif" --regularization "None"  \
--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
--select_args_label "elim.vif5.sp0.7" --select_drivers "P04637;P42336;P01116;O75874" 

## RUN DYSCOVR PER-CANCER ON ALL GENE TARGETS (EACH CANCER TYPE INDIVIDUALLY) ##
#Rscript dyscovr.R --dataset "TCGA" --cancerType "PanCancer" --clinical_df "clinical.csv" \
#--run_query_genes_jointly "T" \
#--specificTypes "ACC;BLCA;BRCA;CESC;COAD;ESCA;HNSC;KICH;KIRC;KIRP;LGG;LIHC;LUAD;LUSC;MESO;PAAD;PCPG;PRAD;SARC;STAD;THCA;UCEC" \
#--test "F" --randomize "F" --expression_df "expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/argo/users/scamilli/Dyscovr/input_files/PanCancer/lm_input_tables/top_0.05" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--run_name "top_0.05" --target_df "allgene_targets.csv" --targets_name "allGenes" \
#--cna_bucketing "rawCNA" --meth_bucketing "F" --meth_type "M" --num_pcs 3 \
#--debug "F" --model_type "linear" --collinearity_diagn "F" \
#--collinearity_corr_method "eliminate_vif" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "elim.vif5.sp0.7" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--adjust_cna_rel "F" --driver_file_prefix "iprotein_protein_ids_df_gr0.05Freq_drivers_nonsynonymous_uniformHypermutRm_" \
#--path_to_driver_lists "/Genomics/argo/users/scamilli/Dyscovr/input_files/PanCancer" 

## For cancer types with >5 drivers mutated at >5% frequency, optionally limit to the top 5
#Rscript dyscovr.R --dataset "TCGA" --cancerType "PanCancer" --clinical_df "clinical.csv" \
#--run_query_genes_jointly "T" --specificTypes "BLCA;COAD;ESCA;HNSC;LGG;LUSC;UCEC" \
#--test "F" --randomize "F" \
#--expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/argo/users/scamilli/Dyscovr/input_files/PanCancer/lm_input_tables/top_0.05" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--run_name "top_0.05" --target_df "allgene_targets_pt1.csv" --targets_name "allGenespt1" \
#--cna_bucketing "rawCNA" --meth_bucketing "F" --meth_type "M" --num_pcs 3 \
#--debug "F" --model_type "linear" --collinearity_diagn "F" \
#--collinearity_corr_method "eliminate_vif" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "Top5.elim.vif5.sp0.7" \
#--cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" 
#--adjust_cna_rel "F" --driver_file_prefix "iprotein_protein_ids_df_gr0.05Freq_vogelstein_top5_drivers_nonsynonymous_uniformHypermutRm_" \
#--path_to_driver_lists "/Genomics/argo/users/scamilli/Dyscovr/input_files/PanCancer" 

# Run only on sex-specific cancers
#Rscript dyscovr.R --dataset "TCGA" --cancerType "PanCancer" --clinical_df "clinical.csv" \
#--run_query_genes_jointly "T" --test "F" --randomize "F" \
#--expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--specificTypes "BRCA;CESC;PRAD;UCEC" \
#--input_lm_filepath "/Genomics/argo/users/scamilli/Dyscovr/input_files/PanCancer/lm_input_tables/top_0.05" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--run_name "top_0.05" --target_df "allgene_targets_pt1.csv" --targets_name "allGenespt1" \
#--cna_bucketing "rawCNA" --meth_bucketing "F" --meth_type "M" --num_pcs 3 \
#--debug "F" --model_type "linear" --collinearity_diagn "F" \
#--collinearity_corr_method "eliminate_vif" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "elim.vif.5" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--adjust_cna_rel "F" --driver_file_prefix "iprotein_protein_ids_df_gr0.05Freq_vogelstein_drivers_nonsynonymous_uniformHypermutRm_" \
#--path_to_driver_lists "/Genomics/argo/users/scamilli/Dyscovr/input_files/PanCancer" 

# Sex specific cancers: BRCA;CESC;PRAD;UCEC

# Run only on NON-sex specific cancers
#Rscript dyscovr.R --dataset "TCGA" --cancerType "PanCancer" --clinical_df "clinical.csv" \
#--run_query_genes_jointly "T" --patientsOfInterest "female_patients.txt" \
#--patientsOfInterestLabel "female" \
#--specificTypes "ACC;BLCA;COAD;ESCA;HNSC;KICH;KIRC;KIRP;LGG;LIHC;LUAD;LUSC;MESO;PAAD;PCPG;SARC;STAD;THCA" \
#--test "F" --randomize "F" --expression_df "expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/argo/users/scamilli/Dyscovr/input_files/PanCancer/lm_input_tables/top_0.05" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv"\
#--run_name "top_0.05" --target_df "metabolic_targets.csv" --targets_name "metabolicTargs" \
#--cna_bucketing "rawCNA" --meth_bucketing "F" --meth_type "M"  --num_pcs 3 \
#--debug "F" --model_type "linear" --collinearity_diagn "F" \
#--collinearity_corr_method "eliminate_vif" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "Nonsyn.Drivers.Vogel.elim.vif5.sp0.7" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" --adjust_cna_rel "F" \
#--driver_file_prefix "iprotein_protein_ids_df_gr0.05Freq_vogelstein_drivers_nonsynonymous_uniformHypermutRm_" \
#--path_to_driver_lists "/Genomics/argo/users/scamilli/Dyscovr/input_files/PanCancer" 

# To limit to either male or female patients, include these lines:
#--patientsOfInterest "male_patients.txt" --patientsOfInterestLabel "male" \
#--patientsOfInterest "female_patients.txt" --patientsOfInterestLabel "female" \
# And ensure you remove the "Gender" argument
