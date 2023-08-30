#!/bin/bash
#SBATCH --qos=1day
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64GB
#SBATCH --job-name=lm_topDrivers
#SBATCH --mail-user=scamilli@princeton.edu
#SBATCH --mail-type=fail,time_limit

module add R/4.0.3


# SAMPLE RUN: Pan-Cancer and All Gene/ Metabolic Targets (NON-TNM)
#Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "ALL" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
#--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_0.05/eQTL" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" --run_name "top_0.05" \
#--target_df "metabolic_targets.csv" --targets_name "metabolicTargs" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "Nonsyn.Drivers.Vogel" --select_drivers "P04637;P42336;P01116;O75874" 
# ;Q8WZ42 TTN

# SAMPLE RUN: Pan-Cancer and All Gene Targets, broken into multiple parts for memory purposes (NON-TNM)
#Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "ALL" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
#--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_0.05/eQTL" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" --run_name "top_0.05" \
#--target_df "allgene_targets_pt5.csv" --targets_name "allGenespt5" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --collinearity_corr_method "eliminate_vif" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "Nonsyn.Drivers.Vogel.elim.vif5.sp0.7" --select_drivers "P04637;P42336;P01116;O75874" 


#--select_drivers "Q9NZR2;P04637;P42336"
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \


# Per-Cancer and All Gene/ Metabolic Targets (0.05)
#Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
#--specificTypes "BLCA;COAD;ESCA;HNSC;LGG;LUSC" \
#--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_0.05/eQTL" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" --run_name "top_0.05" \
#--target_df "allgene_targets_pt5.csv" --targets_name "allGenespt5" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --collinearity_corr_method "eliminate_vif" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "Nonsyn.Drivers.Vogel.Top5.elim.vif5.sp0.7" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" --adjust_cna_rel "FALSE" \
#--driver_file_prefix "iprotein_protein_ids_df_gr0.05Freq_vogelstein_top5_drivers_nonsynonymous_uniformHypermutRm_" \
#--path_to_driver_lists "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/" 

# ;UCEC   

#Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
#--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--specificTypes "BRCA;CESC;PRAD;UCEC" \
#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_0.05/eQTL" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" --run_name "top_0.05" \
#--target_df "allgene_targets_pt1.csv" --targets_name "allGenespt1" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --collinearity_corr_method "eliminate_vif" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "Nonsyn.Drivers.Vogel.elim.vif.5" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" --adjust_cna_rel "FALSE" \
#--driver_file_prefix "iprotein_protein_ids_df_gr0.05Freq_vogelstein_drivers_nonsynonymous_uniformHypermutRm_" \
#--path_to_driver_lists "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/" 

# Sex specific cancers: BRCA;CESC;PRAD;UCEC
#--patientsOfInterest "male_patients.txt" --patientsOfInterestLabel "male" \
#--specificTypes "ACC;BLCA;COAD;ESCA;HNSC;KICH;KIRC;KIRP;LGG;LIHC;LUAD;LUSC;MESO;PAAD;PCPG;SARC;STAD;THCA" \

# --specificTypes "ACC;BLCA;BRCA;CESC;COAD;ESCA;HNSC;KICH;KIRC;KIRP;LGG;LIHC;LUAD;LUSC;MESO;PAAD;PRAD;PCPG;SARC;STAD;THCA;UCEC" \

# Per-Cancer and All Gene/ Metabolic Targets (0.05)
Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
--patientsOfInterest "female_patients.txt" --patientsOfInterestLabel "female" \
--specificTypes "ACC;BLCA;COAD;ESCA;HNSC;KICH;KIRC;KIRP;LGG;LIHC;LUAD;LUSC;MESO;PAAD;PCPG;SARC;STAD;THCA" \
--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv" \
--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_0.05/eQTL" \
--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" --run_name "top_0.05" \
--target_df "metabolic_targets.csv" --targets_name "metabolicTargs" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
--select_args_label "Nonsyn.Drivers.Vogel.elim.vif5.sp0.7" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" --adjust_cna_rel "FALSE" \
--driver_file_prefix "iprotein_protein_ids_df_gr0.05Freq_vogelstein_drivers_nonsynonymous_uniformHypermutRm_" \
--path_to_driver_lists "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/" 

#Gender;

# Per-Cancer and All Gene/ Metabolic Targets (0.05)
#Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
#--specificTypes "CESC;COAD;ESCA;GBM;HNSC;STAD;BLCA;LUAD;LUSC;UCEC;READ" \
#--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_0.05/eQTL" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" --run_name "top_0.05" \
#--target_df "allgene_targets.csv" --targets_name "allGenes" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "Nonsyn.Drivers.Top5" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" --adjust_cna_rel "FALSE" \
#--driver_file_prefix "iprotein_protein_ids_df_gr0.05Freq_top5_drivers_nonsynonymous_uniformHypermutRm_" \
#--path_to_driver_lists "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/" 

# Per-Cancer and All Gene/ Metabolic Targets (0.05)
#Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
#--specificTypes "ESCA;COAD;STAD;BLCA" \
#--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_0.05/eQTL" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" --run_name "top_0.05" \
#--target_df "allgene_targets.csv" --targets_name "allGenes" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "Nonsyn.Drivers.Top5" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" --adjust_cna_rel "FALSE" \
#--driver_file_prefix "iprotein_protein_ids_df_gr0.05Freq_top5_drivers_nonsynonymous_uniformHypermutRm_" \
#--path_to_driver_lists "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/" 

#--select_drivers "P04637;P42336;Q8NEZ4" \

# Per-Cancer and All Gene/ Metabolic Targets (0.1)
#Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
#--specificTypes "BRCA;CESC;COAD;ESCA;GBM;HNSC;KIRC;LGG;LIHC;STAD" \
#--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_0.1/eQTL" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" --run_name "top_0.05" \
#--target_df "allgene_targets.csv" --targets_name "allGenes" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "Nonsyn.Drivers.Top5" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" --adjust_cna_rel "FALSE" \
#--driver_file_prefix "iprotein_protein_ids_df_gr0.05Freq_top5_drivers_nonsynonymous_uniformHypermutRm_" \
#--path_to_driver_lists "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/" 
#--select_drivers "P04637;P42336;Q8NEZ4" \

# Per-Cancer and All Gene/ Metabolic Targets (0.15)
#Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
#--specificTypes "BLCA;COAD;GBM;HNSC;LUAD;LUSC;STAD;UCEC" \
#--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_0.15/eQTL" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" --run_name "top_0.05" \
#--target_df "metabolic_targets.csv" --targets_name "metabolicTargs" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "Nonsyn.Drivers.Top5" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" --adjust_cna_rel "FALSE" \
#--driver_file_prefix "iprotein_protein_ids_df_gr0.05Freq_top5_drivers_nonsynonymous_uniformHypermutRm_" \
#--path_to_driver_lists "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/" 
#--select_drivers "P04637;P42336;Q8NEZ4" \


#--mutation_regprot_df "iprotein_results_nonsynonymous_IntersectPatientsWashU.csv" \

# SAMPLE RUN: Per-Cancer and All Gene/ Metabolic Targets (0.1)
#Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
#--specificTypes "BLCA;BRCA;LUSC;LUAD;STAD;COAD;KIRC;GBM;CESC;LGG;PAAD;KIRP;LIHC;HNSC;ESCA;UCEC;DLBC" \
#--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_noHypermut_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_drivers_0.1/eQTL" \
#--patient_df "combined_patient_sample_inclSubtypes_noHypermut_IntersectPatientsWashU.csv" --run_name "top_drivers_0.1" \
#--target_df "allgene_targets.csv" --targets_name "allGenes" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut" \
#--select_args_label "inclCT.sklearn.Nonsyn.0.1.WashU" --path_to_driver_lists "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/" \
#--driver_file_prefix "iprotein_protein_ids_df_gr0.1Freq_drivers_nonsynonymous_"

#--mutation_regprot_df "iprotein_results_nonsynonymous_IntersectPatientsWashU.csv" \



# SAMPLE RUN: BRCA/BLCA/HNSC and All Gene/ Metabolic Targets (NON-TNM)
#Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "ALL" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
#--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_BRCA_BLCA_HNSC_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_drivers_0.05_BRCABLCAHNSC/eQTL" \
#--patient_df "combined_patient_sample_cibersort_total_frac_BRCA_BLCA_HNSC_inclBRCA.BLCA.HNSCsubtypes_IntersectPatientsWashU.csv" --run_name "top_drivers_0.05_BRCABLCAHNSC" \
#--target_df "metabolic_targets.csv" --targets_name "metabolicTargs" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "inclCT.sklearn.Nonsyn.0.05.WashU" 


# SAMPLE RUN: COAD/READ and All Gene/ Metabolic Targets (NON-TNM)
#Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "ALL" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
#--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_COAD_READ_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_drivers_0.1_COADREAD/eQTL/subtypes/msi" \
#--patient_df "combined_patient_sample_cibersort_total_frac_COAD_READ_IntersectPatientsWashU.csv" --run_name "top_drivers_0.1_COADREAD" \
#--target_df "metabolic_targets.csv" --targets_name "metabolicTargs" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "inclCT.sklearn.Nonsyn.0.1.WashU.MSI" 

#Rscript linear_model3.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "ALL" --clinical_df "clinical.csv" --run_query_genes_jointly "TRUE" \
#--test "FALSE" --randomize "FALSE" --expression_df "expression_quantile_norm_COAD_READ_IntersectPatientsWashU.csv" \
#--input_lm_filepath "/Genomics/grid/users/scamilli/thesis_work/run-model-R/input_files/PanCancer/tumor_only/lm_input_tables/top_all_0.15_COADREAD/eQTL/" \
#--patient_df "combined_patient_sample_cibersort_total_frac_COAD_READ_IntersectPatientsWashU.csv" --run_name "top_all_0.15_COADREAD" \
#--target_df "metabolic_targets.csv" --targets_name "metabolicTargs" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M"  \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --model_type "linear" --collinearity_diagn "FALSE" --regularization "None"  \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Gender;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "inclCT.sklearn.Nonsyn.0.15.WashU.MSI" 


#--mutation_regprot_df "iprotein_results_nonsynonymous_COAD_READ_IntersectPatientsWashU.csv" \
# cms




### ARCHIVAL RUNS

# SAMPLE RUN: TP53 and all targets/ ChIP-eat targets (NON-TNM)
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "HNSC" --clinical_df "clinical.csv" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "TRUE" --expression_df "expression_tmm_filtByExpr_SDGr1_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatients.csv" --mutation_targ_df "mut_count_matrix_missense_nonsense_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_missense_nonsense_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_MN_IntersectPatients.csv" --target_df "tca_cycle_metabolic_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --targets_name "tcaCycleTargs" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 0 --mut_pc_vect "0" --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;PIK3CA_MutStat_i;PIK3CA_AmplStat_i;PIK3CAP53_MutStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut" \
#--select_args_label "PIK3CA_Covs" --removeCis "TRUE"

# Version using BRCA/BLCA/HNSC-specific files
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "ALL" --clinical_df "clinical.csv" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "TRUE" --incl_nextMutDriver "TRUE" --expression_df "expression_tmm_filtByExpr_SDGr1_BRCA_BLCA_HNSC_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_BRCA_BLCA_HNSC_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_BRCA_BLCA_HNSC_CancerOnly_IntersectPatients.csv" --mutation_targ_df "mut_count_matrix_missense_nonsense_BRCA_BLCA_HNSC_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_missense_nonsense_BRCA_BLCA_HNSC_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_BRCA_BLCA_HNSC_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_MN_BRCA_BLCA_HNSC_IntersectPatients.csv" --target_df "metabolic_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --targets_name "metabolicTargs" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 0 --mut_pc_vect "0" --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;PIK3CA_MutStat_i;PIK3CA_AmplStat_i;PIK3CAP53_MutStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "inclCT.PIK3CA_Covs.BRCA.BLCA.HNSC" --removeCis "TRUE"


# SAMPLE RUN: PIK3CA (Version using BRCA/BLCA/HNSC-specific files)
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "ALL" --clinical_df "clinical.csv" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "TRUE" --tester_name "PIK3CA" --tester_uniprot_id "P42336" --tester_ensg_id "ENSG00000121879" --incl_nextMutDriver "TRUE" \
#--expression_df "expression_tmm_filtByExpr_SDGr1_BRCA_BLCA_HNSC_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_BRCA_BLCA_HNSC_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_BRCA_BLCA_HNSC_CancerOnly_IntersectPatients.csv" --mutation_targ_df "mut_count_matrix_missense_nonsense_BRCA_BLCA_HNSC_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_missense_nonsense_BRCA_BLCA_HNSC_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_BRCA_BLCA_HNSC_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_MN_BRCA_BLCA_HNSC_IntersectPatients.csv" --target_df "allgene_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --targets_name "allGenes" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 0 --mut_pc_vect "0" --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;TP53_MutStat_i;TP53_AmplStat_i;TP53PIK3CA_MutStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "inclCT.TP53_Covs.BRCA.BLCA.HNSC" --removeCis "TRUE"

# SAMPLE RUN: TTN, Negative Control (Version using BRCA/BLCA/HNSC-specific files)
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "ALL" --clinical_df "clinical.csv" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "TRUE" --tester_name "TTN" --tester_uniprot_id "Q8WZ42" --tester_ensg_id "ENSG00000155657" --incl_nextMutDriver "TRUE" \
#--expression_df "expression_tmm_filtByExpr_SDGr1_BRCA_BLCA_HNSC_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_BRCA_BLCA_HNSC_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_BRCA_BLCA_HNSC_CancerOnly_IntersectPatients.csv" --mutation_targ_df "mut_count_matrix_missense_nonsense_BRCA_BLCA_HNSC_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_missense_nonsense_BRCA_BLCA_HNSC_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_BRCA_BLCA_HNSC_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_MN_BRCA_BLCA_HNSC_IntersectPatients.csv" --target_df "allgene_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --targets_name "allGenes" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 0 --mut_pc_vect "0" --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;TP53_MutStat_i;TP53_AmplStat_i;TP53PIK3CA_MutStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "TP53_Covs.BRCA.BLCA.HNSC" --removeCis "TRUE"


# SAMPLE RUN: Tester Protein and all targets/ ChIP-eat targets (NON-TNM) -- PER-CANCER
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "UCEC" --clinical_df "clinical.csv" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "TRUE" --tester_name "MUC16" --tester_uniprot_id "Q8WXI7" --tester_ensg_id "ENSG00000181143"  \
#--expression_df "expression_tmm_filtByExpr_SDGr1_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatients.csv" --mutation_targ_df "mut_count_matrix_missense_nonsense_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_missense_nonsense_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_MN_IntersectPatients.csv" --target_df "metabolic_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --targets_name "metabolicTargs" --num_PEER 0 --useNumFunctCopies "FALSE" \
#--num_pcs 0 --mut_pc_vect "0" --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut" \
#--select_args_label "" --removeCis "TRUE"

# Other cancer types of interest: "PAAD" (pancreatic); "COAD" (colon)
# --test "TRUE" --tester_name "IDH1" --tester_uniprot_id "O75874" --tester_ensg_id "ENSG00000138413"  \
# --test "TRUE" --tester_name "KRAS" --tester_uniprot_id "190070" --tester_ensg_id "ENSG00000133703"  \


# SAMPLE RUN: USH2A and all targets (NON-TNM) -- NEGATIVE CONTROL
#Rscript linear_model.R --QTLtype "meQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --randomize "FALSE" \
#--test "TRUE" --tester_name "USH2A" --tester_uniprot_id "O75445" --tester_ensg_id "ENSG00000042781" \
#--expression_df "expression_tmm_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_bucketed_Beta_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_Beta_CancerOnly_IntersectPatients.csv" --mutation_targ_df "mut_count_matrix_missense_CancerOnly_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_missense_CancerOnly_IntersectPatients.csv" --cna_df "CNA_AllGenes_Bucketed_InclAmp_CancerOnly_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_tmm_IntersectPatients.csv" --target_df "allgene_targets.csv" \
#--cna_bucketing "bucket_inclAmp" --meth_bucketing "TRUE" --meth_type "Beta" --targets_name "allGenes" --num_PEER 0 --num_pcs 0 \
#--debug "FALSE" --collinearity_diagn "TRUE" --regularization "None" --select_args "MutStat_i;ExpStat_k"


# SAMPLE RUN: Cancer Related Genes/ Highly Deleted TFs and Metabolic Targets (NON-TNM)
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "BLCA" \
#--patientsOfInterest "BLCA.NotTP53andPIK3CAMNmut_patient_ids.txt" --patientsOfInterestLabel "BLCA.NoMutOverlap.MN" \
#--clinical_df "clinical.csv" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "FALSE" --protein_ids_df "iprotein_protein_ids_df_sixMostMut.csv" \
#--expression_df "expression_tmm_filtByExpr_SDGr1_CancerOnly_IntersectPatients.csv" --methylation_df "methylation_M_CancerOnly_IntersectPatients.csv" \
#--methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatients.csv" --mutation_targ_df "mut_count_matrix_missense_IntersectPatients.csv" \
#--mutation_regprot_df "iprotein_results_missense_nonsense_IntersectPatients.csv" --cna_df "CNA_AllGenes_CancerOnly_IntersectPatients.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_MN_IntersectPatients.csv" --target_df "allgene_targets.csv" \
#--cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "sixMostMut" --targets_name "allGenes" --num_PEER 0 \
#--num_pcs 0 --mut_pc_vect "0" --debug "FALSE" --collinearity_diagn "FALSE" --regularization "None" --useNumFunctCopies "FALSE" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;Prior_malig;Tot_Mut" \
#--select_args_label "allButCancerType" --removeCis "TRUE"


# SAMPLE RUN: Cancer Related Genes/ Highly Deleted TFs and Metabolic Targets (NON-TNM); with regularization
#Rscript linear_model.R --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "ALL" --randomize "FALSE" --removeMetastatic "TRUE" \
#--test "FALSE" --protein_ids_df "iprotein_protein_ids_df_gr0.025Freq_drivers_nonsynonymous_BRCA_BLCA_HNSC.csv" --incl_nextMutDriver "FALSE" --clinical_df "clinical.csv" \
#--expression_df "expression_quantile_norm_BRCA_BLCA_HNSC_IntersectPatientsWashU.csv" \
#--methylation_df "methylation_M_CancerOnly_BRCA_BLCA_HNSC_IntersectPatientsWashU.csv" --methylation_df_meQTL "methylation_M_CancerOnly_BRCA_BLCA_HNSC_IntersectPatientsWashU.csv" \
#--mutation_targ_df "mut_count_matrix_nonsynonymous_BRCA_BLCA_HNSC_IntersectPatientsWashU.csv" --mutation_regprot_df "iprotein_results_nonsynonymous_BRCA_BLCA_HNSC_IntersectPatientsWashU.csv" \
#--cna_df "CNA_AllGenes_CancerOnly_BRCA_BLCA_HNSC_IntersectPatientsWashU.csv" --patient_df "combined_patient_sample_cibersort_total_frac_BRCA_BLCA_HNSC_inclBRCA.BLCA.HNSCsubtypes_IntersectPatientsWashU.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "top_drivers" \
#--targets_name "allGenes" --num_PEER 0 --num_pcs 3 --debug "FALSE" --collinearity_diagn "FALSE" --regularization "L1" --signif_eval_type "randomization" --useNumFunctCopies "FALSE" \
#--select_args "ExpStat_k;MutStat_i;CNAStat_i;MethStat_i;CNAStat_k;MutStat_k;MethStat_k;Age;Tumor_purity;Treatment_rad;Treatment_pharm;Tot_IC_Frac;PC;Prior_malig;Tot_Mut;Cancer_type" \
#--select_args_label "inclCT.sklearn.Nonsyn.0.025.WashU" --removeCis "TRUE" --dataset "TCGA" --run_query_genes_jointly "TRUE" --model_type "linear" 




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
# FPKM: "expression_fpkm_CancerOnly_IntersectPatientsFPKM.csv"  # fpkm-normalized and filtered
# Rank-Normalized: "expression_rank_norm_CancerOnly_IntersectPatientsRN.csv"

###########################################################################################
# Methylation DF
###########################################################################################
# Beta Values, Cancer Only, FPKM-Norm Exp: "methylation_Beta_CancerOnly_IntersectPatientsFPKM.csv"
# Beta Values, Cancer Only, Rank-Norm Exp: "methylation_Beta_CancerOnly_IntersectPatientsRN.csv"

# M-Values, Cancer Only, FPKM-Norm Exp: "methylation_M_CancerOnly_IntersectPatientsFPKM.csv"
# M-Values, Cancer Only, Rank-Norm Exp: "methylation_M_CancerOnly_IntersectPatientsRN.csv"

# Beta (Bucketed, 3 Buckets), Cancer Only, FPKM-Norm Exp: "methylation_bucketed_Beta_CancerOnly_IntersectPatientsFPKM.csv"
# Beta (Bucketed, 3 Buckets), Cancer Only, Rank-Norm Exp: "methylation_bucketed_Beta_CancerOnly_IntersectPatientsRN.csv"

# M (Bucketed, 3 Buckets), Cancer Only, FPKM-Norm Exp: "methylation_bucketed_M_CancerOnly_IntersectPatientsFPKM.csv"
# M (Bucketed, 3 Buckets), Cancer Only, Rank-Norm Exp: "methylation_bucketed_M_CancerOnly_IntersectPatientsRN.csv"

# Threshold (Beta > 0.8), Cancer Only, FPKM-Norm Exp: "methylation_Beta_0.8_CancerOnly_IntersectPatientsFPKM.csv"
# Threshold (Beta > 0.8), Cancer Only, Rank-Norm Exp: "methylation_Beta_0.8_CancerOnly_IntersectPatientsRN.csv"

###########################################################################################
# Mutation Target DF
###########################################################################################
# Missense, FPKM-Norm Exp: "mut_count_matrix_missense_IntersectPatientsFPKM.csv"
# Missense, Rank-Norm Exp: "mut_count_matrix_missense_IntersectPatientsRN.csv"

###########################################################################################
# Regulatory Protein Mutation DF
###########################################################################################
# I-Protein, FPKM-Norm Exp: "iprotein_results_missense_IntersectPatientsFPKM.csv"
# I-Protein, Rank-Norm Exp: "iprotein_results_missense_IntersectPatientsRN.csv"

# TO BE CREATED:
# I-Protein (Nucleic Acids): "iprotein_results_missense_nucacids_IntersectPatients.csv"
# I-Domain: "idomain_results_missense_IntersectPatients.csv"
# I-Domain (Nucleic Acids): "idomain_results_missense_nucacids_IntersectPatients.csv"
# I-Binding Position: "ibindingpos_results_missense_IntersectPatients.csv"
# I-Binding Position (Nucleic Acids): "ibindingpos_results_missense_nucacids_IntersectPatients.csv"

###########################################################################################
# CNA DF
###########################################################################################
# All Genes (Cancer Samples Only), FPKM-Norm Exp: "CNA_AllGenes_CancerOnly_IntersectPatientsFPKM.csv"
# All Genes (Cancer Samples Only), Rank-Norm Exp: "CNA_AllGenes_CancerOnly_IntersectPatientsRN.csv"

# All Genes (Cancer Samples Only), Bucketed, Including Amplifications, FPKM-Norm Exp: "CNA_AllGenes_Bucketed_InclAmp_CancerOnly_IntersectPatientsFPKM.csv"
# All Genes (Cancer Samples Only), Bucketed, Including Amplifications, Rank-Norm Exp: "CNA_AllGenes_Bucketed_InclAmp_CancerOnly_IntersectPatientsRN.csv"

# All Genes (Cancer Samples Only), Bucketed, Excluding Amplifications, FPKM-Norm Exp: "CNA_AllGenes_Bucketed_ExclAmp_CancerOnly_IntersectPatientsFPKM.csv"
# All Genes (Cancer Samples Only), Bucketed, Excluding Amplifications, Rank-Norm Exp: "CNA_AllGenes_Bucketed_ExclAmp_CancerOnly_IntersectPatientsRN.csv"

# TO BE CREATED:
# TFs Only: "CNA_TFsOnly_IntersectPatients.csv"

###########################################################################################
# Patient DF
###########################################################################################
# CIBERSORT Abs I.C.D.: 
# 1. "combined_patient_sample_cibersort_abs_fpkm_IntersectPatientsFPKM.csv"
# 2. "combined_patient_sample_cibersort_abs_rank_norm_IntersectPatientsRN.csv"

# TIMER I.C.D.: 
# 1. "combined_patient_sample_timer_fpkm_IntersectPatientsFPKM.csv"
# 2. "combined_patient_sample_timer_rn_IntersectPatientsRN.csv"

# CIBERSORT Total Frac. I.C.D.: 
# 1. "combined_patient_sample_cibersort_total_frac_fpkm_IntersectPatientsFPKM.csv"
# 2. "combined_patient_sample_cibersort_total_frac_rn_IntersectPatientsRN.csv"

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
