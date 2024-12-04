#!/bin/bash
#SBATCH --qos=1wk
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --job-name=create_dyscovr_input_tables
#SBATCH --mail-user=scamilli@princeton.edu
#SBATCH --mail-type=fail,time_limit

module add R/4.0.5

# If needed:
# install.packages(c("dplyr", "snow", "argparse", "data.table", "doParallel", "dqrng", "foreach", "parallel", "fastmatch"), dependencies = TRUE, INSTALL_opts = '--no-lock')

# Create linear model (lm) input tables for all gene targets, breast cancer (BRCA)
#Rscript create_dyscovr_input_tables.R --dataset "TCGA" --cancerType "BRCA" --removeMetastatic "T" \
#--randomize "F" --test "F" --run_query_genes_jointly "T" \
#--protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_vogelstein_drivers_nonsynonymous.csv" \
#--expression_df "expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--methylation_df "methylation_M_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--mutation_df "mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "F" \
#--meth_type "M" --run_name "top_0.05" --num_pcs 3 --debug "F" --removeCis "F" \

# Add line to add a particular subtype of interest, e.g. Luminal B
#--patientsOfInterest "Luminal.B_patient_ids.txt" --patientsOfInterestLabel "LumB" \
# and adjust the protein IDs DF to include only drivers highly mutated in this subtype
#--protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_all_nonsynonymous_lumB.csv" \


# Create LM input tables for all gene targets, METABRIC
#Rscript create_dyscovr_input_tables.R --dataset "METABRIC" --cancerType "BRCA" \
#--removeMetastatic "F" --randomize "F" --test "F" \
#--run_query_genes_jointly "T" --protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_vogelstein_drivers_nonsynonymous.csv" \
#--expression_df "expression_df_quantile_norm_sklearn_IntersectPatients.csv" \
#--methylation_df "methylation_df_RRBS_logit_IntersectPatients.csv" \
#--mutation_df "mutation_count_df_nonsynonymous_IntersectPatients.csv" \
#--cna_df "cna_df_allgenes_IntersectPatients.csv" \
#--patient_df "patient_sample_df_cibersort_total_frac_IntersectPatients.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "F" \
#--meth_type "M" --run_name "top_0.05" --num_pcs 0 --debug "F" --removeCis "F" 


# Create LM input tables for all gene targets, pan-cancer (all cancer types together)
Rscript create_dyscovr_input_tables.R --dataset "TCGA" --cancerType "PanCancer" --specificTypes "ALL" \
--removeMetastatic "T" --randomize "F" --test "F" --run_query_genes_jointly "T" \
--protein_ids_df "protein_ids_df_gr0.05Freq_vogelstein_all_nonsynonymous_uniformHypermutRm.csv" \
--expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU_orig.csv" \
--methylation_df "methylation_M_CancerOnly_uniformHypermutRm_IntersectPatientsWashU_orig.csv" \
--mutation_df "mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU_orig.csv" \
--cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "F" \
--meth_type "M" --run_name "top_0.05_orig_wttn" --num_pcs 3 --debug "F" --removeCis "F" 


# Create LM input tables for all gene targets, per-cancer (all cancer types individually)
#Rscript create_dyscovr_input_tables.R --dataset "TCGA" --cancerType "PanCancer" \
#--specificTypes "ACC;BLCA;BRCA;CESC;COAD;ESCA;HNSC;KICH;KIRC;KIRP;LGG;LIHC;LUAD;LUSC;MESO;PAAD;PCPG;PRAD;READ;SARC;STAD;THCA;UCEC" \
#--removeMetastatic "T" --randomize "F" --test "F" --run_query_genes_jointly "T" \
#--protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_vogelstein_drivers_nonsynonymous_uniformHypermutRm_" \
#--expression_df "expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--methylation_df "methylation_M_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--mutation_df "mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--patient_df "_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "F" \
#--meth_type "M" --run_name "top_0.05" --num_pcs 3 --debug "F" --removeCis "F" 

# To run on silent mutations, change the mutation_df
#--mutation_df "mut_count_matrix_silent_uniformHypermutRm_IntersectPatientsWashU.csv" \
# And the associated run_name (where output files will be stored)
#--run_name "top_0.05_silent"

# Excluding cancer types with no drivers that exceed 0.05 mutational freq threshold: 
# LAML;UVM;SKCM;THYM;TGCT 