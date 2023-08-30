#!/bin/bash
#SBATCH --qos=1wk
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64GB
#SBATCH --job-name=create_lm_input_tables
#SBATCH --mail-user=scamilli@princeton.edu
#SBATCH --mail-type=fail,time_limit

module add R/4.0.3

# Create LM input tables for all gene targets, BRCA
#Rscript create_lm_input_tables.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --removeMetastatic "TRUE" --randomize "FALSE" \
#--test "FALSE" --incl_nextMutDriver "FALSE" --run_query_genes_jointly "TRUE" --protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_all_nonsynonymous.csv" \
#--expression_df "expression_quantile_norm_edgeRfilt_sklearn_CancerOnly_IntersectPatientsWashU.csv" \
#--methylation_df "methylation_M_CancerOnly_IntersectPatientsWashU.csv" --methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatientsWashU.csv" \
#--mutation_targ_df "mut_count_matrix_nonsynonymous_IntersectPatientsWashU.csv"  --cna_df "CNA_AllGenes_CancerOnly_IntersectPatientsWashU.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_IntersectPatientsWashU.csv" --target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" \
#--run_name "top_0.05" --num_PEER 0 --num_pcs 3 --debug "FALSE" --useNumFunctCopies "FALSE" --removeCis "FALSE" \
   
# Create LM input tables for all gene targets, within a BRCA subtype 
#Rscript create_lm_input_tables.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --removeMetastatic "TRUE" --randomize "FALSE" \
#--test "FALSE" --incl_nextMutDriver "FALSE" --run_query_genes_jointly "TRUE" --protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_all_nonsynonymous_lumB.csv" \
#--patientsOfInterest "Luminal.B_patient_ids.txt" --patientsOfInterestLabel "LumB" \
#--expression_df "expression_quantile_norm_edgeRfilt_sklearn_CancerOnly_IntersectPatientsWashU.csv" \
#--methylation_df "methylation_M_CancerOnly_IntersectPatientsWashU.csv" --methylation_df_meQTL "methylation_M_CancerOnly_IntersectPatientsWashU.csv" \
#--mutation_targ_df "mut_count_matrix_nonsynonymous_IntersectPatientsWashU.csv" --cna_df "CNA_AllGenes_CancerOnly_IntersectPatientsWashU.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_IntersectPatientsWashU.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "top_0.05" \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --useNumFunctCopies "FALSE" --removeCis "FALSE" \


# Create LM input tables for all gene targets, METABRIC
#Rscript create_lm_input_tables.R --dataset "METABRIC" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --removeMetastatic "FALSE" --randomize "FALSE" \
#--test "FALSE" --incl_nextMutDriver "FALSE" --run_query_genes_jointly "TRUE" --protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_all_nonsynonymous.csv" \
#--expression_df "expression_df_log2intensity_filtBySD0.2_IntersectPatients.csv" \
#--methylation_df "methylation_df_RRBS_logit_IntersectPatients.csv" --methylation_df_meQTL "methylation_df_RRBS_logit_IntersectPatients.csv" \
#--mutation_targ_df "mutation_count_df_nonsynonymous_IntersectPatients.csv" --cna_df "cna_df_allgenes_IntersectPatients.csv" \
#--patient_df "patient_sample_df_cibersort_total_frac_IntersectPatients.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "top_0.05" \
#--num_PEER 0 --num_pcs 0 --debug "FALSE" --useNumFunctCopies "FALSE" --removeCis "FALSE" \

# Create LM input tables for all gene targets, Chinese TN
#Rscript create_lm_input_tables.R --dataset "Chinese_TN" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "BRCA" --removeMetastatic "FALSE" --randomize "FALSE" \
#--test "FALSE" --incl_nextMutDriver "FALSE" --run_query_genes_jointly "TRUE" --protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_all_nonsynonymous_uniformHypermutRm.csv" \
#--expression_df "expression_df_FPKM_filtBySD1_uniformHypermutRm_IntersectPatients.csv" --methylation_df "NA" --methylation_df_meQTL "NA" \
#--mutation_targ_df "mutation_count_df_nonsynonymous_uniformHypermutRm_IntersectPatients.csv" --cna_df "cna_df_allgenes_uniformHypermutRm_IntersectPatients.csv" \
#--patient_df "patient_sample_df_uniformHypermutRm_IntersectPatients.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "NA" --run_name "top_0.05" \
#--num_PEER 0 --num_pcs 0 --debug "TRUE" --useNumFunctCopies "FALSE" --removeCis "FALSE" \

# Create LM input tables for all gene targets, pan-cancer
#Rscript create_lm_input_tables.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "ALL" \
#--removeMetastatic "TRUE" --randomize "FALSE" --test "FALSE" --incl_nextMutDriver "FALSE" --run_query_genes_jointly "TRUE" \
#--protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_all_nonsynonymous_uniformHypermutRm.csv" \
#--expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--methylation_df "methylation_M_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" --methylation_df_meQTL "methylation_M_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--mutation_targ_df "mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--patient_df "combined_patient_sample_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "top_0.05" \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --useNumFunctCopies "FALSE" --removeCis "FALSE" \

# Create LM input tables for all gene targets, per-cancer (0.05)
Rscript create_lm_input_tables.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" \
--specificTypes "ACC;BLCA;BRCA;CESC;COAD;ESCA;HNSC;KICH;KIRC;KIRP;LGG;LIHC;LUAD;LUSC;MESO;PAAD;PCPG;PRAD;READ;SARC;STAD;THCA;UCEC" \
--removeMetastatic "TRUE" --randomize "FALSE" --test "FALSE" --incl_nextMutDriver "FALSE" --run_query_genes_jointly "TRUE" \
--protein_ids_df "iprotein_protein_ids_df_gr0.05Freq_all_nonsynonymous_uniformHypermutRm_" \
--expression_df "expression_quantile_norm_sklearn_uniformHypermutRm_IntersectPatientsWashU.csv" \
--methylation_df "methylation_M_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" --methylation_df_meQTL "methylation_M_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
--mutation_targ_df "mut_count_matrix_silent_uniformHypermutRm_IntersectPatientsWashU.csv" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
--patient_df "_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "top_0.05_silent" \
--num_PEER 0 --num_pcs 3 --debug "FALSE" --useNumFunctCopies "FALSE" --removeCis "FALSE" \

# Create LM input tables for all gene targets, per-cancer (0.1)
#Rscript create_lm_input_tables.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" \
#--specificTypes "BRCA;CESC;COAD;ESCA;GBM;HNSC;KIRC;LGG;LIHC;STAD" \
#--removeMetastatic "TRUE" --randomize "FALSE" --test "FALSE" --incl_nextMutDriver "FALSE" --run_query_genes_jointly "TRUE" \
#--protein_ids_df "iprotein_protein_ids_df_gr0.1Freq_all_nonsynonymous_uniformHypermutRm_" \
#--expression_df "expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--methylation_df "methylation_M_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" --methylation_df_meQTL "methylation_M_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--mutation_targ_df "mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--patient_df "_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "top_0.1" \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --useNumFunctCopies "FALSE" --removeCis "FALSE" \

# Create LM input tables for all gene targets, per-cancer (0.15)
#Rscript create_lm_input_tables.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" \
#--specificTypes "BLCA;COAD;GBM;HNSC;LUAD;LUSC;STAD;UCEC" \
#--removeMetastatic "TRUE" --randomize "FALSE" --test "FALSE" --incl_nextMutDriver "FALSE" --run_query_genes_jointly "TRUE" \
#--protein_ids_df "iprotein_protein_ids_df_gr0.15Freq_all_nonsynonymous_uniformHypermutRm_" \
#--expression_df "expression_quantile_norm_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--methylation_df "methylation_M_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" --methylation_df_meQTL "methylation_M_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--mutation_targ_df "mut_count_matrix_nonsynonymous_uniformHypermutRm_IntersectPatientsWashU.csv" --cna_df "CNA_AllGenes_mergedIsoforms_CancerOnly_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--patient_df "_inclSubtypes_uniformHypermutRm_IntersectPatientsWashU.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "top_0.15" \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --useNumFunctCopies "FALSE" --removeCis "FALSE" \

# Excluding cancer types with no drivers that exceed 0.05 mutational freq threshold: LAML;UVM;SKCM;THYM;TGCT 
 
# Create LM input tables for all gene targets, BRCA-BLCA-HNSC with subtypes
#Rscript create_lm_input_tables.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "ALL" \
#--removeMetastatic "TRUE" --randomize "FALSE" --test "FALSE" --incl_nextMutDriver "FALSE" --run_query_genes_jointly "TRUE" \
#--protein_ids_df "iprotein_protein_ids_df_gr0.1Freq_drivers_nonsynonymous_BRCA_BLCA_HNSC.csv" \
#--expression_df "expression_quantile_norm_BRCA_BLCA_HNSC_IntersectPatientsWashU.csv" \
#--methylation_df "methylation_M_CancerOnly_BRCA_BLCA_HNSC_IntersectPatientsWashU.csv" --methylation_df_meQTL "methylation_M_CancerOnly_BRCA_BLCA_HNSC_IntersectPatientsWashU.csv" \
#--mutation_targ_df "mut_count_matrix_nonsynonymous_BRCA_BLCA_HNSC_IntersectPatientsWashU.csv" --cna_df "CNA_AllGenes_CancerOnly_BRCA_BLCA_HNSC_IntersectPatientsWashU.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_BRCA_BLCA_HNSC_inclBRCA.BLCA.HNSCsubtypes_IntersectPatientsWashU.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "top_drivers_0.1_BRCABLCAHNSC" \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --useNumFunctCopies "FALSE" --removeCis "FALSE" \

# Create LM input tables for all gene targets, COAD-READ
#Rscript create_lm_input_tables.R --dataset "TCGA" --QTLtype "eQTL" --tumNormMatched "FALSE" --cancerType "PanCancer" --specificTypes "ALL" \
#--removeMetastatic "TRUE" --randomize "FALSE" --test "FALSE" --incl_nextMutDriver "FALSE" --run_query_genes_jointly "TRUE" \
#--protein_ids_df "iprotein_protein_ids_df_gr0.1Freq_drivers_nonsynonymous_COAD_READ.csv" \
#--expression_df "expression_quantile_norm_COAD_READ_IntersectPatientsWashU.csv" \
#--methylation_df "methylation_M_CancerOnly_COAD_READ_IntersectPatientsWashU.csv" --methylation_df_meQTL "methylation_M_CancerOnly_COAD_READ_IntersectPatientsWashU.csv" \
#--mutation_targ_df "mut_count_matrix_nonsynonymous_COAD_READ_IntersectPatientsWashU.csv"  --cna_df "CNA_AllGenes_CancerOnly_COAD_READ_IntersectPatientsWashU.csv" \
#--patient_df "combined_patient_sample_cibersort_total_frac_COAD_READ_MSIsubtypes_IntersectPatientsWashU.csv" \
#--target_df "allgene_targets.csv" --cna_bucketing "rawCNA" --meth_bucketing "FALSE" --meth_type "M" --run_name "top_drivers_0.1_COADREAD" \
#--num_PEER 0 --num_pcs 3 --debug "FALSE" --useNumFunctCopies "FALSE" --removeCis "FALSE" \

