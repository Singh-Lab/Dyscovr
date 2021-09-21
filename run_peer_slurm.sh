#!/bin/bash
#SBATCH --mem=8192
#SBATCH --qos=1wk
#SBATCH --nodes=1
#SBATCH --job-name=peer_run_fpkm_nocov_noMean_allGenes
#SBATCH --mail-user=scamilli@princeton.edu
#SBATCH --mail-type=fail,time_limit

Rscript run_peer.R -e "expression_fpkm_gn_DF.csv" --incl_CovMat FALSE --inclMean FALSE --top10genes FALSE --cancerType "BRCA"\
--expType "FPKM" 