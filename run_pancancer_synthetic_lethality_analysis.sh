#!/bin/bash
#SBATCH --qos=1day
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH --job-name=synLeth_analysis
#SBATCH --mail-user=scamilli@princeton.edu
#SBATCH --mail-type=fail,time_limit

module add R/4.3.1

# Run the the synthetic lethality analysis
Rscript pancancer_synthetic_lethality_analysis.R