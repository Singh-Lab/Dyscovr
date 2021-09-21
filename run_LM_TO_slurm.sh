#!/bin/bash
#SBATCH --mem=8192
#SBATCH --qos=1wk
#SBATCH --nodes=1
#SBATCH --job-name=linear_model_TP53test
#SBATCH --mail-user=scamilli@princeton.edu
#SBATCH --mail-type=fail,time_limit

Rscript linear_model_TO.R 