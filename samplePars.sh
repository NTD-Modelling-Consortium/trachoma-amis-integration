#!/bin/bash

#SBATCH --output log/sample_pars.out-%A_%a
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00

module purge
module load R/4.3.2-gfbf-2023a
stdbuf -i0 -o0 -e0 command

Rscript preprocess_for_projections.R



