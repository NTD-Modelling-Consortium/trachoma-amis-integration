#!/bin/bash

#SBATCH --output log/mtp-trachoma.out-%A_%a
#SBATCH --array=1-564
#SBATCH -p short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=30:00:00

#1-564

# Change directory
cd trachoma-amis-integration

# Load modules
module purge
module load R/4.3.2-gfbf-2023a
#module load Python/3.11.3-GCCcore-12.3.0
source ../trachoma-venv/bin/activate
unset RETICULATE_PYTHON

stdbuf -i0 -o0 -e0 command

# Run R script
Rscript trachoma_fitting.R 12
