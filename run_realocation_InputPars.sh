#!/bin/bash

#SBATCH --output log/run_realoc_InputPars_MTP.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00


module purge
module load R/4.3.2-gfbf-2023a
stdbuf -i0 -o0 -e0 command

Rscript realocate_InputPars_MTP.R


