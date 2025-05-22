#!/bin/bash

#SBATCH --output log/proj_trachoma.out-%A_%a
#SBATCH --array=0-2788
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=20:00

#2789 IUs total (indexes from 0 to 2788)

module purge
module load Python/3.11.3-GCCcore-12.3.0 
stdbuf -i0 -o0 -e0 command

cd ntd-model-trachoma/
source ../trachoma-venv/bin/activate

python ../trachoma-amis-integration/RunProjectionsTo2026.py ${SLURM_ARRAY_TASK_ID}
