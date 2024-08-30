#!/bin/bash
#SBATCH -p shared
#SBATCH -t 0-01:00
#SBATCH -c 10
#SBATCH --mem=10G
#SBATCH --mail-type=END
#SBATCH -e errors_outputs/%x_%A_%a.sim.err
#SBATCH -o errors_outputs/%x_%A_%a.sim.out

module load R/4.3.3-fasrc01
export R_LIBS_USER=$HOME/apps/R-4.3.3:$R_LIBS_USER
Rscript code/02_extract_wp_Link_cluster_batch.R ${SLURM_ARRAY_TASK_ID}