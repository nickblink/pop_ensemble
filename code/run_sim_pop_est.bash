#!/bin/bash
#SBATCH -p shared
#SBATCH -t 0-01:00
#SBATCH -c 10
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH -e errors_outputs/%x_%A_%a.sim.err
#SBATCH -o errors_outputs/%x_%A_%a.sim.out

module load R/4.3.3-fasrc01
export R_LIBS_USER=$HOME/apps/R_version:$R_LIBS_USER
Rscript code/simulation_main_CAR.R $1 ${SLURM_ARRAY_TASK_ID}