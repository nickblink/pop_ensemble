#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -t 0-00:15
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH -e errors_outputs/%x_%A_%a.sim.err
#SBATCH -o errors_outputs/%x_%A_%a.sim.out

module load R/4.3.3-fasrc01
export R_LIBS_USER=$HOME/apps/R_version:$R_LIBS_USER
Rscript code/simulation_main_onerun.R