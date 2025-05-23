#!/bin/bash
#SBATCH -p shared,hsph,serial_requeue
#SBATCH -t 0-5:00
#SBATCH -c 10
#SBATCH --mem=184G
#SBATCH --mail-type=END
#SBATCH -e errors_outputs/%x_%A_%a.sim.err
#SBATCH -o errors_outputs/%x_%A_%a.sim.out

module load R/4.3.3-fasrc01
export R_LIBS_USER=$HOME/apps/R-4.3.3:$R_LIBS_USER
Rscript code/real_data_main_CAR.R $1