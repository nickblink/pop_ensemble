#!/bin/bash
#SBATCH -p shared
#SBATCH -t 0-10:00
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH -e errors_outputs/%A.sim.err
#SBATCH -o errors_outputs/%A.sim.out

module load R/4.3.3-fasrc01
export R_LIBS_USER=$HOME/apps/R-4.3.3:$R_LIBS_USER
Rscript code/02_extract_wp_Link_cluster.R