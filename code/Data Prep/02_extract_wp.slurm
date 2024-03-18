#!/bin/bash
#SBATCH -J rjob
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-08:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p test             # Partition to submit to
#SBATCH --mem=20000          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o rjob_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e rjob_%j.err  # File to which STDERR will be written, %j inserts jobid

module load gdal/3.2.2-fasrc01
module load gcc/9.3.0-fasrc01  R/4.0.5-fasrc02
module load udunits/2.2.26-fasrc01
module load geos/3.9.1-fasrc01

export R_LIBS_USER=$HOME/apps/R_4.0.5:$R_LIBS_USER


R CMD BATCH --quiet --no-restore --no-save sim20dp.R

