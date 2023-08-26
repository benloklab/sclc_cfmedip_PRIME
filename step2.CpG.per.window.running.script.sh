#!/bin/bash
# the following are sbatch parameters
#SBATCH -t 60:00:00
#SBATCH -p veryhimem
#SBATCH --mem=180G
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o %x-%j.out


module load R/3.6.1

echo 'Starting R driver script (R version 3.6.1)'

Rscript CpG.per.window.R ${1}

echo 'Finished running driver script'
