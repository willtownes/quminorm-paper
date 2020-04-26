#!/bin/bash
#SBATCH -n 4 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 180 # Runtime in minutes
#SBATCH --open-mode=append
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[YOUR EMAIL]

#first command line arg is input RDS file containing SingleCellExperiment
#second arg is output txt file with all distances

Rscript ./real/util/qumi_dist_compare_script.R $1 $2
