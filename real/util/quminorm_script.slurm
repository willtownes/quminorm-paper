#!/bin/bash
#SBATCH -n 16 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 360 # Runtime in minutes
#SBATCH --mem=96000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[YOUR EMAIL]

#module load R_packages

#first command line arg is dataset directory eg ./real/macosko_2015
#second arg is input RDS file relative to the dataset directory (arg 1)
#third arg is output RDS file relative to the dataset dir
#both 2nd and 3rd args are RDS files containing a SingleCellExperiment object

Rscript ./real/util/quminorm_script.R $1 $2 $3
