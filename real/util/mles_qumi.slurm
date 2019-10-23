#!/bin/bash
#SBATCH -n 32 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 360 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=128000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[YOUREMAIL]

module load R_packages

#pass base path as first command line argument eg ./real/macosko_2015
Rscript ./real/util/mles_qumi.R $1
