#!/bin/bash
#SBATCH -n 8 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 360 # Runtime in minutes
#SBATCH --mem=48000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[YOUR EMAIL]

#pass base path as first command line argument eg ./real/macosko_2015
Rscript ./real/util/05_bus2genecount.R $1
