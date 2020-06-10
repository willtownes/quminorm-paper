#!/bin/bash
#SBATCH -n 16 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 90 # Runtime in minutes
#SBATCH --mem=128000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[YOUREMAIL]

#pass base path as first command line argument eg ./real/macosko_2015
#pass SingleCellExperiment file relative to base path as second arg
#eg, data/01_sce_all_genes_all_cells.rds
Rscript ./real/util/mles_qumi.R $1 $2
