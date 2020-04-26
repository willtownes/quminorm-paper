#!/bin/bash
#SBATCH -J kallisto # A single job name for the array
#SBATCH -N 1 # Force all computation on a single node
#SBATCH -n 1 # Number of tasks requested (-N to control number of nodes)
#SBATCH -c 8 # Number of cores per task
#SBATCH -t 60 # Runtime in minutes
#SBATCH --mem=64G # Memory per node (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o hostname_%A_%a.out # Standard out goes to this file
#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[YOUR EMAIL]

sed -n ${SLURM_ARRAY_TASK_ID}p extdata/kallisto_cmd.txt | bash
