#!/bin/bash
#SBATCH -J fasterq-dump # A single job name for the array
#SBATCH -n 16 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 01:01:00 # HH:MM:SS
#SBATCH --open-mode=append
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[YOUR EMAIL]

#load conda environment that includes sra-tools, pigz, etc
module load anaconda3
conda activate fwt

#create scratch directory for fast storage of intermediate files
export scratch=/tmp/${srr}
mkdir -p $scratch

#read the SRA ID from file. SLURM_ARRAY_TASK_ID is a number 1,2,3,...
export srr=$(sed -n ${SLURM_ARRAY_TASK_ID}p extdata/sraFiles.txt)

fasterq-dump $srr -O data/original/fastq -t $scratch -e 16
pigz data/original/fastq/${srr}*.fastq
