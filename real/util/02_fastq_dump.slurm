#!/bin/bash
#SBATCH -n 8 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 360 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=24000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[YOUR EMAIL]
module load sratoolkit
module load parallel

#cat extdata/sraFiles.txt | parallel -j 8 fastq-dump -I --split-files -O data/original/fastq --gzip data/original/sra/{}
#or to avoid useless cat
parallel -j 8 fastq-dump -I --split-files -O data/original/fastq --gzip data/original/sra/{} < extdata/sraFiles.txt
