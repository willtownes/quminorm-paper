#!/bin/bash
#SBATCH -n 8 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 60 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=24000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[YOUR EMAIL]

wget -c -nc -i extdata/sraFilesPath.txt -P data/original/sra

#module load sratoolkit
#time prefetch --transport ascp --max-size 50G -O data/original/sra --option-file extdata/sraFiles.txt

