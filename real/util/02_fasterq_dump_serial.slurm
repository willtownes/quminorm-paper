#!/bin/bash
#SBATCH -J fasterq-dump # A single job name for the array
#SBATCH -n 16 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 02:00:00 # HH:MM:SS
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

export pfx=data/original/fastq

while read srr; do
  echo "$srr"
  #check if the gzipped fastq is already in the output directory
  if [[ -f ${pfx}/${srr}_1.fastq.gz && -f ${pfx}/${srr}_2.fastq.gz ]]
  then
    echo "output files already in directory so no processing done"
  else
    fasterq-dump $srr -e 16 -O $pfx -t $scratch
    pigz ${pfx}/${srr}*.fastq
  fi
done < extdata/sraFiles.txt
