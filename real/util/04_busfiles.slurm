#!/bin/bash
#SBATCH -J busfiles # job name
#SBATCH -n 4 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 60 # Runtime in minutes
#SBATCH --mem=24G # Total memory (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this file
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[YOUR EMAIL]

bp=./data/original/kallisto
while read srr; do
  pth=$bp/$srr
  if [[ -f extdata/barcodes.txt ]]
  then
    bustools correct -w extdata/barcodes.txt -o $pth/output_corrected.bus $pth/output.bus
    bustools sort -t4 -o $pth/output_sort.bus $pth/output_corrected.bus
  else
    bustools sort -t4 -o $pth/output_sort.bus $pth/output.bus
  fi
  bustools text -o $pth/output_sort.txt $pth/output_sort.bus
done < ./extdata/sraFiles.txt

#tar everything except the .bus files for convenient download/upload
tar -zc -f data/original/kallisto.tar.gz --exclude="$bp/*.bus" $bp/*

