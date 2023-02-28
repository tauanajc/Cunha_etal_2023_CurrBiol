#!/bin/bash

#SBATCH -J flye                                 # Job name 
#SBATCH -o flye.%a.%A.out                       # File to which stdout will be written
#SBATCH -e flye.%a.%A.err                       # File to which stderr will be written
#SBATCH -n 32                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 3-00:00                              # Runtime in DD-HH:MM
#SBATCH --mem 184G                              # Memory for all cores in Mbytes
#SBATCH -p shared,giribet                       # Partition

source activate genomes

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../spp-names)
outfolder=flye-5kbp

cd $spp_dir

flye --nano-raw ../$1/*.fastq.gz --out-dir $outfolder --threads $SLURM_NTASKS
# --resume

mv ../flye.$SLURM_ARRAY_TASK_ID.* $outfolder/

