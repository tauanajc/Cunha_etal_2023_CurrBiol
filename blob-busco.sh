#!/bin/bash

#SBATCH -J blob-busco                           # Job name 
#SBATCH -o blob-busco.%a.%A.out                 # File to which stdout will be written
#SBATCH -e blob-busco.%a.%A.err                 # File to which stderr will be written
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 10G                               # Memory for all cores in Mbytes
#SBATCH -p test,shared,giribet                  # Partition

module load Anaconda3/2020.11
source activate blob

mkdir -p blobtools
cd blobtools

busco \
    -i ../hypo/whole_genome.h.fa \
    -o assembly_$1 \
    -l /n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/blobtoolkit/blob_databases/busco_2022_02/$1 \
    -m geno \
    -c $SLURM_NTASKS

mv ../blob-busco.$SLURM_ARRAY_TASK_ID.$SLURM_JOB_ID.* ./

