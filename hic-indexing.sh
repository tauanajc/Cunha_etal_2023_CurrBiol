#!/bin/bash

#SBATCH -J HiC-indexing                         # Job name 
#SBATCH -e indexing.%A.err                      # File to which stderr will be written
#SBATCH -o indexing.%A.out                      # File to which stdout will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 1                                    # Number of cores/cpus
#SBATCH -t 00-00:30                             # Runtime in DD-HH:MM
#SBATCH --mem 1000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p test                                 # Partition shared, serial_requeue, unrestricted, test

module load Anaconda3/2020.11
source activate hic

mkdir -p hic/mapping

# Create assembly.fai index
samtools faidx $1
rsync -a blobtools/decontaminated.fasta.fai hic/mapping/

# Create bwa index
bwa index -a bwtsw -p hic-bwa-index $1
mv hic-bwa-index.* hic/mapping/

rm indexing.$SLURM_JOB_ID.out
mv indexing.$SLURM_JOB_ID.err hic/mapping/

