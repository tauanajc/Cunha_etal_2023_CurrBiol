#!/bin/bash

#SBATCH -J blob-mapONT                          # Job name 
#SBATCH -o blob-mapONT.%a.%A.out                # File to which stdout will be written
#SBATCH -e blob-mapONT.%a.%A.err                # File to which stderr will be written
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 25G                               # Memory for all cores in Mbytes
#SBATCH -p test,shared,giribet                  # Partition

module load Anaconda3/2020.11
source activate blob

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../../spp-names)

mkdir -p blobtools

ont=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_ont.fastq.gz

minimap2 -ax map-ont \
         -t $SLURM_NTASKS hypo/whole_genome.h.fa \
         $ont > blobtools/ont.reads.sam
cat blobtools/ont.reads.sam | samtools sort -O BAM -o blobtools/ont.reads.bam -        
 
mv ./blob-mapONT.$SLURM_ARRAY_TASK_ID.* blobtools/

