#!/bin/bash

#SBATCH -J blob-mapILL                          # Job name 
#SBATCH -o blob-mapILL.%a.%A.out                # File to which stdout will be written
#SBATCH -e blob-mapILL.%a.%A.err                # File to which stderr will be written
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 25G                               # Memory for all cores in Mbytes
#SBATCH -p test,shared,giribet                  # Partition

module load Anaconda3/2020.11
source activate blob

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../../spp-names)

mkdir -p blobtools

illumina1=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_illumina1.fastq.gz
illumina2=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_illumina2.fastq.gz

minimap2 -ax sr \
         -t $SLURM_NTASKS hypo/whole_genome.h.fa \
         $illumina1 $illumina2 \
         | samtools sort -@$SLURM_NTASKS -O BAM -o blobtools/illumina.reads.bam -
         
mv ./blob-mapILL.$SLURM_ARRAY_TASK_ID.* blobtools/
