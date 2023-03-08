#!/bin/bash

#SBATCH -J hisat2                                # Job name
#SBATCH -o hisat2.out                            # File to which stdout will be written
#SBATCH -e hisat2.err                            # File to which stderr will be written
#SBATCH -N 1                                     # Ensure that all cores are on one machine
#SBATCH -n 8                                     # Number of cores/cpus
#SBATCH -t 07-00:00                              # Runtime in DD-HH:MM
#SBATCH --mem 16000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p shared,giribet                        # Partition

module load hisat2/2.1.0-fasrc01
module load samtools/1.10-fasrc01

read1=/path/to/rnaseq_R1.fq.gz
read2=/path/to/rnaseq_R2.fq.gz

## build index
hisat2-build /path/to/masked/fasta masked.index

## alignment
hisat2 -p 8 -x masked.index -1 $read1 -2 $read2 -S rnamapped.sam

## convert to bam and sort
samtools view --threads 8 -b rnamapped.sam | samtools sort -o rnamapped.bam --threads 8