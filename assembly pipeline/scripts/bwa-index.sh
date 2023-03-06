#!/bin/bash

#SBATCH -J bwaIndex                             # Job name 
#SBATCH -o bwaIndex.%A.out                      # File to which stdout will be written
#SBATCH -e bwaIndex.%A.err                      # File to which stderr will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 1                                    # Number of cores/cpus
#SBATCH -t 00-00:30                             # Runtime in DD-HH:MM
#SBATCH --mem 1000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p test                                 # Partition shared, serial_requeue, unrestricted, test

module load bwa/0.7.17-fasrc01

mkdir -p bwa-index

bwa index -p bwa-index $1
mv bwa-index.* bwa-index/

mv bwaIndex.$SLURM_JOB_ID.* bwa-index/

