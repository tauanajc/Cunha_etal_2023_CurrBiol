#!/bin/bash

#SBATCH -J bwaMap                               # Job name 
#SBATCH -o bwaMap.%A.out                        # File to which stdout will be written
#SBATCH -e bwaMap.%A.err                        # File to which stderr will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 16                                   # Number of cores/cpus
#SBATCH -t 07-00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 10000                             # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p shared,giribet                       # Partition shared, serial_requeue, unrestricted, test

module load bwa/0.7.17-fasrc01

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../../spp-names)
ont=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_ont.fastq.gz

mkdir -p bwa-mapping

bwa mem -t $SLURM_NTASKS -x ont2d \
bwa-index/bwa-index $ont > bwa-mapping/mapping.sam

mv bwaMap.$SLURM_JOB_ID.* bwa-mapping/

