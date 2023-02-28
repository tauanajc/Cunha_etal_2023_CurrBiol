#!/bin/bash

#SBATCH -J Medaka                               # Job name 
#SBATCH -o Medaka.%a.%A.out                     # File to which stdout will be written
#SBATCH -e Medaka.%a.%A.err                     # File to which stderr will be written
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 100G                              # Memory for all cores in Mbytes
#SBATCH -p test,shared,giribet                  # Partition

source activate medaka

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../../spp-names)

ont=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_ont.fastq.gz
DRAFT=$1
OUTDIR=medaka

medaka_consensus -i $ont -d ${DRAFT} -o ${OUTDIR} -t $SLURM_NTASKS -m r941_min_high_g360

mv Medaka.$SLURM_ARRAY_TASK_ID.* medaka/

