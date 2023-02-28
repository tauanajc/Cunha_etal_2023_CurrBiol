#!/bin/bash

#SBATCH -J Racon                                # Job name 
#SBATCH -o Racon.%a.%A.out                      # File to which stdout will be written
#SBATCH -e Racon.%a.%A.err                      # File to which stderr will be written
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 70G                               # Memory for all cores in Mbytes
#SBATCH -p test,shared,giribet                  # Partition

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../../spp-names)

ont=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_ont.fastq.gz

mkdir -p racon

/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/racon/build/bin/racon -t $SLURM_NTASKS \
-m 8 -x -6 -g -8 -w 500 \
$ont bwa-mapping/mapping.sam $1 > racon/racon.fasta

mv racon.$SLURM_ARRAY_TASK_ID.* racon/

