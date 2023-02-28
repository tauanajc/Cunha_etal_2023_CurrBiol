#!/bin/bash

#BATCH -J Quast                                # Job name 
#SBATCH -o Quast.%A.out                         # File to which stdout will be written
#SBATCH -e Quast.%A.err                         # File to which stderr will be written
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 1000                              # Memory for all cores in Mbytes
#SBATCH -p test,shared,giribet                  # Partition

module load Anaconda3
source activate quast

quast.py -t $SLURM_NTASKS -o quast "$@"
