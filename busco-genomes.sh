#!/bin/bash

#SBATCH -J Busco                                # Job name 
#SBATCH -o Busco.%A.out                         # File to which stdout will be written
#SBATCH -e Busco.%A.err                         # File to which stderr will be written
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 184000                            # Memory for all cores in Mbytes
#SBATCH -p test,shared,giribet                  # Partition

module load Anaconda3
source activate busco545

#busco -i $1 -l metazoa_odb10 -o busco_metazoa -m genome --cpu $SLURM_NTASKS
busco -i $1 -l $2 -o out_$2 -m genome --cpu $SLURM_NTASKS
#If interactive, append to end of command: 2>&1 | tee busco-metazoa.log

