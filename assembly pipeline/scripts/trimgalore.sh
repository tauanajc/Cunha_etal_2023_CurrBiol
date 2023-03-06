#!/bin/bash

#SBATCH -J TrimG                                # Job name 
#SBATCH -o TrimG.%a.%A.out                      # File to which stdout will be written
#SBATCH -e TrimG.%a.%A.err                      # File to which stderr will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 4                                    # Number of cores/cpus
#SBATCH -t 00-08:00                             # Runtime in DD-HH:MM
#SBATCH --mem 200                               # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p shared,giribet,test                  # Partition general, serial_requeue, unrestricted, interact


module load cutadapt/1.8.1-fasrc01 parallel


# Create variable to hold the directory name (which is the species name) and move inside folder
spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../spp-names)
mkdir -p $spp_dir/trimG                         # Create directory for output if it is not already there
cd $spp_dir/trimG


parallel -j $SLURM_NTASKS ../../../../scripts/Trim_Galore/trim_galore --phred33 --gzip --length 50 --retain_unpaired --paired \
${1}{}_R1* ${1}{}_R2* ::: 1 2 3 4


# Move .err and .out files to trimG folder with rest of output
mv ../../TrimG.$SLURM_ARRAY_TASK_ID.* ./

