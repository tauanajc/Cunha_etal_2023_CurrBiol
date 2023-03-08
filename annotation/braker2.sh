#!/bin/bash

#SBATCH -J braker2                              # Job name
#SBATCH -o braker2.%A.out                       # File to which stdout will be written
#SBATCH -e braker2.%A.err                       # File to which stderr will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 8                                    # Number of cores/cpus
#SBATCH -t 07-00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 30000                             # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p shared                       		# Partition

singularity exec --no-home \
                --home /opt/gm_key \
                --cleanenv \
                --env AUGUSTUS_CONFIG_PATH=/n/holylfs04/LABS/giribet_lab/Lab/adlord/programs/Augustus/config \
                /n/singularity_images/informatics/braker/braker_2.1.6_5-2022-05-04.sif \
                braker.pl --species=species_rna --genome=/path/to/genome/ --bam=/path/to/mapped/rna-seq --softmasking
