#!/bin/bash

#SBATCH -J Repmasker                             # Job name
#SBATCH -o Repmasker.%A.out                      # File to which out will be written
#SBATCH -e Repmasker.%A.err                      # File to which err will be written
#SBATCH -N 1                                     # Ensure that all cores are on one machine
#SBATCH -n 16                                    # Number of cores/cpus
#SBATCH -t 07-00:00                              # Runtime in DD-HH:MM
#SBATCH --mem 20000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p shared                                # Partition shared, serial_requeue, unrestricted, test


singularity exec --cleanenv /n/singularity_images/informatics/dfam-tetools/dfam-tetools_1.4.sif RepeatMasker -lib /path/to/RepeatModel_library /path/to/genome_fasta
