#!/bin/bash

#SBATCH -J RepModler3                            # Job name
#SBATCH -o RepModler3.%A.out                     # File to which out will be written
#SBATCH -e RepModler3.%A.err                     # File to which err will be written
#SBATCH -N 1                                     # Ensure that all cores are on one machine
#SBATCH -n 16                                    # Number of cores/cpus
#SBATCH -t 07-00:00                              # Runtime in DD-HH:MM
#SBATCH --mem 48000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p shared                                # Partition shared, serial_requeue, unrestricted, test

 
singularity exec --cleanenv /n/singularity_images/informatics/dfam-tetools/dfam-tetools_1.4.sif BuildDatabase \
-name species.DB -engine rmblast /path/to/genome_assembly

singularity exec --cleanenv /n/singularity_images/informatics/dfam-tetools/dfam-tetools_1.4.sif RepeatModeler \
-database species.DB -pa 16 -LTRStruct


