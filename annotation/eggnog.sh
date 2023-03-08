#!/bin/bash

#SBATCH -J eggNOG                                # Job name
#SBATCH -o eggNOG.out                            # File to which out will be written
#SBATCH -e eggNOG.err                            # File to which err will be written
#SBATCH -N 1                                     # Ensure that all cores are on one machine
#SBATCH -n 8                                     # Number of cores/cpus
#SBATCH -t 07-00:00                              # Runtime in DD-HH:MM
#SBATCH --mem 16000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p shared                                # Partition shared, serial_requeue, unrestricted, test

module load Anaconda3/2020.11
source activate eggNOG

emapper.py -i /path/to/predicted_gene.fasta \
--itype proteins \
-o taxa_func_anno \
--cpu 8 \
--tax_scope Bilateria \
--excel \
--data_dir /path/to/eggnog_database \
--pident 40 \
--query_cover 20 \
--subject_cover 20
