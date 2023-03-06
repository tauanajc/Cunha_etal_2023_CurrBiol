#!/bin/bash

#SBATCH -J blob-blast                           # Job name 
#SBATCH -o blob-blast.%a.%A.out                 # File to which stdout will be written
#SBATCH -e blob-blast.%a.%A.err                 # File to which stderr will be written
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 3-00:00                              # Runtime in DD-HH:MM
#SBATCH --mem 100G                              # Memory for all cores in Mbytes
#SBATCH -p shared,giribet                       # Partition

module load Anaconda3/2020.11
source activate blob

mkdir -p blobtools

blastn -db /n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/blobtoolkit/blob_databases/nt_2021_07/nt \
       -num_threads $SLURM_NTASKS \
       -query hypo/whole_genome.h.fa \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -out blobtools/blast.out

mv ./blob-blast.$SLURM_ARRAY_TASK_ID.* blobtools/
