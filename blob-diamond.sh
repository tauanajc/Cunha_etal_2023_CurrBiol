#!/bin/bash

#SBATCH -J blob-diamond                         # Job name 
#SBATCH -o blob-diamond.%a.%A.out               # File to which stdout will be written
#SBATCH -e blob-diamond.%a.%A.err               # File to which stderr will be written
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 15G                               # Memory for all cores in Mbytes
#SBATCH -p test,shared,giribet                  # Partition

module load Anaconda3/2020.11
source activate blob

mkdir -p blobtools

diamond blastx \
        --query hypo/whole_genome.h.fa \
        --db /n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/blobtoolkit/blob_databases/uniprot_2021_07/reference_proteomes.dmnd \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-25 \
        --threads $SLURM_NTASKS \
        > blobtools/diamond.out

mv ./blob-diamond.$SLURM_ARRAY_TASK_ID.* blobtools/
