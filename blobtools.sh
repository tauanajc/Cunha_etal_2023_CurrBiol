#!/bin/bash

#SBATCH -J blobtools                            # Job name 
#SBATCH -o blobtools.%a.%A.out                  # File to which stdout will be written
#SBATCH -e blobtools.%a.%A.err                  # File to which stderr will be written
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 10000                             # Memory for all cores in Mbytes
#SBATCH -p test,shared,giribet                  # Partition

module load Anaconda3/2020.11
source activate blob

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../../spp-names)

# merge the two alignment bam files from illumina and ont: out.bam in1.bam in2.bam
#samtools merge -@$SLURM_NTASKS blobtools/${spp_dir}.reads.bam blobtools/*.bam

# blobtools
blobtools create \
    --threads $SLURM_NTASKS \
    --fasta hypo/whole_genome.h.fa \
    --taxid $1 \
    --taxdump /n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/blobtoolkit/blob_databases/taxdump_2021_07 \
    --hits blobtools/blast.out \
    --hits blobtools/diamond.out \
    --taxrule bestsumorder \
    --cov blobtools/illumina.reads.bam \
    --cov blobtools/ont.reads.bam \
    --busco blobtools/assembly_metazoa_odb10/run_metazoa_odb10/full_table.tsv \
    --busco blobtools/assembly_eukaryota_odb10/run_eukaryota_odb10/full_table.tsv \
    --busco blobtools/assembly_nematoda_odb10/run_nematoda_odb10/full_table.tsv \
    --busco blobtools/assembly_arthropoda_odb10/run_arthropoda_odb10/full_table.tsv \
    blobtools/BlobDir

mv ./blobtools.$SLURM_ARRAY_TASK_ID.* blobtools/

