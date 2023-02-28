#!/bin/bash

#SBATCH -J HyPo                                 # Job name 
#SBATCH -o HyPo.%a.%A.out                       # File to which stdout will be written
#SBATCH -e HyPo.%a.%A.err                       # File to which stderr will be written
#SBATCH -n 32                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 184G                              # Memory for all cores in Mbytes
#SBATCH -p test,shared,giribet                  # Partition

module load Anaconda3
source activate hypo

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../../spp-names)

illumina1=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_illumina1.fastq.gz
illumina2=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_illumina2.fastq.gz
ont=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_ont.fastq.gz
DRAFT=purge_dups/purged.fa
genomesize=$1

mkdir -p hypo
cd hypo

##### Mapping the short reads to contigs #####
minimap2 --secondary=no --MD -ax sr -t $SLURM_NTASKS \
../$DRAFT $illumina1 $illumina2 | samtools view -Sb - > mapped-sr.bam

samtools sort -@$SLURM_NTASKS -o mapped-sr.sorted.bam mapped-sr.bam
samtools index mapped-sr.sorted.bam
rm mapped-sr.bam

##### HyPo #####
echo -e "$illumina1\n$illumina2" > illumina_file_names.txt
samtools depth -a mapped-sr.sorted.bam > coverage.txt
coverage=$(cat coverage.txt | awk '{ sum+=$3 } END { printf "%.0f", sum/NR }')

hypo -d ../$DRAFT -r @illumina_file_names.txt -s $genomesize -c $coverage \
-b mapped-sr.sorted.bam -t $SLURM_NTASKS -o whole_genome.h.fa

mv ../HyPo.$SLURM_ARRAY_TASK_ID.* ./

