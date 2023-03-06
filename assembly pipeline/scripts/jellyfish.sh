#!/bin/bash

#SBATCH -J Jellyfish                            # Job name 
#SBATCH -o Jellyfish.%a.%A.out                  # File to which stdout will be written
#SBATCH -e Jellyfish.%a.%A.err                  # File to which stderr will be written
#SBATCH -t 07-00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 15000                             # Memory for all cores in Mbytes
#SBATCH -p giribet,shared                       # Partition

module load jellyfish/2.2.5-fasrc01

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../spp-names)

mkdir -p $spp_dir/jellyfish
cd $spp_dir/jellyfish

# Count all k-mers
ls ../fastp/*.clean.fq.gz | xargs -n 1 echo zcat > generators
jellyfish count -C -s 2G -t 4 --disk -m 21 -Q "5" -o out21kmer.jf -g generators -G 4

# Make histogram
jellyfish histo -t 8 out21kmer.jf > out21kmer-${spp_dir}.histo

mv ../../Jellyfish.$SLURM_ARRAY_TASK_ID.* ./

