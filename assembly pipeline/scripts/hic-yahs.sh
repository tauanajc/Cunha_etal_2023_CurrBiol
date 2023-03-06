#!/bin/bash

#BATCH -J yahs                                 # Job name 
#SBATCH -o yahs.%A.out                          # File to which stdout will be written
#SBATCH -e yahs.%A.err                          # File to which stderr will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 8                                    # Number of cores/cpus
#SBATCH -t 07-00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 184000                            # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p shared,giribet                       # Partition shared, serial_requeue, unrestricted, test

module load Java/1.8
module load parallel

LABEL=$(head -n1 hic/fastq-names-hic.txt | rev | cut -d_ -f3- | rev)
YAHS="/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/yahs/"
JUICER="/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/juicer_tools_1.22.01.jar"

mkdir -p hic/yahs

### convert bam to bed file and sort by read name
bamToBed -i hic/mapping/final-deduplicated-bam/$LABEL.final.bam > hic/$LABEL.bed
sort -k 4 hic/$LABEL.bed > tmp && mv tmp hic/$LABEL.bed

### run yahs
$YAHS/yahs blobtools/decontaminated.fasta hic/$LABEL.bed \
-o hic/yahs/yahs.out

### create .hic contact map file for visualization
cd hic/yahs

# juicer pre to create intermediary file
($YAHS/juicer pre yahs.out.bin yahs.out_scaffolds_final.agp ../mapping/decontaminated.fasta.fai | \
 sort -k2,2d -k6,6d -T ./ --parallel=$SLURM_NTASKS -S32G | \
 awk 'NF' > alignments_sorted.txt.part) 2>&1 | tee juicer1.log && (mv alignments_sorted.txt.part alignments_sorted.txt)

# create file for scaffold sizes
cat juicer1.log | grep "PRE_C_SIZE:" | cut -d' ' -f 2,3 > scaffolds_final.chrom.sizes

# juicer pre to create .hic file - longest and more memory intensive step of the script
(java -jar -Xmx184G $JUICER pre --threads 8 alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes) 2>&1 | tee juicer2.log && (mv out.hic.part out.hic)

mv ../../yahs.$SLURM_JOB_ID.* ./
