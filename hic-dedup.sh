#!/bin/bash

#SBATCH -J HiC-dedup                            # Job name 
#SBATCH -o dedup.%A.out                         # File to which stdout will be written
#SBATCH -e dedup.%A.err                         # File to which stderr will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 16                                   # Number of cores/cpus
#SBATCH -t 00-08:00                             # Runtime in DD-HH:MM
#SBATCH --mem 50000                             # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p test,giribet                         # Partition shared, serial_requeue, unrestricted, test

module load Java/1.8
module load Anaconda3/2020.11
source activate hic

LABEL=$(head -n1 hic/fastq-names-hic.txt | rev | cut -d_ -f3- | rev)
STATS='/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/arima_mapping/get_stats.pl'
PICARD='/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/picard.jar'

mkdir -p hic/mapping/combined-bam
mkdir -p hic/mapping/final-deduplicated-bam

cd hic/mapping

echo "### Step 1: combine technical replicates"
cd paired-bam
java -Xmx8G -Djava.io.tmpdir=../temp/ -jar $PICARD MergeSamFiles \
INPUT=${LABEL}_S6_L001.bam INPUT=${LABEL}_S6_L002.bam INPUT=${LABEL}_S6_L003.bam INPUT=${LABEL}_S6_L004.bam \
OUTPUT=../combined-bam/$LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

cd ..
echo "### Step 2: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates \
INPUT=combined-bam/$LABEL.bam OUTPUT=final-deduplicated-bam/$LABEL.final.bam \
METRICS_FILE=final-deduplicated-bam/metrics.$LABEL.txt TMP_DIR=tmp-bam \
ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

echo "### Step 3: index bam"
samtools index final-deduplicated-bam/$LABEL.final.bam

echo "### Step 4: bam stats"
perl $STATS final-deduplicated-bam/$LABEL.final.bam > final-deduplicated-bam/$LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"

rm ../../dedup.$SLURM_JOB_ID.out
mv ../../dedup.$SLURM_JOB_ID.err ./

