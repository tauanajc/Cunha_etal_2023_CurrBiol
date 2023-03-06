#!/bin/bash

#SBATCH -J HiC-map                              # Job name 
#SBATCH -o mapping.%a.%A.out                    # File to which stdout will be written
#SBATCH -e mapping.%a.%A.err                    # File to which stderr will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 16                                   # Number of cores/cpus
#SBATCH -t 00-08:00                             # Runtime in DD-HH:MM
#SBATCH --mem 10000                             # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p test,giribet                         # Partition shared, serial_requeue, unrestricted, test

module load Java/1.8
module load Anaconda3/2020.11
source activate hic

NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" hic/fastq-names-hic.txt)
LABEL=$(head -n1 hic/fastq-names-hic.txt | rev | cut -d_ -f3- | rev)
IN_DIR='/n/Giribet_Lab/tauanajc/GASTROPODA/rawdata/2020_03-NovaSeqS4/Tauana'
FILTER='/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/arima_mapping/filter_five_end.pl'
COMBINER='/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/arima_mapping/two_read_bam_combiner.pl'
PICARD='/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/picard.jar'
MAPQ_FILTER=10

mkdir -p hic/mapping/raw-bam
mkdir -p hic/mapping/filtered-bam
mkdir -p hic/mapping/tmp-bam
mkdir -p hic/mapping/paired-bam

cd hic/mapping

echo "### Step 1.A: FASTQ to BAM (1st)"
bwa mem -t $SLURM_NTASKS hic-bwa-index $IN_DIR/$NAME\_R1_001.fastq.gz | samtools view -@ $SLURM_NTASKS -Sb - > raw-bam/$NAME\_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -t $SLURM_NTASKS hic-bwa-index $IN_DIR/$NAME\_R2_001.fastq.gz | samtools view -@ $SLURM_NTASKS -Sb - > raw-bam/$NAME\_2.bam

echo "### Step 2.A: Filter 5' end (1st)"
samtools view -h raw-bam/$NAME\_1.bam | perl $FILTER | samtools view -Sb - > filtered-bam/$NAME\_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
samtools view -h raw-bam/$NAME\_2.bam | perl $FILTER | samtools view -Sb - > filtered-bam/$NAME\_2.bam

echo "### Step 3: Pair reads & mapping quality filter"
perl $COMBINER filtered-bam/$NAME\_1.bam filtered-bam/$NAME\_2.bam samtools $MAPQ_FILTER | samtools view -bS -t decontaminated.fasta.fai - | samtools sort -@ $SLURM_NTASKS -o tmp-bam/$NAME.bam -

echo "### Step 4: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=tmp-bam/$NAME.bam OUTPUT=paired-bam/$NAME.bam ID=$NAME LB=$NAME SM=$LABEL PL=ILLUMINA PU=none

rm ../../mapping.$SLURM_ARRAY_TASK_ID.$SLURM_JOB_ID.out
mv ../../mapping.$SLURM_ARRAY_TASK_ID.$SLURM_JOB_ID.err ./

