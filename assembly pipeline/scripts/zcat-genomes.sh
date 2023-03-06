#!/bin/bash

#SBATCH -J zcat_reads                           # Job name 
#SBATCH -o zcat_reads.%a.%A.out                 # File to which stdout will be written
#SBATCH -N 1
#SBATCH -n 1                                    # number of ranks
#SBATCH -t 1-00:00:00                           # Runtime in DD-HH:MM
#SBATCH --mem 10000                             # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p giribet,shared                       # Partition general, serial_requeue, unrestricted, interact

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../spp-names)

illumina1=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_illumina1.fastq.gz
illumina2=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_illumina2.fastq.gz
ont=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_ont.fastq.gz

if [ ! -f $illumina1 ]; then
zcat ../illumina/$spp_dir/trimG/*val_1.fq.gz | gzip > $illumina1
echo $(date +"%T") $illumina1
else echo $illumina1 already existed- skip.
fi

if [ ! -f $illumina2 ]; then
zcat ../illumina/$spp_dir/trimG/*val_2.fq.gz | gzip > $illumina2
echo $(date +"%T") $illumina2
else echo $illumina2 already existed- skip.
fi

if [ ! -f $ont ]; then
zcat ../basecalled-ont/$spp_dir/fastq/*.fastq.gz | gzip > $ont
echo $(date +"%T") $ont
else echo $ont already existed- skip.
fi

