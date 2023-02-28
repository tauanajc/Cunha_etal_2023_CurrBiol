#!/bin/bash

#SBATCH -J Guppy                                # Job name 
#SBATCH -o Guppy.%a.%A.out                      # File to which stdout will be written
#SBATCH -e Guppy.%a.%A.err                      # File to which stderr will be written
#SBATCH --gres=gpu:1                            # Number of GPUs
#SBATCH -t 01-00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 3000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p gpu                                  # Partition gpu_test (for testing), gpu

module load cuda/10.1.243-fasrc01

# Create variable to hold the directory name (which is the species name and voucher)
spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../spp-names)

mkdir -p $spp_dir/fastq

fast5folder=$1

/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/ont-guppy-4.5.2/bin/guppy_basecaller \
-i $fast5folder -s $spp_dir/fastq --recursive \
--flowcell FLO-MIN106 --kit SQK-LSK109 --compress_fastq --disable_qscore_filtering \
--device "auto" --num_callers 8 --gpu_runners_per_device 16

mv Guppy.$SLURM_ARRAY_TASK_ID.$SLURM_ARRAY_JOB_ID.* $spp_dir/

