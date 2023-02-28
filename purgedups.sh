#!/bin/bash

#SBATCH -J PurgeDups                            # Job name 
#SBATCH -o PurgeDups.%a.%A.out                  # File to which stdout will be written
#SBATCH -e PurgeDups.%a.%A.err                  # File to which stderr will be written
#SBATCH -n 32                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-02:00                              # Runtime in DD-HH:MM
#SBATCH --mem 10000                             # Memory for all cores in Mbytes
#SBATCH -p test,shared,giribet                  # Partition

source activate genomes

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../../spp-names)

ont=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_ont.fastq.gz
ASSEMBLY=$1

mkdir -p purge_dups
cd purge_dups

##### STEP 1 #####
minimap2 -x map-ont ../$ASSEMBLY $ont -t $SLURM_NTASKS | gzip -c - > mini2-reads.paf.gz

/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/purge_dups/bin/pbcstat *.paf.gz
    # produces PB.base.cov and PB.stat files
/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/purge_dups/bin/calcuts \
PB.stat > cutoffs 2>calcults.log

##### STEP 1 #####
/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/purge_dups/bin/split_fa ../$ASSEMBLY > assembly.split

minimap2 -x asm5 -DP assembly.split assembly.split -t $SLURM_NTASKS | gzip -c - > assembly.split.self.paf.gz

##### STEP 2 #####
/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/purge_dups/bin/purge_dups -2 \
-T cutoffs -c PB.base.cov assembly.split.self.paf.gz > dups.bed 2> purge_dups.log

##### STEP 3 #####
/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/purge_dups/bin/get_seqs -e dups.bed ../$ASSEMBLY


mv ../PurgeDups.$SLURM_ARRAY_TASK_ID.* ./

