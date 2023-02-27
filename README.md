This repository contains the scripts for reproducing assembly, analyses and figures from Cunha TJ, De Medeiros BAS, Lord A, Giribet G. 2023. Rampant loss of universal metazoan genes revealed by a chromosome-level genome assembly of the parasitic Nematomorpha.

# Genome Assembly

### Folder structure

- Inside project folder, I keep a simple text file (spp\_names) with sample names in one column. I separate names by \_ and voucher by - (this is important for scripts that cut the folder name to find the correct folders with raw data):

```
Acutogordius_australiensis-MCZ152393
Nectonema_munidae-DNA05827
```

For certain jobs below, the array is the line number for the species in the file spp-names. The name will be used to create new folders and file names.

>- spp_names
>- fast5_symlinks/
>- basecalled-ont/
>    - Species1_MCZ/
>        - fastq/
>        - nanoplot/
>    - Species2_MCZ/
>        - fastq/
>        - nanoplot/
>- assemblies/
>- illumina/

### Python environment

Create an enviroment for genome-related tools. At least some of the tools used below are for python 3.<br>
Not everything is python-based, but some steps are, so activate enviroment in the scripts or before running programs in the command line.
```bash
conda create -n genomes -c bioconda python=3 flye nanoplot
source activate genomes
# conda info --envs
```

### SymLinks

For multiple ONT flow cells, put all fast5 files in the same folder as symlinks of the original files. Pass and fail fast5 folders in the flow cell rawdata will have identical file names, therefore create symlinks in separate pass/fail folders per species to avoid losing/overwriting files.

- In **fast5_symlinks**:
```bash
mkdir -p spp_name/fast5_fail
mkdir -p spp_name/fast5_pass
cd spp_name/fast5_fail
ln -s /FULL/PATH/TO/FAST5/FAIL/*.fast5 ./
cd ../fast5_pass
ln -s /FULL/PATH/TO/FAST5/PASS/*.fast5 ./
```

## Basecalling ONT with Guppy

Basecall ONT long reads.

- From folder **basecalled-ont**:

```bash
sbatch --array=1-2 ../../scripts/basecall.sh ../fast5-symlinks/SPECIES-FOLDER
```

This structure allows for multiple flow cells with similar but not identical names to be basecalled into one species folder.

```bash
%%bash
#!/bin/bash

#SBATCH -J Guppy                                # Job name 
#SBATCH --gres=gpu:1                            # Number of GPUs
#SBATCH -t 03-00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 3000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)

module load cuda/10.1.243-fasrc01

# Create variable to hold the directory name (which is the species name and voucher)
spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../spp-names)

mkdir -p $spp_dir/fastq

fast5folder=$1

/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/ont-guppy-4.5.2/bin/guppy_basecaller \
-i $fast5folder -s $spp_dir/fastq --recursive \
--flowcell FLO-MIN106 --kit SQK-LSK109 --compress_fastq --disable_qscore_filtering \
--device "auto" --num_callers 8 --gpu_runners_per_device 16
```

- Basic options:

```bash
# Can see which config file is called by kit/flowcell
guppy_basecaller --print_workflows

# Basic command:
guppy_basecaller -i path/fast5folder -s path/fastqfolder
--flowcell FLO-MIN106 --kit SQK-LSK109 #OR
--config dna_r9.4.1_450bps_hac.cfg # HAC: high accuray mode vs. fast mode

-r # Recursive for folders inside the -i folder

# CPU:
--num_callers # Number of files basecalling (ONT: one at a time and dedicate all processors to it)
--cpu_threads_per_caller 8 # Threads/processors available
# GPU: for our cluster Tesla V100, according to guppy manual
--device "auto"
--num_callers 8
--gpu_runners_per_device 16

--resume # If crushed and want to pick up where it left
```

## NanoPlot visualization of ONT raw reads

- Plots for quality scores, lengths etc. of fastq files.

https://github.com/wdecoster/nanoplot<br>
https://academic.oup.com/bioinformatics/article/34/15/2666/4934939

- From specific **taxa_folder** with basecalled data:

```bash
salloc --gres=gpu:1 -p gpu_test --mem 1500 -t 2:00:00
source activate genomes

mkdir -p nanoplot
cd fastq
NanoPlot --summary sequencing_summary.txt -o ../nanoplot -p $i- -f pdf --N50 -t 8 # Max plots
```

- Then download .html files in nanoplot folder to open in browser

## Remove Illumina adapters with TrimGalore

Remove adaptors and low quality reads from Illumina short reads.

http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/<br>
https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

- From **illumina** folder:

```bash
sbatch --array=1-2 ../../scripts/trimgalore.sh /wholePATHtoRAWDATA/PREFIXofFILES
# example prefix /PATH/DNA-Acutogordius_australiensis_152393_S1_L00
```

```bash
%%bash
#!/bin/bash

#SBATCH -J TrimG                                # Job name 
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 4                                    # Number of cores/cpus
#SBATCH -t 00-08:00                             # Runtime in DD-HH:MM
#SBATCH --mem 200                               # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)

module load cutadapt/1.8.1-fasrc01 parallel

# Create variable to hold the directory name (which is the species name) and move inside folder
spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../spp-names)
mkdir -p $spp_dir/trimG
cd $spp_dir/trimG

parallel -j $SLURM_NTASKS ../../../../scripts/Trim_Galore/trim_galore --phred33 --gzip --length 50 --retain_unpaired --paired \
${1}{}_R1* ${1}{}_R2* ::: 1 2 3 4
```

## Estimate genome size with Jellyfish

Esimate genome size based on short reads.

https://github.com/gmarcais/Jellyfish/tree/master/doc : Jellyfish

http://qb.cshl.edu/genomescope/ : Web service to estimate genome size from histogram of hellyfish

https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/# : Nice tutorial

- From **illumina** folder:

```bash
sbatch --array=1-2 ../../scripts/jellyfish.sh
```

```bash
%%bash
#!/bin/bash

#SBATCH -J Jellyfish                            # Job name 
#SBATCH -t 00-08:00                             # Runtime in DD-HH:MM
#SBATCH --mem 15000                             # Memory for all cores in Mbytes

module load jellyfish/2.2.5-fasrc01

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../spp-names)

mkdir -p $spp_dir/jellyfish
cd $spp_dir/jellyfish

# Count all k-mers
ls ../trimG/*val_*.fq.gz | xargs -n 1 echo zcat > generators
jellyfish count -C -s 2G -t 4 --disk -m 21 -Q "20" -o out21kmer.jf -g generators -G 4

# Make histogram
jellyfish histo -t 8 out21kmer.jf > out21kmer-${spp_dir}.histo
```

```bash
jellyfish mem -m 21 -s 1G # To see how much memory is required
jellyfish info out21kmer.jf # To get info on the run

-m kmer size
-Q (Phred score) Any base with quality below this character is changed to N
-C Count both strand, canonical representation (false) (--canonical)
-g File of commands generating fast[aq] (--generator=path)
-G Number of generators run simultaneously (1) (--Generators=uint32) #number of threads
```

## Assembly of ONT data with Flye

https://github.com/fenderglass/Flye

Flye was primarily developed to run on base-called raw reads and does not require any prior error correction. Flye automatically detects chimeric reads or reads with low quality ends, so you do not need to curate them before the assembly.

Input reads can be in FASTA or FASTQ format. You may specify multiple files with reads (separated by spaces).

Polishing is performed as the final assembly stage. By default, Flye runs one polishing iteration. Additional iterations might correct a small number of extra errors (due to improvements on how reads may align to the corrected assembly).

- In **assemblies** folder: (output folder name can be edited inside script)

```bash
sbatch --array=1 ../../scripts/flye.sh PATH_TO_FASTQ_FOLDER
```

```bash
%%bash
#!/bin/bash

#SBATCH -J flye                                 # Job name 
#SBATCH -n 32                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 2-00:00                              # Runtime in DD-HH:MM
#SBATCH --mem 184G                              # Memory for all cores in Mbytes

source activate genomes

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../spp-names)
outfolder=flye

cd $spp_dir

flye --nano-raw ../$1/*.fastq.gz --out-dir $outfolder --threads $SLURM_NTASKS
```

Main output files are:

- assembly.fasta - Final assembly. Contains contigs and possibly scaffolds (see below).
- assembly_graph.{gfa|gv} - Final repeat graph. Note that the edge sequences might be different (shorter) than contig sequences, because contigs might include multiple graph edges (see below).
- assembly_info.txt - Tab-delimited table with extra information about contigs. Columns as follows:

    - Contig/scaffold id
    - Length
    - Coverage
    - Is circular, (Y)es or (N)o
    - Is repetitive, (Y)es or (N)o
    - Multiplicity (based on coverage)
    - Alternative group
    - Graph path (graph path corresponding to this contig/scaffold).
    - Scaffold gaps are marked with ?? symbols, and * symbol denotes a terminal graph node.

Alternative contigs (representing alternative haplotypes) will have the same alt. group ID. Primary contigs are marked by *

## Polishing with long reads

### BWA mapping

Mapping of long reads to the assembly. Mapping information is used for polishing by racon and pilon.

https://github.com/lh3/bwa#type

Make an index file from the assembly, then run mapping.

- From specific assembly folder (e.g. SPP_FOLDER/flye):
```
sbatch PATH_SCRIPTS/bwa-index.sh <assembly-file>
```

```bash
%%bash
#!/bin/bash

#SBATCH -J bwaIndex                             # Job name 
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 1                                    # Number of cores/cpus
#SBATCH -t 00-00:30                             # Runtime in DD-HH:MM
#SBATCH --mem 200                               # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)

module load bwa/0.7.17-fasrc01

mkdir -p bwa-index

bwa index -p bwa-index $1
mv bwa-index.* bwa-index/
```

- From specific assembly folder (e.g. SPP_FOLDER/flye):
```
sbatch --array=8 PATH_SCRIPTS/bwa-map.sh
```

```bash
%%bash
#!/bin/bash

#SBATCH -J bwaMap                               # Job name 
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 16                                   # Number of cores/cpus
#SBATCH -t 00-08:00                             # Runtime in DD-HH:MM
#SBATCH --mem 1000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)

module load bwa/0.7.17-fasrc01

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../spp-names)
ont=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_ont.fastq.gz

mkdir -p bwa-mapping

bwa mem -t $SLURM_NTASKS -x ont2d \
bwa-index/bwa-index $ont > bwa-mapping/mapping.sam
```

### Racon

Intended as a standalone consensus module to correct raw contigs generated by rapid assembly methods which do not include a consensus step.

https://github.com/lbcb-sci/racon

Racon takes as input only three files: contigs in FASTA/FASTQ format, reads in FASTA/FASTQ format, and alignments between the reads and the contigs in MHAP/PAF/SAM format. Output is a set of polished contigs in FASTA format printed to stdout. All input files can be compressed with gzip (which will have impact on parsing time).

The step below, medaka, was trained with non-default parameters of racon, so we set the same parameters to have similar error profiles. The effect of deviating from this prescription has not been explored with recent basecallers, and it may well be the case that medaka is not overly sensitive to changes to these values ([reference]((https://nanoporetech.github.io/medaka/draft_origin.html#recommendations)).

- In specific assembly folder:

```bash
sbatch --array=1 PATH_TO_scripts/racon.sh <assembly-file>
# Less than 10 min on Acutogordius flye assembly
```

```bash
%%bash
#!/bin/bash

#SBATCH -J Racon                                # Job name 
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-03:00                              # Runtime in DD-HH:MM
#SBATCH --mem 60G                               # Memory for all cores in Mbytes

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../../spp-names)

ont=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_ont.fastq.gz

mkdir -p racon

/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/racon/build/bin/racon -t $SLURM_NTASKS \
-m 8 -x -6 -g -8 -w 500 \ # non-default parameter values used to train medaka
$ont bwa-mapping/mapping.sam $1 > racon/racon.fasta
```

### Medaka

Tool to create consensus sequences and variant calls from nanopore data. This task is performed using neural networks applied a pileup of individual sequencing reads against a draft assembly.

https://github.com/nanoporetech/medaka

```bash
# create new environment with latest medaka
conda create -n medaka -c conda-forge -c bioconda "python=3.8" "medaka=1.3.2"
```

As input medaka accepts reads in either .fasta or .fastq. It requires a draft assembly as a .fasta (e.g. directly from assembler or after racon). The model parameter should be updated as needed (see [medaka details](https://github.com/nanoporetech/medaka#models)), especially note if there is a more recent guppy version specification.

- In specific assembly folder:
```bash
sbatch --array=1 ../../../../scripts/medaka.sh racon/racon.fasta
```

```bash
%%bash
#!/bin/bash

#SBATCH -J Medaka                               # Job name 
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 100G                              # Memory for all cores in Mbytes

source activate medaka

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../../spp-names)

ont=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_ont.fastq.gz
DRAFT=$1
OUTDIR=medaka

medaka_consensus -i $ont -d ${DRAFT} -o ${OUTDIR} -t $SLURM_NTASKS -m r941_min_high_g360
```

## purge_dups

Filter to eliminate haplotigs and contig overlaps from a de novo assembly based on read depth.

https://github.com/dfguan/purge_dups

Install purge_dups from the github, and minimap2 with bioconda:
```bash
conda install -c bioconda minimap2 
```

There are 4 steps described in this [Pipeline Guide](https://github.com/dfguan/purge_dups#--pipeline-guide): "Steps with same number can be run simultaneously. Among all the steps, although step 4 is optional, we highly recommend our users to do so, because assemblers may produce overrepresented sequeences. In such a case, the final step 4 can be applied to remove those seqeuences."

Step 1: Run minimap2 to align ont data and generate paf files, then calculate read depth histogram and base-level read depth.

Step 1: Split an assembly and do a self-self alignment.

Step 2: Purge haplotigs and overlaps.

Step 3: Get purged primary and haplotig sequences from draft assembly.
<br>Notice this command will only remove haplotypic duplications at the ends of the contigs. If you also want to remove the duplications in the middle, please remove -e option at your own risk, it may delete false positive duplications.

(Step 4: Merge hap.fa and \$hap_asm and redo the above steps to get a decent haplotig set. - These instructions are not clear, $hap_asm is an alternative assembly I think, but what would that be? And how exactly to proceed? I'm considering I don't have an alternative assembly. If the BUSCO of the final purged.fa output from step 3 is good, then don't even worry about this)

- In specific assembly folder:
```bash
sbatch --array=1 ../../../../scripts/purgedups.sh medaka/consensus.fasta
```

```bash
%%bash
#!/bin/bash

#SBATCH -J PurgeDups                            # Job name 
#SBATCH -n 32                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-02:00                              # Runtime in DD-HH:MM
#SBATCH --mem 10000                             # Memory for all cores in Mbytes

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

##### coverage distribution #####
/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/purge_dups/scripts/hist_plot.py \
-c cutoffs PB.stat hist.pdf
```

purge_dups limitation:
<br>Read depth cutoffs calculation: the coverage cutoffs can be larger for a low heterozygosity species, which causes the purged assembly size smaller than expected. In such a case, please use script/hist_plot.py to make the histogram plot and set coverage cutoffs manually.

The last command of the script produces a coverage plot that can be used to validate the cutoff values.

## Polishing with short reads - HyPo

Hybrid Polisher, uses short (as well as long reads) within a single run to polish a long reads assembly of small and large genomes.

https://github.com/kensung-lab/hypo

```bash
# install from bioconda in new python env
conda create -n hypo -c conda-forge -c bioconda python=3 hypo minimap2 samtools bedtools
```

"short reads" doesn't necessarily have to be NGS short reads; HiFi genomic reads (e.g. CCS) like those generated from PacBio SequelII could also be used instead. The requirement is that those reads should be highly accurate (>98% accuracy).

Required input:

Short reads (in FASTA/FASTQ format; can be compressed)
<br>Draft contigs (in FASTA/FASTQ format; can be compressed)
<br>Alignments between short reads and the draft (in sam/bam format; should contain CIGAR)
<br>If long (noisy) reads are also to be used for polishing, then alignments between long reads and the draft (in sam/bam format; should contain CIGAR)
<br>Expected mean coverage of short reads and approximate size of the genome

Hypo (conceptually) divides a draft (uncorrected) contig into two types of regions (segments): Strong and Weak. Strong regions are those which have strong evidence (support) of their correctness and thus do not need polishing. Weak regions, on the other hand, will be polished using POA. Each weak region will be polished using either short reads or long reads; short reads taking precedence over long reads. To identify strong regions, we make use of solid kmers (expected unique genomic kmers). Strong regions also play a role in selecting the read-segments to polish their neighbouring weak regions. Furthermore, the approach takes into account that the long reads and thus the assemblies generated from them are prone to homopolymer errors.

From Ritu-Kundu in a github issue: "Hypo currently doesn't use paired-end info for short reads and treat them as individual reads. It means that HiSeq2500 and linked reads can be combined together as short reads input with their coverage set as the combined coverage. As for cleaning, it is usually a good idea to clean the Illumina reads before using them in any downstream analysis. If one has error-corrected reads available, we would recommend using them but it is not a pre-requisite as such.
There is absolutely no need to merge the paired-end reads for Hypo."

- From specific assembly folder:
```bash
sbatch --array=1 ../../../../scripts/hypo.sh 173m
```

```bash
%%bash
#!/bin/bash

#SBATCH -J HyPo                                 # Job name 
#SBATCH -n 32                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 184G                              # Memory for all cores in Mbytes

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
```

Alternative way to calculate coverage:
<br>(result different by very little compared to samtools command above, don't know exactly why)
```bash
bedtools genomecov -ibam mapped-sr.sorted.bam -bga -split > CoverageTotal.bedgraph
coverage=$(cat CoverageTotal.bedgraph | awk '{ $5 = $3 - $2; $6 = $4 * $5; size+=$5; cov+=$6} END { print cov/size }')
```

## Check contamination with BlobTools

Sujai Kumar gave a short workshop on using the BlobTools viewer for the genomics-mgig group (slack).
<br>Slides: https://tiny.cc/btk-molluscs
<br>Viewer: https://blobtoolkit.genomehubs.org
<br>Recording: https://www.youtube.com/watch?v=Lw3WaGpSEFM
<br>Then on how to build blobtools for your data - tutorial: https://github.com/blobtoolkit/tutorials/blob/main/2021-04-28_blobtoolkit_commandline_mollusc.md <br>Workshop video: https://www.youtube.com/watch?v=SuKBKEH0cMA
<br>But the snakemake pipeline he uses did not work for me, probably some broken path in the minimap commands when I give both Illumina and ONT reads. So in the end, I am following the original BlobTools2 tutorial: https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/getting-started-with-blobtools2/

### Installation BlobTools2

Installation commands are mostly the same as in Kumar's tutorial above, except that I create one environment for mamba and then use it to install the other packages in the same env:

- In home folder, install python dependencies:

```bash
module load tools/python/3.7 #module load Anaconda3/2020.11
conda create -n blob -c conda-forge python=3.8 mamba
source activate blob

mamba install -c conda-forge -c bioconda -c tolkit aria2==1.34.0 busco==5.1.2 defusedxml==0.7.1 diamond==2.0.8 docopt==0.6.2 geckodriver==0.29.0 minimap2==2.17 mosdepth==0.2.9 nodejs==14.14.0 parallel pip==21.0.1 psutil==5.8.0 pysam==0.16.0.1 pyvirtualdisplay==2.1 pyyaml==5.4.1 samtools==1.10 selenium==3.141.0 seqtk==1.3 snakemake==6.0.5 tolkein==0.2.6 tqdm==4.59.0 ujson==4.0.2 urllib3==1.26.3
```

- In holylfs scripts folder, install blobtoolkit:

```bash
VERSION=release/v2.5.0
mkdir -p blobtoolkit
cd blobtoolkit
git clone -b $VERSION https://github.com/blobtoolkit/blobtools2
git clone -b $VERSION https://github.com/blobtoolkit/specification
git clone -b $VERSION https://github.com/blobtoolkit/pipeline
git clone -b $VERSION https://github.com/blobtoolkit/viewer

#either run this before using blobtools functions, or put in .bashrc
#but I have not been doing this
export PATH=/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/blobtoolkit/blobtools2:/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/blobtoolkit/specification:/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/blobtoolkit/insdc-pipeline/scripts:$PATH
#export PATH=~/scripts/blobtoolkit/blobtools2:~/scripts/blobtoolkit/specification:~/scripts/blobtoolkit/insdc-pipeline/scripts:$PATH
```

- In blobtools folder, install databases:


1- Download the NCBI taxdump:

```bash
mkdir blob_databases
cd blob_databases

TAXDUMP=taxdump_2021_07
mkdir -p $TAXDUMP;
cd $TAXDUMP;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd -
```

2- Download and extract UniProt reference proteomes:

```bash
UNIPROT=uniprot_2021_07
mkdir -p $UNIPROT

#This actually takes a long time, better to put in a wrap or script with 3 days:
sbatch -t 4-00 -p shared,giribet --mem 10000 -o wget_uniprot.out -n 2 --wrap "wget -O reference_proteomes.tar.gz \
 ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
     -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
     awk '/tar.gz/ {print $9}')"

cd $UNIPROT
tar xf reference_proteomes.tar.gz

touch reference_proteomes.fasta.gz
find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz

printf "accession\taccession.version\ttaxid\tgi\n" > reference_proteomes.taxid_map
zcat */*/*.idmapping.gz | grep "NCBI_TaxID" | awk '{print $1 "\t" $1 "\t" $3 "\t" 0}' >> reference_proteomes.taxid_map

diamond makedb -p 16 --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map --taxonnodes $TAXDUMP/nodes.dmp --taxonnames $TAXDUMP/names.dmp -d reference_proteomes.dmnd
cd -
```

3- Download NCBI nt database:

```bash
NT=nt_2021_07
wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz' -P $NT/
cd $NT
for file in *.tar.gz; do
    tar xf $file && rm $file;
done
```

4- Download BUSCO data and lineages to allow BUSCO to run in offline mode:

```bash
BUSCO=busco_2022_02
mkdir -p $BUSCO
#cd $BUSCO
#wget -r https://busco-data.ezlab.org/v5/data
#find busco-data.ezlab.org -name "*.tar.gz" | parallel "cd {//}; tar -xzf {/}"
wget -q -O eukaryota_odb10.gz "https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2020-09-10.tar.gz" \
        && tar xf eukaryota_odb10.gz -C busco_2022_02
        
wget -q -O metazoa_odb10.gz "https://busco-data.ezlab.org/v5/data/lineages/metazoa_odb10.2021-02-24.tar.gz" \
        && tar xf metazoa_odb10.gz -C busco_2022_02
        
wget -q -O arthropoda_odb10.gz "https://busco-data.ezlab.org/v5/data/lineages/arthropoda_odb10.2020-09-10.tar.gz" \
        && tar xf arthropoda_odb10.gz -C busco_2022_02
        
wget -q -O nematoda_odb10.gz "https://busco-data.ezlab.org/v5/data/lineages/nematoda_odb10.2020-08-05.tar.gz" \
        && tar xf nematoda_odb10.gz -C busco_2022_02
```

### Config file

(This was required in the snakemake tutorial, but not using it in my current scripts)

- Create yaml config file. It has 6 sessions: `assembly`, `busco`, `reads`, `settings`, `similarity`, `taxon` and `keep_intermediates`.

To see what are the available [busco](https://busco.ezlab.org/busco_userguide.html) datasets:
```bash
busco --list-datasets
```

### Build up BlobDir

Before running blobtools, need to map reads to assembly, get blast hits and diamond hits. The following scripts can be run in parallel, and they create the necessary input files to run blobtools right after, adding all these data to the BlobDir.

From specific assembly folder, e.g. Nectonema flye folder:
```bash
sbatch --array=1 ../../../../scripts/blob-blast.sh
sbatch --array=1 ../../../../scripts/blob-diamond.sh
sbatch --array=1 ../../../../scripts/blob-mapILL.sh
sbatch --array=1 ../../../../scripts/blob-mapONT.sh
sbatch --array=1 ../../../../scripts/blob-busco.sh nematoda_odb10
#metazoa_odb10 arthropoda_odb10 eukaryota_odb10
```

```bash
%%bash
#!/bin/bash

#SBATCH -J blob-blast                           # Job name 
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 2-00:00                              # Runtime in DD-HH:MM
#SBATCH --mem 100G                              # Memory for all cores in Mbytes

module load Anaconda3/2020.11
source activate blob

mkdir -p blobtools

blastn -db /n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/blobtoolkit/blob_databases/nt_2021_07/nt \
       -num_threads $SLURM_NTASKS \
       -query hypo/whole_genome.h.fa \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -out blobtools/blast.out
```

```bash
%%bash
#!/bin/bash

#SBATCH -J blob-diamond                         # Job name 
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 15G                               # Memory for all cores in Mbytes

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
```

```bash
%%bash
#!/bin/bash

#SBATCH -J blob-mapILL                          # Job name 
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 25G                               # Memory for all cores in Mbytes

module load Anaconda3/2020.11
source activate blob

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../../spp-names)

mkdir -p blobtools

illumina1=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_illumina1.fastq.gz
illumina2=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_illumina2.fastq.gz

minimap2 -ax sr \
         -t $SLURM_NTASKS hypo/whole_genome.h.fa \
         $illumina1 $illumina2 \
         | samtools sort -@$SLURM_NTASKS -O BAM -o blobtools/illumina.reads.bam -
```

```bash
%%bash
#!/bin/bash

#SBATCH -J blob-mapONT                          # Job name 
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 25G                               # Memory for all cores in Mbytes

module load Anaconda3/2020.11
source activate blob

spp_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../../spp-names)

mkdir -p blobtools

ont=$SCRATCH/giribet_lab/Users/tauanajc/${spp_dir}_ont.fastq.gz

minimap2 -ax map-ont \
         -t $SLURM_NTASKS hypo/whole_genome.h.fa \
         $ont > blobtools/ont.reads.sam
cat blobtools/ont.reads.sam | samtools sort -O BAM -o blobtools/ont.reads.bam -
#command fails if output of minimap is piped to sort (maybe too long reads causing truncated files)
#but works just fine if saves sam and then cat pipe
```

```bash
%%bash
#!/bin/bash

#SBATCH -J blob-busco                           # Job name 
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 10G                               # Memory for all cores in Mbytes

module load Anaconda3/2020.11
source activate blob

mkdir -p blobtools
cd blobtools

busco \
    -i ../hypo/whole_genome.h.fa \
    -o assembly_$1 \
    -l /n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/blobtoolkit/blob_databases/busco_2022_02/$1 \
    -m geno \
    -c $SLURM_NTASKS

mv ../blob-busco.$SLURM_ARRAY_TASK_ID.$SLURM_JOB_ID.* ./
```

### Run BlobTools

From specific assembly folder, e.g. Nectonema flye folder:
```bash
sbatch --array=2 ../../../../scripts/blobtools.sh taxid
```

--taxrule default is bestsumorder, option is best sum. This relates to how the taxonomy of blast hits is applied to contigs, but I don't fully understand the difference, so keep an eye, and if weird, try other option: https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/adding-data-to-a-dataset/adding-hits/


```bash
%%bash
#!/bin/bash

#SBATCH -J blobtools                            # Job name 
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 10000                             # Memory for all cores in Mbytes

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
    --busco assembly_metazoa_odb10/run_metazoa_odb10/full_table.tsv \
    --busco assembly_eukaryota_odb10/run_eukaryota_odb10/full_table.tsv \
    --busco assembly_nematoda_odb10/run_nematoda_odb10/full_table.tsv \
    --busco assembly_arthropoda_odb10/run_arthropoda_odb10/full_table.tsv \
    blobtools/BlobDir
#--meta blobtools/config.yaml \
```

### View blob

https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/opening-a-dataset-in-the-viewer/

View step was giving error in my installation of blobtools2, something about npm. But docker version worked without any issues!

From **blobtools** folder inside a specific assembly:
```bash
salloc -n 16 -p test -t 0-04:00 --mem=10000
singularity exec --cleanenv ../../../../../scripts/blobtoolkit_2.6.5.sif blobtools view --remote BlobDir
```

Then open in terminal in local machine:
```bash
ssh -N -L 8001:localhost:8001 -L 8000:localhost:8000 -J tauanajc@odyssey.fas.harvard.edu holy7c18312
```

Then in browser:
http://localhost:8001/view/BlobDir/dataset/BlobDir/blob

### bioawk

Use `bioawk` to get sequences of specific contigs and blast them, and also get length and other metrics directly from fasta assembly.

https://bioinformaticsworkbook.org/Appendix/Unix/bioawk-basics.html#gsc.tab=0

```bash
source activate blob
conda install -c bioconda bioawk
```

```bash
# convert fasta to table of contig name and sequence, then grep the contig of interest
bioawk -t -c fastx '{ print $name, $seq }' input.fasta | grep "contig_2422_1"
```

### Nectonema Blobs

Highlights from blobtools view (illumina coverage x GC). Images and csv saved in a folder. Need to develop a system to organize blobtools results and which contigs to remove.

Arthropod hits

- TOO LOW GC:<br>
contig_2422_1 & contig_447_1: blasts to nothing (1-4 bad hits of ~20 bp)

- TOO HIGH COVERAGE:<br>
contig_623_1: many blasts but bad score (~75) on tiny pieces
<br>contig_853_1: can see big chunk of repetitive sequence in this contig, so understandable to have high coverage; blasts to 
<br>contig_874_1:
<br>contig_912_1:
<br>contig_968_1:

No hits

- TOO LOW GC:<br>
contig_3243_1
contig_481_1
contig_768_1

- TOO LOW COVERAGE:<br>
contig_3155_1 (lowest coverage on both illumina and ont)

- TOO HIGH GC:
contig_1011_1
contig_2201_1
contig_282_1
contig_517_1

- TOO HIGH COVERAGE:<br>
contig_1351_1
contig_1768_1
contig_3107_1
contig_3108_1
contig_3443_1
contig_478_1

- TOO HIGH COVERAGE AND GC:<br>
contig_547_1

Annelida (/Cnidaria/Chordata) hits

contig_3341_1: higher GC; 

Mucoromycota hits

- TOO HIGH GC:<br>
contig_1450_1
contig_935_1

- regular blob position:
contig_1093_1
contig_1516_1

Other hits

- TOO HIGH COVERAGE:<br>
contig_1844_1: hit to Nematomorpha! Best hits are =<9% of contig, high scores, all Nectonema, Paragordius, Gordius 18S rRNA, so total sense to be high coverage.

- TOO HIGH GC:<br>
contig_3408_1: hit to Negative-strand RNA virus (Negarnaviricota)! RNA virus, so wouldn't have gotten it unless it was incorporated in the genome of the Nectonema? Apparently virus in this group (Peropuvirus) infect parasitoid wasps, barnacles, pillworms, woodlice odonates, or copepods! https://talk.ictvonline.org/ictv-reports/ictv_online_report/negative-sense-rna-viruses/w/artoviridae/1127/genus-peropuvirus

There are also some other hits with regular blob position, but best hits to Proteobacteria, Oomycota, Streptophyta, and Apicompexa. These worry me more, shoud we leave or remove? A couple of them are "larger", 150kb, 70kb, while others are small 4kb

### Acutogordius Blobs

Highlights from blobtools view (illumina coverage x GC). Images and csv saved in a folder.

Arthropod hits

- SLIGHTLY HIGHER GC, but perfect coverage:<br>
contig_2510_1, contig_712_1 (longer chunk to ant), contig_7651_1 (longer chunk to scropion), contig_1232_1, contig_3365_1: is it possible that these hits are to hosts? I'm probably keeping the contigs, as they are just a bit higher GC, and perfect fit to coverage.

No hits - all very off in metrics, small contigs, will remove all 6

- TOO LOW GC:<br>
contig_2528_1

- TOO LOW COVERAGE:<br>
contig_5382_1<br>
contig_6030_1<br>
contig_1382_1<br>
contig_6629_1<br>
contig_8573_1

- SLIGHTLY LOWER COVERAGE AND GC:<br>
contig_3668_1, contig_4869_1: 10-20Kb - will probably remove both

Nematoda

- SLIGHTLY HIGH GC:<br>
contig_5557_1, contig_4559_1 (larger chunk hit): will probably keep, as perfect coverage

Nematomorpha

- TOO HIGH COVERAGE AND TOO HIGH GC:<br>
contig_2666_1: hit to Nematomorpha! Best hits are =10% of contig, 95% identity, all to different genera of Nematomorpha 18S rRNA, so total sense to be high coverage. Will keep.

Other

There are a bunch of contigs that have perfect coverage and gc content, but best hits to bacteria, protists or fungi (also other invert groups, but I'm considering these fine). The hit graph shows very small hit in the contigs, so might be nothing. The ones with slightly longer hits (blob graphs, I did not blast yet) are:<br>
contig_3913_1 (200Kb): Candidatus Entotheonella bacteria, associated to marine sponges https://www.pnas.org/content/114/3/E347<br>
contig_235_1 (75Kb): Eukaryota-undef / Guillardia theta, protist algae https://en.wikipedia.org/wiki/Guillardia<br>
contig_2712_1, contig_4538_1, contig_3095_1 (15-30Kb): Smittium, from an order of fungi found in the guts of insect larvae (most often aquatic flies) https://en.wikipedia.org/wiki/Smittium<br>
contig_7414_1 (20Kb): Dasytricha ruminantium, apicomplexa well studies from ruminants https://lmgtfy.app/?q=Dasytricha+ruminantium<br>
contig_7_1 (15Kb): Winogradskyella sp. PG-2, small contig, bacteria from surface seawater https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4038882/

Mixed Hits

- TOO LOW COVERAGE:<br>
These are not so low as others described above, but still illumina coverage is ~5-22Kb, compared to >=100Kb of most other contigs. I'm already keeping others around 40Kb because Acutogordius has high heterozygosity (estimated 3.72% vs. 0.71% Nectonema), so 50% coverage is expected in haploid segments. But these might be too low, even if some are relatively large contigs of 185Kb and have normal cg content. Some hit to Arthropoda, some Nematoda, some Chordata, and the rest are NoHits. Maybe they are fine, but to be on the safe side, I will remove them at the moment. They are in file Blob-5-22-illumina-coverage.csv. Together with other contigs listed above, the coverage filter ends up being removing all contigs with illumina coverage < 25.

### Check Blasts

For some, might want to explore the hits in more detail. From assembly folder, use contig name to get sequence and paste into blast browser:

```bash
bioawk -t -c fastx '{ print $name, $seq }' hypo/whole_genome.h.fa | grep "contig_2712_1"
```

### Remove contaminants

After inspecting blob plots, decide on which to remove from assembly as possible contaminands, and put contig names in a **text file inside blobtools** folder to remove with `bioawk`.

```bash
bioawk -cfastx 'BEGIN{while((getline k <"contigs-to-exclude.txt")>0)i[k]=1}{if(!(i[$name]))print ">"$name"\n"$seq}' ../hypo/whole_genome.h.fa > decontaminated.fasta

# check the correct number of contigs was removed:
grep -c ">" ../hypo/whole_genome.h.fa
grep -c ">" decontaminated.fasta
```

## HiC

Arima protocol for mapping HiC reads to an assembly: https://github.com/ArimaGenomics/mapping_pipeline/blob/master/Arima_Mapping_UserGuide_A160156_v02.pdf
<br>Then use SALSA2 to create maps.

(Arima's manual to use Juicer: https://arimagenomics.com/wp-content/files/Bioinformatics-User-Guide-Arima-HiC-and-Arima-High-Coverage-HiC.pdf - not used in this pipeline, but saving link for reference. Points to a paper and an online resource listing all HiC-related tools)

Download Arima scripts in **scripts** folder:
```bash
git clone https://github.com/ArimaGenomics/mapping_pipeline.git
mv mapping_pipeline arima_mapping
# below I install picard in conda, but did not work, because arima's command uses the java version, so:
wget https://github.com/broadinstitute/picard/releases/download/2.25.7/picard.jar
```

Create conda environment with required programs: bwa, samtools, picard:
```bash
module load Anaconda3
conda create -n hic -c conda-forge -c bioconda python=3.8 samtools==1.10 picard==2.26.10 bwa=0.7.17
```

There are 4 parts to the mapping pipeline: indexing the assembly, mapping HiC reads to assembly, combining replicate reads from different lanes, and removing duplicates.

### Mapping - Indexing

Create bwa index of the assembly.

- From specific assembly folder:
```
sbatch ../../../../scripts/hic-indexing.sh <post-blobtools-assembly-file>
```

```bash
%%bash
#!/bin/bash

#SBATCH -J HiC-indexing                         # Job name 
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 1                                    # Number of cores/cpus
#SBATCH -t 00-00:30                             # Runtime in DD-HH:MM
#SBATCH --mem 1000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)

module load Anaconda3/2020.11
source activate hic

mkdir -p hic/mapping

# Create assembly.fai index
samtools faidx $1
rsync -a blobtools/decontaminated.fasta.fai hic/mapping/

# Create bwa index
bwa index -a bwtsw -p hic-bwa-index $1
mv hic-bwa-index.* hic/mapping/
```

### Mapping

From Arima manual:<br>
STEP 1: Use BWA-MEM to align the Hi-C paired-end reads to the reference sequences. Because Hi-C captures conformation via proximity-ligated fragments, paired-end reads are first mapped independently (assingle-ends) usingBWA-MEM and are subsequently paired in a later step.<br>
STEP 2: Subsequent to mapping as single-ends, some of these single-end mapped reads can manifest a ligation junction and are therefore considered 'chimeric' (i.e. they do not originate from a contiguous piece of DNA). When BWA-MEM maps these chimeric reads, there can be high quality mapping on both the 5’- side and 3’-side of the ligation junction within a given read. In such cases, only the 5’-side should be retained because the 3’-side can originate from the same contiguous DNA as the 5’-side of the reads mate-pair. Therefore, we retain only the portion of the chimeric read that maps in the 5’-orientation in relation to its read orientation. This is accomplished using the script `filter_five_end.pl`.<br>
STEP 3: After filtering, we pair the filtered single-end Hi-C reads using `two_read_bam_combiner.pl`, which outputs a sorted, mapping quality filtered, paired-end BAM file.
STEP 4: We then add read groups to this BAM file using Picard Tools.

- From specific assembly folder:

```bash
# get list of HiC fastq file names and save in file inside hic folder - CHANGE NAME OF TAXA
\ls /n/Giribet_Lab/tauanajc/GASTROPODA/rawdata/2020_03-NovaSeqS4/Tauana/HiC-Acutogordius_australiensis_152393_S6_L00* | rev | cut -d/ -f 1 | cut -d_ -f3- | rev | uniq > hic/fastq-names-hic.txt

# run
sbatch --array=1-2 ../../../../scripts/hic-map.sh # array is number of lines in fastq-names-hic.txt
```

```bash
%%bash
#!/bin/bash

#SBATCH -J HiC-map                              # Job name 
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 16                                   # Number of cores/cpus
#SBATCH -t 00-08:00                             # Runtime in DD-HH:MM
#SBATCH --mem 5000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)

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
java -Xmx4G -Djava.io.tmpdir=temp/ -jar picard AddOrReplaceReadGroups INPUT=tmp-bam/$NAME.bam OUTPUT=paired-bam/$NAME.bam ID=$NAME LB=$NAME SM=$LABEL PL=ILLUMINA PU=none
```

### Combine lanes and remove duplicates

From Arima's manual:<br>
We also use Picard Tools to discard any PCR duplicates present in the paired-end BAM file generated above. If applicable, we require that you merge paired-end BAM files that were sequenced via multiple Illumina lanes from the same library (i.e. technical replicates) before removing PCR duplicates.

From specific assembly folder:
```bash
sbatch ../../../../scripts/hic-dedup.sh
```

```bash
%%bash
#!/bin/bash

#SBATCH -J HiC-dedup                            # Job name 
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 16                                   # Number of cores/cpus
#SBATCH -t 00-08:00                             # Runtime in DD-HH:MM
#SBATCH --mem 15000                             # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)

module load Java/1.8
module load Anaconda3/2020.11
source activate hic

LABEL=$(head -n1 hic/fastq-names-hic.txt | rev | cut -d_ -f3- | rev)
STATS='/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/arima_mapping/get_stats.pl'
PICARD='/n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/picard.jar'
# This does not expand variable because of single quotes:
#INPUTS_TECH_REPS=("INPUT=${LABEL}_S6_L001.bam" "INPUT=${LABEL}_S6_L002.bam" "INPUT=${LABEL}_S6_L003.bam" "INPUT=${LABEL}_S6_L004.bam")
# This does not read all elements of array in the java command, just the first…
#INPUTS_TECH_REPS=('INPUT='"${LABEL}"'_S6_L002.bam' 'INPUT='"${LABEL}"'_S6_L003.bam' 'INPUT='"${LABEL}"'_S6_L004.bam')
# So decided to list all elements directly in command, instead of using bash array

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
```

From Arima's manual:<br>
The pipeline is complete. The final output of this pipeline is a single BAM file that contains the paired, 5’-filtered, and duplicate-removed Hi-C reads mapped to the reference sequences of choice. The resulting statistics file has a breakdown of the total number of intra-contig read-pairs, long-range intra-contig read-pairs, and inter-contig read-pairs in the final processed BAM file.

The outputs from the Arima Hi-C Mapping Pipeline are the inputs for salsa2.

Salsa2 requires a bed file which is converted from bam using Bedtools, an assembly in fasta format, an index for that assembly in .fai (we already made this at the beginning of the mapping pipeline)

```
# convert bam to bed file and sort by read name
bamToBed -i hic/mapping/final-deduplicated-bam/$LABEL.final.bam > hic/$LABEL.bed
sort -k 4 hic/$LABEL.bed > tmp && mv tmp hic/$LABEL.bed
```

## YaHS

https://github.com/c-zhou/yahs
<br>https://github.com/sanger-tol/yahs

Input file is the same bed file used for salsa.

YaHS runs pretty fast on nematomorph genomes. After that, juicer tools to create the file used for visualization of the contact map. I followed instructions on c-zhou's github for the hole thing (yahs and juicer). One of the commands is to create the file for scaffold sizes (should contain two columns - scaffold name and scaffold size): the github says to create from first two columns of the index file from the yahs scaffolds fasta, but this index is not needed otherwise, so it doesn't exist at this point. In one of the github issues, however, YaHS author Chenxy says to create this file by grepping the juicer log file - it is the preferred method anyway (deals with >2G chromosomes, https://github.com/c-zhou/yahs/issues/4#issuecomment-1159017424)

From specific assembly folder:
```bash
sbatch ../../../../scripts/hic-yahs.sh
# last juicer step is long and memory intensive (12h, >100G for Acutogordius)
```


```bash
%%bash
#!/bin/bash

#SBATCH -J yahs                                 # Job name 
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 8                                    # Number of cores/cpus
#SBATCH -t 07-00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 184000                            # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)

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
#awk '{print $1,$2}' ../mapping/yahs.out_scaffolds_final.fa.fai > scaffolds_final.chrom.sizes # but have to create index first!
cat juicer1.log | grep "PRE_C_SIZE:" | cut -d' ' -f 2,3 > scaffolds_final.chrom.sizes
#asmlen=$(cat juicer1.log | grep "PRE_C_SIZE:" | cut -d' ' -f 2,3) # same content but in a variable instead of file

# juicer pre to create .hic file
(java -jar -Xmx184G $JUICER pre --threads 8 alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes) 2>&1 | tee juicer2.log && (mv out.hic.part out.hic)
#(java -jar -Xmx32G $JUICER pre --threads 8 alignments_sorted.txt out.hic.part <(echo "assembly ${asmlen}")) && (mv out.hic.part out.hic)
```

### Visualize with Juicebox

https://github.com/aidenlab/Juicebox/wiki


Download out.hic file and open with Juicebox (online: https://aidenlab.org/juicebox/ or java or Desktop: https://github.com/aidenlab/Juicebox/wiki/Download). The online version allows to share a link to put in the paper, so that people can load and play with your contact map.

![Acutogordius_HiC_YaHS_nogrid.png](attachment:a2b31c8f-50fb-4d24-850a-8f01ca759eb3.png)

## BUSCO: check assembly quality

Database of orthologs from different large sets of organisms. Based on evolutionarily-informed expectations of gene content of near-universal single-copy orthologs.

https://busco.ezlab.org<br>
https://gitlab.com/ezlab/busco

Interpretation of results: https://busco.ezlab.org/busco_userguide.html#interpreting-the-results

```bash
# new environment for busco dependencies
conda create -n busco545 -c conda-forge -c bioconda busco=5.4.5
# To see which datasets are available
busco --list-datasets
```

[I was getting error "ERROR:busco.BuscoRunner Metaeuk did not recognize any genes matching the dataset" with HiC assembly only. Apparently issue with memory, even though log files do not indicate that. Based on seeing other people mention memory issues and trying to accomodate with metaeuk flag, I simply gave a lot of memory to the job (184GB), and it ran correctly.]

Of interest: eukaryota_odb10, metazoa_odb10, nematoda_odb10, arthropoda_odb10

- From **assemblies/busco/** folder:
```bash
sbatch busco-genomes.sh <assembly.file/folder> <lineage>
#Example:
sbatch busco-genomes.sh assemblies eukaryota_odb10
#for i in assemblies/all_assemblies_fasta/*; do sbatch busco-metazoa.sh $i; done
```

```bash
%%bash
#!/bin/bash

#SBATCH -J Busco                                # Job name 
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 184000                            # Memory for all cores in Mbytes

module load Anaconda3
source activate busco

#busco -i $1 -l metazoa_odb10 -o busco_metazoa -m genome --cpu $SLURM_NTASKS
busco -i $1 -l $2 -o out_$2 -m genome --cpu $SLURM_NTASKS
```

Important output files

Inside folder structure:<br>
busco/out_metazoa_odb10/run_metazoa_odb10/

    - short_summary.txt
    - full_table.tsv
    - missing_busco_list.tsv

## Assembly stats with QUAST

QUality ASsessment Tool evaluates genome/metagenome assemblies by computing various metrics.

https://github.com/ablab/quast

Works both with and without reference genomes. However, it is much more informative if at least a close reference genome is provided along with the assemblies. The tool accepts multiple assemblies, thus is suitable for comparison.

```bash
# install quast
conda create -n quast -c bioconda -c conda-forge quast
```

- From specific assembly folder:
```bash
sbatch ../../scripts/quast.sh <assembly.file> [multiple files possible]
```

```bash
%%bash
#!/bin/bash

#SBATCH -J Quast                                # Job name 
#SBATCH -n 16                                   # Number of cores
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00                              # Runtime in DD-HH:MM
#SBATCH --mem 10000                             # Memory for all cores in Mbytes

module load Anaconda3
source activate quast

quast.py -t $SLURM_NTASKS -o quast "$@"
```
