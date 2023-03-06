This repository contains the scripts used to assemble and analyze the genomes of two nematomorph species. The assemblies were built with ONT long reads, and Illumina whole-genome and Hi-C short reads. If you use this repository, please cite:

Cunha TJ, De Medeiros BAS, Lord A, Sørensen MV, Giribet G. 2023. **Rampant loss of universal metazoan genes revealed by a chromosome-level genome assembly of the parasitic Nematomorpha**.

Below is the pipeline used for genome assembly, with reference links to software and information. Additional folders in this repository contain the final assemblies, and the R script for analyses of functional enrichment.


# Genome Assembly

## Setup

### Folder structure

- Inside the project folder, I keep a simple text file (spp\_names) with sample names in one column, which are used in the scripts to call input files and name new output files/folders:

```
# spp_names file example
Acutogordius_australiensis-MCZ152393
Nectonema_munidae-MCZ153622
```

Example folder organization:

>- spp_names
>- fast5_symlinks/
>- basecalled-ont/
>- illumina/
>- assemblies/

### Python environment

Create a conda environment for genome-related tools. Activate the enviroment in the scripts or before running programs in the command line.
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


## Basecall ONT with Guppy

Basecall ONT long reads.

- From folder **basecalled-ont**:

```bash
sbatch --array=1 ../../scripts/basecall.sh ../fast5-symlinks/SPECIES-FOLDER
```

This structure allows for multiple flow cells with similar but not identical names to be basecalled into one species folder.


## Visualization of ONT reads with NanoPlot

Plots for quality scores, lengths etc. of fastq files.

https://github.com/wdecoster/nanoplot<br>
https://academic.oup.com/bioinformatics/article/34/15/2666/4934939

- From specific **taxon_folder** within **basecalled-ont**:

```bash
salloc --gres=gpu:1 -p gpu_test --mem 1500 -t 2:00:00
source activate genomes

mkdir -p nanoplot
cd fastq
NanoPlot --summary sequencing_summary.txt -o ../nanoplot -p $i- -f pdf --N50 -t 8 # Max plots
```

Then download .html files in nanoplot folder to open in browser


## Remove Illumina adapters with TrimGalore

Remove adaptors and low quality reads from Illumina short reads.

http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/<br>
https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

- From **illumina** folder:

```bash
sbatch --array=1 ../../scripts/trimgalore.sh /wholePATHtoRAWDATA/PREFIXofFILES
# example prefix /PATH/Acutogordius_australiensis_MCZ152393_S1_L00
```


## Estimate genome size with Jellyfish

Esimate genome size based on short reads.

https://github.com/gmarcais/Jellyfish/tree/master/doc : Jellyfish<br>
http://qb.cshl.edu/genomescope/ : Web service to estimate genome size from histogram of jellyfish<br>
https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/# : Nice tutorial

- From **illumina** folder:

```bash
sbatch --array=1 ../../scripts/jellyfish.sh
```


## Genome assembly of ONT data with Flye

https://github.com/fenderglass/Flye

Flye was primarily developed to run on base-called raw reads and does not require any prior error correction. Flye automatically detects chimeric reads or reads with low quality ends, so you do not need to curate them before the assembly.

Input reads can be in FASTA or FASTQ format. You may specify multiple files (separated by spaces).

Polishing is performed as the final assembly stage. By default, Flye runs one polishing iteration.

- In **assemblies** folder: (output folder name can be edited inside script)

```bash
sbatch --array=1 ../../scripts/flye.sh PATH_TO_FASTQ_FOLDER
```

Main output files are:

- assembly.fasta - Final assembly. Contains contigs and possibly scaffolds (see below).
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


## Polish with long reads

### BWA mapping

Mapping of long reads to the assembly. Mapping information is used for polishing by racon and pilon.

https://github.com/lh3/bwa#type

Make an index file from the assembly, then run mapping.

- From specific assembly folder (**assemblies/spp_folder/flye**):
```
sbatch PATH_SCRIPTS/bwa-index.sh <assembly-file>
```

Then:
```
sbatch --array=1 PATH_SCRIPTS/bwa-map.sh
```


### Racon

https://github.com/lbcb-sci/racon

Intended as a standalone consensus module to correct raw contigs generated by rapid assembly methods which do not include a consensus step.<br>
Racon takes as input only three files: contigs in FASTA/FASTQ format, reads in FASTA/FASTQ format, and alignments between the reads and the contigs in MHAP/PAF/SAM format. Output is a set of polished contigs in FASTA format printed to stdout. All input files can be compressed with gzip (which will have impact on parsing time).

The step below, medaka, was trained with non-default parameters of racon, so we set the same parameters to have similar error profiles. The effect of deviating from this prescription has not been explored with recent basecallers, and it may well be the case that medaka is not overly sensitive to changes to these values ([reference](https://nanoporetech.github.io/medaka/draft_origin.html#recommendations)).

- In specific assembly folder (**assemblies/spp_folder/flye**):

```bash
sbatch --array=1 PATH_TO_scripts/racon.sh <assembly-file>
# Less than 10 min on Acutogordius flye assembly
```


### Medaka

https://github.com/nanoporetech/medaka

Tool to create consensus sequences and variant calls from nanopore data. This task is performed using neural networks applied a pileup of individual sequencing reads against a draft assembly.

```bash
# create new environment with latest medaka
conda create -n medaka -c conda-forge -c bioconda "python=3.8" "medaka=1.3.2"
```

As input medaka accepts reads in either .fasta or .fastq. It requires a draft assembly as a .fasta (e.g. directly from assembler or after racon). The model parameter should be updated as needed (see [medaka details](https://github.com/nanoporetech/medaka#models)), especially note if there is a more recent guppy version specification.

- In specific assembly folder (**assemblies/spp_folder/flye**):
```bash
sbatch --array=1 ../../../../scripts/medaka.sh racon/racon.fasta
```


## Remove haplotigs with purge_dups

Filter to eliminate haplotigs and contig overlaps from a de novo assembly based on read depth.

https://github.com/dfguan/purge_dups

Install purge_dups from the github, and minimap2 with bioconda:
```bash
source activate genomes
conda install -c bioconda minimap2 
```

[Pipeline Guide](https://github.com/dfguan/purge_dups#--pipeline-guide):

Step 1: Run minimap2 to align ont data and generate paf files, then calculate read depth histogram and base-level read depth.

Step 1: Split an assembly and do a self-self alignment.

Step 2: Purge haplotigs and overlaps.

Step 3: Get purged primary and haplotig sequences from draft assembly.
<br>Notice this command will only remove haplotypic duplications at the ends of the contigs. If you also want to remove the duplications in the middle, please remove -e option at your own risk, it may delete false positive duplications.

purge_dups limitation:
<br>Read depth cutoffs calculation: the coverage cutoffs can be larger for a low heterozygosity species, which causes the purged assembly size smaller than expected. In such a case, please use script/hist_plot.py to make the histogram plot and set coverage cutoffs manually.

The last command of the script produces a coverage plot that can be used to validate the cutoff values.

- In specific assembly folder (**assemblies/spp_folder/flye**):
```bash
sbatch --array=1 ../../../../scripts/purgedups.sh medaka/consensus.fasta
```


## Polish with short reads - HyPo

Hybrid Polisher, uses short reads within a single run to polish a long reads assembly of small and large genomes. The requirement is that those reads should be highly accurate (>98% accuracy).

https://github.com/kensung-lab/hypo

```bash
# install from bioconda in new python env
conda create -n hypo -c conda-forge -c bioconda python=3 hypo minimap2 samtools bedtools
```

Required input:

Short reads (in FASTA/FASTQ format; can be compressed)
<br>Draft contigs (in FASTA/FASTQ format; can be compressed)
<br>Alignments between short reads and the draft (in sam/bam format; should contain CIGAR)
<br>If long (noisy) reads are also to be used for polishing, then alignments between long reads and the draft (in sam/bam format; should contain CIGAR)
<br>Expected mean coverage of short reads and approximate size of the genome

From Ritu-Kundu in a github issue: "Hypo currently doesn't use paired-end info for short reads and treat them as individual reads. It means that HiSeq2500 and linked reads can be combined together as short reads input with their coverage set as the combined coverage. As for cleaning, it is usually a good idea to clean the Illumina reads before using them in any downstream analysis. If one has error-corrected reads available, we would recommend using them but it is not a pre-requisite as such.
There is absolutely no need to merge the paired-end reads for Hypo."

- From specific assembly folder (**assemblies/spp_folder/flye**):
```bash
sbatch --array=1 ../../../../scripts/hypo.sh 173m # estimated genome size from GenomeScope
```


## Check contamination with BlobTools

BlobTools2 tutorial: https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/getting-started-with-blobtools2/

### Installation BlobTools2

I found it easier/faster to create one environment for mamba and then use it to install the other packages in the same env:

- In home folder, install python dependencies:

```bash
module load Anaconda3/2020.11
conda create -n blob -c conda-forge python=3.8 mamba
source activate blob

mamba install -c conda-forge -c bioconda -c tolkit aria2==1.34.0 busco==5.1.2 defusedxml==0.7.1 diamond==2.0.8 docopt==0.6.2 geckodriver==0.29.0 minimap2==2.17 mosdepth==0.2.9 nodejs==14.14.0 parallel pip==21.0.1 psutil==5.8.0 pysam==0.16.0.1 pyvirtualdisplay==2.1 pyyaml==5.4.1 samtools==1.10 selenium==3.141.0 seqtk==1.3 snakemake==6.0.5 tolkein==0.2.6 tqdm==4.59.0 ujson==4.0.2 urllib3==1.26.3
```

- In my usual scripts folder, install blobtoolkit:

```bash
VERSION=release/v2.5.0
mkdir -p blobtoolkit
cd blobtoolkit
git clone -b $VERSION https://github.com/blobtoolkit/blobtools2
git clone -b $VERSION https://github.com/blobtoolkit/specification
git clone -b $VERSION https://github.com/blobtoolkit/pipeline
git clone -b $VERSION https://github.com/blobtoolkit/viewer
```

- In scripts/blobtoolkit folder, install databases:

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

### Installation with Docker

Both the full instalation above and the docker version work when running the regular blobtool commands. I only used the docker for the View Blob step below. In any case, here are the instructions of how to install and use an example command.

On Harvard Cannon cluster, get container with:
```bash
singularity pull --disable-cache docker://genomehubs/blobtoolkit:2.6.5
# Test simple commands:
singularity exec --cleanenv blobtoolkit_2.6.5.sif blobtools -h
singularity exec --cleanenv blobtoolkit_2.6.5.sif blobtools -v
```

Example of how to run the docker version:
```bash
# This runs fine! meta.json has all info in config file, but is not using paths. When I try the same, but without the --taxdump argument, it fails because doesn't find taxdump…
singularity exec --cleanenv blobtoolkit_2.6.5.sif blobtools create \
    --threads $THREADS \
    --fasta $DATA_DIR/Nectonema_munidae-DNA05827.fasta \
    --meta $DATA_DIR/config.yaml \
    --taxid 190569 \
    --taxdump /n/holylfs04/LABS/giribet_lab/Lab/tauanajc/scripts/blobtoolkit/blob_databases/taxdump_2021_07 \
    $DATA_DIR/BlobDir
```

### Build up BlobDir

Before running blobtools, need to map reads to assembly, get blast hits and diamond hits. The following scripts can be run in parallel, and they create the necessary input files to run blobtools right after, adding all these data to the BlobDir.

From specific assembly folder (**assemblies/spp_folder/flye**):
```bash
sbatch --array=1 ../../../../scripts/blob-blast.sh
sbatch --array=1 ../../../../scripts/blob-diamond.sh
sbatch --array=1 ../../../../scripts/blob-mapILL.sh
sbatch --array=1 ../../../../scripts/blob-mapONT.sh
sbatch --array=1 ../../../../scripts/blob-busco.sh metazoa_odb10
# eukaryota_odb10
```

### Run BlobTools

From specific assembly folder (**assemblies/spp_folder/flye**):
```bash
sbatch --array=1 ../../../../scripts/blobtools.sh taxid
```

### View blob

https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/opening-a-dataset-in-the-viewer/

View step was giving error in my installation of blobtools2, something about npm. But docker version worked without any issues.

From **blobtools** folder inside a specific assembly:
```bash
salloc -n 16 -p test -t 0-04:00 --mem=10000
singularity exec --cleanenv ../../../../../scripts/blobtoolkit_2.6.5.sif blobtools view --remote BlobDir
```

Then open in terminal in local machine:
```bash
ssh -N -L 8001:localhost:8001 -L 8000:localhost:8000 -J <myuser>@odyssey.fas.harvard.edu holy7c18312
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


## Scaffold assembly with Hi-C data

Arima protocol for mapping Hi-C reads to an assembly: https://github.com/ArimaGenomics/mapping_pipeline/blob/master/Arima_Mapping_UserGuide_A160156_v02.pdf
<br>Then use YaHS to scaffold.

### Map Hi-C reads to assembly

Download Arima scripts in **scripts** folder:
```bash
git clone https://github.com/ArimaGenomics/mapping_pipeline.git
mv mapping_pipeline arima_mapping
# picard in conda did not work, because arima's command uses the java version, so:
wget https://github.com/broadinstitute/picard/releases/download/2.25.7/picard.jar
```

Create conda environment with required programs: bwa, samtools, picard:
```bash
module load Anaconda3
conda create -n hic -c conda-forge -c bioconda python=3.8 samtools==1.10 picard==2.26.10 bwa=0.7.17
```

There are 4 parts to the mapping pipeline: indexing the assembly, mapping HiC reads to assembly, combining replicate reads from different lanes, and removing duplicates.

1- Indexing

Create bwa index of the assembly.

- From specific assembly folder (**assemblies/spp_folder/flye**):
```
sbatch ../../../../scripts/hic-indexing.sh <post-blobtools-assembly-file>
```

2- Mapping

From Arima manual:<br>
STEP 1: Use BWA-MEM to align the Hi-C paired-end reads to the reference sequences. Because Hi-C captures conformation via proximity-ligated fragments, paired-end reads are first mapped independently (assingle-ends) usingBWA-MEM and are subsequently paired in a later step.<br>
STEP 2: Subsequent to mapping as single-ends, some of these single-end mapped reads can manifest a ligation junction and are therefore considered 'chimeric' (i.e. they do not originate from a contiguous piece of DNA). When BWA-MEM maps these chimeric reads, there can be high quality mapping on both the 5’- side and 3’-side of the ligation junction within a given read. In such cases, only the 5’-side should be retained because the 3’-side can originate from the same contiguous DNA as the 5’-side of the reads mate-pair. Therefore, we retain only the portion of the chimeric read that maps in the 5’-orientation in relation to its read orientation. This is accomplished using the script `filter_five_end.pl`.<br>
STEP 3: After filtering, we pair the filtered single-end Hi-C reads using `two_read_bam_combiner.pl`, which outputs a sorted, mapping quality filtered, paired-end BAM file.
STEP 4: We then add read groups to this BAM file using Picard Tools.

- From specific assembly folder (**assemblies/spp_folder/flye**):

```bash
# get list of HiC fastq file names and save in file inside hic folder - CHANGE NAME OF TAXA
\ls /n/Giribet_Lab/tauanajc/GASTROPODA/rawdata/2020_03-NovaSeqS4/Tauana/HiC-Acutogordius_australiensis_152393_S6_L00* | rev | cut -d/ -f 1 | cut -d_ -f3- | rev | uniq > hic/fastq-names-hic.txt

# run
sbatch --array=1 ../../../../scripts/hic-map.sh
```

3/4- Combine lanes and remove duplicates

From Arima's manual:<br>
We also use Picard Tools to discard any PCR duplicates present in the paired-end BAM file generated above. If applicable, we require that you merge paired-end BAM files that were sequenced via multiple Illumina lanes from the same library (i.e. technical replicates) before removing PCR duplicates.

From specific assembly folder:
```bash
sbatch ../../../../scripts/hic-dedup.sh
```

From Arima's manual:<br>
The pipeline is complete. The final output of this pipeline is a single BAM file that contains the paired, 5’-filtered, and duplicate-removed Hi-C reads mapped to the reference sequences of choice. The resulting statistics file has a breakdown of the total number of intra-contig read-pairs, long-range intra-contig read-pairs, and inter-contig read-pairs in the final processed BAM file.


### Scaffold with YaHS

https://github.com/c-zhou/yahs
<br>https://github.com/sanger-tol/yahs

The outputs from the Arima Hi-C Mapping Pipeline are the inputs for YaHS.

YaHS requires a bed file which is converted from bam using Bedtools, an assembly in fasta format, an index for that assembly in .fai (we already made this at the beginning of the mapping pipeline)

```
# convert bam to bed file and sort by read name
bamToBed -i hic/mapping/final-deduplicated-bam/$LABEL.final.bam > hic/$LABEL.bed
sort -k 4 hic/$LABEL.bed > tmp && mv tmp hic/$LABEL.bed
```

YaHS runs pretty fast on nematomorph genomes. After that, juicer tools to create the file used for visualization of the contact map. I followed instructions on c-zhou's github for the hole thing (yahs and juicer). One of the commands is to create the file for scaffold sizes (should contain two columns - scaffold name and scaffold size): the github says to create from first two columns of the index file from the yahs scaffolds fasta, but this index is not needed otherwise, so it doesn't exist at this point. In one of the github issues, however, YaHS author Chenxy says to create this file by grepping the juicer log file - it is the preferred method anyway (deals with >2G chromosomes, https://github.com/c-zhou/yahs/issues/4#issuecomment-1159017424)

From specific assembly folder (**assemblies/spp_folder/flye**):
```bash
sbatch ../../../../scripts/hic-yahs.sh
# last juicer step is long and memory intensive (12h, >100G for Acutogordius)
```

### Visualize contact map with Juicebox

https://github.com/aidenlab/Juicebox/wiki

Download out.hic file and open with Juicebox (online: https://aidenlab.org/juicebox/ or java or Desktop: https://github.com/aidenlab/Juicebox/wiki/Download). The online version allows to share a link to put in the paper, so that people can load and play with your contact map.


## Check assembly quality with BUSCO

Database of universal orthologs from different sets of organisms. Based on evolutionarily-informed expectations of gene content of near-universal single-copy orthologs.

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

Of interest: metazoa_odb10. eukaryota_odb10

- From **assemblies** folder:
```bash
mkdir -p QC/busco/assemblies
#copy final assemblies to be evaluated into new QC/busco/assemblies folder, then
cd QC/busco
sbatch ../../../../scripts/busco-genomes.sh <assembly.file/folder> <lineage>
#Example:
sbatch busco-genomes.sh assemblies metazoa_odb10
```

Important output files

Inside folder:<br>
busco/out_metazoa_odb10/run_metazoa_odb10/

    - short_summary.txt
    - full_table.tsv
    - missing_busco_list.tsv
