# Synteny

Code to create Oxford dot plot of _Acutogordius_ genome against the ancestral linkage groups inferred for Metazoa.

odp: https://github.com/conchoecia/odp

Install odp with instructions in github, and also create environment with dependencies:

```bash
conda create -n odp -c conda-forge -c bioconda -c anaconda snakemake matplotlib networkx scipy pandas numpy seaborn blast diamond gawk
```

Required inputs:
- config.yaml
- Genome assembly in .fasta format
- Protein sequences in .fasta format
- A file which details where the proteins are located in the genome, in .chrom format

Create the .chrom file with a script provided by odp using a GFF input from NCBI. Our annotation from BRAKER is in GTF, so first I convert the GTF into GFF, then convert into .chrom.

From inside **synteny/data/Acutogordius_australiensis**:

```bash
python ../../../scripts/gtf_to_gff.py augustus.hints.gtf
# creates augustus.hints.gff
```

odp provides a script called NCBIgff2chrom.py, which uses lines with type "CDS" and containing "protein_id=" in the description to create the .chrom file. Our annotations do not have the protein id, just "transcript_id ", and the description has lots of quotation marks. So I modified their script to remove the quotations and use the transcript_id to generate the .chrom file:

```bash
python ../../../scripts/MYgff2chrom.py augustus.hints.gff > Acutogordius_australiensis-MCZ152393.chrom
```

Create config.yaml inside **synteny** folder. Then run odp with snakemake:
```bash
source activate odp
# vim config.yaml
snakemake --cores 16 --snakefile ../../scripts/odp/scripts/odp > synteny.log 2>&1 &
```

This creates an output folder, inside which there is step2-figures folder with PDFs of the dotplots.
