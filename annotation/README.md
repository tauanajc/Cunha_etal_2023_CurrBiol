# Genome Annotation

Scripts used for genome annotation, by Arianna Lord.


## Create library of repetitive elements with RepeatModeler2

https://github.com/Dfam-consortium/RepeatModeler

Input: final genome assembly.

`repeatmodeler.sh`


## Mask repeats with RepeatMasker

https://github.com/rmhubley/RepeatMasker

Input: output library from RepeatModeler2 and final genome assembly.

`repeatmasker.sh`


## Map RNAseq to assembly with HISAT2

https://github.com/DaehwanKimLab/hisat2

Input: masked assembly and RNAseq reads.

`hisat2.sh`


## Gene prediction with BRAKER2

https://github.com/Gaius-Augustus/BRAKER

Input: masked assembly and mapped RNAseq.

`braker2.sh`


## Functional annotation with eggNOG-mapper

https://github.com/eggnogdb/eggnog-mapper

Input: fasta of predicted genes.

`eggnog.sh`
