# BUSCO Comparison and Functional Enrichment of GO terms

Code for comparison of missing BUSCOs between nematomorph species, enrichment analyses of GO terms, and related figures.


## Folder contents:

- `buscos_enrichment.Rmd`: R script for BUSCO comparison and enrichment analyses
- `buscos_to_GO.py`: Python script to get Gene Ontology terms for a list of BUSCO genes using the OrthoDB v10 database
- `full_table.tsv`: Output files from [BUSCO](https://github.com/tauanajc/Cunha_etal_2023_CurrBiol/tree/main/assembly%20pipeline#check-assembly-quality-with-busco) used as inputs for the R script
- `metazoa_buscos.csv`: List of all 954 metazoan BUSCO genes, created in the R script and used as input for the python script
- `metazoa_buscos_GOterms.csv`: Output of the python script, used as input in the R script
- `tableS1_EnrichedGOterms.csv`: Output of the R script with significant enriched GO terms (TableS2)

