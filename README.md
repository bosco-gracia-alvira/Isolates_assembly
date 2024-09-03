# Isolates assembly pipeline

This project contains scripts for processing pools of short reads and assemble them into whole genomes.

## Description

1. `01.Reads_cleaning_local.sh` and `01.Reads_cleaning_vetgrid.sh`: The pipeline starts with the cleaning and trimming of the reads. In some pools the genomic information is in cram format and it needs to be converted into fastq.
2. `02.Removing_Drosophila.sh`: Drosophila, Canis familiaris and Homo sapiens are common sources of contamination in our reads sets. This script maps the reads sets against a concatenated fasta with common Eukaryotic contaminats and removes those reads that map to the reference. Can only be run in the vetgrid.
3. `03.Assembly_iMac.sh` and `03.Assembly_vetgrid.sh`: Next step is assembling them into contigs with SPAdes.
4. `04.CheckM_iMac.sh` and `04.CheckM_vetgrid.sh`: The completeness and contamination of the contigs is evaluated, determining which reads sets contain a unique taxon or are co-cultures.
5. `05.Genome_cleaning_iMac.sh` and `05.Genome_cleaning_vetgrid.sh`: Co-cultures assemblies are binned manually based on coverage, tetranucleotide frequency and taxonomy of the contigs.
6. The commands `05.Genome_display_iMac.sh` and `05.Genome_display_vetgrid.sh` are used to evaluate the contigs of each assembly.
7. `06.GTDB-Tk_iMac.sh` and `06.GTDB-Tk_vetgrid.sh`: Once there is only 1 species/bin, the taxonomy of each assembly is assessed.
8. `07.Scaffolding.sh`: For some species there are closed publicly available reference genomes, that can be used to scaffold the short reads assembly. This script is not working currently.
