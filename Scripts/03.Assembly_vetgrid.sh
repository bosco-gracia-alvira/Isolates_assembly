#!/bin/bash
# This script assembles genomes from the desired pool of isolates
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
POOL="Pool_$1"
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"

READS="$WORKDIR"/"$POOL"/02.Rm_adapters/fastq_clean
METADATA="$WORKDIR"/"$POOL"/metadata.tsv
ASSEMBLY="$WORKDIR"/"$POOL"/03.Assembly

### COMMANDS
IFS="
"

if [[ -z "$POOL" ]]
then
    echo "You need to choose a pool, idiot! "
    exit
fi

if [[ ! -d "$WORKDIR"/"$POOL" ]]
then
    echo "You, liar! This pool does not exist. Try again."
    exit
else
    cd "$WORKDIR"/"$POOL"
fi

if [[ ! -d "$ASSEMBLY"/Assembled_genomes ]]
then
    mkdir -p "$ASSEMBLY"/Assembled_genomes
fi

for i in $(cut -f2 "$METADATA" | grep -v "sample")
do
    # This if statement allows to re-start the script if a sample failed
    if [[ ! -f "$ASSEMBLY"/Assembled_genomes/${i}.fasta ]]
    then
        # Assemble each reads set with SPAdes, my favourite assembler
        spades.py \
            -1 "$READS"/${i}.clean_1.fq.gz \
            -2 "$READS"/${i}.clean_2.fq.gz \
            -o "$ASSEMBLY"/SPAdes/${i} \
            --meta \
            -k 21,33,55,77,99,127 \
            --threads 12
        
        # Remove contigs <2000 bp, change the contigs name and save the file to "Assembled genomes"
        cat "$ASSEMBLY"/SPAdes/${i}/scaffolds.fasta |\
            seqkit seq -m 2000 | seqkit replace -p .+ -r "${i}_{nr}" --nr-width 5 > "$ASSEMBLY"/Assembled_genomes/${i}.fasta
    fi
done

#Let's now look at N50 and L50. I use the program "assembly-stats". 
# Since it is not implemented in Homebrew I have just installed it in conda base environment...

conda activate base

rm "$ASSEMBLY"/*.col

for i in $(cut -f2 $METADATA | grep -v "sample")
do 
    echo "${i}" >> "$ASSEMBLY"/Name.col
    assembly-stats "$ASSEMBLY"/Assembled_genomes/${i}.fasta | grep "N50" >> "$ASSEMBLY"/N50.col
    assembly-stats "$ASSEMBLY"/Assembled_genomes/${i}.fasta | grep "sum" >> "$ASSEMBLY"/Size.col
done

paste \
    "$ASSEMBLY"/Name.col \
    "$ASSEMBLY"/N50.col \
    "$ASSEMBLY"/Size.col > "$ASSEMBLY"/Assembly_stats.txt

rm "$ASSEMBLY"/*.col

conda deactivate
