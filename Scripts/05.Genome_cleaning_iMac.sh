#!/bin/bash
# This script assesses the taxonomy using CheckM
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
POOL="Pool_$1"
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"

RM="$WORKDIR"/"$POOL"/02.Rm_adapters
ASSEMBLY="$WORKDIR"/"$POOL"/03.Assembly
CHECKM="$WORKDIR"/"$POOL"/04.CheckM2
CLEANING="$WORKDIR"/"$POOL"/05.Genome_cleaning
METADATA="$WORKDIR"/"$POOL"/metadata.tsv
LOCAL="/Users/bgracia/PhD_local/Isolates_assembly/centrifuge"

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

if [[ ! -d "$CLEANING" ]]
then
    mkdir -p "$CLEANING"
fi

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

for i in $(cut -f3 "$METADATA" | grep -v "sample")
do      
        mkdir "$CLEANING"/${i}

        anvi-gen-contigs-database \
                -f "$ASSEMBLY"/Assembled_genomes/${i}.fa \
                -o "$CLEANING"/${i}/contigs_${i}.db \
                -n ${i} \
                --num-threads 16

        anvi-run-hmms \
                -c "$CLEANING"/${i}/contigs_${i}.db \
                --num-threads 16

        anvi-run-scg-taxonomy \
                -c "$CLEANING"/${i}/contigs_${i}.db \
                --num-threads 16

        mkdir "$CLEANING"/${i}/Bowtie2_alignment
        
        bowtie2-build \
                -f "$CLEANING"/${i}/contigs_${i}.fasta \
                "$CLEANING"/${i}/Bowtie2_alignment/contigs_${i}
        
        bowtie2 \
                -q \
                -x "$CLEANING"/${i}/Bowtie2_alignment/contigs_${i} \
                -1 "$RM"/fastq_clean/${i}.clean_1.fq.gz \
                -2 "$RM"/fastq_clean/${i}.clean_2.fq.gz \
                -S "$CLEANING"/${i}/contigs_${i}.sam
        
        samtools view -S -b "$CLEANING"/${i}/contigs_${i}.sam > "$CLEANING"/${i}/contigs_${i}.bam
        anvi-init-bam "$CLEANING"/${i}/contigs_${i}.bam -o "$CLEANING"/${i}/contigs_${i}.sorted.bam

        rm "$CLEANING"/${i}/contigs_${i}.sam

        anvi-profile \
                -i "$CLEANING"/${i}/contigs_${i}.sorted.bam \
                -c "$CLEANING"/${i}/contigs_${i}.db \
                --min-contig-length 2000 \
                -o "$CLEANING"/${i}/Profile \
                --cluster-contigs
done

# Some of the assemblies are co-cultures of two or three genomes. I will bin them using centrifuge
DIRTY=$(awk -F'\t' '$3 > 5 {print $1}' "$CHECKM"/quality_report.tsv | grep -v "Name")

# This link has the instructions to install the centrifuge database. I have downloaded the whole nt index
# http://www.ccb.jhu.edu/software/centrifuge/
CENTRIFUGE_BASE="/Users/bgracia/PhD_local/db/centrifuge"

# This stupid software does not allow to work from Dropbox because of the stupid parentheses in the path
if [[ ! -d "$LOCAL" ]]
then  
  mkdir -p "$LOCAL"
fi

for i in $DIRTY
do
    # Extract sequences
    anvi-get-sequences-for-gene-calls \
        -c "$CLEANING"/${i}/contigs_${i}.db \
        -o "$LOCAL"/${i}_sequences.fa

    # Run Centrifuge
    centrifuge \
        -f \
        -x "$CENTRIFUGE_BASE"/p+h+v \
        "$LOCAL"/${i}_sequences.fa \
        -S "$LOCAL"/${i}_centrifuge_results.tsv \
        --report-file "$LOCAL"/${i}_centrifuge_report.tsv \
        --seed 999 \
        --verbose \
        -p 1

    cp "$LOCAL"/${i}_centrifuge*.tsv "$CLEANING"/${i}/

    anvi-import-taxonomy-for-genes \
        -c "$CLEANING"/${i}/contigs_${i}.db \
        -i "$CLEANING"/${i}/${i}_centrifuge_report.tsv \
        "$CLEANING"/${i}/${i}_centrifuge_results.tsv \
        -p centrifuge
done

# After this, you can bin the dirty genomes interactively making use of the centrifuge results.