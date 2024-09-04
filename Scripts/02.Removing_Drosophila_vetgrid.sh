#!/bin/bash
# This script cleans the raw reads from the desired pool of isolates
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
POOL="Pool_$1"
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"
EUK=""$WORKDIR"/euk/contaminants"

READS="$WORKDIR"/"$POOL"/02.Rm_adapters/fastq_clean
METADATA="$WORKDIR"/"$POOL"/metadata.tsv

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
    cd "$WORKDIR"/$POOL
fi

if [[ ! -d "$WORKDIR"/euk ]]
then
    mkdir -p "$WORKDIR"/euk
fi

# We create a reference genome by concatenating D. melanogaster (v6.52) and D. simulans (v2.02) genomes from flybase.org, as well as other eukaryotic genomes that could potentially be contaminating our samples based on previous experiences.
# These eukaryotes are: H. sapiens, M. musculus, A. thaliana, S. cerevisiae and C. lupus familiaris. We download the latest reference genome from RefSeq.
if [[ ! -f "$WORKDIR"/euk/euk.fa ]]
then
    # D. simulans
    wget -O "$WORKDIR"/euk/dsim202.fa.gz \
        http://ftp.flybase.net/genomes/Drosophila_simulans/dsim_r2.02_FB2020_03/fasta/dsim-all-chromosome-r2.02.fasta.gz
    
    # D. melanogaster
    wget -O "$WORKDIR"/euk/dmel652.fa.gz \
        http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.52_FB2023_03/fasta/dmel-all-chromosome-r6.52.fasta.gz
    
    # H. sapiens
    wget -O "$WORKDIR"/euk/Hsapiens.fa.gz \
        ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz

    # M. musculus
    wget -O "$WORKDIR"/euk/Mmusculus.fa.gz \
        https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/reference/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz

    # A. thaliana
    wget -O "$WORKDIR"/euk/Athaliana.fa.gz \
        https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Arabidopsis_thaliana/reference/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz

    # S. cerevisiae
    wget -O "$WORKDIR"/euk/Scerevisiae.fa.gz \
        https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/reference/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

    # C. lupus
    wget -O "$WORKDIR"/euk/Clupus.fa.gz \
        https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Canis_lupus_familiaris/representative/GCF_011100685.1_UU_Cfam_GSD_1.0/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna.gz

    # We concatenate the 7 genomes into a single file and index it.
    gzcat "$WORKDIR"/euk/dsim202.fa.gz \
          "$WORKDIR"/euk/dmel652.fa.gz \
          "$WORKDIR"/euk/Hsapiens.fa.gz \
          "$WORKDIR"/euk/Mmusculus.fa.gz \
          "$WORKDIR"/euk/Athaliana.fa.gz \
          "$WORKDIR"/euk/Scerevisiae.fa.gz\
          "$WORKDIR"/euk/Clupus.fa.gz > "$WORKDIR"/euk/contaminants.fa

    # We index the reference
    bwa index -p "$EUK" "$WORKDIR"/euk/contaminants.fa

    # We don't want the genomes wasting space anymore
    rm "$WORKDIR"/euk/*.fa.gz

fi

# We align them to the reference file

for i in $(cut -f3 "$METADATA" | grep -v "sample")
do

    bwa mem \
        -M \
        -t 24 \
        "$EUK" \
        "$READS"/${i}.clean_1.fq.gz \
        "$READS"/${i}.clean_2.fq.gz |\
    samtools view \
        -Sbh -F 0x100 -@ 24 -f 0xc - |\
    samtools collate -Ou -@ 24 - |\
    samtools fixmate -u -@ 24 - - |\
    samtools view -u -@ 24 -f 0x1 - |\
    samtools fastq -@ 24 -N -0 /dev/null -s /dev/null \
        -1 "$READS"/${i}.noDro_1.fq.gz \
        -2 "$READS"/${i}.noDro_2.fq.gz \
            -
done

for i in $(cut -f3 "$METADATA" | grep -v "sample")
do
    mv "$READS"/${i}.noDro_1.fq.gz "$READS"/${i}.clean_1.fq.gz
    mv "$READS"/${i}.noDro_2.fq.gz "$READS"/${i}.clean_2.fq.gz  
done