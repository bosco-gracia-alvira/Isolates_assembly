#!/bin/bash
# This script bins the contaminated genomes
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
POOL="Pool_$1"
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"

CHECKM="$WORKDIR"/"$POOL"/04.CheckM2
CLEANING="$WORKDIR"/"$POOL"/05.Genome_cleaning
BINNING="$WORKDIR"/"$POOL"/06.Binning
METADATA="$WORKDIR"/"$POOL"/metadata.tsv

# Get the IDs of the isolates that have >5% contamination according to CheckM
DIRTY=$(awk -F'\t' '$3 > 5 {print $1}' "$CHECKM"/quality_report.tsv | grep -v "Name")

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

if [[ ! -d "$BINNING" ]]
then
    mkdir -p "$BINNING"
fi

#Although Anvi'o is very nice to verify that there is contamination. But if I want to extract the bins in a correct way, I need to use an automatic, reproducible, metagenomic binner.
#I will use METABAT2 because is among the best softwares (CAMI 2022) and has been the easiest to install because it has a conda environment.
eval "$(conda shell.bash hook)"
conda activate anvio-7.1

cd "$BINNING"

# Run the binning software
for i in $DIRTY
do
    runMetaBat.sh \
    --seed 999 \
    "../03.Assembly/Assembled_genomes/${i}.fasta" \
    "../05.Genome_cleaning/${i}/contigs_${i}.sorted.bam" 

done

# Sort the output in two new folders
mkdir -p "$BINNING/depths"
mkdir -p "$BINNING/bins"
mkdir -p "$BINNING/Anvio_profile"

# Rename the new bins based on their isolate ID and move them to the new bins folder
for i in $(basename -a * | grep "meta")
do  
    bin=1
    isolate=$(echo "${i}" | cut -d . -f1)
    for j in $(basename -a "${i}"/*.fa )
    do
        mv "${i}"/"${j}" "$BINNING/bins/${isolate}-${bin}.fa"
        bin=$((bin+1))
    done
    rm -r "${i}"
done

for i in $(basename -a *.fasta.depth.txt)
do
    mv "${i}" "$BINNING/depths"
done

# Add binning results back to anvio
# For each original isolate we create file that contains which contigs belong to which bin
# Then, we include the new information to anvio
for i in $(basename -a "$BINNING"/bins/*.fa)
do  
    grep ">" "$BINNING/bins/${i}" |\
    sed "s/^>//g" |\
    sed "s/$/\t${i%.fa}/g" > "$BINNING"/${i%.fa}.tmp
done

for i in $DIRTY
do
    cat "$BINNING"/${i}-*.tmp > "$BINNING/Anvio_profile/${i}.txt"
done

for i in $DIRTY
do
    anvi-import-collection \
        "$BINNING/Anvio_profile/${i}.txt" \
        --contigs-mode \
        -p "$CLEANING/${i}/Profile/PROFILE.db" \
        -c "$CLEANING/${i}/contigs_${i}.db" \
        -C default
done

rm "$BINNING"/*.tmp
