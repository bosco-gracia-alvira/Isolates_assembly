#!/bin/bash
# This script creates a graph pangenome from all the isolates of a species
# Bosco Gracia Alvira, 2024

### VARIABLES
if [ $# -ne 2 ]; then
    echo "USAGE: $0 <POOL> <ISOLATE>"
    exit 1
fi

POOL="Pool_$1"
GENOME="$2"
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"
CLEANING="$WORKDIR"/"$POOL"/05.Genome_cleaning

### COMMANDS
IFS="
"

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

anvi-interactive -p "$CLEANING/$GENOME/Profile/PROFILE.db" -c "$CLEANING/$GENOME/contigs_$GENOME.db"

#anvi-display-contigs-stats "$CLEANING/$GENOME/contigs_$GENOME.db"
