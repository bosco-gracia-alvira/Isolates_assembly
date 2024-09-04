#!/bin/bash
# This script calls Anvi'o to visualise the isolates' assembly.
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
POOL="Pool_$1"
GENOME="$2"
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"
CLEANING="$WORKDIR"/"$POOL"/05.Genome_cleaning


### COMMANDS
IFS="
"

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

anvi-interactive -p "$CLEANING/$GENOME/Profile/PROFILE.db" -c "$CLEANING/$GENOME/contigs_$GENOME.db"

# This command displays the assembly statistics instead of the profile
#anvi-display-contigs-stats "$CLEANING/$GENOME/contigs_$GENOME.db"
