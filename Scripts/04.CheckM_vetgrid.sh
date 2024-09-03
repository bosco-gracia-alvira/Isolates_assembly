#!/bin/bash
# This script assesses the completeness and contamination using CheckM2.
# It requires connexion to the VetLinux server since my local computer does not support CheckM2.
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
POOL="Pool_$1"
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"

ASSEMBLY="$WORKDIR"/"$POOL"/03.Assembly
CHECKM="$WORKDIR"/"$POOL"/04.CheckM2

### COMMANDS
IFS="
"

ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

if [[ ! -d ~/Bosco/Isolates_assembly/$POOL/04.CheckM2 ]]
then  
    mkdir -p ~/Bosco/Isolates_assembly/$POOL/04.CheckM2
fi

if [[ ! -d ~/Bosco/Isolates_assembly/$POOL/Assembled_genomes ]]
then  
    mkdir -p ~/Bosco/Isolates_assembly/$POOL/Assembled_genomes
fi

FOO

# Copy the genomes to the server
scp -r "$ASSEMBLY"/Assembled_genomes/*.fasta vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Isolates_assembly/$POOL/Assembled_genomes

# We connect to the server again to run CheckM2
ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

cd ~/Bosco/Isolates_assembly/$POOL

eval \$(conda shell.bash hook)
conda activate checkm2

CHECKM2DB="/home/vetlinux05/Bosco/db/CheckM2_database/uniref100.KO.1.dmnd"

checkm2 predict \
    --threads 20 \
    --input ./Assembled_genomes \
    -x fasta \
    --force \
    --output-directory ./04.CheckM2

FOO

# I copy the results back to my computer
scp -r vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Isolates_assembly/$POOL/04.CheckM2 "$WORKDIR"/"$POOL"/
