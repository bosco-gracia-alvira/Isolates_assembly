#!/bin/bash
# This script assesses the taxonomy and completeness of the final bins.
# It requires connexion to the VetLinux server since my local computer does not support CheckM2 nor GTDB-TK.
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
POOL="Pool_$1"
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"

ASSEMBLY="$WORKDIR"/"$POOL"/03.Assembly
CHECKM="$WORKDIR"/"$POOL"/04.CheckM2
CLEANING="$WORKDIR"/"$POOL"/05.Genome_cleaning
GTDBTK="$WORKDIR"/"$POOL"/07.GTDB-Tk

### COMMANDS
IFS="
"

ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

if [[ ! -d ~/Bosco/Isolates_assembly/$POOL/07.GTDB-Tk/Genomes ]]
then  
    mkdir -p ~/Bosco/Isolates_assembly/$POOL/07.GTDB-Tk/Genomes
fi

FOO

DIRTY=$(awk -F'\t' '$3 > 5 {print $1}' "$CHECKM"/quality_report.tsv | grep -v "Name")
CLEAN=$(awk -F'\t' '$3 <= 5 {print $1}' "$CHECKM"/quality_report.tsv | grep -v "Name")

# Copy the genomes to the server
for i in $CLEAN
do
        scp -r "$ASSEMBLY/Assembled_genomes/${i}.fasta" \
        vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Isolates_assembly/$POOL/07.GTDB-Tk/Genomes
done

for i in $DIRTY
do
        scp -r "$CLEANING"/${i}/Profile/SUMMARY_default/bin_by_bin/${i}_?/${i}_?-contigs.fa \
        vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Isolates_assembly/$POOL/07.GTDB-Tk/Genomes
done

# Bins cleaned from contamination are in fa format. I change the extension to fasta
ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO
POOL=$POOL
for i in ~/Bosco/Isolates_assembly/\$POOL/07.GTDB-Tk/Genomes/*.fa
do
    mv "\${i}" "\${i%.fa}.fasta"
done
FOO

# We connect to the server again to run GTDB-TK2
ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

cd ~/Bosco/Isolates_assembly/$POOL

eval \$(conda shell.bash hook)
conda activate gtdbtk-2.1.1

export GTDBTK_DATA_PATH="/home/vetlinux05/Bosco/db/gtdbtk_r214_database"

gtdbtk classify_wf \
        --genome_dir 07.GTDB-Tk/Genomes \
        -x fasta \
        --out_dir 07.GTDB-Tk/output

FOO

# I copy the results and the genomes back to my computer
scp -r vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Isolates_assembly/$POOL/07.GTDB-Tk/output "$GTDBTK"
mkdir "$GTDBTK"/Genomes
scp -r vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Isolates_assembly/$POOL/Assembled_genomes/*.fasta \
        "$GTDBTK"/Genomes

# #We have finally recovered as many clean isolates as possible. Now I will gather them together and run GTDB-Tk in order to know which reference genome I should use for the scaffolding
# eval "$(conda shell.bash hook)"
# conda activate gtdbtk-2.0
# export GTDBTK_DATA_PATH=/Users/bgracia/PhD_local/db/GTDB/release207

# if [[ ! -d "$LOCAL"/Genomes ]]
# then
#     mkdir -p "$LOCAL"/Genomes
# fi

# for i in $CLEAN;
# do 
#         cp "$ASSEMBLY"/Assembled_genomes/${i}.fasta "$LOCAL"/Genomes
# done

# for i in $DIRTY
# do 
#         cp "$CLEANING"/${i}/Profile/SUMMARY_default/bin_by_bin/${i}_*/${i}_*-contigs.fa "$LOCAL"/Genomes
# done

# for i in $(basename -a "$LOCAL"/Genomes/*fa)
# do
#         mv "$LOCAL"/Genomes/${i} "$LOCAL"/Genomes/${i}sta
# done

# for i in $(basename "$LOCAL"/Genomes/*.fasta | cut -d . -f1)
# do      
#         mkdir "$LOCAL"/${i}
#         cp "$LOCAL"/Genomes/${i}.fasta "$LOCAL"/${i}
#         gtdbtk classify_wf \
#                 --genome_dir "$LOCAL"/${i} \
#                 -x fasta \
#                 --out_dir "$LOCAL"/${i} \
#                 --cpus 16
# done

# cat \
#     $(
#     for i in $(basename "$LOCAL"/Genomes/*.fasta | cut -d . -f1)
#     do 
#         echo "$LOCAL/${i}/classify/gtdbtk.bac120.summary.tsv";
#     done
#     ) | sort -r | uniq -i > "$LOCAL"/summary.tsv
