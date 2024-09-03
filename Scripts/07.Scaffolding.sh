#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate ragtag-2.1

# mkdir -p 08.Scaffolding/Acetobacter_indonesiensis

echo 'Which bacteria do you want to scaffold? Type in the terminal the species with the format: "Genus_species"'
read
SPECIES="Acetobacter indonesiensis"
SPECIES=$(echo $REPLY | sed 's/_/ /')

ncbi-genome-download bacteria \
    -g "$SPECIES" \
    -s refseq \
    -F fasta \
    -l "all"\
    -P \
    --flat-output \
    -o . 
gunzip *.gz

#08.Scaffolding/$REPLY

for i in $(basename *);
    do echo ${i} >> Name.col;
    assembly-stats ${i} | grep "N50" >> N50.col;
    assembly-stats ${i} | grep "sum" >> Size.col;
done

paste -d "\t" Name.col N50.col Size.col > Contig_stats.txt

rm *.col

echo 'Which isolate do you want to scaffold? These are the available assemblies'
read ISOLATE

cp ../Isolates_assembly/Pool_$(echo $ISOLATE | cut -f1 -d "_")/07.GTDB-Tk/Genomes/$ISOLATE.fa .

# scaffold a query assembly
ragtag.py scaffold -o single GCF_000963945.1_ASM96394v1_genomic.fna 503_G175_R4_T1_MRS_4.fa

# scaffold with multiple references/maps
x=1
for i in $(basename *.fna)
do  ragtag.py scaffold -o multiple_$x $i $ISOLATE.fa;
    x=$x+1;
done

ragtag.py merge $ISOLATE.fa multiple_*/*.agp -o final_multiple

declare -i x=1

for i in $(basename *.fna)
do  echo $i;
    x=$x+1;
    echo $x;
done