#!/bin/bash

# 

mkdir 06.Binning



#These are the IDs of the clean and dirty isolates based on CheckM results
CLEAN=($(awk -F'\t' 'BEGIN {ORS = " "} $2 !~ /Marker lineage/ && $2 !~ /k__/ && $2 !~ /root/ && $13 < 5 {print $1}' 04.CheckM/Completion.txt))
DIRTY=($(awk -F'\t' 'BEGIN {ORS = " "} $2 !~ /Marker lineage/ && ($2 ~ /k__/ || $2 ~ /root/ || $13 > 5) {print $1}' 04.CheckM/Completion.txt))


#Although Anvi'o is very nice to verify that there is contamination. But if I want to extract the bins in a correct way, I need to use an automatic, reproducible, metagenomic binner.
#I will use METABAT2 because is among the best softwares (CAMI 2022) and has been the easiest to install because it has a conda environment.
eval "$(conda shell.bash hook)"
conda activate metabat-2

 cd 06.Binning

for i in $DIRTY;
do runMetaBat.sh \
    --seed 999 \
    ../05.Genome_cleaning/Anvio_${i}/contigs_${i}.fasta \
    ../05.Genome_cleaning/Anvio_${i}/contigs_${i}.sorted.bam;
done

conda deactivate

#Here I change the name of the bins to include the original assembly from which they come from:
for i in $(ls | grep "meta")
do mv "${i}" $(echo "${i}" | cut -d . -f1);
done

rm *.fasta*

for i in $DIRTY;
do
    for x in $(seq 1 $(ls -l contigs_${i}/ | grep -v ^total | wc -l)); 
        do mv $(basename -a contigs_${i})/bin.${x}.fa $(basename -a contigs_${i})/${i}_bin_${x}.fa;
    done
done

#I put all the bins in the same folder (inside 06.Binning directory still) in order to calculate their completeness, contamination, N50 and L50
mkdir Bins
for i in $DIRTY;
do cp $(basename -a contigs_${i})/*.fa Bins;
done

#Completeness and contamination with CheckM
mkdir CheckM

conda activate checkm-1.2
export CHECKM_DATA_PATH=/Users/bgracia/PhD/db/CheckM

checkm lineage_wf -t 8 -x fa Bins CheckM --tab_table
checkm qa CheckM/lineage.ms CheckM --tab_table > CheckM/Completion.tmp
sed '1,7d' CheckM/Completion.tmp | sed '$d' > CheckM/Completion.txt
rm CheckM/Completion.tmp
conda deactivate

#N50, L50 and final size with assembly-stats
conda activate base
for i in $(ls Bins/*.fa);
    do echo "${i}" | cut -d "/" -f2 | cut -d . -f1 >> Bins/Name.col;
    assembly-stats ${i} | grep "N50" >> Bins/N50.col;
    assembly-stats ${i} | grep "sum" >> Bins/Size.col;
done

paste -d "\t" Bins/Name.col Bins/N50.col Bins/Size.col > Bins/Contig_stats.txt
rm Bins/*.col
conda deactivate

cd ..