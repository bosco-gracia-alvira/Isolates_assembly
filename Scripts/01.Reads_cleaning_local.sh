#!/bin/bash
# This script cleans the raw reads from the desired pool of isolates
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
POOL="Pool_$1"
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"
LOCAL="/Users/bgracia/PhD_local/Isolates_assembly"
ADAPTERS="/Users/bgracia/PhD_local/db/Adapters"

RAW="$LOCAL"/"$POOL"/01.Raw_data/Demultiplexed
RM="$LOCAL"/"$POOL"/02.Rm_adapters
METADATA="$WORKDIR"/"$POOL"/metadata.tsv

### COMMANDS
IFS="
"

if [[ -z "$POOL" ]]
then
    echo "You need to choose a pool, idiot! "
    exit
fi

if [[ ! -d "$LOCAL"/"$POOL" ]]
then
    echo "You, liar! This pool does not exist. Try again."
    exit
else
    cd "$LOCAL"/$POOL
fi

if [[ ! -d "$RM"/fastq_clean ]]
then
    mkdir -p "$RM"/fastq
    mkdir -p "$RM"/fastq_wo_adapt
    mkdir -p "$RM"/fastq_clean
    mkdir -p "$RM"/fastqc_clean
fi

if [[ "$POOL" = "Pool_503" ]]
then
    # In this pool the raw reads were stored in bam format. We have to sort the bam files and convert them to fastq.
    while IFS=$'\t' read -r pool_name pool sample _rest
    do
        # Skip the header line
        if [[ "$pool_name" != "pool_name" ]]
        then
            # Sort the BAM file by read names
            samtools sort -n "$RAW/${pool_name}" |\
            
            # Convert the sorted BAM file to paired-end FASTQ files
            samtools fastq -1 "$RM/fastq/${sample}_1.fq.gz" -2 "$RM/fastq/${sample}_2.fq.gz"
        fi
    done < "$METADATA"

elif [[ "$POOL" = "Pool_591" ]]
then
    while IFS=$'\t' read -r pool_name pool sample _rest
    do  
        # Skip the header line
        if [[ "$pool_name" != "pool_name" ]] 
        then
            samtools collate \
                --reference "$RAW/dsimM252v1.2+microbiome.fa" \
                -o "$RAW/${pool_name}.sorted.cram" "$RAW/${pool_name}.cram"

            samtools fastq \
                --reference $RAW/dsimM252v1.2+microbiome.fa \
                -1 "$RM/fastq/${sample}_1.fq.gz" \
                -2 "$RM/fastq/${sample}_2.fq.gz" \
                "$RAW/${pool_name}.sorted.cram"
        fi
    done < "$METADATA"

elif [[ "$POOL" =~ Pool_6[0-9][0-9] ]]
then
    # In this pool the raw reads were stored in fastq format. We link them to save space.
    rm "$RM"/fastq/*_?.fq.gz
    for i in $(cut -f3 $METADATA | grep -v "sample")
    do
        ln -s "$RAW"/${i}/*_1.fq.gz "$RM"/fastq/${i}_1.fq.gz
        ln -s "$RAW"/${i}/*_2.fq.gz "$RM"/fastq/${i}_2.fq.gz
    done

fi

# We loop over the samples to clean the reads
for i in $(cut -f3 $METADATA | grep -v "sample")
do  
    #We use bbduk, from bbtools, to trim the reads and remove the Illumina adapters.
    bbduk.sh \
        -Xmx24g \
        in1="$RM/fastq/${i}_1.fq.gz" \
        in2="$RM/fastq/${i}_2.fq.gz" \
        out1="$RM/fastq_wo_adapt/${i}.RmAdp_1.fq.gz" \
        out2="$RM/fastq_wo_adapt/${i}.RmAdp_2.fq.gz" \
        ref="$ADAPTERS/adapters.fa" \
        ktrim=r k=23 mink=11 hdist=1 tbo tpe

    # Trim any PhiX174 sequence and remove low quality reads.
    bbduk.sh \
        -Xmx24g \
        in1="$RM/fastq_wo_adapt/${i}.RmAdp_1.fq.gz" \
        in2="$RM/fastq_wo_adapt/${i}.RmAdp_2.fq.gz" \
        out1="$RM/fastq_clean/${i}.clean_1.fq.gz" \
        out2="$RM/fastq_clean/${i}.clean_2.fq.gz" \
        ref="$ADAPTERS/phix174_ill.ref.fa" \
        k=31 hdist=1 qtrim=rl trimq=20

    # Assess the quality of the reads after trimming
    fastqc \
        "$RM"/fastq_clean/${i}.clean_1.fq.gz \
        "$RM"/fastq_clean/${i}.clean_2.fq.gz \
        -o "$RM"/fastqc_clean
done

# We move everything to the Dropbox
if [[ ! -d "$WORKDIR"/"$POOL"/02.Rm_adapters ]]
then
    mkdir "$WORKDIR"/"$POOL"/02.Rm_adapters
fi

mv "$RM"/fastq "$WORKDIR"/"$POOL"/02.Rm_adapters/
mv "$RM"/fastq_clean "$WORKDIR"/"$POOL"/02.Rm_adapters/
mv "$RM"/fastqc_clean "$WORKDIR"/"$POOL"/02.Rm_adapters/

#rm -r "$RM"