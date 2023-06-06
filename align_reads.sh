#! /bin/bash

#Script Romain Coppee
#Creation data: 29/10/2022
#Last modification: 29/10/2022

###############################################################################################################

#Location of the Reference genome
REF_GEN=HIV1_ref.fasta

#FILES contains the fastq file for each sample
FILES_FASTQ=*.fastq.gz

#FILES contains the SAM file for each sample
FILES_SAM=*.sam

#FILES contains the BAM file for each sample
FILES_BAM=*.bam

###############################################################################################################

for f in `ls $FILES_FASTQ`
do
    current_name=$(basename $f .fastq.gz)
    echo "$current_name"
    minimap2 -ax map-ont $REF_GEN $f > $current_name.sam
    echo "MINIMAP2 PROCESSED $f"
    rm $f
done

for f in `ls $FILES_SAM`
do
    current_name=$(basename $f .sam)
    samtools view -b -S $f > $current_name.bam
    echo "SAM to BAM PROCESSED $f"
    rm $f
done

for f in `ls $FILES_BAM`
do
    current_name=$(basename $f .bam)
    samtools sort $f -o $current_name.sorted.bam
    echo "SORTING PROCESSED $f"
    rm $f
done

for f in `ls $FILES_BAM`
do
    samtools index $f
    echo "INDEXATION PROCESSED $f"
done
