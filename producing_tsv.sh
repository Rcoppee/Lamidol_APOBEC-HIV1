#! /bin/bash

#Script Romain Coppee
#Creation data: 29/10/2022
#Last modification: 29/10/2022

###############################################################################################################

#Location of the Reference genome
REF_GEN=HIV1_ref.fasta

#FILES contains the fixed, indexed BAM file for each sample
FILES_FIX=*.fix

#FILES contains the BAM file for each sample
FILES_BAM=*.bam

#FILES contains the indexed BAM file for each sample
FILES_BAI=*.bai

#FILES contains the fixed, indexed BAM file for each sample
FILES_FIX=*.fix

FILES_PILEUP=*.pileup

FILES_TSV=*.tsv

###############################################################################################################

#we assume the name of the reference sequence is HXB2
#here is the example for the RT prot
#for vif gene, replace the coordinates by 5041-5619
for f in `ls $FILES_BAM`
do
    current_name=$(basename $f .sorted.bam)
    samtools mpileup -d 1000000 -f $REF_GEN -r HXB2:2550-3326 $f -o $current_name.pileup
done

for f in `ls $FILES_PILEUP`
do
    current_name=$(basename $f .pileup)
    echo "$current_name"
    perl /home/virologie/Documents/scripts/extract_info_from_pileup.pl -p $f -o ./
    mv info_from_pileup.tsv $current_name.tsv
    awk -v new="$current_name" '$1=new' < $current_name.tsv > $current_name.tsv.fix
done

rm $FILES_PILEUP
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "fixing PROCESSED $f"
done

for f in `ls $FILES_TSV`
do
    current_name=$(basename $f .tsv)
    echo "$f"
    sed -i '1 i\Id_sample Position n_A n_C n_G n_T n_indel major_nucl complement qual' $f
    awk -v OFS="\t" '$1=$1' $f > $f.fix
done

rm $FILES_TSV
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "fixing PROCESSED $f"
done
