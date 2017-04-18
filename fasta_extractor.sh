#!/bin/bash

echo "Chambeando..."

input_file=$1
genome="/run/media/winter/Winter_HDrive_1/Databases/hg19/hg19.fa"
samtoolsDir="/run/media/winter/Winter_HDrive_2/Programs/samtools-1.3.1"
lenght=50

head -n1 $input_file > tmp/extended_fastas.build
tail -n+2 $input_file > tmp/original_lines

while read p
do
	#echo $p | gawk 'BEGIN {FS=","; OFS=","} {print $1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12}' - >> tmp/extended_fastas.build
	precollumns=$(echo $p | cut -d"," -f1-2)
	postcollumns=$(echo $p | cut -d"," -f4-)
	
	ALLELES=$(echo $p | cut -d"," -f3 | grep "\[.*\]" -m1 -o)
	CHROM=$(echo $p | cut -d"," -f4)
	POS=$(echo $p | cut -d"," -f5)
	pos1=$((POS-1))
	pos2=$((POS+1))
	posLen1=$((pos1-lenght))                                                                                                  
	posLen2=$((pos2+lenght))
	
	upfasta=$($samtoolsDir/samtools faidx $genome "chr"$CHROM:$posLen1-$pos1 | tail -n 1 | tr '[:lower:]' '[:upper:]')
	downfasta=$($samtoolsDir/samtools faidx $genome "chr"$CHROM:$pos2-$posLen2 | tail -n 1 | tr '[:lower:]' '[:upper:]')
	echo $precollumns,$upfasta$ALLELES$downfasta,$postcollumns >> tmp/extended_fastas.build
done < tmp/original_lines

