#!/bin/bash

#cd /global/cscratch1/sd/rwint/subread_align/align_bam


for lib in CGATGT TTAGGC TGACCA
do
        touch "$lib"_multi_reads.stats.txt
	while read -r read_id
	do
		#check if read_id is multimapped
		if  samtools view -@ 16 $lib.mapped.sorted.bam | grep -w $read_id | grep -qvw NH:i:1
		then
			#extract the chromosomes it mapped to#
			echo 'first success'
			chrs=$(samtools view -@ 16 $lib.mapped.sorted.bam | grep -w $read_id | cut -f3 | uniq)  #first alignment
			chrs+=$(samtools view -@ 16 $lib.mapped.sorted.bam | grep -w $read_id | cut -f7 | sed '/FBtr/!d') #2 aligns
			#extract cigar info
			cigs=$(samtools view -@ 16 $lib.mapped.sorted.bam | grep -w $read_id | cut -f6)
			#write to file
			echo -e "$read_id\t[$chrs]\t[$cigs]" >>./"$lib"_multi_reads.stats.txt
		fi
	done < "$lib"_uniq_reads.txt
done
