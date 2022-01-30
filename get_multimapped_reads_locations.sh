#!/bin/bash

libs=("CGATGT" "TTAGGC" "TGACCA")
for lib in ${libs[@]}
#create file to store read info for current library
touch $lib_multi_reads.stats.txt
while read -r read_id;
do
		#check if read is multimapped
		if samtools view -@ 16 $lib.mapped.sorted.bam | grep -w "$read" | grep -qvw "NH:i:1"
		then
			#extract the chromosomes it mapped to
			chrs=$(samtools view -@ 16 $lib.mapped.sorted.bam | grep "$read" | cut -f3 | uniq)  #first alignment
			chrs+=$(samtools view -@ 16 $lib.mapped.sorted.bam | grep "$read" | cut -f7 | sed '/FBtr/!d') #2 aligns
			#extract cigar info
			cigs=$(samtools view -@ 16 $lib.mapped.sorted.bam | grep "$read" | cut -f6) 
			#write to file
			echo -e "$read\t[$chrs]\t[$cigs]" >>$lib_multi_reads.stats.txt
		fi
	done < "lib"_uniq_reads.txt
done


  
  