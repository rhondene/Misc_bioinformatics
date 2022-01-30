#!/bin/bash
genes=($(cut -f1 mature_trnas_clus_oligo.gtf| head -n 107| tail -n 106))
files=($(ls *.sorted.bam))
echo ${files[*]}

for file in ${files[@]}
do	
	touch ${file:0:6}_map_stats.txt
	echo "outer for loop"
	for tRNA in ${genes[@]}
	do
		echo "inner for loop"
		#count the unique alignments#
		unique=($(samtools view -@ 16 $file $tRNA | grep -wc 'NH:i:1'))
		echo "$unique "
		#count the multimappe#
		multi=($(samtools view -@ 16 $file $tRNA | grep -wvc 'NH:i:1'))
		echo -e "$tRNA\t$unique\t$multi" >>${file:0:6}_map_stats.txt
	done
done
