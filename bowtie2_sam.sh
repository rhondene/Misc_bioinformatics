#!/bin/bash

# Author: RWint
#convenient script to align to with bowtie to sort sam files to bam then and index files


read -p 'name of input fastq': file; read -p 'name of reference': ref; read -p 'name of output prefix': out; bowtie2 --threads 4 --quiet -x /home/alchemis/Desktop/trnaseq_data/tRNAseq_metrics/refsets/$ref -U $file -S $out.sam; samtools view -S -b $out.sam > ./$out.bam; samtools sort -o ./$out.sorted.bam ./$out.bam; samtools index ./$out.sorted.bam ; samtools flagstat ./$out.sorted.bam; rm ./$out.sam


