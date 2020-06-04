#!/bin/bash

# Author: RWint
#convenient script to pull unmapped reads from a bam file into a fastq output for re-alignment

read -p 'name of input sorted bam file': file; read -p 'name of output': out; read -p 'name of fastq file':fastq; samtools view -f4 $file.sorted.bam > $out.unmapped.bam; cut -f1 $out.unmapped.bam |sort| uniq > unmapped_IDs.txt; seqtk subseq $fastq.fq.gz unmapped_IDs.txt > $out.unmapped.fq

