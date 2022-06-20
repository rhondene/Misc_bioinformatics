#!/bin/bash
#Author: Rwint
# Parses blastn output and stores the top 10 hit species+E-value+Max_Score when the options -num_alignments and -num_descriptions were used
blast_out=$1
grep -A2 'Sequences' $blast_out |grep -v 'Sequences' | tr -s '[:blank:]' |tr [:blank:] \\t | grep . | grep -v '\--' | head > top_10_blast_hits.txt
cat top_10_blast_hits.txt
echo -e ' \n \n The top unique species are:\n'
cut -f1 top_10_blast_hits.txt | sort | uniq -d