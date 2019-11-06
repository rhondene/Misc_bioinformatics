## Author Rhondene Wint
## Purpose: To extract and make custom reference fasta file of pre-tRNA genes using gtf and fasta annotation files.
## Alternatively, You can use gffutils package but it takes too long to build its database and throws errors on dmel release 6.27 so I had to bypass that
## You can replace 'tRNA' feature with any other feature of your choice once it is annotated in the gtf file
## gtf and fasta must be same release

import pickle
import pandas as pd
import numpy as np

# load gtf annotation file 

gtf = pd.read_table('dmel-all-r6.27.gtf', header=None, names=['Chr','Database','Feature','Start','End','?', 'Strand','?', 'Attributes'])
gtf.head()

# extract tRNA gene entries
tRNA_gtf = gtf[gtf['Feature']=='tRNA']


chrom = ['>2L', '>2R', '>3L', '>3R', '>X', '>mitochondrion_genome']  #d.melanogaster chromosomes
#reformat fasta to match chromosomes 
fasta = []
with open('dmel-all-chromosome-r6.27.fasta', 'r') as f:
    for line in f:
        if line.startswith('>'):
            line = line[:3].replace(" ","")
        fasta.append(line.rstrip())

##write to a new file
with open('mod_all_chr_r6.27.fasta','w') as f:
    for line in fasta:
        f.write(line+'\n')

## add extra 40 bases to flank tRNAs
tRNA_flank = tRNA_gtf.drop(['Database', '?', '?.1'], axis=1)
tRNA_flank['Start'] = tRNA_flank['Start']-40
tRNA_flank['End'] = tRNA_flank['End']+40


##build a dictionary {chrom: sequence}from fasta file for faster look-up
fas_dict = dict()
for c in chrom:
    start = fasta.index(c)+1 ##sequence is one psotion to righ tof chr in fasta list
    s= []
    #get all the sequence for the chromosome
    for line in fasta[start:]:
        if line.startswith('>'):
            break
        s.append(line)
    seq = "".join(s)
    fas_dict[c]=seq

fas_dict.keys()

#save as binary file
with open('chrom_fasta_dict_r6.27.pickle','wb') as f:
    pickle.dump(fas_dict,f)
	
	
## now extract tRNA + flanking sequences from fasta dict based on gtf coordinates

chrom = tRNA_flank['Chr'].unique()
chrom_groups = tRNA_flank.groupby('Chr')

df_list = []
for ch in chrom:
    df = chrom_groups.get_group(ch)
    f = []
    for t in range(len(df['Start'].values)):
        begin = df['Start'].values[t]
        end  = df['End'].values[t]
        ##look up the sequence
        flank_seq = fas_dict['>'+ch][begin:end+1] #python is zero-based
        header = '>'+df['Attributes'].values[t]+'\n'
        f.append([header+flank_seq])
    df['Fasta_Seq'] = f
    df_list.append(df)
tRNA_flank2 = pd.concat(df_list, axis=0)

##write tRNAs to  a fasta file
with open('pre-tRNAs_r6.27.fasta', 'w') as f:
    for seq in tRNA_flank2['Fasta_Seq'].values:
        f.write(seq[0]+'\n')


## Extrac tRNA genes and make custom mature sequences by adding CCA at 3' end 		
		
chrom = tRNA['Chr'].unique()
chrom_groups = tRNA.groupby('Chr')

df_list = []
for ch in chrom:
    df = chrom_groups.get_group(ch)
    f = []
    for t in range(len(df['Start'].values)):
        begin = df['Start'].values[t]
        end  = df['End'].values[t]
        ##look up the sequence
        seq = fas_dict['>'+ch][begin:end+1]+'CCA' #python is zero-based
        header = '>'+df['Attributes'].values[t]+'\n'
        f.append([header+seq])
    df['Fasta_Seq'] = f
    df_list.append(df)
tRNA_2 = pd.concat(df_list, axis=0)

##make a fasta file
with open('mature_tRNAs_r6.27.fasta', 'w') as f:
    for seq in tRNA_2['Fasta_Seq'].values:
        f.write(seq[0]+'\n')