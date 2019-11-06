"""Author: Rhondene Wint
Purpose: Using intron-only bed file and the genome gtf file to extract the transcript witht the longest 
intron for each gene in D.melanogaster in order to compute total intronic composition in D.melanogaster
You can also use this output to download these intronic sequences from flybase or ensembl"""


from interval import interval, inf, imath  ## interval requires running python on  Linux or Mac
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn 

## to get total introns, I will  merge all overlapping introns in the bed file of a transcript
## then map the transcript ID to the gene IDs


all_introns = pd.read_table('dmel_introns_ensembl.bed', sep='\t',header=None, names=['Chr','Start','End','ID', 'Phase', 'Strand'])
##split columns to get transcript IDs 
IDs = []
for name in all_introns['ID'].values:
    n = name.split('_')[0]
    IDs.append(n)
all_introns['Transcript_ID']=IDs

# compute intron lengths
all_introns['Length']= all_introns['End']- all_introns['Start']
all_introns.head()

transc_ids = all_introns['Transcript_ID'].unique()
transcripts = all_introns.groupby('Transcript_ID')

"""merge overlapping introons for each gene and compute length"""
merged_regions = dict()
## iterate over each gene
for ID in transc_ids:
    df = transcripts.get_group(ID).sort_values(by='Start', ascending=True) ##sort by ascending order
    #set the initial value to start of the earliest intron
    consol_region = interval[df['Start'].values[0]] 
    ##iterate over exons of the gene
    for i in range(df.shape[0]):
        #create an interval of an individual exon region
        intron_size = interval[df['Start'].values[i],df['End'].values[i]] 
        ##consolidate overlapping the intron region
        consol_region= consol_region | intron_size  
    ##finally store  a list of non-overlapping intronic intervals of a gene
    merged_regions[ID]=consol_region

##store total_intron_size for each transcript
transc_total_introns=dict()
for ID in transc_ids:
    consol_region = merged_regions[ID]
    total=0
    for region in consol_region:
        total+= region[1]-region[0]
    
    transc_total_introns[ID] = total  ##for my own use when I want to look intron sise distbtion
	
##store the total intron lenghths for transcript for each gene in a table
total_introns_transc = pd.DataFrame.from_dict(transc_total_introns,orient='index').reset_index()
total_introns_transc.columns=['Transcript ID', 'Total Intron Size']


##select transcript entries for mRNA and ncRNA fromt the gtf file

all_transc= fb_gtf.query('Feature!="gene" and Feature!="exon" and Feature!= "5UTR" and Feature!="3UTR" and Feature!="stop_codon" and Feature!="start_codon" and Feature!="CDS"')

##parse Attributes  to obtain transcript IDs
trans_id = []
for attr in all_transc['Attributes'].values:
    ID = attr.split(";")[2].split(" ")[2].replace('"',"")
    trans_id.append(ID)
all_transc['Transcript ID'] =  trans_id

#### Okay so 83 transcript ID in the ensembl intron annotation based on  
##Flybase release 6.22 is missing in the flybase release 6.27 gtf, so those 83 gonna get dropped

""" identify missing transcript entries """
present = []
for i in range(total_introns_transc.shape[0]):
    transc_id = total_introns_transc['Transcript ID'].values[i]
    if transc_id not in all_transc['Transcript ID'].values:
        present.append('Missing')
    else:
        present.append('Yes')
total_introns_transc['Present']=present

""" select entries that are present in both annotation, i.e. filter out the 83 transcitps"""
total_introns2 = total_introns_transc.query('Present=="Yes"')

"""update the code for mapping transcript ID to gen ID"""
parent_gene = []
for ID in total_introns2['Transcript ID'].values:
    for gene in list(gene_dict.keys()):
        if ID in gene_dict[gene]:
            parent_gene.append(gene)
            break
        else:
            continue


""" compute total intron size"""
parent_genes = total_introns2['Gene ID'].unique()
total_intron_size = 0
genes =  total_introns2.groupby('Gene ID')

for gene in parent_genes:
    ##get all transcripts for each gene
    df = genes.get_group(gene)
    ##update the total genomic intron with the max intron length
    total_intron_size+=df['Total Intron Size'].max()
total_intron_size


