import numpy as np 
import pandas as pd

##Author: Rhondene Wint
## An example of how to obtain transcript ID of the longest CDS of a gene in order download multi-fasta from flybase


##load gtf file
fb_gtf = pd.read_table('dmel-all-r6.27.gtf',sep='\t',header=None, names=['Chr','Database','Feature','Start','End','?', 'Strand','?.1', 'Attributes'])


##parse Attributes column  to obtain gene IDs
gene_id = []
for attr in fb_gtf['Attributes'].values:
    ID = attr.split(";")[0].split(" ")[1].replace('"',"")
    gene_id.append(ID)
fb_gtf['Gene ID'] =  gene_id
fb_gtf.head()

## add  a length column to the gtf file 
if np.all(fb_gtf['Start']<= fb_gtf['End']):
	fb_gtf['Length']=fb_gtf['End'].values-fb_gtf['Start'].values
else:
	fb_gtf['Length']= np.abs(fb_gtf['End'].values-fb_gtf['Start'].values)
	
## to get transcript ID of protein-coding genes, extract all  mRNA features
all_mRNAs = fb_gtf[fb_gtf['Feature']=='mRNA']

trans_id = []
for attr in all_mRNAs['Attributes'].values:
    ID = attr.split(";")[2].split(" ")[2].replace('"',"")
    trans_id.append(ID)
all_mRNAs['Transcript ID'] =  trans_id

##save new gtf 
all_mRNAs.to_csv('all_mRNAs_fb_r6.27.gtf',sep='\t',index=False)

gene_ids = all_mRNAs['Gene ID'].unique()
genes = all_mRNAs.groupby('Gene ID')

long_list_2 = []
for gene in gene_ids:
    df = genes.get_group(gene)
    max_ = df['Length'].max()
    df[df['Length']==max_]  ##look up transcript that matches length
    long_list_2.append( df[df['Length']==max_].head(1))  ##take the first isoform, if equal lenght

longest_transc = pd.concat(long_list_2,axis=0)
##shape should be equal total number number of protein-coding genes
assert longest_transc.shape[0] == len(gene_ids)

## save transcript IDs to upload directly to Flybase 

longest_trans['Transcript ID'].to_csv('longest_transcript_IDs.txt',sep='\t',index=False, header=None)



	
					