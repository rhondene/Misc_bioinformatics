##An example when a newline separates nucleotides of within the same sequence,
#for example Flybase sometimes writes newline within the same sequence when downloading bulk sequences from gene ID
def fix_fasta(filename):
	"filename: name of fasta formatted file" 

    fasta = []
    file=[]
    with open('../Transcripts/{}'.format(filename), 'r') as f:
        for line in f:
            file.append(line.strip())

        for i in range(len(file)):
            line = file[i]
            if line.startswith('>') is True:
                header =(line.rstrip())
                seq = ""
                j=i  #the index for the sequence
                while True:
                    j+=1  #the sequence is one position right of the header line
                    if j==len(file):
                        break
                    if file[j].startswith('>') is True:
                        break  #when it encounters another header
                    else:
                        s = file[j]
                        seq+=s
                fasta.append(header)
                fasta.append(seq)
            else:
                continue

    ### re-rewrite the file as a proper fasta such that newline only between sequences and headers
    with open('../Transcripts_fixed/file_name', 'w') as f:
        for line in fasta:
            f.write(line+'\n')
