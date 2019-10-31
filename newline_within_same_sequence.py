# An example when a newline separates nucleotides of the same sequence

fasta = []
file = []
with open('host_transc_introns.fasta', 'r') as f:
    for line in f:
        file.append(line.strip())

    for i in range(len(file)):
        line = file[i]
        if line.startswith('>'):
            header = (line.strip().split(">")[1])  # collect just the header without '>'
            seq = ""
            j = i  # the index for the sequence
            while True:
                j += 1  # the sequence is one position right of the header line
                if j == len(file):
                    break
                if file[j].startswith('>'):
                    break  # when it encounters another header
                else:
                    s = file[j]
                    seq += s
            fasta.append(header)
            fasta.append(seq)
        else:
            continue

# re-rewrite the file as a proper fasta such that only newline between sequences and
# headers
with open('host_transc_introns_V2.fasta', 'w') as f:
    for line in fasta:
        f.write("{}\n".format(line))
