import numpy as np


Nucleotides_with_N= {'A', 'C', 'G', 'N', 'T','a','c','g','n','t'}
def check_accepted_characters(seq,acepted_characters=Nucleotides_with_N):
    if not (set(seq) - acepted_characters == set()):
        raise Exception("I found something else than the accepted characters"
                        " in this line in this fasta file. "
                        f"Line:\n{line}\naccepted_cahracters: {acepted_characters}"
                       )

def get_genome_stats(fasta_file):
    """Small function to calculate total genome length and N50
    """

    Lengths=[]

    with open(fasta_file) as fasta:
        seq=''

        for line in fasta:
            if line[0]=='>':
                Lengths.append(0)
            else:
                seq= line.strip()
                check_accepted_characters(seq)
                Lengths[-1]+= len(seq)


    Lengths= np.array(Lengths)


    Total_length= Lengths.sum()

    Lengths.sort()
    Lengths= Lengths[-1:0:-1] # sort length to shortest

    N50= Lengths[Lengths.cumsum() >= Total_length/2][0] # get first contig with cumsum larger than 50% of Total length

    return Total_length, N50
