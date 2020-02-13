import pyfastx
from multiprocessing import Pool
import pandas as pd
import os,sys
from .io import simplify_path
from itertools import groupby
import numpy as np


def genome_stats(fasta_file,remove_index=True):
    """Uses pyfastx to get genome stats from a fasta file. Outputs a tuple with:
       name,Length, n_seq,N50
    """

    name = simplify_path(fasta_file)



    # from Bio import SeqIO

    lengths = []

    with open(fasta_file) as fasta:
        ## parse each sequence by header: groupby(data, key)
        faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))

        for record in faiter:
            ## join sequence lines
            seq = "".join(s.strip() for s in faiter.__next__())
            lengths.append(len(seq))

    ## sort contigs longest>shortest
    all_len=sorted(lengths, reverse=True)
    csum=np.cumsum(all_len)

    Length = int(sum(lengths))
    n_seq = len(lengths)

    n2=int(Length/2)

    # get index for cumsum >= N/2
    csumn2=min(csum[csum >= n2])
    ind=int(np.where(csum == csumn2)[0][0])

    N50 = all_len[ind]
    ## N90
    #nx90=int(sum(lengths)*0.90)

    ## index for csumsum >= 0.9*N
    #csumn90=min(csum[csum >= nx90])
    #ind90=numpy.where(csum == csumn90)
    #n90 = all_len[int(ind90[0])]



    return name,Length, n_seq,N50


def get_many_genome_stats(filenames,output_filename,threads=1):
    """Small function to calculate total genome length and N50
    """

    pool = Pool(threads)


    results= pool.map(genome_stats,filenames)
    Stats= pd.DataFrame(results,columns=["Genome","Length", "Nseqs","N50"])
    Stats.to_csv(output_filename,sep='\t',index=False)
