import pyfastx
from multiprocessing import Pool
import pandas as pd
from numpy import log
import os
from .io import simplify_path


def genome_stats(fasta_file,remove_index=True):
    """Uses pyfastx to get genome stats from a fasta file. Outputs a tuple with:
       name,Length, n_seq,N50,L50
    """

    name = simplify_path(fasta_file)

    fa = pyfastx.Fasta(fasta_file)
    N50,L50 =fa.nl(50)
    n_seq= len(fa)
    Length= fa.size

    if remove_index:
        os.remove(fa.file_name+'.fxi')

    return name,Length, n_seq,N50,L50


def get_many_genome_stats(filenames,output_filename,threads=1):
    """Small function to calculate total genome length and N50
    """

    pool = Pool(threads)


    results= pool.map(genome_stats,filenames)
    Stats= pd.DataFrame(results,columns=["Genome","Length", "Nseqs","N50","L50"])
    Stats['logL50']=log(Stats.L50)
    Stats.to_csv(output_filename,sep='\t',index=False)
