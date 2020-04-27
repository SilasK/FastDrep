from multiprocessing import Pool
import pandas as pd
import os, sys
from .io import simplify_path, simply_open
from itertools import groupby
import numpy as np
import gzip as gz


def get_stats_from_lengths(lengths):

    sorted_lengths = sorted(lengths, reverse=True)
    csum = np.cumsum(sorted_lengths)

    Total_length = int(sum(lengths))
    N = len(lengths)

    n2 = int(Total_length / 2)

    # get index for cumsum >= N/2
    csumn2 = min(csum[csum >= n2])
    ind = int(np.where(csum == csumn2)[0][0])

    N50 = sorted_lengths[ind]

    return Total_length, N, N50


def genome_stats(fasta_file):
    """Get genome stats from a fasta file. Outputs a tuple with:
       name,Length, n_seq,N50
    """

    try:

        name = simplify_path(fasta_file)

        scaffold_lengths = []
        contig_lengths = []

        with simply_open(fasta_file, "r") as fasta:
            ## parse each sequence by header: groupby(data, key)
            faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))

            for record in faiter:
                ## join sequence lines
                sequence = sum(s.strip() for s in faiter.__next__())
                scaffold_lengths.append(len(sequence))
                contig_lengths += [
                    len(contig) for contig in sequence.replace("N", " ").split()
                ]

        Length_scaffolds, N_scaffolds, N50 = get_stats_from_lengths(scaffold_lengths)

        Length_contigs, N_contigs, _ = get_stats_from_lengths(contig_lengths)

    except Exception as e:
        raise Exception(f"Error in calculating stats of {fasta_file}") from e

    return name, Length_scaffolds, N_scaffolds, N50, Length_contigs, N_contigs


def get_many_genome_stats(filenames, output_filename, threads=1):
    """Small function to calculate total genome length and N50
    """

    pool = Pool(threads)

    results = pool.map(genome_stats, filenames)
    Stats = pd.DataFrame(
        results,
        columns=[
            "Genome",
            "Length",
            "N_scaffolds",
            "N50",
            "Length_contigs",
            "N_contigs",
        ],
    )
    Stats.to_csv(output_filename, sep="\t", index=False)
