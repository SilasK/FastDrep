import sys, os
import pandas as pd

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = open(snakemake.log[0], "a")

from snakemake.shell import shell
from common.alignments import parse_paf_files


def run_minimap(
    ref,
    query,
    paf,
    options=snakemake.params.minimap_extra,
    threads=snakemake.threads,
    log=snakemake.log[0],
):

    os.makedirs(os.path.dirname(paf), exist_ok=True)

    shell(
        "minimap2 {options} -t {threads} {ref} {query}  > {paf}.tmp 2> {log} ;"
        "mv {paf}.tmp {paf} 2> {log}"
    )


def many_minimap( alignment_list, genome_folder, paf_folder, genome_stats_file, extension ):

    stats = pd.read_csv(genome_stats_file, index_col=0, sep="\t")

    os.makedirs(paf_folder, exist_ok=True)

    paf_files = []

    with open(alignment_list) as f:
        for line in f:

            # get larger genome as ref and smaler as query
            genome_pair = line.strip().split()
            genome_query, genome_ref = (
                stats.loc[genome_pair, "Length"].sort_values().index
            )

            paf_file = os.path.join(paf_folder, genome_ref, genome_query + ".paf")

            if not os.path.exists(paf_file):

                # run minimap
                fasta_ref = os.path.join(genome_folder, genome_ref + extension)
                fasta_query = os.path.join(genome_folder, genome_query + extension)
                run_minimap(fasta_ref, fasta_query, paf_file)

            paf_files.append(paf_file)

    return paf_files


paf_files = many_minimap(
    snakemake.input.alignment_list,
    snakemake.input.genome_folder,
    snakemake.params.paf_folder,
    snakemake.input.genome_stats,
    snakemake.params.extension,
)


parse_paf_files(
    paf_files, snakemake.input.genome_stats, snakemake.output.alignments_stats
)
