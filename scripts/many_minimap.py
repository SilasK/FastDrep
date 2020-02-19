import sys,os


sys.stdout = open(snakemake.log[0], 'w')

from snakemake.shell import shell
from common.alignments import parse_paf_files


def run_minimap(ref,query,paf,options=snakemake.params.minimap_extra,
                threads=snakemake.threads,log=snakemake.log[0]):
    shell("minimap2 {options} -t {threads} {ref} {query}   > {paf} 2> {log}")

def many_minimap(alignment_list,genome_folder,paf_folder,extension='.fasta'):

    os.makedirs(paf_folder,exist_ok=True)

    paf_files=[]

    with open(alignment_list) as f:
        for line in f:
            genome1,genome2= line.strip().split()

            paf= os.path.join(paf_folder,f"{genome1}/{genome2}.paf")

            if not os.path.exists(paf):
                os.makedirs(os.path.join(paf_folder,genome1),exist_ok=True)
                ref= os.path.join(genome_folder,genome1+extension)
                query =os.path.join(genome_folder,genome2+extension)
                run_minimap(ref,query,paf)
            paf_files.append(paf)

    return paf_files



paf_files = many_minimap(snakemake.input.alignment_list,
                       snakemake.input.genome_folder,
                       snakemake.params.paf_folder,
                       snakemake.params.extension
                       )




parse_paf_files(paf_files,snakemake.input.genome_stats,snakemake.output.alignments_stats)
