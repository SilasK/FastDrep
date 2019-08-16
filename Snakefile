import os,sys
from glob import glob



input_genome_folder=config['genome_folder']

#genomes= glob(os.path.join(genome_folder,"*"))

sys.path.append(os.path.join(os.path.dirname(workflow.snakefile),'scripts'))
from common import genome_pdist as gd





include: "rules/filter.smk"
include: "rules/fastani.smk"
include: "rules/mash.smk"
include: "rules/pyani.smk"
include: "rules/strains.smk"


rule all:
    input:
        # "ANI.tsv",
        #expand("pyani/{method}",method=['ANIm','ANIb']),
        "mash_dists.txt",
        "tables/genome_stats.tsv"
