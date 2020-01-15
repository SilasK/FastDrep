import os,sys
from glob import glob

genome_folder='genomes'

input_genome_folder=config['genome_folder']

#genomes= glob(os.path.join(genome_folder,"*"))

sys.path.append(os.path.join(os.path.dirname(workflow.snakefile),'scripts'))
from common import genome_pdist as gd
from common import io




include: "rules/checkm.smk"
include: "rules/filter.smk"


rule all:
    input:
        "tables/mummer_dist.tsv",
