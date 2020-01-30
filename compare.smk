#this snakefile is for comparison of a given genome set


import os,sys
from glob import glob

genome_folder='genomes'

input_genome_folder=config['genome_folder']

#genomes= glob(os.path.join(genome_folder,"*"))

sys.path.append(os.path.join(os.path.dirname(workflow.snakefile),'scripts'))
from common import genome_pdist as gd
from common import io



include: "rules/fastani.smk"
include: "rules/pyani.smk"
include: "rules/bbsketch.smk"

include: "rules/sketch.smk"
include: "rules/alignments.smk"
include: "rules/cluster.smk"

rule all:
    input:
        #"tables/fastANI_dists.tsv",
        #expand("pyani/{method}",method=['ANIm','ANIb']),
        f"tables/{config['aligner']}_dists.tsv",
        "representatives/strains",
        "representatives/species"


rule all_bbsketch:
    input:
        "tables/bbsketch_aa.tsv",
        "tables/bbsketch_nt.tsv"

rule all_species:
    input:
        "tables/refseq_mapping_species.tsv",
        "representatives/species",
