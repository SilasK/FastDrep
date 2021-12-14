#this snakefile is for comparison of a given genome set


import os,sys
from glob import glob

genome_folder='genomes'

input_genome_folder=config['genome_folder']

#genomes= glob(os.path.join(genome_folder,"*"))

sys.path.append(os.path.join(os.path.dirname(workflow.snakefile),'scripts'))
from common import genome_pdist as gd
from common import io






include: "rules/sketch.smk"
include: "rules/alignments.smk"
include: "rules/cluster.smk"

rule all:
    input:
        #"tables/fastANI_dists.tsv",
        #"representatives/strains",
        "representatives/species"




for r in workflow.rules:
    if not "mem" in r.resources:
        r.resources["mem_mb"]=config["mem"]['default'] *1000
    if not "time" in r.resources:
        r.resources["time_min"]=config["runtime"]["default"] *60

#
