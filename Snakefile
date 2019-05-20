import os,sys
from glob import glob



genome_folder=config['genome_folder']
#genomes= glob(os.path.join(genome_folder,"*"))

sys.path.append(os.path.join(os.path.dirname(workflow.snakefile),'scripts'))


#include: "rules/fastani.smk"
include: "rules/mash.smk"




rule all:
    input:
        "ANI.tsv"
