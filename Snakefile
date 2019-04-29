import os
from glob import glob



genome_folder=config['genome_folder']
genomes= glob(os.path.join(genome_folder,"*"))



include: "rules/mash.smk"
#include: "rules/fastani.smk"

rule all:
    input:
        #"ANIs.txt",
        "mash_dists.txt"
